////////////////////////////////////////////////////////////////////////
// File:        getophitMCPE_module.cc
// Purpose:     Save per-OpHit info: MC-candidate flag, channel, hit PE,
//              and total PE of the parent flash.
//              (Optional) estimate PhotonBackTracker Delay by scanning.
// Author:      Gabriel
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "larcore/Geometry/WireReadout.h"

#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"

#include "TTree.h"

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {
  struct HitForDelay 
  {
    int opdet = -1;
    double t_hit_ns = 0.0;  // PeakTime * 1000
    double w_ns     = 0.0;  // Width    * 1000
  };

  inline bool HasKeyInRange(std::vector<double> const& v, double a, double b)
  {
    if (v.empty()) return false;
    if (a > b) std::swap(a, b);
    auto it = std::lower_bound(v.begin(), v.end(), a);
    return (it != v.end() && *it <= b);
  }
}

class getophitMCPE : public art::EDAnalyzer
{
public:
  explicit getophitMCPE(fhicl::ParameterSet const& p);

  void beginJob() override;
  void analyze(art::Event const& e) override;
  void endJob() override;

private:
  art::InputTag fFlashLabel;

  // TTree
  TTree* fTree = nullptr;
  int   fMonteCarlo = 0;   // 0/1
  int   fChannel    = -1;  // OpChannel
  float fHitPE      = 0.f; // PE of OpHit
  float fFlashPE    = 0.f; // TotalPE of parent flash

  // Optional delay estimation (scan)
  bool   fEstimateDelay = true;
  double fDelayMinNs = -1.0e6;   // ns
  double fDelayMaxNs = +1.0e6;   // ns
  double fDelayStepNs = 5000.0;  // ns
  float  fDelayPEMin = 20.0f;    // only use hits with PE > this in delay scan
  unsigned fDelayMaxHits = 2000; // cap per event

  std::vector<double> fDelayGrid; // ns
  std::vector<long long> fDelayScoreGlobal;
};

getophitMCPE::getophitMCPE(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
  fFlashLabel = p.get<art::InputTag>("FlashLabel");

  // optional delay scan knobs
  fEstimateDelay = p.get<bool>("EstimateDelay", true);
  fDelayMinNs    = p.get<double>("DelayScanMinNs",  -1.0e6);
  fDelayMaxNs    = p.get<double>("DelayScanMaxNs",  +1.0e6);
  fDelayStepNs   = p.get<double>("DelayScanStepNs",  5000.0);
  fDelayPEMin    = p.get<float>("DelayPEMin", 20.0f);
  fDelayMaxHits  = p.get<unsigned>("DelayMaxHits", 2000);

  if (fEstimateDelay) {
    if (fDelayStepNs <= 0) throw cet::exception("getophitMCPE") << "DelayScanStepNs must be > 0\n";
    if (fDelayMinNs > fDelayMaxNs) std::swap(fDelayMinNs, fDelayMaxNs);

    for (double d = fDelayMinNs; d <= fDelayMaxNs + 0.5*fDelayStepNs; d += fDelayStepNs)
      fDelayGrid.push_back(d);

    fDelayScoreGlobal.assign(fDelayGrid.size(), 0);
  }
}

void getophitMCPE::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("hit_tree", "Per-OpHit MC candidate info");
  fTree->Branch("montecarlo",   &fMonteCarlo, "montecarlo/I");
  fTree->Branch("channel",      &fChannel,    "channel/I");
  fTree->Branch("hitPE",        &fHitPE,      "hitPE/F");
  fTree->Branch("flashTotalPE", &fFlashPE,    "flashTotalPE/F");
}

void getophitMCPE::analyze(art::Event const& e)
{
  auto& pbts = *art::ServiceHandle<cheat::PhotonBackTrackerService>();
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

  // ---- flashes ----
  auto flash_h = e.getHandle<std::vector<recob::OpFlash>>(fFlashLabel);
  if (!flash_h) {
    mf::LogWarning("getophitMCPE") << "Cannot load OpFlash: " << fFlashLabel;
    return;
  }

  art::FindManyP<recob::OpHit> fmOpHits(flash_h, e, fFlashLabel);
  if (!fmOpHits.isValid()) {
    mf::LogWarning("getophitMCPE") << "No OpFlash<->OpHit assns for " << fFlashLabel;
    return;
  }

  // If delay scan is enabled, build opdet -> sorted key-times (timePDclockSDPsMap keys)
  std::unordered_map<int, std::vector<double>> opdetTimes;
  std::vector<HitForDelay> hitsForDelay;

  if (fEstimateDelay) {
    auto const& btrs = pbts.OpDetBTRs();
    opdetTimes.reserve(btrs.size());

    for (auto const& btr : btrs) {
      if (!btr) continue;
      int opdet = btr->OpDetNum();
      auto const& m = btr->timePDclockSDPsMap();
      auto& v = opdetTimes[opdet];
      v.reserve(m.size());
      for (auto const& kv : m) v.push_back(kv.first);
      std::sort(v.begin(), v.end());
    }

    hitsForDelay.reserve(std::min<unsigned>(fDelayMaxHits, 2000u));
  }

  std::vector<art::Ptr<recob::OpFlash>> flashes;
  art::fill_ptr_vector(flashes, flash_h);

  for (auto const& fl : flashes)
  {
    const float flashTotalPE = (float)fl->TotalPE();
    auto ophits = fmOpHits.at(fl.key());

    for (auto const& oph : ophits)
    {
      // --- MC flag: require trackID != 0 (SDPs OR TrackIds)
      bool hasMC = false;

      auto sdps = pbts.OpHitToSimSDPs_Ps(oph);
      for (auto const* sdp : sdps) {
        if (!sdp) continue;
        if (std::abs(sdp->trackID) != 0) { hasMC = true; break; }
      }

      if (!hasMC) {
        auto tids = pbts.OpHitToTrackIds(oph);
        for (int tid : tids) {
          if (std::abs(tid) != 0) { hasMC = true; break; }
        }
      }

      // --- fill tree
      fMonteCarlo = hasMC ? 1 : 0;
      fChannel    = (int)oph->OpChannel();
      fHitPE      = (float)oph->PE();
      fFlashPE    = flashTotalPE;
      fTree->Fill();

      // --- collect hits for delay scan (use ONLY "good-ish" hits)
      if (fEstimateDelay && hasMC && oph->PE() > fDelayPEMin && hitsForDelay.size() < fDelayMaxHits) {
        HitForDelay h;
        h.opdet     = wireReadout.OpDetFromOpChannel(oph->OpChannel());
        h.t_hit_ns  = oph->PeakTime() * 1000.0;
        h.w_ns      = oph->Width()    * 1000.0;
        hitsForDelay.push_back(h);
      }
    }
  }

  // --- delay scan scoring (per event -> add to global)
  if (fEstimateDelay && !hitsForDelay.empty() && !opdetTimes.empty())
  {
    for (size_t i = 0; i < fDelayGrid.size(); ++i) {
      double delay = fDelayGrid[i];
      long long sc = 0;

      for (auto const& h : hitsForDelay) {
        auto it = opdetTimes.find(h.opdet);
        if (it == opdetTimes.end() || it->second.empty()) continue;

        // service uses: start=(pTime-w)*1000 - Delay  =>  start=t_hit_ns - w_ns - Delay
        double start = h.t_hit_ns - h.w_ns - delay;
        double end   = h.t_hit_ns + h.w_ns - delay;

        if (HasKeyInRange(it->second, start, end)) sc++;
      }

      fDelayScoreGlobal[i] += sc;
    }
  }
}

void getophitMCPE::endJob()
{
  if (!fEstimateDelay || fDelayGrid.empty()) return;

  // best delay = max score
  size_t ibest = std::distance(
    fDelayScoreGlobal.begin(),
    std::max_element(fDelayScoreGlobal.begin(), fDelayScoreGlobal.end())
  );

  std::cout << "\n[getophitMCPE] Delay scan summary (global over job)\n"
            << "  best Delay = " << fDelayGrid[ibest] << " ns"
            << "  (score=" << fDelayScoreGlobal[ibest] << ")\n";

  // print a few neighbors for sanity
  auto printOne = [&](size_t i){
    std::cout << "    Delay " << fDelayGrid[i] << " ns -> score " << fDelayScoreGlobal[i] << "\n";
  };

  std::cout << "  around best:\n";
  if (ibest > 0) printOne(ibest-1);
  printOne(ibest);
  if (ibest + 1 < fDelayGrid.size()) printOne(ibest+1);
}

DEFINE_ART_MODULE(getophitMCPE)
