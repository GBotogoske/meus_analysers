////////////////////////////////////////////////////////////////////////
// File: CompareSimPhotonsOpHits_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "TTree.h"

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <utility>

class CompareSimPhotonsOpHits : public art::EDAnalyzer
{
public:
  explicit CompareSimPhotonsOpHits(fhicl::ParameterSet const& p);

  void beginJob() override;
  void analyze(art::Event const& evt) override;

private:
  art::InputTag fOpHitTag;
  art::InputTag fSimPhotonsTag;
  bool fUseStartTime;
  double fTimeOffset;

  TTree* fTree = nullptr;

  int    b_run = 0;
  int    b_subrun = 0;
  int    b_event = 0;
  int    b_opchannel = -1;

  int    b_hit_index = -1;
  double b_hit_peak = -9999.0;
  double b_hit_start = -9999.0;
  double b_hit_width = 0.0;
  double b_hit_pe = 0.0;

  double    b_n_simphotons = 0.0;
  double b_t_first_sim = -9999.0;
  double b_t_last_sim = -9999.0;

  int    b_has_hit = 0;
  int    b_has_sim = 0;

   std::vector<double> fPDEVector;
   double Xtalk,Kdup;
};

CompareSimPhotonsOpHits::CompareSimPhotonsOpHits(fhicl::ParameterSet const& p)
  : art::EDAnalyzer(p)
  , fOpHitTag(p.get<art::InputTag>("OpHitTag",
                                   art::InputTag("ophitspe", "", "myFlash")))
  , fSimPhotonsTag(p.get<art::InputTag>("SimPhotonsTag",
                                        art::InputTag("PDFastSim", "", "G4")))
  , fUseStartTime(p.get<bool>("UseStartTime", false))
  , fTimeOffset(p.get<double>("TimeOffset", 250.0))
{
  fPDEVector = p.get<std::vector<double>>("PDEvector",std::vector<double>(180,0.03));
  Xtalk = p.get<double>("XTalk",0.01);
  Kdup = Xtalk/(1-Xtalk);
}

void CompareSimPhotonsOpHits::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "OpHit <-> SimPhotonsLite association");

  fTree->Branch("run",       &b_run,       "run/I");
  fTree->Branch("subrun",    &b_subrun,    "subrun/I");
  fTree->Branch("event",     &b_event,     "event/I");
  fTree->Branch("opchannel", &b_opchannel, "opchannel/I");

  fTree->Branch("hit_index", &b_hit_index, "hit_index/I");
  fTree->Branch("hit_peak",  &b_hit_peak,  "hit_peak/D");
  fTree->Branch("hit_start", &b_hit_start, "hit_start/D");
  fTree->Branch("hit_width", &b_hit_width, "hit_width/D");
  fTree->Branch("hit_pe",    &b_hit_pe,    "hit_pe/D");

  fTree->Branch("n_simphotons", &b_n_simphotons, "n_simphotons/D");
  fTree->Branch("t_first_sim",  &b_t_first_sim,  "t_first_sim/D");
  fTree->Branch("t_last_sim",   &b_t_last_sim,   "t_last_sim/D");

  fTree->Branch("has_hit", &b_has_hit, "has_hit/I");
  fTree->Branch("has_sim", &b_has_sim, "has_sim/I");
}

void CompareSimPhotonsOpHits::analyze(art::Event const& evt)
{
  b_run    = evt.run();
  b_subrun = evt.subRun();
  b_event  = evt.id().event();

  auto simHandle = evt.getValidHandle<std::vector<sim::SimPhotonsLite>>(fSimPhotonsTag);
  auto hitHandle = evt.getValidHandle<std::vector<recob::OpHit>>(fOpHitTag);

  // Por canal: vetor de (tempo, nPhotons)
  std::unordered_map<int, std::vector<std::pair<int,double>>> photonsByCh;

  for (auto const& simphot : *simHandle) {
    int ch = simphot.OpChannel;
    auto& v = photonsByCh[ch];

    for (auto const& kv : simphot.DetectedPhotons) {
    
      int time  = kv.first/1000;
      double nphot = kv.second*fPDEVector[ch]*(1+Kdup);
      v.emplace_back(time, nphot);
    }
  }

  for (auto& kv : photonsByCh) {
    auto& v = kv.second;
    std::sort(v.begin(), v.end(),
              [](auto const& a, auto const& b) { return a.first < b.first; });
  }

  // Hits por canal
  std::unordered_map<int, std::vector<std::pair<int, const recob::OpHit*>>> hitsByCh;
  for (size_t i = 0; i < hitHandle->size(); ++i) {
    auto const& hit = hitHandle->at(i);
    hitsByCh[hit.OpChannel()].push_back({static_cast<int>(i), &hit});
  }

  // União dos canais
  std::unordered_map<int, bool> allChannels;
  for (auto const& kv : photonsByCh) allChannels[kv.first] = true;
  for (auto const& kv : hitsByCh)    allChannels[kv.first] = true;

  for (auto const& kv : allChannels) {
    int ch = kv.first;

    auto& photons = photonsByCh[ch]; // vetor de (tempo, multiplicidade)
    auto& hits    = hitsByCh[ch];

    // Para cada hit, guardar quais bins de tempo foram associados
    std::vector<std::vector<int>> hitToPhotonBins(hits.size());
    std::vector<int> photonBinAssigned(photons.size(), -1);

    // Associa cada bin de tempo ao melhor hit
    for (size_t ip = 0; ip < photons.size(); ++ip) {
      double t = static_cast<double>(photons[ip].first);

      int bestHit = -1;
      double bestDist = std::numeric_limits<double>::max();

      for (size_t ih = 0; ih < hits.size(); ++ih) {
        auto const* hit = hits[ih].second;

        double t0 = fUseStartTime
          ? (hit->StartTime() + fTimeOffset)
          : (hit->PeakTime()  + fTimeOffset - 0.5 * hit->Width());

        double t1 = fUseStartTime
          ? (hit->StartTime() + fTimeOffset + hit->Width())
          : (hit->PeakTime()  + fTimeOffset + 0.5 * hit->Width());

        if (t >= t0 && t <= t1) {
          double dist = std::abs(t - (hit->PeakTime() + fTimeOffset));
          if (dist < bestDist) {
            bestDist = dist;
            bestHit = static_cast<int>(ih);
          }
        }
      }

      if (bestHit >= 0) {
        photonBinAssigned[ip] = bestHit;
        hitToPhotonBins[bestHit].push_back(static_cast<int>(ip));
      }
    }

    // 1) Todos os hits: com ou sem sim
    for (size_t ih = 0; ih < hits.size(); ++ih) {
      auto [hitIndex, hit] = hits[ih];

      b_opchannel = ch;
      b_hit_index = hitIndex;
      b_hit_peak  = hit->PeakTime() + fTimeOffset;
      b_hit_start = hit->StartTime() + fTimeOffset;
      b_hit_width = hit->Width();
      b_hit_pe    = hit->PE();

      b_has_hit = 1;
      b_has_sim = hitToPhotonBins[ih].empty() ? 0 : 1;

      b_n_simphotons = 0;
      if (b_has_sim) {
        for (int idx : hitToPhotonBins[ih]) {
          b_n_simphotons += photons[idx].second;
        }
        b_t_first_sim = static_cast<double>(photons[hitToPhotonBins[ih].front()].first);
        b_t_last_sim  = static_cast<double>(photons[hitToPhotonBins[ih].back()].first);
      }
      else {
        b_t_first_sim = -9999.0;
        b_t_last_sim  = -9999.0;
      }

      fTree->Fill();
    }

    // 2) Bins de sim sem hit -> PE = 0
    for (size_t ip = 0; ip < photons.size(); ++ip) 
    {
      if (photonBinAssigned[ip] >= 0) continue;

      b_opchannel = ch;
      b_hit_index = -1;
      b_hit_peak  = -9999.0;
      b_hit_start = -9999.0;
      b_hit_width = 0.0;
      b_hit_pe    = 0.0;

      b_has_hit = 0;
      b_has_sim = 1;

      b_n_simphotons = photons[ip].second;
      b_t_first_sim  = static_cast<double>(photons[ip].first);
      b_t_last_sim   = static_cast<double>(photons[ip].first);

      fTree->Fill();
    }
  }
}

DEFINE_ART_MODULE(CompareSimPhotonsOpHits)