#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TTree.h"

#include <iostream>
#include <unordered_set>
#include <vector>
#include <string>

class GetMyWireData : public art::EDAnalyzer
{
public:
  explicit GetMyWireData(fhicl::ParameterSet const& p);

  void analyze(art::Event const& e) override;
  void beginJob() override;

private:
  // --- Config ---
  art::InputTag fRawDigitLabel;

  art::InputTag fTrackLabel;
  art::InputTag fShowerLabel;

  art::InputTag fSliceLabel;
  art::InputTag fPFPLabel;

  // Assns (coloque explicitamente no fhicl!)
  art::InputTag fTrackHitAssnLabel;
  art::InputTag fShowerHitAssnLabel;

  art::InputTag fSliceToPFPAssnLabel;
  art::InputTag fPFPToTrackAssnLabel;
  art::InputTag fPFPToShowerAssnLabel;

  bool getShowers = false;
  std::string ClusterType; // "Track" ou "Slice"

  // --- Trees/branches ---
  int fRun = 0, fEvent = 0;

  TTree* fTree  = nullptr; // wire
  TTree* fTreeT = nullptr; // object (track/shower/slice) -> hits

  int fCh = -1;
  int fTrackID = -1;
  double fTime = 0.0;

  std::vector<int> fTPC, fWire, fPlane;
  std::vector<short> fadc;

  double fLength = 0.0;
  std::vector<float> fPeakTime;
  std::vector<float> fAmplitude;
  std::vector<float> fIntegral;
  std::vector<int>   fChHit;
  std::vector<int>   fTPCHit, fWireHit, fPlaneHit;
  int fType = -1; // 0=track, 1=shower, 3=slice

private:
  void ResetHitVectors()
  {
    fPeakTime.clear();
    fAmplitude.clear();
    fIntegral.clear();
    fChHit.clear();
    fTPCHit.clear();
    fWireHit.clear();
    fPlaneHit.clear();
  }

  void FillHitVectors(std::vector<art::Ptr<recob::Hit>> const& hits)
  {
    for (auto const& h : hits)
    {
      if (h.isNull()) continue;
      fPeakTime.push_back(h->PeakTime());
      fAmplitude.push_back(h->PeakAmplitude());
      fIntegral.push_back(h->Integral());
      fChHit.push_back(h->Channel());

      auto wid = h->WireID();
      fTPCHit.push_back(wid.TPC);
      fPlaneHit.push_back(wid.Plane);
      fWireHit.push_back(wid.Wire);
    }
  }
};

GetMyWireData::GetMyWireData(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
  // modo
  ClusterType = p.get<std::string>("ClusterType", "Track");
  getShowers  = p.get<bool>("getShowers", false);

  // labels
  fRawDigitLabel = p.get<art::InputTag>("RawDigitLabel",
                                       art::InputTag("tpcrawdecoder","daq","Detsim"));

  fTrackLabel  = p.get<art::InputTag>("TrackLabel",  art::InputTag("pandoraTrack"));
  fShowerLabel = p.get<art::InputTag>("ShowerLabel", art::InputTag("pandoraShower"));

  fSliceLabel  = p.get<art::InputTag>("SliceLabel",  art::InputTag("pandoraSlice"));
  fPFPLabel    = p.get<art::InputTag>("PFParticleLabel", art::InputTag("pandora"));

  // assns (defaults “seguros” = o próprio produtor do objeto)
  fTrackHitAssnLabel  = p.get<art::InputTag>("TrackHitAssnLabel",  fTrackLabel);
  fShowerHitAssnLabel = p.get<art::InputTag>("ShowerHitAssnLabel", fShowerLabel);

  fSliceToPFPAssnLabel   = p.get<art::InputTag>("SliceToPFPAssnLabel",  fSliceLabel);
  fPFPToTrackAssnLabel   = p.get<art::InputTag>("PFPToTrackAssnLabel",  fTrackLabel);
  fPFPToShowerAssnLabel  = p.get<art::InputTag>("PFPToShowerAssnLabel", fShowerLabel);

  std::cout << "[GetMyWireData] ClusterType=" << ClusterType
            << " getShowers=" << (getShowers ? "true" : "false") << "\n";
}

void GetMyWireData::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("wire_tree", "wire");
  fTree->Branch("run", &fRun);
  fTree->Branch("event", &fEvent);
  fTree->Branch("channel", &fCh);
  fTree->Branch("adc", &fadc);
  fTree->Branch("time", &fTime);
  fTree->Branch("tpc", &fTPC);
  fTree->Branch("wire", &fWire);
  fTree->Branch("plane", &fPlane);

  fTreeT = tfs->make<TTree>("track_tree", "track/slice/shower -> hits");
  fTreeT->Branch("run", &fRun);
  fTreeT->Branch("event", &fEvent);
  fTreeT->Branch("trackID", &fTrackID);
  fTreeT->Branch("type", &fType);
  fTreeT->Branch("length", &fLength);

  fTreeT->Branch("hitTime", &fPeakTime);
  fTreeT->Branch("hitAmplitude", &fAmplitude);
  fTreeT->Branch("hitIntegral", &fIntegral);
  fTreeT->Branch("hitCh", &fChHit);
  fTreeT->Branch("hitTPC", &fTPCHit);
  fTreeT->Branch("hitWire", &fWireHit);
  fTreeT->Branch("hitPlane", &fPlaneHit);
}

void GetMyWireData::analyze(art::Event const& e)
{
  fRun   = e.run();
  fEvent = e.event();

  auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  // ------------------------------
  // 1) RAW DIGITS -> wire_tree
  // ------------------------------
  auto wireHandle = e.getHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
  if (!wireHandle || wireHandle->empty()) {
    std::cout << "[GetMyWireData] No RawDigit in run " << fRun << " event " << fEvent << "\n";
    return;
  }

  for (const auto& wire : *wireHandle)
  {
    fCh = wire.Channel();
    auto wids = wireReadout.ChannelToWire(fCh);

    fTPC.clear(); fPlane.clear(); fWire.clear();
    fadc = wire.ADCs();

    fTime = clockData.TPCTick2TrigTime(0);

    for (auto const& wid : wids)
    {
      fTPC.push_back(wid.TPC);
      fPlane.push_back(wid.Plane);
      fWire.push_back(wid.Wire);
    }
    fTree->Fill();
  }

  // ----------------------------------------
  // 2) OBJETO -> Hits  (track_tree)
  // ----------------------------------------

  // common handles (pra Track/ Slice precisar)
  auto track_h = e.getHandle<std::vector<recob::Track>>(fTrackLabel);

  // =========
  // MODO TRACK
  // =========
  if (ClusterType == "Track")
  {
    if (!track_h) return;

    std::vector<art::Ptr<recob::Track>> tracks;
    art::fill_ptr_vector(tracks, track_h);

    art::FindManyP<recob::Hit> fmHits(track_h, e, fTrackHitAssnLabel);
    if (!fmHits.isValid()) {
      mf::LogWarning("GetMyWireData") << "No Track<->Hit assns for " << fTrackHitAssnLabel;
      return;
    }

    for (auto const& trk : tracks)
    {
      if (trk.isNull()) continue;

      auto hits = fmHits.at(trk.key());

      fTrackID = trk->ID();
      fLength  = trk->Length();
      fType    = 0;

      ResetHitVectors();
      FillHitVectors(hits);
      fTreeT->Fill();
    }

    if (getShowers)
    {
      auto shower_h = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
      if (!shower_h) return;

      std::vector<art::Ptr<recob::Shower>> showers;
      art::fill_ptr_vector(showers, shower_h);

      art::FindManyP<recob::Hit> fmHitsShower(shower_h, e, fShowerHitAssnLabel);
      if (!fmHitsShower.isValid()) {
        mf::LogWarning("GetMyWireData") << "No Shower<->Hit assns for " << fShowerHitAssnLabel;
        return;
      }

      for (auto const& shw : showers)
      {
        if (shw.isNull()) continue;
        auto hits = fmHitsShower.at(shw.key());

        fTrackID = shw->ID();
        fLength  = shw->Length();
        fType    = 1;

        ResetHitVectors();
        FillHitVectors(hits);
        fTreeT->Fill();
      }
    }

    return;
  }

  // =========
  // MODO SLICE
  // =========
  if (ClusterType == "Slice")
  {
    auto slice_h = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);
    auto pfp_h   = e.getHandle<std::vector<recob::PFParticle>>(fPFPLabel);

    if (!slice_h || !pfp_h || !track_h) {
      mf::LogWarning("GetMyWireData") << "Missing Slice/PFParticle/Track collections.";
      return;
    }

    // Track -> Hit (pra puxar hits das tracks dentro do slice)
    art::FindManyP<recob::Hit> fmHits(track_h, e, fTrackHitAssnLabel);
    if (!fmHits.isValid()) {
      mf::LogWarning("GetMyWireData") << "No Track<->Hit assns for " << fTrackHitAssnLabel;
      return;
    }

    // Slice -> PFParticle
    art::FindManyP<recob::PFParticle> slice_to_pfps(slice_h, e, fSliceToPFPAssnLabel);
    if (!slice_to_pfps.isValid()) {
      mf::LogWarning("GetMyWireData") << "No Slice<->PFParticle assns for " << fSliceToPFPAssnLabel;
      return;
    }

    // PFParticle -> Track
    art::FindManyP<recob::Track> pfp_to_tracks(pfp_h, e, fPFPToTrackAssnLabel);
    if (!pfp_to_tracks.isValid()) {
      mf::LogWarning("GetMyWireData") << "No PFParticle<->Track assns for " << fPFPToTrackAssnLabel;
      return;
    }

    // (opcional) PFParticle -> Shower + Shower->Hit
    std::unique_ptr<art::FindManyP<recob::Shower>> pfp_to_showers;
    std::unique_ptr<art::FindManyP<recob::Hit>>    fmHitsShower;
    art::Handle<std::vector<recob::Shower>> shower_h;

    if (getShowers)
    {
      shower_h = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
      if (shower_h)
      {
        auto tmp1 = std::make_unique<art::FindManyP<recob::Shower>>(pfp_h, e, fPFPToShowerAssnLabel);
        auto tmp2 = std::make_unique<art::FindManyP<recob::Hit>>(shower_h, e, fShowerHitAssnLabel);

        if (tmp1->isValid() && tmp2->isValid()) {
          pfp_to_showers = std::move(tmp1);
          fmHitsShower   = std::move(tmp2);
        } else {
          mf::LogWarning("GetMyWireData") << "Showers requested but PFParticle<->Shower or Shower<->Hit assns missing.";
        }
      }
    }

    std::vector<art::Ptr<recob::Slice>> slices;
    art::fill_ptr_vector(slices, slice_h);

    for (auto const& sl : slices)
    {
      if (sl.isNull()) continue;

      // Hits acumulados do slice (tracks + showers)
      std::vector<art::Ptr<recob::Hit>> hits_sel;
      std::unordered_set<size_t> seen;

      double sumLen = 0.0;
      bool anyObj = false;

      auto pfps = slice_to_pfps.at(sl.key());
      for (auto const& pfp : pfps)
      {
        if (pfp.isNull()) continue;

        // tracks do PFP
        auto trks = pfp_to_tracks.at(pfp.key());
        for (auto const& trk : trks)
        {
          if (trk.isNull()) continue;
          anyObj = true;
          sumLen += trk->Length();

          auto trkHits = fmHits.at(trk.key());
          for (auto const& h : trkHits)
          {
            if (h.isNull()) continue;
            if (!seen.insert(h.key()).second) continue;
            hits_sel.push_back(h);
          }
        }

        // showers do PFP (se habilitado e válido)
        if (pfp_to_showers && fmHitsShower)
        {
          auto shws = pfp_to_showers->at(pfp.key());
          for (auto const& shw : shws)
          {
            if (shw.isNull()) continue;
            anyObj = true;
            sumLen += shw->Length();

            auto shwHits = fmHitsShower->at(shw.key());
            for (auto const& h : shwHits)
            {
              if (h.isNull()) continue;
              if (!seen.insert(h.key()).second) continue;
              hits_sel.push_back(h);
            }
          }
        }
      }

      if (!anyObj) continue;
      if (hits_sel.empty()) continue;

      fTrackID = sl->ID();   // “ID lógico” do slice (seu uso anterior)
      fType    = 3;          // slice
      fLength  = sumLen;     // não existe Length do Slice -> uso soma das tracks/showers

      ResetHitVectors();
      FillHitVectors(hits_sel);
      fTreeT->Fill();
    }

    return;
  }

  mf::LogWarning("GetMyWireData") << "Unknown ClusterType='" << ClusterType << "'. Use 'Track' or 'Slice'.";
}

DEFINE_ART_MODULE(GetMyWireData)
