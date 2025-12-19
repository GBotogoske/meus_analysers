////////////////////////////////////////////////////////////////////////
// File:        MyFlashMatchingMC_module.cc
// Purpose:     Flash↔Track truth-matching via MC (weighted overlap)
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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "duneopdet/OpticalDetector/OpFlashSort.h"

#include "TTree.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>

static void Normalize(std::unordered_map<int,double>& w)
{
  double s = 0.0;
  for (auto const& [tid, ww] : w) s += ww;
  if (s <= 0.0) return;
  for (auto& [tid, ww] : w) ww /= s;
}


namespace {
  inline int AbsTID(int tid) { return std::abs(tid); }

  // overlap = sum min(wF, wT)
  double OverlapMin(std::unordered_map<int,double> const& F,
                    std::unordered_map<int,double> const& T)
{
    auto const* small = &F;
    auto const* big   = &T;
    if (F.size() > T.size()) 
    {
        std::swap(small, big);
    }

    double s = 0.0;
    for (auto const& [tid, w] : *small) 
    {
        auto it = big->find(tid);
        if (it == big->end()) continue;
        s += std::min(w, it->second);
    }
    return s;
}

  // jaccard ponderado = sum min / sum max  (0..1)
  double JaccardWeighted(std::unordered_map<int,double> const& F,
                         std::unordered_map<int,double> const& T)
  {
    double smin = 0.0, smax = 0.0;

    for (auto const& [tid, wf] : F) 
    {
        double wt = 0.0;
        auto it = T.find(tid);
        if (it != T.end()) wt = it->second;
        smin += std::min(wf, wt);
        smax += std::max(wf, wt);
    }
    for (auto const& [tid, wt] : T) 
    {
        if (F.count(tid)) continue;
        smax += wt;
    }
    return (smax > 0.0) ? (smin / smax) : 0.0;
  }

  // dominante (maior peso)
  int DominantID(std::unordered_map<int,double> const& w, double& wmax)
  {
        int best = 0;
        double bestW = 0.0;
        for (auto const& [tid, ww] : w)
        {
            if (ww > bestW) 
            {  
                bestW = ww; 
                best = tid; 
            }
        }
        wmax = bestW;
        return best;
  }

  // converte mapa em vetores (IDs, PDGs, pesos) ordenados por peso desc (top N)
  void MapToVectors(std::unordered_map<int,double> const& w,
                    cheat::ParticleInventoryService const& pis,
                    std::vector<int>& ids,
                    std::vector<int>& pdgs,
                    std::vector<float>& weights,
                    int maxStore)
  {
        std::vector<std::pair<int,double>> tmp;
        tmp.reserve(w.size());
        for (auto const& kv : w) tmp.push_back(kv);

        std::sort(tmp.begin(), tmp.end(),
                [](auto const& a, auto const& b){ return a.second > b.second; });

        int n = (maxStore > 0) ? std::min<int>(maxStore, tmp.size()) : (int)tmp.size();
        ids.clear(); pdgs.clear(); weights.clear();
        ids.reserve(n); pdgs.reserve(n); weights.reserve(n);

        for (int i = 0; i < n; ++i)
        {
            int tid = tmp[i].first;
            double ww = tmp[i].second;
            ids.push_back(tid);
            weights.push_back((float)ww);

            if (tid == 0) 
            {
                pdgs.push_back(0);
                continue;
            }
            auto const* p = pis.TrackIdToParticle_P(tid);
            pdgs.push_back(p ? p->PdgCode() : 0);
            
        }
  }
}

class MyFlashMatchingMC : public art::EDAnalyzer
{
    public:
        explicit MyFlashMatchingMC(fhicl::ParameterSet const& p);

        void beginJob() override;
        void analyze(art::Event const& e) override;

    private:
        art::InputTag fFlashLabel;
        art::InputTag fTrackLabel; // (você pode continuar chamando de SliceLabel no FHiCL se quiser)

        double fMinJaccard;   // corte pra salvar par
        int    fMaxStore;     // truncar vetores no TTree

        // Trees
        TTree* fTreeF  = nullptr; // flash
        TTree* fTreeT  = nullptr; // track
        TTree* fTreeFT = nullptr; // flash-track pairs

        // event
        int fRun=0, fEvent=0;

        // --- flash branches ---
        int fFlashKey=-1;
        double fFlashTime=0.0;
        double fFlashTotalPE=0.0;
        int fFlashDomTID=0, fFlashDomPDG=0;
        float fFlashDomW=0.f;
        std::vector<int>   fFlashTIDs;
        std::vector<int>   fFlashPDGs;
        std::vector<float> fFlashW;

        // --- track branches ---
        int fTrackKey=-1;
        int fTrackRecoID=-1; // Track::ID()
        int fTrackDomTID=0, fTrackDomPDG=0;
        float fTrackDomW=0.f;
        std::vector<int>   fTrackTIDs;
        std::vector<int>   fTrackPDGs;
        std::vector<float> fTrackW;


        // --- pair branches ---
        int fPairFlashKey=-1;
        int fPairTrackKey=-1;
        int fPairTrackRecoID=-1;

        float fOverlap=0.f;
        float fJaccard=0.f;
        int   fDomEqual=0;

        std::vector<int> fCommonTIDs; // só IDs (pra debug)
        std::vector<int> fCommonPDGs;

        // helpers
        std::unordered_map<int,double> BuildFlashMap(std::vector<art::Ptr<recob::OpHit>> const& ophits,
                                                    cheat::PhotonBackTrackerService& pbts) const;

        std::unordered_map<int,double> BuildTrackMap(detinfo::DetectorClocksData const& clockData,
                                                    std::vector<art::Ptr<recob::Hit>> const& hits,
                                                    cheat::BackTrackerService const& bts) const;
};

MyFlashMatchingMC::MyFlashMatchingMC(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
    fFlashLabel = p.get<art::InputTag>("FlashLabel");
    // Para não quebrar teu FHiCL antigo:
    // use "TrackLabel" se existir, senão usa "SliceLabel".
    fTrackLabel = p.get<art::InputTag>("TrackLabel", p.get<art::InputTag>("SliceLabel"));

    fMinJaccard = p.get<double>("MinJaccard", 0.0);
    fMaxStore   = p.get<int>("MaxStore", 50);
}

void MyFlashMatchingMC::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTreeF = tfs->make<TTree>("treeF", "Flash truth content (MC)");
    fTreeF->Branch("run", &fRun);
    fTreeF->Branch("event", &fEvent);
    fTreeF->Branch("flashKey", &fFlashKey);
    fTreeF->Branch("time", &fFlashTime);
    fTreeF->Branch("totalPE", &fFlashTotalPE);
    fTreeF->Branch("domTID", &fFlashDomTID);
    fTreeF->Branch("domPDG", &fFlashDomPDG);
    fTreeF->Branch("domW", &fFlashDomW);
    fTreeF->Branch("tids", &fFlashTIDs);
    fTreeF->Branch("pdgs", &fFlashPDGs);
    fTreeF->Branch("w", &fFlashW);

    fTreeT = tfs->make<TTree>("treeT", "Track truth content (MC)");
    fTreeT->Branch("run", &fRun);
    fTreeT->Branch("event", &fEvent);
    fTreeT->Branch("trackKey", &fTrackKey);
    fTreeT->Branch("trackRecoID", &fTrackRecoID);
    fTreeT->Branch("domTID", &fTrackDomTID);
    fTreeT->Branch("domPDG", &fTrackDomPDG);
    fTreeT->Branch("domW", &fTrackDomW);
    fTreeT->Branch("tids", &fTrackTIDs);
    fTreeT->Branch("pdgs", &fTrackPDGs);
    fTreeT->Branch("w", &fTrackW);

    fTreeFT = tfs->make<TTree>("treeFT", "Flash-Track candidates (weighted truth overlap)");
    fTreeFT->Branch("run", &fRun);
    fTreeFT->Branch("event", &fEvent);
    fTreeFT->Branch("flashKey", &fPairFlashKey);
    fTreeFT->Branch("trackKey", &fPairTrackKey);
    fTreeFT->Branch("trackRecoID", &fPairTrackRecoID);
    fTreeFT->Branch("overlapMin", &fOverlap);
    fTreeFT->Branch("jaccardW", &fJaccard);
    fTreeFT->Branch("domEqual", &fDomEqual);
    fTreeFT->Branch("commonTIDs", &fCommonTIDs);
    fTreeFT->Branch("commonPDGs",  &fCommonPDGs);
    
}

std::unordered_map<int,double>
MyFlashMatchingMC::BuildFlashMap(std::vector<art::Ptr<recob::OpHit>> const& ophits,
                                 cheat::PhotonBackTrackerService& pbts) const
{
  std::unordered_map<int,double> w;

  for (auto const& oph : ophits) {
    // SDPs “crus” que contribuíram para esse OpHit
    auto sdps = pbts.OpHitToSimSDPs_Ps(oph); // vector<const sim::SDP*> :contentReference[oaicite:2]{index=2}
    if (sdps.empty()) continue;

    // soma fotões totais (dentro desse OpHit)
    double totPhot = 0.0;
    for (auto const* sdp : sdps) {
      if (!sdp) continue;
      if (sdp->trackID == 0) continue;
      totPhot += std::max(0.f, sdp->numPhotons); // numPhotons existe em sim::SDP :contentReference[oaicite:3]{index=3}
    }

    // fallback: se não tem fotões, volta pro teu método antigo (igualitário)
    if (totPhot <= 0.0) {
      auto tids = pbts.OpHitToTrackIds(oph);
      if (tids.empty()) continue;
      std::unordered_set<int> uniq;
      for (int tidRaw : tids) {
        int tid = AbsTID(tidRaw);
        if (tid) uniq.insert(tid);
      }
      if (uniq.empty()) continue;

      double share = oph->PE() / (double)uniq.size();
      for (int tid : uniq) w[tid] += share;
      continue;
    }

    // distribui o PE do hit proporcional ao numPhotons por trackID
    for (auto const* sdp : sdps) {
      if (!sdp) continue;
      int tid = AbsTID(sdp->trackID);
      if (tid == 0) continue;

      double frac = std::max(0.f, sdp->numPhotons) / totPhot;
      w[tid] += oph->PE() * frac;
    }
  }

  return w;
}


std::unordered_map<int,double>
MyFlashMatchingMC::BuildTrackMap(detinfo::DetectorClocksData const& clockData,
                                 std::vector<art::Ptr<recob::Hit>> const& hits,
                                 cheat::BackTrackerService const& bts) const
{
    std::unordered_map<int,double> w;
    for (auto const& h : hits)
    {
        auto ides = bts.HitToTrackIDEs(clockData, h);
        for (auto const& ide : ides)
        {
            int tid = AbsTID(ide.trackID);
            if (tid == 0) continue;
            w[tid] += ide.energy; // peso físico no TPC
        }
    }
    return w;
}

void MyFlashMatchingMC::analyze(art::Event const& e)
{
    fRun   = e.run();
    fEvent = e.event();

    auto const clockData =
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    auto& pbts = *art::ServiceHandle<cheat::PhotonBackTrackerService>();
    auto const& bts  = *art::ServiceHandle<cheat::BackTrackerService const>();
    auto const& pis  = *art::ServiceHandle<cheat::ParticleInventoryService const>();

    // ---- flashes ----
    auto flash_h = e.getHandle<std::vector<recob::OpFlash>>(fFlashLabel);
    if (!flash_h) 
    {
        mf::LogWarning("MyFlashMatchingMC") << "Cannot load OpFlash: " << fFlashLabel;
        return;
    }

    std::vector<art::Ptr<recob::OpFlash>> flashes;
    art::fill_ptr_vector(flashes, flash_h);
    std::sort(flashes.begin(), flashes.end(), recob::OpFlashPtrSortByPE);

    art::FindManyP<recob::OpHit> fmOpHits(flash_h, e, fFlashLabel);
    if (!fmOpHits.isValid())
    {
        mf::LogWarning("MyFlashMatchingMC") << "No OpFlash<->OpHit assns for " << fFlashLabel;
        return;
    }

    const int nF = (int)flashes.size();
    std::vector<std::unordered_map<int,double>> flashMaps(nF);

    for (int f = 0; f < nF; ++f) 
    {
        auto const& fl = flashes[f];
        auto ophits = fmOpHits.at(fl.key());

        auto wF = BuildFlashMap(ophits, pbts);
        Normalize(wF);
        flashMaps[f] = wF;

        fFlashKey     = fl.key();
        fFlashTime    = fl->Time();
        fFlashTotalPE = fl->TotalPE();

        double wmax=0.0;
        fFlashDomTID = DominantID(wF, wmax);
        fFlashDomW   = (float)wmax;
        if (fFlashDomTID != 0) 
        {
            auto const* p = pis.TrackIdToParticle_P(fFlashDomTID);
            fFlashDomPDG = p ? p->PdgCode() : 0;
        } 
        else 
        {
            fFlashDomPDG = 0;
        }


        MapToVectors(wF, pis, fFlashTIDs, fFlashPDGs, fFlashW, fMaxStore);
        fTreeF->Fill();
    }

    // ---- tracks ----
    auto track_h = e.getHandle<std::vector<recob::Track>>(fTrackLabel);
    if (!track_h) 
    {
        mf::LogWarning("MyFlashMatchingMC") << "Cannot load Track: " << fTrackLabel;
        return;
    }
    art::FindManyP<recob::Hit> fmHits(track_h, e, fTrackLabel);
    if (!fmHits.isValid()) 
    {
        mf::LogWarning("MyFlashMatchingMC") << "No Track<->Hit assns for " << fTrackLabel;
        return;
    }

    std::vector<art::Ptr<recob::Track>> tracks;
    art::fill_ptr_vector(tracks, track_h);
    const int nT = (int)tracks.size();

    std::vector<std::unordered_map<int,double>> trackMaps(nT);

    for (int t = 0; t < nT; ++t)
    {
        auto const& trk = tracks[t];
        auto hits = fmHits.at(trk.key());

        auto wT = BuildTrackMap(clockData, hits, bts);
        Normalize(wT);
        trackMaps[t] = wT;

        fTrackKey    = trk.key();
        fTrackRecoID = trk->ID();

        double wmax=0.0;
        fTrackDomTID = DominantID(wT, wmax);
        fTrackDomW   = (float)wmax;
        if (fTrackDomTID != 0)
        {
            auto const* p = pis.TrackIdToParticle_P(fTrackDomTID);
            fTrackDomPDG = p ? p->PdgCode() : 0;
        } else 
        {
            fTrackDomPDG = 0;
        }

        MapToVectors(wT, pis, fTrackTIDs, fTrackPDGs, fTrackW, fMaxStore);
        fTreeT->Fill();
    }

    // ---- candidates flash ↔ track ----
    for (int f = 0; f < nF; ++f) 
    {
        for (int t = 0; t < nT; ++t) 
        {
            double ov  = OverlapMin(flashMaps[f], trackMaps[t]);
            if (ov <= 0.0) continue;

            double jac = JaccardWeighted(flashMaps[f], trackMaps[t]);
            if (jac < fMinJaccard) continue;

            fPairFlashKey   = flashes[f].key();
            fPairTrackKey   = tracks[t].key();
            fPairTrackRecoID= tracks[t]->ID();
            fOverlap        = (float)ov;
            fJaccard        = (float)jac;

            double tmp=0.0;
            int domF = DominantID(flashMaps[f], tmp);
            int domT = DominantID(trackMaps[t], tmp);
            fDomEqual = (AbsTID(domF) == AbsTID(domT)) ? 1 : 0;

            // comuns (só pra debug, top por min(wF,wT))
            fCommonTIDs.clear();
            fCommonPDGs.clear();
            {
                std::vector<std::pair<int,double>> common;
                for (auto const& [tid, wf] : flashMaps[f]) 
                {
                    auto it = trackMaps[t].find(tid);
                    if (it == trackMaps[t].end()) continue;
                    common.push_back({tid, std::min(wf, it->second)});
                }

                std::sort(common.begin(), common.end(),
                            [](auto const& a, auto const& b){ return a.second > b.second; });

                int n = std::min<int>(fMaxStore, (int)common.size());
                fCommonTIDs.reserve(n);
                fCommonPDGs.reserve(n);

                for (int i = 0; i < n; ++i) 
                {
                    int tid = common[i].first;

                    // PDG do MCParticle (protege nullptr)
                    int pdg = 0;
                    if (tid != 0) 
                    {
                        auto const* p = pis.TrackIdToParticle_P(tid);
                        pdg = p ? p->PdgCode() : 0;
                    }
                    fCommonTIDs.push_back(tid);
                    fCommonPDGs.push_back(pdg);
                }
            }


            fTreeFT->Fill();
        }
    }

    mf::LogInfo("MyFlashMatchingMC")
    << "Processed run " << fRun << " event " << fEvent
    << " (flashes=" << nF << ", tracks=" << nT << ")";
}

DEFINE_ART_MODULE(MyFlashMatchingMC)
