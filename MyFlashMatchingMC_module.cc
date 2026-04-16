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
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larcorealg/Geometry/WireReadoutGeom.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
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
#include <optional>
#include <functional>

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
        art::InputTag fTrackLabel; 
        art::InputTag fShowerLabel;//
        art::InputTag fSliceLabel;
        art::InputTag fPFPLabel;

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
        int fType;
        std::vector<int>   fTrackTIDs;
        std::vector<int>   fTrackPDGs;
        std::vector<float> fTrackW;

        bool getShowers = false;

        // --- pair branches ---
        int fPairFlashKey=-1;
        int fPairTrackKey=-1;
        int fPairTrackRecoID=-1;
        int fPairType=-1;

        float fOverlap=0.f;
        float fJaccard=0.f;
        int   fDomEqual=0;
        int fDomPDG=0;

        std::vector<int> fCommonTIDs; // só IDs (pra debug)
        std::vector<int> fCommonPDGs;

        // helpers
        std::unordered_map<int,double> BuildFlashMap(std::vector<art::Ptr<recob::OpHit>> const& ophits,
                                                    cheat::PhotonBackTrackerService& pbts,
                                                    cheat::ParticleInventoryService const& pis) const;

        std::unordered_map<int,double> BuildTrackMap(detinfo::DetectorClocksData const& clockData,
                                                     detinfo::DetectorPropertiesData const& detProp,
                                                    std::vector<art::Ptr<recob::Hit>> const& hits,
                                                    cheat::BackTrackerService const& bts,
                                                    cheat::ParticleInventoryService const& pis) const;


        double limitMinFlash = 0.0;
        double limitMinFlashSide = 0.0;
        double trackLength = 0.0;
        std::string DetectorZone;
        std::string ClusterType; //Slice,Track

        int nTtotal = 0;
        int nFTotal = 0 ;
};

MyFlashMatchingMC::MyFlashMatchingMC(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
    fFlashLabel = p.get<art::InputTag>("FlashLabel");
    fSliceLabel = p.get<art::InputTag>("SliceLabel");  
    fTrackLabel = p.get<art::InputTag>("TrackLabel");
    fShowerLabel = p.get<art::InputTag>("ShowerLabel");

    fPFPLabel = p.get<art::InputTag>("PFParticleLabel");

    fMinJaccard = p.get<double>("MinJaccard", 0.0);
    fMaxStore   = p.get<int>("MaxStore", 50);

    ClusterType=p.get<std::string>("ClusterType","Track");
    DetectorZone=p.get<std::string>("DetectorZone","All");
    limitMinFlash = p.get<double>("limitMinFlash",0.0);
    limitMinFlashSide = p.get<double>("limitMinFlashSide",0.0);
    trackLength = p.get<double>("trackLengthMin",0.0);

    getShowers = p.get<bool>("getShowers",false);
    if(getShowers)
        std::cout << "Getting Showers!!! " << std::endl;
    else
        std::cout << "Not Getting Showers!!! " << std::endl;

}

void MyFlashMatchingMC::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    //nessa tree é salva as informaceos do flash
    fTreeF = tfs->make<TTree>("treeF", "Flash truth content (MC)");
    fTreeF->Branch("run", &fRun); //run
    fTreeF->Branch("event", &fEvent); //event
    fTreeF->Branch("flashKey", &fFlashKey); //flashkey
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
    fTreeT->Branch("type", &fType);
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
    fTreeFT->Branch("type", &fPairType);
    fTreeFT->Branch("trackRecoID", &fPairTrackRecoID);
    fTreeFT->Branch("overlapMin", &fOverlap);
    fTreeFT->Branch("jaccardW", &fJaccard);
    fTreeFT->Branch("domEqual", &fDomEqual);
    fTreeFT->Branch("domPDG", &fDomPDG);
    fTreeFT->Branch("commonTIDs", &fCommonTIDs);
    fTreeFT->Branch("commonPDGs",  &fCommonPDGs);
    
}

std::unordered_map<int,double> MyFlashMatchingMC::BuildFlashMap(std::vector<art::Ptr<recob::OpHit>> const& ophits,
                                 cheat::PhotonBackTrackerService& pbts,cheat::ParticleInventoryService const& pis) const
{
    std::unordered_map<int,double> w;
    for (auto const& oph : ophits) //vare os hits
    {
        if(DetectorZone == "Positive" && oph->OpChannel() >= 80) continue;
        if(DetectorZone == "Negative" && oph->OpChannel() < 80) continue;
        // SDPs “crus” que contribuíram para esse OpHit
        auto sdps = pbts.OpHitToSimSDPs_Ps(oph); // vector<const sim::SDP*> --> retorna todo ponto de cintilicao simulado que contribuiu para esse ophit
        if (sdps.empty()) continue;

        // soma fotons totais (dentro desse OpHit)
        double totPhot = 0.0; // isso vai conter o numero de fotons criados que contribuiram para esse ophit
        for (auto const* sdp : sdps) //varre todos os pontos de cintilacao
        {
            if (!sdp) continue;
            if (sdp->trackID == 0) continue;
            totPhot += std::max(0.f, sdp->numPhotons); 
        }

        // distribui o PE do hit proporcional ao numPhotons por trackID
        for (auto const* sdp : sdps) //varre todos os pontos de cintilacao
        {
            int tid;
            if (!sdp) continue;
            if(sdp->trackID!=0)
            {
                tid = abs(pis.TrackIdToEveTrackId(abs(sdp->trackID)));   //AbsTID(sdp->trackID);
                if (tid < 0) continue;
            }
            else
            {
                continue;
            }
            
            double frac=0.0;
            if(totPhot>0) frac = std::max(0.f, sdp->numPhotons) / totPhot;
        
            w[tid] += oph->PE() * frac; //para cada track MC G4 id construi o a porcentagem do numero de photon eletron responsveis por esse ophit
        }
    }
    return w;
}


std::unordered_map<int,double> MyFlashMatchingMC::BuildTrackMap(detinfo::DetectorClocksData const& clockData,
                                detinfo::DetectorPropertiesData const& detProp,
                                std::vector<art::Ptr<recob::Hit>> const& hits,
                                cheat::BackTrackerService const& bts,cheat::ParticleInventoryService const& pis) const
{

    std::unordered_map<int,double> w;
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

    for (auto const& h : hits) // varre os hits
    {
        auto wireIDs = wireReadout.ChannelToWire(h->Channel());
        if (wireIDs.empty()) continue;
        auto const& wid = wireIDs.front();

        geo::PlaneID pid{wid.Cryostat, wid.TPC, wid.Plane};
        double xWire = wireReadout.Plane(pid).GetCenter().X();
        if (DetectorZone == "Positive" && xWire<0) continue;
        if (DetectorZone == "Negative" && xWire>0) continue;
       
        auto ides = bts.HitToTrackIDEs(clockData, h); // le todos os tracks G4 ID que contribuiram para esse hit
        for (auto const& ide : ides)
        {
            int tid;
            if(ide.trackID!=0)
            {
                tid = abs(pis.TrackIdToEveTrackId(abs(ide.trackID)));   //AbsTID(sdp->trackID);
                if (tid == 0) continue;
            }
            else
            {
                continue;
            }
            w[tid] += ide.energy; // peso físico no TPC
        }
    }
    return w;
}


void MyFlashMatchingMC::analyze(art::Event const& e)
{
    nFTotal = 0; //numero de flashs totais
    nTtotal = 0; //numero de tracks totais

    fRun   = e.run(); //indice da run
    fEvent = e.event(); //indice do evento

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e); // informacao do clock
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData); //informacao do detector

    auto& pbts = *art::ServiceHandle<cheat::PhotonBackTrackerService>(); //responsavel por buscar fotons da simulacao com base nos ophits
    auto const& bts  = *art::ServiceHandle<cheat::BackTrackerService const>(); //responsavel por buscar trajetorias MC com base nos hits
    auto const& pis  = *art::ServiceHandle<cheat::ParticleInventoryService const>(); //responsavel por obter informacoes da particula MC 

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

    art::FindManyP<recob::OpHit> fmOpHits(flash_h, e, fFlashLabel); // estrutura para vuscar os hits de um flash
    if (!fmOpHits.isValid())
    {
        mf::LogWarning("MyFlashMatchingMC") << "No OpFlash<->OpHit assns for " << fFlashLabel;
        return;
    }

    const int nF = (int)flashes.size(); // carrega o numero de flashs
    std::vector<std::unordered_map<int,double>> flashMaps(nF); // um vetor ( ... para cada flash ... ) de mapas.
                                                            //cada mapa associa um trackID(track aqui eh na sim MC G4) um valor de contribuicao para esse flash

    for (int f = 0; f < nF; ++f) // varre os flashes
    {
        auto const& fl = flashes[f];
        auto ophits = fmOpHits.at(fl.key()); // le os ophits desse flash

        //secao de cut de lado do detector---------------------------------------------------------------------
        auto flash1 = 0.0; 
        auto flash2 = 0.0;
        for(auto const& op: ophits)
        {
            int ch = op->OpChannel();
            if(ch<80)
            {
                flash1+=op->PE();          
            }
            else
            {
                flash2+=op->PE();     
            }
        }
        auto flashT=flash1+flash2;
        bool pass =  (flashT>limitMinFlash);
        if(DetectorZone=="Positive") pass = pass && (flash1>limitMinFlashSide);
        if(DetectorZone=="Negative") pass = pass && (flash2>limitMinFlashSide);
        if(!pass) continue;
        //-------------------------------------------------------------------------------------------------------
        auto wF = BuildFlashMap(ophits, pbts , pis); // construi o mapa para esse flash
        if (wF.empty()) continue;
        nFTotal++;
        Normalize(wF); // normaliza para soma 1
        flashMaps[f] = wF;

        fFlashKey     = fl.key(); // id
        fFlashTime    = fl->AbsTime(); // time
        fFlashTotalPE = fl->TotalPE(); // numero de photo-electrons

        double wmax=0.0;
        fFlashDomTID = DominantID(wF, wmax); // determina track id MC G4 com maior contribuicao
        fFlashDomW   = (float)wmax;
        auto const* p = pis.TrackIdToParticle_P(fFlashDomTID);
        fFlashDomPDG = p ? p->PdgCode() : 0;

        MapToVectors(wF, pis, fFlashTIDs, fFlashPDGs, fFlashW, fMaxStore); //preenche as variavias para salvar na tree
        fTreeF->Fill();
    }

    // ---- tracks ----
    //hora de varrer os tracks
    art::Handle<std::vector<recob::Track>> track_h;
    art::Handle<std::vector<recob::Shower>> shower_h;
    art::Handle<std::vector<recob::Slice>> slice_h;

    std::vector<art::Ptr<recob::Track>> tracks;
    std::vector<art::Ptr<recob::Slice>> slices;
    std::vector<art::Ptr<recob::Shower>> showers;
    std::optional<art::FindManyP<recob::Hit>> fmHits;
    std::optional<art::FindManyP<recob::Hit>> fmHitsShower;
    
    art::InputTag assnLabel;

    std::optional<art::FindManyP<recob::PFParticle>> slice_to_pfps;
    std::optional<art::FindManyP<recob::Track>> pfp_to_tracks;
    std::optional<art::FindManyP<recob::Shower>> pfp_to_showers;

    std::vector<art::Ptr<recob::PFParticle>> allpfps;
    std::unordered_map<size_t, art::Ptr<recob::PFParticle>> pfpMap;
    art::Handle<std::vector<recob::PFParticle>> pfp_h;


    // 1 secao setandos os parametros e produtos
    int nT,nS=0;
    if (ClusterType == "Track") // se estamos no tipo track
    {
        track_h = e.getHandle<std::vector<recob::Track>>(fTrackLabel);
        if (!track_h) 
        {
            mf::LogWarning("MyFlashMatchingMC") << "Cannot load Track: " << fTrackLabel;
            return;
        }
        assnLabel = fTrackLabel;
        fmHits.emplace(track_h, e, assnLabel);
        if (!fmHits->isValid()) 
        {
            mf::LogWarning("MyFlashMatchingMC") << "No <obj><->Hit assns for " << assnLabel;
            return;
        }
        art::fill_ptr_vector(tracks, track_h);
        nT = (int)tracks.size();
        
        if(getShowers) // se queremos showers tambem
        {
            shower_h = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
            if (!shower_h) {
                mf::LogWarning("MyFlashMatchingMC") << "Cannot load Shower: " << fShowerLabel;
                return;
            }
            fmHitsShower.emplace(shower_h, e, fShowerLabel);
            if (!fmHitsShower->isValid()) 
            {
                mf::LogWarning("MyFlashMatchingMC") << "No <obj><->Hit assns for " << fShowerLabel;
                return;
            }
            art::fill_ptr_vector(showers, shower_h);
            nS = (int)showers.size();
        }
    }
    else // se estamos no tipo SLICE
    {
        track_h = e.getHandle<std::vector<recob::Track>>(fTrackLabel);
        if (!track_h)
        {
            mf::LogWarning("MyFlashMatchingMC") << "Cannot load Track: " << fTrackLabel;
            return;
        }

        slice_h = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);
        if (!slice_h) {
            mf::LogWarning("MyFlashMatchingMC") << "Cannot load Slice: " << fSliceLabel;
            return;
        }
        assnLabel = fSliceLabel;
        fmHits.emplace(track_h, e, fTrackLabel);
        if (!fmHits->isValid())
        {
            mf::LogWarning("MyFlashMatchingMC") << "No Track<->Hit assns for " << fTrackLabel;
            return;
        }

        art::fill_ptr_vector(slices, slice_h);
        nT = (int)slices.size();
        pfp_h = e.getHandle<std::vector<recob::PFParticle>>(fPFPLabel);
        if (!pfp_h)
        {
            mf::LogWarning("MyFlashMatchingMC") << "Cannot load PFParticle: " << fPFPLabel;
            return;
        }

        art::fill_ptr_vector(allpfps, pfp_h);
        for (auto const& pfp : allpfps)
        {
            if (!pfp.isNull()) pfpMap[pfp->Self()] = pfp;
        }

        slice_to_pfps.emplace(slice_h, e, fSliceLabel);
        pfp_to_tracks.emplace(pfp_h, e, fTrackLabel);
        if (!slice_to_pfps->isValid() || !pfp_to_tracks->isValid())
        {
            mf::LogWarning("MyFlashMatchingMC") << "Missing Slice<->PFParticle or PFParticle<->Track associations.";
            return;
        }

        if (getShowers) // se queremos showers
        {
            shower_h = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
            if (!shower_h) 
            {
                mf::LogWarning("MyFlashMatchingMC") << "Cannot load Shower: " << fShowerLabel;
                return;
            }

            fmHitsShower.emplace(shower_h, e, fShowerLabel);
            if (!fmHitsShower->isValid()) 
            {
                mf::LogWarning("MyFlashMatchingMC") << "No Shower<->Hit assns for " << fShowerLabel;
                return;
            }

            // PFParticle -> Shower (assns produzidas pelo módulo de shower)
            pfp_to_showers.emplace(pfp_h, e, fShowerLabel);
            if (!pfp_to_showers->isValid()) 
            {
                mf::LogWarning("MyFlashMatchingMC") << "No PFParticle<->Shower assns for " << fShowerLabel;
                return;
            }
        }
    }

    // 2 secao lendo os dados
    int nTotal = nT + nS;
    std::vector<std::unordered_map<int,double>> trackMaps(nTotal); // um vetor ( ... para cada entidade(track/shower/slice) ... ) de mapas.
                                                            //cada mapa associa um trackID(track aqui eh na sim MC G4) um valor de contribuicao para esse cluster

    for (int t = 0; t < (nTotal); ++t) // varre os cluster
    {
        std::vector<art::Ptr<recob::Hit>> hits;
        if(ClusterType == "Track") // tipo track/shower
        {
            if(t<nT)
            {
                auto const& trk = tracks[t];
                if (trk->Length() < trackLength) continue;
                fTrackKey = trk.key();
                fTrackRecoID = trk->ID();
                hits = fmHits->at(trk.key());  
                fType = 0 ;
            }
            else
            {
                if(getShowers)
                {
                    int index = t-nT;
                    auto const& shw = showers[index];
                    if (shw->Length() < trackLength) continue;
                    fTrackKey = shw.key();
                    fTrackRecoID = shw->ID();
                    hits = fmHitsShower->at(shw.key());  
                    fType = 1 ;
                }
            }
            
        }
        else // tipo slice
        {
            auto const& slice = slices[t];
            fTrackKey    = slice.key();
            fTrackRecoID = slice->ID();
            fType        = 3;

            double myLength = 0.0;

            std::vector<art::Ptr<recob::Hit>> hits_sel;
            bool anyGoodObj = false;

            std::unordered_set<size_t> seenPFP;
            std::unordered_set<size_t> seenHits;

            auto seed_pfps = slice_to_pfps->at(slice.key());

            std::function<void(const art::Ptr<recob::PFParticle>&)> visitPFP;
            visitPFP = [&](const art::Ptr<recob::PFParticle>& pfp)
            {
                if (pfp.isNull()) return;

                const size_t self = pfp->Self();
                if (!seenPFP.insert(self).second) return;

                // Tracks associados a este PFP
                auto trks = pfp_to_tracks->at(pfp.key());
                for (auto const& trk : trks)
                {
                    if (trk.isNull()) continue;

                    myLength += trk->Length();
                    anyGoodObj = true;

                    auto trkHits = fmHits->at(trk.key());
                    for (auto const& h : trkHits)
                    {
                        if (h.isNull()) continue;
                        if (!seenHits.insert(h.key()).second) continue;
                        hits_sel.push_back(h);
                    }
                }
                // Showers associados a este PFP
                if (getShowers && pfp_to_showers && fmHitsShower)
                {
                    auto shws = pfp_to_showers->at(pfp.key());
                    for (auto const& shw : shws)
                    {
                        if (shw.isNull()) continue;

                        myLength += shw->Length();
                        anyGoodObj = true;

                        auto shwHits = fmHitsShower->at(shw.key());
                        for (auto const& h : shwHits)
                        {
                            if (h.isNull()) continue;
                            if (!seenHits.insert(h.key()).second) continue;
                            hits_sel.push_back(h);
                        }
                    }
                }

                // Desce para as filhas
                for (size_t dauID : pfp->Daughters())
                {
                    auto it = pfpMap.find(dauID);
                    if (it != pfpMap.end())
                    {
                        visitPFP(it->second);
                    }
                }
            };

            for (auto const& pfp : seed_pfps)
            {
                visitPFP(pfp);
            }

            if (!anyGoodObj) continue;
            if (hits_sel.empty()) continue;
            if (myLength < trackLength) continue;

            hits = std::move(hits_sel);
        }
     
        auto wT = BuildTrackMap(clockData, detProp, hits, bts ,pis); //consturi o mapa
        if (wT.empty()) continue;
        nTtotal++;     
        Normalize(wT); // normaliza para soma 1
        trackMaps[t] = wT;

        double wmax=0.0;
        fTrackDomTID = DominantID(wT, wmax); //determina track com contribuicao mais forte
        fTrackDomW   = (float)wmax;
        auto const* p = pis.TrackIdToParticle_P(fTrackDomTID);
        fTrackDomPDG = p ? p->PdgCode() : 0;

        MapToVectors(wT, pis, fTrackTIDs, fTrackPDGs, fTrackW, fMaxStore); // preenche os vetores para salvar na tree
        fTreeT->Fill();
    }

    // ---- candidates flash ↔ track ----
    //vamos varrer cada possivel match
    for (int f = 0; f < nF; ++f)  
    {
        if (flashMaps[f].empty()) continue; 
        for (int t = 0; t < nTotal; ++t) 
        {
            if (trackMaps[t].empty()) continue;
            double ov  = OverlapMin(flashMaps[f], trackMaps[t]); // calcula o overlap min
            if (ov <= 0.0) continue;

            double jac = JaccardWeighted(flashMaps[f], trackMaps[t]); // calcula o peso jaccard
            if (jac <= fMinJaccard) continue;

            fPairFlashKey   = flashes[f].key();
            if(ClusterType == "Track")
            {
                if(t<nT)
                {   
                    fPairTrackKey   = tracks[t].key();
                    fPairTrackRecoID= tracks[t]->ID();
                    fPairType = 0;
                }
                else
                {
                    if(getShowers)
                    {
                        int index = t-nT;
                        fPairTrackKey   = showers[index].key();
                        fPairTrackRecoID= showers[index]->ID();
                        fPairType = 1;
                    }
                }
            }
            else
            {
                fPairTrackKey   = slices[t].key();
                fPairTrackRecoID= slices[t]->ID();
                fPairType = 3;
            }
            fOverlap        = (float)ov;
            fJaccard        = (float)jac;

            double tmp=0.0;
            int domF = DominantID(flashMaps[f], tmp); // determina o  trackid geant4 dominante para flash
            int domT = DominantID(trackMaps[t], tmp); // determina o trackid geant dominatne para track
            fDomEqual = (AbsTID(domF) == AbsTID(domT)) ? 1 : 0; // determina se eh o mesmo ou 
            if (fDomEqual==1) 
            {
                auto const* p = pis.TrackIdToParticle_P(domF);
                fDomPDG = p ? p->PdgCode() : 0;
            } 
            else 
            {
                fDomPDG = -10000;
            }

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
                    auto const* p = pis.TrackIdToParticle_P(tid);
                    pdg = p ? p->PdgCode() : 0;
                    fCommonTIDs.push_back(tid);
                    fCommonPDGs.push_back(pdg);
                }
            }

            fTreeFT->Fill();
        }
    }

    std::string type_s  = (ClusterType == "Track") ? "tracks" : "slices";

    mf::LogInfo("MyFlashMatchingMC")
    << "Processed run " << fRun << " event " << fEvent
    << " (flashes=" << nFTotal << ", "<< type_s << "=" << nTtotal << " )";
}

DEFINE_ART_MODULE(MyFlashMatchingMC)
