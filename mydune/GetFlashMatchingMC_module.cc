#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/MCCheater/ParticleInventoryService.h" 
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"

#include <iostream> // Para debug opcional

#include "dunecore/DuneObj/OpDetDivRec.h" //novo

#include "duneopdet/OpticalDetector/OpFlashSort.h"


class GetFlashMatchingMC : public art::EDAnalyzer 
{
  public:
    explicit GetFlashMatchingMC(fhicl::ParameterSet const& p);

    void analyze(art::Event const& e) override;
    void beginJob() override;

  private:

    TTree* fTreeF = nullptr;
    TTree* fTreeS = nullptr;
  
    int fRun, fEvent;
    int fflashID,fsliceID;
    std::vector<int> ftracksID;
    std::vector<int> fPDGsID;
    art::InputTag fFlashLabel;  
    art::InputTag fSliceLabel; 

    TTree* fTreeFS = nullptr;

    int fMatchedFlashID;
    int fMatchedSliceID;
    std::vector<int> fCommonTracks;
    std::vector<int> fCommonPDGs;

};

GetFlashMatchingMC::GetFlashMatchingMC(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
    fFlashLabel= p.get<art::InputTag>("FlashLabel");
    fSliceLabel = p.get<art::InputTag>("SliceLabel");  

}

void GetFlashMatchingMC::beginJob() 
{
        art::ServiceHandle<art::TFileService> tfs;
        fTreeF = tfs->make<TTree>("treeF", "my_tree_F");
        fTreeF->Branch("run" , &fRun);
        fTreeF->Branch("event" , &fEvent);
        fTreeF->Branch("flashID" , &fflashID);
        fTreeF->Branch("trackID" , &ftracksID);
        fTreeF->Branch("pdgID" , &fPDGsID);

        fTreeS = tfs->make<TTree>("treeS", "my_tree_S");
        fTreeS->Branch("run" , &fRun);
        fTreeS->Branch("event" , &fEvent);
        fTreeS->Branch("sliceID" , &fsliceID);
        fTreeS->Branch("trackID" , &ftracksID);
        fTreeS->Branch("pdgID" , &fPDGsID);

        fTreeFS = tfs->make<TTree>("treeFS", "Flash-Slice Matching");
        fTreeFS->Branch("run"  , &fRun);
        fTreeFS->Branch("event", &fEvent);
        fTreeFS->Branch("flashID", &fMatchedFlashID);
        fTreeFS->Branch("sliceID", &fMatchedSliceID);
        fTreeFS->Branch("trackID", &fCommonTracks);
        fTreeFS->Branch("pdgID",   &fCommonPDGs);
}

void GetFlashMatchingMC::analyze(art::Event const& e)
{
    fRun = e.run();
    fEvent = e.event();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    art::ServiceHandle<cheat::PhotonBackTrackerService > pbts;
    art::ServiceHandle<cheat::BackTrackerService> bts;
    art::ServiceHandle<cheat::ParticleInventoryService> pis;

    //pegar os flashes do evento
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    auto FlashHandle = e.getHandle< std::vector< recob::OpFlash > >(fFlashLabel);
    if (FlashHandle) 
    {
      art::fill_ptr_vector(flashlist, FlashHandle);
      std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
    }
    else 
    {
      mf::LogWarning("GetFlashMatchingMC") << "Cannot load any flashes. Failing";
      return;
    }

    int number_flashs = flashlist.size();
    std::vector<std::vector<int>> flashTracks(number_flashs);

    //std::cout << "Number of flashs: " << number_flashs << std::endl;
    for(int j=0;j<number_flashs;j++)
    {
        std::vector< art::Ptr<recob::OpHit> > hitFromFlash = pbts->OpFlashToOpHits_Ps(flashlist[j]);
        int number_hits=hitFromFlash.size();
        //std::cout << "Number of Ophits: " << number_hits << std::endl;
        fflashID=j;
        std::unordered_set<int> usedTrackIDs;
        ftracksID.clear();
        fPDGsID.clear();
        usedTrackIDs.clear();
        
        for(int i=0;i<number_hits;i++)
        {
            art::Ptr<recob::OpHit> this_hit = hitFromFlash[i];
            //std::cout<< "#OPHIT: " <<i << std::endl;
            std::vector<int> eveIds = pbts->OpHitToTrackIds(this_hit);
            for (auto const& trackId : eveIds) 
            {
                if (usedTrackIDs.count(trackId))
                {
                    continue;  
                }
                
                const simb::MCParticle* particle = pis->TrackIdToParticle_P(trackId);
                if(!particle) 
                {
                    continue;
                }
                
                usedTrackIDs.insert(trackId);
                ftracksID.push_back(trackId);
                fPDGsID.push_back(particle->PdgCode());
                //std::cout << "Flash TrackID = " << trackId << " --- Pdg Id =  "  << particle->PdgCode() <<"\n";
            
            }
        }
        flashTracks[j] = ftracksID;
        fTreeF->Fill();
    }

    auto slice_h = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
    std::vector<art::Ptr<recob::Slice>> slicelist;
    if (!slice_h) 
    {
      mf::LogWarning("GetFlashMatchingMC") << "Cannot load any slices. Failing";
      return;
    }

    size_t n_slice = slice_h->size();
    std::vector<std::vector<int>> sliceTracks(n_slice);
    //std::cout << "Number of Slices: " << n_slice << std::endl;

    art::FindManyP<recob::Hit> Hits_from_Slice(slice_h, e, fSliceLabel);
    for (size_t j = 0; j < n_slice; j++)
    {
        auto hits = Hits_from_Slice.at(j); 
        int number_hits=hits.size();
        //std::cout << "Number of hits: " << number_hits << std::endl;

        fsliceID=j;
        std::unordered_set<int> usedTrackIDs;
        ftracksID.clear();
        fPDGsID.clear();
        usedTrackIDs.clear();

        for(int i=0;i<number_hits;i++)
        {
            //std::cout<< "#WireHIT: " << i << std::endl;
            art::Ptr<recob::Hit> this_hit = hits[i];
            auto trackIDEs = bts->HitToTrackIDEs(clockData, this_hit);   
            for (auto const& ide : trackIDEs) 
            {
                int trackId = ide.trackID;
                if (usedTrackIDs.count(trackId))
                {
                    continue;  
                }
                const simb::MCParticle* particle = pis->TrackIdToParticle_P(trackId);
                if(!particle) 
                {
                    continue;
                }
                usedTrackIDs.insert(trackId);
                ftracksID.push_back(trackId);
                fPDGsID.push_back(particle->PdgCode());
            }
        }
        sliceTracks[j] = ftracksID;
        fTreeS->Fill();
    }

    for (int f = 0; f < number_flashs; f++)
    {
        for (size_t s = 0; s < n_slice; s++)
        {
            // limpar buffers
            fCommonTracks.clear();
            fCommonPDGs.clear();

            // transformar sliceTracks em set para busca rápida
            std::unordered_set<int> sliceSet(sliceTracks[s].begin(), sliceTracks[s].end());

            std::unordered_set<int> usedCommon;
            // intersect
            for (int trk : flashTracks[f])
            {
                if (usedCommon.count(trk))
                {
                    continue;
                }
                    
                if (sliceSet.count(trk))
                {
                    usedCommon.insert(trk);
                    fCommonTracks.push_back(trk);

                    const simb::MCParticle* p = pis->TrackIdToParticle_P(trk);
                    fCommonPDGs.push_back(p ? p->PdgCode() : 0);
                }
            }

            // se existe interseção → matching é real
            if (!fCommonTracks.empty())
            {
                fMatchedFlashID = f;
                fMatchedSliceID = s;
                fTreeFS->Fill();
            }
        }
    }


    std::cout << "Processed event: run " << fRun << ", event " << fEvent << std::endl;
}

DEFINE_ART_MODULE(GetFlashMatchingMC)
