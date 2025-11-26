#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TTree.h"

#include <iostream> 

#include "dunecore/DuneObj/OpDetDivRec.h" 
#include "larsim/MCCheater/ParticleInventoryService.h" 
#include "larsim/MCCheater/PhotonBackTrackerService.h"

#include "canvas/Persistency/Common/Ptr.h"          // art::Ptr
#include "canvas/Persistency/Common/PtrVector.h"    // art::fill_ptr_vector
#include "duneopdet/OpticalDetector/OpFlashSort.h"
#include <algorithm>    


class GetMCInfoPDSAnalyzer : public art::EDAnalyzer 
{
  public:
    explicit GetMCInfoPDSAnalyzer(fhicl::ParameterSet const& p);

    void analyze(art::Event const& e) override;
    void beginJob() override;

  private:
    art::InputTag fInputTag;

    TTree* fTree = nullptr;

    int fRun, fEvent;
    int fChannel;
    double fTime;
    int fFlashNumber;
    int fHitNumber;
    double fPeak;

    std::vector<int> fTrackID;
    std::vector<int> fPDGIDs;
    std::vector<double> fEnergy;
    std::vector<double> fPosIniX, fPosIniY, fPosIniZ;
    std::vector<double> fMomIniX, fMomIniY, fMomIniZ;

    std::string fOpFlashModuleLabel; 
};

GetMCInfoPDSAnalyzer::GetMCInfoPDSAnalyzer(fhicl::ParameterSet const& p):EDAnalyzer(p)
{
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("my_tree","my_tree");

    //fTree->Branch("run" , &fRun); --> como eh MC, nao preciso disso ( eu acho ...)
    fTree->Branch("event" , &fEvent);
    fTree->Branch("channel" , &fChannel);
    fTree->Branch("time" , &fTime);
    fTree->Branch("flash_number" , &fFlashNumber);
    fTree->Branch("hit_number" , &fHitNumber);
    fTree->Branch("peak" , &fPeak);
   
    fTree->Branch("track", &fTrackID);
    fTree->Branch("pdg", &fPDGIDs);
    
    // Posição inicial
    fTree->Branch("x", &fPosIniX);
    fTree->Branch("y", &fPosIniY);
    fTree->Branch("z", &fPosIniZ);

    // Momento inicial
    fTree->Branch("px", &fMomIniX);
    fTree->Branch("py", &fMomIniY);
    fTree->Branch("pz", &fMomIniZ);

    // Energia_inicial
    fTree->Branch("e", &fEnergy);

    fOpFlashModuleLabel = p.get<std::string>("OpFlashModuleLabel");
}


void GetMCInfoPDSAnalyzer::beginJob()
{
 
}

void GetMCInfoPDSAnalyzer::analyze(const art::Event &evt)
{

    fEvent = evt.event();

    art::ServiceHandle<cheat::PhotonBackTrackerService > pbts;
    art::ServiceHandle<cheat::ParticleInventoryService> pis;

    //pegar os flashes do evento
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    auto FlashHandle = evt.getHandle< std::vector< recob::OpFlash > >(fOpFlashModuleLabel);
    if (FlashHandle) 
    {
      art::fill_ptr_vector(flashlist, FlashHandle);
      std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
    }
    else 
    {
      mf::LogWarning("GetMCInfoPDSAnalyzer") << "Cannot load any flashes. Failing";
      return;
    }

    int number_flashs = flashlist.size();
    std::cout << "Number of flashs: " << number_flashs << std::endl;
    for(int j=0;j<number_flashs;j++)
    {

        fFlashNumber = j;

        std::vector< art::Ptr<recob::OpHit> > hitFromFlash = pbts->OpFlashToOpHits_Ps(flashlist[j]);

        int number_hits=hitFromFlash.size();

        std::cout << "Number of hits: " << number_hits << std::endl;
        
        for(int i=0;i<number_hits;i++)
        {

            fHitNumber = i;

            art::Ptr<recob::OpHit> this_hit = hitFromFlash[i];
            //std::cout << "Channel: " << this_hit->OpChannel() << std::endl;

            fChannel = this_hit->OpChannel();
            fTime =  this_hit->PeakTime();
            fPeak= this_hit->Amplitude();

            std::vector<int> eveIds = pbts->OpHitToTrackIds(this_hit); // antes era: pbts->OpHitToEveTrackIds(this_hit);

            /*  //vamos testar a diferenca:
            std::cout << "Channel: " << fChannel << std::endl;
            std::cout << "Start Time: " << this_hit->StartTime() << std::endl;
            std::cout << "Peak Time: " << this_hit->PeakTime() << std::endl;
            std::cout << "Peak Time Abs: " << this_hit->PeakTimeAbs() << std::endl; */

            fTrackID.clear();
            fPDGIDs.clear();
            fEnergy.clear();
            fPosIniX.clear(); fPosIniY.clear(); fPosIniZ.clear();
            fMomIniX.clear(); fMomIniY.clear(); fMomIniZ.clear();

            for (int trackId : eveIds)
            {
                fTrackID.push_back(trackId);

                const simb::MCParticle* particle = pis->TrackIdToParticle_P(trackId);
                if(particle) 
                {
                    fPDGIDs.push_back(particle->PdgCode());

                    fPosIniX.push_back(particle->Vx());
                    fPosIniY.push_back(particle->Vy());
                    fPosIniZ.push_back(particle->Vz());

                    fMomIniX.push_back(particle->Px());
                    fMomIniY.push_back(particle->Py());
                    fMomIniZ.push_back(particle->Pz());

                    fEnergy.push_back(particle->E());

                     /*  std::cout << "PDG Code: " << mc->PdgCode()
                            << " -- Energy: " << mc->E()
                            << " -- Origin: " << mc->Process() << std::endl; */
                }

            }
            fTree->Fill();
        }
    }

}

DEFINE_ART_MODULE(GetMCInfoPDSAnalyzer)
