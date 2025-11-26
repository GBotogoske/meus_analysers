#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"


#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TTree.h"

#include <iostream> 
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "dunecore/DuneObj/OpDetDivRec.h" 
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/ParticleInventoryService.h" 


class GetAllTimesAnalyzer : public art::EDAnalyzer 
{
  public:
    explicit GetAllTimesAnalyzer(fhicl::ParameterSet const& p);

    void analyze(art::Event const& e) override;
    void beginJob() override;

  private:
    art::InputTag fInputTag;

    TTree* fTreeP = nullptr; // Photons
    TTree* fTreeW = nullptr; // Waveforms
    TTree* fTreeH = nullptr; // Hits
    TTree* fTreeHw = nullptr; // HitsWire

  
    int fRun, fEvent;
    int fOfflineChannel;
    double fTimePhoton;
    double fTimeWaveform;
    double fTimeHit;
    double fTimeHitWire;
};

GetAllTimesAnalyzer::GetAllTimesAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer(p),
    fInputTag(p.get<std::string>("input_tag"))
{
}

void GetAllTimesAnalyzer::beginJob() 
{
        art::ServiceHandle<art::TFileService> tfs;

        fTreeP = tfs->make<TTree>("photon_tree", "Photons");
        fTreeP->Branch("run" , &fRun);
        fTreeP->Branch("event" , &fEvent);
        fTreeP->Branch("timestamp" , &fTimePhoton);

        fTreeW = tfs->make<TTree>("waveform_tree", "Waveforms");
        fTreeW->Branch("run" , &fRun);
        fTreeW->Branch("event" , &fEvent);
        fTreeW->Branch("timestamp" , &fTimeWaveform);

        fTreeH = tfs->make<TTree>("hit_tree", "Hits");
        fTreeH->Branch("run" , &fRun);
        fTreeH->Branch("event" , &fEvent);
        fTreeH->Branch("timestamp" , &fTimeHit);

      /*   fTreeHw = tfs->make<TTree>("wirehit_tree", "HitsWire");
        fTreeHw->Branch("run" , &fRun);
        fTreeHw->Branch("event" , &fEvent);
        fTreeHw->Branch("timestamp" , &fTimeHitWire); */
    
}

void GetAllTimesAnalyzer::analyze(art::Event const& e)
{
        fRun = e.run();
        fEvent = e.event();

        auto wfHandle = e.getHandle<std::vector<raw::OpDetWaveform>>(fInputTag);

        auto PDFastSimHandle = e.getHandle<std::vector<sim::SimPhotonsLite>>(art::InputTag("PDFastSim::G4")); // Ajuste a tag

        auto HitHandle = e.getHandle<std::vector<recob::OpHit>>(art::InputTag("ophitspe::myFlash"));

        //auto HitWireHandle = e.getHandle<std::vector<recob::Hit>>(art::InputTag("gaushit::myCharge"));


        if (!wfHandle || wfHandle->empty()) 
        {
            std::cout << "No pds waveform data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        }

        if (!PDFastSimHandle || PDFastSimHandle->empty()) 
        {
            std::cout << "No PDFast data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        }

        if (!HitHandle || HitHandle->empty()) 
        {
            std::cout << "No Hit data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        }

        /*if (!HitWireHandle || HitWireHandle->empty()) 
        {
            std::cout << "No Wire Hit data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        } */

        int maxChannel = 80;


        for (const auto& wf : *wfHandle) 
        {   
            if (static_cast<int>(wf.ChannelNumber()) <= maxChannel)
            {
                fTimeWaveform = wf.TimeStamp();
                fTreeW->Fill();
            }
        }

        for (const auto& this_photon: *PDFastSimHandle)
        {
            if (this_photon.OpChannel <= maxChannel)
            {
                for (const auto& [tick, nphotons] : this_photon.DetectedPhotons)
                {
                    fTimePhoton = tick;
                    fTreeP->Fill();
                }
            }
            
        }
        
        for (const auto& this_hit: *HitHandle)
        {
            if (this_hit.OpChannel() <= maxChannel)
            {
                fTimeHit =this_hit.StartTime();
                fTreeH->Fill();
            }
        }
        
        /* for (const auto& this_hit_wire: *HitWireHandle)
        {
            
            fTimeHit =this_hit_wire.PeakTime();
            fTreeHw->Fill();
            
        }
         */
        std::cout << "Processed event: run " << fRun << ", event " << fEvent << std::endl;
}

DEFINE_ART_MODULE(GetAllTimesAnalyzer)
