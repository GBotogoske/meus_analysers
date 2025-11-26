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
#include "lardataobj/RecoBase/OpWaveform.h"

#include "dunecore/DuneObj/OpDetDivRec.h" 
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larsim/MCCheater/ParticleInventoryService.h" 


class GetWfHitsAnalyzer : public art::EDAnalyzer 
{
  public:
    explicit GetWfHitsAnalyzer(fhicl::ParameterSet const& p);

    void analyze(art::Event const& e) override;
    void beginJob() override;

  private:
    art::InputTag fInputTag;

    TTree* fTreeW = nullptr; // Waveforms
    TTree* fTreeH = nullptr; // Hits
  
    int fRun, fEvent;
    int fOfflineChannel;
    double fTimestamp;
    std::vector<float> fADCValue;

    int fOpChannelHit;
    double 	fPeakTimeAbs;
    double 	fPeakTime;
    double 	fStartTime;
    double 	fRiseTime;
    unsigned short 	fFrame;
    double 	fWidth ;
    double 	fArea;
    double 	fAmplitude;
    double 	fPE;
    double 	fFastToTotal ;


};

GetWfHitsAnalyzer::GetWfHitsAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
}

void GetWfHitsAnalyzer::beginJob() 
{
    art::ServiceHandle<art::TFileService> tfs;

    fTreeW = tfs->make<TTree>("waveform_tree", "Waveforms");
    fTreeW->Branch("run" , &fRun);
    fTreeW->Branch("event" , &fEvent);
    fTreeW->Branch("offline_channel" , &fOfflineChannel);
    fTreeW->Branch("timestamp" , &fTimestamp);
    fTreeW->Branch("adc" , &fADCValue);      


    fTreeH = tfs->make<TTree>("hit_tree", "Hits");
    fTreeH->Branch("run" , &fRun);
    fTreeH->Branch("event" , &fEvent);
    fTreeH->Branch("OpChannel" , &fOpChannelHit);
    fTreeH->Branch("PeakTimeAbs" , &fPeakTimeAbs);
    fTreeH->Branch("PeakTime" , &fPeakTime);
    fTreeH->Branch("StartTime" , &fStartTime);
    fTreeH->Branch("RiseTime" , &fRiseTime);
    fTreeH->Branch("Frame" , &fFrame);
    fTreeH->Branch("Width" , &fWidth);
    fTreeH->Branch("Area" , &fArea);
    fTreeH->Branch("Amplitude" , &fAmplitude);
    fTreeH->Branch("PE" , &fPE);
    fTreeH->Branch("FastToTotal" , &fFastToTotal);
   
}

void GetWfHitsAnalyzer::analyze(art::Event const& e)
{
        fRun = e.run();
        fEvent = e.event();

        auto wfHandle = e.getHandle<std::vector<recob::OpWaveform>>(art::InputTag("opdec", "", "myFlash"));
        if (!wfHandle || wfHandle->empty()) 
        {
            std::cout << "No OpWaveform data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        }
    
        
        auto HitHandle = e.getHandle<std::vector<recob::OpHit>>(art::InputTag("ophitspe::myFlash"));
        if (!HitHandle || HitHandle->empty()) 
        {
            std::cout << "No OpHit data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        }

        for (const auto& wf : *wfHandle) 
        {
            fOfflineChannel = wf.Channel();      
            fTimestamp      = wf.TimeStamp();    
            std::vector<float> signal = wf.Signal();
            fADCValue = signal;
            
            fTreeW->Fill();
        }

        for (const auto& hit : *HitHandle) 
        {
            fOpChannelHit = hit.OpChannel(); 
            fPeakTime = hit.PeakTime(); 
            fPeakTimeAbs = hit.PeakTimeAbs();
            fStartTime=hit.StartTime();
            fRiseTime=hit.RiseTime();
            fFrame=hit.Frame();
            fWidth=hit.Width();
            fArea=hit.Area();
            fAmplitude=hit.Amplitude();
            fPE=hit.PE();
            fFastToTotal=hit.FastToTotal();
            fTreeH->Fill();
        }
        
}

DEFINE_ART_MODULE(GetWfHitsAnalyzer)
