#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RawData/RawDigit.h"

#include "TTree.h"

#include <iostream> 


class GetMyWireData : public art::EDAnalyzer 
{
  public:
    explicit GetMyWireData(fhicl::ParameterSet const& p);

    void analyze(art::Event const& e) override;
    void beginJob() override;

  private:
    art::InputTag fInputTag;
    int fRun, fEvent;
    TTree* fTree = nullptr;
    int fCh;
    std::vector<short> fadc;
};

GetMyWireData::GetMyWireData(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
{
}

void GetMyWireData::beginJob() 
{
        std::cout  <<"#######################\n[GetMyWireData] beginJob() called\n#####################\n";
        art::ServiceHandle<art::TFileService> tfs;
        fTree = tfs->make<TTree>("wire_tree", "wire");
        fTree->Branch("run" , &fRun);
        fTree->Branch("event" , &fEvent);
        fTree->Branch("channel" , &fCh);
        fTree->Branch("adc" , &fadc);

}

void GetMyWireData::analyze(art::Event const& e)
{
        
        fRun = e.run();
        fEvent = e.event();
 
        auto wireHandle = e.getHandle<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq", "Detsim"));

        if (!wireHandle || wireHandle->empty())
        {
                std::cout << "No Wire data in event: run " << fRun << ", event " << fEvent << std::endl;
                return;
        }
        for (const auto& wire : *wireHandle) 
        {
            fCh = wire.Channel();
            fadc = wire.ADCs();
            //std::cout << fCh << std::endl;     
            fTree->Fill();
        }

}

DEFINE_ART_MODULE(GetMyWireData)
