#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"

#include <iostream> 

#include "lardata/DetectorInfoServices/DetectorClocksService.h"


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
    TTree* fTreeT = nullptr;
    int fCh;
    int fTrackID;
    double fTime;
    std::vector<int> fTPC,fWire,fPlane;
    std::vector<short> fadc;

    double fLength;
    std::vector<float> fPeakTime;
    std::vector <float> fAmplitude;
    std::vector <float> fIntegral;
    std::vector <int> fChHit;
    std::vector<int> fTPCHit,fWireHit,fPlaneHit;


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
        fTree->Branch("time" , &fTime);
        fTree->Branch("tpc" , &fTPC);
        fTree->Branch("wire" , &fWire);
        fTree->Branch("plane" , &fPlane);


        fTreeT = tfs->make<TTree>("track_tree", "track");
        fTreeT->Branch("run" , &fRun);
        fTreeT->Branch("event" , &fEvent);
        fTreeT->Branch("trackID" , &fTrackID);
        fTreeT->Branch("length" , &fLength);
        fTreeT->Branch("hitTime" , &fPeakTime);
        fTreeT->Branch("hitAmplitude" , &fAmplitude);
        fTreeT->Branch("hitIntegral" , &fIntegral);
        fTreeT->Branch("hitCh" , &fChHit);
        fTreeT->Branch("hitTPC" , &fTPCHit);
        fTreeT->Branch("hitWire" , &fWireHit);
        fTreeT->Branch("hitPlane" , &fPlaneHit);

        
}

void GetMyWireData::analyze(art::Event const& e)
{
        
    fRun = e.run();
    fEvent = e.event();

    auto wireHandle = e.getHandle<std::vector<raw::RawDigit>>(art::InputTag("tpcrawdecoder", "daq", "Detsim"));
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);


    if (!wireHandle || wireHandle->empty())
    {
            std::cout << "No Wire data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
    }
    for (const auto& wire : *wireHandle) 
    {
        fCh = wire.Channel();   
        auto wids = wireReadout.ChannelToWire(fCh);
        fTPC.clear();
        fPlane.clear();
        fWire.clear();
        fadc.clear();
        fadc = wire.ADCs();
        fTime = clockData.TPCTick2TrigTime(0);
        //
        for (auto const& wid : wids) 
        {
        // wid tem: Cryostat, TPC, Plane, Wire
            fTPC.push_back(wid.TPC);
            fPlane.push_back(wid.Plane);
            fWire.push_back(wid.Wire);
        }
        fTree->Fill();
    }

    //Preenchendo a segunda ttre com informacoes do track --> hits
    art::Handle<std::vector<recob::Track>> track_h = e.getHandle<std::vector<recob::Track>>("pandoraTrack");
    std::vector<art::Ptr<recob::Track>> tracks;
    art::fill_ptr_vector(tracks, track_h);
    int nT = (int)tracks.size();
    art::FindManyP<recob::Hit> fmHits(track_h, e, "pandoraTrack");
    for (int t = 0; t < nT; ++t)
    {
        std::vector<art::Ptr<recob::Hit>> hits;
        auto const& trk = tracks[t];
        fTrackID = trk->ID();
        hits = fmHits.at(trk.key());  
        fLength = trk->Length();

        fPeakTime.clear();
        fAmplitude.clear();
        fIntegral.clear();
        fChHit.clear();
        fTPCHit.clear();
        fWireHit.clear();
        fPlaneHit.clear();

        for (auto const& h : hits)
        {
            fPeakTime.push_back(h->PeakTime());
            fAmplitude.push_back(h->PeakAmplitude());
            fIntegral.push_back(h->Integral());
            fChHit.push_back(h->Channel());
            auto whit=h->WireID();
            fTPCHit.push_back(whit.TPC);
            fPlaneHit.push_back(whit.Plane);
            fWireHit.push_back(whit.Wire);
        }
        
        fTreeT->Fill();
    }

    art::Handle<std::vector<recob::Shower>> shower_h = e.getHandle<std::vector<recob::Shower>>("pandoraShower");
    std::vector<art::Ptr<recob::Shower>> showers;
    art::fill_ptr_vector(showers, shower_h);
    int nS = (int)showers.size();
    art::FindManyP<recob::Hit> fmHitsShower(shower_h, e, "pandoraShower");

    for (int t = 0; t < nS; ++t)
    {
        std::vector<art::Ptr<recob::Hit>> hits;
        auto const& trk = showers[t];
        fTrackID = trk->ID();
        hits = fmHitsShower.at(trk.key());  
        fLength = trk->Length();

        fPeakTime.clear();
        fAmplitude.clear();
        fIntegral.clear();
        fChHit.clear();
        fTPCHit.clear();
        fWireHit.clear();
        fPlaneHit.clear();

        for (auto const& h : hits)
        {
            fPeakTime.push_back(h->PeakTime());
            fAmplitude.push_back(h->PeakAmplitude());
            fIntegral.push_back(h->Integral());
            fChHit.push_back(h->Channel());
            auto whit=h->WireID();
            fTPCHit.push_back(whit.TPC);
            fPlaneHit.push_back(whit.Plane);
            fWireHit.push_back(whit.Wire);
        }
        
        fTreeT->Fill();
    }

}

DEFINE_ART_MODULE(GetMyWireData)
