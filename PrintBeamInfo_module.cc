////////////////////////////////////////////////////////////////////////
// File: PrintBeamInfo_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"

// Se teu release for antigo e der erro no include acima, troca por:
// #include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include "TTree.h"

#include <vector>
#include <iostream>

class PrintBeamInfo : public art::EDAnalyzer {
public:
    explicit PrintBeamInfo(fhicl::ParameterSet const& p);

    void analyze(art::Event const& evt) override;

private:
    void reset();

    art::InputTag fBeamTag;

    TTree* fTree;

    unsigned int fRun;
    unsigned int fSubRun;
    unsigned long long fEvent;

    bool fHasBeamProduct;
    int fNBeamEvents;

    int fTimingTrigger;

    double fTOF;
    int fTOFChan;

    std::vector<double> fMultipleTOFs;
    std::vector<int> fMultipleTOFChans;
    int fNMultipleTOFs;
    int fNMultipleTOFChans;

    bool fIsMatched;
    int fNBeamTracks;

    int fHighPressureCKov;
    int fLowPressureCKov;

    int fNRecoBeamMomenta;
    std::vector<double> fRecoBeamMomenta;
};

PrintBeamInfo::PrintBeamInfo(fhicl::ParameterSet const& p)
    : art::EDAnalyzer(p)
    , fBeamTag(p.get<art::InputTag>("BeamTag", art::InputTag("beamevent")))
    , fTree(nullptr)
{
    art::ServiceHandle<art::TFileService> tfs;

    fTree = tfs->make<TTree>("beam_tree", "ProtoDUNE beam information");

    fTree->Branch("run",    &fRun);
    fTree->Branch("subrun", &fSubRun);
    fTree->Branch("event",  &fEvent);

    fTree->Branch("has_beam_product", &fHasBeamProduct);
    fTree->Branch("n_beam_events",    &fNBeamEvents);

    fTree->Branch("timing_trigger", &fTimingTrigger);

    fTree->Branch("tof",      &fTOF);
    fTree->Branch("tof_chan", &fTOFChan);

    fTree->Branch("n_multiple_tofs",      &fNMultipleTOFs);
    fTree->Branch("multiple_tofs",        &fMultipleTOFs);
    fTree->Branch("n_multiple_tof_chans", &fNMultipleTOFChans);
    fTree->Branch("multiple_tof_chans",   &fMultipleTOFChans);

    fTree->Branch("is_matched",   &fIsMatched);
    fTree->Branch("n_beam_tracks", &fNBeamTracks);

    fTree->Branch("high_pressure_ckov", &fHighPressureCKov);
    fTree->Branch("low_pressure_ckov",  &fLowPressureCKov);

    fTree->Branch("n_reco_beam_momenta", &fNRecoBeamMomenta);
    fTree->Branch("reco_beam_momenta",   &fRecoBeamMomenta);
}

void PrintBeamInfo::reset()
{
    fRun = 0;
    fSubRun = 0;
    fEvent = 0;

    fHasBeamProduct = false;
    fNBeamEvents = 0;

    fTimingTrigger = -999;

    fTOF = -999.0;
    fTOFChan = -999;

    fMultipleTOFs.clear();
    fMultipleTOFChans.clear();
    fNMultipleTOFs = 0;
    fNMultipleTOFChans = 0;

    fIsMatched = false;
    fNBeamTracks = 0;

    fHighPressureCKov = -999;
    fLowPressureCKov  = -999;

    fNRecoBeamMomenta = 0;
    fRecoBeamMomenta.clear();
}

void PrintBeamInfo::analyze(art::Event const& evt)
{
    reset();

    fRun = evt.run();
    fSubRun = evt.subRun();
    fEvent = evt.event();

    auto beamHandle =
        evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamTag);

    if (!beamHandle) {
        std::cout << "Run " << evt.run()
                  << ", SubRun " << evt.subRun()
                  << ", Event " << evt.event()
                  << " : no beam::ProtoDUNEBeamEvent with tag "
                  << fBeamTag.encode() << "\n";

        fTree->Fill();
        return;
    }

    fHasBeamProduct = true;
    fNBeamEvents = beamHandle->size();

    if (beamHandle->empty()) {
        std::cout << "Run " << evt.run()
                  << ", SubRun " << evt.subRun()
                  << ", Event " << evt.event()
                  << " : empty beam::ProtoDUNEBeamEvent vector\n";

        fTree->Fill();
        return;
    }

    auto prod = beamHandle->at(0);

    fTimingTrigger = prod.GetTimingTrigger();

    fTOF = prod.GetTOF();
    fTOFChan = prod.GetTOFChan();

    fMultipleTOFs = prod.GetTOFs();
    fMultipleTOFChans = prod.GetTOFChans();

    fNMultipleTOFs = fMultipleTOFs.size();
    fNMultipleTOFChans = fMultipleTOFChans.size();

    fIsMatched = prod.CheckIsMatched();
    fNBeamTracks = prod.GetNBeamTracks();

    fHighPressureCKov = prod.GetCKov0Status();
    fLowPressureCKov  = prod.GetCKov1Status();

    fRecoBeamMomenta = prod.GetRecoBeamMomenta();
    fNRecoBeamMomenta = fRecoBeamMomenta.size();

    fTree->Fill();
}

DEFINE_ART_MODULE(PrintBeamInfo)