#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/Track.h"

#include "TTree.h"

#include <vector>

class GetPandoraTracksAnalyzer : public art::EDAnalyzer {
public:
    explicit GetPandoraTracksAnalyzer(fhicl::ParameterSet const& p);

    void beginJob() override;
    void analyze(art::Event const& e) override;

private:
    art::InputTag fTrackTag;

    TTree* fTree;

    int   run, event;
    int   track_id;

    std::vector<float> vx, vy, vz;  // posições do traço
};


GetPandoraTracksAnalyzer::GetPandoraTracksAnalyzer(fhicl::ParameterSet const& p)
    : EDAnalyzer(p),
      fTrackTag(p.get<std::string>("TrackTag", "pandoraTrack"))
{}


void GetPandoraTracksAnalyzer::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTree = tfs->make<TTree>("track_tree", "Pandora Tracks");

    fTree->Branch("run",     &run);
    fTree->Branch("event",   &event);
    fTree->Branch("track_id", &track_id);

    fTree->Branch("vx", &vx);
    fTree->Branch("vy", &vy);
    fTree->Branch("vz", &vz);
}


void GetPandoraTracksAnalyzer::analyze(art::Event const& e)
{
    run   = e.run();
    event = e.event();

    vx.clear();
    vy.clear();
    vz.clear();

    // pega os tracks
    auto trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackTag);

    track_id = 0;

    for (auto const& trk : *trackHandle) {

        vx.clear();
        vy.clear();
        vz.clear();

        size_t N = trk.NumberTrajectoryPoints();

        for (size_t i = 0; i < N; ++i) {
            if (!trk.HasValidPoint(i)) continue;

            auto const& p = trk.LocationAtPoint(i);

            vx.push_back(p.X());
            vy.push_back(p.Y());
            vz.push_back(p.Z());
        }

        fTree->Fill();
        track_id++;
    }
}

DEFINE_ART_MODULE(GetPandoraTracksAnalyzer)
