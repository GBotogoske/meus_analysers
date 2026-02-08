////////////////////////////////////////////////////////////////////////
// File:        DuneSimpleSliceReader_module.cc
// Author: Gabriel
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcore/CoreUtils/ServiceUtil.h" 

#include "dunecore/DuneObj/OpDetDivRec.h" 
#include "duneopdet/OpticalDetector/OpFlashSort.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "art/Utilities/make_tool.h"

#include "TTree.h"

#include <vector>
#include <iostream>

#include "mydune/utils/my_utils.hh"
#include "mydune/utils/MyMatch.hh"

#include <optional>

class getVis : public art::EDAnalyzer
 {
    public:
        explicit getVis(fhicl::ParameterSet const& p);

        void beginJob() override;
        void analyze(art::Event const& e) override;

    private:

        //declaracao da TTree
        TTree* fTreeV;
        TTree* fTreeA;

        art::ServiceHandle<geo::Geometry> geo;
        phot::SemiAnalyticalModel const* fSAM;

        double x,y,z;
        int channel;
        double vis;
};


getVis::getVis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{


    // carrega o servico de visibilidade otica

    fSAM = new phot::SemiAnalyticalModel(
                p.get<fhicl::ParameterSet>("vuvhitspars"), 
                p.get<fhicl::ParameterSet>("vishitspars"),
                std::shared_ptr<phot::OpticalPath>(art::make_tool<phot::OpticalPath>(p.get<fhicl::ParameterSet>("OpticalPathTool"))), 
                p.get<bool>("do_refl", false), 
                p.get<bool>("do_include_anode_refl", false),
                p.get<bool>("do_include_xe_absorption", false)
                ); 

    
    /*
   for (unsigned int ch = 0; ch < fPVS->NOpChannels(); ch++)
    {
        auto const& detGeo = geo->OpDetGeoFromOpDet(ch);
        auto const& pos = detGeo.GetCenter();

        std::cout << "CH " << ch
                << " at (" << pos.X()
                << ", " << pos.Y()
                << ", " << pos.Z() << ")\n";
    }

    geo::Point_t p1{300.0, 0.0, 100.0};
    for (unsigned int ch = 0; ch < fPVS->NOpChannels(); ch++) 
    {
        float vis = fPVS->GetVisibility(p1, ch);
        if (vis > 0) std::cout << "Channel " << ch << " sees vis = " << vis << std::endl;
        else std::cout << "Channel " << ch << " sees vis (-1) = " << vis << std::endl;
    } */  
    
}

void getVis::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTreeV = tfs->make<TTree>("treeVis", "");
    fTreeV->Branch("x", &x);
    fTreeV->Branch("y", &y);
    fTreeV->Branch("z", &z);
    fTreeV->Branch("vis", &vis);
    fTreeV->Branch("ch", &channel);

    double XMin= -391.564606;
    double XMax= 464.195190;
    double YMin= -65.479393;
    double YMax= 725.724670;
    double ZMin= -199.330917;
    double ZMax= 661.863403;

    double step=10;
    for (x=XMin;x<=XMax;x=x+step)
    {
        for (y=YMin;y<=YMax;y=y+step)
        {
            for (z=ZMin;z<=ZMax;z=z+step)
            {
                std::vector<double> direct_visibilities(160, 0.0);
                direct_visibilities.reserve(160);
                geo::Point_t point{x, y, z};
                fSAM->detectedDirectVisibilities(direct_visibilities, point);
       
                for(channel=0;channel<160;++channel)
                {   
                    vis = direct_visibilities[channel];
                    fTreeV->Fill();
                }     
     
            }    
        }        
    }
    
    fTreeA = tfs->make<TTree>("treePos", "");
    fTreeA->Branch("x", &x);
    fTreeA->Branch("y", &y);
    fTreeA->Branch("z", &z);
    fTreeA->Branch("ch", &channel);
    
    for(channel=0;channel<160;++channel)
    {   
        auto opDet = geo->OpDetGeoFromOpDet(channel);
        auto const c = opDet.GetCenter(); 
        x=c.X();
        y=c.Y();
        z=c.Z();
        fTreeA->Fill();
    }  
}

void getVis::analyze(art::Event const& e)
{

    
}



DEFINE_ART_MODULE(getVis)
