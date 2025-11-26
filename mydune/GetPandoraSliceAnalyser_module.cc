////////////////////////////////////////////////////////////////////////
// File:        DuneSimpleSliceReader_module.cc
// Author: Gabriel
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/Simulation/LArG4Parameters.h"

#include "dunecore/DuneObj/OpDetDivRec.h" 
#include "duneopdet/OpticalDetector/OpFlashSort.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"

#include <vector>
#include <iostream>

#include "mydune/utils/my_utils.hh"

class GetPandoraSliceAnalyser : public art::EDAnalyzer
 {
    public:
        explicit GetPandoraSliceAnalyser(fhicl::ParameterSet const& p);

        void beginJob() override;
        void analyze(art::Event const& e) override;

    private:
        //declarao dos nomes dos produtos
        art::InputTag fSliceLabel; 
        art::InputTag fPFPLabel;
        art::InputTag fTrackLabel;
        art::InputTag fCaloLabel;
        art::InputTag fFlashLabel;

        //este aqui eh para converter ADC em e-
        std::vector<float> _cal_area_const;

        //declaracao da TTree
        TTree* fTree;
        int run, event;
        int sliceID;
        int nSlices;
        int pfpID;
        int nPFPs;
        int nTracksPFP;
        int trackID;
        std::vector<float> v_dEdx;
        std::vector<float> v_dQdx;
        std::vector<float> v_pitch;
        std::vector<float> v_x;
        std::vector<float> v_y;
        std::vector<float> v_z;

        ::art::ServiceHandle<geo::Geometry> geo;
        double drift_length;
        double drift_speed;
        double electronlife;
        double W_LAr;


        std::vector<QCluster> getQClusters(art::Event const& e);
        void getFlashs(art::Event const& e);
};


GetPandoraSliceAnalyser::GetPandoraSliceAnalyser(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
    //seta os labels para ler os produtos
    fSliceLabel = p.get<art::InputTag>("SliceLabel");  
    fPFPLabel = p.get<art::InputTag>("PFParticleLabel");
    fTrackLabel = p.get<art::InputTag>("TrackLabel");
    fCaloLabel = p.get<art::InputTag>("CalorimetryLabel");
    fFlashLabel = p.get<art::InputTag>("FlashLabel");

    _cal_area_const    = p.get<std::vector<float>>("CalAreaConstants"); // PEGANDO OS VALORES PADRAO, TA CERTO??

    std::cout << "CalAreaConstants = "
          << _cal_area_const[0] << ", "
          << _cal_area_const[1] << ", "
          << _cal_area_const[2] << std::endl;


    //geometria do detector
    int nTPCs = geo->TotalNTPC();
    int nCrio = geo->Ncryostats();
   
    std::cout << "nCryo: " << nCrio << std::endl;
    std::cout << "nTPCs: " << nTPCs << std::endl;

    geo::TPCID id(0,1);
    auto const& tpc = geo->TPC(id);
    this->drift_length = tpc.DriftDistance(); // ~ 360 cm

    std::cout << "Drift Distance: " << this->drift_length << std::endl;

}

void GetPandoraSliceAnalyser::beginJob()
{
    /* art::ServiceHandle<art::TFileService> fs;

    fTree = fs->make<TTree>("sliceTree","Simple slice info");
    fTree->Branch("run",      &run,      "run/I");
    fTree->Branch("event",    &event,    "event/I");
    fTree->Branch("sliceID",  &sliceID,  "sliceID/I");
    fTree->Branch("nSlices",  &nSlices,  "nSlices/I");
    fTree->Branch("pfpID",  &pfpID,  "pfpID/I");
    fTree->Branch("nPFPs",  &nPFPs,  "nPFPs/I");
    fTree->Branch("nTracksPFP", &nTracksPFP, "nTracksPFP/I");
    fTree->Branch("trackID",    &trackID,    "trackID/I");
    fTree->Branch("dEdx", &v_dEdx);
    fTree->Branch("dQdx", &v_dQdx);
    fTree->Branch("pitch", &v_pitch);
    fTree->Branch("x", &v_x);
    fTree->Branch("y", &v_y);
    fTree->Branch("z", &v_z); */

}

void GetPandoraSliceAnalyser::analyze(art::Event const& e)
{
    run    = e.id().run();
    event  = e.id().event();

    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
    art::ServiceHandle<sim::LArG4Parameters const> g4param;

    drift_speed = det_prop.DriftVelocity(); //~ 0.15 cm/us
    electronlife = det_prop.ElectronLifetime(); //35000 us
    W_LAr=g4param->Wph(); //19.5 eV

    std::cout << "electron vd := " << drift_speed << std::endl;
    std::cout << "electron lifetime:=  " << electronlife << std::endl;
    std::cout << "W LAr:=  " << W_LAr << std::endl;

    auto QLigths = getQClusters(e);
    getFlashs(e);
    
    std::cout << "Slice ID := " << QLigths[0].sliceID << std::endl; // para o compilador nao reclamar que nao esta sendo usado
}


std::vector<QCluster> GetPandoraSliceAnalyser::getQClusters(art::Event const& e)
{
   

    // pegar produtos importantes
    auto slice_h = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel); // SLICE
    auto pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel); //PFParticles
    auto track_h = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel); // tracks
    auto calo_h = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloLabel); //calorimetria

    art::FindManyP<recob::PFParticle> slice_to_pfps(slice_h, e, fSliceLabel); //pega a associacao de PFParticles nos slices
    art::FindManyP<recob::Track> pfp_to_tracks(pfp_h, e, fTrackLabel); // pega a associacao de tracks das PFParticles
    art::FindManyP<anab::Calorimetry> trk_to_calo(track_h, e, fCaloLabel); // pega as info de calorimetria dos tracks

    //numero de slices
    nSlices = slice_h->size();
    //mf::LogInfo("DuneSimpleSliceReader") << "Evento " << event << " tem " << nSlices << " slices.";

    //cria um vetor para encher os slices
    std::vector<art::Ptr<recob::Slice>> slices;
    art::fill_ptr_vector(slices, slice_h);

    std::vector<QCluster> QLigths;

    //varre os slices
    for (size_t i = 0; i < slices.size(); i++)
    {
        QCluster this_qlight;

        //para cada slice ...
        auto sl = slices[i];
        //... pega as pfParticels
        std::vector<art::Ptr<recob::PFParticle>> pfps = slice_to_pfps.at(i);

        sliceID = sl->ID(); // ID do slice
        nPFPs = pfps.size(); // quantOS PFParticels tem o slice
        //mf::LogInfo("DuneSimpleSliceReader") << " Slice " << sliceID << " tem " << nPFPs << " PFParticles.";
        
        //varre todos os pfps deste slice
        for (auto const& pfp : pfps) 
        {
            pfpID = pfp->Self(); // ID do PFParticle
            // pegar tracks associados ao PFP
            std::vector<art::Ptr<recob::Track>> tracks = pfp_to_tracks.at(pfp.key());
            nTracksPFP = tracks.size();

            //varre todos tracks deste PFParticle
            for (auto const& trk : tracks)
            {
                trackID = trk->ID(); //ID do track

                //pega as info de calorimetria
                std::vector<art::Ptr<anab::Calorimetry>> calos = trk_to_calo.at(trk.key());
                //vetores para salvar as informacoes
                v_dEdx.clear();
                v_dQdx.clear();
                v_pitch.clear();
                v_x.clear();
                v_y.clear();
                v_z.clear();

                //No codigo original tem uma secao de verificar o a posicao do track e se cruzou o plano de fios, e se o track eh uncontained
                //------------------------------------------------------------------------------------
                //############### BLA BLA BLA BLA BLA BLA ############################
                //------------------------------------------------------------------------------------

                // ---------------------------------------------------------------------------------------------
                //sao 3 planos --> Determinar o melhor plano para fazer a associaco com os flashs

                std::vector<art::Ptr<anab::Calorimetry>> calo_plane(3, art::Ptr<anab::Calorimetry>());
                for (auto const& calo : calos)
                {
                    int plane = calo->PlaneID().Plane;

                    if (plane < 0 || plane > 2) {
                        //std::cout << "Calo inválido (plane = " << plane << ")" << std::endl;
                        continue;
                    }

                    calo_plane[plane] = calo;
                }
                
                int bestPlane = -1;
                size_t maxSize = 0;

                for (int pl = 0; pl < 3; ++pl)
                {
                    if (calo_plane[pl].isNull()) continue;
                    size_t nhits = calo_plane[pl]->dEdx().size();

                    if (nhits > maxSize)
                    {
                        maxSize = nhits;
                        bestPlane = pl;
                    }
                }
                if (bestPlane < 0) 
                {
                    std::cout << "Track sem calorimetria válida em nenhum plano." << std::endl;
                    continue;
                }

                int plane = bestPlane;
                //for (int plane = 0; plane < 3; ++plane){

                auto const& calo = calo_plane[plane];
                
                if (calo.isNull())
                {
                    std::cout << "Plano " << plane << " não tem calorimetria." << std::endl;
                    continue;
                }
                //------------------------- termino de buscar o plano -------------------------------------------------------------------
                
                auto const& dEdx_v  = calo->dEdx();
                auto const& dADCdx_v = calo->dQdx();
                auto const& pitch_v = calo->TrkPitchVec();
                auto const& pos_v   = calo->XYZ();

                // create vector of e- instead of ADC units
                std::vector<float> dQdx_v(dADCdx_v.size(),0);
                for (size_t s = 0; s < dADCdx_v.size(); s++)
                {
                    dQdx_v[s] = dADCdx_v[s]*(1/_cal_area_const.at(plane));
                }

                //std::cout << "calos : " << plane << " - " << calo->PlaneID().Plane << std::endl;

                this_qlight.sliceID=sliceID;
                //varre todas as posicoes/energia depositadas
                for (size_t s = 0; s < dEdx_v.size(); s++)
                {
                    float x = pos_v[s].X();
                    float y = pos_v[s].Y();
                    float z = pos_v[s].Z();
                    if(x<0) //por hora vamos so pegar as posicoes com x positivo
                    {
                        continue;
                    }
                    double drift_time = (drift_length - abs(x))/(drift_speed); //
                    double atten_corr = std::exp(drift_time/electronlife); //

                    float pitch;
                    float dE,dQ;
                    float nphotons;

                    //NESTA PARTE O CODIGO DO SBND SEPARA EM 2 PARTES (VALORES NORMAIS E ESTRANHOS)
                    if(true)//valores normais (depois tenho que fazer o outro caso)
                    {
                        pitch = (s < pitch_v.size()) ? pitch_v[s] : -1;
                        dQ = dQdx_v[s] * pitch * atten_corr; // corigido pelo drift
                        dE = dEdx_v[s] * pitch; // talvez precise corrigir pelo drfit, de uma olhada na fcl de reconstrucao depois ...
                        nphotons = dE/(W_LAr*1e-6) - dQ;
                        
                    }
                    else
                    {
                        //aqui depois colocamos os valores estranhos
                    }

                    this_qlight.push_back(QPoint(x,y,z,nphotons));
                    

                    //mf::LogInfo("DuneSimpleSliceReader")<< "   Calo step (plane " << plane << "): dEdx=" << dEdx<< " pitch=" << pitch<< " pos=(" << x << "," << y << "," << z << ")";
                    /* 
                    v_dEdx.push_back(dEdx);
                    v_pitch.push_back(pitch);
                    v_x.push_back(x);
                    v_y.push_back(y);
                    v_z.push_back(z); */
                }
                //}
                //fTree->Fill();

            }
        }
        QLigths.push_back(this_qlight);
    }
    return QLigths;
}

void GetPandoraSliceAnalyser::getFlashs(art::Event const& e)
{
    //aqui vamos pegar os flashs
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    auto FlashHandle = e.getHandle< std::vector< recob::OpFlash >> (fFlashLabel);
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
    //numero de flashs
    int number_flashs = flashlist.size();
    std::cout << "N flashs " << number_flashs << std::endl;
}




DEFINE_ART_MODULE(GetPandoraSliceAnalyser)
