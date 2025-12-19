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

class MyFlashMatching : public art::EDAnalyzer
 {
    public:
        explicit MyFlashMatching(fhicl::ParameterSet const& p);

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
        TTree* fTreeF;
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

        art::ServiceHandle<geo::Geometry> geo;
        
        double drift_length;
        double drift_speed;
        double electronlife;
        double W_LAr;

        int nOPdet;

        void returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, int const& thisId);
        std::vector<QCluster> getQClustersSlices(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersPFPs(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersTracks(art::Event const& e); // esse aqui eh por track
        std::vector<QFlash> getFlashs(art::Event const& e);

        phot::PhotonVisibilityService const* fPVS;
        phot::SemiAnalyticalModel const* fSAM;

        bool norm=false;
        std::string ClusterType; //Slice,PFP,Track
        std::string DetectorZone; //Positive,Negative,All

        int fflashID,fclusterID,fclusterTrivialID;
        double fScore,fxoffset;
        double fScoreTrivial,fxoffsetTrivial;
};


MyFlashMatching::MyFlashMatching(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
    //seta os labels para ler os produtos
    fSliceLabel = p.get<art::InputTag>("SliceLabel");  
    fPFPLabel = p.get<art::InputTag>("PFParticleLabel");
    fTrackLabel = p.get<art::InputTag>("TrackLabel");
    fCaloLabel = p.get<art::InputTag>("CalorimetryLabel");
    fFlashLabel = p.get<art::InputTag>("FlashLabel");

    norm =p.get<bool>("NormFlash");
    ClusterType=p.get<std::string>("ClusterType");
    DetectorZone=p.get<std::string>("DetectorZone");

    _cal_area_const    = p.get<std::vector<float>>("CalAreaConstants"); // PEGANDO OS VALORES PADRAO, TA CERTO??

    std::cout << "CalAreaConstants = "
          << _cal_area_const[0] << ", "
          << _cal_area_const[1] << ", "
          << _cal_area_const[2] << std::endl;


    //geometria do detector
    int nTPCs = geo->TotalNTPC();
    int nCrio = geo->Ncryostats();
    nOPdet=geo->NOpDets();
   
    std::cout << "nCryo: " << nCrio << std::endl;
    std::cout << "nTPCs: " << nTPCs << std::endl;
    std::cout << "nARAPUCAs: " << nOPdet << std::endl;

    geo::TPCID id(0,1);
    auto const& tpc = geo->TPC(id);
    this->drift_length = tpc.DriftDistance(); // ~ 360 cm 

    std::cout << "Drift Distance: " << this->drift_length << std::endl;

    // carrega o servico de visibilidade otica
    fPVS = art::ServiceHandle<phot::PhotonVisibilityService>().get();
    std::cout << "Loaded VS. NOpChannels = " << fPVS->NOpChannels() << std::endl;

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

void MyFlashMatching::beginJob()
{
     art::ServiceHandle<art::TFileService> tfs;

    fTreeF = tfs->make<TTree>("treeMatch", "Flash Match with likelihood fit using poisson and Hungarian algorithm to assing  light to charge");
    fTreeF->Branch("run", &run);
    fTreeF->Branch("event", &event);
    fTreeF->Branch("clusterID", &fclusterID);
    fTreeF->Branch("clusterID_trivial", &fclusterTrivialID);
    fTreeF->Branch("flashID", &fflashID);
    fTreeF->Branch("score", &fScore);
    fTreeF->Branch("x", &fxoffset);
    fTreeF->Branch("score_trivial", &fScoreTrivial);
    fTreeF->Branch("x_trivial", &fxoffsetTrivial);
    
}

void MyFlashMatching::analyze(art::Event const& e)
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

    std::vector<QCluster> QClusters;
    
    if(ClusterType=="Slice") QClusters = getQClustersSlices(e);
    if(ClusterType=="PFP")  QClusters = getQClustersPFPs(e);
    if(ClusterType=="Track") QClusters = getQClustersTracks(e);
    
    auto QFlashs = getFlashs(e);
    
    std::cout <<"Final number Qcluster: " << QClusters.size() << std::endl;
    std::cout <<"Final number Flashs: " << QFlashs.size() << std::endl;

    myMatch* match_operator = new myMatch(QClusters,QFlashs,drift_length,drift_speed,fPVS,fSAM,norm,DetectorZone);
    std::cout << "pudim" << std::endl;

    auto& HR = match_operator->HR;
    auto& MYScore = match_operator->MYScore;
    auto& MYOffset = match_operator->MYOffset;

    auto& qqs = QClusters;
    auto& qfs = QFlashs;
    const int nRows = std::min((int)HR.assign.size(), match_operator->Nf);

    for (int nf = 0; nf < nRows; ++nf) 
    {
        auto const& row = MYScore[nf];
        auto it = std::min_element(row.begin(), row.begin() + match_operator->Nc); // só clusters reais
        int ncTrivial = (int)std::distance(row.begin(), it);
        fclusterTrivialID = qqs[ncTrivial].objID; 
        fScoreTrivial = MYScore[nf][ncTrivial];
        fxoffsetTrivial = MYOffset[nf][ncTrivial];

        const int nc = HR.assign[nf];
        fflashID = qfs[nf].flashID;
        // se caiu em dummy
        if (nc < 0 || nc >= match_operator->Nc)
        {
            fclusterID = -1;
            fScore     = -1000;
            fxoffset   = 0.0;
        } 
        else
        {
            fclusterID = qqs[nc].objID;
            fScore     = MYScore[nf][nc];
            fxoffset   = MYOffset[nf][nc];
        }

        fTreeF->Fill();
    }


    delete match_operator;
}

void MyFlashMatching::returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, int const& thisId)
{
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
        return;
    }

    int plane = bestPlane;
    //for (int plane = 0; plane < 3; ++plane){

    auto const& calo = calo_plane[plane];
    
    if (calo.isNull())
    {
        std::cout << "Plano " << plane << " não tem calorimetria." << std::endl;
        return;
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

    //varre todas as posicoes/energia depositadas
    for (size_t s = 0; s < dEdx_v.size(); s++)
    {
        float x = pos_v[s].X();
        float y = pos_v[s].Y();
        float z = pos_v[s].Z();
        if(DetectorZone == "Positive" && x<0) //o pegar as posicoes com x positivo
        {
            continue;
        }
        else if(DetectorZone == "Negative" && x>0) // so pegar as posicoes com x negativo
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
    }
    this_qlight.objID=thisId;
}


//eu vi que para cosmics eh melhor fazer match diretamente com track/PFParticle em vez de slices ( muito quebrado para cosmics )

std::vector<QCluster> MyFlashMatching::getQClustersSlices(art::Event const& e)
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

    std::cout << "N Slices " << nSlices << std::endl;

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
                //trackID = trk->ID(); //ID do track
                //pega as info de calorimetria
                returnQCluster(this_qlight,trk,trk_to_calo,sl->ID()); 
            }
        }
        if(this_qlight.size()>0)
        {
            QLigths.push_back(this_qlight);
        }  
    }
    return QLigths;
}


std::vector<QCluster> MyFlashMatching::getQClustersPFPs(art::Event const& e)
{
   
    // pegar produtos importantes
    
    auto pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel); //PFParticles
    auto track_h = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel); // tracks
    auto calo_h = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloLabel); //calorimetria

    art::FindManyP<recob::Track> pfp_to_tracks(pfp_h, e, fTrackLabel); // pega a associacao de tracks das PFParticles
    art::FindManyP<anab::Calorimetry> trk_to_calo(track_h, e, fCaloLabel); // pega as info de calorimetria dos tracks

    std::vector<QCluster> QLigths;
    
    std::vector<art::Ptr<recob::PFParticle>> pfps;
    art::fill_ptr_vector(pfps, pfp_h);
    nPFPs = pfps.size();

    std::cout << "N PFPs: " << nPFPs << std::endl; 

    //varre todos os pfps deste slice
    for (auto const& pfp : pfps) 
    {
        QCluster this_qlight;
        pfpID = pfp->Self(); // ID do PFParticle
        // pegar tracks associados ao PFP
        std::vector<art::Ptr<recob::Track>> tracks = pfp_to_tracks.at(pfp.key());
        nTracksPFP = tracks.size();

        //varre todos tracks deste PFParticle
        for (auto const& trk : tracks)
        {
            returnQCluster(this_qlight,trk,trk_to_calo, pfp->Self()); 
        }
        if(this_qlight.size()>0)
        {
            QLigths.push_back(this_qlight);
        }  
    }
    return QLigths;
}

std::vector<QCluster> MyFlashMatching::getQClustersTracks(art::Event const& e)
{ 
    // pegar produtos importantes
    auto track_h = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel); // tracks
    auto calo_h = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloLabel); //calorimetria

    art::FindManyP<anab::Calorimetry> trk_to_calo(track_h, e, fCaloLabel); // pega as info de calorimetria dos tracks

    std::vector<QCluster> QLigths;
    //varre os tracks

    std::vector<art::Ptr<recob::Track>> tracks;
    art::fill_ptr_vector(tracks, track_h);
    auto nTracks = tracks.size();

    std::cout << "N TRACKS: " << nTracks << std::endl; 
    //varre todos tracks deste PFParticle
    for (auto const& trk : tracks)
    {
        QCluster this_qlight;
        returnQCluster(this_qlight,trk,trk_to_calo,trk->ID()); 
       
        if(this_qlight.size()>0)
        {
            QLigths.push_back(this_qlight);
        }  
    }
    return QLigths;
}

std::vector<QFlash> MyFlashMatching::getFlashs(art::Event const& e)
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
      mf::LogWarning("MyFlashMatching") << "Cannot load any flashes. Failing";
      return {};
    }
    art::FindManyP<recob::OpHit> OpHits_from_Flashs(FlashHandle, e, fFlashLabel);

    //numero de flashs
    int number_flashs = flashlist.size();
    std::cout << "N flashs " << number_flashs << std::endl;

    std::vector<QFlash> QFlashs;
    //varrer os flashs
    for(int i=0;i<number_flashs;i++)
    {
        QFlash this_qflash;
        
        auto const& flash = flashlist[i];
        auto hits = OpHits_from_Flashs.at(flash.key());
        int number_hits=hits.size();
        this_qflash.flashID=flash.key();
        this_qflash.PE_CH.resize(this->nOPdet);

        for(int j=0; j<this->nOPdet ; j++)
        {
            this_qflash.PE_CH[j]=0.0;
        }
        for(int j=0;j<number_hits;j++)
        {
            int this_ch=hits[j]->OpChannel();
            if(DetectorZone=="Positive" && this_ch>=80) continue;
            if(DetectorZone=="Negative" && this_ch<80) continue;
            this_qflash.PE_CH[this_ch]+=hits[j]->PE();
        }

        auto flash1 = 0.0;
        auto flash2 = 0.0;

        for(int i=0;i<this->nOPdet;i++)
        {
            if(i<80)
            {
                flash1+=this_qflash.PE_CH[i];
            }
            else
            {
                flash2+=this_qflash.PE_CH[i];
            }
        }
        if(norm) this_qflash.norm_this_flash();
        /* std::cout <<"(x,y,z) : (" << x << " , " << y <<" , " << z << ")"<<std::endl;
        std::cout <<"flash1 : "<< flash1 << "   flash 2 : " << flash2 <<std::endl;
        std::cout << "Press ENTER." << std::endl;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); */
        if(true)//flash1>=flash2) // apenas APA 1 E 2
        {
            this_qflash.y = flash->YCenter();
            this_qflash.z = flash->ZCenter();
            this_qflash.y_err = flash->YWidth();
            this_qflash.z_err = flash->ZWidth();
            this_qflash.time = flash->Time(); 
            this_qflash.time_err = flash->TimeWidth(); 
            QFlashs.push_back(this_qflash);
        }
     
    }
    
    return QFlashs;
}


DEFINE_ART_MODULE(MyFlashMatching)
