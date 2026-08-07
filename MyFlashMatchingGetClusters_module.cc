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
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include "dunecore/DuneObj/OpDetDivRec.h" 
#include "duneopdet/OpticalDetector/OpFlashSort.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "art/Utilities/make_tool.h"

#include "TTree.h"

#include <vector>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <functional>

#include "mydune/utils/my_utils.hh"
#include "mydune/utils/MyMatch.hh"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


class MyFlashMatchingGetClusters : public art::EDAnalyzer
 {
    public:
        explicit MyFlashMatchingGetClusters(fhicl::ParameterSet const& p);

        void beginJob() override;
        void analyze(art::Event const& e) override;

    private:
        //declarao dos nomes dos produtos
        art::InputTag fSliceLabel; 
        art::InputTag fPFPLabel;
        art::InputTag fTrackLabel;
        art::InputTag fShowerLabel;
        art::InputTag fCaloLabel;
        art::InputTag fCaloShowerLabel;
        art::InputTag fFlashLabel;
        art::InputTag fHitLabel;
        art::InputTag fTriggerLabel;

        //este aqui eh para converter ADC em e-
        std::vector<float> _cal_area_const;

        //declaracao da TTree
        TTree* fTreeF;
        TTree* fTreeT;
        TTree* fTreeFT;
        int run, event;
        int sliceID;
        int nSlices;
        int pfpID;
        int nPFPs;
        int nTracksPFP;
        int nShowersPFP;
        int trackID;
        std::vector<float> v_dEdx;
        std::vector<float> v_dQdx;
        std::vector<float> v_pitch;
        std::vector<float> v_x;
        std::vector<float> v_y;
        std::vector<float> v_z;

        art::ServiceHandle<geo::Geometry> geo;

        QFlash fFlash;
        QCluster fCluster;
        QFlash fFlash_fit;
        QCluster fCluster_fit;

        int ftrackID,fflashID,ftrackType;

        int nOPdet;
        double fvis,fvislight;

        double drift_length;
        double drift_speed;
        double electronlife;
        double W_LAr;
        double density;
        double Efield;

        double gain_factor=1.0;
        double vis_map_factor=1.0;

        void returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, int const& thisId, int const& typeObj );
        //void returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<recob::Hit> const& shw_to_hit, art::FindManyP<recob::SpacePoint> const& hit_to_sp, int const& thisId, int const& typeObj );
        void returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<anab::Calorimetry> const& shw_to_calo, int const& thisId,int const& typeObj );

        std::vector<QCluster> getQClustersSlices(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersPFPs(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersTracks(art::Event const& e); // esse aqui eh por track
        std::vector<QFlash> getFlashs(art::Event const& e);

        bool getTrackDirFromCaloPoint(art::Ptr<recob::Track> const& trk,double x,double y,double z,float& dirx,float& diry,float& dirz);


        phot::PhotonVisibilityService const* fPVS;
        phot::SemiAnalyticalModel const* fSAM;

        std::vector<double> fPDEVector;
        std::vector<double> fXTalkVector;
        std::vector<int> fCHActiveVector;

        std::vector<double> xch,ych,zch,dch;

        std::vector<double> LYCH;

        bool norm=false;
        bool getShowers = false;
        bool fitMode = true;
        bool isMC = true;
        bool useSCE = false;
        std::string ClusterType; //Slice,PFP,Track
        std::string DetectorZone; //Positive,Negative,All

        double limitMinFlash = 0 ;
        double limitMinFlashSide = 0;
        double trackLength = 0.0;

        double ftriggerTime = 0.0;

};

MyFlashMatchingGetClusters::MyFlashMatchingGetClusters(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
    //seta os labels para ler os produtos
    fSliceLabel = p.get<art::InputTag>("SliceLabel");  
    fPFPLabel = p.get<art::InputTag>("PFParticleLabel");
    fTrackLabel = p.get<art::InputTag>("TrackLabel");
    fShowerLabel = p.get<art::InputTag>("ShowerLabel");
    fCaloLabel = p.get<art::InputTag>("CalorimetryLabel");
    fCaloShowerLabel = p.get<art::InputTag>("CalorimetryShowerLabel");
    fFlashLabel = p.get<art::InputTag>("FlashLabel");
    fHitLabel = p.get<art::InputTag>("HitLabel");

    norm = p.get<bool>("NormFlash",false);
    if(norm)
        std::cout << "Normalizing the flash!!! " << std::endl;
    else
        std::cout << "Not Normalizing the flash!!! " << std::endl;

    isMC = p.get<bool>("isMC",true);
    if(isMC)
        std::cout << "Analysing MC!!! " << std::endl;
    else
        std::cout << "Analysing Real Data!!! " << std::endl;

    useSCE = p.get<bool>("useSCE",false);
    if(useSCE)
        std::cout << "Using SCE corrections " << std::endl;
    else
        std::cout << " NOTTTT Using SCE corrections " << std::endl;

    ClusterType=p.get<std::string>("ClusterType","Track");
    DetectorZone=p.get<std::string>("DetectorZone","All");;

    limitMinFlash = p.get<double>("limitMinFlash",0.0);
    limitMinFlashSide = p.get<double>("limitMinFlashSide",0.0);

    trackLength = p.get<double>("trackLengthMin",0.0);

    gain_factor = p.get<double>("gain_factor",1.0);
    vis_map_factor = p.get<double>("vis_map_factor",1.0);

    std::cout << "gain_factor: " << gain_factor << std::endl;
    std::cout << "vis_map_factor: " << vis_map_factor << std::endl;

    std::cout << "Track Length min: " << trackLength << std::endl;
    std::cout << "limitMinFlash: " << limitMinFlash << std::endl;
    std::cout << "limitMinFlashSide: " << limitMinFlashSide << std::endl;

    getShowers = p.get<bool>("getShowers",false);
    if(getShowers)
        std::cout << "Getting Showers!!! " << std::endl;
    else
        std::cout << "Not Getting Showers!!! " << std::endl;

    //geometria do detector
    int nTPCs = geo->TotalNTPC();
    int nCrio = geo->Ncryostats();
    nOPdet=geo->NOpDets();
   
    std::cout << "nCryo: " << nCrio << std::endl;
    std::cout << "nTPCs: " << nTPCs << std::endl;
    std::cout << "nARAPUCAs: " << nOPdet << std::endl;

    fTriggerLabel = p.get<std::string>("TriggerTag","daq:trigger:pdhdkeepupstage1");

    ftriggerTime = p.get<double>("triggerTime",0.0);
    std::cout << ftriggerTime << " : Trigger Time loaded (only used when MONTE CARLO)" << std::endl;

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

    fPDEVector = p.get<std::vector<double>>("PDEvector",std::vector<double>(nOPdet,0.03));
    fXTalkVector = p.get<std::vector<double>>("XTalkvector",std::vector<double>(nOPdet,0.09));
    fCHActiveVector = p.get<std::vector<int>>("CHActive",std::vector<int>(nOPdet,1));

    std::cout << "PDE loaded " << fPDEVector[0] << std::endl;
    std::cout << "Xtalk loaded " << fXTalkVector[0] << std::endl;
    std::cout << "Active Vector loaded " << fCHActiveVector[0] << std::endl;

    _cal_area_const    = p.get<std::vector<float>>("CalAreaConstants"); // PEGANDO OS VALORES PADRAO, TA CERTO??

    std::cout << "CalAreaConstants = "
          << _cal_area_const[0] << ", "
          << _cal_area_const[1] << ", "
          << _cal_area_const[2] << std::endl;


    xch=std::vector<double>(nOPdet,0.0);
    ych=std::vector<double>(nOPdet,0.0);
    zch=std::vector<double>(nOPdet,0.0);
    dch=std::vector<double>(nOPdet,0.0);
    for(int channel=0;channel<nOPdet;++channel)
    {   
        auto opDet = geo->OpDetGeoFromOpDet(channel);
        auto const c = opDet.GetCenter(); 
        xch[channel]=c.X();
        ych[channel]=c.Y();
        zch[channel]=c.Z();
    }  
    
}

void MyFlashMatchingGetClusters::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTreeF = tfs->make<TTree>("treeF", "");
    
    fTreeF->Branch("run", &run);
    fTreeF->Branch("event", &event);
    fTreeF->Branch("flash", &fFlash);

    fTreeT = tfs->make<TTree>("treeT", "");

    fTreeT->Branch("run", &run);
    fTreeT->Branch("event", &event);
    fTreeT->Branch("cluster", &fCluster);

    fTreeFT = tfs->make<TTree>("treeFT", "");

    fTreeFT->Branch("run", &run);
    fTreeFT->Branch("event", &event);
    fTreeFT->Branch("flashID", &fflashID);
    fTreeFT->Branch("trackID", &ftrackID);
    fTreeFT->Branch("trackType", &ftrackType);
    fTreeFT->Branch("vis", &fvis);
    fTreeFT->Branch("vis_light", &fvislight);
    fTreeFT->Branch("flash", &fFlash_fit);
    fTreeFT->Branch("cluster", &fCluster_fit);
    fTreeFT->Branch("LYCH", &LYCH);
    fTreeFT->Branch("dCH", &dch);
    
}

void MyFlashMatchingGetClusters::analyze(art::Event const& e)
{
    run    = e.id().run();
    event  = e.id().event();

    auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);
    art::ServiceHandle<sim::LArG4Parameters const> g4param;

    drift_speed = det_prop.DriftVelocity(); //~ 0.15 cm/us
    electronlife = det_prop.ElectronLifetime(); //35000 us
    W_LAr=g4param->Wph(); //19.5 eV
    this->density = det_prop.Density();
    this->Efield = det_prop.Efield();

    std::cout << "electron vd := " << drift_speed << std::endl;
    std::cout << "electron lifetime:=  " << electronlife << std::endl;
    std::cout << "W LAr :=  " << W_LAr << std::endl;
    std::cout << "Density :=  " << density << std::endl;
    std::cout << "Efield :=  " << Efield << std::endl;

    std::vector<QCluster> QClusters;
    
    if(ClusterType=="Slice") QClusters = getQClustersSlices(e);
    if(ClusterType=="PFP")  QClusters = getQClustersPFPs(e);
    if(ClusterType=="Track") QClusters = getQClustersTracks(e);
    
    auto QFlashs = getFlashs(e);

    if (QClusters.empty() || QFlashs.empty()) 
    {
        mf::LogWarning("MyFlashMatchingGetClusters")
            << "Skipping event " << e.id()
            << " because clusters or flashes are missing/empty.";
        return;
    }

    const int Nc = QClusters.size();
    const int Nf = QFlashs.size();

    for (int nf = 0; nf < Nf; ++nf) 
    {
        fFlash = QFlashs[nf];
        fTreeF->Fill();
    }
   
    for (int nc = 0; nc < Nc; ++nc) 
    {
        fCluster = QClusters[nc];
        fTreeT->Fill();
    }

    auto& qqs = QClusters;
    auto& qfs = QFlashs;
    myMatch* match_operator = new myMatch(drift_length,drift_speed,electronlife,density,Efield,fPVS,fSAM,fPDEVector,fXTalkVector,fCHActiveVector,useSCE,sce);
    match_operator->vis_map_factor = vis_map_factor;
    for (int nf = 0; nf < Nf; ++nf) 
    {
        match_operator->flash_actual = qfs[nf]; 
        fflashID = qfs[nf].flashID;
        for (int nc = 0; nc < Nc; ++nc) 
        {
            match_operator->cluster_actual = qqs[nc]; 
            ftrackID = qqs[nc].objID;
            ftrackType =  qqs[nc].type;
           if(match_operator->checkPossibility(&(match_operator->cluster_actual),&(match_operator->flash_actual),50))
           {
                double deltax_T0 =  match_operator->flash_actual.time*drift_speed;
                match_operator->ChargeHypothesis(deltax_T0);
                fFlash_fit=match_operator->flash_fit;
                fCluster_fit = match_operator->cluster_fit;
                fCluster_fit.Charge = match_operator->save_charge;
                fCluster_fit.Length = fCluster_fit.calcLength();
                fCluster_fit.Energy = fCluster_fit.TotalEnergy();
                fFlash_fit.Light = fFlash_fit.TotalLight();
                fvis = match_operator->returnVisEff();
                fvislight = match_operator->returnVisEffLight();
                LYCH =  match_operator->returnVisEffCh();
                dch = match_operator->returndCh(xch,ych,zch);
                fTreeFT->Fill();
           }

        }
    }
    delete match_operator;
}

void MyFlashMatchingGetClusters::returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, int const& thisId,int const& typeObj )
{
    std::vector<art::Ptr<anab::Calorimetry>> calos = trk_to_calo.at(trk.key());
   

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
    int ob_APA=0;

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
    ob_APA= calo->PlaneID().TPC;
    //std::cout << ob_APA << std::endl;
    
    if(ob_APA != 1 && ob_APA != 2 && ob_APA != 5 && ob_APA != 6) return;

    //if(ob_APA==1) return; // to tirando o apa1, pode colocar dpeois se quiser

    //------------------------- termino de buscar o plano -------------------------------------------------------------------
    
    auto const& dEdx_v  = calo->dEdx();
    auto const& dADCdx_v = calo->dQdx();
    auto const& pitch_v = calo->TrkPitchVec();
    auto const& pos_v   = calo->XYZ();
    auto const& indexpoints = calo->TpIndices();

    // create vector of e- instead of ADC units
    std::vector<float> dQdx_v(dADCdx_v.size(),0);
    for (size_t s = 0; s < dADCdx_v.size(); s++)
    {
        dQdx_v[s] = dADCdx_v[s]*(1/_cal_area_const.at(plane))*gain_factor;
    }

    //std::cout << "calos : " << plane << " - " << calo->PlaneID().Plane << std::endl;
    //varre todas as posicoes/energia depositadas
    
    for (size_t s = 0; s < dEdx_v.size(); s++)
    {
        float x = pos_v[s].X();
        float y = pos_v[s].Y();
        float z = pos_v[s].Z();

        bool has_dir = false;
        float dirx = -10;
        float diry = -10;
        float dirz = -10;
        // 1) tenta o índice oficial da calorimetria
        if (s < indexpoints.size())
        {
            int indp = indexpoints[s];

            if (indp >= 0 && trk->HasValidPoint(indp))
            {
                auto const& dir_v = trk->DirectionAtPoint(indp);

                double norm = std::sqrt(
                    dir_v.X()*dir_v.X() +
                    dir_v.Y()*dir_v.Y() +
                    dir_v.Z()*dir_v.Z()
                );

                if (norm > 0.0)
                {
                    dirx = dir_v.X() / norm;
                    diry = dir_v.Y() / norm;
                    dirz = dir_v.Z() / norm;
                    has_dir = true;
                }
            }
        }

        // 2) fallback: usa o ponto do calo e outro ponto válido da track
        if (!has_dir)
        {
            has_dir = getTrackDirFromCaloPoint(trk, x, y, z, dirx, diry, dirz);
        }

        if (!has_dir) //ultimo fallback
        {
        auto const& sd = trk->StartDirection();

        double norm = std::sqrt(
            sd.X()*sd.X() +
            sd.Y()*sd.Y() +
            sd.Z()*sd.Z()
        );

        if (norm > 1.0e-12)
        {
            dirx = sd.X() / norm;
            diry = sd.Y() / norm;
            dirz = sd.Z() / norm;
            has_dir = true;
        }
        }

        if(DetectorZone == "Positive" && x<0) //o pegar as posicoes com x positivo
        {
            continue;
        }
        else if(DetectorZone == "Negative" && x>0) // so pegar as posicoes com x negativo
        {
            continue;
        }  
        //double drift_time = (drift_length - abs(x))/(drift_speed); //
        //double atten_corr = std::exp(drift_time/electronlife); //

        float pitch;
        float dQ;//, dE;
        //float nphotons;

        //NESTA PARTE O CODIGO DO SBND SEPARA EM 2 PARTES (VALORES NORMAIS E ESTRANHOS)
        if(true)//valores normais (depois tenho que fazer o outro caso)
        {
            pitch = (s < pitch_v.size()) ? pitch_v[s] : -1;
            dQ = dQdx_v[s]; // * pitch * atten_corr; // corigido pelo drift
            /* dE = dEdx_v[s] * pitch; // talvez precise corrigir pelo drfit, de uma olhada na fcl de reconstrucao depois ...
            nphotons = dE/(W_LAr*1e-6) - dQ;
            nphotons = std::max(0.0f, nphotons); */    
        }
        else
        {
            //aqui depois colocamos os valores estranhos
        }
        this_qlight.push_back(QPoint(x,y,z,dQ,dirx,diry,dirz,pitch,ob_APA));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
}

void MyFlashMatchingGetClusters::returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<anab::Calorimetry> const& shw_to_calo, int const& thisId,int const& typeObj )
{
    std::vector<art::Ptr<anab::Calorimetry>> calos = shw_to_calo.at(shw.key());
    //vetores para salvar as informacoes

    //No codigo original tem uma secao de verificar o a posicao do track e se cruzou o plano de fios, e se o track eh uncontained
    //------------------------------------------------------------------------------------
    //############### BLA BLA BLA BLA BLA BLA ############################
    //------------------------------------------------------------------------------------

    
    // ---------------------------------------------------------------------------------------------
    //sao 3 planos --> Determinar o melhor plano para fazer a associaco com os flashs

    std::vector<art::Ptr<anab::Calorimetry>> calo_plane(3, art::Ptr<anab::Calorimetry>());
    //std::cout << calos.size() <<std::endl;
    for (auto const& calo : calos)
    {
        int plane = calo->PlaneID().Plane;

        if (plane < 0 || plane > 2)
        {
            continue;
        }

        calo_plane[plane] = calo;
    }
    
    int bestPlane = -1;
    size_t maxSize = 0;
    int ob_APA=0;

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
    ob_APA= calo->PlaneID().TPC;
    //std::cout << ob_APA << std::endl;
    //std::cout << "#NHITSCALOR_SHOWERS: " <<  maxSize << " --- " << ob_APA << std::endl ;
    if(ob_APA != 1 && ob_APA != 2 && ob_APA != 5 && ob_APA != 6) return;
    //if(ob_APA==1) return; // to tirando o apa1, pode colocar dpeois se quiser
    //std::cout << "Passou "<<std::endl << "------------------" <<std::endl;
    
    //------------------------- termino de buscar o plano -------------------------------------------------------------------
    
    auto const& dEdx_v  = calo->dEdx();
    auto const& dADCdx_v = calo->dQdx();
    auto const& pitch_v = calo->TrkPitchVec();
    auto const& pos_v   = calo->XYZ();

    // create vector of e- instead of ADC units
    std::vector<float> dQdx_v(dADCdx_v.size(),0);
    for (size_t s = 0; s < dADCdx_v.size(); s++)
    {
        dQdx_v[s] = dADCdx_v[s]*(1/_cal_area_const.at(plane))*gain_factor;
    }

    //varre todas as posicoes/energia depositadas
    for (size_t s = 0; s < dEdx_v.size(); s++)
    {
        float x = pos_v[s].X();
        float y = pos_v[s].Y();
        float z = pos_v[s].Z();

        const auto& dir = shw->Direction();
        float dirx = dir.X();
        float diry = dir.Y();
        float dirz = dir.Z();
        
        if(DetectorZone == "Positive" && x<0) //o pegar as posicoes com x positivo
        {
            continue;
        }
        else if(DetectorZone == "Negative" && x>0) // so pegar as posicoes com x negativo
        {
            continue;
        }  
        //double drift_time = (drift_length - abs(x))/(drift_speed); //
        //double atten_corr = std::exp(drift_time/electronlife); //

        float pitch;
        float dQ;//, dE;
        //float nphotons;

        //NESTA PARTE O CODIGO DO SBND SEPARA EM 2 PARTES (VALORES NORMAIS E ESTRANHOS)
        if(true)//valores normais (depois tenho que fazer o outro caso)
        {
            pitch = (s < pitch_v.size()) ? pitch_v[s] : -1;
            dQ = dQdx_v[s];// * pitch * atten_corr; // corigido pelo drift
            /* dE = dEdx_v[s] * pitch; // talvez precise corrigir pelo drfit, de uma olhada na fcl de reconstrucao depois ...
            nphotons = dE/(W_LAr*1e-6) - dQ;
            nphotons = std::max(0.0f, nphotons); */
            
        }
        else
        {
            //aqui depois colocamos os valores estranhos
        }
        this_qlight.push_back(QPoint(x,y,z,dQ,dirx,diry,dirz,pitch,ob_APA));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
        
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
}

//eu vi que para cosmics eh melhor fazer match diretamente com track/PFParticle em vez de slices ( muito quebrado para cosmics )

std::vector<QCluster> MyFlashMatchingGetClusters::getQClustersSlices(art::Event const& e)
{
    auto slice_h = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);
    if (!slice_h) {
        mf::LogWarning("MyFlashMatchingGetClusters")
            << "Slice product not found: " << fSliceLabel
            << " ; skipping event " << e.id();
        return {};
    }
    auto pfp_h   = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel);
    auto track_h = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
    auto calo_h  = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloLabel);

    art::FindManyP<recob::PFParticle>  slice_to_pfps(slice_h, e, fSliceLabel);
    art::FindManyP<recob::Track>       pfp_to_tracks(pfp_h, e, fTrackLabel);
    art::FindManyP<anab::Calorimetry>  trk_to_calo(track_h, e, fCaloLabel);

    art::Handle<std::vector<recob::Shower>> shower_h;
    std::optional<art::FindManyP<recob::Shower>> pfp_to_showers;
    std::optional<art::FindManyP<anab::Calorimetry>> shw_to_calo;
    bool haveShowers = false;

    if (getShowers)
    {
        shower_h = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
        if (shower_h.isValid())
        {
            pfp_to_showers.emplace(pfp_h, e, fShowerLabel);
            shw_to_calo.emplace(shower_h, e, fCaloShowerLabel);
            haveShowers = true;
        }
        else
        {
            mf::LogWarning("MyFlashMatchingGetClusters")
                << "getShowers=true, but Shower handle is invalid for label " << fShowerLabel;
        }
    }

    nSlices = slice_h->size();

    std::vector<art::Ptr<recob::Slice>> slices;
    art::fill_ptr_vector(slices, slice_h);

    std::vector<art::Ptr<recob::PFParticle>> allpfps;
    art::fill_ptr_vector(allpfps, pfp_h);

    std::unordered_map<int, art::Ptr<recob::PFParticle>> pfpMap;
    for (auto const& pfp : allpfps)
    {
        if (!pfp.isNull()) pfpMap[pfp->Self()] = pfp;
    }

    std::vector<QCluster> QLights;

    std::cout << "N Slices " << nSlices << std::endl;

    for (auto const& sl : slices)
    {
        QCluster this_qlight;
        double totalLength = 0.0;

        auto seed_pfps = slice_to_pfps.at(sl.key());

        sliceID = sl->ID();
        nPFPs   = seed_pfps.size();

        std::unordered_set<int> visited;

        std::function<void(const art::Ptr<recob::PFParticle>&)> visitPFP;
        int Nt=0;
        int Ns=0;
        visitPFP = [&](const art::Ptr<recob::PFParticle>& pfp)
        {
            if (pfp.isNull()) return;

            const int self = pfp->Self();
            if (!visited.insert(self).second) return; // evita contar duas vezes

            pfpID = self;

            // Tracks deste PFP
            auto tracks = pfp_to_tracks.at(pfp.key());
            nTracksPFP = tracks.size();

            for (auto const& trk : tracks)
            {
                if (trk.isNull()) continue;
                if (trk->Length() < 0.0) continue;

                size_t nBefore = this_qlight.size();
                returnQCluster(this_qlight, trk, trk_to_calo, sl->ID(), 3);
                
                if (this_qlight.size() > nBefore)
                {
                    totalLength += trk->Length();
                    Nt +=1;
                }
                
            }

            // Showers deste PFP
            if (haveShowers && pfp_to_showers && shw_to_calo)
            {
                auto showers = pfp_to_showers->at(pfp.key());
                nShowersPFP = showers.size();

                for (auto const& shw : showers)
                {
                    if (shw.isNull()) continue;
                    if (shw->Length() < 0.0) continue;

                    size_t nBefore = this_qlight.size();
                    returnQClusterShower(this_qlight, shw, *shw_to_calo, sl->ID(), 3);

                    if (this_qlight.size() > nBefore)
                    {
                        totalLength += shw->Length();
                        Ns +=1;
                    }
                   
                }
            }

            // Desce para as filhas
            for (int dauID : pfp->Daughters())
            {
                auto it = pfpMap.find(dauID);
                if (it != pfpMap.end())
                {
                    visitPFP(it->second);
                }
            }
        };

        // começa dos PFPs associados ao slice
        for (auto const& pfp : seed_pfps)
        {
            visitPFP(pfp);
        }

        if (!this_qlight.empty())
        {
            this_qlight.objID  = sl->ID();
            this_qlight.type   = 3;     // Slice
            this_qlight.APA    = -1;    // slice pode misturar APAs
            this_qlight.Length = totalLength;
            this_qlight.Nt = Nt;
            this_qlight.Ns = Ns;

            if (totalLength >= trackLength)
            {
                this_qlight.Charge = this_qlight.TotalCharge();
                QLights.push_back(this_qlight);
            }
        }
    }

    return QLights;
}


std::vector<QCluster> MyFlashMatchingGetClusters::getQClustersPFPs(art::Event const& e)
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
            if(trk->Length() >= trackLength)
            {
                returnQCluster(this_qlight,trk,trk_to_calo, pfp->Self(),2); 
            }
        }
        if(this_qlight.size()>0)
        {
            QLigths.push_back(this_qlight);
        }  
    }
    return QLigths;
}

std::vector<QCluster> MyFlashMatchingGetClusters::getQClustersTracks(art::Event const& e)
{ 
    // pegar produtos importantes
    auto track_h = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel); // tracks
    auto calo_h = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloLabel); //calorimetria

    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const det_prop   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

    art::FindManyP<anab::Calorimetry> trk_to_calo(track_h, e, fCaloLabel); // pega as info de calorimetria dos tracks

    std::vector<QCluster> QLigths;
    //varre os tracks

    std::vector<art::Ptr<recob::Track>> tracks;
    art::fill_ptr_vector(tracks, track_h);
    auto nTracks = tracks.size();

    art::FindManyP<recob::Hit> fmHits(track_h, e, fTrackLabel);
    if (!fmHits.isValid()) 
    {
      mf::LogWarning("GetMyWireData") << "No Track<->Hit assns for " << fTrackLabel << std::endl;  
    }

    std::cout << "N TRACKS: " << nTracks << std::endl; 
    //varre todos tracks deste PFParticle
    for (auto const& trk : tracks)
    {
        QCluster this_qlight;
        auto hits = fmHits.at(trk.key());
        double xmint =  100000;
        double xmaxt = -100000;
       
        geo::WireID wid;
        for (auto const& h : hits) 
        {
            if (h.isNull()) continue;

            wid = h->WireID();

            // X do TPC *no mesmo frame do reco*
            auto TPC = wid.TPC;
            if(TPC==0 || TPC==3 || TPC==4 || TPC==7) continue;
            double x_hit = det_prop.ConvertTicksToX(h->PeakTime(), wid.Plane, wid.TPC, wid.Cryostat);

            xmint = std::min(xmint, x_hit);
            xmaxt = std::max(xmaxt, x_hit);
        }

        if(xmint==100000) continue;
        if(xmaxt==-100000) continue;

        //std::cout << det_prop.ConvertTicksToX(0, wid.Plane, 2, wid.Cryostat) << "   " << det_prop.ConvertTicksToX(5999, wid.Plane, 2, wid.Cryostat) << std::endl;
        //std::cout << det_prop.ConvertTicksToX(0, wid.Plane, 1, wid.Cryostat) << "   " << det_prop.ConvertTicksToX(5999, wid.Plane, 1, wid.Cryostat) << std::endl;

        if(trk->Length() >= trackLength)
        {
            returnQCluster(this_qlight,trk,trk_to_calo,trk->ID(),0);
            this_qlight.Length = trk->Length();
            this_qlight.Charge = this_qlight.TotalCharge();
            this_qlight.Nt = 1;
            this_qlight.Ns = 0;
        }
        
        if(this_qlight.size()>0)
        {
            QLigths.push_back(this_qlight);
        }
        
    }

    if(getShowers)
    {
        auto shower_h = e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel);
        auto calo_sh_h = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloShowerLabel); // shower
        art::FindManyP<anab::Calorimetry> shw_to_calo(shower_h, e, fCaloShowerLabel);
        //auto hit_h = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
        //art::FindManyP<recob::Hit> shower_to_hits (shower_h, e, fShowerLabel);
        //art::FindManyP<recob::SpacePoint> hit_to_sps(hit_h, e, fHitLabel);
        

        art::FindManyP<recob::Hit> fmHits(shower_h, e, fShowerLabel);
        if (!fmHits.isValid()) 
        {
            mf::LogWarning("GetMyWireData") << "No Track<->Hit assns for " << fShowerLabel << std::endl;  
        }
    
        std::vector<art::Ptr<recob::Shower>> showers;
        art::fill_ptr_vector(showers, shower_h);
        auto nShowers = showers.size();

        std::cout << "N SHOWERS: " << nShowers << std::endl; 
        for (auto const& shw : showers)
        {   
            double xmint =  100000;
            double xmaxt = -100000;
           // int nHits=0;
            auto hits = fmHits.at(shw.key());
            for (auto const& h : hits) 
            {
                if (h.isNull()) continue;
                auto const wid = h->WireID();

                // X do TPC *no mesmo frame do reco*
                auto TPC = wid.TPC;
                if(TPC==0 || TPC==3 || TPC==4 || TPC==7) continue;
                double x_hit = det_prop.ConvertTicksToX(h->PeakTime(), wid.Plane, wid.TPC, wid.Cryostat);

                xmint = std::min(xmint, x_hit);
                xmaxt = std::max(xmaxt, x_hit);
                //nHits++;
            }

            //checa se tem pontos no volume ativo
            if(xmint==100000) continue;
            if(xmaxt==-100000) continue;
               
            //std::cout << "#NHITS_SHOWERS: " <<  nHits << "  " << shw->Length() << std::endl;
            QCluster this_qlight;
            if(shw->Length() >= trackLength)
            {
                //std::cout << "tenho track length o suficiente" << std::endl;
                //returnQClusterShower(this_qlight,shw,shower_to_hits,hit_to_sps,shw->ID(),1); 
                returnQClusterShower(this_qlight,shw,shw_to_calo,shw->ID(),1); 
                this_qlight.Length = shw->Length();
                this_qlight.Charge = this_qlight.TotalCharge();
                this_qlight.Nt = 0;
                this_qlight.Ns = 1;
            }
            
            if(this_qlight.size()>0)
            {
                QLigths.push_back(this_qlight);
            } 
        }
    }
    
    return QLigths;
}

std::vector<QFlash> MyFlashMatchingGetClusters::getFlashs(art::Event const& e)
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
      mf::LogWarning("MyFlashMatchingGetClusters") << "Cannot load any flashes. Failing";
      return {};
    }
    art::FindManyP<recob::OpHit> OpHits_from_Flashs(FlashHandle, e, fFlashLabel);

    //numero de flashs
    int number_flashs = flashlist.size();
    std::cout << "N flashs " << number_flashs << std::endl;

    std::vector<QFlash> QFlashs;

    //se nao for MC temos que ler o triiger time
    Long64_t triggerTime=0;
    if(!isMC)
    {
        auto hTrigger = e.getHandle<raw::RDTimeStamp>(fTriggerLabel);
        if(!hTrigger)
        {
            std::cout << "No Trigger Time found" << std::endl; 
            return {};
        }
        triggerTime =  static_cast<Long64_t>((*hTrigger).GetTimeStamp());
    }


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
            this_qflash.PE_CH[this_ch]+=hits[j]->PE()*fCHActiveVector[this_ch];
        }

        auto flash1 = 0.0;
        auto flash2 = 0.0;
      
        int NCh_active=0;
        for(int i=0;i<this->nOPdet;i++)
        {
            if(i<80)
            {
                flash1+=this_qflash.PE_CH[i];
                if(DetectorZone=="Negative") this_qflash.PE_CH[i]=0;              
            }
            else
            {
                flash2+=this_qflash.PE_CH[i];
                if(DetectorZone=="Positive") this_qflash.PE_CH[i]=0;   
            }
            if(this_qflash.PE_CH[i]>0.0)
            {
                NCh_active++;
            }
        }
        auto flashT=flash1+flash2;

        if(norm) this_qflash.norm_this_flash();

        bool pass =  (flashT>limitMinFlash);
        if(DetectorZone=="Positive") pass = pass && (flash1>limitMinFlashSide);
        if(DetectorZone=="Negative") pass = pass && (flash2>limitMinFlashSide);

        if(pass)//flash1>=flash2) // apenas APA 1 E 2
        {
            this_qflash.y = flash->YCenter();
            this_qflash.z = flash->ZCenter();
            this_qflash.y_err = flash->YWidth();
            this_qflash.z_err = flash->ZWidth();
            this_qflash.NCh_active = NCh_active;

            if(isMC)
            {
                this_qflash.time = flash->AbsTime()-ftriggerTime; //flash->Time(); 
            }
            else
            {
                Long64_t timelong = flash->AbsTime();
                this_qflash.time = ((double)(timelong-triggerTime))*16.0/1000.0; //isso estava em ticks de 16ns --> converter para us
            }
             
            this_qflash.time_err = flash->TimeWidth(); 
            this_qflash.Light = this_qflash.TotalLight();
            QFlashs.push_back(this_qflash);
        }
     
    }
    
    return QFlashs;
}


bool MyFlashMatchingGetClusters::getTrackDirFromCaloPoint(
    art::Ptr<recob::Track> const& trk,
    double x, double y, double z,
    float& dirx, float& diry, float& dirz)
{
    const size_t npts = trk->NumberTrajectoryPoints();
    if (npts == 0) return false;

    // ------------------------------------------------------------
    // 1) acha o ponto da track mais próximo geometricamente
    // ------------------------------------------------------------
    double best_dist2 = std::numeric_limits<double>::max();
    int best_idx = -1;

    for (size_t ip = 0; ip < npts; ++ip)
    {
        if (!trk->HasValidPoint(ip)) continue;

        auto const& p = trk->LocationAtPoint(ip);

        double dx = p.X() - x;
        double dy = p.Y() - y;
        double dz = p.Z() - z;

        double dist2 = dx*dx + dy*dy + dz*dz;

        if (dist2 < best_dist2)
        {
            best_dist2 = dist2;
            best_idx = static_cast<int>(ip);
        }
    }

    if (best_idx < 0) return false;

    // ------------------------------------------------------------
    // 2) primeira tentativa: usa a direção da própria track
    // no ponto geometricamente mais próximo
    // ------------------------------------------------------------
    {
        auto const& d = trk->DirectionAtPoint(static_cast<size_t>(best_idx));

        double norm = std::sqrt(
            d.X()*d.X() +
            d.Y()*d.Y() +
            d.Z()*d.Z()
        );

        if (norm > 1.0e-12)
        {
            dirx = d.X() / norm;
            diry = d.Y() / norm;
            dirz = d.Z() / norm;
            return true;
        }
    }

    // ------------------------------------------------------------
    // 3) fallback: PCA local com os pontos geometricamente
    // mais próximos do ponto de calo.
    //
    // Isso NÃO assume que vizinhos no vetor são vizinhos geométricos.
    // ------------------------------------------------------------
    struct NearPoint
    {
        double dist2;
        double x;
        double y;
        double z;
    };

    std::vector<NearPoint> pts;
    pts.reserve(npts);

    for (size_t ip = 0; ip < npts; ++ip)
    {
        if (!trk->HasValidPoint(ip)) continue;

        auto const& p = trk->LocationAtPoint(ip);

        double dx = p.X() - x;
        double dy = p.Y() - y;
        double dz = p.Z() - z;

        double dist2 = dx*dx + dy*dy + dz*dz;

        pts.push_back({dist2, p.X(), p.Y(), p.Z()});
    }

    if (pts.size() < 2) return false;

    std::sort(
        pts.begin(),
        pts.end(),
        [](NearPoint const& a, NearPoint const& b)
        {
            return a.dist2 < b.dist2;
        }
    );

    // Usa no máximo os N pontos geometricamente mais próximos
    const size_t NLOCAL = std::min<size_t>(pts.size(), 8);

    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;

    for (size_t i = 0; i < NLOCAL; ++i)
    {
        mx += pts[i].x;
        my += pts[i].y;
        mz += pts[i].z;
    }

    mx /= static_cast<double>(NLOCAL);
    my /= static_cast<double>(NLOCAL);
    mz /= static_cast<double>(NLOCAL);

    double cxx = 0.0;
    double cxy = 0.0;
    double cxz = 0.0;
    double cyy = 0.0;
    double cyz = 0.0;
    double czz = 0.0;

    for (size_t i = 0; i < NLOCAL; ++i)
    {
        double dx = pts[i].x - mx;
        double dy = pts[i].y - my;
        double dz = pts[i].z - mz;

        cxx += dx*dx;
        cxy += dx*dy;
        cxz += dx*dz;
        cyy += dy*dy;
        cyz += dy*dz;
        czz += dz*dz;
    }

    // ------------------------------------------------------------
    // Power iteration para pegar o maior autovetor da matriz
    // de covariância. Esse vetor é a direção principal local.
    // ------------------------------------------------------------
    double vx = 1.0;
    double vy = 1.0;
    double vz = 1.0;

    for (int it = 0; it < 20; ++it)
    {
        double nx = cxx*vx + cxy*vy + cxz*vz;
        double ny = cxy*vx + cyy*vy + cyz*vz;
        double nz = cxz*vx + cyz*vy + czz*vz;

        double norm = std::sqrt(nx*nx + ny*ny + nz*nz);

        if (norm <= 1.0e-12) return false;

        vx = nx / norm;
        vy = ny / norm;
        vz = nz / norm;
    }

    // ------------------------------------------------------------
    // 4) orienta o sinal usando StartDirection(), se possível.
    // PCA dá eixo, não sentido; então o sinal precisa ser escolhido.
    // ------------------------------------------------------------
    {
        auto const& sd = trk->StartDirection();

        double sd_norm = std::sqrt(
            sd.X()*sd.X() +
            sd.Y()*sd.Y() +
            sd.Z()*sd.Z()
        );

        if (sd_norm > 1.0e-12)
        {
            double dot = vx*sd.X() + vy*sd.Y() + vz*sd.Z();

            if (dot < 0.0)
            {
                vx *= -1.0;
                vy *= -1.0;
                vz *= -1.0;
            }
        }
    }

    dirx = vx;
    diry = vy;
    dirz = vz;

    return true;
}


DEFINE_ART_MODULE(MyFlashMatchingGetClusters)
