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
#include "lardataobj/RecoBase/TrackHitMeta.h"

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


class MyFlashMatchingGetClustersEfield : public art::EDAnalyzer
 {
    public:
        explicit MyFlashMatchingGetClustersEfield(fhicl::ParameterSet const& p);

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
        double fvis;

        double drift_length;
        double drift_speed,drift_speed_std;
        double electronlife;
        double W_LAr;
        double density;
        double Efield,Efield_std;

        double gain_factor=1.0;

        void returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, art::FindManyP<recob::Hit, recob::TrackHitMeta> const& trk_to_hit_meta,int const& thisId, int const& typeObj );
        //void returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<recob::Hit> const& shw_to_hit, art::FindManyP<recob::SpacePoint> const& hit_to_sp, int const& thisId, int const& typeObj );
        void returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<anab::Calorimetry> const& shw_to_calo, int const& thisId,int const& typeObj );


        bool getCorrectedDirectionFromIndex(art::Ptr<recob::Track> const& trk,int indp,double scalev,float& dirx,float& diry,float& dirz,
            double& stretch,double& ux,double& uy,double& uz) const;

        int findClosestValidTrackPoint(art::Ptr<recob::Track> const& trk,double x,double y,double z) const;

        std::vector<QCluster> getQClustersSlices(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersPFPs(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersTracks(art::Event const& e); // esse aqui eh por track
        std::vector<QFlash> getFlashs(art::Event const& e);

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

        std::vector<int> run_list;
        std::vector<double> E_list;
        std::vector<double> v_list;


};

MyFlashMatchingGetClustersEfield::MyFlashMatchingGetClustersEfield(fhicl::ParameterSet const& p)
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
    
    run_list = p.get<std::vector<int>>("run_list",std::vector<int>(1,1));
    E_list = p.get<std::vector<double>>("E_list",std::vector<double>(1,0.5));
    v_list = p.get<std::vector<double>>("v_list",std::vector<double>(1,0.16));

}

void MyFlashMatchingGetClustersEfield::beginJob()
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
    fTreeFT->Branch("flash", &fFlash_fit);
    fTreeFT->Branch("cluster", &fCluster_fit);
    fTreeFT->Branch("LYCH", &LYCH);
    fTreeFT->Branch("dCH", &dch);
    
}

void MyFlashMatchingGetClustersEfield::analyze(art::Event const& e)
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

    this->Efield_std=this->Efield;
    this->drift_speed_std=this->drift_speed;

    for (size_t i = 0; i < run_list.size(); ++i)
    {
        if(run_list[i]==run)
        {
            drift_speed=v_list[i];
            this->Efield=E_list[i];
        }
    }
    std::cout << "run := " << run << std::endl;
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
        mf::LogWarning("MyFlashMatchingGetClustersEfield")
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
                fFlash_fit.Light = fFlash_fit.TotalLight();
                fCluster_fit.Charge = match_operator->save_charge;
                fCluster_fit.Length = fCluster_fit.calcLength();
                fvis = match_operator->returnVisEff();
                LYCH =  match_operator->returnVisEffCh();
                dch = match_operator->returndCh(xch,ych,zch);
                fTreeFT->Fill();
           }

        }
    }
    delete match_operator;
}

void MyFlashMatchingGetClustersEfield::returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, art::FindManyP<recob::Hit, recob::TrackHitMeta> const& trk_to_hit_meta, int const& thisId,int const& typeObj )
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

    if(ob_APA==1) return; // to tirando o apa1, pode colocar dpeois se quiser

    //------------------------- termino de buscar o plano -------------------------------------------------------------------
    
    auto const& dEdx_v     = calo->dEdx();
    auto const& dADCdx_v   = calo->dQdx();
    auto const& pitch_v    = calo->TrkPitchVec();
    auto const& pos_v      = calo->XYZ();
    auto const& indexpoints = calo->TpIndices();

    const double scalev = drift_speed / drift_speed_std;

    // Usa só o tamanho comum para não acessar vetor fora do range
    size_t npts = dEdx_v.size();
    npts = std::min(npts, dADCdx_v.size());
    npts = std::min(npts, pitch_v.size());
    npts = std::min(npts, pos_v.size());

    if (npts == 0) return;

    std::unordered_map<int, int> hitKeyToTrajIndex;

    if (trk_to_hit_meta.isValid())
    {
        auto trkHits  = trk_to_hit_meta.at(trk.key());
        auto trkMetas = trk_to_hit_meta.data(trk.key());

        size_t nAssoc = std::min(trkHits.size(), trkMetas.size());

        for (size_t ih = 0; ih < nAssoc; ++ih)
        {
            if (trkHits[ih].isNull()) continue;
            if (trkMetas[ih] == nullptr) continue;

            int hitKey = static_cast<int>(trkHits[ih].key());

            auto idx = trkMetas[ih]->Index();

            if (idx < trk->NumberTrajectoryPoints() && trk->HasValidPoint(idx))
            {
                hitKeyToTrajIndex[hitKey] = static_cast<int>(idx);
            }
        }
    }

    // Define o plano de anodo/cátodo desse APA
    double xAPA = 0.0;

    if (ob_APA == 2 || ob_APA == 6) { xAPA = +drift_length;}
    else if (ob_APA == 1 || ob_APA == 5) { xAPA = -drift_length; }
    else {return;}

    // create vector of e- instead of ADC units
    std::vector<float> dQdx_v(npts, 0.0);

    for (size_t s = 0; s < npts; s++)
    {
        dQdx_v[s] = dADCdx_v[s] * (1.0 / _cal_area_const.at(plane)) * gain_factor;
    }

    //std::cout << "calos : " << plane << " - " << calo->PlaneID().Plane << std::endl;
    //varre todas as posicoes/energia depositadas
    for (size_t s = 0; s < npts; s++)
    {
        float x = pos_v[s].X();
        float y = pos_v[s].Y();
        float z = pos_v[s].Z();

        x = xAPA + scalev * (x - xAPA);
        // ----------------------------
        // Direção
        // ----------------------------
        float dirx = -10.0;
        float diry = -10.0;
        float dirz = -10.0;

        double stretch = 1.0;
        bool validDir = false;

        double ux = 0.0;
        double uy = 0.0;
        double uz = 0.0;

        int indp = -1;

        // 1) Primeiro usa TpIndices como HIT KEY
        if (s < indexpoints.size())
        {
            int hitKey = indexpoints[s];

            auto it = hitKeyToTrajIndex.find(hitKey);

            if (it != hitKeyToTrajIndex.end())
            {
                indp = it->second;
            }
        }

        // 2) Se não achou via TrackHitMeta, usa fallback geométrico
        if (indp < 0)
        {
            indp = findClosestValidTrackPoint(
                trk,
                pos_v[s].X(),
                pos_v[s].Y(),
                pos_v[s].Z()
            );
        }

        // 3) Calcula direção corrigida
        validDir = getCorrectedDirectionFromIndex(trk,indp,scalev,dirx,diry,dirz,stretch,ux,uy,uz);
        // ----------------------------
        // Corte por lado depois da correção de x
        // ----------------------------
        if (DetectorZone == "Positive" && x < 0) {
            continue;
        }
        else if (DetectorZone == "Negative" && x > 0) {
            continue;
        }

        // ----------------------------
        // Pitch e carga
        // ----------------------------
        float pitch = -1.0;
        float dQdx_corr = -1.0;

        if (s < pitch_v.size() && s < dQdx_v.size())
        {
            float pitch_old = pitch_v[s];

            if (pitch_old > 0.0 && validDir)
            {
                /* std::cout << "######### TRACK Entrei  #########"<<std::endl; */
                pitch = pitch_old * stretch;

                // Se dQdx_v[s] é carga por unidade de comprimento,
                // então a carga total do segmento deve ser preservada:
                //
                // dQ_old = dQdx_old * pitch_old
                // dQdx_new = dQ_old / pitch_new
                //
                dQdx_corr = dQdx_v[s] * pitch_old / pitch;
            }
            else
            {
                // fallback: sem correção geométrica
                pitch = pitch_old;
                dQdx_corr = dQdx_v[s];
            }
        }
        /* std::cout << "######### TRACK  #########"<<std::endl;
        std::cout << ux << " " << uy << " " << uz << " " << scalev << std::endl; 
        std::cout << dirx << " " << diry << " " << dirz << " " << stretch << std::endl;
        std::cout << "##################"<<std::endl; */
        this_qlight.push_back(QPoint(x, y, z, dQdx_corr, dirx, diry, dirz, pitch, ob_APA));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
}

void MyFlashMatchingGetClustersEfield::returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<anab::Calorimetry> const& shw_to_calo, int const& thisId,int const& typeObj )
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
    if(ob_APA==1) return; // to tirando o apa1, pode colocar dpeois se quiser
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
        const double scalev = drift_speed / drift_speed_std;

        // ----------------------------
        // Corrige coordenada x
        // ----------------------------
        double xAPA = 0.0;
        bool validAPA = true;

        if (ob_APA == 2 || ob_APA == 6) {
            xAPA = +drift_length;
        }
        else if (ob_APA == 1 || ob_APA == 5) {
            xAPA = -drift_length;
        }
        else {
            validAPA = false;
        }

        if (validAPA) {
            x = xAPA + scalev * (x - xAPA);
        }

        // Direção do shower
        const auto& dir = shw->Direction();

        double ux = dir.X();
        double uy = dir.Y();
        double uz = dir.Z();

        float dirx = -10.0;
        float diry = -10.0;
        float dirz = -10.0;

        double stretch = 1.0;
        bool validDir = false;

        // Aplica a deformação em x
        double vx = scalev * ux;
        double vy = uy;
        double vz = uz;

        stretch = std::sqrt(vx*vx + vy*vy + vz*vz);

        if (stretch > 0.0)
        {
            dirx = vx / stretch;
            diry = vy / stretch;
            dirz = vz / stretch;

            validDir = true;
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

        float pitch = -1.0;
        float q_corr = -1.0;

        if (s < pitch_v.size() && s < dQdx_v.size())
        {
            float pitch_old = pitch_v[s];

            if (pitch_old > 0.0 && validDir)
            {
                //std::cout << "######### Shower Entrei  #########"<<std::endl;
                pitch = pitch_old * stretch;

                // Se dQdx_v[s] é dQ/dx, preserve a carga total:
                // dQ_total = dQdx_old * pitch_old
                // dQdx_new = dQ_total / pitch_new
                q_corr = dQdx_v[s] * pitch_old / pitch;
            }
            else
            {
                pitch = pitch_old;
                q_corr = dQdx_v[s];
            }
        }
       /*  std::cout << "######### SHOWER  #########"<<std::endl;
        std::cout << ux << " " << uy << " " << uz << " " << scalev << std::endl; 
        std::cout << dirx << " " << diry << " " << dirz << " " << stretch << std::endl;
        std::cout << "##################"<<std::endl; */
        this_qlight.push_back(QPoint(x,y,z,q_corr,dirx,diry,dirz,pitch,ob_APA));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
}

//eu vi que para cosmics eh melhor fazer match diretamente com track/PFParticle em vez de slices ( muito quebrado para cosmics )

std::vector<QCluster> MyFlashMatchingGetClustersEfield::getQClustersSlices(art::Event const& e)
{
    auto slice_h = e.getHandle<std::vector<recob::Slice>>(fSliceLabel);
    if (!slice_h) {
        mf::LogWarning("MyFlashMatchingGetClustersEfield")
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
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(track_h, e, fTrackLabel);

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
            mf::LogWarning("MyFlashMatchingGetClustersEfield")
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
                returnQCluster(this_qlight, trk, trk_to_calo, fmHits, sl->ID(), 3);

                if (this_qlight.size() > nBefore)
                {
                    totalLength += trk->Length();
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

            if (totalLength >= trackLength)
            {
                this_qlight.Charge = this_qlight.TotalCharge();
                QLights.push_back(this_qlight);
            }
        }
    }

    return QLights;
}


std::vector<QCluster> MyFlashMatchingGetClustersEfield::getQClustersPFPs(art::Event const& e)
{
   
    // pegar produtos importantes
    
    auto pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel); //PFParticles
    auto track_h = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel); // tracks
    auto calo_h = e.getValidHandle<std::vector<anab::Calorimetry>>(fCaloLabel); //calorimetria

    art::FindManyP<recob::Track> pfp_to_tracks(pfp_h, e, fTrackLabel); // pega a associacao de tracks das PFParticles
    art::FindManyP<anab::Calorimetry> trk_to_calo(track_h, e, fCaloLabel); // pega as info de calorimetria dos tracks

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(track_h, e, fTrackLabel);
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
                returnQCluster(this_qlight, trk, trk_to_calo, fmHits, pfp->Self(), 2);
            }
        }
        if(this_qlight.size()>0)
        {
            QLigths.push_back(this_qlight);
        }  
    }
    return QLigths;
}

std::vector<QCluster> MyFlashMatchingGetClustersEfield::getQClustersTracks(art::Event const& e)
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

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(track_h, e, fTrackLabel);
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
            returnQCluster(this_qlight, trk, trk_to_calo, fmHits, trk->ID(), 0);
            this_qlight.Length = trk->Length();
            this_qlight.Charge = this_qlight.TotalCharge();
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
            }
            
            if(this_qlight.size()>0)
            {
                QLigths.push_back(this_qlight);
            } 
        }
    }
    
    return QLigths;
}

std::vector<QFlash> MyFlashMatchingGetClustersEfield::getFlashs(art::Event const& e)
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
      mf::LogWarning("MyFlashMatchingGetClustersEfield") << "Cannot load any flashes. Failing";
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

bool MyFlashMatchingGetClustersEfield::getCorrectedDirectionFromIndex(art::Ptr<recob::Track> const& trk, int indp,double scalev, float& dirx, float& diry,float& dirz,
    double& stretch, double& ux,double& uy,double& uz) const
{
    if (indp < 0) return false;
    if (!trk->HasValidPoint(indp)) return false;

    auto const& dir_v = trk->DirectionAtPoint(indp);

    ux = dir_v.X();
    uy = dir_v.Y();
    uz = dir_v.Z();

    double norm = std::sqrt(ux*ux + uy*uy + uz*uz);

    if (norm <= 0.0 || !std::isfinite(norm)) return false;

    // Normaliza a direção original por segurança
    ux /= norm;
    uy /= norm;
    uz /= norm;

    stretch = std::sqrt(scalev*scalev*ux*ux + uy*uy + uz*uz);

    if (stretch <= 0.0 || !std::isfinite(stretch)) return false;

    dirx = scalev * ux / stretch;
    diry = uy / stretch;
    dirz = uz / stretch;

    return true;
}

int MyFlashMatchingGetClustersEfield::findClosestValidTrackPoint( art::Ptr<recob::Track> const& trk,double x,double y,double z) const
{
    int bestIndex = -1;
    double bestDist2 = 1.0e30;

    for (size_t i = 0; i < trk->NumberTrajectoryPoints(); ++i)
    {
        if (!trk->HasValidPoint(i)) continue;

        auto const& p = trk->LocationAtPoint(i);

        double dx = p.X() - x;
        double dy = p.Y() - y;
        double dz = p.Z() - z;

        double dist2 = dx*dx + dy*dy + dz*dz;

        if (dist2 < bestDist2)
        {
            bestDist2 = dist2;
            bestIndex = static_cast<int>(i);
        }
    }

    return bestIndex;
}


DEFINE_ART_MODULE(MyFlashMatchingGetClustersEfield)
