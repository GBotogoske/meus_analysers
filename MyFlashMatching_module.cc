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
        art::InputTag fShowerLabel;
        art::InputTag fCaloLabel;
        art::InputTag fCaloShowerLabel;
        art::InputTag fFlashLabel;
        art::InputTag fHitLabel;

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
        int nShowersPFP;
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
        double density;
        double Efield;

        int nOPdet;

        void returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, int const& thisId, int const& typeObj );
        //void returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<recob::Hit> const& shw_to_hit, art::FindManyP<recob::SpacePoint> const& hit_to_sp, int const& thisId, int const& typeObj );
        void returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<anab::Calorimetry> const& shw_to_calo, int const& thisId,int const& typeObj );

        std::vector<QCluster> getQClustersSlices(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersPFPs(art::Event const& e); // esse aqui eh por slice
        std::vector<QCluster> getQClustersTracks(art::Event const& e); // esse aqui eh por track
        std::vector<QFlash> getFlashs(art::Event const& e);

        phot::PhotonVisibilityService const* fPVS;
        phot::SemiAnalyticalModel const* fSAM;

        bool norm=false;
        bool getShowers = false;
        std::string ClusterType; //Slice,PFP,Track
        std::string DetectorZone; //Positive,Negative,All
        std::string type_order;

        int fflashID,fclusterID,fclusterTrivialID, fclusterType,fclusterTypeTrivial;
        double fflashTime;
        double fScore,fxoffset;
        double fScoreTrivial,fxoffsetTrivial;
        int fflashIDTrivial;
        double fflashTimeTrivial;
        double fdeltat0Trivial, fdeltat0;

        double limitMinFlash = 0 ;
        double limitMinFlashSide = 0;

        double trackLength = 0.0;
};


MyFlashMatching::MyFlashMatching(fhicl::ParameterSet const& p)
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
    ClusterType=p.get<std::string>("ClusterType","Track");
    DetectorZone=p.get<std::string>("DetectorZone","All");
    type_order=p.get<std::string>("type_order","flash");

    limitMinFlash = p.get<double>("limitMinFlash",0.0);
    limitMinFlashSide = p.get<double>("limitMinFlashSide",0.0);

    trackLength = p.get<double>("trackLengthMin",0.0);

    if(norm)
        std::cout << "Normalizing the flash!!! " << std::endl;
    else
        std::cout << "Not Normalizing the flash!!! " << std::endl;

    std::cout << "Track Lenght min: " << trackLength << std::endl;
    std::cout << "limitMinFlash: " << limitMinFlash << std::endl;
    std::cout << "limitMinFlashSide: " << limitMinFlashSide << std::endl;

    getShowers = p.get<bool>("getShowers",false);
    if(getShowers)
        std::cout << "Getting Showers!!! " << std::endl;
    else
        std::cout << "Not Getting Showers!!! " << std::endl;

    std::cout << "type_order: " << type_order<< std::endl;


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
    
}

void MyFlashMatching::beginJob()
{
     art::ServiceHandle<art::TFileService> tfs;

    fTreeF = tfs->make<TTree>("treeMatch", "Flash Match with likelihood fit using poisson and Hungarian algorithm to assing  light to charge");
    fTreeF->Branch("run", &run);
    fTreeF->Branch("event", &event);
    if(type_order=="flash")
    {
        fTreeF->Branch("clusterID", &fclusterID);
        fTreeF->Branch("clusterType", &fclusterType);
        fTreeF->Branch("clusterType_trivial", &fclusterTypeTrivial);
        fTreeF->Branch("clusterID_trivial", &fclusterTrivialID);
        fTreeF->Branch("flashID", &fflashID);
        fTreeF->Branch("flashTime", &fflashTime); 
    }
    if(type_order=="track")
    {
        fTreeF->Branch("clusterID", &fclusterID);
        fTreeF->Branch("clusterType", &fclusterType);
        fTreeF->Branch("flashID", &fflashID);
        fTreeF->Branch("flashTime", &fflashTime);
        fTreeF->Branch("flashID_trivial", &fflashIDTrivial);
        fTreeF->Branch("flashTime_trivial", &fflashTimeTrivial);
    }
    fTreeF->Branch("score", &fScore);
    fTreeF->Branch("x", &fxoffset);
    fTreeF->Branch("score_trivial", &fScoreTrivial);
    fTreeF->Branch("x_trivial", &fxoffsetTrivial);   
    fTreeF->Branch("deltat0_trivial", &fdeltat0Trivial);
    fTreeF->Branch("deltat0", &fdeltat0);   

}

void MyFlashMatching::analyze(art::Event const& e)
{

    /* for (geo::TPCGeo const& tpc : geo->Iterate<geo::TPCGeo>()) {

            geo::Point_t  const cath = tpc.GetCathodeCenter(); // cm
            geo::Vector_t const dir  = tpc.DriftDir();         // unit vector -> anode side
            double        const L    = tpc.DriftDistance();    // cm

            geo::Point_t const anode{
            cath.X() + dir.X()*L,
            cath.Y() + dir.Y()*L,
            cath.Z() + dir.Z()*L
            };

            auto const id = tpc.ID(); // geo::TPCID

            std::cout
            << "Cryo " << id.Cryostat << "  TPC " << id.TPC << "\n"
            << "  Cathode center (cm): (" << cath.X()  << ", " << cath.Y()  << ", " << cath.Z()  << ")\n"
            << "  Anode   center (cm): (" << anode.X() << ", " << anode.Y() << ", " << anode.Z() << ")\n"
            << "  Drift dir: (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")  L=" << L << " cm\n";
        } */

    run    = e.id().run();
    event  = e.id().event();

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
    
    std::cout <<"Final number Qcluster: " << QClusters.size() << std::endl;
    std::cout <<"Final number Flashs: " << QFlashs.size() << std::endl;

    myMatch* match_operator = new myMatch(QClusters,QFlashs,drift_length,drift_speed,electronlife,density,Efield,fPVS,fSAM,norm,DetectorZone,type_order);
    std::cout << "pudim" << std::endl;

    auto& HR = match_operator->HR;
    auto& MYScore = match_operator->MYScore;
    auto& MYOffset = match_operator->MYOffset;
    auto& MYdeltaT0 = match_operator->MYdeltaT0;

    auto& qqs = QClusters;
    auto& qfs = QFlashs;

    if(type_order=="flash")
    {
        const int nRows = std::min((int)HR.assign.size(), match_operator->Nf);

        const bool haveClusters = (!QClusters.empty() && match_operator->Nc > 0);
        for (int nf = 0; nf < nRows; ++nf) 
        {
            
            fflashID = qfs[nf].flashID;
            fflashTime = qfs[nf].time;

            if (haveClusters && (int)MYScore[nf].size() >= match_operator->Nc) 
            {
                auto const& row = MYScore[nf];
                auto it = std::min_element(row.begin(), row.begin() + match_operator->Nc); // só clusters reais
                int ncTrivial = (int)std::distance(row.begin(), it);
                fclusterTrivialID = qqs[ncTrivial].objID; 
                fScoreTrivial = MYScore[nf][ncTrivial];
                fxoffsetTrivial = MYOffset[nf][ncTrivial];
                fclusterTypeTrivial = qqs[ncTrivial].type; 
                fdeltat0Trivial = MYdeltaT0[nf][ncTrivial];
            }
            else 
            {
                fclusterTrivialID = -1;
                fScoreTrivial     = 1e6;
                fxoffsetTrivial   = -1e6;
                fclusterTypeTrivial = -1;
                fdeltat0Trivial = -1e6;
            }
            // se caiu em dummy

            const int nc = HR.assign[nf];
            if (!haveClusters || nc < 0 || nc >= match_operator->Nc)
            {
                fclusterID = -1;
                fScore     = 1e6;
                fxoffset   = -1e6;
                fclusterType = -1;
                fdeltat0 = -1e6;

            } 
            else
            {
                fclusterID = qqs[nc].objID;
                fScore     = MYScore[nf][nc];
                fxoffset   = MYOffset[nf][nc];
                fclusterType = qqs[nc].type;
                fdeltat0 = MYdeltaT0[nf][nc];
            }
            fTreeF->Fill();
        }
    }
    else if(type_order=="track")
    {
        const int Nc = match_operator->Nc;
        const int Nf = match_operator->Nf;

        const int nRows = std::min((int)HR.assign.size(), Nc);
        const bool haveFlashes = (Nf > 0 && !qfs.empty());

        for (int nc = 0; nc < nRows; ++nc)
        {
            // info do cluster (sempre existe aqui)
            fclusterID   = qqs[nc].objID;
            fclusterType = qqs[nc].type;

            // ------------------------
            // TRIVIAL: melhor FLASH p/ esse cluster (mínimo na linha nc)
            // ------------------------
            
            if (haveFlashes && (int)MYScore[nc].size() >= Nf)
            {
                auto const& row = MYScore[nc]; // size Nf
                auto it = std::min_element(row.begin(), row.begin() + Nf);
                int nfTrivial = (int)std::distance(row.begin(), it);

                fflashIDTrivial   = qfs[nfTrivial].flashID;
                fflashTimeTrivial = qfs[nfTrivial].time;

                fScoreTrivial   = MYScore[nc][nfTrivial];
                fxoffsetTrivial = MYOffset[nc][nfTrivial];
                fdeltat0Trivial = MYdeltaT0[nc][nfTrivial];
            }
            else
            {
                fflashIDTrivial   = -1;
                fflashTimeTrivial = -9e6;
                fScoreTrivial     = 1e6;
                fxoffsetTrivial   = -5000.0;
                fdeltat0Trivial = -1e6;
            }

            // ------------------------
            // HUNGARIAN: flash atribuído ao cluster
            // ------------------------
            const int nf = HR.assign[nc]; // agora assign[row=cluster] = col=flash

            if (!haveFlashes || nf < 0 || nf >= Nf)
            {
                fflashID   = -1;
                fflashTime = -9e6;
                fScore     = 1e6;
                fxoffset   = -5000.0;
                fdeltat0 = -1e6;
            }
            else
            {
                fflashID   = qfs[nf].flashID;
                fflashTime = qfs[nf].time;

                fScore   = MYScore[nc][nf];
                fxoffset = MYOffset[nc][nf];
                fdeltat0 = MYdeltaT0[nc][nf];
            }

            fTreeF->Fill();
        }
    }

   
    delete match_operator;
}

void MyFlashMatching::returnQCluster(QCluster& this_qlight, art::Ptr<recob::Track> const& trk, art::FindManyP<anab::Calorimetry> const& trk_to_calo, int const& thisId,int const& typeObj )
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
    
    if(ob_APA != 1 && ob_APA != 2 && ob_APA != 5 && ob_APA != 6) return;
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
        //double drift_time = (drift_length - abs(x))/(drift_speed); //
        //double atten_corr = std::exp(drift_time/electronlife); //

        float pitch;
        float dQ;//, dE;
        //float nphotons;

        //NESTA PARTE O CODIGO DO SBND SEPARA EM 2 PARTES (VALORES NORMAIS E ESTRANHOS)
        if(true)//valores normais (depois tenho que fazer o outro caso)
        {
            pitch = (s < pitch_v.size()) ? pitch_v[s] : -1;
            dQ = dQdx_v[s] ; // * pitch * atten_corr; // corigido pelo drift
            /* dE = dEdx_v[s] * pitch; // talvez precise corrigir pelo drfit, de uma olhada na fcl de reconstrucao depois ...
            nphotons = dE/(W_LAr*1e-6) - dQ;
            nphotons = std::max(0.0f, nphotons); */    
        }
        else
        {
            //aqui depois colocamos os valores estranhos
        }
        this_qlight.push_back(QPoint(x,y,z,dQ,pitch));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
}

void MyFlashMatching::returnQClusterShower(QCluster& this_qlight, art::Ptr<recob::Shower> const& shw, art::FindManyP<anab::Calorimetry> const& shw_to_calo, int const& thisId,int const& typeObj )
{
    std::vector<art::Ptr<anab::Calorimetry>> calos = shw_to_calo.at(shw.key());
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
    if(ob_APA != 1 && ob_APA != 2 && ob_APA != 5 && ob_APA != 6) return;
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
        //double drift_time = (drift_length - abs(x))/(drift_speed); //
        //double atten_corr = std::exp(drift_time/electronlife); //

        float pitch;
        float dQ;//, dE;
        //float nphotons;

        //NESTA PARTE O CODIGO DO SBND SEPARA EM 2 PARTES (VALORES NORMAIS E ESTRANHOS)
        if(true)//valores normais (depois tenho que fazer o outro caso)
        {
            pitch = (s < pitch_v.size()) ? pitch_v[s] : -1;
            dQ = dQdx_v[s] ;// * pitch * atten_corr; // corigido pelo drift
            /* dE = dEdx_v[s] * pitch; // talvez precise corrigir pelo drfit, de uma olhada na fcl de reconstrucao depois ...
            nphotons = dE/(W_LAr*1e-6) - dQ;
            nphotons = std::max(0.0f, nphotons); */
            
        }
        else
        {
            //aqui depois colocamos os valores estranhos
        }
        this_qlight.push_back(QPoint(x,y,z,dQ,pitch));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
        
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
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

    art::Handle<std::vector<recob::Shower>> shower_h;
    //art::Handle<std::vector<recob::Hit>> hit_h; 
    std::optional<art::FindManyP<recob::Shower>> pfp_to_showers;
    std::optional<art::FindManyP<anab::Calorimetry>> shw_to_calo;
    //std::optional<art::FindManyP<recob::Hit>> shower_to_hits ;
    //std::optional<art::FindManyP<recob::SpacePoint>> hit_to_sps ;

    if(getShowers)
    {
        shower_h = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
        //hit_h = e.getHandle<std::vector<recob::Hit>>(fHitLabel);
        pfp_to_showers.emplace(pfp_h, e, fShowerLabel);
        //shower_to_hits.emplace(shower_h, e, fShowerLabel);
        //hit_to_sps.emplace(hit_h, e, fHitLabel);
        shw_to_calo.emplace(shower_h,e,fCaloShowerLabel);
    } 

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
        std::vector<art::Ptr<recob::PFParticle>> pfps = slice_to_pfps.at(sl.key());//i);

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
                if(trk->Length() >= trackLength)
                {
                    returnQCluster(this_qlight,trk,trk_to_calo,sl->ID(),3); 
                }
            }
            if(getShowers)
            {
                std::vector<art::Ptr<recob::Shower>> showers = pfp_to_showers->at(pfp.key());
                nShowersPFP = showers.size();
                for (auto const& shw : showers)
                {
                    //trackID = trk->ID(); //ID do track
                    //pega as info de calorimetria
                    if(shw->Length() >= trackLength)
                    {
                        //returnQClusterShower(this_qlight,trk,trk_to_calo,sl->ID(),3); 
                        //returnQClusterShower(this_qlight,shw,*shower_to_hits,*hit_to_sps,sl->ID(),3); 
                        returnQClusterShower(this_qlight,shw,*shw_to_calo,sl->ID(),3);    
                    }
                }
            }
            

        }
        if(this_qlight.size()>0)
        {
            this_qlight.APA=-1;
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

std::vector<QCluster> MyFlashMatching::getQClustersTracks(art::Event const& e)
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
            }

            //checa se tem pontos no volume ativo
            if(xmint==100000) continue;
            if(xmaxt==-100000) continue;
               
            QCluster this_qlight;
            if(shw->Length() >= trackLength)
            {
                //returnQClusterShower(this_qlight,shw,shower_to_hits,hit_to_sps,shw->ID(),1); 
                returnQClusterShower(this_qlight,shw,shw_to_calo,shw->ID(),1); 
            }
            
            if(this_qlight.size()>0)
            {
                QLigths.push_back(this_qlight);
            } 
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
            this_qflash.PE_CH[this_ch]+=hits[j]->PE();
        }

        auto flash1 = 0.0;
        auto flash2 = 0.0;
      
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
        }
        auto flashT=flash1+flash2;
        /* std::cout <<"(x,y,z) : (" << flash->XCenter() << " , " << flash->YCenter() <<" , " << flash->ZCenter() << ")"<<std::endl;
        std::cout <<"flash1 : "<< flash1 << "   flash 2 : " << flash2 <<std::endl;
        std::cout << flash1/(flash1+flash2) << "   " << flash2/(flash1+flash2) << std::endl;
        std::cout<<"###########"<<std::endl; */
        /* std::cout << "Press ENTER." << std::endl;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); */


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
            this_qflash.time = flash->Time(); 
            this_qflash.time_err = flash->TimeWidth(); 
            QFlashs.push_back(this_qflash);
        }
     
    }
    
    return QFlashs;
}


DEFINE_ART_MODULE(MyFlashMatching)
