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
        art::InputTag fTriggerLabel;

        //este aqui eh para converter ADC em e-
        std::vector<float> _cal_area_const;

        //declaracao da TTree
        TTree* fTreeF;
        TTree* fTreeT;
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

        std::vector<double> fPDEVector;
        double XTalk;

        bool norm=false;
        bool getShowers = false;
        bool fitMode = true;
        bool isMC = true;
        std::string ClusterType; //Slice,PFP,Track
        std::string DetectorZone; //Positive,Negative,All
        std::string type_fit;

        int fflashID,fclusterID,fclusterTrivialID, fclusterType,fclusterTypeTrivial;
        double fflashTime;
        double fScore,fxoffset;
        double fScoreTrivial,fxoffsetTrivial;
        int fflashIDTrivial;
        double fflashTimeTrivial;
        double fdeltat0Trivial, fdeltat0;
        double fcloseAnode, fcloseAnodeTrivial;

        double limitMinFlash = 0 ;
        double limitMinFlashSide = 0;

        double trackLength = 0.0;
        double flength=0.0,flengthTrivial=0.0;
        double fTotalCharge=0.0 , fTotalChargeTrivial=0.0;
        double fTotalLight=0.0 , fTotalLightTrivial=0.0;

        double fvisEff, fvisEffTrivial;
        int fNCh_Flash,fNCh_FlashTrivial;

        std::vector<double> flash_pe,flash_peTrivial;
        std::vector<double> charge_pe,charge_peTrivial;

        double ftriggerTime = 0.0;
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
    if(norm)
        std::cout << "Normalizing the flash!!! " << std::endl;
    else
        std::cout << "Not Normalizing the flash!!! " << std::endl;

    fitMode = p.get<bool>("FitMode",true);
    if(fitMode)
        std::cout << "Fitting mode!!! " << std::endl;
    else
        std::cout << "NOTTTT Fitting mode!!! " << std::endl;

    isMC = p.get<bool>("isMC",true);
    if(isMC)
        std::cout << "Analysing MC!!! " << std::endl;
    else
        std::cout << "Analysing Real Data!!! " << std::endl;

    ClusterType=p.get<std::string>("ClusterType","Track");
    DetectorZone=p.get<std::string>("DetectorZone","All");
    type_fit=p.get<std::string>("type_fit","flash");

    limitMinFlash = p.get<double>("limitMinFlash",0.0);
    limitMinFlashSide = p.get<double>("limitMinFlashSide",0.0);

    trackLength = p.get<double>("trackLengthMin",0.0);

    std::cout << "Track Length min: " << trackLength << std::endl;
    std::cout << "limitMinFlash: " << limitMinFlash << std::endl;
    std::cout << "limitMinFlashSide: " << limitMinFlashSide << std::endl;

    getShowers = p.get<bool>("getShowers",false);
    if(getShowers)
        std::cout << "Getting Showers!!! " << std::endl;
    else
        std::cout << "Not Getting Showers!!! " << std::endl;

    std::cout << "type_fit: " << type_fit<< std::endl;

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
    
    fTriggerLabel = p.get<std::string>("TriggerTag","daq:trigger:pdhdkeepupstage1");

    fPDEVector = p.get<std::vector<double>>("PDEvector",std::vector<double>(nOPdet,0.03));
    XTalk =p.get<double>("XTalk",0.01);

    std::cout << "PDE loaded " << fPDEVector[0] << std::endl;
    std::cout << "Xtalk loaded " << XTalk << std::endl;

    flash_pe=std::vector<double>(nOPdet);
    flash_peTrivial=std::vector<double>(nOPdet);
    charge_pe=std::vector<double>(nOPdet);
    charge_peTrivial=std::vector<double>(nOPdet);

    ftriggerTime = p.get<double>("triggerTime",0.0);
    std::cout << ftriggerTime << " : Trigger Time loaded (only used when MONTE CARLO)" << std::endl;
}

void MyFlashMatching::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTreeF = tfs->make<TTree>("treeMatchF", "Flash Match with likelihood fit using poisson and Hungarian algorithm to assing  light to charge");
    
    fTreeF->Branch("run", &run);
    fTreeF->Branch("event", &event);
    fTreeF->Branch("clusterID", &fclusterID);
    fTreeF->Branch("clusterType", &fclusterType);
    fTreeF->Branch("clusterType_trivial", &fclusterTypeTrivial);
    fTreeF->Branch("clusterID_trivial", &fclusterTrivialID);
    fTreeF->Branch("flashID", &fflashID);
    fTreeF->Branch("flashTime", &fflashTime);
    fTreeF->Branch("NCh_flash", &fNCh_Flash);
    fTreeF->Branch("score", &fScore);
    fTreeF->Branch("x", &fxoffset);
    fTreeF->Branch("score_trivial", &fScoreTrivial);
    fTreeF->Branch("x_trivial", &fxoffsetTrivial);   
    fTreeF->Branch("deltat0_trivial", &fdeltat0Trivial);
    fTreeF->Branch("deltat0", &fdeltat0); 
    fTreeF->Branch("length_trivial", &flengthTrivial);
    fTreeF->Branch("length", &flength); 
    fTreeF->Branch("light", &fTotalLight); 
    fTreeF->Branch("charge_trivial", &fTotalChargeTrivial);
    fTreeF->Branch("charge", &fTotalCharge); 
    fTreeF->Branch("distanceAnode", &fcloseAnode);
    fTreeF->Branch("distanceAnode_trivial", &fcloseAnodeTrivial); 
    fTreeF->Branch("visEff", &fvisEff);
    fTreeF->Branch("visEff_trivial", &fvisEffTrivial); 
    fTreeF->Branch("flashPE",&flash_pe);
    fTreeF->Branch("chargePE",&charge_pe);
    fTreeF->Branch("chargePE_trivial",&charge_peTrivial);

    fTreeT = tfs->make<TTree>("treeMatchT", "Flash Match 2 with likelihood fit using poisson and Hungarian algorithm to assing  light to charge");

    fTreeT->Branch("run", &run);
    fTreeT->Branch("event", &event);
    fTreeT->Branch("clusterID", &fclusterID);
    fTreeT->Branch("clusterType", &fclusterType);
    fTreeT->Branch("flashID", &fflashID);
    fTreeT->Branch("flashTime", &fflashTime);
    fTreeT->Branch("flashID_trivial", &fflashIDTrivial);
    fTreeT->Branch("flashTime_trivial", &fflashTimeTrivial);
    fTreeT->Branch("NCh_flash", &fNCh_Flash);
    fTreeT->Branch("NCh_flash_trivial", &fNCh_FlashTrivial);
    fTreeT->Branch("score", &fScore);
    fTreeT->Branch("x", &fxoffset);
    fTreeT->Branch("score_trivial", &fScoreTrivial);
    fTreeT->Branch("x_trivial", &fxoffsetTrivial);   
    fTreeT->Branch("deltat0_trivial", &fdeltat0Trivial);
    fTreeT->Branch("deltat0", &fdeltat0);
    fTreeT->Branch("length", &flength); 
    fTreeT->Branch("light", &fTotalLight); 
    fTreeT->Branch("light_trivial", &fTotalLightTrivial);
    fTreeT->Branch("charge", &fTotalCharge); 
    fTreeT->Branch("distanceAnode", &fcloseAnode);
    fTreeT->Branch("distanceAnode_trivial", &fcloseAnodeTrivial);
    fTreeT->Branch("visEff", &fvisEff);
    fTreeT->Branch("visEff_trivial", &fvisEffTrivial);
    fTreeT->Branch("flashPE",&flash_pe);
    fTreeT->Branch("flashPE_trivial",&flash_peTrivial);
    fTreeT->Branch("chargePE",&charge_pe);
    fTreeT->Branch("chargePE_trivial",&charge_peTrivial);
}

void MyFlashMatching::analyze(art::Event const& e)
{
    /* for (geo::TPCGeo const& tpc : geo->Iterate<geo::TPCGeo>()) 
    {
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

    myMatch* match_operator = new myMatch(QClusters,QFlashs,drift_length,drift_speed,electronlife,density,Efield,fPVS,fSAM,
        fPDEVector,XTalk,norm,DetectorZone,type_fit,fitMode);
    std::cout << "pudim" << std::endl;

    auto& HR = match_operator->HR;
    auto& MYScore = match_operator->MYScore;
    auto& MYOffset = match_operator->MYOffset;
    auto& MYdeltaT0 = match_operator->MYdeltaT0;
    auto& MYcloseAnode = match_operator->MYcloseAnode;
    auto& MYvisEf = match_operator->MYvisEf;

    auto& qqs = QClusters;
    auto& qfs = QFlashs;

    const int Nc = match_operator->Nc;
    const int Nf = match_operator->Nf;

    const int nRows = std::min((int)HR.row2col.size(), Nf);
    const bool haveClusters = (!QClusters.empty() && Nc > 0);
    for (int nf = 0; nf < nRows; ++nf) 
    {
        fflashID = qfs[nf].flashID;
        fflashTime = qfs[nf].time;
        fTotalLight = qfs[nf].TotalLight();
        fNCh_Flash = qfs[nf].NCh_active;

        flash_pe = qfs[nf].PE_CH;

        fclusterTrivialID = -1;
        fScoreTrivial     = 1e12;
        fxoffsetTrivial   = -1e6;
        fclusterTypeTrivial = -1;
        fdeltat0Trivial = -1e6;
        flength=-100.0;
        flengthTrivial=-100.0;
        fTotalCharge=-100.0;
        fTotalChargeTrivial=-100.0;
        fcloseAnodeTrivial = 1e8;
        fvisEffTrivial = -1e6;
        charge_peTrivial = std::vector<double>(nOPdet,-1.0);
        
        if (haveClusters && (int)MYScore[nf].size() >= Nc) 
        {
            auto const& row = MYScore[nf];
            auto it = std::min_element(row.begin(), row.begin() + Nc); // só clusters reais
            int ncTrivial = (int)std::distance(row.begin(), it);
            fclusterTrivialID = qqs[ncTrivial].objID; 
            fScoreTrivial = MYScore[nf][ncTrivial];
            fxoffsetTrivial = MYOffset[nf][ncTrivial];
            fclusterTypeTrivial = qqs[ncTrivial].type; 
            fdeltat0Trivial = MYdeltaT0[nf][ncTrivial];
            flengthTrivial = qqs[ncTrivial].Length; 
            fTotalChargeTrivial =  qqs[ncTrivial].TotalCharge();
            fcloseAnodeTrivial = MYcloseAnode[nf][ncTrivial];
            fvisEffTrivial = MYvisEf[nf][ncTrivial];

            match_operator->cluster_actual = qqs[ncTrivial]; //setar cluster
            match_operator->ChargeHypothesis(fdeltat0Trivial);//calcular flash hip
            charge_peTrivial= match_operator->flash_fit.PE_CH;//pegar flash hip

        }
        // se caiu em dummy

        const int nc = HR.row2col[nf];
        if (!haveClusters || nc < 0 || nc >= Nc)
        {
            fclusterID = -1;
            fScore     = 1e12;
            fxoffset   = -1e6;
            fclusterType = -1;
            fdeltat0 = -1e6;
            fcloseAnode = 1e8;
            fvisEff = -1e6;
            charge_pe = std::vector<double>(nOPdet,-1.0);
        } 
        else
        {
            fclusterID = qqs[nc].objID;
            fScore     = MYScore[nf][nc];
            fxoffset   = MYOffset[nf][nc];
            fclusterType = qqs[nc].type;
            fdeltat0 = MYdeltaT0[nf][nc];
            flength = qqs[nc].Length;
            fTotalCharge = qqs[nc].TotalCharge();
            fcloseAnode = MYcloseAnode[nf][nc];
            fvisEff = MYvisEf[nf][nc];

            match_operator->cluster_actual = qqs[nc]; //setar cluster
            match_operator->ChargeHypothesis(fdeltat0);//calcular flash hip
            charge_pe= match_operator->flash_fit.PE_CH;//pegar flash hip

        }
        fTreeF->Fill();
    }
    const int nCols = std::min((int)HR.col2row.size(), Nc);
    const bool haveFlashes = (Nf > 0 && !qfs.empty());
    for (int nc = 0; nc < nCols; ++nc)
    {
        // info do cluster (sempre existe aqui)
        fclusterID   = qqs[nc].objID;
        fclusterType = qqs[nc].type;
        flength = qqs[nc].Length;
        fTotalCharge = qqs[nc].TotalCharge();

        // ------------------------
        // TRIVIAL: melhor FLASH p/ esse cluster (mínimo na linha nc)
        // ------------------------
        fflashIDTrivial   = -1;
        fflashTimeTrivial = -9e6;
        fScoreTrivial     = 1e12;
        fxoffsetTrivial   = -5000.0;
        fdeltat0Trivial   = -1e6;
        fTotalLight=-100.0;
        fTotalLightTrivial=-100.0;
        fcloseAnodeTrivial = 1e8;
        fvisEffTrivial = -1e6;
        fNCh_FlashTrivial = -5000;
        charge_peTrivial = std::vector<double>(nOPdet,-1.0);
        flash_peTrivial = std::vector<double>(nOPdet,-1.0);

        fflashID   = -1;
        fflashTime = -9e6;
        fScore     = 1e12;
        fxoffset   = -5000.0;
        fdeltat0   = -1e6;
        fcloseAnode = 1e8;
        fvisEff = -1e6;
        fNCh_Flash = -5000;
        charge_pe = std::vector<double>(nOPdet,-1.0);
        flash_pe = std::vector<double>(nOPdet,-1.0);

        match_operator->cluster_actual = qqs[nc]; //setar cluster

        if (haveFlashes)
        {
            int nfTrivial = -1;
            double best = 1e100;

            for (int nf = 0; nf < Nf; ++nf)
            {
                if ((int)MYScore[nf].size() <= nc) continue;
                double s = MYScore[nf][nc];
                if (s < best) { best = s; nfTrivial = nf; }
            }
            if (nfTrivial >= 0) 
            {
                fflashIDTrivial   = qfs[nfTrivial].flashID;
                fflashTimeTrivial = qfs[nfTrivial].time;
                fNCh_FlashTrivial = qfs[nfTrivial].NCh_active;
                fScoreTrivial     = MYScore[nfTrivial][nc];
                fxoffsetTrivial   = MYOffset[nfTrivial][nc];
                fdeltat0Trivial   = MYdeltaT0[nfTrivial][nc];
                fTotalLightTrivial = qfs[nfTrivial].TotalLight();
                fcloseAnodeTrivial = MYcloseAnode[nfTrivial][nc];
                fvisEffTrivial = MYvisEf[nfTrivial][nc];

                flash_peTrivial = qfs[nfTrivial].PE_CH;
                match_operator->ChargeHypothesis(fdeltat0Trivial);//calcular flash hip
                charge_peTrivial= match_operator->flash_fit.PE_CH;//pegar flash hip

            }
        }

        // ------------------------
        // HUNGARIAN: flash atribuído ao cluster
        // ------------------------
        const int nf = HR.col2row[nc]; // agora row2col[row=cluster] = col=flash
        if (haveFlashes && nf >= 0 && nf < Nf)
        {
            fflashID   = qfs[nf].flashID;
            fflashTime = qfs[nf].time;
            fNCh_Flash = qfs[nf].NCh_active;
            fScore     = MYScore[nf][nc];
            fxoffset   = MYOffset[nf][nc];
            fdeltat0   = MYdeltaT0[nf][nc];
            fTotalLight = qfs[nf].TotalLight();
            fcloseAnode = MYcloseAnode[nf][nc];
            fvisEff = MYvisEf[nf][nc];

            flash_pe = qfs[nf].PE_CH;
            match_operator->ChargeHypothesis(fdeltat0);//calcular flash hip
            charge_pe= match_operator->flash_fit.PE_CH;//pegar flash hip
        }
        fTreeT->Fill();
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
    //std::cout << ob_APA << std::endl;
    
    if(ob_APA != 1 && ob_APA != 2 && ob_APA != 5 && ob_APA != 6) return;

    if(ob_APA==1) return; // to tirando o apa1, pode colocar dpeois se quiser

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
        this_qlight.push_back(QPoint(x,y,z,dQ,-1,-1,-1,pitch,ob_APA));
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
        this_qlight.push_back(QPoint(x,y,z,dQ,-1,-1,-1,pitch,ob_APA));
        //this_qlight.push_back(QPoint(x,y,z,nphotons));
        
    }
    this_qlight.objID = thisId;
    this_qlight.type = typeObj;
    this_qlight.APA = ob_APA;
}

//eu vi que para cosmics eh melhor fazer match diretamente com track/PFParticle em vez de slices ( muito quebrado para cosmics )

std::vector<QCluster> MyFlashMatching::getQClustersSlices(art::Event const& e)
{
    auto slice_h = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
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
            mf::LogWarning("MyFlashMatching")
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
                returnQCluster(this_qlight, trk, trk_to_calo, sl->ID(), 3);

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
                QLights.push_back(this_qlight);
            }
        }
    }

    return QLights;
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
            this_qlight.Length = trk->Length();
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
            this_qflash.PE_CH[this_ch]+=hits[j]->PE();
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
            QFlashs.push_back(this_qflash);
        }
     
    }
    
    return QFlashs;
}


DEFINE_ART_MODULE(MyFlashMatching)
