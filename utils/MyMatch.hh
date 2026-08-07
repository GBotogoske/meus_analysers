#ifndef MYMATCH_HH
#define MYMATCH_HH

#include "my_utils.hh"

#include <TMinuit.h>
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"


#include <map>

namespace spacecharge {
  class SpaceCharge;
}

class myMatch
{
    public:
        myMatch();
        myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, 
            double drift_speed, double elec_atenuation, double density, double Efield, phot::PhotonVisibilityService const* PVS,  phot::SemiAnalyticalModel const* SAM,
            const std::vector<double>& eff, std::vector<double>& XTalk, std::vector<int>& CHActive, bool fuseSCE, spacecharge::SpaceCharge const* sce_service, 
            bool norm=false, std::string DetectorZone = "All", std::string typeFit = "flash", bool fit_mode=true);
        myMatch(double drift_length, double drift_speed, double elec_atenuation, double density, double Efield, 
            phot::PhotonVisibilityService const* PVS,  phot::SemiAnalyticalModel const* SAM,
            const std::vector<double>& eff, std::vector<double>& XTalk, std::vector<int>& CHActive, 
            bool fuseSCE, spacecharge::SpaceCharge const* sce_service);
        ~myMatch();

        int Nc;
        int Nf;
        int nc,nf;

        int Nline,Ncol;

        bool fit_mode=true;

        bool checkPossibility(const QCluster* qs,const  QFlash* qf, double time_buffer=300);
        bool startFlash(const QCluster* qs,const  QFlash* qf);

        double returnVisEff(const QCluster* qs, double xoffset);
        double returnVisEff(const QCluster* qs, const QFlash* qf, const double xoffset);
        std::vector<double> returnVisEffCh(const QCluster* qs, const QFlash* qf, const double xoffset);
        std::vector<double> returndCh(const QCluster* qs,const double xoffset,const std::vector<double>& xch,const std::vector<double>& ych,const std::vector<double>& zch);
        double returnVisEff();
        double returnVisEffLight();
        std::vector<double> returnVisEffCh();
        std::vector<double> returndCh(const std::vector<double>& xch,const std::vector<double>& ych,const std::vector<double>& zch);

        void fixPositionSce();
        double getLocalEfield(double x, double y, double z, int tpcid) const;

        double drift_length;
        double drift_speed;
        double elec_atenuation;
        double density;
        double Efield;

        double save_charge=0.0;

        int CH_MAX=80;
        int APA=0;

        TMinuit* MyMinuit = nullptr;
        int num_var = 1;
        void ChargeHypothesis(const double xoffset);
        double NLL(); //Negative Log-Likelihood
        static myMatch* s_me;
        static void FCN(Int_t&, Double_t*, Double_t& f, Double_t* x, Int_t); //FUNCAO DE MINIMAZAO PARA MINUIT

        std::vector<std::vector<double>> MYScore;
        std::vector<std::vector<double>> MYOffset;
        std::vector<std::vector<double>> MYdeltaT0;
        std::vector<std::vector<double>> MYcloseAnode;
        std::vector<std::vector<double>> MYvisEf;

        std::string type_fit = "flash";

        QFlash flash_actual;
        QCluster cluster_actual;

        QFlash flash_fit;
        QCluster cluster_fit;
        phot::PhotonVisibilityService const* fPVS;
        phot::SemiAnalyticalModel const* fSAM;
        const spacecharge::SpaceCharge* fSCE = nullptr;

        struct HungarianResult HR;

        double eff = 0.03;
        bool normPE = false;    
        std::string fDetectorZone = "All";
        std::string fActualDetectorZone = "All";

        bool pause_now=false;

        std::map<int,int> ch_min_map = {{1,120}, {2,40}, {5,80}, {6,0}};
        std::map<int,int> ch_max_map = {{1,159}, {2,79}, {5,119}, {6,39}};

        std::vector<double> direct_visibilities;
        std::vector<double> effVector;
        std::vector<double> XtalkVector;
        std::vector<int> CHActiveVector;

        std::vector<double> ratio_sce;

        double Xtalk;
        double Kdup;
        bool useSCE = false;
        double vis_map_factor = 1.0;

};


#endif