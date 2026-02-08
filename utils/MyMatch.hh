#ifndef MYMATCH_HH
#define MYMATCH_HH

#include "my_utils.hh"

#include <TMinuit.h>
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include <map>

class myMatch
{
    public:
        myMatch();
        myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, 
            double drift_speed, double elec_atenuation, double density, double Efield, phot::PhotonVisibilityService const* PVS,  phot::SemiAnalyticalModel const* SAM,
             bool norm=false, std::string DetectorZone = "All", std::string Order = "flash", bool fit_mode=true);
        ~myMatch();

        int Nc;
        int Nf;
        int nc,nf;

        int Nline,Ncol;

        bool fit_mode=true;

        bool checkPossibility(const QCluster* qs,const  QFlash* qf);
        bool startFlash(const QCluster* qs,const  QFlash* qf);
        
        double drift_length;
        double drift_speed;
        double elec_atenuation;
        double density;
        double Efield;

        int CH_MAX=80;
        int APA=0;

        TMinuit* MyMinuit = nullptr;
        int num_var = 1;
        void ChargeHypothesis(const double xoffset);
        void ChargeHypothesis_2(const double xoffset);
        double NLL(); //Negative Log-Likelihood
        static myMatch* s_me;
        static void FCN(Int_t&, Double_t*, Double_t& f, Double_t* x, Int_t); //FUNCAO DE MINIMAZAO PARA MINUIT

        std::vector<std::vector<double>> MYScore;
        std::vector<std::vector<double>> MYOffset;
        std::vector<std::vector<double>> MYdeltaT0;

        std::string type_order = "flash";

        QFlash flash_actual;
        QCluster cluster_actual;

        QFlash flash_fit;
        QCluster cluster_fit;
        phot::PhotonVisibilityService const* fPVS;
        phot::SemiAnalyticalModel const* fSAM;

        struct HungarianResult HR;

        double eff = 0.03;
        bool normPE = false;    
        std::string fDetectorZone = "All";
        std::string fActualDetectorZone = "All";

        bool pause_now=false;

        std::map<int,int> ch_min_map = {{1,120}, {2,40}, {5,80}, {6,0}};
        std::map<int,int> ch_max_map = {{1,159}, {2,79}, {5,119}, {6,39}};

};


#endif