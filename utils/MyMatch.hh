#ifndef MYMATCH_HH
#define MYMATCH_HH

#include "my_utils.hh"

#include <TMinuit.h>
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

class myMatch
{
    public:
        myMatch();
        myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, 
            double drift_speed, phot::PhotonVisibilityService const* PVS,  phot::SemiAnalyticalModel const* SAM);
        ~myMatch();

        int Nc;
        int Nf;
        int nc,nf;

        bool checkPossibility(const QCluster* qs,const  QFlash* qf);
        bool startFlash(const QCluster* qs,const  QFlash* qf);
        
        double drift_length;
        double drift_speed;

        int CH_MAX=80;

        TMinuit* MyMinuit = nullptr;
        int num_var = 1;
        void ChargeHypothesis(const double xoffset);
        double NLL(); //Negative Log-Likelihood
        static myMatch* s_me;
        static void FCN(Int_t&, Double_t*, Double_t& f, Double_t* x, Int_t); //FUNCAO DE MINIMAZAO PARA MINUIT

        std::vector<std::vector<double>> MYScore;
        std::vector<std::vector<double>> MYOffset;

        QFlash flash_actual;
        QCluster cluster_actual;

        QFlash flash_fit;
        QCluster cluster_fit;
        phot::PhotonVisibilityService const* fPVS;
        phot::SemiAnalyticalModel const* fSAM;

        double eff = 0.03;

        bool pause_now=false;

        


};


#endif