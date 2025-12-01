#ifndef MYMATCH_HH
#define MYMATCH_HH

#include "my_utils.hh"

#include <TMinuit.h>

class myMatch
{
    public:
        myMatch();
        myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, double drift_speed);
        ~myMatch();

        bool checkPossibility(const QCluster* qs,const  QFlash* qf);
        bool startFlash(const QCluster* qs,const  QFlash* qf);
        
        double drift_length;
        double drift_speed;

        TMinuit* MyMinuit;
        int num_var = 4;
        const QFlash ChargeHypothesis(const double xoffset);

        QFlash flash_actual;
        QCluster cluster_actual;

        QFlash flash_fit;
        QCluster cluster_fit;

};


#endif