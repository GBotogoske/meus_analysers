#ifndef MYMATCH_HH
#define MYMATCH_HH

#include "my_utils.hh"

class myMatch
{
    public:
        myMatch();
        myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs);
        ~myMatch();

        bool checkPossibility(const QCluster* qs,const  QFlash* qf);
};


#endif