#include "MyMatch.hh"
#include "iostream"

myMatch::myMatch()
{}

myMatch::~myMatch()
{}

bool myMatch::checkPossibility(const QCluster* qs,const  QFlash* qf)
{
    auto flash_time = qf->time;
    double x_max = -10e3;
    double x_min = 10e3;
    for (auto const& pt : *qs)
    {
      if (pt.x > x_max) { x_max = pt.x; }
      if (pt.x < x_min) { x_min = pt.x; }
    }
    std::cout << flash_time << std::endl;
    return true;
}

myMatch::myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs)
{
    int Nc = qqs.size();
    int Nf = qfs.size();

    //talvez colocar um algoritmo de filtro para filtrar cluster e flashs? --> tipo isso usado no codigo no SBND _alg_tpc_filter->Filter(_tpc_object_v);
    

    for(int nc=0;nc<Nc;nc++)
    {
        for(int nf=0;nf<Nf;nf++)
        {
            //AQUI testamos se essa dupla eh possivel
            this->checkPossibility(&qqs[nc],&qfs[nf]);
            
        }
    }
   


}