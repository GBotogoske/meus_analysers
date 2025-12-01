#include "MyMatch.hh"
#include "iostream"


myMatch::myMatch()
{}

myMatch::~myMatch()
{}

bool myMatch::checkPossibility(const QCluster* qs,const  QFlash* qf)
{
    auto flash_time = qf->time;
    double x_max = -10e9;
    double x_min = 10e9;
    for (auto const& pt : *qs)
    {
      if (pt.x > x_max) { x_max = pt.x; }
      if (pt.x < x_min) { x_min = pt.x; }
    }
    //std::cout << "Flash Time :=> " << flash_time << std::endl;

    //DEPOIS PARA DADOS TALVEZ TENHA QUE USAR OFFSETS ...
    double clus_t_min = (x_max - this->drift_length) /  this->drift_speed ;
    double clus_t_max = x_min /  this->drift_speed ;
    //std::cout << clus_t_min << " - - - " << clus_t_max << std::endl;
    //std::cout  << " - - - - - - - - - - - - - - - - - - - - - - - - "<< std::endl;

    return (clus_t_min  < flash_time && flash_time < clus_t_max );
}

const QFlash myMatch::ChargeHypothesis(const double xoffset)
{
    for (size_t pt_index = 0; pt_index < cluster_actual.size(); ++pt_index) 
    {
      cluster_fit[pt_index].x = cluster_actual[pt_index].x + xoffset;
      cluster_fit[pt_index].y = cluster_actual[pt_index].y;
      cluster_fit[pt_index].z = cluster_actual[pt_index].z;
      cluster_fit[pt_index].q = cluster_actual[pt_index].q;
    }
    return QFlash();
}

bool myMatch::startFlash(const QCluster* qs,const  QFlash* qf)
{
    double x_max = -10e9;
    double x_min = 10e9;

    QCluster track = *qs;
    for (auto const& pt : *qs)
    {
        if (pt.x > x_max) 
        { 
            x_max = pt.x; 
        }
        if (pt.x < x_min) 
        {
            x_min = pt.x; 
        }
    }

    for (auto &pt : track)
    {
        pt.x -= x_min;
    } 
    this->flash_actual = *qf;
    this->cluster_actual = *qs;
    //talvez aqui normalizar o flash
    //-------

    if (!MyMinuit)
    {
        MyMinuit = new TMinuit(this->num_var);
    } 
    // fazer configuracoes do minuit
    // 


    return true;
}

myMatch::myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, double drift_speed)
{
    this->drift_length = drift_length;
    this->drift_speed = drift_speed;

    int Nc = qqs.size();
    int Nf = qfs.size();

    //talvez colocar um algoritmo de filtro para filtrar cluster e flashs? --> tipo isso usado no codigo no SBND _alg_tpc_filter->Filter(_tpc_object_v);
    
    for(int nc=0;nc<Nc;nc++)
    {
        for(int nf=0;nf<Nf;nf++)
        {
            //AQUI testamos se essa dupla eh possivel
            if(this->checkPossibility(&qqs[nc],&qfs[nf]))
            {
                //testa o par
                startFlash(&qqs[nc],&qfs[nf]);
            }
            
        }
    }
   


}