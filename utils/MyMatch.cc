#include "MyMatch.hh"
#include "iostream"

myMatch::myMatch()
{}

myMatch::~myMatch()
{
    delete MyMinuit;
    MyMinuit = nullptr;
}

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

myMatch* myMatch::s_me = nullptr;

void myMatch::FCN(Int_t&, Double_t*, Double_t& f, Double_t* x, Int_t)
{
  // x[0] = xoffset
  s_me->ChargeHypothesis(x[0]);
  f = s_me->NLL();
}


double myMatch::NLL()
{
    
    if(pause_now)
    {
        std::cout << "REAL FLASH" << std::endl;
        for(int ch=0;ch<CH_MAX;++ch)
        {
            std::cout << "CH: " << ch << " => " << this->flash_actual.PE_CH[ch] << std::endl;
        }   
        std::cout << "Press ENTER." << std::endl;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    double nll = 0.0;
    double eps = 1e-6;
    if(normPE) flash_fit.norm_this_flash();
    for (int ch = 0; ch < CH_MAX; ++ch)
    {
        double O = flash_actual.PE_CH[ch]; //OBSERVADO
        double H = flash_fit.PE_CH[ch]; //HIPOTESE

        // thresholds pra não explodir log(0)
        if (H < eps) H = eps;
        // Poisson NLL (com gamma pra O double)
        nll += (H - O * std::log(H));// + std::lgamma(O + 1.0));
    }
    return nll;

}

void myMatch::ChargeHypothesis(const double xoffset)
{
    //faz o deslocamento
    cluster_fit.resize(cluster_actual.size());
    for (size_t pt_index = 0; pt_index < cluster_actual.size(); ++pt_index) 
    {
        cluster_fit[pt_index].x = cluster_actual[pt_index].x + xoffset;
        cluster_fit[pt_index].y = cluster_actual[pt_index].y;
        cluster_fit[pt_index].z = cluster_actual[pt_index].z;
        cluster_fit[pt_index].q = cluster_actual[pt_index].q;
    }

    //vamos trabalahos com apenas os primerios 2 APAS por hora (canal maximo 0-79)
    
    for(int ch=0;ch<CH_MAX;++ch)
    {
        this->flash_fit.PE_CH[ch]=0.0;
    }  

    std::vector<double> direct_visibilities(fPVS->NOpChannels(), 0.0);
    direct_visibilities.reserve(fPVS->NOpChannels());
    //varre o cluster
    for(size_t i=0;i<this->cluster_fit.size();i++)
    {   
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double q = 0.0;

        x = this->cluster_fit[i].x ;
        y = this->cluster_fit[i].y ;
        z = this->cluster_fit[i].z ;
        q = this->cluster_fit[i].q ;
        geo::Point_t point{x, y, z};

        //calcula a visibilidade
        
        fSAM->detectedDirectVisibilities(direct_visibilities, point);
       
        int chi=0;
        int chf=CH_MAX;

        if(fDetectorZone=="Positive")
        {
            chf=80;
        }
        if(fDetectorZone=="Negative")
        {
            chi=80;
        }
        for(int ch=chi;ch<chf;++ch)
        {   
            //std::cout << "i := " << i << std::endl;  
            //std::vector<double> reflected_visibilities;
            //fSAM->detectedDirectVisibilities(reflected_visibilities, point);
            //double v1 = fPVS->GetVisibility(point, ch);
            //double v2 = fPVS->GetVisibility(point, ch , true);
            //double v3 = direct_visibilities[ch];
            //double v4 = reflected_visibilities[ch];
            double vis = direct_visibilities[ch];
            double n1 = q*vis*eff;
            this->flash_fit.PE_CH[ch]+=n1;
        }     
    }

    if(pause_now)
    {   
         std::cout << "HYPOTESIS FLASH" << std::endl;
        for(int ch=0;ch<CH_MAX;++ch)
        {
            std::cout << "CH: " << ch << " => " << this->flash_fit.PE_CH[ch] << std::endl;
        }   
        std::cout << "Press ENTER." << std::endl;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    //return QFlash();
}

bool myMatch::startFlash(const QCluster* qs,const  QFlash* qf)
{
    s_me = this;

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
    this->cluster_actual = track;
    //talvez aqui normalizar o flash
    //-------

    if (!MyMinuit)
    {
        MyMinuit = new TMinuit(this->num_var);
    } 
    // fazer configuracoes do minuit
    // 

    const double deltaX = x_max - x_min;
    const double step = 0.5;
    double xfitmin = -this->drift_length;
    double xfitmax = this->drift_length-deltaX;
    if(fDetectorZone=="Positive")
    {
        xfitmin = 0;
    }
     if(fDetectorZone=="Negative")
    {
        xfitmax = 0-deltaX;
    }
    
    const double x0 = (xfitmax-xfitmin)/2.0;
    int ierr = 0;
    MyMinuit->mnexcm("CLEAR", nullptr, 0, ierr);
    MyMinuit->SetPrintLevel(-1);
    MyMinuit->SetFCN(myMatch::FCN);
    
    MyMinuit->DefineParameter(0, "Xoffset", x0, step, xfitmin, xfitmax);
    double arglist[2] = {5000, 0.01};
    MyMinuit->mnexcm("MIGRAD", arglist, 2, ierr);
    
    double bestx=0, bestxerr=0;
    MyMinuit->GetParameter(0, bestx, bestxerr);
    ChargeHypothesis(bestx);          // hipótese final
    this->MYScore[nf][nc] = NLL();
    this->MYOffset[nf][nc] = bestx; 

    return true;
}

myMatch::myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, 
    double drift_speed,phot::PhotonVisibilityService const* PVS, phot::SemiAnalyticalModel const* SAM, bool norm, std::string DetectorZone)
{
    this->drift_length = drift_length;
    this->drift_speed = drift_speed;
    this->fPVS = PVS;
    this->fSAM = SAM;

    this->normPE = norm; 
    this->fDetectorZone = DetectorZone;
    
    //ZERA o vetor de fit
    flash_fit.PE_CH.clear();
    flash_fit.PE_CH.resize(fPVS->NOpChannels());

    this->Nc = qqs.size();
    this->Nf = qfs.size();

    const double BIG = 1e15; 
    MYScore.assign(Nf, std::vector<double>(Nc, BIG));
    MYOffset.assign(Nf, std::vector<double>(Nc, -1.0));

    this->CH_MAX = fPVS->NOpChannels();

    //talvez colocar um algoritmo de filtro para filtrar cluster e flashs? --> tipo isso usado no codigo no SBND _alg_tpc_filter->Filter(_tpc_object_v);
    for(nf=0;nf<Nf;nf++)
    {
        for(nc=0;nc<Nc;nc++)
        {
            //AQUI testamos se essa dupla eh possivel
            if(true)//this->checkPossibility(&qqs[nc],&qfs[nf]))
            {
                //testa o par
                std::cout << nf << " -- " << nc << std::endl; 
                startFlash(&qqs[nc],&qfs[nf]);
            }
            
        }
    }

    this->HR =  hungarian_min(MYScore);
    std::cout << "Custo minimo = " << this->HR.cost << "\n";

    const int nRows = std::min((int)this->HR.assign.size(), this->Nf);

    for (int nf = 0; nf < nRows; ++nf) 
    {
        const int nc = this->HR.assign[nf];

        // se caiu em dummy / não atribuído
        if (nc < 0 || nc >= this->Nc) 
        {
            std::cout << "Flash row " << nf
                    << " (FlashID=" << qfs[nf].flashID << ")"
                    << " -> unassigned/dummy (col=" << nc << ")\n";
            continue;
        }

        std::cout << "Flash row " << nf
                << " (FlashID=" << qfs[nf].flashID << ")"
                << " -> Cluster col " << nc
                << " (ClusterID=" << qqs[nc].objID << ")"
                << " | cost=" << this->MYScore[nf][nc]
                << " | xoffset=" << this->MYOffset[nf][nc] << " cm\n";
    }

   
}