#include "MyMatch.hh"
#include "iostream"
#include <cmath>
#include <limits>
#include <algorithm>

myMatch::myMatch()
{}

myMatch::~myMatch()
{
    delete MyMinuit;
    MyMinuit = nullptr;
}

bool myMatch::checkPossibility(const QCluster* qs, const QFlash* qf)
{
    //if(fActualDetectorZone=="All") return true;
    if (!qs || !qf) return false;
    if (qs->empty()) return false;
    if (drift_speed <= 0.0 || drift_length <= 0.0) return true;

    // --- Janela de t0 usando x local do lado ---

    // --- Decide o lado do flash por PE (melhor que "xflashmedio" em ±L) ---
    double PE1=0.0,PE2=0.0,PE5=0.0,PE6=0.0;
    double peNeg = 0.0, pePos = 0.0;

    for (int i = 0; i < CH_MAX; ++i) 
    {
        const double pe = qf->PE_CH[i];
        if (pe <= 0.0) continue;
        if(i>=0 && i <40) PE6+=pe;
        else if(i>=40 && i<80) PE2+=pe;
        else if(i>=80 && i<120) PE5+=pe;
        else if(i>=120 && i<160) PE1+=pe;
    }

    pePos = PE6+PE2;
    peNeg = PE5+PE1;
    const double peTot = peNeg + pePos;
    if (peTot <= 0.0) return false;

    const int side = (pePos > peNeg) ? +1 : -1; // +1 -> +x, -1 -> -x

    int hits1=0,hits2=0,hits5=0,hits6=0;
    for (auto const& pt : *qs) 
    {
        const int thisAPA = (pt.APA);
        if(thisAPA==1) hits1++; 
        if(thisAPA==2) hits2++; 
        if(thisAPA==5) hits5++; 
        if(thisAPA==6) hits6++; 
    }

    if((hits6+hits2)>(hits1+hits5)) fActualDetectorZone="Positive";
    else if ((hits6 + hits2)<(hits1+hits5)) fActualDetectorZone="Negative";
    else fActualDetectorZone="Equal";

    // --- Consistência de lado (heurística) ---
    if (side > 0 && fActualDetectorZone == "Negative") return false;
    if (side < 0 && fActualDetectorZone == "Positive") return false;

    // --- Tempo (deixa explícito: t0 = time - offset_total) ---
    const double ttrigger     = -0; // aqui funciona para MC
    const double tflash_shift = 0; //1 --> -250.0; 2 --> 0 ; 3--> -500
    const double flash_time = qf->time - (ttrigger + tflash_shift); // ajuste sinais como você define offsets

    double xminP =  100000;
    double xmaxP = -100000;
    double xminN =  100000;
    double xmaxN = -100000;

    bool hasNeg=false;
    bool hasPos=false;

    for (auto const& pt : *qs) 
    {
        const double xlocal = (pt.x);   // se lado -, espelha
        if(pt.APA==6 || pt.APA==2)
        {
            xminP = std::min(xminP, xlocal);
            xmaxP = std::max(xmaxP, xlocal);
            hasPos=true;
        }
        if(pt.APA==1 || pt.APA==5)
        {
            xminN = std::min(xminN, xlocal);
            xmaxN = std::max(xmaxN, xlocal);
            hasNeg=true;
        }
        
    }

    const double v = drift_speed;
    const double L = drift_length;

    // 0 <= x_real = x_pandora + v*t0 <= L  (teu comentário)
    double tminP = 0.0;
    double tmaxP = 0.0;
    double tminN = 0.0;
    double tmaxN = 0.0;
    
    bool checkP=true;
    bool checkN=true;
    const double pad = 300.0; //300.0; // opcional: tolerância (ex.: 5 us)
    if(hasPos)
    {
        tminP = -xminP / v;
        tmaxP = (L - xmaxP) / v;
        if (tminP > tmaxP) return false;
        checkP = (flash_time >= tminP - pad) && (flash_time <= tmaxP + pad);
    }
    if(hasNeg)
    {
        tmaxN = (xminN+L)/v;
        tminN = xmaxN / v;
        if (tminN > tmaxN) return false;
        checkN = (flash_time >= tminN - pad) && (flash_time <= tmaxN + pad);
    }
    //std::cout << tmin << " < " << flash_time << " < " << tmax << std::endl;
    return checkP && checkN;
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

    double nll = 0.0;
    double eps = 1e-12;
    int chi=0;
    int chf=CH_MAX-1;

    if(normPE) flash_fit.norm_this_flash();

    if(type_fit=="flash")
    {
        for(int ch=chi;ch<=chf;++ch)
        {   
            double O = flash_actual.PE_CH[ch]; //OBSERVADO
            double H = flash_fit.PE_CH[ch]; //HIPOTESE
            if (H < eps) H = eps;

            // Poisson NLL (com gamma pra O double) --> se normalizado vira (cross entropy) (sum(H)=1)
            nll += (H - O * std::log(H)); //+ std::lgamma(O + 1.0); remover o gamma porque eh cte 
        }
    }
    else
    {
        for(int ch=chi;ch<=chf;++ch)
        {   
            double O = flash_actual.PE_CH[ch]; //OBSERVADO
            double H = flash_fit.PE_CH[ch]; //HIPOTESE
            if (H < eps) H = eps;

            // Poisson NLL (com gamma pra O double) --> se normalizado vira (cross entropy) (sum(H)=1)
            if(!normPE)
            {
                if (O < 0) O = 0;
                if (H < eps) H = eps;

                if (O > 0) nll += 2.0 * (H - O + O * std::log(O / H));
                else       nll += 2.0 * H;
            }         
        
                //nll += (H - O * std::log(H));// + std::lgamma(O + 1.0); //remover o gamma porque eh cte 
            else
                nll += (- O * std::log(H));
        }
    }
    
    return nll;
}


void myMatch::ChargeHypothesis(const double xoffset)
{
    //faz o deslocamento
    cluster_fit.resize(cluster_actual.size());
    double drift_time = 0.0;
    double atten_corr = 0.0;
    double dQdxCorr=0.0;
    double dQ=0.0;
    double dE=0.0;
    double dEdx = 0.0;
    double nphotons = 0.0;

    double W_ion = 23.6e-6;
    double W_Lar = 19.6e-6;
    double C1=0.212; // g/(MeV cm²)*kV/cm. 
    double B=C1/(density*Efield);
    double A=0.930; //
    for (size_t pt_index = 0; pt_index < cluster_actual.size(); ++pt_index) 
    {
        if(cluster_actual[pt_index].APA==2 ||cluster_actual[pt_index].APA==6 ) cluster_fit[pt_index].x = cluster_actual[pt_index].x + xoffset;
        else if(cluster_actual[pt_index].APA==1 ||cluster_actual[pt_index].APA==5 ) cluster_fit[pt_index].x = cluster_actual[pt_index].x - xoffset;

        cluster_fit[pt_index].y = cluster_actual[pt_index].y;
        cluster_fit[pt_index].z = cluster_actual[pt_index].z;

        //ESTIMAR A ENERGIA
        drift_time = (drift_length - std::abs(cluster_fit[pt_index].x)) / drift_speed;
        if (drift_time < 0.0) drift_time = 0.0; // fora de geometria / proteção
        atten_corr = std::exp(drift_time/elec_atenuation);
        dQdxCorr=cluster_actual[pt_index].q*atten_corr;
        
        dEdx = (std::exp(B * W_ion * dQdxCorr) - A) / B;

        dQ = dQdxCorr * cluster_actual[pt_index].pitch;
        dE = dEdx * cluster_actual[pt_index].pitch;
        nphotons = dE/W_Lar-dQ;
        cluster_fit[pt_index].q = std::max(0.0f, static_cast<float>(nphotons));

    }

    //------------------------------------------------------------
     //ESTIMAR A LUZ
    for(int ch=0;ch<CH_MAX;++ch)
    {
        this->flash_fit.PE_CH[ch]=0.0;
    }  

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
        int chf=CH_MAX-1;
   
        for(int ch=chi;ch<=chf;++ch)
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
}

bool myMatch::startFlash(const QCluster* qs,const  QFlash* qf)
{
    s_me = this;

    double x_maxP = -10e9;
    double x_minP = 10e9;
    double x_maxN = -10e9;
    double x_minN = 10e9;

    //pega o t0 candidato
    double t0=qf->time;
    double shiftt0=t0*drift_speed;

    bool hasP = false;
    bool hasN = false;

    //determina o valor minino e maximo em coordenadas do pandora
    QCluster track = *qs;
    for (auto const& pt : *qs)
    {   
        if(pt.APA==2 || pt.APA==6 )
        {
            hasP=true;
            if (pt.x > x_maxP) 
            { 
                x_maxP = pt.x; 
            }
            if (pt.x < x_minP) 
            {
                x_minP = pt.x; 
            }
        }
        if(pt.APA==1 || pt.APA==5 )
        {
            hasN=true;
            if (pt.x > x_maxN) 
            { 
                x_maxN = pt.x; 
            }
            if (pt.x < x_minN) 
            {
                x_minN = pt.x; 
            }
        }      
    }

    int side_shifted=0;
    //desloca para comecar em x=0;
    if(hasN==false)
    {
        for (auto &pt : track)
        {
            pt.x -= x_minP;
        }
    }
    else if(hasP==false)
    {
        for (auto &pt : track)
        {
            pt.x -= x_maxN;
        }
    }
    else // ambos positivos
    {
        if(x_maxN+x_minP<0)
        {
            side_shifted=+1;
            for (auto &pt : track)
            {
                if(pt.APA==2 ||pt.APA==6 )  pt.x -= x_minP;
                if(pt.APA==1 ||pt.APA==5 )  pt.x += x_minP;   
            }
        }
        else
        {
            side_shifted=-1;
            for (auto &pt : track)
            {
                if(pt.APA==2 ||pt.APA==6 )  pt.x += x_maxN;
                if(pt.APA==1 ||pt.APA==5 )  pt.x -= x_maxN;  
            }
        }
    }
    
    double x_maxAbs = -10e9;
    for (auto const& pt : track)
    {   
        if (abs(pt.x) > x_maxAbs) 
        { 
            x_maxAbs = abs(pt.x); 
        }    
    }

    this->flash_actual = *qf;
    this->cluster_actual = track;
    //-------
    //arruma minuit
    if (!MyMinuit)
    {
        MyMinuit = new TMinuit(this->num_var);
    } 
    // fazer configuracoes do minuit
    // 

    //const double deltaXP = x_maxP - x_minP; //tamanho do track
    //const double deltaXN = x_maxP - x_minP; 
    const double step = 0.1;
    double xfitmin = 0;
    double xfitmax = this->drift_length-x_maxAbs;

    double x0=0.0;
    double x_tolerance = 5;
  
    if(hasN==false)
    {
        x0 = x_minP + shiftt0;
    }
    else if(hasP==false)
    {
        x0 = x_maxN - shiftt0;
    }
    else
    {
        if(side_shifted==+1)
        {
            x0 = x_minP + shiftt0;
        }
        if(side_shifted==-1)
        {
            x0 = x_maxN - shiftt0;
        }   
    }

    x0=abs(x0);

    //-----------------------------------
    double bestx=x0, bestxerr=0;
    if(fit_mode)
    {
        //x0 = (xfitmax+xfitmin)/2.0;
        xfitmin -= x_tolerance; 
        xfitmax += x_tolerance; 

        int ierr = 0;
        MyMinuit->mnexcm("CLEAR", nullptr, 0, ierr);
        MyMinuit->SetPrintLevel(-1);
        MyMinuit->SetFCN(myMatch::FCN);
        
        MyMinuit->DefineParameter(0, "Xoffset", x0, step, xfitmin, xfitmax);
        double arglist[2] = {5000, 0.01};//{5000, 0.01};
        MyMinuit->mnexcm("MIGRAD", arglist, 2, ierr);
        MyMinuit->GetParameter(0, bestx, bestxerr);

        //if(ierr!=0) return false;
    }
     
    if(hasN==false || side_shifted==+1 ) 
    {
        shiftt0 = bestx - x_minP;
    }
    else if(hasP==false || side_shifted==-1 ) 
    {
        shiftt0 = x_maxN - (-bestx);
    }

    ChargeHypothesis(bestx);    

    //buscar o ponto mais proximo do anodo...
    double x_close_anode;
    if(hasN==false) 
    {
        x_close_anode = x_maxP + shiftt0;
    }
    else if(hasP==false)
    {
        x_close_anode = x_minN - shiftt0;
    }
    else
    {
        double xcp = x_maxP + shiftt0;
        double xcn = x_minN - shiftt0;
        if(xcp>abs(xcn))
        {
            x_close_anode = xcp;
        }
        else
        {
            x_close_anode = xcn;
        }
    }
    
    this->MYScore[nf][nc] = NLL();
    this->MYOffset[nf][nc] = bestx;
    this->MYdeltaT0[nf][nc] = shiftt0;  
    this->MYcloseAnode[nf][nc] = abs(x_close_anode);  
    this->MYvisEf[nf][nc] = returnVisEff(qs,shiftt0);

    return true;
}

myMatch::myMatch(std::vector<QCluster> qqs ,std::vector<QFlash> qfs, double drift_length, 
    double drift_speed, double elec_atenuation, double density, double Efield,
     phot::PhotonVisibilityService const* PVS, phot::SemiAnalyticalModel const* SAM, bool norm, std::string DetectorZone, std::string typeFit, bool fit_mode)
{
    this->drift_length = drift_length;
    this->drift_speed = drift_speed;
    this->elec_atenuation = elec_atenuation;
    this->density = density;
    this->Efield = Efield;

    this->fPVS = PVS;
    this->fSAM = SAM;

    this->normPE = norm; 
    this->fDetectorZone = DetectorZone;
    
    this->type_fit = typeFit;
    this->fit_mode = fit_mode;

    //ZERA o vetor de fit
    flash_fit.PE_CH.clear();
    flash_fit.PE_CH.resize(fPVS->NOpChannels());

    this->Nc = qqs.size();
    this->Nf = qfs.size();

    this->Nline = this->Nf;
    this->Ncol = this->Nc;
   
    const double BIG = 1e8; 
    MYScore.assign(Nline, std::vector<double>(Ncol, BIG));
    MYOffset.assign(Nline, std::vector<double>(Ncol, -100000.0));
    MYdeltaT0.assign(Nline, std::vector<double>(Ncol, -100000.0));
    MYcloseAnode.assign(Nline, std::vector<double>(Ncol, -100000.0));
    MYvisEf.assign(Nline, std::vector<double>(Ncol, -100000.0));

    this->CH_MAX = fPVS->NOpChannels();
    direct_visibilities.assign(CH_MAX, 0.0);

    //talvez colocar um algoritmo de filtro para filtrar cluster e flashs? --> tipo isso usado no codigo no SBND _alg_tpc_filter->Filter(_tpc_object_v);
    for(nf=0;nf<Nf;nf++)
    {
        for(nc=0;nc<Nc;nc++)
        {
            //AQUI testamos se essa dupla eh possivel
            if(this->checkPossibility(&qqs[nc],&qfs[nf]))
            {
                //testa o par
                /* std::cout << nf << " -- " << nc << std::endl;  */
                startFlash(&qqs[nc],&qfs[nf]);
            }
            
        }
    }

    this->HR =  hungarian_min(MYScore); // aqui tava comentado, acho que era por isso
  
}

double myMatch::returnVisEff(const QCluster* qs, const double xoffset)
{
    if (!qs) return 100000.0;

    double total = 0.0;
    double L = qs->Length;
    if (L <= 0.0) return 100000.0;

    for (size_t pt_index = 0; pt_index < qs->size(); ++pt_index)
    {
        const auto& pt = (*qs)[pt_index];

        double x = pt.x;
        double y = pt.y;
        double z = pt.z;
        double pitch = (pt.pitch >= 0.0) ? pt.pitch : 0.0;

        if (pt.APA == 2 || pt.APA == 6)
            x = pt.x + xoffset;
        else if (pt.APA == 1 || pt.APA == 5)
            x = pt.x - xoffset;
        else
            continue;

        geo::Point_t point{x, y, z};
        fSAM->detectedDirectVisibilities(direct_visibilities, point);

        for (int ch = 0; ch < CH_MAX; ++ch)
        {
            total += direct_visibilities[ch] * pitch;
        }
    }

    return total / L;
}