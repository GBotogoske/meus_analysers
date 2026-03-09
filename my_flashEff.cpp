#include <string>
#include <vector>
#include <limits>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <cstddef>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>

struct TrackKey
{
    int run=0, event=0, track=0, type=0;
    bool operator==(TrackKey const& o) const { return run==o.run && event==o.event && track==o.track && type==o.type; }
};

struct FlashKey
{
    int run=0, event=0, flash=0;
    bool operator==(FlashKey const& o) const { return run==o.run && event==o.event && flash==o.flash; }
};

static inline void hash_combine(std::size_t& seed, std::size_t v) noexcept {
    // versão padrão (boost-like)
    seed ^= v + 0x9e3779b97f4a7c15ULL + (seed<<6) + (seed>>2);
}

struct TrackKeyHash 
{
    std::size_t operator()(TrackKey const& k) const noexcept 
    {
        std::size_t seed = 0;
        hash_combine(seed, std::hash<int>{}(k.run));
        hash_combine(seed, std::hash<int>{}(k.event));
        hash_combine(seed, std::hash<int>{}(k.track));
        hash_combine(seed, std::hash<int>{}(k.type));
        return seed;
    }
};

struct FlashKeyHash 
{
    std::size_t operator()(FlashKey const& k) const noexcept
     {
        std::size_t seed = 0;
        hash_combine(seed, std::hash<int>{}(k.run));
        hash_combine(seed, std::hash<int>{}(k.event));
        hash_combine(seed, std::hash<int>{}(k.flash));
        return seed;
    }
};

struct MCInfo
{
    double score = std::numeric_limits<double>::quiet_NaN();
    bool hasScore = false;
    int matchedPDG = -1;
};

struct MatchInfo
{
    double score = std::numeric_limits<double>::quiet_NaN();
    bool hasScore = false;
    double flashTime;
    double deltaT0;
    double Length=-50e3;
    bool isCorrect = false;
    bool isTrustable = false;
    double X=-50e3;
    double light=-50e3;
    double charge=-50e3;

    double scoreMC=-1; 

};

void fill_MCMap(std::unordered_map<FlashKey,std::unordered_map<TrackKey, MCInfo, TrackKeyHash>,FlashKeyHash> &flashMapMC, std::unordered_map<TrackKey,std::unordered_map<FlashKey, MCInfo, FlashKeyHash>,TrackKeyHash> &trackMapMC,
                TFile* fMC  ,std::string treeMC_name , bool only_maximum_MC_score)
{
    int run, event;
    int TrackType, TrackID;
    int FlashID;
    float jaccard;
    int domPDG;

    TTree* treeMC = dynamic_cast<TTree*>(fMC->Get(treeMC_name.c_str()));
    treeMC->SetBranchAddress("run",&run);
    treeMC->SetBranchAddress("event",&event);
    treeMC->SetBranchAddress("trackKey",&TrackID);
    treeMC->SetBranchAddress("type",&TrackType);
    treeMC->SetBranchAddress("flashKey",&FlashID);
    treeMC->SetBranchAddress("jaccardW",&jaccard);
    treeMC->SetBranchAddress("domPDG",&domPDG);
    int nMC = treeMC->GetEntries();

    for(int i=0;i<nMC;i++)
    {
        treeMC->GetEntry(i);
        FlashKey fk{run, event, FlashID};
        TrackKey tk{run, event, TrackID, TrackType};
        MCInfo info;
        info.score = jaccard;
        info.hasScore = true;
        info.matchedPDG=domPDG;
    
        flashMapMC[fk][tk]=info;
        trackMapMC[tk][fk]=info; 
    }

    //------------------------------------------------------------------------------------------
    if(only_maximum_MC_score)
    {
        // Mantém apenas o TrackKey com maior score dentro de cada flashMapMC[fk]
        for (auto & [fk, inner] : flashMapMC) 
        {
            if (inner.size() <= 1) continue;

            auto bestIt = inner.end();
            double best = -std::numeric_limits<double>::infinity();

            for (auto it = inner.begin(); it != inner.end(); ++it)
            {
                const MCInfo &info = it->second;
                if (!info.hasScore) continue;
                double s = info.score;
                if (bestIt == inner.end() || s > best) 
                {
                    best = s;
                    bestIt = it;
                }
            }
            if (bestIt == inner.end()) 
            {
                inner.clear(); // não achou nenhum score válido
            } 
            else 
            {
                auto keep = *bestIt;  // copia (TrackKey, MCInfo)
                inner.clear();
                inner.emplace(keep.first, keep.second);
            }
        }
        // Mantém apenas o FlashKey com maior score dentro de cada trackMapMC[tk]
        for (auto & [tk, inner] : trackMapMC) 
        {
            if (inner.size() <= 1) continue;
            auto bestIt = inner.end();
            double best = -std::numeric_limits<double>::infinity();
            for (auto it = inner.begin(); it != inner.end(); ++it) 
            {
                const MCInfo &info = it->second;
                if (!info.hasScore) continue;
                double s = info.score;
                if (bestIt == inner.end() || s > best) 
                {
                    best = s;
                    bestIt = it;
                }
            }
            if (bestIt == inner.end()) 
            {
                inner.clear(); // nenhum score válido
            } 
            else 
            {
                auto keep = *bestIt; // copia (FlashKey, MCInfo)
                inner.clear();
                inner.emplace(keep.first, keep.second);
            }
        }
    }
}

void fill_DataFMap(std::unordered_map<FlashKey,std::unordered_map<TrackKey, MatchInfo, TrackKeyHash>,FlashKeyHash>& flashMapData_Hungarian, std::unordered_map<FlashKey,std::unordered_map<TrackKey, MatchInfo, TrackKeyHash>,FlashKeyHash>& flashMapData_Hungry,
                    TFile* fData,  std::string treeF_name,std::unordered_map<FlashKey,std::unordered_map<TrackKey, MCInfo, TrackKeyHash>,FlashKeyHash> flashMapMC)
{
    int run,event;
    int FlashID_Ftype;
    int clusterID_Hungarian, clusterID_Hungry;
    int clusterType_Hungarian, clusterType_Hungry;
    double flashTime_Ftype;
    double score_Hungarian,score_Hungry;
    double deltaT0_Hungarian,deltaT0_Hungry;
    double length_Hungarian,length_Hungry;
    double x_Hungarian , x_Hungry;
    double charge_Hungarian,charge_Hungry;
    double light;

    TTree* treeF = dynamic_cast<TTree*>(fData->Get(treeF_name.c_str()));
    treeF->SetBranchAddress("run",&run);
    treeF->SetBranchAddress("event",&event);
    treeF->SetBranchAddress("flashID",&FlashID_Ftype);
    treeF->SetBranchAddress("clusterID",&clusterID_Hungarian);
    treeF->SetBranchAddress("clusterID_trivial",&clusterID_Hungry);
    treeF->SetBranchAddress("clusterType",&clusterType_Hungarian);
    treeF->SetBranchAddress("clusterType_trivial",&clusterType_Hungry);
    treeF->SetBranchAddress("flashTime",&flashTime_Ftype);
    treeF->SetBranchAddress("score",&score_Hungarian);
    treeF->SetBranchAddress("score_trivial",&score_Hungry);
    treeF->SetBranchAddress("deltat0",&deltaT0_Hungarian);
    treeF->SetBranchAddress("deltat0_trivial",&deltaT0_Hungry);
    treeF->SetBranchAddress("length",&length_Hungarian);
    treeF->SetBranchAddress("length_trivial",&length_Hungry);
    treeF->SetBranchAddress("x_trivial",&x_Hungry);
    treeF->SetBranchAddress("x",&x_Hungarian);
    treeF->SetBranchAddress("charge_trivial",&charge_Hungry);
    treeF->SetBranchAddress("charge",&charge_Hungarian);
    treeF->SetBranchAddress("light",&light);
    int nF = treeF->GetEntries();

    int i1f=0,i2f=0,i3f=0;
    for(int i=0;i<nF;i++)
    {
        treeF->GetEntry(i);
        FlashKey fk{run, event, FlashID_Ftype};
        TrackKey tk_Hungarian{run, event, clusterID_Hungarian, clusterType_Hungarian};
        TrackKey tk_Hungry{run, event, clusterID_Hungry, clusterType_Hungry};
        MatchInfo infoHungarian;
        MatchInfo infoHungry;

        infoHungarian.score = score_Hungarian;
        infoHungarian.deltaT0 = deltaT0_Hungarian;
        infoHungarian.flashTime = flashTime_Ftype;
        infoHungarian.hasScore = true;
        infoHungarian.Length = length_Hungarian;
        infoHungarian.X = x_Hungarian;
        infoHungarian.charge = charge_Hungarian;
        infoHungarian.light = light;

        infoHungry.score = score_Hungry;
        infoHungry.deltaT0 = deltaT0_Hungry;
        infoHungry.flashTime = flashTime_Ftype;
        infoHungry.hasScore = true;
        infoHungry.Length = length_Hungry;
        infoHungry.X = x_Hungry;
        infoHungry.charge = charge_Hungry;
        infoHungry.light = light;

        auto itF = flashMapMC.find(fk);
        if (itF != flashMapMC.end()) 
        {
            infoHungarian.isTrustable = true;
            auto itT = itF->second.find(tk_Hungarian);
            infoHungarian.isCorrect   = itT != itF->second.end();
            if (itT != itF->second.end())
            {
                infoHungarian.isCorrect = true;
                infoHungarian.scoreMC   = itT->second.score;   
            }

            infoHungry.isTrustable = true;
            itT = itF->second.find(tk_Hungry);
            infoHungry.isCorrect   = (itT != itF->second.end());
            if (itT != itF->second.end())
            {
                infoHungry.isCorrect = true;
                infoHungry.scoreMC   = itT->second.score;   
            }

            if(infoHungarian.isCorrect ) i1f++;
            if(infoHungry.isCorrect ) i2f++;
            i3f++;
        } 
        else 
        {
            infoHungarian.isTrustable = false;
            infoHungarian.isCorrect   = (infoHungarian.score >= 1e6);

            infoHungry.isTrustable = false;
            infoHungry.isCorrect   = (infoHungry.score >= 1e6);
        }
        flashMapData_Hungarian[fk][tk_Hungarian] = infoHungarian;
        flashMapData_Hungry[fk][tk_Hungry] = infoHungry;
    } 
    std::cout << i1f << " " << i2f << " " << i3f <<std::endl;
    std::cout << (double) i1f/i3f << "  " << (double) i2f/i3f << std::endl;
}

void fill_DataTMap(std::unordered_map<TrackKey,std::unordered_map<FlashKey, MatchInfo, FlashKeyHash>,TrackKeyHash>& trackMapData_Hungarian, std::unordered_map<TrackKey,std::unordered_map<FlashKey, MatchInfo, FlashKeyHash>,TrackKeyHash>& trackMapData_Hungry,
                    TFile* fData,  std::string treeT_name, std::unordered_map<TrackKey,std::unordered_map<FlashKey, MCInfo, FlashKeyHash>,TrackKeyHash> trackMapMC)
{
    int run,event;
    int clusterID_Ttype;
    int clusterType_Ttype;
    int FlashID_Hungry,FlashID_Hungarian;
    double flashTime_Hungry,flashTime_Hungarian;
    double score_Hungarian,score_Hungry;
    double deltaT0_Hungarian,deltaT0_Hungry;
    double x_Hungarian , x_Hungry;
    double length;
    double light_Hungarian,light_Hungry;
    double charge;
   
    TTree* treeT = dynamic_cast<TTree*>(fData->Get(treeT_name.c_str()));
    treeT->SetBranchAddress("run",&run);
    treeT->SetBranchAddress("event",&event);
    treeT->SetBranchAddress("clusterID",&clusterID_Ttype);
    treeT->SetBranchAddress("clusterType",&clusterType_Ttype);
    treeT->SetBranchAddress("flashID",&FlashID_Hungarian);
    treeT->SetBranchAddress("flashID_trivial",&FlashID_Hungry);
    treeT->SetBranchAddress("flashTime",&flashTime_Hungarian);
    treeT->SetBranchAddress("flashTime_trivial",&flashTime_Hungry);
    treeT->SetBranchAddress("score",&score_Hungarian);
    treeT->SetBranchAddress("score_trivial",&score_Hungry);
    treeT->SetBranchAddress("deltat0",&deltaT0_Hungarian);
    treeT->SetBranchAddress("deltat0_trivial",&deltaT0_Hungry);
    treeT->SetBranchAddress("length",&length);
    treeT->SetBranchAddress("x_trivial",&x_Hungry);
    treeT->SetBranchAddress("x",&x_Hungarian);
    treeT->SetBranchAddress("light_trivial",&light_Hungry);
    treeT->SetBranchAddress("light",&light_Hungarian);
    treeT->SetBranchAddress("charge",&charge);

    int nT = treeT->GetEntries();

    int i1t=0,i2t=0,i3t=0;
    for(int i=0;i<nT;i++)
    {
        treeT->GetEntry(i);
        FlashKey fk_Hungarian{run, event, FlashID_Hungarian};
        FlashKey fk_Hungry{run, event, FlashID_Hungry};
        TrackKey tk{run, event, clusterID_Ttype, clusterType_Ttype};

        MatchInfo infoHungarian;
        MatchInfo infoHungry;

        infoHungarian.score = score_Hungarian;
        infoHungarian.deltaT0 = deltaT0_Hungarian;
        infoHungarian.flashTime = flashTime_Hungarian;
        infoHungarian.hasScore = true;
        infoHungarian.Length = length;
        infoHungarian.X = x_Hungarian;
        infoHungarian.charge = charge;
        infoHungarian.light = light_Hungarian;

        infoHungry.score = score_Hungry;
        infoHungry.deltaT0 = deltaT0_Hungry;
        infoHungry.flashTime = flashTime_Hungry;
        infoHungry.hasScore = true;
        infoHungry.Length = length;
        infoHungry.X = x_Hungry;
        infoHungry.charge = charge;
        infoHungry.light = light_Hungry;

        auto itT = trackMapMC.find(tk);
        if (itT != trackMapMC.end()) 
        {
            infoHungarian.isTrustable = true;
            auto itF = itT->second.find(fk_Hungarian);
            if (itF != itT->second.end())
            {
                infoHungarian.isCorrect = true;
                infoHungarian.scoreMC   = itF->second.score;   
            }
            infoHungry.isTrustable = true;
            itF = itT->second.find(fk_Hungry);
            if (itF != itT->second.end())
            {
                infoHungry.isCorrect = true;
                infoHungry.scoreMC   = itF->second.score;   
            }

            if(infoHungarian.isCorrect ) i1t++;
            if(infoHungry.isCorrect ) i2t++;
            i3t++;
        } 
        else 
        {
            infoHungarian.isTrustable = false;
            infoHungarian.isCorrect   = (infoHungarian.score >= 1e6);

            infoHungry.isTrustable = false;
            infoHungry.isCorrect   = (infoHungry.score >= 1e6);
        }
        trackMapData_Hungarian[tk][fk_Hungarian] = infoHungarian;
        trackMapData_Hungry[tk][fk_Hungry] = infoHungry;
    
    }
    std::cout << i1t << " " << i2t << " " << i3t <<std::endl;
    std::cout << (double) i1t/i3t << "  " << (double) i2t/i3t << std::endl;
}

void save_FMap(std::unordered_map<FlashKey,std::unordered_map<TrackKey, MatchInfo, TrackKeyHash>,FlashKeyHash> fmap,std::string treeName)
{  
    int run_o,event_o;
    int trackID_o,type_o, flashID_o;
    double score_o,deltato_o,flashTime_o,scoreMC_o;
    int isCorrect_o;
    double length;
    double x;
    double charge,light;

    TTree tout(treeName.c_str(),treeName.c_str());
    tout.Branch("run",&run_o,"run/I");
    tout.Branch("event",&event_o,"event/I");
    tout.Branch("trackID",&trackID_o,"trackID/I");
    tout.Branch("type",&type_o,"type/I");
    tout.Branch("flashID",&flashID_o,"flashID/I");
    tout.Branch("scoreData",&score_o,"scoreData/D");
    tout.Branch("deltaT0",&deltato_o,"deltaT0/D");
    tout.Branch("flashTime",&flashTime_o,"flashTime/D");
    tout.Branch("scoreMC",&scoreMC_o,"scoreMC/D");
    tout.Branch("isCorrect",&isCorrect_o,"isCorrect/I");
    tout.Branch("length",&length,"length/D");
    tout.Branch("x",&x,"x/D");
    tout.Branch("charge",&charge,"charge/D");
    tout.Branch("light",&light,"light/D");
    
    for(auto thisInfo:fmap)
    {
        run_o=thisInfo.first.run;
        event_o=thisInfo.first.event;
        flashID_o=thisInfo.first.flash;
        for(auto thisInfo2:thisInfo.second)
        {
            trackID_o=thisInfo2.first.track;
            type_o=thisInfo2.first.type;
            score_o=thisInfo2.second.score;
            deltato_o=thisInfo2.second.deltaT0;
            flashTime_o=thisInfo2.second.flashTime;
            scoreMC_o=thisInfo2.second.scoreMC;
            length = thisInfo2.second.Length;
            x = thisInfo2.second.X;
            charge = thisInfo2.second.charge;
            light = thisInfo2.second.light;
            if(thisInfo2.second.isCorrect && thisInfo2.second.isTrustable) isCorrect_o=1; //correto
            else if(!thisInfo2.second.isCorrect && thisInfo2.second.isTrustable) isCorrect_o=0; // errado
            else if(thisInfo2.second.isCorrect && !thisInfo2.second.isTrustable) isCorrect_o=-1; //correto nao tem match real
            else if(!thisInfo2.second.isCorrect && !thisInfo2.second.isTrustable) isCorrect_o=-2; //erado nao tem match real

        }
        tout.Fill();
    }
    tout.Write();
}

void save_TMap(std::unordered_map<TrackKey,std::unordered_map<FlashKey, MatchInfo, FlashKeyHash>,TrackKeyHash> tmap,std::string treeName)
{  
    int run_o,event_o;
    int trackID_o,type_o, flashID_o;
    double score_o,deltato_o,flashTime_o,scoreMC_o;
    int isCorrect_o;
    double length;
    double x;
    double charge,light;

    TTree tout(treeName.c_str(),treeName.c_str());
    tout.Branch("run",&run_o,"run/I");
    tout.Branch("event",&event_o,"event/I");
    tout.Branch("trackID",&trackID_o,"trackID/I");
    tout.Branch("type",&type_o,"type/I");
    tout.Branch("flashID",&flashID_o,"flashID/I");
    tout.Branch("scoreData",&score_o,"scoreData/D");
    tout.Branch("deltaT0",&deltato_o,"deltaT0/D");
    tout.Branch("flashTime",&flashTime_o,"flashTime/D");
    tout.Branch("scoreMC",&scoreMC_o,"scoreMC/D");
    tout.Branch("isCorrect",&isCorrect_o,"isCorrect/I");
    tout.Branch("length",&length,"length/D");
    tout.Branch("x",&x,"x/D");
    tout.Branch("charge",&charge,"charge/D");
    tout.Branch("light",&light,"light/D");

    for(auto thisInfo:tmap)
    {
        run_o=thisInfo.first.run;
        event_o=thisInfo.first.event;
        trackID_o=thisInfo.first.track;
        type_o=thisInfo.first.type;
        for(auto thisInfo2:thisInfo.second)
        {
            flashID_o = thisInfo2.first.flash;
            score_o = thisInfo2.second.score;
            deltato_o = thisInfo2.second.deltaT0;
            flashTime_o = thisInfo2.second.flashTime;
            scoreMC_o = thisInfo2.second.scoreMC;
            length = thisInfo2.second.Length;
            x = thisInfo2.second.X;
            if(thisInfo2.second.isCorrect && thisInfo2.second.isTrustable) isCorrect_o=1;
            else if(!thisInfo2.second.isCorrect && thisInfo2.second.isTrustable) isCorrect_o=0;
            else if(thisInfo2.second.isCorrect && !thisInfo2.second.isTrustable) isCorrect_o=-1;
            else if(!thisInfo2.second.isCorrect && !thisInfo2.second.isTrustable) isCorrect_o=-2;
        }
        tout.Fill();
    }
    tout.Write();
}


void my_flashEff()
{
    //arquivo com dados de match com data
    std::string DataFileName = "data/track_25cm_ALL_20PE_50PETotal_N10_newfit5/cosmics_detsim_000010_pandora_reco_flash_reco_FlashMatchingDataAnalyser.root";
    std::string MCFileName = "data/track_25cm_ALL_20PE_50PETotal_N10_newfit3/cosmics_detsim_000010_pandora_reco_flash_reco_FlashMatchingMCAnalyser.root";

    std::string treeMC_name = "flashinfo/treeFT";
    std::string treeF_name = "waveAna/treeMatchF";
    std::string treeT_name = "waveAna/treeMatchT";

    bool only_maximum_MC_score=false;

    TFile* fMC = TFile::Open(MCFileName.c_str(), "READ");
    if (!fMC || fMC->IsZombie())
    {
        std::cerr << "Erro abrindo MC file: " << MCFileName << "\n";
        return;
    }
   
    std::unordered_map<FlashKey,std::unordered_map<TrackKey, MCInfo, TrackKeyHash>,FlashKeyHash> flashMapMC;
    std::unordered_map<TrackKey,std::unordered_map<FlashKey, MCInfo, FlashKeyHash>,TrackKeyHash> trackMapMC;
    
    fill_MCMap(flashMapMC,trackMapMC,fMC,treeMC_name,only_maximum_MC_score);

    //--------------------------------------------------------------------------

    TFile* fData = TFile::Open(DataFileName.c_str(), "READ");
    if (!fData || fData->IsZombie())
    {
        std::cerr << "Erro abrindo truth file: " << DataFileName << "\n";
        return;
    }

    std::unordered_map<FlashKey,std::unordered_map<TrackKey, MatchInfo, TrackKeyHash>,FlashKeyHash> flashMapData_Hungarian;
    std::unordered_map<FlashKey,std::unordered_map<TrackKey, MatchInfo, TrackKeyHash>,FlashKeyHash> flashMapData_Hungry;
    std::unordered_map<TrackKey,std::unordered_map<FlashKey, MatchInfo, FlashKeyHash>,TrackKeyHash> trackMapData_Hungarian;
    std::unordered_map<TrackKey,std::unordered_map<FlashKey, MatchInfo, FlashKeyHash>,TrackKeyHash> trackMapData_Hungry;

    //---------------------------------------------------------------------------

    fill_DataFMap(flashMapData_Hungarian,flashMapData_Hungry,fData,treeF_name,flashMapMC);
    std::cout << "############################################" <<std::endl;
    fill_DataTMap(trackMapData_Hungarian,trackMapData_Hungry,fData,treeT_name,trackMapMC);
   

// =================== WRITE OUTPUT ROOT (per map) ===================
    namespace fs = std::filesystem;
    fs::path outPath = fs::path(DataFileName).parent_path() /
                       (fs::path(DataFileName).stem().string() + "_track_match_flags.root");
    if (outPath.parent_path().empty()) outPath = fs::path("track_match_flags.root");

    int run_o,event_o;
    int trackID_o,type_o, flashID_o;
    double score_o,deltato_o,flashTime_o,scoreMC_o;

    TFile fout(outPath.string().c_str(), "RECREATE");
    fout.cd();
    
    save_FMap(flashMapData_Hungry,"treeFHungry");
    save_FMap(flashMapData_Hungarian,"treeFHungarian");
    save_TMap(trackMapData_Hungry,"treeTHungry");
    save_TMap(trackMapData_Hungarian,"treeTHungarian");

    fout.Write();
    fout.Close();

}