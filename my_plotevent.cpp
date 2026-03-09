#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLegend.h"
#include "Rtypes.h" 

#include <filesystem>
#include <limits>

//1 TPC tick = 512 ns
//1 PDS tick = 16 ns

int correct_color = kGreen+3;
int wrong_color = kRed;
int nothing_color = kBlack;
int shift_color = kViolet;

std::map<int,  int> colorMap = {{1, correct_color},{0, wrong_color},{-2, wrong_color},{-1,correct_color}};

struct TrackKey 
{
    int run=0, event=0, track=0, type=0;
    bool operator==(TrackKey const& o) const 
    {
        return run==o.run && event==o.event && track==o.track && type==o.type;
    }
};

struct TrackKeyHash
{
    std::size_t operator()(TrackKey const& k) const noexcept 
    {
        std::size_t h = 1469598103934665603ull;
        auto mix = [&](std::size_t x){ h ^= x + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2); };
        mix((std::size_t)k.run);
        mix((std::size_t)k.event);
        mix((std::size_t)k.track);
        mix((std::size_t)k.type);
        return h;
    }
};

struct FlagInfo 
{
    int correct = 0;
    double FlashTime=0.0;
};

struct TrackHits 
{
    double length = 0.0;
    std::vector<float> x[8][3]; // wire
    std::vector<float> y[8][3]; // tick
};

namespace fs = std::filesystem;

void my_plotevents(int runSel = 1, int evtSel = 1)
{
    gROOT->SetBatch(kTRUE);//nao abrir janelas

    std::string fileName  = "data/trackshower_20cm_ALL_20PE_50PETotal_Norm_N10/1/cosmics_detsim_000001_pandora_reco_flash_reco_myAnalyserWire.root";
    std::string treeNameWire  = "waveAna/wire_tree";
    std::string treeTrackWire = "waveAna/track_tree";

    std::string dataMatchName = "data/track_25cm_ALL_20PE_50PETotal_N10_newfit5/cosmics_detsim_000001_pandora_reco_flash_reco_FlashMatchingDataAnalyser_track_match_flags.root";
    std::string treeMatchName  = "treeFHungry";


    fs::path inPath(fileName);
    // pasta de saída base
    fs::path outBase = inPath.has_parent_path() ? (inPath.parent_path() / "plots"): fs::path("plots");

    float thrPlane[3] = {50.0, 50.0, 18.0};
    float lengthCut = 0;

    TFile f(fileName.c_str(), "READ");
    if (f.IsZombie()) 
    {
        std::cerr << "ERROR: cannot open " << fileName << "\n";
        return;
    }

    auto* tw = dynamic_cast<TTree*>(f.Get(treeNameWire.c_str()));
    auto* tt = dynamic_cast<TTree*>(f.Get(treeTrackWire.c_str()));
    if (!tw) { std::cerr << "ERROR: no wire_tree\n"; return; }
    if (!tt) { std::cerr << "ERROR: no track_tree\n"; return; }

    std::cout << "Plotting run=" << runSel << " event=" << evtSel << "\n";
    std::cout << "Length Cut=" << lengthCut << "\n";

    std::unordered_map<TrackKey, FlagInfo, TrackKeyHash> flagsByTrack;

    {
        TFile fFlags(dataMatchName.c_str(), "READ");
        if (fFlags.IsZombie()) 
        {
            std::cerr << "WARN: nao consegui abrir flags file: " << dataMatchName << std::endl;
        }  
        else 
        {
            auto* tDump = dynamic_cast<TTree*>(fFlags.Get(treeMatchName.c_str()));
            if (!tDump) 
            {
                std::cerr << "WARN: nao achei TTree "<< treeMatchName.c_str() << " em " << dataMatchName << std::endl;
            } 
            else 
            {
                int run=0, event=0;
                int trackID=0, typeID =0;
                int isCorrect = 0;
                double flashTime;
                tDump->SetBranchAddress("run", &run);
                tDump->SetBranchAddress("event", &event);
                tDump->SetBranchAddress("trackID", &trackID);
                tDump->SetBranchAddress("type", &typeID);
                tDump->SetBranchAddress("flashTime", &flashTime);
                tDump->SetBranchAddress("isCorrect", &isCorrect);

                const Long64_t n = tDump->GetEntries();
                for (Long64_t i=0; i<n; ++i)
                {
                    tDump->GetEntry(i);
                    if (run != runSel || event != evtSel) continue;
                    TrackKey track_k{run, event, trackID, typeID};
                    flagsByTrack[track_k] = FlagInfo{isCorrect,flashTime};
                }
            }
        }
    }

    //LENDO OS DADOS DA TPC -- CHARGE
    //primeiro procurar os maximos indices

    int maxWire[8][3];
    int maxTick[8][3];
    for (int T=0; T<8; ++T)
        for (int P=0; P<3; ++P) 
            { maxWire[T][P] = -1; maxTick[T][P] = -1; }

    {
        TTreeReader r1(tw);
        TTreeReaderValue<int> run(r1, "run");
        TTreeReaderValue<int> evt(r1, "event");
        TTreeReaderValue<std::vector<int>> tpc(r1, "tpc");
        TTreeReaderValue<std::vector<int>> plane(r1, "plane");
        TTreeReaderValue<std::vector<int>> wire(r1, "wire");
        TTreeReaderValue<std::vector<short>> adc(r1, "adc");

        while (r1.Next()) 
        {
            if (*run != runSel || *evt != evtSel) continue;

            int nTicks = (int)adc->size();
            int nMap = std::min({ (int)tpc->size(), (int)plane->size(), (int)wire->size() });

            for (int i=0; i<nMap; ++i) 
            {
                int T = (*tpc)[i];
                int P = (*plane)[i];
                int W = (*wire)[i];
                if (T<0 || T>=8 || P<0 || P>=3) continue;
                maxWire[T][P] = std::max(maxWire[T][P], W);
                maxTick[T][P] = std::max(maxTick[T][P], nTicks);
            }
        }
    }

    //crias os histogramas e ler os dados
    TH2F* h[8][3] = {{nullptr}};
    for (int T=0; T<8; ++T) 
    {
        for (int P=0; P<3; ++P)
        {
            if (maxWire[T][P] < 0 || maxTick[T][P] <= 0) continue;
            int nX = maxWire[T][P] + 1;
            int nY = maxTick[T][P];
            std::string name  = "h_tpc" + std::to_string(T) + "_p" + std::to_string(P)+std::to_string(runSel)+std::to_string(evtSel);
            std::string title = "Run " + std::to_string(runSel) + " Event " + std::to_string(evtSel) +
                                " - TPC " + std::to_string(T) + " Plane " + std::to_string(P) +
                                ";wire;tick;ADC";
            h[T][P] = new TH2F(name.c_str(), title.c_str(), nX, 0, nX, nY, 0, nY);
        }
    }

    // ===================== preencher imagem =====================
    {
        TTreeReader r2(tw);
        TTreeReaderValue<int> run(r2, "run");
        TTreeReaderValue<int> evt(r2, "event");
        TTreeReaderValue<std::vector<int>> tpc(r2, "tpc");
        TTreeReaderValue<std::vector<int>> plane(r2, "plane");
        TTreeReaderValue<std::vector<int>> wire(r2, "wire");
        TTreeReaderValue<std::vector<short>> adc(r2, "adc");

        while (r2.Next()) 
        {
            if (*run != runSel || *evt != evtSel) continue;

            int nTicks = (int)adc->size();
            int nMap = std::min({ (int)tpc->size(), (int)plane->size(), (int)wire->size() });

            double sum = 0.0;
            for (short s : *adc) sum += (double)s;
            float mean = (float)(sum / nTicks);

            for (int i=0; i<nMap; ++i) 
            {
                int T = (*tpc)[i];
                int P = (*plane)[i];
                int W = (*wire)[i];
                if (T<0 || T>=8 || P<0 || P>=3) continue;
                if (!h[T][P]) continue;
                if (W < 0 || W > maxWire[T][P]) continue;

                float thr = thrPlane[P];

                for (int tick=0; tick<nTicks; ++tick) 
                {
                    float v = (float)(*adc)[tick] - mean;
                    float bit = (std::abs(v) > thr) ? 1.0f : 0.2f;

                    int xbin = W + 1;
                    int ybin = tick + 1;

                    float old = h[T][P]->GetBinContent(xbin, ybin);
                    if (bit > old) h[T][P]->SetBinContent(xbin, ybin, bit);
                }
            }
        }
    }

    // ===================== HITS POR TRACK (para colorir) =====================

     std::unordered_map<TrackKey, TrackHits, TrackKeyHash> hitsByTrack;// trackID -> hits arrays

    {
        TTreeReader rH(tt);
        TTreeReaderValue<int>    runH(rH, "run");
        TTreeReaderValue<int>    evtH(rH, "event");
        TTreeReaderValue<double> lengthH(rH, "length");
        TTreeReaderValue<int> trackID(rH, "trackID");
        TTreeReaderValue<int> tracktype(rH, "type");

        TTreeReaderValue<std::vector<float>> hitTime(rH, "hitTime");
        TTreeReaderValue<std::vector<int>>   hitTPC(rH, "hitTPC");
        TTreeReaderValue<std::vector<int>>   hitPlane(rH, "hitPlane");
        TTreeReaderValue<std::vector<int>>   hitWire(rH, "hitWire");
      
        while (rH.Next()) 
        {
            if (*runH != runSel || *evtH != evtSel) continue;
            if (*lengthH < lengthCut) continue;

            TrackKey k{*runH, *evtH, *trackID, *tracktype};
            auto& th = hitsByTrack[k];
            th.length = *lengthH;
            int n = std::min({ (int)hitTime->size(), (int)hitTPC->size(),(int)hitPlane->size(), (int)hitWire->size() });

            for (int i=0; i<n; ++i) 
            {
                int T = (*hitTPC)[i];
                int P = (*hitPlane)[i];
                int W = (*hitWire)[i];
                float ttick = (*hitTime)[i];

                if (T < 0 || T >= 8 || P < 0 || P >= 3) continue;

                th.x[T][P].push_back((float)W);
                th.y[T][P].push_back(ttick);
            }
        }
    }
    std::cout << "Tracks com hits no evento: " << hitsByTrack.size() << "\n";

    // ===================== PLOT =====================
    gStyle->SetOptStat(0);

    std::vector<std::string> name_plane = {"U","V","C"};
    const char* tag = treeMatchName.c_str();

    fs::path outDirPath = outBase / tag;
    std::string outDir = outDirPath.string();
    gSystem->mkdir(outDir.c_str(), true);

    for (int T=0; T<8; ++T)
    {
        bool hasAny = false;
        for (int P=0; P<3; ++P) if (h[T][P]) hasAny = true;
        if(T==0 || T==3 || T== 4 || T==7) continue;
        if (!hasAny) continue;

        TCanvas c(Form("c_tpc%d_%s", T, tag), Form("TPC %d (%s)", T, tag), 1400, 1200);
        c.SetBatch(kTRUE); //nao abrir janela
        c.Divide(1, 3, 0.001, 0.001);
        std::vector<TGraph*> graphsToDelete;
        std::vector<TLegend*> legendsToDelete;

        for (int P=0; P<3; ++P)
        {
            int k;
            if (P==0) k=3;
            if (P==1) k=2;
            if (P==2) k=1;

            c.cd(k);

            gPad->SetLeftMargin(0.08);
            gPad->SetRightMargin(0.12);
            gPad->SetBottomMargin(0.10);
            gPad->SetTopMargin(0.08);

            if (h[T][P]) 
            {
                h[T][P]->SetTitle(Form("Run %d Event %d - TPC %d Plane %s  (color=%s);wire;tick;ADC",
                                        runSel, evtSel, T, name_plane[P].c_str(), tag));
                h[T][P]->SetMinimum(0);
                h[T][P]->SetMaximum(1);
                h[T][P]->Draw("COLZ");
            } 
            else 
            {
                TH2F empty("empty","No data;wire;tick;ADC", 10,0,10, 10,0,10);
                empty.Draw();
            }

            for (auto const& kv : hitsByTrack)
            {
                    TrackKey k = kv.first;
                    TrackHits const& th = kv.second;

                    auto const& vx = th.x[T][P];
                    auto const& vy = th.y[T][P];
                    if (vx.empty()) continue;

                    auto itF = flagsByTrack.find(k);
                    const FlagInfo* fptr = (itF != flagsByTrack.end()) ? &itF->second : nullptr;
            }

            //#####################plot dos hits####################################
             // overlay hits por track, com cor
            for (auto const& kv : hitsByTrack)
            {
                TrackKey k = kv.first;
                TrackHits const& th = kv.second;

                auto const& vx = th.x[T][P];
                auto const& vy = th.y[T][P];
                if (vx.empty()) continue;

                // achar flag (se não existir -> preto)
                auto itF = flagsByTrack.find(k);
                const FlagInfo* fptr = (itF != flagsByTrack.end()) ? &itF->second : nullptr;

                int col;
                if(fptr==nullptr) col=nothing_color;
                else col=colorMap[fptr->correct];
                auto* gr = new TGraph((int)vx.size(), vx.data(), vy.data());

                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.35);
                gr->SetMarkerColor(col);
                gr->Draw("P SAME");

                if(col==correct_color)
                {
                    //find min time
                    auto t_vector = vy;
                    auto adc_vector = vx;

                    const size_t n = std::min(t_vector.size(), adc_vector.size());
                    t_vector.resize(n);
                    adc_vector.resize(n);

                    auto t_min= t_vector[0];
                    for (size_t i = 0; i < vy.size(); ++i)
                    {
                        const auto ti = t_vector[i];
                    }

                    double drift_speed = 0.16; // cm/us
                    double tick_tpc = 0.512; // us

                    //(hit*tick_tpc)-250=t0+delta_x/v
                    //delta_x = {[(hit*tick_tpc)-250] - t0}/tick_tpc=hit-250/tick-t0/tick

                    float deltat = 250+fptr->FlashTime;

                    for(auto& tj : t_vector)
                    {
                        tj -= deltat;
                    }

                    auto* gr_shift = new TGraph((int)adc_vector.size(), adc_vector.data(), t_vector.data());

                    gr_shift->SetMarkerStyle(20);
                    gr_shift->SetMarkerSize(0.35);
                    gr_shift->SetMarkerColor(shift_color);
                    gr_shift->Draw("P SAME");

                    graphsToDelete.push_back(gr_shift);
                }

                graphsToDelete.push_back(gr);
            }

            // legenda simples (uma vez por pad)
            auto* leg = new TLegend(0.14, 0.80, 0.45, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            // dummy graphs pra legenda
            double xx[1] = {0}, yy[1] = {0};
            auto* gOk  = new TGraph(1, xx, yy); gOk->SetMarkerStyle(20); gOk->SetMarkerColor(correct_color); gOk->SetMarkerSize(0.6);
            auto* gBad = new TGraph(1, xx, yy); gBad->SetMarkerStyle(20); gBad->SetMarkerColor(wrong_color);     gBad->SetMarkerSize(0.6);
            auto* gBlk = new TGraph(1, xx, yy); gBlk->SetMarkerStyle(20); gBlk->SetMarkerColor(nothing_color);   gBlk->SetMarkerSize(0.6);
            auto* gShift = new TGraph(1, xx, yy); gShift->SetMarkerStyle(20); gShift->SetMarkerColor(shift_color);   gShift->SetMarkerSize(0.6);


            leg->AddEntry(gOk,  "truthable + acerto", "P");
            leg->AddEntry(gBad, "truthable + erro",   "P");
            leg->AddEntry(gBlk, "outros / sem flag",  "P");
            leg->AddEntry(gShift, "shifted ",  "P");
            leg->Draw();

            graphsToDelete.push_back(gOk);
            graphsToDelete.push_back(gBad);
            graphsToDelete.push_back(gBlk);
            graphsToDelete.push_back(gShift);
            legendsToDelete.push_back(leg);

        }
        c.SaveAs(Form("%s/TPC_%s_run_%d_event_%d_tpc_%d.png", outDir.c_str(),tag,runSel,evtSel,T));
    }
  
    std::cout << "Saved plots (" << tag << ") in: " << outDir << "\n";

}

void my_plotevent()
{
    for(int i=1;i<=10;i++)
        my_plotevents(1,i);
}