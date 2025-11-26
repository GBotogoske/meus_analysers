#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TTree.h"

#include <iostream> // Para debug opcional

#include "dunecore/DuneObj/OpDetDivRec.h" //novo
#include "larsim/MCCheater/ParticleInventoryService.h" //novo


class MyWaveformAnalyzerDetSimMC : public art::EDAnalyzer 
{
  public:
    explicit MyWaveformAnalyzerDetSimMC(fhicl::ParameterSet const& p);

    void analyze(art::Event const& e) override;
    void beginJob() override;

  private:
    art::InputTag fInputTag;

    TTree* fTree = nullptr;
  
    int fRun, fEvent;
    int fOfflineChannel;
    double fTimestamp;
    std::vector<raw::ADC_Count_t> fADCValue;

    std::vector<int> fTrackIDs;
    std::vector<int> fPDGIDs;
    std::vector<int> fMother;

    std::vector<double> fPosIniX, fPosIniY, fPosIniZ;
    std::vector<double> fPosEndX, fPosEndY, fPosEndZ;
    std::vector<double> fMomIniX, fMomIniY, fMomIniZ;
    std::vector<double> fMomEndX, fMomEndY, fMomEndZ;
    std::vector<double> fEneIni, fEneEnd;
 
};

MyWaveformAnalyzerDetSimMC::MyWaveformAnalyzerDetSimMC(fhicl::ParameterSet const& p)
  : EDAnalyzer(p),
    fInputTag(p.get<std::string>("input_tag"))
{
}

void MyWaveformAnalyzerDetSimMC::beginJob() 
{
        art::ServiceHandle<art::TFileService> tfs;
        fTree = tfs->make<TTree>("waveform_tree", "Waveforms");
        fTree->Branch("run" , &fRun);
        fTree->Branch("event" , &fEvent);
        fTree->Branch("offline_channel" , &fOfflineChannel);
        fTree->Branch("timestamp" , &fTimestamp);
        fTree->Branch("adc" , &fADCValue);
        fTree->Branch("track", &fTrackIDs);
        fTree->Branch("pdg", &fPDGIDs);
        fTree->Branch("mother", &fMother);

        // Posição inicial
        fTree->Branch("vx", &fPosIniX);
        fTree->Branch("vy", &fPosIniY);
        fTree->Branch("vz", &fPosIniZ);

        // Posição final
        fTree->Branch("ex", &fPosEndX);
        fTree->Branch("ey", &fPosEndY);
        fTree->Branch("ez", &fPosEndZ);

        // Momento inicial
        fTree->Branch("px", &fMomIniX);
        fTree->Branch("py", &fMomIniY);
        fTree->Branch("pz", &fMomIniZ);

        // Momento final
        fTree->Branch("epx", &fMomEndX);
        fTree->Branch("epy", &fMomEndY);
        fTree->Branch("epz", &fMomEndZ);

        // Energia
        fTree->Branch("e", &fEneIni);
        fTree->Branch("ee", &fEneEnd);

}

void MyWaveformAnalyzerDetSimMC::analyze(art::Event const& e)
{
        fRun = e.run();
        fEvent = e.event();

        auto wfHandle = e.getHandle<std::vector<raw::OpDetWaveform>>(fInputTag);

        auto divrecHandle = e.getHandle<std::vector<sim::OpDetDivRec>>(art::InputTag("opdigi")); // Ajuste a tag

        if (!divrecHandle || divrecHandle->empty()) 
        {
            std::cout << "No DivRec data in event: run " << fRun << ", event " << fEvent << std::endl;
            return;
        }

        std::unordered_map<int, const sim::OpDetDivRec*> divrec_map;
        for (const auto& divrec : *divrecHandle) 
        {
                divrec_map[divrec.OpDetNum()] = &divrec;
        }

        if (!wfHandle || wfHandle->empty())
        {
                std::cout << "No PDS data in event: run " << fRun << ", event " << fEvent << std::endl;
                return;
        }

        // Agora converte para vetor e busca PDGIDs
        art::ServiceHandle<cheat::ParticleInventoryService const> pis;                        

        for (const auto& wf : *wfHandle) 
        {
                fOfflineChannel = wf.ChannelNumber();
                fTimestamp = wf.TimeStamp();

                fADCValue.assign(wf.begin(), wf.end());
                if (fADCValue.empty()) continue;

                fTrackIDs.clear();
                fPDGIDs.clear();

                fMother.clear();

                fPosIniX.clear(); fPosIniY.clear(); fPosIniZ.clear();
                fPosEndX.clear(); fPosEndY.clear(); fPosEndZ.clear();
                fMomIniX.clear(); fMomIniY.clear(); fMomIniZ.clear();
                fMomEndX.clear(); fMomEndY.clear(); fMomEndZ.clear();
                fEneIni.clear(); fEneEnd.clear();


                if(fOfflineChannel >= 80)
                {
                    continue;  
                }
                
               std::set<int> uniqueTrackIDs; // para evitar duplicados temporários

                auto it = divrec_map.find(fOfflineChannel);
                if (it != divrec_map.end())
                {
                        const auto& divrec = *(it->second);
                        uniqueTrackIDs.clear();
                        
                         auto timeChans = divrec.GetTimeChans();

                        /* std::sort(timeChans.begin(), timeChans.end(), 
                                [](const sim::OpDet_Time_Chans& a, const sim::OpDet_Time_Chans& b) {
                                    return a.time < b.time;
                                }); */

                        for (const auto& timeChan : timeChans) 
                        {
                            double time_us = timeChan.time / 1000.0;

                            if (time_us < fTimestamp) continue;
                            if (time_us > fTimestamp + 16.0*1024/1000.0) continue;

                            for (const auto& phot : timeChan.phots) 
                            {
                                uniqueTrackIDs.insert(phot.trackID);
                            }
                        }

                        const auto& particleList = pis->ParticleList();

                        for (const int trackID : uniqueTrackIDs) 
                        {
                                fTrackIDs.push_back(trackID);

                                const simb::MCParticle* particle = particleList.Particle(trackID);
                                if (particle) 
                                {
                                        fPDGIDs.push_back(particle->PdgCode());
                                        fMother.push_back(particle->Mother());

                                        fPosIniX.push_back(particle->Vx());
                                        fPosIniY.push_back(particle->Vy());
                                        fPosIniZ.push_back(particle->Vz());

                                        fPosEndX.push_back(particle->EndX());
                                        fPosEndY.push_back(particle->EndY());
                                        fPosEndZ.push_back(particle->EndZ());

                                        fMomIniX.push_back(particle->Px());
                                        fMomIniY.push_back(particle->Py());
                                        fMomIniZ.push_back(particle->Pz());

                                        fMomEndX.push_back(particle->EndPx());
                                        fMomEndY.push_back(particle->EndPy());
                                        fMomEndZ.push_back(particle->EndPz());

                                        fEneIni.push_back(particle->E());
                                        fEneEnd.push_back(particle->EndE());
                                } 
                                else 
                                {
                                        fPDGIDs.push_back(-5e3);
                                        fMother.push_back(-5e3);
                                        fPosIniX.push_back(-9999.);
                                        fPosIniY.push_back(-9999.);
                                        fPosIniZ.push_back(-9999.);

                                        fPosEndX.push_back(-9999.);
                                        fPosEndY.push_back(-9999.);
                                        fPosEndZ.push_back(-9999.);

                                        fMomIniX.push_back(-9999.);
                                        fMomIniY.push_back(-9999.);
                                        fMomIniZ.push_back(-9999.);

                                        fMomEndX.push_back(-9999.);
                                        fMomEndY.push_back(-9999.);
                                        fMomEndZ.push_back(-9999.);

                                        fEneIni.push_back(-9999.);
                                        fEneEnd.push_back(-9999.);
                                }
                        }
                        
                } 
                fTree->Fill();
        }
        std::cout << "Processed event: run " << fRun << ", event " << fEvent << std::endl;
}

DEFINE_ART_MODULE(MyWaveformAnalyzerDetSimMC)
