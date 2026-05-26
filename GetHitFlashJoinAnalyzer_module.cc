#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "TTree.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

class GetHitFlashJoinAnalyzer : public art::EDAnalyzer {
public:
    explicit GetHitFlashJoinAnalyzer(fhicl::ParameterSet const& p);

    void beginJob() override;
    void analyze(art::Event const& e) override;

private:

    struct HitInfo {
        double time;
        double pe;
        int channel;
    };

    struct HitGroup {
        int id = -1;

        double start_time = 0.0;
        double end_time   = 0.0;

        double sum_pe  = 0.0;
        double sum_tpe = 0.0;
        double sum_t   = 0.0;

        std::vector<int> channels;
        std::vector<double> hit_times;
        std::vector<double> hit_pes;

        void add(HitInfo const& h)
        {
            if (channels.empty()) {
                start_time = h.time;
                end_time   = h.time;
            }
            else {
                if (h.time < start_time) start_time = h.time;
                if (h.time > end_time)   end_time   = h.time;
            }

            sum_pe  += h.pe;
            sum_tpe += h.time * h.pe;
            sum_t   += h.time;

            channels.push_back(h.channel);
            hit_times.push_back(h.time);
            hit_pes.push_back(h.pe);
        }

        int nHits() const
        {
            return static_cast<int>(channels.size());
        }

        double time() const
        {
            if (sum_pe > 0.0) {
                return sum_tpe / sum_pe; // tempo médio pesado por PE
            }

            if (!channels.empty()) {
                return sum_t / static_cast<double>(channels.size());
            }

            return 0.0;
        }

        double width() const
        {
            return end_time - start_time;
        }
    };

    struct FlashInfo {
        int id = -1;
        double time = 0.0;
        double pe = 0.0;
        double width = 0.0;
    };

private:

    art::InputTag fOpHitTag;
    art::InputTag fOpFlashTag;

    double fDeltaTHitGroup;
    double fDeltaTMatch;
    double fMinHitPE;

    bool fUseAbsTime;
    bool fSaveUnmatched;
    bool fRequireUniqueFlash;

    TTree* fTree = nullptr;

    int fRun;
    int fEvent;

    int fMatched;

    int fGroupID;
    int fFlashID;

    int fGroupNHits;

    double fHitGroupTime;
    double fHitGroupStartTime;
    double fHitGroupEndTime;
    double fHitGroupWidth;
    double fHitGroupPE;

    double fFlashTime;
    double fFlashWidth;
    double fFlashPE;

    double fDeltaT;

    std::vector<int> fGroupChannels;
    std::vector<double> fGroupHitTimes;
    std::vector<double> fGroupHitPEs;
};

GetHitFlashJoinAnalyzer::GetHitFlashJoinAnalyzer(fhicl::ParameterSet const& p)
    : EDAnalyzer(p)
{
    fOpHitTag  = p.get<art::InputTag>("OpHitTag",  art::InputTag("ophitspe", "", "myFlash"));
    fOpFlashTag = p.get<art::InputTag>("OpFlashTag", art::InputTag("opflash",  "", "myFlash"));

    // Delta de tempo para juntar OpHits no mesmo grupo.
    // A unidade é a mesma de OpHit::PeakTime().
    fDeltaTHitGroup = p.get<double>("DeltaTHitGroup", 0.50);

    // Delta máximo para associar um grupo de hits a um OpFlash.
    fDeltaTMatch = p.get<double>("DeltaTMatch", 0.50);

    // Corte mínimo em PE por hit.
    fMinHitPE = p.get<double>("MinHitPE", 0.0);

    // false: usa OpHit::PeakTime() e OpFlash::Time()
    // true : usa OpHit::PeakTimeAbs() e OpFlash::AbsTime()
    fUseAbsTime = p.get<bool>("UseAbsTime", false);

    // Se true, salva também grupos de hits sem flash associado.
    fSaveUnmatched = p.get<bool>("SaveUnmatched", true);

    // Se true, cada flash só pode ser usado uma vez.
    fRequireUniqueFlash = p.get<bool>("RequireUniqueFlash", false);
}

void GetHitFlashJoinAnalyzer::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fTree = tfs->make<TTree>("hit_flash_tree", "OpHit temporal groups matched to OpFlash");

    fTree->Branch("run",   &fRun,   "run/I");
    fTree->Branch("event", &fEvent, "event/I");

    fTree->Branch("matched", &fMatched, "matched/I");

    fTree->Branch("group_id", &fGroupID, "group_id/I");
    fTree->Branch("flash_id", &fFlashID, "flash_id/I");

    fTree->Branch("group_nhits", &fGroupNHits, "group_nhits/I");

    fTree->Branch("hit_group_time",       &fHitGroupTime,      "hit_group_time/D");
    fTree->Branch("hit_group_start_time", &fHitGroupStartTime, "hit_group_start_time/D");
    fTree->Branch("hit_group_end_time",   &fHitGroupEndTime,   "hit_group_end_time/D");
    fTree->Branch("hit_group_width",      &fHitGroupWidth,     "hit_group_width/D");
    fTree->Branch("hit_group_pe",         &fHitGroupPE,        "hit_group_pe/D");

    fTree->Branch("flash_time",  &fFlashTime,  "flash_time/D");
    fTree->Branch("flash_width", &fFlashWidth, "flash_width/D");
    fTree->Branch("flash_pe",    &fFlashPE,    "flash_pe/D");

    fTree->Branch("delta_t", &fDeltaT, "delta_t/D");

    fTree->Branch("group_channels",  &fGroupChannels);
    fTree->Branch("group_hit_times", &fGroupHitTimes);
    fTree->Branch("group_hit_pes",   &fGroupHitPEs);
}

void GetHitFlashJoinAnalyzer::analyze(art::Event const& e)
{
    fRun   = e.run();
    fEvent = e.event();

    auto hitHandle = e.getHandle<std::vector<recob::OpHit>>(fOpHitTag);
    if (!hitHandle || hitHandle->empty()) {
        std::cout << "No OpHit data in run " << fRun
                  << ", event " << fEvent << std::endl;
        return;
    }

    auto flashHandle = e.getHandle<std::vector<recob::OpFlash>>(fOpFlashTag);

    // -------------------------------------------------------
    // 1. Lê e ordena os OpHits por tempo
    // -------------------------------------------------------

    std::vector<HitInfo> hits;
    hits.reserve(hitHandle->size());

    for (auto const& hit : *hitHandle) {
        double pe = hit.PE();

        if (pe < fMinHitPE) continue;

        double time = fUseAbsTime ? hit.PeakTimeAbs() : hit.PeakTime();

        HitInfo h;
        h.time    = time;
        h.pe      = pe;
        h.channel = static_cast<int>(hit.OpChannel());

        hits.push_back(h);
    }

    if (hits.empty()) return;

    std::sort(
        hits.begin(),
        hits.end(),
        [](HitInfo const& a, HitInfo const& b) {
            return a.time < b.time;
        }
    );

    // -------------------------------------------------------
    // 2. Agrupa OpHits próximos em tempo
    // -------------------------------------------------------

    std::vector<HitGroup> groups;

    HitGroup current;
    current.id = 0;
    current.add(hits.front());

    for (size_t i = 1; i < hits.size(); ++i) {
        double dt = hits[i].time - current.end_time;

        if (dt <= fDeltaTHitGroup) {
            current.add(hits[i]);
        }
        else {
            groups.push_back(current);

            current = HitGroup();
            current.id = static_cast<int>(groups.size());
            current.add(hits[i]);
        }
    }

    groups.push_back(current);

    // -------------------------------------------------------
    // 3. Lê e ordena os OpFlash por tempo
    // -------------------------------------------------------

    std::vector<FlashInfo> flashes;

    if (flashHandle && !flashHandle->empty()) {
        flashes.reserve(flashHandle->size());

        for (size_t i = 0; i < flashHandle->size(); ++i) {
            auto const& flash = flashHandle->at(i);

            FlashInfo f;
            f.id    = static_cast<int>(i);
            f.time  = fUseAbsTime ? flash.AbsTime() : flash.Time();
            f.pe    = flash.TotalPE();
            f.width = flash.TimeWidth();

            flashes.push_back(f);
        }

        std::sort(
            flashes.begin(),
            flashes.end(),
            [](FlashInfo const& a, FlashInfo const& b) {
                return a.time < b.time;
            }
        );
    }

    std::vector<bool> flashUsed(flashes.size(), false);

    // -------------------------------------------------------
    // 4. Para cada grupo de hits, acha o flash mais próximo
    // -------------------------------------------------------

    for (auto const& group : groups) {

        double groupTime = group.time();

        int bestFlashIndex = -1;
        double bestAbsDT = std::numeric_limits<double>::max();

        for (size_t i = 0; i < flashes.size(); ++i) {
            if (fRequireUniqueFlash && flashUsed[i]) continue;

            double dt = groupTime - flashes[i].time;
            double absdt = std::abs(dt);

            if (absdt < bestAbsDT) {
                bestAbsDT = absdt;
                bestFlashIndex = static_cast<int>(i);
            }
        }

        bool matched = false;

        if (bestFlashIndex >= 0 && bestAbsDT <= fDeltaTMatch) {
            matched = true;
        }

        if (!matched && !fSaveUnmatched) {
            continue;
        }

        if (matched && fRequireUniqueFlash) {
            flashUsed[bestFlashIndex] = true;
        }

        // ---------------------------------------------------
        // 5. Preenche a árvore
        // ---------------------------------------------------

        fMatched = matched ? 1 : 0;

        fGroupID = group.id;
        fGroupNHits = group.nHits();

        fHitGroupTime      = groupTime;
        fHitGroupStartTime = group.start_time;
        fHitGroupEndTime   = group.end_time;
        fHitGroupWidth     = group.width();
        fHitGroupPE        = group.sum_pe;

        fGroupChannels = group.channels;
        fGroupHitTimes = group.hit_times;
        fGroupHitPEs   = group.hit_pes;

        if (matched) {
            auto const& flash = flashes[bestFlashIndex];

            fFlashID    = flash.id;
            fFlashTime  = flash.time;
            fFlashWidth = flash.width;
            fFlashPE    = flash.pe;

            fDeltaT = fHitGroupTime - fFlashTime;
        }
        else {
            fFlashID    = -1;
            fFlashTime  = std::numeric_limits<double>::quiet_NaN();
            fFlashWidth = std::numeric_limits<double>::quiet_NaN();
            fFlashPE    = std::numeric_limits<double>::quiet_NaN();

            fDeltaT = std::numeric_limits<double>::quiet_NaN();
        }

        fTree->Fill();
    }
}

DEFINE_ART_MODULE(GetHitFlashJoinAnalyzer)