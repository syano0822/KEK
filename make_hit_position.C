/**
 * make_hit_position.C
 *
 * Reads per-event TGraph waveforms stored by make_event_waveforms.C and
 * computes the amplitude-weighted hit position for each strip detector group.
 *
 * Input file structure (from make_event_waveforms.C):
 *   /event_000000/
 *     DUT sensor/   strip1..strip7   (7 strips, voltage in V pedestal-subtracted)
 *     MCP-PMT/      wf               (single channel — skipped for position)
 *     Tracking_1/   strip1..strip8
 *     Tracking_2/   strip1..strip8
 *     Tracking_4/   strip1..strip8
 *     Tracking_5/   strip1..strip8
 *   /event_000001/ ...
 *
 * Hit-position estimators:
 *   Linear  :  pos = Σ(w_i  · x_i) / Σ(w_i)     where w_i = amplitude_i [mV]
 *   Squared :  pos = Σ(w_i² · x_i) / Σ(w_i²)
 *
 * Strip positions:   x_i = x_first + (i−1) · pitch   [i = 1-based strip index]
 *   Default: strip1 = 14.5 mm, pitch = 0.5 mm
 *     → DUT strip7 = 17.5 mm,  Tracking strip8 = 18.0 mm
 *
 * Four position estimators are studied in parallel:
 *   top2 : 2 strips with highest amplitude  (from those > pos_thr_mV)
 *   top3 : 3 strips with highest amplitude  (from those > pos_thr_mV)
 *   thr  : all strips with amplitude > pos_thr_mV
 *   gaus : Gaussian fit to amplitude vs. strip position (all strips > pos_thr_mV);
 *          hit position = Gaussian mean;  requires >= 3 strips above threshold;
 *          fit range = [leftmost − pitch/2, rightmost + pitch/2];
 *          sigma constrained to [0.1, 5] × pitch;  failed fits are discarded.
 *
 * 1D output histograms per detector-group directory:
 *   h_pos_top2_lin,  h_pos_top2_sq
 *   h_pos_top3_lin,  h_pos_top3_sq
 *   h_pos_thr_lin,   h_pos_thr_sq
 *   h_pos_gaus
 *
 * 2D correlation histograms (directory "Correlations/"), 10 pairs
 * (strip detectors only — MCP-PMT excluded):
 *   DUT sensor × Tracking_1,2,4,5
 *   Tracking_1 × Tracking_2,4,5
 *   Tracking_2 × Tracking_4,5
 *   Tracking_4 × Tracking_5
 *   Each pair: one TH2D per method tag → 10 × 6 = 60 histograms total.
 *   Filled only when both groups have a valid hit in the same event.
 *
 * 1D position-difference histograms (directory "Differences/"):
 *   Tracking_1 − Tracking_5,  Tracking_1 − DUT sensor,
 *   DUT sensor − Tracking_5,  Tracking_2 − Tracking_4
 *   Each pair: one TH1D per method tag → 4 × 6 = 24 histograms total.
 *
 * 1D strip-multiplicity histograms (per detector-group directory):
 *   h_nstrips_<group>  : number of strips with amplitude > pos_thr_mV
 *   Range: 0..8 (9 bins, −0.5 to 8.5); filled for every event including zero.
 *
 * Channel name → strip index:  "strip3" → 3,  "wf" → 1
 *
 * ── Usage ───────────────────────────────────────────────────────────────────
 *  Interpreted:
 *    root 'make_hit_position.C("waveforms.root","pos.root")'
 *    root 'make_hit_position.C("waveforms.root","pos.root",8.,14.5,0.5,200,13.,19.)'
 *
 *  ACLiC compiled (recommended for large files):
 *    root 'make_hit_position.C+("waveforms.root","pos.root")'
 *
 * ── Strip-column modes ───────────────────────────────────────────────────────
 *  2-cols mode (use_4cols = false, default):
 *    All strips treated as one column; strip index used directly for position.
 *    One position histogram per group.
 *
 *  4-cols mode (use_4cols = true):
 *    DUT sensor channels are split into two physical columns:
 *      Column-1: Ch3→sub-strip1, Ch4→sub-strip2, Ch5→sub-strip3, Ch6→sub-strip4
 *      Column-2: Ch1→sub-strip1, Ch2→sub-strip2, Ch7→sub-strip3 (no Ch8 on DUT)
 *    Position within each column: x_first_mm + (sub_strip − 1) × x_pitch_mm.
 *    Two independent position histograms for DUT (_col1 / _col2).
 *    Tracking planes are always treated as 2-cols (unchanged).
 *    Correlations: DUT_col1 × Tracking, DUT_col2 × Tracking, Tracking × Tracking
 *                  (14 pairs total per tag).
 *    Differences: Trk1−Trk5, Trk1−DUT_col1, Trk1−DUT_col2,
 *                 DUT_col1−Trk5, DUT_col2−Trk5, Trk2−Trk4 (6 pairs per tag).
 *
 * ── Parameters ──────────────────────────────────────────────────────────────
 *  input_path  : waveform ROOT file (output of make_event_waveforms.C)
 *  output_path : output ROOT file
 *  pos_thr_mV  : min amplitude to include a strip [mV]    (default   8.0)
 *  x_first_mm  : position of sub-strip 1 [mm]             (default  14.5)
 *  x_pitch_mm  : strip pitch [mm]                         (default   0.5)
 *                In 4-cols mode, DUT pitch is automatically set to 1.0 mm
 *                and Tracking pitch remains 0.5 mm.
 *  nbins       : histogram bins (1D and each 2D axis)      (default 600)
 *  pos_min_mm  : histogram x-axis minimum [mm]            (default  13.0)
 *  pos_max_mm  : histogram x-axis maximum [mm]            (default  19.0)
 *  amp_min_mV  : amplitude axis minimum for pos-vs-amp [mV] (default   0.0)
 *  amp_max_mV  : amplitude axis maximum for pos-vs-amp [mV] (default 300.0)
 *  skip_edge   : if true, skip events where the leading strip is at the edge
 *                (strip1 or strip7/8 in 2-cols; sub-strip1 or sub-strip4 in 4-cols)
 *                (default false)
 *  diff_min_mm : position-difference histogram minimum [mm](default  -6.0)
 *  diff_max_mm : position-difference histogram maximum [mm](default  +6.0)
 *  use_4cols   : false = 2-cols mode (default); true = 4-cols mode
 * ────────────────────────────────────────────────────────────────────────────
 */

#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TF1.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

// ── Method tags (used as histogram suffixes and map keys) ─────────────────────
static const char* kTags[]  = {
    "top2_lin", "top2_sq",
    "top3_lin", "top3_sq",
    "thr_lin",  "thr_sq",
    "gaus"
};
static const int kNTags = 7;

// ── Extract 1-based strip index from a TGraph name ───────────────────────────
// "strip3" → 3,  "wf" → 1 (no trailing digits)
static int ChIndex(const std::string& name)
{
    int i = (int)name.size() - 1;
    while (i >= 0 && std::isdigit((unsigned char)name[i])) --i;
    ++i;
    if (i < (int)name.size()) return std::stoi(name.substr(i));
    return 1;
}

// ── 4-cols channel mapping ────────────────────────────────────────────────────
// Maps a 1-based strip index to {column (1|2), sub-strip index (1..4)}.
// Col-1: strip3,4,5,6 → sub 1,2,3,4
// Col-2: strip1,2,7,8 → sub 1,2,3,4
// Returns {0,0} for unrecognised indices (e.g. strip8 on a 7-strip DUT).
static std::pair<int,int> ColMap(int s)
{
    switch (s) {
        case 3: return {1, 1};
        case 4: return {1, 2};
        case 5: return {1, 3};
        case 6: return {1, 4};
        case 1: return {2, 1};
        case 2: return {2, 2};
        case 7: return {2, 3};
        case 8: return {2, 4};
        default: return {0, 0};
    }
}

// ── Amplitude-weighted centroid ───────────────────────────────────────────────
// pw    : vector of (position [mm], amplitude [mV])
// scheme: 0 = linear weight w,  1 = squared weight w²
// Returns −999 if total weight is zero.
static double Centroid(const std::vector<std::pair<double,double>>& pw, int scheme)
{
    double sumW = 0.0, sumWX = 0.0;
    for (const auto& p : pw) {
        double w = (scheme == 0) ? p.second : p.second * p.second;
        sumW  += w;
        sumWX += w * p.first;
    }
    return (sumW > 0.0) ? sumWX / sumW : -999.0;
}

void make_hit_position(
    const char* input_path,
    const char* output_path = "hit_position.root",
    Double_t    pos_thr_mV  =   8.0,
    Double_t    x_first_mm  =   14.5,
    Double_t    x_pitch_mm  =    0.5,
    Int_t       nbins       =   600,
    Double_t    pos_min_mm  =   13.0,
    Double_t    pos_max_mm  =   19.0,
    Double_t    amp_min_mV  =    0.0,
    Double_t    amp_max_mV  =  300.0,
    Bool_t      skip_edge   = false,   // if true, skip events where ch1 or ch8 is the leading strip
    Double_t    diff_min_mm =  -6.0,  // position-difference histogram minimum [mm]
    Double_t    diff_max_mm =   6.0,  // position-difference histogram maximum [mm]
    Bool_t      use_4cols   = false)  // false = 2-cols mode;  true = 4-cols mode
{
    // ── Open input file ───────────────────────────────────────────────────────
    TFile* fin = TFile::Open(input_path, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open input: " << input_path << "\n";
        return;
    }
    printf("[INFO] Input  file      : %s\n",  input_path);
    printf("[INFO] Amplitude thr    : > %.1f mV\n", pos_thr_mV);
    // Pitch per group: in 4-cols mode DUT sub-strip pitch is 2× the base pitch.
    const double dut_pitch_mm = use_4cols ? 2.0 * x_pitch_mm : x_pitch_mm;
    const double trk_pitch_mm = x_pitch_mm;

    printf("[INFO] Strip positions  : ch1 = %.2f mm,  DUT pitch = %.2f mm,  Trk pitch = %.2f mm\n",
           x_first_mm, dut_pitch_mm, trk_pitch_mm);
    printf("[INFO] Position range   : [%.2f, %.2f] mm  (%d bins)\n",
           pos_min_mm, pos_max_mm, (int)nbins);

    // ── Collect and sort event directories ───────────────────────────────────
    std::vector<std::string> evdirs;
    {
        TIter it(fin->GetListOfKeys());
        TKey* k;
        while ((k = (TKey*)it())) {
            std::string nm = k->GetName();
            if (nm.size() >= 6 && nm.compare(0, 6, "event_") == 0)
                evdirs.push_back(nm);
        }
    }
    std::sort(evdirs.begin(), evdirs.end());

    if (evdirs.empty()) {
        std::cerr << "[ERROR] No event_ directories found in " << input_path << "\n";
        fin->Close();
        return;
    }
    printf("[INFO] Events found     : %zu\n\n", evdirs.size());

    // ── Gaussian TF1 (created once, reused every event/group) ────────────────
    // Parameters: [0]=amplitude, [1]=mean (position), [2]=sigma
    TF1* fg_gaus = new TF1("_fg_gaus_pos", "gaus", pos_min_mm, pos_max_mm);

    // ── 1D histogram registry ─────────────────────────────────────────────────
    // key = "group/tag"
    std::map<std::string, TH1D*> h1;

    auto getH1 = [&](const std::string& grp,
                     const std::string& tag) -> TH1D* {
        std::string key = grp + "/" + tag;
        auto it = h1.find(key);
        if (it != h1.end()) return it->second;

        std::string hname = "h_pos_" + grp + "_" + tag;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        std::string title = grp + " hit pos [" + tag + "];Position [mm];Entries";
        TH1D* h = new TH1D(hname.c_str(), title.c_str(),
                           nbins, pos_min_mm, pos_max_mm);
        h->SetDirectory(nullptr);
        h1[key] = h;
        return h;
    };

    // ── 2D histogram registry ─────────────────────────────────────────────────
    // Strip detectors only (MCP-PMT has no strip position → excluded).
    // 5 strip groups → C(5,2) = 10 pairs × 6 tags = 60 TH2D total.
    // key = "corrLabel/tag"
    std::map<std::string, TH2D*> h2;

    struct CorrDef { const char* label; const char* xgrp; const char* ygrp; };
    const CorrDef kCorr[] = {
        // DUT sensor vs Tracking planes
        { "DUT_vs_Trk1",   "DUT sensor",  "Tracking_1" },
        { "DUT_vs_Trk2",   "DUT sensor",  "Tracking_2" },
        { "DUT_vs_Trk4",   "DUT sensor",  "Tracking_4" },
        { "DUT_vs_Trk5",   "DUT sensor",  "Tracking_5" },
        // Tracking_1 vs Tracking planes
        { "Trk1_vs_Trk2",  "Tracking_1",  "Tracking_2" },
        { "Trk1_vs_Trk4",  "Tracking_1",  "Tracking_4" },
        { "Trk1_vs_Trk5",  "Tracking_1",  "Tracking_5" },
        // Tracking_2 vs Tracking_4,5
        { "Trk2_vs_Trk4",  "Tracking_2",  "Tracking_4" },
        { "Trk2_vs_Trk5",  "Tracking_2",  "Tracking_5" },
        // Tracking_4 vs Tracking_5
        { "Trk4_vs_Trk5",  "Tracking_4",  "Tracking_5" },
    };
    const int kNCorr = 10;

    auto getH2 = [&](const std::string& label,
                     const std::string& tag,
                     const char* xgrp,
                     const char* ygrp) -> TH2D* {
        std::string key = label + "/" + tag;
        auto it = h2.find(key);
        if (it != h2.end()) return it->second;

        std::string hname = "h2_" + label + "_" + tag;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        std::string title = std::string(xgrp) + " vs " + ygrp +
                            " [" + tag + "];" +
                            xgrp + " pos [mm];" +
                            ygrp + " pos [mm]";
        TH2D* h = new TH2D(hname.c_str(), title.c_str(),
                           nbins, pos_min_mm, pos_max_mm,
                           nbins, pos_min_mm, pos_max_mm);
        h->SetDirectory(nullptr);
        h2[key] = h;
        return h;
    };

    // ── Position-difference histogram registry ────────────────────────────────
    // key = "diffLabel/tag"
    std::map<std::string, TH1D*> h1diff;

    struct DiffDef { const char* label; const char* agrp; const char* bgrp; };
    const DiffDef kDiff[] = {
        { "Trk1_minus_Trk5", "Tracking_1", "Tracking_5" },
        { "Trk1_minus_DUT",  "Tracking_1", "DUT sensor" },
        { "DUT_minus_Trk5",  "DUT sensor", "Tracking_5" },
        { "Trk2_minus_Trk4", "Tracking_2", "Tracking_4" },
    };
    const int kNDiff = 4;

    // ── 4-cols correlation pairs (14 = DUT_col×4 Tracking×2 + C(4,2) Tracking) ─
    // Only DUT is split into columns; Tracking planes stay as single groups.
    const CorrDef kCorr4[] = {
        // DUT col-1 vs each Tracking plane
        { "DUT_col1_vs_Trk1",  "DUT sensor_col1", "Tracking_1" },
        { "DUT_col1_vs_Trk2",  "DUT sensor_col1", "Tracking_2" },
        { "DUT_col1_vs_Trk4",  "DUT sensor_col1", "Tracking_4" },
        { "DUT_col1_vs_Trk5",  "DUT sensor_col1", "Tracking_5" },
        // DUT col-2 vs each Tracking plane
        { "DUT_col2_vs_Trk1",  "DUT sensor_col2", "Tracking_1" },
        { "DUT_col2_vs_Trk2",  "DUT sensor_col2", "Tracking_2" },
        { "DUT_col2_vs_Trk4",  "DUT sensor_col2", "Tracking_4" },
        { "DUT_col2_vs_Trk5",  "DUT sensor_col2", "Tracking_5" },
        // Tracking plane pairs (identical to 2-cols)
        { "Trk1_vs_Trk2",      "Tracking_1",      "Tracking_2" },
        { "Trk1_vs_Trk4",      "Tracking_1",      "Tracking_4" },
        { "Trk1_vs_Trk5",      "Tracking_1",      "Tracking_5" },
        { "Trk2_vs_Trk4",      "Tracking_2",      "Tracking_4" },
        { "Trk2_vs_Trk5",      "Tracking_2",      "Tracking_5" },
        { "Trk4_vs_Trk5",      "Tracking_4",      "Tracking_5" },
    };
    const int kNCorr4 = 14;

    // ── 4-cols difference pairs (6) ───────────────────────────────────────────
    const DiffDef kDiff4[] = {
        { "Trk1_minus_Trk5",     "Tracking_1",      "Tracking_5"      },
        { "Trk1_minus_DUT_col1", "Tracking_1",      "DUT sensor_col1" },
        { "Trk1_minus_DUT_col2", "Tracking_1",      "DUT sensor_col2" },
        { "DUT_col1_minus_Trk5", "DUT sensor_col1", "Tracking_5"      },
        { "DUT_col2_minus_Trk5", "DUT sensor_col2", "Tracking_5"      },
        { "Trk2_minus_Trk4",     "Tracking_2",      "Tracking_4"      },
    };
    const int kNDiff4 = 6;

    // Active pair arrays — selected once based on mode
    const CorrDef* corrArr = use_4cols ? kCorr4 : kCorr;
    const DiffDef* diffArr = use_4cols ? kDiff4 : kDiff;
    const int      nCorr   = use_4cols ? kNCorr4 : kNCorr;
    const int      nDiff   = use_4cols ? kNDiff4 : kNDiff;

    auto getH1Diff = [&](const std::string& label,
                         const std::string& tag,
                         const char* agrp,
                         const char* bgrp) -> TH1D* {
        std::string key = label + "/" + tag;
        auto it = h1diff.find(key);
        if (it != h1diff.end()) return it->second;

        std::string hname = "hdiff_" + label + "_" + tag;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        std::string title = std::string(agrp) + " - " + bgrp +
                            " [" + tag + "];Position difference [mm];Entries";
        TH1D* h = new TH1D(hname.c_str(), title.c_str(),
                           nbins, diff_min_mm, diff_max_mm);
        h->SetDirectory(nullptr);
        h1diff[key] = h;
        return h;
    };

    // ── Strip-multiplicity histogram registry ─────────────────────────────────
    // key = group name
    // DUT sensor: 0..7 strips (8 bins, −0.5 to 7.5)
    // Tracking planes: 0..8 strips (9 bins, −0.5 to 8.5)
    std::map<std::string, TH1D*> h_nstrips;

    auto getNStripsH = [&](const std::string& grp) -> TH1D* {
        auto it = h_nstrips.find(grp);
        if (it != h_nstrips.end()) return it->second;

        std::string hname = "h_nstrips_" + grp;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        // _col1 / _col2 sub-groups always have 4 sub-strips
        bool is_col = (grp.size() > 5 &&
                       (grp.substr(grp.size() - 5) == "_col1" ||
                        grp.substr(grp.size() - 5) == "_col2"));
        int nmax = is_col ? 4 : (grp == "DUT sensor") ? 7 : 8;
        std::string title = grp + " strips above threshold;"
                            "N strips above threshold;Entries";
        TH1D* h = new TH1D(hname.c_str(), title.c_str(), nmax + 1, -0.5, nmax + 0.5);
        h->SetDirectory(nullptr);
        h_nstrips[grp] = h;
        return h;
    };

    // ── Position vs leading-strip amplitude registries ────────────────────────
    // h2pa : TH2D   X = hit position [mm], Y = max-amplitude [mV]
    // hppa : TProfile X = hit position [mm], Y-mean = max-amplitude [mV]
    // key = "group/tag"  (same structure as h1)
    std::map<std::string, TH2D*>     h2pa;
    std::map<std::string, TProfile*> hppa;

    auto getH2PA = [&](const std::string& grp,
                       const std::string& tag) -> TH2D* {
        std::string key = grp + "/" + tag;
        auto it = h2pa.find(key);
        if (it != h2pa.end()) return it->second;

        std::string hname = "h2pa_" + grp + "_" + tag;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        std::string title = grp + " pos vs max-amp [" + tag + "];"
                            "Position [mm];Max-strip amplitude [mV]";
        TH2D* h = new TH2D(hname.c_str(), title.c_str(),
                           nbins, pos_min_mm, pos_max_mm,
                           nbins, amp_min_mV, amp_max_mV);
        h->SetDirectory(nullptr);
        h2pa[key] = h;
        return h;
    };

    auto getProfPA = [&](const std::string& grp,
                         const std::string& tag) -> TProfile* {
        std::string key = grp + "/" + tag;
        auto it = hppa.find(key);
        if (it != hppa.end()) return it->second;

        std::string hname = "hppa_" + grp + "_" + tag;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        std::string title = grp + " mean max-amp vs pos [" + tag + "];"
                            "Position [mm];Mean max-strip amplitude [mV]";
        TProfile* h = new TProfile(hname.c_str(), title.c_str(),
                                   nbins, pos_min_mm, pos_max_mm,
                                   amp_min_mV, amp_max_mV);
        h->SetDirectory(nullptr);
        hppa[key] = h;
        return h;
    };

    // ── Event loop ────────────────────────────────────────────────────────────
    for (size_t iev = 0; iev < evdirs.size(); ++iev) {

        TDirectory* evdir = nullptr;
        fin->GetObject(evdirs[iev].c_str(), evdir);
        if (!evdir) {
            fprintf(stderr, "[WARN] %s not found – skipping\n",
                    evdirs[iev].c_str());
            continue;
        }

        // Per-event position results: event_pos[group][tag] = position [mm]
        // Used after the group loop to fill 2D correlation histograms.
        std::map<std::string, std::map<std::string, double>> event_pos;

        // ── Group loop ────────────────────────────────────────────────────────
        TIter itg(evdir->GetListOfKeys());
        TKey* gk;
        while ((gk = (TKey*)itg())) {
            if (std::string(gk->GetClassName()) != "TDirectoryFile") continue;

            std::string grpname = gk->GetName();

            // MCP-PMT has only one channel ("wf") — no strip position to compute.
            if (grpname == "MCP-PMT") continue;

            TDirectory* grpdir  = evdir->GetDirectory(grpname.c_str());
            if (!grpdir) continue;

            // Per-group pitch: DUT uses dut_pitch_mm, Tracking uses trk_pitch_mm.
            const double pitch_mm = (grpname == "DUT sensor") ? dut_pitch_mm : trk_pitch_mm;

            // Collect (position [mm], amplitude [mV]) for channels above threshold.
            std::vector<std::pair<double,double>> pw_all;

            TIter itc(grpdir->GetListOfKeys());
            TKey* ck;
            while ((ck = (TKey*)itc())) {
                TObject* obj = ck->ReadObj();
                TGraph* gr   = dynamic_cast<TGraph*>(obj);
                if (!gr) { delete obj; continue; }

                double amp_mV = 0.0;
                int n = gr->GetN();
                if (n > 0) {
                    const double* Y = gr->GetY();
                    amp_mV = -(*std::min_element(Y, Y + n)) * 1e3;
                }
                delete gr;

                if (amp_mV <= pos_thr_mV) continue;

                int    idx = ChIndex(std::string(ck->GetName()));
                double x   = x_first_mm + (idx - 1) * pitch_mm;
                pw_all.push_back({x, amp_mV});
            }

            // ── 2-cols mode: single position per group ────────────────────────
            // Always used for Tracking planes; used for DUT when !use_4cols.
            if (!use_4cols || grpname != "DUT sensor") {
                // Fill multiplicity histogram for every event (including zero).
                getNStripsH(grpname)->Fill((int)pw_all.size());

                if (pw_all.empty()) continue;

                // Sort by amplitude descending for top-N selection.
                std::sort(pw_all.begin(), pw_all.end(),
                          [](const std::pair<double,double>& a,
                             const std::pair<double,double>& b) {
                              return a.second > b.second;
                          });

                // Edge-strip veto.
                // DUT sensor: edges = strip1, strip7.
                // Tracking planes: edges = strip1, strip8.
                if (skip_edge) {
                    int leading_idx = (int)std::round(
                        (pw_all[0].first - x_first_mm) / pitch_mm) + 1;
                    int edge_max = (grpname == "DUT sensor") ? 7 : 8;
                    if (leading_idx == 1 || leading_idx == edge_max) continue;
                }

                const double max_amp = pw_all[0].second;

                auto compute = [&](const std::vector<std::pair<double,double>>& pw,
                                   const std::string& tag_lin,
                                   const std::string& tag_sq) {
                    double pos_lin = Centroid(pw, 0);
                    double pos_sq  = Centroid(pw, 1);
                    if (pos_lin > -998.) {
                        getH1(grpname, tag_lin)->Fill(pos_lin);
                        event_pos[grpname][tag_lin] = pos_lin;
                        getH2PA(grpname, tag_lin)->Fill(pos_lin, max_amp);
                        getProfPA(grpname, tag_lin)->Fill(pos_lin, max_amp);
                    }
                    if (pos_sq > -998.) {
                        getH1(grpname, tag_sq)->Fill(pos_sq);
                        event_pos[grpname][tag_sq] = pos_sq;
                        getH2PA(grpname, tag_sq)->Fill(pos_sq, max_amp);
                        getProfPA(grpname, tag_sq)->Fill(pos_sq, max_amp);
                    }
                };

                compute({ pw_all.begin(),
                          pw_all.begin() + std::min((size_t)2, pw_all.size()) },
                        "top2_lin", "top2_sq");
                compute({ pw_all.begin(),
                          pw_all.begin() + std::min((size_t)3, pw_all.size()) },
                        "top3_lin", "top3_sq");
                compute(pw_all, "thr_lin", "thr_sq");

                if (pw_all.size() >= 3) {
                    std::vector<std::pair<double,double>> pw_bypos = pw_all;
                    std::sort(pw_bypos.begin(), pw_bypos.end());
                    double x_lo = pw_bypos.front().first - 0.5 * pitch_mm;
                    double x_hi = pw_bypos.back().first  + 0.5 * pitch_mm;
                    std::vector<double> gx, gy;
                    for (const auto& p : pw_bypos) { gx.push_back(p.first); gy.push_back(p.second); }
                    fg_gaus->SetParLimits(2, 0.1 * pitch_mm, 5.0 * pitch_mm);
                    fg_gaus->SetRange(x_lo, x_hi);
                    fg_gaus->SetParameters(max_amp, Centroid(pw_all, 0), pitch_mm);
                    TGraph tgfit((int)gx.size(), gx.data(), gy.data());
                    if (tgfit.Fit(fg_gaus, "QNS") == 0) {
                        double pos_gaus = fg_gaus->GetParameter(1);
                        if (pos_gaus > pos_min_mm && pos_gaus < pos_max_mm) {
                            getH1(grpname, "gaus")->Fill(pos_gaus);
                            event_pos[grpname]["gaus"] = pos_gaus;
                            getH2PA(grpname, "gaus")->Fill(pos_gaus, max_amp);
                            getProfPA(grpname, "gaus")->Fill(pos_gaus, max_amp);
                        }
                    }
                }
            }  // end !use_4cols

            // ── 4-cols mode: DUT sensor only ──────────────────────────────────
            // Col-1: strips 3,4,5,6 → sub-strip 1,2,3,4
            // Col-2: strips 1,2,7   → sub-strip 1,2,3  (no strip8 on DUT)
            if (use_4cols && grpname == "DUT sensor") {
                std::vector<std::pair<double,double>> pw_col[2];  // [0]=col1, [1]=col2

                // Re-use the already-read amplitudes in pw_all (keyed by strip idx).
                // But pw_all stores (position-in-2cols, amp); we need to rebuild
                // from scratch using sub-strip position.  Re-iterate channel keys.
                TIter itc2(grpdir->GetListOfKeys());
                TKey* ck2;
                while ((ck2 = (TKey*)itc2())) {
                    TObject* obj2 = ck2->ReadObj();
                    TGraph*  gr2  = dynamic_cast<TGraph*>(obj2);
                    if (!gr2) { delete obj2; continue; }

                    double amp_mV2 = 0.0;
                    int n2 = gr2->GetN();
                    if (n2 > 0) {
                        const double* Y2 = gr2->GetY();
                        amp_mV2 = -(*std::min_element(Y2, Y2 + n2)) * 1e3;
                    }
                    delete gr2;

                    if (amp_mV2 <= pos_thr_mV) continue;

                    int strip_idx = ChIndex(std::string(ck2->GetName()));
                    std::pair<int,int> cm = ColMap(strip_idx);
                    if (cm.first == 0) continue;   // unmapped (e.g. strip8 on DUT)

                    double x_sub = x_first_mm + (cm.second - 1) * pitch_mm;
                    pw_col[cm.first - 1].push_back({x_sub, amp_mV2});
                }

                // Process each column independently
                for (int col = 0; col < 2; ++col) {
                    std::string colname = grpname + "_col" + std::to_string(col + 1);

                    getNStripsH(colname)->Fill((int)pw_col[col].size());
                    if (pw_col[col].empty()) continue;

                    // Sort by amplitude descending
                    std::sort(pw_col[col].begin(), pw_col[col].end(),
                              [](const std::pair<double,double>& a,
                                 const std::pair<double,double>& b) {
                                  return a.second > b.second;
                              });

                    // Edge veto: sub-strip 1 and sub-strip 4 are the edges
                    if (skip_edge) {
                        int lead = (int)std::round(
                            (pw_col[col][0].first - x_first_mm) / pitch_mm) + 1;
                        if (lead == 1 || lead == 4) continue;
                    }

                    const double col_max_amp = pw_col[col][0].second;

                    auto compute_col = [&](const std::vector<std::pair<double,double>>& pw,
                                           const std::string& tag_lin,
                                           const std::string& tag_sq) {
                        double pos_lin = Centroid(pw, 0);
                        double pos_sq  = Centroid(pw, 1);
                        if (pos_lin > -998.) {
                            getH1(colname, tag_lin)->Fill(pos_lin);
                            event_pos[colname][tag_lin] = pos_lin;
                            getH2PA(colname, tag_lin)->Fill(pos_lin, col_max_amp);
                            getProfPA(colname, tag_lin)->Fill(pos_lin, col_max_amp);
                        }
                        if (pos_sq > -998.) {
                            getH1(colname, tag_sq)->Fill(pos_sq);
                            event_pos[colname][tag_sq] = pos_sq;
                            getH2PA(colname, tag_sq)->Fill(pos_sq, col_max_amp);
                            getProfPA(colname, tag_sq)->Fill(pos_sq, col_max_amp);
                        }
                    };

                    compute_col({ pw_col[col].begin(),
                                  pw_col[col].begin() + std::min((size_t)2, pw_col[col].size()) },
                                "top2_lin", "top2_sq");
                    compute_col({ pw_col[col].begin(),
                                  pw_col[col].begin() + std::min((size_t)3, pw_col[col].size()) },
                                "top3_lin", "top3_sq");
                    compute_col(pw_col[col], "thr_lin", "thr_sq");

                    if (pw_col[col].size() >= 3) {
                        std::vector<std::pair<double,double>> pw_bypos = pw_col[col];
                        std::sort(pw_bypos.begin(), pw_bypos.end());
                        double x_lo = pw_bypos.front().first - 0.5 * pitch_mm;
                        double x_hi = pw_bypos.back().first  + 0.5 * pitch_mm;
                        std::vector<double> gx, gy;
                        for (const auto& p : pw_bypos) { gx.push_back(p.first); gy.push_back(p.second); }
                        fg_gaus->SetParLimits(2, 0.1 * pitch_mm, 5.0 * pitch_mm);
                        fg_gaus->SetRange(x_lo, x_hi);
                        fg_gaus->SetParameters(col_max_amp, Centroid(pw_col[col], 0), pitch_mm);
                        TGraph tgfit((int)gx.size(), gx.data(), gy.data());
                        if (tgfit.Fit(fg_gaus, "QNS") == 0) {
                            double pos_gaus = fg_gaus->GetParameter(1);
                            if (pos_gaus > pos_min_mm && pos_gaus < pos_max_mm) {
                                getH1(colname, "gaus")->Fill(pos_gaus);
                                event_pos[colname]["gaus"] = pos_gaus;
                                getH2PA(colname, "gaus")->Fill(pos_gaus, col_max_amp);
                                getProfPA(colname, "gaus")->Fill(pos_gaus, col_max_amp);
                            }
                        }
                    }
                }  // end column loop
            }  // end use_4cols
        }

        // ── 2D correlations ───────────────────────────────────────────────────
        // Fill only when both groups have a valid hit in this event.
        for (int ic = 0; ic < nCorr; ++ic) {
            const CorrDef& cd = corrArr[ic];
            auto itx = event_pos.find(cd.xgrp);
            auto ity = event_pos.find(cd.ygrp);
            if (itx == event_pos.end() || ity == event_pos.end()) continue;

            for (int it = 0; it < kNTags; ++it) {
                std::string tag = kTags[it];
                auto px = itx->second.find(tag);
                auto py = ity->second.find(tag);
                if (px == itx->second.end() || py == ity->second.end()) continue;
                getH2(cd.label, tag, cd.xgrp, cd.ygrp)->Fill(px->second, py->second);
            }
        }

        // ── Position differences ──────────────────────────────────────────────
        for (int id = 0; id < nDiff; ++id) {
            const DiffDef& dd = diffArr[id];
            auto ita = event_pos.find(dd.agrp);
            auto itb = event_pos.find(dd.bgrp);
            if (ita == event_pos.end() || itb == event_pos.end()) continue;

            for (int it = 0; it < kNTags; ++it) {
                std::string tag = kTags[it];
                auto pa = ita->second.find(tag);
                auto pb = itb->second.find(tag);
                if (pa == ita->second.end() || pb == itb->second.end()) continue;
                double diff = pa->second - pb->second;
                if (diff == 0.0) continue;
                getH1Diff(dd.label, tag, dd.agrp, dd.bgrp)->Fill(diff);
            }
        }

        if ((iev + 1) % 200 == 0 || iev + 1 == evdirs.size()) {
            printf("[INFO] %zu / %zu events processed\r", iev + 1, evdirs.size());
            fflush(stdout);
        }
    }
    putchar('\n');

    delete fg_gaus;
    fin->Close();

    // ── Write histograms to output file ───────────────────────────────────────
    TFile* fout = TFile::Open(output_path, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output: " << output_path << "\n";
        return;
    }

    // One sub-directory per detector group for 1D histograms.
    std::map<std::string, TDirectory*> outdirs;
    auto getOutDir = [&](const std::string& grp) -> TDirectory* {
        auto it = outdirs.find(grp);
        if (it != outdirs.end()) return it->second;
        fout->cd();
        TDirectory* d = fout->mkdir(grp.c_str());
        outdirs[grp] = d;
        return d;
    };

    for (auto& kv : h1) {
        size_t slash = kv.first.find('/');
        std::string grp = (slash != std::string::npos)
                          ? kv.first.substr(0, slash) : kv.first;
        getOutDir(grp)->cd();
        kv.second->Write();
    }
    for (auto& kv : h_nstrips) {
        getOutDir(kv.first)->cd();
        kv.second->Write();
    }
    for (auto& kv : h2pa) {
        size_t slash = kv.first.find('/');
        std::string grp = (slash != std::string::npos)
                          ? kv.first.substr(0, slash) : kv.first;
        getOutDir(grp)->cd();
        kv.second->Write();
    }
    for (auto& kv : hppa) {
        size_t slash = kv.first.find('/');
        std::string grp = (slash != std::string::npos)
                          ? kv.first.substr(0, slash) : kv.first;
        getOutDir(grp)->cd();
        kv.second->Write();
    }

    // 2D histograms go into "Correlations/<corrLabel>/" sub-directories.
    fout->cd();
    TDirectory* corrTop = fout->mkdir("Correlations");
    std::map<std::string, TDirectory*> corrdirs;
    for (auto& kv : h2) {
        size_t slash = kv.first.find('/');
        std::string label = (slash != std::string::npos)
                            ? kv.first.substr(0, slash) : kv.first;
        if (corrdirs.find(label) == corrdirs.end()) {
            corrTop->cd();
            corrdirs[label] = corrTop->mkdir(label.c_str());
        }
        corrdirs[label]->cd();
        kv.second->Write();
    }

    // 1D difference histograms go into "Differences/<diffLabel>/" sub-directories.
    fout->cd();
    TDirectory* diffTop = fout->mkdir("Differences");
    std::map<std::string, TDirectory*> diffdirs;
    for (auto& kv : h1diff) {
        size_t slash = kv.first.find('/');
        std::string label = (slash != std::string::npos)
                            ? kv.first.substr(0, slash) : kv.first;
        if (diffdirs.find(label) == diffdirs.end()) {
            diffTop->cd();
            diffdirs[label] = diffTop->mkdir(label.c_str());
        }
        diffdirs[label]->cd();
        kv.second->Write();
    }

    fout->Close();
    printf("[INFO] Output file      : %s\n", output_path);
    printf("[INFO] 1D pos           : %zu\n", h1.size());
    printf("[INFO] N-strips hists   : %zu\n", h_nstrips.size());
    printf("[INFO] 2D pos×amp       : %zu\n", h2pa.size());
    printf("[INFO] Profile amp(pos) : %zu\n", hppa.size());
    printf("[INFO] 2D pos corr      : %zu\n", h2.size());
    printf("[INFO] 1D pos diff      : %zu\n", h1diff.size());

    for (auto& kv : h1)        delete kv.second;
    for (auto& kv : h_nstrips) delete kv.second;
    for (auto& kv : h2pa)      delete kv.second;
    for (auto& kv : hppa)      delete kv.second;
    for (auto& kv : h2)        delete kv.second;
    for (auto& kv : h1diff)    delete kv.second;
}
