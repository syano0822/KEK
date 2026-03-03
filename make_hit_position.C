/**
 * make_hit_position.C
 *
 * Reads per-event TGraph waveforms stored by make_event_waveforms.C and
 * computes the amplitude-weighted hit position for each detector group.
 *
 * Hit-position estimators:
 *   Linear  :  pos = Σ(w_i  · x_i) / Σ(w_i)     where w_i = |amplitude_i|
 *   Squared :  pos = Σ(w_i² · x_i) / Σ(w_i²)
 *
 * Strip positions:   x_i = x_first + (i−1) · pitch   [i = 1-based channel index]
 *   Default: ch1 = 14.5 mm, pitch = 0.5 mm  →  ch8 = 18.0 mm
 *
 * Three channel-selection methods are studied in parallel:
 *   top2 : 2 strips with highest amplitude  (from those > pos_thr_mV)
 *   top3 : 3 strips with highest amplitude  (from those > pos_thr_mV)
 *   thr  : all strips with amplitude > pos_thr_mV
 *
 * 1D output histograms per detector-group directory:
 *   h_pos_top2_lin,  h_pos_top2_sq
 *   h_pos_top3_lin,  h_pos_top3_sq
 *   h_pos_thr_lin,   h_pos_thr_sq
 *
 * 2D correlation histograms (directory "Correlations/"), all 15 pairs:
 *   DUT sensor × MCP-PMT, DUT sensor × Tracking_0..3,
 *   MCP-PMT    × Tracking_0..3,
 *   Tracking_0 × Tracking_1..3,
 *   Tracking_1 × Tracking_2, Tracking_1 × Tracking_3,
 *   Tracking_2 × Tracking_3
 *   Each pair: one TH2D per method tag → 15 × 6 = 90 histograms total.
 *   Filled only when both groups have a valid hit in the same event.
 *
 * 1D position-difference histograms (directory "Differences/"):
 *   DUT sensor − Tracking_0,  DUT sensor − Tracking_3,
 *   Tracking_0 − Tracking_3,  Tracking_1 − Tracking_2
 *   Each pair: one TH1D per method tag → 4 × 6 = 24 histograms total.
 *
 * 1D strip-multiplicity histograms (per detector-group directory):
 *   h_nstrips_<group>  : number of strips with amplitude > pos_thr_mV
 *   Range: 0..8 (9 bins, −0.5 to 8.5); filled for every event including zero.
 *
 * Channel name → index mapping (trailing digits, no digits → 1):
 *   "strip3" → 3,   "wf_C5" → 5,   "wf" → 1
 *
 * ── Usage ───────────────────────────────────────────────────────────────────
 *  Interpreted:
 *    root 'make_hit_position.C("waveforms.root")'
 *    root 'make_hit_position.C("waveforms.root","pos.root",8.,14.5,0.5,200,13.,19.)'
 *
 *  ACLiC compiled (recommended for large files):
 *    root 'make_hit_position.C+("waveforms.root")'
 *
 * ── Parameters ──────────────────────────────────────────────────────────────
 *  input_path  : waveform ROOT file (output of make_event_waveforms.C)
 *  output_path : output ROOT file                         (default "hit_position.root")
 *  pos_thr_mV  : min amplitude to include a channel [mV]  (default   8.0)
 *  x_first_mm  : position of channel 1 [mm]               (default  14.5)
 *  x_pitch_mm  : strip pitch [mm]                         (default   0.5)
 *  nbins       : histogram bins (1D and each 2D axis)      (default 200)
 *  pos_min_mm  : histogram x-axis minimum [mm]            (default  13.0)
 *  pos_max_mm  : histogram x-axis maximum [mm]            (default  19.0)
 *  amp_min_mV  : amplitude axis minimum for pos-vs-amp [mV] (default   0.0)
 *  amp_max_mV  : amplitude axis maximum for pos-vs-amp [mV] (default 300.0)
 *  skip_edge   : if true, skip events where ch1 or ch8 is the leading strip
 *                (edge-strip veto, applied per group)      (default false)
 *  diff_min_mm : position-difference histogram minimum [mm](default  -6.0)
 *  diff_max_mm : position-difference histogram maximum [mm](default  +6.0)
 * ────────────────────────────────────────────────────────────────────────────
 */

#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>

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
    "thr_lin",  "thr_sq"
};
static const int kNTags = 6;

// ── Extract 1-based channel index from a TGraph name ─────────────────────────
// "strip3" → 3,  "wf_C5" → 5,  "wf" → 1 (no trailing digits)
static int ChIndex(const std::string& name)
{
    int i = (int)name.size() - 1;
    while (i >= 0 && std::isdigit((unsigned char)name[i])) --i;
    ++i;
    if (i < (int)name.size()) return std::stoi(name.substr(i));
    return 1;
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
    Double_t    diff_max_mm =   6.0)  // position-difference histogram maximum [mm]
{
    // ── Open input file ───────────────────────────────────────────────────────
    TFile* fin = TFile::Open(input_path, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open input: " << input_path << "\n";
        return;
    }
    printf("[INFO] Input  file      : %s\n",  input_path);
    printf("[INFO] Amplitude thr    : > %.1f mV\n", pos_thr_mV);
    printf("[INFO] Strip positions  : ch1 = %.2f mm,  pitch = %.2f mm\n",
           x_first_mm, x_pitch_mm);
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
    // Two correlations, six tags each → 12 TH2D total.
    // key = "corrLabel/tag"  e.g. "EUDAQ5_vs_EUDAQ2/top2_lin"
    std::map<std::string, TH2D*> h2;

    // corrLabel, xGroup, yGroup  — all C(6,2) = 15 pairs
    struct CorrDef { const char* label; const char* xgrp; const char* ygrp; };
    const CorrDef kCorr[] = {
        // DUT sensor vs others
        { "DUT_vs_MCP",    "DUT sensor",  "MCP-PMT"    },
        { "DUT_vs_Trk0",   "DUT sensor",  "Tracking_0" },
        { "DUT_vs_Trk1",   "DUT sensor",  "Tracking_1" },
        { "DUT_vs_Trk2",   "DUT sensor",  "Tracking_2" },
        { "DUT_vs_Trk3",   "DUT sensor",  "Tracking_3" },
        // MCP-PMT vs Tracking planes
        { "MCP_vs_Trk0",   "MCP-PMT",     "Tracking_0" },
        { "MCP_vs_Trk1",   "MCP-PMT",     "Tracking_1" },
        { "MCP_vs_Trk2",   "MCP-PMT",     "Tracking_2" },
        { "MCP_vs_Trk3",   "MCP-PMT",     "Tracking_3" },
        // Tracking_0 vs Tracking planes
        { "Trk0_vs_Trk1",  "Tracking_0",  "Tracking_1" },
        { "Trk0_vs_Trk2",  "Tracking_0",  "Tracking_2" },
        { "Trk0_vs_Trk3",  "Tracking_0",  "Tracking_3" },
        // Tracking_1 vs Tracking_2,3
        { "Trk1_vs_Trk2",  "Tracking_1",  "Tracking_2" },
        { "Trk1_vs_Trk3",  "Tracking_1",  "Tracking_3" },
        // Tracking_2 vs Tracking_3
        { "Trk2_vs_Trk3",  "Tracking_2",  "Tracking_3" },
    };
    const int kNCorr = 15;

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
        { "DUT_minus_Trk0",  "DUT sensor", "Tracking_0" },
        { "DUT_minus_Trk3",  "DUT sensor", "Tracking_3" },
        { "Trk0_minus_Trk3", "Tracking_0", "Tracking_3" },
        { "Trk1_minus_Trk2", "Tracking_1", "Tracking_2" },
    };
    const int kNDiff = 4;

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
    // Bins: 0..8 strips above threshold (9 integer bins, −0.5 to 8.5)
    std::map<std::string, TH1D*> h_nstrips;

    auto getNStripsH = [&](const std::string& grp) -> TH1D* {
        auto it = h_nstrips.find(grp);
        if (it != h_nstrips.end()) return it->second;

        std::string hname = "h_nstrips_" + grp;
        for (char& c : hname) if (c == ' ' || c == '-') c = '_';

        std::string title = grp + " strips above threshold;"
                            "N strips above threshold;Entries";
        TH1D* h = new TH1D(hname.c_str(), title.c_str(), 9, -0.5, 8.5);
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
            TDirectory* grpdir  = evdir->GetDirectory(grpname.c_str());
            if (!grpdir) continue;
	    
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
                double x   = x_first_mm + (idx - 1) * x_pitch_mm;
                pw_all.push_back({x, amp_mV});
            }

            // Fill multiplicity histogram for every event (including zero strips).
            getNStripsH(grpname)->Fill((int)pw_all.size());

            if (pw_all.empty()) continue;

            // Sort by amplitude descending for top-N selection.
            std::sort(pw_all.begin(), pw_all.end(),
                      [](const std::pair<double,double>& a,
                         const std::pair<double,double>& b) {
                          return a.second > b.second;
                      });

            // Edge-strip veto: skip this group/event when the leading strip
            // (highest amplitude) is ch1 or ch8.
            if (skip_edge) {
                int leading_idx = (int)std::round(
                    (pw_all[0].first - x_first_mm) / x_pitch_mm) + 1;
                if (grpname == "DUT sensor") {
		  if (leading_idx == 1 || leading_idx == 7) continue;
		} else {
		  if (leading_idx == 1 || leading_idx == 8) continue;
		}		
            }

            // Leading-strip amplitude (highest-amplitude strip after sorting).
            // Used as Y-axis in pos-vs-amp histograms — same value for all methods.
            const double max_amp = pw_all[0].second;

            // Compute and store positions for both weighting schemes.
            // Also fills the 1D, 2D, and profile histograms immediately.
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
        }

        // ── 2D correlations ───────────────────────────────────────────────────
        // Fill only when both groups have a valid hit in this event.
        for (int ic = 0; ic < kNCorr; ++ic) {
            const CorrDef& cd = kCorr[ic];
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
        for (int id = 0; id < kNDiff; ++id) {
            const DiffDef& dd = kDiff[id];
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
