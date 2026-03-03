/**
 * draw_hit_position.C
 *
 * Reads the hit-position ROOT file produced by make_hit_position.C
 * and opens interactive canvases for visual inspection.
 *
 * Canvas layout
 * ─────────────
 *  c_pos_<group>   : 1D position distributions per group (all 6 methods overlaid)
 *  c_nstrips       : strip multiplicity — all groups (3×2 pads)
 *  c_posamp_<group>: 2D position vs max-strip amplitude per group [method_tag]
 *  c_corr_1, _2,…  : 2D position correlations, 4 pairs per canvas [method_tag]
 *  c_diff          : 1D position differences — all 4 pairs (2×2 pads) [method_tag]
 *
 * ── Usage ───────────────────────────────────────────────────────────────────
 *  root -l 'draw_hit_position.C("91232923")'
 *  root -l 'draw_hit_position.C("91232923","thr_lin")'
 *  root -l 'draw_hit_position.C("91232923","top2_lin","/path/to/files")'
 *
 * ── Parameters ──────────────────────────────────────────────────────────────
 *  run_number : run number string (XXXXXX in pos_XXXXXX.root)
 *  method_tag : method for 2D and difference plots     (default "top2_lin")
 *               one of: top2_lin top2_sq top3_lin top3_sq thr_lin thr_sq
 *  input_dir  : directory containing pos_XXXXXX.root   (default ".")
 * ────────────────────────────────────────────────────────────────────────────
 */

#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPad.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

void draw_hit_position(
    const char* run_number,
    const char* method_tag = "top2_lin",
    const char* input_dir  = ".")
{
    // ── Style ─────────────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPalette(kBird);
    gROOT->SetBatch(kFALSE);

    // ── Open input file ───────────────────────────────────────────────────────
    std::string in_path = std::string(input_dir) + "/pos_" + run_number + ".root";

    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open: " << in_path << "\n";
        return;
    }
    printf("[INFO] Input  : %s\n", in_path.c_str());
    printf("[INFO] Method : %s  (used for 2D and difference plots)\n\n", method_tag);

    // ── Detector groups ───────────────────────────────────────────────────────
    const char* kGroups[] = {
        "DUT sensor", "MCP-PMT",
        "Tracking_0", "Tracking_1", "Tracking_2", "Tracking_3"
    };
    const int kNGroups = 6;

    // ── Method tags ───────────────────────────────────────────────────────────
    const char* kTags[] = {
        "top2_lin", "top2_sq",
        "top3_lin", "top3_sq",
        "thr_lin",  "thr_sq"
    };
    const int kNTags = 6;

    // ── Correlation pairs ─────────────────────────────────────────────────────
    struct CorrDef { const char* label; const char* xtitle; const char* ytitle; };
    const CorrDef kCorr[] = {
        { "DUT_vs_Trk0",   "DUT sensor [mm]",  "Tracking_0 [mm]"  },
        { "DUT_vs_Trk1",   "DUT sensor [mm]",  "Tracking_1 [mm]"  },
        { "DUT_vs_Trk2",   "DUT sensor [mm]",  "Tracking_2 [mm]"  },
        { "DUT_vs_Trk3",   "DUT sensor [mm]",  "Tracking_3 [mm]"  },
        { "Trk0_vs_Trk1",  "Tracking_0 [mm]",  "Tracking_1 [mm]"  },
        { "Trk0_vs_Trk2",  "Tracking_0 [mm]",  "Tracking_2 [mm]"  },
        { "Trk0_vs_Trk3",  "Tracking_0 [mm]",  "Tracking_3 [mm]"  },
        { "Trk1_vs_Trk2",  "Tracking_1 [mm]",  "Tracking_2 [mm]"  },
        { "Trk1_vs_Trk3",  "Tracking_1 [mm]",  "Tracking_3 [mm]"  },
        { "Trk2_vs_Trk3",  "Tracking_2 [mm]",  "Tracking_3 [mm]"  },
    };
    const int kNCorr = 10;

    // ── Difference pairs ──────────────────────────────────────────────────────
    struct DiffDef { const char* label; const char* title; };
    const DiffDef kDiff[] = {
        { "DUT_minus_Trk0",  "DUT sensor #minus Tracking_0"  },
        { "DUT_minus_Trk3",  "DUT sensor #minus Tracking_3"  },
        { "Trk0_minus_Trk3", "Tracking_0 #minus Tracking_3"  },
        { "Trk1_minus_Trk2", "Tracking_1 #minus Tracking_2"  },
    };
    const int kNDiff = 4;

    // ── Colour / style for 6 methods ──────────────────────────────────────────
    const Color_t kMethodCol[] = {
        kBlue, kBlue+2,
        kRed,  kRed+2,
        kGreen+2, kGreen+4
    };
    const Style_t kMethodSty[] = { 1, 2, 1, 2, 1, 2 };

    int n_canvas = 0;

    // ── §1  1D position histograms — one canvas per group ────────────────────
    printf("[INFO] §1  1D position histograms...\n");

    for (int g = 0; g < kNGroups; ++g) {
        TDirectory* grpdir = nullptr;
        fin->GetObject(kGroups[g], grpdir);
        if (!grpdir) continue;

        std::string grp_norm = kGroups[g];
        for (char& ch : grp_norm) if (ch == ' ' || ch == '-') ch = '_';

        std::string cname  = "c_pos_" + grp_norm;
        std::string ctitle = "Run " + std::string(run_number)
                           + "  " + kGroups[g] + "  — position";
        TCanvas* cv = new TCanvas(cname.c_str(), ctitle.c_str(), 1000, 600);
        cv->SetLogy(0);
        cv->SetLeftMargin(0.12);
        cv->SetRightMargin(0.20);
        cv->SetBottomMargin(0.13);
        cv->SetTopMargin(0.10);
        cv->SetGrid();

        TLegend* leg = new TLegend(0.81, 0.15, 0.99, 0.90);
        leg->SetTextSize(0.038);
        leg->SetBorderSize(1);

        double gmax = 0.0;
        std::vector<TH1D*> hv;

        for (int it = 0; it < kNTags; ++it) {
            std::string hname = "h_pos_" + grp_norm + "_" + kTags[it];
            TH1D* h = nullptr;
            grpdir->GetObject(hname.c_str(), h);
            if (!h) { hv.push_back(nullptr); continue; }

            TH1D* hc = (TH1D*)h->Clone();
            hc->SetDirectory(nullptr);
            hc->SetLineColor(kMethodCol[it]);
            hc->SetLineStyle(kMethodSty[it]);
            hc->SetLineWidth(2);
            if (hc->GetMaximum() > gmax) gmax = hc->GetMaximum();
            hv.push_back(hc);
        }

        bool first = true;
        for (int it = 0; it < kNTags; ++it) {
            if (!hv[it]) continue;
            hv[it]->SetMaximum(gmax * 3.0);
            hv[it]->SetMinimum(0.5);
            if (first) {
                hv[it]->SetTitle(
                    Form("Run %s  %s;Position [mm];Entries", run_number, kGroups[g]));
                hv[it]->GetXaxis()->SetTitleSize(0.05);
                hv[it]->GetXaxis()->SetLabelSize(0.045);
                hv[it]->GetYaxis()->SetTitleSize(0.05);
                hv[it]->GetYaxis()->SetLabelSize(0.045);
                hv[it]->Draw("HIST");
                first = false;
            } else {
                hv[it]->Draw("HIST SAME");
            }
            leg->AddEntry(hv[it], kTags[it], "l");
        }
        if (!first) {
            leg->Draw();
            cv->Modified(); cv->Update();
            printf("[INFO]   %s\n", ctitle.c_str());
            ++n_canvas;
        }
    }

    // ── §2  Strip multiplicity — all groups, one canvas (3×2) ────────────────
    printf("[INFO] §2  Strip multiplicities...\n");
    {
        TCanvas* cv = new TCanvas(
            "c_nstrips",
            Form("Run %s — strip multiplicity", run_number),
            1200, 700);
        cv->Divide(3, 2, 0.01, 0.01);

        bool any = false;
        for (int g = 0; g < kNGroups; ++g) {
            TDirectory* grpdir = nullptr;
            fin->GetObject(kGroups[g], grpdir);
            if (!grpdir) continue;

            std::string grp_norm = kGroups[g];
            for (char& ch : grp_norm) if (ch == ' ' || ch == '-') ch = '_';

            TH1D* h = nullptr;
            grpdir->GetObject(("h_nstrips_" + grp_norm).c_str(), h);
            if (!h) continue;

            TVirtualPad* pad = cv->cd(g + 1);
            pad->SetLeftMargin(0.14);
            pad->SetBottomMargin(0.15);
            pad->SetTopMargin(0.12);
            pad->SetGrid();

            TH1D* hc = (TH1D*)h->Clone();
            hc->SetDirectory(nullptr);
            hc->SetFillColor(kAzure + 7);
            hc->SetFillStyle(1001);
            hc->SetLineColor(kAzure - 3);
            hc->SetTitle(Form("%s;N strips > threshold;Entries", kGroups[g]));
            hc->GetXaxis()->SetTitleSize(0.055);
            hc->GetXaxis()->SetLabelSize(0.05);
            hc->GetYaxis()->SetTitleSize(0.055);
            hc->GetYaxis()->SetLabelSize(0.05);
            hc->Draw("HIST");
            any = true;
        }
        if (any) { cv->Modified(); cv->Update(); ++n_canvas; }
    }

    // ── §3  2D pos vs amplitude — one canvas per group ───────────────────────
    printf("[INFO] §3  2D pos vs amplitude (method: %s)...\n", method_tag);

    for (int g = 0; g < kNGroups; ++g) {
        TDirectory* grpdir = nullptr;
        fin->GetObject(kGroups[g], grpdir);
        if (!grpdir) continue;

        std::string grp_norm = kGroups[g];
        for (char& ch : grp_norm) if (ch == ' ' || ch == '-') ch = '_';

        TH2D* h = nullptr;
        grpdir->GetObject(("h2pa_" + grp_norm + "_" + method_tag).c_str(), h);
        if (!h) {
            printf("[WARN]   h2pa not found for '%s' [%s]\n", kGroups[g], method_tag);
            continue;
        }

        std::string cname  = "c_posamp_" + grp_norm;
        std::string ctitle = "Run " + std::string(run_number)
                           + "  " + kGroups[g] + "  pos vs amp  [" + method_tag + "]";
        TCanvas* cv = new TCanvas(cname.c_str(), ctitle.c_str(), 900, 700);
        cv->SetLeftMargin(0.12);
        cv->SetRightMargin(0.14);
        cv->SetBottomMargin(0.13);
        cv->SetTopMargin(0.10);

        TH2D* hc = (TH2D*)h->Clone();
        hc->SetDirectory(nullptr);
        hc->SetTitle(
            Form("Run %s  %s  [%s];Position [mm];Max-strip amplitude [mV]",
                 run_number, kGroups[g], method_tag));
        hc->GetXaxis()->SetTitleSize(0.05);
        hc->GetXaxis()->SetLabelSize(0.045);
        hc->GetYaxis()->SetTitleSize(0.05);
        hc->GetYaxis()->SetLabelSize(0.045);
        hc->Draw("COLZ");
        cv->Modified(); cv->Update();
        printf("[INFO]   %s\n", ctitle.c_str());
        ++n_canvas;
    }

    // ── §4  2D correlations — 4 pairs per canvas ─────────────────────────────
    printf("[INFO] §4  2D correlations (method: %s)...\n", method_tag);
    {
        TDirectory* corrTop = nullptr;
        fin->GetObject("Correlations", corrTop);

        if (!corrTop) {
            printf("[WARN]   'Correlations' directory not found\n");
        } else {
            const int kPerPage = 4;
            int cpage = 0;
            for (int ic = 0; ic < kNCorr; ic += kPerPage) {
                ++cpage;
                std::string cname  = Form("c_corr_%d", cpage);
                std::string ctitle = Form("Run %s — correlations [%s] (%d)",
                                         run_number, method_tag, cpage);
                TCanvas* cv = new TCanvas(cname.c_str(), ctitle.c_str(), 1200, 900);
                cv->Divide(2, 2, 0.01, 0.01);

                bool any = false;
                for (int ip = 0; ip < kPerPage && (ic + ip) < kNCorr; ++ip) {
                    const CorrDef& cd = kCorr[ic + ip];
                    TDirectory* pairdir = nullptr;
                    corrTop->GetObject(cd.label, pairdir);
                    if (!pairdir) continue;

                    TH2D* h = nullptr;
                    pairdir->GetObject(
                        (std::string("h2_") + cd.label + "_" + method_tag).c_str(), h);
                    if (!h) continue;

                    TVirtualPad* pad = cv->cd(ip + 1);
                    pad->SetLeftMargin(0.14);
                    pad->SetRightMargin(0.13);
                    pad->SetBottomMargin(0.15);
                    pad->SetTopMargin(0.12);
                    pad->SetLogz(0);

                    TH2D* hc = (TH2D*)h->Clone();
                    hc->SetDirectory(nullptr);
                    hc->SetTitle(
                        Form("Run %s  %s  [%s];%s;%s",
                             run_number, cd.label, method_tag,
                             cd.xtitle, cd.ytitle));
                    hc->GetXaxis()->SetTitleSize(0.055);
                    hc->GetXaxis()->SetLabelSize(0.05);
                    hc->GetYaxis()->SetTitleSize(0.055);
                    hc->GetYaxis()->SetLabelSize(0.05);
                    hc->Draw("COLZ");
                    any = true;
                }
                if (any) { cv->Modified(); cv->Update(); ++n_canvas; }
            }
        }
    }

    // ── §5  Position differences — all 4 pairs, one canvas ───────────────────
    printf("[INFO] §5  Position differences (method: %s)...\n", method_tag);
    {
        TDirectory* diffTop = nullptr;
        fin->GetObject("Differences", diffTop);

        if (!diffTop) {
            printf("[WARN]   'Differences' directory not found\n");
        } else {
            TCanvas* cv = new TCanvas(
                "c_diff",
                Form("Run %s — position differences [%s]", run_number, method_tag),
                1200, 900);
            cv->Divide(2, 2, 0.01, 0.01);

            bool any = false;
            for (int id = 0; id < kNDiff; ++id) {
                const DiffDef& dd = kDiff[id];
                TDirectory* pairdir = nullptr;
                diffTop->GetObject(dd.label, pairdir);
                if (!pairdir) continue;

                TH1D* h = nullptr;
                pairdir->GetObject(
                    (std::string("hdiff_") + dd.label + "_" + method_tag).c_str(), h);
                if (!h) continue;

                TVirtualPad* pad = cv->cd(id + 1);
                pad->SetLeftMargin(0.14);
                pad->SetBottomMargin(0.15);
                pad->SetTopMargin(0.12);
                pad->SetGrid();

                TH1D* hc = (TH1D*)h->Clone();
                hc->SetDirectory(nullptr);
                hc->SetFillColor(kOrange - 3);
                hc->SetFillStyle(1001);
                hc->SetLineColor(kOrange + 3);
                hc->SetTitle(
                    Form("Run %s  %s  [%s];#Delta pos [mm];Entries",
                         run_number, dd.title, method_tag));
                hc->GetXaxis()->SetTitleSize(0.055);
                hc->GetXaxis()->SetLabelSize(0.05);
                hc->GetYaxis()->SetTitleSize(0.055);
                hc->GetYaxis()->SetLabelSize(0.05);
                hc->Draw("HIST");
                any = true;
            }
            if (any) { cv->Modified(); cv->Update(); ++n_canvas; }
        }
    }

    // ── Done ──────────────────────────────────────────────────────────────────
    fin->Close();
    printf("\n[INFO] Done — %d canvas(es) open.\n", n_canvas);
    // All canvases remain alive for interactive inspection.
    // Close them manually or type .q to exit ROOT.
}
