/**
 * draw_event_waveforms.C
 *
 * Reads per-event TGraph waveforms from a waveform ROOT file produced by
 * make_event_waveforms.C and displays all channels for a chosen event.
 *
 * One canvas per detector group (= one oscilloscope),
 * all channels overlaid with different colours and a legend:
 *
 *   c_DUT_sensor   :  DUT sensor   strip1..strip7  (7 overlaid)
 *   c_MCP_PMT      :  MCP-PMT      wf              (1 channel)
 *   c_Tracking_1   :  Tracking_1   strip1..strip8  (8 overlaid)
 *   c_Tracking_2   :  Tracking_2   strip1..strip8  (8 overlaid)
 *   c_Tracking_4   :  Tracking_4   strip1..strip8  (8 overlaid)
 *   c_Tracking_5   :  Tracking_5   strip1..strip8  (8 overlaid)
 *
 * EUDAQ_ID → detector mapping (from make_event_waveforms.C):
 *   EUDAQ_ID=4 → Tracking_1   C1..C8 → strip1..strip8
 *   EUDAQ_ID=1 → Tracking_2   C1..C8 → strip1..strip8
 *   EUDAQ_ID=3 → DUT sensor   C1..C7 → strip1..strip7
 *              → MCP-PMT      C8     → wf
 *   EUDAQ_ID=2 → Tracking_4   C1..C8 → strip1..strip8
 *   EUDAQ_ID=0 → Tracking_5   C1..C8 → strip1..strip8
 *
 * X axis : time [ns]  (TGraph X × t_scale)
 * Y axis : voltage [mV] (TGraph Y × 1000)
 * A dashed zero line and a legend are drawn on each canvas.
 *
 * ── Usage ────────────────────────────────────────────────────────────────────
 *  root -l 'draw_event_waveforms.C("waveform_102134314_260303134320.root")'
 *  root -l 'draw_event_waveforms.C("waveform_102134314_260303134320.root",5)'
 *  root -l 'draw_event_waveforms.C("waveform_102134314_260303134320.root",5,1.0)'
 *
 * ── Parameters ───────────────────────────────────────────────────────────────
 *  input_path   : waveform ROOT file (output of make_event_waveforms.C)
 *  event_number : event index to display (0-based, default 0)
 *  t_scale      : multiply TGraph X to get ns  (default 1e9, s→ns)
 *                 set to 1.0 if TGraph X is already in ns
 * ────────────────────────────────────────────────────────────────────────────
 */

#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void draw_event_waveforms(
    const char* input_path,
    Int_t       event_number = 0,
    Double_t    t_scale      = 1e9)   // multiply TGraph X to get ns (1e9: s→ns)
{
    // ── Style ─────────────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gROOT->SetBatch(kFALSE);

    // ── Open input file ───────────────────────────────────────────────────────
    TFile* fin = TFile::Open(input_path, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open: " << input_path << "\n";
        return;
    }

    // ── Locate event directory ────────────────────────────────────────────────
    char evname[32];
    snprintf(evname, sizeof(evname), "event_%06d", (int)event_number);

    TDirectory* evdir = nullptr;
    fin->GetObject(evname, evdir);
    if (!evdir) {
        fprintf(stderr, "[ERROR] '%s' not found in %s\n", evname, input_path);
        fin->Close();
        return;
    }

    printf("[INFO] File       : %s\n", input_path);
    printf("[INFO] Event      : %s\n", evname);
    printf("[INFO] Time scale : %.4g  (TGraph X × scale = ns)\n\n", t_scale);

    // ── Detector groups ───────────────────────────────────────────────────────
    const char* kGroups[] = {
        "DUT sensor", "MCP-PMT",
        "Tracking_1", "Tracking_2", "Tracking_4", "Tracking_5"
    };
    const int kNGroups = 6;

    // ── Colour palette (one colour per channel) ───────────────────────────────
    const Color_t kCols[] = {
        kBlue, kRed, kGreen+2, kMagenta+1,
        kCyan+2, kOrange+1, kViolet+1, kTeal+1
    };
    const int kNCols = 8;

    int n_canvas = 0;

    // ── Loop over groups ──────────────────────────────────────────────────────
    for (int g = 0; g < kNGroups; ++g) {
        const char* grpname = kGroups[g];

        TDirectory* grpdir = evdir->GetDirectory(grpname);
        if (!grpdir) {
            printf("[WARN]   '%s' not found in %s — skipping\n", grpname, evname);
            continue;
        }

        // Collect channels as (name, TGraph*), sorted by name
        std::vector<std::pair<std::string, TGraph*>> channels;
        {
            TIter it(grpdir->GetListOfKeys());
            TKey* k;
            while ((k = (TKey*)it())) {
                TObject* obj = k->ReadObj();
                TGraph* gr = dynamic_cast<TGraph*>(obj);
                if (!gr) { delete obj; continue; }
                channels.push_back({ std::string(k->GetName()), gr });
            }
        }
        if (channels.empty()) continue;

        std::sort(channels.begin(), channels.end(),
                  [](const std::pair<std::string,TGraph*>& a,
                     const std::pair<std::string,TGraph*>& b) {
                      return a.first < b.first;
                  });

        // ── Convert all channels to ns / mV, find global range ────────────────
        struct ChData { std::string name; std::vector<double> xns, ymv; };
        std::vector<ChData> cdata;
        double tmin =  1e18, tmax = -1e18;
        double ymin_mV =  1e18, ymax_mV = -1e18;

        for (auto& ch : channels) {
            TGraph* gr = ch.second;
            const int     n = gr->GetN();
            const double* X = gr->GetX();
            const double* Y = gr->GetY();
            if (n <= 0) { delete gr; continue; }

            ChData cd;
            cd.name = ch.first;
            cd.xns.resize(n);
            cd.ymv.resize(n);
            for (int i = 0; i < n; ++i) {
                cd.xns[i] = X[i] * t_scale;
                cd.ymv[i] = Y[i] * 1e3;
                if (cd.xns[i] < tmin) tmin = cd.xns[i];
                if (cd.xns[i] > tmax) tmax = cd.xns[i];
                if (cd.ymv[i] < ymin_mV) ymin_mV = cd.ymv[i];
                if (cd.ymv[i] > ymax_mV) ymax_mV = cd.ymv[i];
            }
            cdata.push_back(std::move(cd));
            delete gr;
        }
        if (cdata.empty()) continue;

        // Y range with 15 % margin; keep zero visible
        double yrange = ymax_mV - ymin_mV;
        if (yrange < 1.0) yrange = 1.0;
        double ylo = ymin_mV - 0.15 * yrange;
        double yhi = ymax_mV + 0.15 * yrange;
        if (yhi < 0.0) yhi =  0.1 * yrange;
        if (ylo > 0.0) ylo = -0.1 * yrange;

        // ── Create canvas ─────────────────────────────────────────────────────
        std::string grp_norm = grpname;
        for (char& c : grp_norm) if (c == ' ' || c == '-') c = '_';

        std::string cname  = "c_" + grp_norm;
        std::string ctitle = Form("Event %d  —  %s  (%d channels)",
                                  (int)event_number, grpname,
                                  (int)cdata.size());

        TCanvas* cv = new TCanvas(cname.c_str(), ctitle.c_str(), 1000, 600);
        cv->SetLeftMargin(0.12);
        cv->SetRightMargin(0.18);
        cv->SetBottomMargin(0.13);
        cv->SetTopMargin(0.10);
        cv->SetGrid();

        // Legend on the right
        TLegend* leg = new TLegend(0.83, 0.15, 0.99, 0.90);
        leg->SetTextSize(0.038);
        leg->SetBorderSize(1);

        // ── Draw all channels overlaid ────────────────────────────────────────
        for (int ci = 0; ci < (int)cdata.size(); ++ci) {
            ChData& cd = cdata[ci];
            const int n = (int)cd.xns.size();

            TGraph* gplot = new TGraph(n, cd.xns.data(), cd.ymv.data());
            gplot->SetName(cd.name.c_str());
            gplot->SetLineColor(kCols[ci % kNCols]);
            gplot->SetLineWidth(2);

            if (ci == 0) {
                gplot->SetTitle(Form("Event %d  %s;Time [ns];Voltage [mV]",
                                     (int)event_number, grpname));
                gplot->SetMinimum(ylo);
                gplot->SetMaximum(yhi);
                gplot->GetXaxis()->SetLimits(tmin, tmax);
                gplot->GetXaxis()->SetTitleSize(0.05);
                gplot->GetXaxis()->SetLabelSize(0.045);
                gplot->GetYaxis()->SetTitleSize(0.05);
                gplot->GetYaxis()->SetLabelSize(0.045);
                gplot->Draw("AL");
            } else {
                gplot->Draw("L SAME");
            }
            leg->AddEntry(gplot, cd.name.c_str(), "l");
        }

        // Dashed zero reference line
        TLine* lz = new TLine(tmin, 0., tmax, 0.);
        lz->SetLineColor(kGray + 1);
        lz->SetLineStyle(2);
        lz->SetLineWidth(1);
        lz->Draw();

        leg->Draw();
        cv->Modified();
        cv->Update();
        printf("[INFO]   %-12s  %d channel(s)\n", grpname, (int)cdata.size());
        ++n_canvas;
    }

    // ── Done ──────────────────────────────────────────────────────────────────
    fin->Close();
    printf("\n[INFO] Done — %d canvas(es) open for event %d.\n",
           n_canvas, (int)event_number);
    // Canvases remain alive for interactive inspection.
    // Type .q to exit ROOT or close canvases manually.
}
