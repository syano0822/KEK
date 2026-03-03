/**
 * make_event_waveforms.C
 *
 * Reads oscilloscope waveforms from a TTree, applies pedestal subtraction,
 * and writes per-event TGraphs to an output ROOT file.
 *
 * TTree branch format:
 *   EventN/L, EUDAQ_ID/L
 *   nC#_volt/I,  C#_volt[nC#_volt]/D    (# = 1..8)
 *   nC#_time/I,  C#_time[nC#_time]/D
 *
 * EUDAQ_ID → detector mapping:
 *   EUDAQ_ID=0 → Tracking_5   C1..C8 → strip1..strip8
 *   EUDAQ_ID=1 → Tracking_2   C1..C8 → strip1..strip8
 *   EUDAQ_ID=2 → Tracking_4   C1..C8 → strip1..strip8
 *   EUDAQ_ID=3 → DUT sensor   C1..C7 → strip1..strip7
 *              → MCP-PMT      C8     → wf
 *   EUDAQ_ID=4 → Tracking_1   C1..C8 → strip1..strip8
 *
 * Output structure:
 *   /event_000000/
 *     DUT sensor/   strip1..strip7
 *     MCP-PMT/      wf
 *     Tracking_1/   strip1..strip8
 *     Tracking_2/   strip1..strip8
 *     Tracking_4/   strip1..strip8
 *     Tracking_5/   strip1..strip8
 *   /event_000001/ ...
 *
 * ── Usage ────────────────────────────────────────────────────────────────────
 *  Interpreted:
 *    root 'make_event_waveforms.C("in.root","out.root")'
 *  ACLiC compiled (recommended for large files):
 *    root 'make_event_waveforms.C+("in.root","out.root")'
 *
 * ── Parameters ───────────────────────────────────────────────────────────────
 *  input_path   : input ROOT file
 *  output_path  : output ROOT file
 *  tree_name    : TTree name                              (default "events")
 *  max_entries  : max TTree entries to read, -1 = all    (default -1)
 *  thr_mV       : signal threshold [mV], negative pulse  (default -10.0)
 *  thr_mcp_mV   : MCP-PMT threshold [mV] (EUDAQ_ID=3 C8) (default -200.0)
 *  ped_n        : samples used for pedestal estimate      (default 20)
 * ─────────────────────────────────────────────────────────────────────────────
 */

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TDirectory.h>

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  CONFIGURATION                                                              ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

/// Maximum waveform samples per channel.
const Int_t kMaxSamp = 1000;

/// Number of channels per EUDAQ scope.
const Int_t kNCh = 8;

/// EUDAQ_ID of the DUT scope (contains DUT strips C1..C7 and MCP-PMT C8).
const Long64_t kDUT_ID = 3LL;

/// Tracking plane definitions: EUDAQ_ID → output directory name.
struct TrkDef { Long64_t eid; const char* name; };
const TrkDef kTracking[] = {
    {4LL, "Tracking_1"},
    {1LL, "Tracking_2"},
    {2LL, "Tracking_4"},
    {0LL, "Tracking_5"},
};
const int kNTracking = 4;

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  BRANCH BUFFERS  (global: avoids stack overflow for large arrays)           ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

Long64_t g_EventN   = 0;
Long64_t g_EUDAQ_ID = 0;
Int_t    g_nvolt[kNCh + 1] = {};   // index 1..kNCh, index 0 unused
Int_t    g_ntime[kNCh + 1] = {};
Double_t g_volt[kNCh + 1][kMaxSamp];
Double_t g_time[kNCh + 1][kMaxSamp];

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  BRANCH SETUP                                                               ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

void SetupBranches(TTree* tree)
{
    tree->SetBranchStatus("*", false);

    tree->SetBranchStatus("EventN",   true);
    tree->SetBranchAddress("EventN",   &g_EventN);
    tree->SetBranchStatus("EUDAQ_ID", true);
    tree->SetBranchAddress("EUDAQ_ID", &g_EUDAQ_ID);

    for (int ch = 1; ch <= kNCh; ++ch) {
        char bnt[32], bnv[32], bnn_t[32], bnn_v[32];
        snprintf(bnt,   sizeof(bnt),   "C%d_time",  ch);
        snprintf(bnv,   sizeof(bnv),   "C%d_volt",  ch);
        snprintf(bnn_t, sizeof(bnn_t), "nC%d_time", ch);
        snprintf(bnn_v, sizeof(bnn_v), "nC%d_volt", ch);

        if (tree->GetBranch(bnn_t)) {
            tree->SetBranchStatus(bnn_t, true);
            tree->SetBranchAddress(bnn_t, &g_ntime[ch]);
        } else fprintf(stderr, "[WARN] Missing branch: %s\n", bnn_t);

        if (tree->GetBranch(bnn_v)) {
            tree->SetBranchStatus(bnn_v, true);
            tree->SetBranchAddress(bnn_v, &g_nvolt[ch]);
        } else fprintf(stderr, "[WARN] Missing branch: %s\n", bnn_v);

        if (tree->GetBranch(bnt)) {
            tree->SetBranchStatus(bnt, true);
            tree->SetBranchAddress(bnt, g_time[ch]);
        } else fprintf(stderr, "[WARN] Missing branch: %s\n", bnt);

        if (tree->GetBranch(bnv)) {
            tree->SetBranchStatus(bnv, true);
            tree->SetBranchAddress(bnv, g_volt[ch]);
        } else fprintf(stderr, "[WARN] Missing branch: %s\n", bnv);
    }
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  PER-EVENT DATA STRUCTURES                                                  ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

/// Raw waveform for one channel.
struct ChRaw {
    int n = 0;
    std::vector<double> time;
    std::vector<double> volt;
};

/// All channels from one TTree entry (one scope reading).
struct ScopeRaw {
    Long64_t eudaq_id = -1;
    std::map<int, ChRaw> ch;   // key = 1-based channel number
};

/// Pedestal-corrected result for one channel.
struct ChProc {
    double pedestal      = 0.0;
    double min_volt_corr = 0.0;   // most negative sample after pedestal subtraction
    std::vector<double> time;
    std::vector<double> volt_corr;
};

/// Copy current global branch buffers into a ScopeRaw.
ScopeRaw CaptureEntry()
{
    ScopeRaw sr;
    sr.eudaq_id = g_EUDAQ_ID;
    for (int ch = 1; ch <= kNCh; ++ch) {
        int n = std::min(g_ntime[ch], g_nvolt[ch]);
        if (n <= 0) continue;
        if (n > kMaxSamp) {
            fprintf(stderr, "[WARN] EUDAQ_ID=%lld C%d: n=%d > kMaxSamp=%d, truncating\n",
                    (long long)g_EUDAQ_ID, ch, n, (int)kMaxSamp);
            n = kMaxSamp;
        }
        ChRaw cr;
        cr.n = n;
        cr.time.assign(g_time[ch], g_time[ch] + n);
        cr.volt.assign(g_volt[ch], g_volt[ch] + n);
        sr.ch[ch] = std::move(cr);
    }
    return sr;
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  STATISTICS HELPERS                                                         ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

std::string FmtCP(Long64_t cnt, Long64_t total)
{
    double p   = (total > 0) ? (double)cnt / (double)total : 0.0;
    double pct = 100.0 * p;
    double err = (total > 0) ? 100.0 * std::sqrt(p * (1.0 - p) / (double)total) : 0.0;
    std::ostringstream o;
    o << std::right << std::setw(10) << (long long)cnt
      << "   (" << std::fixed << std::setprecision(1)
      << std::setw(5) << pct << " ± " << std::setw(5) << err << "%)";
    return o.str();
}

std::string EventDirName(Long64_t idx)
{
    char buf[32];
    snprintf(buf, sizeof(buf), "event_%06lld", (long long)idx);
    return buf;
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  STATISTICS STRUCT                                                          ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

struct Stats {
    Long64_t n_events             = 0;
    Long64_t count_MCP            = 0;   // MCP-PMT below thr_mcp
    Long64_t count_MCP_DUT        = 0;   // MCP-PMT and DUT
    Long64_t count_ANY            = 0;   // any detector below thr (saved)
    Long64_t count_trk[4]         = {};  // one per kTracking[] entry
    Long64_t count_Trk1_Trk2      = 0;   // Tracking_1 and Tracking_2
    Long64_t count_Trk4_Trk5      = 0;   // Tracking_4 and Tracking_5
    Long64_t count_Trk1_Trk5      = 0;   // Tracking_1 and Tracking_5
    Long64_t count_Trk2_Trk4      = 0;   // Tracking_2 and Tracking_4
    Long64_t count_Trk2_DUT_Trk4  = 0;   // Tracking_2 and DUT and Tracking_4
    Long64_t count_MCP_DUT_allTrk = 0;   // MCP-PMT and DUT and all 4 tracking planes
    Long64_t graphs_total         = 0;
    Long64_t saved_idx            = 0;
};

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  PROCESS ONE PHYSICS EVENT                                                  ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

void ProcessEvent(const std::vector<ScopeRaw>& scopes,
                  Long64_t evIdx,
                  TFile*   fout,
                  double   thr_V,
                  double   thr_mcp_V,
                  int      ped_n,
                  Stats&   st)
{
    // ── Step 1: pedestal subtraction ─────────────────────────────────────────
    // proc[eudaq_id][ch] = ChProc
    std::map<Long64_t, std::map<int, ChProc>> proc;

    for (size_t si = 0; si < scopes.size(); ++si) {
        const ScopeRaw& sr = scopes[si];
        for (auto cit = sr.ch.begin(); cit != sr.ch.end(); ++cit) {
            const ChRaw& raw = cit->second;
            if (raw.n <= 0) continue;

            int nped = std::min(ped_n, raw.n);
            double ped = 0.0;
            for (int j = 0; j < nped; ++j) ped += raw.volt[j];
            if (nped > 0) ped /= nped;

            ChProc cp;
            cp.pedestal = ped;
            cp.time = raw.time;
            cp.volt_corr.resize(raw.n);
            cp.min_volt_corr = std::numeric_limits<double>::max();
            for (int j = 0; j < raw.n; ++j) {
                cp.volt_corr[j] = raw.volt[j] - ped;
                if (cp.volt_corr[j] < cp.min_volt_corr)
                    cp.min_volt_corr = cp.volt_corr[j];
            }
            proc[sr.eudaq_id][cit->first] = std::move(cp);
        }
    }

    // ── Step 2a: MCP-PMT — EUDAQ_ID=3, C8 ───────────────────────────────────
    bool cond_MCP = false;
    {
        auto sit = proc.find(kDUT_ID);
        if (sit != proc.end()) {
            auto cit = sit->second.find(8);
            if (cit != sit->second.end() &&
                cit->second.min_volt_corr < thr_mcp_V)
                cond_MCP = true;
        }
    }
    if (cond_MCP) ++st.count_MCP;

    // ── Step 2b: build S — detector groups with any channel below thr_V ──────
    // DUT sensor: EUDAQ_ID=3, C1..C7
    // Tracking planes: EUDAQ_ID=0,1,2,4, all channels
    // MCP-PMT (EUDAQ_ID=3 C8) excluded from S.
    std::set<Long64_t> S;

    // DUT sensor
    {
        auto sit = proc.find(kDUT_ID);
        if (sit != proc.end()) {
            for (int ch = 1; ch <= 7; ++ch) {
                auto cit = sit->second.find(ch);
                if (cit != sit->second.end() &&
                    cit->second.min_volt_corr < thr_V) {
                    S.insert(kDUT_ID);
                    break;
                }
            }
        }
    }

    // Tracking planes
    for (int t = 0; t < kNTracking; ++t) {
        Long64_t eid = kTracking[t].eid;
        auto sit = proc.find(eid);
        if (sit == proc.end()) continue;
        for (auto cit = sit->second.begin(); cit != sit->second.end(); ++cit) {
            if (cit->second.min_volt_corr < thr_V) {
                S.insert(eid);
                break;
            }
        }
    }

    const bool keep_event = !S.empty();
    if (keep_event) ++st.count_ANY;
    if (cond_MCP && S.count(kDUT_ID)) ++st.count_MCP_DUT;

    // Per-tracking-plane counts and combined condition
    // kTracking[0]=Tracking_1, [1]=Tracking_2, [2]=Tracking_4, [3]=Tracking_5
    bool cond_allTrk = true;
    for (int t = 0; t < kNTracking; ++t) {
        bool hit = S.count(kTracking[t].eid) > 0;
        if (hit) ++st.count_trk[t];
        if (!hit) cond_allTrk = false;
    }

    // Pairwise and triple coincidences
    bool h1 = S.count(kTracking[0].eid) > 0;   // Tracking_1
    bool h2 = S.count(kTracking[1].eid) > 0;   // Tracking_2
    bool h4 = S.count(kTracking[2].eid) > 0;   // Tracking_4
    bool h5 = S.count(kTracking[3].eid) > 0;   // Tracking_5
    bool hD = S.count(kDUT_ID)          > 0;   // DUT sensor

    if (h1 && h2)      ++st.count_Trk1_Trk2;
    if (h4 && h5)      ++st.count_Trk4_Trk5;
    if (h1 && h5)      ++st.count_Trk1_Trk5;
    if (h2 && h4)      ++st.count_Trk2_Trk4;
    if (h2 && hD && h4) ++st.count_Trk2_DUT_Trk4;
    if (cond_MCP && hD && cond_allTrk) ++st.count_MCP_DUT_allTrk;

    if (!keep_event) return;

    // ── Step 3: write TGraphs ────────────────────────────────────────────────
    fout->cd();
    TDirectory* evdir = fout->mkdir(EventDirName(st.saved_idx).c_str());
    if (!evdir) {
        fprintf(stderr, "[ERROR] mkdir failed for %s\n", EventDirName(st.saved_idx).c_str());
        return;
    }
    ++st.saved_idx;

    // Helper: look up a ChProc by (eudaq_id, channel).
    auto getProc = [&](Long64_t eid, int ch) -> const ChProc* {
        auto sit = proc.find(eid);
        if (sit == proc.end()) return nullptr;
        auto cit = sit->second.find(ch);
        if (cit == sit->second.end()) return nullptr;
        return &cit->second;
    };

    // Helper: write one TGraph into dir.
    auto writeGraph = [&](TDirectory* dir, const char* gname, const char* gtitle,
                          const ChProc* cp) {
        if (!cp || cp->volt_corr.empty()) return;
        int n = (int)cp->volt_corr.size();
        TGraph* gr = new TGraph(n, cp->time.data(), cp->volt_corr.data());
        gr->SetName(gname);
        gr->SetTitle(gtitle);
        dir->cd();
        gr->Write();
        delete gr;
        ++st.graphs_total;
    };

    // ── 3a: DUT sensor  (EUDAQ_ID=3, C1..C7 → strip1..strip7) ───────────────
    evdir->cd();
    TDirectory* dutdir = evdir->mkdir("DUT sensor");
    if (dutdir) {
        for (int ch = 1; ch <= 7; ++ch) {
            char gname[32], gtitle[256];
            snprintf(gname, sizeof(gname), "strip%d", ch);
            snprintf(gtitle, sizeof(gtitle),
                     "evt=%lld DUT strip%d (EUDAQ_ID=3 C%d) thr=%.1fmV ped=%dpts",
                     (long long)evIdx, ch, ch, thr_V * 1e3, ped_n);
            writeGraph(dutdir, gname, gtitle, getProc(kDUT_ID, ch));
        }
    }

    // ── 3b: MCP-PMT  (EUDAQ_ID=3, C8 → wf) ──────────────────────────────────
    evdir->cd();
    TDirectory* mcpdir = evdir->mkdir("MCP-PMT");
    if (mcpdir) {
        char gtitle[256];
        snprintf(gtitle, sizeof(gtitle),
                 "evt=%lld MCP-PMT (EUDAQ_ID=3 C8) thr_mcp=%.1fmV ped=%dpts",
                 (long long)evIdx, thr_mcp_V * 1e3, ped_n);
        writeGraph(mcpdir, "wf", gtitle, getProc(kDUT_ID, 8));
    }

    // ── 3c: Tracking planes  (8 channels each → strip1..strip8) ─────────────
    for (int t = 0; t < kNTracking; ++t) {
        evdir->cd();
        TDirectory* trkdir = evdir->mkdir(kTracking[t].name);
        if (!trkdir) continue;
        for (int ch = 1; ch <= kNCh; ++ch) {
            char gname[32], gtitle[256];
            snprintf(gname, sizeof(gname), "strip%d", ch);
            snprintf(gtitle, sizeof(gtitle),
                     "evt=%lld %s strip%d (EUDAQ_ID=%lld C%d) thr=%.1fmV ped=%dpts",
                     (long long)evIdx, kTracking[t].name, ch,
                     (long long)kTracking[t].eid, ch, thr_V * 1e3, ped_n);
            writeGraph(trkdir, gname, gtitle, getProc(kTracking[t].eid, ch));
        }
    }
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  MACRO ENTRY POINT                                                          ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

void make_event_waveforms(const char* input_path,
                          const char* output_path,
                          const char* tree_name   = "events",
                          Long64_t    max_entries = -1,
                          Double_t    thr_mV      = -10.0,
                          Double_t    thr_mcp_mV  = -200.0,
                          Int_t       ped_n       = 20)
{
    const double thr_V     = thr_mV     * 1e-3;
    const double thr_mcp_V = thr_mcp_mV * 1e-3;

    printf("[INFO] ── Configuration ──────────────────────────────────────────\n");
    printf("[INFO]   Global threshold      : %.1f mV\n",    thr_mV);
    printf("[INFO]   MCP-PMT threshold     : %.1f mV\n",    thr_mcp_mV);
    printf("[INFO]   Pedestal window       : %d samples\n", (int)ped_n);

    TFile* fin = TFile::Open(input_path, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open input: " << input_path << "\n";
        return;
    }
    printf("[INFO]   Input  file           : %s\n", input_path);

    TTree* tree = nullptr;
    fin->GetObject(tree_name, tree);
    if (!tree) {
        fprintf(stderr, "[ERROR] TTree \"%s\" not found\n", tree_name);
        fin->Close();
        return;
    }

    Long64_t n_entries = tree->GetEntries();
    printf("[INFO]   TTree \"%s\"          : %lld entries\n",
           tree_name, (long long)n_entries);

    if (max_entries > 0 && max_entries < n_entries) {
        n_entries = max_entries;
        printf("[INFO]   Capped to             : %lld entries\n", (long long)n_entries);
    }
    SetupBranches(tree);

    TFile* fout = TFile::Open(output_path, "RECREATE");
    if (!fout || fout->IsZombie()) {
        fprintf(stderr, "[ERROR] Cannot create output: %s\n", output_path);
        fin->Close();
        return;
    }
    printf("[INFO]   Output file           : %s\n", output_path);
    printf("[INFO] ─────────────────────────────────────────────────────────\n");

    // ── Statistics ────────────────────────────────────────────────────────────
    Stats st;

    // ── Event loop with EventN-based grouping ─────────────────────────────────
    // Entries sharing the same EventN are buffered into one group, then
    // processed together (one group = one physics event).

    Long64_t currentEventN = std::numeric_limits<Long64_t>::min();
    Long64_t evIdx         = -1;    // sequential physics-event counter
    std::vector<ScopeRaw> group;

    auto flush_group = [&]() {
        if (group.empty()) return;
        ProcessEvent(group, evIdx, fout, thr_V, thr_mcp_V, ped_n, st);
        ++st.n_events;
        if (st.n_events % 500 == 0) fout->Flush();
        group.clear();
    };

    for (Long64_t i = 0; i < n_entries; ++i) {
        if (i == 0 || (i + 1) % 1000 == 0 || i + 1 == n_entries) {
            printf("[INFO] Entry %7lld / %7lld  "
                   "events: %lld  saved: %lld  graphs: %lld\r",
                   (long long)(i + 1), (long long)n_entries,
                   (long long)st.n_events, (long long)st.count_ANY,
                   (long long)st.graphs_total);
            fflush(stdout);
            if (i + 1 == n_entries) putchar('\n');
        }

        tree->GetEntry(i);

        if (g_EventN != currentEventN) {
            flush_group();
            currentEventN = g_EventN;
            ++evIdx;
        }
        group.push_back(CaptureEntry());
    }
    flush_group();

    fout->Close();
    fin->Close();

    // ── Summary statistics ────────────────────────────────────────────────────
    const int LW = 36;
    std::string bar(LW + 32, '=');
    printf("\n  %s\n", bar.c_str());
    printf("  %-*s\n", (int)bar.size() - 2, "  Summary Statistics");
    printf("  %s\n", bar.c_str());

    auto KvLine = [&](const std::string& label, const std::string& val) {
        printf("  %-*s : %s\n", LW, label.c_str(), val.c_str());
    };
    auto KvCount = [&](const std::string& label, Long64_t cnt) {
        std::ostringstream o; o << cnt;
        KvLine(label, o.str());
    };

    KvCount("Physics events (unique EventN)",    st.n_events);
    KvLine("MCP-PMT signal",                     FmtCP(st.count_MCP,             st.n_events));
    KvLine("MCP-PMT and DUT",                    FmtCP(st.count_MCP_DUT,         st.n_events));
    KvLine("Any detector crossed thr (saved)",   FmtCP(st.count_ANY,             st.n_events));
    for (int t = 0; t < kNTracking; ++t)
        KvLine(std::string(kTracking[t].name),   FmtCP(st.count_trk[t],          st.n_events));
    KvLine("Tracking_1 and Tracking_2",          FmtCP(st.count_Trk1_Trk2,       st.n_events));
    KvLine("Tracking_4 and Tracking_5",          FmtCP(st.count_Trk4_Trk5,       st.n_events));
    KvLine("Tracking_1 and Tracking_5",          FmtCP(st.count_Trk1_Trk5,       st.n_events));
    KvLine("Tracking_2 and Tracking_4",          FmtCP(st.count_Trk2_Trk4,       st.n_events));
    KvLine("Tracking_2 and DUT and Tracking_4",  FmtCP(st.count_Trk2_DUT_Trk4,   st.n_events));
    KvLine("MCP-PMT and DUT and all Tracking",   FmtCP(st.count_MCP_DUT_allTrk,  st.n_events));
    KvCount("TGraphs written",                   st.graphs_total);
    printf("  %s\n\n", bar.c_str());
}
