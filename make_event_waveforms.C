/**
 * make_event_waveforms.C
 *
 * Reads oscilloscope waveforms from a TTree, applies pedestal subtraction,
 * performs threshold-based event selection, and writes per-event TGraphs.
 *
 * TTree branch format (confirmed from Print() output):
 *   EventN/L, EUDAQ_ID/L
 *   nC#_time/I,  C#_time[nC#_time]/D     (# = 1..8)
 *   nC#_volt/I,  C#_volt[nC#_volt]/D
 *
 * Output structure:
 *   Directories are numbered sequentially for saved events only
 *   (event_000000 = 1st kept event, event_000001 = 2nd kept event, …).
 *   /event_000000/
 *     DUT sensor/   strip1..strip7  ← EUDAQ=5 C1..C7 (pedestal-subtracted)
 *     MCP-PMT/      wf              ← EUDAQ=5 C8
 *     Tracking_1/   strip1..strip8  ← EUDAQ=6 C1..C8
 *     Tracking_2/   strip1..strip8  ← EUDAQ=3
 *     Tracking_4/   strip1..strip8  ← EUDAQ=4
 *     Tracking_5/   strip1..strip8  ← EUDAQ=2
 *   /event_000001/…
 *
 * Statistics counting rule for S (set of detector groups crossing global threshold):
 *   DUT sensor  (key=0)  : any of DUT strip1..strip7 (EUDAQ=5 C1..C7) below thr_mV.
 *   Tracking_1  (key=10) : any of EUDAQ=6 channels below thr_mV.
 *   Tracking_2,4,5 (keys 3,4,2) : any channel in EUDAQ=3,4,2 below thr_mV.
 *   MCP-PMT (EUDAQ=5 C8) : excluded from S; tracked separately as count_MCP.
 *
 * ── Usage ────────────────────────────────────────────────────────────────────
 *  Interpreted:
 *    root 'make_event_waveforms.C("in.root","out.root")'
 *    root 'make_event_waveforms.C("in.root","out.root","events",-1,-10.,-200.,20)'
 *
 *  ACLiC compiled (recommended for large files):
 *    root 'make_event_waveforms.C+("in.root","out.root")'
 *
 * ── Parameters ───────────────────────────────────────────────────────────────
 *  input_path    : input ROOT file
 *  output_path   : output ROOT file
 *  tree_name     : TTree name                              (default "events")
 *  max_entries   : max TTree entries to read, -1 = all    (default -1)
 *  thr_mV        : global threshold [mV], negative pulse  (default -10.0)
 *  thr_mcp_mV    : special threshold for MCP-PMT (EUDAQ=5 C8) [mV](default -200.0)
 *  ped_n         : samples used for pedestal estimate      (default 20)
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
// ║  USER-EDITABLE CONFIGURATION                                                ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

/// Maximum waveform samples per channel (must be >= largest nC#_volt in file).
const Int_t kMaxSamp = 10000;

/// Channels to read per EUDAQ_ID.
const std::map<Long64_t, int> kChannelsPerScope = {
    {2, 8}, {3, 8}, {4, 8}, {5, 8}, {6, 8},
};

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  DUT SENSOR CHANNEL MAPPING                                                 ║
// ║                                                                             ║
// ║  DUT sensor comes from EUDAQ=5.  C1..C7 are the seven DUT strip channels;  ║
// ║  C8 of EUDAQ=5 is the MCP-PMT signal and is written separately.            ║
// ║                                                                             ║
// ║  Output DUT strip  →  source (EUDAQ_ID, channel):                          ║
// ║    DUT strip1 ← EUDAQ=5 C1    DUT strip5 ← EUDAQ=5 C5                     ║
// ║    DUT strip2 ← EUDAQ=5 C2    DUT strip6 ← EUDAQ=5 C6                     ║
// ║    DUT strip3 ← EUDAQ=5 C3    DUT strip7 ← EUDAQ=5 C7                     ║
// ║    DUT strip4 ← EUDAQ=5 C4                                                 ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

struct DutSrc { Long64_t eid; int src_ch; };

// Index 0 is unused; valid DUT strip numbers are 1..kDutNCh.
const DutSrc kDutChannels[8] = {
    {-1, -1},   // [0] unused
    { 5,  1},   // DUT strip1 ← EUDAQ=5 C1
    { 5,  2},   // DUT strip2 ← EUDAQ=5 C2
    { 5,  3},   // DUT strip3 ← EUDAQ=5 C3
    { 5,  4},   // DUT strip4 ← EUDAQ=5 C4
    { 5,  5},   // DUT strip5 ← EUDAQ=5 C5
    { 5,  6},   // DUT strip6 ← EUDAQ=5 C6
    { 5,  7},   // DUT strip7 ← EUDAQ=5 C7
};
const int kDutNCh = 7;   // DUT strip channels: 1..kDutNCh

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  TRACKING_1 CHANNEL MAPPING                                                 ║
// ║                                                                             ║
// ║  Tracking_1 comes from EUDAQ=6 (8 channels).                               ║
// ║                                                                             ║
// ║  Output Tracking_1 strip  →  source (EUDAQ_ID, channel):                   ║
// ║    strip1 ← EUDAQ=6 C1    strip5 ← EUDAQ=6 C5                             ║
// ║    strip2 ← EUDAQ=6 C2    strip6 ← EUDAQ=6 C6                             ║
// ║    strip3 ← EUDAQ=6 C3    strip7 ← EUDAQ=6 C7                             ║
// ║    strip4 ← EUDAQ=6 C4    strip8 ← EUDAQ=6 C8                             ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

// Index 0 is unused; valid Tracking_1 strip numbers are 1..kTrackNCh.
const DutSrc kTrackChannels[9] = {
    {-1, -1},   // [0] unused
    { 6,  1},   // Tracking strip1 ← EUDAQ=6 C1
    { 6,  2},   // Tracking strip2 ← EUDAQ=6 C2
    { 6,  3},   // Tracking strip3 ← EUDAQ=6 C3
    { 6,  4},   // Tracking strip4 ← EUDAQ=6 C4
    { 6,  5},   // Tracking strip5 ← EUDAQ=6 C5
    { 6,  6},   // Tracking strip6 ← EUDAQ=6 C6
    { 6,  7},   // Tracking strip7 ← EUDAQ=6 C7
    { 6,  8},   // Tracking strip8 ← EUDAQ=6 C8
};
const int kTrackNCh = 8;   // Tracking_1 strip channels: 1..kTrackNCh

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  BRANCH BUFFERS  (global: avoids stack overflow from large arrays)          ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

Long64_t g_EventN   = 0;
Long64_t g_EUDAQ_ID = 0;
Int_t    g_ntime[9] = {};   // index 1..8, index 0 unused
Int_t    g_nvolt[9] = {};
Double_t g_time[9][kMaxSamp];
Double_t g_volt[9][kMaxSamp];

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  BRANCH SETUP                                                               ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

void SetupBranches(TTree* tree)
{
    tree->SetBranchStatus("*", false);

    // metadata
    tree->SetBranchStatus("EventN",   true);
    tree->SetBranchStatus("EUDAQ_ID", true);
    tree->SetBranchAddress("EventN",   &g_EventN);
    tree->SetBranchAddress("EUDAQ_ID", &g_EUDAQ_ID);

    for (int ch = 1; ch <= 8; ++ch) {
        char bnt[32], bnv[32], bnn_t[32], bnn_v[32];
        snprintf(bnt,   sizeof(bnt),   "C%d_time",  ch);
        snprintf(bnv,   sizeof(bnv),   "C%d_volt",  ch);
        snprintf(bnn_t, sizeof(bnn_t), "nC%d_time", ch);
        snprintf(bnn_v, sizeof(bnn_v), "nC%d_volt", ch);

        // count branches  (nC#_time / nC#_volt)
        if (tree->GetBranch(bnn_t)) {
            tree->SetBranchStatus(bnn_t, true);
            tree->SetBranchAddress(bnn_t, &g_ntime[ch]);
        } else std::cerr << "[WARN] Missing: " << bnn_t << "\n";

        if (tree->GetBranch(bnn_v)) {
            tree->SetBranchStatus(bnn_v, true);
            tree->SetBranchAddress(bnn_v, &g_nvolt[ch]);
        } else std::cerr << "[WARN] Missing: " << bnn_v << "\n";

        // array branches  (C#_time / C#_volt)
        if (tree->GetBranch(bnt)) {
            tree->SetBranchStatus(bnt, true);
            tree->SetBranchAddress(bnt, g_time[ch]);
        } else std::cerr << "[WARN] Missing: " << bnt << "\n";

        if (tree->GetBranch(bnv)) {
            tree->SetBranchStatus(bnv, true);
            tree->SetBranchAddress(bnv, g_volt[ch]);
        } else std::cerr << "[WARN] Missing: " << bnv << "\n";
    }
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  PER-EVENT DATA STRUCTURES                                                  ║
// ║                                                                             ║
// ║  Each TTree entry belongs to one EUDAQ_ID.  Several entries share the same  ║
// ║  EventN (one per scope).  We buffer a full EventN group before deciding     ║
// ║  whether to keep the event and what to write.                               ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

/// Raw waveform for one channel (copied from global branch buffers).
struct ChRaw {
    int n = 0;
    std::vector<double> time;
    std::vector<double> volt;
};

/// All channels contributed by one TTree entry (one scope reading).
struct ScopeRaw {
    Long64_t eudaq_id = -1;
    std::map<int, ChRaw> ch;   // key = channel number (1-based)
};

/// Pedestal-corrected result for one channel (computed from ChRaw).
struct ChProc {
    double pedestal      = 0.0;
    double min_volt_corr = 0.0;
    std::vector<double> time;
    std::vector<double> volt_corr;
};

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  HELPERS                                                                    ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

int GetNch(Long64_t eid) {
    auto it = kChannelsPerScope.find(eid);
    return it != kChannelsPerScope.end() ? it->second : 8;
}

std::string EventDirName(Long64_t evIdx) {
    char buf[32];
    snprintf(buf, sizeof(buf), "event_%06lld", (long long)evIdx);
    return buf;
}

std::string ScopeDirName(Long64_t eid) {
    char buf[32];
    snprintf(buf, sizeof(buf), "EUDAQ_%lld", (long long)eid);
    return buf;
}

// Output directory name for EUDAQ_2,3,4 (renamed to Tracking planes).
// EUDAQ=2 → Tracking_5,  EUDAQ=3 → Tracking_2,  EUDAQ=4 → Tracking_4
std::string TrackingDirName(Long64_t eid) {
    if (eid == 2) return "Tracking_5";
    if (eid == 3) return "Tracking_2";
    if (eid == 4) return "Tracking_4";
    return ScopeDirName(eid);   // fallback (should not be reached)
}

/// Copy the current global branch buffers into a ScopeRaw struct.
ScopeRaw CaptureEntry()
{
    ScopeRaw sr;
    sr.eudaq_id = g_EUDAQ_ID;
    const int nch = GetNch(g_EUDAQ_ID);
    for (int ch = 1; ch <= nch; ++ch) {
        int n = (int)std::min(g_ntime[ch], g_nvolt[ch]);
        if (n <= 0) continue;
        if (n > (int)kMaxSamp) {
            fprintf(stderr,
                "[WARN] EUDAQ=%lld C%d: n=%d > kMaxSamp=%d, truncating.\n",
                (long long)g_EUDAQ_ID, ch, n, (int)kMaxSamp);
            n = (int)kMaxSamp;
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
// ║  PROCESS ONE PHYSICS EVENT                                                  ║
// ║                                                                             ║
// ║  1. Compute pedestal-subtracted waveforms for all scope/channel pairs.      ║
// ║  2. Evaluate threshold conditions and accumulate statistics.                ║
// ║  3. If keep_event: create directories and write TGraphs.                   ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

void ProcessEvent(const std::vector<ScopeRaw>& scopes,
                  Long64_t    evIdx,          // original sequential event index (for titles)
                  TFile*      fout,
                  double      thr_V,          // global threshold [V]
                  double      thr_mcp_V,      // MCP-PMT threshold [V]  (EUDAQ=5 C8)
                  int         ped_n,          // pedestal window [samples]
                  Long64_t&   count_MCP,      // MCP-PMT crosses its threshold
                  Long64_t&   count_MCP_DUT,  // MCP-PMT AND DUT sensor both cross
                  Long64_t&   count_ANY,      // any detector group crosses global thr
                  std::map<int,Long64_t>& count_k,   // distribution of |S|, k=0..5
                  Long64_t&   graphs_written,
                  Long64_t&   saved_idx,      // sequential index of saved events (directory name)
                  // ── Coincidence monitors: MCP & DUT & tracking plane(s) ──────
                  Long64_t&   cnt_MDT1,       // MCP & DUT & Trk1
                  Long64_t&   cnt_MDT2,       // MCP & DUT & Trk2
                  Long64_t&   cnt_MDT4,       // MCP & DUT & Trk4
                  Long64_t&   cnt_MDT5,       // MCP & DUT & Trk5
                  Long64_t&   cnt_MDT12,      // MCP & DUT & Trk1 & Trk2
                  Long64_t&   cnt_MDT24,      // MCP & DUT & Trk2 & Trk4
                  Long64_t&   cnt_MDT45,      // MCP & DUT & Trk4 & Trk5
                  Long64_t&   cnt_MDT1245)    // MCP & DUT & Trk1 & Trk2 & Trk4 & Trk5
{
    // ── Step 1: pedestal subtraction ─────────────────────────────────────────
    // proc[eudaq_id][ch] = ChProc
    std::map<Long64_t, std::map<int, ChProc>> proc;

    for (size_t si = 0; si < scopes.size(); ++si) {
        const ScopeRaw& sr = scopes[si];
        for (auto cit = sr.ch.begin(); cit != sr.ch.end(); ++cit) {
            int ch = cit->first;
            const ChRaw& raw = cit->second;
            if (raw.n <= 0) continue;

            // pedestal = mean of first ped_n samples
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

            proc[sr.eudaq_id][ch] = std::move(cp);
        }
    }

    // ── Step 2a: MCP-PMT  – EUDAQ=5 C8 below its own threshold ─────────────
    // cond_MCP is kept in function scope so it can be combined with S below.
    // Does NOT affect keep_event.
    bool cond_MCP = false;
    {
        auto sit = proc.find(5LL);
        if (sit != proc.end()) {
            auto cit = sit->second.find(8);
            if (cit != sit->second.end() &&
                cit->second.min_volt_corr < thr_mcp_V)
                cond_MCP = true;
        }
    }
    if (cond_MCP) ++count_MCP;

    // ── Step 2b: build S  – logical detector groups with any channel below thr_V
    //
    // Groups and their integer keys used in S:
    //    0 = "DUT sensor"  (kDutChannels[1..kDutNCh]: EUDAQ=5 C1..C7)
    //   10 = "Tracking_1"  (kTrackChannels[1..kTrackNCh]: EUDAQ=6)
    //    2 = "Tracking_5"  (proc[2] = EUDAQ=2, all channels)
    //    3 = "Tracking_2"  (proc[3] = EUDAQ=3, all channels)
    //    4 = "Tracking_4"  (proc[4] = EUDAQ=4, all channels)
    // MCP-PMT (EUDAQ=5 C8) is excluded from S entirely.
    // |S| is in [0..5], matching the five logical detector groups.
    std::set<int> S;

    // DUT sensor: EUDAQ=5 C1..C7 (via kDutChannels)
    for (int dut_ch = 1; dut_ch <= kDutNCh; ++dut_ch) {
        Long64_t src_eid = kDutChannels[dut_ch].eid;
        int      src_ch  = kDutChannels[dut_ch].src_ch;
        auto sit = proc.find(src_eid);
        if (sit == proc.end()) continue;
        auto cit = sit->second.find(src_ch);
        if (cit == sit->second.end()) continue;
        if (cit->second.min_volt_corr < thr_V) {
            S.insert(0);    // DUT sensor → key 0
            break;
        }
    }

    // Tracking_1: EUDAQ=6 (via kTrackChannels)
    for (int tr_ch = 1; tr_ch <= kTrackNCh; ++tr_ch) {
        Long64_t src_eid = kTrackChannels[tr_ch].eid;
        int      src_ch  = kTrackChannels[tr_ch].src_ch;
        auto sit = proc.find(src_eid);
        if (sit == proc.end()) continue;
        auto cit = sit->second.find(src_ch);
        if (cit == sit->second.end()) continue;
        if (cit->second.min_volt_corr < thr_V) {
            S.insert(10);   // Tracking_1 → key 10
            break;
        }
    }

    // Tracking_2,4,5 (EUDAQ=3,4,2): any channel crossing threshold
    for (auto sit = proc.begin(); sit != proc.end(); ++sit) {
        Long64_t eid = sit->first;
        if (eid == 5 || eid == 6) continue;   // EUDAQ=5 → DUT+MCP; EUDAQ=6 → Tracking_1
        for (auto cit = sit->second.begin(); cit != sit->second.end(); ++cit) {
            if (cit->second.min_volt_corr < thr_V) {
                S.insert((int)eid);
                break;
            }
        }
    }

    const bool keep_event = !S.empty();
    if (keep_event) ++count_ANY;

    // Convenience flags for coincidence evaluation.
    const bool has_DUT  = (S.count(0)  > 0);
    const bool has_Trk1 = (S.count(10) > 0);  // Tracking_1 (EUDAQ=6)
    const bool has_Trk2 = (S.count(3)  > 0);  // Tracking_2 (EUDAQ=3)
    const bool has_Trk4 = (S.count(4)  > 0);  // Tracking_4 (EUDAQ=4)
    const bool has_Trk5 = (S.count(2)  > 0);  // Tracking_5 (EUDAQ=2)

    // Basic MCP+DUT coincidence.
    if (cond_MCP && has_DUT) ++count_MCP_DUT;

    // MCP & DUT & single tracking plane.
    if (cond_MCP && has_DUT && has_Trk1) ++cnt_MDT1;
    if (cond_MCP && has_DUT && has_Trk2) ++cnt_MDT2;
    if (cond_MCP && has_DUT && has_Trk4) ++cnt_MDT4;
    if (cond_MCP && has_DUT && has_Trk5) ++cnt_MDT5;

    // MCP & DUT & adjacent tracking plane pairs.
    if (cond_MCP && has_DUT && has_Trk1 && has_Trk2) ++cnt_MDT12;
    if (cond_MCP && has_DUT && has_Trk2 && has_Trk4) ++cnt_MDT24;
    if (cond_MCP && has_DUT && has_Trk4 && has_Trk5) ++cnt_MDT45;

    // MCP & DUT & all four tracking planes.
    if (cond_MCP && has_DUT && has_Trk1 && has_Trk2 && has_Trk4 && has_Trk5)
        ++cnt_MDT1245;

    count_k[(int)S.size()]++;

    if (!keep_event) return;   // ← discard: no directories, no TGraphs

    // ── Step 3: write TGraphs ────────────────────────────────────────────────
    fout->cd();
    TDirectory* evdir = fout->mkdir(EventDirName(saved_idx).c_str());
    if (!evdir) {
        fprintf(stderr, "[ERROR] mkdir failed for %s – skipping.\n",
                EventDirName(saved_idx).c_str());
        return;
    }
    ++saved_idx;

    // Helper lambda: look up a ChProc by original (eudaq_id, channel).
    // Returns nullptr if not found.
    auto getProc = [&](Long64_t eid, int ch) -> const ChProc* {
        auto sit = proc.find(eid);
        if (sit == proc.end()) return nullptr;
        auto cit = sit->second.find(ch);
        if (cit == sit->second.end()) return nullptr;
        return &cit->second;
    };

    // ── 3a: "DUT sensor"  (EUDAQ=5 C1..C7 → strip1..strip7) ─────────────────
    evdir->cd();
    TDirectory* dutdir = evdir->mkdir("DUT sensor");
    if (dutdir) {
        dutdir->cd();
        for (int dut_ch = 1; dut_ch <= kDutNCh; ++dut_ch) {
            Long64_t src_eid = kDutChannels[dut_ch].eid;
            int      src_ch  = kDutChannels[dut_ch].src_ch;

            const ChProc* cp = getProc(src_eid, src_ch);
            if (!cp) continue;
            int n = (int)cp->volt_corr.size();
            if (n <= 0) continue;

            char gname[32], gtitle[256];
            snprintf(gname, sizeof(gname), "strip%d", dut_ch);
            snprintf(gtitle, sizeof(gtitle),
                     "event=%lld DUT_strip%d (src EUDAQ=%lld C%d) "
                     "thr=%.1fmV ped=%dpts",
                     (long long)evIdx, dut_ch,
                     (long long)src_eid, src_ch,
                     thr_V * 1e3, ped_n);

            TGraph* gr = new TGraph(n, cp->time.data(), cp->volt_corr.data());
            gr->SetName(gname);
            gr->SetTitle(gtitle);
            dutdir->cd();
            gr->Write();
            delete gr;
            ++graphs_written;
        }
    }

    // ── 3b: "MCP-PMT"  (EUDAQ=5 C8) ─────────────────────────────────────────
    evdir->cd();
    TDirectory* mcpdir = evdir->mkdir("MCP-PMT");
    if (mcpdir) {
        const ChProc* cp = getProc(5LL, 8);
        if (cp && (int)cp->volt_corr.size() > 0) {
            int n = (int)cp->volt_corr.size();
            char gtitle[256];
            snprintf(gtitle, sizeof(gtitle),
                     "event=%lld MCP-PMT (src EUDAQ=5 C8) "
                     "thr_mcp=%.1fmV ped=%dpts",
                     (long long)evIdx, thr_mcp_V * 1e3, ped_n);

            TGraph* gr = new TGraph(n, cp->time.data(), cp->volt_corr.data());
            gr->SetName("wf");
            gr->SetTitle(gtitle);
            mcpdir->cd();
            gr->Write();
            delete gr;
            ++graphs_written;
        }
    }

    // ── 3c: "Tracking_1"  (EUDAQ=6 C1..C8 → strip1..strip8) ────────────────
    evdir->cd();
    TDirectory* trkdir = evdir->mkdir("Tracking_1");
    if (trkdir) {
        trkdir->cd();
        for (int tr_ch = 1; tr_ch <= kTrackNCh; ++tr_ch) {
            Long64_t src_eid = kTrackChannels[tr_ch].eid;
            int      src_ch  = kTrackChannels[tr_ch].src_ch;

            const ChProc* cp = getProc(src_eid, src_ch);
            if (!cp) continue;
            int n = (int)cp->volt_corr.size();
            if (n <= 0) continue;

            char gname[32], gtitle[256];
            snprintf(gname, sizeof(gname), "strip%d", tr_ch);
            snprintf(gtitle, sizeof(gtitle),
                     "event=%lld Tracking_1 strip%d (src EUDAQ=%lld C%d) "
                     "thr=%.1fmV ped=%dpts",
                     (long long)evIdx, tr_ch,
                     (long long)src_eid, src_ch,
                     thr_V * 1e3, ped_n);

            TGraph* gr = new TGraph(n, cp->time.data(), cp->volt_corr.data());
            gr->SetName(gname);
            gr->SetTitle(gtitle);
            trkdir->cd();
            gr->Write();
            delete gr;
            ++graphs_written;
        }
    }

    // ── 3d: Tracking_2,4,5  (EUDAQ=2,3,4 renamed; 8 channels each) ──────────
    for (auto sit = proc.begin(); sit != proc.end(); ++sit) {
        Long64_t eid = sit->first;
        if (eid == 5 || eid == 6) continue;   // EUDAQ=5 → DUT+MCP; EUDAQ=6 → Tracking_1

        const std::string dirName = TrackingDirName(eid);
        evdir->cd();
        TDirectory* scopedir = evdir->mkdir(dirName.c_str());
        if (!scopedir) continue;
        scopedir->cd();

        for (auto cit = sit->second.begin(); cit != sit->second.end(); ++cit) {
            int ch = cit->first;
            const ChProc& cp = cit->second;
            int n = (int)cp.volt_corr.size();
            if (n <= 0) continue;

            char gname[32], gtitle[256];
            snprintf(gname, sizeof(gname), "strip%d", ch);
            snprintf(gtitle, sizeof(gtitle),
                     "event=%lld %s strip%d (src EUDAQ=%lld) thr=%.1fmV ped=%dpts",
                     (long long)evIdx, dirName.c_str(), ch,
                     (long long)eid, thr_V * 1e3, ped_n);

            TGraph* gr = new TGraph(n, cp.time.data(), cp.volt_corr.data());
            gr->SetName(gname);
            gr->SetTitle(gtitle);
            scopedir->cd();
            gr->Write();
            delete gr;
            ++graphs_written;
        }
    }
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  STATISTICS BOX HELPERS                                                     ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

/// Repeat a UTF-8 multi-byte character n times (e.g. RepeatChar("═", 40)).
std::string RepeatChar(const char* utf8ch, int n) {
    std::string s;
    for (int i = 0; i < n; ++i) s += utf8ch;
    return s;
}

/// Format count + percentage with binomial uncertainty (normal approximation).
/// Right-align count in 10 cols, pct/err in 7 cols each.
/// Example: "      1234   ( 12.345 ±  0.678%)"
std::string FmtCP(Long64_t cnt, Long64_t total) {
    // Binomial proportion (normal approximation) uncertainty:
    //   p = cnt/total,  sigma_p = sqrt(p(1-p)/total)
    // Display in percent points: 100*sigma_p
    double p   = (total > 0) ? (double)cnt / (double)total : 0.0;
    double pct = 100.0 * p;
    double err = (total > 0) ? 100.0 * std::sqrt(p * (1.0 - p) / (double)total) : 0.0;

    std::ostringstream o;
    o << std::right << std::setw(10) << (long long)cnt
      << "   (" << std::fixed << std::setprecision(1)
      << std::setw(5) << pct
      << " ± " << std::setw(5) << err
      << "%)";
    return o.str();
}

/// Left-align label in lw chars, then ": ", then val.
std::string KvLine(const std::string& label, const std::string& val, int lw) {
    std::ostringstream o;
    o << std::left << std::setw(lw) << label << ": " << val;
    return o.str();
}

/// Three-column row: c1 left-padded to w1, c2 right-padded to w2, "  ", c3 left-padded to w3.
std::string Row3(const std::string& c1, const std::string& c2,
                 const std::string& c3, int w1, int w2, int w3) {
    std::ostringstream o;
    o << std::left  << std::setw(w1) << c1
      << std::right << std::setw(w2) << c2
      << "  "
      << std::left  << std::setw(w3) << c3;
    return o.str();
}

/// Print a dynamic-width Unicode box around the content lines.
/// Lines equal to "__HR__" are rendered as ╠═══╣ separator rows.
void PrintBoxed(const std::vector<std::string>& lines, const std::string& title) {
    // Compute inner width = max(content line lengths + 2 padding, title + 4)
    int inner = (int)title.size() + 4;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i] == "__HR__") continue;
        int w = (int)lines[i].size() + 2;   // 1 space on each side
        if (w > inner) inner = w;
    }

    const std::string bar = RepeatChar("═", inner);

    // Top border + centred title
    printf("╔%s╗\n", bar.c_str());
    int tl = (int)title.size();
    int pl = (inner - tl) / 2;
    int pr = inner - tl - pl;
    printf("║%*s%s%*s║\n", pl, "", title.c_str(), pr, "");

    // Content lines
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i] == "__HR__")
            printf("╠%s╣\n", bar.c_str());
        else
            printf("║ %-*s║\n", inner - 1, lines[i].c_str());
    }

    printf("╚%s╝\n", bar.c_str());
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
    const double thr_V     = thr_mV    * 1e-3;
    const double thr_mcp_V = thr_mcp_mV * 1e-3;

    printf("[INFO] ── Configuration ─────────────────────────────────\n");
    printf("[INFO]   Global threshold        : %.1f mV\n",  thr_mV);
    printf("[INFO]   MCP-PMT threshold       : %.1f mV\n",  thr_mcp_mV);
    printf("[INFO]   Pedestal window         : %d samples\n", (int)ped_n);

    // ── Open input ────────────────────────────────────────────────────────────
    TFile* fin = TFile::Open(input_path, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[ERROR] Cannot open input: " << input_path << "\n";
        return;
    }
    printf("[INFO]   Input  file             : %s\n", input_path);

    TTree* tree = nullptr;
    fin->GetObject(tree_name, tree);
    if (!tree) {
        fprintf(stderr, "[ERROR] TTree \"%s\" not found.\n", tree_name);
        fin->Close();
        return;
    }

    Long64_t n_entries = tree->GetEntries();
    printf("[INFO]   TTree \"%s\"            : %lld entries\n",
           tree_name, (long long)n_entries);

    if (max_entries > 0 && max_entries < n_entries) {
        n_entries = max_entries;
        printf("[INFO]   Capped to               : %lld entries\n",
               (long long)n_entries);
    }
    
    SetupBranches(tree);

    // ── Open output ───────────────────────────────────────────────────────────
    TFile* fout = TFile::Open(output_path, "RECREATE");
    if (!fout || fout->IsZombie()) {
        fprintf(stderr, "[ERROR] Cannot create output: %s\n", output_path);
        fin->Close();
        return;
    }
    printf("[INFO]   Output file             : %s\n", output_path);
    printf("[INFO] ─────────────────────────────────────────────────\n");

    // ── Statistics ────────────────────────────────────────────────────────────
    Long64_t count_MCP     = 0;  // MCP-PMT crosses its threshold
    Long64_t count_MCP_DUT = 0;  // MCP-PMT AND DUT sensor both cross
    Long64_t count_ANY     = 0;  // any detector group crosses global thr
    Long64_t graphs_total  = 0;
    Long64_t n_events      = 0;  // physics events processed (unique EventN groups)
    // count_k[k] = events where exactly k logical detector groups crossed thr_V
    // Groups: {DUT sensor, Tracking_1, Tracking_2, Tracking_4, Tracking_5} → k in [0..5]
    std::map<int, Long64_t> count_k;
    for (int k = 0; k <= 5; ++k) count_k[k] = 0;

    // ── Coincidence monitors ──────────────────────────────────────────────────
    Long64_t cnt_MDT1    = 0;  // MCP & DUT & Trk1
    Long64_t cnt_MDT2    = 0;  // MCP & DUT & Trk2
    Long64_t cnt_MDT4    = 0;  // MCP & DUT & Trk4
    Long64_t cnt_MDT5    = 0;  // MCP & DUT & Trk5
    Long64_t cnt_MDT12   = 0;  // MCP & DUT & Trk1 & Trk2
    Long64_t cnt_MDT24   = 0;  // MCP & DUT & Trk2 & Trk4
    Long64_t cnt_MDT45   = 0;  // MCP & DUT & Trk4 & Trk5
    Long64_t cnt_MDT1245 = 0;  // MCP & DUT & Trk1 & Trk2 & Trk4 & Trk5

    // ── Event loop with EventN-based grouping ─────────────────────────────────
    // Entries with the same EventN are collected into one group and processed
    // together.  This assumes entries are sorted by EventN (standard for EUDAQ).
    //
    // evIdx is a sequential 0-based physics-event counter used for directory
    // naming.  It advances for every unique EventN, regardless of whether the
    // event is kept, so event_000042 always refers to the 42nd physics event
    // in the file, not the 42nd *kept* event.

    Long64_t currentEventN = std::numeric_limits<Long64_t>::min();
    Long64_t evIdx         = -1;   // sequential index of all physics events
    Long64_t saved_idx     = 0;    // sequential index of saved (kept) events
    std::vector<ScopeRaw> group;

    auto flush_group = [&]() {
        if (group.empty()) return;
        ProcessEvent(group, evIdx, fout,
                     thr_V, thr_mcp_V, ped_n,
                     count_MCP, count_MCP_DUT, count_ANY, count_k, graphs_total,
                     saved_idx,
                     cnt_MDT1, cnt_MDT2, cnt_MDT4, cnt_MDT5,
                     cnt_MDT12, cnt_MDT24, cnt_MDT45, cnt_MDT1245);
        ++n_events;
        if (n_events % 500 == 0) fout->Flush();
        group.clear();
    };

    for (Long64_t i = 0; i < n_entries; ++i) {

        // Progress line (overwrites itself)
        if (i == 0 || (i + 1) % 1000 == 0 || i + 1 == n_entries) {
            printf("[INFO] Entry %7lld / %7lld  "
                   "events: %lld  kept: %lld  graphs: %lld\r",
                   (long long)(i + 1), (long long)n_entries,
                   (long long)n_events, (long long)count_ANY,
                   (long long)graphs_total);
            fflush(stdout);
            if (i + 1 == n_entries) putchar('\n');
        }

        tree->GetEntry(i);

        // EventN changed → flush the accumulated group and start a new one
        if (g_EventN != currentEventN) {
            flush_group();
            currentEventN = g_EventN;
            ++evIdx;
        }

        group.push_back(CaptureEntry());
    }

    flush_group();   // handle the final group

    // ── Finalise ──────────────────────────────────────────────────────────────
    fout->Close();
    fin->Close();

    // ── Print statistics ──────────────────────────────────────────────────────
    {
        const Long64_t total = n_events;    // denominator: physics events (unique EventN)
        const int LW = 36;                 // label column width for KvLine
        const int DW = 22, CW = 12, PW = 12;  // distribution table column widths

        std::vector<std::string> lines;

        // ── Overview ──────────────────────────────────────────────────────────
        {
            std::ostringstream o;
            o << std::right << std::setw(10) << (long long)n_entries;
            lines.push_back(KvLine("TTree entries scanned", o.str(), LW));
        }
        {
            std::ostringstream o;
            o << std::right << std::setw(10) << (long long)n_events;
            lines.push_back(KvLine("Physics events (unique EventN)", o.str(), LW));
        }

        lines.push_back("__HR__");

        // ── Thresholds ────────────────────────────────────────────────────────
        {
            std::ostringstream o;
            o << std::fixed << std::setprecision(1) << thr_mV << " mV";
            lines.push_back(KvLine("Global thr (DUT + Tracking_1,2,4,5)", o.str(), LW));
        }
        {
            std::ostringstream o;
            o << std::fixed << std::setprecision(1) << thr_mcp_mV << " mV";
            lines.push_back(KvLine("MCP-PMT thr (EUDAQ=5 C8)", o.str(), LW));
        }
        {
            std::ostringstream o;
            o << (int)ped_n << " samples";
            lines.push_back(KvLine("Pedestal window", o.str(), LW));
        }

        lines.push_back("__HR__");

        // ── Counters ──────────────────────────────────────────────────────────
        lines.push_back(KvLine("At least MCP-PMT",         FmtCP(count_MCP,     total), LW));
        lines.push_back(KvLine("MCP-PMT and DUT",          FmtCP(count_MCP_DUT, total), LW));
        lines.push_back(KvLine("Any detector (Count_ANY)", FmtCP(count_ANY,     total), LW));
        {
            std::ostringstream o;
            o << std::right << std::setw(10) << (long long)graphs_total;
            lines.push_back(KvLine("TGraphs written", o.str(), LW));
        }

        lines.push_back("__HR__");

        // ── Coincidence monitors ──────────────────────────────────────────────
        lines.push_back("Coincidence monitors  (MCP & DUT & ...)");
        lines.push_back(KvLine("  + Trk1",                        FmtCP(cnt_MDT1,    total), LW));
        lines.push_back(KvLine("  + Trk2",                        FmtCP(cnt_MDT2,    total), LW));
        lines.push_back(KvLine("  + Trk4",                        FmtCP(cnt_MDT4,    total), LW));
        lines.push_back(KvLine("  + Trk5",                        FmtCP(cnt_MDT5,    total), LW));
        lines.push_back(KvLine("  + Trk1 & Trk2",                 FmtCP(cnt_MDT12,   total), LW));
        lines.push_back(KvLine("  + Trk2 & Trk4",                 FmtCP(cnt_MDT24,   total), LW));
        lines.push_back(KvLine("  + Trk4 & Trk5",                 FmtCP(cnt_MDT45,   total), LW));
        lines.push_back(KvLine("  + Trk1 & Trk2 & Trk4 & Trk5",  FmtCP(cnt_MDT1245, total), LW));


        printf("\n");
        PrintBoxed(lines, "Summary Statistics");
    }
}
