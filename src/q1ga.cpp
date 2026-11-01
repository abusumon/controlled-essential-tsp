#include <bits/stdc++.h>
using namespace std;
namespace fs = std::filesystem;

string DIR = "inputs", OUT = "output", INIT = "random", KSPEC = "10", XOVER = "cx", SEL = "roulette", PARTNER = "uniform";
int POP = 300, ELITE = 60, GENS = 2000, GK = 0, NRUNS = 30, TK = 3, KNN = 10, FILEORDER = 0;
double PM = 0.005, ESSR = 0.1, SA = 4.0, SB = 10.5;
vector<string> INSTS, VARS = {"EM_swap","C_swap","EM_2opt","C_2opt","EM_rev","C_rev","EM_ins","C_ins"};
vector<int> CPS;

enum Move { SWAP, OPT2, REV, INS, SCHED, PROBE, MIX, STAG };
struct Variant { string name; bool em; Move move; };

int n; vector<double> X, Y; mt19937_64 rng;
vector<vector<int>> pop_, kids; vector<double> len_;
struct Edge { double sum; int cnt; };
vector<char> ess; vector<Edge> edges; vector<int> touched, essKeys;
vector<int> posInB, posT; vector<char> usedPos, oxUsed;
vector<vector<int>> nbr;
long long gAtt, gAcc; double gGain;
Move curMove; int phase; double probeCur, probeNext;
double LOG1M; double STAGC = 1.0, STAGC2 = 0.3; long long genAtt, genAcc; vector<double> accHist;

double rnd() { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
double dist(int a, int b) { double dx = X[a]-X[b], dy = Y[a]-Y[b]; return sqrt(dx*dx+dy*dy); }
double tourLength(const vector<int>& t) { double s = 0; for (int k = 0; k + 1 < n; k++) s += dist(t[k], t[k+1]); return s + dist(t[n-1], t[0]); }

bool parseLine(const string& line, double v[3]) {
    const char* p = line.c_str(); int count = 0;
    while (true) { char* end; double x = strtod(p, &end); if (end == p) break; if (count < 3) v[count] = x; count++; p = end; }
    return count == 3;
}
bool loadCities(const string& path) {
    ifstream fin(path); if (!fin) return false; X.clear(); Y.clear(); string line; double v[3];
    while (getline(fin, line)) { if (!parseLine(line, v)) continue; X.push_back(v[1]); Y.push_back(v[2]); }
    n = X.size(); return true;
}

struct Cand { double avg; int key; };
bool betterCand(const Cand& a, const Cand& b) { if (a.avg != b.avg) return a.avg > b.avg; return a.key < b.key; }
void computeEssential(bool symmetric) {
    vector<double> fit(POP); for (int i = 0; i < POP; i++) fit[i] = 1.0 / tourLength(kids[i]);
    touched.clear();
    for (int pos = 0; pos < n; pos++) { int np = (pos + 1 == n) ? 0 : pos + 1;
        for (int i = 0; i < POP; i++) { int key = kids[i][pos] + kids[i][np] * n; Edge& e = edges[key];
            if (e.cnt == 0) touched.push_back(key); e.sum += fit[i]; e.cnt++; } }
    int P = touched.size(); vector<Cand> c(P);
    for (int k = 0; k < P; k++) { Edge& e = edges[touched[k]]; c[k].avg = e.sum / e.cnt; c[k].key = touched[k]; e.sum = 0; e.cnt = 0; }
    int top = (int)llround(P * ESSR); essKeys.clear();
    if (top > 0) { if (top < P) nth_element(c.begin(), c.begin() + top, c.end(), betterCand);
        for (int k = 0; k < top; k++) { int key = c[k].key; ess[key] = 1; essKeys.push_back(key);
            if (symmetric) { int f = key % n, to = key / n; int rv = to + f * n; ess[rv] = 1; essKeys.push_back(rv); } } }
}
inline bool essential(int a, int b) { return ess[a + b * n]; }

void cyclicCrossover(const vector<int>& a, const vector<int>& b, vector<int>& child) {
    for (int k = 0; k < n; k++) posInB[b[k]] = k; fill(usedPos.begin(), usedPos.end(), 0);
    bool fromB = false;
    for (int s = 0; s < n; s++) { if (usedPos[s]) continue; const vector<int>& src = fromB ? b : a; int p = s;
        do { child[p] = src[p]; usedPos[p] = 1; p = posInB[a[p]]; } while (p != s); fromB = !fromB; }
}
void orderCrossover(const vector<int>& a, const vector<int>& b, vector<int>& child) {
    fill(oxUsed.begin(), oxUsed.end(), 0);
    int c1 = (int)(rnd() * n), c2 = (int)(rnd() * n); if (c1 > c2) swap(c1, c2);
    for (int k = c1; k <= c2; k++) { child[k] = a[k]; oxUsed[a[k]] = 1; }
    int pos = (c2 + 1) % n;
    for (int t = 0; t < n; t++) { int c = b[(c2 + 1 + t) % n]; if (oxUsed[c]) continue; child[pos] = c; oxUsed[c] = 1; pos = (pos + 1) % n; }
}

bool evalMove(const vector<int>& t, Move m, int x, int y, bool needEss, double& gain, bool& essOk, int& I, int& J) {
    if (x == y) return false;
    int i, j; if (m == INS) { i = x; j = y; } else { i = min(x, y); j = max(x, y); }
    int pI = (i == 0) ? n-1 : i-1, nI = (i == n-1) ? 0 : i+1, pJ = (j == 0) ? n-1 : j-1, nJ = (j == n-1) ? 0 : j+1;
    gain = 0; essOk = true; I = i; J = j;
    if (m == SWAP) {
        int ci = t[i], cj = t[j], pi = t[pI], ni = t[nI], pj = t[pJ], nj = t[nJ];
        if (needEss) essOk = !essential(pi,ci) && !essential(ci,ni) && !essential(pj,cj) && !essential(cj,nj);
        double o, w;
        if (j - i == 1) { o = dist(pi,ci) + dist(cj,nj); w = dist(pi,cj) + dist(ci,nj); }
        else if (i == 0 && j == n-1) { o = dist(pj,cj) + dist(ci,ni); w = dist(pj,ci) + dist(cj,ni); }
        else { o = dist(pi,ci)+dist(ci,ni)+dist(pj,cj)+dist(cj,nj); w = dist(pi,cj)+dist(cj,ni)+dist(pj,ci)+dist(ci,nj); }
        gain = o - w; return true;
    }
    if (m == OPT2) {
        if (j < i + 2 || (i == 0 && j == n-1)) return false;
        int ci = t[i], ci1 = t[nI], cj = t[j], cj1 = t[nJ];
        if (needEss) essOk = !essential(ci,ci1) && !essential(cj,cj1);
        gain = dist(ci,ci1) + dist(cj,cj1) - dist(ci,cj) - dist(ci1,cj1); return true;
    }
    if (m == REV) {
        if (i == 0 && j == n-1) return false;
        int c0 = t[pI], ci = t[i], cj = t[j], cj1 = t[nJ];
        if (needEss) essOk = !essential(c0,ci) && !essential(cj,cj1);
        gain = dist(c0,ci) + dist(cj,cj1) - dist(c0,cj) - dist(ci,cj1); return true;
    }

    int ci = t[i], pi = t[pI], ni = t[nI], cj = t[j];
    if (i < j) { if (i == 0 && j == n-1) return false; int nj = t[nJ];
        if (needEss) essOk = !essential(pi,ci) && !essential(ci,ni) && !essential(cj,nj);
        gain = dist(pi,ci)+dist(ci,ni)+dist(cj,nj) - (dist(pi,ni)+dist(cj,ci)+dist(ci,nj)); return true; }
    if (pJ == i) return false; int pj = t[pJ];
    if (needEss) essOk = !essential(pj,cj) && !essential(pi,ci) && !essential(ci,ni);
    gain = dist(pj,cj)+dist(pi,ci)+dist(ci,ni) - (dist(pj,ci)+dist(ci,cj)+dist(pi,ni)); return true;
}
void applyMove(vector<int>& t, Move m, int i, int j, bool trackPos) {
    int lo, hi;
    if (m == SWAP) { swap(t[i], t[j]); if (trackPos) { posT[t[i]] = i; posT[t[j]] = j; } return; }
    if (m == OPT2) { reverse(t.begin()+i+1, t.begin()+j+1); lo = i+1; hi = j; }
    else if (m == REV) { reverse(t.begin()+i, t.begin()+j+1); lo = i; hi = j; }
    else { int ci = t[i];
        if (i < j) { for (int k = i; k < j; k++) t[k] = t[k+1]; t[j] = ci; lo = i; hi = j; }
        else { for (int k = i; k > j; k--) t[k] = t[k-1]; t[j] = ci; lo = j; hi = i; } }
    if (trackPos) for (int k = lo; k <= hi; k++) posT[t[k]] = k;
}
Move nextOf(Move m) { return m == SWAP ? INS : OPT2; }

void mutate(vector<int>& t, const Variant& v) {
    bool knn = (PARTNER == "knn");
    if (knn) for (int k = 0; k < n; k++) posT[t[k]] = k;
    double a = -1;
    while (true) {
        a += 1 + floor(log(1.0 - rnd()) / LOG1M);
        if (a >= n) break;
        int x = (int)a;
        int y = knn ? posT[nbr[t[x]][(int)(rnd() * nbr[t[x]].size())]] : (int)(rnd() * n);
        Move m = v.move;
        if (m == SCHED || m == PROBE || m == STAG) m = curMove;
        else if (m == MIX) { double u = rnd(); m = u < 1.0/3 ? SWAP : (u < 2.0/3 ? INS : OPT2); }
        double gain; bool essOk; int I, J;
        if (!evalMove(t, m, x, y, v.em, gain, essOk, I, J)) continue;
        if (v.move == PROBE && phase < 2) {
            double g2; bool e2; int I2, J2;
            probeCur += max(gain, 0.0);
            if (evalMove(t, nextOf(m), x, y, false, g2, e2, I2, J2)) probeNext += max(g2, 0.0);
        }
        bool ok = v.em ? essOk : (gain > 1e-9);
        if (ok) applyMove(t, m, I, J, knn);
        gAtt++; genAtt++; if (ok) { gAcc++; genAcc++; gGain += gain; }
    }
}

void kmeansTour(vector<int>& t, int K) {
    K = max(1, min(K, n)); vector<double> cx(K), cy(K);
    int f = (int)(rnd() * n); cx[0] = X[f]; cy[0] = Y[f]; vector<double> d2(n, INFINITY);
    for (int c = 1; c < K; c++) { double tot = 0;
        for (int i = 0; i < n; i++) { double dx = X[i]-cx[c-1], dy = Y[i]-cy[c-1]; d2[i] = min(d2[i], dx*dx+dy*dy); tot += d2[i]; }
        double r = rnd() * tot; int pick = n-1; for (int i = 0; i < n; i++) { r -= d2[i]; if (r <= 0) { pick = i; break; } }
        cx[c] = X[pick]; cy[c] = Y[pick]; }
    vector<int> lab(n, 0);
    for (int it = 0; it < 30; it++) { bool ch = false;
        for (int i = 0; i < n; i++) { int b = 0; double bd = INFINITY;
            for (int c = 0; c < K; c++) { double dx = X[i]-cx[c], dy = Y[i]-cy[c], q = dx*dx+dy*dy; if (q < bd) { bd = q; b = c; } }
            if (lab[i] != b) ch = true; lab[i] = b; }
        vector<double> sx(K,0), sy(K,0); vector<int> cn(K,0);
        for (int i = 0; i < n; i++) { sx[lab[i]] += X[i]; sy[lab[i]] += Y[i]; cn[lab[i]]++; }
        for (int c = 0; c < K; c++) if (cn[c]) { cx[c] = sx[c]/cn[c]; cy[c] = sy[c]/cn[c]; }
        if (!ch && it > 0) break; }
    int p = 0;
    for (int c = 0; c < K; c++) { int s0 = p; for (int i = 0; i < n; i++) if (lab[i] == c) t[p++] = i;
        if (!FILEORDER) for (int k = p-1; k > s0; k--) swap(t[k], t[s0 + (int)(rnd() * (k - s0 + 1))]); }
}
int kclusters() {
    if (KSPEC == "sqrt") return (int)llround(sqrt((double)n));
    if (KSPEC == "n10") return max(2, (int)llround(n/10.0));
    if (KSPEC == "n30") return max(2, (int)llround(n/30.0));
    return atoi(KSPEC.c_str());
}
void initPopulation() {
    int K = kclusters();
    for (int i = 0; i < POP; i++) {
        if (INIT == "kmeans") { kmeansTour(pop_[i], K); continue; }
        if (INIT == "nn") { vector<char> used(n, 0); int cur = (int)(rnd() * n); pop_[i][0] = cur; used[cur] = 1;
            for (int k = 1; k < n; k++) { int best = -1; double bd = INFINITY;
                for (int c = 0; c < n; c++) if (!used[c]) { double d = dist(cur, c); if (d < bd) { bd = d; best = c; } }
                pop_[i][k] = best; used[best] = 1; cur = best; }
            continue; }
        iota(pop_[i].begin(), pop_[i].end(), 0);
        for (int k = n-1; k > 0; k--) swap(pop_[i][k], pop_[i][(int)(rnd() * (k+1))]);
    }
}

void nextGeneration(const Variant& v) {
    vector<int> order(POP); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [](int a, int b) { if (len_[a] != len_[b]) return len_[a] < len_[b]; return a < b; });
    vector<int> sel(POP);
    for (int k = 0; k < ELITE; k++) sel[k] = order[k];
    if (SEL == "tour") {
        for (int k = ELITE; k < POP; k++) { int b = (int)(rnd() * POP);
            for (int q = 1; q < TK; q++) { int c = (int)(rnd() * POP); if (len_[c] < len_[b]) b = c; } sel[k] = b; }
    } else {
        vector<double> cum(POP); double tot = 0;
        for (int k = 0; k < POP; k++) { tot += 1.0 / len_[order[k]]; cum[k] = tot; }
        for (int k = 0; k < POP; k++) cum[k] = 100.0 * cum[k] / tot;
        for (int k = ELITE; k < POP; k++) { double pick = 100.0 * rnd(); int j = lower_bound(cum.begin(), cum.end(), pick) - cum.begin();
            if (j == POP) j = POP-1; sel[k] = order[j]; }
    }
    vector<int> sh(POP); iota(sh.begin(), sh.end(), 0);
    for (int k = POP-1; k > 0; k--) swap(sh[k], sh[(int)(rnd() * (k+1))]);
    for (int i = 0; i < ELITE; i++) kids[i] = pop_[sel[i]];
    for (int i = 0; i < POP - ELITE; i++) {
        const vector<int>& p1 = pop_[sel[sh[i]]]; const vector<int>& p2 = pop_[sel[sh[POP-1-i]]];
        if (XOVER == "ox") orderCrossover(p1, p2, kids[ELITE+i]); else cyclicCrossover(p1, p2, kids[ELITE+i]);
    }
    if (v.em) computeEssential(v.move != SWAP);
    for (int i = 0; i < POP; i++) mutate(kids[i], v);
    if (v.em) for (int key : essKeys) ess[key] = 0;
    swap(pop_, kids);
}

struct CP { int gen; double best; long long att, acc; double gain; int phase; };
void runGA(const Variant& v, bool curve, double& bestLen, int& convGen, vector<int>& bestTour, vector<double>& fr, vector<CP>& cps, int& sw1, int& sw2) {
    initPopulation(); bestLen = INFINITY; convGen = 1; fr.clear(); cps.clear(); gAtt = gAcc = 0; gGain = 0;
    phase = 0; curMove = SWAP; sw1 = sw2 = -1; accHist.clear();
    for (int g = 1; g <= GENS; g++) {
        double mn = INFINITY; int bi = 0;
        for (int i = 0; i < POP; i++) { len_[i] = tourLength(pop_[i]); if (len_[i] < mn) { mn = len_[i]; bi = i; } }
        double x = round(mn * 10.0) / 10.0;
        if (x < bestLen) { bestLen = x; bestTour = pop_[bi]; convGen = g; }
        if (curve) fr.push_back(x);
        if (g == GENS || binary_search(CPS.begin(), CPS.end(), g)) { cps.push_back({g, bestLen, gAtt, gAcc, gGain, phase}); gAtt = gAcc = 0; gGain = 0; }
        if (g < GENS) {

            if (v.move == SCHED) { double gpc = (double)g / n; phase = gpc < SA ? 0 : (gpc < SB ? 1 : 2); curMove = phase == 0 ? SWAP : (phase == 1 ? INS : OPT2);
                if (phase >= 1 && sw1 < 0) sw1 = g; if (phase == 2 && sw2 < 0) sw2 = g; }
            probeCur = probeNext = 0; genAtt = genAcc = 0;
            nextGeneration(v);
            if (v.move == STAG && phase < 2) { accHist.push_back(genAtt ? (double)genAcc / genAtt : 0);
                if (accHist.size() >= 3) { double m3 = (accHist[accHist.size()-1] + accHist[accHist.size()-2] + accHist[accHist.size()-3]) / 3;
                    if (m3 < (phase == 0 ? STAGC : STAGC2) / n) { phase++; curMove = phase == 1 ? INS : OPT2; if (phase == 1) sw1 = g; else sw2 = g; accHist.clear(); } } }
            if (v.move == PROBE && phase < 2 && probeNext > probeCur) { phase++; curMove = phase == 1 ? INS : OPT2; if (phase == 1) sw1 = g; else sw2 = g; }
        }
    }
}

double meanv(const vector<double>& a) { double s = 0; for (double x : a) s += x; return s / a.size(); }
double sdv(const vector<double>& a) { double m = meanv(a), s = 0; for (double x : a) s += (x-m)*(x-m); return a.size() > 1 ? sqrt(s / (a.size()-1)) : 0; }
vector<string> split(const string& s) { vector<string> r; string c; stringstream ss(s); while (getline(ss, c, ',')) if (!c.empty()) r.push_back(c); return r; }

int main(int argc, char** argv) {
    string cpArg;
    for (int k = 1; k < argc; k++) { string a = argv[k]; auto nx = [&]() { return string(argv[++k]); };
        if (a == "--dir") DIR = nx(); else if (a == "--out") OUT = nx(); else if (a == "--inst") INSTS = split(nx());
        else if (a == "--runs") NRUNS = stoi(nx()); else if (a == "--gens") GENS = stoi(nx()); else if (a == "--gk") GK = stoi(nx());
        else if (a == "--pop") POP = stoi(nx()); else if (a == "--elite") ELITE = stoi(nx()); else if (a == "--pm") PM = stod(nx());
        else if (a == "--ess") ESSR = stod(nx()); else if (a == "--init") INIT = nx(); else if (a == "--kspec") KSPEC = nx();
        else if (a == "--fileorder") FILEORDER = stoi(nx()); else if (a == "--xover") XOVER = nx(); else if (a == "--sel") SEL = nx();
        else if (a == "--tk") TK = stoi(nx()); else if (a == "--partner") PARTNER = nx(); else if (a == "--knn") KNN = stoi(nx());
        else if (a == "--variants") VARS = split(nx()); else if (a == "--sa") SA = stod(nx()); else if (a == "--sb") SB = stod(nx());
        else if (a == "--cp") cpArg = nx(); else if (a == "--stagc") STAGC = stod(nx()); else if (a == "--stagc2") STAGC2 = stod(nx());
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; } }
    LOG1M = log(1.0 - PM);
    map<string, Move> mv = {{"swap",SWAP},{"2opt",OPT2},{"rev",REV},{"ins",INS},{"sched",SCHED},{"probe",PROBE},{"mix",MIX},{"stag",STAG}};
    vector<Variant> V;
    for (auto& s : VARS) { auto p = s.find('_'); string r = s.substr(0, p), m = s.substr(p+1);
        if (!mv.count(m) || (r != "C" && r != "EM")) { fprintf(stderr, "bad variant %s\n", s.c_str()); return 2; }
        V.push_back({s, r == "EM", mv[m]}); }
    fs::create_directories(OUT);
    FILE* diag = fopen((OUT + "/DIAG.csv").c_str(), "a");
    FILE* summ = fopen((OUT + "/SUMMARY.csv").c_str(), "a");
    for (const string& inst : INSTS) {
        if (!loadCities(DIR + "/" + inst + ".tsp")) { fprintf(stderr, "cannot open %s\n", inst.c_str()); return 1; }
        if (GK) GENS = GK * n;
        CPS.clear();
        if (!cpArg.empty()) for (auto& s : split(cpArg)) CPS.push_back(stoi(s));
        else { for (double f : {0.005,0.0125,0.025,0.05,0.1,0.15,0.25,0.375,0.5,0.75,1.0}) CPS.push_back(max(1, (int)llround(f * GENS))); CPS.push_back(1); }
        sort(CPS.begin(), CPS.end()); CPS.erase(unique(CPS.begin(), CPS.end()), CPS.end());
        pop_.assign(POP, vector<int>(n)); kids.assign(POP, vector<int>(n)); len_.assign(POP, 0);
        ess.assign((size_t)n*n, 0); edges.assign((size_t)n*n, Edge{0,0}); posInB.assign(n,0); usedPos.assign(n,0); oxUsed.assign(n,0); posT.assign(n,0);
        if (PARTNER == "knn") { nbr.assign(n, {}); for (int c = 0; c < n; c++) { vector<pair<double,int>> d;
            for (int o = 0; o < n; o++) if (o != c) d.push_back({dist(c,o), o});
            int K = min(KNN, n-1); partial_sort(d.begin(), d.begin()+K, d.end()); for (int q = 0; q < K; q++) nbr[c].push_back(d[q].second); } }
        for (const Variant& v : V) {
            string tag = inst + "_" + v.name;
            if (fs::exists(OUT + "/" + tag + "_runs.csv")) { printf("skip %s (done)\n", tag.c_str()); continue; }
            vector<double> best, tm, conv; string rows; vector<double> fr1;
            for (int r = 1; r <= NRUNS; r++) {
                rng.seed(1 + r); double bl; int cg, s1, s2; vector<int> bt; vector<double> fr; vector<CP> cps;
                auto t0 = chrono::steady_clock::now();
                runGA(v, r == 1, bl, cg, bt, fr, cps, s1, s2);
                double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
                vector<char> seen(n, 0); for (int c : bt) { if (seen[c]) { fprintf(stderr, "invalid tour\n"); return 1; } seen[c] = 1; }
                best.push_back(bl); tm.push_back(sec); conv.push_back(cg); if (r == 1) fr1 = fr;
                char buf[256]; snprintf(buf, sizeof buf, "%d,%.1f,%.4f,%d,%d,%d\n", r, bl, sec, cg, s1, s2); rows += buf;
                for (auto& c : cps) fprintf(diag, "%s,%d,%s,%d,%d,%d,%.1f,%lld,%lld,%.4f,%d\n", inst.c_str(), n, v.name.c_str(), r, GENS, c.gen, c.best, c.att, c.acc, c.gain, c.phase);
                printf("%s %s run %d: %.1f (%.1fs)\n", inst.c_str(), v.name.c_str(), r, bl, sec); fflush(stdout);
            }
            fflush(diag);
            FILE* f = fopen((OUT + "/" + tag + "_runs.tmp").c_str(), "w"); fprintf(f, "Run,Best,Time,ConvGen,Switch1,Switch2\n%s", rows.c_str()); fclose(f);
            fs::rename(OUT + "/" + tag + "_runs.tmp", OUT + "/" + tag + "_runs.csv");
            FILE* g = fopen((OUT + "/" + tag + "_FR.csv").c_str(), "w"); fprintf(g, "Best\n"); for (double x : fr1) fprintf(g, "%.1f\n", x); fclose(g);
            fprintf(summ, "%s,%d,%s,%d,%.1f,%.4f,%.4f,%.2f,%.4f\n", inst.c_str(), n, v.name.c_str(), GENS, *min_element(best.begin(), best.end()), meanv(best), sdv(best), meanv(conv), meanv(tm));
            fflush(summ);
        }
    }
    fclose(diag); fclose(summ);
    return 0;
}
