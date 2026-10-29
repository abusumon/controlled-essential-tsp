#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

const string TSP_DIR = "inputs";
const string OUT_DIR = "output";

const int    POP_SIZE    = 300;
const int    ELITE_SIZE  = 60;
const int    GENERATIONS = 2000;
const double MUT_RATE    = 0.005;
const double ESS_RATE    = 0.1;

const int N_RUNS    = 30;
const int BASE_SEED = 1;

const vector<string> INST_ONLY = {};

enum Move { SWAP, OPT2, REV, INS };

struct Variant {
    string name;
    bool useEM;
    Move move;
};

const vector<Variant> VARIANTS = {
    {"EM_swap", true,  SWAP}, {"C_swap", false, SWAP},
    {"EM_2opt", true,  OPT2}, {"C_2opt", false, OPT2},
    {"EM_rev",  true,  REV},  {"C_rev",  false, REV},
    {"EM_ins",  true,  INS},  {"C_ins",  false, INS}
};

int n;
vector<double> X, Y;
vector<string> cityId;

mt19937_64 rng;

vector<vector<int>> pop, children;
vector<double> len;

struct Edge {
    double sum;
    int cnt;
};
vector<char> ess;
vector<Edge> edges;
vector<int> touched;
vector<int> essKeys;

vector<int> posInB;
vector<char> usedPos;

double rnd() {
    return (rng() >> 11) * (1.0 / 9007199254740992.0);
}

double dist(int a, int b) {
    double dx = X[a] - X[b];
    double dy = Y[a] - Y[b];
    return sqrt(dx * dx + dy * dy);
}

double tourLength(const vector<int>& t) {
    double s = 0;
    for (int k = 0; k + 1 < n; k++) s += dist(t[k], t[k + 1]);
    return s + dist(t[n - 1], t[0]);
}

bool parseLine(const string& line, double v[3]) {
    const char* p = line.c_str();
    int count = 0;
    while (true) {
        char* end;
        double x = strtod(p, &end);
        if (end == p) break;
        if (count < 3) v[count] = x;
        count++;
        p = end;
    }
    return count == 3;
}

bool loadCities(const string& path) {
    ifstream fin(path);
    if (!fin) return false;
    X.clear(); Y.clear(); cityId.clear();
    string line;
    double v[3];
    while (getline(fin, line)) {
        if (!parseLine(line, v)) continue;
        char buf[40];
        snprintf(buf, sizeof(buf), "%.15g", v[0]);
        cityId.push_back(buf);
        X.push_back(v[1]);
        Y.push_back(v[2]);
    }
    n = X.size();
    return true;
}

struct Cand {
    double avg;
    int key;
};

bool betterCand(const Cand& a, const Cand& b) {
    if (a.avg != b.avg) return a.avg > b.avg;
    return a.key < b.key;
}

void computeEssential(bool symmetric) {
    vector<double> fit(POP_SIZE);
    for (int i = 0; i < POP_SIZE; i++) fit[i] = 1.0 / tourLength(children[i]);

    touched.clear();
    for (int pos = 0; pos < n; pos++) {
        int nextPos = (pos + 1 == n) ? 0 : pos + 1;
        for (int i = 0; i < POP_SIZE; i++) {
            int key = children[i][pos] + children[i][nextPos] * n;
            Edge& e = edges[key];
            if (e.cnt == 0) touched.push_back(key);
            e.sum += fit[i];
            e.cnt++;
        }
    }

    int poolSize = touched.size();
    vector<Cand> cand(poolSize);
    for (int k = 0; k < poolSize; k++) {
        Edge& e = edges[touched[k]];
        cand[k].avg = e.sum / e.cnt;
        cand[k].key = touched[k];
        e.sum = 0;
        e.cnt = 0;
    }

    int topN = (int)llround(poolSize * ESS_RATE);
    essKeys.clear();
    if (topN > 0) {
        if (topN < poolSize)
            nth_element(cand.begin(), cand.begin() + topN, cand.end(), betterCand);
        for (int k = 0; k < topN; k++) {
            int key = cand[k].key;
            ess[key] = 1;
            essKeys.push_back(key);
            if (symmetric) {
                int from = key % n, to = key / n;
                int rev = to + from * n;
                ess[rev] = 1;
                essKeys.push_back(rev);
            }
        }
    }
}

bool essential(int from, int to) {
    return ess[from + to * n];
}

void cyclicCrossover(const vector<int>& a, const vector<int>& b, vector<int>& child) {
    for (int k = 0; k < n; k++) posInB[b[k]] = k;
    fill(usedPos.begin(), usedPos.end(), 0);

    bool takeFromB = false;
    for (int start = 0; start < n; start++) {
        if (usedPos[start]) continue;
        const vector<int>& src = takeFromB ? b : a;
        int p = start;
        do {
            child[p] = src[p];
            usedPos[p] = 1;
            p = posInB[a[p]];
        } while (p != start);
        takeFromB = !takeFromB;
    }
}

const double LOG_1M = log(1.0 - MUT_RATE);

void mutate(vector<int>& t, bool useEM, Move move) {
    double a = -1;
    while (true) {
        a += 1 + floor(log(1.0 - rnd()) / LOG_1M);
        if (a >= n) break;

        int x = (int)a;
        int y = (int)(rnd() * n);
        if (x == y) continue;

        int i, j;
        if (move == INS) { i = x; j = y; }
        else             { i = min(x, y); j = max(x, y); }

        int prevI = (i == 0) ? n - 1 : i - 1;
        int nextI = (i == n - 1) ? 0 : i + 1;
        int prevJ = (j == 0) ? n - 1 : j - 1;
        int nextJ = (j == n - 1) ? 0 : j + 1;

        bool ok = false;

        if (move == SWAP) {
            int ci = t[i], cj = t[j];
            int pi = t[prevI], ni = t[nextI];
            int pj = t[prevJ], nj = t[nextJ];

            if (useEM) {
                ok = !essential(pi, ci) && !essential(ci, ni) &&
                     !essential(pj, cj) && !essential(cj, nj);
            } else {
                double oldCost, newCost;
                if (j - i == 1) {
                    oldCost = dist(pi, ci) + dist(cj, nj);
                    newCost = dist(pi, cj) + dist(ci, nj);
                } else if (i == 0 && j == n - 1) {
                    oldCost = dist(pj, cj) + dist(ci, ni);
                    newCost = dist(pj, ci) + dist(cj, ni);
                } else {
                    oldCost = dist(pi, ci) + dist(ci, ni) + dist(pj, cj) + dist(cj, nj);
                    newCost = dist(pi, cj) + dist(cj, ni) + dist(pj, ci) + dist(ci, nj);
                }
                ok = (oldCost - newCost) > 1e-9;
            }
            if (ok) swap(t[i], t[j]);
        }
        else if (move == OPT2) {
            if (j < i + 2 || (i == 0 && j == n - 1)) continue;
            int ci = t[i], ci1 = t[nextI];
            int cj = t[j], cj1 = t[nextJ];

            if (useEM) {
                ok = !essential(ci, ci1) && !essential(cj, cj1);
            } else {
                double gain = dist(ci, ci1) + dist(cj, cj1) - dist(ci, cj) - dist(ci1, cj1);
                ok = gain > 1e-9;
            }
            if (ok) reverse(t.begin() + i + 1, t.begin() + j + 1);
        }
        else if (move == REV) {
            if (i == 0 && j == n - 1) continue;
            int c0 = t[prevI], ci = t[i];
            int cj = t[j], cj1 = t[nextJ];

            if (useEM) {
                ok = !essential(c0, ci) && !essential(cj, cj1);
            } else {
                double gain = dist(c0, ci) + dist(cj, cj1) - dist(c0, cj) - dist(ci, cj1);
                ok = gain > 1e-9;
            }
            if (ok) reverse(t.begin() + i, t.begin() + j + 1);
        }
        else {
            int ci = t[i], pi = t[prevI], ni = t[nextI], cj = t[j];

            if (i < j) {
                if (i == 0 && j == n - 1) continue;
                int nj = t[nextJ];
                if (useEM) {
                    ok = !essential(pi, ci) && !essential(ci, ni) && !essential(cj, nj);
                } else {
                    double oldCost = dist(pi, ci) + dist(ci, ni) + dist(cj, nj);
                    double newCost = dist(pi, ni) + dist(cj, ci) + dist(ci, nj);
                    ok = (oldCost - newCost) > 1e-9;
                }
                if (ok) {
                    for (int k = i; k < j; k++) t[k] = t[k + 1];
                    t[j] = ci;
                }
            } else {
                if (prevJ == i) continue;
                int pj = t[prevJ];
                if (useEM) {
                    ok = !essential(pj, cj) && !essential(pi, ci) && !essential(ci, ni);
                } else {
                    double oldCost = dist(pj, cj) + dist(pi, ci) + dist(ci, ni);
                    double newCost = dist(pj, ci) + dist(ci, cj) + dist(pi, ni);
                    ok = (oldCost - newCost) > 1e-9;
                }
                if (ok) {
                    for (int k = i; k > j; k--) t[k] = t[k - 1];
                    t[j] = ci;
                }
            }
        }
    }
}

void initPopulation() {
    for (int i = 0; i < POP_SIZE; i++) {
        iota(pop[i].begin(), pop[i].end(), 0);
        for (int k = n - 1; k > 0; k--) swap(pop[i][k], pop[i][(int)(rnd() * (k + 1))]);
    }
}

void nextGeneration(const Variant& v) {

    vector<int> order(POP_SIZE);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [](int a, int b) {
        if (len[a] != len[b]) return len[a] < len[b];
        return a < b;
    });

    vector<double> cum(POP_SIZE);
    double total = 0;
    for (int k = 0; k < POP_SIZE; k++) {
        total += 1.0 / len[order[k]];
        cum[k] = total;
    }
    for (int k = 0; k < POP_SIZE; k++) cum[k] = 100.0 * cum[k] / total;

    vector<int> selected(POP_SIZE);
    for (int k = 0; k < ELITE_SIZE; k++) selected[k] = order[k];
    for (int k = ELITE_SIZE; k < POP_SIZE; k++) {
        double pick = 100.0 * rnd();
        int j = lower_bound(cum.begin(), cum.end(), pick) - cum.begin();
        if (j == POP_SIZE) j = POP_SIZE - 1;
        selected[k] = order[j];
    }

    vector<int> shuf(POP_SIZE);
    iota(shuf.begin(), shuf.end(), 0);
    for (int k = POP_SIZE - 1; k > 0; k--) swap(shuf[k], shuf[(int)(rnd() * (k + 1))]);

    for (int i = 0; i < ELITE_SIZE; i++) children[i] = pop[selected[i]];
    for (int i = 0; i < POP_SIZE - ELITE_SIZE; i++) {
        const vector<int>& p1 = pop[selected[shuf[i]]];
        const vector<int>& p2 = pop[selected[shuf[POP_SIZE - 1 - i]]];
        cyclicCrossover(p1, p2, children[ELITE_SIZE + i]);
    }

    if (v.useEM) {
        bool symmetric = (v.move != SWAP);
        computeEssential(symmetric);
    }

    for (int i = 0; i < POP_SIZE; i++) mutate(children[i], v.useEM, v.move);

    if (v.useEM) {
        for (int key : essKeys) ess[key] = 0;
    }

    swap(pop, children);
}

void runGA(const Variant& v, bool recordCurve,
           double& bestLen, int& convGen, vector<int>& bestTour, vector<double>& curve) {
    initPopulation();
    bestLen = INFINITY;
    convGen = 1;
    curve.clear();

    for (int g = 1; g <= GENERATIONS; g++) {
        double mn = INFINITY;
        int bi = 0;
        for (int i = 0; i < POP_SIZE; i++) {
            len[i] = tourLength(pop[i]);
            if (len[i] < mn) { mn = len[i]; bi = i; }
        }
        double x = round(mn * 10.0) / 10.0;
        if (x < bestLen) {
            bestLen = x;
            bestTour = pop[bi];
            convGen = g;
        }
        if (recordCurve) curve.push_back(x);

        if (g < GENERATIONS) nextGeneration(v);
    }
}

double mean(const vector<double>& a) {
    double s = 0;
    for (double x : a) s += x;
    return s / a.size();
}

double stdDev(const vector<double>& a) {
    double m = mean(a), s = 0;
    for (double x : a) s += (x - m) * (x - m);
    return sqrt(s / (a.size() - 1));
}

double ranksumP(const vector<double>& x, const vector<double>& y) {
    int n1 = x.size(), n2 = y.size(), N = n1 + n2;
    vector<pair<double, int>> all;
    for (double v : x) all.push_back({v, 0});
    for (double v : y) all.push_back({v, 1});
    sort(all.begin(), all.end());

    double W = 0;
    double tieSum = 0;
    for (int i = 0; i < N;) {
        int j = i;
        while (j < N && all[j].first == all[i].first) j++;
        double rank = (i + 1 + j) / 2.0;
        for (int k = i; k < j; k++)
            if (all[k].second == 0) W += rank;
        double t = j - i;
        tieSum += t * t * t - t;
        i = j;
    }

    double mu = n1 * (N + 1) / 2.0;
    double var = (double)n1 * n2 / 12.0 * ((N + 1) - tieSum / ((double)N * (N - 1)));
    double diff = W - mu;
    if (var <= 0 || diff == 0) return 1.0;
    double z = (fabs(diff) - 0.5) / sqrt(var);
    return erfc(z / sqrt(2.0));
}

string toLower(string s) {
    for (char& c : s) c = tolower((unsigned char)c);
    return s;
}

int main() {
    fs::create_directories(OUT_DIR);

    char stamp[32];
    time_t now = time(nullptr);
    strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", localtime(&now));

    vector<fs::path> files;
    for (auto& entry : fs::directory_iterator(TSP_DIR))
        if (entry.is_regular_file()) files.push_back(entry.path());
    sort(files.begin(), files.end());

    vector<string> only;
    for (const string& s : INST_ONLY) only.push_back(toLower(s));

    string summaryPath = OUT_DIR + "/ALL_EMvsC_summary_" + stamp + ".csv";
    string summaryText = "Instance,Variant,Best,Avg,Std,AvgConvGen,StdConvGen,AvgTime,p_vs_twin\n";

    int nV = VARIANTS.size();

    for (const fs::path& file : files) {
        string instName = file.stem().string();
        if (!only.empty() && find(only.begin(), only.end(), toLower(instName)) == only.end())
            continue;

        if (!loadCities(file.string())) {
            printf("Cannot open %s\n", file.string().c_str());
            return 1;
        }
        if (n < 5) continue;

        pop.assign(POP_SIZE, vector<int>(n));
        children.assign(POP_SIZE, vector<int>(n));
        len.assign(POP_SIZE, 0.0);
        ess.assign((size_t)n * n, 0);
        edges.assign((size_t)n * n, Edge{0.0, 0});
        posInB.assign(n, 0);
        usedPos.assign(n, 0);

        vector<vector<double>> allBest(nV, vector<double>(N_RUNS));
        vector<vector<double>> allTime(nV, vector<double>(N_RUNS));
        vector<vector<double>> allConv(nV, vector<double>(N_RUNS));

        for (int vi = 0; vi < nV; vi++) {
            const Variant& v = VARIANTS[vi];
            printf("\n=== %s | %s ===\n", instName.c_str(), v.name.c_str());

            vector<string> routes(N_RUNS);
            vector<double> firstRunCurve;

            for (int r = 1; r <= N_RUNS; r++) {
                rng.seed(BASE_SEED + r);

                double bestLen;
                int convGen;
                vector<int> bestTour;
                vector<double> curve;

                auto t0 = chrono::steady_clock::now();
                runGA(v, r == 1, bestLen, convGen, bestTour, curve);
                double seconds = chrono::duration<double>(chrono::steady_clock::now() - t0).count();

                vector<char> seen(n, 0);
                for (int c : bestTour) {
                    if (seen[c]) {
                        printf("Invalid tour produced by %s (run %d)\n", v.name.c_str(), r);
                        return 1;
                    }
                    seen[c] = 1;
                }

                allBest[vi][r - 1] = bestLen;
                allTime[vi][r - 1] = seconds;
                allConv[vi][r - 1] = convGen;

                string route = "[";
                for (int k = 0; k < n; k++) {
                    if (k > 0) route += ";";
                    route += cityId[bestTour[k]];
                }
                routes[r - 1] = route + "]";

                if (r == 1) firstRunCurve = curve;
                printf("RUN %2d: %.1f  (%.1f s, converged at gen %d)\n", r - 1, bestLen, seconds, convGen);
                fflush(stdout);
            }

            string tag = instName + "_" + v.name;

            FILE* f1 = fopen((OUT_DIR + "/" + tag + "_runs_" + stamp + ".csv").c_str(), "w");
            fprintf(f1, "Best,ExecutionTime,ConvGen,Plot\n");
            for (int r = 0; r < N_RUNS; r++)
                fprintf(f1, "%.1f,%.6f,%d,%s\n", allBest[vi][r], allTime[vi][r],
                        (int)allConv[vi][r], routes[r].c_str());
            fclose(f1);

            FILE* f2 = fopen((OUT_DIR + "/" + tag + "_FR_" + stamp + ".csv").c_str(), "w");
            fprintf(f2, "Best\n");
            for (double x : firstRunCurve) fprintf(f2, "%.1f\n", x);
            fclose(f2);
        }

        printf("\n--- %s summary over %d runs ---\n", instName.c_str(), N_RUNS);
        printf("%-8s %10s %10s %8s %9s %9s %8s %9s\n",
               "Variant", "Best", "Avg", "Std", "AvgConv", "StdConv", "AvgTime", "p_twin");

        for (int vi = 0; vi < nV; vi++) {
            double best = *min_element(allBest[vi].begin(), allBest[vi].end());
            double avg = mean(allBest[vi]);
            double sd = stdDev(allBest[vi]);
            double avgConv = mean(allConv[vi]);
            double sdConv = stdDev(allConv[vi]);
            double avgTime = mean(allTime[vi]);
            double p = ranksumP(allBest[vi], allBest[vi ^ 1]);

            printf("%-8s %10.1f %10.2f %8.2f %9.1f %9.1f %8.2f %9.4f\n",
                   VARIANTS[vi].name.c_str(), best, avg, sd, avgConv, sdConv, avgTime, p);

            char row[256];
            snprintf(row, sizeof(row), "%s,%s,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6g\n",
                     instName.c_str(), VARIANTS[vi].name.c_str(), best, avg, sd,
                     avgConv, sdConv, avgTime, p);
            summaryText += row;
        }

        FILE* fs_ = fopen(summaryPath.c_str(), "w");
        fputs(summaryText.c_str(), fs_);
        fclose(fs_);
    }

    printf("\nCombined summary: %s\n", summaryPath.c_str());
    return 0;
}
