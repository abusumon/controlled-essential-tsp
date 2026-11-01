# Controlled versus Essential Mutation in a Genetic Algorithm for the TSP

Code, raw results and analysis scripts of an undergraduate thesis at the Mathematics Discipline, Khulna University, Bangladesh.

**Author:** Abu Bakkar Siddik · **Supervisor:** Md. Azizur Rahman

The study looks at mutation in a simple genetic algorithm (GA) for the travelling salesman problem (TSP). It uses four moves (swap, insertion, 2-opt and reversion) on 30 TSPLIB instances with 51 to 1432 cities, and it answers four questions.

| Question | Short answer |
|:--|:--|
| **RQ1.** Is controlled mutation (keep a move only if it shortens the tour) better than essential mutation (EM: keep a move only if it removes no edge that is common in good tours)? | Yes. Controlled mutation wins 111 of 120 per-instance comparisons and loses none. It stays ahead with 20,000 generations, from clustered starts and with order crossover and tournament selection. |
| **RQ2.** Does each controlled move have its own range of instances? | At 2000 generations, yes: 2-opt or reversion up to about 150 cities, insertion from about 200 to 440, swap from about 500 (three exceptions). The ranges shift with the budget, the partner rule, the starting tours and the GA. |
| **RQ3.** Why do the moves behave differently? | Each move converges to its own local optimum, and speed and depth are ordered in opposite ways (swap fastest and shallowest, 2-opt slowest and deepest). With a uniformly random partner, the hand-overs between moves come after a number of generations proportional to the number of cities. An acceptance-triggered schedule (C_stag) uses this. |
| **RQ4.** What if the initial population is clustered? | Final tours become about four times closer to the best-known tours, and the best move shifts towards 2-opt. Keeping the TSPLIB file order inside clusters adds a head start that comes from the files. |

---

## Repository layout

```
controlled-essential-tsp/
├── src/
│   ├── em_vs_c.cpp            EM and controlled variants                              E1, E4
│   ├── em_vs_c_diag.cpp       same GA with checkpoints and move statistics            E2, E3
│   ├── em_vs_c_cem.cpp        controlled essential mutation                           CEM, E5
│   ├── q1ga.cpp               unified GA with all options of E7 to E14
│   └── localopt.cpp           descent of one move to a local optimum                  E6
├── scripts/
│   └── run_experiments.py     runs E7 to E14 in parallel with q1ga
├── analysis/
│   ├── common.py              paths, best-known tour lengths, helpers
│   ├── summarize.py           summaries and tests of E7 to E14  ->  results/
│   ├── figures.py             Figures 13 to 16                  ->  figures/
│   └── city_groups.py         effective number of city groups   ->  results/city_groups.csv
├── notebooks/
│   └── analysis.ipynb         tables and figures of E1 to E5    ->  results/, figures/
├── inputs/                    30 TSPLIB instances (.tsp)
│   └── synthetic/             30 synthetic instances with k city groups (E9)
├── data/                      raw results of every experiment (see below)
├── results/                   summary tables (CSV, LaTeX)
├── figures/                   all figures (PNG)
├── requirements.txt
└── LICENSE
```

---

## Experiments

All experiments use the essential rate ρ = 0.1, the mutation rate p_m = 0.005 (unless stated) and the seed `1 + r` for run `r`, so runs are paired across variants and experiments.

| Exp. | Program | Instances | Start | Generations | Runs | N / E | What changes | Folder in `data/` |
|:--|:--|:--|:--|:--|--:|:--|:--|:--|
| E1 | `em_vs_c` | 30 TSPLIB | random | 2000 | 30 | 300 / 60 | 4 moves × {EM, controlled} | `E1_benchmark` |
| E2 | `em_vs_c_diag` | 6 TSPLIB | random | 2000 | 20 | 200 / 50 | checkpoints, move statistics | `E2_trajectories` |
| E3 | `em_vs_c_diag` | d493, u574 | random | 10,000 | 20 | 200 / 50 | long budget | `E3_long_budget` |
| E4 | `em_vs_c` | d493, u574 | random | 20,000 | 50 | 500 / 100 | long budget, large population | `E4_large_population` |
| CEM | `em_vs_c_cem` | 30 TSPLIB | random | 2000 | 30 | 300 / 60 | controlled essential mutation | `CEM` |
| E5 | `em_vs_c_cem` | 30 TSPLIB | random | 2000 | 30 | 200 / 50 | CEM with p_m = 0.05 | `E5_CEM_high_rate` |
| E6 | `localopt` | 30 TSPLIB | random | descent | 3 | – | local optimum of each move | `E6_local_optima` |
| E7 | `q1ga` | 30 TSPLIB | random | 2n, 7n | 30 | 300 / 60 | budget proportional to n | `E7_budget_2n`, `E7_budget_7n` |
| E8 | `q1ga` | 30 TSPLIB | nearest neighbour | 2000 | 30 | 300 / 60 | good starting tours | `E8_nn_start` |
| E9 | `q1ga` | 30 synthetic | random | 2000 | 30 | 300 / 60 | number of city groups | `E9_synthetic` |
| E10 | `q1ga` | 30 TSPLIB | k-means clusters | 2000 | 30 | 300 / 60 | k = 10, √n, n/30, n/10, file order | `E10_clustered_*` |
| E11 | `q1ga` | 30 TSPLIB | clustered, k = √n | 2000 | 30 | 300 / 60 | EM from clustered starts | `E11_clustered_em` |
| E12 | `q1ga` | 30 TSPLIB | random | 2000 | 30 | 300 / 60 | order crossover, tournament selection (size 3) | `E12_ox_tournament` |
| E13 | `q1ga` | 30 TSPLIB | random | 2n, 2000, 7n | 30 | 300 / 60 | move schedules C_stag, C_sched, C_mix | `E13_schedules_*` |
| E14 | `q1ga` | 30 TSPLIB | random | 2n, 2000 | 30 | 300 / 60 | partner from the 10 nearest cities | `E14_knn_partner_*` |

`q1ga` reproduces the E1 runs exactly when it is run with the E1 settings.

---

## Requirements

- A C++17 compiler (tested with g++ on Windows/MinGW-w64 and Linux). No external libraries.
- Python 3.9 or newer: `pip install -r requirements.txt`

---

## Building

From the repository root:

```bash
g++ -O2 -std=c++17 src/em_vs_c.cpp      -o em_vs_c
g++ -O2 -std=c++17 src/em_vs_c_diag.cpp -o em_vs_c_diag
g++ -O2 -std=c++17 src/em_vs_c_cem.cpp  -o em_vs_c_cem
g++ -O2 -std=c++17 src/q1ga.cpp         -o q1ga
g++ -O2 -std=c++17 src/localopt.cpp     -o localopt
```

On Windows the same commands produce `.exe` files. Older MinGW versions may need `-lstdc++fs`.

---

## Running the experiments

### E1 to E5

The three `em_vs_c*` programs read every `.tsp` file in `inputs/` (or only the instances in `INST_ONLY`) and write time-stamped files to `output/`. Run them from the repository root. The settings are constants at the top of each source file; the files ship with the settings of E1, E3 and CEM.

| Exp. | Program | `POP_SIZE` | `ELITE_SIZE` | `GENERATIONS` | `MUT_RATE` | `N_RUNS` | `INST_ONLY` |
|:--|:--|--:|--:|--:|--:|--:|:--|
| E1 | `em_vs_c.cpp` | 300 | 60 | 2000 | 0.005 | 30 | all |
| E2 | `em_vs_c_diag.cpp` | 200 | 50 | 2000 | 0.005 | 20 | kroA100, ch150, a280, d493, u574, fl1400 |
| E3 | `em_vs_c_diag.cpp` | 200 | 50 | 10000 | 0.005 | 20 | d493, u574 |
| E4 | `em_vs_c.cpp` | 500 | 100 | 20000 | 0.005 | 50 | d493, u574 |
| CEM | `em_vs_c_cem.cpp` | 300 | 60 | 2000 | 0.005 | 30 | all |
| E5 | `em_vs_c_cem.cpp` | 200 | 50 | 2000 | 0.05 | 30 | all |

Checkpoints of the diagnostic program: `{1, 10, 25, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000}` for E2 and `{1, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000}` for E3.

### E6

```bash
./localopt inputs/a280.tsp 3
```

prints `n,start,move,length` for a random tour and for the swap, 2-opt and insertion local optima of three random starts. Setting the environment variable `DUMP=<prefix>` also writes each local-optimum tour to a text file.

### E7 to E14

```bash
python scripts/run_experiments.py                 # all 480 jobs, all CPU cores but one
python scripts/run_experiments.py --workers 4     # fewer parallel jobs
python scripts/run_experiments.py --only E12,E13  # selected experiments
python scripts/run_experiments.py --inst eil51 --runs 1   # quick test
python scripts/run_experiments.py --list
```

Results go to `data/<experiment>/<instance>/`. Finished variants are skipped, so the script can be stopped and restarted. A single configuration can also be run directly, for example:

```bash
./q1ga --dir inputs --out out --inst a280 --runs 30 --xover ox --sel tour --variants C_2opt,EM_2opt
```

| Option | Values | Meaning |
|:--|:--|:--|
| `--gens`, `--gk` | integer | fixed generations, or `gk × n` generations |
| `--pop`, `--elite`, `--pm`, `--ess` | numbers | population, elite size, mutation rate, essential rate |
| `--init` | `random`, `nn`, `kmeans` | initial population |
| `--kspec`, `--fileorder` | `10`, `sqrt`, `n30`, `n10`; `0/1` | number of clusters; keep file order inside clusters |
| `--xover`, `--sel`, `--tk` | `cx`, `ox`; `roulette`, `tour`; integer | crossover, selection, tournament size |
| `--partner`, `--knn` | `uniform`, `knn`; integer | how the partner position is drawn |
| `--variants` | `EM_<move>`, `C_<move>`, `C_sched`, `C_stag`, `C_mix` | moves: `swap`, `2opt`, `rev`, `ins` |
| `--stagc`, `--stagc2`, `--sa`, `--sb` | numbers | C_stag thresholds (×1/n) and C_sched switch times (×n) |

### Run time

All experiments were run on an HP EliteBook x360 830 G6 (Intel Core i5-8265U, 4 cores / 8 threads, 8 GB RAM, Windows). E1 takes about 10 hours on one thread and E7 to E14 about 70 CPU hours in total. Some E7 to E14 runs were interrupted by the computer going to sleep, so their recorded times are not meaningful; the tours are unaffected because every run is seeded.

---

## Output files

| File | Content |
|:--|:--|
| `ALL_EMvsC_summary_<stamp>.csv`, `ALL_CEM_summary_<stamp>.csv` | one row per instance and variant: best, mean, standard deviation, convergence generation, time, Mann–Whitney p-value against the twin variant |
| `<instance>_<variant>_runs[_<stamp>].csv` | one row per run: best length, time, convergence generation (q1ga also logs the switch generations of the schedules) |
| `<instance>_<variant>_FR[_<stamp>].csv` | best length in every generation of run 1 |
| `DIAG[_<stamp>].csv` | per run and checkpoint: best-so-far length, attempted and accepted moves, total gain |
| `SUMMARY.csv` (q1ga) | one row per variant for that instance |
| `E6_local_optima/local_optima_raw.csv` | instance, n, start, move, length |
| `E6_local_optima/cross_group_share.csv` | share of the local-optimum length spent on edges between city groups |

Tour lengths use unrounded Euclidean distances and are rounded to one decimal.

---

## Reproducing the tables and figures

```bash
cd notebooks && jupyter nbconvert --to notebook --execute analysis.ipynb   # E1 to E5, Figures 1 to 12
cd ../analysis
python summarize.py      # E7 to E14 summaries and tests
python figures.py        # Figures 13 to 16
python city_groups.py    # effective number of city groups
```

The flowchart is drawn with Graphviz: `dot -Tpng figures/fig1_flowchart.dot -o figures/fig1_flowchart.png`.

---

## References

- Abu Doush, I., Al-Khatib, R. M., Alotaibi, B., & Alhatamleh, S. (2026). Essential mutation framework for solving travelling salesman problem. *Cluster Computing, 29*, 49. https://doi.org/10.1007/s10586-025-05789-w
- Rahman, M. A., & Ma, J. (2024). Two-stage probe-based search optimization algorithm for the traveling salesman problems. *Mathematics, 12*(9), 1340. https://doi.org/10.3390/math12091340
- Rahman, M. A., Nazib, K. M., Ali, L. E., & Islam, M. R. (2025). Solving traveling salesman problem through genetic algorithm with clustering. *International Journal of Intelligent Systems and Applications, 17*(3). https://doi.org/10.5815/ijisa.2025.03.02
- Reinelt, G. (1991). TSPLIB—A traveling salesman problem library. *ORSA Journal on Computing, 3*(4), 376–384.

## License

MIT, see `LICENSE`.
