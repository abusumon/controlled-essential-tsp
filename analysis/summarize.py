import os

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from common import (BKS, C3, C4, DATA, GROUPED, INSTANCES, N, RESULTS, TRAIN, e1, gap, mw, runs,
                    synthetic_k, winner)

MOVES = ['swap', '2opt', 'rev', 'ins']


def wtl(better, worse, p):
    w = sum(1 for b, x, q in zip(better, worse, p) if b < x and q < .05)
    l = sum(1 for b, x, q in zip(better, worse, p) if b > x and q < .05)
    return w, len(p) - w - l, l


def e7_budget():
    print('E7  budget proportional to n: winning move')
    rows = []
    for exp in ['E7_budget_2n', 'E7_budget_7n']:
        for i in INSTANCES:
            w, marg, p, m = winner(exp, i, C3)
            rows.append(dict(experiment=exp, instance=i, n=N[i], winner=w, margin_pct=marg, p=p,
                             **{f'gap_{v}': gap(m[v], i) for v in C3}))
    d = pd.DataFrame(rows)
    for exp, g in d.groupby('experiment'):
        print(' ', exp, g.winner.value_counts().to_dict(), 'significant:', int((g.p < .05).sum()))
    d.to_csv(os.path.join(RESULTS, 'E7_budget_summary.csv'), index=False)


def e8_nn_start():
    print('E8  nearest-neighbour start')
    rows = []
    for i in INSTANCES:
        w, marg, p, m = winner('E8_nn_start', i, C4)
        rank_swap = sorted(m.values()).index(m['C_swap']) + 1
        rows.append(dict(instance=i, n=N[i], winner=w, margin_pct=marg, p=p, swap_rank=rank_swap,
                         **{f'gap_{v}': gap(m[v], i) for v in C4}))
    d = pd.DataFrame(rows)
    print(' ', d.winner.value_counts().to_dict(), 'swap last on', int((d.swap_rank == 4).sum()))
    d.to_csv(os.path.join(RESULTS, 'E8_nn_start_summary.csv'), index=False)


def e9_synthetic():
    print('E9  synthetic instances with k city groups')
    rows = []
    for i in sorted(os.listdir(os.path.join(DATA, 'E9_synthetic'))):
        d = {v: runs('E9_synthetic', i, v) for v in C3}
        m = {v: x.mean() for v, x in d.items()}
        o = sorted(m, key=m.get)
        n = len(open(os.path.join(DATA, '..', 'inputs', 'synthetic', i + '.tsp')).read()
                .split('NODE_COORD_SECTION')[1].split('\n')) - 3
        rows.append(dict(instance=i, n=n, k=synthetic_k(i), ratio_2opt_swap=m['C_2opt'] / m['C_swap'],
                         winner=o[0], p=mw(d[o[0]], d[o[1]]), **{f'mean_{v}': m[v] for v in C3}))
    d = pd.DataFrame(rows).sort_values(['k', 'instance'])
    few, many = d[d.k <= 12], d[(d.k >= 20) & (d.k < 999)]
    print('  k<=12: 2-opt wins', int((few.winner == 'C_2opt').sum()), 'of', len(few),
          '| k>=20: swap wins', int((many.winner == 'C_swap').sum()), 'of', len(many))
    print('  mean 2-opt/swap ratio by k:', d.groupby('k').ratio_2opt_swap.mean().round(3).to_dict())
    d.to_csv(os.path.join(RESULTS, 'E9_synthetic_summary.csv'), index=False)


def start_gap(exp, inst):
    d = pd.read_csv(os.path.join(DATA, exp, inst, 'DIAG.csv'), header=None)
    d = d[d[5] == 1]
    return gap(d[6].mean(), inst)


def e10_clustered():
    print('E10 clustered initial populations')
    rows = []
    for exp in ['E10_clustered_10', 'E10_clustered_sqrt', 'E10_clustered_n30', 'E10_clustered_n10',
                'E10_clustered_10_fileorder']:
        wins, gaps, starts, better = [], [], [], 0
        for i in INSTANCES:
            w, _, _, m = winner(exp, i, C3)
            wins.append(w)
            best = min(m.values())
            gaps.append(gap(best, i))
            starts.append(start_gap(exp, i))
            better += best < min(e1(i, v).mean() for v in C4)
        vc = pd.Series(wins).value_counts()
        rows.append(dict(setting=exp, swap=vc.get('C_swap', 0), ins=vc.get('C_ins', 0), opt2=vc.get('C_2opt', 0),
                         median_start_gap=np.median(starts), median_best_gap=np.median(gaps), better_than_E1=better))
    d = pd.DataFrame(rows)
    print(d.round(1).to_string(index=False))
    d.to_csv(os.path.join(RESULTS, 'E10_clustered_summary.csv'), index=False)
    ratio = [min(runs('E10_clustered_10', i, v).mean() for v in C3) /
             min(runs('E10_clustered_10_fileorder', i, v).mean() for v in C3) for i in INSTANCES]
    print('  file order better on', sum(r > 1 for r in ratio), 'instances; median shuffled/file ratio',
          round(float(np.median(ratio)), 2))


def e11_clustered_em():
    print('E11 EM from clustered starts (k = sqrt n) against controlled')
    for m in ['swap', 'ins', '2opt']:
        em = [runs('E11_clustered_em', i, 'EM_' + m) for i in INSTANCES]
        c = [runs('E10_clustered_sqrt', i, 'C_' + m) for i in INSTANCES]
        w, t, l = wtl([x.mean() for x in c], [x.mean() for x in em], [mw(a, b) for a, b in zip(em, c)])
        ratio = np.median([a.mean() / b.mean() for a, b in zip(em, c)])
        print(f'  {m:5s} controlled W/T/L {w}/{t}/{l}  median EM/C {ratio:.2f}')


def e12_ox_tournament():
    print('E12 order crossover and tournament selection')
    rows = []
    for i in INSTANCES:
        r = dict(instance=i, n=N[i])
        for m in MOVES:
            em, c = runs('E12_ox_tournament', i, 'EM_' + m), runs('E12_ox_tournament', i, 'C_' + m)
            r['EM_' + m], r['C_' + m], r['p_' + m] = em.mean(), c.mean(), mw(em, c)
        for v in ['C_stag', 'C_mix']:
            r[v] = runs('E12_ox_tournament', i, v).mean()
        rows.append(r)
    a = pd.DataFrame(rows)
    for m in MOVES:
        w, t, l = wtl(a['C_' + m], a['EM_' + m], a['p_' + m])
        ratio = a['EM_' + m] / a['C_' + m]
        print(f'  {m:5s} controlled W/T/L {w}/{t}/{l}  median EM/C {ratio.median():.2f}  '
              f'signed-rank p {wilcoxon(a["EM_" + m], a["C_" + m]).pvalue:.1e}')
    single = a[['C_swap', 'C_2opt', 'C_rev', 'C_ins']]
    a['best_single'] = single.min(axis=1)
    a['winner'] = single.idxmin(axis=1)
    a['gap_best'] = [gap(b, i) for b, i in zip(a.best_single, a.instance)]
    print('  winners:', a.winner.value_counts().to_dict())
    a.to_csv(os.path.join(RESULTS, 'E12_ox_tournament_summary.csv'), index=False)


def e13_schedules():
    print('E13 move schedules against the best single move at the same budget')
    rows = []
    for i in INSTANCES:
        singles = {'2000': {v: e1(i, v) for v in C4},
                   '2n': {v: runs('E7_budget_2n', i, v) for v in C3},
                   '7n': {v: runs('E7_budget_7n', i, v) for v in C3}}
        for b in ['2000', '2n', '7n']:
            sg = singles[b]
            bv = min(sg, key=lambda v: sg[v].mean())
            best = sg[bv]
            for v in ['C_stag', 'C_sched', 'C_mix']:
                x = runs('E13_schedules_' + b, i, v)
                p = mw(x, best)
                res = 'tie' if p >= .05 else ('win' if x.mean() < best.mean() else 'loss')
                rows.append(dict(instance=i, n=N[i], budget=b, schedule=v, best_single=bv,
                                 rel_pct=100 * (x.mean() / best.mean() - 1), p=p, result=res))
    g = pd.DataFrame(rows)
    for b in ['2n', '2000', '7n']:
        for v in ['C_stag', 'C_sched', 'C_mix']:
            s = g[(g.budget == b) & (g.schedule == v)]
            t = s[~s.instance.isin(TRAIN)]
            cnt = lambda x: '/'.join(str(int((x.result == r).sum())) for r in ['win', 'tie', 'loss'])
            print(f'  {b:5s} {v:8s} all {cnt(s):8s} test {cnt(t):8s} median {s.rel_pct.median():+.1f}%')
    g.to_csv(os.path.join(RESULTS, 'E13_schedules_summary.csv'), index=False)


def e14_knn_partner():
    print('E14 partner drawn from the 10 nearest cities')
    rows = []
    for i in INSTANCES:
        r = dict(instance=i, n=N[i])
        for b in ['2000', '2n']:
            r['winner_' + b] = winner('E14_knn_partner_' + b, i, C3)[0]
        knn = min(runs('E14_knn_partner_2000', i, v).mean() for v in C3)
        uni = min(e1(i, v).mean() for v in C4)
        r.update(knn_best_gap=gap(knn, i), uniform_best_gap=gap(uni, i),
                 stag_vs_knn_best_pct=100 * (runs('E14_knn_partner_2000', i, 'C_stag').mean() / knn - 1))
        rows.append(r)
    d = pd.DataFrame(rows)
    print('  2000 gen.:', d.winner_2000.value_counts().to_dict(), '| 2n:', d.winner_2n.value_counts().to_dict())
    d.to_csv(os.path.join(RESULTS, 'E14_knn_partner_summary.csv'), index=False)


if __name__ == '__main__':
    os.makedirs(RESULTS, exist_ok=True)
    for f in [e7_budget, e8_nn_start, e9_synthetic, e10_clustered, e11_clustered_em, e12_ox_tournament,
              e13_schedules, e14_knn_partner]:
        f()
