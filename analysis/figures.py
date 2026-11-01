import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from common import C3, C4, DATA, FIGURES, GROUPED, INPUTS, INSTANCES, N, RESULTS, e1_winner, runs, synthetic_k, winner

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': '#e6e6e6', 'grid.linewidth': 0.7, 'axes.axisbelow': True})
COL = {'swap': '#D55E00', 'ins': '#0072B2', '2opt': '#009E73', 'rev': '#CC79A7'}
MARK = {'swap': 'o', 'ins': 's', '2opt': '^', 'rev': 'D'}
LABEL = {'swap': 'Swap', 'ins': 'Insertion', '2opt': '2-opt', 'rev': 'Reversion'}


def winners(exp, variants):
    return [winner(exp, i, variants)[0][2:] for i in INSTANCES]


def winner_strips():
    e1w = [e1_winner(i)[2:] for i in INSTANCES]
    panels = [
        ('(a) Uniform partner, random start, three budgets',
         [('$G = 2n$', winners('E7_budget_2n', C3)), ('$G = 2000$ (E1)', e1w), ('$G = 7n$', winners('E7_budget_7n', C3))]),
        ('(b) Partner from the 10 nearest cities',
         [('$G = 2n$', winners('E14_knn_partner_2n', C3)), ('$G = 2000$', winners('E14_knn_partner_2000', C3))]),
        ('(c) Order crossover and tournament selection, 2000 generations',
         [('OX + tournament', winners('E12_ox_tournament', C4))]),
        ('(d) 2000 generations, different initial populations',
         [('Random (E1)', e1w)] +
         [(lab, winners(exp, C3)) for exp, lab in [('E10_clustered_10', 'Clustered, $k = 10$'),
                                                   ('E10_clustered_sqrt', 'Clustered, $k = \\sqrt{n}$'),
                                                   ('E10_clustered_n30', 'Clustered, $k = n/30$'),
                                                   ('E10_clustered_n10', 'Clustered, $k = n/10$')]] +
         [('Nearest neighbour', winners('E8_nn_start', C4))]),
    ]
    fig, ax = plt.subplots(4, 1, figsize=(10, 8.4), gridspec_kw={'height_ratios': [3, 2, 1, 6]})
    for a, (title, rows) in zip(ax, panels):
        for r, (lab, w) in enumerate(rows):
            y = len(rows) - 1 - r
            for x, m in enumerate(w):
                a.plot(x, y, MARK[m], color=COL[m], ms=6.5, mec='white', mew=0.6)
        a.set_yticks(range(len(rows)))
        a.set_yticklabels([l for l, _ in rows][::-1])
        a.set_xticks(range(len(INSTANCES)))
        a.set_xticklabels(INSTANCES if a is ax[3] else [''] * len(INSTANCES), rotation=90, fontsize=7)
        a.set_xlim(-0.6, len(INSTANCES) - 0.4)
        a.set_ylim(-0.6, len(rows) - 0.4)
        a.grid(axis='y', visible=False)
        a.set_title(title, fontsize=10, loc='left')
    for i in GROUPED:
        ax[3].get_xticklabels()[INSTANCES.index(i)].set_fontweight('bold')
    handles = [Line2D([], [], marker=MARK[m], color=COL[m], ls='none', ms=7, label=LABEL[m]) for m in COL]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(FIGURES, 'fig14_winner_strips.png'), dpi=300)
    plt.close()


def city_groups():
    rows = []
    for i in sorted(os.listdir(os.path.join(DATA, 'E9_synthetic'))):
        m = {v: runs('E9_synthetic', i, v).mean() for v in C3}
        n = len(open(os.path.join(INPUTS, 'synthetic', i + '.tsp')).read().split('NODE_COORD_SECTION')[1].split('\n')) - 3
        rows.append([i, n, synthetic_k(i), m['C_2opt'] / m['C_swap']])
    d = pd.DataFrame(rows, columns=['i', 'n', 'k', 'ratio'])
    x = pd.read_csv(os.path.join(DATA, 'E6_local_optima', 'cross_group_share.csv')).set_index('inst')
    x = x.loc[['synk_05_0', 'syn_s05_0', 'pr226', 'p654', 'fl1400']]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
    a = ax[0]
    s = d[d.n == 1000]
    kx = s.k.where(s.k < 999, 150)
    tight = ~s.i.str.startswith(('syn_s10', 'syn_s20'))
    a.plot(kx[tight], s.ratio[tight], 'o', color='#0072B2', ms=6, label='$n = 1000$, tight groups')
    a.plot(kx[~tight], s.ratio[~tight], 'o', mfc='white', color='#0072B2', ms=6, label='$n = 1000$, looser groups ($k = 20$)')
    for n, mk in [(300, 'v'), (2000, 'P')]:
        t = d[d.n == n]
        a.plot(t.k.where(t.k < 999, 150), t.ratio, mk, color='#555555', ms=6, label=f'$n = {n}$')
    a.axhline(1, color='#888888', lw=1, ls='--')
    a.set_xscale('log')
    a.set_xticks([3, 5, 8, 12, 20, 30, 50, 150])
    a.set_xticklabels(['3', '5', '8', '12', '20', '30', '50', 'uniform'])
    a.set_xlabel('Number of city groups $k$')
    a.set_ylabel('Mean length, 2-opt / swap')
    a.text(60, 0.97, '2-opt better', fontsize=8, va='top', ha='center', color='#444')
    a.text(60, 1.03, 'swap better', fontsize=8, va='bottom', ha='center', color='#444')
    a.set_title('(a) Synthetic instances, 2000 generations, 30 runs', fontsize=10)
    a.legend(frameon=False, fontsize=7.5, loc='lower right')
    a = ax[1]
    names = {'synk_05_0': 'synthetic\n5 groups', 'syn_s05_0': 'synthetic\n20 groups', 'pr226': 'pr226', 'p654': 'p654',
             'fl1400': 'fl1400'}
    w = 0.26
    xs = np.arange(len(x))
    for j, (op, col) in enumerate([('swap', 'swap: share of length'), ('ins', 'ins: share'), ('2opt', '2opt: share')]):
        a.bar(xs + (j - 1) * w, 100 * x[col], w, color=COL[op], label=LABEL[op])
    a.set_xticks(xs)
    a.set_xticklabels([names[i] for i in x.index], fontsize=8.5)
    a.set_ylabel('Tour length on cross-group edges (%)')
    a.set_ylim(0, 118)
    a.set_yticks([0, 20, 40, 60, 80, 100])
    a.set_title('(b) Local optima: length spent between groups', fontsize=10)
    a.legend(frameon=False, fontsize=8.5, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.02))
    a.grid(axis='x', visible=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig15_groups.png'), dpi=300)
    plt.close()


def schedules():
    g = pd.read_csv(os.path.join(RESULTS, 'E13_schedules_summary.csv'))
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.6), sharey=True)
    for a, b, t in zip(ax, ['2n', '2000', '7n'], ['$G = 2n$', '$G = 2000$', '$G = 7n$']):
        for v, col, mk, lab in [('C_mix', '#999999', 'x', 'Random mix'), ('C_sched', '#E69F00', 's', 'Fixed schedule'),
                                ('C_stag', '#0072B2', 'o', 'Acceptance-triggered (C_stag)')]:
            s = g[(g.budget.astype(str) == b) & (g.schedule == v)].set_index('instance').reindex(INSTANCES)
            a.plot([N[i] for i in INSTANCES], s.rel_pct.clip(-30, 60), mk, color=col, ms=5 if mk != 'x' else 6,
                   label=lab, ls='none')
        a.axhline(0, color='#555', lw=1)
        a.set_xscale('log')
        a.set_title(t, fontsize=10)
        a.set_xlabel('Number of cities $n$')
    ax[0].set_ylabel('Mean length vs best single move (%)')
    ax[0].set_ylim(-32, 62)
    ax[2].text(1400, 58, 'fl1400, p654 off scale', fontsize=7, ha='right', va='top', color='#444')
    ax[0].legend(frameon=False, fontsize=7.5, loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig16_schedule.png'), dpi=300)
    plt.close()


def local_optima():
    raw = pd.read_csv(os.path.join(DATA, 'E6_local_optima', 'local_optima_raw.csv'))
    return raw.groupby(['instance', 'move']).length.mean().unstack()


def depth_speed():
    lo = local_optima()
    from common import BKS
    m = lo.loc[[i for i in INSTANCES]]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    a = ax[0]
    for op in ['swap', 'ins', '2opt']:
        g = 100 * (m[op] / [BKS[i] for i in m.index] - 1)
        a.plot([N[i] for i in m.index], g, MARK[op], color=COL[op], ms=5.5, label=LABEL[op], ls='none')
        if op == 'swap':
            for i in ['fl1400', 'p654']:
                a.annotate(i, (N[i], g[i]), xytext=(-7, 0), textcoords='offset points', fontsize=8, ha='right', va='center')
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlabel('Number of cities $n$')
    a.set_ylabel('Gap of the local optimum (%)')
    a.set_title('(a) Depth: local optimum of each move', fontsize=10)
    a.legend(frameon=False, fontsize=9, loc='upper left')
    diag = pd.read_csv([os.path.join(DATA, 'E2_trajectories', f) for f in os.listdir(os.path.join(DATA, 'E2_trajectories'))
                        if f.startswith('DIAG')][0])
    curve = diag.groupby(['Instance', 'Variant', 'Gen']).BestSoFar.mean().reset_index()
    curve = curve[curve.Gen > 1]
    a = ax[1]
    for v, op in {'C_swap': 'swap', 'C_ins': 'ins', 'C_2opt': '2opt'}.items():
        for k, (inst, g) in enumerate(curve[curve.Variant == v].groupby('Instance')):
            f, r = lo.loc[inst, op], lo.loc[inst, 'random']
            e = ((g.BestSoFar - f) / (r - f)).clip(lower=1e-3)
            a.plot(g.Gen / N[inst], e, '-', marker=MARK[op], ms=3, color=COL[op], lw=1.2, label=LABEL[op] if k == 0 else None)
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_ylim(5e-3, 1.2)
    a.set_xlabel('Generations per city $g/n$')
    a.set_ylabel('Remaining distance to own local optimum')
    a.set_title('(b) Speed: approach to that optimum (E2)', fontsize=10)
    a.legend(frameon=False, fontsize=9, loc='lower left')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig13_depth_speed.png'), dpi=300)
    plt.close()


if __name__ == '__main__':
    os.makedirs(FIGURES, exist_ok=True)
    depth_speed()
    winner_strips()
    city_groups()
    schedules()
