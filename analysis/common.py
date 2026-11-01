import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
INPUTS = os.path.join(ROOT, 'inputs')
RESULTS = os.path.join(ROOT, 'results')
FIGURES = os.path.join(ROOT, 'figures')

BKS = dict(eil51=426, berlin52=7542, st70=675, kroA100=21282, bier127=118282, ch150=6528, pr152=73682, rat195=2323,
           kroA200=29368, pr226=80369, a280=2579, pr299=48191, lin318=42029, pr439=107217, d493=35002, u574=36905,
           rat575=6773, p654=34643, u724=41910, rat783=8806, pr1002=259045, u1060=224094, vm1084=239297,
           pcb1173=56892, d1291=50801, rl1304=252948, rl1323=270199, nrw1379=56638, fl1400=20127, u1432=152970)


def load(path):
    pts = []
    for line in open(path):
        t = line.split()
        if len(t) != 3:
            continue
        try:
            pts.append([float(x) for x in t])
        except ValueError:
            continue
    return np.array(pts)[:, 1:]


def coords(inst):
    return load(os.path.join(INPUTS, inst + '.tsp'))


INSTANCES = sorted(BKS, key=lambda i: len(coords(i)))
N = {i: len(coords(i)) for i in INSTANCES}
C3 = ['C_swap', 'C_ins', 'C_2opt']
C4 = ['C_swap', 'C_ins', 'C_2opt', 'C_rev']
GROUPED = ['pr226', 'p654', 'fl1400']
TRAIN = ['kroA100', 'ch150', 'a280', 'd493', 'u574', 'fl1400']


def runs(exp, inst, var):
    f = os.path.join(DATA, exp, inst, f'{inst}_{var}_runs.csv')
    return pd.read_csv(f).Best.values if os.path.exists(f) else None


def e1(inst, var):
    f = glob.glob(os.path.join(DATA, 'E1_benchmark', f'{inst}_{var}_runs_*.csv'))[0]
    return pd.read_csv(f).Best.values


def mw(a, b):
    return mannwhitneyu(a, b).pvalue


def gap(length, inst):
    return 100 * (length / BKS[inst] - 1)


def winner(exp, inst, variants):
    d = {v: runs(exp, inst, v) for v in variants}
    m = {v: x.mean() for v, x in d.items()}
    o = sorted(m, key=m.get)
    return o[0], 100 * (m[o[1]] / m[o[0]] - 1), mw(d[o[0]], d[o[1]]), m


def e1_winner(inst):
    m = {v: e1(inst, v).mean() for v in C4}
    return min(m, key=m.get)


def synthetic_k(name):
    if '_u' in name:
        return 999
    if name.startswith('syn_s'):
        return 20
    import re
    return int(re.search(r'k_?0?(\d+)', name).group(1))
