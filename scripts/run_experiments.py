import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, 'q1ga.exe' if os.name == 'nt' else 'q1ga')
TSPLIB = os.path.join(ROOT, 'inputs')
SYNTHETIC = os.path.join(ROOT, 'inputs', 'synthetic')
C3 = 'C_swap,C_ins,C_2opt'


def names(folder):
    return sorted(f[:-4] for f in os.listdir(folder) if f.endswith('.tsp'))


def size(folder, inst):
    with open(os.path.join(folder, inst + '.tsp')) as f:
        return sum(1 for line in f if len(line.split()) == 3 and line.split()[0].replace('.', '').isdigit())


def jobs():
    T, S = names(TSPLIB), names(SYNTHETIC)
    J = []

    def add(exp, prio, inst, folder, args):
        J.append(dict(exp=exp, prio=prio, inst=inst, folder=folder, args=args))

    for i in T:
        add('E13_schedules_2000', 1, i, TSPLIB, ['--variants', 'C_stag,C_sched,C_mix'])
        add('E13_schedules_2n', 1, i, TSPLIB, ['--gk', '2', '--variants', 'C_stag,C_sched,C_mix'])
        add('E13_schedules_7n', 5, i, TSPLIB, ['--gk', '7', '--variants', 'C_stag,C_sched,C_mix'])
        add('E14_knn_partner_2000', 2, i, TSPLIB, ['--partner', 'knn', '--knn', '10', '--variants', C3 + ',C_stag'])
        add('E14_knn_partner_2n', 2, i, TSPLIB, ['--partner', 'knn', '--knn', '10', '--gk', '2', '--variants', C3 + ',C_stag'])
        add('E12_ox_tournament', 3, i, TSPLIB, ['--xover', 'ox', '--sel', 'tour', '--variants',
                                               'EM_swap,C_swap,EM_2opt,C_2opt,EM_rev,C_rev,EM_ins,C_ins,C_stag,C_mix'])
        add('E11_clustered_em', 4, i, TSPLIB, ['--init', 'kmeans', '--kspec', 'sqrt', '--variants', 'EM_swap,EM_ins,EM_2opt'])
        add('E7_budget_2n', 6, i, TSPLIB, ['--gk', '2', '--variants', C3])
        add('E7_budget_7n', 7, i, TSPLIB, ['--gk', '7', '--variants', C3])
        add('E8_nn_start', 6, i, TSPLIB, ['--init', 'nn', '--variants', 'C_swap,C_2opt,C_rev,C_ins'])
        for k in ['sqrt', 'n30', '10', 'n10']:
            add('E10_clustered_' + k, 9, i, TSPLIB, ['--init', 'kmeans', '--kspec', k, '--variants', C3])
        add('E10_clustered_10_fileorder', 9, i, TSPLIB, ['--init', 'kmeans', '--kspec', '10', '--fileorder', '1', '--variants', C3])
    for i in S:
        add('E9_synthetic', 8, i, SYNTHETIC, ['--variants', C3])
    for j in J:
        j['n'] = size(j['folder'], j['inst'])
        j['out'] = os.path.join(ROOT, 'data', j['exp'], j['inst'])
    J.sort(key=lambda j: (j['prio'], -j['n']))
    return J


def run(j, n_runs):
    os.makedirs(j['out'], exist_ok=True)
    cmd = [EXE, '--dir', j['folder'], '--out', j['out'], '--inst', j['inst'], '--runs', str(n_runs)] + j['args']
    t0 = time.time()
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    return j, r.returncode, r.stderr, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description='Run experiments E7 to E14 with q1ga.')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--only', default='', help='comma-separated experiment prefixes, e.g. E12,E13')
    ap.add_argument('--inst', default='', help='comma-separated instance names')
    ap.add_argument('--runs', type=int, default=30)
    ap.add_argument('--list', action='store_true')
    a = ap.parse_args()
    J = jobs()
    if a.inst:
        J = [j for j in J if j['inst'] in a.inst.split(',')]
    if a.only:
        J = [j for j in J if j['exp'].startswith(tuple(x.strip() for x in a.only.split(',')))]
    if a.list:
        for j in J:
            print(f"{j['exp']:28s} {j['inst']}")
        return
    if not os.path.exists(EXE):
        sys.exit('Compile first: g++ -O2 -std=c++17 src/q1ga.cpp -o q1ga.exe')
    print(f'{len(J)} jobs, {a.workers} in parallel. Finished variants are skipped on restart.')
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futures = [ex.submit(run, j, a.runs) for j in J]
        for k, f in enumerate(as_completed(futures), 1):
            j, rc, err, dt = f.result()
            status = 'ok' if rc == 0 else f'FAILED ({rc}): {err.strip()[:200]}'
            print(f"[{k}/{len(J)}] {j['exp']:28s} {j['inst']:10s} {status} {dt / 60:.1f} min "
                  f"(elapsed {(time.time() - t0) / 3600:.2f} h)", flush=True)


if __name__ == '__main__':
    main()
