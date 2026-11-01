import glob
import os

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform

from common import INPUTS, INSTANCES, N, RESULTS, coords, load


def k_eff(xy, factor=10):
    tree = minimum_spanning_tree(squareform(pdist(xy))).toarray()
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    tree[tree > factor * np.median(nn)] = 0
    _, label = connected_components(tree, directed=False)
    p = np.bincount(label) / len(label)
    return 1 / (p ** 2).sum()


if __name__ == '__main__':
    rows = [dict(instance=i, n=N[i], k_eff=k_eff(coords(i))) for i in INSTANCES]
    for f in sorted(glob.glob(os.path.join(INPUTS, 'synthetic', '*.tsp'))):
        xy = load(f)
        rows.append(dict(instance=os.path.basename(f)[:-4], n=len(xy), k_eff=k_eff(xy)))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(RESULTS, 'city_groups.csv'), index=False)
    print(d.round(1).to_string(index=False))
