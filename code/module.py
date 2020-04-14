from __future__ import division
from functools import reduce
import pandas as pd, numpy as np, scipy.io as sio, multiprocessing as mp, itertools as it, networkx as nx #, igraph as ig
import os, sys, time, string as s, copy, glob

sys.path.append('/home/aa/git/lcn_fmri')
import aani

sys.path.append('/home/aa/git/aa_env/aa_python')
import aa_misc

sys.path.append('/home/aa/git/lcn_fmri/modules')
import modules_louvain

get_part_0 = aa_misc.get_part(0)
get_part_1 = aa_misc.get_part(1)
def outer_equal(a_list):
    return map(lambda x: x[0]==x[1], it.product(a_list, a_list))

print(__file__)
euclidean_dist = np.loadtxt(os.path.abspath(os.path.join(
    __file__, 
    '../../utils',
    'willardROI_euclidean-distance_asMatrix'
)))
print(euclidean_dist.shape)

filenames_RT = sorted(
    glob.glob(os.path.join(
        '/media/aa/aa_10TB/external_LARGE/LMP_MC/DATA_analysis',
        'LMP_MC_*_v*/00*/all_3mm_mean/3mm_rs_RT_30_30.mat')))
print(len(filenames_RT))
print(filenames_RT)

def filenames_RT2louvain(filename_RT):
    if filename_RT[-4:] == '.mat':
        return filename_RT.replace('.mat', '_louvain-plus1')
    return filename_RT + '_louvain-plus1'
    
filenames_louvain = map(lambda filename_RT: filenames_RT2louvain(filename_RT), filenames_RT)
print(np.unique(list(map(
    os.path.isfile, 
    filenames_louvain)), return_counts=True))
# list(map(lambda pair: map(os.path.isfile, pair), 
for pair in zip(filenames_RT, filenames_louvain):
    print([os.path.isfile(filename) for filename in pair])
    
for filename_RT, filename_louvain in zip(filenames_RT, filenames_louvain):
    # filename_RT, filename_louvain = map(get_part_0, [filenames_RT, filenames_louvain])
    filename_louvain_instab_grouped = filename_louvain + '_instab-grouped-byEuclidean'
    print filename_louvain_instab_grouped
    RT = aani.load_snapshots(filename_RT)
    RT_flat = map(lambda RT_snapshot: RT_snapshot.flatten(), RT)
    louvain = np.loadtxt(filename_louvain)
    n_nodes = len(louvain[snapshot_prev])
    print map(np.shape, [RT, louvain])
    # diff_by_moduleSwitch_all = []
    df_snapshots = []
    for snapshot_prev, snapshot_next in zip(range(len(RT_flat)-1),np.arange(1,len(RT_flat))):
        outer_noSwitchAndDistance_flat = zip(
            # np.repeat(snapshot_prev, len(louvain[snapshot_prev])**2),
            outer_equal(louvain[snapshot_prev]),
            outer_equal(louvain[snapshot_next]),
            np.around(map(
                lambda pair: euclidean_dist[pair[0],pair[1]], 
                it.product(range(n_nodes),range(n_nodes))
            ), -1)
        )
        indices_grouped_by_parts = get_indices_grouped_by_parts(outer_noSwitchAndDistance_flat, [0,1,2])
        mean_instab_grouped = map(
            lambda elem: (
                [snapshot_next] + elem[0] + [np.mean(RT_flat[snapshot_next][elem[1]] - RT_flat[snapshot_prev][elem[1]])]
            ), 
            indices_grouped_by_parts
        )
        df_snapshot = pd.DataFrame.from_records(mean_instab_grouped, index=[0,1,2,3])
        df_snapshots.append(df_snapshot)
    df = pd.concat(df_snapshots)
    df.to_csv(filename_louvain_instab_grouped)    
