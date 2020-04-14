# -*- coding: utf-8 -*-
"""
This script takes pre-processed ROI space time-series and calculates brain
network instabilities for whole brain and for subnetworks shown below:

Subnetwork labels:
whole: whole brain
unlabeled: unlabeled
Motor: motor
Auditory: auditory
Basal: basal ganglia
dDMN: dorsal default mode network
high: high visual
Language: language
LECN: left executive control network
post: posterior salience
Precuneus: precuneus
prim: primary visual
RECN: right executive control network
Salience: salience
vDMN: ventral default mode networ
Visuospatial: visuospatial

"""

import os
import sys
import numpy as np
import pandas as pd
import itertools
from nilearn.connectome import ConnectivityMeasure
import networkx as nx, community
import itertools as it

sys.path.append('/home/aa/git/aa_env/aa_python')
import aa_misc


# =============================================================================
# Setup
# =============================================================================

# Define filepaths
HOMEDIR = os.path.abspath(os.path.join(__file__, "../..")) + "/"
WORKDIR = HOMEDIR + "data/time_series/"
OUTDIR = HOMEDIR + "data/instabilities/"

# Important inputs
TOTAL_LENGTH = 720  # Total length of the timeseries (number of frames)
WINDOW_LENGTH = 30  # Length of one time window (number of frames)
NUM_ROI_TOTAL = 498  # Number of ROIs

# Get run identifiers
subjects = [int(subid) for subid in sys.argv[1:]]  # Subject IDs from bash script
sessions = ["BHB", "GLC"]
tasks = ["rest"]
runs = ["1", "2"]

# Items to be analyzed
items = list(itertools.product(subjects, sessions, tasks, runs))

# Preallocate
merged_csv_data = []

# =============================================================================
# Get annotated subregions
# =============================================================================

func_label_data = pd.read_csv(HOMEDIR + \
                              "utils/functional_anatomical_willard.csv")
func_labels = func_label_data["network"].unique()  # Get func region labels
func_subnet_ixs = {}  # Dictionary for storing ROI indexes of each func region

# First create item for whole brain
func_labels = np.insert(func_labels, 0, "whole")
ids = np.arange(NUM_ROI_TOTAL)
func_subnet_ixs["whole"] = np.array(list(itertools.product(ids, ids))).T

# Then extract indexes and take their product for each func region label
for label in func_labels[1:]:
    ids = func_label_data.loc[func_label_data["network"] == label, "index"] \
        .to_numpy()
    func_subnet_ixs[f"{label}"] = \
        np.array(list(itertools.product(ids, ids))).T

#
#
#

euclidean_dist = np.loadtxt(os.path.abspath(os.path.join(
    __file__, 
    '../../utils',
    'willardROI_euclidean-distance_asMatrix'
)))
print("euclidean_dist.shape=", euclidean_dist.shape)


# =============================================================================
# Read in time-series and calculate instabilities
# =============================================================================

# Function for calculating instability values for a run which is already
# preprocessed and converted into ROI space

def outer_equal(a_list):
    return [x[0]==x[1] for x in it.product(a_list, a_list)]


def get_indices_grouped_by_parts(rectangle, parts):
    f_get_part_0 = aa_misc.get_part(0)
    '''
    f_get_parts = map(aa_misc.get_part, parts)
    print f_get_parts
    '''
    return [
        (group[0], [f_get_part_0(elem) for elem in group[1]])
        for group 
        in it.groupby(
            sorted(
                enumerate(rectangle), 
                key=lambda element: [element[1][part] for part in parts]
            ), 
            key=lambda element: [element[1][part] for part in parts]
        )]
    

def calculate_instabilities(item):

    # Load time-series into the pipeline
    roi_time_series = pd.read_csv(WORKDIR + "sub-{0:0>3}_ses-{1}_task-{2}_run-" \
                    "{3}.csv".format(item[0], item[1].lower(), item[2], item[3]),
                                     header=None, index_col=False)

    # Convert into numpy array
    roi_time_series = np.array(roi_time_series)

    # Check length
    assert roi_time_series.shape == (TOTAL_LENGTH, NUM_ROI_TOTAL)

    # Divide timeseries up into time windows
    # --------------------------------------

    # Number of windows-1
    num = int((roi_time_series.shape[0] - 1 - \
               (roi_time_series.shape[0] - 1) % WINDOW_LENGTH) /
              WINDOW_LENGTH)
    # Indicies at which slicing should occur
    indicies = np.linspace(WINDOW_LENGTH, num*WINDOW_LENGTH, num).astype(int)

     # 30 frame long time windows (last window might be shorter)
    windows = np.split(roi_time_series, indicies, axis=0)

    # Remove last window if shorter than window_size
    windows = windows[:-1] if roi_time_series.shape[0] % WINDOW_LENGTH > 0 \
        else windows

    # Calculate correlations
    # ----------------------

    # Take pairwise correlation of ROI tracers separately within all windows
    correlation_measure = ConnectivityMeasure(kind='correlation')
    corr_matrices = correlation_measure.fit_transform(windows)
    corr_matrices_flat = [corr_matrix.flatten() for corr_matrix in corr_matrices]
    corr_matrices_shape = corr_matrices.shape

    print(corr_matrices_shape)
    print([corr_matrix_flat.shape for corr_matrix_flat in corr_matrices_flat])

    partitions = np.zeros(corr_matrices_shape[:2])
    print(partitions.shape)
    for i, corr_matrix in enumerate(corr_matrices):
        # + 1. to eliminate negative numbers in the matrix, which cause problems
        g = nx.Graph(corr_matrix + 1.)
        communities = community.best_partition(g)
        # print(communities)
        partition = np.array(list(communities.values()))
        # print(partition)
        partitions[i] = partition

    print(partitions)
    df_partitions = pd.DataFrame(partitions)
    df_partitions.to_csv(OUTDIR + "module_instabilities_" + str(item[0]) + ".csv")
    

    df_snapshots = []
    for snapshot_prev, snapshot_next in zip(
        range(len(corr_matrices_flat)-1),
        np.arange(1, len(corr_matrices_flat))
    ):
        outer_noSwitchAndDistance_flat = zip(
            # np.repeat(snapshot_prev, len(louvain[snapshot_prev])**2),
            outer_equal(partitions[snapshot_prev]),
            outer_equal(partitions[snapshot_next]),
            np.around([
                euclidean_dist[pair[0],pair[1]] 
                for pair 
                in it.product(range(498), range(498))
            ], -1)
        )
        # TODO find out where the def of get_indices_grouped_by_parts went
        indices_grouped_by_parts = get_indices_grouped_by_parts(
            outer_noSwitchAndDistance_flat, 
            [0,1,2])

        mean_instab_grouped = [
            (
                [snapshot_next] + elem[0] + [np.mean(corr_matrices_flat[snapshot_next][elem[1]] - corr_matrices_flat[snapshot_prev][elem[1]])]
            ) 
            for elem 
            in indices_grouped_by_parts
        ]
        df_snapshot = pd.DataFrame.from_records(mean_instab_grouped, index=[0,1,2,3])
        df_snapshots.append(df_snapshot)
    df_snapshots_all = pd.concat(df_snapshots)
    print(df_snapshots_all)
    df_snapshots_all.to_csv(WORKDIR + "sub-{0:0>3}_ses-{1}_task-{2}_run-" \
        "{3}_louvain_values.csv".format(item[0], item[1].lower(), item[2], item[3]))
        
        
        

    """
    # Calculate instabilities
    # -----------------------
    # Calculate instabilities for every tau in [1, len(corr_matrices)] interval
    # method: elementwise difference of correlations of two time windows
    # separated by tau, then l2norm of the differences

    tau_vals = np.arange(1, len(corr_matrices))  # Range of tau values
    instabilities_all = [None] * len(tau_vals)  # Array for storing
    # instability values

    # Loop through all taus
    for i, tau in enumerate(tau_vals):

        # Pairs of windows which the elementwise difference of correlations
        # is taken between
        window_pairs = [[start, start + tau] \
                   for start in range(len(corr_matrices) - tau)]

        # Calculate and concatenate elementwise differences
        diffs = np.concatenate(
            [np.diff((corr_matrices[win1], corr_matrices[win2]), axis=0) \
             for win1, win2 in window_pairs])

        # Preallocate for storing, filled with nans
        subnets_instab = np.full((len(func_subnet_ixs), len(tau_vals)), np.nan)

        # Loop through all subnetworks
        for j, ixs in enumerate(func_subnet_ixs.values()):

            # Extract diffs corresponding to (sub)network
            subnet_diffs = diffs[:, ixs[0], ixs[1]]

            # Take l2norm and normalize with number of nodes - 1(for diagonal)
            subnets_instab[j][:len(diffs)] = np.linalg.norm(subnet_diffs, axis=1) \
                /np.sqrt(len(ixs.T)-1)

        # Store instabilities
        instabilities_all[i] = list(subnets_instab)

    # Make csv data
    # -------------

    # Labels
    dflabels_A = pd.DataFrame(
        [
        ["sub{0:0>3}".format(item[0]), item[1], item[2], item[3]]
                    for tau in tau_vals
                    for label in func_labels
                    for time in tau_vals
         ]
                                )

    dflabels_B = pd.DataFrame(itertools.product(tau_vals, func_labels, tau_vals))
    dflabels = pd.concat((dflabels_A, dflabels_B), axis=1)
    mtx = pd.MultiIndex.from_frame(
        dflabels,
        names=["subject", "bolus", "task", "run", "tau", "funclabel", "time"]
        )

    # Values
    csv_data = pd \
        .DataFrame(instabilities_all) \
        .T.melt() \
        .drop(labels="variable", axis=1) \
        .explode(column="value") \
        .set_index(mtx) \
        .unstack(["funclabel"]) \
        .xs('value', axis=1, drop_level=True) \
        [func_labels]

    return csv_data
    """
    return 0

# Loop through all items
for item in items:
    print("Computing instabilities for:", item)

    # Call the function to calculate intabilities for current run
    calculated_instabs = calculate_instabilities(item)

    # Add calculated instabilities to output data
    merged_csv_data.append(calculated_instabs)

# =============================================================================
# Write results into csv
# =============================================================================
"""
pd \
    .concat(merged_csv_data, axis=0) \
    .to_csv(OUTDIR + f"instabilities_N{len(subjects)}.csv")
"""

