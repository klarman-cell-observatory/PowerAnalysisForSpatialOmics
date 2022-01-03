import neighborhoods.permutationtest as perm_test
import numpy as np
import networkx as nx
import multiprocessing as mp
from datetime import datetime
import errno
from joblib import Parallel, delayed
import os
from glob import glob


def build_assignment_matrix(attribute_dict, n_cell_types):
    data = list(attribute_dict.items())
    data = np.array(data)  # Assignment matrix

    B = np.zeros((data.shape[0], n_cell_types))  # Empty matrix

    for i in range(0, data.shape[0]):
        t = data[i, 1]
        B[i, t] = 1

    return B


def parse_results(results, size, out_dir):
    print("Writing results...")
    print(str(out_dir))
    print(len(results))
    print(results[0])
    for i in range(0, len(results)):
        arr = results[i]
        np.save(str(out_dir) + str(size) + 'cells_shuffle' + str(i), arr)
    return


def run_test(results_path, A, B, H_gt, size, n_jobs, trials, plot, graph, graph_id, threshold):
    '''
    Runs the permutation test, and calculates signficant interaction pairs.

    Parameters
    ----------
        results_path: str, the root results dir
        size : int, size of graph to calculate.
        n_jobs: int, number of parallel jobs to spawn
        trials: int, number of shuffles in empirical distribution
        plot : bool, generate histogram of each pairwise relation if True.

    Returns
    -------
        None
    '''
    # Make results dir
    try:
        os.mkdir(results_path)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    # Perform calculations.
    results = []
    if graph_id == None:
        out_dir = results_path + str(size) + '_cells/'
    else:
        out_dir = results_path + str(size) + '_cells_' + str(graph_id) + '/'

    try:
        os.mkdir(out_dir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    
    n_cell_types = B.shape[1]
    args = (A, B, size, graph, n_cell_types)
    arg_list = [args for i in range(0, trials)]
    results = Parallel(n_jobs=n_jobs, verbose=1, backend="sequential")(
        delayed(perm_test.permutation_test_trial_wrapper)(args) for args in arg_list)
    #parse_results(results, size, out_dir)

    # Process results

    '''# size_list = []
    result_list = []

    file_list = glob(out_dir + '*.npy')
    for f in file_list:
        arr = np.load(f)
        # size_list.append(size)
        result_list.append(arr)'''

    arr = np.dstack(results)  # stack into a 3-D array
    n_types = arr.shape[0]

    enriched_pairs = []
    depleted_pairs = []

    for i in range(0, n_types):
        for j in range(0, n_types):
            ground_truth_score = H_gt[i, j]
            emp_dist = arr[i, j, :]
            indices, = np.where(emp_dist < ground_truth_score)
            p = (len(emp_dist) - len(indices) + 1) / (len(emp_dist) + 1)
            if p <= threshold:
                enriched_pairs.append([i, j, p])
            elif p >= 1 - threshold:
                depleted_pairs.append([i, j, p])

            # Visualize empirical distribution
            if plot == True:
                plt.clf()
                # sns.set(style = 'white')
                plt.hist(arr[2, 2, :], color='k')
                plt.xlim(0, 1)
                plt.xlabel("Probability of Interaction between " + str(i) + " and " + str(j))
                plt.ylabel("Count")
                plt.savefig(out_dir + "distplot_" + str(i) + "_" + str(j) + ".pdf")

    # Write results matrix.
    np.save(out_dir + "enriched_pairs.npy", np.array(enriched_pairs))
    np.save(out_dir + "depleted_pairs.npy", np.array(depleted_pairs))

    return

def run_test_nosave(A, B, H_gt, size, n_jobs, trials, graph, threshold):
    '''
    Runs the permutation test, and calculates signficant interaction pairs.

    Parameters
    ----------
        size : int, size of graph to calculate.
        n_jobs: int, number of parallel jobs to spawn
        trials: int, number of shuffles in empirical distribution
        plot : bool, generate histogram of each pairwise relation if True.

    Returns
    -------
        enriched_pairs  :   array-like
        depleted_pairs  :   array-like
    '''
   
    n_cell_types = B.shape[1]
    args = (A, B, size, graph, n_cell_types)
    arg_list = [args for i in range(0, trials)]
    results = Parallel(n_jobs=n_jobs, verbose=1, backend="sequential")(
        delayed(perm_test.permutation_test_trial_wrapper)(args) for args in arg_list)
    #parse_results(results, size, out_dir)

    arr = np.dstack(results)  # stack into a 3-D array
    n_types = arr.shape[0]

    enriched_pairs = []
    depleted_pairs = []

    for i in range(0, n_types):
        for j in range(0, n_types):
            ground_truth_score = H_gt[i, j]
            emp_dist = arr[i, j, :]
            indices, = np.where(emp_dist < ground_truth_score)
            p = (len(emp_dist) - len(indices) + 1) / (len(emp_dist) + 1)
            if p <= threshold:
                enriched_pairs.append([i, j, p])
            elif p >= 1 - threshold:
                depleted_pairs.append([i, j, p])

    # Write results matrix.
    np.save(out_dir + "enriched_pairs.npy", np.array(enriched_pairs))
    np.save(out_dir + "depleted_pairs.npy", np.array(depleted_pairs))

    return enriched_pairs, depleted_pairs