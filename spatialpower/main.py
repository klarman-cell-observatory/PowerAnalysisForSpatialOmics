import numpy as np
import networkx as nx
import neighborhoods.permutationtest as perm_test 
import neighborhoods.neighborhoods as nbr 
from scipy import sparse
from datetime import datetime
#import yappi

if __name__ == '__main__':

    # Data Loading

    A = np.load('../results/sample_adjmat_20200601.npy')
    B = np.load('../results/sample_ass_matrix_4types_20200601.npy')

    # Build required structures
    graph = nx.from_numpy_matrix(A)
    node_id_list = list(graph.nodes)
    attribute_dict = dict(zip(node_id_list, [-1 for i in graph.nodes]))
    n_cell_types = B.shape[1]
    S = sparse.coo_matrix(A)
    S_csc = S.tocsc() # Sparse Adj Mat

    # Configure trial parameters 
    trials = 1000
    max_size = S_csc.shape[0]
    size_list = [i for i in range(0, max_size, 50)]

    if max(size_list) != max_size:
        size_list.append(max_size)

    # Configure results
    out_dir = '../results'
    now = datetime.now()
    datestamp = date_time = now.strftime("%m%d%Y")
    results_path = str(out_dir) + '/neighborhood_perm_test_' + str(datestamp) + '/'

    size = max_size
    H_gt = perm_test.calculate_neighborhood_distribution_sparse(S_csc, B)

    nbr.run_test(results_path, S_csc, B, H_gt, size, n_jobs=-1, trials=750, plot=False, graph=graph, graph_id=None, threshold=0.1)


    frac_list = np.linspace(0, 1, 11)  # Measure every 10%
    size_list = np.round(frac_list * max_size).astype(int)

    same_size_trials = 1
    size = size_list[1]

    '''for size in size_list[1:-1]:
        print(size)
        for j in range(0, same_size_trials):
            # Generate the subgraph
            sg = perm_test.create_subgraph(graph, size, 1)[0]
            sg_A, sg_B = perm_test.parse_subgraph(sg, graph, B)

            # Get the ground truth distributions for this subgraph.
            sg_gt = perm_test.calculate_neighborhood_distribution_sparse(sg_A, sg_B)

            # Now permute this graph
            nbr.run_test(results_path, sg_A, sg_B, sg_gt, size, n_jobs=-1, trials=750, plot=False, graph = graph, graph_id=j,
                     threshold=0.1)'''