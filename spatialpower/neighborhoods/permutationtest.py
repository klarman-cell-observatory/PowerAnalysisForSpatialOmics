import numpy as np
import networkx as nx
from scipy import sparse

def create_subgraph(graph, subgraph_size, n=1):
    '''
    Returns the nodes in a breadth-first subgraph.

    Parameters
    ----------
        graph   :   NetworkX.Graph
            The graph object of the full tissue. 
        subgraph_size   :   int
            the size of the subgraph to generate (number of nodes)
        n   :   int
            the number of subgraphs to generate. Default = 1
    
    Returns
    -------
        subgraphs   :   Array-like
            list of n lists of nodes contained in subgraph.
    '''
    # subgraph_size = 500
    counter = 0
    subgraphs = []

    while counter < n:
        subgraph_nodes = []
        searched = []
        node_list = list(graph.nodes)

        start_node = np.random.randint(0, max(node_list))
        subgraph_nodes.append(start_node)

        for i in graph.neighbors(start_node):
            subgraph_nodes.append(i)

        searched.append(start_node)

        while len(subgraph_nodes) < subgraph_size:
            for node in subgraph_nodes:
                if node not in searched:
                    searched.append(node)  # Now we're searching that node.
                    for i in graph.neighbors(node):
                        if i not in subgraph_nodes:
                            subgraph_nodes.append(i)
                            if len(subgraph_nodes) >= subgraph_size:
                                break
                            else:
                                continue
                        else:
                            continue
                else:
                    continue

                if len(subgraph_nodes) >= subgraph_size:
                    break
                else:
                    continue

        subgraphs.append(subgraph_nodes)
        counter += 1

    return subgraphs

def shuffle_labels(ass_matrix, n_cell_types):
    '''
    Shuffles the cell assignment matrix.

    Parameters
    ----------
        ass_matrix: array_like
            The assignment matrix, B.

    Returns
    -------
        B: array_like
            An assignment matrix with shuffled assignmments.
    '''
    node_ids = [i for i in range(0, ass_matrix.shape[0])]  # This has to be for the original B.

    assignments = []

    for i in node_ids:
        counter = 0
        for j in range(0, n_cell_types):
            if ass_matrix[i, j] == 1:
                break
            else:
                counter += 1
        assignments.append(j)

    assignments = np.array(assignments)
    np.random.shuffle(assignments)  # in place shuffling

    shuffled_data = np.vstack((node_ids, assignments)).T
    B = np.zeros((shuffled_data.shape[0], n_cell_types))

    for i in range(0, shuffled_data.shape[0]):
        t = int(shuffled_data[i, 1])
        # print(i,t)
        B[i, t] = 1

    return B


def calculate_neighborhood_distribution(adj_matrix, ass_matrix):
    '''
    Calculates the probabilities of cell type adjacencies.

    Parameters
    ----------
        adj_matrix: array_like
            The N x N adjacency matrix for the graph.
        ass_matrix: array-like
            The assignment matrix, B.

    Returns
    -------
        H: array_like
            A K x K matrix where element (i, j) is the fraction
            of neighors of type i that are of type j.
    '''
    A = adj_matrix
    B = ass_matrix

    AB = np.matmul(A, B)  # Number of neighbors by type (cols), per node (row)

    edge_count_by_node = np.matmul(A, A.T)  # The diagonal of this matrix is the edge count for each node
    edge_count_by_node = np.diag(edge_count_by_node)

    aux1 = 1 / edge_count_by_node
    aux1[aux1 == np.inf] = 0
    aux1 = aux1 * np.identity(A.shape[0])

    cell_type_counts = np.matmul(B.T, B)  # The diagonal of this matrix is the number of cells of each type
    cell_type_counts = np.diag(cell_type_counts)

    aux2 = 1 / cell_type_counts
    aux2[aux2 == np.inf] = 0
    aux2 = aux2 * np.identity(B.shape[1])
    # diag(B.T*B)^-1

    H = np.matmul(np.matmul(aux2, B.T), np.matmul(aux1, AB))

    return H

def calculate_neighborhood_distribution_sparse(adj_matrix, ass_matrix):
    '''
    Calculates the probabilities of cell type adjacencies.

    Parameters
    ----------
        adj_matrix: sparse matrix
            The N x N adjacency matrix for the graph in scipy sparse format.
        ass_matrix: array-like
            The assignment matrix, B.

    Returns
    -------
        H: array_like
            A K x K matrix where element (i, j) is the fraction
            of neighors of type i that are of type j.
    '''
    A = adj_matrix
    B = ass_matrix

    edge_count_by_node = np.multiply(A, A.T)
    edge_count_by_node = edge_count_by_node.diagonal()
    aux1 = 1/edge_count_by_node
    aux1[aux1 == np.inf] = 0
    aux1 = sparse.identity(A.shape[0]).multiply(aux1)

    cell_type_counts = np.matmul(B.T, B) # The diagonal of this matrix is the number of cells of each type
    cell_type_counts = np.diag(cell_type_counts)

    aux2 = 1/cell_type_counts
    aux2[aux2 == np.inf] = 0
    aux2 = aux2 * np.identity(B.shape[1])

    AB = A*B

    aux3 = np.matmul(aux2, B.T)
    aux4 = aux1*AB

    H = np.matmul(aux3, aux4)
    
    return H

def calculate_enrichment_statistic(adj_matrix, ass_matrix, type_a, type_b):
    '''
    Calculates the probabilities of cell type adjacencies.

    Parameters
    ----------
        adj_matrix: sparse matrix
            The N x N adjacency matrix for the graph in scipy sparse format.
        ass_matrix: array-like
            The assignment matrix, B.
        type_a  : int
            index of first type in the interaction pair 
        type_b  : int
            index of first type in the interaction pair 

    Returns
    -------
        X: float
            the interaction enrichment statistic, centered at 0. 
    '''
    A = adj_matrix
    B = ass_matrix

    AB = A @ B
    C = B.T @ AB

    E = np.sum(np.sum(A, axis=1))/2 # Number of edges
    n = A.shape[0] #Number of cells
    i = type_a
    j = type_b

    f_a = np.sum(B, axis=0)[i]/n #in full graph
    f_b = np.sum(B, axis=0)[j]/n

    if i != j:
        N_ab = C[i,j]
    else:
        N_ab = C[i,j] / 2 # When i==j, you're on the diagonal and there's a double count.

    X_ab = N_ab/(2*f_a*f_b*E) - 1
    
    return X_ab

def perform_z_test(X1, X2):
    '''
    Tests the difference between the distributions of 
    enrichment statistics.

    Paramters
    ---------
        X1  :   np.array
            The array of x statistics from tissue 1
        X2  :   np.array
            The array of x statistics from tissue 2

    Returns
    -------
        z   : float
            The z statistics
        p   : float
            the p-value of the z statistic (1 sided)

    '''
    from scipy import stats 
    X1_bar = np.mean(X1)
    X2_bar = np.mean(X2)

    sigma_1 = np.std(X1)
    sigma_2 = np.std(X2)

    n_1 = len(X1)
    n_2 = len(X2)

    sem_1 = sigma_1/np.sqrt(n_1)
    sem_2 = sigma_2/np.sqrt(n_2)

    z = (np.abs(X1_bar-X2_bar)) / np.sqrt(sem_1 + sem_2)
    p = 1 - stats.norm.cdf(z)

    return z, p
    

def parse_subgraph(subgraph_nodes, graph, ass_matrix):
    """
    Convert list of nodes to subgraph induced by that list.

    Parameters
    ----------
        subgraph_nodes: array_like
            The list of nodes on which to induce the subgraph
        graph: nx.Graph
            The networkx graph of the full network.
        ass_matrix: array-like
            The assignment matrix, B, of the full network.

    Returns
    -------
        sg_adj: CSC Sparse Matrix
            The adjacency matrix of the subgraph.
        sg_ass: array_like
            The assignment matrix of the subgraph.
    """

    sg = graph.subgraph(subgraph_nodes)
    sg_adj = nx.to_scipy_sparse_matrix(sg, format='csc')  # New adjacency matrix.
    sg_ass = ass_matrix[list(sg.nodes)]

    return sg_adj, sg_ass


def permutation_test_trial(adj_matrix, ass_matrix, size, graph, n_cell_types):
    """
    Conducts a trial of the permutation test.
    Builds subgraph, shuffles network, then recalculates H. 
    
    Parameters
    ----------
        adj_matrix: array_like
            The N x N adjacency matrix for the graph. 
        ass_matrix: array-like
            The assignment matrix, B.
        size: int
            The size of the subgraph to calculate
        graph: nx.Graph
            The networkx graph of the full network.   
    
    Returns
    -------
        H: array_like
            A K x K matrix where element (i, j) is the fraction 
            of neighbors of type i that are of type j.
    """

    if size == adj_matrix.shape[0]:
        shuffled_graph = shuffle_labels(ass_matrix, n_cell_types)
        H = calculate_neighborhood_distribution_sparse(adj_matrix, shuffled_graph)

    else:
        subgraph_nodes = create_subgraph(graph, size, 1)[0]
        sg_adj, sg_ass = parse_subgraph(subgraph_nodes, graph, ass_matrix)
        shuffled_graph = shuffle_labels(sg_ass, n_cell_types)
        H = calculate_neighborhood_distribution_sparse(sg_adj, shuffled_graph)

    return H


def permutation_test_trial_wrapper(args):
    """
    Parallelization wrapper for the permutation test trial. 
    
    Parameters
    ----------
        args: tuple
            Format (adj_matrix, ass_matrix, size, graph, n_cell_types)
    
    Returns
    -------
        H: array_like
            the result of the permutation test trial
    """
    # print("starting " + str(mp.current_process()))
    adj_matrix = args[0]
    ass_matrix = args[1]
    size = args[2]
    graph = args[3]
    n_cell_types = args[4]
    # H = permumation_test_trial
    H = permutation_test_trial(adj_matrix, ass_matrix, size, graph, n_cell_types)
    # print("ending " + str(mp.current_process()))

    return H
