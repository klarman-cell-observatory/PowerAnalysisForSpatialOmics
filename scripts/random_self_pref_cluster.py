from glob import glob
import numpy as np
import scipy.sparse as sparse
import matplotlib.pyplot as plt
import networkx as nx
import operator
from spatialpower.tissue_generation import assign_labels
from spatialpower.tissue_generation import visualization

results_dir = './results/motif_detection/'
adj_mat_list = np.sort(glob(results_dir + 'blank_graph_network*.npy'))
pos_mat_list = np.sort(glob(results_dir + 'blank_graph_positions*.npy'))

dim = 300

##RANDOM##
cell_type_probabilities = np.ones(10) * 0.1
neighborhood_probabilities = np.ones((10,10)) * 0.1 
n_cell_types = len(cell_type_probabilities)

for ii in range(0, len(adj_mat_list)):
    A = np.load(adj_mat_list[ii])
    C = np.load(pos_mat_list[ii])
    
    j = adj_mat_list[ii].split('_')[-1].split('.')[0]
    # Blank assignment structure
    n_cell_types = len(cell_type_probabilities)
    position_dict = dict()
    for i in range(0, C.shape[0]):
        position_dict[i] = C[i, :]

    graph = nx.from_numpy_matrix(A)
    node_id_list = list(graph.nodes)
    attribute_dict = dict(zip(node_id_list, [-1 for i in graph.nodes]))
    
    attribute_dict = assign_labels.heuristic_assignment(graph, cell_type_probabilities, neighborhood_probabilities, 'region', dim, position_dict)
    observed_cell_type_dist, kl = assign_labels.check_cell_type_dist(n_cell_types, attribute_dict, cell_type_probabilities)
    observed_neighborhood_dist, kl_neighbor = assign_labels.check_neighborhood_dist(n_cell_types, attribute_dict, neighborhood_probabilities, graph, 1)
    B = assign_labels.build_assignment_matrix(attribute_dict, n_cell_types)
    np.save(results_dir + 'random_B_' + str(j) + '.npy', B)

    visualization.make_vor(dim, attribute_dict, position_dict, n_cell_types, results_dir, 'random_B_' + str(j), node_id_list)

## High Self Preference ##
'''cell_type_probabilities = [0.03, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.10, 0.11, 0.10]
neighborhood_probabilities = np.array([[0.50, 0.06, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05],
                                       [0.06, 0.11, 0.11, 0.11, 0.11, 0.10, 0.10, 0.10, 0.10, 0.10],
                                       [0.06, 0.11, 0.11, 0.11, 0.11, 0.10, 0.10, 0.10, 0.10, 0.10],
                                       [0.06, 0.11, 0.11, 0.11, 0.11, 0.10, 0.10, 0.10, 0.10, 0.10],
                                       [0.06, 0.11, 0.11, 0.11, 0.11, 0.10, 0.10, 0.10, 0.10, 0.10],
                                       [0.06, 0.10, 0.10, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11],
                                       [0.05, 0.10, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11],
                                       [0.05, 0.10, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11],
                                       [0.05, 0.10, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11],
                                       [0.05, 0.10, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11]])
n_cell_types = len(cell_type_probabilities)

for ii in range(0, len(adj_mat_list)):
    A = np.load(adj_mat_list[ii])
    C = np.load(pos_mat_list[ii])
    j = adj_mat_list[ii].split('_')[-1].split('.')[0]

    # Blank assignment structure
    n_cell_types = len(cell_type_probabilities)
    position_dict = dict()
    for i in range(0, C.shape[0]):
        position_dict[i] = C[i, :]

    graph = nx.from_numpy_matrix(A)
    node_id_list = list(graph.nodes)
    attribute_dict = dict(zip(node_id_list, [-1 for i in graph.nodes]))
    
    attribute_dict = assign_labels.heuristic_assignment(graph, cell_type_probabilities, neighborhood_probabilities, 'region', dim, position_dict)

    preferred_node_type = 0 
    for i in list(graph.nodes):
        if attribute_dict[i] == preferred_node_type:
            #print(i)
            graph_distance = 1
            neighborhood = nx.ego_graph(graph, i, radius = graph_distance)
            neighborhood_nodes = list(neighborhood.nodes)

            # Now set the remaining probabilities in the region. 

            for node in neighborhood_nodes:
                if node != i:
                    attribute_dict[node] = assign_labels.sample_cell_type(neighborhood_probabilities[preferred_node_type])
                else:
                    continue

    observed_cell_type_dist, kl = assign_labels.check_cell_type_dist(n_cell_types, attribute_dict, cell_type_probabilities)
    observed_neighborhood_dist, kl_neighbor = assign_labels.check_neighborhood_dist(n_cell_types, attribute_dict, neighborhood_probabilities, graph, 1)
    B = assign_labels.build_assignment_matrix(attribute_dict, n_cell_types)
    np.save(results_dir + 'selfpref_B_' + str(j) + '.npy', B)

    visualization.make_vor(dim, attribute_dict, position_dict, n_cell_types, results_dir, 'selfpref_B_' + str(j), node_id_list)'''

## 3 Cell Motif ##
cell_type_probabilities = [0.04, 0.04, 0.04, 0.13, 0.13, 0.13, 0.12, 0.12, 0.13, 0.12]
neighborhood_probabilities = np.array([[0.15, 0.40, 0.15, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04],
                                       [0.40, 0.06, 0.40, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02],
                                       [0.15, 0.40, 0.15, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04],
                                       [0.05, 0.02, 0.05, 0.13, 0.12, 0.13, 0.13, 0.13, 0.12, 0.12],
                                       [0.05, 0.02, 0.05, 0.12, 0.13, 0.13, 0.12, 0.12, 0.13, 0.13],
                                       [0.04, 0.02, 0.04, 0.13, 0.13, 0.13, 0.12, 0.13, 0.13, 0.13],
                                       [0.04, 0.02, 0.04, 0.13, 0.12, 0.12, 0.13, 0.13, 0.14, 0.13],
                                       [0.04, 0.02, 0.04, 0.13, 0.12, 0.13, 0.13, 0.12, 0.14, 0.13],
                                       [0.04, 0.02, 0.04, 0.12, 0.13, 0.13, 0.14, 0.14, 0.12, 0.12],
                                       [0.04, 0.02, 0.04, 0.12, 0.13, 0.13, 0.13, 0.13, 0.12, 0.14]])
n_cell_types = len(cell_type_probabilities)

for ii in range(0, len(adj_mat_list)):
    A = np.load(adj_mat_list[ii])
    C = np.load(pos_mat_list[ii])
    j = adj_mat_list[ii].split('_')[-1].split('.')[0]

    # Blank assignment structure
    n_cell_types = len(cell_type_probabilities)
    position_dict = dict()
    for i in range(0, C.shape[0]):
        position_dict[i] = C[i, :]

    graph = nx.from_numpy_matrix(A)
    node_id_list = list(graph.nodes)
    attribute_dict = dict(zip(node_id_list, [-1 for i in graph.nodes]))
    
    attribute_dict = assign_labels.heuristic_assignment(graph, cell_type_probabilities, neighborhood_probabilities, 'region', dim, position_dict)

    #preferred_node_type = 0 
    for i in list(graph.nodes):
        if ((attribute_dict[i] == 0) or (attribute_dict[i] == 1) or (attribute_dict[i] == 2)):
            #print(i)
            graph_distance = 1
            neighborhood = nx.ego_graph(graph, i, radius = graph_distance)
            neighborhood_nodes = list(neighborhood.nodes)

            # Now set the remaining probabilities in the region. 

            for node in neighborhood_nodes:
                if node != i:
                    attribute_dict[node] = assign_labels.sample_cell_type(neighborhood_probabilities[attribute_dict[i]])
                else:
                    continue

    observed_cell_type_dist, kl = assign_labels.check_cell_type_dist(n_cell_types, attribute_dict, cell_type_probabilities)
    observed_neighborhood_dist, kl_neighbor = assign_labels.check_neighborhood_dist(n_cell_types, attribute_dict, neighborhood_probabilities, graph, 1)
    B = assign_labels.build_assignment_matrix(attribute_dict, n_cell_types)
    np.save(results_dir + '3cellmotif_B_' + str(j) + '.npy', B)

    visualization.make_vor(dim, attribute_dict, position_dict, n_cell_types, results_dir, '3cellmotif_B_' + str(j), node_id_list)