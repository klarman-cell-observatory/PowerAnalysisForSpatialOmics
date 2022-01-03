##Used to generate all of the tiles for the tile shuffling experiement.
import numpy as np
import networkx as nx
from spatialpower.tissue_generation import assign_labels


def build_assignment_matrix(attribute_dict, n_cell_types):
    data = list(attribute_dict.items())
    data = np.array(data)  # Assignment matrix

    B = np.zeros((data.shape[0], n_cell_types))  # Empty matrix

    for i in range(0, data.shape[0]):
        t = data[i, 1]
        B[i, t] = 1

    return B

#Set parameters
n_images = 20
b_follicle_count = 51 * n_images
pals_count = 94 * n_images 
red_pulp_count = 174 * n_images
marginal_zone_count = 69 * n_images

#Load relevant data
adjacency_matrix = np.load('./adj_mat_239_cell_tissue_scaffold_tile.npy')
C = np.load('./C_239_cell_tissue_scaffold_tile.npy')
R = np.load('./R_239_cell_tissue_scaffold_tile.npy')
graph = nx.from_numpy_matrix(adjacency_matrix)



#B Follicle
p = np.load('./spleen_data/for_paper/p_bfollicle_balbc1.npy')
H = np.load('./spleen_data/for_paper/H_bfollicle_balbc1.npy')

position_dict = dict()
for i in range(0, C.shape[0]):
    position_dict[i] = C[i, :]

i=0
while i <= b_follicle_count:
    attribute_dict = assign_labels.heuristic_assignment(graph, p, H, 'region', 350, position_dict, 100)
    heuristic_B = build_assignment_matrix(attribute_dict, 27)
    np.save('./spleen_data/for_paper/tiles/239_cell_tiles/shuffling_experiment/B_hueristic_bfollicle_' + str(i) + '.npy', heuristic_B)
    if i % 10 == 0:
        print(i)
    i += 1

#PALS Zone
p = np.load('./spleen_data/for_paper/p_pals_balbc1.npy')
H = np.load('./spleen_data/for_paper/H_pals_balbc1.npy')

n_cell_types = 27
position_dict = dict()
for i in range(0, C.shape[0]):
    position_dict[i] = C[i, :]

i=0
while i <= pals_count:
    attribute_dict = assign_labels.heuristic_assignment(graph, p, H, 'region', 350, position_dict, 50)
    heuristic_B = build_assignment_matrix(attribute_dict, 27)
    np.save('./spleen_data/for_paper/tiles/239_cell_tiles/shuffling_experiment/B_hueristic_pals_' + str(i) + '.npy', heuristic_B)
    if i % 10 == 0:
        print(i)
    i += 1

#Red Pulp

p = np.load('./spleen_data/for_paper/p_redpulp_balbc1.npy')
H = np.load('./spleen_data/for_paper/H_redpulp_balbc1.npy')

n_cell_types = 27
position_dict = dict()
for i in range(0, C.shape[0]):
    position_dict[i] = C[i, :]

i=0
while i <= red_pulp_count:
    attribute_dict = assign_labels.heuristic_assignment(graph, p, H, 'region', 350, position_dict, 50)
    heuristic_B = build_assignment_matrix(attribute_dict, 27)
    np.save('./spleen_data/for_paper/tiles/239_cell_tiles/shuffling_experiment/B_hueristic_redpulp_' + str(i) + '.npy', heuristic_B)
    if i % 10 == 0:
        print(i)
    i += 1

#Marginal Zone
p = np.load('./spleen_data/for_paper/p_marginalzone_balbc1.npy')
H = np.load('./spleen_data/for_paper/H_marginalzone_balbc1.npy')


n_cell_types =27
position_dict = dict()
for i in range(0, C.shape[0]):
    position_dict[i] = C[i, :]

i=0
while i <= marginal_zone_count:
    attribute_dict = assign_labels.heuristic_assignment(graph, p, H, 'region', 350, position_dict, 50)
    heuristic_B = build_assignment_matrix(attribute_dict, 27)
    np.save('./spleen_data/for_paper/tiles/239_cell_tiles/shuffling_experiment/B_hueristic_marginalzone_' + str(i) + '.npy', heuristic_B)
    if i % 10 == 0:
        print(i)
    i += 1

