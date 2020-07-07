import numpy as np
import networkx as nx
import operator

def get_unassigned_nodes():
    unassigned = []
    for n in node_id_list:
        if attribute_dict[n] == -1:
            unassigned.append(n)
    return unassigned

def sample_cell_type(dist):
    x_i = np.random.multinomial(1, dist)
    for i in range(0, len(x_i)):
        if x_i[i] == 1:
            return i

def kl_divergence(p, q):
    '''
    Calculate the KL Divergence between two distrbitions

    Parameters
    ----------
    p : Array-like, distribution 1
    q : Array-like, distribution 2

    Returns
    -------
    div : float, KL divergence
    '''
    p = np.array(p)
    q = np.array(q)
    
    e = 0.00000001 # this is just to avoid div 0
    
    p = p + e
    q = q + e
    
    div = np.sum(p*np.log(p/q))
    return div

def get_type_proportions(n_types, observed):
    try:
        unique, counts = np.unique(observed, return_counts = True)
    except:
        unique, counts = onp.unique(observed, return_counts = True)
        
    if len(unique) == n_types:
            return counts/len(observed)
    
    else:
        proportions = []
        for i in range(0, n_types):
            counter = 0
            for j in observed:
                if j == i:
                    counter += 1
            proportions.append(counter/len(observed))
        return proportions

def check_cell_type_dist(n_cell_types, attribute_dict, cell_type_probabilities, silent=False):
    n_cells = len(attribute_dict)
    values = list(attribute_dict.values())
    
    type_proportions = get_type_proportions(n_cell_types, values)
    
    try:
        type_proportions = np.round(get_type_proportions(n_cell_types, list(attribute_dict.values())), decimals=3)
    except:
        type_proportions = onp.round(get_type_proportions(n_cell_types, list(attribute_dict.values())), decimals=3)
    
    
    kl = kl_divergence(type_proportions, cell_type_probabilities)
    
    if silent == False:
        print("Expected Cell Type Dist: " + str(cell_type_probabilities))
        print("Observed Cell Type Dist: " + str(type_proportions))
        print("K-L Divergence: " + str(kl))
    
    return type_proportions, kl

def check_neighborhood_dist(n_cell_types, attribute_dict, neighborhood_probabilities, graph, d, silent=False):
    values = np.array(list(attribute_dict.values()))
    nodes_of_type = []
    results = [[] for i in range(0, n_cell_types)]

    for i in range(0, n_cell_types): 
        nodes_of_type.append(np.where(values == i))

    for i in range(0, n_cell_types):
        node_arr = nodes_of_type[i][0]

        for node in node_arr:
            distance = d

            neighborhood = nx.ego_graph(graph, node, radius = distance)
            neighborhood_nodes = list(neighborhood.nodes)
            if len(neighborhood_nodes) > 1:
                neighborhood_nodes = [n for n in neighborhood_nodes if n != node]
                neighborhood_types = [attribute_dict[n] for n in neighborhood_nodes]
                neighborhood_proportions = get_type_proportions(n_cell_types, neighborhood_types)
                results[i].append(neighborhood_proportions)

    mean_distributions = [] 

    for i in range(0, n_cell_types):
        trials = results[i]
        trials = np.array(trials)
        sums = np.sum(trials, axis = 0)
        mean = np.divide(sums, len(nodes_of_type[i][0]))
        mean = np.round(mean, decimals = 3)
        mean_distributions.append(mean)

    mean_distributions = np.array(mean_distributions)
    kl = kl_divergence(mean_distributions, neighborhood_probabilities)
    
    if silent == False:
        print("Observed Neighborhood Distribution:")
        print(str(mean_distributions))
        print("Expected Neighborhood Distribution:")
        print(str(neighborhood_probabilities))
        print("K-L Divergence: " + str(kl))
    
    return mean_distributions, kl

def swap_types(observed_neighborhood, expected_neighborhood, attribute_dict, region_nodes):
    dif = observed_neighborhood - expected_neighborhood
    max_type = np.argmax(dif) #max difference
    min_type = np.argmin(dif)
    n_in_region = len(region_nodes)

    n_max_type = n_in_region * np.abs(observed_neighborhood[max_type])
    n_min_type = n_in_region * np.abs(observed_neighborhood[min_type])

    expected_max_type = np.rint(expected_neighborhood[max_type] * n_in_region)
    expected_min_type = np.rint(expected_neighborhood[min_type] * n_in_region)
   
    if n_min_type <= n_max_type:
        expected_to_swap = n_max_type - expected_max_type

        if expected_to_swap > expected_min_type:
            n_to_swap = expected_min_type
        else:
            n_to_swap = expected_to_swap

        # select nodes to swap
        max_type_nodes = []

        for node in region_nodes:
            if attribute_dict[node] == max_type:
                max_type_nodes.append(node)

        swap_count = 0
        if len(max_type_nodes) == 0:
            pass
        else:
            #print('here')
            while swap_count < n_to_swap:
                node_to_swap = np.random.choice(max_type_nodes)
                attribute_dict[node_to_swap] = min_type
                swap_count += 1

    return attribute_dict

def heuristic_assignment(graph, cell_type_probabilities, neighborhood_probabilities, mode, dim, position_dict):
    n_cell_types = neighborhood_probabilities.shape[0]
    node_id_list = list(graph.nodes)
    attribute_dict = dict(zip(node_id_list, [-1 for i in graph.nodes]))
    
    x_lims = [z for z in range(0, dim+1, 50)]
    y_lims = [z for z in range(0, dim+1, 50)]

    j = 0
    while j < len(x_lims) - 1:
        x_min = x_lims[j]
        x_max = x_lims[j + 1]
        k = 0
        while k < len(y_lims) - 1:
            y_min = y_lims[k]
            y_max = y_lims[k + 1]

            region_nodes = []
            for i in node_id_list:
                x = position_dict[i][0]
                if x_min <= x < x_max or (x_max == dim and x_min <= x <= x_max):
                    y = position_dict[i][1]
                    if y_min <= y < y_max or (y_max == dim and y_min <= y <= y_max):
                        region_nodes.append(i)
                else:
                    continue
            if mode == 'region':
                # Pick a starting point by eigenvector centrality and sample its type
                subgraph = graph.subgraph(region_nodes)
                subgraph_centrality = nx.eigenvector_centrality(subgraph, max_iter=1000)
                start_node_id = max(subgraph_centrality.items(), key = operator.itemgetter(1))[0]

                start_node_type = sample_cell_type(cell_type_probabilities)
                attribute_dict[start_node_id] = start_node_type # set it in the attributes dict

                # Now set the remaining probabilities in the region. 

                for node in region_nodes:
                    if node != start_node_id:
                        attribute_dict[node] = sample_cell_type(neighborhood_probabilities[start_node_type])
                    else:
                        continue
                        
            elif mode == 'graph':
                 # Pick a starting point by at random.
                start_node_id = np.random.choice(region_nodes)
                start_node_type = sample_cell_type(cell_type_probabilities)
                attribute_dict[start_node_id] = start_node_type # set it in the attributes dict

                # Calculate the neighborhood of the start node. 
                graph_distance = 1
                neighborhood = nx.ego_graph(graph, start_node_id, radius = graph_distance)
                neighborhood_nodes = list(neighborhood.nodes)

                # Now set the remaining probabilities in the region. 

                for node in neighborhood_nodes:
                    if node != start_node_id:
                        attribute_dict[node] = sample_cell_type(neighborhood_probabilities[start_node_type])
                    else:
                        continue
            k += 1
        j += 1

    # Shift and Recalculate
    x_lims = [z for z in range(25, dim + 1, 50)]
    y_lims = [z for z in range(25, dim + 1, 50)]

    j = 0
    while j < len(x_lims) - 1:
        x_min = x_lims[j]
        x_max = x_lims[j + 1]
        k = 0
        while k < len(y_lims) - 1:
            y_min = y_lims[k]
            y_max = y_lims[k + 1]

            region_nodes = []
            for i in node_id_list:
                x = position_dict[i][0]
                if x_min <= x < x_max or (x_max == dim and x_min <= x <= x_max):
                    y = position_dict[i][1]
                    if y_min <= y < y_max or (y_max == dim and y_min <= y <= y_max):
                        region_nodes.append(i)
                else:
                    continue

            if mode == 'region':
                # Pick a starting point by eigenvector centrality and sample its type
                subgraph = graph.subgraph(region_nodes)
                subgraph_centrality = nx.eigenvector_centrality(subgraph, max_iter = 1000)
                start_node_id = max(subgraph_centrality.items(), key = operator.itemgetter(1))[0]
                start_node_type = attribute_dict[start_node_id]

                # Calculate the observed and expected distributions of cell type in the neighborhood
                expected_neighborhood = neighborhood_probabilities[start_node_type]
                observed_neighborhood = [attribute_dict[x] for x in region_nodes]
                observed_neighborhood = get_type_proportions(n_cell_types, observed_neighborhood)

                divergence = kl_divergence(observed_neighborhood, expected_neighborhood)
                # print("First divergence: ", divergence)
                if divergence > 0.25:
                    #print(divergence)
                    div = 100
                    attempts = 0
                    while div > 0.25 and attempts < 50:
                        attribute_dict = swap_types(observed_neighborhood, expected_neighborhood, 
                                                attribute_dict, region_nodes)
                        expected_neighborhood = neighborhood_probabilities[start_node_type]
                        observed_neighborhood = [attribute_dict[x] for x in region_nodes]
                        observed_neighborhood = get_type_proportions(n_cell_types, observed_neighborhood)
                        div = kl_divergence(observed_neighborhood, expected_neighborhood)
                        attempts += 1
                   # print(div)
            elif mode == 'graph':
                # Pick a starting point by at random.
                start_node_id = np.random.choice(region_nodes)
                start_node_type = sample_cell_type(cell_type_probabilities)
                attribute_dict[start_node_id] = start_node_type # set it in the attributes dict

                # Calculate the neighborhood of the start node. 
                graph_distance = 1
                neighborhood = nx.ego_graph(graph, start_node_id, radius = graph_distance)
                neighborhood_nodes = list(neighborhood.nodes)

                # Now set the remaining probabilities in the region. 

                for node in neighborhood_nodes:
                    if node != start_node_id:
                        attribute_dict[node] = sample_cell_type(neighborhood_probabilities[start_node_type])
                    else:
                        continue
            k += 1
        j += 1
    
    if mode == 'graph':
        unassigned = get_unassigned_nodes(node_id_list, attribute_dict)

        while len(unassigned) > 0:
            start_node_id = np.random.choice(unassigned)
            start_node_type = sample_cell_type(cell_type_probabilities)
            attribute_dict[start_node_id] = start_node_type # set it in the attributes dict

             # Calculate the neighborhood of the start node. 
            graph_distance = 1
            neighborhood = nx.ego_graph(graph, start_node_id, radius = graph_distance)
            neighborhood_nodes = list(neighborhood.nodes)

            # Now set the remaining probabilities in the region. 

            for node in neighborhood_nodes:
                if node != start_node_id:
                    attribute_dict[node] = sample_cell_type(neighborhood_probabilities[start_node_type])
                else:
                    continue

            unassigned = get_unassigned_nodes(node_id_list, attribute_dict)

    return attribute_dict

def build_assignment_matrix(attribute_dict, n_cell_types):
    data = list(attribute_dict.items())
    data = np.array(data) # Assignment matrix
    
    B = np.zeros((data.shape[0],n_cell_types)) # Empty matrix
    
    for i in range(0, data.shape[0]):
        t = data[i,1]
        B[i,t] = 1
    
    return B 