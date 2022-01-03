from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Adapted from: https://stackoverflow.com/questions/36063533/clipping-a-voronoi-diagram-python/43023639#43023639
    
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp(axis=0).max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all([v >= 0 for v in vertices]):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

def get_cmap(n, name='Spectral'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def make_vor(dim, attribute_dict, position_dict, n_cell_types, results_dir, graph_id, node_id_list):

    c = get_cmap(n_cell_types)
    colors = [c(attribute_dict[i]) for i in node_id_list]

    ab = [dim, dim]
    plt.clf()
    fig, ax2 = plt.subplots(1,1,figsize=(8,8))

    '''# Plot the network 
    nx.draw(graph, pos = position_dict, ax=ax1, node_color = colors, node_size=17, with_labels = False)
    plt.axis('on')
    ax1.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    '''
    # Now plot the Voronoi Diagram 

    point_list = [position_dict[i] for i in node_id_list]
    vor = Voronoi(point_list)
    regions, vertices = voronoi_finite_polygons_2d(vor)

    for r in range(0, len(regions)):
        region = regions[r]
        polygon = vertices[region]
        ax2.fill(*zip(*polygon), alpha=0.8, facecolor=colors[r], edgecolor = 'k')
        
    # plt.plot(C[:,0], C[:,1], 'o', markersize=1, color = 'k')
    ax2.axis('equal')
    ax2.set(xlim=(0 - 0.01*ab[0], ab[0] + 0.01 * ab[0]), ylim=(0 - 0.01*ab[1], ab[1] + 0.01 * ab[1]))
    plt.savefig(results_dir + 'vor_' + str(graph_id) + '.pdf')
    plt.close()
    return 
