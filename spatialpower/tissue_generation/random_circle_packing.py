import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import random
import argparse
import time
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

def sample_line_segment(x_a, x_b, rc_a, rc_b, r_min, r_max):
    r = x_b - x_a
    l = np.linalg.norm(r, ord=2)  # L2 norm
    r = r / np.linalg.norm(l)
    dr = r_max - r_min
    center = x_a
    radius = rc_a
    while True:
        R_new = dr * np.random.rand() + r_min
        C_new = center[-1] + r * (radius[-1] + R_new + r_max * np.random.rand())
        d = l - np.linalg.norm(C_new + r * (R_new + rc_b) - x_a, ord=2)
        if d < 2 * (r_min + 1e-12):
            center = np.vstack((center, x_b))
            radius = np.vstack((radius, rc_b))
            break
        else:
            center = np.vstack((center, C_new))
            radius = np.vstack((radius, R_new))

    return center, radius


def update_grid(X_new, R_new, G, r_min):
    _D = (G - X_new)
    _D = np.power(_D, 2)
    D = np.sum(_D, axis=1)
    thresh = np.power((R_new + r_min + 1e-12), 2)
    mask = D < thresh
    mask = np.invert(mask)
    G = G[mask]
    return G, mask


def set_color(r_list, max_r, min_r, n_types):
    #color_list = ['#24557D', '#156944', '#5EC358', '#EB9420', '#E54116', '#F4BF1F', '#BC1826', '#CCB98E',
                #  '#90B3B7', '#C4B6C1', '#E6ED40', '#5ECCB7', '#BB4E64', '#8C4FC9', '#E05B7C', '#26A8BD']
    color_list = [(.14,.33,.49), (.08,.41,.27), (.37,.76,.35), (.92,.58,.13), (.90,.25,.09), (.96,.75,.12), (.74,.09,.15),
                  (0.8,0.76,0.56), (0.56,0.7,0.72), (0.77,0.71,.76)]
    random.shuffle(color_list)
    assert n_types <= len(color_list)
    boundary_list = []
    increment = (max_r - min_r) / n_types
    assigned_colors = []

    for j in range(0, n_types):
        b = min_r + j * increment
        boundary_list.append(b)

    for r in r_list:

        upper_bound_idx = -1  # void value
        for i in range(0, len(boundary_list)):
            if r < boundary_list[i]:
                upper_bound_idx = i
                break
        assigned_colors.append(color_list[upper_bound_idx])

    return assigned_colors


def circle_intersect_check(r_0, r_1, x_0, x_1, y_0, y_1):
    """
    Check if two circles intersect. 

    Parameters
    ----------
        r_0 :   float
            radius of search circle
        r_1 :   float
            radius of check circle
        x_0 :   float
            x coord of search circle center
        x_1 :   float
            x coord of check circle center
        y_0 :   
            y coord of search circle center
        y_1 :   
            x coord of check circle center
    Returns
    --------
        indices :   np.Array
            Numpy array of the circle indicies that intersect.
    """

    distance_squared = np.power(x_0 - x_1, 2) + np.power(y_0 - y_1, 2)
    dr_squared = np.power(r_0 - r_1, 2)
    r_sum_squared = np.power(r_0 + r_1, 2)

    distance_squared = distance_squared.reshape(distance_squared.size, 1)
    dr_squared = dr_squared.reshape(dr_squared.size, 1)
    r_sum_squared = r_sum_squared.reshape(r_sum_squared.size, 1)

    # print(distance_squared.shape, dr_squared.shape, r_sum_squared.shape)

    mask_a = dr_squared <= distance_squared
    mask_b = distance_squared <= r_sum_squared
    mask_a = mask_a.astype(int)  # 1 if true
    mask_b = mask_b.astype(int)  # 1 if True
    mask = mask_a + mask_b
    mask = mask == 2  # Meets both conditions
    mask = mask.reshape(mask.size)
    mask = mask.astype(int)
    indices = mask.nonzero()

    return indices


def voronoi_finite_polygons_2d(vor, radius=None):
    """
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
        radius = vor.points.ptp(axis=0).max() * 2

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

            t = vor.points[p2] - vor.points[p1]  # tangent
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
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Configuration parameters for circle packing.")
    parser.add_argument('-x', '--width', type=float, default=500, help="Width of constraining rectangle")
    parser.add_argument('-y', '--height', type=float, default=500, help="Height of the contstraining rectangle")
    parser.add_argument('--rmax', type=float, default=25, help="Maximum circle radius.")
    parser.add_argument('--rmin', type=float, default=4, help="Minimum circle radius.")
    parser.add_argument('--visualization', action="store_true",
                        help="If present, generate and save a visualization of the circle packing.")
    # parser.add_argument('--constrained', action="store_true",
    #  help="If present, do not allow circle edges past boundary.")
    parser.add_argument('-o', '--outdir', type=str, default='./', help="The folder in which to store results.")
    parser.add_argument('-n', '--n_colors', type=int, default=10, help='The number of bins to color cell size.')
    parser.add_argument('-g', '--graph', default=True,
                        help='If present, draw the graph representation of the circle packing.', action="store_true")
    parser.add_argument('-v', '--voronoi', default=True, action="store_true", help='If present, draw the voronoi.')
    args = parser.parse_args()

    ab = np.array([args.width, args.height])

    r_min = args.rmin
    r_max = args.rmax
    visualization = args.visualization
    constrained = False

    if constrained:
        pass
    else:
        constrained = False

    outdir = args.outdir

    dx = max(min(ab) / 2e3, r_min / 50)  # Sets up the increment.
    x = np.arange(0, ab[0] + dx, dx)
    y = np.arange(0, ab[1] + dx, dx)

    x, y = np.meshgrid(x, y)

    x_col = np.reshape(x, (np.size(x), 1))
    y_col = np.reshape(y, (np.size(y), 1))
    G = np.concatenate((x_col, y_col), axis=1)

    # Start by placing circles around the edges if constrained is false.

    dr = r_max - r_min
    corners = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    corner_vertices = np.multiply(ab, corners)

    if not constrained:
        x_a = corner_vertices
        x_b = np.roll(corner_vertices, [-1, -1])  # This command empirically seems equivalent to circshift(X, [-1 0])

        rc = np.multiply(dr, np.random.rand(4, 1)) + r_min
        rc_a = rc
        rc_b = np.roll(rc, [-1])

        C = []
        R = []

        for i in range(0, 4):
            Ci, Ri = sample_line_segment(x_a[i, :], x_b[i, :], rc_a[i], rc_b[i], r_min, r_max)
            Ci = Ci[:-1, :]
            # C[i] = Ci
            C.append(Ci)
            Ri = Ri[:-1, :]
            # R[i] = Ri
            R.append(Ri)

        R = np.vstack((R[0], R[1], R[2], R[3]))
        C = np.vstack((C[0], C[1], C[2], C[3]))

        for i in range(0, C.shape[0]):
            G, _ = update_grid(C[i, :], R[i], G, r_min)

    else:
        G_max = G + r_min + 1e-12
        G_min = G - r_min - 1e-12
        chk_in = (G_max <= ab) & (G_min >= [0, 0])
        chk_in = chk_in.astype(int)
        chk_in = np.sum(chk_in, axis=1) == 2
        # keep [T, T]
        G = G[chk_in]
        C = []
        R = []

    circle_list = []

    t = np.linspace(0, 2 * np.pi, 100)[np.newaxis]
    t = t.T
    sin = np.sin(t)
    cos = np.cos(t)
    P = np.hstack((sin, cos))
    if np.count_nonzero(C) != 0:  # Not empty
        for i in range(0, C.shape[0]):
            Pm = R[i] * P + C[i, :]
            circle_list.append(Pm)

    # Rejection sampling to populate interior of rectangle
    flag = True
    n = 0
    cnt = 0
    m = 0
    Ng = G.shape[0]

    while (G.shape[0] != 0) and cnt < 3e5:
        n += 1

        # New circle
        if flag and (cnt > 500 or (G.shape[0] < 0.95 * Ng)):
            # print("CONDITION 1")
            flag = False
            Rg = r_max * np.ones((G.shape[0], 1))

        i = []
        if cnt <= 500 and flag:
            # print("CONDITION 2")
            X_new = np.multiply(ab, np.random.rand(1, 2))
        else:
            # print("CONDITION 3")
            i = np.random.randint(0, G.shape[0])
            X_new = G[i, :] + (dx / 2) * (2 * np.random.rand(1, 2) - 1)
            X_new = np.minimum(np.maximum(X_new, np.array([0, 0])), ab)
            if cnt > 1e3:
                Rg[:] = np.maximum(0.95 * Rg, r_min)

        if np.count_nonzero(i) == 0:
            R_new = dr * np.random.rand() + r_min  # radius
        else:
            R_new = (Rg[i] - r_min) * np.random.rand() + r_min

        if constrained:
            # print("CONSTRAINED CONDITION")
            X_new_max = X_new + R_new + 1e-12
            X_new_min = X_new - R_new - 1e-12
            mask_max = X_new_max <= ab
            mask_min = X_new_min >= np.array([0, 0])
            mask = mask_max & mask_min
            mask = mask.astype(int)
            if np.sum(mask, axis=1) < 2:
                cnt += 1
                continue
        if np.count_nonzero(C) != 0:  # not empty
            # print("NON EMPTY C")
            _d_in = C - X_new
            _d_in = np.power(_d_in, 2)
            d_in = np.sqrt(np.sum(_d_in, axis=1))
            v = d_in
            v = v.reshape(v.size, 1)
            mask = v < (R + R_new)
            mask = mask.astype(int)
            if np.sum(mask) > 0:
                cnt += 1
                if np.count_nonzero(i) != 0:  # Not empty
                    Rg[i] = min([0.95 * Rg[i], min([0.99 * (R_new + dx / 2), r_max])])
                    Rg[i] = max([Rg[i], r_min])
                continue

        # Accept new circle
        cnt = 0
        m += 1
        C = np.vstack((C, X_new))
        R = np.vstack((R, R_new))
        G, mask = update_grid(X_new, R_new, G, r_min)

        if not flag:
            Rg = Rg[mask]  # make sure this doesn't need to be flipped.

        if visualization:
            Pm = R_new * P + X_new
            circle_list.append(Pm)

    time_stamp = str(time.strftime("%Y-%m-%d-%H%M%S"))
    assigned_colors = set_color(R, r_max, r_min, args.n_colors)

    if visualization:
        plt.clf()
        #fig = Figure()
        #canvas = FigureCanvas(fig)

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set(xlim=(0 - 0.05 * ab[0], ab[0] + 0.05 * ab[0]), ylim=(0 - 0.05 * ab[1], ab[1] + 0.05 * ab[1]))
        ax.axis('equal')
        ax.axis('off')

        for i in range(0, R.shape[0]):
            radius = R[i]
            circle = circle_list[i]
            color = assigned_colors[i]
            plt.fill(circle[:, 0], circle[:, 1], facecolor=color, antialiased=False)

        plt.savefig(str(outdir) + 'circle_packing_' + time_stamp + '.png', dpi=350)

        arr = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        arr = arr.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        np.save(str(outdir) + 'circle_packing_' + time_stamp + '.npy', arr)
    if args.graph:
        import networkx as nx

        # For each circle, check all other circles to see if there is overlap.
        # Will need to track which are which, somehow use the indicies for this.

        search_radii = R + r_min
        search_array = np.hstack((C, R))

        empty = np.array([np.nan for i in range(0, C.shape[0])])
        empty = empty.reshape(empty.size, 1)
        search_array = np.hstack((search_array, empty))

        neighbors_list = []  # This is a list of arrays of the neighbors of circles at index i

        for i in range(0, search_array.shape[0]):
            search_array[:, 3] = np.nan  # Clear the column
            search_x = C[i, 0]
            search_y = C[i, 1]
            search_radius = search_radii[i]

            neighbors = circle_intersect_check(search_radius, search_array[:, 2], search_x, search_array[:, 0],
                                               search_y, search_array[:, 1])
            neighbors_list.append(neighbors)

        adjacency_matrix = np.zeros((search_array.shape[0], search_array.shape[0]))

        for i in range(0, len(neighbors_list)):
            neighbors = neighbors_list[i][0]

            for n in neighbors:
                if i == n:
                    continue

                else:
                    adjacency_matrix[i, n] = 1
                    adjacency_matrix[n, i] = 1

        position_dict = dict()

        for i in range(0, C.shape[0]):
            position_dict[i] = C[i, :]

        graph = nx.from_numpy_matrix(adjacency_matrix)

        plt.clf()
        fig, ax = plt.subplots(figsize=(8, 8))

        nx.draw(graph, pos=position_dict, node_size=17, node_color=assigned_colors, with_labels=False)
        plt.savefig(str(outdir) + 'graph_' + time_stamp + '.png', dpi=350)
        # Save the adjacency matrix and the positions of the nodes.
        np.save(str(outdir) + 'C_' + time_stamp + '.npy', C)
        np.save(str(outdir) + 'A_' + time_stamp + '.npy', adjacency_matrix)

    if args.voronoi:
        from scipy.spatial import Voronoi, voronoi_plot_2d

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))

        vor = Voronoi(C)
        regions, vertices = voronoi_finite_polygons_2d(vor)
        # colorize
        for r in range(0, len(regions)):
            region = regions[r]
            polygon = vertices[region]
            plt.fill(*zip(*polygon), alpha=0.7, facecolor=assigned_colors[r], edgecolor='k')

        ax.axis('equal')
        ax.axis('off')
        ax.set(xlim=(0 - 0.05 * ab[0], ab[0] + 0.05 * ab[0]), ylim=(0 - 0.05 * ab[1], ab[1] + 0.05 * ab[1]))
        plt.savefig(str(outdir) + 'voronoi_' + time_stamp + '.png', dpi=350)
