import numpy as np
import networkx as nx
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
# Create N delaunays
seed = 43
np.random.seed(seed)
totit = 1

for it in range(totit):

    nnodes = 100 # number of nodes in the graph

    # x,y coords of points
    clust = np.random.uniform(-1, 1, size=(nnodes, 2))

    # for each Delaunay triangle
    delaunay = nx.Graph()
    n_ls = list(range(0, nnodes))
    delaunay.add_nodes_from(n_ls)

    # make a graph based on the Delaunay triangulation edges

    # dictionary of node:position
    pos = dict(zip(range(len(clust)), clust))
    with open('/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/'+\
            f'KSP_cpp/100n_delaunays/pos_delaunay_{it}.txt', 'w') as f:
        for key, value in pos.items():
            f.write(f"{key},{value[0]},{value[1]}\n")
    

    nx.set_node_attributes(delaunay, pos, name='pos')

    tri2 = Delaunay(list(pos.values()))

    for path in tri2.simplices:
        nx.add_path(delaunay, path)

    edge_dists = {edge : np.linalg.norm(np.array(pos[edge[1]]) - 
             np.array(pos[edge[0]])) for edge in delaunay.edges}

    nx.set_edge_attributes(delaunay, edge_dists, name = 'weight')
    # save the graph   
    nx.write_weighted_edgelist(delaunay, '/Users/robertbenassai/Documents/UOC/k_shortest_path_betweenness/'+\
                               f'KSP_cpp/100n_delaunays/delaunay_{it}.csv', 
                               delimiter=',')#, header='vertex1,vertex2,edge_weight')
    nx.draw_networkx(delaunay, pos=pos, with_labels=False, 
                     node_size=10, font_size=8, width=0.1)
    # plt.scatter([pos[47][0], pos[91][0]],[pos[47][1], pos[91][1]], color='green', label='Parella 1')
    # plt.scatter([pos[24][0], pos[180][0]], [pos[24][1], pos[180][1]], color='red', label='Parella 2')
    plt.legend()
    plt.show()

# x = np.linspace(-1, 1, 200)
# y = np.linspace(-1, 1, 200)

# X, Y = np.meshgrid(x, y)

# def number_of_paths(X, Y, s, t, gamma1, eps=0.1):
#     """Calculate the number of paths based on a Gaussian-like decay."""
#     d = np.linalg.norm(t - s)  # Scalar distance between s and t

#     # Compute distance from every (X, Y) point to s and t
#     d_sv = np.sqrt((X - s[0])**2 + (Y - s[1])**2)
#     d_tv = np.sqrt((X - t[0])**2 + (Y - t[1])**2)

#     return np.exp(d * (1 + eps) - (d_sv + d_tv))

# Z = number_of_paths(X, Y, np.array([-0.25, -0.25]), np.array([0.25, 0.25]), eps=0.1)

# plt.figure(figsize=(8, 6))
# plt.contourf(X, Y, Z, levels=50, cmap='viridis')
# plt.colorbar(label='Number of Paths')
# plt.scatter(-0.25, -0.25, color='red', label='Source (s)')
# plt.scatter(0.25, 0.25, color='blue', label='Target (t)')
# plt.show()