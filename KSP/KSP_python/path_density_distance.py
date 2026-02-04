import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import networkx as nx
import os
from scipy.optimize import curve_fit

# read the graph from the file
def straight_path_dist(X, Y, s, t):
        # Compute distance from every (X, Y) point to s and t
    d_sv = np.sqrt((X - s[0])**2 + (Y - s[1])**2)
    d_tv = np.sqrt((X - t[0])**2 + (Y - t[1])**2)
    
    return d_sv + d_tv

def point_to_line_distance(point, line_start, line_end):
    """Calculate the distance from a point to a line segment defined by two points."""
    line_vec = np.array(line_end) - np.array(line_start)
    point_vec = np.array(point) - np.array(line_start)
    line_len = np.linalg.norm(line_vec)
    
    if line_len == 0:
        return np.linalg.norm(point_vec)  # Line segment is a point
    
    line_unitvec = line_vec / line_len
    projection_length = np.dot(point_vec, line_unitvec)
    
    if projection_length < 0:
        closest_point = line_start
    elif projection_length > line_len:
        closest_point = line_end
    else:
        closest_point = line_start + projection_length * line_unitvec
    
    return np.linalg.norm(np.array(point) - closest_point)

def points_to_line_distance(x, y, line_start, line_end):
    """
    Calculate the distance from points (x, y) to a line segment defined by line_start and line_end.
    Supports x, y as scalars or NumPy arrays.
    """
    # Convert inputs to arrays
    x = np.asarray(x)
    y = np.asarray(y)
    points = np.stack((x, y), axis=-1)  # shape: (..., 2)

    # Vectors
    a = np.array(line_start)
    b = np.array(line_end)
    ab = b - a
    ab_len_sq = np.dot(ab, ab)

    # Vector from line_start to each point
    ap = points - a

    # Projection factor t of each point onto the line
    t = np.clip(np.dot(ap, ab) / ab_len_sq, 0.0, 1.0)
    closest = a + np.outer(t, ab) if t.ndim > 0 else a + t * ab

    # Distance from point to closest point on segment
    dist = np.linalg.norm(points - closest, axis=-1)
    return dist

def number_of_paths(x, gamma):#, gamma2):
    """
    coords: flattened input of shape (N, 4) where each row is [x, y, s, t]
    gamma: fitting parameter
    """
    sp_s_tot, sp_t_tot, d = x[0], x[1], x[2]
    eps=0.1

    # d = np.linalg.norm(t - s)  # shape: (N,)
    #d_sv = np.sqrt((x - s[0])**2 + (y - s[1])**2)
    d_sv = np.array([sp_s_tot[v] for v in range(len(sp_s_tot))])
    #d_tv = np.sqrt((x - t[0])**2 + (y - t[1])**2)
    d_tv = np.array([sp_t_tot[v] for v in range(len(sp_t_tot))])
    diff = d * (1 + eps) - (d_sv + d_tv)
    # norm = (gamma*d*eps+1)*np.exp(gamma*d*eps)
    norm = np.exp(gamma*d*eps)
    # res = (gamma*diff+1)*np.exp(gamma*diff)
    res = (gamma*diff-1)*np.exp(gamma*diff)+1
    out = res/norm
    out[diff < 0] = 0  # Set negative differences to zero
    return out


def gamma_dist(x, a):
    gamma = a*x**0.7
    return gamma

def number_of_paths_2(coords, gamma):
    eps=0.6
    """
    coords: flattened input of shape (N, 4) where each row is [x, y, s, t]
    gamma: fitting parameter
    """
    # x, y, s, t = coords[:, 0], coords[:, 1], np.array([coords[:, 2][0], coords[:, 3][0]]), np.array([coords[:, 4][0], coords[:, 5][0]])
    # d = np.linalg.norm(t - s)  # shape: (N,)
    d = sp_s[t_ind]
    n = len(sp_s)
    d_prim = d * (1 + eps)
    #d_sv = np.sqrt((x - s[0])**2 + (y - s[1])**2)
    d_sv = np.array([sp_s[v] for v in coords])
    #d_tv = np.sqrt((x - t[0])**2 + (y - t[1])**2)
    d_tv = np.array([sp_t[v] for v in coords])
    diff = d_prim - (d_sv + d_tv)
    norm = np.exp(gamma*d_prim)*(d*eps-1/gamma)/gamma + np.exp(gamma*d)/gamma**2
    res = np.exp(gamma*d_prim)*(diff-1/gamma)/gamma + np.exp(gamma*(d_sv+d_tv))/gamma**2
 
    out = res/norm
    out[diff < 0] = 0  # Set negative differences to zero
    return out


G = nx.Graph()
with open(f'../KSP_cpp/100n_delaunays/delaunay_0.csv', 'r') as f:
    lines = f.readlines()
    # read the path and sum the instances of each node
    tot_lines = len(lines)
    for line in lines:
        path = line.strip().split(',')
        G.add_edge(int(path[0]), int(path[1]), weight=float(path[2]))

N = len(G) # Assuming 200 nodes as per the original code
global sp_s, sp_t
distances = []
fig, ax = plt.subplots(figsize=(10, 6))

# Arrays recording data from all pairs for final fit
dists_s = np.ones(0)
dists_t = np.ones(0)
node_occurr_tot = np.ones(0)
dists_st = np.ones(0)
eps_ls = np.ones(0)

node_x = np.zeros(N)
node_y = np.zeros(N)
with open('../KSP_cpp/100n_delaunays/pos_delaunay_0.txt', 'r') as f:
    it = 0
    lines = f.readlines()
    for line in lines:
        parts = line.strip().split(',')
        node_id = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        # Store the position of each node
        node_pos = [x, y]
        node_x[node_id] = x
        node_y[node_id] = y

for file in os.listdir('../KSP_cpp/100n_delaunays/'):
    node_occurr = np.zeros(N)  
    if file.startswith('paths'):
        it += 1
        with open(f'../KSP_cpp/100n_delaunays/'+file, 'r') as f:
            lines = f.readlines()
            # read the path and sum the instances of each node
            tot_lines = len(lines)
            for line in lines:
                path = line.strip().split(',')
                s_ind = int(path[0])
                t_ind = int(path[-1])
                for node in path:
                    node_occurr[int(node)] += 1
        
        with open(f'../KSP_cpp/100n_delaunays/dists_{s_ind}_{t_ind}.csv') as f:
            lines = f.readlines()
            # read the path and sum the instances of each node
            for i,line in enumerate(lines):
                d_sp = float(line)
                d_max = float(line)
            

    else:
        continue
# s_ind = 24
# t_ind = 180

# s_ind = 47
# t_ind = 91

    sp_s = nx.shortest_path_length(G, source=s_ind, target=None, weight='weight')
    sp_t = nx.shortest_path_length(G, source=t_ind, target=None, weight='weight')

    sp_s_arr = np.zeros(N)
    sp_t_arr = np.zeros(N)
    for node_id in range(N):
        sp_s_arr[node_id] = sp_s[node_id]
        sp_t_arr[node_id] = sp_t[node_id]

    dists_s = np.append(dists_s, np.copy(sp_s_arr))
    dists_t = np.append(dists_t, np.copy(sp_t_arr))
    dists_st = np.append(dists_st, np.ones(N)*sp_s[t_ind])
    eps_ls = np.append(eps_ls, np.ones(N)*d_sp/d_max)


# with open(f'./100n_delaunays/paths_{s_ind}_{t_ind}.csv', 'r') as f:
#     lines = f.readlines()
#     # read the path and sum the instances of each node
#     tot_lines = len(lines)
#     for line in lines:
#         path = line.strip().split(',')
#         for node in path:
#             node_occurr[int(node)] += 1

    # Normalize the occurrences
    node_occurr /= tot_lines
    node_occurr_tot= np.append(node_occurr_tot, np.copy(node_occurr))


    # Calculate distances for each node
    node_dists = np.zeros(N)
    for node_id in range(N):
        node_dists[node_id] = sp_s[node_id] + sp_t[node_id]  # path distance from s to t through the node
        if node_id == s_ind:
            s = np.array([node_x[node_id],node_y[node_id]])
        elif node_id == t_ind:
            t = np.array([node_x[node_id],node_y[node_id]])

        # node_dists[node_id] = straight_path_dist(node_pos[0], node_pos[1], s, t)

    # return np.exp(gamma1*(d * (1 + eps) - (d_sv + d_tv)))-np.exp(gamma2*(d * (1 + eps) - (d_sv + d_tv)))




# coords = np.arange(N)

    # Fit the model

# distances = straight_path_dist(coords[:, 0],coords[:, 1], sources[0]["center"], sources[1]["center"])
    distances += [sp_s[i]+ sp_t[i] - sp_s[t_ind] for i in range(N)]
    if it == 1:
        ax.scatter(node_dists-sp_s[t_ind], node_occurr, marker='o', c = 'orange', s=3,
            label = 'Observed')
    else:
        ax.scatter(node_dists-sp_s[t_ind], node_occurr, marker='o', c = 'orange', s=3)

popt, pcov = curve_fit(number_of_paths, [dists_s, dists_t, dists_st], node_occurr_tot, [15]) 

ax.set_xlabel(r"Distance $d_{svt} - d_{st}$", fontsize = 20)
ax.set_ylabel(r"$\sigma_{st}^{\epsilon}(v)/\sigma_{st}^{\epsilon}$", fontsize = 20)
print("Fitted gamma:", popt[0], 'exp. gamma', np.log(len(G.edges)/(N-1)))#, "Fitted gamma 2:", popt[1])

pred = number_of_paths([dists_s, dists_t, dists_st], popt[0])  # Use the fitted gamma value
r2_score = r2_score(node_occurr_tot, pred)
ax.plot(distances, pred, marker='o', markersize=1, linestyle='None', 
        label = r'Predicted, $R^2$ = '+f'{r2_score:.2f}')#*popt

ax.set_xlim(0,2)
# ax.set_yscale('log')
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()


# fig, ax = plt.subplots(figsize=(10, 6))
# ax.set_xlabel(r"Predicted fraction of paths", fontsize = 20)
# ax.set_ylabel(r"Experimental fraction of paths", fontsize = 20)

# ax.scatter(pred, node_occurr, marker='o', c = 'blue', s=10)
# ax.plot(sorted(pred), sorted(pred), c = 'black', ls = '--',
#            label = 'y = x')
# #ax.set_yscale('log')
# plt.legend(fontsize=16)
# plt.tight_layout()
# plt.show()