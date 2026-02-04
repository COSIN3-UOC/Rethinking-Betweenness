import igraph
import numpy as np
import warnings
from tqdm import tqdm
import os
import networkx as nx
import argparse
import pandas as pd
import re
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.io import savemat

import igraph as ig
import math

def all_sp_paths_from_s(g, s, weights=None, mode="OUT"):
    # 1) One-source distances
    dist = g.shortest_paths(source=s, weights=weights, mode=mode)[0]

    # 2) Build predecessor list: pred[v] = {u s.t. u→v is on some shortest path}
    pred = [[] for _ in range(g.vcount())]
    directed = g.is_directed()
    for e in g.es:
        u, v = e.tuple
        w = e[weights] if (weights and weights in e.attributes()) else 1.0

        def consider(a, b):  # edge a->b
            if math.isfinite(dist[a]) and math.isfinite(dist[b]) and abs(dist[a] + w - dist[b]) < 1e-12:
                pred[b].append(a)

        if directed:
            consider(u, v)
        else:
            # undirected edge acts both ways
            consider(u, v)
            consider(v, u)

    # 3) Backtrack all shortest paths s→t using the predecessor DAG
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def enumerate_to(t):
        if t == s:
            return [[s]]
        paths = []
        for u in pred[t]:
            for p in enumerate_to(u):
                paths.append(p + [t])
        return paths

    return enumerate_to  # call enumerate_to(t) to get all paths s→t



# def route_cost_nx(G, route_nx, attribute):
#     return sum(G[u][v][0][attribute] for u, v in zip(route_nx[:-1], route_nx[1:]))
    

def node_to_edge_list_ig(G, route):
    return [G.get_eid(route[i], route[i + 1]) for i in range(len(route) - 1)]
    

def route_cost_ig(G, route, attribute):

    path_edges = node_to_edge_list_ig(G, route)

    # Sum the weights of these edges
    route_cost = sum(G.es[path_edges][attribute])

    return route_cost


def check_if_connected(G, node_src, node_dest):

    # dict_nx_to_ig = G["info"]["node_nx_to_ig"]

    # node_src_ig = dict_nx_to_ig[node_src]
    # node_dest_ig = dict_nx_to_ig[node_dest]

    return len(G.get_shortest_paths(node_src, node_dest, mode="OUT", output="vpath")[0])>0


def path_penalization(G_ig, node_src, node_dest, k, p, attribute, all_distinct=True, 
                      remove_tmp_attribute=True, max_iter=1e3, dict_ig_to_nx=None, eps=0.0, 
                      shortest_paths_st=None):

    # dict_ig_to_nx dict to map IG ids back into the "original" ones (that is NetworkX)
    
    # Initialize iteration counter and result containers
    it = 0
    result_list, path_list = [], []
    
    # Create a temporary copy of the edge attribute to penalize
    G_ig.es[f"tmp_{attribute}"] = G_ig.es[attribute]

    # dict_ig_to_nx = G_ig["info"]["node_ig_to_nx"]
    # dict_nx_to_ig = G_ig["info"]["node_nx_to_ig"]

    # node_src_ig = dict_nx_to_ig[node_src]
    # node_dest_ig = dict_nx_to_ig[node_dest]
    
    # Iterate to find k distinct paths or until max_iter is reached
    # sp = G_ig.get_shortest_paths(node_src, node_dest, weights=attribute, mode="OUT", output="vpath")[0];
    if shortest_paths_st is not None:
        if np.array(shortest_paths_st).ndim == 1:
            sp = shortest_paths_st
        elif np.array(shortest_paths_st).ndim == 2:
            sp = shortest_paths_st[0]
    else:
        sp = G_ig.get_shortest_paths(node_src, node_dest, weights=attribute, mode="OUT", output="vpath")[0];
    
    sp_cost = route_cost_ig(G_ig, sp, attribute)
    original_cost = sp_cost if sp_cost > 0 else 1e-10  # avoid division by zero
    while len(result_list) < k and it < max_iter:# and original_cost <= sp_cost* (1 + eps):
        
        # Compute the shortest path using the temporary penalized attribute
        if len(result_list) == 0 or len(result_list) < len(shortest_paths_st): #do not compute the shortest path again if already computed

            if shortest_paths_st is not None:
                sp_k = shortest_paths_st[it] 
            else:
                sp_k = G_ig.get_shortest_paths(node_src, node_dest, weights=f"tmp_{attribute}", mode="OUT", output="vpath")[0];
        else:
            sp_k = G_ig.get_shortest_paths(node_src, node_dest, weights=f"tmp_{attribute}", mode="OUT", output="vpath")[0];
        # If the path is distinct (or if duplicates are allowed), add it to the result list
        if not (all_distinct and sp_k in path_list):
            
            # Compute the original cost
            original_cost = route_cost_ig(G_ig, sp_k, attribute)
        
            # Compute the penalized cost
            penalized_cost = route_cost_ig(G_ig, sp_k, f"tmp_{attribute}")
                
            # Store the path details in the result list
            # store the route as sequence of edges
            result_list.append({"node_list_nx": [dict_ig_to_nx[v] for v in sp_k], # translate into networkx vertices
                                    "original_cost": original_cost,
                                    "penalized_cost": penalized_cost,
                                    "iteration": it})
            
            if len(result_list) == 1:
                sp_path_edges = node_to_edge_list_ig(G_ig, sp_k)
            elif len(result_list) > 1:
                sp_k_path_edges = node_to_edge_list_ig(G_ig, sp_k)

                # compute similarity with the first path
                similarity = len(set(sp_k_path_edges) & set(sp_path_edges)) / len(set(sp_path_edges))
                result_list[-1]["similarity_to_sp"] = similarity
    
            # Keep track of the path to ensure uniqueness
            path_list.append(sp_k)
    
        # apply the penalization to the current sp's edges
        edge_list_sp_k = node_to_edge_list_ig(G_ig, sp_k)
        for e in G_ig.es(edge_list_sp_k):
            e[f"tmp_{attribute}"] *= (1 + p)
    
        # Increment the iteration counter
        it += 1
    
    # Check if the maximum number of iterations was reached
    # if it == max_iter:
    #    print(f'Iterations limit reached, returned {len(result_list)} distinct paths (instead of {k}).')
    #    warnings.warn(f'Iterations limit reached, returned {len(result_list)} distinct paths (instead of {k}).', RuntimeWarning)
        
    # Remove the temporary attribute from the graph if required
    if remove_tmp_attribute:
        del(G_ig.es[f"tmp_{attribute}"])
    
    result_list = sorted(result_list, key=lambda x: x["original_cost"])
    # throw away paths that have original cost longer than (1+eps)*SP
    filtered_result_list = []
    for path_info in result_list:
        if path_info["original_cost"] <= sp_cost * (1 + eps):
            filtered_result_list.append(path_info)


    return filtered_result_list


def compute_penalized_bwss(G_ig, k, p, attribute, all_distinct=True, remove_tmp_attribute=True,
                           max_iter=1e3, directed=False, norm = False, eps=0.0, dict_ig_to_nx=None, list_paths = False):
    
    if eps > 1:
        eps = eps/100  # convert percentage to fraction
    elif eps < 0:
        raise ValueError("eps must be non-negative.")
    
    N = len(G_ig.vs)
    max_len_paths = 2*int(np.sqrt(N))
    bwss_penalized = np.zeros(N)
    total = 0
    no_paths_found = 0
    tolerance_ls = []
    similarity_ls = []
    distance_ls = []
    n_paths = []
    sp_len_ls = []

    path_list_all = []  # to store all paths if list_paths is True
    for _ in range(N):
        path_list_all.append([])
        for _ in range(N):
            path_list_all[-1].append([])
    
    with tqdm(total=N*(N-1)/2) as pbar:  # or set a known upper bound if you have one
        for node_src_ig in range(len(G_ig.vs)):
            all_sh_paths = all_sp_paths_from_s(G_ig, node_src_ig, weights=attribute, mode="OUT")
            for node_dest_ig in range(len(G_ig.vs)) if directed else range(node_src_ig+1, len(G_ig.vs)):
                if node_src_ig != node_dest_ig:
                    shortest_paths_st = all_sh_paths(node_dest_ig)  # to cache paths
                    res_pp = path_penalization(G_ig, node_src_ig, node_dest_ig,
                                                k, p, attribute, all_distinct=all_distinct, 
                                                remove_tmp_attribute=remove_tmp_attribute, 
                                                max_iter=max_iter, dict_ig_to_nx=dict_ig_to_nx, eps=eps, shortest_paths_st=shortest_paths_st)
                    # accumulate the contribution to each edge in each path
                    sp_dist = res_pp[0]["original_cost"] if len(res_pp) > 0 else None

                    valid_paths = 0
                    valid_p_arr = [1 if path_info["original_cost"]<= sp_dist*(1+eps) else 0 for path_info in res_pp]

                    tot_valid_paths = sum(valid_p_arr)
                    
                    for path_info in res_pp:
                        node_list = path_info["node_list_nx"]
                        distance = path_info["original_cost"]
                        
                        if distance > sp_dist:
                            similarity = path_info["similarity_to_sp"]

                        if distance <= sp_dist*(1+eps):

                            for v in node_list:
                                # ignore source and target nodes
                                if v != dict_ig_to_nx[node_src_ig] and v != dict_ig_to_nx[node_dest_ig]:
                                    bwss_penalized[v] += 1 / tot_valid_paths

                            if distance > sp_dist:
                                tolerance_ls.append((distance - sp_dist)/sp_dist * 100)  # in percentage
                                similarity_ls.append(similarity)
                                distance_ls.append(distance)
                                        
                            valid_paths += 1

                            if list_paths:
                                if len(node_list) == 0:
                                    raise ValueError("Empty path encountered.")
                                final_list = node_list+np.zeros(max_len_paths-len(node_list)).tolist()+[distance]  # pad to avoid issues
                                path_list_all[dict_ig_to_nx[node_src_ig]][dict_ig_to_nx[node_dest_ig]].append(final_list)
                                if directed == False:
                                    final_list = node_list[::-1]+np.zeros(max_len_paths-len(node_list)).tolist()+[distance]  # pad to avoid issues
                                    path_list_all[dict_ig_to_nx[node_dest_ig]][dict_ig_to_nx[node_src_ig]].append(final_list)

                    if valid_paths == 1:
                        if valid_paths == tot_valid_paths:
                            no_paths_found += 1
                        else: 
                            raise ValueError("Inconsistent number of valid paths found.")
                    elif valid_paths == 0:
                        raise ValueError("No valid paths found, but nodes are connected.")
                    
                    n_paths.append(valid_paths)
                    sp_len_ls.append(sp_dist)
                    
                total += 1
                pbar.update(1)
        pbar.total = total 

    if norm:
        if directed:
            bwss_penalized = bwss_penalized / ((N-1)*(N-2))
            print("fraction of pairs without QSPs", no_paths_found/((N)*(N-1)))
            print("average tolerance (%)", np.mean(tolerance_ls), np.std(tolerance_ls), np.any(np.array(tolerance_ls)>eps*100))


        else:
            bwss_penalized = 2*bwss_penalized / ((N-1)*(N-2))
            print("fraction of pairs without QSPs", 2*no_paths_found/((N)*(N-1)))
            print("average tolerance (%)", np.mean(tolerance_ls), np.std(tolerance_ls), np.any(np.array(tolerance_ls)>eps*100))

    if list_paths:
        return(bwss_penalized, similarity_ls, tolerance_ls, distance_ls, n_paths, sp_len_ls, path_list_all)
    else:
        return(bwss_penalized, similarity_ls, tolerance_ls, distance_ls, n_paths, sp_len_ls)

def get_graph_number(fname):
    num = ''
    for char in fname:
        if char.isdigit():
            num += char
    return int(num)


def build_upper_tri_cell_from_dict(idx2mat, M, N_expected=None, varname="X", out_path="upper_tri_cell.mat"):
    """
    Build an MxM MATLAB cell array where cells (i,j) for i<j hold K_ij x N matrices,
    and diagonal/lower cells are empty (0 x N). Save to .mat.

    Parameters
    ----------
    idx2mat : dict[(int,int) -> np.ndarray]
        Keys are (i,j) with 0 <= i < j < M. Values are K_ij x N numeric arrays.
    M : int
        Size of the square cell array.
    N_expected : int or None
        If given, enforce that every matrix has N columns. If None, infer from first.
    varname : str
        Variable name in the MAT-file.
    out_path : str or path-like
        Output .mat path.
    """
    # Infer N if needed
    if N_expected is None:
        for (i, j), A in idx2mat.items():
            A = np.asarray(A)
            if A.ndim != 2 or A.size == 0:
                continue
            N_expected = A.shape[1]
            break
        if N_expected is None:
            raise ValueError("Cannot infer N (no non-empty 2D arrays); pass N_expected.")

    # Create object array (MATLAB cell)
    cell = np.empty((M, M), dtype=object)

    # Fill with empties first (0 x N) so MATLAB gets [] with correct width
    empty_row = np.empty((0, N_expected), dtype=float)
    for i in range(M):
        for j in range(M):
            cell[i, j] = empty_row

    # Place upper-triangle matrices
    for (i, j), A in idx2mat.items():
        if not (0 <= i < j < M):
            raise ValueError(f"Invalid key {(i,j)}; require 0 <= i < j < {M}")
        A = np.asarray(A, dtype=float)
        if A.ndim != 2 or A.shape[1] != N_expected:
            raise ValueError(f"(i={i}, j={j}) must be Kx{N_expected}, got {A.shape}")
        cell[i, j] = A

    sio.savemat(out_path, {varname: cell}, oned_as="row")
    return cell



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Penalized Betweenness-Weighted Shortest Paths (BWSS) for graphs.")
    parser.add_argument("graph_type", type=str, help="ring, lattice or star)")
    parser.add_argument("speed_ratio", type=str, help="the speed ratio between fast and slow roads (1, 2, 3)")
    parser.add_argument("max_k", type=int, default=10, help="maximum number of paths to compute (default: 10)")
    parser.add_argument("p", type=float, default=10, help="penalization ([0,1] will be times 1+p)")
    parser.add_argument("eps", type=float, default=0.3, help="tolerance for path length increase (default: 0.3)")
    args = parser.parse_args()

    graph_type = args.graph_type
    speed_ratio = args.speed_ratio
    max_k = args.max_k  # number of paths to compute
    p = args.p
    eps = args.eps

    if graph_type == 'lattice':
        folder_path = f'../../graphs_qspbwss/lattice_roads/vx{int(speed_ratio)}/graph_variations/'
    elif graph_type == 'star':
        folder_path = f'../../graphs_qspbwss/star_roads/vx{int(speed_ratio)}/graph_variations/'
    elif graph_type == 'ring':
        folder_path = f'../../graphs_qspbwss/ring_roads/vx{int(speed_ratio)}/graph_variations/'
    elif graph_type == 'other':
        folder_path = input("Enter the folder path containing the graphs: ")
        nickname = input("Enter a nickname for this set of graphs (used for saving results): ")

    iter = 0
    for graph_fname in os.listdir(folder_path):
        if graph_fname.endswith('.csv'):
            if graph_type != 'other':
                graph_num = get_graph_number(graph_fname)
                
            graph_path = os.path.join(folder_path, graph_fname)
            df = pd.read_csv(graph_path)  # header handled automatically
            G = nx.from_pandas_edgelist(
                df, source="vertex1", target="vertex2",
                edge_attr="edge_weight", create_using=nx.Graph())

            G_ig = igraph.Graph.from_networkx(G)

            dict_ig_to_nx = {v.index: v["_nx_name"] for v in G_ig.vs}
            
            G_ig.es["time"] = G_ig.es["edge_weight"]
            print(f'Computing penalized bwss for graph {graph_fname}...')

            bwss_penalized, simil, tol_ls, ds_ls, n_paths, sp_len_ls, path_mat =\
                  compute_penalized_bwss(G_ig, k=max_k, p=p, attribute='time', all_distinct=True, 
                                        remove_tmp_attribute=True, max_iter=10000, directed=False, 
                                        norm=True, eps = eps, dict_ig_to_nx=dict_ig_to_nx, list_paths=True)
            # save to file

            if iter == 0:
                simil_tot = simil
                tol_ls_tot = tol_ls
                ds_ls_tot = ds_ls
                n_paths_tot = n_paths
                sp_len_ls_tot = sp_len_ls
            else:
                simil_tot += simil
                tol_ls_tot += tol_ls
                ds_ls_tot += ds_ls
                n_paths_tot += n_paths
                sp_len_ls_tot += sp_len_ls
            iter += 1

            p_str = str(p)
            if p<1:
                p_str = p_str[0]+p_str[2:]


            e_str = str(eps)
            if eps<1:
                e_str = e_str[0]+e_str[2:]

            if graph_type == 'other':
                out_path = f'./penalized_bwss/other/{nickname}/{e_str}_tol/penbwss_k{max_k}_p{p_str}.csv'
            else:
                out_path = f'./penalized_bwss/{graph_type}/vx{speed_ratio}/{e_str}_tol/k{max_k}/p{p_str}/{graph_num}_penbwss_k{max_k}.csv'

            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            np.savetxt(out_path, bwss_penalized, delimiter=',')

            idx2mat = {(i,j): path_mat[i][j] for i in range(len(path_mat)) for j in range(i+1, len(path_mat)) if len(path_mat[i][j])>0}
            # print(len(idx2mat) == len(G_ig.vs)*(len(G_ig.vs)-1)/2, len(G_ig.vs), len(idx2mat), "pairs have penalized paths.")
            N = len(G_ig.vs)+1
            M = N-1

            if graph_type == 'other':
                out_path_paths = f'../../../MCM/mcm/other_graphs/{nickname}/pen_path_ls_t_{e_str}_k_{max_k}_p_{p_str}.mat'
                os.makedirs(os.path.dirname(out_path_paths), exist_ok=True)
            else:
                out_path_paths = f'../../../MCM/mcm/{graph_type}_road_p_structure_vx{speed_ratio}/ring_road_{graph_num}/pen_path_ls_t_{e_str}_k_{max_k}_p_{p_str}.mat'
            try:
                cell = build_upper_tri_cell_from_dict(idx2mat, M, N_expected=N,
                                      varname="mycell", out_path=out_path_paths)
            except FileNotFoundError:
                continue

    # # plot tolerance vs similarity
    # plt.figure(figsize=(8,6))
    # plt.scatter(tol_ls_tot, simil_tot, alpha=0.5, s=10)
    # plt.xlabel('Tolerance (%)', fontsize=14)
    # plt.ylabel('Similarity to SP', fontsize=14)
    # plt.title(f'Penalized BWSS paths similarity vs tolerance (p={p}, eps={eps})', fontsize=16)
    # plt.tight_layout()
    # plt.show()

    # #plot distance vs similarity
    # plt.figure(figsize=(8,6))
    # plt.scatter(ds_ls_tot, simil_tot, alpha=0.5, s=10)
    # plt.xlabel('Path Distance', fontsize=14)
    # plt.ylabel('Similarity to SP', fontsize=14)
    # plt.title(f'Penalized BWSS paths similarity vs distance (p={p}, eps={eps})', fontsize=16)
    # plt.tight_layout()
    # plt.show()

    # plt.figure(figsize=(8,6))
    # plt.scatter(sp_len_ls_tot, n_paths_tot, alpha=0.3, s=10)
    # plt.xlabel('SP Length', fontsize=14)
    # plt.ylabel('Number of Penalized BWSS paths', fontsize=14)
    # plt.title(f'Number of Penalized BWSS paths vs SP length (p={p}, eps={eps})', fontsize=16)
    # plt.tight_layout()
    # plt.show()  

    #plot histogram of similarity
    # plt.figure(figsize=(8,6))
    # plt.hist(simil_tot, bins=30, alpha=0.7, color='blue', edgecolor='black')
    # plt.xlabel('Similarity to SP', fontsize=14)
    # plt.ylabel('Frequency', fontsize=14)
    # plt.title(f'Histogram of Penalized BWSS paths similarity to SP (p={p}, eps={eps})', fontsize=16)
    # plt.tight_layout()
    # plt.show()  

