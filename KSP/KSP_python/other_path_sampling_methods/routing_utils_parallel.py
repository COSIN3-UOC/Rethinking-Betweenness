import argparse
import math
import os
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import igraph
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
# scipy imports removed (MATLAB writers no longer used)
from tqdm import tqdm
import h5py


def node_to_edge_list_ig(G, route):
    return [G.get_eid(route[i], route[i + 1]) for i in range(len(route) - 1)]


def route_cost_ig(G, route, attribute):
    path_edges = node_to_edge_list_ig(G, route)
    return sum(G.es[path_edges][attribute])


def check_if_connected(G, node_src, node_dest):
    return len(G.get_shortest_paths(node_src, node_dest, mode="OUT", output="vpath")[0]) > 0


def path_penalization(
    G_ig,
    node_src,
    node_dest,
    k,
    p,
    attribute,
    all_distinct=True,
    remove_tmp_attribute=True,
    max_iter=1e3,
    dict_ig_to_nx=None,
    eps=0.0,
):
    it = 0
    result_list, path_list = [], []

    G_ig.es[f"tmp_{attribute}"] = G_ig.es[attribute]

    sp = G_ig.get_shortest_paths(node_src, node_dest, weights=attribute, mode="OUT", output="vpath")[0]
    sp_cost = route_cost_ig(G_ig, sp, attribute)
    original_cost = sp_cost if sp_cost > 0 else 1e-10

    sp_k = sp
    while len(result_list) < k and it < max_iter:# and original_cost <= sp_cost * (1 + eps):

        if not (all_distinct and sp_k in path_list):
            penalized_cost = route_cost_ig(G_ig, sp_k, f"tmp_{attribute}")

            if original_cost > 100000:
                Exception(f"Original cost too high, {original_cost}, something is wrong.")

            result_list.append(
                {
                    "node_list_nx": [dict_ig_to_nx[v] for v in sp_k],
                    "original_cost": original_cost,
                    "penalized_cost": penalized_cost,
                    "iteration": it,
                }
            )

            if len(result_list) == 1:
                sp_path_edges = node_to_edge_list_ig(G_ig, sp_k)
            elif len(result_list) > 1:
                sp_k_path_edges = node_to_edge_list_ig(G_ig, sp_k)
                similarity = len(set(sp_k_path_edges) & set(sp_path_edges)) / len(set(sp_path_edges))
                result_list[-1]["similarity_to_sp"] = similarity

            path_list.append(sp_k)

        edge_list_sp_k = node_to_edge_list_ig(G_ig, sp_k)
        for e in G_ig.es(edge_list_sp_k):
            e[f"tmp_{attribute}"] *= 1 + p

        sp_k = G_ig.get_shortest_paths(
            node_src, node_dest, weights=f"tmp_{attribute}", mode="OUT", output="vpath")[0]
        if not (all_distinct and sp_k in path_list):
            original_cost = route_cost_ig(G_ig, sp_k, attribute)
        it += 1

    if remove_tmp_attribute:
        del G_ig.es[f"tmp_{attribute}"]

    result_list = sorted(result_list, key=lambda x: x["original_cost"])

    # throw away paths that have original cost longer than (1+eps)*SP
    filtered_result_list = []
    for path_info in result_list:
        if path_info["original_cost"] <= sp_cost * (1 + eps):
            filtered_result_list.append(path_info)


    return filtered_result_list


def _process_source_batch(
    source_indices,
    graph,
    k,
    p,
    attribute,
    all_distinct,
    remove_tmp_attribute,
    max_iter,
    directed,
    eps,
    dict_ig_to_nx,
    list_paths,
):
    N = len(graph.vs)
    max_len_path_ls = 2*np.sqrt(N).astype(int)
    bwss_partial = np.zeros(N)
    tolerance_ls = []
    similarity_ls = []
    distance_ls = []
    n_paths = []
    sp_len_ls = []
    no_paths_found = 0
    pairs_processed = 0
    path_dict = {}

    for node_src_ig in source_indices:
        iterable = range(N) if directed else range(node_src_ig + 1, N)
        for node_dest_ig in iterable:
            if node_src_ig == node_dest_ig:
                continue

            res_pp = path_penalization(
                graph,
                node_src_ig,
                node_dest_ig,
                k,
                p,
                attribute,
                all_distinct=all_distinct,
                remove_tmp_attribute=remove_tmp_attribute,
                max_iter=max_iter,
                dict_ig_to_nx=dict_ig_to_nx,
                eps=eps,
            )

            if not res_pp:
                raise ValueError("No valid paths returned, but nodes are connected.")

            sp_dist = res_pp[0]["original_cost"]

            tot_valid_paths = len(res_pp)

            if tot_valid_paths == 0:
                raise ValueError("No valid paths found, but nodes are connected.")

            src_nx = dict_ig_to_nx[node_src_ig]
            dest_nx = dict_ig_to_nx[node_dest_ig]
            valid_paths = 0

            for path_info in res_pp:
                node_list = path_info["node_list_nx"]
                distance = path_info["original_cost"]

                similarity = path_info.get("similarity_to_sp")

                for v in node_list:
                    if v != dict_ig_to_nx[node_src_ig] and v != dict_ig_to_nx[node_dest_ig]:
                        bwss_partial[v] += 1 / tot_valid_paths

                if distance > sp_dist:
                    tolerance_ls.append((distance - sp_dist) / sp_dist * 100)
                    if similarity is not None:
                        similarity_ls.append(similarity)
                    distance_ls.append(distance)

                valid_paths += 1

                if list_paths:
                    if len(node_list) == 0:
                        raise ValueError("Empty path encountered.")
                    
                    # if undirected, store paths in (min, max) order to avoid duplicates
                    if not directed and src_nx > dest_nx:
                        path_dict.setdefault((dest_nx, src_nx), []).append(
                                            node_list[::-1] + [np.float32(distance)])
                    else:
                        path_dict.setdefault((src_nx, dest_nx), []).append(
                                            node_list + [np.float32(distance)])

            if valid_paths == 1:
                if valid_paths == tot_valid_paths:
                    no_paths_found += 1
                else:
                    raise ValueError("Inconsistent number of valid paths found.")
            elif valid_paths == 0:
                raise ValueError("No valid paths found, but nodes are connected.")

            n_paths.append(valid_paths)
            sp_len_ls.append(sp_dist)
            pairs_processed += 1

    return {
        "bwss_partial": bwss_partial,
        "tolerance_ls": tolerance_ls,
        "similarity_ls": similarity_ls,
        "distance_ls": distance_ls,
        "n_paths": n_paths,
        "sp_len_ls": sp_len_ls,
        "no_paths_found": no_paths_found,
        "pairs_processed": pairs_processed,
        "path_dict": path_dict,
    }


def compute_penalized_bwss_parallel(
    G_ig,
    k,
    p,
    attribute,
    all_distinct=True,
    remove_tmp_attribute=True,
    max_iter=1e3,
    directed=False,
    norm=False,
    eps=0.0,
    dict_ig_to_nx=None,
    list_paths=False,
    num_workers=None,
    chunk_size=1,
    path_consumer=None,
):
    if eps > 1:
        eps = eps / 100
    elif eps < 0:
        raise ValueError("eps must be non-negative.")

    if dict_ig_to_nx is None:
        raise ValueError("dict_ig_to_nx mapping must be provided.")

    N = len(G_ig.vs)
    if isinstance(dict_ig_to_nx, dict):
        dict_ig_to_nx = [dict_ig_to_nx[i] for i in range(N)]

    total_pairs = N * (N - 1) if directed else (N * (N - 1)) // 2

    bwss_penalized = np.zeros(N)
    tolerance_ls = []
    similarity_ls = []
    distance_ls = []
    n_paths = []
    sp_len_ls = []
    no_paths_found = 0
    if list_paths and path_consumer is not None and not callable(path_consumer):
        raise ValueError("path_consumer must be callable when provided.")

    path_list_all = defaultdict(list) if list_paths and path_consumer is None else None

    source_indices = list(range(N))
    if not directed:
        # ignore last node when chunking to avoid empty work units
        source_indices = source_indices[:-1]

    chunks = [source_indices[i : i + chunk_size] for i in range(0, len(source_indices), chunk_size)]

    with ProcessPoolExecutor(max_workers=num_workers) as executor, tqdm(total=total_pairs) as pbar:
        futures = [
            executor.submit(
                _process_source_batch,
                chunk,
                G_ig,
                k,
                p,
                attribute,
                all_distinct,
                remove_tmp_attribute,
                max_iter,
                directed,
                eps,
                dict_ig_to_nx,
                list_paths,
            )
            for chunk in chunks
        ]

        for future in as_completed(futures):
            result = future.result()
            bwss_penalized += result["bwss_partial"]
            tolerance_ls.extend(result["tolerance_ls"])
            similarity_ls.extend(result["similarity_ls"])
            distance_ls.extend(result["distance_ls"])
            n_paths.extend(result["n_paths"])
            sp_len_ls.extend(result["sp_len_ls"])
            no_paths_found += result["no_paths_found"]
            pbar.update(result["pairs_processed"])

            if list_paths and result["path_dict"]:
                if path_consumer is not None:
                    path_consumer(result["path_dict"])
                else:
                    for key, entries in result["path_dict"].items():
                        path_list_all[key].extend(entries)

    if norm:
        if directed:
            bwss_penalized = bwss_penalized / ((N - 1) * (N - 2))
            print("fraction of pairs without QSPs", no_paths_found / (N * (N - 1)))
            print("average tolerance (%)",
                np.mean(tolerance_ls),
                np.std(tolerance_ls),
                np.any(np.array(tolerance_ls) > eps * 100),
            )
        else:
            bwss_penalized = 2 * bwss_penalized / ((N - 1) * (N - 2))
            print("fraction of pairs without QSPs", 2 * no_paths_found / (N * (N - 1)))
            print("average tolerance (%)",
                np.mean(tolerance_ls),
                np.std(tolerance_ls),
                np.any(np.array(tolerance_ls) > eps * 100),
            )

    if list_paths:
        path_payload = dict(path_list_all) if path_list_all is not None else None
        return (
            bwss_penalized,
            similarity_ls,
            tolerance_ls,
            distance_ls,
            n_paths,
            sp_len_ls,
            path_payload,
        )
    return bwss_penalized, similarity_ls, tolerance_ls, distance_ls, n_paths, sp_len_ls


def get_graph_number(fname):
    num = ""
    for char in fname:
        if char.isdigit():
            num += char
    return int(num)
class IncrementalPathWriter:
    """Incrementally append OD-paths into the ragged HDF5 layout."""

    def __init__(self, out_path, N, varname="paths", dtype=np.float32):
        self._file = h5py.File(out_path, "w")
        self._dtype = dtype
        self._N = N
        self._pair_counts = defaultdict(int)

        root = self._file.create_group(varname)
        root.attrs["N"] = N

        self._paths_ds = root.create_dataset("paths", shape=(0,), maxshape=(None,), dtype=dtype, chunks=True)
        int_dtype = np.int32
        self._offsets_ds = root.create_dataset("offsets", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=True)
        self._lengths_ds = root.create_dataset("lengths", shape=(0,), maxshape=(None,), dtype=int_dtype, chunks=True)
        self._i_ds = root.create_dataset("i", shape=(0,), maxshape=(None,), dtype=int_dtype, chunks=True)
        self._j_ds = root.create_dataset("j", shape=(0,), maxshape=(None,), dtype=int_dtype, chunks=True)
        self._k_ds = root.create_dataset("k", shape=(0,), maxshape=(None,), dtype=int_dtype, chunks=True)

    def close(self):
        if self._file is not None:
            self._file.close()
            self._file = None

    def _append_meta(self, start, length, i_val, j_val, k_val):
        cur_size = self._offsets_ds.shape[0]
        new_size = cur_size + 1
        self._offsets_ds.resize((new_size,))
        self._lengths_ds.resize((new_size,))
        self._i_ds.resize((new_size,))
        self._j_ds.resize((new_size,))
        self._k_ds.resize((new_size,))
        self._offsets_ds[cur_size] = start
        self._lengths_ds[cur_size] = length
        self._i_ds[cur_size] = i_val
        self._j_ds[cur_size] = j_val
        self._k_ds[cur_size] = k_val

    def append_pair(self, i, j, paths):
        if not (0 <= i < j < self._N):
            raise ValueError(f"Invalid key {(i, j)}; require 0 <= i < j < {self._N}")
        if not paths:
            return

        k_base = self._pair_counts[(i, j)]
        for offset, path in enumerate(paths):
            arr = np.asarray(path, dtype=self._dtype)
            if arr.ndim != 1:
                raise ValueError("Each path must be a 1D sequence.")
            start = self._paths_ds.shape[0]
            new_total = start + len(arr)
            self._paths_ds.resize((new_total,))
            self._paths_ds[start:new_total] = arr
            self._append_meta(start, len(arr), i, j, k_base + offset)

        self._pair_counts[(i, j)] = k_base + len(paths)


def main():

    # input arguments
    parser = argparse.ArgumentParser(
        description="Compute Penalized Betweenness-Weighted Shortest Paths (BWSS) for graphs (parallel version)."
    )
    parser.add_argument("graph_type", type=str, help="ring, lattice or star)")
    parser.add_argument("speed_ratio", type=str, help="the speed ratio between fast and slow roads (1, 2, 3)")
    parser.add_argument("max_k", type=int, default=10, help="maximum number of paths to compute (default: 10)")
    parser.add_argument("p", type=float, default=10, help="penalization ([0,1] will be times 1+p)")
    parser.add_argument("eps", type=float, default=0.3, help="tolerance for path length increase (default: 0.3)")
    parser.add_argument("--workers", type=int, default=8, help="number of worker processes (default: os.cpu_count())")
    parser.add_argument("--chunk-size", type=int, default=1, help="number of sources per worker task")
    parser.add_argument("--list-paths", action="store_true", default = True, help="store all penalised paths, default: True")
    args = parser.parse_args()

    graph_type = args.graph_type
    speed_ratio = args.speed_ratio
    max_k = args.max_k
    p = args.p
    eps = args.eps
    num_workers = args.workers
    chunk_size = max(1, args.chunk_size)
    list_paths_flag = args.list_paths

    # determine directory of the graphs based on graph type
    if graph_type == "lattice":
        folder_path = f"../../graphs_qspbwss/lattice_roads/vx{int(speed_ratio)}/graph_variations/"
    elif graph_type == "star":
        folder_path = f"../../graphs_qspbwss/star_roads/vx{int(speed_ratio)}/graph_variations/"
    elif graph_type == "ring":
        folder_path = f"../../graphs_qspbwss/ring_roads/vx{int(speed_ratio)}/graph_variations/"
    elif graph_type == "other":
        folder_path = input("Enter the folder path containing the graphs: ")
        nickname = input("Enter a nickname for this set of graphs (used for saving results): ")
    else:
        raise ValueError("Invalid graph type specified.")

    # iterate over all graph files in the folder
    iter_idx = 0
    for graph_fname in os.listdir(folder_path):
        if not graph_fname.endswith(".csv"):
            continue

        if graph_type != "other":
            graph_num = get_graph_number(graph_fname)

        graph_path = os.path.join(folder_path, graph_fname)
        df = pd.read_csv(graph_path)
        G = nx.from_pandas_edgelist(
            df, source="vertex1", target="vertex2", edge_attr="edge_weight", create_using=nx.Graph()
        )
        # check if node labels are in default format (0, 1, 2, ..., N-1)
        if set(G.nodes) != set(range(len(G.nodes))):
            raise ValueError("Node labels in the graph must be integers from 0 to N-1.")
        
        G_ig = igraph.Graph.from_networkx(G)
        dict_ig_to_nx = {v.index: v["_nx_name"] for v in G_ig.vs}

        G_ig.es["time"] = G_ig.es["edge_weight"]
        print(f"Computing penalized bwss (parallel) for graph {graph_fname}...")
        print(f"maximum edge weight (time): {max(G_ig.es['time'])}")

        p_str = str(p)
        if p < 1:
            p_str = p_str.replace(".", "")

        e_str = str(eps)
        if eps < 1:
            e_str = e_str.replace(".", "")

        path_writer = None
        path_consumer_fn = None
        out_path_paths = None
        if list_paths_flag:
            num_nodes = len(G_ig.vs)
            if graph_type == "other":
                out_path_paths = (
                    f"../../../MCM/mcm/other_graphs/{nickname}/pen_path_ls_t_{e_str}_k_{max_k}_p_{p_str}.h5"
                )
            else:
                out_path_paths = (
                    f"../../../MCM/mcm/{graph_type}_road_p_structure_vx{speed_ratio}_test/ring_road_{graph_num}/"
                    f"pen_path_ls_t_{e_str}_k_{max_k}_p_{p_str}.h5"
                )
            os.makedirs(os.path.dirname(out_path_paths), exist_ok=True)
            path_writer = IncrementalPathWriter(out_path_paths, num_nodes, varname="paths", dtype=np.float32)

            def path_consumer_fn(chunk_dict, writer=path_writer):
                if not chunk_dict:
                    return
                for (i, j), entries in chunk_dict.items():
                    if not entries:
                        continue

                    if i == j:
                        continue
                    if i > j:
                        i, j = j, i
                    writer.append_pair(i, j, entries)

        try:
            results = compute_penalized_bwss_parallel(
                G_ig,
                k=max_k,
                p=p,
                attribute="time",
                all_distinct=True,
                remove_tmp_attribute=True,
                max_iter=10000,
                directed=False,
                norm=True,
                eps=eps,
                dict_ig_to_nx=dict_ig_to_nx,
                list_paths=list_paths_flag,
                num_workers=num_workers,
                chunk_size=chunk_size,
                path_consumer=path_consumer_fn,
            )
        finally:
            if path_writer is not None:
                path_writer.close()

        if list_paths_flag:
            (
                bwss_penalized,
                simil,
                tol_ls,
                ds_ls,
                n_paths,
                sp_len_ls,
                _,
            ) = results
        else:
            (
                bwss_penalized,
                simil,
                tol_ls,
                ds_ls,
                n_paths,
                sp_len_ls,
            ) = results

        if iter_idx == 0:
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
        iter_idx += 1

        if graph_type == "other":
            out_path = f"./penalized_bwss/other/{nickname}/{e_str}_tol/penbwss_k{max_k}_p{p_str}.csv"
        else:
            out_path = (
                f"./penalized_bwss/{graph_type}/vx{speed_ratio}/{e_str}_tol/k{max_k}/p{p_str}/"
                f"{graph_num}_penbwss_k{max_k}.csv"
            )

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        np.savetxt(out_path, bwss_penalized, delimiter=",")


if __name__ == "__main__":
    main()
