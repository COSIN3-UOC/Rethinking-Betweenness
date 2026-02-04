import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap

# -----------------------------
# PLOT THE CONGESTION IN DIFFERENT PARTS OF THE GRAPH
# -----------------------------

# --- find ring folders and extract numeric IDs robustly
def cmap_with_alpha_ramp(cmap, alpha_min=0.0, alpha_max=1.0, ramp_fraction=0.3):
    """
    ramp_fraction: fraction of the colormap (from bottom) over which alpha ramps up
    """
    colors = cmap(np.linspace(0, 1, cmap.N))
    alphas = np.ones(cmap.N)

    n_ramp = int(cmap.N * ramp_fraction)
    print('n_ramp:', n_ramp)  
    alphas[:n_ramp] = np.linspace(alpha_min, alpha_max, n_ramp)

    colors[:, -1] = alphas
    return mcolors.ListedColormap(colors)

import numpy as np
import colorsys
from matplotlib.colors import ListedColormap

def desaturate_cmap(cmap, factor=0.6):
    """
    factor < 1  → less saturated
    factor = 1  → unchanged
    """
    colors = cmap(np.linspace(0, 1, cmap.N))
    new_colors = []

    for r, g, b, a in colors:
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        s *= factor
        new_colors.append(colorsys.hsv_to_rgb(h, s, v) + (a,))

    return ListedColormap(new_colors)


base_dir = "../../MCM/mcm/ring_road_p_structure_vx2/"
folders = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
ring_folders = [d for d in folders if "ring_road_" in d]

nnodes = 58
graph_nums = []
for name in ring_folders:
    m = re.findall(r"\d+", name)
    if m:
        graph_nums.append(int(m[0]))

graph_nums = np.array(graph_nums, dtype=int)

# --- you need possGenRates defined somewhere (MATLAB had it in workspace)
possGenRates = np.logspace(-2, np.log10(3), 100)
# -----------------------------
# buffers
# -----------------------------
len_fold = len(graph_nums)

accum_x = np.zeros((len_fold, nnodes), dtype=float)
accum_y = np.zeros((len_fold, nnodes), dtype=float)
accum_congestion0  = np.zeros((len_fold, nnodes), dtype=float)
accum_congestion10 = np.zeros((len_fold, nnodes), dtype=float)
accum_congestion20 = np.zeros((len_fold, nnodes), dtype=float)
accum_congestion30 = np.zeros((len_fold, nnodes), dtype=float)

genrate_ind = 40

# -----------------------------
# loop over graphs
# -----------------------------
for gi, g in enumerate(graph_nums):
    file1 = f"../../MCM/sim_results/all_graphs/ring/vx2/tol0_{g}_vx2.mat"
    file2 = f"../../MCM/sim_results/all_graphs/ring/vx2/tol10_{g}_vx2.mat"
    file3 = f"../../MCM/sim_results/all_graphs/ring/vx2/tol20_{g}_vx2.mat"
    file4 = f"../../MCM/sim_results/all_graphs/ring/vx2/tol30_{g}_vx2.mat"

    need = [file1, file2, file3]  # (MATLAB had file4 optional but loads it later)
    if any(not os.path.isfile(f) for f in need):
        print(f"Warning: Missing file(s) for graph {g}. Skipping.")
        continue

    S0  = loadmat(file1)
    S10 = loadmat(file2)
    S20 = loadmat(file3)
    S30 = loadmat(file4) if os.path.isfile(file4) else None

    # expects DeltaNExp in each file
    DeltaNExp0  = np.asarray(S0["DeltaNExp"])
    DeltaNExp10 = np.asarray(S10["DeltaNExp"])
    DeltaNExp20 = np.asarray(S20["DeltaNExp"])
    DeltaNExp30 = np.asarray(S30["DeltaNExp"]) if S30 is not None and "DeltaNExp" in S30 else None

    # MATLAB: y0 = S0.DeltaNExp./ possGenRates(genrate_ind);
    # NOTE: possGenRates(genrate_ind) is scalar, divides the whole matrix
    denom = possGenRates[genrate_ind]
    y0  = DeltaNExp0 / denom
    y10 = DeltaNExp10 / denom
    y20 = DeltaNExp20 / denom
    y30 = (DeltaNExp30 / denom) if DeltaNExp30 is not None else None

    # positions
    pos_path = f"../graphs_qspbwss/ring_roads/vx2/node_pos/pos_ring{g}.csv"
    T = pd.read_csv(pos_path)

    # MATLAB columns: vertex1, x, y
    x = T["x"].to_numpy(dtype=float)
    y = T["y"].to_numpy(dtype=float)

    # store
    accum_x[gi, :] = x
    accum_y[gi, :] = y

    # MATLAB: accum_congestion10(gi,:) = y10(:,genrate_ind);
    # so y10 is (nodes × rho). Take column genrate_ind.
    accum_congestion0[gi, :]  = y0[:, genrate_ind]
    accum_congestion10[gi, :] = y10[:, genrate_ind]
    accum_congestion20[gi, :] = y20[:, genrate_ind]
    if y30 is not None:
        accum_congestion30[gi, :] = y30[:, genrate_ind]

# -----------------------------
# remove rows that are all zero in congestion20
# -----------------------------
idx_all = ~np.all(accum_congestion20 == 0, axis=1)

accum_bw_noz = accum_congestion20[idx_all, :]
x_nnz = accum_x[idx_all, :]
y_nnz = accum_y[idx_all, :]

# flatten to vectors
z = accum_bw_noz.reshape(-1)
xv = x_nnz.reshape(-1)
yv = y_nnz.reshape(-1)


# Base magma colormap
base_cmap = plt.colormaps["inferno"]

# Sample only the first 80% of the colormap
colors = base_cmap(np.linspace(0.0, 0.95, 256))

# Create a new ListedColormap
cmap = ListedColormap(colors)

# Apply your alpha ramp
cmap = cmap_with_alpha_ramp(
    cmap,
    alpha_min=0.1,
    alpha_max=1.0,
    ramp_fraction=0.5
)

#cmap = desaturate_cmap(cmap, factor=0.5)

# ---- filter finite values ----
m = np.isfinite(xv) & np.isfinite(yv) & np.isfinite(z)
x = xv[m]
y = yv[m]
zz = z[m]

fig, ax = plt.subplots(figsize=(8, 6))

# ---- hexbin with mean aggregation ----
hb = ax.hexbin(
    x, y,
    C=zz,
    reduce_C_function=np.mean,   # mean(z) per hex
    gridsize=40,                 # roughly comparable to 40x40 bins
    mincnt=1,                    # ignore empty hexes
    cmap=cmap,
    linewidths=0,                # avoid seam artifacts
)

# match your color scaling
hb.set_clim(0, 1)

ax.set_aspect("equal")
ax.tick_params(labelsize=16)

# # ---- highlight the max hex (red outline) ----
arr = hb.get_array()  # aggregated values per hex (length = #hexes)
# if arr.size:
#     i_max = int(np.argmax(arr))
#     print(i_max, arr.shape, np.shape(hb.get_paths()))
#     verts = hb.get_paths()[i_max].vertices  # hex polygon in data coords
#     ax.add_patch(Polygon(verts, fill=False, edgecolor="red", linewidth=2))

# ---- colorbar ----
cbar = plt.colorbar(hb, ax=ax)
cbar_lab = rf"$\eta (\varepsilon = 20\%$, $\rho = {possGenRates[genrate_ind]:.2f}$)"
cbar.set_label(cbar_lab, fontsize=22)

print("max mean per hex:", float(arr.max()) if arr.size else np.nan)

plt.tight_layout()
plt.show()

