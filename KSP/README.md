As stated in the `README` of the parent folder, this directory contains the graphs and path distributions in `graphs_qspbwss`, the python files generating most plots and the algorithm of the penalised path calculation in `KSP_python` and `KSP_python/other_sampling_methods/` respectively. It also has the scripts to compute the QSPs and the QSP-BW in `KSP_cpp` written in `C++` for efficiency.

In the following, a detailed description of the functionalities of each script is provided.

# KSP_python scripts

## Computing and plotting QSP-BW in the continuum:
1. `compute_amnalyt_bwss_cont_homog.py`: Computes the QSP betweenness in a continuous 2D disk representing an infinitely homogeneous planar graph. 
Inputs: 
- $\gamma$: The exponent of the relationship $\sigma^\varepsilon = \sigma e^{\gamma (\delta-\delta^*)}$
- $\varepsilon$: Tolerance to the shortest path ($\varepsilon\in[0,1]$)
Output:
Files containing the Betweenness and Error in a polar mesh in the files:

`'./bwss_continuum/bwss_polar_{int(args.eps*100)}_gam{int(args.gamma)}.txt'`
`'./bwss_continuum/err_bwss_polar_{int(args.eps*100)}_gam{int(args.gamma)}.txt'`

2. `plot_cont_bwss.ipynb`: Plots the results from 1 giving the manuscript's Figure 5 (a) and (b).

## Building the networks

1. Getting Barcelona's road network: The road network from Barcelona and surroundings is obtained usinf `get_bcn_osm.ipynb`.
2. Toy models: `generate_toy_graphs.ipynb``
3. Delaunay Graphs: `create_delaunays.py`

## Other plots

1. `plot_congestion_hexbin.py` is used to get the insets in Fig. 8 showing the state of congestion for different generation rates.
2. `plot_degeneracy_empiric.ipynb` is used for the plots in Fig. 2.
3. `plot_hexbins_bw.ipynb` is used for Fig. 7.

## Path penalisation algorithm

The path penalisation algorithm computes at most k paths given a s,t pair. At each iteration in computes the shortest path and then multiplies the path's weights by (1+$p$). The algorithm keeps the path if its distance is within (1+$\varepsilon$) of the SP. Then it recomputes the shortest path and so on until it finds k distinct paths or reaches maximum iterations.

The scripts `routing_utils_optim.py` and `routing_utils_parallel.py` both compute the QSP betweenness with the paths obtained with the path penalisation algorithm. The script runs the algorithm in for each s,t pair and stores the path lists in a file to then run the packet congestion simulations. You may need to modify the output paths if using `'other'`as graph type.

input params:
`graph_type`:(str) `'ring'`, `'star'`, `'lattice'`, `'other'`
`speed_ratio`:(str) important for toy models, options = (`'1'`, `'2'`)
`k`:(int) Maximum number of paths
`p`:(float) penalisation (1+`p`)
$\varepsilon$:(float) tolerance

Example Call
`python routing_utils_optim.py 'ring' '2' 20 0.1 0.3`
`python routing_utils_parallel.py 'ring' '2' 20 0.1 0.3`

# KSP_cpp scripts

- `compute_KSPs.cpp`gives K first shortest paths of randomly selected node pairs. The number of paths and node pairs can be modified from within the file, as well as the input and output paths corresponding to the graph and path list respectively. This script was used to compute the path distributions shown in Fig. 2.

- `compute_qsp_bw.cpp` Computes the QSP-BW given a tolerance set within the file (inside main). The path to the input graph and the output path may need to be modified.

The user may compile the scripts using the `Makefile`, although it is set for MacOS. for Linux users, you may use a similar `Makefile` as `Makefile_linux` in the `MCM`folder.

- `qsp_bw.h` contains the functions that compute the paths using a modification of Yen's Algorithm [1].

[1] Chen, B. Y., Chen, X. W., Chen, H. P., & Lam, W. H. (2021). A fast algorithm for finding K shortest paths using generalized spur path reuse technique. Transactions in GIS, 25(1), 516-533.

