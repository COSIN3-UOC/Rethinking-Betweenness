# Rethinking Betweenness Centrality Beyond Shortest Paths

## Table of Contents
1. [Citing](#citing)
2. [Abstract](#abstract)
3. [Description](#description)
4. [Usage](#usage)
5. [Authors and acknowledgment](#authors-and-acknowledgment)
6. [Requirements](#requirements)
7. [License](#license)

## 1. Citing
If you use our approach or code in this repository, please cite our paper: <br>


## 2. Abstract

Centrality measures are fundamental tools for quantifying the importance of nodes and links in networks. Geodesic-based metrics, such as shortest-path betweenness, implicitly rely on specific assumptions about traversal dynamics, most notably that flows systematically follow shortest paths. While widely adopted, this assumption is often overly idealized and may fail to capture realistic dynamics, motivating a reassessment of how such centrality measures are defined and interpreted. Using transportation networks as a well-characterized case study, we show that even small relaxations of the shortest-path assumption can induce qualitative changes in centrality estimates. Empirical evidence indicates that real mobility frequently deviates from optimal routes by approximately 15–20\%, a feature not captured by standard metrics. To account for this, we introduce Quasi-Shortest Path Betweenness, which incorporates all paths lying within a prescribed tolerance of the shortest path. By comparing it with classical shortest-path betweenness, we show that increasing path tolerance systematically amplifies node betweenness and network load, generally degrading performance. When priority routes are present, however, a two-regime behaviour emerges in which performance either improves or deteriorates depending on transport demand. Overall, our results demonstrate that shortest-path betweenness corresponds to a restrictive, idealized limit of more general traversal dynamics, and that moving beyond this limit fundamentally alters centrality-based predictions, calling for a systematic reassessment of geodesic-based centrality measures.

## 3. Description

This repository contains the code needed to compute and store the quasi-shortest paths (QSPs), the Quasi-Shortest Path Betweenness (QSP-BW) and the packet dynamics together with all plots present in the manuscript. The repository is structured as follows:

- MCM contains the scripts and files needed to simulate the packet dynamics routed through QSPs reproducing the congestion phase transitions.

- KSP contains three subfolders
	- **graphs_qspbwss** which has all the graphs, paths and distance distributions used in the work.
	- **KSP_cpp** which contains the scripts to compute the QSPs and the QSP-BW using a prallelised K-shortest path listing algorithm
	- **KSP_cpp** which has most of the scripts dedicated to producing the plots in the manuscript as well as the path penalisation algorithm used to restrict overlap in the QSPs.

Information on the usage of the scripts in each folder is provided in the READMES within them.

## 4. Data and Extraction

The only dataset used in this project is OpenStreetMap (OSM) to access the road network of Barcelona. Details on the construction of the network can be found on the manuscript's Appendix and the script for the construction of the network is located at `./KSP/KSP_python/get_bcn_osm.ipynb`.

## 5. Usage

This repository allows the computation of:

1. The QSPs
2. The QSP-BW
3. The packet routing simulations

### 5.1. Computing the QSPs and the QSP-BW

If the objective of the user is to compute a fixed number of quasi-shortest paths, the user may use the script `./KSP/KSP_cpp/compute_KSPs.cpp`. Alternatively, if the user wants to compute the QSP-BW up to a given tolerance, then the file to run is `./KSP/KSP_cpp/compute_KSPs.cpp`. Please note that the original graph from which the paths need to be computed has to be stored in a `.csv` file containing columns `source_ID`, `target_ID`, `edge_weight`. Additionally, the source and output paths need to be edited in each script. More details on the compilation process can be found on the `README` in the `KSP_cpp`folder.

### 5.2. Computing the packet simulations

In the case the user wants to reproduce the packet simulations to get the congestion transition curves, the corresponding scripts are located in the `./MCM` and `./MCM/mcm` folders. The scripts that run the simulations are written in MATLAB and are located at `./MCM`. They they rely on `C` functions that need to be compiled beforehand and are located at `./MCM/mcm`.

More information on the compilation and usage of these functions can be found in the `README` files present in each folder.

## 6. Requirements

- Python version: 3.11.4

- MATLAB version: MATLAB_R2025b. Previous versions may work, but the minimum is the 2018 one.

### 6.1. C and C++

#### Requirements (macOS)

**Operating System**
- macOS (Apple Silicon assumed)
- Homebrew installed under `/opt/homebrew`

**Compiler**
- C++ compiler with **C++20** support  
  - Apple `clang++` (via Xcode Command Line Tools)

**Package Manager**
- Homebrew

**MATLAB**
- MATLAB **R2024a** or more installed at:
  - `/Applications/MATLAB_R2024a.app`
- MATLAB headers and libraries must exist at:
  - `/Applications/MATLAB_R2024a.app/extern/include`
  - `/Applications/MATLAB_R2024a.app/bin/maca64`
  - `/Applications/MATLAB_R2024a.app/extern/bin/maca64`

**Libraries**
Installed (typically via Homebrew) and available under `/opt/homebrew`:
- **Boost**
  - `boost_filesystem`
  - `boost_graph`
  - `boost_process`
- **OpenMP**
  - `libomp`

Install dependencies:
```bash
brew install boost libomp
```
### 6.2. Requirements (Linux)

**Operating System**
- 64-bit Linux distribution

**Compiler**
- `g++` (default) or `clang++`
- C++ compiler with **C++20** support (`-std=c++2a`)

**MATLAB**
- MATLAB installed on the system  
  - Default expected path:
    ```
    /usr/local/R2024a
    ```
  - Custom location can be provided at build time:
    ```bash
    make MATLABROOT=/path/to/MATLAB
    ```

Required MATLAB directories:
- `${MATLABROOT}/extern/include`
- `${MATLABROOT}/bin/glnxa64`
- `${MATLABROOT}/extern/bin/glnxa64`

**Libraries**
The following libraries must be available on the system:

- **Boost**
  - `boost_filesystem`
  - `boost_graph`
  - `boost_system`
- **POSIX Threads**
  - `pthread`

**OpenMP**
- **GCC**:
  - Uses `-fopenmp` and links `libgomp` automatically
- **Clang**:
  - Requires `libomp`

**Include and Library Paths**
Defaults used by the Makefile:
- Include paths:
	`/usr/include`	
- Library paths:
	`/usr/lib`

Custom paths can be provided at build time:
```bash
make EXTRA_INC=/custom/include EXTRA_LIB=/custom/lib
```

- Python libraries:
  ```
	contextily==1.4.0
	geopandas==0.14.0
	h5py==3.7.0
	matplotlib==3.7.1
	networkx==3.2.1
	numpy==1.24.3
	osmnx==2.0.7
	pandas==1.5.3
	pydot==1.4.2
	python_igraph==1.0.0
	scikit_learn==1.3.0
	scipy==1.10.1
	seaborn==0.13.2
	Shapely==2.1.2
	tqdm==4.65.0
	treelib==1.7.0
  ```


## 7. License
This project is licensed under the [MIT License](https://mit-license.org/).


