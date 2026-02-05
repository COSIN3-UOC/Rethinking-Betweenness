#  MCM: Microscopic congestion model

The code here reproduces partially the experiments published for the Microscopic Congestion Model developed by Solé-Ribalta et al, Royal Society Open Science 2016 (https://doi.org/10.1098/rsos.160098)

## System Requirements

- Matlab R2018a or above
- C99 and C++ compiler
- Boost C++ library (https://www.boost.org/)

## Structure and setup of the modified MCM packet

This is a modification of the code available from the [Microscopic congestion model](https://github.com/COSIN3-UOC/MCM).

The structure is as follows, a main file that runs the experiments `congestion_dyn.m` and a directory `mcm` with auxiliary files. The parallelised version is `congestion_dyn_parallel.m` and the version for computing the dynamics where paths are stored in `.h5` files (for large graphs like the city of Barcelona with thousands of nodes) is `congestion_dyn_h5.m`.

To run, the scripts read the stored paths in `.mat` files or `.h5` for large graphs present in `mcm`. Examples of folders where paths are stored are all ending in `path_structure`. The files storing the paths between all node pairs can be generated compiling and running `generate_p_list.cpp`. Two makefiles exist, on for MacOS and another for Linux. You may need to modify the include paths to suit your computer.

Computationally costly functions are implemented in C and C++ and need to be compiled before running the toy example.

There only one file needs to be compiled depending on the format of the file storing the paths. There are two possibilities:

1. Padded case: Paths are stored as:\
     `s, ..., t, 0, ..., 0, d`\
     `s, ..., t, 0, ..., 0, d`\
     `s, ..., t, 0, ..., 0, d`\
     ...\
    As an MxN padded array. In that case, compile `cSPCongestion_statMem_dir_weighted_local_list_paths.c`.

2. Ragged case: paths are stored in a ragged array with uneven shape\
    `s, ..., t, d`\
    `s, ......, t, d`\
    `s, .........., t, d`\
    ...\
    In this case, compile `cSPCongestion_statMem_dir_weighted_local_list_paths_no_pad.c`.

Compilation line example in the Matlab terminal within the Matlab interface:

`mex -R2018a cSPCongestion_statMem_dir_weighted_local_list_paths.c`

One needs to make sure that Matlab is correctly configured to compile C and C++ files, usually one does that using command `mex -setup c++`, see the documentation [here](https://es.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html) if required.

### Configuration of the C file

Definition `#define MAX_COLA 2000` within the C file can be set to account for maximum node queue length.

We remind the user to recompile the file once modified.


## Usage

1. If you are planning to compute the simulations in a graph that was not in the manuscript, first compute all the paths by running `generate_p_list.cpp`. This file .mat file containing an NxN CellArray. Each cell (i,j) has a list of the computed paths between a the node pair (i,j), and the last element in each row is the distance of the path. To run the file, first modify the input and output paths in the file `generate_p_list.cpp`. Then run `Make` and run the executable.

2. Once the paths are stored, compile the `cSPCongestion_statMem_dir_weighted_local_list_paths_no_pad.c` as detailed above. 

3. In any of the `congestion_dyn` files, modify the input paths of the CellArray outputed in step 1 and the output paths containing the transition curves. Then run the script.

You may plot the results using `plot_phasetr.m`.

# Citations

Solé-Ribalta, A., Gómez, S., & Arenas, A. (2016). A model to identify urban traffic congestion hotspots in complex networks. [Royal Society open science, 3(10), 160098](https://doi.org/10.1007/s11067-017-9349-y).

Solé-Ribalta, A., Gómez, S., & Arenas, A. (2018). Decongestion of urban areas with hotspot pricing. [Networks and Spatial Economics, 18, 33-50](https://link.springer.com/article/10.1007/s11067-017-9349-y).
