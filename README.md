## Build QSP-BWSS structure for MCM

Modify the path of the csv with the graph (edgelist with weights) in `generate_p_structure_cpp.cpp`. 
Modify the output path for the successor matrices (1 successor matrix for each node) and for the probability matrices.

```
# macports include files, for brew, typically in /usr/local
INC_DIRS = -I/opt/homebrew/include -I/Applications/MATLAB_R2023b.app/extern/include
LIB_DIRS = -L/opt/homebrew/lib -L/Applications/MATLAB_R2023b.app/bin/maca64 -L/Applications/MATLAB_R2023b.app/extern/bin/maca64

CPPFLAGS += -std=c++20
CFLAGS += -O3 -ffast-math -Wall -pedantic

# Compiler options for ThreadSanitizer
TSAN_FLAGS = -g


#-fsanitize=thread 

# Link boost and other libraries
LDFLAGS = -lboost_filesystem -lboost_graph -Xpreprocessor -fopenmp -lomp -lmat -lmx -leng -lmex

#Yen_VF: VF_KSP.cpp
Yen_tolSPR: generate_p_structure_cpp.cpp
# amazon_SPR: SPR_tol_amazon.cpp
# Yen_SPR: SPR_KSP.cpp
#Test: test_lattice.cpp
#Yen: yen_ksp.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(TSAN_FLAGS) $(INC_DIRS) $(LIB_DIRS) $(LDFLAGS) $< -o $@
	install_name_tool -add_rpath /System/Volumes/Data/Applications/MATLAB_R2023b.app/bin/maca64 $@


clean:
	rm Yen_SPR
```

