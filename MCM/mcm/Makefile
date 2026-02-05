# # macports include files, for brew, typically in /usr/local
# INC_DIRS = -I/opt/homebrew/include -I/Applications/MATLAB_R2024a.app/extern/include
# LIB_DIRS = -L/opt/homebrew/lib -L//Applications/MATLAB_R2024a.app/bin/maca64 -L/Applications/MATLAB_R2024a.app/extern/bin/maca64
# CPPFLAGS += -std=c++20
# CFLAGS += -O3 -ffast-math -Wall -pedantic

# # Compiler options for ThreadSanitizer
# TSAN_FLAGS = -g


# #-fsanitize=thread 

# # Link boost and other libraries
# LDFLAGS = -lboost_filesystem -lboost_graph -Xpreprocessor -fopenmp -lomp -lmat -lmx -leng -lmex -lboost_process


# #Yen_VF: VF_KSP.cpp
# Yen_tolSPR: generate_p_structure_cpp.cpp
# # amazon_SPR: SPR_tol_amazon.cpp
# # Yen_SPR: SPR_KSP.cpp
# #Test: test_lattice.cpp
# #Yen: yen_ksp.cpp
# 	$(CXX) $(CFLAGS) $(CPPFLAGS) $(TSAN_FLAGS) $(INC_DIRS) $(LIB_DIRS) $(LDFLAGS) $< -o $@
# 	install_name_tool -add_rpath /Applications/MATLAB_R2024a.app/bin/maca64 $@


# clean:
# 	rm Yen_tolSPR

# macports include files, for brew, typically in /usr/local
INC_DIRS = -I/opt/homebrew/include -I/Applications/MATLAB_R2024a.app/extern/include
LIB_DIRS = -L/opt/homebrew/lib -L/Applications/MATLAB_R2024a.app/bin/maca64 -L/Applications/MATLAB_R2024a.app/extern/bin/maca64

CPPFLAGS += -std=c++20
CFLAGS += -O3 -ffast-math -Wall -pedantic

# Compiler options for ThreadSanitizer
TSAN_FLAGS = -g
# -fsanitize=thread  # Uncomment if needed

# Link boost and other libraries + set RPATH at link time
LDFLAGS = -lboost_filesystem -lboost_graph -Xpreprocessor -fopenmp -lomp -lmat -lmx -leng -lmex -lboost_process -Wl,-rpath,/Applications/MATLAB_R2024a.app/bin/maca64

Yen_tolSPR: generate_p_list.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(TSAN_FLAGS) $(INC_DIRS) $(LIB_DIRS) $(LDFLAGS) $< -o $@

clean:
	rm -f Yen_tolSPR
