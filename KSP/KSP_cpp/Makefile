
# macports include files, for brew, typically in /usr/local
INC_DIRS = /opt/homebrew/include
LIB_DIRS = /opt/homebrew/lib

CPPFLAGS += -std=c++20
CFLAGS += -O3 -Wall -pedantic

# Compiler options for ThreadSanitizer
TSAN_FLAGS = -g


#-fsanitize=thread 

# Link boost and other libraries
LDFLAGS = -lboost_filesystem -lboost_graph -Xpreprocessor -fopenmp -lomp -lboost_process

Yen_tolSPR: compute_KSPs.cpp
#amazon_SPR: SPR_tol_amazon_digraph.cpp
#Yen_SPR: compute_qsp_bw.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(TSAN_FLAGS) -I$(INC_DIRS) -L$(LIB_DIRS) $(LDFLAGS) $< -o $@

clean:
	rm Yen_tolSPR
