# Same build as on IBM Polus (http://hpc.cmc.msu.ru/node/241).
# After `module load SpectrumMPI` the XL wrapper is used verbatim:
#   mpixlC_r -O3 -qarch=pwr8 -qsmp=omp -std=c++11
# Without it, mpicxx or g++ gets the matching options:
#   -O3 -std=c++11 -fopenmp
# -fopenmp is the GNU form of -qsmp=omp. -qarch=pwr8 is POWER8-only.

XL := $(shell command -v mpixlC_r 2>/dev/null)
MPICXX := $(shell command -v mpicxx 2>/dev/null)
MPIRUN := $(shell command -v mpirun 2>/dev/null)

ifneq ($(XL),)
  CXX := mpixlC_r
  CXXFLAGS := -O3 -qarch=pwr8 -qsmp=omp -std=c++11
else ifneq ($(MPICXX),)
  CXX := mpicxx
  CXXFLAGS := -O3 -std=c++11 -fopenmp
else
  CXX := g++
  CXXFLAGS := -O3 -std=c++11 -fopenmp
endif

# make run mirrors: mpisubmit.pl -p NP -t THREADS ./wave -- N LX LY LZ T K
NP ?= 1
THREADS ?= 1
N ?= 128
LX ?= 1
LY ?= 1
LZ ?= 1
K ?= 20
# 90% of the CFL limit, so tau = T/K stays below 1/sqrt(1/hx^2+1/hy^2+1/hz^2).
T ?= $(shell awk -v n='$(N)' -v lx='$(LX)' -v ly='$(LY)' -v lz='$(LZ)' -v k='$(K)' 'BEGIN { printf "%.16g", 0.9 * k / sqrt((n/lx)^2 + (n/ly)^2 + (n/lz)^2) }')

wave: main.cpp solution.cpp solution.hpp
	$(CXX) $(CXXFLAGS) -o wave main.cpp solution.cpp

run: wave
ifneq ($(MPIRUN),)
	OMP_NUM_THREADS=$(THREADS) $(MPIRUN) -np $(NP) ./wave $(N) $(LX) $(LY) $(LZ) $(T) $(K)
else
	@if [ "$(NP)" != "1" ]; then echo "mpirun not found: running one process, NP=$(NP) ignored" >&2; fi
	OMP_NUM_THREADS=$(THREADS) ./wave $(N) $(LX) $(LY) $(LZ) $(T) $(K)
endif

clean:
	rm -f wave

.PHONY: clean run
