# Makefile for PEM Simulator with Performance Optimizations

# Compiler selection
FC = ifort
GFORTRAN = gfortran

# Base directories
SRCDIR = src
BINDIR = bin
OBJDIR = obj

# Source files
ORIG_SRC = $(SRCDIR)/pem_simulator.f90
OPT_SRC = $(SRCDIR)/optimized_pem_simulator.f90

# Executable names
ORIG_EXE = $(BINDIR)/pem_simulator
OPT_EXE = $(BINDIR)/pem_optimized
DEBUG_EXE = $(BINDIR)/pem_debug

# Common flags
COMMON_FLAGS = -module $(OBJDIR)

# Optimization flags for Intel Fortran
IFORT_OPT_FLAGS = -O3 -xHost -ipo -qopenmp -fp-model fast=2 -align array64byte
IFORT_AGGRESSIVE_FLAGS = -Ofast -march=native -qopenmp -parallel -unroll-aggressive
IFORT_DEBUG_FLAGS = -O0 -g -check all -traceback -qopenmp
IFORT_PROFILE_FLAGS = -O3 -g -pg -qopenmp

# Optimization flags for GNU Fortran
GFORTRAN_OPT_FLAGS = -O3 -march=native -fopenmp -funroll-loops -ffast-math
GFORTRAN_AGGRESSIVE_FLAGS = -Ofast -march=native -fopenmp -funroll-all-loops
GFORTRAN_DEBUG_FLAGS = -O0 -g -fcheck=all -fbacktrace -fopenmp

# Default target
all: dirs optimize

# Create necessary directories
dirs:
	@mkdir -p $(BINDIR) $(OBJDIR)

# Original version with basic optimization
original: dirs
	$(FC) $(COMMON_FLAGS) -O3 -o $(ORIG_EXE) $(ORIG_SRC)
	@echo "Built original version with -O3 optimization"

# Optimized version with Intel compiler
optimize: dirs
	$(FC) $(COMMON_FLAGS) $(IFORT_OPT_FLAGS) -o $(OPT_EXE) $(OPT_SRC)
	@echo "Built optimized version with Intel compiler"

# Maximum optimization for Intel CPUs
optimize-intel: dirs
	$(FC) $(COMMON_FLAGS) $(IFORT_AGGRESSIVE_FLAGS) -o $(OPT_EXE)_intel $(OPT_SRC)
	@echo "Built with aggressive Intel optimizations"

# Optimized version with GNU compiler
optimize-gnu: dirs
	$(GFORTRAN) $(GFORTRAN_OPT_FLAGS) -o $(OPT_EXE)_gnu $(OPT_SRC)
	@echo "Built optimized version with GNU compiler"

# Debug build with optimization info
debug-opt: dirs
	$(FC) $(COMMON_FLAGS) $(IFORT_DEBUG_FLAGS) -qopt-report=5 -o $(DEBUG_EXE) $(OPT_SRC)
	@echo "Built debug version with optimization reports"

# Profile-guided optimization build
profile: dirs
	@echo "Building for profiling..."
	$(FC) $(COMMON_FLAGS) $(IFORT_PROFILE_FLAGS) -o $(OPT_EXE)_prof $(OPT_SRC)
	@echo "Run the program with typical workload, then run 'make profile-use'"

# Use profile data for optimization
profile-use: dirs
	$(FC) $(COMMON_FLAGS) $(IFORT_OPT_FLAGS) -prof-use -o $(OPT_EXE)_pgo $(OPT_SRC)
	@echo "Built with profile-guided optimization"

# Vectorization report
vec-report: dirs
	$(FC) $(COMMON_FLAGS) $(IFORT_OPT_FLAGS) -qopt-report=5 -qopt-report-phase=vec,loop -c $(OPT_SRC)
	@echo "Vectorization report generated in $(OBJDIR)/"

# Benchmarking targets
benchmark: optimize
	@echo "Running benchmark with optimized version..."
	cd scripts && ./benchmark_optimized.sh

benchmark-compare: original optimize
	@echo "Comparing original vs optimized performance..."
	cd scripts && ./compare_performance.sh

# Testing targets
test: optimize
	@echo "Running validation tests..."
	$(OPT_EXE) input/test_validation.dat

test-all: original optimize optimize-gnu
	@echo "Testing all versions..."
	@$(ORIG_EXE) input/test_validation.dat > test_orig.out
	@$(OPT_EXE) input/test_validation.dat > test_opt.out
	@$(OPT_EXE)_gnu input/test_validation.dat > test_gnu.out
	@echo "Comparing outputs..."
	@diff test_orig.out test_opt.out || echo "Intel optimized differs from original"
	@diff test_orig.out test_gnu.out || echo "GNU optimized differs from original"

# Performance analysis
perf-analysis: optimize
	@echo "Running performance analysis..."
	perf record -g $(OPT_EXE) input/benchmark.dat
	perf report

vtune-analysis: optimize
	@echo "Running Intel VTune analysis..."
	vtune -collect hotspots -r vtune_results $(OPT_EXE) input/benchmark.dat
	vtune -report summary -r vtune_results

# Clean targets
clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod

distclean: clean
	rm -rf $(BINDIR)/* vtune_results perf.data* *.out

# Installation
install: optimize
	@echo "Installing optimized version..."
	cp $(OPT_EXE) /usr/local/bin/pem_simulator
	@echo "Installation complete"

# Help
help:
	@echo "PEM Simulator Makefile targets:"
	@echo "  make              - Build optimized version (default)"
	@echo "  make original     - Build original version with -O3"
	@echo "  make optimize     - Build optimized version with Intel compiler"
	@echo "  make optimize-gnu - Build optimized version with GNU compiler"
	@echo "  make debug-opt    - Build debug version with optimization info"
	@echo "  make profile      - Build for profile-guided optimization"
	@echo "  make vec-report   - Generate vectorization report"
	@echo "  make benchmark    - Run performance benchmark"
	@echo "  make test         - Run validation tests"
	@echo "  make clean        - Remove object files"
	@echo "  make help         - Show this help message"

.PHONY: all dirs original optimize optimize-intel optimize-gnu debug-opt profile profile-use vec-report benchmark benchmark-compare test test-all perf-analysis vtune-analysis clean distclean install help