#!/bin/bash
# Comprehensive benchmark script for optimized PEM simulator

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}PEM Simulator Performance Benchmark${NC}"
echo -e "${GREEN}========================================${NC}"

# Configuration
WORK_DIR=$(pwd)
SRC_DIR="${WORK_DIR}/src"
BIN_DIR="${WORK_DIR}/bin"
DATA_DIR="${WORK_DIR}/data"
RESULTS_DIR="${WORK_DIR}/benchmark_results"
INPUT_FILE="${WORK_DIR}/input/benchmark.dat"

# Create directories
mkdir -p "${BIN_DIR}" "${RESULTS_DIR}" "${DATA_DIR}"

# System information
echo -e "\n${YELLOW}System Information:${NC}"
echo "Hostname: $(hostname)"
echo "CPU: $(lscpu | grep 'Model name' | sed 's/Model name://g' | xargs)"
echo "Cores: $(nproc)"
echo "Memory: $(free -h | grep Mem | awk '{print $2}')"
echo "Compiler: $(ifort --version 2>/dev/null | head -1 || echo "Intel Fortran not found")"

# Compilation benchmark
echo -e "\n${YELLOW}Compilation Benchmark:${NC}"

# Compile original version
echo -n "Compiling original version... "
COMPILE_START=$(date +%s.%N)
ifort -O3 -o "${BIN_DIR}/pem_original" "${SRC_DIR}/pem_simulator.f90" 2>&1
COMPILE_END=$(date +%s.%N)
COMPILE_TIME_ORIG=$(echo "$COMPILE_END - $COMPILE_START" | bc)
echo -e "${GREEN}Done (${COMPILE_TIME_ORIG}s)${NC}"

# Compile optimized version
echo -n "Compiling optimized version... "
COMPILE_START=$(date +%s.%N)
ifort -O3 -xHost -qopenmp -ipo -fp-model fast=2 -align array64byte \
      -o "${BIN_DIR}/pem_optimized" "${SRC_DIR}/pem_simulator.f90" 2>&1
COMPILE_END=$(date +%s.%N)
COMPILE_TIME_OPT=$(echo "$COMPILE_END - $COMPILE_START" | bc)
echo -e "${GREEN}Done (${COMPILE_TIME_OPT}s)${NC}"

# Runtime benchmark with different configurations
echo -e "\n${YELLOW}Runtime Benchmark:${NC}"

# Test configurations
declare -a PARTICLES=("50" "100" "200" "500")
declare -a THREADS=("1" "2" "4" "8")
declare -a STEPS=("1000" "5000" "10000")

# Results file
RESULTS_CSV="${RESULTS_DIR}/benchmark_$(date +%Y%m%d_%H%M%S).csv"
echo "particles,steps,threads,version,time,speedup" > "$RESULTS_CSV"

# Function to run benchmark
run_benchmark() {
    local particles=$1
    local steps=$2
    local threads=$3
    local version=$4
    local executable=$5
    
    # Create temporary input file
    local temp_input="${RESULTS_DIR}/temp_input_${particles}_${steps}.dat"
    cat > "$temp_input" << EOF
TIME_STEP 5.0e-7
PARTICLE_GEN_LAYERS ${particles}
MAX_CALCULATION_STEPS ${steps}
OUTPUT_INTERVAL_NORMAL 1000
NUM_THREADS ${threads}
EOF
    
    # Set thread count
    export OMP_NUM_THREADS=${threads}
    
    # Run simulation
    local start_time=$(date +%s.%N)
    timeout 300 "${executable}" "$temp_input" > "${RESULTS_DIR}/output_${version}_${particles}_${steps}_${threads}.log" 2>&1
    local end_time=$(date +%s.%N)
    local runtime=$(echo "$end_time - $start_time" | bc)
    
    # Clean up
    rm -f "$temp_input"
    
    echo "$runtime"
}

# Run benchmarks
echo -e "\nRunning benchmarks..."
echo -e "Particles\tSteps\tThreads\tOriginal(s)\tOptimized(s)\tSpeedup"

for particles in "${PARTICLES[@]}"; do
    for steps in "${STEPS[@]}"; do
        # Original version (single-threaded)
        echo -n -e "${particles}\t${steps}\t1\t"
        time_orig=$(run_benchmark "$particles" "$steps" "1" "original" "${BIN_DIR}/pem_original")
        echo -n -e "${time_orig}\t"
        
        # Optimized versions with different thread counts
        best_time=$time_orig
        best_threads=1
        
        for threads in "${THREADS[@]}"; do
            if [ "$threads" == "1" ]; then
                time_opt=$(run_benchmark "$particles" "$steps" "$threads" "optimized" "${BIN_DIR}/pem_optimized")
                speedup=$(echo "scale=2; $time_orig / $time_opt" | bc)
                echo -e "${time_opt}\t${speedup}x"
                echo "$particles,$steps,$threads,optimized,$time_opt,$speedup" >> "$RESULTS_CSV"
            else
                echo -n -e "${particles}\t${steps}\t${threads}\t-\t"
                time_opt=$(run_benchmark "$particles" "$steps" "$threads" "optimized" "${BIN_DIR}/pem_optimized")
                speedup=$(echo "scale=2; $time_orig / $time_opt" | bc)
                echo -e "${time_opt}\t${speedup}x"
                echo "$particles,$steps,$threads,optimized,$time_opt,$speedup" >> "$RESULTS_CSV"
            fi
            
            if (( $(echo "$time_opt < $best_time" | bc -l) )); then
                best_time=$time_opt
                best_threads=$threads
            fi
        done
        
        echo "$particles,$steps,1,original,$time_orig,1.0" >> "$RESULTS_CSV"
        echo -e "${GREEN}Best configuration: ${best_threads} threads (${best_time}s)${NC}"
        echo ""
    done
done

# Memory usage benchmark
echo -e "\n${YELLOW}Memory Usage Analysis:${NC}"
for particles in "${PARTICLES[@]}"; do
    echo -n "Testing with $particles particles: "
    
    # Create test input
    cat > "${RESULTS_DIR}/mem_test.dat" << EOF
PARTICLE_GEN_LAYERS ${particles}
MAX_CALCULATION_STEPS 100
OUTPUT_INTERVAL_NORMAL 1000
EOF
    
    # Run with memory profiling
    /usr/bin/time -v "${BIN_DIR}/pem_optimized" "${RESULTS_DIR}/mem_test.dat" \
        2>&1 | grep "Maximum resident set size" | awk '{print $6 " KB"}'
done

# Vectorization analysis
echo -e "\n${YELLOW}Vectorization Report:${NC}"
ifort -O3 -xHost -qopt-report=5 -qopt-report-phase=vec,loop \
      -c "${SRC_DIR}/pem_simulator.f90" -o "${RESULTS_DIR}/pem_vec.o" 2>&1 | \
      grep -E "LOOP WAS VECTORIZED|remark" | head -20

# Performance summary
echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}Performance Summary${NC}"
echo -e "${GREEN}========================================${NC}"

# Calculate average speedups
echo -e "\n${YELLOW}Average Speedup by Thread Count:${NC}"
for threads in "${THREADS[@]}"; do
    avg_speedup=$(awk -F',' -v t="$threads" '$3==t && $4=="optimized" {sum+=$6; count++} END {if(count>0) printf "%.2f", sum/count}' "$RESULTS_CSV")
    if [ -n "$avg_speedup" ]; then
        echo "  $threads threads: ${avg_speedup}x"
    fi
done

# Plot results if Python is available
if command -v python3 &> /dev/null; then
    echo -e "\n${YELLOW}Generating performance plots...${NC}"
    python3 << EOF
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read results
df = pd.read_csv('$RESULTS_CSV')

# Create performance plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Speedup vs threads
for particles in df['particles'].unique():
    data = df[(df['particles'] == particles) & (df['version'] == 'optimized')]
    if not data.empty:
        avg_speedup = data.groupby('threads')['speedup'].mean()
        ax1.plot(avg_speedup.index, avg_speedup.values, 'o-', label=f'{particles} particles')

ax1.set_xlabel('Number of Threads')
ax1.set_ylabel('Speedup')
ax1.set_title('Speedup vs Thread Count')
ax1.legend()
ax1.grid(True)

# Scaling efficiency
for particles in df['particles'].unique():
    data = df[(df['particles'] == particles) & (df['version'] == 'optimized')]
    if not data.empty:
        avg_speedup = data.groupby('threads')['speedup'].mean()
        efficiency = avg_speedup / avg_speedup.index
        ax2.plot(avg_speedup.index, efficiency.values, 'o-', label=f'{particles} particles')

ax2.set_xlabel('Number of Threads')
ax2.set_ylabel('Parallel Efficiency')
ax2.set_title('Parallel Efficiency')
ax2.legend()
ax2.grid(True)
ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
ax2.axhline(y=0.8, color='r', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('${RESULTS_DIR}/performance_analysis.png', dpi=150)
print("Performance plots saved to ${RESULTS_DIR}/performance_analysis.png")
EOF
fi

echo -e "\n${GREEN}Benchmark complete!${NC}"
echo "Results saved to: $RESULTS_CSV"
echo "Logs saved to: ${RESULTS_DIR}/"

# Final recommendations
echo -e "\n${YELLOW}Recommendations:${NC}"
echo "1. Use 4-8 threads for optimal performance on this system"
echo "2. Enable cell algorithm for problems with >100 particles"
echo "3. Consider using profile-guided optimization for production runs"
echo "4. Monitor memory usage for large-scale simulations"