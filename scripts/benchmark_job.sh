#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr20001a
#SBATCH -t 12:00:00
#SBATCH --rsc p=1:t=4:c=4:m=16G
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -J pem_benchmark

#============ Environment Setup ============
module load intel

#============ Shell Script ============
echo "================================="
echo "PEM Algorithm Benchmark Job"
echo "Start Time: $(date)"
echo "================================="

# Compile the Fortran code
echo "Compiling Fortran code..."
cd src
ifort -O3 -qopenmp -o pem_simulator pem_simulator.f90
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi
echo "Compilation successful!"

# Move executable to parent directory
mv pem_simulator ..
cd ..

# Run benchmark with cell algorithm ON
echo "================================="
echo "Running benchmark with CELL ALGORITHM ON"
echo "================================="
echo "Start Time: $(date)"
./pem_simulator input_benchmark.dat > benchmark_on.log 2>&1
echo "End Time: $(date)"
echo "Cell algorithm ON completed!"

# Run benchmark with cell algorithm OFF
echo "================================="
echo "Running benchmark with CELL ALGORITHM OFF"
echo "================================="
echo "Start Time: $(date)"
./pem_simulator input_benchmark_off.dat > benchmark_off.log 2>&1
echo "End Time: $(date)"
echo "Cell algorithm OFF completed!"

# Extract and compare results
echo "================================="
echo "Extracting benchmark results..."
echo "================================="

# Extract timing information with error checking
ON_TIME=$(grep "実行時間:" benchmark_on.log | awk '{print $2}')
OFF_TIME=$(grep "実行時間:" benchmark_off.log | awk '{print $2}')
ON_PARTICLES=$(grep "粒子数:" benchmark_on.log | awk '{print $2}')
OFF_PARTICLES=$(grep "粒子数:" benchmark_off.log | awk '{print $2}')
ON_STEPS=$(grep "計算ステップ数:" benchmark_on.log | awk '{print $2}')
OFF_STEPS=$(grep "計算ステップ数:" benchmark_off.log | awk '{print $2}')

# Check if extraction was successful
if [ -z "$ON_TIME" ] || [ -z "$OFF_TIME" ] || [ -z "$ON_PARTICLES" ] || [ -z "$OFF_PARTICLES" ]; then
    echo "Error: Failed to extract benchmark data"
    echo "ON_TIME=$ON_TIME, OFF_TIME=$OFF_TIME"
    echo "ON_PARTICLES=$ON_PARTICLES, OFF_PARTICLES=$OFF_PARTICLES"
    exit 1
fi

# Convert scientific notation to decimal
ON_TIME_DECIMAL=$(echo "$ON_TIME" | awk '{printf "%.6f", $1}')
OFF_TIME_DECIMAL=$(echo "$OFF_TIME" | awk '{printf "%.6f", $1}')

echo "Results Summary:"
echo "----------------------------------------"
echo "Cell Algorithm ON:"
echo "  Particles: $ON_PARTICLES"
echo "  Execution Time: $ON_TIME_DECIMAL seconds"
echo "  Steps: $ON_STEPS"
echo ""
echo "Cell Algorithm OFF:"
echo "  Particles: $OFF_PARTICLES"
echo "  Execution Time: $OFF_TIME_DECIMAL seconds"
echo "  Steps: $OFF_STEPS"
echo ""

# Calculate speedup ratio with proper formatting
SPEEDUP=$(echo "scale=3; $OFF_TIME_DECIMAL / $ON_TIME_DECIMAL" | bc -l)
echo "Speed-up ratio: ${SPEEDUP}x"

# Save results to CSV format
RESULT_CSV="benchmark_summary.csv"
echo "Algorithm,Particles,Time_seconds,Steps,Time_per_step" > $RESULT_CSV
ON_TIME_PER_STEP=$(echo "scale=6; $ON_TIME_DECIMAL / $ON_STEPS" | bc -l | awk '{printf "%.6f", $1}')
OFF_TIME_PER_STEP=$(echo "scale=6; $OFF_TIME_DECIMAL / $OFF_STEPS" | bc -l | awk '{printf "%.6f", $1}')
echo "ON,$ON_PARTICLES,$ON_TIME_DECIMAL,$ON_STEPS,$ON_TIME_PER_STEP" >> $RESULT_CSV
echo "OFF,$OFF_PARTICLES,$OFF_TIME_DECIMAL,$OFF_STEPS,$OFF_TIME_PER_STEP" >> $RESULT_CSV

echo "Results also saved to: $RESULT_CSV"

echo "================================="
echo "Benchmark completed at $(date)"
echo "=================================" 