#!/bin/bash

#============ Complete Benchmark for Missing Data ============
echo "================================="
echo "PEM Complete Benchmark"
echo "Completing missing OFF data..."
echo "Start Time: $(date)"
echo "================================="

# Compile the Fortran code if needed
if [ ! -f "pem_simulator" ]; then
    echo "Compiling Fortran code..."
    cd src
    ifort -O3 -qopenmp -o pem_simulator pem_simulator.f90
    if [ $? -ne 0 ]; then
        echo "Compilation failed!"
        exit 1
    fi
    mv pem_simulator ..
    cd ..
    echo "Compilation successful!"
fi

# Create result file
RESULT_FILE="benchmark/results/complete_benchmark_results.csv"
echo "Particle_Layers,Particles,Algorithm,Time_seconds,Steps,Time_per_step" > $RESULT_FILE

# Missing layer counts for OFF algorithm
LAYER_COUNTS=(10 12)
STEPS=30001  # Same as original benchmark

for LAYERS in "${LAYER_COUNTS[@]}"; do
    echo "================================="
    echo "Running OFF algorithm with $LAYERS layers"
    echo "================================="
    
    # Create input file for this layer count with cell algorithm OFF
    cat > input_complete_off.dat << EOF
# Complete benchmark configuration - Cell Algorithm OFF
TIME_STEP                   5.0e-7
MAX_CALCULATION_STEPS       $STEPS
YOUNG_MODULUS_PARTICLE      4.9e9
YOUNG_MODULUS_WALL          3.9e9
POISSON_RATIO_PARTICLE      0.23
POISSON_RATIO_WALL          0.25
PARTICLE_DENSITY            2480.0
FRICTION_COEFF_PARTICLE     0.25
FRICTION_COEFF_WALL         0.17
VALIDATION_MODE             0
PARTICLE_RADIUS_LARGE       1.0e-2
PARTICLE_RADIUS_SMALL       5.0e-3
CONTAINER_WIDTH             0.5
PARTICLE_GEN_LAYERS         $LAYERS
RANDOM_SEED                 584287
DISABLE_CELL_ALGORITHM      1
CELL_SIZE_OVERRIDE          0.0
OUTPUT_INTERVAL_NORMAL      10000
OUTPUT_INTERVAL_VALIDATION  10
EOF
    
    echo "Running with Cell Algorithm OFF ($LAYERS layers)..."
    ./pem_simulator input_complete_off.dat > benchmark/logs/complete_off_${LAYERS}.log 2>&1
    
    # Extract results with error checking and clean up whitespace/newlines
    OFF_TIME=$(grep "実行時間:" benchmark/logs/complete_off_${LAYERS}.log | awk '{print $2}' | sed 's/[[:space:]]*$//')
    OFF_PARTICLES=$(grep "粒子数:" benchmark/logs/complete_off_${LAYERS}.log | tail -1 | awk '{print $2}' | sed 's/[[:space:]]*$//')
    OFF_STEPS=$(grep "計算ステップ数:" benchmark/logs/complete_off_${LAYERS}.log | awk '{print $2}' | sed 's/[[:space:]]*$//')
    
    # Check if extraction was successful
    if [ -z "$OFF_TIME" ] || [ -z "$OFF_PARTICLES" ] || [ -z "$OFF_STEPS" ]; then
        echo "Error: Failed to extract data from complete_off_${LAYERS}.log"
        echo "OFF_TIME=$OFF_TIME, OFF_PARTICLES=$OFF_PARTICLES, OFF_STEPS=$OFF_STEPS"
    else
        # Convert scientific notation to decimal and calculate time per step
        OFF_TIME_DECIMAL=$(echo "$OFF_TIME" | awk '{printf "%.6f", $1}')
        OFF_TIME_PER_STEP=$(echo "scale=6; $OFF_TIME_DECIMAL / $OFF_STEPS" | bc -l | awk '{printf "%.6f", $1}')
        
        # Output to CSV with explicit formatting
        echo "${LAYERS},${OFF_PARTICLES},OFF,${OFF_TIME_DECIMAL},${OFF_STEPS},${OFF_TIME_PER_STEP}" >> $RESULT_FILE
        
        echo "Completed $LAYERS layers test"
        echo "Time: ${OFF_TIME_DECIMAL}s, Particles: ${OFF_PARTICLES}, Steps: ${OFF_STEPS}"
    fi
done

# Clean up temporary files
rm -f input_complete_off.dat

echo "================================="
echo "Complete Benchmark Results:"
echo "================================="
cat $RESULT_FILE
echo "================================="

# Merge with existing results
echo "Merging with existing results..."
cat benchmark/results/benchmark_results_fixed.csv > benchmark/results/benchmark_results_complete.csv
tail -n +2 $RESULT_FILE >> benchmark/results/benchmark_results_complete.csv

echo "================================="
echo "Final Complete Results:"
echo "================================="
cat benchmark/results/benchmark_results_complete.csv
echo "================================="
echo "Complete benchmark finished at $(date)"
echo "Results saved to: benchmark/results/benchmark_results_complete.csv"
echo "=================================" 