#!/bin/bash

#============ Quick Test Benchmark for CSV Format Verification ============
echo "================================="
echo "PEM Quick Test Benchmark"
echo "Testing CSV output format..."
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
RESULT_FILE="quick_test_results.csv"
echo "Particle_Layers,Particles,Algorithm,Time_seconds,Steps,Time_per_step" > $RESULT_FILE

# Test with minimal particle counts for quick verification
LAYER_COUNTS=(4 6)
STEPS=5000  # Much shorter for quick testing

for LAYERS in "${LAYER_COUNTS[@]}"; do
    echo "================================="
    echo "Quick testing with $LAYERS layers"
    echo "================================="
    
    # Create input files for this layer count
    cat > input_quick_on.dat << EOF
# Quick test configuration - Cell Algorithm ON
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
DISABLE_CELL_ALGORITHM      0
CELL_SIZE_OVERRIDE          0.0
OUTPUT_INTERVAL_NORMAL      10000
OUTPUT_INTERVAL_VALIDATION  10
EOF
    
    cat > input_quick_off.dat << EOF
# Quick test configuration - Cell Algorithm OFF
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
    
    # Run with cell algorithm ON
    echo "Running with Cell Algorithm ON ($LAYERS layers)..."
    ./pem_simulator input_quick_on.dat > quick_on_${LAYERS}.log 2>&1
    
    # Extract results with error checking and clean up whitespace/newlines
    ON_TIME=$(grep "実行時間:" quick_on_${LAYERS}.log | awk '{print $2}' | sed 's/[[:space:]]*$//')
    ON_PARTICLES=$(grep "粒子数:" quick_on_${LAYERS}.log | tail -1 | awk '{print $2}' | sed 's/[[:space:]]*$//')
    ON_STEPS=$(grep "計算ステップ数:" quick_on_${LAYERS}.log | awk '{print $2}' | sed 's/[[:space:]]*$//')
    
    # Check if extraction was successful
    if [ -z "$ON_TIME" ] || [ -z "$ON_PARTICLES" ] || [ -z "$ON_STEPS" ]; then
        echo "Error: Failed to extract data from quick_on_${LAYERS}.log"
        echo "ON_TIME=$ON_TIME, ON_PARTICLES=$ON_PARTICLES, ON_STEPS=$ON_STEPS"
        continue
    fi
    
    # Convert scientific notation to decimal and calculate time per step
    ON_TIME_DECIMAL=$(echo "$ON_TIME" | awk '{printf "%.6f", $1}')
    ON_TIME_PER_STEP=$(echo "scale=6; $ON_TIME_DECIMAL / $ON_STEPS" | bc -l | awk '{printf "%.6f", $1}')
    
    # Output to CSV with explicit formatting
    echo "${LAYERS},${ON_PARTICLES},ON,${ON_TIME_DECIMAL},${ON_STEPS},${ON_TIME_PER_STEP}" >> $RESULT_FILE
    
    echo "Running with Cell Algorithm OFF ($LAYERS layers)..."
    ./pem_simulator input_quick_off.dat > quick_off_${LAYERS}.log 2>&1
    
    # Extract results with error checking and clean up whitespace/newlines
    OFF_TIME=$(grep "実行時間:" quick_off_${LAYERS}.log | awk '{print $2}' | sed 's/[[:space:]]*$//')
    OFF_PARTICLES=$(grep "粒子数:" quick_off_${LAYERS}.log | tail -1 | awk '{print $2}' | sed 's/[[:space:]]*$//')
    OFF_STEPS=$(grep "計算ステップ数:" quick_off_${LAYERS}.log | awk '{print $2}' | sed 's/[[:space:]]*$//')
    
    # Check if extraction was successful
    if [ -z "$OFF_TIME" ] || [ -z "$OFF_PARTICLES" ] || [ -z "$OFF_STEPS" ]; then
        echo "Error: Failed to extract data from quick_off_${LAYERS}.log"
        echo "OFF_TIME=$OFF_TIME, OFF_PARTICLES=$OFF_PARTICLES, OFF_STEPS=$OFF_STEPS"
    else
        # Convert scientific notation to decimal and calculate time per step
        OFF_TIME_DECIMAL=$(echo "$OFF_TIME" | awk '{printf "%.6f", $1}')
        OFF_TIME_PER_STEP=$(echo "scale=6; $OFF_TIME_DECIMAL / $OFF_STEPS" | bc -l | awk '{printf "%.6f", $1}')
        
        # Output to CSV with explicit formatting
        echo "${LAYERS},${OFF_PARTICLES},OFF,${OFF_TIME_DECIMAL},${OFF_STEPS},${OFF_TIME_PER_STEP}" >> $RESULT_FILE
        
        # Calculate speedup
        if [ -n "$ON_TIME_DECIMAL" ] && [ -n "$OFF_TIME_DECIMAL" ]; then
            SPEEDUP=$(echo "scale=2; $OFF_TIME_DECIMAL / $ON_TIME_DECIMAL" | bc -l)
            echo "Speedup for $LAYERS layers: ${SPEEDUP}x"
        fi
    fi
    
    echo "Completed $LAYERS layers test"
done

# Clean up temporary files
rm -f input_quick_on.dat input_quick_off.dat

echo "================================="
echo "Quick Test Results:"
echo "================================="
cat $RESULT_FILE
echo "================================="
echo "Testing CSV format with Python script..."

# Test the CSV with the analysis script
if [ -f "analyze_benchmark.py" ]; then
    python analyze_benchmark.py $RESULT_FILE
else
    echo "analyze_benchmark.py not found, skipping format verification"
fi

echo "================================="
echo "Quick test completed at $(date)"
echo "Results saved to: $RESULT_FILE"
echo "=================================" 