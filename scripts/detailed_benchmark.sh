#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr20001a
#SBATCH -t 24:00:00
#SBATCH --rsc p=1:t=8:c=8:m=32G
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -J pem_detailed_benchmark

#============ Environment Setup ============
module load intel

#============ Shell Script ============
echo "================================="
echo "PEM Detailed Algorithm Benchmark"
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

# Create result file
RESULT_FILE="benchmark_results.csv"
echo "Particle_Layers,Particles,Algorithm,Time_seconds,Steps,Time_per_step" > $RESULT_FILE

# Test different particle numbers
LAYER_COUNTS=(4 6 8 10 12)
STEPS=30000

for LAYERS in "${LAYER_COUNTS[@]}"; do
    echo "================================="
    echo "Testing with $LAYERS layers"
    echo "================================="
    
    # Create input files for this layer count
    cat > input_test_on.dat << EOF
# Test configuration - Cell Algorithm ON
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
    
    cat > input_test_off.dat << EOF
# Test configuration - Cell Algorithm OFF
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
    ./pem_simulator input_test_on.dat > test_on_${LAYERS}.log 2>&1
    
    # Extract results with error checking
    ON_TIME=$(grep "実行時間:" test_on_${LAYERS}.log | awk '{print $2}')
    ON_PARTICLES=$(grep "粒子数:" test_on_${LAYERS}.log | awk '{print $2}')
    ON_STEPS=$(grep "計算ステップ数:" test_on_${LAYERS}.log | awk '{print $2}')
    
    # Check if extraction was successful
    if [ -z "$ON_TIME" ] || [ -z "$ON_PARTICLES" ] || [ -z "$ON_STEPS" ]; then
        echo "Error: Failed to extract data from test_on_${LAYERS}.log"
        echo "ON_TIME=$ON_TIME, ON_PARTICLES=$ON_PARTICLES, ON_STEPS=$ON_STEPS"
        continue
    fi
    
    # Convert scientific notation to decimal and calculate time per step
    ON_TIME_DECIMAL=$(echo "$ON_TIME" | awk '{printf "%.6f", $1}')
    ON_TIME_PER_STEP=$(echo "scale=6; $ON_TIME_DECIMAL / $ON_STEPS" | bc -l | awk '{printf "%.6f", $1}')
    
    # Output to CSV with consistent formatting
    echo "$LAYERS,$ON_PARTICLES,ON,$ON_TIME_DECIMAL,$ON_STEPS,$ON_TIME_PER_STEP" >> $RESULT_FILE
    
    # Run with cell algorithm OFF (only for smaller particle counts)
    if [ $LAYERS -le 8 ]; then
        echo "Running with Cell Algorithm OFF ($LAYERS layers)..."
        ./pem_simulator input_test_off.dat > test_off_${LAYERS}.log 2>&1
        
        # Extract results with error checking
        OFF_TIME=$(grep "実行時間:" test_off_${LAYERS}.log | awk '{print $2}')
        OFF_PARTICLES=$(grep "粒子数:" test_off_${LAYERS}.log | awk '{print $2}')
        OFF_STEPS=$(grep "計算ステップ数:" test_off_${LAYERS}.log | awk '{print $2}')
        
        # Check if extraction was successful
        if [ -z "$OFF_TIME" ] || [ -z "$OFF_PARTICLES" ] || [ -z "$OFF_STEPS" ]; then
            echo "Error: Failed to extract data from test_off_${LAYERS}.log"
            echo "OFF_TIME=$OFF_TIME, OFF_PARTICLES=$OFF_PARTICLES, OFF_STEPS=$OFF_STEPS"
        else
            # Convert scientific notation to decimal and calculate time per step
            OFF_TIME_DECIMAL=$(echo "$OFF_TIME" | awk '{printf "%.6f", $1}')
            OFF_TIME_PER_STEP=$(echo "scale=6; $OFF_TIME_DECIMAL / $OFF_STEPS" | bc -l | awk '{printf "%.6f", $1}')
            
            # Output to CSV with consistent formatting
            echo "$LAYERS,$OFF_PARTICLES,OFF,$OFF_TIME_DECIMAL,$OFF_STEPS,$OFF_TIME_PER_STEP" >> $RESULT_FILE
            
            # Calculate speedup
            if [ -n "$ON_TIME_DECIMAL" ] && [ -n "$OFF_TIME_DECIMAL" ]; then
                SPEEDUP=$(echo "scale=2; $OFF_TIME_DECIMAL / $ON_TIME_DECIMAL" | bc -l)
                echo "Speedup for $LAYERS layers: ${SPEEDUP}x"
            fi
        fi
    else
        echo "Skipping Cell Algorithm OFF for $LAYERS layers (too many particles)"
    fi
    
    echo "Completed $LAYERS layers test"
done

# Clean up temporary files
rm -f input_test_on.dat input_test_off.dat

echo "================================="
echo "Benchmark Results Summary"
echo "================================="
cat $RESULT_FILE
echo "================================="
echo "Detailed benchmark completed at $(date)"
echo "Results saved to: $RESULT_FILE"
echo "=================================" 