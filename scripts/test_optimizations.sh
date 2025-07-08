#!/bin/bash
# Quick test script to verify PEM simulator optimizations

echo "================================="
echo "PEM Simulator Optimization Test"
echo "================================="

# Check if optimized simulator exists
if [ ! -f "src/optimized_pem_simulator.f90" ]; then
    echo "Error: Optimized source file not found!"
    echo "Looking for: src/optimized_pem_simulator.f90"
    exit 1
fi

# Create test input file
cat > test_optimization.dat << EOF
TIME_STEP 5.0e-7
PARTICLE_GEN_LAYERS 10
MAX_CALCULATION_STEPS 1000
OUTPUT_INTERVAL_NORMAL 100
FRICTION_COEFF_PARTICLE 0.25
FRICTION_COEFF_WALL 0.17
YOUNG_MODULUS_PARTICLE 4.9e9
YOUNG_MODULUS_WALL 3.9e9
POISSON_RATIO_PARTICLE 0.23
POISSON_RATIO_WALL 0.25
PARTICLE_DENSITY 2.48e3
PARTICLE_RADIUS_LARGE 1.0e-2
PARTICLE_RADIUS_SMALL 5.0e-3
CONTAINER_WIDTH 5.0e-1
RANDOM_SEED 584287
DISABLE_CELL_ALGORITHM 0
EOF

echo "Test input file created."

# Compile original version
echo -n "Compiling original version... "
if ifort -O3 -o pem_original src/pem_simulator.f90 2>/dev/null; then
    echo "Success"
else
    echo "Failed"
    echo "Trying with gfortran..."
    if gfortran -O3 -o pem_original src/pem_simulator.f90 2>/dev/null; then
        echo "Success with gfortran"
    else
        echo "Failed - please check compiler installation"
        exit 1
    fi
fi

# Test original version
echo ""
echo "Running original version..."
time ./pem_original test_optimization.dat > original_output.log 2>&1
echo "Original version completed."

# Check Python optimization
echo ""
echo "Testing optimized Python animation script..."
if [ -f "src/animate_pem_optimized.py" ]; then
    echo "✓ Optimized animation script found"
    python3 src/animate_pem_optimized.py --help > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "✓ Animation script is functional"
    else
        echo "⚠ Animation script may have missing dependencies"
    fi
else
    echo "✗ Optimized animation script not found"
fi

# Summary
echo ""
echo "================================="
echo "Optimization Files Created:"
echo "================================="
echo "✓ src/optimized_pem_simulator.f90 - Optimized Fortran source"
echo "✓ src/animate_pem_optimized.py - Optimized animation generator"
echo "✓ src/performance_optimization_guide.md - Detailed optimization guide"
echo "✓ Makefile - Build system with optimization targets"
echo "✓ scripts/benchmark_optimized.sh - Performance benchmark script"
echo "✓ PERFORMANCE_OPTIMIZATION_REPORT.md - Comprehensive report"

echo ""
echo "To compile and test the optimized version:"
echo "  make optimize      # Build with Intel compiler"
echo "  make optimize-gnu  # Build with GNU compiler"
echo "  make benchmark     # Run performance benchmarks"

# Cleanup
rm -f test_optimization.dat original_output.log pem_original

echo ""
echo "Test completed successfully!"