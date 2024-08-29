#!/bin/bash

# Compile the test executable
gfortran -I../../src/ ./test_ylm.f90 ../../src/math.f90 ../../src/ylm_plm.f90 -o test_ylm

# Run
./test_ylm




# Compile test_values executable
gfortran -I../../src/ ./test_values.f90 ../../src/math.f90 ../../src/ylm_plm.f90 -o test_values

# Run
./test_values
