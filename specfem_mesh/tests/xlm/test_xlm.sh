#!/bin/bash

# Compile the test executable
gfortran -I../../src/ ./test_xlm.f90 ../../src/params.f90 ../../src/math.f90 ../../src/integrate.f90  ../../src/ylm_plm.f90 -o test_xlm

# Run
./test_xlm
