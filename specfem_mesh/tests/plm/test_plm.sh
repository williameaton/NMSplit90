#!/bin/bash

# Compile the test executable
gfortran -I../../src/ ./test_plm.f90  ../../src/math.f90 ../../src/ylm_plm.f90 -o test_plm

# Run
./test_plm
