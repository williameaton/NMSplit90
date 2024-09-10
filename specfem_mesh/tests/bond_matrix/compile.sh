#!/bin/bash

# Compile the test executable
gfortran -I../../src/ ./bond_matrix.f90  ../../src/math.f90 ../../src/V_ani_matrix.f90 -o bond_matrix

# Run
./bond_matrix

