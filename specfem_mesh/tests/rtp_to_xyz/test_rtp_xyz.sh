#!/bin/bash

# Compile the test executable
gfortran -I../../src/ ./test_rtp_xyz.f90 \
     ../../src/params.f90 \
     ../../src/math.f90 \
    ../../src/allocation.f90 \
     ../../src/mesh_utils.f90 \
     -o test_rtp_xyz

# Run
./test_rtp_xyz
