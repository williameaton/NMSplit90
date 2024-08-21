#!/bin/bash

src_dir="../../src"

# Compile the test executable
gfortran -llapack -lblas -I../../src/ ./test_1_integration.f90 \
    ${src_dir}/params.f90           \
    ${src_dir}/math.f90             \
    ${src_dir}/mineos_model.f90     \
    ${src_dir}/mesh_utils.f90       \
    ${src_dir}/get_mode.f90         \
    ${src_dir}/allocation.f90       \
    ${src_dir}/projection.f90       \
    ${src_dir}/spline.f90           \
    ${src_dir}/integrate.f90        \
    ${src_dir}/gll.f90              \
    ${src_dir}/ylm_plm.f90          \
    -o test_1_integration

# Run
./test_1_integration
