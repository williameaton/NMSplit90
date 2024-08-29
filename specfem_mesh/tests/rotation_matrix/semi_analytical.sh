#!/bin/bash

src_dir="../../src"

# Compile the test executable
gfortran -llapack -lblas -I../../src/ ./rotation_semi_analytical.f90 \
    ${src_dir}/params.f90           \
    ${src_dir}/math.f90             \
    ${src_dir}/mineos_model.f90     \
    ${src_dir}/mesh_utils.f90       \
    ${src_dir}/get_mode.f90         \
    ${src_dir}/allocation.f90       \
    ${src_dir}/projection.f90       \
    ${src_dir}/spline.f90           \
    ${src_dir}/integrate.f90        \
    ${src_dir}/visual.f90           \
    ${src_dir}/ylm_plm.f90          \
    ${src_dir}/rotation_matrix.f90  \
    ${src_dir}/gll.f90              \
    -o rotation_semi_analytical


if [ $? -ne 0 ]; then
    echo "1" > status.txt
    exit 1
fi


./rotation_semi_analytical

# Check if the executable ran successfully
if [ $? -eq 0 ]; then
    echo "0" > status.txt
else
    echo "1" > status.txt
fi

