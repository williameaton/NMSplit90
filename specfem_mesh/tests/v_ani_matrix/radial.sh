#!/bin/bash

src_dir="../../src"

# Compile the SEM executable
gfortran -llapack -lblas -I../../src/ ./tromp_1995.f90 \
    ${src_dir}/params.f90           \
    ${src_dir}/math.f90             \
    ${src_dir}/mineos_model.f90     \
    ${src_dir}/mesh_utils.f90       \
    ${src_dir}/get_mode.f90         \
    ${src_dir}/allocation.f90       \
    ${src_dir}/projection.f90       \
    ${src_dir}/spline.f90           \
    ${src_dir}/integrate.f90        \
    ${src_dir}/visual.f90        \
    ${src_dir}/ylm_plm.f90          \
    ${src_dir}/V_ani_matrix.f90  \
    ${src_dir}/gll.f90  \
    ${src_dir}/output.f90              \
    -o radial

./radial
