
# Look in src for source files 
VPATH = src

# Ensure that it is reading/writing the mod files to mod dir with the -I -J
AM_FCFLAGS = -fopenmp  -Wall -O2 -llapack -lblas -cuda -lcudart

bin_PROGRAMS = split_mesh read_mesh compute_vani

# Basic files for all programs 
SETUP_SRCS = params.f90 allocation.f90 math.f90 mesh_utils.f90 ylm_plm.f90 spline.f90 output.f90 projection.f90 


split_mesh_SOURCES = $(SETUP_SRCS) split_mesh.f90 

read_mesh_SOURCES = $(SETUP_SRCS) gll.f90  mineos_model.f90 visual.f90  read_mesh.f90 


compute_vani_SOURCES = $(SETUP_SRCS) \
                        gll.f90 \
                        integrate.f90 \
                        get_mode.f90 \
                        mineos_model.f90 \
                        V_ani_matrix.f90 \
                        vani_kernel.o\
                        compute_ani_splitting.f90
                        
compute_vani_LDADD = vani_kernel.o

# Specify where to find the object file (if necessary)
EXTRA_DIST = vani_kernel.o

