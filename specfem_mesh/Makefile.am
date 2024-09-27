
# Look in src for source files 
VPATH  = src

# -cpp allows linknig to the .cuf file
# -D WITH_CUDA to set preprocessor WITH_CUDA = true
if USE_CUDA
    AM_FCFLAGS = -cpp -Wall -O4 -llapack -lblas -cuda -lcudart -D WITH_CUDA
else
    AM_FCFLAGS = -cpp -Wall -O4 -llapack -lblas -cuda -lcudart
endif

# Append flag for MPI
if USE_OPENMPI
    AM_FCFLAGS += -D WITH_MPI -Ktrap=inv
endif 




bin_PROGRAMS = split_mesh read_mesh compute_vani linbreak_mesh

# Basic files for all programs 
SETUP_SRCS = params.f90 allocation.f90 math.f90 output.f90 mesh_utils.f90 ylm_plm.f90 spline.f90  projection.f90 


# Specify that vani_kernel.o needs to be built before anything else
BUILT_SOURCES = params.o vani_kernel.o



# Create params.f90 from params.f90.in by replacing nmodes
params: src/params.in
	rm -f src/params.f90 &&  \
    $(AM_V_GEN) sed -e 's/@n_sets@/$(N_SETS)/g' src/params.in > $@ &&\
     sed -i -e 's/@n_modes@/$(N_MODES)/g' params && \
     mv params src/params.f90


# Rule to build params.o from params.f90
params.o: params.f90
	$(FC) -c $(AM_FCFLAGS) src/params.f90 -o params.o


# Rule to build vani_kernel.o from vani_kernel.cuf
vani_kernel.o: src/params.f90 src/vani_kernel.cuf
	$(FC) -c $(AM_FCFLAGS) src/vani_kernel.cuf -o vani_kernel.o


linbreak_mesh_SOURCES = $(SETUP_SRCS) linearly_breakup_mesh.f90 

split_mesh_SOURCES    = $(SETUP_SRCS) split_mesh.f90 

read_mesh_SOURCES     = $(SETUP_SRCS) gll.f90  mineos_model.f90 visual.f90  read_mesh.f90 


compute_vani_SOURCES = $(SETUP_SRCS) \
                        gll.f90 \
                        integrate.f90 \
                        get_mode.f90 \
                        mineos_model.f90 \
                        cuda_proxies.f90\
                        V_ani_matrix.f90 \
                        vani_kernel.o\
                        compute_ani_splitting.f90
                        
compute_vani_LDADD = vani_kernel.o

# Specify where to find the object file (if necessary)
EXTRA_DIST = vani_kernel.o

