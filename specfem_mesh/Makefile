# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -Wall -O2 -llapack -lblas -Wpedantic


# Directories
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

# Source files
SRC = $(wildcard $(SRC_DIR)/*.f90)

# Object files (excluding main program files)
COMMON_OBJ = $(filter-out $(OBJ_DIR)/integrate_mesh.o $(OBJ_DIR)/read_mesh.o $(OBJ_DIR)/split_mesh.o, $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC)))

# Executable names
EXEC_READ = read_mesh
EXEC_SPLIT = split_mesh
EXEC_INTEG = integrate_mesh

# Default target
all: $(EXEC_READ) $(EXEC_SPLIT) $(EXEC_INTEG)

# Link the object files to create the read_mesh executable
$(EXEC_READ): $(COMMON_OBJ) $(OBJ_DIR)/read_mesh.o
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJ) $(OBJ_DIR)/read_mesh.o

# Link the object files to create the split_mesh executable
$(EXEC_SPLIT): $(COMMON_OBJ) $(OBJ_DIR)/split_mesh.o
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJ) $(OBJ_DIR)/split_mesh.o

$(EXEC_INTEG): $(COMMON_OBJ) $(OBJ_DIR)/integrate_mesh.o
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJ) $(OBJ_DIR)/integrate_mesh.o



# Compile params.f90 first to generate params.mod
$(OBJ_DIR)/params.o: $(SRC_DIR)/params.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -c $< -o $@

# Compile math.f90 after params.mod is generated
$(OBJ_DIR)/math.o: $(SRC_DIR)/math.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR)  -J$(MOD_DIR) -c $< -o $@

# Compile mesh_utils.f90 after params.mod is generated
$(OBJ_DIR)/mesh_utils.o: $(SRC_DIR)/mesh_utils.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR)  -J$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/gll.o: $(SRC_DIR)/gll.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR)  -J$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/allocation.o: $(SRC_DIR)/allocation.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

# Compile get_mode.f90 after params.mod is generated
$(OBJ_DIR)/get_mode.o: $(SRC_DIR)/get_mode.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile ylm_plm.f90 after params.mod is generated
$(OBJ_DIR)/ylm_plm.o: $(SRC_DIR)/ylm_plm.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

# Compile spline.f90 after params.mod is generated
$(OBJ_DIR)/spline.o: $(SRC_DIR)/spline.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR)  -c $< -o $@

# Compile projection.f90 after params.mod is generated
$(OBJ_DIR)/projection.o: $(SRC_DIR)/projection.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile mineos_model.f90 after params.mod is generated
$(OBJ_DIR)/mineos_model.o: $(SRC_DIR)/mineos_model.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile visual.f90 after params.mod is generated
$(OBJ_DIR)/visual.o: $(SRC_DIR)/visual.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile integrate.f90 after params.mod is generated
$(OBJ_DIR)/integrate.o: $(SRC_DIR)/integrate.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

# Compile read_mesh.f90 after params.mod is generated
$(OBJ_DIR)/read_mesh.o: $(SRC_DIR)/read_mesh.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile split_mesh.f90 after params.mod is generated
$(OBJ_DIR)/split_mesh.o: $(SRC_DIR)/split_mesh.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile integrate_mesh.f90 after params.mod is generated
$(OBJ_DIR)/integrate_mesh.o: $(SRC_DIR)/integrate_mesh.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@


# Compile integrate_mesh.f90 after params.mod is generated
$(OBJ_DIR)/V_ani_matrix.o: $(SRC_DIR)/V_ani_matrix.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile integrate_mesh.f90 after params.mod is generated
$(OBJ_DIR)/W_matrix.o: $(SRC_DIR)/W_matrix.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile integrate_mesh.f90 after params.mod is generated
$(OBJ_DIR)/output.o: $(SRC_DIR)/output.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@



# Ensure the obj and mod directories exist
$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

# Clean up compiled files
clean:
	rm -f $(OBJ_DIR)/*.o $(EXEC_READ) $(EXEC_SPLIT) $(EXEC_INTEG)
	rm -rf $(OBJ_DIR) $(MOD_DIR)

.PHONY: all clean
