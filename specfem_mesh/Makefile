# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -Wall -O2

# Directories
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

# Source files
SRC = $(wildcard $(SRC_DIR)/*.f90)

# Object files
OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC))

# Executable name
EXEC = read_mesh

# Default target
all: $(EXEC)

# Link the object files to create the executable
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

# Compile params.f90 first to generate params.mod
$(OBJ_DIR)/params.o: $(SRC_DIR)/params.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -c $< -o $@

# Compile mesh_utils.f90 after params.mod is generated
$(OBJ_DIR)/mesh_utils.o: $(SRC_DIR)/mesh_utils.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Compile read_mesh.f90 after params.mod is generated
$(OBJ_DIR)/read_mesh.o: $(SRC_DIR)/read_mesh.f90 $(OBJ_DIR)/params.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c $< -o $@

# Ensure the obj and mod directories exist
$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

# Clean up compiled files
clean:
	rm -f $(OBJ) $(EXEC)
	rm -rf $(OBJ_DIR) $(MOD_DIR)

.PHONY: all clean
