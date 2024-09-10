
# Compile the test executable
gfortran -I../src/ ./extract_mode.f90 \
    ../src/params.f90 \
    ../src/mineos_model.f90                \
    ../src/get_mode.f90 -o extract_mode

