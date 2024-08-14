#!/bin/bash

# Compile the test executable
gfortran -I../src/ ./test_ylm.f90 ../src/ylm_plm.f90 -o test_ylm

# Run
./test_ylm

# Doesnt run yet with but should test correctly 
python3 -m test_ylm.py 

# cleanup 
rm testylm.txt
rm ylm_plm.mod
rm -r __pycache__
rm test_ylm