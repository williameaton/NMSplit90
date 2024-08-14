#!/bin/bash

# Compile the test executable
gfortran -I../src/ ./test_plm.f90 ../src/ylm_plm.f90 -o test_plm

# Run
./test_plm

# Doesnt run yet with but should test correctly 
python3 -m test_plm.py 

# cleanup 
rm testplm.txt
rm ylm_plm.mod
rm test_plm
rm -r __pycache__