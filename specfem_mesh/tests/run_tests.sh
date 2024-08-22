#!/bin/bash

# -------------------------------------------------------
cd plm 
echo "Running test_plm..."
bash test_plm.sh
cd .. 
# -------------------------------------------------------
cd ylm 
echo "Running test_ylm..."
bash test_ylm.sh
cd .. 
# -------------------------------------------------------
cd integration 
echo "Running test_ylm_integration..."
bash test_ylm_integration.sh
bash test_1_integration.sh
cd .. 
# -------------------------------------------------------
cd spline 
echo "Running test_spline..."
bash test_spline.sh
cd .. 

# -------------------------------------------------------
# Run the pytests
pytest

# -------------------------------------------------------
# cleanup 
export test_dir=$(pwd)
source ./plm/cleanup.sh
source ./ylm/cleanup.sh
source ./integration/cleanup.sh
source ./spline/cleanup.sh
