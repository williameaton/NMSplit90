#!/bin/bash


# -------------------------------------------------------
cd v_ani_matrix
echo "Running v ani matrix test..."
bash test_vani.sh
cd .. 
# -------------------------------------------------------
cd rotation_matrix 
echo "Running rotation matrix tests..."
bash test_rotation.sh
cd .. 
# -------------------------------------------------------
cd xlm 
echo "Running test_xlm..."
bash test_xlm.sh
cd .. 
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
echo "Running test_1_integration..."
bash test_1_integration.sh
cd .. 
# -------------------------------------------------------
cd spline 
echo "Running test_spline..."
bash test_spline.sh
cd .. 
# -------------------------------------------------------
cd mode_normalisation 
echo "Running test_norm..."
bash test_norm.sh
cd .. 


# -------------------------------------------------------
# Run the pytests
pytest

# -------------------------------------------------------
# cleanup 
export test_dir=$(pwd)
source ./rotation_matrix/cleanup.sh
source ./v_ani_matrix/cleanup.sh
source ./plm/cleanup.sh
source ./xlm/cleanup.sh
source ./ylm/cleanup.sh
source ./integration/cleanup.sh
source ./spline/cleanup.sh
source ./mode_normalisation/cleanup.sh