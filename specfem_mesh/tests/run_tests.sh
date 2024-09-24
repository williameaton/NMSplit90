#!/bin/bash

# Running with the flags 
# --nc --np so that it does not clean up the directories or run pytest
# since I want it to be a single pytest at the end in this case 

# -------------------------------------------------------
# bash test_vani.sh
bash run_spline.sh --nc --np
bash run_xlm.sh --nc --np
bash run_plm.sh --nc --np
bash run_ylm.sh --nc --np
bash run_rotmat.sh --nc --np
bash run_integration.sh --nc --np
bash run_mode_norm.sh --nc --np

# -------------------------------------------------------
# Run the pytests
pytest

# -------------------------------------------------------
# cleanup 
export this_dir=$(pwd)

#source ./v_ani_matrix/cleanup.sh
source ./plm/cleanup.sh
source ./xlm/cleanup.sh
source ./ylm/cleanup.sh
source ./rot_mat/cleanup.sh
source ./spline/cleanup.sh
source ./integration/cleanup.sh
source ./mode_normalisation/cleanup.sh