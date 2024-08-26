#!/bin/bash

this_dir=${test_dir}/mode_normalisation/

# cleanup 
# Note that you should not delete spline.mod
rm ${this_dir}/*T*.txt
rm ${this_dir}/*spline*.txt
rm ${this_dir}/params.mod
rm ${this_dir}/math.mod
rm ${this_dir}/allocation_module.mod
rm ${this_dir}/integrate.mod
rm ${this_dir}/mesh_utils.mod
rm ${this_dir}/ylm_plm.mod
rm ${this_dir}/test_norm
rm -r ${this_dir}/__pycache__
rm ${this_dir}/status.txt