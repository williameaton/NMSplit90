#!/bin/bash

this_dir=${test_dir}/ylm/

# cleanup 
rm ${this_dir}/ylm_value_error.txt
rm ${this_dir}/testylm.txt
rm ${this_dir}/ylm_plm.mod
rm ${this_dir}/math.mod
rm -r ${this_dir}/__pycache__
rm ${this_dir}/test_ylm
rm ${this_dir}/test_values