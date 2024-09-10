#!/bin/bash

this_dir=${test_dir}/rotation_matrix/

# cleanup 
rm ${this_dir}/*.mod
rm ${this_dir}/rotation_semi_analytical
rm ${this_dir}/test_rotation
rm ${this_dir}/test_stored_Wmat
rm -r ${this_dir}/__pycache__
rm ${this_dir}/*.txt
