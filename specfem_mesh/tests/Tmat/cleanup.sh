#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
else
    my_dir="${this_dir}/Tmat"
fi

# cleanup 
rm -fr ${my_dir}/__pycache__
rm -f  ${my_dir}/*.o
rm -f   ${my_dir}/radial_0S2_0S2.txt
rm -f   ${my_dir}/Tmat_0S2_0S2.txt
rm -f  ${my_dir}/Wmat_stored_0S2_0S2.txt