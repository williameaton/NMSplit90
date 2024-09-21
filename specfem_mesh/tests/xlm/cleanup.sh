#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
else
    my_dir="${this_dir}/xlm"
fi

# cleanup 
rm -f   ${my_dir}/xlm_integral.txt
rm -rf ${my_dir}/__pycache__
rm -f  ${my_dir}/test_xlm.o