#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
else
    my_dir="${this_dir}/ylm"
fi

# cleanup 
rm -f ${my_dir}/ylm_value_error.txt
rm -f ${my_dir}/testylm.txt
rm -f ${my_dir}/*.o
rm -rf ${my_dir}/__pycache__
