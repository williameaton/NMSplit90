#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
elsec
    my_dir="${this_dir}/integration/"
fi

rm -rf ${my_dir}/__pycache__

rm -f ${my_dir}/*.o
rm -f ${my_dir}/test_ylm_int_*.txt
rm -f ${my_dir}/test_1_int.txt
rm -f ${my_dir}/test_r_int.txt