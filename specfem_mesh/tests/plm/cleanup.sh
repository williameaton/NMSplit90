#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
else
    my_dir="${this_dir}/plm"
fi
rm -rf ${my_dir}/__pycache__
rm -f ${my_dir}/testplm.txt
rm -f ${my_dir}/plm_value_error.txt
rm -f ${my_dir}/test_plm.o
rm -f ${my_dir}/test_plm_values.o