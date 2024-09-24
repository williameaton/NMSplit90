#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
else
    my_dir="${this_dir}/mode_normalisation"
fi

rm -f ${my_dir}/*T*.txt
rm -f ${my_dir}/*spline*.txt
rm -rf ${my_dir}/__pycache__
rm -f ${my_dir}/status.txt
rm -f ${my_dir}/*.o