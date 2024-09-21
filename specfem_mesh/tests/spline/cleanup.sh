#!/bin/bash

if [ -z "$this_dir" ]; then
    my_dir='./'
else
    my_dir="${this_dir}/spline"
fi


rm -rf ${my_dir}/__pycache__
rm -f ${my_dir}/*.o
rm -f ${my_dir}/*23S12.txt
