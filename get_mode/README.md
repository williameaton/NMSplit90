
## get_mode

This directory contains a program called `test_modes`. The functionality of this program is as follows: 

- Load a mode from a MINEOS binary database
- Extract the eigenfunctions and frequencies
- Eigenfunctions U, U', V, V' or W, W' are then stored in file `mode.data`
- Tests normalisation of the eigenfunctions to ensure proper amplitudes and outputs norm which should be 1 


### Setup 
- You will need the following outputs from MINEOS: 
    1. Radial eigenfunction _out_plain_file_ and _out_bin_file_ 
    2. Toroidal eigenfunction _out_plain_file_ and _out_bin_file_ 
    3. Spheroidal eigenfunction _out_plain_file_ and _out_bin_file_ 
    4. Model file - the first 3 columns must be knot id, radius, density. The rest does not matter. No header lines can be included. The easiest way to get this is to take it from one of the out_plane)files and delete the eigenfunction content
    5. **note** that you must similarly remove the model rows of the _out_plain_files_ so that the first line starts with the first eigenfunction 

- You can then make the code using the `make` command. This is functioning for gfortran but needs testing with other compilers.

- Execute the code with `$ ./test_modes`. It will ask you to input the n and l values of the mode of interest, as well as the type, either spheroidal (S) or toroidal (T). 


### To-do: 

- Increase flexibility in binary file names/not hard-coded