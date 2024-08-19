# NMSplit90

This respository hosts some F90 programs to compute splitting functions of 
normal modes due to Inner Core anisotropy. Some of the original code
was translated from Jeroen Tromp's old F77 code, hence the emphasis
on '90'. 


### Repository overview and structure

Within this repository are two directories, described below:

- #### `get_mode` 
    - This directory contains the program `test_modes`. The majority of
    this code was translated from F77 with some updates to the MINEOS
    read format to enable compatibility with the most recent version of
    MINEOS. 
    - The program loads a mode and conducts radial integration to test
    the normalisation of the radial eigenfunction in accordance with
    the MINEOS normalisation conventions

- #### `specfem_mesh` 
     This directory contains a number of programs
     - `read_mesh`\
        This is the main driver program, which evidently needs renaming. 
        So far, this code has the following functionality: 
        - Reads MINEOS radial knots
        - Reads SPECFEM mesh 
        - Loads mode eigenfunctions
        - Computes mode displacement at each GLL point in the mesh
        - Computes mode strain at each GLL point in the mesh
     - `split_mesh`\
        The goal of this package is to compute splitting functions due to
        inner core anisotropy, which requires evaluation of the integral 
        $$ \int_{V} \mathbf{\tilde{\epsilon}}_k^{*} : \mathbf{\gamma} : \mathbf{\tilde{\epsilon}}_{k'} \: \mathrm{d} V  $$ 
        between two modes $k$ and $k'$. We therefore need to integrate 
        over the entire inner core.
        - The current problem is that while the SPECFEM\_GLOBE mesher will 
        do the load balancing and mesh generation, the innermost cube
        of the inner core is shared between each process. 
        - Note that this is true when using 6 processes. If you use more 
        processes then only parts of the central cube are shared between
        different processes. 
        - For simplicity, I therefore use 6 here. 
        - Split mesh therefore processes these mesh files so that the central
          cube is only stored on one process, while the other processes contain
          only the parts that are non-overlapping in the original mesh.
        

### Runtime considerations

- The mesh is read in from the directory hard coded in `params.f90`, stored
  in the variable `datadir`. If you run `split_mesh` it will output the 
  new mesh files into `${datadir}/sliced`. Once you have run this program, 
  you would then need to edit your `datadir` to `${datadir}/sliced` so that
  it is then loading the new mesh from the `sliced` directory. 
        

       




