program read_specfem_mesh
use params 
use allocation_module
use mesh_utils
use gll 
use ylm_plm

implicit none 

integer :: iproc             
integer :: nprocs                ! Number of processors used by mesher 
integer :: region                ! Region code

integer :: n, l, m 


! Setup parameters: 
region = 3      ! Inner core
nprocs = 6

n = 10
l = 5
m = 2





! Read mineos model 
call process_mineos_model()

do iproc = 0, 5

        ! Read the mesh info and coordinates
        call read_proc_coordinates(iproc, region)

        ! Load ibool variable: 
        call load_ibool(iproc, region)

        call setup_gll()
        call compute_jacobian()

        ! Get unique mesh radii that are present
        call get_mesh_radii()

        ! We will need the r theta phi coordinates: 
        call compute_rtp_from_xyz()

        ! Compute strain tensor for mode
        allocate(strain1(6,ngllx, nglly, ngllz, nspec), & 
                globalstrain(6, nglob))
        call compute_gll_mode_strain('S', n, l, m, strain1)

        call map_complex_vector(6, strain1, globalstrain, 0)

        ! needed for ensight geo file
        call allocate_if_unallocated(nglob, x_glob)
        call allocate_if_unallocated(nglob, y_glob)
        call allocate_if_unallocated(nglob, z_glob)
        call map_local_global(xstore, x_glob, 0)
        call map_local_global(ystore, y_glob, 0)
        call map_local_global(zstore, z_glob, 0)


        !call compute_rotation_matrix()
        !call allocate_if_unallocated(nglob, theta_glob)
        !call allocate_if_unallocated(nglob, phi_glob)
        !call map_local_global(thetastore, theta_glob, 0)
        !call map_local_global(phistore,   phi_glob,   0)


        call create_ensight_file_prefix(iproc, region)
        call create_proc_case_file()
        call create_proc_geo_file(1)
        call write_complex_symtensor_to_ensight(globalstrain, 'strain', 1)

        ! Compute global mode displacement 
        !call allocate_if_unallocated(3, nglob, globaldisp)
        !call allocate_if_unallocated(3, ngllx, nglly, ngllz, nspec, disp1)

        !call compute_global_mode_displacement('S', n, l, m, disp1)
        !call map_complex_vector(3, disp1, globaldisp, 0)
        !call write_complex_vector_to_ensight(globaldisp, 'disp', 1)


        call cleanup_for_mode()

        write(*,*)'Finished processor ', iproc
enddo 




end program



