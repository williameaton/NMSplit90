program read_specfem_mesh
use params
use allocation_module
use mesh_utils
use gll, only: setup_gll
use ylm_plm

implicit none 

real(kind=8), allocatable :: proc_id(:)
integer :: iproc             
integer :: region                ! Region code
integer :: n, l, m 


region = 3

n = 10
l = 5
m = 2

! Read mineos model 
call process_mineos_model(.true.)



do iproc = 0, nprocs -1 

        ! Read the mesh info and coordinates
        call read_proc_coordinates(iproc, region)
        allocate(proc_id(nglob))

        call load_ibool(iproc, region)

        call setup_gll()
        call compute_jacobian(iproc, .false.)

        ! Get unique mesh radii that are present
        call get_mesh_radii(iproc, .false., NR)
        
        ! We will need the r theta phi coordinates: 
        call compute_rtp_from_xyz(iproc, .false.)

        ! needed for ensight geo file
        call allocate_if_unallocated(nglob, x_glob)
        call allocate_if_unallocated(nglob, y_glob)
        call allocate_if_unallocated(nglob, z_glob)
        call map_local_global(xstore, x_glob, 0)
        call map_local_global(ystore, y_glob, 0)
        call map_local_global(zstore, z_glob, 0)



        proc_id = real(iproc,kind=8) 

        call create_ensight_file_prefix(iproc, region)
        call create_proc_case_file()
        call create_proc_geo_file(1)
        call write_real_scalar_to_ensight(proc_id, 'proc_id', 1)

        call cleanup_for_mode()

        deallocate(proc_id)


        write(*,*)'Finished processor ', iproc
enddo 




end program



