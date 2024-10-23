program read_specfem_mesh
use params, only: nprocs
use allocation_module
use mesh_utils
use ylm_plm
use mineos_model, only: mineos, mineos_ptr
use specfem_mesh, only: SetMesh, create_SetMesh
implicit none 

type(SetMesh) :: sm 
real(kind=8), allocatable :: proc_id(:)
integer :: iset             
integer, parameter :: region = 0

! Read mineos model 
call mineos%process_mineos_model(.false.)


write(*,*)'number of procs: ', nprocs
do iset = 0, nprocs -1 

        sm = create_SetMesh(iset, region) ! ic = 3

        ! Read the mesh info and coordinates
        call sm%read_proc_coordinates()
        call sm%load_ibool()

        call sm%setup_global_coordinate_arrays(.false.)
        call sm%compute_rtp_from_xyz(.false.)

        allocate(proc_id(sm%nglob))

        proc_id = real(iset,kind=8) 

        call create_ensight_file_prefix(iset, region)
        call create_proc_case_file()
        call create_proc_geo_file(sm, 1)
        call write_real_scalar_to_ensight(sm,proc_id, 'proc_id', 1)

        call sm%cleanup()

        deallocate(proc_id)
        write(*,*)'Finished processor ', iset
enddo 


end program



