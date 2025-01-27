program read_specfem_mesh
use params, only: nprocs
use allocation_module
use mesh_utils
use ylm_plm
use mineos_model, only: mineos, mineos_ptr
use specfem_mesh, only: SetMesh, create_SetMesh

use modes, only: Mode, get_mode
implicit none 

type(SetMesh) :: sm 
real(kind=8), allocatable :: globstrain(:)
integer :: iset             , l1 ,n1, im 
integer, parameter :: region = 1
character(len=12)  :: strainname
character   :: t1

type(Mode)    :: mode_1


! Read mineos model 
call mineos%process_mineos_model(.false.)
mineos_ptr => mineos

! Allocate strain matrix
n1 = 13
l1 = 2
t1 = 'S'


! Load model and coordinates 
tree = KdTree(vor_x, vor_y, vor_z) 




write(*,*)'number of procs: ', nprocs
do iset = 0, nprocs -1 

        sm = create_SetMesh(iset, region) ! ic = 3

        ! Read the mesh info and coordinates
        call sm%read_proc_coordinates()
        call sm%load_ibool()
        call sm%setup_global_coordinate_arrays(.false.)
        call sm%compute_rtp_from_xyz(.false.)
        call sm%get_unique_radii(.true.)

        allocate(globstrain(sm%nglob))

        call create_ensight_file_prefix(iset, region)
        call create_proc_case_file()
        call create_proc_geo_file(sm, 1)

        ! Get the mode: 
        mode_1 =  get_mode(n1, t1, l1, mineos_ptr)
        call sm%interp%interpolate_mode_eigenfunctions(mode_1)

        ! Load strain: 
        allocate(sm%strain1(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, 6))
        do im =  -l1, l1 


                ! Compute, for example, the mode strain
                call sm%compute_mode_strain(im, mode_1, sm%strain1)
                !call sm%load_mode_strain_binary(n1, 'S', l1, im, &  
                !                                strain1(:,:,:,:,:))
                ! Map to global: 
                call sm%map_local_global_double_precision(real(sm%strain1(:,:,:,:,1), kind=8), globstrain, 0) 
                
                if(im < 0)then 
                        write(strainname,'(a,i2)')'strain_m', im    
                else 
                        write(strainname,'(a,i1)')'strain_m', im  
                endif   
                call write_real_scalar_to_ensight(sm, globstrain, strainname, 1)
        enddo 


        call sm%cleanup()

        !deallocate(globstrain)
        deallocate(sm%strain1)
        write(*,*)'Finished processor ', iset
enddo 


end program



