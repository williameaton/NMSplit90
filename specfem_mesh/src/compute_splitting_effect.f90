
program compute_splitting_effect
    ! The aim of this program is to be a generic (non-optimised) code that is flexible
    ! for different models, computing a general splitting matrix: 
    use params, only: nprocs
    use allocation_module
    use mesh_utils
    use ylm_plm
    use mineos_model, only: mineos, mineos_ptr
    use specfem_mesh, only: SetMesh, create_SetMesh


    use model3d, only: M3D
    
    use modes, only: Mode, get_mode
    implicit none 
    
    type(SetMesh) :: sm 
    real(kind=8), allocatable :: globalvariable(:)
    integer :: iset             , l1 ,n1, im 
    integer, parameter :: region = 1
    character(len=12)  :: strainname
    character   :: t1
    
    type(Mode)  :: mode_1
    
    type(M3D) :: model3D
    
    real(kind=8), allocatable :: gpsi(:,:,:,:,:)
    real(kind=8) :: ggpsi(3,3)

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos
    
    ! Allocate strain matrix
    n1 = 13
    l1 = 2
    t1 = 'S'
    


    ! Read text 3D model:  
    !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/gladm35/perturbation/dvp.txt"
    !call Model3D%read_model_from_file()

    ! Create a KD tree based on the points in the 3D model
    !call Model3D%create_KDtree()
    
    
    write(*,*)'number of procs: ', nprocs
    do iset = 0, nprocs -1 
    
            sm = create_SetMesh(iset, region) ! ic = 3
    
            ! Read the mesh info and coordinates
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            call sm%setup_gll()
            call sm%compute_jacobian(.false.)
            call sm%setup_global_coordinate_arrays(.false.)
            call sm%compute_rtp_from_xyz(.false.)
            call sm%get_unique_radii(.true.)
    

            allocate(gpsi(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
            call compute_grad_centrifugal(sm, gpsi, ggpsi)

            call create_ensight_file_prefix(iset, region)
            call create_proc_case_file()
            call create_proc_geo_file(sm, 1)
      

            ! Useful for outputting to ensight: 
            allocate(globalvariable(sm%nglob))

            ! Project the 3D model to the mesh: 
            call Model3D%project_to_gll(sm, globalvariable, id=1)
  

            call write_real_scalar_to_ensight(sm, globalvariable, 'example', 1)


    
            call sm%cleanup()

            deallocate(globalvariable)
            write(*,*)'Finished processor ', iset
    enddo 
    
    
end program compute_splitting_effect
 
