
program compute_splitting_effect
    ! The aim of this program is to be a generic (non-optimised) code that is flexible
    ! for different models, computing a general splitting matrix: 
    use params, only: nprocs, Vcen
    use allocation_module
    use mesh_utils
    use ylm_plm
    use mineos_model, only: mineos, mineos_ptr
    use specfem_mesh, only: SetMesh, create_SetMesh
    use m_KdTree, only: KdTree, KdTreeSearch


    use model3d, only: M3D
    
    use modes, only: Mode, get_mode
    implicit none 
    
    type(SetMesh) :: sm 
    real(kind=8), allocatable :: globalvariable(:)
    integer :: iset , l1 ,n1, l2, n2, im 
    integer, parameter :: region = 3
    character(len=12)  :: strainname
    character   :: t1, t2
    
    type(Mode)  :: mode_1
    type(KdTree)           :: tree
    
    type(M3D) :: model3D
    
    real(kind=8), allocatable :: gpsi(:,:,:,:,:)
    real(kind=8) :: ggpsi(3,3)

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos
    
    
    ! Allocate strain matrix
    n1 = 0
    l1 = 2
    t1 = 'S'

    n2 = 0
    l2 = 2
    t2 = 'S'
    
    allocate(Vcen(2*l1+1, 2*l2 +1))


    ! Read text 3D model:  
    !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/gladm35/perturbation/dvp.txt"
    !call Model3D%read_model_from_file()

    Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
    call Model3D%read_model_from_file()
    call Model3D%create_KDtree()
    

    tree = KdTree(Model3D%xcoord, Model3D%ycoord, Model3D%zcoord) 


    Vcen = (0.0d0,0.0d0)

    write(*,*)'NPROCS: ', nprocs

    do iset = 0, nprocs -1 

            sm = create_SetMesh(iset, region) ! ic = 3
    
            ! Read the mesh info and coordinates
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            call sm%setup_gll()

            call sm%compute_jacobian(.false.)
            call sm%compute_wglljac(.false.)
            call sm%setup_global_coordinate_arrays(.false.)
            call sm%compute_rtp_from_xyz(.false.)
            call sm%get_unique_radii(.true.)



            call sm%compute_rotation_matrix()

            !allocate(gpsi(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
            !call compute_grad_centrifugal(sm, gpsi, ggpsi)

            !call compute_Vcen_matrix(sm, sm%interp, gpsi, n1, t1, l1, n2, t2, l2)
            call create_ensight_file_prefix(iset, region)
            call create_proc_case_file()
            call create_proc_geo_file(sm, 1)
      
            ! Useful for outputting to ensight: 
            allocate(globalvariable(sm%nglob))

            ! Project the 3D model to the mesh: 
            call Model3D%project_to_gll(sm, globalvariable, id=2)
  
            write(*,*)'Write scalar ', iset
            call write_real_scalar_to_ensight(sm, globalvariable, 'eta2', 1)
    
            call sm%cleanup()

            deallocate(globalvariable)
            !deallocate(gpsi)
            write(*,*)'Finished processor ', iset
    enddo 
    

    call save_Vcen_matrix(l1, l2, 'testvcen')

    
end program compute_splitting_effect
 
