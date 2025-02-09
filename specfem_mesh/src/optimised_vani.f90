
program optimised_vani
    use params, only: Vani, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                      MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                      glob_eta2,  nmodes, nprocs, all_warnings, datadir, max_tl1
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                     compute_Vani_matrix, compute_vani_matrix_stored, & 
                     convert_imag_to_real, save_Vani_real_matrix
    use splitting_function, only: get_Ssum_bounds, cst_to_H, H_to_cst, write_cst_to_file

    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
    use m_KdTree, only: KdTree, KdTreeSearch
    use voronoi, only: vor_x, vor_y, vor_z, & 
                       vor_A, vor_C, vor_L, vor_N, vor_F, &
                       load_voronoi_model, project_voroni_to_gll
    use specfem_mesh, only: SetMesh, create_SetMesh
    use modes, only: get_mode, Mode 
    use mineos_model, only: mineos, mineos_ptr
    use model3d, only: M3D

    implicit none
    include "constants.h"
    include 'mpif.h'

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2, n2, ierr, & 
               tl2, h, b, cluster_size, sets_per_process, & 
               myset_start, myset_end, i_mode, smin, smax, num_s, ncols, this_tl1, imodeliter
    character(len=2) nstr, lstr
    character(len=5) iterstr
    character(len=250) :: out_name

    real(kind=SPLINE_REAL), allocatable :: Vani_real(:,:)
    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)
    real(kind=SPLINE_REAL), allocatable :: cst(:,:)
    character(len=20) :: model_ti

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  
    type(M3D) :: model3D

    ! Modes: 
    integer, dimension(40), parameter :: modeNs = (/2, 3, 3, 5, 6, 8, 8, 9, 9, 9, 11, 11, 11, 11, 13, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21, 21, 22, 23, 23, 25, 25, 27/)
    integer, dimension(40), parameter :: modeLs = (/3, 1, 2, 2, 3, 1, 5, 2, 3, 4,  1,  4,  5,  6, 1,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7,  8,  1,  4,  5,  1,  2,  2/)


    ! Simulation parameters: 
    integer, parameter    :: region      = 3
    integer, parameter    :: nmodeliter  = 1000
    character, parameter  :: t1          = 'S'
    logical, parameter :: force_VTI      = .false.

    ! Setup MPI 
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, cluster_size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    if(SPLINE_REAL.eq.4)then 
         MPI_SPLINE_REAL    = MPI_REAL
         MPI_SPLINE_COMPLEX = MPI_COMPLEX
    elseif(SPLINE_REAL.eq.8)then
         MPI_SPLINE_REAL    = MPI_DOUBLE_PRECISION 
         MPI_SPLINE_COMPLEX = MPI_DOUBLE_COMPLEX
    endif
    if(CUSTOM_REAL.eq.4)then
         MPI_CUSTOM_REAL = MPI_REAL
    elseif(CUSTOM_REAL.eq.8)then
         MPI_CUSTOM_REAL = MPI_DOUBLE_PRECISION 
    endif 


    ! ASSUMING 1 mesh per proc
    sets_per_process = nprocs/cluster_size
    myset_start      = myrank*sets_per_process
    myset_end        = myset_start + sets_per_process - 1 

    ! Check equal load balance across the processes:
    if(all_warnings)then 
        if( mod(nprocs,cluster_size).ne.0)then 
            write(*,*)'Error: you are using '
            write(*,*)'     -- nprocs ',nprocs 
            write(*,*)'     -- nnodes ',cluster_size 
            write(*,*)'And therefore nnodes is not divisible by nnodes. Stop.'
            stop 
        else 
            if(myrank.eq.0 .and.verbose.ge.1)write(*,*)'Sets for each node:', sets_per_process
            print *, 'Process: ', myrank, 'does sets', myset_start, 'to ', myset_end
        endif 
    endif 
    
    if(force_VTI)then 
        model_ti = '_VTI'
    else 
        model_ti = ''
    endif 

    ! Determine if each task is doing more than one set
    IIN  = myrank
    IOUT = IIN + 2000


    ! Setup normal mode 1D model (may possibly be able to remove with some
    ! edits to the source code)

    call mineos%load_mineos_radial_info_MPI()
    mineos_ptr => mineos


    ! Load mesh data for this proc (1 set per proc)
    ! True false indicates load from disc and dont save to disc
    sm   = create_SetMesh(myset_start, region)
    call sm%setup_mesh_sem_details(.true., .false.)


    ! Compute elastic tensor in cartesian
    ! NOTE - we could compute rotation matrix once and for all and 
    ! load from disc instead of computing each time 
    call sm%compute_rotation_matrix()


    allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))


    ! Allocate Vani matrices
    allocate(Vani(max_tl1, max_tl1))
    if(myrank.eq.0)then 
        allocate(Vani_modesum(max_tl1, max_tl1))
        allocate(Vani_real(max_tl1, max_tl1))
    endif 


    ! (1) READ IN VORONOI MODEL, LOAD GLL MESH AND PROJECT 
    ! Read Hen's model and build K-d tree: 
    Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
    call Model3D%read_model_from_file()
    ! We can probably just do this for the first one since these should be the same each time
    call Model3D%create_KDtree()


    ! Loop through the model iterations: 
    do imodeliter = 38, nmodeliter

        call buffer_int(iterstr, imodeliter)
        Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS//voronoi/colat_adjusted_instances/instance_"//trim(iterstr)
        call Model3D%re_readmodel()


        vor_A = Model3D%valconsts(1)/100.0d0
        vor_C = Model3D%valconsts(2)/100.0d0
        vor_L = Model3D%valconsts(3)/100.0d0
        vor_N = Model3D%valconsts(4)/100.0d0
        vor_F = Model3D%valconsts(5)/100.0d0



        ! Project Voronoi model --> GLL grid 
        ! This allocates 2 angles for each global GLL point and finds the
        ! nearest neighbour values using the KD tree 
        ! Use a more efficient KD tree? 
    
        call Model3D%project_to_gll(sm, glob_eta1, id=1)
        call Model3D%project_to_gll(sm, glob_eta2, id=2)

        if(force_VTI)then 
            ! for benchmark - will delete for real runs 
            glob_eta1 = zero 
            glob_eta2 = zero
        endif 


        ! Convert to CUDA kernel? 
        call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, vor_N, & 
                                            vor_F, glob_eta1, glob_eta2, & 
                                            perturbation_on_prem=.true.)


        ! Loop for each mode to compute the splitting and the Cst value
        do i_mode = 1, nmodes
            n1       = modeNs(i_mode)
            l1       = modeLs(i_mode)
            this_tl1 = 2*l1 +1

            Vani         = SPLINE_iZERO
            if(myrank.eq.0)then
                Vani_modesum = SPLINE_iZERO
            endif
            
            

            ! Compute the Vani matrix
            call cuda_Vani_matrix_stored_selfcoupling(sm, n1, t1, l1)

            !Reduce the matrices across all of the MPI procs
            call MPI_Reduce(Vani, Vani_modesum, max_tl1**2, MPI_SPLINE_COMPLEX, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            ! USE Vani_modesum to output/compute CSTs
            if(myrank.eq.0)then 

                call buffer_int(nstr, n1)
                call buffer_int(lstr, l1)
                if(force_VTI)then 
                    out_name =  './output/sem_fast_'//trim(nstr)// t1//trim(lstr)//'_VTI.txt'
                else 
                    out_name =  './output/sem_fast_'//trim(nstr)// t1//trim(lstr)//'.txt'
                endif 
                Vani = Vani_modesum
                !call save_Vani_matrix(l1, out_name)
    
                call convert_imag_to_real(l1, l1, Vani_modesum(1:this_tl1, 1:this_tl1), Vani_real(1:this_tl1, 1:this_tl1))
                call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
                allocate(cst(num_s, ncols))
                call H_to_cst(Vani_real(1:this_tl1, 1:this_tl1), l1, l1, cst, ncols, num_s, t1, t1, 2)
                out_name = 'output/instance_csts/cst_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'_'//trim(iterstr)
                call write_cst_to_file(out_name, cst, ncols, num_s, smin, 2)
                deallocate(cst)
            endif

            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        enddo ! i_mode 
        if(myrank.eq.0)write(*,*)'Completed iteration: ', imodeliter
    enddo ! model iteration



    ! Cleanup memory 
    deallocate(Vani)
    if(myrank.eq.0) deallocate(Vani_modesum, Vani_real)

    call mpi_finalize(ierr)
end program optimised_vani

