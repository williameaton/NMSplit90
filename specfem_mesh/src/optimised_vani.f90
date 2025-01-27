
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

    implicit none
    include "constants.h"
    include 'mpif.h'

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2, n2, ierr, & 
               tl2, h, b, cluster_size, sets_per_process, & 
               myset_start, myset_end, i_mode, smin, smax, num_s, ncols, this_tl1
    character(len=2) nstr, lstr
    character(len=250) :: out_name

    real(kind=SPLINE_REAL), allocatable :: Vani_real(:,:)
    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)
    real(kind=SPLINE_REAL), allocatable :: cst(:,:)
    character(len=20) :: model_ti

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  

    ! Modes: 
    integer, dimension(27), parameter :: modeNs = (/ 6, 2, 3, 5, 7, 7, 8, 9, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 16, 18, 18, 20, 21, 21, 21, 25, 27/)
    integer, dimension(27), parameter :: modeLs = (/10, 3, 2, 3, 4, 5, 5, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  7,  3,  4,  1,  6,  8,  7,  2,  2/)

    ! Simulation parameters: 
    integer, parameter    :: region  = 3
    character, parameter  :: t1      = 'S'
    logical, parameter :: force_VTI = .true.

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



    ! (1) READ IN VORONOI MODEL, LOAD GLL MESH AND PROJECT 
    ! Read in voronoi model and build K-d tree: 
    call load_voronoi_model()
    ! Benchmark values - delete for real runs 
    vor_A =  0.4d0
    vor_C = -0.2d0
    vor_L =  0.3d0
    vor_N = -0.5d0
    vor_F =  0.1d0

    ! Setup KD tree using x,y,z of voronoi nodes
    tree = KdTree(vor_x, vor_y, vor_z) 

    ! Setup normal mode 1D model (may possibly be able to remove with some
    ! edits to the source code)

    call mineos%load_mineos_radial_info_MPI()
    mineos_ptr => mineos

    ! Load mesh data for this proc (1 set per proc)
    ! True false indicates load from disc and dont save to disc

    sm   = create_SetMesh(myset_start, region)
    call sm%setup_mesh_sem_details(.true., .false.)

    ! Project Voronoi model --> GLL grid 
    ! This allocates 2 angles for each global GLL point and finds the
    ! nearest neighbour values using the KD tree 
    ! Use a more efficient KD tree? 
    allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
    call project_voroni_to_gll(sm, tree)
 
    ! Compute elastic tensor in cartesian
    ! NOTE - we could compute rotation matrix once and for all and 
    ! load from disc instead of computing each time 
    call sm%compute_rotation_matrix()
    if(force_VTI)then 
        ! for benchmark - will delete for real runs 
        glob_eta1 = zero 
        glob_eta2 = zero
    endif 
    ! Convert to CUDA kernel? 

    call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, vor_N, & 
                                           vor_F, glob_eta1, glob_eta2)

    ! Allocate Vani matrices
    allocate(Vani(max_tl1, max_tl1))
    if(myrank.eq.0)then 
        allocate(Vani_modesum(max_tl1, max_tl1))
        allocate(Vani_real(max_tl1, max_tl1))
    endif 

    ! Loop for each mode to compute the splitting and the Cst value
    do i_mode = 1, 1!nmodes
        n1       =  6  !modeNs(i_mode)
        l1       =  10 !modeLs(i_mode)
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
            call save_Vani_matrix(l1, out_name)
   
            call convert_imag_to_real(l1, l1, Vani_modesum(1:this_tl1, 1:this_tl1), Vani_real(1:this_tl1, 1:this_tl1))
            call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
            allocate(cst(num_s, ncols))
            call H_to_cst(Vani_real(1:this_tl1, 1:this_tl1), l1, l1, cst, ncols, num_s, t1, t1, 2)
            out_name = 'output/cst_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'.txt'
            call write_cst_to_file(out_name, cst, ncols, num_s, smin, 2)
            deallocate(cst)
        endif


        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    enddo ! i_mode 



    ! Cleanup memory 
    deallocate(Vani)
    if(myrank.eq.0) deallocate(Vani_modesum, Vani_real)

    call mpi_finalize(ierr)
end program optimised_vani

