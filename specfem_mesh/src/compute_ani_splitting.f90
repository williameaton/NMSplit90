
program compute_vani_splitting
    use params, only: Vani, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                      MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                      glob_eta2,  nmodes, nprocs
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                     compute_Vani_matrix, compute_vani_matrix_stored
#ifdef WITH_CUDA
    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
#endif

    use m_KdTree, only: KdTree, KdTreeSearch
    use voronoi, only: vor_x, vor_y, vor_z, & 
                       vor_A, vor_C, vor_L, vor_N, vor_F, &
                       load_voronoi_model, project_voroni_to_gll
    use specfem_mesh, only: SetMesh, create_SetMesh
    use modes, only: get_mode, Mode 
    use mineos_model, only: mineos, mineos_ptr
    implicit none
    include "constants.h"

#ifdef WITH_MPI
    include 'mpif.h'
#endif 

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2, n2, region, ierr, & 
               tl1, tl2, h, b, cluster_size, sets_per_process, & 
               myset_start, myset_end, i_mode, maxknot
    character ::  t1
    character(len=2) nstr, lstr
    character(len=12) nprocstr, nmodestr, timing_fmt_vals
    character(len=250) :: out_name

    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  
    type(Mode)             :: mode_1 

    ! Switches 
    logical :: ONLY_ONE_TASK_PER_SET
    logical, parameter :: load_from_bin = .false.
    logical, parameter :: save_to_bin   = .true.
    logical, parameter :: force_VTI     = .true.

    ! Modes: 
    integer, dimension(28), parameter :: modeNs = (/5, 6, 7, 8, 21, 7, 9, 2, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21, 25, 27, 21, 21, 16/)
    integer, dimension(28), parameter :: modeLs = (/3, 3, 4, 5,  7, 5, 2, 3, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,  2,  2,  8,  6,  7/)

 
#ifdef WITH_MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, cluster_size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    ! Setup MPI precisions: 
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

    ! Check equal load balance across the processes:
    sets_per_process =  nprocs/cluster_size
    if( mod(nprocs,cluster_size).ne.0)then 
        write(*,*)'Error: you are using '
        write(*,*)'     -- nprocs ',nprocs 
        write(*,*)'     -- nnodes ',cluster_size 
        write(*,*)'And therefore nnodes is not divisible by nnodes. Stop.'
        stop 
    else 
        if(myrank.eq.0 .and.verbose.ge.1)write(*,*)'Sets for each node:', sets_per_process
        myset_start = myrank*sets_per_process
        myset_end = myset_start + sets_per_process - 1 

        print *, 'Process: ', myrank, 'does sets', myset_start, 'to ', myset_end

    endif 
    
    ! Determine if each task is doing more than one set
    ONLY_ONE_TASK_PER_SET = (myset_start.eq.myset_end)

    IIN = myrank
    IOUT = IIN + 2000
#else 
    write(*,*)"Computing without OpenMPI"
    myset_start = 0
    myset_end   = nprocs-1

    IIN  = 1
    IOUT = 101
    ONLY_ONE_TASK_PER_SET = .false.
#endif

region = 3


! Read in voronoi model and build K-d tree: 
call load_voronoi_model()

! Benchmark value
vor_A =  0.4d0
vor_C = -0.2d0
vor_L =  0.3d0
vor_N = -0.5d0
vor_F =  0.1d0


tree = KdTree(vor_x, vor_y, vor_z) 


#ifdef WITH_MPI
    call mineos%load_mineos_radial_info_MPI()
#else
    ! Read mineos model 
    call mineos%process_mineos_model(.true.) 
#endif
mineos_ptr => mineos



if(ONLY_ONE_TASK_PER_SET)then 
    iset = myset_start
    sm = create_SetMesh(iset, region)


    call sm%setup_mesh_sem_details(load_from_bin, save_to_bin)



    allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
    call project_voroni_to_gll(sm, tree)


    call sm%compute_rotation_matrix()
    if(force_VTI)then 
        glob_eta1 = zero 
        glob_eta2 = zero
    endif 
    call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, vor_N, & 
                                           vor_F, glob_eta1, glob_eta2)
endif 



do i_mode = 1, nmodes
    n1      =  modeNs(i_mode)
    t1      = 'S'
    l1      =  modeLs(i_mode)

    mode_1  = get_mode(n1, t1, l1, mineos_ptr)


    allocate(Vani(mode_1%tl1, mode_1%tl1))
    Vani = SPLINE_iZERO


    do iset = myset_start, myset_end
        if(.not.ONLY_ONE_TASK_PER_SET)then
            
            sm = create_SetMesh(iset, region)
            call sm%setup_mesh_sem_details(load_from_bin, save_to_bin)

            allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
            call project_voroni_to_gll(sm, tree)
        
            call sm%compute_rotation_matrix()

            if(force_VTI)then 
                glob_eta1 = zero 
                glob_eta2 = zero
            endif 
            call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, & 
                                                   vor_N, vor_F, glob_eta1, glob_eta2)
        endif 




! Compute the Vani matrix
#ifdef WITH_CUDA
        call cuda_Vani_matrix_stored_selfcoupling(sm, n1, t1, l1)
#else
        if(load_from_bin)then 
            call compute_Vani_matrix_stored(sm, t1, l1, n1, t1, l1, n1) 
        else 
            call compute_Vani_matrix(sm, n1, t1, l1, n1, t1, l1, .true.)
            write(*,*)'Done iset', iset
        endif
#endif

        if(.not.ONLY_ONE_TASK_PER_SET)then 
            deallocate(glob_eta1, glob_eta2)
            call sm%cleanup()
        endif

        
    enddo !iset 


! ---------------------- OUTPUT THE V MATRIX FOR A MODE ----------------------
#ifdef WITH_MPI
    if(myrank.eq.0)then 
        allocate(Vani_modesum(mode_1%tl1, mode_1%tl1))
        Vani_modesum = SPLINE_iZERO  
    endif 

    call MPI_Reduce(Vani, Vani_modesum, mode_1%tl1**2, MPI_SPLINE_COMPLEX, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

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
        deallocate(Vani_modesum)
    endif 

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
    call buffer_int(nstr, n1)
    call buffer_int(lstr, l1)
    if(force_VTI)then 
        out_name =  './output/sem_fast_'//trim(nstr)// t1//trim(lstr)//'_VTI.txt'
    else 
        out_name =  './output/sem_fast_'//trim(nstr)// t1//trim(lstr)//'.txt'
    endif 
    
    call save_Vani_matrix(l1, out_name)
#endif
    deallocate(Vani)
! ----------------- END OF OUTPUT THE V MATRIX FOR A MODE --------------

enddo ! i_mode 




#ifdef WITH_MPI
    call mpi_finalize(ierr)
#endif
end program compute_vani_splitting

