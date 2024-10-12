
program compute_vani_splitting
    use params, only: Vani, nspec, nprocs, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                      MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                      glob_eta2, nglob, nmodes, IC_ID
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: cleanup_for_mode, compute_rotation_matrix
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                     compute_Vani_matrix_stored_selfcoupling, compute_Vani_matrix
#ifdef WITH_CUDA
    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
#endif

    use m_KdTree, only: KdTree, KdTreeSearch
    use voronoi, only: vor_x, vor_y, vor_z, & 
                       vor_A, vor_C, vor_L, vor_N, vor_F, &
                       load_voronoi_model, project_voroni_to_gll

    implicit none
    include "constants.h"

#ifdef WITH_MPI
include 'mpif.h'
#endif 

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2, n2, region, ierr, & 
               tl1, tl2, h, b, cluster_size, sets_per_process, & 
               myset_start, myset_end, i_mode, maxknot
    character ::  type_1, type_2
    character(len=2) nstr, lstr
    character(len=12) nprocstr, nmodestr, timing_fmt_vals
    character(len=250) :: out_name
    real(kind=SPLINE_REAL) :: min_r, min_i, thirty, twone

    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)

    ! KD tree: 
    type(KdTree)           :: tree


    ! Switches 
    logical :: ONLY_ONE_TASK_PER_SET
    logical, parameter :: load_from_bin = .true.
    logical, parameter :: save_to_bin   = .false.
    logical, parameter :: force_VTI     = .true.

    ! Modes: 
    !integer, dimension(20), parameter :: modeNs = (/2, 3, 9, 9, 9, 11, 11, 13,13,13,13,15,15,18,18,20,21,25,27, 6/)
    !integer, dimension(20), parameter :: modeLs = (/3, 2, 2, 3, 4,  4,  5,  1, 2, 3, 6, 3, 4, 3, 4, 1, 6,2,2, 10/)
    
    integer, dimension(16), parameter :: modeNs = (/6, 6, 6, 6, 6,  6,  6,  6, 6, 6, 6, 6, 6, 6, 6, 6 /)
    integer, dimension(16), parameter :: modeLs = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)

    ! Timing: 
    integer :: start_clock, end_clock, count_rate, mid1, mid2
    real(8),allocatable :: elapsed_time(:)


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


call system_clock(count_rate=count_rate)


! Read in voronoi model and build K-d tree: 
call load_voronoi_model()
tree = KdTree(vor_x, vor_y, vor_z) 


#ifdef WITH_MPI
    call load_mineos_radial_info_MPI()
#else
    ! Read mineos model 
    call process_mineos_model(.true.) 
#endif

maxknot = IC_ID


if(ONLY_ONE_TASK_PER_SET)then 
    iset = myset_start
    call setup_mesh_sem_details(iset, region, load_from_bin, save_to_bin, maxknot)

    allocate(glob_eta1(nglob), glob_eta2(nglob))
    call project_voroni_to_gll(tree)


    call compute_rotation_matrix()
    if(force_VTI)then 
        glob_eta1 = zero 
        glob_eta2 = zero
    endif 
    call compute_Cxyz_at_gll_constantACLNF(vor_A, vor_C, vor_L, vor_N, & 
                                           vor_F, glob_eta1, glob_eta2)
endif 

allocate(elapsed_time(nmodes))


do i_mode = 1, nmodes

    ! Start clock count
    call system_clock(start_clock)
    n1      = modeNs(i_mode)
    type_1  = 'S'
    l1      = modeLs(i_mode)

    ! Setup Vani matrix
    tl1 = 2*l1 + 1
    allocate(Vani(tl1, tl1))
    Vani = SPLINE_iZERO




    do iset = myset_start, myset_end
        if(.not.ONLY_ONE_TASK_PER_SET)then
            call setup_mesh_sem_details(iset, region, load_from_bin, save_to_bin, maxknot)

            allocate(glob_eta1(nglob), glob_eta2(nglob))
            call project_voroni_to_gll(tree)
        
            call compute_rotation_matrix()

            if(force_VTI)then 
                glob_eta1 = zero 
                glob_eta2 = zero
            endif 
            call compute_Cxyz_at_gll_constantACLNF(vor_A, vor_C, vor_L, & 
                                                   vor_N, vor_F, glob_eta1, glob_eta2)
        endif 



! Compute the Vani matrix
#ifdef WITH_CUDA
        call cuda_Vani_matrix_stored_selfcoupling(type_1, l1, n1, iset)
#else
        if(load_from_bin)then 
            call compute_Vani_matrix_stored_selfcoupling(type_1, l1, n1, iset)
        else 
            call compute_Vani_matrix(type_1, l1, n1, & 
                                    type_1, l1, n1, & 
                                    .true., iset)
            write(*,*)'Done iset', iset
        endif
#endif

        if(.not.ONLY_ONE_TASK_PER_SET)then 
            deallocate(glob_eta1, glob_eta2)
            call cleanup_for_mode()
        endif
    enddo !iset 


! ---------------------- OUTPUT THE V MATRIX FOR A MODE ----------------------
#ifdef WITH_MPI
    if(myrank.eq.0)then 
        allocate(Vani_modesum(tl1, tl1))
        Vani_modesum = SPLINE_iZERO  
    endif 

    call MPI_Reduce(Vani, Vani_modesum, tl1*tl1, MPI_SPLINE_COMPLEX, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if(myrank.eq.0)then 
        call buffer_int(nstr, n1)
        call buffer_int(lstr, l1)
        out_name =  './output/sem_fast_'//trim(nstr)//type_1//trim(lstr)// '.txt'

        Vani = Vani_modesum
        call save_Vani_matrix(l1, out_name)
        ! Compute run time
        call system_clock(end_clock)
        elapsed_time(i_mode) = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
        deallocate(Vani_modesum)
    endif 

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
    call buffer_int(nstr, n1)
    call buffer_int(lstr, l1)
    out_name =  './output/sem_fast_'//trim(nstr)// type_1//trim(lstr)// '.txt'
    call save_Vani_matrix(l1, out_name)
    call system_clock(end_clock)
    elapsed_time(i_mode) = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
#endif
    deallocate(Vani)
! ----------------- END OF OUTPUT THE V MATRIX FOR A MODE --------------

enddo ! i_mode 



#ifdef WITH_MPI 
#ifdef WITH_CUDA
if(myrank.eq.0)then 

    call buffer_int(nprocstr, nprocs)
    call buffer_int(nmodestr, nmodes)

    open(1, file='output/timings/timing_176_'//trim(nprocstr)//'_'//trim(nmodestr), & 
            status='unknown', action='write', position='append' )

    ! Use the nproc to buffer the number of floats needed for the format of elapsedtime
    call buffer_int(nprocstr, nmodes)
    timing_fmt_vals = "(a,"//trim(nprocstr)//"f12.6)"

    ! Use the nmodestr to write the np 
    call buffer_int(nmodestr, cluster_size)
    write(1,timing_fmt_vals) nmodestr, elapsed_time
    close(1)
endif
#endif
#endif




#ifdef WITH_MPI
    call mpi_finalize(ierr)
#endif
end program compute_vani_splitting













subroutine setup_mesh_sem_details(iset, region, load_from_bin, save_to_bin, maxknot)
    use mesh_utils, only: read_proc_coordinates, load_ibool, & 
                          compute_jacobian, compute_rtp_from_xyz, & 
                          setup_global_coordinate_arrays
    use gll, only: setup_gll, compute_wglljac
    implicit none 
    integer :: iset, region, maxknot
    logical :: load_from_bin, save_to_bin

    ! Things that need to be done for each set  
    call read_proc_coordinates(iset, region)
    call load_ibool(iset, region)
    call setup_gll()


    if(load_from_bin)then 
        call load_jacobian(iset)
        call load_wglljac(iset)
        call load_global_xyz(iset)
        call load_elem_rtp(iset)
        call load_get_mesh_radii_results(iset)
    else
        call compute_jacobian(iset, save_to_bin)
        call compute_wglljac(iset, save_to_bin)
        call setup_global_coordinate_arrays(iset, save_to_bin)
        call compute_rtp_from_xyz(iset, save_to_bin)
        call get_mesh_radii(iset, save_to_bin, maxknot)


    endif 
end subroutine setup_mesh_sem_details