program kdTree_test

    ! KD TREE stuff from coretran 
    use variableKind, only: i32, r64
    use m_allocate, only: allocate
    use m_deallocate, only: deallocate
    use m_random, only: rngNormal
    use m_KdTree, only: KdTree, KdTreeSearch
    use dArgDynamicArray_Class, only: dArgDynamicArray

    use params, only: IIN, IOUT, x_glob, y_glob, z_glob, xstore, ystore, & 
                      zstore, myrank, MPI_SPLINE_REAL, MPI_CUSTOM_REAL,  & 
                      MPI_SPLINE_COMPLEX, verbose, ONLY_ONE_TASK_PER_SET,& 
                      nglob, nprocs, start_clock, end_clock, count_rate, & 
                      elapsed_time
    use allocation_module, only: deallocate_if_allocated, & 
                                 allocate_if_unallocated
    use mesh_utils, only: read_proc_coordinates, cleanup_for_mode , & 
                          map_local_global, load_ibool
    
    implicit none 

include "constants.h"

#ifdef WITH_MPI
include 'mpif.h'
#endif 

    real(r64), allocatable :: x(:), y(:), z(:), D(:,:), v_val(:)
    integer(i32) :: ia, N, ind
    type(KdTree) :: tree
    type(KdTreeSearch) :: search
    type(dArgDynamicArray) :: da
    integer :: i, iset
    character(len=250) :: f_voroni
   
    real(kind=8), allocatable :: proc_id(:)
    integer :: iproc, region, sets_per_process, cluster_size, & 
               myset_start, myset_end, ierr


! SETUP THE MPI 
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


    call system_clock(count_rate=count_rate)
    call system_clock(start_clock)


    ! Read in voronoi x: 
    f_voroni = '/scratch/gpfs/we3822/NMSplit90/specfem_mesh/DATABASES_MPI/voronoi/voronoi_x.txt'
    open(1, file=f_voroni, form='formatted')
    read(1, *)N
    call allocate(x, N)
    call allocate(y, N)
    call allocate(z, N)
    call allocate(v_val, N)
    do i = 1, N
        read(1,*)  x(i), y(i), z(i), v_val(i)
    enddo 
    close(1)

    ! Build the tree
    !write(*,*)'Building tree...'
    tree = KdTree(x, y, z)
    !write(*,*)'Done'

    region = 3
    ! Read mineos model 
    call process_mineos_model(.false.)

    do iset = myset_start, myset_end

        ! Read the mesh info and coordinates
        call read_proc_coordinates(iset, region)
        allocate(proc_id(nglob))

        call load_ibool(iset, region)

        ! needed for ensight geo file
        call allocate_if_unallocated(nglob, x_glob)
        call allocate_if_unallocated(nglob, y_glob)
        call allocate_if_unallocated(nglob, z_glob)
        call map_local_global(xstore, x_glob, 0)
        call map_local_global(ystore, y_glob, 0)
        call map_local_global(zstore, z_glob, 0)

        ! For each global coordinate we want to map its voroni cell variable value 
        do i = 1, nglob
            da = search%kNearest(tree, x, y, z, & 
                                 xQuery = x_glob(i), &
                                 yQuery = y_glob(i), &
                                 zQuery = z_glob(i), & 
                                 k = 1)
            ind =  da%i%values(1)
            proc_id(i) = v_val(ind)
        enddo 

        call create_ensight_file_prefix(iset, region)
        call create_proc_case_file()
        call create_proc_geo_file(1)
        call write_real_scalar_to_ensight(proc_id, 'proc_id', 1)

        call cleanup_for_mode()

        deallocate(proc_id)
    enddo 



    ! Compute run time
    call system_clock(end_clock)
    ! Calculate the elapsed time in seconds
    elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
    ! Print the elapsed time
    write(*,*) 'Wall clock :', elapsed_time, 'seconds'



    
    !Deallocate any tree memory
    call tree%deallocate()
    call deallocate(x)
    call deallocate(y)
    call deallocate(z)
    call deallocate(D)
    end program