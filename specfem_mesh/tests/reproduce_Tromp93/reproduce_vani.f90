
program reproduce_tromp93
    use params, only: Vani, verbose,  MPI_SPLINE_COMPLEX, & 
                      MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                      glob_eta2,  nmodes, nprocs, cluster_size, Arad, Crad, Lrad, Nrad, Frad
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                     compute_Vani_matrix, compute_vani_matrix_stored, & 
                     compute_Cxyz_at_gll_radialACLNF, compute_Cxyz_at_gll_generalVTI

    use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst_8, write_cst_complex_to_file

    use m_KdTree, only: KdTree, KdTreeSearch
    use voronoi, only: vor_x, vor_y, vor_z, & 
                       vor_A, vor_C, vor_L, vor_N, vor_F, &
                       load_voronoi_model, project_voroni_to_gll

    use model3d, only: M3D
    

    use specfem_mesh, only: SetMesh, create_SetMesh
    use modes, only: get_mode, Mode 
    use mineos_model, only: mineos, mineos_ptr
    implicit none
    include "constants.h"

#ifdef WITH_MPI
    include 'mpif.h'
#endif 

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2,  region, ierr, & 
               tl1, tl2, h, b, sets_per_process, & 
               myset_start, myset_end, i_mode, maxknot, smin, smax, ncols, num_s, myrank
    character ::  t1
    character(len=2) n1str, l1str
    character(len=12) nprocstr, nmodestr, timing_fmt_vals
    character(len=250) :: out_name
    complex(kind=SPLINE_REAL), allocatable :: cst(:,:)

    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)

    real(kind=CUSTOM_REAL), allocatable    :: Aspl(:), Cspl(:), Lspl(:), & 
                                              Nspl(:), Fspl(:)
    real(kind=CUSTOM_REAL), allocatable    :: A0(:), vpspl(:), rhospl(:)

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  
    type(Mode)             :: mode_1 , mode_2

    ! 3D model
    type(M3D) :: model3D

    ! Switches 
    logical :: ONLY_ONE_TASK_PER_SET
    logical, parameter :: load_from_bin  = .false.
    logical, parameter :: save_to_bin    = .true.

    ! Added: 
    integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 6, 8, 9, 11, 11, 13, 13, 14, 15, 16, 18, 20, 21, 23, 25, 27/)
    integer, dimension(nmodes), parameter :: modeLs = (/3, 2, 3, 5, 3, 4, 5, 2, 3, 4, 3, 6, 4, 5, 6, 5, 2, 2/)



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


#ifdef WITH_MPI
    call mineos%load_mineos_radial_info_MPI(MPI_COMM_WORLD)
#else
    ! Read mineos model 
    call mineos%process_mineos_model(.true.) 
#endif
mineos_ptr => mineos



Model3D%filename = "/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/3D_MODELS/benchmarks/Tromp93_benchmark_model.txt"


call Model3D%read_model_from_file()
call Model3D%create_KDtree()



if(ONLY_ONE_TASK_PER_SET)then 
    iset = myset_start
    sm = create_SetMesh(iset, region)

    call sm%setup_mesh_sem_details(load_from_bin, save_to_bin)


    ! These are not splines but we can use the arrays: 
    allocate(Aspl(sm%nglob))
    allocate(Cspl(sm%nglob))
    allocate(Lspl(sm%nglob))
    allocate(Nspl(sm%nglob))
    allocate(Fspl(sm%nglob))

    call Model3D%project_to_gll(sm, Aspl, id=1) ! A 
    call Model3D%project_to_gll(sm, Cspl, id=2) ! C
    call Model3D%project_to_gll(sm, Lspl, id=3) ! L
    call Model3D%project_to_gll(sm, Nspl, id=4) ! N
    call Model3D%project_to_gll(sm, Fspl, id=5) ! N

    call sm%compute_rotation_matrix()

    call compute_Cxyz_at_gll_generalVTI(sm, Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
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

            ! These are not splines but we can use the arrays: 
            allocate(Aspl(sm%nglob))
            allocate(Cspl(sm%nglob))
            allocate(Lspl(sm%nglob))
            allocate(Nspl(sm%nglob))
            allocate(Fspl(sm%nglob))

            call Model3D%project_to_gll(sm, Aspl, id=1) ! A 
            call Model3D%project_to_gll(sm, Cspl, id=2) ! C
            call Model3D%project_to_gll(sm, Lspl, id=3) ! L
            call Model3D%project_to_gll(sm, Nspl, id=4) ! N
            call Model3D%project_to_gll(sm, Fspl, id=5) ! F

            call sm%compute_rotation_matrix()

            call compute_Cxyz_at_gll_generalVTI(sm, Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
        endif 



        ! Compute the Vani matrix
        if(load_from_bin)then 
            call compute_Vani_matrix_stored(sm, t1, l1, n1, t1, l1, n1) 
        else 
            call compute_Vani_matrix(sm, n1, t1, l1, n1, t1, l1, .true.)
            write(*,*)'Done iset', iset
        endif


        if(.not.ONLY_ONE_TASK_PER_SET)then 
            deallocate(Aspl, Cspl, Lspl, Nspl, Fspl)
            call sm%cleanup()
        endif

        
    enddo !iset 


! ---------------------- OUTPUT THE V MATRIX FOR A MODE ----------------------
#ifdef WITH_MPI
    if(myrank.eq.0)then 
        allocate(Vani_modesum(mode_1%tl1, mode_1%tl1))
        Vani_modesum = SPLINE_iZERO  
    endif 

    call MPI_Reduce(Vani, Vani_modesum, mode_1%tl1*mode_1%tl1 , MPI_SPLINE_COMPLEX, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if(myrank.eq.0)then 
        call buffer_int(n1str, n1)
        call buffer_int(l1str, l1)

      
        out_name =  './reproduce_Tromp93/matrices//vani'//trim(n1str)//t1//trim(l1str)//'.txt'

        Vani = Vani_modesum
        call save_Vani_matrix(l1, l1, out_name)
        deallocate(Vani_modesum)
    endif 

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
    call buffer_int(n1str, n1)
    call buffer_int(l1str, l1)

    out_name =  './reproduce_Tromp93/matrices//vani'//trim(n1str)//t1//trim(l1str)//'.txt'

    
    call save_Vani_matrix(l1,l1, out_name)

    ! Write as a CST
    !call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
    !allocate(cst(num_s, ncols))
    !call Hcomplex_to_cst_8(Vani, l1, l1, cst, ncols, num_s, t1, t2, 2)
    !out_name = 'output/cst_'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)
    !call write_cst_complex_to_file(out_name, cst, ncols, num_s, smin, 2)
    !deallocate(cst)


#endif
    deallocate(Vani)
! ----------------- END OF OUTPUT THE V MATRIX FOR A MODE --------------

enddo ! i_mode 




#ifdef WITH_MPI
    call mpi_finalize(ierr)
#endif
end program reproduce_tromp93

