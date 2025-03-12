
program compute_vani_cross
    use params, only: Vani, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                      MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                      glob_eta2,  nmodes, nprocs, cluster_size, Arad, Crad, Lrad, Nrad, Frad
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                     compute_Vani_matrix, compute_vani_matrix_stored, & 
                     compute_Cxyz_at_gll_radialACLNF, compute_Cxyz_at_gll_generalVTI
use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst_8, write_cst_complex_to_file

#ifdef WITH_CUDA
    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
#endif

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

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2, n2, region, ierr, & 
               tl1, tl2, h, b, sets_per_process, & 
               myset_start, myset_end, i_mode, maxknot, smin, smax, ncols, num_s
    character ::  t1, t2
    character(len=2) n1str, l1str, n2str, l2str
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
    logical, parameter :: force_VTI      = .false.
    logical, parameter :: tromp93_model  = .false.
    logical, parameter :: benchmark_deuss  = .true.


    ! Added: 
    integer, dimension(1), parameter :: modeN1s = (/16/)
    integer, dimension(1), parameter :: modeL1s = (/5/)

    integer, dimension(1), parameter :: modeN2s = (/17/)
    integer, dimension(1), parameter :: modeL2s = (/4/)

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
    call mineos%load_mineos_radial_info_MPI()
#else
    ! Read mineos model 
    call mineos%process_mineos_model(.true.) 
#endif
mineos_ptr => mineos



if(tromp93_model)then 
    ! Read TROMP ACLNF model with 33 points (mineos for IC)
    call load_ACLNF_from_files('/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/tromp93/ACLNF', 33)
else

    

    ! Cross couple benchmark
    if(benchmark_deuss)then 
        Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/benchmarks/DR_benchmark_model_alt2.txt"
    else 
        ! Read Hen's model and build K-d tree: 
        Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
    endif 

    call Model3D%read_model_from_file()
    call Model3D%create_KDtree()
endif

! Benchmark value
!vor_A =  0.4d0
!vor_C = -0.2d0
!vor_L =  0.3d0
!vor_N = -0.5d0
!vor_F =  0.1d0

! Model values are a % perturbation on PREM so need to divide by 100 
! to get actual value
if(.not.tromp93_model .and. .not.benchmark_deuss)then
    vor_A = Model3D%valconsts(1)/100.0d0
    vor_C = Model3D%valconsts(2)/100.0d0
    vor_L = Model3D%valconsts(3)/100.0d0
    vor_N = Model3D%valconsts(4)/100.0d0
    vor_F = Model3D%valconsts(5)/100.0d0
endif 



if(ONLY_ONE_TASK_PER_SET)then 
    iset = myset_start
    sm = create_SetMesh(iset, region)

    call sm%setup_mesh_sem_details(load_from_bin, save_to_bin)

    if(.not.tromp93_model .and. .not.benchmark_deuss)then 
        allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
        !call project_voroni_to_gll(sm, tree)
        call Model3D%project_to_gll(sm, glob_eta1, id=1)
        call Model3D%project_to_gll(sm, glob_eta2, id=2)
        if(force_VTI)then 
            glob_eta1 = zero 
            glob_eta2 = zero
        endif 
    endif 

    ! Benchmark:
    if(benchmark_deuss)then 
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
    endif 


    call sm%compute_rotation_matrix()


    if(tromp93_model)then 
        call compute_Cxyz_at_gll_radialACLNF(sm, sm%interp%n_radial, &
                                             Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
    endif 

    if(benchmark_deuss)then 
        write(*,*)'DEUSS BENCHMARK...'
        ! VTI model
        call compute_Cxyz_at_gll_generalVTI(sm, Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
    else   
        call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, vor_N, & 
                                            vor_F, glob_eta1, glob_eta2, &
                                            perturbation_on_prem=.true.)
    endif 
endif 



do i_mode = 1, nmodes
    n1      =  modeN1s(i_mode)
    t1      = 'S'
    l1      =  modeL1s(i_mode)

    n2      =  modeN2s(i_mode)
    t2      = 'S'
    l2      =  modeL2s(i_mode)


    mode_1  = get_mode(n1, t1, l1, mineos_ptr)
    mode_2  = get_mode(n2, t2, l2, mineos_ptr)


    allocate(Vani(mode_1%tl1, mode_2%tl1))
    Vani = SPLINE_iZERO

    do iset = myset_start, myset_end
        if(.not.ONLY_ONE_TASK_PER_SET)then
            
            sm = create_SetMesh(iset, region)
            call sm%setup_mesh_sem_details(load_from_bin, save_to_bin)

            if(.not.tromp93_model .and. .not.benchmark_deuss)then
                allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
                !call project_voroni_to_gll(sm, tree)
                call Model3D%project_to_gll(sm, glob_eta1, id=1)
                call Model3D%project_to_gll(sm, glob_eta2, id=2)
                if(force_VTI)then 
                    glob_eta1 = zero 
                    glob_eta2 = zero
                endif 
            endif 

            ! Benchmark:
            if(benchmark_deuss)then 
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
            endif 

            call sm%compute_rotation_matrix()

            if(tromp93_model)then
                if(load_from_bin)call sm%get_unique_radii(save_to_bin)

                ! interpolate the ACLNF to the SM radii: 
                allocate(Aspl(sm%interp%n_radial))
                allocate(Cspl(sm%interp%n_radial))
                allocate(Lspl(sm%interp%n_radial))
                allocate(Nspl(sm%interp%n_radial))
                allocate(Fspl(sm%interp%n_radial))
                call sm%interp%interpolate_mineos_variable(real(Arad, kind=SPLINE_REAL), Aspl)
                call sm%interp%interpolate_mineos_variable(real(Crad, kind=SPLINE_REAL), Cspl)
                call sm%interp%interpolate_mineos_variable(real(Lrad, kind=SPLINE_REAL), Lspl)
                call sm%interp%interpolate_mineos_variable(real(Nrad, kind=SPLINE_REAL), Nspl)
                call sm%interp%interpolate_mineos_variable(real(Frad, kind=SPLINE_REAL), Fspl)
            

                allocate(A0(sm%interp%n_radial))
                allocate(vpspl(sm%interp%n_radial))
                allocate(rhospl(sm%interp%n_radial))

                call sm%interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rhospl)
                call sm%interp%interpolate_mineos_variable(real(mineos%vp_mineos,  kind=SPLINE_REAL), vpspl)
            
                ! Multiply by A0: 
                ! Note that rhospl and vpspl are already non-dimensionalised
                ! So i dont think we need to then non-dimensionalise the Aspl
                
                A0 = rhospl *  vpspl * vpspl 
                Aspl = Aspl * A0
                Cspl = Cspl * A0
                Lspl = Lspl * A0
                Nspl = Nspl * A0
                Fspl = Fspl * A0

                call compute_Cxyz_at_gll_radialACLNF(sm, sm%interp%n_radial, & 
                                                     Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
            elseif(benchmark_deuss) then
                call compute_Cxyz_at_gll_generalVTI(sm, Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
            else 
                call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, & 
                                                       vor_N, vor_F, glob_eta1, glob_eta2, &
                                                       perturbation_on_prem=.true.)
            endif 
        endif 



! Compute the Vani matrix
#ifdef WITH_CUDA
        write(*,*)'CUDA cross coupling not implemented. stop'
        stop
#else
        if(load_from_bin)then 
            call compute_Vani_matrix_stored(sm, t1, l1, n1, t2, l2, n2) 
        else 
            call compute_Vani_matrix(sm, n1, t1, l1, n2, t2, l2, .true.)
            write(*,*)'Done iset', iset
        endif
#endif

        if(.not.ONLY_ONE_TASK_PER_SET)then 
            if(tromp93_model)then 
                deallocate(Aspl, Cspl, Lspl, Nspl, Fspl, rhospl, vpspl, A0)
            elseif(benchmark_deuss)then
                deallocate(Aspl, Cspl, Lspl, Nspl, Fspl)
            else 
                deallocate(glob_eta1, glob_eta2)
            endif 
            call sm%cleanup()
        endif

        
    enddo !iset 


! ---------------------- OUTPUT THE V MATRIX FOR A MODE ----------------------
#ifdef WITH_MPI
    if(myrank.eq.0)then 
        allocate(Vani_modesum(mode_1%tl1, mode_2%tl2))
        Vani_modesum = SPLINE_iZERO  
    endif 

    call MPI_Reduce(Vani, Vani_modesum, mode_1%tl1*mode_2%tl2 , MPI_SPLINE_COMPLEX, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if(myrank.eq.0)then 
        call buffer_int(n1str, n1)
        call buffer_int(n2str, n2)
        call buffer_int(l1str, l1)
        call buffer_int(l2str, l2)
        if(force_VTI)then 
            out_name =  './output/vani'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)//'_VTI.txt'
        else 
            out_name =  './output/vani'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)//'.txt'
        endif 
        Vani = Vani_modesum
        call save_Vani_matrix(l1, l2, out_name)
        deallocate(Vani_modesum)
    endif 

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
    call buffer_int(n1str, n1)
    call buffer_int(n2str, n2)
    call buffer_int(l1str, l1)
    call buffer_int(l2str, l2)
    if(force_VTI)then 
        out_name =  './output/vani'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)//'_VTI.txt'
    else 
        out_name =  './output/vani'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)//'.txt'
    endif 
    
    call save_Vani_matrix(l1,l2, out_name)

    ! Write as a CST
    call get_Ssum_bounds(l1, l2, smin, smax, num_s, ncols)
    allocate(cst(num_s, ncols))
    call Hcomplex_to_cst_8(Vani, l1, l2, cst, ncols, num_s, t1, t2, 2)
    out_name = 'output/cst_'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)
    call write_cst_complex_to_file(out_name, cst, ncols, num_s, smin, 2)
    deallocate(cst)


#endif
    deallocate(Vani)
! ----------------- END OF OUTPUT THE V MATRIX FOR A MODE --------------

enddo ! i_mode 




#ifdef WITH_MPI
    call mpi_finalize(ierr)
#endif
end program compute_vani_cross

