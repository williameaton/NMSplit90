
program compute_vani_splitting
    use params, only: Vani, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                      MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                      glob_eta2,  nmodes, nprocs, cluster_size, Arad, Crad, Lrad, Nrad, Frad
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                     compute_Vani_matrix, compute_vani_matrix_stored, & 
                     compute_Cxyz_at_gll_radialACLNF, compute_Cxyz_at_gll

#ifdef WITH_CUDA
    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
#endif
use splitting_function, only: get_Ssum_bounds, Hreal_to_cst, write_cst_to_file, Hcomplex_to_cst_8, write_cst_complex_to_file

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
               myset_start, myset_end, i_mode, maxknot, ib
    character ::  t1
    character(len=2) nstr, lstr
    character(len=12) nprocstr, nmodestr, timing_fmt_vals
    character(len=250) :: out_name
    integer, allocatable :: updated(:)
    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)

    real(kind=CUSTOM_REAL), allocatable    :: Aspl(:), Cspl(:), Lspl(:), & 
                                              Nspl(:), Fspl(:)
    real(kind=CUSTOM_REAL), allocatable    :: A0(:), vpspl(:), rhospl(:)

    ! CSTS: 
    integer :: smin, smax, num_s, ncols
    complex(kind=SPLINE_REAL), allocatable :: cst_imag(:,:)
    complex(kind=SPLINE_REAL) :: trace 

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  
    type(Mode)             :: mode_1 

    ! 3D model
    type(M3D) :: model3D

    ! Switches 
    logical :: ONLY_ONE_TASK_PER_SET
    logical, parameter :: load_from_bin         = .false.
    logical, parameter :: save_to_bin           = .false.
    logical, parameter :: force_VTI             = .true.

    logical, parameter :: constant_ACLNF        = .true.
    logical, parameter :: use_radial_eta12      = .false.
    logical, parameter :: tromp93_model         = .false.
    logical, parameter :: perturbation_on_prem  = .true.


    logical, parameter :: compute_csts          = .true.
    logical, parameter :: redimensionalise_csts = .false.

    ! Modes: 
    !integer, dimension(29), parameter :: modeNs = (/2, 5, 6, 7, 8, 21, 7, 9, 2, 3, 9, 9, 11, 11, 13, 13, 13, 13, 15, 15, 18, 18, 20, 21, 25, 27, 21, 21, 16/)
    !integer, dimension(29), parameter :: modeLs = (/3, 3, 3, 4, 5,  7, 5, 2, 3, 2, 3, 4,  4,  5,  1,  2,  3,  6,  3,  4,  3,  4,  1,  6,  2,  2,  8,  6,  7/)
    !integer, dimension(40), parameter :: modeNs =  (/2, 3, 3, 5, 6, 8, 8, 9, 9, 9, 11, 11, 11, 11, 13, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21, 21, 22, 23, 23, 25, 25, 27/)
    !integer, dimension(40), parameter :: modeLs =  (/3, 1, 2, 2, 3, 1, 5, 2, 3, 4,  1,  4,  5,  6, 1,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7,  8,  1,  4,  5,  1,  2,  2/)

    ! Added: 
    ! REAL 30: 
    integer, dimension(1), parameter :: modeNs = (/2/)
    integer, dimension(1), parameter :: modeLs = (/3/)


    ! Deuss fig lists: 
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 6, 8, 8, 9, 11, 11, 13, 13, 13, 16, 16, 17, 18, 18, 21, 21, 23, 23/)
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 3, 1, 5, 3, 4,  5,  1,  2,  3,  5,  7,  1,  3,  4,  6,  7,  4,  5/)

    !integer, dimension(5), parameter :: modeNs = (/5, 8, 7, 3, 13/)
    !integer, dimension(5), parameter :: modeLs = (/3, 1, 5, 2, 2/)



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


if(tromp93_model)then 
    ! Read TROMP ACLNF model with 33 points (mineos for IC)
    call load_ACLNF_from_files('/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/3D_MODELS/tromp93/ACLNF', 33)
elseif(constant_ACLNF)then 
    Model3D%nconst = 5  
    allocate(Model3D%valconsts(Model3D%nconst))

    Model3D%valconsts(1) =  0.4d0
    Model3D%valconsts(2) = -0.2d0
    Model3D%valconsts(3) =  0.3d0
    Model3D%valconsts(4) = -0.5d0
    Model3D%valconsts(5) =  0.1d0
else
    if(myrank.eq.0)write(*,*)'Reading 3D model...'
    ! Read Hen's model and build K-d tree: 
      !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
      Model3D%filename = "/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/3D_MODELS/benchmarks/benchmark_TI_perturb_prem.txt"
    ! Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi//MCMC_models/instances/c1_m10900.txt"
    call Model3D%read_model_from_file()
    call Model3D%create_KDtree()
    if(myrank.eq.0)write(*,*)'Done.'
endif



if(ONLY_ONE_TASK_PER_SET)then 

    iset = myset_start
    sm = create_SetMesh(iset, region)

    call sm%setup_mesh_sem_details(load_from_bin, save_to_bin)

    if(tromp93_model)then 
        call sm%compute_rotation_matrix()
        call compute_Cxyz_at_gll_radialACLNF(sm, sm%interp%n_radial, &
                                             Aspl, Cspl, Lspl, Nspl, Fspl, zero, zero)
    else
        allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
        !call project_voroni_to_gll(sm, tree)
        ! call Model3D%project_to_gll(sm, glob_eta1, id=1)
        ! call Model3D%project_to_gll(sm, glob_eta2, id=2)


        allocate(Aspl(sm%nglob), Cspl(sm%nglob), Lspl(sm%nglob), Nspl(sm%nglob), Fspl(sm%nglob))

        if(constant_ACLNF)then 
            Aspl = Model3D%valconsts(1)
            Cspl = Model3D%valconsts(2)
            Lspl = Model3D%valconsts(3)
            Nspl = Model3D%valconsts(4)
            Fspl = Model3D%valconsts(5)
        else 
                    ! If the ACLNF vary at each point: 
                    ! we can use the Aspl arrays: 
        call Model3D%project_to_gll(sm, Aspl, id=1, scaling=one/(RHOAV*SCALE_V*SCALE_V))
        call Model3D%project_to_gll(sm, Cspl, id=2, scaling=one/(RHOAV*SCALE_V*SCALE_V))
        call Model3D%project_to_gll(sm, Lspl, id=3, scaling=one/(RHOAV*SCALE_V*SCALE_V))
        call Model3D%project_to_gll(sm, Nspl, id=4, scaling=one/(RHOAV*SCALE_V*SCALE_V))
        call Model3D%project_to_gll(sm, Fspl, id=5, scaling=one/(RHOAV*SCALE_V*SCALE_V))    
        endif 
    
        if(force_VTI)then 
            glob_eta1 = zero
            glob_eta2 = zero
        endif 
        
        call sm%compute_rotation_matrix()

        if(use_radial_eta12)then 
            call sm%compute_eta12_radial()
        endif 

        if(constant_ACLNF)then 
            call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, vor_N, & 
                                                vor_F, glob_eta1, glob_eta2, &
                                                perturbation_on_prem=perturbation_on_prem)
        else 
            call compute_Cxyz_at_gll(sm, Aspl, Cspl, Lspl, Nspl, Fspl, &
                            glob_eta1, glob_eta2, & 
                            perturbation_on_PREM)
        endif


    endif 
endif 



do i_mode = 1, 1
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

            if(.not.tromp93_model)then
                allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
                !call project_voroni_to_gll(sm, tree)
                !call Model3D%project_to_gll(sm, glob_eta1, id=1)
                !call Model3D%project_to_gll(sm, glob_eta2, id=2)

                ! If the ACLNF vary at each point: 
                ! we can use the Aspl arrays: 
                allocate(Aspl(sm%nglob), Cspl(sm%nglob), Lspl(sm%nglob), Nspl(sm%nglob), Fspl(sm%nglob))

                ! Constant values
                if(constant_ACLNF)then 
                    Aspl = Model3D%valconsts(1)
                    Cspl = Model3D%valconsts(2)
                    Lspl = Model3D%valconsts(3)
                    Nspl = Model3D%valconsts(4)
                    Fspl = Model3D%valconsts(5)   
                else
                    call Model3D%project_to_gll(sm, Aspl, id=1, scaling=one/(RHOAV*SCALE_V*SCALE_V))
                    call Model3D%project_to_gll(sm, Cspl, id=2, scaling=one/(RHOAV*SCALE_V*SCALE_V))
                    call Model3D%project_to_gll(sm, Lspl, id=3, scaling=one/(RHOAV*SCALE_V*SCALE_V))
                    call Model3D%project_to_gll(sm, Nspl, id=4, scaling=one/(RHOAV*SCALE_V*SCALE_V))
                    call Model3D%project_to_gll(sm, Fspl, id=5, scaling=one/(RHOAV*SCALE_V*SCALE_V))                     
                endif 

                if(force_VTI)then 
                    glob_eta1 = zero 
                    glob_eta2 = zero 
                endif 
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
            else 
                if(use_radial_eta12)then 
                    call sm%compute_eta12_radial()
                endif 

                call compute_Cxyz_at_gll(sm, Aspl, Cspl, Lspl, Nspl, Fspl, &
                                         glob_eta1, glob_eta2, & 
                                         perturbation_on_PREM)

                !call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, & 
                !                                    vor_N, vor_F, glob_eta1, glob_eta2, &
                !                                    perturbation_on_prem=perturbation_on_prem)
            endif 

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
            if(tromp93_model)then 
                deallocate(Aspl, Cspl, Lspl, Nspl, Fspl, rhospl, vpspl, A0)
            else
                deallocate(Aspl, Cspl, Lspl, Nspl, Fspl)
                deallocate(glob_eta1, glob_eta2)
            endif 
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
            out_name =  './output/sem_fast_'//trim(nstr)// t1//trim(lstr)//'.txt'
        else 
            out_name =  './output/sem_fast_'//trim(nstr)// t1//trim(lstr)//'.txt'
        endif 
        Vani = Vani_modesum
        call save_Vani_matrix(l1, l1, out_name)

        ! Output CSTS if desired
        if(compute_csts)then 
            call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
            allocate(cst_imag(num_s, ncols))
            call Hcomplex_to_cst_8(Vani, l1, l1, cst_imag, ncols, num_s, t1, t1, 2)
    
            out_name = 'output/csts/'//'/cst_'//trim(nstr)//trim(t1)//trim(lstr)//'.txt'
    
            ! COMPUTE THE TRACE AND SCALE IT ACCORDINGLY 
            ! trace = SPLINE_iZERO
            ! do i = 1, 2*l1+1
            !     trace = trace + Vani(i,i)
            ! enddo 

            ! dimensionalise: 
            ! note the scale_T*mode_1%wcom is non-dim omega
            ! then multiply by 1/SCALE_T for dimensionalisation and the 2pi is --> Hz
            if(redimensionalise_csts)then 
                cst_imag = cst_imag/(four*pi*SCALE_T*SCALE_T * mode_1%wcom)
            endif 

            call write_cst_complex_to_file(trim(out_name), cst_imag, ncols, num_s, smin, 2)
            write(*,*)'Cst written to '//trim(out_name)
            deallocate(cst_imag)
        endif 
        deallocate(Vani_modesum)
    endif 

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
    call buffer_int(nstr, n1)
    call buffer_int(lstr, l1)

    out_name =  './output/sem_fast_'//trim(nstr)//t1//trim(lstr)//'.txt'
    call save_Vani_matrix(l1, l1, out_name)


    ! Output CSTS if desired
    ! NOTE ON NON-DIMENSIONALISATION: 
    ! The Vani matrix is re-dimensionalised in save_Vani_matrix by multiplying
    ! by 1/SCALE_T^2 
    ! It does not include the multiplication of the 1/(2omega) part
    ! that would be done in computing Hmat...
    ! we include that here for the Cst calculation so that the output Cst is the 
    ! true value Cst contribution from this matrix
    ! also convert it to Hz, but not mHz
    if(compute_csts)then 
        call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
        allocate(cst_imag(num_s, ncols))
        call Hcomplex_to_cst_8(Vani, l1, l1, cst_imag, ncols, num_s, t1, t1, 2)

        out_name = 'output/csts/'//'/cst_'//trim(nstr)//trim(t1)//trim(lstr)//'.txt'

        ! dimensionalise: 
        ! note the scale_T*mode_1%wcom is non-dim omega
        ! then multiply by 1/SCALE_T for dimensionalisation and the 2pi is --> Hz
        if(redimensionalise_csts)then 
            cst_imag = cst_imag/(four*pi*SCALE_T*SCALE_T*mode_1%wcom)
        endif 

        call write_cst_complex_to_file(trim(out_name), cst_imag, ncols, num_s, smin, 2)
        write(*,*)'Cst written to '//trim(out_name)
        deallocate(cst_imag)
    endif 

#endif
    deallocate(Vani)
! ----------------- END OF OUTPUT THE V MATRIX FOR A MODE --------------

enddo ! i_mode 




#ifdef WITH_MPI
    call mpi_finalize(ierr)
#endif
end program compute_vani_splitting

