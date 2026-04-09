
program fairhead_optimised_vani
    ! This is a version of optimised vani that couples to the Fairhead codesuite 
    ! for joint inversions
    use params, only: VaniAllModes_4, VaniAllModes_8, verbose, myGlobalrank, MPI_SPLINE_COMPLEX, & 
                        MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT,   &
                         nmodes, nprocs, all_warnings, datadir, max_tl1, & 
                        Cxyz, MaxBrettModelPts, glob_eta1, glob_eta2, compute_cst_smax, timingNEX, Vani
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                        compute_Vani_matrix, compute_vani_matrix_stored, & 
                        convert_imag_to_real, save_Vani_real_matrix
    use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst_4, Hcomplex_to_cst_8, write_cst_complex_to_file, write_cst_complex_to_file_4
    use mesh_utils, only: find_row_col
    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
    use m_KdTree, only: KdTree, KdTreeSearch
    use voronoi, only: vor_x, vor_y, vor_z, & 
                        vor_A, vor_C, vor_L, vor_N, vor_F, &
                        load_voronoi_model, project_voroni_to_gll
    use specfem_mesh, only: SetMesh, create_SetMesh
    use modes, only: get_mode, Mode 
    use mineos_model, only: mineos, mineos_ptr
    use vani_kernel, only: d_allstrains_r, d_allstrains_i, copy_allstrain_to_device, & 
                           d_wglljac_g, d_Cxyz_g, d_vani_imag_G, d_vani_real_G, & 
                           compute_Vani_onemode_allstrains, check_cuda_device_allocations,& 
                           Model3D_dev, cuda_project_eta_to_GLL, Model_ACLNF, & 
                           d_rad, d_eta1, d_eta2, d_xcoord, d_ycoord, d_zcoord, d_NN1_TOTAL, & 
                           d_modeLUT, compute_Vani_all_modes_at_once
    use cudafor 
    
    use iso_c_binding
    use Vanibindings

    
    implicit none
    include "constants.h"
    include 'mpif.h'

    integer :: iset, i,j,k,s, p, ispec, l1, l2, n1, m1,m2, n2, ierr, & 
                tl2, h, b, gcluster_size, f90cluster_size, JFcluster_size, sets_per_process, & 
                myset_start, myset_end, i_mode, smin, smin_theoretical, smax, smax_theoretical, num_s, & 
                ncols, this_tl1, imode, im, imodel_iter, igll, iproc, & 
                success, idx, thisrow, thiscol, is, it, ncstsvals, icst,& 
                thissmax,allscalars,cxyzsize,LUTsize,m3dsize,irank, GPUSharesize
    character(len=2) nstr, lstr
    character(len=5) iterstr
    character(len=3) chainstr
    character(len=450) :: out_name
    real :: gb_per_set


    real(kind=8), allocatable :: Vani_modesum_r(:), Vani_modesum_i(:)
    character(len=20) :: model_ti
    complex(kind=SPLINE_REAL), allocatable :: allstrains(:, :, :, :, :, :, :)

    integer :: max_nn1

    real(kind=8), allocatable :: wglljac_loc(:, :)
    real(kind=8), allocatable :: Cxyz_loc(:, :, :, :)

    real(kind=4), allocatable :: BrettModelToTransfer_4(:,:)
    real(kind=8), allocatable :: BrettModelToTransfer_8(:,:)

    integer, allocatable, target :: modeLUT(:)

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  

    integer :: thisnn1, ival
    integer :: start_clock, end_clock, count_rate, total_nn1, iii, jjj, loop_clock_start
    real(8) :: elapsed_time

    ! Fairhead stuff: 
    INTEGER         :: MPI_COMM_F90, MPI_COMM_PT, MPI_COMM_JF, MPI_COMM_GPU,  MPIF90colour 
    integer         :: myf90rank, myJFrank, myPTrank, GPUSharerank
    integer         :: iter, upID, NodeID
    integer         :: status(MPI_STATUS_SIZE)
    integer(kind=8) :: FHnvoronoi, niterations, model_start, STORE_FHnvoronoi, NPTrungs, JL_MaxBrettModelPts, Nmodelbuffersize, PTexchangeevery
    integer         :: juliaJFrank
    real(kind=8),allocatable  ::  JLModelBuffer(:)
    real(kind=8), target      :: FHupdates(6)
    integer, allocatable :: dummyrecbuffer(:)


    ! Modes: 
    ! REAL 27:
    type(Mode)             :: mode_1 

    integer, dimension(nmodes), parameter :: modeNs=(/ 3,8,13,6,8,21,9,13,13,15,18,18,23,23,2,3,5,3,9,11,11,16 /)
    integer, dimension(nmodes), parameter :: modeLs=(/ 1,1,1,3,5,6,2,2,3,3,3,4,4,5,3,2,3,8,3,4,5,5 /)
    integer, dimension(nmodes), parameter :: dataSmax=(/ 2,2,2,6,6,6,4,4,6,4,6,6,6,6,4,4,6,6,6,6,6,6 /)

    

    INTEGER :: request1, request2
    INTEGER, dimension(2) :: requests  ! Array of requests


    logical :: using_PT
    ! BINDING PARAMETERS TO CPP : 
    integer :: cppprec
    logical :: cppdouble
    integer :: size_of_array
    integer(kind=8) :: strainsize, straingb
    type(C_PTR) :: ta_ptr, eta1_ptr, eta2_ptr, cxyz_ptr, LUT_ptr, strain_r_ptr, strain_i_ptr, Vani_real_ptr, Vani_imag_ptr, wgll_ptr

    real(4), allocatable :: STORE_flat3Dmodel_4(:)
    real(4), pointer     :: flat3Dmodel_4(:)
    real(8), pointer     :: flat3Dmodel_8(:)
    type(C_PTR) :: ptr_m3D, ptr_FHupdates

    integer, parameter :: HANDLE_SIZE = 128
    character(kind=C_CHAR), dimension(HANDLE_SIZE) :: my_ipc_handle

    integer :: MHAcceptance ! 0 = accept, 1 = reject

    real(4), allocatable, target :: flatarray_4(:), flatstrain_r_4(:), flatstrain_i_4(:), Vani_real_4(:), Vani_imag_4(:)
    real(8), allocatable, target :: flatarray_8(:), flatstrain_r_8(:), flatstrain_i_8(:), Vani_real_8(:), Vani_imag_8(:)


    real(4), allocatable :: Vani_real_4_REDUCED(:), Vani_imag_4_REDUCED(:)


    real(4), dimension(nmodes) :: two_nondimomega

    real(kind=4), allocatable :: allcsts_r_4(:), allcsts_i_4(:), allcsts_r_RED_4(:), allcsts_i_RED_4(:)
    real(kind=8), allocatable :: allcsts_r_8(:), allcsts_i_8(:), allcsts_r_RED_8(:), allcsts_i_RED_8(:)

    complex(kind=4), allocatable :: cst_4(:,:)
    complex(kind=8), allocatable :: cst_8(:,:)


    real(kind=4) :: aclnf_4(5)
    real(kind=8) :: aclnf_8(5), STORE_aclnf_8(5)


    ! Simulation parameters: 
    integer               :: model_chain 
    integer, parameter    :: region      = 3
    character, parameter  :: t1          = 'S'
    logical, parameter    :: force_VTI   = .false.

    
    
    smin = 0

    ! Pre initialise CUDA before F90 
    ierr = fh_cuda_preinit()
    write(*,*)'Pre-initialising CUDA ', ierr

    ! Setup Global MPI with Julia 
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, gcluster_size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myGlobalrank, ierr)

    ! IO 
    IIN  = myGlobalrank
    IOUT = IIN + 2000

    ! Julia broadcasts number 1 or 2 telling us what splits are required
    ! If its not 1 then its the number of PT rungs being used by JF nodes
    call MPI_Bcast(NPTrungs, 1, MPI_INT64_T, 0, MPI_COMM_WORLD, ierr)

    ! If this broadcasted value is 1 then no PT. This means we just need
    ! the split that gives us the F90 only communicator separate from Julia 

    if(NPTrungs.eq.1)then 
        ! Split F90 group communicator from Julia 
        ! the rank is used for ordering of the new procs (key)
        ! These are equivalent 
        MPI_COMM_JF     = MPI_COMM_WORLD
        myJFrank        = myGlobalrank
        JFcluster_size  = gcluster_size

        MPIF90colour = 1
        call MPI_Comm_split(MPI_COMM_JF, MPIF90colour, myGlobalrank, MPI_COMM_F90, ierr)
        call MPI_COMM_RANK(MPI_COMM_F90, myF90rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_F90, f90cluster_size, ierr)

        using_PT = .false.
        GPUSharerank = 0 
    else
        ! Two splits need to occur
        ! (1) PT Comm split - none of these nodes should be in this
        !     Julia color = 2, so use value 3 
        ! We will never use this communicator here
        MPIF90colour = 3
        call MPI_Comm_split(MPI_COMM_WORLD, MPIF90colour, myGlobalrank, MPI_COMM_PT, ierr)
        call MPI_COMM_RANK(MPI_COMM_PT, myPTrank, ierr)


        ! (2) JF Comm split - this breaks up each group of 4 F90 CPUS + 1 JF CPU 
        !     into a group of 5. Currently the method is 
        ! integer division of +ve numbers == floor function 
        MPIF90colour = (myGlobalrank - NPTrungs )/ 4

        call MPI_Comm_split(MPI_COMM_WORLD, MPIF90colour, myGlobalrank, MPI_COMM_JF, ierr)
        call MPI_COMM_RANK(MPI_COMM_JF, myJFrank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_JF, JFcluster_size, ierr)


        ! (3) We need to do the F90 split from the JF node of each group of 5 
        MPIF90colour = 1
        call MPI_Comm_split(MPI_COMM_JF, MPIF90colour, myJFrank, MPI_COMM_F90, ierr)
        call MPI_COMM_RANK(MPI_COMM_F90, myF90rank, ierr) 
        call MPI_COMM_SIZE(MPI_COMM_F90, f90cluster_size, ierr)


        ! (4) For strain sharing on the GPU, we create a Comm groups based on the 
        ! rank within COMM_F90 since this is the GPU they will end up on
        ! They key (sorting in the new set) is done by global rank
        ! the colour is their F90 rank 
        call MPI_Comm_split(MPI_COMM_WORLD, myF90rank, myGlobalrank, MPI_COMM_GPU, ierr)
        call MPI_COMM_RANK(MPI_COMM_GPU, GPUSharerank, ierr) 
        call MPI_COMM_SIZE(MPI_COMM_GPU, GPUSharesize, ierr)

        if (GPUSharesize.ne.NPTrungs)then 
            write(*,*)"Error: GPU Share Comm size should be same as Number of rungs but is ", GPUSharesize, NPTrungs
            stop 
        endif 

        using_PT = .true.

    endif 
    

    if(myJFrank.eq.0)then 
        write(*,*)'ERROR: global rank ', myGlobalrank, "has myJFrank = 0"
        stop 
    else 
        juliaJFrank = 0 
    endif   
   

    ! C++ precision
    cppprec = get_cpp_precision()
    if(cppprec==4)then 
        CPPDOUBLE = .false.
    elseif(cppprec==8)then 
        CPPDOUBLE = .true.
    else 
        write(*,*)'Error, cppprec should be 4 or 8 but was ', cppprec
    endif 

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


    ! Forces each F90rank to its accompanyting device number 
    ierr =  force_proc_to_device(nprocs, myf90rank)

    ! ASSUMING 1 mesh per proc for this optimised code
    sets_per_process = nprocs/f90cluster_size
    myset_start      = myf90rank*sets_per_process
    myset_end        = myset_start + sets_per_process - 1 

    ! Check equal load balance across the processes:
    if( mod(nprocs, f90cluster_size).ne.0)then 
        write(*,*)'Error: you are using '
        write(*,*)'     -- nprocs ', nprocs 
        write(*,*)'     -- nnodes ', f90cluster_size 
        write(*,*)'     (Note the rank 0 is reserved for Julia Fairhead)'
        write(*,*)'and therefore nnodes is not divisible by nnodes. Stop.'
        stop 
        !else    
        ! if(myf90rank.eq.0 .and.verbose.ge.1)write(*,*)'Sets for each node:', sets_per_process
        ! print *, 'F90 process: ', myf90rank, 'does sets', myset_start, 'to ', myset_end
    endif 


    ! Setup mineos only on f90 nodes
    ! POSSIBLE ERROR
    call mineos%load_mineos_radial_info_MPI(myf90rank, MPI_COMM_F90)
    mineos_ptr => mineos


    ! Get the mode frequencies for Vani division: 
    do imode = 1, nmodes
        mode_1  = get_mode(modeNs(imode), 'S', modeLs(imode), mineos_ptr)
        two_nondimomega(imode) = mode_1%wcom * two * SCALE_T
    enddo 


    ! Load mesh data for this proc (1 set per proc)
    ! True false indicates load from disc and dont save to disc
    sm   = create_SetMesh(myset_start, region)
    call sm%setup_mesh_sem_details(.true., .false.)
    call sm%compute_rotation_matrix()


    ! Until I can think of a better system, lets setup a mode look up table on the gpu
    ! this is the number of groups launched to deal with different elements of a matrix
    total_nn1 = 0
    do imode  = 1, nmodes
        l1        = modeLs(imode)
        !this_tl1  = 2*l1 +1
        !total_nn1 = total_nn1 +  (l1+1)*(l1) + 1   ! we dont need separate launches for the lower right box because they are computed when some of the upper right ones are 
        total_nn1 = total_nn1 +  ((l1+1)*(l1+2))/2     ! reduced to only needing the triangular red/blue shading and computing the lower rows in the upper row launches 
    enddo  
    
    allocate(modeLUT(total_nn1*4), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating modeLUT on', myGlobalrank
        stop 
    endif 


    ! Sizes of arrays:
    size_of_array = sm%ngllx * sm%nglly * sm%ngllz * sm%nspec   ! standard mesh scalar
    strainsize    = nmodes * max_tl1 * sm%ngllx * sm%nglly * sm%ngllz * sm%nspec * 6 

    ! Safety check: 
    if(max_tl1.ne. maxval(modeLs)*2 + 1)then 
        write(*,*)'max tl1 is not the max of that in the model1 array - should be ', maxval(modeLs)*2 + 1, ' but is ', max_tl1
        stop 
    endif 


    ! Setup Vani all modes - unique for each proc 
    ! Allocate Vani matrices
    allocate(VaniAllModes_4(max_tl1, max_tl1, nmodes), stat=ierr)
    VaniAllModes_4 = SPLINE_iZERO
    if(ierr.ne.0)then 
        write(*,*)'Error allocating VaniAllModes for proc ', myGlobalrank
        stop 
    endif 


    ! Loading strains and determining estimate of the memory cost: 
    ! Dimension of strains would be: ngllx, nglly, ngllz, nspec, 2*l1+1, 6, nmodes
    ! In double precision (8 bytes) for complex numbers (x2): 
    ! Strain for each mode is about 14 Mb for NEX 176 1 of sets 16
    if(myGlobalrank.eq.NPTrungs)then
        write(*,*)
        write(*,*)'-------------------- GPU MEMORY ESTIMATES --------------------' 
        allscalars = size_of_array * cppprec * 3 ! x, y, z 
        cxyzsize   = size_of_array * cppprec * 36 
        LUTsize    = total_nn1*4   * 4
        m3dsize    = MaxBrettModelPts*5*cppprec
        straingb   = strainsize * cppprec * 2 ! imag + real 
        gb_per_set = real(allscalars+cxyzsize+LUTsize+m3dsize+straingb)/1073741824.0

        write(*,*)' C++ precision          :             ', cppprec
        write(*,*)' Number of modes        :             ', nmodes
        write(*,*)' Standard DP mesh scalar:             ', allscalars, ' bytes'
        write(*,*)' Cxyz matrix            :             ', cxyzsize, ' bytes'
        write(*,*)' Look up table          :             ', LUTsize, ' bytes'
        write(*,*)' 3D Voronoi Model       :             ', m3dsize, ' bytes'
        write(*,*)' Strain arrays          : ', straingb, ' bytes'
        write(*,*)' ........................................................... '
        write(*,*)' Total per processor    :       ', gb_per_set, 'Gb'
        write(*,*)' Total for all sets     :       ', gb_per_set*nprocs, 'Gb'
        write(*,*)'-------------------------------------------------------------' 
        write(*,*)
    endif


    
   
    ! -------------- COPY XYZ COORDINATES TO GPU  --------------
    ! Copy over the xcoord, ycoord, zcoord, rstore arrays: 
    ierr=0
    allocate(flatarray_4(size_of_array), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating flatarray for proc ', myGlobalrank
        stop 
    endif 

    ! Copy x coordinates
    iii = 1
    do k = 1, sm%ngllz
        do j = 1, sm%nglly
            do i = 1, sm%ngllx
                do ispec = 1, sm%nspec 
                        flatarray_4(iii) = real(sm%xstore(i,j,k,ispec), kind=4)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ta_ptr = c_loc(flatarray_4)
    ierr = copythisarraytodevice(ta_ptr, size_of_array, 1)

    ! Copy y coordinates
    iii = 1
    do k = 1, sm%ngllz
        do j = 1, sm%nglly
            do i = 1, sm%ngllx
                do ispec = 1, sm%nspec 
                        flatarray_4(iii) = real(sm%ystore(i,j,k,ispec),kind=4)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ta_ptr = c_loc(flatarray_4)
    ierr = copythisarraytodevice(ta_ptr, size_of_array, 2)

    ! Copy z coordinates
    iii = 1
    do k = 1, sm%ngllz
        do j = 1, sm%nglly
            do i = 1, sm%ngllx
                do ispec = 1, sm%nspec 
                        flatarray_4(iii) = real(sm%zstore(i,j,k,ispec),kind=4)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ta_ptr = c_loc(flatarray_4)
    ierr = copythisarraytodevice(ta_ptr, size_of_array, 3)

    ! -------------- END OF COPYING XYZ TO GPUS --------------



    ! ----------------- Determine number of CSTS ----------------------
    ncstsvals = 0
    do imode = 1, nmodes
        ! There are (l+1) s values we need (for self coupling) where 
        ! Now dynamic smax based on real data: 
        ! We now also include s=0 
        thissmax = dataSmax(imode)
        do s = smin, thissmax, 2 
            ncstsvals = ncstsvals + (s+1)
        enddo 
    enddo 


    allocate(allcsts_r_4(ncstsvals), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating allcsts_r on ', myGlobalrank
        stop
    endif 
    allocate(allcsts_i_4(ncstsvals), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating allcsts_r on ', myGlobalrank
        stop
    endif 

    ! ONLY ALLOCATED ON THE LEAD NODE OF ONE OF THE 4 PROCS in F90 
    ! DO WE STILL NEED TO ALLOCATE THIS? I DONT THINK SO
    if(myf90rank.eq.0)then 
        allocate(allcsts_r_RED_4(ncstsvals), stat=ierr)
        if(ierr.ne.0)then 
            write(*,*)'Error allcsts_r_RED allcsts_r'
            stop
        endif 
        allocate(allcsts_i_RED_4(ncstsvals), stat=ierr)
        if(ierr.ne.0)then 
            write(*,*)'Error allcsts_i_RED allcsts_r'
            stop
        endif 
    endif 






    ! -------------- COPY STRAINS TO GPU  --------------
    ! Also goes into this first case if not using PT 
    if (GPUSharerank.eq.0)then 
        ! Allocate mega strain array
        allocate(allstrains(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, max_tl1, 6, nmodes), stat=ierr)
        if(ierr.ne.0)then 
            write(*,*)'Error allocating allstrains for proc ', myGlobalrank
            stop 
        endif 
        allstrains = SPLINE_iZERO

        ! Load all the strains in a loop
        do imode = 1, nmodes
            n1       = modeNs(imode)
            l1       = modeLs(imode)
            this_tl1 = 2*l1 +1
            do im =  -l1, l1 
                call sm%load_mode_strain_binary(n1, t1, l1, im, &  
                                                allstrains(:,:,:,:,l1+im+1,:,imode))
            enddo 
        enddo  

        if(GPUSharerank.eq.0)write(*,*)'Loaded all strain binaries from disc.'

        ! Compute the size of all the strains:
        ! (sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, max_tl1, 6, nmodes)
        ! Now we have collapsed the matrix we can store only the number of tl1s 
        ! that each mode needs 
        ! FOR NOW WE ARE USING THE MORE MEMORY INEFFICIENT VERSION OF ASSUMING THEY ALL 
        ! HAVE max_tl1 - this makes the indexing a bit simpler in the kernel for now.
        ! probs need a LUT otherwise
        allocate(flatstrain_r_4(strainsize), stat=ierr)
        if(ierr.ne.0)then 
            write(*,*)'Error in flatstrain_r on ', myGlobalrank
            stop 
        endif 
        allocate(flatstrain_i_4(strainsize), stat=ierr)
        if(ierr.ne.0)then 
            write(*,*)'Error in flatstrain_i on ', myGlobalrank
            stop 
        endif 

        ! The strain mega-array
        iii = 1
        do imode = 1, nmodes 
            l1  = modeLs(imode)
            do p = 1,6
                do im =  1, max_tl1!-l1, l1 
                    do k = 1, sm%ngllz
                        do j = 1, sm%nglly
                            do i = 1, sm%ngllx
                                do ispec = 1, sm%nspec 
                                    if (im.le. 2*l1 + 1)then 
                                        flatstrain_r_4(iii) =  real(allstrains(i, j, k, ispec, im, p, imode) *   sm%wglljac(i,j,k,ispec)**half ) 
                                        flatstrain_i_4(iii) = aimag(allstrains(i, j, k, ispec, im, p, imode) *   sm%wglljac(i,j,k,ispec)**half ) 
                                    else
                                        ! regions outside of the tl1 that is valid for this mode
                                        flatstrain_r_4(iii) = zero 
                                        flatstrain_i_4(iii) = zero 
                                    endif 
                                    iii = iii + 1
                                enddo 
                            enddo 
                        enddo 
                    enddo
                enddo 
            enddo 
        enddo 
        strain_r_ptr = c_loc(flatstrain_r_4)
        strain_i_ptr = c_loc(flatstrain_i_4)    

        if(using_PT)then 
            ! Need to transfer the mega-strains! 
            ierr = copy_allstrains_fromrank0(strain_r_ptr, strain_i_ptr, strainsize, myf90rank, my_ipc_handle)
            if(ierr.ne.0)then 
                write(*,*)'Error in copy_allstrains on ', myGlobalrank
                stop 
            endif 

            deallocate(allstrains) 
            deallocate(flatstrain_r_4)
            deallocate(flatstrain_i_4)

            ! Broadcast the IPC handle to others in the group 
            call MPI_Bcast(my_ipc_handle, HANDLE_SIZE, MPI_CHARACTER, 0, MPI_COMM_GPU, ierr)
        else 
            ! Still copy strains but we dont need the handles; 
            ! could probabily just use the rank0 function but i dont know
            ! if the shared handles slows it down a bit? 
            ierr = copy_allstrains(strain_r_ptr, strain_i_ptr, strainsize)
        endif 
    else 
        ! I am not the rank 0 on this GPU group so I will be sent the
        ! handles for CUDA IPC 
        call MPI_Bcast(my_ipc_handle, HANDLE_SIZE, MPI_CHARACTER, 0, MPI_COMM_GPU, ierr)

        ierr = copy_allstrains_higherranks(myf90rank, my_ipc_handle)
    endif !GPUSharerank 




    ! Create look up table for how each cuda group maps to a matrix elemet 
    idx = 1
    do imode = 1, nmodes
        l1 = modeLs(imode)

        !this_tl1 = 2*l1 +1
        !thisnn1  =  (l1+1)*(l1) + 1 !+ (l1+1)*(l1)/2
        thisnn1  =  ((l1+1)*(l1+2))/2

        do iii = 1, thisnn1
            modeLUT(idx) = imode - 1 ! mode
            idx = idx + 1
            modeLUT(idx) = iii - 1   ! place in matrix
            idx = idx + 1            ! global index 

            ! Get the row and column of this point in the matrix
            ! ie the m1, m2: 
            ! Note that the find_row_col fucks up the iii value so
            ! copy it to a new integer before parsing -- this took me
            ! way too long to work out
            ival = iii 

            call find_row_col(ival, thisrow, thiscol, l1)

            modeLUT(idx) = thisrow -1  ! row in c++ -1 
            idx = idx + 1

            modeLUT(idx) = thiscol -1  ! col in c++ -1 
            idx = idx + 1
        enddo 
    enddo  
    
    ! Compute the maximum number of matrix entries we need the kernel 
    ! to compute based on the largest mode. We store max_tl1
    ! but the number of elements is now l+1^2 so
    max_nn1 = ((max_tl1 - 1)/2  + 1 )**2 
    LUT_ptr = c_loc(modeLUT)
    ierr = copy_LUT_array(LUT_ptr, total_nn1*4)
    if(ierr.ne.0)then 
        write(*,*)'Error in copy_LUT_array on ', myGlobalrank
        stop
    endif 



    ! Allocate the arrays for the output Vani matrix (flattened)
    ierr = allocate_Vani_arrays(nmodes*max_nn1)
    if(ierr.ne.0)then 
        write(*,*)'Error in allocate_Vani_arrays on', myGlobalrank
        stop 
    endif 
    

    ! Local F90 versions on host cpu 
    allocate(Vani_real_4(nmodes*max_nn1), stat=ierr)
    Vani_real_ptr = c_loc(Vani_real_4)
    allocate(Vani_imag_4(nmodes*max_nn1), stat=ierr)
    Vani_imag_ptr = c_loc(Vani_imag_4)

    if(ierr.ne.0)then 
        write(*,*)'Error in allocating Vani_real on', myGlobalrank
        stop 
    endif 

    ! Allocate the Cxyz array
    ierr = allocate_Cxyz_array(size_of_array*36)
    if(ierr.ne.0)then 
        write(*,*)'Error in function allocate_Cxyz_array', myGlobalrank
        stop 
    endif 


    ! The flat model array that is used on the GPU
    allocate(flat3Dmodel_4(MaxBrettModelPts*5), stat=ierr)
    allocate(STORE_flat3Dmodel_4(MaxBrettModelPts*5), stat=ierr)
    ! Get ptr to host flatmodel array
    ptr_m3D = c_loc(flat3Dmodel_4)

    if(ierr.ne.0)then 
        write(*,*)'Error allocating flat3Dmodel on', myGlobalrank
        stop 
    endif 

    ! Allocate Model array on the device
    ierr = allocate_M3D_array(MaxBrettModelPts*5)
    if(ierr.ne.0)then 
        write(*,*)'Error running allocate_M3D_array on', myGlobalrank
        stop 
    endif 



    ! -------- BROADCASTING MODEL/MCMC PARAMS FROM JULIA BEGINS --------
    
    ! Get number of model iterations: 
    call MPI_Bcast(niterations, 1, MPI_INT64_T, juliaJFrank, MPI_COMM_JF, ierr)
    call MPI_Bcast(model_start, 1, MPI_INT64_T, juliaJFrank, MPI_COMM_JF, ierr)

    ! Ensure maximum number of volumes is consistent 
    call MPI_Bcast(JL_MaxBrettModelPts, 1, MPI_INT64_T, juliaJFrank, MPI_COMM_JF, ierr)
    if(JL_MaxBrettModelPts.ne.MaxBrettModelPts)then 
        write(*,*)"Error: max model points not same for JL/F90:  ", & 
                  JL_MaxBrettModelPts, MaxBrettModelPts
        stop 
    endif 

    ! Allocate full model buffer once and for all 
    ! X, Y, Z, eta1, eta2 for N nodes, + 5 for ACLNF + 1 for number of nodes
    Nmodelbuffersize = 1 + (MaxBrettModelPts+1)*5
    allocate(JLModelBuffer(Nmodelbuffersize))

    ! Recieve the first model 
    call MPI_Bcast(JLModelBuffer, Nmodelbuffersize, MPI_DOUBLE_PRECISION, juliaJFrank, MPI_COMM_JF, ierr)

    FHnvoronoi = JLModelBuffer(1)
    aclnf_8    = JLModelBuffer(2:6)

    ! Copy over the original model to GPU:
    flat3Dmodel_4(:) = zero 
    flat3Dmodel_4(1:int(FHnvoronoi)*5) = real(JLModelBuffer(7:7+FHnvoronoi*5 - 1), kind=4)
    ierr = copy_M3D_array(ptr_m3D, MaxBrettModelPts*5)


    ! Ensure that there is a pointer to the updates array: 
    ptr_FHupdates = c_loc(FHupdates)

    if(myf90rank.eq.0)then 
        write(*,*)
        write(*,*)' NMSPLIT90 has: '
        write(*,*)' Fairhead starts at               : ', model_start, myGlobalrank
        write(*,*)' Fairhead number of iterations    : ', niterations, myGlobalrank
        write(*,*)' Fairhead initial model points    : ', FHnvoronoi,  myGlobalrank
        write(*,*)' Fairhead initial model ACLNF     : ', aclnf_8,     myGlobalrank
        write(*,*)
    endif 



    ! Julia needs to know the number of cst values it will receive from each rank
    ! should be identical but Julia will check this 
    allocate(dummyrecbuffer(JFcluster_size))
    call MPI_Gather(ncstsvals, 1, MPI_INTEGER, dummyrecbuffer, 1, MPI_INTEGER, juliaJFrank, MPI_COMM_JF, ierr)


    ! If using PT get how often the switches occur: 
    if(using_PT)then 
        call MPI_Bcast(PTexchangeevery, 1, MPI_INT64_T, juliaJFrank, MPI_COMM_JF, ierr)
    endif 

    ! For now we will ignore any mode calcluations for the first iteration since its never going to
    ! be used anyway 

    ! ! ------------------------------------------------------------------------------------------
    ! ! THINGS THAT HAPPEN EACH MODEL 
    ! Fairhead is going to generate a perturbation to the system and send that over. 
    ! The options are: 
    !   ID      Action          N vals sent     Ordering
    !  1,2,3    Modify ACLNF    3               ID, index of ACLNF (1-5), value
    !  4,5,6    Modify angles   3               ID, nodeID, new eta1, new eta2
    !  7        Move node       5               ID, nodeID, new x, new y, new z
    !  8        Birth node      7               ID, new x, new y, new z, index of ACLNF, new value * 
    !  9        Kill node       2               ID, nodeID
    ! All values sent as floats so can be sent in single MPI call 

    ! We are about to start iterations so lets get everyhthing up to date
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    

    !call system_clock(count_rate=count_rate)
    ! Iterations start at 0 because we need the 0th run through to 
    ! get the NLL before the julia iterations start. 
    do iter = model_start-1, niterations
        
        ierr = copy_M3D_array(ptr_m3D, MaxBrettModelPts*5)
        ierr = cpp_project_eta_to_gll(int(FHnvoronoi), sm%nspec, sm%ngllx, aclnf_8)     

        ! The rank being parsed is just for debugging
        ierr = launch_vanikernel(sm%ngllx, sm%nspec, total_nn1, max_nn1, max_tl1, nmodes, myGlobalrank)

        ! if (GPUSharerank.eq.0)then 
        !     if (mod(iter, 100) == 0) then  ! Every 100 iterations
        !         ierr = check_gpu_utilization(myF90rank)
        !     endif
        ! endif 

        ierr = copyfromdevice(Vani_real_ptr, max_nn1*nmodes, 8)
        ierr = copyfromdevice(Vani_imag_ptr, max_nn1*nmodes, 9)

        ! Overall the number of cst values we need to compute and 
        ! conduct reduction of is as follows: 
        ! Each mode of degree l = we need the s from 0 to 2l but only using smax based on data (inclusive)
        ! for each s there are s+1 values we need to transfer (only computing the negative)
        ! Due to hermitian nature 
        ! call system_clock(count_rate=count_rate)
        ! call system_clock(start_clock)
        icst = 1
        do imode = 1, nmodes

            l1  = modeLs(imode)
            n1  = modeNs(imode)
            this_tl1 = 2*l1 + 1

            ! Now need all the values: 
            do iii = 1, (l1+1)*(l1+1) ! thisnn1
                ival = iii 

                call find_row_col(ival, thisrow, thiscol, l1)
                        
                if(imode.eq.4.and.myf90rank.eq.0) write(*,*)iii, thisrow, thiscol
                ! if(thisrow.eq.l1+1 .and. thiscol.gt.l1+1)then
                !     ! do nothing for now 
                ! else 
                ! Normal index
                VaniAllModes_4(thisrow, thiscol, imode) = Vani_real_4((imode-1)*max_nn1 + iii) + & 
                                              SPLINE_iONE*Vani_imag_4((imode-1)*max_nn1 + iii)
        
                m1 = thisrow - l1 -1 
                m2 = thiscol - l1 -1 
            
                ! Maps the lower right triangular to the top left triangular
                if(m1.gt.0)then 


                    ! nm1 = 2*l1 + 2 - thiscol   
                    ! nm2 = 2*l1 + 2 - thisrow   

                    VaniAllModes_4(this_tl1 - thiscol + 1, this_tl1 - thisrow + 1, imode) = VaniAllModes_4(thisrow, thiscol, imode) * (-one)**real( (thisrow + thiscol - two*(l1 +1)) ,kind=8)
                    ! if(imode.eq.23 .and.myGlobalrank.eq.10)write(*,*)"2:: ", this_tl1 - thiscol + 1, this_tl1 - thisrow + 1, VaniAllModes_4(this_tl1 - thiscol + 1, this_tl1 - thisrow + 1, imode)
                    if(imode.eq.4.and.myf90rank.eq.0)write(*,*)"  maps to  ", this_tl1 - thiscol + 1, this_tl1 - thisrow + 1
                else 
                    ! if(thisrow+thiscol.ne.this_tl1+1)then 
                    !     ! Avoids the diagonal from centre to top right - others are reflected

                    !     if(imode.eq.4.and.myf90rank.eq.0)write(*,*)"  maps to  ", l1 + 1 -(thiscol - l1 - 1 ), l1 + 1  - (thisrow - l1 - 1)


                    !     !VaniAllModes_4(thiscol, this_tl1-thisrow+1, imode) = VaniAllModes_4(thisrow, thiscol, imode) *  ((-one)**real( l1 - thisrow +1 , kind=8 ))

                    !     m1 = thisrow - l1 -1 
                    !     m2 = thiscol - l1 -1 

                    !     ! Shouldnt always reflect? 
                    !     VaniAllModes_4( l1 + 1 -(thiscol - l1 - 1 ), l1 + 1  - (thisrow - l1 - 1) , imode) = VaniAllModes_4(thisrow, thiscol, imode) * ((-one)**real( m1 + m2 , kind=8 ))
                    ! endif 
                endif 

              
                if(imode.eq.4.and.myf90rank.eq.0)write(*,*)

                !endif
            enddo 

            if(imode.eq.4)then 
                if(myf90rank.eq.0)then 
                    l1 = modeLs(imode)
                    allocate(Vani(l1, l1))
                    Vani = VaniAllModes_4(:, :, imode)
                    call save_Vani_matrix(l1, l1, "ranktest.txt", .false.)
                endif

                stop 
            endif 



            ! smax is now the dataSmax
            call get_Ssum_bounds(l1, l1, smin_theoretical, smax_theoretical, num_s, ncols)
            allocate(cst_4(num_s, ncols))
            call Hcomplex_to_cst_4(VaniAllModes_4(1:this_tl1, 1:this_tl1, imode)/two_nondimomega(imode) , &
                                   l1, l1, cst_4, ncols, num_s, t1, t1, 2)


            ! Starts at smin (s=2)
            thissmax = dataSmax(imode)
            do is = smin+1, thissmax+1, 2
                do it = 1, (is-1) + 1
                    allcsts_r_4(icst) =  real(cst_4(is,it))
                    allcsts_i_4(icst) = aimag(cst_4(is,it))
                    icst = icst + 1
                enddo 
            enddo 
            deallocate(cst_4)
        enddo ! imode  


        ! Let us now reduce the csts directlry on the Julia node 
        call MPI_Reduce(allcsts_r_4, allcsts_r_RED_4, ncstsvals, MPI_REAL, &
                        MPI_SUM, juliaJFrank, MPI_COMM_JF, ierr) 
        call MPI_Reduce(allcsts_i_4, allcsts_i_RED_4, ncstsvals, MPI_REAL, &
                        MPI_SUM, juliaJFrank, MPI_COMM_JF, ierr) 

        ! Listen for decision from the MH accepta/rejectance
        ! For the pre-loop NLL evaluation we send an acceptance
        ! this means that we dont update the model, but will assign the 
        ! model to the store
        call MPI_Bcast(MHAcceptance, 1, MPI_INT64_T, juliaJFrank, MPI_COMM_JF, ierr)


        if(MHAcceptance.eq.1)then 
            ! rejected
            ! need to reset the arrays to stored model:
            aclnf_8(:)        = STORE_aclnf_8(:)
            FHnvoronoi        = STORE_FHnvoronoi
            flat3Dmodel_4(:)  = STORE_flat3Dmodel_4(:)
        endif


        ! This is where PT goes 
        if (iter.ne.model_start-1 .and. MOD(iter, PTexchangeevery).eq.0 .and. using_PT)then 
            ! Swaps happen on JL nodes - we dont care about this 
            ! Regardless of if swap occurs, we broadcast the model
            ! since latency attached to opening the MPI line is the 
            ! main overhead 
            ! Recieve the first model 
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            JLModelBuffer(:) = zero 
            call MPI_Bcast(JLModelBuffer, Nmodelbuffersize, MPI_DOUBLE_PRECISION, juliaJFrank, MPI_COMM_JF, ierr)

            ! Update mode
            FHnvoronoi = JLModelBuffer(1)
            aclnf_8    = JLModelBuffer(2:6)

            flat3Dmodel_4(:) = zero 
            flat3Dmodel_4(1:int(FHnvoronoi)*5) = real(JLModelBuffer(7:7+FHnvoronoi*5 - 1), kind=4)

            ! Store this as our new base model: 
            STORE_aclnf_8(:)    = aclnf_8(:)
            STORE_FHnvoronoi    = FHnvoronoi
            STORE_flat3Dmodel_4 = flat3Dmodel_4(:)
        endif 



        ! Listen for updates to the next model
        if(iter.ne.niterations)then 
            ! Store the current model
            if(MHAcceptance.eq.0)then 
                ! If accepted then we need to update the 'store' 
                ! otherwise the store from last iter is same 
                ! as the current values 
                STORE_aclnf_8(:)    = aclnf_8(:)
                STORE_FHnvoronoi    = FHnvoronoi
                STORE_flat3Dmodel_4 = flat3Dmodel_4(:)
            endif 


            ! Listen for the the next iteration: 
            call MPI_Bcast(FHupdates, 6, MPI_DOUBLE_PRECISION, juliaJFrank, MPI_COMM_JF, ierr)

            ! Process the update
            ! 3rd argument is the size of the M3D array
            ! returns the updated number of voronoi points
            !npts = update_FH_model(ptr_FHupdates, ptr_m3D, MaxBrettModelPts*5)
            upID   = int(FHupdates(1))
            ! Node ID wont always be sent but, if it is, then its in slot 2
            ! index the nodes starting at 0 
            NodeID = int(FHupdates(2))-1 

            if(upID.le.3)then 
                ! Velocity change
                aclnf_8(int(FHupdates(2))) = FHupdates(3)
            elseif(upID.ge.4.and.upID.le.6)then
                ! Symmetry axis update - eta1, eta2 are send as
                flat3Dmodel_4(NodeID*5 +4 : NodeID*5 +5) = FHupdates(3:4)
            elseif(upID.eq.7)then
                ! Move node - update x,y,z
                flat3Dmodel_4(NodeID*5 +1 : NodeID*5 +3) = FHupdates(3:5)
            elseif(upID.eq.8)then
                ! Birth node
                ! Coordinates for new node
                flat3Dmodel_4(FHnvoronoi*5 +1: FHnvoronoi*5 +3) = FHupdates(2:4)
                ! New angles 
                flat3Dmodel_4(FHnvoronoi*5 +4: FHnvoronoi*5 +5) = FHupdates(5:6)
                FHnvoronoi = FHnvoronoi + 1
            elseif(upID.eq.9)then
                ! This is trickier. First we need to identify the node and then shuffle down the nodes values that
                ! are above it in the array so its overwritten
                ! Let s_i and f_i be the start and end indices of each node
                ! s_i = (i-1)*5 + 1   -->  s_i = NodeID*5 + 1
                ! f_i = 5*i           -->  f_i = (NodeID+1)*5
                ! Note that NodeID = i-1 
                ! The replacement is 
                ! Model[ s_i : f_{n-1}]  is given Model[ s_{i+1} : f_{n} ]        
                flat3Dmodel_4(NodeID*5 +1 : (FHnvoronoi-1)*5) = flat3Dmodel_4( (NodeID+1)*5 +1: FHnvoronoi*5)
                ! Now anything above this should be given a 0 
                flat3Dmodel_4((FHnvoronoi-1)*5 +1: MaxBrettModelPts*5) = zero 
                FHnvoronoi = FHnvoronoi - 1
            elseif(upID.eq.10)then
                ! BURNIN proposal 
                ! This is only used in the burn-in window so we can afford a little
                ! inefficiency - the aim here is that we dont want to use a bigger buffer
                ! as default (for cases 1-9) so use a second broadcast here
                ! this is analogous to the PT broadcast
                JLModelBuffer(:) = zero 
                call MPI_Bcast(JLModelBuffer, Nmodelbuffersize, MPI_DOUBLE_PRECISION, juliaJFrank, MPI_COMM_JF, ierr)

                ! Update mode
                FHnvoronoi = JLModelBuffer(1)
                aclnf_8    = JLModelBuffer(2:6)

                flat3Dmodel_4(:) = zero 
                flat3Dmodel_4(1:int(FHnvoronoi)*5) = real(JLModelBuffer(7:7+FHnvoronoi*5 - 1), kind=4)
            else
                !ERROR 
                write(*,*)"Error updating FH model. ID is not 1-9:", upID
                stop
            endif

        endif ! listening for new  
        
    enddo !iterations



    ! ! ------------------------------------------------------------------------------------------
    ! ! STUFF TO DO AFTER ALL THE ITERATIONS ARE COMPLETE: 
    ! ! Cleanup memory 
    ! deallocate(VaniAllModes_4)


    call mpi_finalize(ierr)
end program fairhead_optimised_vani