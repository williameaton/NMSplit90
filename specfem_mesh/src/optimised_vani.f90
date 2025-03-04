
program optimised_vani
    use params, only: VaniAllModes, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                        MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT,   &
                         nmodes, nprocs, all_warnings, datadir, max_tl1, & 
                        Cxyz, MaxBrettModelPts, glob_eta1, glob_eta2, compute_cst_smax, timingNEX
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                        compute_Vani_matrix, compute_vani_matrix_stored, & 
                        convert_imag_to_real, save_Vani_real_matrix
    use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst, write_cst_complex_to_file
    use mesh_utils, only: find_row_col
    use v_ani, only: cuda_Vani_matrix_stored_selfcoupling
    use m_KdTree, only: KdTree, KdTreeSearch
    use voronoi, only: vor_x, vor_y, vor_z, & 
                        vor_A, vor_C, vor_L, vor_N, vor_F, &
                        load_voronoi_model, project_voroni_to_gll
    use specfem_mesh, only: SetMesh, create_SetMesh
    use modes, only: get_mode, Mode 
    use mineos_model, only: mineos, mineos_ptr
    use model3d, only: M3D
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
                tl2, h, b, cluster_size, sets_per_process, & 
                myset_start, myset_end, i_mode, smin, smax, num_s, & 
                ncols, this_tl1, imode, im, imodel_iter, igll, iproc, & 
                success, idx, thisrow, thiscol, is, it, ncstsvals, icst, thissmax
    character(len=2) nstr, lstr
    character(len=5) iterstr
    character(len=3) chainstr
    character(len=450) :: out_name
    real :: gb_per_set

    real(kind=8), allocatable :: allcsts_r(:), allcsts_i(:), allcsts_r_RED(:), allcsts_i_RED(:)

    real(kind=8), allocatable :: Vani_modesum_r(:), Vani_modesum_i(:)
    complex(kind=SPLINE_REAL), allocatable :: cst(:,:)
    character(len=20) :: model_ti
    complex(kind=SPLINE_REAL), allocatable :: allstrains(:, :, :, :, :, :, :)

    integer :: max_nn1

    real(kind=8), allocatable :: wglljac_loc(:, :)
    real(kind=8), allocatable :: Cxyz_loc(:, :, :, :)

    real(kind=8), allocatable :: BrettModelToTransfer(:,:)

    integer, allocatable, target :: modeLUT(:)

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  
    type(M3D) :: model3D

    integer :: thisnn1, ival
    integer :: start_clock, end_clock, count_rate, total_nn1, iii, loop_clock_start
    real(8) :: elapsed_time

    ! Modes: 
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 5, 6, 8, 8, 9, 9, 11, 11, 11, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21, 21, 22, 23, 23, 25, 25, 27/)
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 2, 3, 1, 5, 3, 4,  4,  5,  6,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7,  8,  1,  4,  5,  1,  2,  2/)
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 5, 6, 8, 8, 9, 9, 11, 11, 11, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21/) 
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 2, 3, 1, 5, 3, 4,  4,  5,  6,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7/)
    
    ! Deuss figure modes
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 6, 8, 8, 9, 11, 11, 13, 13, 13,  16, 16, 17,  18, 18, 21, 21, 23 , 23 /) 
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 3, 1, 5, 3,  4,  5,  1,  2,  3,   5,  7,  1,   3,  4,  6,  7,  4 ,  5 /)
    
    ! The 33
    !integer, dimension(33), parameter :: modeNs = (/7, 27, 9, 5, 17, 16, 3, 23, 3, 8, 11, 18, 21, 3, 16, 13, 6, 13, 21, 2,  8,  7, 23, 11, 13, 18, 21,5, 27, 9, 22, 15, 14/)
    !integer, dimension(33), parameter :: modeLs = (/4, 1,  3, 3, 1,  7,  2, 4,  8, 5, 5,   3, 7,  1,  5,  3, 3,  2,  6,  3, 1,  5,  5,  4,  1,  4,  8,2,  2, 2,  1,  3,  4/)


    ! Half the 33
    !integer, dimension(nmodes), parameter :: modeNs = (/7, 27, 9, 5, 17, 16, 3, 23, 3, 8, 11, 18, 21, 3, 16, 13/)
    !integer, dimension(nmodes), parameter :: modeLs = (/4, 1,  3, 3, 1,  7,  2, 4,  8, 5, 5,   3, 7,  1,  5,  3/)
    integer :: sumoftl1s

    integer, dimension(nmodes), parameter :: modeNs = (/ 13, 3, 16, 23/) 
    integer, dimension(nmodes), parameter :: modeLs = (/  1, 2,  7,  4/)

    real(kind=8) :: testval

    ! BINDING PARAMETERS: 
    integer     :: size_of_array, strainsize
    type(C_PTR) :: ta_ptr, eta1_ptr, eta2_ptr, cxyz_ptr, LUT_ptr, strain_r_ptr, strain_i_ptr, Vani_real_ptr, Vani_imag_ptr, wgll_ptr

    real(8), pointer :: flat3Dmodel(:)
    type(C_PTR) :: ptr_m3D


    real(8), allocatable, target :: flatarray(:), eta1(:), eta2(:), localcxyz(:), flatstrain_r(:), flatstrain_i(:), Vani_real(:), Vani_imag(:)
    real(8), allocatable, target :: wgllflat(:)


    ! Simulation parameters: 
    integer, parameter    :: region      = 3
    integer               :: model_chain 
    integer, parameter    :: nmodeliter  = 11000
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


    if(force_VTI)then 
        model_ti = '_VTI'
    else 
        model_ti = ''
    endif 

    ! For some reason Myrank =6 seems to not print...is this an issue? 
    call check_cuda_device_allocations(nprocs)
    

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

    IIN  = myrank
    IOUT = IIN + 2000

    ! Setup mineos: 
    call mineos%load_mineos_radial_info_MPI()
    mineos_ptr => mineos

    ! Load mesh data for this proc (1 set per proc)
    ! True false indicates load from disc and dont save to disc
    sm   = create_SetMesh(myset_start, region)
    call sm%setup_mesh_sem_details(.true., .false.)
    call sm%compute_rotation_matrix()



    ! Safety check: 
    if(max_tl1.ne. maxval(modeLs)*2 + 1)then 
        write(*,*)'max tl1 is not the max of that in the model1 array - should be ', maxval(modeLs)*2 + 1, ' but is ', max_tl1
        stop 
    endif 


    ! Setup global eta1, eta2 arrays
    ! Allocate Vani matrices
    allocate(VaniAllModes(max_tl1, max_tl1, nmodes))
    VaniAllModes = SPLINE_iZERO



    ! Loading strains and determining estimate of the memory cost: 
    ! Dimension of strains would be: ngllx, nglly, ngllz, nspec, 2*l1+1, 6, nmodes
    ! In double precision (8 bytes) for complex numbers (x2): 

    ! Strain for each mode is about 14 Mb for NEX 176 1 of sets 16
    if(myrank.eq.0)then
        ! For one set:
        gb_per_set = sm%ngllx*sm%nglly*sm%ngllz*sm%nspec* SIX * EIGHT * TWO * max_tl1 * nmodes/1.e9
        write(*,*)'Estimated strain array size ', gb_per_set,' Gb' 
        write(*,*)'For all sets:               ', gb_per_set * nprocs, ' Gb' 
    endif


    ! Load all strains once and for all for each m value
    ierr=0
    allocate(allstrains(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, max_tl1, 6, nmodes), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating allstrains for proc ', myrank
        stop 
    endif 


    ! Allocate memory on the GPU device: 
    size_of_array = sm%ngllx * sm%nglly * sm%ngllz * sm%nspec

    ! Copy over the xcoord, ycoord, zcoord, rstore arrays: 
    allocate(flatarray(size_of_array))

    ! Copy x coordinates
    iii = 1
    do ispec = 1, sm%nspec 
        do k = 1, sm%ngllz
            do j = 1, sm%nglly
                do i = 1, sm%ngllx
                    flatarray(iii) = sm%xstore(i,j,k,ispec)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ta_ptr = c_loc(flatarray)
    ierr = copythisarraytodevice(ta_ptr, size_of_array, 1)

    ! Copy y coordinates
    iii = 1
    do ispec = 1, sm%nspec 
        do k = 1, sm%ngllz
            do j = 1, sm%nglly
                do i = 1, sm%ngllx
                    flatarray(iii) = sm%ystore(i,j,k,ispec)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ta_ptr = c_loc(flatarray)
    ierr = copythisarraytodevice(ta_ptr, size_of_array, 2)


    ! Copy z coordinates
    iii = 1
    do ispec = 1, sm%nspec 
        do k = 1, sm%ngllz
            do j = 1, sm%nglly
                do i = 1, sm%ngllx
                    flatarray(iii) = sm%zstore(i,j,k,ispec)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ta_ptr = c_loc(flatarray)
    ierr = copythisarraytodevice(ta_ptr, size_of_array, 3)


    ! Transfer the wgll: 
    allocate(wgllflat(size_of_array))
    iii = 1
    do k = 1, sm%ngllz
        do j = 1, sm%nglly
            do i = 1, sm%ngllx
                do ispec = 1, sm%nspec 
                    wgllflat(iii) = sm%wglljac(i,j,k,ispec)
                    iii = iii + 1
                enddo 
            enddo 
        enddo 
    enddo 



    wgll_ptr = c_loc(wgllflat)
    ierr = copy_wgll_array(wgll_ptr, size_of_array)




    allstrains = SPLINE_iZERO

    ncstsvals = 0 
    do imode = 1, nmodes
        n1       = modeNs(imode)
        l1       = modeLs(imode)
        this_tl1 = 2*l1 +1
        do im =  -l1, l1 
            call sm%load_mode_strain_binary(n1, t1, l1, im, &  
                                            allstrains(:,:,:,:,l1+im+1,:,imode))
        enddo 

        ! There are (l+1) s values we need (for self coupling) where 
        ! Each of those s there are s + 1 
        ! Maxing at s = 6
        if (2*l1.gt.compute_cst_smax)then 
            thissmax = compute_cst_smax
        else 
            thissmax = 2*l1 
        endif 
        do s = 0, thissmax, 2 
            ncstsvals = ncstsvals + (s+1)
        enddo 
        
    enddo  
    if(myrank.eq.0)write(*,*)'Loaded all strain binaries from disc.'


    allocate(allcsts_r(ncstsvals))
    allocate(allcsts_i(ncstsvals))
    
    if(myrank.eq.0)then 
        allocate(allcsts_r_RED(ncstsvals))
        allocate(allcsts_i_RED(ncstsvals))
    endif 
     


    ! Until I can think of a better system, lets setup a mode look up table on the gpu
    total_nn1 = 0
    do imode = 1, nmodes
        l1       = modeLs(imode)
        this_tl1 = 2*l1 +1
        total_nn1 = total_nn1 + (this_tl1*(this_tl1+1)/2 - l1*(l1+1)/2)
    enddo  
    allocate(modeLUT(total_nn1*4))


    idx = 1
    do imode = 1, nmodes
        l1       = modeLs(imode)
        this_tl1 = 2*l1 +1

        thisnn1  = (this_tl1 * (this_tl1+1)/2) - l1*(l1+1)/2

        do iii = 1, (this_tl1 * (this_tl1+1)/2) - l1*(l1+1)/2
            modeLUT(idx) = imode-1 ! mode
            idx = idx + 1
            modeLUT(idx) = iii -1   ! place in matrix
            idx = idx + 1

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

    ! I think we want to keep this without the subtraction since it 
    ! is used for the spacing of the arrays etc 
    max_nn1 = max_tl1*(max_tl1+1)/2


    LUT_ptr = c_loc(modeLUT)
    ierr = copy_LUT_array(LUT_ptr, total_nn1*4)






    ! Compute the size of all the strains:
    ! (sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, max_tl1, 6, nmodes)
    ! Now we have collapsed the matrix we can store only the number of tl1s 
    ! that each mode needs 

    ! FOR NOW WE ARE USING THE MORE MEMORY INEFFICIENT VERSION OF ASSUMING THEY ALL 
    ! HAVE max_tl1 - this makes the indexing a bit simpler in the kernel for now.
    ! probs need a LUT otherwise
    strainsize = nmodes * max_tl1 * sm%ngllx * sm%nglly * sm%ngllz * sm%nspec * 6 
    allocate(flatstrain_r(strainsize))
    allocate(flatstrain_i(strainsize))

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
                                    flatstrain_r(iii) =  real(allstrains(i, j, k, ispec, im, p, imode), kind=8)
                                    flatstrain_i(iii) = aimag(allstrains(i, j, k, ispec, im, p, imode))
                                else
                                    ! regions outside of the tl1 that is valid for this mode
                                    flatstrain_r(iii) = zero 
                                    flatstrain_i(iii) = zero 
                                endif 
                                iii = iii + 1
                            enddo 
                        enddo 
                    enddo 
                enddo
            enddo 
        enddo 
    enddo 

    strain_r_ptr = c_loc(flatstrain_r)
    strain_i_ptr = c_loc(flatstrain_i)

    ! Need to transfer the megastrains! 
    ierr = copy_allstrains(strain_r_ptr, strain_i_ptr, strainsize)


    ierr = allocate_Vani_arrays(nmodes*max_nn1)

    allocate(Vani_real(nmodes*max_nn1))
    Vani_real_ptr = c_loc(Vani_real)

    allocate(Vani_imag(nmodes*max_nn1))
    Vani_imag_ptr = c_loc(Vani_imag)



    ierr = allocate_eta_arrays(size_of_array)
    ierr = allocate_Cxyz_array(size_of_array*36)

    allocate(eta1(size_of_array))
    allocate(eta2(size_of_array))
    eta1_ptr = c_loc(eta1)
    eta2_ptr = c_loc(eta2)

    allocate(localcxyz(size_of_array*36))
    cxyz_ptr = c_loc(localcxyz)

    ! Assuming a maximum number of points of 150 ish 
    allocate(BrettModelToTransfer(MaxBrettModelPts, 5))
    allocate(flat3Dmodel(MaxBrettModelPts*5))


    if(myrank.eq.0)write(*,*)'------------------ BEGIN ALL THE LOOPS ------------------ '

    ! Loop over the models: 
    do imodel_iter = 10900, 10905 !nmodeliter
        call buffer_int(iterstr, imodel_iter)

        call system_clock(count_rate=count_rate)
        call system_clock(loop_clock_start)



        do model_chain = 1, 1
            call buffer_int(chainstr, model_chain)

            if(myrank.eq.0)then
                open(34,file='./output/timing/'//trim(timingNEX)//'it_'//trim(iterstr)//'_ch_'//trim(chainstr), form='formatted')
            endif 
          
            call system_clock(count_rate=count_rate)
            call system_clock(start_clock)    
            ! Load the model: 
            !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
            Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/MCMC_models/instances/c"//trim(chainstr)//"_m"//trim(iterstr)//".txt"
            call Model3D%read_model_from_file()
            ! Flatten 
            iii=1
            do i = 1, Model3D%npts
                flat3Dmodel(iii  ) = Model3D%xcoord(i)
                flat3Dmodel(iii+1) = Model3D%ycoord(i)
                flat3Dmodel(iii+2) = Model3D%zcoord(i)
                flat3Dmodel(iii+3) = Model3D%valspats(i,1)
                flat3Dmodel(iii+4) = Model3D%valspats(i,2)
                iii = iii+5
            enddo 
            call system_clock(end_clock)
            elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
            if(myrank.eq.0)write(34,*) 'Read model + flatten:', elapsed_time*1000, ' ms'
            if(myrank.eq.0)write(*,*) 'Read model + flatten:', elapsed_time*1000, ' ms'
                         

            call system_clock(count_rate=count_rate)
            call system_clock(start_clock)
            ptr_m3D = c_loc(flat3Dmodel)
            ierr = copythisarraytodevice(ptr_m3D, MaxBrettModelPts*5, 4)            
            ierr = cpp_project_eta_to_gll(Model3D%npts, sm%nspec, sm%ngllx, & 
                                          Model3D%valconsts/100.0d0) 
            call system_clock(end_clock)
            elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
            if(myrank.eq.0)write(34,*) 'Copy model to device:', elapsed_time*1000, ' ms'
            if(myrank.eq.0)write(*,*) 'Copy model to device:', elapsed_time*1000, ' ms'
                              
                                    
            call system_clock(count_rate=count_rate)
            call system_clock(start_clock)
            ierr = launch_vanikernel(sm%ngllx, sm%nspec, total_nn1, max_nn1,  max_tl1, nmodes)
            call system_clock(end_clock)
            elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
            if(myrank.eq.0)write(34,*) 'Kernel time         :', elapsed_time*1000, ' ms'
            if(myrank.eq.0)write(*,*) 'Kernel time         :', elapsed_time*1000, ' ms'




            call system_clock(count_rate=count_rate)
            call system_clock(start_clock)
            ierr = copyfromdevice(Vani_real_ptr, max_nn1*nmodes, 8)
            ierr = copyfromdevice(Vani_imag_ptr, max_nn1*nmodes, 9)
            call system_clock(end_clock)
            elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
            if(myrank.eq.0)write(34,*) 'Copy back           :', elapsed_time*1000, ' ms'
            if(myrank.eq.0)write(*,*) 'Copy back           :', elapsed_time*1000, ' ms'




            ! Overall the number of cst values we need to compute and 
            ! conduct reduction of is as follows: 
            ! Each mode of degree l = we need the s from 0 to 2l (inclusive)
            ! for each s there are s+1 values we need to transfer (only computing the negative)
            ! Due to hermitian nature 

            call system_clock(count_rate=count_rate)
            call system_clock(start_clock)

            icst = 1
            do imode = 1, nmodes
                l1  = modeLs(imode)
                n1  = modeNs(imode)
                this_tl1 = 2*l1 + 1
                thisnn1  = this_tl1*(this_tl1+1)/2 - l1*(l1+1)/2


                do iii = 1, thisnn1
                    ival = iii 
                    call find_row_col(ival, thisrow, thiscol, l1)
                    VaniAllModes(thisrow, thiscol, imode) = Vani_real( (imode-1)*max_nn1 + iii) + SPLINE_iONE*Vani_imag((imode-1)*max_nn1 + iii)

                    ! For the values that are in colums < l we can use the symmetries of the matrix to compute
                    ! The relationship between V_mm' and V_-m -m' is in D.153 - we can use this but this will give us 
                    ! the elements  in V (col < m=0) that are on the wrong side of the diagonal, so then we can apply
                    ! the fact it is hermitian - this is why we dont take the conjugate in the end.
                    if(thisrow > l1 )then 
                        VaniAllModes(this_tl1 - thiscol + 1, this_tl1 - thisrow + 1, imode) = VaniAllModes(thisrow, thiscol, imode) * (-one)**(thisrow + thiscol - two*real(l1 +1 ,kind=8) )
                    endif 
                enddo 

         
                call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
                allocate(cst(num_s, ncols))
                call Hcomplex_to_cst(VaniAllModes(1:this_tl1, 1:this_tl1, imode), l1, l1, cst, ncols, num_s, t1, t1, 2)
                
                if (2*l1.gt.compute_cst_smax)then 
                    thissmax = compute_cst_smax
                else 
                    thissmax = 2*l1 
                endif 


                do is = 1, thissmax-smin+1, 2
                    do it = 1, (smin+is-1) +1
                        allcsts_r(icst) =  real(cst(is,it))
                        allcsts_i(icst) =  aimag(cst(is,it))
                        icst = icst + 1
                    enddo 
                enddo 
                    
                deallocate(cst)
            enddo      

            ! Now each process has compute the csts we can reduce them 
            call MPI_Reduce(allcsts_r, allcsts_r_RED, ncstsvals, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            call MPI_Reduce(allcsts_i, allcsts_i_RED, ncstsvals, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)


            if(myrank.eq.0)then 
                icst = 1
                do imode = 1, nmodes
                    l1  = modeLs(imode)
                    n1  = modeNs(imode)

                    if (2*l1.gt.compute_cst_smax)then 
                        thissmax = compute_cst_smax
                    else 
                        thissmax = 2*l1 
                    endif 


                    call buffer_int(nstr, n1)
                    call buffer_int(lstr, l1)
                
                    out_name =  './output/NEX_'//trim(timingNEX)//'cst_'//trim(nstr)//t1//trim(lstr)//'_'//trim(iterstr)//'_'//trim(chainstr)//'.txt'
                    open(1,file=trim(out_name), form='formatted')
                    do s = 0, thissmax, 2
                        do it = 1, s+1
                            write(1,*) s, it-s-1, allcsts_r_RED(icst), allcsts_i_RED(icst)
                            icst = icst + 1
                        enddo !it
                    enddo ! s 
                    close(1) ! close file
                enddo ! loop over modes for writing 

            endif 
            call system_clock(end_clock)
            elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
            if(myrank.eq.0)write(34,*) 'CST computation     :', elapsed_time*1000, ' ms'
            if(myrank.eq.0)write(*,*) 'CST computation     :', elapsed_time*1000, ' ms'
                

            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

            deallocate(Model3D%idspats)
            deallocate(Model3D%valspats)
            deallocate(Model3D%idconsts)
            deallocate(Model3D%valconsts)
            deallocate(Model3D%xcoord)
            deallocate(Model3D%ycoord)
            deallocate(Model3D%zcoord)

        enddo ! chain




        call system_clock(end_clock)
        elapsed_time = real(end_clock - loop_clock_start, kind=8) / real(count_rate, kind=8)
        if(myrank.eq.0)then 
            write(34,*)  'Iteration time      :', elapsed_time*1000, ' ms'
            write(*,*)  'Iteration time      :', elapsed_time*1000, ' ms'
            write(*,*)  'Completed iteration : ', imodel_iter
            write(*,*)
        endif 
        
    enddo ! imodel_iter


    ! Cleanup memory 
    deallocate(VaniAllModes)


    call mpi_finalize(ierr)
end program optimised_vani