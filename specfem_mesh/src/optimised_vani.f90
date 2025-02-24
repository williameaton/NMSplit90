
program optimised_vani
    use params, only: Vani, verbose, myrank, MPI_SPLINE_COMPLEX, & 
                        MPI_SPLINE_REAL, MPI_CUSTOM_REAL, IIN, IOUT, glob_eta1,   &
                        glob_eta2,  nmodes, nprocs, all_warnings, datadir, max_tl1, Cxyz
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use v_ani, only: save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF, & 
                        compute_Vani_matrix, compute_vani_matrix_stored, & 
                        convert_imag_to_real, save_Vani_real_matrix
    use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst, write_cst_complex_to_file

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
                            compute_Vani_onemode_allstrains, check_cuda_device_allocations

    implicit none
    include "constants.h"
    include 'mpif.h'

    integer :: iset, i,j,k,ispec, l1, l2, n1, m1,m2, n2, ierr, & 
                tl2, h, b, cluster_size, sets_per_process, & 
                myset_start, myset_end, i_mode, smin, smax, num_s, & 
                ncols, this_tl1, imode, im, imodel_iter, igll, iproc
    character(len=2) nstr, lstr
    character(len=5) iterstr
    character(len=3) chainstr
    character(len=450) :: out_name
    real :: gb_per_set

    real(kind=SPLINE_REAL), allocatable :: Vani_real(:,:)
    complex(kind=SPLINE_REAL), allocatable :: Vani_modesum(:,:)
    complex(kind=SPLINE_REAL), allocatable :: cst(:,:)
    character(len=20) :: model_ti
    complex(kind=SPLINE_REAL), allocatable :: allstrains(:, :, :, :, :, :, :)


    real(kind=8), allocatable :: wglljac_loc(:, :)
    real(kind=8), allocatable :: Cxyz_loc(:, :, :, :)

    ! KD tree: 
    type(KdTree)           :: tree
    type(SetMesh)          :: sm  
    type(M3D) :: model3D


    ! Modes: 
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 5, 6, 8, 8, 9, 9, 11, 11, 11, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21, 21, 22, 23, 23, 25, 25, 27/)
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 2, 3, 1, 5, 3, 4,  4,  5,  6,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7,  8,  1,  4,  5,  1,  2,  2/)
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 5, 6, 8, 8, 9, 9, 11, 11, 11, 13, 13, 13, 14, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 20, 20, 21, 21/) 
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 2, 3, 1, 5, 3, 4,  4,  5,  6,  2,  3,  6,   4,  3,  4,  5,  6,  7,  1,  8,  3,  4,  6,  1,  5,  6,  7/)
    
    ! Deuss figure modes
    !integer, dimension(nmodes), parameter :: modeNs = (/2, 3, 3, 6, 8, 8, 9, 11, 11, 13, 13, 13,  16, 16, 17,  18, 18, 21, 21, 23 , 23 /) 
    !integer, dimension(nmodes), parameter :: modeLs = (/3, 1, 2, 3, 1, 5, 3,  4,  5,  1,  2,  3,   5,  7,  1,   3,  4,  6,  7,  4 ,  5 /)
    
    integer, dimension(nmodes), parameter :: modeNs = (/13/) 
    integer, dimension(nmodes), parameter :: modeLs = (/ 1/)


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
    allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
    ! Allocate Vani matrices
    allocate(Vani(max_tl1, max_tl1))
    if(myrank.eq.0)then 
        allocate(Vani_modesum(max_tl1, max_tl1))
        allocate(Vani_real(max_tl1, max_tl1))
    endif



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


    if(myrank.eq.0)write(*,*)'Loading mode strains'

    do imode = 1, nmodes
        n1       = modeNs(imode)
        l1       = modeLs(imode)
        this_tl1 = 2*l1 +1
        do im =  -l1, l1 
            call sm%load_mode_strain_binary(n1, t1, l1, im, &  
                                            allstrains(:,:,:,:,l1+im+1,:,imode))
        enddo 
    enddo  
    if(myrank.eq.0)write(*,*)'Done.'


    ! For some reason Myrank =6 seems to not print...is this an issue? 
    call check_cuda_device_allocations(nprocs)
    

    ! Copying this all to the device: 
    ierr=0
    allocate(d_allstrains_r(sm%ngllx*sm%nglly*sm%ngllz, sm%nspec, max_tl1, 6, nmodes), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating d_allstrains_r for proc ', myrank
        stop 
    endif 

    ierr=0
    allocate(d_allstrains_i(sm%ngllx*sm%nglly*sm%ngllz, sm%nspec, max_tl1, 6, nmodes), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating d_allstrains_i for proc ', myrank
        stop 
    endif 


    call copy_allstrain_to_device(allstrains, modeLs, sm%ngllx, sm%nspec)
    if(myrank.eq.0)write(*,*)'Copied over all the strains.'


    ! WGLL on device
    ierr=0
    allocate(d_wglljac_g(sm%ngllx*sm%nglly*sm%ngllz, sm%nspec), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating d_wglljac_g for proc ', myrank
        stop 
    endif 

    ! Allocate local copy to cast over to new format: 
    ierr=0
    allocate(wglljac_loc(sm%ngllx*sm%nglly*sm%ngllz, sm%nspec), stat=ierr)
    if(ierr.ne.0)then 
        write(*,*)'Error allocating wglljac_loc for proc ', myrank
        stop 
    endif 

    do ispec = 1, sm%nspec 
        igll = 1
        do i = 1, sm%ngllx
            do j = 1, sm%nglly 
                do k = 1, sm%ngllz
                    wglljac_loc(igll, ispec) = sm%wglljac(i,j,k,ispec)
                    igll = igll + 1
                enddo 
            enddo 
        enddo 
    enddo 
    ! Copy wglljac to the device:
    d_wglljac_g = wglljac_loc
    if(myrank.eq.0)write(*,*)'Copied wglljac to device.'



    ! Allocate the cxyz for the device: 
    allocate(d_Cxyz_g(sm%ngllx*sm%nglly*sm%ngllz, sm%nspec, 6, 6))
    allocate(Cxyz_loc(sm%ngllx*sm%nglly*sm%ngllz, sm%nspec, 6, 6))
    if(myrank.eq.0)write(*,*)'Allocated Cxyz.'


    ! Allocate the Vmatrices on the GPU: 
    allocate(d_vani_real_G(max_tl1, max_tl1))
    allocate(d_vani_imag_G(max_tl1, max_tl1))


    ! Loop over the models: 
    do imodel_iter = 10900, nmodeliter
        call buffer_int(iterstr, imodel_iter)

        do model_chain = 1, 20 
            call buffer_int(chainstr, model_chain)

            ! Load the model: 
            !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
            Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/MCMC_models/instances/c"//trim(chainstr)//"_m"//trim(iterstr)//".txt"
            call Model3D%read_model_from_file()
            call Model3D%create_KDtree()
            
            !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS//voronoi/instances/instance_"//trim(iterstr)
            !call Model3D%re_readmodel()

            vor_A = Model3D%valconsts(1)/100.0d0
            vor_C = Model3D%valconsts(2)/100.0d0
            vor_L = Model3D%valconsts(3)/100.0d0
            vor_N = Model3D%valconsts(4)/100.0d0
            vor_F = Model3D%valconsts(5)/100.0d0

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


            ! Copy the Cxyz over to GPU once per model: 
            do ispec = 1, sm%nspec 
                igll = 1
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly 
                        do k = 1, sm%ngllz
                            Cxyz_loc(igll, ispec, :, :) = Cxyz(i,j,k,ispec, :, :)
                            igll = igll + 1
                        enddo 
                    enddo 
                enddo 
            enddo 
            ! Copy wglljac to the device:


            d_Cxyz_g = Cxyz_loc                                  



            ! Loop for each mode to compute the splitting and the Cst value
            do i_mode = 1, nmodes

                n1       = modeNs(i_mode)
                l1       = modeLs(i_mode)
                this_tl1 = 2*l1 +1

                ! Reset the matrix:
                Vani         = SPLINE_iZERO
                if(myrank.eq.0)then
                    Vani_modesum = SPLINE_iZERO
                endif
                
                ! Compute the Vani matrix
                call compute_Vani_onemode_allstrains(i_mode, l1, sm%ngllx, sm%nspec)

                !Reduce the matrices across all of the MPI procs
                call MPI_Reduce(Vani, Vani_modesum, max_tl1**2, MPI_SPLINE_COMPLEX, &
                                MPI_SUM, 0, MPI_COMM_WORLD, ierr)

                ! USE Vani_modesum to output/compute CSTs
                if(myrank.eq.0)then 

                    call buffer_int(nstr, n1)
                    call buffer_int(lstr, l1)
                    Vani = Vani_modesum

                    if(force_VTI)then 
                        out_name =  './output/instance_matrices/vani_'//trim(nstr)// t1//trim(lstr)//'_'//trim(iterstr)//'_'//trim(chainstr)//'_VTI.txt'
                    else 
                        out_name =  './output/instance_matrices/vani_'//trim(nstr)// t1//trim(lstr)//'_'//trim(iterstr)//'_'//trim(chainstr)//'.txt'
                    endif 
                    call save_Vani_matrix(l1,l1, out_name)
        
                    !call convert_imag_to_real(l1, l1, Vani_modesum(1:this_tl1, 1:this_tl1), Vani_real(1:this_tl1, 1:this_tl1))
                    call get_Ssum_bounds(l1, l1, smin, smax, num_s, ncols)
                    allocate(cst(num_s, ncols))
                    call Hcomplex_to_cst(Vani(1:this_tl1, 1:this_tl1), l1, l1, cst, ncols, num_s, t1, t1, 2)
                    out_name = 'output/instance_csts/cst_'//trim(nstr)//trim(t1)//trim(lstr)//trim(model_ti)//'_'//trim(iterstr)//'_'//trim(chainstr)
                    call write_cst_complex_to_file(out_name, cst, ncols, num_s, smin, 2)
                    deallocate(cst)
                endif

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            enddo ! i_mode 


            deallocate(Model3D%idspats)
            deallocate(Model3D%valspats)
            deallocate(Model3D%idconsts)
            deallocate(Model3D%valconsts)
            deallocate(Model3D%xcoord)
            deallocate(Model3D%ycoord)
            deallocate(Model3D%zcoord)

        enddo ! chain

        if(myrank.eq.0)write(*,*)'Completed iteration: ', imodel_iter
    enddo ! imodel_iter


    ! Cleanup memory 
    deallocate(Vani)
    if(myrank.eq.0) deallocate(Vani_modesum, Vani_real)

    call mpi_finalize(ierr)
end program optimised_vani