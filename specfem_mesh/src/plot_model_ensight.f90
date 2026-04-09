program plot_vani_to_ensight
    ! program loads in a Vani matrix (complex) and outputs the splitting functions
        use params, only: Vani, nprocs, nmodes, glob_eta1, glob_eta2, Cxyz
        use v_ani, only: load_vani_from_file, convert_imag_to_real, save_Vani_real_matrix, compute_Cxyz_at_gll_constantACLNF
        use splitting_function, only: get_Ssum_bounds, Hreal_to_cst, write_cst_to_file, write_cst_complex_to_file
        use specfem_mesh,       only: SetMesh, create_SetMesh
        use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
        use modes,              only: get_mode, Mode 
        use mineos_model,       only: mineos, mineos_ptr
        use ylm_plm,            only: ylm_real    
        use voronoi, only: vor_x, vor_y, vor_z, & 
                       vor_A, vor_C, vor_L, vor_N, vor_F, &
                       load_voronoi_model, project_voroni_to_gll
        use model3d, only: M3D

        implicit none 
        include "constants.h"
        character(len=3) :: nstr, lstr
        integer   :: l1, n1 , tl1, smin, smax, num_s, ncols, nrows, is, it, & 
                     iproc, i, j, k, ispec, s, t, ilat, ilon, nlat, nlon, i_mode
        character :: t1
        character(len=20) :: model_ti
        character(len=250) :: modefile , out_name
        real(kind=CUSTOM_REAL)  :: dlat, dlon, theta, phi
        real(kind=SPLINE_REAL)  :: sum
        real(kind=SPLINE_REAL), allocatable :: Vani_real(:,:)
        real(kind=SPLINE_REAL), allocatable :: sig_st(:,:)
        complex(kind=SPLINE_REAL), allocatable :: cst_imag(:,:)
        real(kind=CUSTOM_REAL), allocatable :: sigma(:), ylm_global(:)
        type(SetMesh) :: sm

        real(kind=SPLINE_REAL), allocatable :: xvec(:), yvec(:), zvec(:)

        integer, dimension(9), parameter :: Vcont = (/1, 2, 3, 4, 4, 5, 5, 6, 6/)
        integer :: m1, p, q 
        type(Mode)    :: mode_1
        real(kind=CUSTOM_REAL) :: cont

    
        integer, parameter :: region = 3
        type(M3D) :: model3D



        call mineos%process_mineos_model(.false.)
        mineos_ptr => mineos


        ! Read 3D model and build K-d tree: 
        !Model3D%filename = "/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
        Model3D%filename = "/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/3D_MODELS/benchmarks/6node_model.txt"
        !Model3D%filename = "/scratch/gpfs/TROMP/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/MCMC_models/instances/c1_m5000.txt"
        call Model3D%read_model_from_file()
        call Model3D%create_KDtree()

        do iproc = 0, nprocs -1 
            sm = create_SetMesh(iproc, region)
            ! Read the mesh info and coordinates
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            call sm%setup_gll()
            !call sm%load_original_boundaries()

            ! needed for ensight geo file
            call sm%setup_global_coordinate_arrays(.false.)
            call sm%compute_rtp_from_xyz(.false.)

            call sm%compute_jacobian(.false.)
            call sm%compute_wglljac(.false.)
            call sm%get_unique_radii(.false.)

            !call sm%compute_elliptical_boundary_perturbation('./ellipticity/epsi_discontinuities.txt', mineos%ndisc)

            allocate(glob_eta1(sm%nglob))

            call create_ensight_file_prefix(iproc, sm%region)
            call create_proc_case_file()
            call create_proc_geo_file(sm, 1)

            !Test boundaries first: 
            !glob_eta1(:) = zero
            !do ispec = 1, sm%nspec2D_bottom
            !     do i = 1, sm%ngllx
            !         do j = 1, sm%nglly
            !             glob_eta1(sm%ibool(i,j,1, sm%ibelm_bottom(ispec))) = sm%delta_surf_bottom(i,j,ispec)
            !         enddo 
            !     enddo 
            ! enddo 
            ! do ispec = 1, sm%nspec2D_top
            !     do i = 1, sm%ngllx
            !         do j = 1, sm%nglly
            !             !write(*,*)sm%ibelm_top(ispec)
            !             glob_eta1(sm%ibool(i,j,5, sm%ibelm_top(ispec))) = sm%delta_surf_top(i,j,ispec)
            !         enddo 
            !     enddo 
            ! enddo 

            !write(*,*)minval(glob_eta1), maxval(glob_eta1)
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'ellipticity', 1)



            call Model3D%project_to_gll(sm, glob_eta1, id=1)
            call write_real_scalar_to_ensight(sm, glob_eta1, 'longitude', 1)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=2)
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'love_C', 2)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=3)
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'love_L', 3)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=4)
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'love_N', 4)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=5)
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'love_F', 5)


            call sm%cleanup()


            deallocate(glob_eta1)



            write(*,*)'Finished processor ', iproc
        enddo 


    
    end program