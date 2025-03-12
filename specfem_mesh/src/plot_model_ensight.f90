program plot_vani_to_ensight
    ! program loads in a Vani matrix (complex) and outputs the splitting functions
        use params, only: Vani, nprocs, nmodes, glob_eta1, glob_eta2
        use v_ani, only: load_vani_from_file, convert_imag_to_real, save_Vani_real_matrix
        use splitting_function, only: get_Ssum_bounds, Hreal_to_cst, write_cst_to_file, write_cst_complex_to_file
        use specfem_mesh,       only: SetMesh, create_SetMesh
        use modes,              only: get_mode, Mode 
        use mineos_model,       only: mineos, mineos_ptr
        use ylm_plm,            only: ylm_real    
    
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

        integer, parameter :: region = 3
        type(M3D) :: model3D



        call mineos%process_mineos_model(.false.)
        

        ! Read 3D model and build K-d tree: 
        !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
        Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/benchmarks/DR_benchmark_model.txt"
        !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/MCMC_models/instances/c1_m5000.txt"
        call Model3D%read_model_from_file()
        call Model3D%create_KDtree()

        do iproc = 0, nprocs -1 
            sm = create_SetMesh(iproc, region)
            ! Read the mesh info and coordinates
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            ! needed for ensight geo file
            call sm%setup_global_coordinate_arrays(.false.)
            call sm%compute_rtp_from_xyz(.false.)

            
            !allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
            !call Model3D%project_to_gll(sm, glob_eta1, id=1)
            !call Model3D%project_to_gll(sm, glob_eta2, id=2)


            call create_ensight_file_prefix(iproc, 3)
            call create_proc_case_file()

            call create_proc_geo_file(sm, 1)


            allocate(glob_eta1(sm%nglob))
            call Model3D%project_to_gll(sm, glob_eta1, id=1) ! A 
            call write_real_scalar_to_ensight(sm, glob_eta1, 'A', 1)
            
            call Model3D%project_to_gll(sm, glob_eta1, id=2) ! C
            call write_real_scalar_to_ensight(sm, glob_eta1, 'C', 1)

            call Model3D%project_to_gll(sm, glob_eta1, id=3) ! L
            call write_real_scalar_to_ensight(sm, glob_eta1, 'L', 1)

            call Model3D%project_to_gll(sm, glob_eta1, id=4) ! N
            call write_real_scalar_to_ensight(sm, glob_eta1, 'N', 1)

            call Model3D%project_to_gll(sm, glob_eta1, id=5) ! N
            call write_real_scalar_to_ensight(sm, glob_eta1, 'F', 1)

            call sm%cleanup()


            deallocate(glob_eta1)
            !deallocate(glob_eta2)


            write(*,*)'Finished processor ', iproc
        enddo 


    
    end program