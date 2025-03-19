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
        Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/voronoi_model_new_format.txt"
        !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/benchmarks/DR_benchmark_model.txt"
        !Model3D%filename = "/scratch/gpfs/we3822/NMSplit90/specfem_mesh/3D_MODELS/voronoi/MCMC_models/instances/c1_m5000.txt"
        call Model3D%read_model_from_file()
        call Model3D%create_KDtree()

        do iproc = 0, nprocs -1 
            sm = create_SetMesh(iproc, region)
            ! Read the mesh info and coordinates
            call sm%read_proc_coordinates()
            call sm%load_ibool()
            call sm%setup_gll()

            ! needed for ensight geo file
            call sm%setup_global_coordinate_arrays(.false.)
            call sm%compute_rtp_from_xyz(.false.)

            call sm%compute_jacobian(.false.)
            call sm%compute_wglljac(.false.)
            call sm%get_unique_radii(.false.)

            

            allocate(glob_eta1(sm%nglob), glob_eta2(sm%nglob))
            call Model3D%project_to_gll(sm, glob_eta1, id=1)
            call Model3D%project_to_gll(sm, glob_eta2, id=2) ! colatitude

            allocate(xvec(sm%nglob), yvec(sm%nglob), zvec(sm%nglob))
            xvec = sin(glob_eta2) * cos(glob_eta1)
            yvec = sin(glob_eta2) * sin(glob_eta1)
            zvec = cos(glob_eta2)

            call create_ensight_file_prefix(iproc, 3)
            call create_proc_case_file()

            call create_proc_geo_file(sm, 1)


            call write_real_scalar_to_ensight(sm, 180.d0*glob_eta2/PI, 'colatitude', 1)

            call write_real_vector_to_ensight(sm, xvec, yvec, zvec, 'TTI', 1)


            ! coordinate longitude; 
            do ispec = 1, sm%nspec
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        do k = 1, sm%ngllz
                            xvec(sm%ibool(i,j,k,ispec)) = sm%phistore(i,j,k,ispec)*180.d0/PI 
                            if( xvec(sm%ibool(i,j,k,ispec)).gt.180.d0)then 
                                xvec(sm%ibool(i,j,k,ispec)) = xvec(sm%ibool(i,j,k,ispec)) - 360.d0
                            endif 
                        enddo 
                    enddo 
                enddo 
            enddo 
            call write_real_scalar_to_ensight(sm, xvec, 'coordlong', 1)


                        ! coordinate longitude; 
            xvec = zero 
            do ispec = 1, sm%nspec
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        do k = 1, sm%ngllz
                            if(ispec.eq.2201)then
                                xvec(sm%ibool(i,j,k,ispec)) = 99.0d0
                            endif 
                            if(ispec.eq.2202)then
                                xvec(sm%ibool(i,j,k,ispec)) = 99.0d0
                            endif 
                            if(ispec.eq.2221)then
                                xvec(sm%ibool(i,j,k,ispec)) = 99.0d0
                            endif 
                            if(ispec.eq.2222)then
                                xvec(sm%ibool(i,j,k,ispec)) = 99.0d0
                            endif 
                        enddo 
                    enddo 
                enddo 
            enddo 
            call write_real_scalar_to_ensight(sm, xvec, 'spec', 1)




            ! For an example mode, compute the contraction of Eps C Eps 
            n1      =  27
            t1      = 'S'
            l1      =  6
            m1      =  4

            vor_A = Model3D%valconsts(1)/100.0d0
            vor_C = Model3D%valconsts(2)/100.0d0
            vor_L = Model3D%valconsts(3)/100.0d0
            vor_N = Model3D%valconsts(4)/100.0d0
            vor_F = Model3D%valconsts(5)/100.0d0


            mode_1  = get_mode(n1, t1, l1, mineos_ptr)

            call sm%compute_rotation_matrix()
            call compute_Cxyz_at_gll_constantACLNF(sm, vor_A, vor_C, vor_L, & 
                                                    vor_N, vor_F, glob_eta1, glob_eta2, &
                                                    perturbation_on_prem=.true.)
            call sm%interp%interpolate_mode_eigenfunctions(mode_1)
            call deallocate_if_allocated(sm%strain1)
            allocate(sm%strain1(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec, 6))

            call sm%compute_mode_strain(m1, mode_1, sm%strain1)

            call sm%rotate_complex_sym_matrix_rtp_to_xyz(sm%strain1)

            xvec = zero 
            do ispec = 1, sm%nspec 
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        do k = 1, sm%ngllz
                            cont = SPLINE_iZERO
                            do p = 1, 9
                                do q = 1, 9
                                    cont = cont + ( conjg(sm%strain1(i,j,k,ispec,Vcont(p))) * & 
                                                    Cxyz(i,j,k,ispec,Vcont(p),Vcont(q))     * & 
                                                    sm%strain1(i,j,k,ispec,Vcont(q)) ) 
                                enddo
                            enddo 
                            xvec(sm%ibool(i,j,k,ispec)) = cont 
                        enddo 
                    enddo
                enddo
            enddo
            call write_real_scalar_to_ensight(sm, xvec, 'contraction', 1)

            deallocate(sm%strain1)


            ! allocate(glob_eta1(sm%nglob))
            ! call Model3D%project_to_gll(sm, glob_eta1, id=1) ! A 
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'A', 1)
            
            ! call Model3D%project_to_gll(sm, glob_eta1, id=2) ! C
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'C', 1)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=3) ! L
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'L', 1)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=4) ! N
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'N', 1)

            ! call Model3D%project_to_gll(sm, glob_eta1, id=5) ! N
            ! call write_real_scalar_to_ensight(sm, glob_eta1, 'F', 1)



            call sm%cleanup()


            deallocate(glob_eta1)
            deallocate(glob_eta2)
            deallocate(xvec, yvec, zvec)
    


            write(*,*)'Finished processor ', iproc
        enddo 


    
    end program