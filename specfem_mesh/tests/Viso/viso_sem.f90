program viso
    ! Computes the Viso matrix using SEM method to compare against radial method:
    use params, only:  Viso, nprocs
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use specfem_mesh, only: create_SetMesh, SetMesh
    use rho_st_profiles, only: single_r_rho_st
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use ylm_plm, only: ylm_complex
    implicit none
    include "constants.h"

    integer :: iproc, region, knot_lower, knot_upper, n1, n2, l1, l2, kap_smax, & 
               ispec, i,j,k, s, it , t
    character :: t1, t2
    type(SetMesh)         :: sm
    type(InterpPiecewise) :: interp
    character(len=30)     :: out_name
    complex(kind=CUSTOM_REAL), allocatable :: kap_st(:,:), mu_st(:,:), rho_st(:,:)
    real(kind=CUSTOM_REAL),    allocatable :: delta_kap(:,:,:,:), & 
                                              delta_mu(:,:,:,:),  & 
                                              delta_rho(:,:,:,:)                  
    real(kind=CUSTOM_REAL) ::  r_l, th_l, phi_l
    complex(kind=CUSTOM_REAL):: l_ylm


    ! Setup the delta kappa global perturbation
    kap_smax = 3
    allocate(kap_st(kap_smax, 2*kap_smax+1))
    allocate(mu_st(kap_smax, 2*kap_smax+1))
    allocate(rho_st(kap_smax, 2*kap_smax+1))
    
    rho_st(1, 1:3) = (/0.4d0, -0.14d0,  -0.4d0 /)
    rho_st(2, 1:5) = (/ 0.7d0,  -0.43d0, -0.63d0,   0.43d0, 0.7d0/)
    rho_st(3, 1:7) = (/0.22d0,  0.1d0,   0.09d0,   0.0d0, -0.09d0, 0.1d0, -0.22d0/)
    
    kap_st(1, 1:3) = (/0.2d0,   0.55d0,  -0.2d0 /)
    kap_st(2, 1:5) = (/ 0.54d0,  0.0d0, -0.23d0,   0.0d0, 0.54d0/)
    kap_st(3, 1:7) = (/0.43d0,  0.34d0,   -0.09d0,  1.0d0, 0.09d0, 0.34d0, -0.43d0/)
    
    mu_st(1, 1:3) = (/ 0.6d0,  -0.2d0,  -0.6d0 /)
    mu_st(2, 1:5) = (/ 0.14d0,  0.13d0,  0.47d0,   -0.13d0, 0.14d0/)
    mu_st(3, 1:7) = (/0.52d0,  -0.1d0,   0.39d0,   0.2d0, -0.39d0, -0.1d0, -0.52d0/)
    

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos


    n1 = 0
    n2 = 0

    t1 = 'S'
    t2 = 'S'

    l1 = 2
    l2 = 2

    region = 3

    ! The matrix should be 2l + 1 from -m to m 
    allocate(Viso(2*l1+1, 2*l2+1))
    Viso = SPLINE_iZERO

    do iproc = 0, nprocs-1
        write(*,*)'iproc = ', iproc

        ! Things that need to be done for each processor
        sm = create_SetMesh(iproc, region)

        call sm%read_proc_coordinates()

        call sm%load_ibool()
        call sm%setup_gll()
        call sm%compute_jacobian(.true.)
        call sm%compute_rtp_from_xyz(.true.)
        call sm%compute_wglljac(.true.)

        call sm%setup_global_coordinate_arrays(.true.)
        call sm%get_unique_radii(.true.)
        call sm%compute_rotation_matrix()
        call sm%compute_background_g()


        allocate(delta_kap(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
        allocate(delta_mu( sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
        allocate(delta_rho(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))

        delta_kap = zero 
        delta_mu  = zero 
        delta_rho = zero 

        ! Compute delta kappa based on the cst expansion:
        do ispec = 1, sm%nspec 
            do i = 1, sm%ngllx 
                do j = 1, sm%nglly 
                    do k = 1, sm%ngllz 

                        ! Get the r, theta, phi values: 
                        r_l   = sm%rstore(i,j,k,ispec)
                        th_l  = sm%thetastore(i,j,k,ispec)
                        phi_l = sm%phistore(i,j,k,ispec)

                        ! Sum over s and t: 
                        do s = 1, kap_smax           
                            do it = 1, 2*s + 1
                                t = it - 1 - s

                                l_ylm = ylm_complex(s, t, th_l, phi_l)
                                
                                ! Compute the delta kappa value at each GLL 
                                delta_kap(i,j,k,ispec) = delta_kap(i,j,k,ispec)                + &
                                                        (single_r_rho_st(s, kap_st(s,it), r_l) * l_ylm) 

                                delta_mu(i,j,k,ispec) = delta_mu(i,j,k,ispec)                + &
                                                        (single_r_rho_st(s, mu_st(s,it), r_l) * l_ylm) 

                                delta_rho(i,j,k,ispec) = delta_rho(i,j,k,ispec)                + &
                                                        (single_r_rho_st(s, rho_st(s,it), r_l) * l_ylm)              
                            enddo ! it
                        enddo !is 

                    enddo 
                enddo 
            enddo 
        enddo

        call compute_Viso_matrix(sm, sm%interp, delta_kap, &
                                 delta_mu, delta_rho,      & 
                                 n1, t1, l1,               & 
                                 n2, t2, l2,               & 
                                 .true.)
        call sm%cleanup()
        deallocate(delta_kap)
        deallocate(delta_mu)
        deallocate(delta_rho)
    enddo 

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'Viso/Viso_', n1, t1, l1, '_', n2, t2, l2, '.txt'
    call save_Viso_matrix(l1, l2, trim(out_name))

end program viso