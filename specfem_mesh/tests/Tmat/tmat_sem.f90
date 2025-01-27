! Computes the Tmat using SEM method to compare against radial method:
program tmat_sem_test
    use params, only:  Tmat, rho_spl, nprocs
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use specfem_mesh, only: create_SetMesh, SetMesh
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use ylm_plm, only: ylm_complex
    use rho_st_profiles, only: single_r_rho_st
    implicit none
    include "constants.h"

    integer :: iproc, region, knot_lower, knot_upper, n1, n2, l1, l2, rho_smax, & 
               ispec, i,j,k, s, it , t
    character :: t1, t2
    type(SetMesh)         :: sm
    type(InterpPiecewise) :: interp
    character(len=30)     :: out_name
    complex(kind=CUSTOM_REAL), allocatable :: rho_st(:,:)
    real(kind=CUSTOM_REAL),    allocatable :: delta_rho(:,:,:,:), r_l, th_l, phi_l

    ! Setup the delta rho global perturbation
    rho_smax = 3
    allocate(rho_st(rho_smax, 2*rho_smax+1))
    
    rho_st(1, 1:3) = (/ 0.4d0,  -0.14d0,  -0.4d0 /)
    rho_st(2, 1:5) = (/ 0.7d0,  -0.43d0, -0.63d0,   0.43d0, 0.7d0/)
    rho_st(3, 1:7) = (/0.22d0,    0.1d0,   0.09d0,   0.0d0, -0.09d0, 0.1d0, -0.22d0/)
    

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
    allocate(Tmat(2*l1+1, 2*l2+1))
    Tmat = SPLINE_iZERO

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


        allocate(delta_rho(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
        delta_rho = zero 

        ! Compute delta rho based on the cst expansion:
        do ispec = 1, sm%nspec 
            do i = 1, sm%ngllx 
                do j = 1, sm%nglly 
                    do k = 1, sm%ngllz 

                        ! Get the r, theta, phi values: 
                        r_l   = sm%rstore(i,j,k,ispec)
                        th_l  = sm%thetastore(i,j,k,ispec)
                        phi_l = sm%phistore(i,j,k,ispec)


                        ! Sum over s and t: 
                        do s = 1, rho_smax           
                            do it = 1, 2*s + 1
                                t = it - 1 - s
                                ! compute the actual t value:                                 
                                delta_rho(i,j,k,ispec) = delta_rho(i,j,k,ispec)                + &
                                                        (single_r_rho_st(s, rho_st(s,it), r_l) * & 
                                                         ylm_complex(s, t, th_l, phi_l))
                            enddo ! it
                        enddo !is 
                    enddo 
                enddo 
            enddo 
        enddo


        call compute_T_matrix(sm, sm%interp, delta_rho, & 
                              n1, t1, l1,               & 
                              n2, t2, l2,               & 
                              .true.)
        call sm%cleanup()
        deallocate(delta_rho)
    enddo 

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'Tmat/Tmat_', n1, t1, l1, '_', n2, t2, l2, '.txt'
    call save_T_matrix(l1, l2, trim(out_name))







end program tmat_sem_test
