
program test_W_matrix
    use params, only: NL, IC_ID, ndisc, rdisc, disc, disp1, disp2, & 
                     rad_mineos, rho_mineos, u_spl, v_spl, interp_map, & 
                     interp_id_r, unique_r, n_unique_rad, ngllx, nglly, &
                     ngllz, nspec, Wmat, rad_id, wgll, detjac, rho_spl, & 
                      nglob, realfmt
    use spline, only: create_interpolation_radial_map, interpolate_mode_eigenfunctions, & 
                      write_mode_spline, interpolate_mineos_variable, quad_spline_interp_3, & 
                      write_scalar_spline
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: compute_jacobian, compute_rtp_from_xyz, load_ibool, & 
                          read_proc_coordinates, rotate_complex_vector_rtp_to_xyz, & 
                          compute_rotation_matrix, setup_global_coordinate_arrays, & 
                          cleanup_for_mode, map_local_global, map_complex_vector
    use gll, only: setup_gll
    use math, only: sqrtp
    implicit none
    include "constants.h"

    integer :: i, j, k, l, n , ispec, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, component
    character(len=1)   :: char, type
    character(len=10)  :: mode
    character(len=250) :: arg, ensight_nm
    character(len=30)  :: out_name
    character(len=4), parameter :: fmti1 = "(i1)"
    character(len=4), parameter :: fmti2 = "(i2)"
    character(len=4), parameter :: fmti3 = "(i3)"
    character(len=4) :: fmt_n, fmt_l
    real(SPLINE_REAL) :: wcom, qmod
    real(SPLINE_REAL), allocatable :: u(:), v(:), du(:), dv(:)
    integer :: npoints, tl1, nproc
    real(kind=CUSTOM_REAL), allocatable :: radial_vals(:), r_lower, r_upper

    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:)
    complex(SPLINE_REAL), allocatable :: disp1_glob(:,:)
    real(SPLINE_REAL) :: total_integral,  precision,  lf, kf
    complex(SPLINE_REAL) :: sum, val



    ! Read mineos model 
    call process_mineos_model()
    allocate(u(NL), v(NL), du(NL), dv(NL))


    ! Now we can load the mode: 
    ! Choose a mode: 
    type = 'S'
    l      = 9
    n      = 4
    region = 3
    nproc  = 6

    ! Setup W matrix
    tl1 = 2*l + 1
    ! The matrix should be 2l + 1 from -m to m 
    call deallocate_if_allocated(Wmat)
    call allocate_if_unallocated(tl1, tl1, Wmat)
    Wmat = SPLINE_iZERO


    do iproc = 0, nproc-1
        ! Things that need to be done for each processor
        call read_proc_coordinates(iproc, region)

        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian()

        call setup_global_coordinate_arrays()
        call compute_rtp_from_xyz()
        call get_mesh_radii()
        call compute_rotation_matrix()

        ! Get density at each radius 
        allocate(rho_spl(n_unique_rad))
        call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), 1, IC_ID, & 
                                        unique_r, n_unique_rad, rho_spl, interp_id_r)
        !call write_scalar_spline(unique_r, rho_spl, n_unique_rad, 'sem_spline_rho.txt')

        

        ! Load the eigenfunctions 
        call get_mode(type, n, l, wcom, qmod, u, du, v, dv, .true.)
        call interpolate_mode_eigenfunctions(type, u, v, du, dv, 1, IC_ID, &  
                                            unique_r, n_unique_rad, interp_id_r)


        !call setup_ensight_for_proc(iproc, region, 1)

        allocate(disp1(3, ngllx, nglly, ngllz, nspec) )
        allocate(disp2(3, ngllx, nglly, ngllz, nspec) )


        ! Since we are self coupling spheroidal currently we would only need the 
        ! diagonals 
        ! NOT TRUE FOR OTHER CASES
        do m1 = -l,  l
            m2 = m1

                call compute_global_mode_displacement(type, l, m1, disp1)
                call rotate_complex_vector_rtp_to_xyz(disp1)

                ! Save to ensight
                !if (m1.lt.0)then 
                !    write(ensight_nm, '(a,i2)')'disp1_', m1
                !else 
                !    write(ensigh t_nm, '(a,i1)')'disp1_', m1
                !endif 
                !allocate(disp1_glob(3, nglob))
                !call map_complex_vector(3, disp1, disp1_glob, 0)
                !call write_complex_vector_to_ensight(disp1_glob, trim(ensight_nm), 1)
                !deallocate(disp1_glob)


                !if(m1.eq.m2)then 
                disp2(:,:,:,:,:) = disp1(:,:,:,:,:)
                !else
                !    call compute_global_mode_displacement(type, l, m2, disp2)
                !    call rotate_complex_vector_rtp_to_xyz(disp2)
                !endif

                sum = SPLINE_iZERO
                ! Compute D.36 integral
                ! Note that \Omega should only be non-zero in the z-direction
                ! Hence the integrand should be 
                ! rho \tilde{s}^*  \dot (-i\Omega S_y, i\Omega s_x, 0)
                ! Note also that naturally the global mode displacement
                ! is computed in r, theta, phi and needs to be converted
                do ispec = 1, nspec 
                    do i = 1, ngllx
                        do j = 1, nglly
                            do k = 1, ngllz
                                val = rho_spl(rad_id(i,j,k,ispec)) * & 
                                (- conjg(disp1(1,i,j,k,ispec))*disp2(2,i,j,k,ispec) +  & 
                                   conjg(disp1(2,i,j,k,ispec))*disp2(1,i,j,k,ispec)) * & 
                                  detjac(i,j,k,ispec) * wgll(i) * wgll(j) * wgll(k)    
                                  
                                  sum = sum + val

                            enddo 
                        enddo
                    enddo
                enddo

                Wmat(m1+l+1, m2+l+1) = Wmat(m1+l+1, m2+l+1) + sum
        enddo ! m1 

        call cleanup_for_mode()

    enddo 


    ! Constants so multiply after
    Wmat = Wmat * OMEGA * iONE

    write(out_name, '(a,i1,a)')'Wmat_', l, '.txt'
    call save_W_matrix(l, trim(out_name))



end program test_W_matrix