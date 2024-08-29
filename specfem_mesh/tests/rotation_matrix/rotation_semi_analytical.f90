
program semi_analytical_W_matrix
    use params, only: NL, IC_ID, ndisc, rdisc, disc, disp1, disp2, & 
                     rad_mineos, rho_mineos, u_spl, v_spl, interp_map, & 
                     interp_id_r, unique_r, n_unique_rad, ngllx, nglly, &
                     ngllz, nspec, Wmat, rad_id, wgll, detjac, rho_spl, & 
                     nglob, realfmt, u_spl, v_spl, udot_spl, vdot_spl
    use spline, only: create_interpolation_radial_map, interpolate_mode_eigenfunctions, write_mode_spline, interpolate_mineos_variable, write_scalar_spline
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
    real(SPLINE_REAL) :: total_integral,  precision, val, lf, kf
    complex(SPLINE_REAL) :: sum


    ! Read mineos model 
    call process_mineos_model()
    allocate(u(NL), v(NL), du(NL), dv(NL))


    ! Now we can load the mode: 
    ! Choose a mode: 
    type = 'S'
    l      = 9
    n      = 4
    region = 3

    ! Setup W matrix
    tl1 = 2*l + 1
    call allocate_if_unallocated(tl1, tl1, Wmat)



    ! Values for the inner core
    knot_lower = 1
    r_lower    = zero        
    knot_upper =  disc(1)
    r_upper    = rdisc(1)
    npoints    = 10*(knot_upper-knot_lower)


    allocate(radial_vals(npoints))
    allocate(interp_map(npoints), W_s(npoints))

    ! Create radial array where eigenfunction will be interpolated
    do j = 1, npoints
        radial_vals(j) = (r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R
    enddo 

    call get_mode(type, n, l, wcom, qmod, u, du, v, dv, .true.)


    ! Interpolate
    call create_interpolation_radial_map(radial_vals, interp_map, npoints, knot_lower, knot_upper)
    call interpolate_mode_eigenfunctions(type, u, v, du, dv, knot_lower, knot_upper, &  
                                        radial_vals, npoints, interp_map)
    call write_mode_spline(n, type, l, radial_vals, npoints)


    ! We also need the density
    call deallocate_if_allocated(rho_spl)
    allocate(rho_spl(npoints))
    call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), knot_lower, knot_upper, & 
                                    radial_vals, npoints, rho_spl, interp_map)
    call write_scalar_spline(radial_vals, rho_spl, npoints, 'spline_rho.txt')


    lf = real(l, kind=SPLINE_REAL)
    kf = sqrtp(lf*(lf+SPLINE_ONE))
    v_spl    = v_spl/kf
    vdot_spl = vdot_spl/kf



    ! Compute the integrand: 
    if(type.eq.'S')then 
        W_s = (v_spl*v_spl + SPLINE_TWO*u_spl*v_spl ) * rho_spl * radial_vals * radial_vals
    elseif(type.eq.'T')then
        W_s = (u_spl*u_spl) * rho_spl * radial_vals * radial_vals
    else
        write(*,*)'Error in mode type', type
        stop
    endif 


    ! Now we need to integrate 
    sum = zero
    sum =  integrate_r_traps(radial_vals, W_s, npoints)

    Wmat = SPLINE_iZERO
    do m1 = -l, l 
        Wmat(m1+l+1,m1+l+1) = real(m1, kind=SPLINE_REAL)*OMEGA*sum

        ! If you want it in the form D.180
        !Wmat(l+1+m1, l+1-m1) = spline_ione * real(m1, kind=SPLINE_REAL) * OMEGA * sum
    enddo 


    write(out_name, '(a,i1,a)')'semi_analytical', l, '.txt'
    call save_W_matrix(l, trim(out_name))


end program semi_analytical_W_matrix