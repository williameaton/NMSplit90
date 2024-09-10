
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
                          cleanup_for_mode, map_local_global, map_complex_vector, delta_spline
    use gll, only: setup_gll
    use math, only: sqrtp
    implicit none
    include "constants.h"

    integer :: i, j, k, l1, n1, l2, n2, ispec, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, component
    character(len=1)   :: char, type_1, type_2
    character(len=10)  :: mode
    character(len=250) :: arg, ensight_nm
    character(len=30)  :: out_name
    character(len=4), parameter :: fmti1 = "(i1)"
    character(len=4), parameter :: fmti2 = "(i2)"
    character(len=4), parameter :: fmti3 = "(i3)"
    character(len=4) :: fmt_n, fmt_l
    real(SPLINE_REAL) :: wcom, qmod
    real(SPLINE_REAL), allocatable :: u1(:), v1(:), du1(:), dv1(:)
    real(SPLINE_REAL), allocatable :: u2(:), v2(:), du2(:), dv2(:)
    integer :: npoints, tl1, tl2, nproc
    real(kind=CUSTOM_REAL), allocatable :: radial_vals(:), r_lower, r_upper

    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:), W_a(:)
    complex(SPLINE_REAL), allocatable :: disp1_glob(:,:)
    real(SPLINE_REAL) :: total_integral,  precision, val, lf1, lf2, kf1, kf2, kmkm2, kpkm2, kmkp2, Sl1m, Sl2m, mf 
    real(SPLINE_REAL) :: int_Ws, int_Wa


    real(SPLINE_REAL), allocatable ::  u_spl_1(:),  v_spl_1(:), & 
                                      du_spl_1(:), dv_spl_1(:)
    real(SPLINE_REAL), allocatable ::  u_spl_2(:),  v_spl_2(:), & 
                                      du_spl_2(:), dv_spl_2(:)
    ! Timing: 
        integer :: start_clock, end_clock, count_rate
        real(8) :: elapsed_time

        
    ! Read mineos model 
    call process_mineos_model()
    allocate(u1(NL), v1(NL), du1(NL), dv1(NL))
    allocate(u2(NL), v2(NL), du2(NL), dv2(NL))

    region  = 3

    ! Now we can load the mode: 
    ! Choose a mode: 
    type_1 = 'S'
    l1      = 4
    n1      = 0

    type_2 = 'S'
    l2      = 4
    n2      = 0


    ! Setup W matrix
    tl1 = 2*l1 + 1
    tl2 = 2*l2 + 1
    call allocate_if_unallocated(tl1, tl2, Wmat)



    ! Start clock count
    call system_clock(count_rate=count_rate)
    call system_clock(start_clock)

    ! Values for the inner core
    knot_lower = 1
    r_lower    = zero        
    knot_upper =  disc(1)
    r_upper    = rdisc(1)
    npoints    = 10*(knot_upper-knot_lower)


    allocate(radial_vals(npoints))
    allocate(interp_map(npoints))
    allocate(W_a(npoints), W_s(npoints))

    ! Create radial array where eigenfunction will be interpolated
    do j = 1, npoints
        radial_vals(j) = (r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R
    enddo 

    call create_interpolation_radial_map(radial_vals, interp_map, npoints, knot_lower, knot_upper)

    ! Get the first mode: 
    write(*,*)'get 1st mode: ', n1, type_1, l1



    call get_mode(type_1, n1, l1, wcom, qmod, u1, du1, v1, dv1, .false.)


    allocate(u_spl_1(npoints), v_spl_1(npoints), du_spl_1(npoints), dv_spl_1(npoints))
    call interpolate_mode_eigenfunctions(type_1, u1, v1, du1, dv1, knot_lower, knot_upper, &  
                                        radial_vals, npoints, interp_map, & 
                                        u_spl_1, v_spl_1, du_spl_1, dv_spl_1)

    ! Second mode: 
    write(*,*)'get 2nd mode', n2, type_2, l2

    call get_mode(type_2, n2, l2, wcom, qmod, u2, du2, v2, dv2, .false.)
    allocate(u_spl_2(npoints), v_spl_2(npoints), du_spl_2(npoints), dv_spl_2(npoints))
    call interpolate_mode_eigenfunctions(type_2, u2, v2, du2, dv2, knot_lower, knot_upper, &  
                                        radial_vals, npoints, interp_map, & 
                                        u_spl_2, v_spl_2, du_spl_2, dv_spl_2)


    ! We also need the density
    call deallocate_if_allocated(rho_spl)
    allocate(rho_spl(npoints))
    call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), knot_lower, knot_upper, & 
                                    radial_vals, npoints, rho_spl, interp_map)


    lf1 = real(l1, kind=SPLINE_REAL)
    kf1 = sqrtp(lf1*(lf1+SPLINE_ONE))

    lf2 = real(l2, kind=SPLINE_REAL)
    kf2 = sqrtp(lf2*(lf2+SPLINE_ONE))


    if (type_1.eq.'S')then 
        v_spl_1  = v_spl_1/kf1
        dv_spl_1 = dv_spl_1/kf1
    elseif(type_1.eq.'T')then 
        u_spl_1  = u_spl_1/kf1
        du_spl_1 = du_spl_1/kf1
    else
        write(*,*)'type_1 needs to be S or T but is', type_1
    endif


    if (type_2.eq.'S')then 
        v_spl_2  = v_spl_2/kf2
        dv_spl_2 = dv_spl_2/kf2
    elseif(type_2.eq.'T')then 
        u_spl_2  = u_spl_2/kf2
        du_spl_2 = du_spl_2/kf2
    else
        write(*,*)'type_2 needs to be S or T but is', type_2
    endif


    ! Compute Ws (D.70)
    if (type_1.ne.type_2)then 
        W_s = SPLINE_ZERO
        write(*,*)'Ws will be 0'

    elseif(type_1.eq.'S' .and. type_2.eq.'S')then 
        W_s = (v_spl_1*v_spl_1 + u_spl_1*v_spl_2 + u_spl_2*v_spl_1  ) * rho_spl * radial_vals * radial_vals
    elseif(type_1.eq.'T' .and. type_2.eq.'T')then
        W_s = (u_spl_1*u_spl_2) * rho_spl * radial_vals * radial_vals
    else
        write(*,*)'Error in mode type', type_1, type_2
        stop
    endif 



    ! Now we need to integrate for rho Ws r^2 
    int_Ws =  integrate_r_traps(radial_vals, W_s, npoints)



    ! Compute Wa (D.71)

    kmkm2 = kf1*kf1 - kf2*kf2 - two 
    kmkp2 = kf1*kf1 - kf2*kf2 + two
    kpkm2 = kf1*kf1 + kf2*kf2 - two

    if(type_1.eq.'T' .and. type_2.eq.'S')then 
        W_a = (kmkp2 * u_spl_1 * u_spl_1) - &
              (kpkm2 * u_spl_1 * v_spl_1)
    elseif(type_1.eq.'S' .and. type_2.eq.'T')then
        W_a = (kmkm2 * u_spl_1 * u_spl_1) + &
              (kpkm2 * v_spl_1 * u_spl_2)
    else 
        W_a = SPLINE_ZERO
    endif 

    W_a = W_a * SPLINE_HALF * rho_spl * radial_vals * radial_vals

    ! Now we need to integrate for rho Ws r^2 
    int_Wa =  integrate_r_traps(radial_vals, W_a, npoints)


    write(*,*)'integral of W_s: ', int_Ws
    write(*,*)'integral of W_a: ', int_Wa



    ! Only non zero if m1 = m2 
    Wmat = SPLINE_iZERO
    do m1 = -l1, l1
        do m2 = -l2, l2

            if (m1.eq.m2) then 
                mf = real(m1, kind=SPLINE_REAL)

                ! First line of D.68 
                if (l1.eq.l2)then
                    Wmat(m1+l1+1, m2+l2+1) = Wmat(m1+l1+1, m2+l2+1) +   mf * OMEGA * int_Ws
                endif 


                Sl1m = (((lf1 +mf)*(lf1-mf))/((SPLINE_TWO*lf1 + one)*(SPLINE_TWO*lf1 - one)))**SPLINE_HALF
                Sl2m = (((lf2 +mf)*(lf2-mf))/((SPLINE_TWO*lf2 + one)*(SPLINE_TWO*lf2 - one)))**SPLINE_HALF


                ! Second line of D.68 
                Wmat(m1+l1+1, m2+l2+1) = Wmat(m1+l1+1, m2+l2+1) - & 
                                        (SPLINE_iONE * OMEGA * int_Wa *  & 
                                         (delta_spline(l1, l2+1)*Sl1m +  & 
                                          delta_spline(l1, l2-1)*Sl2m))
            endif 



        enddo 
    enddo 

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'semi_analytical_', n1, type_1, l1, '_', n2, type_2, l2, '.txt'
    call save_W_matrix(l1, l2, trim(out_name))



    ! Compute run time
    call system_clock(end_clock)
    ! Calculate the elapsed time in seconds
    elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
    ! Print the elapsed time
    write(*,*) 'Wall clock time taken for full SEM:', elapsed_time, 'seconds'


end program semi_analytical_W_matrix