
program radial_method
    use params, only: NL, IC_ID, rdisc, disc, disp1, disp2, & 
                     rad_mineos, rho_mineos, u_spl, v_spl, interp_map, & 
                     interp_id_r, unique_r, n_unique_rad, ngllx, nglly, &
                     ngllz, nspec, Vani, rad_id, wgll, detjac, rho_spl, & 
                     nglob, realfmt, u_spl, v_spl, udot_spl, vdot_spl, & 
                     Arad, Crad, Lrad, Nrad, Frad, vp, vp_spl, A0
    use spline, only: create_interpolation_radial_map, interpolate_mode_eigenfunctions,& 
     write_mode_spline, interpolate_mineos_variable, write_scalar_spline
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: compute_jacobian, compute_rtp_from_xyz, load_ibool, & 
                          read_proc_coordinates, rotate_complex_vector_rtp_to_xyz, & 
                          compute_rotation_matrix, setup_global_coordinate_arrays, & 
                          cleanup_for_mode, map_local_global, map_complex_vector, delta_spline
    use gll, only: setup_gll
    use math, only: sqrtp
    use V_ani, only: integrate_GNIr2, save_Vani_matrix
    use w3j, only: thrj
    implicit none
    include "constants.h"

    integer :: i, j, k, l, m, n,n1, s, q, ispec, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, component, is
    character(len=1)   :: char, type
    character(len=10)  :: mode
    character(len=250) :: arg, ensight_nm
    character(len=30)  :: out_name
    character(len=4), parameter :: fmti1 = "(i1)"
    character(len=4), parameter :: fmti2 = "(i2)"
    character(len=4), parameter :: fmti3 = "(i3)"
    character(len=4) :: fmt_n, fmt_l
    real(SPLINE_REAL) :: wcom, qmod, sum
    real(SPLINE_REAL), allocatable :: u1(:), v1(:), du1(:), dv1(:)
    integer :: npoints, tl, nproc
    real(kind=CUSTOM_REAL), allocatable :: radial_vals(:), r_lower, r_upper

    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:), W_a(:)
    complex(SPLINE_REAL), allocatable :: disp1_glob(:,:)
    real(SPLINE_REAL) :: total_integral,  precision, val, lf, mf, kf, kmkm2, kpkm2, kmkp2, Sl1m, Sl2m
    real(SPLINE_REAL) :: int_Ws, int_Wa

    integer, dimension(5), parameter :: I_n = (/ 5, 3, 3, 1, 1/)

    real(SPLINE_REAL), allocatable ::  u_spl_1(:),  v_spl_1(:), & 
                                      du_spl_1(:), dv_spl_1(:)

    real(kind=CUSTOM_REAL) :: dA, dC, dL, dN, dF, thirty, twone

    character(len=250) ::  aclnf_dir

    aclnf_dir = '/Users/eaton/Documents/Software/NMSplit90/& 
                specfem_mesh/tests/v_ani_matrix/Tromp_1993_model/ACLNF/radial'


    thirty = three * ten 
    twone =  three * seven 

    ! Choose a mode: 
    n1     = 9
    type   = 'S'
    l      = 3

    write(*,*)SCALE_T
        
    ! Read mineos model 
    call process_mineos_model()
    allocate(u1(NL), v1(NL), du1(NL), dv1(NL))


    ! Setup Vani matrix
    tl = 2*l + 1
    allocate(Vani(tl, tl))
    Vani = SPLINE_iZERO

    ! Values for the inner core
    knot_lower = 1
    r_lower    = zero        
    knot_upper = disc(1)
    r_upper    = rdisc(1)
    npoints    = 50*(knot_upper-knot_lower)

    allocate(radial_vals(npoints))
    allocate(interp_map(npoints))

    ! Create radial array where eigenfunction will be interpolated
    do j = 1, npoints
        radial_vals(j) = (r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R
    enddo 

    ! CAN USE THIS TO SAVE THE RADIUS
    !write(*,*)'Max radius: ', maxval(radial_vals)
    !open(1, file='ACLNF/radial/r_for_radmethod', form='formatted')
    !do i = 1, npoints
    !    write(1,*)radial_vals(i)
    !enddo 
    !close(1)
    
    call create_interpolation_radial_map(radial_vals, interp_map, npoints, knot_lower, knot_upper)

    ! Get the first mode: 
    call get_mode(type, n1, l, wcom, qmod, u1, du1, v1, dv1, .false.)

    write(*,*)'wcom = ', wcom

    allocate(u_spl_1(npoints), v_spl_1(npoints), du_spl_1(npoints), dv_spl_1(npoints))
    call interpolate_mode_eigenfunctions(type, u1, v1, du1, dv1, knot_lower, knot_upper, &  
                                        radial_vals, npoints, interp_map, & 
                                        u_spl_1, v_spl_1, du_spl_1, dv_spl_1)


    lf = real(l, kind=SPLINE_REAL)
    kf = (lf*(lf+SPLINE_ONE))**SPLINE_HALF

    if (type.eq.'S')then 
        v_spl_1  = v_spl_1/kf
        dv_spl_1 = dv_spl_1/kf
    elseif(type.eq.'T')then 
        u_spl_1  = u_spl_1/kf
        du_spl_1 = du_spl_1/kf
    else
        write(*,*)'type_1 needs to be S or T but is', type
    endif
    


    call load_ACLNF_from_files(aclnf_dir, 0, npoints)


    ! Interpolate the vp and rho to the relevant points 
    allocate(rho_spl(npoints))
    call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), knot_lower, knot_upper, & 
                                    radial_vals, npoints, rho_spl, interp_map)
    allocate(vp_spl(npoints))
    call interpolate_mineos_variable(real(vp, kind=SPLINE_REAL), knot_lower, knot_upper, & 
                                    radial_vals, npoints, vp_spl, interp_map)
    call deallocate_if_allocated(A0)
    allocate(A0(npoints))

    ! kappa + 4*mu 
    A0 = (vp_spl*vp_spl)*rho_spl

    Arad = Arad * A0 
    Crad = Crad * A0 
    Lrad = Lrad * A0 
    Nrad = Nrad * A0 
    Frad = Frad * A0 



    do m = -l, l
        mf = real(m, kind=CUSTOM_REAL)

        do s = 0, 4, 2 ! Loop with step of 2 from 0-4
            ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
            sum = zero
            do N = 0, 4
                do I = 1, I_n(N+1)
                    sum = sum + integrate_GNIr2(s, l, N, I, u_spl_1, du_spl_1, v_spl_1, dv_spl_1, npoints, radial_vals, Arad, Crad, Lrad, Nrad, Frad, type)
                enddo
            enddo 

            ! Computing D.208 but not including the (2s + 1 / 4pi)^1/2 term 
            ! since that is already added in to the Gamma_NI via the gammaD1_coeff 
            ! function 
            Vani(m+l+1,m+l+1) = Vani(m+l+1,m+l+1) + (-SPLINE_ONE)**mf * (two*lf + one) * thrj(l, s, l, -m, 0, m) * sum 

        enddo 
    enddo 

    write(out_name, '(a,i1,a,i1,a)')'./matrices/radial_', n1, type, l, '.txt'
    call save_Vani_matrix(l, out_name)

end program radial_method



