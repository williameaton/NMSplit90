! Benchmark for a VTI hemisphere 
program hemisphere
    use params, only: Vani, rho_spl, vp_spl, Arad, Crad, Lrad, Nrad, Frad, nprocs, nmodes
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only:  delta_spline
    use V_ani, only: integrate_gnir2_mochizuki, save_Vani_matrix, integrate_gnir2
    use w3j, only: thrj
    use mineos_model, only: mineos, mineos_ptr
    use ylm_plm, only: XNlm 
    use modes, only: Mode, get_mode
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    type(Mode)            :: mode_1
    type(InterpPiecewise) :: interp

    ! Expanded I_n 
    integer, dimension(5), parameter :: I_nsimple = (/ 5, 3, 3, 1, 1/)
    integer, dimension(9), parameter :: I_n = (/1, 1, 3, 3, 5, 3, 3, 1, 1/)

    real(kind=CUSTOM_REAL) :: Alove, Clove, Llove, Nlove, Flove, phi1, phi2, theta(101), out(101), xval, pol, mf, mf1, mf2 ,& 
                              r_lower, r_upper
    integer :: i,j,k, N, s, m, t, nl, ll, m1, m2, knot_lower, knot_upper, npoints, tl1
    complex(kind=CUSTOM_REAL) :: sum
    ! Set a constant ACLNF values 
    Alove =  0.04d0
    Clove = -0.02d0
    Llove =  0.03d0
    Nlove = -0.05d0
    Flove =  0.01d0

    phi1 = 0.0
    phi2 = two*PI!PI/three


    ! Test the values for XNlm: 

    !XNlm(theta, N, l, m)

    ! ! We want to evaluate the 21 radial integrals that come from the 
    ! ! more generalised Mochizuchi equation 

    ! ! We shouldnt need to go higher than this because the original L does not
    ! ! have anything higher than s = 0, 2, or 

    ! ! SELF COUPLING so l = l' 

    nl = 3
    ll = 1

    tl1 = 2*ll + 1

    allocate(Vani(tl1, tl1))

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

    mode_1 = get_mode(nl, 'S', ll, mineos_ptr)



    ! Values for the inner core
    knot_lower = 1
    r_lower    = zero        
    knot_upper = mineos%disc(2)
    r_upper    = mineos%rdisc(2)
    npoints    = 50*(knot_upper-knot_lower)
    interp = create_PieceInterp(npoints)
    interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 

    call interp%setup()
    call interp%create_interpolation_radial_map()

    ! Interpolate the mode splines
    call interp%interpolate_mode_eigenfunctions(mode_1)

        
    ! Constant value over the radius
    allocate(Arad(npoints))
    allocate(Crad(npoints))
    allocate(Lrad(npoints))
    allocate(Nrad(npoints))
    allocate(Frad(npoints))

    Arad(:) = Alove
    Crad(:) = Clove
    Lrad(:) = Llove
    Nrad(:) = Nlove
    Frad(:) = Flove


    ! Original method
        Vani(:,:) = SPLINE_iZERO

    do m = -mode_1%l, mode_1%l
        mf1 = real(m, kind=CUSTOM_REAL)

        do s = 0, 4, 2 ! Loop with step of 2 from 0-4
            ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
            sum = zero
            do N = 0, 4
                do I = 1, I_nsimple(N+1)
                    sum = sum + integrate_GNIr2(s, mode_1%l, N, I, mode_1%u_spl, & 
                                                mode_1%du_spl, & 
                                                mode_1%v_spl/mode_1%kf, mode_1%dv_spl/mode_1%kf, npoints, & 
                                                interp%radial, Arad, Crad, Lrad, &
                                                Nrad, Frad, mode_1%t)
                enddo
            enddo 

            write(*,*)'sum ', sum
            write(*,*)'Adds', (-SPLINE_ONE)**mf1 * (two*mode_1%lf + one) * thrj(mode_1%l, s, mode_1%l, -m, 0, m) * sum 
            write(*,*)

            ! Computing D.208 but not including the (2s + 1 / 4pi)^1/2 term 
            ! since that is already added in to the Gamma_NI via the gammaD1_coeff 
            ! function 
            Vani(m+mode_1%l+1,m+mode_1%l+1) = Vani(m+mode_1%l+1,m+mode_1%l+1) + & 
                                            (-SPLINE_ONE)**mf1 * (two*mode_1%lf + one) * thrj(mode_1%l, s, mode_1%l, -m, 0, m) * sum 
        enddo 
    enddo 

    call save_Vani_matrix(ll, ll, "hemisphere1.txt")

    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)



    ! General method 
    write(*,*)"Starting integration"

    Vani(:,:) = SPLINE_iZERO


    do m1 = -mode_1%l, mode_1%l
        mf1 = real(m1, kind=CUSTOM_REAL)
        m2 = m1 
        !do m2 = m1, mode_1%l 
            mf2 = real(m2, kind=CUSTOM_REAL)
            write(*,*)mf1, mf2

            do s = 0, 4, 2 
                t = 0
                !do t = -s, s
                    write(*,*)s, t
                    ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
                    sum = SPLINE_iZERO
                    do N = -4, 4
                        do I = 1, I_n(N+5)
                            sum = sum + integrate_GNIr2_Mochizuki(s, t, mode_1%l, N, I, mode_1%u_spl, & 
                                                                    mode_1%du_spl, & 
                                                                    mode_1%v_spl/mode_1%kf, mode_1%dv_spl/mode_1%kf, npoints, & 
                                                                    interp%radial, Arad, Crad, Lrad, &
                                                                    Nrad, Frad, mode_1%t, phi1, phi2)                        
                        enddo
                    enddo 

                    write(*,*)'sum ', sum
                        write(*,*)'Adds', (-SPLINE_ONE)**mf1 * (two*mode_1%lf + one) * thrj(mode_1%l, s, mode_1%l, -m1, t, m2) * sum 
                        write(*,*)

                    ! Computing D.208 but not including the (2s + 1 / 4pi)^1/2 term 
                    ! since that is already added in to the Gamma_NI via the gammaD1_coeff 
                    ! function 
                    Vani(m1+mode_1%l+1, m2+mode_1%l+1) = Vani(m1+mode_1%l+1,m2+mode_1%l+1) + & 
                                            (-SPLINE_ONE)**mf1 * (two*mode_1%lf + one) * thrj(mode_1%l, s, mode_1%l, -m1, t, m2) * sum 
                !enddo ! t
            enddo ! s

        !enddo! m2
    enddo ! m1 


    call save_Vani_matrix(ll, ll, "hemisphere2.txt")


end program hemisphere