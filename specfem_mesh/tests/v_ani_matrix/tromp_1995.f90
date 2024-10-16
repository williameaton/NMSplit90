
program tromp_1995
    use params, only: Vani, rho_spl, Arad, Crad, Lrad, Nrad, Frad, nprocs
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only:  delta_spline
    use V_ani, only: integrate_GNIr2, save_Vani_matrix
    use w3j, only: thrj
    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    integer :: i, j, k, l, m, n, n1, s, q, ispec, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, is
    character(len=1)   :: t1
    character(len=30)  :: out_name
    real(SPLINE_REAL)  :: sum
    integer :: npoints
    real(kind=CUSTOM_REAL) r_lower, r_upper

    real(SPLINE_REAL) :: val, lf, mf, kf
    integer, dimension(5), parameter :: I_n = (/ 5, 3, 3, 1, 1/)
    real(kind=CUSTOM_REAL) :: dA, dC, dL, dN, dF, thirty, twone

    type(Mode)            :: mode_1
    type(InterpPiecewise) :: interp

    thirty = three * ten 
    twone  = three * seven 
    dA =  0.4d0
    dC = -0.2d0
    dL =  0.3d0
    dN = -0.5d0
    dF =  0.1d0
        
    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

    ! Choose a mode: 
    mode_1 = get_mode(6, 'S', 10, mineos_ptr)

    allocate(Vani(mode_1%tl1, mode_1%tl1))
    Vani = SPLINE_iZERO

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


    if (mode_1%t.eq.'S')then 
        mode_1%v_spl  = mode_1%v_spl/mode_1%kf
        mode_1%dv_spl = mode_1%dv_spl/mode_1%kf
    elseif(mode_1%t.eq.'T')then 
        mode_1%w_spl  = mode_1%w_spl/mode_1%kf
        mode_1%dw_spl = mode_1%dw_spl/mode_1%kf
    else
        write(*,*)'type_1 needs to be S or T but is', mode_1%t
    endif
    
    ! Constant value over the radius
    allocate(Arad(npoints))
    allocate(Crad(npoints))
    allocate(Lrad(npoints))
    allocate(Nrad(npoints))
    allocate(Frad(npoints))
    Arad(:) = dA
    Crad(:) = dC
    Lrad(:) = dL
    Nrad(:) = dN
    Frad(:) = dF


    do m = -mode_1%l, mode_1%l
        mf = real(m, kind=CUSTOM_REAL)

        do s = 0, 4, 2 ! Loop with step of 2 from 0-4
            ! Compute sum_N sum_I int \Gamma_{NI} r^2 dr
            sum = zero
            do N = 0, 4
                do I = 1, I_n(N+1)
                    sum = sum + integrate_GNIr2(s, mode_1%l, N, I, mode_1%u_spl, & 
                                                mode_1%du_spl, & 
                                                mode_1%v_spl, mode_1%dv_spl, npoints, & 
                                                interp%radial, Arad, Crad, Lrad, &
                                                Nrad, Frad, mode_1%t)
                enddo
            enddo 

            ! Computing D.208 but not including the (2s + 1 / 4pi)^1/2 term 
            ! since that is already added in to the Gamma_NI via the gammaD1_coeff 
            ! function 
            Vani(m+mode_1%l+1,m+mode_1%l+1) = Vani(m+mode_1%l+1,m+mode_1%l+1) + & 
                                             (-SPLINE_ONE)**mf * (two*mode_1%lf + one) * thrj(mode_1%l, s, mode_1%l, -m, 0, m) * sum 
        enddo 
    enddo 

    write(out_name, '(a,i1,a,i1,a)')'./v_ani_matrix/radial_', mode_1%n, mode_1%t, mode_1%l, '.txt'
    call save_Vani_matrix(mode_1%l, out_name)


end program tromp_1995



