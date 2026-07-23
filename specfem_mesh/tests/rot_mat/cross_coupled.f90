
program test_W_cross_coupled
    use params, only: rho_spl, Wmat
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: delta_spline
    use modes, only: Mode, get_mode
    use mineos_model, only: mineos, mineos_ptr
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    integer :: i, j, k,  m1, m2, n1, n2, l1, l2, npoints, knot_lower, knot_upper, nval, lval, nstart, nstop, lstart, lfinish 
    character(len=80)  :: out_name
    character(len=1)  :: Tval, t1, t2
    real(kind=CUSTOM_REAL), allocatable :: r_lower, r_upper
    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:), W_a(:)
    real(SPLINE_REAL) :: kmkm2, kpkm2, kmkp2, Sl1m, Sl2m, mf 
    real(SPLINE_REAL) :: int_Ws, int_Wa, wcomnondim, bparam

    logical :: toroidal
    type(Mode)            :: mode_1, mode_2    
    type(InterpPiecewise) :: interp

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

    n1 = 2
    t1 = "S"
    l1 = 3

    n2 = 3
    t2 = "T"
    l2 = 2


    knot_lower = 1
    r_lower    = zero        
    knot_upper = mineos%disc(mineos%ndisc)
    r_upper    = mineos%rdisc(mineos%ndisc)
    npoints    = 1000*(knot_upper-knot_lower)

    ! Create interpolator with evenly spaced points in IC 
    interp = create_PieceInterp(npoints)
    interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 
    call interp%setup()
    call interp%create_interpolation_radial_map()


    ! We also need the density: 
    allocate(rho_spl(npoints))
    call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)


    allocate(W_s(npoints))
    allocate(W_a(npoints))

    mode_1 = get_mode(n1, t1, l1, mineos_ptr)
    mode_2 = get_mode(n2, t2, l2, mineos_ptr)

    ! Setup W matrix
    allocate(Wmat(mode_1%tl1, mode_2%tl1))

    ! Interpolate mode splines
    call interp%interpolate_mode_eigenfunctions(mode_1)
    call interp%interpolate_mode_eigenfunctions(mode_2)

    write(*,*)"Interpolated!"

    if (mode_1%t.eq.'S')then  
        mode_1%v_spl  = mode_1%v_spl  / mode_1%kf
        mode_1%dv_spl = mode_1%dv_spl / mode_1%kf
    elseif(mode_1%t.eq.'T')then 
        mode_1%w_spl  = mode_1%w_spl  / mode_1%kf
        mode_1%dw_spl = mode_1%dw_spl / mode_1%kf
    else
        write(*,*)'mode 1 needs to be S or T but is', mode_1%t
    endif


    if (mode_2%t.eq.'S')then  
        mode_2%v_spl  = mode_2%v_spl  / mode_2%kf
        mode_2%dv_spl = mode_2%dv_spl / mode_2%kf
    elseif(mode_2%t.eq.'T')then 
        mode_2%w_spl  = mode_2%w_spl  / mode_2%kf
        mode_2%dw_spl = mode_2%dw_spl / mode_2%kf
    else
        write(*,*)'mode 2 needs to be S or T but is', mode_2%t
    endif


    ! Compute Ws (D.70)
    if (mode_1%t.ne.mode_2%t)then 
        W_s = SPLINE_ZERO
        write(*,*)'Ws will be 0'

    elseif(mode_1%t.eq.'S' .and. mode_2%t.eq.'S')then 
        W_s = mode_1%v_spl * mode_2%v_spl & 
            + mode_1%u_spl * mode_2%v_spl & 
            + mode_2%u_spl * mode_1%v_spl 
        
    elseif(mode_1%t.eq.'T' .and. mode_2%t.eq.'T')then
        W_s = mode_1%w_spl * mode_2%w_spl
    else
        write(*,*)'Error in mode type', mode_1%t, mode_2%t
        stop
    endif 

    
    W_s = W_s * rho_spl * interp%radial * interp%radial


    ! Now we need to integrate for rho Ws r^2 
    int_Ws =  integrate_r_traps(interp%radial, W_s, npoints)


    ! Compute Wa (D.71)
    ! kmkm2 = (mode_1%kf * mode_1%kf) - (mode_2%kf * mode_2%kf) - two 
    ! kmkp2 = (mode_1%kf * mode_1%kf) - (mode_2%kf * mode_2%kf) + two
    ! kpkm2 = (mode_1%kf * mode_1%kf) + (mode_2%kf * mode_2%kf) - two

    kmkm2 = (mode_2%kf * mode_2%kf) - (mode_1%kf * mode_1%kf) + two 
    kmkp2 = (mode_1%kf * mode_1%kf) - (mode_2%kf * mode_2%kf) + two
    kpkm2 = (mode_1%kf * mode_1%kf) + (mode_2%kf * mode_2%kf) - two


    ! if(mode_1%t.eq.'T' .and. mode_2%t.eq.'S')then 
    !     W_a = (kmkp2 * mode_1%w_spl * mode_2%u_spl) - &
    !           (kpkm2 * mode_1%w_spl * mode_2%v_spl)
    ! elseif(mode_1%t.eq.'S' .and. mode_2%t.eq.'T')then
    !     W_a = (kmkm2 * mode_1%u_spl * mode_2%w_spl) + &
    !           (kpkm2 * mode_1%v_spl * mode_2%w_spl)
    ! else 
    !     W_a = SPLINE_ZERO
    ! endif 


    if(mode_1%t.eq.'T' .and. mode_2%t.eq.'S')then 
        W_a = half * mode_1%w_spl *  (kmkp2*mode_2%u_spl  - kpkm2*mode_2%v_spl )
    elseif(mode_1%t.eq.'S' .and. mode_2%t.eq.'T')then
        W_a = - half * mode_2%w_spl *  (kmkm2*mode_1%u_spl  - kpkm2*mode_1%v_spl )
    else 
        W_a = SPLINE_ZERO
    endif 


    W_a = W_a  * rho_spl * interp%radial * interp%radial

    ! Now we need to integrate for rho Ws r^2 
    int_Wa =  integrate_r_traps(interp%radial, W_a, npoints)

    ! Only non zero if m1 = m2 
    Wmat = SPLINE_iZERO
    do m1 = -mode_1%l, mode_1%l
        do m2 = -mode_2%l, mode_2%l

            if (m1.eq.m2) then 
                mf = real(m1, kind=SPLINE_REAL)

                ! First line of D.68 
                if (mode_1%l.eq.mode_2%l)then
                    Wmat(m1+mode_1%l+1, m2+mode_2%l+1) = Wmat(m1+mode_1%l+1, m2+mode_2%l+1) + mf * OMEGA * int_Ws
                endif 


                ! Sl1m = (((mode_1%lf +mf)*(mode_1%lf-mf))/((SPLINE_TWO*mode_1%lf + one)*(SPLINE_TWO*mode_1%lf - one)))**SPLINE_HALF
                ! Sl2m = (((mode_2%lf +mf)*(mode_2%lf-mf))/((SPLINE_TWO*mode_2%lf + one)*(SPLINE_TWO*mode_2%lf - one)))**SPLINE_HALF

                ! ! Second line of D.68 
                ! Wmat(m1+mode_1%l+1, m2+mode_2%l+1) = Wmat(m1+mode_1%l+1, m2+mode_2%l+1) - & 
                !                                     (SPLINE_iONE * OMEGA * int_Wa *  & 
                !                                     (delta_spline(mode_1%l, mode_2%l+1)*Sl1m +  & 
                !                                      delta_spline(mode_1%l, mode_2%l-1)*Sl2m))


                Sl1m = (((mode_2%lf +mf + one)*(mode_2%lf+one-mf))/((SPLINE_TWO*mode_2%lf + three)*(SPLINE_TWO*mode_2%lf + one)))**SPLINE_HALF

                Sl2m = (((mode_1%lf +mf + one)*(mode_1%lf+one-mf))/((SPLINE_TWO*mode_1%lf + three)*(SPLINE_TWO*mode_1%lf + one)))**SPLINE_HALF

                ! Second line of D.68 
                Wmat(m1+mode_1%l+1, m2+mode_2%l+1) = Wmat(m1+mode_1%l+1, m2+mode_2%l+1) + & 
                                                    (SPLINE_iONE * OMEGA * int_Wa *  & 
                                                    (delta_spline(mode_1%l, mode_2%l+1)*Sl1m +  & 
                                                     delta_spline(mode_2%l, mode_1%l+1)*Sl2m))

            endif 
        enddo 
    enddo 




    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'./rot_mat/CrossCoupled/IConly_', mode_1%n, mode_1%t, mode_1%l, '_', mode_2%n, mode_2%t, mode_2%l, '.txt'
    write(*,*)'Saved as '//trim(out_name)
    call save_W_matrix(mode_1%l, mode_2%l, trim(out_name))


    deallocate(Wmat)





end program test_W_cross_coupled