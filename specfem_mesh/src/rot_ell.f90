! Computes the effect of rotation and ellipticity for mode: 


program semi_analytical_W_matrix
    use params, only: rho_spl, Wmat, Tell, Vell, Vcen, mu_spl, kappa_spl
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: delta_spline
    use modes, only: Mode, get_mode
    use w3j, only: thrj
    use mineos_model, only: mineos, mineos_ptr
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use woodhouse_kernels, only: Slm,Rlm, WK_TbarSrho, WK_TcaronSrho, WK_VbarSk, &
    WK_VbarSmu, WK_VbarSrho, WK_VcaronSk, WK_VcaronSmu, WK_VcaronSrho, WK_Vphi, WK_Vphi_dot
    implicit none
    include "constants.h"

    integer :: i, j, k,  m1, m2, npoints, knot_lower, knot_upper 
    character(len=80)  :: out_name
    real(kind=CUSTOM_REAL), allocatable :: r_lower, r_upper
    real(SPLINE_REAL), allocatable :: integrand(:), gravacc(:)
    complex(SPLINE_REAL), allocatable :: W_s(:), W_a(:), TbarSrho(:), TcaronSrho(:), & 
                                         Tellintegrand(:), VbarSk(:), VbarSmu(:), VbarSrho(:), & 
                                         VcaronSk(:), VcaronSmu(:), VcaronSrho(:), Vellintegrand(:), VS2dotphi(:), VS2phi(:), VS2integrand(:)
    real(kind=SPLINE_REAL), allocatable :: eta(:), epsi(:)
    real(SPLINE_REAL) :: kmkm2, kpkm2, kmkp2, Sl1m, Sl2m, mf 
    real(SPLINE_REAL) :: int_Ws, int_Wa

    complex(SPLINE_REAL) :: int_Tell, int_Vell, int_2ell

    type(Mode)            :: mode_1, mode_2    
    type(InterpPiecewise) :: interp

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

    ! ONLY FOR SELF COUPLING CURRENTLY
    mode_1 = get_mode(3, 'S', 2 , mineos_ptr)
    mode_2 = get_mode(3, 'S', 2 , mineos_ptr)


    ! Setup W matrix
    allocate(Wmat(mode_1%tl1, mode_2%tl1))
    allocate(Vell(mode_1%tl1, mode_2%tl1))
    allocate(Tell(mode_1%tl1, mode_2%tl1))
    allocate(Vcen(mode_1%tl1, mode_2%tl1))

    ! Values for the inner core
    knot_lower = 1
    r_lower    = zero        
    knot_upper = mineos%disc(mineos%ndisc)
    r_upper    = mineos%rdisc(mineos%ndisc)
    npoints    = 10*(knot_upper-knot_lower)

    write(*,*)
    write(*,*)
    write(*,*)

    
    ! Create interpolator with evenly spaced points in IC 
    interp = create_PieceInterp(npoints)
    interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 
    call interp%setup()
    call interp%create_interpolation_radial_map()

    ! Interpolate mode splines
    call interp%interpolate_mode_eigenfunctions(mode_1)
    call interp%interpolate_mode_eigenfunctions(mode_2)

    ! We also need the density: 
    allocate(rho_spl(npoints))
    allocate(kappa_spl(npoints))
    allocate(mu_spl(npoints))
    call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)


    ! Compute the kappa and mu from rho, vp, vs 
    allocate(mineos%mu_mineos(mineos%NR))
    allocate(mineos%kappa_mineos(mineos%NR))

    mineos%mu_mineos    = mineos%rho_mineos * (mineos%vs_mineos)**(two) 
    mineos%kappa_mineos = mineos%rho_mineos * (mineos%vp_mineos)**(two)  - four * mineos%mu_mineos / three

    call interp%interpolate_mineos_variable(real(mineos%mu_mineos, kind=SPLINE_REAL), mu_spl)
    call interp%interpolate_mineos_variable(real(mineos%kappa_mineos, kind=SPLINE_REAL), kappa_spl)


    open(1,file='ellipticity/radialpoints.txt', form='formatted')
    do i = 1, npoints
        write(1,*)interp%radial(i)
    enddo 
    close(1)


    open(1,file='ellipticity/rhopoints.txt', form='formatted')
    do i = 1, npoints
        write(1,*)rho_spl(i)*RHOAV
    enddo 
    close(1)

    ! ------------------- COMPUTE ROTATION TERMS   -------------------
    ! Compute Ws (D.70)
    allocate(W_s(npoints))
    if (mode_1%t.ne.mode_2%t)then 
        W_s = SPLINE_ZERO
        write(*,*)'Ws will be 0'

    elseif(mode_1%t.eq.'S' .and. mode_2%t.eq.'S')then 
        W_s = mode_1%v_spl/mode_1%kf * mode_2%v_spl/mode_2%kf & 
            + mode_1%u_spl * mode_2%v_spl/mode_2%kf & 
            + mode_2%u_spl * mode_1%v_spl/mode_1%kf 
        
    elseif(mode_1%t.eq.'T' .and. mode_2%t.eq.'T')then
        W_s = mode_1%w_spl/mode_1%kf * mode_2%w_spl/mode_2%kf
    else
        write(*,*)'Error in mode type', mode_1%t, mode_2%t
        stop
    endif 

    W_s = W_s * rho_spl * interp%radial * interp%radial

    ! Now we need to integrate for rho Ws r^2 
    int_Ws =  integrate_r_traps(interp%radial, W_s, npoints)

    ! Compute Wa (D.71)
    allocate(W_a(npoints))

    kmkm2 = (mode_1%kf * mode_1%kf) - (mode_2%kf * mode_2%kf) - two 
    kmkp2 = (mode_1%kf * mode_1%kf) - (mode_2%kf * mode_2%kf) + two
    kpkm2 = (mode_1%kf * mode_1%kf) + (mode_2%kf * mode_2%kf) - two

    if(mode_1%t.eq.'T' .and. mode_2%t.eq.'S')then 
        W_a = (kmkp2 * mode_1%w_spl/mode_1%kf * mode_2%u_spl) - &
              (kpkm2 * mode_1%w_spl/mode_1%kf * mode_2%v_spl/mode_2%kf)
    elseif(mode_1%t.eq.'S' .and. mode_2%t.eq.'T')then
        W_a = (kmkm2 * mode_1%u_spl * mode_2%w_spl/mode_2%kf) + &
              (kpkm2 * mode_1%v_spl/mode_1%kf * mode_2%w_spl/mode_2%kf)
    else 
        W_a = SPLINE_ZERO
    endif 

    W_a = W_a * SPLINE_HALF * rho_spl * interp%radial * interp%radial

    ! Now we need to integrate for rho Ws r^2 
    int_Wa =  integrate_r_traps(interp%radial, W_a, npoints)

    ! ------------------- COMPUTE ELLIPTICITY TERMS   -------------------

    ! Load the interpolated eta and epsilon values:
    allocate(eta(npoints))
    allocate(epsi(npoints))
    allocate(gravacc(npoints))

    open(1, file = 'ellipticity/eta_interpolated.txt', status = 'old', form='formatted')
    do i = 1, npoints
        read(1,*)eta(i)
    enddo 
    close(1)
    open(1, file = 'ellipticity/epsilon_interpolated.txt', status = 'old', form='formatted')
    do i = 1, npoints
        read(1,*)epsi(i)
    enddo 
    close(1)

    open(1, file = 'ellipticity/gravity.txt', status = 'old', form='formatted')
    do i = 1, npoints
        read(1,*)gravacc(i)
    enddo 
    close(1)
    gravacc = gravacc/(SCALE_V/SCALE_T)   ! non dimensionalise


    ! Compute Tell integral (D.80)
    allocate(Tellintegrand(npoints))
    allocate(TbarSrho(npoints))
    allocate(TcaronSrho(npoints))

    call WK_TbarSrho(mode_1, mode_2,   TbarSrho)
    call WK_TcaronSrho(mode_1, mode_2, TcaronSrho)

    Tellintegrand = (two/three) * epsi * rho_spl  * interp%radial  * interp%radial * & 
                    (TbarSrho - ((eta+three)*TcaronSrho))


    allocate(Vellintegrand(npoints))
    allocate(VbarSk(npoints))
    allocate(VbarSmu(npoints))
    allocate(VbarSrho(npoints))
    allocate(VcaronSk(npoints))
    allocate(VcaronSmu(npoints))
    allocate(VcaronSrho(npoints))


    call WK_VbarSk(mode_1, mode_2, interp%radial, VbarSk)   
    call WK_VbarSmu(mode_1, mode_2, interp%radial, VbarSmu)
    call WK_VbarSrho(mode_1, mode_2, interp%radial, gravacc, rho_spl, VbarSrho)

    call WK_VcaronSk(mode_1, mode_2, interp%radial, VcaronSk)
    call WK_VcaronSmu(mode_1, mode_2, interp%radial, VcaronSmu)
    call WK_VcaronSrho(mode_1, mode_2, interp%radial, gravacc, rho_spl, VcaronSrho)

   

    Vellintegrand = (two/three) * epsi * (  kappa_spl * (VbarSk   - (eta+one)*VcaronSk  )  & 
                                          + mu_spl    * (VbarSmu  - (eta+one)*VcaronSmu )  & 
                                          + rho_spl   * (VbarSrho - (eta+three)*VcaronSrho)) &
                    * interp%radial * interp%radial  




    allocate(VS2phi(npoints))
    allocate(VS2dotphi(npoints))

    call WK_Vphi(mode_1, mode_2, 2, interp%radial, rho_spl,VS2phi)               
    call WK_Vphi_dot(mode_1, mode_2, 2, interp%radial, rho_spl, VS2dotphi)

    ! Vs=2 integrand
    VS2integrand = (interp%radial**three) * (OMEGA**two) * (two*VS2dotphi  + VS2phi*interp%radial )/three 

    ! We only compute the contribtions for the self coupling case here 
    int_Tell =  integrate_r_traps(interp%radial, Tellintegrand, npoints)
    int_Vell =  integrate_r_traps(interp%radial, Vellintegrand, npoints)
    int_2ell =  integrate_r_traps(interp%radial, VS2integrand,  npoints)


    write(*,*) int_Tell
    write(*,*) int_Vell
    write(*,*) int_2ell


    ! Only non zero if m1 = m2 
    Wmat = SPLINE_iZERO
    Vell = SPLINE_iZERO
    Tell = SPLINE_iZERO
    Vcen = SPLINE_iZERO
    do m1 = -mode_1%l, mode_1%l
        do m2 = -mode_2%l, mode_2%l

            if (m1.eq.m2) then 
                mf = real(m1, kind=SPLINE_REAL)

                ! First line of D.68 
                if (mode_1%l.eq.mode_2%l)then
                    Wmat(m1+mode_1%l+1, m2+mode_2%l+1) = Wmat(m1+mode_1%l+1, m2+mode_2%l+1) + mf * OMEGA * int_Ws
                endif 

                Sl1m = Slm(mode_1%l, m1)
                Sl2m = Slm(mode_2%l, m1)

                ! Second line of D.68 
                Wmat(m1+mode_1%l+1, m2+mode_2%l+1) = Wmat(m1+mode_1%l+1, m2+mode_2%l+1) - & 
                                                     (SPLINE_iONE * OMEGA * int_Wa *  & 
                                                     (delta_spline(mode_1%l, mode_2%l+1)*Sl1m +  & 
                                                      delta_spline(mode_1%l, mode_2%l-1)*Sl2m))

                Vell(m1+mode_1%l+1, m2+mode_2%l+1) = Vell(m1+mode_1%l+1, m2+mode_2%l+1) + & 
                                                     Rlm(mode_1%l, m1)*int_Vell         - &
                                                     ((-one)**mf  *                       &
                                                     ((two*mode_1%lf + one)*(two*mode_2%lf + one))**half * &
                                                     thrj(mode_1%l, 2, mode_2%l, -m1, 0, m1) * int_2ell)


                ! Tell contribution for Self coupling: 
                if(mode_1%l.eq.mode_2%l)then 
                    Tell(m1+mode_1%l+1, m2+mode_2%l+1) = Tell(m1+mode_1%l+1, m2+mode_2%l+1) + & 
                                                         Rlm(mode_1%l, m1)*int_Tell

                  
                    if(mode_1%t.eq.mode_2%t .and. mode_1%n.eq.mode_2%n)then 
                        Vcen(m1+mode_1%l+1, m2+mode_2%l+1) = Vcen(m1+mode_1%l+1, m2+mode_2%l+1) &
                                                           + (two/three)*OMEGA*OMEGA
                    endif 
                    Vcen(m1+mode_1%l+1, m2+mode_2%l+1) = Vcen(m1+mode_1%l+1, m2+mode_2%l+1) & 
                                                       - (two/three)*mode_1%kf*mode_1%kf*OMEGA*OMEGA*int_Ws
                endif 

                Vcen(m1+mode_1%l+1, m2+mode_2%l+1) = Vcen(m1+mode_1%l+1, m2+mode_2%l+1) & 
                                                + ((-one)**mf  *                       &
                                                ((two*mode_1%lf + one)*(two*mode_2%lf + one))**half * &
                                                thrj(mode_1%l, 2, mode_2%l, -m1, 0, m1) * int_2ell)

            endif 
        enddo 
    enddo 

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'./ellipticity/matrices/Wmat_', mode_1%n, mode_1%t, mode_1%l, '_', mode_2%n, mode_2%t, mode_2%l, '.txt'

    call save_W_matrix(mode_1%l, mode_2%l, trim(out_name))

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'./ellipticity/matrices/Vell_', mode_1%n, mode_1%t, mode_1%l, '_', mode_2%n, mode_2%t, mode_2%l, '.txt'
    call save_Vell_matrix(mode_1%l, mode_2%l, out_name)

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'./ellipticity/matrices/Tell_', mode_1%n, mode_1%t, mode_1%l, '_', mode_2%n, mode_2%t, mode_2%l, '.txt'
    call save_Tell_matrix(mode_1%l, mode_2%l, out_name)

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'./ellipticity/matrices/Vcen_', mode_1%n, mode_1%t, mode_1%l, '_', mode_2%n, mode_2%t, mode_2%l, '.txt'
    call save_Vcen_matrix(mode_1%l, mode_2%l, out_name)


end program semi_analytical_W_matrix