! Computes the effect of rotation and ellipticity for mode: 
program semi_analytical_W_matrix
    use params, only: rho_spl, Wmat, Tell, Vell, Vcen, mu_spl, kappa_spl, nmodes
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: delta_spline
    use modes, only: Mode, get_mode
    use w3j, only: thrj
    use mineos_model, only: mineos, mineos_ptr
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use woodhouse_kernels, only: Slm,Rlm, WK_TbarSrho, WK_TcaronSrho, WK_VbarSk, &
    WK_VbarSmu, WK_VbarSrho, WK_VcaronSk, WK_VcaronSmu, WK_VcaronSrho, WK_Vphi, WK_Vphi_dot

    use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst_8, write_cst_complex_to_file

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
    real(SPLINE_REAL) :: int_Ws, int_Wa, PiG, tau, nu,wcomnondim, aparam, bparam, cparam, ksq

    complex(SPLINE_REAL) :: int_Tell, int_Vell, int_2ell
    character(len=2)    :: nstr ,lstr
    type(Mode)            :: mode_1, mode_2    
    type(InterpPiecewise) :: interp

    real(kind=CUSTOM_REAL) :: tau_nu_pref

    character(len=1) :: t1

    ! cst
    complex(kind=SPLINE_REAL), allocatable :: Ellmat(:,:)
    complex(kind=SPLINE_REAL), allocatable :: cst(:,:)
    integer :: smin, smax, num_s, ncols

    ! 19 rn
    !integer, dimension(19), parameter :: modeNs = (/16, 2, 3, 6, 8, 9, 11, 11, 13, 13, 14, 15, 16, 18, 20, 21, 23, 25, 27/)
    !integer, dimension(19), parameter :: modeLs = (/ 5, 3, 2, 3, 5, 3, 4, 5, 2, 3, 4, 3, 6, 4, 5, 6, 5, 2, 2/)
    integer, dimension(31), parameter :: modeNs = (/0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 6, 8, 8, 9, 11, 11, 13, 13, 18, 18 /)
    integer, dimension(31), parameter :: modeLs = (/2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 1, 2, 3, 1, 5, 3,  4,  5,  1,  2,  3,  4 /)
    integer :: imode


    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos


    knot_lower = 1
    r_lower    = zero        
    knot_upper = mineos%disc(mineos%ndisc)
    r_upper    = mineos%rdisc(mineos%ndisc)
    npoints    = 10*(knot_upper-knot_lower)

    ! Create interpolator with evenly spaced points in IC 
    interp = create_PieceInterp(npoints)
    interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 
    call interp%setup()
    call interp%create_interpolation_radial_map()


    allocate(W_s(npoints))
    allocate(W_a(npoints))
    
    allocate(eta(npoints))
    allocate(epsi(npoints))
    allocate(gravacc(npoints))


    allocate(Tellintegrand(npoints))
    allocate(TbarSrho(npoints))
    allocate(TcaronSrho(npoints))
    allocate(Vellintegrand(npoints))
    allocate(VbarSk(npoints))
    allocate(VbarSmu(npoints))
    allocate(VbarSrho(npoints))
    allocate(VcaronSk(npoints))
    allocate(VcaronSmu(npoints))
    allocate(VcaronSrho(npoints))



    ! Load the interpolated eta and epsilon values:
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
    !gravacc = gravacc/(SCALE_V/SCALE_T)   ! non dimensionalise
    ! mineos has pi*g = 1 as normalisation, so im assuming that g at surface
    ! needs to be 1 / pi
    gravacc = gravacc/(gravacc(npoints) * PI)   ! non dimensionalise



    ! We also need the density: 
    allocate(rho_spl(npoints))
    allocate(kappa_spl(npoints))
    allocate(mu_spl(npoints))
    call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)


    ! Compute the kappa and mu from rho, vp, vs 
    allocate(mineos%mu_mineos(mineos%NR))
    allocate(mineos%kappa_mineos(mineos%NR))

    mineos%mu_mineos    = mineos%rho_mineos*(mineos%vs_mineos**two)
    mineos%kappa_mineos = mineos%rho_mineos*(mineos%vp_mineos**two)  - (four*mineos%mu_mineos/three)

    call interp%interpolate_mineos_variable(real(mineos%mu_mineos, kind=SPLINE_REAL), mu_spl)
    call interp%interpolate_mineos_variable(real(mineos%kappa_mineos, kind=SPLINE_REAL), kappa_spl)



    do imode = 1, 7
        write(*,*)'Mode ', modeNs(imode), ' S ', modeLs(imode)

        t1 = 'S'
        ! ONLY FOR SELF COUPLING CURRENTLY
        mode_1 = get_mode(modeNs(imode), t1, modeLs(imode) , mineos_ptr)
        mode_2 = get_mode(modeNs(imode), t1, modeLs(imode) , mineos_ptr)

        ! Setup matrix
        allocate(Wmat(mode_1%tl1, mode_2%tl1))
        allocate(Vell(mode_1%tl1, mode_2%tl1))
        allocate(Tell(mode_1%tl1, mode_2%tl1))
        allocate(Vcen(mode_1%tl1, mode_2%tl1))
        allocate(Ellmat(mode_1%tl1, mode_2%tl1))

        ! Interpolate mode splines
        call interp%interpolate_mode_eigenfunctions(mode_1)
        call interp%interpolate_mode_eigenfunctions(mode_2)


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


        open(1,file='ellipticity/eigen_u.txt', form='formatted')
        do i = 1, npoints
            write(1,*)mode_1%u_spl(i)
        enddo 
        close(1)

        open(1,file='ellipticity/eigen_v.txt', form='formatted')
        do i = 1, npoints
            write(1,*)mode_1%v_spl(i)
        enddo 
        close(1)

        ksq = mode_1%kf*mode_1%kf


        ! ------------------- COMPUTE ROTATION TERMS   -------------------
        ! Compute Ws (D.70)
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

        ! Compute Tell integral (D.80)
        !call WK_TbarSrho(mode_1, mode_2,   TbarSrho)
        !call WK_TcaronSrho(mode_1, mode_2, TcaronSrho)

        if(mode_1%t.eq.'S')then 
            TbarSrho    = -six * mode_1%u_spl * mode_1%v_spl/mode_1%kf
            TcaronSrho  = mode_1%u_spl*mode_1%u_spl + (mode_1%kf*mode_1%kf - three)*(mode_1%v_spl*mode_1%v_spl)/(ksq)
        else 
            TbarSrho = zero
            TcaronSrho =  (mode_1%kf*mode_1%kf - three) * (mode_1%w_spl*mode_1%w_spl)/(mode_1%kf*mode_1%kf)
        endif 


        Tellintegrand = (two/three) * epsi * rho_spl  * interp%radial  * interp%radial * & 
                        (TbarSrho - ((eta+three)*TcaronSrho))        


        ! call WK_VbarSk(mode_1, mode_2, interp%radial, VbarSk)   
        ! call WK_VbarSmu(mode_1, mode_2, interp%radial, VbarSmu)
        ! call WK_VbarSrho(mode_1, mode_2, interp%radial, gravacc, rho_spl, VbarSrho)

        ! call WK_VcaronSk(mode_1, mode_2, interp%radial, VcaronSk)
        ! call WK_VcaronSmu(mode_1, mode_2, interp%radial, VcaronSmu)
        ! call WK_VcaronSrho(mode_1, mode_2, interp%radial, gravacc, rho_spl, VcaronSrho)

        PiG = PI*GRAV

        if(mode_1%t.eq.'S')then 
            ! correct DT98 version
            VbarSk = -two * (mode_1%du_spl + mode_1%aux_f)*(mode_1%du_spl + three*mode_1%v_spl/(mode_1%kf*interp%radial))   ! D.184
            ! incorrect woodhuse 80 version
            ! VbarSk = -two * (mode_1%du_spl + mode_1%aux_f)*(mode_1%du_spl + three*mode_1%v_spl/(mode_1%kf*interp%radial))   ! D.184
        else 
            VbarSk = zero
        endif 


        if(mode_1%t.eq.'S')then 
        VbarSmu = -(two/three)*(two*mode_1%du_spl - mode_1%aux_f)*(two*mode_1%du_spl + nine*mode_1%dv_spl/mode_1%kf - two*six*mode_1%v_spl/(mode_1%kf * interp%radial) ) & 
                + two*mode_1%aux_x*( three*mode_1%du_spl - (ksq - three)*mode_1%dv_spl/mode_1%kf - three*mode_1%kf*mode_1%v_spl/interp%radial ) & 
                + nine*two*(ksq - two)*(mode_1%dv_spl*mode_1%v_spl)/(interp%radial*ksq)
        else 
        VbarSmu = nine*two*(mode_1%kf*mode_1%kf - two)*( mode_1%dw_spl*mode_1%w_spl)/(interp%radial*mode_1%kf*mode_1%kf) & 
                - two * (mode_1%kf*mode_1%kf - three) * mode_1%aux_z*mode_1%dw_spl/mode_1%kf
        endif 


        if(mode_1%t.eq.'S')then   
        VbarSrho = two*mode_1%aux_f*(interp%radial*mode_1%dp_spl + four*PiG*rho_spl*interp%radial*mode_1%u_spl + gravacc*mode_1%u_spl) & 
                    - six*gravacc*mode_1%u_spl*mode_1%v_spl/(mode_1%kf*interp%radial)  & 
                    + six*gravacc*mode_1%u_spl*mode_1%u_spl/interp%radial & 
                    + two*mode_1%p_spl * ((ksq - three)*mode_1%v_spl/mode_1%kf  - ksq*mode_1%u_spl)/interp%radial  
        else 
        VbarSrho = zero 
        endif 


        if(mode_1%t.eq.'S')then 
            VcaronSk = -(mode_1%du_spl + mode_1%aux_f)*(mode_1%du_spl - mode_1%aux_f - six*mode_1%v_spl/(mode_1%kf*interp%radial ))
        else 
            VcaronSk = zero
        endif 



        if(mode_1%t.eq.'S')then 
            VcaronSmu = (ksq - six*two)*(ksq - two)*(mode_1%v_spl*mode_1%v_spl)/(ksq *interp%radial * interp%radial) & 
                    + (ksq - three)*(mode_1%aux_x*mode_1%aux_x - two*mode_1%aux_x*mode_1%dv_spl/mode_1%kf)    & 
                    - (two/three)*(two*mode_1%du_spl - mode_1%aux_f)*(mode_1%du_spl + half*mode_1%aux_f - six*mode_1%v_spl/(mode_1%kf * interp%radial ) )
        else 
            VcaronSmu = (mode_1%kf*mode_1%kf - six*two)*(mode_1%kf*mode_1%kf - two)*(mode_1%w_spl*mode_1%w_spl)/(mode_1%kf*mode_1%kf *interp%radial * interp%radial) & 
                    +(mode_1%kf*mode_1%kf - three)*(mode_1%aux_z*mode_1%aux_z - two*mode_1%aux_z*mode_1%dw_spl/mode_1%kf)   
        endif 


        if(mode_1%t.eq.'S')then 
            VcaronSrho = two*(ksq - three)*mode_1%p_spl*mode_1%v_spl/(mode_1%kf*interp%radial) + & 
                        mode_1%u_spl*(two*mode_1%dp_spl + eight*PiG*rho_spl*mode_1%u_spl - six*gravacc*mode_1%v_spl/(mode_1%kf*interp%radial))
        else 
            VcaronSrho = zero
        endif 


        VbarSk(1)     = zero 
        VbarSmu(1)    = zero 
        VbarSrho(1)   = zero 
        VcaronSk(1)   = zero 
        VcaronSmu(1)  = zero 
        VcaronSrho(1) = zero 



        Vellintegrand = (two/three) * epsi * interp%radial * interp%radial * &
                        (  kappa_spl * (VbarSk   - (eta+one)*VcaronSk  )    & 
                        +  mu_spl    * (VbarSmu  - (eta+one)*VcaronSmu )    & 
                        +  rho_spl   * (VbarSrho - (eta+three)*VcaronSrho)) 
                        



        !allocate(VS2phi(npoints))
        !allocate(VS2dotphi(npoints))

        !call WK_Vphi(mode_1, mode_2, 2, interp%radial, rho_spl,VS2phi)               
        !call WK_Vphi_dot(mode_1, mode_2, 2, interp%radial, rho_spl, VS2dotphi)

        ! Vs=2 integrand
        !VS2integrand = (interp%radial**three) * (OMEGA**two) * (two*VS2dotphi  + VS2phi*interp%radial )/three 

        ! We only compute the contribtions for the self coupling case here 
        int_Tell =  integrate_r_traps(interp%radial, Tellintegrand, npoints)
        int_Vell =  integrate_r_traps(interp%radial, Vellintegrand, npoints)
        !int_2ell =  integrate_r_traps(interp%radial, VS2integrand,  npoints)

        tau_nu_pref = (ksq )/( (mode_1%lf*two + three )*(two*mode_1%lf  - one) )




        tau = int_Tell * tau_nu_pref
        nu  = int_Vell * tau_nu_pref


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
                                                        (one- three*m1*m1/(mode_2%kf*mode_1%kf))*nu      



                    Tell(m1+mode_1%l+1, m2+mode_2%l+1) = Tell(m1+mode_1%l+1, m2+mode_2%l+1) + & 
                                                            (one- three*m1*m1/(mode_2%kf*mode_1%kf))*tau

                    
                    Vcen(m1+mode_1%l+1, m2+mode_2%l+1) = Vcen(m1+mode_1%l+1, m2+mode_2%l+1) &
                                                        + (two/three)*OMEGA*OMEGA * (one - mode_1%kf*mode_1%kf*int_Ws)
                endif 
            enddo 
        enddo 



        ! Parameters a b c: 
        wcomnondim = (SCALE_T * mode_1%wcom)
        write(*,*)"Mode wcom dimensional ", mode_1%wcom

 

        if(mode_1%t .eq.'T')then 
            aparam =  half * (nu/(wcomnondim*wcomnondim) - tau)
        else 
            aparam =  ((OMEGA/wcomnondim)**two)*(one - mode_1%kf*mode_1%kf*int_Ws)/three &
                       +  half * (nu - wcomnondim*wcomnondim*tau)/(wcomnondim*wcomnondim)
        endif 
        bparam =  (int_Ws * OMEGA / wcomnondim)  
        cparam = -three*(nu - wcomnondim*wcomnondim*tau)/(two * wcomnondim*wcomnondim * ksq )

        ! note these are unitless parameters
        write(*,*)'A param: ', aparam * 1000.0
        write(*,*)'B param: ', bparam * 1000.0
        write(*,*)'C param: ', cparam * 1000.0
        write(*,*)"nu  is ", tau
        write(*,*)"tau is ", nu


        ! Not scaled 
        call buffer_int(nstr, mode_1%n)
        call buffer_int(lstr,  mode_1%l)

        write(out_name, '(a)')'./paper_benchmarks/Coriolis/abc/'//trim(nstr)//mode_1%t//trim(lstr)//'_'//trim(nstr)//mode_1%t//trim(lstr)//'.txt'
        open(1,file=trim(out_name), form='formatted')
        write(1,*)aparam 
        write(1,*)bparam 
        write(1,*)cparam 
        close(1)


        write(out_name, '(a)')'./ellipticity/matrices/Wmat_'//trim(nstr)//mode_1%t//trim(lstr)//'_'//trim(nstr)//mode_1%t//trim(lstr)//'.txt'
        call save_W_matrix(mode_1%l, mode_2%l, trim(out_name))

        write(out_name, '(a)')'./ellipticity/matrices/Vell_'//trim(nstr)//mode_1%t//trim(lstr)//'_'//trim(nstr)//mode_1%t//trim(lstr)//'.txt'
        call save_Vell_matrix(mode_1%l, mode_2%l, out_name)

        write(out_name, '(a)')'./ellipticity/matrices/Tell_'//trim(nstr)//mode_1%t//trim(lstr)//'_'//trim(nstr)//mode_1%t//trim(lstr)//'.txt'
        call save_Tell_matrix(mode_1%l, mode_2%l, out_name)

        write(out_name, '(a)')'./ellipticity/matrices/Vcen_'//trim(nstr)//mode_1%t//trim(lstr)//'_'//trim(nstr)//mode_1%t//trim(lstr)//'.txt'
        call save_Vcen_matrix(mode_1%l, mode_2%l, out_name)


        ! Total splitting - we are going to store it in Wmat even though it should be
        ! Hmat
        ! 14.84
        !Wmat = Wmat + (Vell + Vcen - wcomnondim*wcomnondim*Tell)/(two*wcomnondim)
        !write(out_name, '(a)')'./ellipticity/matrices/RotEll_'//trim(nstr)//mode_1%t//trim(lstr)//'_'//trim(nstr)//mode_1%t//trim(lstr)//'.txt'
        !call save_W_matrix(mode_1%l, mode_2%l, trim(out_name))
        

        ! Compute ellipticity contribution cst: 
        ! Nondimensional
        Ellmat =  (Vell  - wcomnondim*wcomnondim*Tell)/(two*wcomnondim)

        ! Dimensionalise the Ellmat into micro Hz
        Ellmat = 1.0e6* Ellmat/(SCALE_T*two*PI)

        ! Compute the cst for this matrix: 
        ! Write as a CST
        call get_Ssum_bounds(mode_1%l, mode_1%l, smin, smax, num_s, ncols)
        allocate(cst(num_s, ncols))
        call Hcomplex_to_cst_8(Ellmat, mode_1%l, mode_1%l, cst, ncols, num_s, t1, t1, 2)
        out_name = 'ellipticity/cst/cst_'//trim(nstr)//t1//trim(lstr)//'_'//trim(nstr)//t1//trim(lstr)
        call write_cst_complex_to_file(out_name, cst, ncols, num_s, smin, 2)



        deallocate(Wmat)
        deallocate(Vcen)
        deallocate(Tell)
        deallocate(Vell)
        deallocate(Ellmat)
        deallocate(cst)



        write(*,*)
    enddo ! i mode



end program semi_analytical_W_matrix