module woodhouse_kernels
    use params, only: safety_checks
    use w3j,    only: thrj
    use modes,  only: Mode
    implicit none 
    include "constants.h"

    contains

    real(kind=CUSTOM_REAL) function BNpmlsld(N, pm, l, s, ld)
    ! Computes D.42 of DT98 B^{(N)\pm}_{l s l'}
    ! To indicate plus or minus should be 1 or -1

    implicit none 
    integer                ::  N, pm, l, s, ld
    real(kind=CUSTOM_REAL) :: fN, fl, fld, fpm, fs


    if(safety_checks)then 
        if(abs(pm).ne.1)then 
            write(*,*)"Error in BNpmlsld. pm is not 1 or -1: ", pm
            stop 
        endif
    endif 

    fN   = real(N,  kind=CUSTOM_REAL)
    fl   = real(l,  kind=CUSTOM_REAL)
    fs   = real(s,  kind=CUSTOM_REAL)
    fld  = real(ld, kind=CUSTOM_REAL)
    fpm  = real(pm, kind=CUSTOM_REAL)

    BNpmlsld = half * (-one)**fN * (one + (fpm)*(-one)**(fl+fs+fld)) & 
                    * thrj(l, s, ld, -N, 0, N)                     * &
                      ((gamma(fl+fN+one)*gamma(fld+fN+one))/         & 
                       (gamma(fl-fN+one)*gamma(fld-fN+one)))**half  

    end function BNpmlsld



subroutine WK_Trho(m_1, m_2, s, Tp)
    ! D.46 Woodhouse kernel T_\rho
    implicit none 

    type(Mode) :: m_1
    type(Mode) :: m_2
    integer    :: s

    complex(kind=SPLINE_REAL) :: Tp(m_1%spl_len)

    if(safety_checks)then 
        if(m_1%spl_len.ne.m_2%spl_len)then 
            write(*,*)"Error in WK_Trho. Mode splines arent same length "
            write(*,*)"Mode 1: ", m_1%spl_len
            write(*,*)"Mode 2: ", m_2%spl_len
            stop 
        endif
    endif 


    if(m_1%t.eq.'S')then 
        if(m_2%t.eq.'S')then 
            ! S S

            Tp = m_1%u_spl * m_2%u_spl * BNpmlsld(0,1,m_1%l,s,m_2%l) + &
                 m_1%v_spl/m_1%kf * m_2%v_spl/m_2%kf * BNpmlsld(1,1,m_1%l,s,m_2%l) 
        else 
            ! S T
            Tp = -SPLINE_iONE * m_1%v_spl/m_1%kf * m_2%w_spl/m_2%kf * BNpmlsld(1,-1,m_1%l,s,m_2%l) 
        endif 
    else 
        if(m_2%t.eq.'S')then 
            ! T S
            Tp = SPLINE_iONE * m_1%w_spl/m_1%kf * m_2%v_spl/m_2%kf * BNpmlsld(1,-1,m_1%l,s,m_2%l) 
        else 
            ! T T
            Tp = m_1%w_spl/m_1%kf * m_2%w_spl/m_2%kf * BNpmlsld(1,1,m_1%l,s,m_2%l)
        endif 
    endif ! 1 = S

end subroutine WK_Trho


subroutine WK_Vrho(m_1, m_2, s, spline_rad, rho_spl, g_spl, Vp, addgrav)
    ! D.50 Woodhouse kernel V_\rho
    ! if addgrav is true then it modifies the kernel as D.55
    use integrate, only: integrate_r_traps
    implicit none 

    type(Mode) :: m_1
    type(Mode) :: m_2
    integer    :: s, N, ir 
    real(kind=CUSTOM_REAL) :: fs, rr
    logical    :: addgrav
    complex(kind=SPLINE_REAL) :: integral_below, integral_above

    complex(kind=SPLINE_REAL),dimension(m_1%spl_len) :: & 
     Vp, grav_eff, Vphi, Vdotphi, b_integrand, a_integrand
    real(kind=CUSTOM_REAL), dimension(m_1%spl_len) :: spline_rad, g_spl, rho_spl

    if(safety_checks)then 
        if(m_1%spl_len.ne.m_2%spl_len)then 
            write(*,*)"Error in WK_Trho. Mode splines arent same length "
            write(*,*)"Mode 1: ", m_1%spl_len
            write(*,*)"Mode 2: ", m_2%spl_len
            stop 
        endif
    endif 

    if(m_1%t.eq.'S')then 
        if(m_2%t.eq.'S')then 
            ! S S NOTE HERE I AM ASSUMING THE 8piG non dimensionalised
            ! is just 8 
            Vp = ((m_1%u_spl*m_2%dp_spl + m_1%dp_spl*m_2%u_spl            -  & 
                  half*g_spl*( m_1%aux_f*m_2%u_spl + m_1%u_spl*m_2%aux_f     &
                              + four*m_1%u_spl*m_2%u_spl/spline_rad)      +  & 
                  eight*rho_spl*m_1%u_spl*m_2%u_spl) *                       &
                  BNpmlsld(0,1,m_1%l,s,m_2%l))                               &      
                + ((m_1%p_spl*m_2%v_spl/m_2%kf + m_1%v_spl*m_2%p_spl/m_1%kf + & 
                  half*g_spl*(m_1%u_spl*m_2%v_spl/m_2%kf +                   & 
                              m_1%v_spl*m_2%u_spl/m_1%kf)) *                 &
                  BNpmlsld(1,1,m_1%l,s,m_2%l)/spline_rad)            
        else 
            ! S T
            ! Note i have factored the division by k for the w splines
            Vp = - BNpmlsld(1,-1,m_1%l,s,m_2%l)                          * &
                    SPLINE_iONE*( m_1%p_spl * m_2%w_spl                  + & 
                                half*g_spl * m_1%u_spl * m_2%w_spl         &
                               )/(spline_rad*m_2%kf) 
        endif 
    else 
        if(m_2%t.eq.'S')then 
            ! T S
            ! Note i have factored the division by k for the w splines
            Vp =  BNpmlsld(1,-1,m_1%l,s,m_2%l)                           * &
                   SPLINE_iONE*( m_2%p_spl * m_1%w_spl                   + & 
                                 half*g_spl * m_2%u_spl * m_1%w_spl        &
                               )/(spline_rad*m_1%kf) 
        else 
            ! T T
            Vp = SPLINE_iZERO
        endif 
    endif ! 1 = S

    if(spline_rad(1).eq.zero) Vp(1) = SPLINE_iZERO


    ! Modify the density kernel to incorporate background gravity effects
    if(addgrav)then

        N =  m_1%spl_len
        fs = real(s, kind=CUSTOM_REAL)
        grav_eff = SPLINE_iZERO
        
        ! Get the vphi and vphidot
        call WK_Vphi(m_1, m_2, s, spline_rad, rho_spl, Vphi)
        call WK_Vphi_dot(m_1, m_2, s, spline_rad, rho_spl, Vdotphi)

        b_integrand = spline_rad**(fs+one)*(fs*Vdotphi + spline_rad*Vphi)
        a_integrand = spline_rad**(-fs)*( (fs+one)*Vdotphi - spline_rad*Vphi )

        ! Avoid nan
        if(spline_rad(1).eq.zero)then 
            a_integrand(1) = zero
        endif 

        ! Loop for each radius and compute the integrals: 
        do ir = 1, N
            rr = spline_rad(ir)

            if(ir.eq.1)then
                integral_below = SPLINE_iZERO
            else  
                integral_below = rr**(-fs-one) * integrate_r_traps(spline_rad(1:ir) , & 
                                                                   b_integrand(1:ir), &
                                                                   ir)
            endif 

            if(ir.eq.N)then 
                integral_above = SPLINE_iZERO
            else 
                integral_above = (rr**fs) * integrate_r_traps(spline_rad(ir:N),  & 
                                                            a_integrand(ir:N), &
                                                            N-ir+1)
            endif 
            grav_eff(ir) = four*(integral_above - integral_below)/(2*fs + one)

            if(grav_eff(ir).ne.grav_eff(ir))then 
                write(*,*)'ir = ', ir
                write(*,*)'Nan found'
                stop 
            endif
        enddo 

        ! Combine 
        Vp = Vp + grav_eff
    endif

end subroutine WK_Vrho






subroutine WK_Vkappa(m_1, m_2, s, spline_rad, Vk)
    ! D.48 Woodhouse kernel V_\kappa
    implicit none 

    type(Mode) :: m_1
    type(Mode) :: m_2
    integer    :: s

    complex(kind=SPLINE_REAL) :: Vk(m_1%spl_len)
    real(kind=CUSTOM_REAL)    :: spline_rad(m_1%spl_len)

    if(safety_checks)then 
        if(m_1%spl_len.ne.m_2%spl_len)then 
            write(*,*)"Error in WK_Vkappa. Mode splines arent same length "
            write(*,*)"Mode 1: ", m_1%spl_len
            write(*,*)"Mode 2: ", m_2%spl_len
            stop 
        endif
    endif 

    if(m_1%t.eq.'S' .and. m_2%t.eq.'S')then 
        Vk = (m_1%du_spl + m_1%aux_f) * &
             (m_2%du_spl + m_2%aux_f) * &
             BNpmlsld(0,1,m_1%l,s,m_2%l)
        
        if(spline_rad(1).eq.zero) Vk(1) = SPLINE_iZERO
    else 
        Vk = SPLINE_iZERO
    endif 


end subroutine WK_Vkappa




subroutine WK_Vmu(m_1, m_2, s, spline_rad, Vmu)
    ! D.49 Woodhouse kernel T_\mu
    implicit none 

    type(Mode) :: m_1
    type(Mode) :: m_2
    integer    :: s

    complex(kind=SPLINE_REAL) :: Vmu(m_1%spl_len)
    real(kind=CUSTOM_REAL)    :: spline_rad(m_1%spl_len)

    if(safety_checks)then 
        if(m_1%spl_len.ne.m_2%spl_len)then 
            write(*,*)"Error in WK_Trho. Mode splines arent same length "
            write(*,*)"Mode 1: ", m_1%spl_len
            write(*,*)"Mode 2: ", m_2%spl_len
            stop 
        endif
    endif 


    if(m_1%t.eq.'S')then 
        if(m_2%t.eq.'S')then 
            ! S S
            Vmu = ((SPLINE_TWO*m_1%du_spl - m_1%aux_f)              * & 
                   (SPLINE_TWO*m_2%du_spl - m_2%aux_f)              * & 
                   BNpmlsld(0,1,m_1%l,s,m_2%l))/SPLINE_THREE        + & 
                  (m_1%aux_x*m_2%aux_x*BNpmlsld(1,1,m_1%l,s,m_2%l)) + & 
                  ((m_1%v_spl/m_1%kf) * (m_2%v_spl/m_2%kf )         * & 
                   BNpmlsld(2,1,m_1%l,s,m_2%l)/(spline_rad*spline_rad) )     
        else 
            ! S T
            Vmu = -SPLINE_iONE * (m_1%aux_x * m_2%aux_z * BNpmlsld(1,-1,m_1%l,s,m_2%l) + &
                                  (m_1%v_spl*m_2%w_spl*BNpmlsld(2,-1,m_1%l,s,m_2%l) & 
                                    /(m_1%kf*m_2%kf*spline_rad*spline_rad))  )
        endif 
    else 
        if(m_2%t.eq.'S')then 
            ! T S
            Vmu = SPLINE_iONE*(m_1%aux_z * m_2%aux_x * BNpmlsld(1,-1,m_1%l,s,m_2%l) + &   
                                (m_1%w_spl*m_2%v_spl*BNpmlsld(2,-1,m_1%l,s,m_2%l) & 
                                   /(m_1%kf*m_2%kf*spline_rad*spline_rad))    )
        else 
            ! T T
            Vmu = m_1%aux_z*m_2%aux_z +  (m_1%w_spl*m_2%w_spl*BNpmlsld(2,1,m_1%l,s,m_2%l) & 
                                            /(m_1%kf*m_2%kf*spline_rad*spline_rad))
        endif 
    endif ! 1 = S


    if(spline_rad(1).eq.zero) Vmu(1) = SPLINE_iZERO

end subroutine WK_Vmu






subroutine WK_Vphi(m_1, m_2, s, spline_rad, rho_spl, Vphi)
    ! D.51 Woodhouse kernel V_\Phi
    implicit none 

    type(Mode) :: m_1
    type(Mode) :: m_2
    integer    :: s

    real(kind=CUSTOM_REAL), dimension(m_1%spl_len) :: rho_spl
    complex(kind=SPLINE_REAL) :: Vphi(m_1%spl_len)
    real(kind=CUSTOM_REAL)    :: spline_rad(m_1%spl_len), fs
    integer :: l1, l2

    l1 = m_1%l
    l2 = m_2%l


    if(safety_checks)then 
        if(m_1%spl_len.ne.m_2%spl_len)then 
            write(*,*)"Error in WK_Vkappa. Mode splines arent same length "
            write(*,*)"Mode 1: ", m_1%spl_len
            write(*,*)"Mode 2: ", m_2%spl_len
            stop 
        endif
    endif 

    fs = real(s, kind=CUSTOM_REAL)


    if(m_1%t.eq.'S') then
        if(m_2%t.eq.'S') then
            ! S S: Use the first, second, and third lines of the equation
            Vphi = (fs * (fs + one) * rho_spl * m_1%u_spl * m_2%u_spl * BNpmlsld(0, 1, l1, s, l2)/(spline_rad*spline_rad)) + & 
                   (half * rho_spl * (m_1%u_spl * m_2%dv_spl/m_2%kf - m_1%du_spl * m_2%v_spl/m_2%kf + &
                     (m_1%u_spl * m_2%v_spl/m_2%kf)/spline_rad - two*m_1%aux_f*m_2%v_spl/m_2%kf ) * BNpmlsld(1, 1, l2, l1, s) /(spline_rad)) + & 
                   (half * rho_spl * (m_2%u_spl * m_1%dv_spl/m_2%kf - m_2%du_spl * m_1%v_spl/m_1%kf + &
                     (m_2%u_spl * m_1%v_spl/m_1%kf)/spline_rad - two*m_2%aux_f*m_1%v_spl/m_1%kf ) * BNpmlsld(1, 1, l1, l2, s) /(spline_rad))    
        else
            ! S T: Use only the third line of the equation
            Vphi = half * SPLINE_iONE * rho_spl * (m_1%u_spl * m_2%dw_spl/m_2%kf - m_1%du_spl * m_2%w_spl/m_2%kf + &
                     (m_1%u_spl * m_2%w_spl/m_2%kf)/spline_rad - two * m_1%aux_f * m_2%w_spl/m_2%kf) * BNpmlsld(1, -1, l2, l1, s) / spline_rad
        endif
  
    else
        if(m_2%t.eq.'S') then
            ! T S: Use only the last line of the equation
            Vphi = - half * SPLINE_iONE * rho_spl * (m_2%u_spl * m_1%dw_spl/m_1%kf - m_2%du_spl * m_1%w_spl/m_1%kf + &
                     (m_2%u_spl * m_1%w_spl/m_1%kf)/spline_rad - two * m_2%aux_f * m_1%w_spl/m_1%kf) * BNpmlsld(1, -1, l2, l1, s) / spline_rad
        else
            ! T T: Vphi is zero
            Vphi = SPLINE_iZERO
        endif
    endif

    if(spline_rad(1).eq.zero) Vphi(1) = SPLINE_iZERO


end subroutine WK_Vphi



subroutine WK_Vphi_dot(m_1, m_2, s, spline_rad, rho_spl, Vdotphi)
    ! D.52 Woodhouse kernel V_\dot{\Phi}
    implicit none 

    type(Mode) :: m_1
    type(Mode) :: m_2
    integer    :: s

    real(kind=CUSTOM_REAL), dimension(m_1%spl_len) :: rho_spl
    complex(kind=SPLINE_REAL) :: Vdotphi(m_1%spl_len)
    real(kind=CUSTOM_REAL)    :: spline_rad(m_1%spl_len), fs
    integer :: l1, l2

    l1 = m_1%l
    l2 = m_2%l

    if(safety_checks)then 
        if(m_1%spl_len.ne.m_2%spl_len)then 
            write(*,*)"Error in WK_Vkappa. Mode splines arent same length "
            write(*,*)"Mode 1: ", m_1%spl_len
            write(*,*)"Mode 2: ", m_2%spl_len
            stop 
        endif
    endif 

    if(m_1%t.eq.'S') then
        if(m_2%t.eq.'S') then
            ! S S:
            Vdotphi = -rho_spl*(m_1%aux_f*m_2%u_spl + m_1%u_spl*m_2%aux_f)*BNpmlsld(0,1,l1,s,l2) + & 
                      (half * rho_spl * ( m_1%u_spl*(m_2%v_spl/m_2%kf)*BNpmlsld(1,1,l2,l1,s) &
                                        + m_2%u_spl*(m_1%v_spl/m_1%kf)*BNpmlsld(1,1,l1,l2,s)) / spline_rad)
        else
            ! S T: 
            Vdotphi = half*SPLINE_iONE*rho_spl*m_1%u_spl*m_2%w_spl*BNpmlsld(1,-1,l2,l1,s) /(m_2%kf*spline_rad)
        endif
  
    else
        if(m_2%t.eq.'S') then
            ! T S: 
            Vdotphi = -half*SPLINE_iONE*rho_spl*m_2%u_spl*m_1%w_spl*BNpmlsld(1,-1,l1,l2,s) /(m_1%kf*spline_rad)
        else
            ! T T: 
            Vdotphi = SPLINE_iONE
        endif
    endif

    if(spline_rad(1).eq.zero) Vdotphi(1) = SPLINE_iZERO


end subroutine WK_Vphi_dot





end module woodhouse_kernels