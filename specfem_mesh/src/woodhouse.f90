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
            Tp = -SPLINE_iONE * m_1%v_spl * m_2%w_spl * BNpmlsld(1,-1,m_1%l,s,m_2%l) 
        endif 
    else 
        if(m_2%t.eq.'S')then 
            ! T S

            Tp = SPLINE_iONE * m_1%w_spl * m_2%v_spl * BNpmlsld(1,-1,m_1%l,s,m_2%l) 
        else 
            ! T T
            Tp = m_1%w_spl * m_2%w_spl * BNpmlsld(1,1,m_1%l,s,m_2%l)
        endif 
    endif ! 1 = S

end subroutine WK_Trho



end module woodhouse_kernels