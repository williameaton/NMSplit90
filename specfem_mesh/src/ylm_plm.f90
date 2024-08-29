module ylm_plm

use math, only: cosp, sinp, tanp, sqrtp
implicit none 
include "constants.h"

contains
   
    real(kind=CUSTOM_REAL) function Plm(x,l,mm)
    ! Function Plm Calculates the Associate Legendre Polinomial
    ! computes Plm at point x 
    use params 
    implicit none 

    ! IO variables
    integer, intent(in) :: l,mm
    real(kind=CUSTOM_REAL):: x

    ! Local variables
    real(kind=CUSTOM_REAL):: mf, lf, y, Pmin2, Pmin1, rll
    integer :: q, ll, m
    logical :: neg_m

    ! We only use positive m - if negative then convert at the end
    if (mm .lt. 0)then 
        neg_m = .true.
    else
        neg_m = .false.
    endif 
    m = abs(mm)

    ! Float versions of m and l (positive m)
    mf = abs(real(m, kind=CUSTOM_REAL))
    lf = abs(real(l, kind=CUSTOM_REAL))

    Plm = 1.0
    if (m>l)then
        !write(*,*)'Warning abs(m) > l for Plm: returning 0'
        Plm= 0.0
    elseif (m.eq.0 .and. l.eq.0)then 
        return ! Plm = 1
    end if


    ! (1 - x^2)^0.5
    y = sqrtp(1-x**2)
    
    ! Increase degree and order from P_0^0 to P_q^q until q = desired m 
    do q = 0, m-1
        Plm = (TWO*real(q, kind=CUSTOM_REAL)+ONE) * y * Plm
    enddo 

    if (l.eq.m)then 
        call convert_neg(neg_m, Plm, mf, lf)
        return
    endif


    ! Use P_m^m to compute P_{l=m+1}^{m}
    ! This is stored in Pmin 1
    ! Pmin_2 is holding P_m^m 
    Pmin2 = Plm
    Pmin1 = (2*mf + 1)*x*Pmin2

    ! If what we actually want is P_{m, l=m+1} then we can return Pmin1
    if (l == m+1)then
        Plm = Pmin1
        call convert_neg(neg_m, Plm, mf, lf)
        return
    endif

    do ll = m+1, l-1
        rll = real(ll,kind=CUSTOM_REAL)
        Plm = ((TWO*rll + ONE) * x * Pmin1 - (rll + mf) * Pmin2) / (rll - mf + ONE)
        Pmin2 = Pmin1
        Pmin1 = Plm
    enddo 
    call convert_neg(neg_m, Plm, mf, lf)

    return 
    end function Plm
    
    ! ------------------------------------------------------------------


    subroutine convert_neg(neg_m, Plm, mf, lf)
        implicit none 
        !include "constants.h"


        logical ::  neg_m
        real(kind=CUSTOM_REAL) :: Plm, mf, lf

        !write(*,*)lf,mf,Plm

        if (neg_m) then 

            !write(*,*)"(-ONE)**mf", (-ONE)**mf
            !write(*,*)"Gamma(lf-mf+ONE)", Gamma(lf-mf+ONE)
            !write(*,*)"Gamma(lf+mf+ONE)", Gamma(lf+mf+ONE)
            !write(*,*)"Plm", Plm


            Plm = (-ONE)**mf * (Gamma(lf-mf+ONE)/Gamma(lf+mf+ONE)) * Plm
        else
            Plm = Plm
        endif 

        return 
    end subroutine convert_neg

    ! ------------------------------------------------------------------
    real(kind=CUSTOM_REAL) function xlm(l,m,theta)
        implicit none 
        include "constants.h"
    
        ! IO variables: 
        integer :: l, m 
        real(kind=CUSTOM_REAL) :: theta
        ! Local: 
        real(kind=CUSTOM_REAL) :: lf, mf, blm, abs_m, bl0, pith

        ! Convert l and m to floats for computation: 
        lf = real(l, kind=CUSTOM_REAL)
        mf = real(m, kind=CUSTOM_REAL)

        ! Gamma function isnt defined in these cases
        ! But the plm is 0 so should xlm be 
        ! Technically plm isnt defined for |m| > l
        if (abs(m).gt.l)then
            xlm = zero
            return
        endif        

    
        ! Xlm has limiting behaviour at the poles: 
        ! It is 0 for any m != 0 terms 
        ! X_{l0}(theta = 0)  = (2l+1 / 4pi)**0.5     
        ! X_{l0}(theta = pi) = (-1)**l * (2l+1 / 4pi)**0.5     
        abs_m = abs(mf)
        bl0   = ((two*lf + one)/(four*PI))**half 
        blm   = (((-1)**mf)/((TWO**abs_m)*gamma(abs_m+1))) * & 
                bl0 * & 
                (gamma(lf + abs_m + one)/gamma(lf - abs_m + one))**half


        if(abs(theta).le.pole_tolerance) then 
            ! theta close to 0 -- North Pole 
            if (m.ge.-l .and. m.lt.0) then
                xlm = ((-1)**mf) * blm * theta**abs_m
            elseif (m.eq.0) then
                xlm = bl0 * (one - lf*(lf+one)*theta*theta/four)
            elseif (m.lt.l .and. m.gt.0) then
                xlm = blm * (theta**mf)
            else 
                xlm = zero
            endif
        elseif (abs(PI - theta) .le. pole_tolerance)then 
            ! theta close to pi -- South Pole 
            pith = PI - theta
            if (m.ge.-l .and. m.lt.0) then
                xlm = ((-1)**lf) * blm * pith**abs_m
            elseif (m.eq.0) then
                xlm = ((-1)**lf) * bl0 * (one - lf*(lf+one)*pith*pith/four)
            elseif (m.lt.l .and. m.gt.0) then
                xlm = blm * (pith**mf) * (-1)**(lf+mf) 
            else 
                xlm = zero
            endif
        elseif (theta.gt.pole_tolerance .and. theta.lt.(PI-pole_tolerance)) then 
            xlm = ((-ONE)**mf) *                         & 
                ( ((TWO*lf)+ONE)*gamma(lf-mf+ONE)/         & 
                    (FOUR*PI*gamma(lf+mf+ONE))  )**HALF  & 
                *Plm(cosp(theta),l,m)   
        else
            write(*,*)'Error in Xlm: theta value is', theta
        endif
            


    end function xlm

    ! ------------------------------------------------------------------

    real(kind=CUSTOM_REAL) function Fact(n)
    ! Computes a factorial for n
    ! this code was taken from: 
    ! https://gitlab.surrey.ac.uk/phs3ps/Sky3D/-/blob/v1.2/Code/ylm.f90
        implicit none 
        integer:: n,i
        Fact=1.0d0
        do i =1,n
            Fact = Fact*i
        end do
    end function Fact

    ! ------------------------------------------------------------------

    real(kind=CUSTOM_REAL) function Fact2(n)
    ! Calculates the double Factorial
    ! this code was taken from: 
    ! https://gitlab.surrey.ac.uk/phs3ps/Sky3D/-/blob/v1.2/Code/ylm.f90
        implicit none
        integer :: n
        Fact2=1.0d0
        do while (n>0)
            Fact2 = Fact2*n
        n=n-2
        end do
    end function Fact2

    ! ------------------------------------------------------------------

    real(kind=CUSTOM_REAL) function ylm_real(l, m, theta, phi)
    ! Calculates the real value of Ylm(theta, phi) following definition in 
    ! Dahlen and Tromp 1998 Eqn. B.58: 
    ! Theta should be between 0 and pi
    ! Phi should be between 0 and 2 pi
    ! Y_lm(theta, phi) = X_{lm}(theta) e^{i m phi}
        implicit none
        !include "constants.h"

        ! IO variables: 
        integer :: l, m 
        real(kind=CUSTOM_REAL) :: theta, phi 
        ! Local variables
        real(kind=CUSTOM_REAL) :: lf, mf

        ! Sanity checks on inputs: 
        if (theta.lt.0.0d0 .or. theta.gt.PI+PI_TOL)then 
            write(*,*)'ERROR: Theta should be between 0 and pi but value is ', theta
            stop
        endif 
        if (phi.lt.0.0d0 .or. phi.gt.TWO_PI+PI_TOL)then 
            write(*,*)'ERROR: Phi should be between 0 and 2pi but value is ', phi
            stop
        endif 
        if (l.lt.0)then 
            write(*,*)'ERROR: l must be >= 0 but has value ', l
            stop
        endif 
        if (abs(m).gt.l)then 
            write(*,*)'ERROR: m must be <= |l| but has value ', m
            stop
        endif 

        ! Convert l and m to floats for computation: 
        lf = real(l, kind=CUSTOM_REAL)
        mf = real(m, kind=CUSTOM_REAL)

        ! Real component is the cos phi part
        ylm_real = cosp(mf*phi)*xlm(l,m,theta)

    end function ylm_real

    ! ------------------------------------------------------------------

    complex(kind=CUSTOM_REAL) function ylm_complex(l, m, theta, phi)
    ! Calculates the complex Ylm(theta, phi) following definition in 
    ! Dahlen and Tromp 1998 Eqn. B.58: 
    ! Theta should be between 0 and pi
    ! Phi should be between 0 and 2 pi
    ! Y_lm(theta, phi) = X_{lm}(theta) e^{i m phi}
        implicit none

        ! IO variables: 
        integer :: l, m 
        real(kind=CUSTOM_REAL) :: theta, phi 
        ! Local variables
        real(kind=CUSTOM_REAL) :: lf, mf
        ! Sanity checks on inputs: 
        if (theta.lt.0.0d0 .or. theta.gt.PI+PI_TOL)then 
            write(*,*)'ERROR: Theta should be between 0 and pi but value is ', theta
            stop
        endif 
        if (phi.lt.0.0d0 .or. phi.gt.TWO_PI+PI_TOL)then 
            write(*,*)'ERROR: Phi should be between 0 and 2pi but value is ', phi
            stop
        endif 
        if (l.lt.0)then 
            write(*,*)'ERROR: l must be >= 0 but has value ', l
            stop
        endif 
        if (abs(m).gt.l)then 
            write(*,*)'ERROR: m must be <= |l| but has value ', m
            stop
        endif 

        ! Convert l and m to floats for computation: 
        lf = real(l, kind=CUSTOM_REAL)
        mf = real(m, kind=CUSTOM_REAL)

        ! Assign complex Ylm
        ylm_complex =  cmplx(cosp(mf*phi), sinp(mf*phi)) * xlm(l,m,theta)

        

    end function ylm_complex

    ! ------------------------------------------------------------------

    subroutine ylm_deriv(l, m, theta, phi, dylm_dth, dylm_dphi)
        ! Computes partial derivatives of Ylm at theta, phi with
        ! respect to phi (dylm_dphi) and theta (dylm_dth)
        implicit none

        ! IO variables: 
        integer :: l, m 
        real(kind=CUSTOM_REAL) :: theta, phi
        complex(kind=CUSTOM_REAL) :: dylm_dphi, dylm_dth
        ! Local variables
        real(kind=CUSTOM_REAL) :: lf, mf, dxlm_dt


        ! Sanity checks on inputs: 
        if (theta.lt.0.0d0 .or. theta.gt.PI+PI_TOL)then 
            write(*,*)'ERROR: Theta should be between 0 and pi but value is ', theta
            stop
        endif 
        if (phi.lt.0.0d0 .or. phi.gt.TWO_PI+PI_TOL)then 
            write(*,*)'ERROR: Phi should be between 0 and 2pi but value is ', phi
            stop
        endif 
        if (l.lt.0)then 
            write(*,*)'ERROR: l must be >= 0 but has value ', l
            stop
        endif 
        if (abs(m).gt.l)then 
            write(*,*)'ERROR: m must be <= |l| but has value ', m
            stop
        endif 

        ! Convert l and m to floats for computation: 
        lf = real(l, kind=CUSTOM_REAL)
        mf = real(m, kind=CUSTOM_REAL)

        ! d Ylm / d phi:
        ! ylm_complex   = cmplx(cos(mf*phi)*xlm, sin(mf*phi)*xlm)
        ! d Ylm / d phi = Xlm(theta) * (-msin(m phi) + im cos(m phi) )
        dylm_dphi = cmplx(-mf*sinp(mf*phi), mf*cosp(mf*phi)) * xlm(l,m,theta)

        ! d Ylm / d theta:
        ! DT98 B.120
        !dxlm_dt = HALF*(  xlm(l,m+1,theta)*((lf-mf)*(lf+mf+1))**HALF   & 
        !                - xlm(l,m-1,theta)* ((lf+mf)*(lf-mf+1))**HALF  )
        !dylm_dth = cmplx(cos(mf*phi), sin(mf*phi)) * dxlm_dt
        ! Using the above option produces the wrong sign for dylm_dth when testing for Ylm with l = m
        ! But they define their ylm with the -1
        ! https://math.stackexchange.com/questions/3256898/partial-derivatives-of-m-l-spherical-harmonics

        dxlm_dt = half * ( xlm(l, m+1, theta) * ((lf-mf)*(lf+mf+one))**half  - & 
                           xlm(l, m-1, theta) * ((lf+mf)*(lf-mf+one))**half    )
  

        dylm_dth = cmplx(cosp(mf*phi), sinp(mf*phi)) * dxlm_dt

        !write(*,*)'dylm_dth :', dylm_dth
        ! Wolfram alpha solution but note their Ylm doesnt contain the (-1)**m 
        !write(*,*)'dylm_dth :', m*ylm_complex(l,m,theta,phi)/(tan(theta)*(-one)**mf) + cmplx(cosp(-phi), sinp(-phi))* (ylm_complex(l,m+1,theta, phi)/((-one)**(mf+1))) *((lf-mf)*(lf+mf+one))**half



    end subroutine ylm_deriv

end module ylm_plm