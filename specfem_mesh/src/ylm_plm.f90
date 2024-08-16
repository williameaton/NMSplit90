module ylm_plm

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
    y = sqrt(1-x**2)
    
    do q = 0, m-1
        Plm = -(TWO*real(q, kind=CUSTOM_REAL)+ONE) * y * Plm
    enddo 

    if (l.eq.m)then 
        call convert_neg(neg_m, Plm, mf, lf)
        return
    endif

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
        real(kind=CUSTOM_REAL) :: lf, mf, tmp, pole_tolerance

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

        pole_tolerance = 1e-4 ! ~ 600 metres of pole

        
        ! Xlm has limiting behaviour at the poles: 
        ! It is 0 for any m != 0 terms 
        ! X_{l0}(theta = 0)  = (2l+1 / 4pi)**0.5     
        ! X_{l0}(theta = pi) = (-1)**l * (2l+1 / 4pi)**0.5     

        if (theta.gt.zero .and. theta.le.pole_tolerance) then 
            ! theta approx 0 -- North Pole 
            if (m.eq.0)then 
                xlm = sqrt((TWO*lf + ONE)/(4*PI))
            else
                xlm = zero
            endif 

        elseif (abs(PI - theta) .le. pole_tolerance)then 
            if (m.eq.0)then 
                xlm = (lf**(-ONE))*sqrt((TWO*lf + ONE)/(4*PI))
            else
                xlm = zero
            endif 
        else 
            xlm = ((-ONE)**mf) *                         & 
                ( (TWO*lf+ONE)*gamma(lf-mf+ONE)/        & 
                    (FOUR*PI*gamma(lf+mf+ONE))  )**HALF       & 
                *Plm(cos(theta),l,m)      
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
        ylm_real = cos(mf*phi)*xlm(l,m,theta)

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
        ylm_complex = cmplx(cos(mf*phi), sin(mf*phi))*xlm(l,m,theta)

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
        real(kind=CUSTOM_REAL) :: lf, mf, dplm_dt, dxlm_dt

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
        dylm_dphi = cmplx(-mf*sin(mf*phi), mf*cos(mf*phi)) * xlm(l,m,theta)

        ! d Ylm / d theta:
        ! DT98 B.120
        !dxlm_dt = HALF*(  xlm(l,m+1,theta)*((lf-mf)*(lf+mf+1))**HALF   & 
        !                - xlm(l,m-1,theta)* ((lf+mf)*(lf-mf+1))**HALF  )
        !dylm_dth = cmplx(cos(mf*phi), sin(mf*phi)) * dxlm_dt
        ! DT98 B.116
        ! Using the above option produces the wrong sign for dylm_dth when testing for Ylm with l = m
        ! https://math.stackexchange.com/questions/3256898/partial-derivatives-of-m-l-spherical-harmonics
        dplm_dt = lf*Plm(cos(theta),l,m)/tan(theta) - (lf+mf)*Plm(cos(theta),l-1,m)/sin(theta)     
        dylm_dth = cmplx(cos(mf*phi), sin(mf*phi)) * dplm_dt * ((-ONE)**mf)*((TWO*lf+ONE)*gamma(lf-mf+ONE)/(FOUR*PI*gamma(lf+mf+ONE)))**HALF

    end subroutine ylm_deriv

end module ylm_plm