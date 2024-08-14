module ylm_plm

include "precision.h"
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
        write(*,*)'Warning abs(m) > l for Plm: returning 0'
        Plm= 0
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
        use params, only: ONE

        implicit none 
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
        use params, only: ONE
        implicit none

        ! IO variables: 
        integer :: l, m 
        real(kind=CUSTOM_REAL) :: theta, phi, pp 
        ! Local variables
        real(kind=CUSTOM_REAL) :: lf, mf, xlm

        ! Sanity checks on inputs: 
        if (theta.lt.0.0d0 .or. theta.gt.PI)then 
            write(*,*)'ERROR: Theta should be between 0 and pi but value is ', theta
            stop
        endif 
        if (phi.lt.0.0d0 .or. phi.gt.TWO_PI)then 
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

        ! Note gamma(n+1) = fact(n)
        xlm = ((-ONE)**mf) *                         & 
               ( (2.0d0*lf+ONE)*gamma(lf-mf+ONE)/        & 
                (4.0d0*PI*gamma(lf+mf+ONE))  )**0.5       & 
              *Plm(cos(theta),l,m)

        ! Real component is the cos phi part
        ylm_real = cos(mf*phi)*xlm

    end function ylm_real

end module ylm_plm