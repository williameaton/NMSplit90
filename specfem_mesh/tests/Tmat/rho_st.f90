! For this benchmark test we want to prescribe coefficients of delta rho st 
! that vary radially. In the test we are using smax = 3 so need radial 
! profiles for each rho_st of s = 1, 2, 3
! For simplicity lets assume that the t values of each s have the same radial 
! profile, which is multiplied by some constant
module rho_st_profiles 
    use math, only: sinp
    implicit none 
    include "constants.h"


    contains 

    subroutine radial_rho_st(s, npts, coeff, arr, radius)
        implicit none 
        integer :: s, npts, i
        complex(kind=SPLINE_REAL) :: coeff
        complex(kind=SPLINE_REAL) :: arr(npts)
        real(kind=CUSTOM_REAL)    :: radius(npts)


        if(s.eq.1)then 
            ! When s = 1 
            ! It will be quadratic, going through a value of +0.6 at each end and -0.4 at the centre 
            arr  =  coeff * ((six/ten) + (radius - half)**two )
        elseif(s.eq.2)then 
            do i = 1, npts
                arr(i) =  coeff * (sinp(TWO_PI*radius(i)))
            enddo 
        elseif(s.eq.3)then
            arr  =  coeff * ((seven/ten) - TWO*radius  )
        else
            write(*,*)'Only supporting s=1,3 but s = ', s
            stop 
        endif

    end subroutine



    ! Same as above but for a single r value: 
    complex(kind=SPLINE_REAL) function single_r_rho_st(s, coeff, radius)
        implicit none 
        integer :: s
        complex(kind=SPLINE_REAL) :: coeff
        real(kind=CUSTOM_REAL)    :: radius


        if(s.eq.1)then 
            ! When s = 1 
            ! It will be quadratic, going through a value of +0.6 at each end and -0.4 at the centre 
            single_r_rho_st =  coeff * ((six/ten) + (radius - half)**two )
        elseif(s.eq.2)then 
            single_r_rho_st =  coeff * (sinp(TWO_PI*radius)) 
        elseif(s.eq.3)then
            single_r_rho_st =  coeff * ((seven/ten) -TWO*radius  )
        else
            write(*,*)'Only supporting s=1,3 but s = ', s
            stop 
        endif

    end function single_r_rho_st


end module 