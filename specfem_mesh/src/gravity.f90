module gravitation 

    implicit none 
    include "constants.h"

    contains


    subroutine compute_background_g(r, p, N, g)
        ! Computes background MAGNITUDE of gravitational acceleration
        ! r is radius
        ! p is density
        ! N is length of array
        ! g is output
        use integrate, only: integrate_r_traps
        implicit none 
        integer                :: N, i
        real(kind=CUSTOM_REAL) :: r(N), p(N), g(N), integrand(N), integral

        g = zero 
        do i = 2, N 
            !integrand = zero 
            !integrand(1:N) = p(1:N) * r(1:N) * r(1:N)
            integral = integrate_r_traps(r(1:N), p(1:N)*r(1:N)*r(1:N), i)
            g(i) = FOUR * integral /(r(i)*r(i))

            if(r(i).eq.zero)g(i)=zero

        enddo 



    end subroutine

end module gravitation