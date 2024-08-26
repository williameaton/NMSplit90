module integrate 
    implicit none 
    include "constants.h"


    interface integrate_over_mesh
        module procedure integrate_real_mesh_scalar     
        module procedure integrate_complex_mesh_scalar
    end interface

    contains

     real(kind=CUSTOM_REAL) function integrate_real_mesh_scalar(scalar)
        ! Integrates a real scalar, defined at each GLL point, over the mesh
        use params, only: ngllx, nglly, ngllz, nspec, detjac, wgll
        implicit none 
        ! IO variables: 
        real(kind=CUSTOM_REAL) :: scalar(ngllx, nglly, ngllz, nspec)
        ! Local 
        real(kind=CUSTOM_REAL) :: sum
        integer :: i, j, k, ispec

        sum = zero 
        do ispec = 1, nspec 
            do i = 1, ngllx 
                do j = 1, nglly
                    do k = 1, ngllz 
                        sum = sum + scalar(i, j, k, ispec) * detjac(i, j, k, ispec) * wgll(i) * wgll(j) * wgll(k)
                    enddo 
                enddo 
            enddo 
        enddo 
        integrate_real_mesh_scalar = sum 
        return
    end function  integrate_real_mesh_scalar

    
    complex(kind=CUSTOM_REAL) function integrate_complex_mesh_scalar(compl_scal)
        ! Integrates a complex scalar, defined at each GLL point, over the mesh
        use params, only: ngllx, nglly, ngllz, nspec, detjac, wgll
        implicit none 

        complex(kind=CUSTOM_REAL) :: compl_scal(ngllx, nglly, ngllz, nspec) 
        complex(kind=CUSTOM_REAL) :: sum
        integer :: i, j, k, ispec

        sum = (zero, zero) 
        do ispec = 1, nspec 
            do i = 1, ngllx 
                do j = 1, nglly
                    do k = 1, ngllz 
                        sum = sum + compl_scal(i, j, k, ispec) * detjac(i, j, k, ispec) * wgll(i) * wgll(j) * wgll(k)
                    enddo 
                enddo 
            enddo 
        enddo 

        integrate_complex_mesh_scalar = sum 
        return
    end function  integrate_complex_mesh_scalar


    real(kind=SPLINE_REAL) function integrate_r_traps(r, f, n)
        ! integrates a radial function f(r) using a trapesoid method

        implicit none
        include "constants.h"

        ! IO variables
        integer :: n
        real(kind=CUSTOM_REAL) r(n)
        real(kind=SPLINE_REAL) f(n)

        ! Local variables
        integer :: i 
        real(kind=SPLINE_REAL) sum, dr(n-1)
        
        if(n .lt. 2)then 
            write(*,*)'Error: must be at least 2 points to integrate'
        elseif(n.eq.2)then 
            sum = real(r(2)-r(1),  kind=SPLINE_REAL) *  (f(1)+f(2))/SPLINE_TWO 
        else
            dr(:) = real(r(2:n)-r(1:n-1), kind=SPLINE_REAL)

            sum = zero
            do i = 1, n-1
                sum = sum +  dr(i)*(f(i) + f(i+1))/SPLINE_TWO
            enddo 
        endif

        integrate_r_traps = sum
end function integrate_r_traps

end module integrate 