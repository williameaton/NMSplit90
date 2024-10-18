module integrate 
    use specfem_mesh, only: SetMesh
    implicit none 
    include "constants.h"


    interface integrate_over_mesh
        module procedure integrate_real_mesh_scalar     
        module procedure integrate_complex_mesh_scalar
    end interface


    interface integrate_r_traps
        module procedure integrate_r_traps_real_4
        module procedure integrate_r_traps_real_8
        module procedure integrate_r_traps_complex_4
        module procedure integrate_r_traps_complex_8
    end interface integrate_r_traps


    contains

     real(kind=CUSTOM_REAL) function integrate_real_mesh_scalar(SM, scalar)
        ! Integrates a real scalar, defined at each GLL point, over the mesh
        implicit none 
        ! IO variables: 
        type(SetMesh) :: SM 
        real(kind=CUSTOM_REAL) :: scalar(SM%ngllx, SM%nglly, SM%ngllz, SM%nspec)
        ! Local 
        real(kind=CUSTOM_REAL) :: sum
        integer :: i, j, k, ispec

        sum = zero 
        do ispec = 1, SM%nspec 
            do i = 1, SM%ngllx 
                do j = 1, SM%nglly
                    do k = 1, SM%ngllz 
                        sum = sum + scalar(i, j, k, ispec) * SM%detjac(i, j, k, ispec) * SM%wgll(i) * SM%wgll(j) *SM%wgll(k)
                    enddo 
                enddo 
            enddo 
        enddo 
        integrate_real_mesh_scalar = sum 
        return
    end function  integrate_real_mesh_scalar

    
    complex(kind=CUSTOM_REAL) function integrate_complex_mesh_scalar(SM, compl_scal)
        ! Integrates a complex scalar, defined at each GLL point, over the mesh
        implicit none 
        type(SetMesh) :: SM 
        complex(kind=CUSTOM_REAL) :: compl_scal(SM%ngllx, SM%nglly, SM%ngllz, SM%nspec) 
        complex(kind=CUSTOM_REAL) :: sum
        integer :: i, j, k, ispec

        sum = (zero, zero) 
        do ispec = 1, SM%nspec 
            do i = 1, SM%ngllx 
                do j = 1, SM%nglly
                    do k = 1, SM%ngllz 
                        sum = sum + compl_scal(i, j, k, ispec) * SM%detjac(i, j, k, ispec) * SM%wgll(i) * SM%wgll(j) * SM%wgll(k)
                    enddo 
                enddo 
            enddo 
        enddo 
        integrate_complex_mesh_scalar = sum 
        return
    end function  integrate_complex_mesh_scalar


    real(kind=4) function integrate_r_traps_real_4(r, f, n)
        ! integrates a radial function f(r) using a trapesoid method

        implicit none
        include "constants.h"

        ! IO variables
        integer :: n
        real(kind=8) r(n)
        real(kind=4) f(n)

        ! Local variables
        integer :: i 
        real(kind=4) sum, dr(n-1)
        
        if(n .lt. 2)then 
            write(*,*)'Error: must be at least 2 points to integrate'
        elseif(n.eq.2)then 
            sum = real(r(2)-r(1),  kind=SPLINE_REAL) *  (f(1)+f(2))/SPLINE_TWO 
        else
            dr(:) = real(r(2:n)-r(1:n-1), kind=SPLINE_REAL)

            sum = SPLINE_ZERO
            do i = 1, n-1
                sum = sum +  dr(i)*(f(i) + f(i+1))/SPLINE_TWO
            enddo 
        endif

        integrate_r_traps_real_4 = sum
end function integrate_r_traps_real_4



real(kind=8) function integrate_r_traps_real_8(r, f, n)
    ! integrates a radial function f(r) using a trapesoid method

    implicit none
    include "constants.h"

    ! IO variables
    integer :: n
    real(kind=8) r(n)
    real(kind=8) f(n)

    ! Local variables
    integer :: i 
    real(kind=8) sum, dr(n-1)

    if(n .lt. 2)then 
        write(*,*)'Error: must be at least 2 points to integrate'
    elseif(n.eq.2)then 
        sum = real(r(2)-r(1),  kind=8) *  (f(1)+f(2))/TWO 
    else
        dr(:) = real(r(2:n)-r(1:n-1), kind=8)

        sum = zero
        do i = 1, n-1
            sum = sum +  dr(i)*(f(i) + f(i+1))/TWO
        enddo 
    endif

    integrate_r_traps_real_8 = sum
end function integrate_r_traps_real_8



complex(kind=4) function integrate_r_traps_complex_4(r, f, n)
! integrates a radial function f(r) using a trapesoid method

implicit none
include "constants.h"

! IO variables
integer :: n
real(kind=8) r(n)
complex(kind=4) f(n)

! Local variables
integer :: i 
real(kind=4) dr(n-1)
complex(kind=4) sum

if(n .lt. 2)then 
    write(*,*)'Error: must be at least 2 points to integrate'
elseif(n.eq.2)then 
    sum = real(r(2)-r(1),  kind=4) *  (f(1)+f(2))/SPLINE_TWO 
else
    dr(:) = real(r(2:n)-r(1:n-1), kind=4)

    sum = zero
    do i = 1, n-1
        sum = sum +  dr(i)*(f(i) + f(i+1))/SPLINE_TWO
    enddo 
endif

integrate_r_traps_complex_4 = sum
end function integrate_r_traps_complex_4




complex(kind=8) function integrate_r_traps_complex_8(r, f, n)
! integrates a radial function f(r) using a trapesoid method

implicit none
include "constants.h"

! IO variables
integer :: n
real(kind=8) r(n)
complex(kind=8) f(n)

! Local variables
integer :: i 
real(kind=8) dr(n-1)
complex(kind=8) sum

if(n .lt. 2)then 
    write(*,*)'Error: must be at least 2 points to integrate'
elseif(n.eq.2)then 
    sum = (r(2)-r(1)) *  (f(1)+f(2))/TWO 
else
    dr(:) = real(r(2:n)-r(1:n-1), kind=8)

    sum = zero
    do i = 1, n-1
        sum = sum +  dr(i)*(f(i) + f(i+1))/TWO
    enddo 
endif

integrate_r_traps_complex_8 = sum
end function integrate_r_traps_complex_8


end module integrate 