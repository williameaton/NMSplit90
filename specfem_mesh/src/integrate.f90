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






end module integrate 