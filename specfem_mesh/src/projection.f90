! subroutine gll_strain_from_disp(disp, strain)
!     ! Uses xyz displacement (d) to compute strain 
!     use params, only: ngllx, nglly, ngllz, nspec, dgll, jacinv
!     implicit none 
!     include "constants.h"

!     complex(SPLINE_REAL) :: disp(3, ngllx, nglly, ngllz, nspec)
!     complex(SPLINE_REAL) :: strain(6, ngllx, nglly, ngllz, nspec)

!     ! Local: 
!     integer :: ispec, i, j, k, p, q, g, h
!     complex(SPLINE_REAL) :: sum_c(3), strn_tmp(3,3)

!     strain = SPLINE_iZERO


!         ! Compute gradient of the displacement
!     do ispec = 1, nspec
!         do i = 1, ngllx 
!             do j = 1, nglly
!                 do k = 1, ngllz

!                     ! Loop over strain elementsa
!                     do p = 1, 3
!                         do q = 1, 3

!                             sum_c(1) = SPLINE_iZERO
!                             sum_c(2) = SPLINE_iZERO
!                             sum_c(3) = SPLINE_iZERO

!                             do g = 1, ngllx
!                                 sum_c(1) = sum_c(1) + disp(q, g, j, k, ispec) * real(dgll(g, i), kind=SPLINE_REAL)
!                                 sum_c(2) = sum_c(2) + disp(q, i, g, k, ispec) * real(dgll(g, j), kind=SPLINE_REAL)
!                                 sum_c(3) = sum_c(3) + disp(q, i, j, g, ispec) * real(dgll(g, k), kind=SPLINE_REAL)
!                             enddo 

!                             strn_tmp(q,p) = sum_c(1) * jacinv(1,p,i,j,k,ispec) +   &  
!                                             sum_c(2) * jacinv(2,p,i,j,k,ispec) +   & 
!                                             sum_c(3) * jacinv(3,p,i,j,k,ispec) 
!                         enddo !q
!                     enddo !p

!                     ! Symmetric part of tensor
!                     strn_tmp(:,:) =  SPLINE_HALF * (strn_tmp + transpose(strn_tmp) )

!                     ! Save in voigt notation 
!                     do h = 1, 3
!                         strain(h, i, j, k, ispec) = strn_tmp(h,h)
!                     enddo 
!                     strain(4, i, j, k, ispec) = strn_tmp(2,3)
!                     strain(5, i, j, k, ispec) = strn_tmp(1,3)
!                     strain(6, i, j, k, ispec) = strn_tmp(1,2)


!                 enddo 
!             enddo 
!         enddo 
!     enddo




! end subroutine gll_strain_from_disp
