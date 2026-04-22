program vcen_woodhouse 
! Computes V centrifugal matrix using woodhouse kernel for testing against
! SEM method. 

use params, only: rho_spl, Vcen
use Integrate, only: integrate_r_traps
use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
use modes, only: Mode, get_mode
use mineos_model, only: mineos, mineos_ptr
use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
use w3j, only: thrj
use woodhouse_kernels, only: WK_Vphi, WK_Vphi_dot, BNpmlsld


implicit none
include "constants.h"

type(Mode)            :: mode_1
type(InterpPiecewise) :: interp

character :: t1, t2 
integer :: l1,  n1, mmin, im 

complex(kind=SPLINE_REAL), allocatable :: integrand(:), Vphi(:), Vdotphi(:)
complex(SPLINE_REAL), allocatable :: W_s(:)
integer :: i, j, k,  m1,  npoints, knot_lower, knot_upper 
real(kind=CUSTOM_REAL) :: r_lower, r_upper, sameone, mf, om2, tl1
complex(SPLINE_REAL) :: int_Ws, intvphis, sum



! Read mineos model 
call mineos%process_mineos_model(.false.)
mineos_ptr => mineos

n1 = 0
t1 = 'S'
l1 = 2


mode_1 = get_mode(n1, t1, l1 , mineos_ptr)


! Allocate the matrix: 
allocate(Vcen(mode_1%tl1, mode_1%tl1))

Vcen = SPLINE_iZERO

! We want to interpolate the radial eigenfunctions to make things a bit 
! more accurate 

! Values for the inner core
! knot_lower = 1
! r_lower    = zero        
! knot_upper = mineos%disc(2)
! r_upper    = mineos%rdisc(2)
! npoints    = 1000*(knot_upper-knot_lower)


! Values for the whole Earth
knot_lower = 1
r_lower    = zero        
knot_upper = mineos%disc(2)
r_upper    = mineos%rdisc(2)
npoints    = 10*(knot_upper-knot_lower)


! Create interpolator with evenly spaced points in IC 
interp = create_PieceInterp(npoints)
interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 
call interp%setup()
call interp%create_interpolation_radial_map()


! Interpolate mode splines
call interp%interpolate_mode_eigenfunctions(mode_1)

! We also need the density: 
allocate(rho_spl(npoints))
call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)


! Compute Ws (D.70)
allocate(W_s(npoints))

W_s = (mode_1%v_spl/mode_1%kf)**two  + two * mode_1%u_spl * mode_1%v_spl/mode_1%kf 

integrand = W_s * rho_spl * interp%radial * interp%radial
! Now we need to integrate for rho Ws r^2 
int_Ws =  integrate_r_traps(interp%radial, integrand, npoints)


! Compute integral of the two VPhi woodhouse kernels for s=2
allocate(Vdotphi(npoints))
allocate(Vphi(npoints))



! Compute integrand: 
integrand = rho_spl * (interp%radial)**three * & 
            (two* BNpmlsld(0, 1, l1, 2, l1) * (three * mode_1%u_spl * mode_1%u_spl/interp%radial - two*mode_1%u_spl * mode_1%aux_f ) +        &
             BNpmlsld(1, 1, l1, 2, l1)*( mode_1%u_spl*mode_1%dv_spl/mode_1%kf - mode_1%du_spl*mode_1%v_spl/mode_1%kf +                        & 
                                       three*mode_1%u_spl*mode_1%v_spl/(mode_1%kf*interp%radial) - two*mode_1%aux_f*mode_1%v_spl/mode_1%kf) & 
             )

integrand(1) = SPLINE_iZERO

intVphis =  integrate_r_traps(interp%radial, integrand, npoints) * OMEGA * OMEGA / three

! Find the smallest of the two l's: 
mmin = l1 




do im = -mmin, mmin

    sum = SPLINE_iZERO
    mf = real(im, kind=SPLINE_REAL)

    sum = sum + OMEGA*OMEGA*(TWO/THREE)*(one - (mode_1%kf**two)*int_Ws)

    sum = sum + ((-one)**(mf)) * thrj(l1, 2, l1, -im, 0, im) * mode_1%tl1 * intVphis

    Vcen(im+l1+1, im+l1+1) = Vcen(im+l1+1, im+l1+1) +  sum 

enddo !im 


call save_Vcen_matrix(l1, l1, './Vcen/vcen_woodhouse')



end program vcen_woodhouse