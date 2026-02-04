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
use woodhouse_kernels, only: WK_Vphi, WK_Vphi_dot


implicit none
include "constants.h"

type(Mode)            :: mode_1, mode_2    
type(InterpPiecewise) :: interp

character :: t1, t2 
integer :: l1, l2, n1, n2, mmin, im 

complex(kind=SPLINE_REAL), allocatable :: integrand(:), Vphi(:), Vdotphi(:)
complex(SPLINE_REAL), allocatable :: W_s(:)
integer :: i, j, k,  m1, m2, npoints, knot_lower, knot_upper 
real(kind=CUSTOM_REAL), allocatable :: r_lower, r_upper, sameone, mf, om2, half_tl12
complex(SPLINE_REAL) :: int_Ws, intvphis, sum



! Read mineos model 
call mineos%process_mineos_model(.false.)
mineos_ptr => mineos

n1 = 0
t1 = 'S'
l1 = 2

n2 = 0
t2 = 'S'
l2 = 2

mode_1 = get_mode(n1, t1, l1 , mineos_ptr)
mode_2 = get_mode(n2, t2, l2,  mineos_ptr)


! Allocate the matrix: 
allocate(Vcen(mode_1%tl1, mode_2%tl1))

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
knot_upper = mineos%NR
r_upper    = scale_R
npoints    = 100*(knot_upper-knot_lower)


! Create interpolator with evenly spaced points in IC 
interp = create_PieceInterp(npoints)
interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R, j = 1, npoints)] 
call interp%setup()
call interp%create_interpolation_radial_map()


! Interpolate mode splines
call interp%interpolate_mode_eigenfunctions(mode_1)
call interp%interpolate_mode_eigenfunctions(mode_2)

! We also need the density: 
allocate(rho_spl(npoints))
call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)


! Compute Ws (D.70)
allocate(W_s(npoints))
if (mode_1%t.ne.mode_2%t)then 
    W_s = SPLINE_ZERO
elseif(mode_1%t.eq.'S' .and. mode_2%t.eq.'S')then 
    W_s = mode_1%v_spl/mode_1%kf * mode_2%v_spl/mode_2%kf & 
        + mode_1%u_spl * mode_2%v_spl/mode_2%kf & 
        + mode_2%u_spl * mode_1%v_spl/mode_1%kf
    
elseif(mode_1%t.eq.'T' .and. mode_2%t.eq.'T')then
    W_s = mode_1%w_spl/mode_1%kf * mode_2%w_spl/mode_2%kf
else
    write(*,*)'Error in mode type', mode_1%t, mode_2%t
    stop
endif 

integrand = W_s * rho_spl * interp%radial * interp%radial
! Now we need to integrate for rho Ws r^2 
int_Ws =  integrate_r_traps(interp%radial, integrand, npoints)


! Compute integral of the two VPhi woodhouse kernels for s=2
allocate(Vdotphi(npoints))
allocate(Vphi(npoints))

call WK_Vphi(mode_1, mode_2, 2, interp%radial, rho_spl, Vphi)
call WK_Vphi_dot(mode_1, mode_2, 2, interp%radial, rho_spl, Vdotphi)

! Compute integrand: 
integrand = (Vphi * interp%radial  + two * Vdotphi )* (interp%radial**three) / three

intVphis =  integrate_r_traps(interp%radial, integrand, npoints) * OMEGA * OMEGA 

! Find the smallest of the two l's: 
mmin = min(l1, l2) 

! because of the delta_mm' we only need the diagonal 

! Delta_sigmasigma' Delta_nn'
if(t1.eq.t2 .and. n1.eq.n1)then
    sameone = one
else 
    sameone = zero
endif

!write(*,*)'sameone: ', sameone
!write(*,*)'int_Ws: ', int_Ws
!write(*,*)'intVphis: ', intVphis
!write(*,*)'2/3 OMEGA^2: ', TWO/THREE* OMEGA*OMEGA

half_tl12 = (real(mode_1%tl1)*real(mode_2%tl1))**(half)


do im = -mmin, mmin

    sum = SPLINE_iZERO
    mf = real(im, kind=SPLINE_REAL)

    if(l1.eq.l2)then 
        sum = sum + OMEGA*OMEGA*(TWO/THREE)*( sameone - (mode_1%kf**two)*int_Ws)
    endif 
    sum = sum + ((-one)**(mf)) * thrj(l1, 2, l2, -im, 0, im)  * half_tl12 * intVphis

    Vcen(im+l1+1, im+l2+1) = Vcen(im+l1+1, im+l2+1) +  sum 

enddo !im 


call save_Vcen_matrix(l1, l2, './Vcen/vcen_woodhouse')



end program vcen_woodhouse