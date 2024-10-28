! Computes Viso from the radial integral expressions of D.44: 
program viso_radial 
use params, only: Viso, rho_spl, g_spl
use w3j, only: thrj
use modes, only: Mode, get_mode
use rho_st_profiles, only: radial_rho_st, phist_from_rhost, Gradphist_from_rhost
use mineos_model, only: mineos, mineos_ptr
use splitting_function, only: xi_kkst
use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
use integrate, only: integrate_r_traps
use woodhouse_kernels, only: WK_Vkappa, WK_Vmu, WK_Vrho, WK_Vphi, WK_Vphi_dot
use gravitation, only: compute_background_g
implicit none 
include "constants.h"

integer   :: kap_smax, mu_smax, rho_smax, n1, l1, n2, l2, s, it, t, m1, m2, im1, im2, j, i 
integer   :: knot_lower, knot_upper, npoints
real(kind=CUSTOM_REAL) :: r_lower, r_upper
character :: t1, t2 
character(len=300) :: outfile
complex(kind=CUSTOM_REAL) :: rad_int
complex(kind=CUSTOM_REAL), allocatable :: kap_st(:,:), kapst_radarray(:)
complex(kind=CUSTOM_REAL), allocatable :: mu_st(:,:), must_radarray(:)
complex(kind=CUSTOM_REAL), allocatable :: rho_st(:,:), rhost_radarray(:)
complex(kind=CUSTOM_REAL), allocatable :: Phist_radarray(:)
complex(kind=CUSTOM_REAL), allocatable :: GradPhist_radarray(:)
type(Mode)                :: mode_1, mode_2    
type(InterpPiecewise)     :: interp
complex(kind=SPLINE_REAL), allocatable :: Vk(:), Vmu(:), Vrho(:), integrand(:), Vphi(:), Vdotphi(:)

logical, parameter :: only_IC = .true.
logical, parameter :: USE_PHI_IN_RHO = .true.


! First mode: 
n1 = 0
t1 = 'S'
l1 = 2

n2 = 0
t2 = 'S'
l2 = 2

! Read mineos model 
call mineos%process_mineos_model(.false.)
mineos_ptr => mineos


! Values for the inner core
knot_lower = 1
r_lower    = zero    

if(only_IC)then 
    knot_upper = mineos%disc(2)
    r_upper    = mineos%rdisc(2)
else 
    knot_upper = mineos%NR !mineos%disc(2)
    r_upper    = SCALE_R !mineos%rdisc(2)
endif

npoints    = 100*(knot_upper-knot_lower)

! Create interpolator with evenly spaced points in IC 
interp = create_PieceInterp(npoints)
interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/SCALE_R, j = 1, npoints)] 

write(*,*)'Max radial value: ', maxval(interp%radial)

! Setup the kap_st values for this example: 
! Note here that s should start at s = 0 
! FOR NOW IT WILL BE A CONSTANT VALUE WITH RADIUS
kap_smax  = 3
mu_smax   = 3
rho_smax  = 3
allocate(kap_st(kap_smax, 2*kap_smax+1))
allocate(mu_st(mu_smax, 2*mu_smax+1))
allocate(rho_st(rho_smax, 2*rho_smax+1))

rho_st(1, 1:3) = (/0.4d0, -0.14d0,  -0.4d0 /)
rho_st(2, 1:5) = (/ 0.7d0,  -0.43d0, -0.63d0,   0.43d0, 0.7d0/)
rho_st(3, 1:7) = (/0.22d0,  0.1d0,   0.09d0,   0.0d0, -0.09d0, 0.1d0, -0.22d0/)

kap_st(1, 1:3) = (/0.2d0,   0.55d0,  -0.2d0 /)
kap_st(2, 1:5) = (/ 0.54d0,  0.0d0, -0.23d0,   0.0d0, 0.54d0/)
kap_st(3, 1:7) = (/0.43d0,  0.34d0,   -0.09d0,  1.0d0, 0.09d0, 0.34d0, -0.43d0/)

mu_st(1, 1:3) = (/ 0.6d0,  -0.2d0,  -0.6d0 /)
mu_st(2, 1:5) = (/ 0.14d0,  0.13d0,  0.47d0,   -0.13d0, 0.14d0/)
mu_st(3, 1:7) = (/0.52d0,  -0.1d0,   0.39d0,   0.2d0, -0.39d0, -0.1d0, -0.52d0/)


call interp%setup()
call interp%create_interpolation_radial_map()

! Load the first and second modes
mode_1 = get_mode(n1, t1, l1, mineos_ptr)
mode_2 = get_mode(n2, t2, l2, mineos_ptr)

! Interpolate mode splines
call interp%interpolate_mode_eigenfunctions(mode_1)
call interp%interpolate_mode_eigenfunctions(mode_2)

! Need to get the density spline and the gravity magnitude: 
allocate(rho_spl(npoints))
call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)

allocate(g_spl(npoints))
call compute_background_g(interp%radial, rho_spl, npoints, g_spl)

! Allocate Viso 
allocate(Viso(mode_1%tl1, mode_2%tl1))
Viso = SPLINE_iZERO

allocate(Vk(npoints))
allocate(Vmu(npoints))
allocate(Vrho(npoints))
if(.not.USE_PHI_IN_RHO)then
    allocate(Vdotphi(npoints))
    allocate(Vphi(npoints))
endif

allocate(integrand(npoints))

allocate(kapst_radarray(npoints))
allocate(must_radarray(npoints))
allocate(rhost_radarray(npoints))
allocate(Phist_radarray(npoints))
allocate(GradPhist_radarray(npoints))

do m1 = -l1, l1 
    do m2 = -l2, l2 
        
        ! Matrix indices
        im1 = m1 + l1 + 1
        im2 = m2 + l2 + 1

        ! Loop over s and t: 
        ! sum over the s
        do s = 1, kap_smax 
            ! Compute the kernel Vkappa and Vmu for this s 
            Vk   = zero
            Vmu  = zero
            Vrho = zero
            call WK_Vkappa(mode_1, mode_2, s, interp%radial, Vk)
            call WK_Vmu(mode_1, mode_2, s, interp%radial, Vmu)
            call WK_Vrho(mode_1, mode_2, s, interp%radial, rho_spl, g_spl, Vrho, USE_PHI_IN_RHO)


            if(.not.USE_PHI_IN_RHO)then
                call WK_Vphi(mode_1, mode_2, s,  interp%radial, rho_spl, Vphi)
                call WK_Vphi_dot(mode_1, mode_2, s, interp%radial, rho_spl, Vdotphi)
            endif 

            do it = 1, 2*s + 1
                t = it - 1 - s
                
                ! Compute a radial profile for this kap_st: 
                kapst_radarray = SPLINE_iZERO
                call radial_rho_st(s, npoints, kap_st(s, it), kapst_radarray, interp%radial)
                
                ! Compute for mu_st - for now use the rho profiles 
                must_radarray = SPLINE_iZERO
                call radial_rho_st(s, npoints, mu_st(s, it), must_radarray, interp%radial)

                ! Compute for mu_st - for now use the rho profiles 
                rhost_radarray = SPLINE_iZERO
                call radial_rho_st(s, npoints, rho_st(s, it), rhost_radarray, interp%radial)

                if(USE_PHI_IN_RHO)then
                    integrand =  ( Vmu  * must_radarray  +  & 
                                Vk   * kapst_radarray +  & 
                                Vrho * rhost_radarray)   &
                                * interp%radial * interp%radial
                else 
                    ! Explicitly have phi and phidot kernels 
                    Phist_radarray     = SPLINE_iZERO
                    GradPhist_radarray = SPLINE_iZERO
                    call phist_from_rhost(s,     interp%radial, rhost_radarray, Phist_radarray, npoints )
                    call Gradphist_from_rhost(s, interp%radial, rhost_radarray, GradPhist_radarray, npoints )


                    integrand =  ( Vmu  * must_radarray  +  & 
                                   Vk   * kapst_radarray +  & 
                                   Vrho * rhost_radarray +  &
                                   Vphi * phist_radarray +  & 
                                   Vdotphi * GradPhist_radarray) &
                                   * interp%radial * interp%radial
                endif 

                rad_int   =  integrate_r_traps(interp%radial, integrand, npoints)

                Viso(im1,im2) = Viso(im1,im2) + & 
                                rad_int * xi_kkst(m1,m2,s,t,l1,l2)
            enddo ! it
        enddo !is 

    enddo ! m2 
enddo ! m1

write(outfile, '(a,i1,a,i1,a,i1,a,i1,a)')'Viso/radial_', n1, t1, l1, '_', n2, t2, l2, '.txt'
call save_Viso_matrix(l1, l2,  outfile)



end program viso_radial