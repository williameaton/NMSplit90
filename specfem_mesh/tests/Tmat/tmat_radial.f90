! Computes Tmat from the radial integral expressions of D.43: 
program tmat_radial 
use params, only: Tmat
use w3j, only: thrj
use modes, only: Mode, get_mode
use rho_st_profiles, only: radial_rho_st
use mineos_model, only: mineos, mineos_ptr
use splitting_function, only: xi_kkst
use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
use integrate, only: integrate_r_traps
use woodhouse_kernels, only: WK_Trho
use splitting_function, only: get_Ssum_bounds, Hcomplex_to_cst, write_cst_complex_to_file

implicit none 
include "constants.h"

integer   :: rho_smax, n1, l1, n2, l2, s, it, t, m1, m2, im1, im2, j, i , thjval, num_s, ncols
integer   :: knot_lower, knot_upper, npoints, smin, smax
real(kind=CUSTOM_REAL) :: r_lower, r_upper
character :: t1, t2 
character(len=2) :: n1str,n2str, l1str, l2str
character(len=300) :: outfile, out_name
complex(kind=SPLINE_REAL), allocatable :: cst(:,:)

complex(kind=CUSTOM_REAL) :: rad_int
complex(kind=CUSTOM_REAL), allocatable :: rho_st(:,:), rhost_radarray(:)
type(Mode)                :: mode_1, mode_2    
type(InterpPiecewise)     :: interp
complex(kind=SPLINE_REAL), allocatable :: Tp(:), integrand(:)

logical, parameter :: only_IC = .true.

! First mode: 
n1 = 0
t1 = 'S'
l1 = 3

n2 = 0
t2 = 'S'
l2 = 2

if(t1.eq.t2)then
    thjval = 0
else 
    thjval = 1
endif 



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
    r_upper    = SCALE_R   !mineos%rdisc(2)
endif

npoints    = 100*(knot_upper-knot_lower)


! Create interpolator with evenly spaced points in IC 
interp = create_PieceInterp(npoints)
interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/SCALE_R, j = 1, npoints)] 


write(*,*)'Max radial value: ', maxval(interp%radial)


! Setup the rho_st values for this example: 
! Note here that s should start at s = 0 
! FOR NOW IT WILL BE A CONSTANT VALUE WITH RADIUS
rho_smax = 5
allocate(rho_st(rho_smax, 2*rho_smax+1))

rho_st(1, 1:1) = (/0.0d0/)    ! s=0
rho_st(2, 1:3) = (/(0.4d0, -0.2d0), (-0.14d0, 0.0d0),  (-0.4d0, 0.2d0) /)    ! s=1
rho_st(3, 1:5) = (/ 0.7d0,  -0.43d0, -0.63d0,   0.43d0, 0.7d0/)  ! s=2
rho_st(4, 1:7) = (/0.22d0,  0.1d0,   0.09d0,   0.0d0, -0.09d0, 0.1d0, -0.22d0/)  ! s=3
rho_st(5, 1:9) = (/0.22d0, 0.02d0,  0.1d0,   0.09d0,   0.0d0, -0.09d0, 0.1d0, -0.02d0, 0.22d0/)  ! s=4


write(*,*)rho_st(2,:)

call interp%setup()
call interp%create_interpolation_radial_map()

! Load the first and second modes
mode_1 = get_mode(n1, t1, l1, mineos_ptr)
mode_2 = get_mode(n2, t2, l2, mineos_ptr)

! Interpolate mode splines
call interp%interpolate_mode_eigenfunctions(mode_1)
call interp%interpolate_mode_eigenfunctions(mode_2)


! Allocate Tmat 
allocate(Tmat(mode_1%tl1, mode_2%tl1))
Tmat = SPLINE_iZERO

allocate(Tp(npoints))
allocate(integrand(npoints))
allocate(rhost_radarray(npoints))

smin = abs(l2-l1)
smax = l2+l1

if(smin.eq.0)then
    smin=2
endif 


do m1 = -l1, l1 
    do m2 = -l2, l2 
        
        ! Matrix indices
        im1 = m1 + l1 + 1
        im2 = m2 + l2 + 1

        ! Loop over s and t: 
        ! sum over the s
        do s = smin, smax, 2
            ! Compute the kernel Tp for this s 
            Tp = zero
            call WK_Trho(mode_1, mode_2, s, Tp)

            do it = 1, 2*s + 1
                t = it - 1 - s
                
                ! Compute a radial profile for this rho_st: 
                rhost_radarray = SPLINE_iZERO
                call radial_rho_st(s, npoints, rho_st(s+1, it), rhost_radarray, interp%radial)

                ! Save the rho_st profile (real part only)
                if(t.lt.0)then
                    write(outfile,'(a,i1,a,i2)')'./Tmat/rho_st_profiles/rhost_',s,'_',t
                else 
                    write(outfile,'(a,i1,a,i1)')'./Tmat/rho_st_profiles/rhost_',s,'_',t
                endif

                open(1,file=trim(outfile))
                write(1,*) real(rho_st(s, it), kind=SPLINE_REAL), aimag(rho_st(s, it))
                do i =1, npoints
                    write(1,*)interp%radial(i), real(rhost_radarray(i), kind=SPLINE_REAL)
                enddo 
                close(1)

                integrand =  Tp * rhost_radarray * interp%radial * interp%radial
                rad_int   =  integrate_r_traps(interp%radial, integrand, npoints)

                Tmat(im1,im2) = Tmat(im1,im2) + & 
                                rad_int * xi_kkst(m1,m2,s,t,l1,l2)
                
                ! Independent of m so write onces
                if(m1.eq.0.and.m2.eq.0) write(*,*)'cst: ', s, t, rad_int/thrj(l1+thjval, s+thjval, l2+thjval, 0, 0, 0)

            enddo ! it
        enddo !is 

    enddo ! m2 
enddo ! m1

write(outfile, '(a,i1,a,i1,a,i1,a,i1,a)')'Tmat/radial_', n1, t1, l1, '_', n2, t2, l2, '.txt'
call save_T_matrix(l1, l2,  outfile)


! Convert the t matrix to Csts separately: 
! Write as a CST
call buffer_int(n1str, n1)
call buffer_int(n2str, n2)
call buffer_int(l1str, l1)
call buffer_int(l2str, l2)

call get_Ssum_bounds(l1, l2, smin, smax, num_s, ncols)
allocate(cst(num_s, ncols))
call Hcomplex_to_cst(Tmat, l1, l2, cst, ncols, num_s, t1, t2, 2)
out_name = 'Tmat/cst_from_mat_'//trim(n1str)//t1//trim(l1str)//'_'//trim(n2str)//t2//trim(l2str)
call write_cst_complex_to_file(out_name, cst, ncols, num_s, smin, 2)
deallocate(cst)


end program tmat_radial