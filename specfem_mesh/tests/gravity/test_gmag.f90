program test_g

    use params, only: rho_spl
    use mineos_model, only: mineos, mineos_ptr
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use gravitation, only: compute_background_g
    implicit none 
    include "constants.h"

    type(InterpPiecewise) :: interp
    integer :: knot_lower, knot_upper, npoints, j
    real(kind=CUSTOM_REAL) :: r_lower, r_upper
    real(kind=CUSTOM_REAL), allocatable :: gmag(:)

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos

    knot_lower = 1
    r_lower    = zero   
    knot_upper = mineos%NR !mineos%disc(2)
    r_upper    = SCALE_R !mineos%rdisc(2)
    npoints    = 100*(knot_upper-knot_lower)
    interp = create_PieceInterp(npoints)
    interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/SCALE_R, j = 1, npoints)] 

    call interp%setup()
    call interp%create_interpolation_radial_map()

    allocate(rho_spl(npoints))
    call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)


    allocate(gmag(npoints))

    call compute_background_g(interp%radial, rho_spl, npoints, gmag)

    open(1,file='gravity/gmag.txt', form='formatted')
    do j = 1, npoints
        write(1,*)interp%radial(j)*SCALE_R, gmag(j)*SCALE_R/(SCALE_T*SCALE_T)
    enddo 

    close(1)

end program test_g
