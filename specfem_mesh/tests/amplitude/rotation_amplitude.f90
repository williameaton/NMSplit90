
program test_amplitude
    use params, only: rho_spl, Wmat
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: delta_spline
    use modes, only: Mode, get_mode
    use mineos_model, only: mineos, mineos_ptr
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"

    integer :: i, j, k,  m1, m2, npoints, knot_lower, knot_upper 
    character(len=80)  :: out_name
    real(kind=CUSTOM_REAL), allocatable :: r_lower, r_upper
    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:), W_a(:)
    real(SPLINE_REAL) :: kmkm2, kpkm2, kmkp2, Sl1m, Sl2m, mf 
    real(SPLINE_REAL) :: chi, bparam,  scaling98
    character(len=6) :: modename

    type(Mode)            :: mode_1
    type(InterpPiecewise) :: interp

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos


    
    mode_1 = get_mode(1, 'S', 9, mineos_ptr)

    write(modename, '(i1, a, i1)')mode_1%n, mode_1%t, mode_1%l


    if (mode_1%t.eq.'S')then
        out_name = 'amplitude/eigens/'//trim(modename)//'_U'
        open(1,file=trim(out_name), form='formatted')
        write(*,*)'writing to '//trim(out_name)
        ! Save the eigenfunctions: 
        do i = 1, mineos%NR
            write(1,*) mineos%rad_mineos(i),  mode_1%u(i)
        enddo
        close(1)
        ! V eigenfunction
        out_name = 'amplitude/eigens/'//trim(modename)//'_V'
        open(1,file=trim(out_name), form='formatted')
        write(*,*)'writing to '//trim(out_name)
        ! Save the eigenfunctions: 
        do i = 1, mineos%NR
            write(1,*) mineos%rad_mineos(i),  mode_1%v(i)
        enddo 
        close(1)
    else 
        ! W eigenfunction
        out_name = 'amplitude/eigens/'//trim(modename)//'_W'
        open(1,file=trim(out_name), form='formatted')
        write(*,*)'writing to '//trim(out_name)
        ! Save the eigenfunctions: 
        do i = 1, mineos%NR
            write(1,*) mineos%rad_mineos(i),  mode_1%w(i)
        enddo 
        close(1)
    endif

    ! Density 
    ! V eigenfunction
    out_name = 'amplitude/eigens/rho'
    open(1,file=trim(out_name), form='formatted')
    write(*,*)'writing to '//trim(out_name)
    ! Save the eigenfunctions: 
    do i = 1, mineos%NR
        write(1,*) mineos%rad_mineos(i),  mineos%rho_mineos(i)
    enddo 


    close(1)



    ! Values for the whole Earth
    knot_lower = 1
    r_lower    = zero        
    knot_upper = mineos%NR
    r_upper    = one
    npoints    = 1000*(knot_upper-knot_lower)


    ! Create interpolator with evenly spaced points in IC 
    interp = create_PieceInterp(npoints)
    interp%radial = [((r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower)), j = 1, npoints)] 


    call interp%setup()
    call interp%create_interpolation_radial_map()

    ! Interpolate mode splines
    call interp%interpolate_mode_eigenfunctions(mode_1)

    ! We also need the density: 
    allocate(rho_spl(npoints))
    call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)

 
    ! compute chi (splititng parameter): 

    allocate(W_s(npoints))
    W_s = SPLINE_ZERO

    if(mode_1%t.eq.'S')then 
        W_s = (mode_1%v_spl/mode_1%kf)**two & 
            + two * (mode_1%u_spl) * (mode_1%v_spl/mode_1%kf) 
    
    else
        W_s = (mode_1%w_spl/mode_1%kf ) * (mode_1%w_spl/mode_1%kf )
    endif 

    W_s = W_s * rho_spl * interp%radial * interp%radial    
    ! Now we need to integrate for rho Ws r^2 
    chi =  integrate_r_traps(interp%radial, W_s, npoints)

    write(*,*)chi


    ! b parameter (dt98 14.56) is chi * omega / w0 
    ! but the OMEGA needs to be re-dimensionalised(?) and we non-dimensionalised it
    ! with time so need to now divide by time 

    ! Need to multiply by omega^2 (non dimensionalised to account for difference in DT98 norm with mineos?)
    scaling98 = (SCALE_T * mode_1%wcom)**two
    
    bparam =scaling98 * chi * OMEGA / (SCALE_T * mode_1%wcom)

    write(*,*)'Mode ', mode_1%n, mode_1%t, mode_1%l
    write(*,*)'-----------------------------------'
    write(*,*)'Angular freq: ', mode_1%wcom
    write(*,*)
    write(*,*)'b           : ', bparam * 1000.0d0


end program test_amplitude