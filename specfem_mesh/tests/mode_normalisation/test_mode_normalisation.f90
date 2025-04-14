
program test_mode_norm
    ! Tests the mode normalisation following the convention of DT98
    ! NMSplit almost always uses the notation of DT98, but they have different
    ! normalisation to MINEOS, where the Eigenfunctions are derived from: 
    !
    ! MINEOS - see manual section 3.1:
    !   - Toroidal  :  \omega^2 \int_0^R \rho W^2(r) r^2 dr = 1 
    !   - Spheroidal:  \omega^2 \int_0^R \rho [U^2 + V^2] r^2 dr = 1
    !
    ! DT98: - Equations 8.107 - 8.108
    !   - Toroidal   : \int_0^R \rho W^2 r^2 dr = 1 
    !   - Spheroidal : \int_0^R \rho ( U^2 + V^2 ) r^2 dr = 1 
    ! In the function Get_Mode we multiply the MINEOS eigenfunctions 
    ! by omega each so that any combination would be omega^2
    use mineos_model, only: MineosModel, mineos, mineos_ptr
    use modes, only: Mode, get_mode
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    use params, only: rho_spl
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    implicit none
    include "constants.h"

    integer :: l, n, npoints, j, itype
    character(len=1)   :: type
    character(len=80)  :: out_name
    real(SPLINE_REAL), allocatable :: integrand(:)
    real(SPLINE_REAL) :: total_integral, precision

    type(InterpPiecewise) :: interp
    type(Mode)            :: mode_1

    ! Read mineos model
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos
    npoints = mineos%NR*100

    ! Create and setup the interpolation object
    interp = create_PieceInterp(npoints)
    ! Fill up the array
    interp%radial = [(mineos%min_r + (real(j-1)/real(npoints-1))*(mineos%max_r -mineos%min_r), j = 1, npoints)] 
    call interp%setup()
    call interp%create_interpolation_radial_map()

    ! Interpolate the density for this region: 
    allocate(rho_spl(npoints))
    call interp%interpolate_mineos_variable(real(mineos%rho_mineos, kind=SPLINE_REAL), rho_spl)
    ! save rhospline
    !write(out_name,'(a, i1, a)')'mode_normalisation/rhospline.txt'
    !open(1,file=trim(out_name))
    !do j =1, npoints
    !    write(1,'(2e18.8)')interp%radial(j), rho_spl(j)
    !enddo 
    !close(1)

    allocate(integrand(npoints))

    do itype = 1, 2 

        ! Setup
        do n = 1,15
        !n       = 20
            if(itype.eq.1)then
                type    = 'T'
            elseif(itype.eq.2)then 
                type    = 'S'
            endif 
            l       = 6

            ! Load the mode
            mode_1 = get_mode(n, type, l, mineos_ptr)

            call interp%interpolate_mode_eigenfunctions(mode_1)                       

            ! Compute mode integrand a la MINEOS
            ! Eigenfunction normalisation from MINEOS: 
            !  \int_0^a rho(r) W^2 r^2    
            if (type.eq.'S')then
                integrand = mode_1%u_spl*mode_1%u_spl + mode_1%v_spl * mode_1%v_spl 
            elseif(type.eq.'T')then
                integrand = mode_1%w_spl * mode_1%w_spl
            else
                write(*,*)'Only for T,S'
                stop
            endif 
            integrand = rho_spl * integrand *  interp%radial *  interp%radial

            ! Integrate and dimensionalise
            total_integral =  integrate_r_traps(interp%radial, integrand, npoints)    
            !write(*,*)'Mode normalisation total integral: ', total_integral
            precision = 1e-2
            if(abs(total_integral-SPLINE_ONE).gt.precision .or. total_integral.ne.total_integral)then 
                write(*,*)'Integral not close enough to 1'
                write(*,'(a,i2, a, i2)')'Mode : ', n, type, l
                write(*,*)'Total integral: ', total_integral
                ! -1 status output
                open(2,file=trim('./mode_normalisation/status.txt'))
                write(2,'(i2)')-1
                close(2)

                stop
            endif
        enddo ! n 
    enddo 

    ! If we made it this far give a 0 status output
    open(2,file=trim('./mode_normalisation/status.txt'))
    write(2,'(i1)')0
    close(2)


end program test_mode_norm