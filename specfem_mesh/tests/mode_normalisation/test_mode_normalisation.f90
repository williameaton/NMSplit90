
program test_mode_norm
    use params, only: NL, IC_ID, ndisc, rdisc, disc, rad_mineos, rho_mineos, u_spl, v_spl,udot_spl, vdot_spl, interp_map, rho_spl
    use spline, only: create_interpolation_radial_map, interpolate_mode_eigenfunctions, write_mode_spline, interpolate_mineos_variable, quad_spline_interp_3
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    implicit none
    include "constants.h"

    integer :: i, j, l, n , lentrim, lenl, lenn
    character(len=1)   :: char, type
    character(len=10)  :: mode
    character(len=250) :: arg
    character(len=30)  :: out_name
    character(len=4), parameter :: fmti1 = "(i1)"
    character(len=4), parameter :: fmti2 = "(i2)"
    character(len=4), parameter :: fmti3 = "(i3)"
    character(len=4) :: fmt_n, fmt_l
    real(SPLINE_REAL) :: wcom, qmod
    real(SPLINE_REAL), allocatable :: u(:), v(:), du(:), dv(:)
    integer :: npoints, knot_lower, knot_upper
    logical :: found
    real(kind=CUSTOM_REAL), allocatable :: radial_vals(:), r_lower, r_upper

    real(SPLINE_REAL), allocatable :: integrand(:)
    real(SPLINE_REAL) :: total_integral, sum,  precision


    
    ! Read mineos model 
    call process_mineos_model()
    allocate(u(NL), v(NL), du(NL), dv(NL))

    ! Now we can load the mode: 
    ! Choose a mode: 
    
    type = 'T'
    do l = 1, 12
        do n = 1, 12


            call get_mode(type, n, l,wcom,qmod, u, du, v, dv, .true.)

            total_integral = SPLINE_ZERO

            do i  = 1, ndisc + 1
                ! Sort out the first and last sections (above, below the first/last continuity)
                if (i == 1) then 
                    knot_lower = 1
                    r_lower = zero        
                else
                    knot_lower =  disc(i-1)+1
                    r_lower    = rdisc(i-1)
                endif 

                if (i == ndisc+1) then 
                    knot_upper =  NL
                    r_upper    = rad_mineos(NL)*SCALE_R
                else
                    knot_upper =  disc(i)
                    r_upper    = rdisc(i)
                endif

                ! Spline has 2x the number of knots in each section 
                npoints = 6*(knot_upper-knot_lower)
                allocate(radial_vals(npoints), interp_map(npoints), integrand(npoints))


                call deallocate_if_allocated(u_spl)
                call allocate_if_unallocated(npoints, u_spl)
                if(type.eq.'S')then
                    call deallocate_if_allocated(v_spl)
                    call allocate_if_unallocated(npoints, v_spl)
                endif 

                ! Create radial array where eigenfunction will be interpolated
                do j = 1, npoints
                    radial_vals(j) = (r_lower +  (real(j-1)/real(npoints-1))*(r_upper-r_lower))/scale_R
                enddo 

                ! Interpolate
                call deallocate_if_allocated(u_spl)
                call allocate_if_unallocated(npoints, u_spl)
                call deallocate_if_allocated(udot_spl)
                call allocate_if_unallocated(npoints, udot_spl)
                call deallocate_if_allocated(v_spl)
                call allocate_if_unallocated(npoints, v_spl)
                call deallocate_if_allocated(vdot_spl)
                call allocate_if_unallocated(npoints, vdot_spl)

                call create_interpolation_radial_map(radial_vals, interp_map, npoints, knot_lower, knot_upper)
                call interpolate_mode_eigenfunctions(type, u, v, du, dv, knot_lower, knot_upper, &  
                                                    radial_vals, npoints, interp_map, & 
                                                    u_spl, v_spl, udot_spl, vdot_spl)

                ! write with custom name
                if (i .lt.10)then
                    write(out_name,'(a, i1, a)')'spline_', i,'.txt'
                else 
                    write(out_name,'(a, i2, a)')'spline_', i,'.txt'
                endif
                !call write_mode_spline(n, type, l, radial_vals, npoints, out_name)


                ! Interpolate the density for this region: 
                allocate(rho_spl(npoints))
                call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), knot_lower, knot_upper, & 
                                                radial_vals, npoints, rho_spl, interp_map)



                ! Save the rho spline
                if (i .lt.10)then
                    write(out_name,'(a, i1, a)')'rhospline_', i,'.txt'
                else 
                    write(out_name,'(a, i2, a)')'rhospline_', i,'.txt'
                endif
                open(1,file=trim(out_name))
                do j =1, npoints
                    write(1,*)radial_vals(j), rho_spl(j)
                enddo 
                close(1)
                

                ! Compute mode integrand a la MINEOS
                ! Eigenfunction normalisation from MINEOS: 
                !  \int_0^a rho(r) W^2 r^2    1 / omega^2 
                if (type.eq.'S')then
                    integrand = rho_spl * (u_spl*u_spl + v_spl*v_spl) * radial_vals * radial_vals
                elseif(type.eq.'T')then
                    integrand = rho_spl * (u_spl*u_spl) * radial_vals * radial_vals
                else
                    write(*,*)'Only for T,S rn.'
                    stop
                endif 
                sum =  integrate_r_traps(radial_vals, integrand, npoints)
                total_integral = total_integral + sum

                deallocate(radial_vals, interp_map, rho_spl, integrand)                      
            enddo 

            total_integral = total_integral *  wcom*SCALE_T * wcom*SCALE_T
            

            precision = 1e-2
            if(abs(total_integral-SPLINE_ONE).gt.precision)then 
                write(*,*)'Integral not close enough to 1'
                write(*,'(a,i2, a, i2)')'Mode : ', n, type, l
                write(*,*)'Total integral: ', total_integral
                stop
            endif

        enddo ! n 
    enddo ! l


end program test_mode_norm