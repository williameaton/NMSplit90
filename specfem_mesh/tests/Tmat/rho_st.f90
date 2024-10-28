! For this benchmark test we want to prescribe coefficients of delta rho st 
! that vary radially. In the test we are using smax = 3 so need radial 
! profiles for each rho_st of s = 1, 2, 3
! For simplicity lets assume that the t values of each s have the same radial 
! profile, which is multiplied by some constant
module rho_st_profiles 
    use math, only: sinp
    implicit none 
    include "constants.h"


    contains 

    subroutine radial_rho_st(s, npts, coeff, arr, radius)
        implicit none 
        integer :: s, npts, i
        complex(kind=SPLINE_REAL) :: coeff
        complex(kind=SPLINE_REAL) :: arr(npts)
        real(kind=CUSTOM_REAL)    :: radius(npts)


        if(s.eq.1)then 
            ! When s = 1 
            ! It will be quadratic, going through a value of +0.6 at each end and -0.4 at the centre 
            arr  =  coeff * ((six/ten) + (radius - half)**two )
        elseif(s.eq.2)then 
            do i = 1, npts
                arr(i) =  coeff * (sinp(TWO_PI*radius(i)))
            enddo 
        elseif(s.eq.3)then
            arr  =  coeff * ((seven/ten) - TWO*radius  )
        else
            write(*,*)'Only supporting s=1,3 but s = ', s
            stop 
        endif

    end subroutine



    ! Same as above but for a single r value: 
    complex(kind=SPLINE_REAL) function single_r_rho_st(s, coeff, radius)
        implicit none 
        integer :: s
        complex(kind=SPLINE_REAL) :: coeff
        real(kind=CUSTOM_REAL)    :: radius


        if(s.eq.1)then 
            ! When s = 1 
            ! It will be quadratic, going through a value of +0.6 at each end and -0.4 at the centre 
            single_r_rho_st =  coeff * ((six/ten) + (radius - half)**two )
        elseif(s.eq.2)then 
            single_r_rho_st =  coeff * (sinp(TWO_PI*radius)) 
        elseif(s.eq.3)then
            single_r_rho_st =  coeff * ((seven/ten) -TWO*radius  )
        else
            write(*,*)'Only supporting s=1,3 but s = ', s
            stop 
        endif

    end function single_r_rho_st


    subroutine phist_from_rhost(s, phi_rad, rhost_radarray, phist_radarray, npoints)
        ! DT98 Eqn. D.54 - computes phi st from rho st
        use integrate, only: integrate_r_traps

        implicit none 
        integer :: s, npoints, ir 
        real(kind=CUSTOM_REAL) :: phi_rad(npoints), fs, rr
        complex(kind=CUSTOM_REAL) ::  b_integrand(npoints), a_integrand(npoints),&
                                  rhost_radarray(npoints), & 
                                  phist_radarray(npoints)
        complex(kind=SPLINE_REAL) :: integral_below, integral_above

        fs = real(s, kind=CUSTOM_REAL)
                            
        b_integrand = phi_rad**(fs+two)*rhost_radarray
        a_integrand = phi_rad**(-fs+one)*rhost_radarray
                          

        do ir = 1, npoints
                rr = phi_rad(ir)

            if(ir.eq.1)then
                integral_below = SPLINE_iZERO
            else  
                integral_below = rr**(-fs-one) * integrate_r_traps(phi_rad(1:ir) , & 
                                                                b_integrand(1:ir), &
                                                                ir)
            endif 

            if(ir.eq.npoints)then 
                integral_above = SPLINE_iZERO
            else 
                integral_above = (rr**fs) * integrate_r_traps(phi_rad(ir:npoints),  & 
                                                              a_integrand(ir:npoints), &
                                                              npoints-ir+1)
            endif 
            phist_radarray(ir) = -four*(integral_above + integral_below)/(2*fs + one)


            if(phi_rad(ir).eq.zero)phist_radarray(ir) = zero

            if(phist_radarray(ir).ne.phist_radarray(ir))then 
                write(*,*)'ir = ', ir
                write(*,*)'Nan found'
                stop 
            endif
        enddo 


    end subroutine phist_from_rhost




    subroutine Gradphist_from_rhost(s, phi_rad, rhost_radarray, Gradphist_radarray, npoints)
        ! Radial gradient of DT98 Eqn. D.54 - computes phi st from rho st
        use integrate, only: integrate_r_traps

        implicit none 
        integer :: s, npoints, ir 
        real(kind=CUSTOM_REAL) :: phi_rad(npoints), fs, rr

        complex(kind=CUSTOM_REAL) ::  b_integrand(npoints), a_integrand(npoints),&
                                        rhost_radarray(npoints), & 
                                        Gradphist_radarray(npoints)
        complex(kind=SPLINE_REAL) :: integral_below, integral_above


        fs = real(s, kind=CUSTOM_REAL)
                            
        b_integrand = phi_rad**(fs+two)*rhost_radarray
        a_integrand = phi_rad**(-fs+one)*rhost_radarray
                          

        do ir = 1, npoints
                rr = phi_rad(ir)

            if(ir.eq.1)then
                integral_below = SPLINE_iZERO
            else  
                integral_below = (fs+one)*(rr**(-fs-two)) * integrate_r_traps(phi_rad(1:ir) , & 
                                                                b_integrand(1:ir), &
                                                                ir)
            endif 

            if(ir.eq.npoints)then 
                integral_above = SPLINE_iZERO
            else 
                integral_above = fs*(rr**(fs-one)) * integrate_r_traps(phi_rad(ir:npoints),  & 
                                                                      a_integrand(ir:npoints), &
                                                                      npoints-ir+1)
            endif 
            Gradphist_radarray(ir) = -four*(integral_above - integral_below)/(2*fs + one)

            if(phi_rad(ir).eq.zero)Gradphist_radarray(ir) = zero

            if(Gradphist_radarray(ir).ne.Gradphist_radarray(ir))then 
                write(*,*)'ir = ', ir
                write(*,*)'Nan found'
                stop 
            endif
        enddo 

    end subroutine Gradphist_from_rhost

end module 