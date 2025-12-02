module PREMModel
    implicit none 
    include "constants.h"

    contains 
    


    subroutine get_PREM_ACLNF_at_radius(radius, A, C, L, N, F)
        ! Gets NON-DIMENSIONALISED ACLNF for PREM model at radius
        ! Currently only for IC
        ! Note that IC is isotropic in prem 
        real(kind=CUSTOM_REAL) :: radius, A, C, L, N, F
        real(kind=CUSTOM_REAL) :: RIC = 1221.5d0/6371d0



        ! Local 
        real(kind=CUSTOM_REAL) :: eta ! eta = F / (A - 2L) and is 1 for isotropic solid
        real(kind=CUSTOM_REAL) :: rho, vp, vs


        if(radius.gt.RIC+ 1.0e-5)then 
            write(*,'(a, f12.6, a, f12.6)')'Error in get_PREM_ACLNF_at_radius: Only setup for IC but using radius ', radius, ' when maximum is ', RIC
            stop
        else 

            eta = ONE 
            
            ! Originally g/cm^3 --> kg/m^3 -> nondimensionalised
            rho = (13.0885d0 - 8.8381d0*radius*radius)*1000.d0/RHOAV
            ! Originally km/s --> m/s 
            vp  = (11.2622d0 - 6.3640d0*radius*radius)*1000.d0/SCALE_V
            vs  = (3.6678d0 - 4.4475d0*radius*radius)*1000.d0/SCALE_V
  
            A = rho * vp * vp 
            C = A
            L = rho * vs * vs
            N = L
            F = eta*(A - TWO*L)

        endif 
    end subroutine get_PREM_ACLNF_at_radius


    subroutine get_radial_density_derivative_at_radius(x, drhodr, k)
        implicit none 
        real(kind=CUSTOM_REAL) :: x, drhodr
        integer :: k ! The z gll value
        
        if(x.lt.1221500.0d0/SCALE_R)then 
            ! Inner core: 
            drhodr =          -two*8.8381d0*x 
        elseif(x.gt.1221500.0d0/SCALE_R .and. x.lt.3480000.0d0/SCALE_R)then 
            ! Outer core
            drhodr = -1.2638  -two*3.6426d0*x  - three*5.5281d0*x*x
        elseif(x.gt.3480000.0d0/SCALE_R .and. x.lt.5701000.0d0/SCALE_R)then 
            ! Lower mantle
            drhodr = -6.4761  +two*5.5283d0*x  - three*3.0807d0*x*x
        elseif(x.gt.5701000.0d0/SCALE_R .and. x.lt.5771000.0d0/SCALE_R)then 
            ! Transition zone 1
            drhodr = -1.4836d0 
        elseif(x.gt.5771000.0d0/SCALE_R .and. x.lt.5971000.0d0/SCALE_R)then 
            ! Transition zone 2
            drhodr = -8.0298d0
        elseif(x.gt.5971000.0d0/SCALE_R.and.x.lt.6151000.0d0/SCALE_R)then 
            ! Transition zone 3
            drhodr = -3.8045
        elseif(x.gt.6151000.0d0/SCALE_R .and. x.lt.6346600.0d0/SCALE_R)then 
            ! LVZ + LID
            drhodr = 0.6924d0
        elseif(x.gt.6346600.0d0/SCALE_R)then 
            ! Crust
            drhodr = 0.0d0
        !! -------------------------------------
        !! BOUNDARY CASES
        elseif(x.eq.3480000.0d0/SCALE_R)then 
            if(k.eq.1)then 
                ! Above the boundary
                drhodr = -6.4761  +two*5.5283d0*x  - three*3.0807d0*x*x
            elseif(k.eq.5)then 
                ! Below the boundary
                drhodr = -1.2638  -two*3.6426d0*x  - three*5.5281d0*x*x
            else 
                write(*,*)'Error - on an boundary gll but on PREM boundary.'
            endif

        elseif(x.eq.5701000.0d0/SCALE_R)then 
            if(k.eq.1)then 
                ! Above the boundary
                drhodr = -8.0298d0
            elseif(k.eq.5)then 
                ! Below the boundary
                drhodr = -6.4761  +two*5.5283d0*x  - three*3.0807d0*x*x
            else 
                write(*,*)'Error - on an boundary gll but on PREM boundary.'
            endif
        elseif(x.eq.5771000.0d0/SCALE_R)then 
            if(k.eq.1)then 
                ! Above the boundary
                drhodr = -8.0298d0
            elseif(k.eq.5)then 
                ! Below the boundary
                drhodr = -1.4836d0 
            else 
                write(*,*)'Error - on an boundary gll but on PREM boundary.'
            endif

        elseif(x.eq.5971000.0d0/SCALE_R)then 
            if(k.eq.1)then 
                ! Above the boundary
                drhodr = -3.8045
            elseif(k.eq.5)then 
                ! Below the boundary
                drhodr = -8.0298d0
            else 
                write(*,*)'Error - on an boundary gll but on PREM boundary.'
            endif

        elseif(x.eq.6151000.0d0/SCALE_R)then 
            if(k.eq.1)then 
                ! Above the boundary
                drhodr = 0.6924d0
            elseif(k.eq.5)then 
                ! Below the boundary
                drhodr = -3.8045
            else 
                write(*,*)'Error - on an boundary gll but on PREM boundary.'
            endif

        elseif(x.eq.6346600.0d0/SCALE_R)then 
            if(k.eq.1)then 
                ! Above the boundary
                drhodr = 0.0d0
            elseif(k.eq.5)then 
                ! Below the boundary
                drhodr = 0.6924d0
            else 
                write(*,*)'Error - on an boundary gll but on PREM boundary.'
            endif

        else 
            write(*,*)'PREM radius not covered. Stop.', x, x*SCALE_R
            stop
        endif 

        ! Note the extra SCALE_R accounts for the derivative
        drhodr = drhodr *1000.d0/(RHOAV*SCALE_R)
    end subroutine

end module PREMModel