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


end module PREMModel