  real*8, parameter :: PI = 3.14159265358979d0

  integer, parameter :: NR = 185
  integer, parameter :: NL = 185
  
  real*8, parameter :: RHOAV = 5514.3d0
  real*8, parameter :: GRAV = 6.6723d-11
  real*8, parameter :: RA = 6371000.0d0

  real*8, parameter :: ACCENORM = PI * GRAV * RHOAV * RA
  real*8, parameter :: MOMENORM = PI * GRAV * (RHOAV ** 2) * (RA ** 5)


!    for reference, won't compile
!  real*8, parameter :: TIMENORM = 1.0 / dsqrt (PI * GRAV * RHOAV)
!  real*8, parameter :: VELONORM = dsqrt (PI * GRAV * RHOAV) * RA
!  real*8, parameter :: ANGUNORM = dsqrt (PI * GRAV * RHOAV)



  
