include "precision.h"               
                                    
double precision, parameter :: PI = 3.141592653589793d0
double precision, parameter :: TWO_PI = 2.d0 * PI
double precision, parameter :: PI_OVER_FOUR = PI / 4.d0
double precision, parameter :: PI_OVER_TWO = PI / 2.0d0
double precision, parameter :: PI_TOL = 0.0000002d0 


! Numbers: 
real(kind=CUSTOM_REAL), parameter :: ZERO  = 0.0000000d0
real(kind=CUSTOM_REAL), parameter :: HALF  = 0.5000000d0
real(kind=CUSTOM_REAL), parameter :: ONE   = 1.0000000d0
real(kind=CUSTOM_REAL), parameter :: TWO   = 2.0000000d0
real(kind=CUSTOM_REAL), parameter :: THREE = 3.0000000d0
real(kind=CUSTOM_REAL), parameter :: FOUR  = 4.0000000d0

! Imaginary i: 
complex(kind=CUSTOM_REAL), parameter :: iZERO  = (ZERO, ZERO)
complex(kind=CUSTOM_REAL), parameter :: iONE   = (ZERO, ONE)

! Scaling factors
real(8), parameter :: RHOAV     = 5514.3d0
real(8), parameter :: GRAV      = 6.6723d-11
real(8), parameter :: SCALE_R   = 6371000.0d0
real(8), parameter :: ACCENORM  = PI * GRAV * RHOAV * SCALE_R
real(8), parameter :: MOMENORM  = PI * GRAV * (RHOAV ** 2) * (SCALE_R ** 5)
