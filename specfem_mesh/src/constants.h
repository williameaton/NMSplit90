include "precision.h"               
                                    
real(kind=CUSTOM_REAL), parameter :: PI = 3.141592653589793_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: TWO_PI = 2._CUSTOM_REAL * PI
real(kind=CUSTOM_REAL), parameter :: PI_OVER_FOUR = PI / 4._CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: PI_OVER_TWO = PI / 2.0_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: PI_TOL = 0.0000002_CUSTOM_REAL

! Tolerances (not really constants!)
double precision, parameter :: pole_tolerance = 1e-4 ! ~600 metres from pole

! Numbers: 
real(kind=CUSTOM_REAL), parameter :: ZERO  = 0.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: HALF  = 0.5000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: ONE   = 1.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: TWO   = 2.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: THREE = 3.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: FOUR  = 4.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: FIVE  = 5.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: SIX   = 6.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: SEVEN = 7.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: EIGHT = 8.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: NINE  = 9.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: TEN   = 10.0000000_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: RT2   = (TWO)**HALF

! Numbers with spline (eigenfunction) precision:
real(kind=SPLINE_REAL), parameter :: SPLINE_PI = 3.141592653589793_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_ZERO  = 0.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_HALF  = 0.5000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_ONE   = 1.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_TWO   = 2.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_THREE = 3.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_FOUR  = 4.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_FIVE  = 5.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_SIX   = 6.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_SEVEN = 7.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_EIGHT = 8.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_NINE  = 9.0000000_SPLINE_REAL
real(kind=SPLINE_REAL), parameter :: SPLINE_TEN   = 10.0000000_SPLINE_REAL
real(kind=CUSTOM_REAL), parameter :: SPLINE_RT2   = (SPLINE_TWO)**SPLINE_HALF

! Imaginary i: 
complex(kind=CUSTOM_REAL), parameter :: iZERO  = (ZERO, ZERO)
complex(kind=CUSTOM_REAL), parameter :: iONE   = (ZERO, ONE)

complex(kind=SPLINE_REAL), parameter :: SPLINE_iZERO  = (ZERO, ZERO)
complex(kind=SPLINE_REAL), parameter :: SPLINE_iONE   = (ZERO, ONE)


! Scaling factors

real(kind=CUSTOM_REAL), parameter :: RHOAV     = 5514.3_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: GRAV      = 6.6723e-11
real(kind=CUSTOM_REAL), parameter :: SCALE_R   = 6371000.0_CUSTOM_REAL
real(kind=CUSTOM_REAL), parameter :: SCALE_T  =  (PI * GRAV * RHOAV)**(-HALF)
real(kind=CUSTOM_REAL), parameter :: ACCENORM  = PI * GRAV * RHOAV * SCALE_R
real(kind=CUSTOM_REAL), parameter :: MOMENORM  = PI * GRAV * (RHOAV ** 2) * (SCALE_R ** 5)

! Rate of rotation
real(kind=CUSTOM_REAL), parameter :: OMEGA    = (TWO_PI/86400.0_CUSTOM_REAL) * SCALE_T
