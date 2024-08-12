!       NR: number of radial knots in eigenfunctions and number of
!       layers in PREM
!       NOC: layer index of ocean floor in PREM (- side)
!       NMOHO: layer index of moho in PREM (+ side)
!       NLVZ: layer index of top of low velocity zone in PREM (- side)
!       NTZ: layer index of top of transition zone in PREM (- side)
!       NCMB: layer index of CMB in PREM (+ side)
!       NICB: layer index of ICB in PREM (- side)

! ---  use new prem layer model and catalogues in the datalib
integer, parameter :: NL      =     185 
integer, parameter :: NR      =     185
integer, parameter :: NOC     =     999 ! no ocean in this model 
integer, parameter :: NMOHO   =     169 
integer, parameter :: NLVZ    =     999
integer, parameter :: NTZ     =     999
integer, parameter :: NCMB    =     67
integer, parameter :: NICB    =     33

!       PI: pi
!       PI2: 2 * pi
real, parameter :: PI   = 3.14159265358979 
real, parameter :: PI2  = 6.28318530717958


!
!       RA: radius of Earth (m)
!       GRAV: gravitational constant (N m^2 / kg^2)
!       RHOAV: average density of the Earth (kg / m^3)
!       SCALE_FAC: scaling factor for seismograms; one needs to divide
!                  data and synthetics by SCALE_FAC to get true
!                  ground accelleration
!
real, parameter :: RA         = 6371000.0
real, parameter :: GRAV       = 6.6723e-11
real, parameter :: RHOAV      = 5515.0
real, parameter :: SCALE_FAC  = 1.0E+10


!       DMDL: perturbation in model parameters used in sb
!       DMDL_FD: perturbation in model parameters used in sb_fd
real, parameter :: DMDL       = 1.0E-03
real, parameter :: DMDL_FD    = 0.5E-02


!       SCALE_KAPPA: scale factor to convert perturbations in shear velocity
!                    to perturbations in bulk modulus
!       SCALE_MU: scale factor to convert perturbations in shear velocity
!                 to perturbations in shear modulus
!       SCALE_ALPHA: scale factor to convert perturbations in shear velocity
!                    to perturbations in compressional velocity
!       SCALE_RHO: scale factor to convert perturbations in shear velocity
!                  to perturbations in density
!
real, parameter :: SCALE_KAPPA = 0.5
real, parameter :: SCALE_MU    = 2.4
real, parameter :: SCALE_ALPHA = 0.55
real, parameter :: SCALE_RHO   = 0.4

!       NK: number of radial eigenfunctions for mantle model
!       NS: number of angular degrees for mantle model
!       ND: number of discontinuities for mantle model
integer, parameter :: NK = 13
integer, parameter :: NS = 12
integer, parameter :: ND = 1

!       LMAX: largest angular degree of modes included in inversion
integer, parameter :: LMAX = 24

!       NMMAX: maximum length of the model vector
integer, parameter :: NMMAX = 10000

!       NSMAX: maximum number of spectra per component
!       NWIN:  maximum number of time windows

integer, parameter :: NSMAX   = 100
integer, parameter :: NWIN    = 1


!       EMAX: maximum number of events
!       CMAX: maximum number of clusters per component per event
!       MMAX: maximum number of modes per cluster
!       GMAX: maximum number of glitches per record

integer, parameter :: EMAX = 10
integer, parameter :: CMAX = 150
integer, parameter :: MMAX = 10
integer, parameter :: GMAX = 3


!       NTFMAX: maximum number samples in the fast time series
!       NTSMAX: maximum number samples in the slow time series
!       NWMAX:  maximum number samples in spectrum

integer, parameter :: NTFMAX  = 1000000
integer, parameter :: NTSMAX  = NTFMAX
integer, parameter :: NWMAX   = 5000

!     SSMAX : maximum number of point sources
integer, parameter :: SSMAX=1500
