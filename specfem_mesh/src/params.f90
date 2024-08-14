module params 
    
include "precision.h"


! Numbers: 
real(KIND=CUSTOM_REAL), parameter :: HALF  = 0.5000000d0
real(KIND=CUSTOM_REAL), parameter :: ONE   = 1.0000000d0
real(KIND=CUSTOM_REAL), parameter :: TWO   = 2.0000000d0
real(KIND=CUSTOM_REAL), parameter :: FOUR  = 2.0000000d0





! Data: 
character(len=250) :: datadir = 'DATABASES_MPI/inner_core'

! Mesh coordinates 
double precision, allocatable :: xstore(:,:,:,:), & 
                                 ystore(:,:,:,:), & 
                                 zstore(:,:,:,:)
! GLL values
integer :: nspec
integer :: nglob
integer :: ngllx
integer :: nglly
integer :: ngllz

! Local mesh variables:
integer, allocatable                :: ibool(:,:,:,:)
real(kind=CUSTOM_REAL), allocatable :: rho(:,:,:,:)

! Global mesh variables 
real(kind=CUSTOM_REAL), allocatable :: globalrho(:)


end module params