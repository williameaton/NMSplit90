module params 
    
! Data: 
character(len=250) :: datadir = 'DATABASES_MPI/inner_core'

! Precision - see specfem3d_globe/shared/constants.h (not .h.in): 
integer, parameter :: CUSTOM_REAL = 4

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