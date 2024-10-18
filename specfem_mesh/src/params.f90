module params 
    
include "precision.h"

! Increasing number = increasing verbosity
! 0  == no printing
! 1  == minimal updates for fast run
! 2  == while coding
! 3  == while debugging
integer, parameter :: verbose       = 0
logical, parameter :: all_warnings  = .false.
logical, parameter :: safety_checks = .false.

integer, parameter :: nprocs       = 8
integer, parameter :: nmodes       = 28

! Specfem mesh files: 
character(len=250) :: datadir = '/scratch/gpfs/we3822/NMSplit90/specfem_mesh/DATABASES_MPI/NEX176/sliced/linear/sets8/'

! Mineos model parameters: 
character(len=250), parameter  :: ddir = '/scratch/gpfs/we3822/NMSplit90/databases/prem_ani_att_database/'
character(len=60),  parameter  :: model_fname = 'model'

! Local mesh variables:
real(kind=CUSTOM_REAL),    allocatable :: rho(:,:,:,:)
complex(kind=SPLINE_REAL), allocatable :: strain1(:,:,:,:,:), strain2(:,:,:,:,:)
complex(kind=SPLINE_REAL), allocatable :: strains1(:,:,:,:,:,:), strains2(:,:,:,:,:,:)
complex(kind=SPLINE_REAL), allocatable :: disp1(:,:,:,:,:), disp2(:,:,:,:,:)

! Matrices: 
complex(kind=SPLINE_REAL), allocatable :: Wmat(:,:)
complex(kind=SPLINE_REAL), allocatable :: Vani(:,:)

! Angles of rotation for VTI 
real(kind=CUSTOM_REAL), allocatable    :: eta1(:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: eta2(:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: glob_eta1(:)
real(kind=CUSTOM_REAL), allocatable    :: glob_eta2(:)



! Perturbed elastic tensor in xyz at each GLL (6 x 6) voigt notation
real(kind=SPLINE_REAL), allocatable    :: Cxyz(:,:,:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: Arad(:), Crad(:), Lrad(:), & 
                                          Nrad(:), Frad(:)

! Spline arrays: 
real(kind=SPLINE_REAL), allocatable :: rho_spl(:),  &  
                                       vp_spl(:),   & 
                                       A0(:)



! Visual: 
character(len=250),parameter :: en_dir='./ensight/'      ! Ensight prefix file name
character(len=250) :: en_fname                     ! Ensight prefix file name
integer, parameter :: CASEUNIT         = 11
integer, parameter :: GEOUNIT          = 12
integer, parameter :: TENSORSYMOUT_I   = 14
integer, parameter :: TENSORSYMOUT_R   = 15
integer, parameter :: VECOUT_R         = 16
integer, parameter :: VECOUT_I         = 17
integer, parameter :: REALSCALOUT      = 18
character(len=5), parameter  :: intfmt  = "(i10)"
character(len=7), parameter  :: realfmt = "(e12.5)"



! Voronoi model: 
character(len=250) :: fname_voronoi = '/scratch/gpfs/we3822/NMSplit90/specfem_mesh/DATABASES_MPI/voronoi/voronoi_model.txt'



! MPI 
integer :: myrank
integer :: MPI_CUSTOM_REAL 
integer :: MPI_SPLINE_REAL
integer :: MPI_SPLINE_COMPLEX
integer :: IIN, IOUT
logical :: ONLY_ONE_TASK_PER_SET



! Timing: 
integer :: start_clock, end_clock, count_rate
real(8) :: elapsed_time


end module params

