module params 
    
include "precision.h"

! Increasing number = increasing verbosity
! 0  == no printing
! 1  == minimal updates for fast run
! 2  == while coding
! 3  == while debugging
integer, parameter :: verbose      = 0
logical, parameter :: all_warnings = .false.

! Specfem mesh files: 
character(len=250) :: datadir = '/scratch/gpfs/we3822/NMSplit90/specfem_mesh/DATABASES_MPI/NEX176/sliced/'

! Mineos model parameters: 
character(len=250), parameter  :: ddir = '/scratch/gpfs/we3822/NMSplit90/databases/prem_ani_att_database/'
character(len=60),  parameter  :: model_fname = 'model'
integer                        :: NR
integer                        :: NL
integer                        :: IC_ID     ! IC side of ICB
integer                        :: CMB_ID    ! OC side of CMB
real(CUSTOM_REAL)              :: RA  

real(kind=CUSTOM_REAL), allocatable :: rad_mineos(:)  ! Non-dimensionalised
real(kind=CUSTOM_REAL), allocatable :: radius(:)      ! Dimensionalised
real(kind=CUSTOM_REAL), allocatable :: rho_mineos(:)  ! Density (non-dimensionalised by RHOAV)
real(kind=CUSTOM_REAL), allocatable :: vp(:)         ! Non Dimensionalised


! Radial discontinuities
!   ndisc: number of discontinuities
!   rdisc: radius of discontinuities
!   disc: knot ID of discontinuity (-ve side)
integer :: ndisc
integer, allocatable                :: disc(:)  
real(kind=CUSTOM_REAL), allocatable :: rdisc(:)  


! Mesh coordinates 
real(kind=CUSTOM_REAL), allocatable, target ::  xstore(:,:,:,:),     & 
                                                ystore(:,:,:,:),     & 
                                                zstore(:,:,:,:),     &
                                                rstore(:,:,:,:),     & 
                                                thetastore(:,:,:,:), & 
                                                phistore(:,:,:,:)
                                        

real(kind=CUSTOM_REAL), allocatable :: unique_r(:)
integer :: n_unique_rad
integer, allocatable :: rad_id(:, :, :, :)  ! Points to the unique radius
integer, allocatable :: interp_id_r(:)
integer, allocatable :: interp_map(:)


! GLL values
integer :: nspec
integer :: nglob
integer :: ngllx
integer :: nglly
integer :: ngllz
real(kind=CUSTOM_REAL), allocatable :: xi(:), wgll(:), dgll(:,:)
real(kind=SPLINE_REAL), allocatable :: wglljac(:,:,:,:)



! Local mesh variables:
integer, allocatable                   :: ibool(:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: rho(:,:,:,:)
complex(kind=SPLINE_REAL), allocatable :: strain1(:,:,:,:,:), strain2(:,:,:,:,:)
complex(kind=SPLINE_REAL), allocatable :: strains(:,:,:,:,:,:)
complex(kind=SPLINE_REAL), allocatable :: disp1(:,:,:,:,:), disp2(:,:,:,:,:)

real(kind=CUSTOM_REAL), allocatable    :: jac(:,:,:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: detjac(:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: jacinv(:,:,:,:,:,:)

! Global mesh variables 
real(kind=CUSTOM_REAL), allocatable    :: globalrho(:)
real(kind=CUSTOM_REAL), allocatable    :: x_glob(:), y_glob(:), z_glob(:)
real(kind=CUSTOM_REAL), allocatable    :: theta_glob(:), phi_glob(:)
complex(kind=SPLINE_REAL), allocatable :: globalstrain(:,:)
complex(kind=CUSTOM_REAL), allocatable :: globaldisp(:,:)
real(kind=CUSTOM_REAL), allocatable :: Rmat(:,:,:)


! Matrices: 
complex(kind=SPLINE_REAL), allocatable :: Wmat(:,:)
complex(kind=SPLINE_REAL), allocatable :: Vani(:,:)

! Angles of rotation for VTI 
real(kind=CUSTOM_REAL), allocatable    :: eta1(:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: eta2(:,:,:,:)

! Perturbed elastic tensor in xyz at each GLL (6 x 6) voigt notation
real(kind=SPLINE_REAL), allocatable    :: Cxyz(:,:,:,:,:,:)
real(kind=CUSTOM_REAL), allocatable    :: Arad(:), Crad(:), Lrad(:), & 
                                          Nrad(:), Frad(:)

! Spline arrays: 
real(kind=SPLINE_REAL), allocatable :: u_spl(:),    &
                                       udot_spl(:), &
                                       v_spl(:),    &
                                       vdot_spl(:), &
                                       xx(:), zz(:),&
                                       rho_spl(:),  &  
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

end module params
