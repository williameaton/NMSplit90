program read_specfem_mesh
    use params, only: NL, u_spl, v_spl, udot_spl, vdot_spl, n_unique_rad, unique_r, interp_id_r, IC_ID
    use mesh_utils, only: read_proc_coordinates, load_ibool
    use spline, only: interpolate_mode_eigenfunctions, write_mode_spline, create_interpolation_radial_map
    use allocation_module, only: allocate_if_unallocated
    implicit none
    include "constants.h"

    integer :: iproc             
    integer :: nprocs                ! Number of processors used by mesher 
    integer :: region                ! Region code

    ! Variables for splines/modes: 
    real(kind=SPLINE_REAL)  :: omega  ! normalized anglar frequency   
    real(kind=SPLINE_REAL)  :: qval   ! normalized Q factor        
    real(kind=SPLINE_REAL), allocatable :: u(:), du(:)
    real(kind=SPLINE_REAL), allocatable :: v(:), dv(:)
    character(len=1) :: mode_type
    integer :: n, l, m, i
    
    ! Setup parameters: 
    region = 3      ! Inner core
    nprocs = 6
    iproc  = 0
    mode_type = 'S' 
    n = 23
    l = 12

    ! Read mineos model 
    call process_mineos_model()

    ! Read the mesh info and coordinates
    call read_proc_coordinates(iproc, region)

    ! Load ibool variable: 
    call load_ibool(iproc, region)

    ! Get unique mesh radii that are present
    call get_mesh_radii()

    allocate(u(NL), v(NL), du(NL), dv(NL))


    ! Load the mode from the database and save it (.true.)
    call get_mode(mode_type, n, l, omega, qval, u, du, v, dv, .true.)

    call interpolate_mode_eigenfunctions(mode_type, u, v, du, dv, 1, IC_ID, &  
                                         unique_r, n_unique_rad, interp_id_r)



    call write_mode_spline(n, mode_type, l, unique_r, n_unique_rad)

    end program
    
    
    
    