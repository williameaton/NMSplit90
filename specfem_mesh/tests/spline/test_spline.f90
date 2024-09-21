program test_spline
    use params, only: NL, u_spl, v_spl, udot_spl, vdot_spl, n_unique_rad, unique_r, interp_id_r, IC_ID, interp_map
    use mesh_utils, only: read_proc_coordinates, load_ibool
    use spline, only: interpolate_mode_eigenfunctions, write_mode_spline, create_interpolation_radial_map
    use allocation_module, only: allocate_if_unallocated
    implicit none
    include "constants.h"

    integer :: iproc             
    integer :: nprocs                ! Number of processors used by mesher 
    integer :: region                ! Region code

    ! Variables for splines/modes: 
    real(kind=SPLINE_REAL)  :: om  ! normalized anglar frequency   
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
    call get_mode(mode_type, n, l, om, qval, u, du, v, dv, .true., './spline/')


    call allocate_if_unallocated(n_unique_rad, u_spl)
    call allocate_if_unallocated(n_unique_rad, v_spl)
    call allocate_if_unallocated(n_unique_rad, udot_spl)
    call allocate_if_unallocated(n_unique_rad, vdot_spl)



    allocate(interp_map(n_unique_rad))
    call create_interpolation_radial_map(unique_r, interp_map, n_unique_rad, 1, IC_ID)
    call interpolate_mode_eigenfunctions(mode_type, u, v, du, dv, 1, IC_ID, &  
                                        unique_r, n_unique_rad, interp_id_r, & 
                                        u_spl, v_spl, udot_spl, vdot_spl)
    call write_mode_spline(n, mode_type, l, unique_r, n_unique_rad)

    end program
    
    
    
    