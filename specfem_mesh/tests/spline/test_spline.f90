program test_spline
    use params, only: nprocs
    use allocation_module, only: allocate_if_unallocated
    use mineos_model, only: mineos, mineos_ptr
    use specfem_mesh, only: SetMesh, create_SetMesh
    use modes, only: get_mode, Mode
    use piecewise_interpolation, only: InterpPiecewise, create_PieceInterp
    implicit none
    include "constants.h"
    integer :: iproc          
    type(SetMesh) :: sm   
    type(Mode)    :: M 
    type(InterpPiecewise) :: IP 

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos 

    ! Read the mesh info and coordinates for set 0 region 3
    sm = create_SetMesh(0, 3)
    call sm%read_proc_coordinates()
    call sm%load_ibool()
    call sm%get_unique_radii(.false.)

    ! Load the mode from the database and save it (.true.)
    M = get_mode(23, 'S', 12, mineos_ptr, './spline/')


    IP = create_PieceInterp(sm%n_unique_rad)
    IP%radial = sm%unique_r
    call IP%setup()
    call IP%create_interpolation_radial_map()
    call IP%interpolate_mode_eigenfunctions(M)

    call M%write_spline(SM%unique_r)

    end program
    
    
    
    