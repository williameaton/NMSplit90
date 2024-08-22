program read_specfem_mesh
    use params, only: NL, u_spl, v_spl, udot_spl, vdot_spl, n_unique_rad, unique_r
    use mesh_utils, only: read_proc_coordinates, load_ibool
    use spline, only: interpolate_mode_eigenfunctions
    
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
    character(len=30) :: eigstring
    integer :: n, l, m, i
    
    ! Setup parameters: 
    region = 3      ! Inner core
    nprocs = 6
    iproc  = 0
    mode_type = 'S' 
    n = 23
    l = 12
    m = 7

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

    ! Spline interpolation: 
    !   - Computes the value of u, v, du, dv at each of the unique radial values
    call interpolate_mode_eigenfunctions(mode_type, u, v, du, dv)
    
    ! Output the spline values: 
    ! Save eigenfunctions to text file in column format 

    write(eigstring,'(a,i2,a, i2, a)') 'spline_', n, mode_type, l, '.txt'
    open(1,file=trim(eigstring))
    do i =1, n_unique_rad
        if(mode_type=='S')then 
            write(1,*)unique_r(i), u_spl(i), udot_spl(i), v_spl(i), vdot_spl(i)
        elseif (mode_type=='T')then 
            write(1,*)unique_r(i), u_spl(i), udot_spl(i)
        endif 
    enddo 
    close(1)


    end program
    
    
    
    