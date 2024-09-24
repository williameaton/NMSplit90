
program test_W_matrix
    use params, only: NL, IC_ID, ndisc, rdisc, disc, disp1, disp2, & 
                     rad_mineos, rho_mineos, u_spl, v_spl, interp_map, & 
                     interp_id_r, unique_r, n_unique_rad, ngllx, nglly, &
                     ngllz, nspec, Wmat, rad_id, wgll, detjac, rho_spl, & 
                      nglob, realfmt
    use spline, only: create_interpolation_radial_map, interpolate_mode_eigenfunctions, & 
                      write_mode_spline, interpolate_mineos_variable, quad_spline_interp_3, & 
                      write_scalar_spline
    use Integrate, only: integrate_r_traps
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: compute_jacobian, compute_rtp_from_xyz, load_ibool, & 
                          read_proc_coordinates, rotate_complex_vector_rtp_to_xyz, & 
                          compute_rotation_matrix, setup_global_coordinate_arrays, & 
                          cleanup_for_mode, map_local_global, map_complex_vector
    use gll, only: setup_gll
    use math, only: sqrtp
    implicit none
    include "constants.h"

    integer :: i, j, k, l1, n1, l2, n2 , ispec, lentrim, iproc, region, m1, m2, knot_lower, knot_upper, component
    character(len=1)   :: char, type_1, type_2
    character(len=10)  :: mode
    character(len=250) :: arg, ensight_nm
    character(len=30)  :: out_name
    character(len=4), parameter :: fmti1 = "(i1)"
    character(len=4), parameter :: fmti2 = "(i2)"
    character(len=4), parameter :: fmti3 = "(i3)"
    character(len=4) :: fmt_n, fmt_l
    real(SPLINE_REAL) :: wcom, qmod
    integer :: npoints, tl1, tl2, nproc
    real(kind=CUSTOM_REAL), allocatable :: radial_vals(:), r_lower, r_upper

    real(SPLINE_REAL), allocatable :: integrand(:)
    complex(SPLINE_REAL), allocatable :: W_s(:)
    complex(SPLINE_REAL), allocatable :: disp1_glob(:,:)
    real(SPLINE_REAL) :: total_integral,  precision,  lf, kf

    ! Timing: 
    integer :: start_clock, end_clock, count_rate
    real(8) :: elapsed_time


    ! Read mineos model 
    call process_mineos_model(.false.)


    ! Now we can load the mode: 
    ! Choose a mode: 
    type_1 = 'S'
    l1      = 2
    n1      = 0

    type_2 = 'S'
    l2      = 2
    n2      = 0


    region = 3
    nproc  = 6

    ! Setup W matrix
    tl1 = 2*l1 + 1
    tl2 = 2*l2 + 1

    ! The matrix should be 2l + 1 from -m to m 
    call deallocate_if_allocated(Wmat)
    call allocate_if_unallocated(tl1, tl2, Wmat)
    Wmat = SPLINE_iZERO


    ! Start clock count
    call system_clock(count_rate=count_rate)
    call system_clock(start_clock)


    do iproc = 0, nproc-1
        ! Things that need to be done for each processor
        call read_proc_coordinates(iproc, region)

        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian(iproc,.false.)

        call setup_global_coordinate_arrays(iproc, .false.)
        call compute_rtp_from_xyz(iproc, .false.)
        call get_mesh_radii(iproc, .false.)
        call compute_rotation_matrix()

        call compute_W_matrix(type_1, l1, n1, type_2, l2, n2, .true., iproc)

        call cleanup_for_mode()
    enddo 


    ! Constants so multiply after
    Wmat = Wmat * OMEGA * iONE

    write(out_name, '(a,i1,a,i1,a,i1,a,i1,a)')'rot_mat/Wmat_', n1, type_1, l1, '_', n2, type_2, l2, '.txt'
    call save_W_matrix(l1, l2, trim(out_name))




    ! Compute run time
    call system_clock(end_clock)
    ! Calculate the elapsed time in seconds
    elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
    ! Print the elapsed time
    write(*,*) 'Wall clock time taken for full SEM:', elapsed_time, 'seconds'




end program test_W_matrix