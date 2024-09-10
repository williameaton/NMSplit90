
program test_Vani_matrix
    use params, only: Vani
    use spline, only: 
    use Integrate, only: 
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: cleanup_for_mode, compute_jacobian, compute_rotation_matrix, & 
                          compute_rtp_from_xyz, load_ibool, read_proc_coordinates, & 
                          setup_global_coordinate_arrays
    use gll, only: setup_gll
    use v_ani, only: compute_Vani_matrix, save_Vani_matrix, compute_Cxyz_at_gll

    implicit none
    include "constants.h"

    real(kind=CUSTOM_REAL) :: A, C, L, N, F
    integer :: iproc, i,j,k,ispec, l1, l2, n1, n2, nproc, region, tl1, tl2, h, b
    character ::  type_1, type_2
    character(len=250) :: out_name
    real(kind=SPLINE_REAL) :: min_r, min_i, thirty, twone

    ! Read mineos model 
    call process_mineos_model()

    ! Choose modes: 
    n1      = 6
    type_1 = 'S'
    l1      = 3

    n2      = 6
    type_2 = 'S'
    l2      = 3

    A =  0.4d0
    C = -0.2d0
    L =  0.3d0
    N = -0.5d0
    F =  0.1d0

    region = 3
    nproc  = 6


    ! Setup Vani matrix
    tl1 = 2*l1 + 1
    tl2 = 2*l2 + 1

    ! The matrix should be 2l + 1 from -m to m 
    call allocate_if_unallocated(tl1, tl2, Vani)
    Vani = SPLINE_iZERO

    do iproc = 0, nproc-1
        write(*,*)"Processor: ", iproc
        ! Things that need to be done for each processor
        call read_proc_coordinates(iproc, region)

        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian()

        call setup_global_coordinate_arrays()
        call compute_rtp_from_xyz()
        call get_mesh_radii()
        call compute_rotation_matrix()

        call compute_Cxyz_at_gll(A, C, L, N, F)

        call compute_Vani_matrix(type_1, l1, n1, type_2, l2, n2, .false., iproc)

        call cleanup_for_mode()
    enddo 


    write(out_name, '(a,i1,a,i1,a)')'sem_', n1, type_1, l1, '.txt'
    call save_Vani_matrix(l1, out_name)


end program test_Vani_matrix