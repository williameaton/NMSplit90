
program test_fast_Vani_matrix
    use params, only: Vani
    use spline, only: 
    use Integrate, only: 
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: cleanup_for_mode, compute_jacobian, compute_rotation_matrix, & 
                          compute_rtp_from_xyz, load_ibool, read_proc_coordinates, & 
                          setup_global_coordinate_arrays
    use gll, only: setup_gll
    use v_ani, only: compute_Vani_matrix_stored_selfcoupling, save_Vani_matrix, compute_Cxyz_at_gll_constantACLNF
    
    implicit none
    include "constants.h"

    real(kind=CUSTOM_REAL) :: A, C, L, N, F
    integer :: iproc, i,j,k,ispec, l1, l2, n1, m1,m2, n2, nproc, region, tl1, tl2, h, b
    character ::  type_1, type_2
    character(len=250) :: out_name
    real(kind=SPLINE_REAL) :: min_r, min_i, thirty, twone

   ! Timing: 
    integer :: start_clock, end_clock, count_rate
    real(8) :: elapsed_time


    ! Start clock count
    call system_clock(count_rate=count_rate)
    call system_clock(start_clock)





    ! Read mineos model 
    call load_mineos_radial_info()

    ! Choose modes: 
    n1      = 6
    type_1 = 'S'
    l1      = 5

    A =  0.4d0
    C = -0.2d0
    L =  0.3d0
    N = -0.5d0
    F =  0.1d0

    region = 3
    nproc  = 6

    ! Setup Vani matrix
    tl1 = 2*l1 + 1
    allocate(Vani(tl1, tl1))
    Vani = SPLINE_iZERO

    do iproc = 0, nproc-1
        write(*,*)"Processor: ", iproc
        ! Things that need to be done for each processor
        call read_proc_coordinates(iproc, region)

        call load_ibool(iproc, region)
        call setup_gll()
        call compute_jacobian(iproc, .false.)

        call setup_global_coordinate_arrays(iproc, .false.)
        call compute_rtp_from_xyz(iproc, .false.)
        call get_mesh_radii(iproc, .false.)
        call compute_rotation_matrix()
        call compute_Cxyz_at_gll_constantACLNF(A, C, L, N, F, zero, zero)

        call compute_Vani_matrix_stored_selfcoupling(type_1, l1, n1, iproc)

        call cleanup_for_mode()
    enddo 


    ! Symmetry D.153
    !do m1 = -l1+1, l1
     !   do m2 = -m1+1, l1 
     !       write(*,*)m1, m2, Vani(-m1+l1+1, -m2+l1+1), (-SPLINE_ONE)**real(m1+m2, kind=SPLINE_REAL)
     !        Vani(m1+l1+1, m2+l1+1) =  (-SPLINE_ONE)**real(m1+m2, kind=SPLINE_REAL) * conjg(Vani(-m1+l1+1, -m2+l1+1))
     !    enddo 
     !enddo 


    write(out_name, '(a,i1,a,i1,a)')'./v_ani_matrix/sem_fast_', n1, type_1, l1, '.txt'
    call save_Vani_matrix(l1, out_name)



    ! Compute run time
    call system_clock(end_clock)
    ! Calculate the elapsed time in seconds
    elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
    ! Print the elapsed time
    write(*,*) 'Wall clock time taken for test_fast_Vani_matrix:', elapsed_time, 'seconds'




end program test_fast_Vani_matrix