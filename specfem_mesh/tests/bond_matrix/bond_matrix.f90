program bond_matrix 

    implicit none 
    include "constants.h"
    integer :: i, j, iters
    real(kind=CUSTOM_REAL) :: n1, n2, M(6,6), Q(6,6), Cnat(6,6), Cxyz(6,6), A, C, L, N, F

    ! Timing: 
    integer :: start_clock, end_clock, count_rate
    real(8) :: elapsed_time

    n1 = 0.0
    n2 = 0.0
    iters = 1000000

    ! We want to time the bond matrix to see which is faster to compute 
    !! Start clock count
    !call system_clock(count_rate=count_rate)
    !call system_clock(start_clock)

    !do i = 1, iters
       call compute_bond_matrix(PI/iters, n2, M)
    !enddo 

    call print_bond_matrix(M)

    !! Compute run time
    !call system_clock(end_clock)
    !! Calculate the elapsed time in seconds
    !elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
    !! Print the elapsed time
    !write(*,*) 'Wall clock time taken with matrix:', elapsed_time, 'seconds'

    !! Start clock count
    !call system_clock(count_rate=count_rate)
    !call system_clock(start_clock)

    !do i = 1, iters
    !    call compute_bond_matrix_explicit(PI/iters, n2, Q)
    !enddo 
    !! Compute run time
    !call system_clock(end_clock)
    !! Calculate the elapsed time in seconds
    !elapsed_time = real(end_clock - start_clock, kind=8) / real(count_rate, kind=8)
    !! Print the elapsed time
    !write(*,*) 'Wall clock time taken for explicit:', elapsed_time, 'seconds'


end program bond_matrix
