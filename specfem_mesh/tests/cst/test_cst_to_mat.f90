program test_cst_to_mat
    use splitting_function, only: get_Ssum_bounds, cst_to_H, H_to_cst
    use V_ani, only: save_Vani_matrix, save_Vani_real_matrix
    implicit none 

    include "constants.h"

    character :: t1, t2
    integer :: n1, n2, l1, l2, tl1, tl2, smin, smax, num_s, is, it, ncols, s
    real(kind=SPLINE_REAL), allocatable :: cst(:,:), Hmat(:,:)

    ! Example 
    t1 = 'T'
    t2 = 'S'
    

    l1  = 8 ! ld
    n1  = 0 
    tl1 = 2*l1 + 1

    ! Self coupling
    l2  = 3 ! l
    n2  = n1
    tl2 = tl1
    
    allocate(Hmat(tl1, tl2))

    call get_Ssum_bounds(l1, l2, smin, smax, num_s, ncols)

    ! Store the Csts in a 2D array. Each row is an s value 
    ! first row is smin and last row is smax
    ! columns go from t = -s to +s so that the max number of columns
    ! is 2*smax + 1
    allocate(cst(num_s, ncols))
    cst = ZERO

    write(*,*)'Min s: ', smin
    write(*,*)'Max s: ', smax
    write(*,*)'Num s: ', num_s
    write(*,*)'Num c: ', ncols


    ! Lets fill up the Cst array: 
    do is = 1, num_s 
        s = smin + is - 1
        write(*,*)'s is', s
        do it = 1, 2*s + 1
            cst(is,it) = real(-s + it -1, kind=CUSTOM_REAL)/TEN
        enddo 
    enddo 

    write(*,*)'------ Original cst: -------'
     !Print Cst matrix:
    do is = 1, num_s 
        do it = 1, ncols
            write(*,*) real(cst(is,it))
        enddo 
        write(*,*)
    enddo 
    write(*,*)'----------------------------'

    ! compute M from cst
    call cst_to_H(Hmat, l1, l2, cst, ncols, num_s, t1, t2)


    call save_Vani_real_matrix(l1, Hmat, './cst/Hmat.txt')


    call H_to_cst(Hmat, l1, l2, cst, ncols, num_s, t1, t2, 1)
    
    !Print Cst matrix:
    do is = 1, num_s 
        do it = 1, ncols
            write(*,*) real(cst(is,it))
        enddo 
        write(*,*)
    enddo 

end program 
