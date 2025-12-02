program test_cst_to_mat
    use params, only: Vani
    use splitting_function, only: get_Ssum_bounds, cst_to_H, Hcomplex_to_cst
    use V_ani, only: save_Vani_matrix, save_Vani_real_matrix
    implicit none 

    include "constants.h"

    character :: t1, t2
    real(kind=8) :: Ast, Bst
    integer  :: tval, irow, icol
    integer :: n1, n2, l1, l2, tl1, tl2, smin, smax, num_s, is, it, ncols, s
    complex(kind=SPLINE_REAL), allocatable :: cst(:,:)

    ! Example 
    t1 = 'S'
    t2 = 'T'
    

    l1  = 5 ! ld
    n1  = 0 
    tl1 = 2*l1 + 1

    l2  = 2
    n2  = 0
    tl2 = 2*l2 + 1
    
    allocate(Vani(tl1, tl2))

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
        do it = 1, 2*s + 1
            tval = it - s - 1 

            Ast = abs(real(tval))/TEN
            Bst = dsin(real(abs(tval), kind=8))

            ! Ignoring the 2pi^1/2 nad 4pi^1/2
            if (tval.eq.0)then
                cst(is,it) = Ast 
            elseif(tval.gt.0)then
                cst(is,it) =  (-1.0d0)**(real(tval)) * (Ast  - SPLINE_iONE*Bst)
            else 
                cst(is,it) =  Ast  + SPLINE_iONE*Bst
            endif 
            
        enddo 
    enddo 

    write(*,*)'------ Original cst: -------'
     !Print Cst matrix:
    do is = 1, num_s 
        s = smin + is - 1

        write(*,*)'s is', s

        do it = 1, 2*s + 1
            tval = it - s - 1 

            write(*,*) cst(is,it)
        enddo 
        write(*,*)
    enddo 
    write(*,*)'----------------------------'

    ! compute M from cst
    call cst_to_H(Vani, l1, l2, cst, ncols, num_s, t1, t2)

    
    ! open(1,file='cst/vani.txt')

    ! do irow = 1,tl1 
    !     do icol = 1, tl2 
    !         if (icol.eq.tl2)then
    !             write(1,'(E15.6)', advance='yes') real(Vani(irow,icol))
    !         else
    !             write(1,'(E15.6)', advance='no') real(Vani(irow,icol))
    !         endif
    !     enddo
    ! enddo 

    ! do irow = 1,tl1 
    !     do icol = 1, tl2 
    !         if (icol.eq.tl2)then
    !             write(1,'(E15.6)', advance='yes') aimag(Vani(irow,icol))
    !         else
    !             write(1,'(E15.6)', advance='no') aimag(Vani(irow,icol))
    !         endif
    !     enddo
    ! enddo 
    ! close(1)





    !call save_Vani_matrix(l1,l2, Vani, './cst/Hmat.txt')
    call Hcomplex_to_cst(Vani, l1, l2, cst, ncols, num_s, t1, t2, 1)


    !Print Cst matrix:
    do is = 1, num_s
        s = smin + is - 1

        write(*,*)'s = ', s
        do it = 1, 2*s + 1
            tval = it - s - 1 

            write(*,*) cst(is,it)
        enddo 
        write(*,*)
    enddo 

end program 
