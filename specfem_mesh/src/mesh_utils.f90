module mesh_utils
    
    use math, only: cosp, atan2p, acosp

    implicit none 
    include "constants.h"


    contains 
    
    real(kind=SPLINE_REAL) function delta_spline(i,j)
    ! Delta function with same precision as eigenfunctions
    integer :: i,j
        if (i.eq.j)then 
            delta_spline = ONE
            return  
        else
            delta_spline = ZERO
            return  
        endif 
    end function delta_spline


    integer function delta_int(i,j)
    integer :: i,j
        if (i.eq.j)then 
            delta_int = 1
            return  
        else
            delta_int = 0
            return  
        endif 
    end function delta_int


subroutine find_row_col(i, row, col, l)
    ! Find row and column related to index, i, for a subset of the upper triangular
    ! matrix (2l+1) x (2l+1), where only columns >= l are included.
    integer :: i, row, col, l
    integer :: r, num_elements_in_row, col_start
    
    ! Iterate through rows to find the corresponding row and column
    do r = 1, 2*l+1
        col_start = max(l + 1, r)  ! The first column to consider in this row
        num_elements_in_row = 2*l + 1  - col_start + 1
        
        if (i <= num_elements_in_row) then
            row = r
            col = i + col_start - 1
            return
        else
            i = i - num_elements_in_row
        end if
    end do
    
end subroutine find_row_col


end module mesh_utils
