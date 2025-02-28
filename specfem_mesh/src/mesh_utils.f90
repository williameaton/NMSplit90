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
        ! Find row and column related to index, i, of an upper triangular
        ! of square matrix (2l+1) x (2l+1) where the index is contiguous
        ! across a row 
        integer :: i, row, col, l
        integer :: r, num_elements_in_row
        
        ! Find the row for the given index i
        num_elements_in_row = 0
        do r = 1, 2*l+1
            num_elements_in_row = 2*l + 2 - r
            if (i <= num_elements_in_row) then
                row = r
                col = i + r - 1
                return
            else
                i = i - num_elements_in_row
            end if
        end do
        end subroutine find_row_col


end module mesh_utils
