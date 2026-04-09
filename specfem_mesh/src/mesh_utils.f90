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


! subroutine find_row_col(i, row, col, l)
!     ! Find row and column related to index, i, for a subset of the upper triangular
!     ! matrix (2l+1) x (2l+1), where only columns >= l are included.
!     integer :: i, row, col, l, origi, inew

    
!     ! Store i
!     origi = i 

!     if(i.le. l*(l+1) +1)then
!         !normal case for upper right square (l x l+1) + central element
!         ! Iterate through rows to find the corresponding row and column
!         call myrowcol(i, row, col, l)
!     else
!         inew = origi - (l*(l+1) +1)
!         call upper_triangular_index_to_ij(inew, row, col, l)
!         row = row + l+1
!         col = col + l+1
!     endif     


! end subroutine find_row_col


subroutine find_row_col(i, row, col, l)
    ! Updated column/row search for the triangular (ignoring Yellow shading)
    ! from Eaton thesis figure. 
    ! Ordering for a 2 l +1 matrix: 
    ! example below l = 3 so 16 elements 
    ! .  .  .  .  .  .  1
    ! .  .  .  .  .  2  3
    ! .  .  .  .  4  5  6
    ! .  .  .  7  8  9 10
    ! .  .  .  . 13 12 11
    ! .  .  .  .  . 15 14
    ! .  .  .  .  .  . 16
    ! note the reflection in the lower indices: this is so we dont
    ! need a look up table - can be deterministic 
    implicit none
    integer, intent(in)  :: i, l
    integer, intent(out) :: row, col

    integer  :: n, top_count, r, k, pos, i_bot, i_bot_prime
    real(8)  :: disc

    n         = 2*l + 1
    top_count = (l + 1) * (l + 2) / 2   ! number of elements in first l+1 rows

    if (i <= top_count) then
        !--- First l+1 rows: filled left-to-right ---
        ! Row r contains elements r(r-1)/2 + 1 ... r(r+1)/2
        ! Invert: r = ceil( (-1 + sqrt(1 + 8i)) / 2 )
        disc = (-1.0d0 + sqrt(1.0d0 + 8.0d0 * real(i, 8))) / 2.0d0
        r    = ceiling(disc)
        pos  = i - r * (r - 1) / 2          ! 1-based position within the row
        row  = r
        col  = n - r + pos                   ! row r starts at column n-r+1

    else
        !--- Final l rows: filled right-to-left ---
        ! Map to a reversed local index so it mirrors the top-section pattern
        i_bot       = i - top_count           ! 1-based index within the bottom block
        i_bot_prime = l * (l + 1) / 2 - i_bot + 1   ! reverse to get triangle index

        disc = (-1.0d0 + sqrt(1.0d0 + 8.0d0 * real(i_bot_prime, 8))) / 2.0d0
        k    = ceiling(disc)
        pos  = i_bot_prime - k * (k - 1) / 2  ! 1-based position within the group
        row  = n + 1 - k                       ! row k from the bottom
        col  = n - k + pos                     ! same column formula as top
    end if
end subroutine find_row_col




subroutine myrowcol(i, row, col, l)
    integer :: i, row, col, l
    integer :: r, num_elements_in_row, col_start

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
endsubroutine 


subroutine upper_triangular_index_to_ij(k, i, j, l)
    implicit none
    integer :: l       ! size of the matrix (l x l)
    integer :: k       ! index in the flattened upper triangular
    integer :: i, j   ! row and column indices (1-based)
    integer :: count, row, col

    count = 0
    do row = 1, l
        do col = row, l
            count = count + 1
            if (count == k) then
                i = row
                j = col
                return
            end if
        end do
    end do

    ! If k is invalid
    print *, 'Error: Index out of bounds in upper triangular matrix.'
    i = -1
    j = -1
end subroutine upper_triangular_index_to_ij



end module mesh_utils
