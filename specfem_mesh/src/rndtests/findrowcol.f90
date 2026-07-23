program test_index
    implicit none
    integer :: i, l, row, col

    do l = 1, 3
        write(*,'(A,I0)') "=== l = ", l
        do i = 1, (l+1)**2  ! total elements = (l+1)^2
            call find_row_col(i, row, col,  l)
            write(*,'(A,I3,A,I3,A,I3)') "  i=", i, "  row=", row, "  col=", col
        end do
    end do

end program test_index


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

