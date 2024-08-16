subroutine cubic_spline_interp(n, x, y, interp)
  ! Interpolate a radial eigenfunction (y) of length n specified at 
  ! radial points x, to the unique_r radii of the mesh
  use params, only: interp_id_r, n_unique_rad, unique_r, rad_mineos

  implicit none 
  include "constants.h"

  ! I/O variables
  integer :: n 
  real(kind=CUSTOM_REAL) :: x(n), y(n)
  real(kind=CUSTOM_REAL) :: interp(n_unique_rad)

  ! Local variables
  real(kind=CUSTOM_REAL):: dx(n-1),dy(n-1), slope(n-1), A(3, n), b(n), d, & 
                           du(n-1), dd(n), dl(n-1), bb(n), t(n-1), c(4,n-1), xp
  integer :: i, info, j, m, k 

  ! Compute dx
  do i = 1, n-1
    dx(i) = x(i+1)-x(i)
    dy(i) = y(i+1)-y(i)
  enddo

  slope(:) = dy(:)/dx(:)

  ! Compute A 
  A(:,:) = zero
  do i = 1, n-2
    A(1, 2+i) = dx(i)                     ! The upper diagonal
    A(2, i+1) = two * (dx(i) + dx(i+1) )  ! The diagonal
    A(3, i) = dx(i+1)                     ! The lower diagonal
  enddo 

  ! Compute b
  b(:) = zero
  do i = 2, n-1
    b(i) = THREE * (dx(i) * slope(i-1) + dx(i-1) * slope(i))
  enddo 

  ! LHS not a knot BC: 
  A(2, 1) = dx(2)
  A(1, 2) = x(3) - x(1)
  d = x(3) - x(1)
  b(1) = ((dx(1) + TWO * d) * dx(2) * slope(1) + (dx(1)**TWO) * slope(2))/d

  ! RHS not a knot BC: 
  A(2, n) = dx(n-2)
  A(3, n-1) = x(n) - x(n-2)
  d = x(n) - x(n-2)
  b(n) = (dx(n-1)**TWO * slope(n-2) + &
            (TWO*d + dx(n-1))*dx(n-2)*slope(n-1))/d

  ! Now we need to solve this Ax = b where A is tridiagonal
  ! dimension, num RHS, subdiagonal (n-1), diagonal(n), superdiag(n-1)

  ! Safe copies that are over-written
  dl(:) = A(3, 1:n-1)
  du(:) = A(1, 2:n)
  dd(:) = A(2,:)
  bb(:) = b(:)

  ! Solves tridiagional A eqn Ax = b for x and stores it in bb
  call sgtsv(n, 1, dl, dd, du, bb, n, info )

  !bb now holds what is referred to as dyxy in the superinit
  ! This is starting to lose some of its accuracy relative to the scipy
  t(:) = (bb(1:n-1) + bb(2:n) - 2*slope(:))/dx(:)

  ! This is the c held in the SciPy cubic spline object
  c(1,:) = t/dx
  c(2,:) = ((slope - bb(1:n-1))/dx) - t
  c(3,:) = bb(1:n-1)
  c(4,:) = y(1:n-1)

  ! We now can evaluate the polynomial as follows: 
  ! Value of spline at point xp is given by 
  !  sum_m [ c[m, i] * (xp - x[i])**(k-m) ] for m = 0 to k
  ! where x[i] is the knot below the point xp 
  k = 3
  interp(:) = zero
  do i = 1, n_unique_rad
    j = interp_id_r(i)
      do m = 0, k
          interp(i) = interp(i) + c(m+1,j)*(unique_r(i) - rad_mineos(j))**(k-m) 
      enddo
  enddo 

end subroutine cubic_spline_interp