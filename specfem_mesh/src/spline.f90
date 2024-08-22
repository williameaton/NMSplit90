module spline 
  implicit none
  include "precision.h" 
  contains 


  subroutine interpolate_mode_eigenfunctions(mode_type, u, v, du, dv)
    use params, only: n_unique_rad, u_spl, v_spl, udot_spl, vdot_spl, & 
                      IC_ID, rad_mineos, NL
    use allocation_module, only: allocate_if_unallocated
    implicit none

    ! IO variables
    character(len=1)        :: mode_type     ! mode type; S or T
    real(kind=SPLINE_REAL)            :: u(NL), du(NL)
    real(kind=SPLINE_REAL)            :: v(NL), dv(NL)

    ! Interpolate u, du (or w, dw if toroidal)
    call allocate_if_unallocated(n_unique_rad, u_spl)
    call allocate_if_unallocated(n_unique_rad, udot_spl)
 
    call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID),  u(1:IC_ID),    u_spl)
    call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID), du(1:IC_ID), udot_spl)

    ! Interpolate v, dv if spheroidal
    if(mode_type.eq.'S')then 
      call allocate_if_unallocated(n_unique_rad, v_spl)
      call allocate_if_unallocated(n_unique_rad, vdot_spl)
      call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID),  v(1:IC_ID),    v_spl)
      call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID), dv(1:IC_ID), vdot_spl)
    endif 

  end subroutine interpolate_mode_eigenfunctions




  subroutine cubic_spline_interp(n, x, y, interp)
    ! Interpolate a radial eigenfunction (y) of length n specified at 
    ! radial points x, to the unique_r radii of the mesh
    use params, only: interp_id_r, n_unique_rad, unique_r, rad_mineos

    implicit none 
    include "constants.h"

    ! I/O variables
    integer :: n 
    real(kind=CUSTOM_REAL) :: x(n)
    real(kind=SPLINE_REAL) :: y(n)
    real(kind=SPLINE_REAL) :: interp(n_unique_rad)

    ! Local variables
    real(kind=CUSTOM_REAL) :: dx(n-1)
    real(kind=SPLINE_REAL) :: dy(n-1), slope(n-1), dxs(n-1), A(3, n), & 
                              b(n), d, du(n-1), dd(n), dl(n-1), bb(n), & 
                              xs(n), t(n-1), c(4,n-1)
    integer :: i, info, j, m, k 


    ! Compute dx
    do i = 1, n-1
      dx(i) = x(i+1)-x(i)
      dy(i) = y(i+1)-y(i)
    enddo


    xs     =  real(x(:),kind=SPLINE_REAL)
    dxs(:) = real(dx(:),kind=SPLINE_REAL)
    slope(:) = dy(:)/dxs

    ! Compute A 
    A(:,:) = SPLINE_ZERO
    do i = 1, n-2
      A(1, 2+i) = dxs(i)                             ! The upper diagonal
      A(2, i+1) = SPLINE_TWO * (dxs(i) + dxs(i+1) )  ! The diagonal
      A(3, i)   = dxs(i+1)                           ! The lower diagonal
    enddo 

    ! Compute b
    b(:) = SPLINE_ZERO
    do i = 2, n-1
      b(i) = SPLINE_THREE * (dxs(i) * slope(i-1) + dxs(i-1) * slope(i))
    enddo 

    ! LHS not a knot BC: 
    A(2, 1) = dxs(2)
    A(1, 2) = xs(3) - xs(1)
    d = xs(3) - xs(1)
    b(1) = ((dxs(1) + SPLINE_TWO * d) * dxs(2) * slope(1) + (dxs(1)**SPLINE_TWO) * slope(2))/d

    ! RHS not a knot BC: 
    A(2, n) = dxs(n-2)
    A(3, n-1) = xs(n) - xs(n-2)
    d = xs(n) - xs(n-2)
    b(n) = (dxs(n-1)**SPLINE_TWO * slope(n-2) + &
              (SPLINE_TWO*d + dxs(n-1))*dxs(n-2)*slope(n-1))/d

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
    t(:) = (bb(1:n-1) + bb(2:n) - 2*slope(:))/dxs(:)

    ! This is the c held in the SciPy cubic spline object
    c(1,:) = t/dxs
    c(2,:) = ((slope - bb(1:n-1))/dxs) - t
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

            !interp(i) = interp(i) + c(m+1,j)*(unique_r(i) - rad_mineos(j))**(k-m) 
            interp(i) = interp(i) +  c(m+1,j)*(real(unique_r(i) - rad_mineos(j), kind=SPLINE_REAL))**(real(k-m, kind=SPLINE_REAL)) 
        enddo
    enddo 


  end subroutine cubic_spline_interp

end module spline