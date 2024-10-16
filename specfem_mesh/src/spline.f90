module spline 
  implicit none
  include "constants.h" 


  contains 
 

  subroutine quad_spline_interp_3(x_in, y_in, n_in, x_out, y_out, n_out)
      ! Very basic spline interpolation if n = 3
      implicit none 
      
      integer :: n_in, n_out
      real(kind=CUSTOM_REAL) :: x_in(n_in), x_out(n_out)
      real(kind=SPLINE_REAL) :: y_in(n_in), y_out(n_out)

      real(kind=SPLINE_REAL) ::  A(6,6), b(6)
      real(kind=SPLINE_REAL) :: x0, x1, x2, y0, y1, y2, xx
      integer :: junk(6), info, i

      if (n_in.ne.3)then 
        write(*,*)'Using quad_spline when n not =  3. Stop.'
        stop
      endif 

      A = SPLINE_ZERO
      b = SPLINE_ZERO

      x0 = real(x_in(1), kind=SPLINE_REAL)
      x1 = real(x_in(2), kind=SPLINE_REAL)
      x2 = real(x_in(3), kind=SPLINE_REAL)
      y0 = y_in(1)
      y1 = y_in(2)
      y2 = y_in(3)

      A(1,1) = x0**SPLINE_TWO
      A(1,2) = x0
      A(1,3) = SPLINE_ONE
      
      A(2,1) = x1**SPLINE_TWO
      A(2,2) = x1
      A(2,3) = SPLINE_ONE
      
      A(3,1) = SPLINE_TWO*x1
      A(3,2) = SPLINE_ONE
      A(3,4) = -SPLINE_TWO*x1
      A(3,5) = -SPLINE_ONE
      
      A(4,4) = x1**SPLINE_TWO
      A(4,5) = x1
      A(4,6) = SPLINE_ONE

      A(5,4) = x2**SPLINE_TWO
      A(5,5) = x2
      A(5,6) = SPLINE_ONE

      A(6,2) = SPLINE_ONE

      b(1) = y0
      b(2) = y1
      b(4) = y1
      b(5) = y2


      call sgesv(6, 1, A, 6, junk, b, 6, info)
      if(info.ne.0)then
        write(*,*)'Error in quad_spline sgesv. stop'
        stop
      endif

      y_out = zero
      do i = 1, n_out
        xx = real(x_out(i), kind=SPLINE_REAL)
        if (xx.ge.x0 .and. xx.lt.x1) then 
          y_out(i) = b(1)*xx*xx + b(2)*xx + b(3)
        else if (xx.ge.x1 .and. xx.le.x2) then 
          y_out(i) = b(4)*xx*xx + b(5)*xx + b(6)
        else
          write(*,*)'Error in quad_spline: x_out trying to extrapolate.'
          stop
        endif
      enddo 


  end subroutine quad_spline_interp_3







  
  subroutine cubic_spline_interp(n, x, y, nout, outradial, interp, interp_map, min_knot)
    ! Interpolate a radial eigenfunction (y) of length n specified at 
    ! radial points x, to the unique_r radii of the mesh
    ! n = number of knots in input model 
    ! x = radial/x values for function 
    ! y = y values for function to be interpolated
    ! nout = number of knots/radial points in the output (interpolated) function
    ! interp = array with interpolated values at the nout radial points
    ! interp_map: a map which, for each radial value in out radial, gives the index
    ! of the mineos knot below the point
    ! min_knot = the first knot iD (e.g. in mineos_radius)

    use params, only: all_warnings

    implicit none 
    include "constants.h"

    ! I/O variables
    integer :: n , nout , interp_map(nout), min_knot
    real(kind=CUSTOM_REAL) :: x(n), outradial(nout)
    real(kind=SPLINE_REAL) :: y(n)
    real(kind=SPLINE_REAL) :: interp(nout)

    ! Local variables
    real(kind=CUSTOM_REAL) :: dx(n-1)
    real(kind=SPLINE_REAL) :: dy(n-1), slope(n-1), dxs(n-1), A(3, n), & 
                              b(n), d, du(n-1), dd(n), dl(n-1), bb(n), & 
                              xs(n), t(n-1), c(4,n-1)
    integer :: i, info, j, m, k 

    if(n.lt.2)then
      write(*,*)'Error in cubic_spline_inter: trying to interpolate given x and y but they have length n < 2'
      stop
    else if (n.eq.2) then 
      if (all_warnings) write(*,*)'Warning in cubic_spline_interp: n = 2 - only linearly interpolating'
      call linear_interpolate(x, y, nout, outradial, interp)
      return
    else if (n.eq.3) then 
      call quad_spline_interp_3(x, y, n, outradial, interp, nout)
      return
    endif

    ! Compute dx
    do i = 1, n-1
      dx(i) = x(i+1)-x(i)
      dy(i) = y(i+1)-y(i)
    enddo

    xs       = real(x(:),  kind=SPLINE_REAL)
    dxs(:)   = real(dx(:), kind=SPLINE_REAL)
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
    call sgtsv(n, 1, dl, dd, du, bb, n, info)

    if (info.ne.0)then 
      write(*,*)'Error in sgtsv, info =', info
      stop
    endif


    ! bb now holds what is referred to as dyxy in the superinit
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
    do i = 1, nout
       j = interp_map(i) - min_knot +1
      
        do m = 0, k
            interp(i) = interp(i) +  c(m+1,j) * & 
              (real(outradial(i) - x(j), & 
                  kind=SPLINE_REAL))**(real(k-m, kind=SPLINE_REAL)) 
        enddo
    enddo 


  end subroutine cubic_spline_interp



  subroutine write_scalar_spline(rad_arr, spline, npoints, eigstring)
    ! Output the spline values: 
    use params, only: verbose
    implicit none 

    ! IO variables: 
    real(kind=CUSTOM_REAL) :: rad_arr(npoints)
    real(kind=SPLINE_REAL) :: spline(npoints)
    integer :: npoints
    character(len=*) :: eigstring

    integer :: i 

    if (verbose.ge.3)then
      write(*,'(/,a)')'Saving spline to '//trim(eigstring)
    endif 
    

    open(1,file=trim(eigstring))
    do i =1, npoints
        write(1,*)rad_arr(i), spline(i)
    enddo 
    close(1)
  end subroutine write_scalar_spline





  subroutine linear_interpolate(x, y, nout, out_x, out_y)
    ! Linearly interpolates function y(x) specified at [x1, y1] and [x2,y2] 
    ! at the outradial x values, stored in interp
    implicit none

    ! IO variables
    integer :: nout
    real(kind=CUSTOM_REAL) :: x(2), out_x(nout)
    real(kind=SPLINE_REAL) :: y(2), out_y(nout)
    
    ! local: 
    integer :: i
    real(kind=SPLINE_REAL) :: slope, x1, x2, xtmp, y1

    ! Checks
    if (x(2).le.x(1)) then 
      write(*,*)'Error in linear_interpolate: x(2) is larger than x(1): '
      write(*,*)'   x = : ', x
      stop
    endif 
    
    y1 = real(y(1), kind=SPLINE_REAL)
    x1 = real(x(1), kind=SPLINE_REAL)
    x2 = real(x(2), kind=SPLINE_REAL)

    slope = (y(2)-y(1))/(x2-x1)

    do i = 1, nout

      xtmp = real(out_x(i), kind=SPLINE_REAL)
      if (xtmp.lt.x1 .or. xtmp.gt.x2) then
        write(*,*)'Error in linear_interpolate: trying to interpolate at'
        write(*,*)'x =', xtmp
        write(*,*)'but function specified at ', x
        stop
      else 
        out_y(i) = y1 + slope * (xtmp - x1)
      endif 

    enddo 

  end subroutine linear_interpolate


end module spline
