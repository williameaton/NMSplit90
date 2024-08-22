module spline 
  implicit none
  include "precision.h" 
  contains 


  subroutine interpolate_mode_eigenfunctions(mode_type, u, v, du, dv, new_rad, nlength, interp_map)
    use params, only: u_spl, v_spl, udot_spl, vdot_spl, & 
                      IC_ID, rad_mineos, NL
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    implicit none
    ! IO variables
    integer :: nlength
    character(len=1)        :: mode_type     ! mode type; S or T
    real(kind=SPLINE_REAL)            :: u(NL), du(NL)
    real(kind=SPLINE_REAL)            :: v(NL), dv(NL)
    real(kind=CUSTOM_REAL)            :: new_rad(nlength)
    integer                           :: interp_map(nlength)

    ! Interpolate u, du (or w, dw if toroidal)
    call deallocate_if_allocated(u_spl)
    call allocate_if_unallocated(nlength, u_spl)
    call deallocate_if_allocated(udot_spl)
    call allocate_if_unallocated(nlength, udot_spl)
 
    call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID),  u(1:IC_ID), nlength, new_rad, u_spl,    interp_map)
    call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID), du(1:IC_ID), nlength, new_rad, udot_spl, interp_map)

    ! Interpolate v, dv if spheroidal
    if(mode_type.eq.'S')then 
      call allocate_if_unallocated(nlength, v_spl)
      call allocate_if_unallocated(nlength, vdot_spl)
      call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID),  v(1:IC_ID), nlength, new_rad,    v_spl, interp_map)
      call cubic_spline_interp(IC_ID, rad_mineos(1:IC_ID), dv(1:IC_ID), nlength, new_rad, vdot_spl, interp_map)
    endif 
  end subroutine interpolate_mode_eigenfunctions





  subroutine create_interpolation_radial_map(radial_arr, map_arr, arr_len, min_knot_id, max_knot_id)
    use params, only: rad_mineos
    implicit none 

    ! IO variables
    integer :: arr_len, map_arr(arr_len), min_knot_id, max_knot_id
    real(kind=CUSTOM_REAL) :: radial_arr(arr_len)

    ! Local variables
    integer :: i_unq, i_knot

    ! Find the knot id's of the radii just below the points we will want to interpolate to 
    map_arr(:) = 0
    do i_unq = 1, arr_len
        ! For each radius we want to find the last knot that it is larger than
        do i_knot = min_knot_id, max_knot_id
            if (radial_arr(i_unq) .le. rad_mineos(i_knot))then 
                ! This should be the first time that the knot is above the
                ! radius in question and so we want the i_knot - 1 to be stored
                ! Note that in some cases the radius will equal a knot (e.g. the IC radius)
                ! In this case set it equal to that knot with the assumption the spline interpolation
                ! will hold at the the boundaries of each knot
                map_arr(i_unq) = i_knot - 1 
                exit 
            else
            endif
        enddo 
    enddo 

    ! Check the interp_id_r 
    if(minval(map_arr).lt.min_knot_id .or. maxval(map_arr).gt.max_knot_id)then
        write(*,*)'Error in assigning values to map'
        write(*,*)' -- min value: ', minval(map_arr)
        write(*,*)' -- max value: ', maxval(map_arr)

        write(*,*)' ------ some useful debugging things -----'
        write(*,*)' -- min value of your radial array  : ', minval(radial_arr)
        write(*,*)' -- min value of mineos radial array: ', rad_mineos(min_knot_id)

        write(*,*)' -- max value of your radial array  : ', maxval(radial_arr)
        write(*,*)' -- max value of mineos radial array: ', rad_mineos(max_knot_id)
    endif 



  end subroutine
  
  






  subroutine cubic_spline_interp(n, x, y, nout, outradial, interp, interp_map)
    ! Interpolate a radial eigenfunction (y) of length n specified at 
    ! radial points x, to the unique_r radii of the mesh
    ! n = number of knots in input model 
    ! x = radial/x values for function 
    ! y = y values for function to be interpolated
    ! nout = number of knots/radial points in the output (interpolated) function
    ! interp = array with interpolated values at the nout radial points
    ! interp_map: a map which, for each radial value in out radial, gives the index
    ! of the mineos knot below the point

    use params, only: n_unique_rad, rad_mineos

    implicit none 
    include "constants.h"

    ! I/O variables
    integer :: n , nout , interp_map(nout)
    real(kind=CUSTOM_REAL) :: x(n), outradial(nout)
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
    do i = 1, nout
      j = interp_map(i)
        do m = 0, k
            interp(i) = interp(i) +  c(m+1,j)*(real(outradial(i) - rad_mineos(j), kind=SPLINE_REAL))**(real(k-m, kind=SPLINE_REAL)) 
        enddo
    enddo 


  end subroutine cubic_spline_interp



  subroutine write_mode_spline(n, mode_type, l, rad_arr, npoints)
    ! Output the spline values: 
    ! Save eigenfunctions to text file in column format 
    use params, only: u_spl, udot_spl, v_spl, vdot_spl
    implicit none 
    ! IO variables: 
    integer :: n, l 
    character :: mode_type
    real(kind=CUSTOM_REAL) :: rad_arr(npoints)
    integer :: npoints

    ! Local :
    integer :: i
    character(len=30) :: eigstring


    write(eigstring,'(a,i2,a, i2, a)') 'spline_', n, mode_type, l, '.txt'
    open(1,file=trim(eigstring))
    do i =1, npoints
        if(mode_type=='S')then 
            write(1,*)rad_arr(i), u_spl(i), udot_spl(i), v_spl(i), vdot_spl(i)
        elseif (mode_type=='T')then 
            write(1,*)rad_arr(i), u_spl(i), udot_spl(i)
        endif 
    enddo 
    close(1)



  end subroutine write_mode_spline



end module spline