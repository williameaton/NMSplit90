module spline 
  implicit none
  include "constants.h" 


  contains 



  subroutine load_rho_spline(rvals, spline, length, iproc)
    use params, only: datadir
    implicit none 
    integer :: length, iproc
    real(kind=CUSTOM_REAL) :: rvals(length)
    real(kind=SPLINE_REAL) :: spline(length)

    ! Local 
    integer :: stored_len
    character(len=250) :: fname 

    call create_rhospline_fname(iproc, fname)

    open(1, file=trim(datadir)//'/store/rho/'//trim(fname), form='unformatted')

    ! Read the data from the binary file
    read(1) stored_len

    if(stored_len.ne.length)then 
      write(*,*)'Error: the length parsed to load_rho_spline was different from that stored in the binary you are trying to load: '
      write(*,*)'length parsed: ', length
      write(*,*)'stored value : ', stored_len
      stop 
    endif

    read(1)rvals
    read(1)spline
    close(1)
  end subroutine load_rho_spline





  subroutine interpolate_mode_eigenfunctions(mode_type, u, v, du, dv, min_knot, max_knot, & 
                                             new_rad, nout, interp_map, & 
                                             u_spline, v_spline, du_spline, dv_spline)
    use params, only: IC_ID, rad_mineos, NL
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    implicit none

    ! IO variables
    integer :: nout, min_knot, max_knot
    character(len=1)        :: mode_type     ! mode type; S or T
    real(kind=SPLINE_REAL)            :: u(NL), du(NL)
    real(kind=SPLINE_REAL)            :: v(NL), dv(NL)
    real(kind=CUSTOM_REAL)            :: new_rad(nout)
    integer                           :: interp_map(nout)
    real(kind=SPLINE_REAL) :: u_spline(nout), v_spline(nout), du_spline(nout), dv_spline(nout)

    integer :: nin 

    ! Length of the input array 
    nin = max_knot-min_knot+1

    u_spline  = SPLINE_ZERO
    du_spline = SPLINE_ZERO
 
    ! If only 1 or 2 points to 'interpolate' then just return them 
    select case (nout)
    case(1)
      u_spline(1) = u(min_knot)
    case(2)
      u_spline(1) = u(min_knot)
      u_spline(2) = u(max_knot)
    case default

      call cubic_spline_interp(nin, rad_mineos(min_knot:max_knot), u(min_knot:max_knot), & 
                              nout, new_rad, u_spline, interp_map, min_knot)



      call cubic_spline_interp(nin, rad_mineos(min_knot:max_knot), du(min_knot:max_knot), & 
                              nout, new_rad, du_spline, interp_map, min_knot) 

    end select




    ! Interpolate v, dv if spheroidal
    if(mode_type.eq.'S')then 
      v_spline    = SPLINE_ZERO
      dv_spline   = SPLINE_ZERO
      
      select case (nout)
      case(1)
        v_spline(1) = v(min_knot)
      case(2)
        v_spline(1) = v(min_knot)
        v_spline(2) = v(max_knot)
      case default
        call cubic_spline_interp(nin, rad_mineos(min_knot:max_knot), v(min_knot:max_knot), & 
        nout, new_rad, v_spline, interp_map, min_knot)

        call cubic_spline_interp(nin, rad_mineos(min_knot:max_knot), dv(min_knot:max_knot), & 
        nout, new_rad, dv_spline, interp_map, min_knot)
      end select
    endif 

  end subroutine interpolate_mode_eigenfunctions








  


  subroutine interpolate_mineos_variable(in_variable, min_knot, max_knot, new_rad, nout, spl_out, interp_map)
    use params, only: rad_mineos
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    implicit none

    ! IO variables
    integer :: nout, min_knot, max_knot
    real(kind=CUSTOM_REAL)  :: new_rad(nout)
    integer                 :: interp_map(nout)
    real(kind=SPLINE_REAL) :: in_variable(max_knot-min_knot+1)
    real(kind=SPLINE_REAL) :: spl_out(nout)

    integer :: nin

    spl_out    = SPLINE_ZERO
    ! Length of the input array 
    nin = max_knot-min_knot+1

    select case(nout)
      case(1)
        spl_out(1) = in_variable(min_knot)
      case(2)
        spl_out(1) = in_variable(min_knot)
        spl_out(2) = in_variable(max_knot)
      case default
        call cubic_spline_interp(nin, rad_mineos(min_knot:max_knot), in_variable(min_knot:max_knot), & 
                                nout, new_rad, spl_out, interp_map, min_knot)
      end select 
  end subroutine interpolate_mineos_variable




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
            if (radial_arr(i_unq) .lt. rad_mineos(i_knot))then 
                ! This should be the first time that the knot is above the
                ! radius in question and so we want the i_knot - 1 to be stored
                map_arr(i_unq) = i_knot - 1 

                exit 
            else if (radial_arr(i_unq) .eq. rad_mineos(i_knot))then 
                ! It will equal a knot at the boundaries 
                ! If at lower boundary (iknot = min_knot_id) then 
                ! we need to use the min_knot_id value, else can use the 
                ! lower end (i.e. the point is at the far (upper) end)
                ! of this piecewise continuous section 
                ! In this case set it equal to that knot with the assumption the spline interpolation
                ! will hold at the the boundaries of each knot
                if (i_knot .eq. min_knot_id) then 
                  map_arr(i_unq) = i_knot 
                else
                  map_arr(i_unq) = i_knot - 1 
                endif
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

    use params, only: rad_mineos, all_warnings

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
       j = interp_map(i)
      
        do m = 0, k
            interp(i) = interp(i) +  c(m+1,j- min_knot +1)*(real(outradial(i) - rad_mineos(j), kind=SPLINE_REAL))**(real(k-m, kind=SPLINE_REAL)) 
        enddo
    enddo 


  end subroutine cubic_spline_interp



  subroutine write_mode_spline(n, mode_type, l, rad_arr, npoints, name)
    ! Output the spline values: 
    ! Save eigenfunctions to text file in column format 
    use params, only: u_spl, udot_spl, v_spl, vdot_spl, verbose
    implicit none 
    ! IO variables: 
    integer :: n, l 
    character :: mode_type
    real(kind=CUSTOM_REAL) :: rad_arr(npoints)
    integer :: npoints
    character(len=*), optional :: name

    ! Local :
    integer :: i
    character(len=30) :: eigstring
    character(len=2)  :: nfmt, lfmt
    character(len=13)  :: fmtstring


    if (len_trim(name) .eq. 0)then
      if(n.ge.0 .and. n.lt.10)then
        nfmt = 'i1'
      elseif(n.ge.10 .and. n.lt.100)then
        nfmt = 'i2'
      else
        write(*,*)'Format not set for n = ', n
      endif
      if(l.ge.0 .and. l.lt.10)then
        lfmt = 'i1'
      elseif(l.ge.10 .and. l.lt.100)then
        lfmt = 'i2'
      else
        write(*,*)'Format not set for l = ', l
      endif
      write(fmtstring, '(a)') '(a,' // nfmt // ',a,' // lfmt // ',a)'
      write(eigstring,fmtstring) 'spline_', n, mode_type, l, '.txt'
    else 
      write(eigstring,'(a)')name
    endif 

    if (verbose.ge.3)then
      write(*,'(/,a)')'Saving spline to '//trim(eigstring)
    endif 
    

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
