
subroutine lagrange_any(xi,NGLL,xigll,h,hprime)
    ! subroutine to compute the Lagrange interpolants based upon the interpolation points
    ! and their first derivatives at any point xi in [-1,1]
    
      implicit none
    
      double precision,intent(in) :: xi
    
      integer,intent(in) :: NGLL
      double precision,dimension(NGLL),intent(in) :: xigll
      double precision,dimension(NGLL),intent(out) :: h, hprime
    
      ! local parameters
      integer :: dgr,i,j
      double precision :: prod1,prod2,prod3
      double precision :: prod2_inv
      double precision :: sum
      double precision :: x0,x
      double precision :: x_1,x_2,x_3,x_4,x_5
    
    ! note: this routine is hit pretty hard by the mesher, optimizing the loops here will be beneficial
    
      select case(NGLL)
      case (5)
        ! NGLL == 5
    
        ! loads positions
        x_1 = xigll(1)
        x_2 = xigll(2)
        x_3 = xigll(3)
        x_4 = xigll(4)
        x_5 = xigll(5)
    
        ! for dgr == 1
        ! lagrangian interpolants
        ! do i = 1,5 - skips i == 1
        prod1 = (xi - x_2) * (xi - x_3) * (xi - x_4) * (xi - x_5)
        prod2 = (x_1 - x_2) * (x_1 - x_3) * (x_1 - x_4) * (x_1 - x_5)
        ! takes inverse to avoid additional divisions
        ! (multiplications are cheaper than divisions)
        prod2_inv = 1.d0 / prod2
        h(1) = prod1 * prod2_inv
        ! first derivatives
        ! do i = 1,5 - skips i == 1
        !   do j = 1,5 - skips j == 1 and j == 2
        sum = (xi - x_3) * (xi - x_4) * (xi - x_5) &
            + (xi - x_2) * (xi - x_4) * (xi - x_5) &
            + (xi - x_2) * (xi - x_3) * (xi - x_5) &
            + (xi - x_2) * (xi - x_3) * (xi - x_4)
        ! hprime
        hprime(1) = sum * prod2_inv
    
        ! for dgr == 2:
        ! lagrangian interpolants
        ! do i = 1,5 - skips i == 2
        prod1 = (xi - x_1) * (xi - x_3) * (xi - x_4) * (xi - x_5)
        prod2 = (x_2 - x_1) * (x_2 - x_3) * (x_2 - x_4) * (x_2 - x_5)
        ! takes inverse to avoid additional divisions
        ! (multiplications are cheaper than divisions)
        prod2_inv = 1.d0 / prod2
        h(2) = prod1 * prod2_inv
        ! first derivatives
        ! do i = 1,5 - skips i == 2
        !   do j = 1,5 - skips j == 2 and j == 1
        sum = (xi - x_3) * (xi - x_4) * (xi - x_5) &
            + (xi - x_1) * (xi - x_4) * (xi - x_5) &
            + (xi - x_1) * (xi - x_3) * (xi - x_5) &
            + (xi - x_1) * (xi - x_3) * (xi - x_4)
        ! hprime
        hprime(2) = sum * prod2_inv
    
        ! for dgr == 3
        ! lagrangian interpolants
        ! do i = 1,5 - skips i == 3
        prod1 = (xi - x_1) * (xi - x_2) * (xi - x_4) * (xi - x_5)
        prod2 = (x_3 - x_1) * (x_3 - x_2) * (x_3 - x_4) * (x_3 - x_5)
        ! takes inverse to avoid additional divisions
        ! (multiplications are cheaper than divisions)
        prod2_inv = 1.d0/prod2
        h(3) = prod1 * prod2_inv
        ! first derivatives
        ! do i = 1,5 - skips i == 3
        !   do j = 1,5 - skips j == 3 and j == 1
        sum = (xi - x_2) * (xi - x_4) * (xi - x_5) &
            + (xi - x_1) * (xi - x_4) * (xi - x_5) &
            + (xi - x_1) * (xi - x_2) * (xi - x_5) &
            + (xi - x_1) * (xi - x_2) * (xi - x_4)
        ! hprime
        hprime(3) = sum * prod2_inv
    
        ! for dgr == 4
        ! lagrangian interpolants
        ! do i = 1,5 - skips i == 4
        prod1 = (xi - x_1) * (xi - x_2) * (xi - x_3) * (xi - x_5)
        prod2 = (x_4 - x_1) * (x_4 - x_2) * (x_4 - x_3) * (x_4 - x_5)
        ! takes inverse to avoid additional divisions
        ! (multiplications are cheaper than divisions)
        prod2_inv = 1.d0 / prod2
        h(4) = prod1 * prod2_inv
        ! first derivatives
        ! do i = 1,5 - skips i == 4
        !   do j = 1,5 - skips j == 4 and j == 1
        sum = (xi - x_2) * (xi - x_3) * (xi - x_5) &
            + (xi - x_1) * (xi - x_3) * (xi - x_5) &
            + (xi - x_1) * (xi - x_2) * (xi - x_5) &
            + (xi - x_1) * (xi - x_2) * (xi - x_3)
        ! hprime
        hprime(4) = sum * prod2_inv
    
        ! for dgr == 5
        ! lagrangian interpolants
        ! do i = 1,5 - skips i == 5
        prod1 = (xi - x_1) * (xi - x_2) * (xi - x_3) * (xi - x_4)
        prod2 = (x_5 - x_1) * (x_5 - x_2) * (x_5 - x_3) * (x_5 - x_4)
        ! takes inverse to avoid additional divisions
        ! (multiplications are cheaper than divisions)
        prod2_inv = 1.d0 / prod2
        h(5) = prod1 * prod2_inv
        ! first derivatives
        ! do i = 1,5 - skips i == 5
        !   do j = 1,5 - skips j == 5 and j == 1
        sum = (xi - x_2) * (xi - x_3) * (xi - x_4) &
            + (xi - x_1) * (xi - x_3) * (xi - x_4) &
            + (xi - x_1) * (xi - x_2) * (xi - x_4) &
            + (xi - x_1) * (xi - x_2) * (xi - x_3)
        ! hprime
        hprime(5) = sum * prod2_inv
    
      case default
        ! general NGLL
        do dgr = 1,NGLL
    
          prod1 = 1.0d0
          prod2 = 1.0d0
    
          ! lagrangian interpolants
          x0 = xigll(dgr)
          do i = 1,NGLL
            if (i /= dgr) then
              x = xigll(i)
              prod1 = prod1*(xi-x)
              prod2 = prod2*(x0-x)
            endif
          enddo
    
          ! takes inverse to avoid additional divisions
          ! (multiplications are cheaper than divisions)
          prod2_inv = 1.d0/prod2
    
          h(dgr) = prod1 * prod2_inv
    
          ! first derivatives
          sum = 0.0d0
          do i = 1,NGLL
            if (i /= dgr) then
              prod3 = 1.0d0
              do j = 1,NGLL
                if (j /= dgr .and. j /= i) prod3 = prod3*(xi-xigll(j))
              enddo
              sum = sum + prod3
            endif
          enddo
    
          hprime(dgr) = sum * prod2_inv
    
        enddo
      end select
    
      end subroutine lagrange_any