program test_xlm 
    ! Tests accuracy of Ylm computed
    use ylm_plm
    use integrate, only: integrate_r_traps
    use math, only: sinp
    implicit none 

    integer :: l1, l2, m, npoints, i
    real(kind=CUSTOM_REAL) :: integral
    real(kind=CUSTOM_REAL), allocatable :: theta(:), x(:)
   
    m  = 1
    npoints = 10000
    allocate(theta(npoints))
    allocate(x(npoints))

    ! Write integral to file for pytest 
    open(unit=1,file='xlm/xlm_integral.txt', &
    status='unknown',form='formatted',action='write')

    do l1 = 1, 12
        do l2 = 3, 10
            do i = 1, npoints
                ! Create an array between 0 and pi 
                theta(i) = (real((i-1),kind=CUSTOM_REAL)/real((npoints-1),kind=CUSTOM_REAL)) * PI  
                
                ! Compute integrand (B.59)
                x(i) = xlm(l1, m, theta(i)) * xlm(l2, m, theta(i))  * sinp(theta(i))
            enddo 

            ! Output
            integral = integrate_r_traps(theta, x, npoints)
            write(1,*)l1, l2, integral
        enddo 
    enddo

    close(1)

end program test_xlm