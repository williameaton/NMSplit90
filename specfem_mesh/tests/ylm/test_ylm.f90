program test_ylm 
    ! Tests accuracy of Ylm computed
    use ylm_plm
    implicit none 

    integer :: l, m 
    real(kind=CUSTOM_REAL) :: theta, phi , rr 


    theta = 0.5 
    phi   = 0.124

    ! Write Plms for testing
    open(unit=1,file='ylm/testylm.txt', &
    status='unknown',form='formatted',action='write')
    do l = 0,10
        do m = - l, + l
            rr = ylm_real(l, m, theta, phi)
            write(1,*)l, m, rr
        enddo 
    enddo
    close(1)


end program test_ylm