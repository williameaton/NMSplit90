program test_ylm_values
    ! Tests accuracy of Ylm computed against the table B.1 in DT98
    use ylm_plm
    use math, only: cosp, sinp
    implicit none 

    integer :: l, m 
    real(kind=CUSTOM_REAL) :: theta, phi , sinth, costh
    complex(kind=CUSTOM_REAL) ::  y, ei1p, ei1m, ei2p, ei2m, ei3p, ei3m, err


    theta = 0.5 
    phi   = 0.124

    ei1p = cmplx(cosp(phi), sinp(phi))
    ei1m = cmplx(cosp(-phi), sinp(-phi))
    ei2p = cmplx(cosp(two*phi), sinp(phi*two))
    ei2m = cmplx(cosp(-phi*two), sinp(-phi*two))
    ei3p = cmplx(cosp(three*phi), sinp(phi*three))
    ei3m = cmplx(cosp(-phi*three), sinp(-phi*three))

    sinth = sinp(theta)
    costh = cosp(theta)


    open(unit=1,file='ylm_value_error.txt', &
    status='unknown',form='formatted',action='write')
    
    
    ! l = 0    m = 0
    y = ylm_complex(0, 0, theta, phi)
    err =  y -  one/(two*(PI**half))
    write(1,*) real(err), aimag(err)

    ! l = 1    m = -1
    y = ylm_complex(1, -1, theta, phi)
    err =   y - (sinth *  half * (three/TWO_PI)**half * ei1m)
    write(1,*) real(err), aimag(err)

    ! l = 1    m = 0
    y = ylm_complex(1, 0, theta, phi)
    err =  y - (half * costh * (three/pi)**half )
    write(1,*) real(err), aimag(err)

    ! l = 1    m = 1
    y = ylm_complex(1, 1, theta, phi)
    err =  y - ( ei1p * sinth *  -half * (three/TWO_PI)**half )
    write(1,*) real(err), aimag(err)

  
    ! l = 2    m = -2
    y = ylm_complex(2, -2, theta, phi)
    err =   y - (((five*three/TWO_PI)**half) * sinth * sinth * ei2m  / four)
    write(1,*) real(err), aimag(err)

    ! l = 2    m = -1
    y = ylm_complex(2, -1, theta, phi)
    err =  y - (half * ((five*three/TWO_PI)**half) * sinth * costh * ei1m)
    write(1,*) real(err), aimag(err)

    ! l = 2    m = 0
    y = ylm_complex(2, 0, theta, phi)
    err =  y -  ((five/pi)**half) * (three*(cosp(theta)**two) - one )/four
    write(1,*) real(err), aimag(err)

    ! l = 2    m = +1
    y = ylm_complex(2, 1, theta, phi)
    err =  y - ( -half * ((five*three/TWO_PI)**half) * sinth * costh * ei1p)
    write(1,*) real(err), aimag(err)

    
    ! l = 2    m = +2
    y = ylm_complex(2, +2, theta, phi)
    err =   y - (((five*three/TWO_PI)**half) * sinth * sinth * ei2p  / four)
    write(1,*) real(err), aimag(err)


   ! l = 3    m = -3
    y = ylm_complex(3, -3, theta, phi)
    err =  y - ( (seven*five/PI)**half  *  sinth**three * ei3m / eight)
    write(1,*) real(err), aimag(err)


   ! l = 3    m = -2
    y = ylm_complex(3, -2, theta, phi)
    err =  y - ((seven*three*five/TWO_PI)**half  * costh * sinth*sinth * ei2m  / four)
    write(1,*) real(err), aimag(err)


   ! l = 3    m = -1
    y = ylm_complex(3, -1, theta, phi)
    err =  y - (((seven*three/PI)**half) * sinth * (five*costh*costh - one) * ei1m / eight)
    write(1,*) real(err), aimag(err)


   ! l = 3    m = 0
    y = ylm_complex(3, 0, theta, phi)
    err =  y - ((seven/pi)**half  *costh*(five*costh*costh - three)/four)
    write(1,*) real(err), aimag(err)


    ! l = 3    m = -1
    y = ylm_complex(3, 1, theta, phi)
    err =  y - ( -((seven*three/PI)**half) * sinth * (five*costh*costh - one) * ei1p / eight)
    write(1,*) real(err), aimag(err)


   ! l = 3    m = 2
    y = ylm_complex(3, 2, theta, phi)
    err =  y - ((seven*three*five/TWO_PI)**half  * costh * sinth*sinth * ei2p  / four)
    write(1,*) real(err), aimag(err)


    ! l = 3    m =  3
    y = ylm_complex(3, 3, theta, phi)
    err =  y - ( -(seven*five/PI)**half  *  sinth**three * ei3p / eight)
    write(1,*) real(err), aimag(err)


    close(1)


end program test_ylm_values