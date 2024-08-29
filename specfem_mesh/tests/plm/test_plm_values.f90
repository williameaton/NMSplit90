program test_values
    ! Tests accuracy of Plm
    ! NOTE THAT THESE ARE NOT THE SAME AS DEFINED on WIKI etc
    ! DT 98 defines the Plm to not include the (-1)**m term 
    ! See DT 98 B.48
    use ylm_plm
    use math, only: sinp, cosp
    implicit none 

    integer :: l, m 
    real(kind=CUSTOM_REAL) :: theta, phi , rr, costh, sinth

    theta = pi/three
    sinth = sinp(theta)
    costh = cosp(theta)

    open(unit=1,file='plm_value_error.txt', &
    status='unknown',form='formatted',action='write')
    
    !(1) 0 0 
    rr = Plm(costh,0,0)
    write(1,*)rr - one 


    !(2) 1 -1
    rr = Plm(costh,1,-1)
    write(1,*)rr - ( -half*sinth )

    !(3) 1 0 
    rr = Plm(costh,1,0)
    write(1,*)rr - costh 

    !(4) 1 1
    rr = Plm(costh,1,1)
    write(1,*)rr - sinth 



    !(5) 2 -2 
    rr = Plm(costh, 2, -2)
    write(1,*)rr - (sinth*sinth/eight )

    !(6) 2 -1 
    rr = Plm(costh, 2, -1)
    write(1,*)rr - (- half * sinth * costh)

    !(7) 2 0 
    rr = Plm(costh, 2, 0)
    write(1,*)rr - (half *(three*costh*costh - one) )

    !(8) 2 1
    rr = Plm(costh, 2, 1)
    write(1,*)rr - (three * sinth * costh )

    !(9) 2 2
    rr = Plm(costh, 2, 2)
    write(1,*)rr - (three * sinth * sinth )



    !(10) 3 -3
    rr = Plm(costh, 3, -3)
    write(1,*)rr - (- sinth**3 /(four*three*four) )

    !(11) 3 -2
    rr = Plm(costh, 3, -2)
    write(1,*)rr - (costh*sinth*sinth/eight)


    !(12) 3 -1 
    rr = Plm(costh, 3, -1)
    write(1,*)rr - (sinth*(one - five*costh*costh)/eight )


    !(13) 3 0 
    rr = Plm(costh, 3, 0)
    write(1,*)rr - (half*costh* (five*costh*costh - three) )

    !(14) 3 1 
    rr = Plm(costh, 3, 1)
    write(1,*)rr - ( sinth * (five*costh*costh - one ) * three/two)

    !(15) 3 2
    rr = Plm(costh, 3, 2)
    write(1,*)rr - (five*three*costh*sinth*sinth)

    !(16) 3 3
    rr = Plm(costh, 3, 3)
    write(1,*)rr - (three*five * sinth**three )

    close(1)


end program test_values
