module gll

    include "constants.h"

    contains 

    real(kind=CUSTOM_REAL) function lagrange(N, i, x)
        ! Function to calculate  Lagrange polynomial for order N and polynomial
        ! i[0, N] at location x (not necessarily GLL point)
        implicit none 

        ! IO variables: 
        integer :: N, i
        real(kind=CUSTOM_REAL) :: x
        ! Local: 
        integer :: j, j1, i2
        real(kind=CUSTOM_REAL) ::  xi(N+1), weights(N+1)

        call get_gll(N, xi, weights)
        lagrange = one


        do j = 0, N
            if (j .ne. i+1) then
                j1 = j+1
                i2 = i+2
                lagrange = lagrange * ((x - xi(j1)) / (xi(i2) - xi(j1)))
            endif 
        enddo 
    end function lagrange
    ! ------------------------------------------------------------------

    subroutine compute_wglljac()
        use params, only: wgll, detjac, ngllx, nglly, ngllz, nspec, & 
                          wglljac
        use allocation_module, only: deallocate_if_allocated
        implicit none
        integer :: i, j, k, ispec

        call deallocate_if_allocated(wglljac)
        allocate(wglljac(ngllx,nglly,ngllz,nspec))

        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz 
                        wglljac(i,j,k,ispec) =real(wgll(i) * & 
                                                   wgll(j) * & 
                                                   wgll(k) * & 
                                                   detjac(i,j,k,ispec), & 
                                                   kind=SPLINE_REAL)
                    enddo 
                enddo
            enddo 
        enddo 


    end subroutine compute_wglljac



    subroutine setup_gll()
        use params, only: xi, wgll, ngllx, nglly, ngllz, verbose
        use allocation_module, only: allocate_if_unallocated
        implicit none 


        if (ngllx.ne.nglly .or. ngllx.ne.ngllz .or. nglly.ne.ngllz)then
            write(*,*)'Error: currently only working for case of ngllx = nglly = ngllz'
            stop 
        endif 

        call allocate_if_unallocated(ngllx, xi)
        call allocate_if_unallocated(ngllx, wgll)

        call get_gll(ngllx-1, xi, wgll)

        if (verbose.ge.5)then 
            write(*,'(/,a)')'• Setup GLL points'
            write(*,'(/,a, i1)')'  --> ngll: ', ngllx
            write(*,'(/,a)')'  -->  gll: '
            write(*,*) xi
            write(*,'(/,a)')'  -->  wgll: '
            write(*,*)wgll
        endif 

        ! Get derivative of lagrange polynomials
        call lagrange1st(ngllx-1)
        
    end subroutine setup_gll


    subroutine lagrange1st(N)
    ! Calculation of 1st derivatives of Lagrange polynomials
    ! at GLL collocation points
    ! dgll = legendre1st(N)
    ! dgll is a matrix with columns -> GLL nodes
    !                        rows -> order
        use allocation_module, only: allocate_if_unallocated
        use params, only: dgll,  xi
        implicit none 

        integer :: N

        integer :: i, j, m
        real(kind=CUSTOM_REAL) :: sum, d(N+1,N+1)

        call allocate_if_unallocated(N+1, N+1, dgll)

        d = zero 

        do i = 0, N 
            do j = 0, N
                if (i .ne. j) then
                    d(i+1, j+1) = legendre(N, xi(i + 1))/legendre(N, xi(j + 1)) * ONE/(xi(i + 1) - xi(j + 1))
                endif 

                if (i .eq. 0) then
                    if (j .eq. 0) then
                        d(i+1, j+1) = - N * (N + 1) / FOUR
                    endif
                endif

                if (i .eq. N) then
                    if (j .eq. N) then
                        d(i+1, j+1) =  N * (N + 1) / FOUR
                    endif
                endif
            enddo
        enddo 

        dgll = zero
        ! Calculate matrix with 1st derivatives of Lagrange polynomials
        do m = 0, N
            do i = 0, N 
                sum = zero
                do j = 0, N
                    sum = sum + d(i+1, j+1) * lagrange(N, m-1, xi(j+1))
                enddo 
                dgll(m+1,i+1) = sum
            enddo 
        enddo 

    end subroutine lagrange1st



    subroutine get_gll(N, xi, weights)
        ! Returns GLL (Gauss Lobato Legendre module with collocation points and
        ! weights)
        ! Initialization of integration weights and collocation points
        ! (/xi, weights/), kind=CUSTOM_REAL) =  gll(N)
        ! Values taken from Diploma Thesis Bernhard Schuberth
        ! Code from Lion Krischner's github
        ! https://github.com/krischer/2014_AdvancedSeismologySeminar/blob/master/MESS/DiscontinuousGalerkin/Python

        implicit none 

        ! IO variables: 
        integer :: N 
        real(kind=CUSTOM_REAL) :: xi(N+1), weights(N+1)

        if (N .eq. 2) then
            xi = real((/-1.0, 0.0, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.33333333, 1.33333333, 0.33333333/), kind=CUSTOM_REAL)
        elseif (N .eq. 3) then
            xi = real((/-1.0, -0.447213595499957, 0.447213595499957, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.1666666667, 0.833333333, 0.833333333, 0.1666666666/), kind=CUSTOM_REAL)
        elseif (N .eq. 4) then
            xi = real((/-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.1, 0.544444444, 0.711111111, 0.544444444, 0.1/), kind=CUSTOM_REAL)
        elseif (N .eq. 5) then
            xi = real((/-1.0, -0.7650553239294647, -0.285231516480645, 0.285231516480645,&
                0.7650553239294647, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0666666666666667,  0.3784749562978470,&
                    0.5548583770354862, 0.5548583770354862, 0.3784749562978470,&
                    0.0666666666666667/), kind=CUSTOM_REAL)
        elseif (N .eq. 6) then
            xi = real((/-1.0, -0.8302238962785670, -0.4688487934707142, 0.0,&
                0.4688487934707142, 0.8302238962785670, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0476190476190476, 0.2768260473615659, 0.4317453812098627,&
                    0.4876190476190476, 0.4317453812098627, 0.2768260473615659,&
                    0.0476190476190476/), kind=CUSTOM_REAL)
        elseif (N .eq. 7) then
            xi = real((/-1.0, -0.8717401485096066, -0.5917001814331423, &
                -0.2092992179024789, 0.2092992179024789, 0.5917001814331423, &
                0.8717401485096066, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0357142857142857, 0.2107042271435061, 0.3411226924835044, &
                    0.4124587946587038, 0.4124587946587038, 0.3411226924835044,&
                    0.2107042271435061, 0.0357142857142857/), kind=CUSTOM_REAL)
        elseif (N .eq. 8) then
            xi = real((/-1.0, -0.8997579954114602, -0.6771862795107377, &
                -0.3631174638261782, 0.0, 0.3631174638261782, &
                0.6771862795107377, 0.8997579954114602, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0277777777777778, 0.1654953615608055, 0.2745387125001617, &
                    0.3464285109730463, 0.3715192743764172, 0.3464285109730463,  &
                    0.2745387125001617, 0.1654953615608055, 0.0277777777777778/), kind=CUSTOM_REAL)
        elseif (N .eq. 9) then
            xi = real((/-1.0, -0.9195339081664589, -0.7387738651055050, &
                -0.4779249498104445, -0.1652789576663870, 0.1652789576663870, &
                0.4779249498104445, 0.7387738651055050, 0.9195339081664589, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0222222222222222, 0.1333059908510701, 0.2248893420631264, &
                    0.2920426836796838, 0.3275397611838976, 0.3275397611838976, &
                    0.2920426836796838, 0.2248893420631264, 0.1333059908510701, &
                    0.0222222222222222/), kind=CUSTOM_REAL)
        elseif (N .eq. 10) then
            xi = real((/-1.0, -0.9340014304080592, -0.7844834736631444, &
                -0.5652353269962050, -0.2957581355869394, 0.0, & 
                0.2957581355869394, 0.5652353269962050, 0.7844834736631444, &
                0.9340014304080592, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0181818181818182, 0.1096122732669949, 0.1871698817803052,&
                    0.2480481042640284, 0.2868791247790080, 0.3002175954556907,&
                    0.2868791247790080, 0.2480481042640284, 0.1871698817803052,&
                    0.1096122732669949, 0.0181818181818182/), kind=CUSTOM_REAL)
        elseif (N .eq. 11) then
            xi = real((/-1.0, -0.9448992722228822, -0.8192793216440067,&
                -0.6328761530318606, -0.3995309409653489, -0.1365529328549276,&
                0.1365529328549276, 0.3995309409653489, 0.6328761530318606,&
                0.8192793216440067, 0.9448992722228822, 1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0151515151515152, 0.0916845174131962, 0.1579747055643701,&
                    0.2125084177610211, 0.2512756031992013, 0.2714052409106962,&
                    0.2714052409106962, 0.2512756031992013, 0.2125084177610211,&
                    0.1579747055643701, 0.0916845174131962, 0.0151515151515152/), kind=CUSTOM_REAL)
        elseif (N .eq. 12) then
            xi = real((/-1.0, -0.9533098466421639, -0.8463475646518723,&
                -0.6861884690817575, -0.4829098210913362, -0.2492869301062400,&
                0.0, 0.2492869301062400, 0.4829098210913362,&
                0.6861884690817575, 0.8463475646518723, 0.9533098466421639,&
                1.0/), kind=CUSTOM_REAL)
            weights = real((/0.0128205128205128, 0.0778016867468189, 0.1349819266896083,&
                    0.1836468652035501, 0.2207677935661101, 0.2440157903066763,&
                    0.2519308493334467, 0.2440157903066763, 0.2207677935661101,&
                    0.1836468652035501, 0.1349819266896083, 0.0778016867468189,&
                    0.0128205128205128/), kind=CUSTOM_REAL)
        else
            write(*,*)'GLL: value N not implemented'
            stop
        endif 

    end subroutine get_gll 


    real(kind=CUSTOM_REAL) function legendre(N, x)
        ! Returns the value of Legendre Polynomial P_N(x) at position x[-1, 1].
        implicit none 

        ! IO variables: 
        real(kind=CUSTOM_REAL) :: x
        integer :: N

        ! Local: 
        integer :: i
        real(kind=CUSTOM_REAL) :: P(2 * N), ii 

        P(:) = zero
        if (N .eq. 0)then
            P(1) = one
        elseif (N .eq. 1) then
            P(2) = x
        else
            P(1) = 1
            P(2) = x
        endif 
        do i = 2, N 
            ii = real(i, kind=CUSTOM_REAL)
            P(i+1) = (one/ii)  * ((2*ii - 1) * x * P(i) - (ii - 1)*P(i - 1))
        enddo 
        legendre = P(N+1)

    end function legendre
    ! ------------------------------------------------------------------

end module gll