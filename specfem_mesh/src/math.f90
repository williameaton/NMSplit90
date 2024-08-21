module math 
    ! Used for playing around with the precision of the code

    implicit none 


    interface cosp
        module procedure cosine_real
        module procedure cosine_double
    end interface cosp
    
    interface sinp
        module procedure sin_real
        module procedure sin_double
    end interface sinp
     
    interface tanp
        module procedure tan_real
        module procedure tan_double
    end interface tanp   


    interface acosp
        module procedure inv_cosine_real
        module procedure inv_cosine_double
    end interface acosp
    
    interface asinp
        module procedure inv_sin_real
        module procedure inv_sin_double
    end interface asinp
     
    interface atan2p
        module procedure inv_tan2_real
        module procedure inv_tan2_double
    end interface atan2p  

    interface sqrtp
        module procedure sqrt_real
        module procedure sqrt_double
    end interface sqrtp
    
    contains 


    ! Square root 
    real(kind=4) function sqrt_real(x)
        implicit none 
        real(kind=4) :: x
        sqrt_real = sqrt(x)
    end function sqrt_real

    double precision function sqrt_double(x)
        implicit none 
        real(kind=8) :: x
        sqrt_double = dsqrt(x)
    end function sqrt_double




    ! Cos 
    real(kind=4) function cosine_real(x)
        implicit none 
        real(kind=4) :: x
        cosine_real = cos(x)
    end function cosine_real

    double precision function cosine_double(x)
        implicit none 
        real(kind=8) :: x
        cosine_double = dcos(x)
    end function cosine_double


    ! Sin
    real(kind=4) function sin_real(x)
        implicit none 
        real(kind=4) :: x
        sin_real = sin(x)
    end function sin_real

    double precision function sin_double(x)
        implicit none 
        real(kind=8) :: x
        sin_double = dsin(x)
    end function sin_double


    ! Tan
    real(kind=4) function tan_real(x)
        implicit none 
        real(kind=4) :: x
        tan_real = tan(x)
    end function tan_real

    double precision function tan_double(x)
        implicit none 
        real(kind=8) :: x
        tan_double = dtan(x)
    end function tan_double


    ! INVERSE TRIGONOMETRIC FUNCTIONS 


    ! arccos 
    real(kind=4) function inv_cosine_real(x)
        implicit none 
        real(kind=4) :: x
        inv_cosine_real = acos(x)
    end function inv_cosine_real

    double precision function inv_cosine_double(x)
        implicit none 
        real(kind=8) :: x
        inv_cosine_double = dcos(x)
    end function inv_cosine_double


    ! arcsin
    real(kind=4) function inv_sin_real(x)
        implicit none 
        real(kind=4) :: x
        inv_sin_real = asin(x)
    end function inv_sin_real

    double precision function inv_sin_double(x)
        implicit none 
        real(kind=8) :: x
        inv_sin_double = dasin(x)
    end function inv_sin_double


    ! arctan
    real(kind=4) function inv_tan2_real(y,x)
        implicit none 
        real(kind=4) :: y,x
        inv_tan2_real = atan2(y,x)
    end function inv_tan2_real

    double precision function inv_tan2_double(y,x)
        implicit none 
        real(kind=8) :: y,x
        inv_tan2_double = datan2(y,x)
    end function inv_tan2_double


end module math