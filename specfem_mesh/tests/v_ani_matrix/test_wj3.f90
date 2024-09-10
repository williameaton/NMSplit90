program test_wj3 
    use Vani, only: thrj
    implicit none 
    include "constants.h"

    real(kind=CUSTOM_REAL) :: out

    out = thrj(22, 2, 20, 0, 0, 0)

    write(*,*)out
end program