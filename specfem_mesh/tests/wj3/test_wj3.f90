program test_wj3 
    use w3j, only: thrj
    implicit none 
    include "constants.h"

    real(kind=CUSTOM_REAL) :: out

    open(unit=1,file='w3j.txt', &
        status='unknown',form='formatted',action='write')

    out = thrj(22, 2, 20, 0, 0, 0)
    write(1,*)out


end program