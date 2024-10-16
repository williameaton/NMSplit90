module cuda_proxies
        include "constants.h"
        contains 
        subroutine compute_vani_sc_cuda(l1, n, p, q)
                integer :: l1, n, p
                real(kind=CUSTOM_REAL)    :: q(n,n,n,p)
        end subroutine compute_vani_sc_cuda

end module cuda_proxies