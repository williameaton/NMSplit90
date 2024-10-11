module w3j
    include "constants.h"

    contains 
    real(CUSTOM_REAL) function thrj(j1, j2, j3, m1, m2, m3)
    ! Arguments
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    integer :: jc, ja, jb, ma, mb, mc, lm1, lm2, m, n, mx, my
    real(CUSTOM_REAL) :: ss, alpha, beta, gamma, y(100)

    !-----
    ! Evaluation of Wigner 3-j coefficients
    thrj = zero

    ! Early returns based on selection criteria
    if ((j1 + j2 - j3 < 0) .or. (j2 + j3 - j1 < 0) .or. (j3 + j1 - j2 < 0)) return
    if (j1 - abs(m1) < 0) return
    if (j2 - abs(m2) < 0) return
    if (j3 - abs(m3) < 0) return
    if (m1 + m2 + m3 /= 0) return

    !-----
    ! Use symmetries to make j3 the largest of j1, j2, j3
    jc = max(j1, j2, j3)

    if(jc.eq.j3)then
            ja = j1; jb = j2; jc = j3
            ma = m1; mb = m2; mc = m3
    elseif (jc.eq.j2)then
            ja = j3; jb = j1; jc = j2
            ma = m3; mb = m1; mc = m2
    else 
            ja = j2; jb = j3; jc = j1
            ma = m2; mb = m3; mc = m1
    endif

    lm2 = -jb
    if (ja + mc - jb < 0) lm2 = -mc - ja

    lm1 = -mc - lm2
    m = lm1 + jb + mc + 1

    if (ja - jb - mc < 0) m = lm1 + ja + 1

    ! Initialize array y
    y(1) = zero
    y(2) = one
    ss   = one

    ! Main calculation loop
    if (m > 1) then
        do n = 2, m
            mx = lm1 - n + 2
            my = lm2 + n - 2
            alpha = sqrt(real((ja - mx + 1) * (ja + mx) * (jb + my + 1) * (jb - my), CUSTOM_REAL))
            beta = real(jc * (jc + 1) - ja * (ja + 1) - jb * (jb + 1) - 2 * mx * my, CUSTOM_REAL)
            gamma = sqrt(real((ja + mx + 1) * (ja - mx) * (jb - my + 1) * (jb + my), CUSTOM_REAL))
            y(n+1) = (beta * y(n) - gamma * y(n-1)) / alpha
            ss = ss + y(n+1)**two
        end do
    end if

    ! Compute final result
    n = lm1 - ma + 2
    thrj = y(n) * ((-one)**real(-ja + jb + mc,kind=CUSTOM_REAL)) / sqrt(ss * real(2 * jc + 1, CUSTOM_REAL))

    return
    end function thrj

end module 
