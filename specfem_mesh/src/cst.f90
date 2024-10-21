module splitting_function

    implicit none 
    include "constants.h"
    
    contains 



    real(kind=CUSTOM_REAL) function xi_kkst(md, m, s, t, ld, l)
        ! I made this one up. Its the same as Gilbert and Masters 2000
        ! Eqn 5 but without the second Wj3 symbol. If you look at 
        ! the woodhouse kernels in DT98 they have things written using 
        ! just the first wj3 so its easier to have this also 
        ! E.g. Apx eqn D.43 
        use w3j, only: thrj
        implicit none
        integer :: m, md, s, t, l, ld, j
        real(kind=CUSTOM_REAL) :: f_md, f_s, f_l, f_ld

        ! Convert integers to floats: 
        f_md = real(md, kind=CUSTOM_REAL)
        f_s  = real(s,  kind=CUSTOM_REAL)
        f_ld = real(ld, kind=CUSTOM_REAL)
        f_l  = real(l,  kind=CUSTOM_REAL)

        xi_kkst = (-ONE)**f_md                      * & 
                    ((TWO*f_ld+ONE)*(TWO*f_s+ONE)   * & 
                    (TWO*f_l+ONE)/(FOUR*PI))**HALF  * &
                    thrj(ld, s, l, -md, t, m)       
    end function



    real(kind=CUSTOM_REAL) function gamma_kkst_alt(md, m, s, t, ld, l, j)
        ! Gilbert & Woodhouse 2000 
        ! Determination of structure coefcients from splitting matrices
        ! Eqn 5
        ! NOTE THE SECOND w3j compared to the function above 
        use w3j, only: thrj
        implicit none
        integer :: m, md, s, t, l, ld, j

        gamma_kkst_alt = xi_kkst(md, m, s, t, ld, l) * &
                         thrj(ld+j, s+j, l+j, 0, 0, 0) 
    end function



    real(kind=CUSTOM_REAL) function gamma_kkst(md, m, s, t, ld, l, j)
        ! Gilbert & Woodhouse 2000 
        ! Determination of structure coefcients from splitting matrices
        ! Eqn 5
        ! NOTE THE SECOND w3j compared to the function above 
        use w3j, only: thrj
        implicit none
        integer :: m, md, s, t, l, ld, j
        real(kind=CUSTOM_REAL) :: f_md, f_s, f_l, f_ld

        ! Convert integers to floats: 
        f_md = real(md, kind=CUSTOM_REAL)
        f_s  = real(s,  kind=CUSTOM_REAL)
        f_ld = real(ld, kind=CUSTOM_REAL)
        f_l  = real(l,  kind=CUSTOM_REAL)

        gamma_kkst = (-ONE)**f_md                     * & 
                     ((TWO*f_ld+ONE)*(TWO*f_s+ONE)    * & 
                      (TWO*f_l+ONE)/(FOUR*PI))**HALF  * &
                     thrj(ld, s, l, -md, t, m)        * &
                     thrj(ld+j, s+j, l+j, 0, 0, 0) 
    end function


    real(kind=CUSTOM_REAL) function F_mst(m, s, t, l, ld, j)
        use w3j, only: thrj
        implicit none
        integer :: m, s, t, l, ld, j
        real(kind=CUSTOM_REAL) :: f_t, f_s, f_ld, f_l, f_m
 
        f_t  = real(t,  kind=CUSTOM_REAL)
        f_s  = real(s,  kind=CUSTOM_REAL)
        f_ld = real(ld, kind=CUSTOM_REAL)
        f_l  = real(l,  kind=CUSTOM_REAL)
        f_m  = real(m,  kind=CUSTOM_REAL)

        F_mst = (-ONE)**(f_m + f_t) * (FOUR*PI*(TWO*f_s + ONE)  / & 
                       ((TWO*f_ld + ONE)*(TWO*f_l +ONE)))**HALF * & 
                        thrj(ld, s, l, -m-t, t, m)/thrj(ld+j, s+j, l+j, 0, 0, 0)
    end function F_mst




    subroutine get_Ssum_bounds(l1, l2, smin, smax, num_s, max_num_t)
        implicit none 
        integer :: l1, l2, smin, smax, num_s, max_num_t

        ! Compute bounds of summation over s and t: 
        smin = abs(l1-l2)           ! minimum s value
        smax = l1+l2                ! maximum s value
        num_s = smax - smin + 1     ! number of s to sum
        max_num_t = 2*smax+1        ! -t to t for smax
    end subroutine



    subroutine cst_to_H(H, ld, l, cst, ncols, nrows, type1, type2)
        ! Compute elements of the splitting matrix H from Csts
        ! Gilbert & Woodhouse 2000 Eqn. 4
        ! m, md index the element of the matrix to be computed
        ! l and ld indicate the order of the modes involved
        implicit none
        integer :: l, ld, ncols, nrows
        real(kind=SPLINE_REAL)   :: cst(nrows, ncols)
        real(kind=SPLINE_REAL):: H(2*ld+1, 2*l+1)
        character :: type1, type2
        ! Local: 
        integer :: smin, smax, num_s, max_num_t, md, m, t, s,j,is,it
        real(kind=CUSTOM_REAL) :: sum 

        ! j = 0 for S-S or T-T coupling and 1 for mixed
        if(type1.eq.type2)then
            j = 0
        else
            j = 1
        endif

        ! Get the summation bounds
        ! Note here that ncols and nrows should be the same as num_s
        ! max_num_t
        call get_Ssum_bounds(ld, l, smin, smax, num_s, max_num_t)
        if(nrows.ne.num_s .or. ncols.ne.max_num_t)then
            write(*,*)'Error: discrepency in nrow/num_s or ncol/max_num_t: '
            write(*,*)'nrows      = ', nrows
            write(*,*)'num_s      = ', num_s
            write(*,*)'ncols      = ', ncols
            write(*,*)'max_num_t  = ', max_num_t
            stop
        endif
        
        ! Compute for each element of the matrix: 
        do md = -ld, ld ! row
            do m = -l, l ! col
                sum = zero

                ! sum over the s
                do is = 1, num_s 
                    ! compute the actual S value
                    s = smin + is - 1
                    do it = 1, 2*s + 1
                        ! compute the actual t value: 
                        t = it - 1 - s
                        sum = sum + cst(is,it) * & 
                                    gamma_kkst(md, m, s, t, l, ld, j)
                    enddo ! it
                enddo !is 

                H(md+ld+1, m+l+1)  = sum
            enddo ! m
        enddo ! md

    end subroutine cst_to_H



    subroutine H_to_cst(H, ld, l, cst, ncols, nrows, type1, type2)
        ! Computes splitting coefficients from the matrix H
        implicit none 

        integer :: l, ld, ncols, nrows
        real(kind=SPLINE_REAL)   :: cst(nrows, ncols)
        real(kind=SPLINE_REAL):: H(2*ld+1, 2*l+1)
        character :: type1, type2


        ! Local: 
        integer :: smin, smax, num_s, max_num_t, is, s, it, t, m, & 
                   r0, rs, Nd, rt, ct, N, R, im, j, md
        
        

        if(ld.lt.l)then
            write(*,*)"ERROR: l' (ld) must be greater than l: "
            write(*,*)"ERROR: l' = ", ld
            write(*,*)"ERROR: l  = ", l  
            stop
        endif

            

        ! j = 0 for S-S or T-T coupling and 1 for mixed
        if(type1.eq.type2)then
            j = 0
        else
            j = 1
        endif
        
        cst = ZERO

        call get_Ssum_bounds(ld, l, smin, smax, num_s, max_num_t)
        if(nrows.ne.num_s .or. ncols.ne.max_num_t)then
            write(*,*)'Error: discrepency in nrow/num_s or ncol/max_num_t: '
            write(*,*)'nrows      = ', nrows
            write(*,*)'num_s      = ', num_s
            write(*,*)'ncols      = ', ncols
            write(*,*)'max_num_t  = ', max_num_t
            stop
        endif
        
        r0 = 1 + ld - l       ! starting row for t = 0
        Nd = 2*l +1           ! Maximum diagonal length 
  
        do is = 1, num_s 
            s = smin + is - 1
            do it = 1, 2*s + 1
                t = it - 1 - s
                ! For each s, t, we need to sum over the (sub)diagonal 
                ! defined by the t value 

                rs = r0 + t           ! Theoretical starting row
                if (rs.le.0)then 
                    rt = 1            ! actual starting row
                    ct = 2 - rs       ! actual starting column 
                    N  = Nd + rs - 1  ! number of elements in diagonal 
                else
                    rt = rs 
                    ct = 1 
                    R  = 2*(ld-l) + 1 - rs
                    if(R.lt.0)then 
                        N = Nd + R
                    else 
                        N = Nd 
                    endif 
                endif

                do im = 1, N
                    md = -(ld) + (rt+im-1) -1 
                    m  = -(l)  + (ct+im-1) -1
                    cst(is,it) = cst(is,it) + F_mst(m, s, t, l, ld, j) * & 
                                              H(md+ld+1, m+l+1)
                enddo 

            enddo 
        enddo 



    end subroutine H_to_cst



    subroutine write_cst_to_file(fname, cst, ncols, nrows, smin, jump)
        ! Jump of 2 for only even stuff, 1 for all 
        implicit none 
        character(len=*)         :: fname
        real(kind=SPLINE_REAL)   :: cst(nrows, ncols)
        integer                  :: nrows, ncols, smin, jump, is ,it
        
        write(*,*)'Writing to '//trim(fname)

        open(1,file=trim(fname), form='formatted')
        do is = 1, nrows, jump
            do it = 1, 2*(smin + is-1)+1
                write(1,*)smin+is-1, it - (smin+is-1) - 1, real(cst(is,it))
            enddo 
        enddo 

    end subroutine write_cst_to_file






end module splitting_function