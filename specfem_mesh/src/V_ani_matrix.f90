

module V_ani 
    include "constants.h"
contains 

    subroutine compute_bond_matrix(n1, n2, M)
        ! Computes bond matrix for a rotation defined by angles n1 and n2 
        ! M should be a 6 x 6 matrix 
        ! M is built for a rotation matrix A which is defined by 
        !       | c1c2   -s1   s2c1 |   c1 = cos n1   ;   s1 = sin n1
        ! A =   | s1c2    c1   s1s2 |   c2 = cos n2   ;   s2 = sin n2
        !       | -s2     0     c2  |
        ! See Auld 1973. Eqn 3.34 
        use math, only: cosp, sinp 
        implicit none   
        include "constants.h"

        real(kind=CUSTOM_REAL) :: n1, n2 
        real(kind=CUSTOM_REAL) :: M(6, 6)

        real(kind=CUSTOM_REAL) :: A(3, 3), c1, c2, s1, s2, c11, c22, s11, s22
        integer :: i, j 

        c1 = cosp(n1)
        c2 = cosp(n2)
        s1 = sinp(n1)
        s2 = sinp(n2)

        c11 = c1*c1
        c22 = c2*c2
        s11 = s1*s1
        s22 = s2*s2

        ! ------------ Build rotation matrix:  ------------
        A(1,1) = c1 * c2 
        A(2,1) = s1 * c2 
        A(3,1) = -s2

        A(1,2) = -s1
        A(2,2) = c1
        A(3,2) = zero

        A(1,3) = s2 * c1
        A(2,3) = s1 * s2
        A(3,3) = c2

        ! -------------- Build bond matrix:  --------------
        ! Upper left submatrix
        do i = 1,3 
            do j = 1,3
                M(i,j) = A(i,j)**2
            enddo
        enddo 

        ! Lower left submatrix
        do i = 1, 3
            M(4,i) = A(2,i) * A(3,i)
            M(5,i) = A(3,i) * A(1,i)
            M(6,i) = A(1,i) * A(2,i)
        enddo 

        ! Upper right submatrix 
        do i = 1, 3
            M(i,4) = two * A(i,2) * A(i,3)
            M(i,5) = two * A(i,3) * A(i,1)
            M(i,6) = two * A(i,1) * A(i,2)
        enddo 

        ! Lower right submatrix
        M(4,4) = A(2,2) * A(3,3)   +   A(2,3) * A(3,2)
        M(5,4) = A(1,2) * A(3,3)   +   A(1,3) * A(3,2)
        M(6,4) = A(1,2) * A(2,3)   +   A(1,3) * A(2,2)

        M(4,5) = A(2,1) * A(3,3)   +   A(2,3) * A(3,1)
        M(5,5) = A(1,3) * A(3,1)   +   A(1,1) * A(3,3)
        M(6,5) = A(1,3) * A(2,1)   +   A(1,1) * A(2,3)

        M(4,6) = A(2,2) * A(3,1)   +   A(2,1) * A(3,2)
        M(5,6) = A(1,1) * A(3,2)   +   A(1,2) * A(3,1)
        M(6,6) = A(1,1) * A(2,2)   +   A(1,2) * A(2,1)
        return 
    end subroutine compute_bond_matrix




    subroutine compute_bond_matrix_explicit(n1, n2, Q)
        ! Computes bond matrix for a rotation defined by angles n1 and n2 
        ! M should be a 6 x 6 matrix 
        ! M is built for a rotation matrix A which is defined by 
        !       | c1c2   -s1   s2c1 |   c1 = cos n1   ;   s1 = sin n1
        ! A =   | s1c2    c1   s1s2 |   c2 = cos n2   ;   s2 = sin n2
        !       | -s2     0     c2  |
        ! See Auld 1973. Eqn 3.34 
        use math, only: cosp, sinp 
        implicit none   
        include "constants.h"

        real(kind=CUSTOM_REAL) :: n1, n2 
        real(kind=CUSTOM_REAL) :: Q(6, 6)

        real(kind=CUSTOM_REAL) :: c1, c2, s1, s2, c11, c22, s11, s22

        c1 = cosp(n1)
        c2 = cosp(n2)
        s1 = sinp(n1)
        s2 = sinp(n2)

        c11 = c1*c1
        c22 = c2*c2
        s11 = s1*s1
        s22 = s2*s2
        ! -------------- Manual bond matrix:  --------------
        Q(1,1) = c11 * c22 
        Q(2,1) = s11 * c22 
        Q(3,1) = s22
        Q(4,1) = -s1 * s2 * c2
        Q(5,1) = -c1 * c2 * s2
        Q(6,1) = c1 * c22 * s1 

        Q(1,2) = s11
        Q(2,2) = c11 
        Q(3,2) = zero
        Q(4,2) = zero
        Q(5,2) = zero
        Q(6,2) = -s1 * c1

        Q(1,3) = s22 * c11
        Q(2,3) = s11 * s22
        Q(3,3) = c22 
        Q(4,3) = s1*s2*c2
        Q(5,3) = c1*c2*s2
        Q(6,3) = s1*s22*c1

        Q(1,4) = -two * s1 * s2 * c1
        Q(2,4) = two * s1 * s2 * c1 
        Q(3,4) = zero 
        Q(4,4) = c1*c2 
        Q(5,4) = -s1*c2
        Q(6,4) = -s2*(c11 + s11*c2)

        Q(1,5) = two * c11 * s2 * c2
        Q(2,5) = two * s11 * c2 * s2
        Q(3,5) = - two * s2 * c2
        Q(4,5) = s1 * (c22 - s22)
        Q(5,5) = c1 * (c22 - s22)
        Q(6,5) = two * c1 * c2 * s1 * s2

        Q(1,6) = - two * c1 * c2 * s1
        Q(2,6) =   two * s1 * c1 * c2
        Q(3,6) =  zero 
        Q(4,6) = - c1 * s2
        Q(5,6) = s1 * s2
        Q(6,6) = c2*(c11 - s11)
        return 
    end subroutine compute_bond_matrix_explicit


    subroutine print_bond_matrix(M)
        implicit none 
        include "constants.h"
        real(kind=CUSTOM_REAL) :: M(6,6)

        integer :: i 

        do i = 1,6 
            write(*,*)M(i,:)
        enddo 

    end subroutine print_bond_matrix

    subroutine setup_Cnatural(Mat, A, C, L, N, F)
        ! Sets VTI matrix given A, C, L, N, F
        implicit none 
        real(kind=CUSTOM_REAL) :: Mat(6,6)
        real(kind=CUSTOM_REAL) :: A, C, L, N, F

        Mat(:,:) = zero
        Mat(1,1) = A 
        Mat(2,2) = A 
        Mat(3,3) = C 
        Mat(4,4) = L 
        Mat(5,5) = L 
        Mat(6,6) = N 

        Mat(1,2) = A - two*N
        Mat(2,1) = A - two*N

        Mat(1,3) = F
        Mat(3,1) = F
        Mat(2,3) = F
        Mat(3,2) = F
    end subroutine setup_Cnatural



    subroutine compute_Cxyz_at_gll_radialACLNF(n1, n2, radmap)
        ! Nlen is the length of the array of unique A,C,L,N,F values
        ! probably = unique_r
        ! radmap is the map from the GLL to the unique radius
        use params, only: ngllx, nglly, ngllz, nspec, Cxyz, Arad, Crad, Lrad, Nrad, Frad
        use allocation_module, only: deallocate_if_allocated
        implicit none 
        include "constants.h"

        real(kind=CUSTOM_REAL) :: n1, n2
        real(kind=CUSTOM_REAL) :: M(6,6), Cnat(6,6)
        integer :: i, j, k, ispec, radmap(ngllx, nglly, ngllz, nspec), r_id

        call deallocate_if_allocated(Cxyz)
        allocate(Cxyz(6, 6, ngllx, nglly, ngllz, nspec))


        ! Compute rotation matrix per GLL in xyz
        do ispec = 1, nspec 
            do i = 1, ngllx
                do j = 1, nglly 
                    do k = 1, ngllz 

                        ! Get C matrix in natural orientation for this gll
                        r_id = radmap(i,j,k,ispec)
                        call setup_Cnatural(Cnat, Arad(r_id), Crad(r_id), Lrad(r_id), &
                                             Nrad(r_id), Frad(r_id))

                        call compute_bond_matrix_explicit(n1, n2, M)

                        Cxyz(:,:,i,j,k,ispec) = real(matmul(matmul(M, Cnat), transpose(M)), kind=SPLINE_REAL)

                    enddo
                enddo
            enddo   
        enddo


    end subroutine compute_Cxyz_at_gll_radialACLNF



    subroutine compute_Cxyz_at_gll_constantACLNF(A, C, L, N, F, n1, n2)
        use params, only: ngllx, nglly, ngllz, nspec, Cxyz
        use allocation_module, only: deallocate_if_allocated
        implicit none 
        include "constants.h"

        real(kind=CUSTOM_REAL) :: A, C, L, N, F, n1, n2
        real(kind=CUSTOM_REAL) :: M(6,6), Cnat(6,6)
        integer :: i, j, k, ispec 


        call deallocate_if_allocated(Cxyz)
        allocate(Cxyz(6, 6, ngllx, nglly, ngllz, nspec))

        ! Transversly isotropic matrix 
        call setup_Cnatural(Cnat, A, C, L, N, F)

        ! Compute rotation matrix per GLL in xyz
        do ispec = 1, nspec 
            do i = 1, ngllx
                do j = 1, nglly 
                    do k = 1, ngllz 
                        call compute_bond_matrix_explicit(n1, n2, M)

                        Cxyz(:,:,i,j,k,ispec) = real(matmul(matmul(M, Cnat), transpose(M)), kind=SPLINE_REAL)

                    enddo
                enddo
            enddo   
        enddo

    end subroutine compute_Cxyz_at_gll_constantACLNF





    subroutine compute_Vani_matrix(t1, l1, n1, t2, l2, n2, store, iproc)
        ! Computes V^{ani}_{kk'} between l and l'
        ! This assumes that the Cxyz matrix has already been computed for each GLL point 
        use params, only: NL, IC_ID, n_unique_rad, unique_r, interp_id_r, strain1, ngllx, & 
                        nglly, ngllz, nspec, strain2, ngllx, nglly, ngllz, & 
                        nspec, Cxyz, Vani, wglljac
        use allocation_module, only: deallocate_if_allocated
        use spline, only: interpolate_mode_eigenfunctions
        use mesh_utils, only: rotate_complex_sym_matrix_rtp_to_xyz
        implicit none 
        include "constants.h"

        ! IO variables
        integer   :: l1, n1, l2, n2
        character :: t1, t2
        logical   :: store
        integer, optional :: iproc

        ! Local variables 
        logical :: self_coupling
        integer :: m1, m2
        real(SPLINE_REAL) :: wcom, qmod
        real(SPLINE_REAL):: u1(NL), v1(NL), du1(NL), dv1(NL)
        real(SPLINE_REAL):: u2(NL), v2(NL), du2(NL), dv2(NL)

        real(SPLINE_REAL)::  u_spl_1(n_unique_rad),  v_spl_1(n_unique_rad),& 
                            du_spl_1(n_unique_rad), dv_spl_1(n_unique_rad)

        real(SPLINE_REAL), allocatable :: u_spl_2(:), v_spl_2(:), du_spl_2(:), dv_spl_2(:)

        integer :: i, j, k, ispec, p, q
        complex(SPLINE_REAL) :: sum, cont

        integer, dimension(9), parameter :: Vcont = (/1, 2, 3, 4, 4, 5, 5, 6, 6/)

        ! Check if self_coupling because it determines if we need to get
        ! strain for two modes or just one
        if (t1.eq.t2 .and. l1.eq.l2 .and. n1.eq.n2)then 
            self_coupling = .true.
        else
            self_coupling = .false.
        endif 

        ! Load modes and interpolates eigenfunctions:  
        call get_mode(t1, n1, l1, wcom, qmod, u1, du1, v1, dv1, .false.)
        call interpolate_mode_eigenfunctions(t1, u1, v1, du1, dv1, 1, IC_ID,      &  
                                            unique_r, n_unique_rad, interp_id_r, & 
                                            u_spl_1, v_spl_1, du_spl_1, dv_spl_1)

        if(.not.self_coupling)then 
            call get_mode(t2, n2, l2, wcom, qmod, u2, du2, v2, dv2, .false.)
            allocate( u_spl_2(n_unique_rad),  v_spl_2(n_unique_rad))
            allocate(du_spl_2(n_unique_rad), dv_spl_2(n_unique_rad))
            call interpolate_mode_eigenfunctions(t2, u2, v2, du2, dv2, 1, IC_ID,      &  
                                                unique_r, n_unique_rad, interp_id_r, & 
                                                u_spl_2, v_spl_2, du_spl_2, dv_spl_2)
            ! We will need this strain also 
            call deallocate_if_allocated(strain2)
            allocate(strain2(6, ngllx, nglly, ngllz, nspec))
        endif
        
        ! Allocate the strain for this proc. 
        call deallocate_if_allocated(strain1)
        allocate(strain1(6, ngllx, nglly, ngllz, nspec))

        do m1 = -l1, l1 
            ! Compute strain and store if desired
            call compute_gll_mode_strain(t1, n1, l1, m1, u_spl_1, v_spl_1, du_spl_1, & 
                                         dv_spl_1, n_unique_rad, strain1)
            call rotate_complex_sym_matrix_rtp_to_xyz(strain1)
            if(store)then
                call save_mode_strain_binary(n1, t1, l1, m1, strain1, iproc)
            endif

            if(self_coupling) then 
                ! Only need 1 strain 
                sum = SPLINE_iZERO
                do ispec = 1, nspec 
                    do i = 1, ngllx
                        do j = 1, nglly
                            do k = 1, ngllz
                                cont = SPLINE_iZERO
                                do p = 1, 9
                                    do q = 1, 9
                                        cont = cont + ( conjg(strain1(Vcont(p),i,j,k,ispec)) * & 
                                                        Cxyz(Vcont(p),Vcont(q),i,j,k,ispec)  * & 
                                                        strain1(Vcont(q),i,j,k,ispec) ) 
                                    enddo
                                enddo 
                                sum = sum  +  wglljac(i,j,k,ispec) * cont 
                            enddo 
                        enddo
                    enddo
                enddo
                Vani(m1+l1+1, m1+l2+1) = Vani(m1+l1+1, m1+l2+1) + sum

            else 
                ! Not self coupling so need 2 strains 
                do m2 = -l2, l2
                    call compute_gll_mode_strain(t2, n2, l2, m2, u_spl_2, v_spl_2, & 
                                                 du_spl_2, dv_spl_2, n_unique_rad, strain2)
                    call rotate_complex_sym_matrix_rtp_to_xyz(strain2)
                    if(store)then
                        call save_mode_strain_binary(n2, t2, l2, m2, strain2, iproc)
                    endif

                    sum = SPLINE_iZERO
                    do ispec = 1, nspec 
                        do i = 1, ngllx
                            do j = 1, nglly
                                do k = 1, ngllz
                                    cont = SPLINE_iZERO
                                    do p = 1, 9
                                        do q = 1, 9
                                            cont = cont + ( conjg(strain1(Vcont(p),i,j,k,ispec)) * & 
                                                            Cxyz(Vcont(p),Vcont(q),i,j,k,ispec) * & 
                                                            strain2(Vcont(q),i,j,k,ispec) ) 
                                        enddo
                                    enddo 
                                    sum = sum  +  wglljac(i,j,k,ispec) * cont 
                                enddo 
                            enddo
                        enddo
                    enddo
                    Vani(m1+l1+1, m2+l2+1) = Vani(m1+l1+1, m2+l2+1) + sum

                    write(*,*)' done m2 ', m2
                enddo ! m2 
            endif 
            write(*,*)' done m1 ', m1

        enddo ! m1

    end subroutine compute_Vani_matrix




    subroutine compute_Vani_matrix_stored_selfcoupling(t1, l1, n1, iproc)
        ! A specific case of compute_Vani_matrix used for fast self coupling computation
        ! by loading stored strain values for modes
        use params, only: strains, ngllx, nglly, ngllz, nspec, wglljac, Cxyz, Vani
        use allocation_module, only: deallocate_if_allocated
        implicit none 
        include "constants.h"

        ! IO variables
        integer           :: l1, n1
        character         :: t1
        integer, optional :: iproc

        ! Local variables 
        integer :: m1, m2, im, i, j, k, ispec, p, q
        complex(SPLINE_REAL) :: sum, cont

        integer, dimension(9), parameter :: Vcont = (/1, 2, 3, 4, 4, 5, 5, 6, 6/)

               
        ! Allocate the strain for this proc. 
        call deallocate_if_allocated(strains)
        allocate(strains(6, ngllx, nglly, ngllz, nspec, 2*l1+1))

        ! Load all strains once and for all for each m value
        ! m is an integer from -l to +l
        ! In this case l is some integer between 1 and ~20
        do im =  -l1, l1 
            call load_mode_strain_binary(n1, t1, l1, im, strains(:,:,:,:,:,l1+im+1), iproc)
        enddo 


        ! Loop over row (m1) and column (m2) of matrix
        do m1 = -l1, l1 
            do m2 = -l1, l1 !-m1
                    sum = SPLINE_iZERO
                    do ispec = 1, nspec 
                        do i = 1, ngllx
                            do j = 1, nglly
                                do k = 1, ngllz

                                    ! Contraction of strain with elastic tensor
                                    ! E_ij c_ijkl E_k^*
                                    cont = SPLINE_iZERO
                                    do p = 1, 9
                                        do q = 1, 9
                                            cont = cont +                                          &
                                                   (conjg(strains(Vcont(p),i,j,k,ispec,m1+l1+1)) * & 
                                                    Cxyz(Vcont(p),Vcont(q),i,j,k,ispec)          * & 
                                                    strains(Vcont(q),i,j,k,ispec,m2+l1+1) ) 
                                        enddo ! q
                                    enddo ! p 
                    
                                    sum = sum  +  wglljac(i,j,k,ispec) * cont 
                                enddo ! k
                            enddo ! j
                        enddo ! i
                    enddo ! ispec

                    ! Add contribution to V_ani matrix
                    Vani(m1+l1+1, m2+l1+1) = Vani(m1+l1+1, m2+l1+1) + sum
            enddo !m2
        enddo !m1 

    end subroutine compute_Vani_matrix_stored_selfcoupling




    subroutine cuda_Vani_matrix_stored_selfcoupling(t1, l1, n1, iproc)
        !use omp_lib
        use params, only: strains, ngllx, nglly, ngllz, nspec, wglljac, Cxyz, Vani
        use allocation_module, only: deallocate_if_allocated
        use vani_kernel, only: compute_vani_sc_cuda
        implicit none 
        include "constants.h"

        ! IO variables
        integer           :: l1, n1
        character         :: t1
        integer, optional :: iproc

        ! Local variables 
        integer :: m1, m2, im, i, j, k, ispec, p, q
        complex(SPLINE_REAL) :: sum, cont

        integer, dimension(9), parameter :: Vcont = (/1, 2, 3, 4, 4, 5, 5, 6, 6/)

        ! Timing: 
        integer ::  loop_begin, loop_end, load_begin, load_end, count_rate
        real(8) :: elapsed_time
    
               
        ! Allocate the strain for this proc. 
        call deallocate_if_allocated(strains)
        allocate(strains(6, ngllx, nglly, ngllz, nspec, 2*l1+1))

        ! Load all strains once and for all for each m value

        do im =  -l1, l1 
            call load_mode_strain_binary(n1, t1, l1, im, strains(:,:,:,:,:,l1+im+1), iproc)
        enddo 

        call compute_vani_sc_cuda(l1)

 
    end subroutine cuda_Vani_matrix_stored_selfcoupling

















    subroutine compute_Vani_matrix_stored(t1, l1, n1, t2, l2, n2, iproc)
        ! A specific case of compute_Vani_matrix used for fast computation
        ! by loading stored strain values for modes
        use params, only: strain1, ngllx, nglly, ngllz, nspec, &   
                          strain2, wglljac, Cxyz, Vani
        use allocation_module, only: deallocate_if_allocated
        implicit none 
        include "constants.h"

        ! IO variables
        integer   :: l1, n1, l2, n2
        character :: t1, t2
        integer, optional :: iproc

        ! Local variables 
        integer :: m1, m2
        integer :: i, j, k, ispec, p, q
        complex(SPLINE_REAL) :: sum, cont

        integer, dimension(9), parameter :: Vcont = (/1, 2, 3, 4, 4, 5, 5, 6, 6/)

               
        ! Allocate the strain for this proc. 
        call deallocate_if_allocated(strain1)
        call deallocate_if_allocated(strain2)
        allocate(strain1(6, ngllx, nglly, ngllz, nspec))
        allocate(strain2(6, ngllx, nglly, ngllz, nspec))

        do m1 = -l1, l1 
            call load_mode_strain_binary(n1, t1, l1, m1, strain1, iproc)
            do m2 = -l2, l2
                call load_mode_strain_binary(n2, t2, l2, m2, strain2, iproc)
                    sum = SPLINE_iZERO
                    do ispec = 1, nspec 
                        do i = 1, ngllx
                            do j = 1, nglly
                                do k = 1, ngllz
                                    cont = SPLINE_iZERO
                                    do p = 1, 9
                                        do q = 1, 9
                                            cont = cont + ( conjg(strain1(Vcont(p),i,j,k,ispec)) * & 
                                                            Cxyz(Vcont(p),Vcont(q),i,j,k,ispec)  * & 
                                                            strain2(Vcont(q),i,j,k,ispec) ) 
                                        enddo
                                    enddo 
                                    sum = sum  +  wglljac(i,j,k,ispec) * cont 
                                enddo 
                            enddo
                        enddo
                    enddo
                    Vani(m1+l1+1, m2+l2+1) = Vani(m1+l1+1, m2+l2+1) + sum
            enddo !m2
        enddo !m1 

    end subroutine compute_Vani_matrix_stored









    
    real(kind=SPLINE_REAL) function OmNL(N, l)
        ! Computes C.146 of Dahlen and Tromp

        integer :: N, l
        real(kind=CUSTOM_REAL) :: lf, Nf, O

        lf = real(l, kind=CUSTOM_REAL)
        Nf = real(abs(N), kind=CUSTOM_REAL)

        if(N.ge.0)then 
            O = (half*(lf+Nf)*(lf-Nf+one))**half 
        else 
            O = (half*(lf-Nf)*(lf+Nf+one))**half
        endif 

        OmNL = real(O, kind=SPLINE_REAL)
    end function OmNL



    integer function ij_to_voigt(i, j)

        integer :: i, j, indx 

        indx = 10*i + j 

        select case(indx)
            case(11)
                ij_to_voigt = 1
            case(12)
                ij_to_voigt = 6
            case(13)
                ij_to_voigt = 5
            case(21)
                ij_to_voigt = 6
            case(22)
                ij_to_voigt = 2
            case(23)
                ij_to_voigt = 4
            case(31)
                ij_to_voigt = 5
            case(32)
                ij_to_voigt = 4
            case(33)
                ij_to_voigt = 3
            case default
                write(*,*)'Error in ij_to_voigt with values of ij', i, j
                stop
        end select 

    end function ij_to_voigt



    subroutine compute_Gamma_NI()
        ! Computes the radial Gamma_NI values (D.209 - D.221) of DT98
    end subroutine compute_Gamma_NI




    real(kind=CUSTOM_REAL) function gammaD1_coeff(row, s, a, c, l, n, f)
    ! Computes value of table D.2 of DT98 
    integer :: row
    integer :: s 
    real(kind=CUSTOM_REAL) :: a, c, l, n , f

    real(kind=CUSTOM_REAL) :: lam(5), coeff, fifteen, twone, i_rt3, i_rt6, two_35, two_7rt10, two_rt70

    fifteen = three * five
    twone   = three*seven
    i_rt3   = one/(three**half)
    i_rt6   = one/(six**half)

    two_35    = two/(five*seven)
    two_7rt10 = two/(seven * (ten**half) )
    two_rt70  = two/((seven*ten)**half)

    ! compute Lambda values 
    lam(1) = c + six * a - four  * l - ten       * n + eight * f
    lam(2) = c +       a + six   * l + five      * n - two   * f
    lam(3) = c - six * a - four  * l + seven*two * n + five  * f
    lam(4) = c +       a + three * l - seven     * n - two   * f
    lam(5) = c +       a - four  * l                 - two   * f


    if(row.lt.1 .or. row.gt.13)then
        write(*,*)'Error in gammaD1_coeff. row must be 1-13 but was ', row
        stop 
    endif 

    select case(s)
        ! S = 0 
        case(0)
            select case(row)
                case(1)
                    coeff = (lam(1) + two*lam(2))/fifteen
                case(2)
                    coeff = lam(2)*two/fifteen
                case(3)
                    coeff = (lam(1) + lam(2))/fifteen
                case(4)
                    coeff = -lam(1)/fifteen
                case(5)
                    coeff = -lam(2)/fifteen
                case default
                    coeff = zero 
            end select 
        ! S = 2
        case(2)
            select case(row)
                case(1)
                    coeff = four*(lam(3)+ two*lam(4))/twone
                case(2)
                    coeff = -four*lam(4)/twone
                case(3)
                    coeff = -two*(lam(3)+lam(4))/twone
                case(4)
                    coeff = -lam(3)/twone
                case(5)
                    coeff = -lam(4)/twone
                case(6)
                    coeff = (lam(3) + two*lam(4))*i_rt3/seven
                case(7)
                    coeff = -two * lam(4) * i_rt3 / seven
                case(8)
                    coeff = -(lam(3)+lam(4))*i_rt3/seven
                case(9)
                    coeff = two*lam(3)*i_rt6/seven
                case(10)
                    coeff = two*lam(4)*i_rt6/seven    
                case(11)
                    coeff = -two*(lam(3)+ two*lam(4)) * i_rt6/seven
                case default
                    coeff = zero 
            end select 
        ! S = 4
        case(4)
            select case(row)
                case(1)
                    coeff = two_35 * four
                case(2)
                    coeff = two_35 
                case(3)
                    coeff = two_35 
                case(4)
                    coeff = two_35 * two 
                case(5)
                    coeff = two_35 * two 
                case(6)
                    coeff = two_7rt10 * two 
                case(7)
                    coeff = two_7rt10 
                case(8)
                    coeff = two_7rt10 
                case(9)
                    coeff = two_7rt10 * two 
                case(10)
                    coeff = two_7rt10 * two 
                case(11)
                    coeff = two_7rt10 
                case(12)
                    coeff = two_rt70 
                case(13)
                    coeff = two_rt70 * two 
            end select

            coeff = coeff * lam(5)

        ! Error s != 0, 2, 4
        case default
            write(*,*)'Error in gammaD1_coeff. S must be 0, 2, 4 but was ', s
            stop
    end select 

    ! Compute coefficient including the s term
    gammaD1_coeff =  coeff
    end function gammaD1_coeff




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




real(kind=SPLINE_REAL) function integrate_GNIr2(s, l, N, I, u, du, v, dv, nlen, rvals, dA, dC, dL, dN, dF, type)
    ! Computes the integral in D.208
    use integrate, only: integrate_r_traps
    implicit none 
    include 'constants.h'
    integer :: N, I, nlen, s, l, coeff_ind, w1, w2, w3, h
    real(kind=CUSTOM_REAL) :: dA(nlen), dC(nlen), dL(nlen), dN(nlen), dF(nlen)
    real(kind=CUSTOM_REAL) :: rvals(nlen)
    real(kind=SPLINE_REAL) :: u(nlen), v(nlen), du(nlen), dv(nlen), lf, kf, om2, om0, om2sq, om0sq
    real(kind=SPLINE_REAL) :: gam(nlen), r(nlen), r2(nlen), coeff, xi(nlen), chi(nlen), phi(nlen), integrand(nlen)
    character :: type

    ! l and k in floats
    lf = real(l, kind=SPLINE_REAL)
    kf = (lf*(lf+SPLINE_ONE))**SPLINE_HALF

    ! Compute r^2
    r  = real(rvals, kind=SPLINE_REAL)
    r2 = r*r

    om0 = OmNL(0,l)
    om2 = OmNL(2,l)


    om0sq = om0*om0
    om2sq = om2*om2

    if (type.eq.'S')then 
        chi = dv*r + u - v
        phi = (SPLINE_TWO*u - kf*kf*v)
        xi  = SPLINE_ZERO
    elseif(type.eq.'T')then 
        chi = SPLINE_ZERO
        phi = SPLINE_ZERO
        xi  = du*r - u
    else
        write(*,*)'Error: type must be S, or T: ', type
    endif 


    ! Generate integral: 
    select case(N)
        case(0) ! N = 0
            coeff_ind = I
            select case(I)
                case(1)
                    gam = du * du  * r2
                    w1 = 0 
                    w2 = 0
                    w3 = 0
                case(2)
                    gam = SPLINE_TWO * om0sq * om2sq  
                    if (type.eq.'S')then 
                        gam = gam * v * v
                    elseif(type.eq.'T')then 
                        gam = gam * u * u
                    endif 
                    w1 = -2 
                    w2 = 0
                    w3 = 2

                case(3)
                    gam = phi*phi 
                    w1 = 0 
                    w2 = 0
                    w3 = 0
                case(4)
                    gam = -SPLINE_TWO * phi * du  * r
                    w1 = 0 
                    w2 = 0
                    w3 = 0
                case(5)
                    gam = SPLINE_TWO * om0sq * (chi*chi + xi*xi)
                    w1 = -1
                    w2 = 0
                    w3 = 1
                case default
                    write(*,*)'Error in I value', I
                    stop
            end select

        case(1) ! N = 1
            coeff_ind = 5 + I 
            select case (I)
                case(1)
                    gam = - SPLINE_FOUR * om0 * chi * du * r 
                    w1 = -1
                    w2 = 1
                    w3 = 0
                case(2)
                    gam = - SPLINE_FOUR * om0sq * om2 
                    if (type.eq.'S')then 
                        gam = gam * v * chi
                    elseif(type.eq.'T')then 
                        gam = gam * u * xi
                    endif 
                    w1 = -2
                    w2 = 1
                    w3 = 1

                case(3)
                    gam = SPLINE_FOUR * om0 * chi * phi 
                    w1 = 0
                    w2 = 1
                    w3 = -1
                case default
                    write(*,*)'Error in I value', I
                    stop
            end select 
            

        case(2) ! N = 2
            coeff_ind = 8 + I 
            select case (I)
                case(1)
                    gam = SPLINE_FOUR * om0 * om2 * du * v * r 
                    w1 = -2 
                    w2 = 2
                    w3 = 0
                case(2)
                    gam = SPLINE_TWO * om0sq * (chi*chi - xi*xi) 
                    w1 = -1
                    w2 = 2
                    w3 = -1
                case(3)
                    gam = - SPLINE_FOUR * om0 * om2 * v * phi 
                    w1 = -2
                    w2 = 2
                    w3 = 0
                case default
                    write(*,*)'Error in I value', I
                    stop
            end select 


        case(3) ! N = 3
            if(I.eq.1)then 
                coeff_ind = 12 
                gam = -SPLINE_FOUR * om0sq * om2  
                if (type.eq.'S')then 
                    gam = gam * v * chi
                elseif(type.eq.'T')then 
                    gam = -gam * u * xi 
                endif 
                w1 = -2
                w2 = 3
                w3 = -1
            else 
                write(*,*)'Error in I value', I
                stop
            endif

        case(4) ! N = 4
            if(I.eq.1)then 
                coeff_ind = 13 
                gam = SPLINE_TWO * om0sq * om2sq 
                if (type.eq.'S')then 
                    gam = gam * v * v
                elseif(type.eq.'T')then 
                    gam = -gam * u * u
                endif 
                w1 = -2
                w2 = 4
                w3 = -2
            else 
                write(*,*)'Error in I value', I
                stop
            endif
        case default
            write(*,*)'Error in N value', N
            stop
    end select 

    do h = 1, nlen
        integrand(h) =  real(gammaD1_coeff(coeff_ind, s, dA(h), dC(h), dL(h), dN(h), dF(h)), kind=SPLINE_REAL)
    enddo 

    integrand = integrand * gam * real(thrj(l,s,l,w1,w2,w3), kind=SPLINE_REAL)
    integrate_GNIr2 =  integrate_r_traps(rvals, integrand, nlen)

    return
end function integrate_GNIr2





subroutine save_Vani_matrix(l, fname)
    use params, only: Vani
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l

    integer :: row

    complex(kind=SPLINE_REAL) :: Vani_i(2*l + 1)

    write(*,*)'Writing to ', trim(fname)

    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l + 1
        Vani_i =  Vani(row,:)
        write(1,*)real(Vani_i)
    enddo 
    ! Write the imaginary part
    do row =1, 2*l + 1
        Vani_i =  Vani(row,:)
        write(1,*)aimag(Vani_i)
    enddo ! row
    close(1)

end subroutine save_Vani_matrix



end module V_ani


