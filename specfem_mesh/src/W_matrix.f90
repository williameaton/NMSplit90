subroutine compute_W_matrix(t1, l1, n1, t2, l2, n2, store, iproc)
    use params, only: Wmat, rho_spl, n_unique_rad, rho_mineos, IC_ID, & 
     unique_r, interp_id_r, disp1, disp2, nglly, ngllx, ngllz, nspec, NL, &
     rad_id, detjac, wgll, datadir
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use spline, only: interpolate_mineos_variable, interpolate_mode_eigenfunctions
    use mesh_utils, only: rotate_complex_vector_rtp_to_xyz
    implicit none 
    include "constants.h"

    integer, optional :: iproc
    integer :: l1, n1, l2, n2
    character :: t1, t2
    character(len=250):: out_name

    integer ::  m, l_loop
    logical :: self_coupling
    real(SPLINE_REAL) :: wcom, qmod

    real(SPLINE_REAL):: u1(NL), v1(NL), du1(NL), dv1(NL)
    real(SPLINE_REAL):: u2(NL), v2(NL), du2(NL), dv2(NL)

    real(SPLINE_REAL):: u_spl_1(n_unique_rad), v_spl_1(n_unique_rad), & 
                        du_spl_1(n_unique_rad), dv_spl_1(n_unique_rad)
    real(SPLINE_REAL), allocatable :: u_spl_2(:), v_spl_2(:), du_spl_2(:), dv_spl_2(:)
    complex(SPLINE_REAL) :: sum
    integer :: i, j, k, ispec
    logical :: store

    
    ! Check if self_coupling  
    if (t1.eq.t2 .and. l1.eq.l2 .and. n1.eq.n2)then 
        self_coupling = .true.
    else
        self_coupling = .false.
    endif 


    ! Get density at each radius 
    call deallocate_if_allocated(rho_spl)
    allocate(rho_spl(n_unique_rad))
    call interpolate_mineos_variable(real(rho_mineos, kind=SPLINE_REAL), 1, IC_ID, & 
                                     unique_r, n_unique_rad, rho_spl, interp_id_r)
    
    if(store)then
        call save_rhospline_binary(unique_r, rho_spl, n_unique_rad, iproc)
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
    endif


    allocate(disp1(3, ngllx, nglly, ngllz, nspec) )
    allocate(disp2(3, ngllx, nglly, ngllz, nspec) )


    ! In general I believe there are only non-zero values when m1 = m2 
    ! We can therefore work out which is the smaller of the two l values and loop over that range
    if (l2 .gt. l1) then 
        l_loop = l1 
    else 
        l_loop = l2
    endif


    do m = -l_loop,  l_loop

            call compute_global_mode_displacement(t1, l1, m, u_spl_1, v_spl_1, n_unique_rad, disp1)
            call rotate_complex_vector_rtp_to_xyz(disp1)
            
            if(store)then
                call save_mode_disp_binary(n1, t1, l1, m, disp1, iproc)
            endif

            if(self_coupling)then 
                disp2(:,:,:,:,:) = disp1(:,:,:,:,:)
            else
                call compute_global_mode_displacement(t2, l2, m, u_spl_2, v_spl_2, n_unique_rad, disp2)
                call rotate_complex_vector_rtp_to_xyz(disp2)

                if(store)then
                    call save_mode_disp_binary(n2, t2, l2, m, disp2, iproc)
                endif

            endif


            sum = SPLINE_iZERO
            ! Compute D.36 integral
            ! Note that \Omega should only be non-zero in the z-direction
            ! Hence the integrand should be 
            ! rho \tilde{s}^*  \dot (-i\Omega S_y, i\Omega s_x, 0)
            ! Note also that naturally the global mode displacement
            ! is computed in r, theta, phi and needs to be converted
            do ispec = 1, nspec 
                do i = 1, ngllx
                    do j = 1, nglly
                        do k = 1, ngllz

                            sum = sum + rho_spl(rad_id(i,j,k,ispec)) * & 
                            (- conjg(disp1(1,i,j,k,ispec))*disp2(2,i,j,k,ispec) +  & 
                               conjg(disp1(2,i,j,k,ispec))*disp2(1,i,j,k,ispec)) * & 
                              detjac(i,j,k,ispec) * wgll(i) * wgll(j) * wgll(k)    
                        enddo 
                    enddo
                enddo
            enddo

            Wmat(m+l1+1, m+l2+1) = Wmat(m+l1+1, m+l2+1) + sum

    enddo ! m
end subroutine compute_W_matrix




subroutine compute_W_matrix_with_stored(t1, l1, n1, t2, l2, n2, iproc)
    use params, only: Wmat, rho_spl, n_unique_rad, rho_mineos, IC_ID, & 
     unique_r, interp_id_r, disp1, disp2, nglly, ngllx, ngllz, nspec, NL, &
     rad_id, detjac, wgll
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use spline, only: interpolate_mineos_variable, interpolate_mode_eigenfunctions, load_rho_spline
    use mesh_utils, only: rotate_complex_vector_rtp_to_xyz
    implicit none 
    include "constants.h"

    integer :: l1, n1, l2, n2, iproc
    character :: t1, t2

    integer :: m, l_loop
    logical :: self_coupling
    real(SPLINE_REAL) :: wcom, qmod

    real(SPLINE_REAL):: u1(NL), v1(NL), du1(NL), dv1(NL)
    real(SPLINE_REAL):: u2(NL), v2(NL), du2(NL), dv2(NL)

    real(SPLINE_REAL):: u_spl_1(n_unique_rad), v_spl_1(n_unique_rad), & 
                        du_spl_1(n_unique_rad), dv_spl_1(n_unique_rad)
    real(SPLINE_REAL), allocatable :: u_spl_2(:), v_spl_2(:), du_spl_2(:), dv_spl_2(:)
    complex(SPLINE_REAL) :: sum
    integer :: i, j, k, ispec

    
    ! Check if self_coupling  
    if (t1.eq.t2 .and. l1.eq.l2 .and. n1.eq.n2)then 
        self_coupling = .true.
    else
        self_coupling = .false.
    endif 


    ! Get density at each radius 
    call deallocate_if_allocated(rho_spl)
    allocate(rho_spl(n_unique_rad))

    call load_rho_spline(unique_r, rho_spl, n_unique_rad, iproc)


    ! In general I believe there are only non-zero values when m1 = m2 
    ! We can therefore work out which is the smaller of the two l values and loop over that range
    if (l2 .gt. l1) then 
        l_loop = l1 
    else 
        l_loop = l2
    endif


    allocate(disp1(3,ngllx,nglly,ngllz, nspec))


    do m = -l_loop,  l_loop

            call  load_mode_disp_binary(n1, t1, l1, m, disp1, iproc)
            !call compute_global_mode_displacement(t1, l1, m, u_spl_1, v_spl_1, n_unique_rad, disp1)
            !call rotate_complex_vector_rtp_to_xyz(disp1)

            sum = SPLINE_iZERO

            if(.not. self_coupling)then 
                allocate(disp2(3,ngllx,nglly,ngllz, nspec))

                call  load_mode_disp_binary(n2, t2, l2, m, disp2, iproc)
                !call compute_global_mode_displacement(t2, l2, m, u_spl_2, v_spl_2, n_unique_rad, disp2)
                !call rotate_complex_vector_rtp_to_xyz(disp2)
                do ispec = 1, nspec 
                    do i = 1, ngllx
                        do j = 1, nglly
                            do k = 1, ngllz

                                sum = sum + rho_spl(rad_id(i,j,k,ispec)) * & 
                                (- conjg(disp1(1,i,j,k,ispec))*disp2(2,i,j,k,ispec) +  & 
                                    conjg(disp1(2,i,j,k,ispec))*disp2(1,i,j,k,ispec)) * & 
                                    detjac(i,j,k,ispec) * wgll(i) * wgll(j) * wgll(k)    
                            enddo 
                        enddo
                    enddo
                enddo
            else
                ! Self coupling so avoiding the second load and copying to disp2
                do ispec = 1, nspec 
                    do i = 1, ngllx
                        do j = 1, nglly
                            do k = 1, ngllz
                                sum = sum + rho_spl(rad_id(i,j,k,ispec)) * & 
                                (- conjg(disp1(1,i,j,k,ispec))*disp1(2,i,j,k,ispec) +  & 
                                    conjg(disp1(2,i,j,k,ispec))*disp1(1,i,j,k,ispec)) * & 
                                    detjac(i,j,k,ispec) * wgll(i) * wgll(j) * wgll(k)    
                            enddo 
                        enddo
                    enddo
                enddo
            endif ! if self coupling

            Wmat(m+l1+1, m+l2+1) = Wmat(m+l1+1, m+l2+1) + sum
    enddo ! m



    deallocate(disp1)
    if(allocated(disp2))deallocate(disp2)

end subroutine compute_W_matrix_with_stored



subroutine save_W_matrix(l1, l2, fname)
    use params, only: Wmat
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l1, l2

    integer :: col, row

    complex(kind=SPLINE_REAL) :: Wmat_i(2*l2 + 1)

    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')real(Wmat(row,col))
            else 
                write(1,'(E15.6)', advance='yes')real(Wmat(row,col))
            endif
        enddo 
    enddo 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')aimag(Wmat(row,col))
            else 
                write(1,'(E15.6)', advance='yes')aimag(Wmat(row,col))
            endif
        enddo 
    enddo 
    close(1)

end subroutine save_W_matrix