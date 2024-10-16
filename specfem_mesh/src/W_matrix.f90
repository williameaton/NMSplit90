subroutine compute_W_matrix(SM, interp, n1, t1, l1, n2, t2, l2, store)
    use params, only: Wmat, datadir, rho_spl
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use modes, only: Mode, get_mode
    use specfem_mesh, only: SetMesh 
    use piecewise_interpolation, only: InterpPiecewise
    use mineos_model, only: mineos_ptr
    implicit none 
    include "constants.h"

    ! IO variables 
    character             :: t1, t2 
    integer               :: n1, n2, l1, l2 
    type(SetMesh)         :: sm
    type(InterpPiecewise) :: interp
    logical               :: store

    ! Local: 
    integer               :: i, j, k, ispec, m, l_loop
    logical               :: self_coupling
    complex(SPLINE_REAL)  :: sum
    type(Mode)            :: mode_1, mode_2


    
    ! Check if self_coupling  
    if (t1.eq.t2 .and. l1.eq.l2 .and. n1.eq.n2)then 
        self_coupling = .true.
    else
        self_coupling = .false.
    endif 

    ! Get density at each radius 
    call deallocate_if_allocated(rho_spl)
    allocate(rho_spl(interp%n_radial))
    call interp%interpolate_mineos_variable(real(interp%m%rho_mineos, kind=SPLINE_REAL), rho_spl)

    write(*,*)'minmax of rho:', minval(rho_spl), maxval(rho_spl)
    
    if(store)then
        call save_rhospline_binary(sm%unique_r, rho_spl, sm%n_unique_rad, sm%iset)
    endif


    ! Load modes and interpolates eigenfunctions:  
    mode_1 =  get_mode(n1, t1, l1, mineos_ptr)
    call interp%interpolate_mode_eigenfunctions(mode_1)

    if(.not.self_coupling)then 
        mode_2 =  get_mode(n2, t2, l2, mineos_ptr)
        call interp%interpolate_mode_eigenfunctions(mode_2)    
    endif


    allocate(sm%disp1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%disp2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )


    ! In general I believe there are only non-zero values when m1 = m2 
    ! THIS IS NOT TRUE BUT CAN BE USED TO GENERATE THE BINARIES WE NEED 
    if (l2 .gt. l1) then 
        l_loop = l1 
    else 
        l_loop = l2
    endif

    do m = -l_loop,  l_loop

            ! Get 1st displacement 
            call sm%compute_mode_displacement(m, mode_1, sm%disp1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%disp1)
            if(store)then
                call sm%save_mode_disp_binary(n1, t1, l1, m, 1)
            endif

            ! Get 2nd displacement
            if(self_coupling)then 
                sm%disp2(:,:,:,:,:) = sm%disp1(:,:,:,:,:)
            else
                call sm%compute_mode_displacement(m, mode_2, sm%disp2)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%disp2)

                if(store)then
                    call sm%save_mode_disp_binary(n2, t2, l2, m, 2)
                endif
            endif

            sum = SPLINE_iZERO
            ! Compute D.36 integral
            ! Note that \Omega should only be non-zero in the z-direction
            ! Hence the integrand should be 
            ! rho \tilde{s}^*  \dot (-i\Omega S_y, i\Omega s_x, 0)
            ! Note also that naturally the global mode displacement
            ! is computed in r, theta, phi and needs to be converted
            do ispec = 1, sm%nspec 
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        do k = 1, sm%ngllz

                            sum = sum + rho_spl(sm%rad_id(i,j,k,ispec)) * & 
                            (- conjg(sm%disp1(1,i,j,k,ispec))*sm%disp2(2,i,j,k,ispec) +  & 
                               conjg(sm%disp1(2,i,j,k,ispec))*sm%disp2(1,i,j,k,ispec)) * & 
                              sm%detjac(i,j,k,ispec) * sm%wgll(i) * sm%wgll(j) * sm%wgll(k)    
                        enddo 
                    enddo
                enddo
            enddo

            Wmat(m+l1+1, m+l2+1) = Wmat(m+l1+1, m+l2+1) + sum

    enddo ! m
end subroutine compute_W_matrix




subroutine compute_W_matrix_with_stored(sm, t1, l1, n1, l2, t2, n2)
    use params, only: Wmat, rho_spl
    use specfem_mesh, only: SetMesh
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    implicit none 
    include "constants.h"

    type(SetMesh)         :: sm

    integer   :: l1, n1, l2, n2
    character :: t1, t2

    integer :: m, l_loop
    logical :: self_coupling

    complex(SPLINE_REAL) :: sum
    integer :: i, j, k, ispec

    
    ! Check if self_coupling  
    if (t1.eq.t2 .and. l1.eq.l2 .and. n1.eq.n2)then 
        self_coupling = .true.
    else
        self_coupling = .false.
        allocate(sm%disp2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
    endif 
    allocate(sm%disp1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))


    ! Get density at each radius 
    call sm%load_rho_spline()

    ! In general I believe there are only non-zero values when m1 = m2 
    ! We can therefore work out which is the smaller of the two l values and loop over that range
    if (l2 .gt. l1) then 
        l_loop = l1 
    else 
        l_loop = l2
    endif



    do m = -l_loop,  l_loop
            call  sm%load_mode_disp_binary(n1, t1, l1, m, 1)


            sum = SPLINE_iZERO

            if(.not. self_coupling)then 
                call sm%load_mode_disp_binary(n2, t2, l2, m, 2)

                do ispec = 1, sm%nspec 
                    do i = 1, sm%ngllx
                        do j = 1, sm%nglly
                            do k = 1, sm%ngllz

                                sum = sum + rho_spl(sm%rad_id(i,j,k,ispec)) * & 
                                (- conjg(sm%disp1(1,i,j,k,ispec))*sm%disp2(2,i,j,k,ispec) +  & 
                                    conjg(sm%disp1(2,i,j,k,ispec))*sm%disp2(1,i,j,k,ispec)) * & 
                                    sm%detjac(i,j,k,ispec) * sm%wgll(i) * sm%wgll(j) *sm%wgll(k)    
                            enddo 
                        enddo
                    enddo
                enddo
            else
                ! Self coupling so avoiding the second load and copying to disp2
                do ispec = 1, sm%nspec 
                    do i = 1, sm%ngllx
                        do j = 1, sm%nglly
                            do k = 1, sm%ngllz
                                sum = sum + rho_spl(sm%rad_id(i,j,k,ispec)) * & 
                                (- conjg(sm%disp1(1,i,j,k,ispec))*sm%disp1(2,i,j,k,ispec) +  & 
                                    conjg(sm%disp1(2,i,j,k,ispec))*sm%disp1(1,i,j,k,ispec)) * & 
                                    sm%detjac(i,j,k,ispec) * sm%wgll(i) * sm%wgll(j) * sm%wgll(k)    
                            enddo 
                        enddo
                    enddo
                enddo
            endif ! if self coupling
            Wmat(m+l1+1, m+l2+1) = Wmat(m+l1+1, m+l2+1) + sum
    enddo ! m

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