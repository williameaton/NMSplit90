subroutine compute_T_matrix(SM, interp, delta_rho, n1, t1, l1, n2, t2, l2, store)
    use params, only: Tmat, datadir, rho_spl
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
    real(kind=CUSTOM_REAL):: delta_rho(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec)


    ! Local: 
    integer               :: i, j, k, ispec, m1, m2, l_loop
    logical               :: self_coupling
    complex(SPLINE_REAL)  :: sum
    type(Mode)            :: mode_1, mode_2

    
    ! Check if self_coupling  
    if (t1.eq.t2 .and. l1.eq.l2 .and. n1.eq.n2)then 
        self_coupling = .true.
    else
        self_coupling = .false.
    endif 


    ! Load modes and interpolates eigenfunctions:  
    mode_1 =  get_mode(n1, t1, l1, mineos_ptr)
    call interp%interpolate_mode_eigenfunctions(mode_1)

    !if(self_coupling)then 
    !    mode_2 = mode_1
    !else 
    mode_2 =  get_mode(n2, t2, l2, mineos_ptr)
    call interp%interpolate_mode_eigenfunctions(mode_2)    
    !endif

    allocate(sm%disp1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%disp2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )


    do m1 = -l1,  l1
            write(*,*)m1
            ! Get 1st displacement 
            call sm%compute_mode_displacement(m1, mode_1, sm%disp1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%disp1)
            if(store)then
                call sm%save_mode_disp_binary(n1, t1, l1, m1, 1)
            endif

            do m2 = -l2,  l2
                write(*,*)'   ', m2
                call sm%compute_mode_displacement(m2, mode_2, sm%disp2)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%disp2)

                if(store)then
                    call sm%save_mode_disp_binary(n2, t2, l2, m2, 2)
                endif

                sum = SPLINE_iZERO
                ! Here we are ignoring perturbations to the boundaries such that 
                ! the kernel is just \delta \rho (\tilde{\bfs}_k^* \cdot \tilde{\bfs}_{k}) dV
                do ispec = 1, sm%nspec 
                    do i = 1, sm%ngllx
                        do j = 1, sm%nglly
                            do k = 1, sm%ngllz
                                sum = sum + (delta_rho(i,j,k,ispec) * & 
                                    (conjg(sm%disp1(1,i,j,k,ispec))*sm%disp2(1,i,j,k,ispec) +  & 
                                     conjg(sm%disp1(2,i,j,k,ispec))*sm%disp2(2,i,j,k,ispec) +  & 
                                     conjg(sm%disp1(3,i,j,k,ispec))*sm%disp2(3,i,j,k,ispec))   & 
                                    * sm%detjac(i,j,k,ispec) * sm%wgll(i) * sm%wgll(j) * sm%wgll(k)  )
                            enddo 
                        enddo
                    enddo
                enddo

                Tmat(m1+l1+1, m2+l2+1) = Tmat(m1+l1+1, m2+l2+1) + sum
            enddo ! m2
    enddo ! m1

end subroutine compute_T_matrix






subroutine save_T_matrix(l1, l2, fname)
    use params, only: Tmat
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l1, l2

    integer :: col, row

    complex(kind=SPLINE_REAL) :: Tmat_i(2*l2 + 1)

    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')real(Tmat(row,col))
            else 
                write(1,'(E15.6)', advance='yes')real(Tmat(row,col))
            endif
        enddo 
    enddo 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')aimag(Tmat(row,col))
            else 
                write(1,'(E15.6)', advance='yes')aimag(Tmat(row,col))
            endif
        enddo 
    enddo 
    close(1)

end subroutine save_T_matrix