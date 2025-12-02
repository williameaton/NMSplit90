subroutine compute_Vcen_matrix(SM, interp, gpsi, n1, t1, l1, n2, t2, l2)
    use params, only: Vcen, datadir, rho_spl, Z_AXIS_EARTH_ROTATION
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

    real(kind=CUSTOM_REAL) :: gpsi(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec), gpsi_l(3)

    ! Local: 
    integer               :: i, j, k, p, q,ispec, m1, m2, l_loop, pp, qq, mm, rid
    logical               :: self_coupling
    complex(SPLINE_REAL)  :: sum, kap_div1, kap_div2, mucont, rho_cont, ups_cont, & 
                             gp_cont, gp_cont_t1, gp_cont_t2,gp_cont_t3,gp_cont_t4, g2p_cont, & 
                             tmpgrad(3), gradgradphi_ij, cs1(3), s2(3), cgs1(3,3), gs2(3,3), cdiv_s1, &
                             div_s2, cont
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


      ! Load modes and interpolates eigenfunctions:  
    mode_1 =  get_mode(n1, t1, l1, mineos_ptr)
    call interp%interpolate_mode_eigenfunctions(mode_1)

    mode_2 =  get_mode(n2, t2, l2, mineos_ptr)
    call interp%interpolate_mode_eigenfunctions(mode_2)    


    ! Displacement
    allocate(sm%disp1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%disp2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

     ! Grad displacement
    allocate(sm%gradS_1(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%gradS_2(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

    write(*,*)'Min and max of gpsi: ', minval(gpsi), maxval(gpsi)
    write(*,*)'Min and max of rhospline: ', minval(rho_spl), maxval(rho_spl)

    write(*,*)'mode 1', n1, t1, l1
    write(*,*)'mode 2', n2, t2, l2

    do m1 = -l1,  l1
            write(*,*)'M1 value: ', m1

            ! Get 1st mode's displacement & store the radial 
            ! component before rotating
            call sm%compute_mode_displacement(m1, mode_1, sm%disp1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%disp1)

            ! Get 1st mode's gradient of displacement 
            call sm%compute_mode_gradS(m1, mode_1, sm%gradS_1)
            call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_1)

            do m2 = -l2,  l2
                ! Displacement: store the radial component before rotating
                call sm%compute_mode_displacement(m2, mode_2, sm%disp2)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%disp2)

                ! Grad disp and strain deviator
                call sm%compute_mode_gradS(m2, mode_2, sm%gradS_2)
                call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_2)

                gpsi_l(3) = zero 

                sum = SPLINE_iZERO
                do ispec = 1, sm%nspec 
                    do i = 1, sm%ngllx
                        do j = 1, sm%nglly
                            do k = 1, sm%ngllz

                                !gpsi_l(:) = gpsi(:,i,j,k,ispec)
                                gpsi_l(1) = - OMEGA * OMEGA * sm%xstore(i,j,k,ispec)
                                gpsi_l(2) = - OMEGA * OMEGA * sm%ystore(i,j,k,ispec)

                                cs1 = conjg(sm%disp1(:,i,j,k,ispec))
                                s2  =       sm%disp2(:,i,j,k,ispec)


                                cgs1 = conjg(sm%gradS_1(:,:,i,j,k,ispec))
                                gs2  = sm%gradS_2(:,:,i,j,k,ispec)


                                ! Compute divergence of the mode displacements: 
                                cdiv_s1 = cgs1(1,1) + cgs1(2,2) + cgs1(3,3)
                                 div_s2 =  gs2(1,1) +  gs2(2,2) +  gs2(3,3)


                                cont = SPLINE_iZERO
                                do p = 1,3
                                    do q = 1,3 
                                        cont = cont +  gpsi_l(p) * cs1(q) *  gs2(p,q)  + &
                                                       gpsi_l(p) *  s2(q) * cgs1(p,q)
                                    enddo 
                                    cont = cont - gpsi_l(p)*( cs1(p)*div_s2 + s2(p)*cdiv_s1)
                                enddo 
                                cont = cont * half  
                                
                                ! Adding the second term, which is really simple for the case of alignment with the 
                                ! z axis since \nabla\nabla\psi is - omega^2 [1 0 0 // 0 1 0 // 0 0 0 ]
                                if(Z_AXIS_EARTH_ROTATION)then
                                    cont = cont - OMEGA*OMEGA*(cs1(1)*s2(1) +  cs1(2)*s2(2) )
                                else
                                    ! in this case just multiply with the tensor 
                                    write(*,*)'Error in compute_Vcen_matrix: Z_AXIS_EARTH_ROTATION = False nor implemented yet.'
                                    stop 
                                endif 



                                ! Background rho                             
                                sum = sum + (rho_spl(sm%rad_id(i,j,k,ispec)) * cont * sm%wglljac(i,j,k,ispec))

                            enddo 
                        enddo
                    enddo
                enddo

                Vcen(m1+l1+1, m2+l2+1) = Vcen(m1+l1+1, m2+l2+1) + sum
            enddo ! m2
    enddo ! m1

end subroutine compute_Vcen_matrix






subroutine save_Vcen_matrix(l1, l2, fname)
    use params, only: Vcen
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l1, l2
    integer :: col, row

    !Dimensionalisation: 
    ! To contribute to H it is 1/2w0 * Vcen 
    ! Hence the dimensionalised units of Vcen should be ang freq ^2 
    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')real(Vcen(row,col)/(SCALE_T*SCALE_T))
            else 
                write(1,'(E15.6)', advance='yes')real(Vcen(row,col)/(SCALE_T*SCALE_T))
            endif
        enddo 
    enddo 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')aimag(Vcen(row,col)/(SCALE_T*SCALE_T))
            else 
                write(1,'(E15.6)', advance='yes')aimag(Vcen(row,col)/(SCALE_T*SCALE_T))
            endif
        enddo 
    enddo 
    close(1)

end subroutine save_Vcen_matrix


subroutine save_Vell_matrix(l1, l2, fname)
    use params, only: Vell
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l1, l2
    integer :: col, row


    !Dimensionalisation: 
    ! To contribute to H it is 1/2w0 * Vell 
    ! Hence the dimensionalised units of Vell should be ang freq ^2 
    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')real(Vell(row,col)/(SCALE_T*SCALE_T))
            else 
                write(1,'(E15.6)', advance='yes')real(Vell(row,col)/(SCALE_T*SCALE_T))
            endif
        enddo 
    enddo 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')aimag(Vell(row,col)/(SCALE_T*SCALE_T))
            else 
                write(1,'(E15.6)', advance='yes')aimag(Vell(row,col)/(SCALE_T*SCALE_T))
            endif
        enddo 
    enddo 
    close(1)

end subroutine save_Vell_matrix



subroutine save_Tell_matrix(l1, l2, fname)
    use params, only: Tell
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l1, l2
    integer :: col, row

    !Dimensionalisation: 
    ! To contribute to H it is w0/2 * Vell 
    ! Hence the dimensionalised units of Vell should be unitless
    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')real(Tell(row,col))
            else 
                write(1,'(E15.6)', advance='yes')real(Tell(row,col))
            endif
        enddo 
    enddo 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')aimag(Tell(row,col))
            else 
                write(1,'(E15.6)', advance='yes')aimag(Tell(row,col))
            endif
        enddo 
    enddo 
    close(1)

end subroutine save_Tell_matrix