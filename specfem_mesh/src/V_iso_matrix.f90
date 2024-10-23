subroutine compute_Viso_matrix(SM, interp, delta_kap, delta_mu, delta_rho, & 
                               n1, t1, l1, n2, t2, l2, store)
    use params, only: Viso, datadir, rho_spl
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use modes, only: Mode, get_mode
    use specfem_mesh, only: SetMesh 
    use piecewise_interpolation, only: InterpPiecewise
    use mineos_model, only: mineos_ptr
    use gravitation, only: compute_background_g
    implicit none 
    include "constants.h"

    ! IO variables 
    character             :: t1, t2 
    integer               :: n1, n2, l1, l2
    type(SetMesh)         :: sm
    type(InterpPiecewise) :: interp
    logical               :: store
    real(kind=CUSTOM_REAL), dimension(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec):: & 
                             delta_mu, delta_kap, delta_rho

    ! Local: 
    integer               :: i, j, k, ispec, m1, m2, l_loop, pp, qq, rid
    logical               :: self_coupling
    complex(SPLINE_REAL)  :: sum, kap_div1, kap_div2, mucont, rho_cont, ups_cont
    type(Mode)            :: mode_1, mode_2
    complex(SPLINE_REAL), allocatable :: s_rad_1(:,:,:,:), s_rad_2(:,:,:,:)

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

    ! Stores the radial part of the displacement
    allocate(s_rad_1(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(s_rad_2(sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

    ! Grad displacement
    allocate(sm%gradS_1(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%gradS_2(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

    ! Strain deviator
    allocate(sm%devE_1(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec)  )
    allocate(sm%devE_2(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec)  )

    ! Grad phi 
    allocate(sm%gradphi_1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%gradphi_2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

    ! Gradient of S_r 
    allocate(sm%gSR_1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))
    allocate(sm%gSR_2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec))


    do m1 = -l1,  l1
            write(*,*)m1

            ! Get 1st mode's displacement & store the radial 
            ! component before rotating
            call sm%compute_mode_displacement(m1, mode_1, sm%disp1)
            s_rad_1(:,:,:,:) = sm%disp1(1,:,:,:,:)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%disp1)

            ! Get 1st mode's gradient of displacement 
            call sm%compute_mode_gradS(m1, mode_1, sm%gradS_1)
            call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_1)
            call sm%compute_strain_deviator(sm%gradS_1, sm%devE_1)

            ! Get 1st mode's gradient of the gravitational potential 
            call sm%compute_mode_gradphi(m1, mode_1, sm%gradphi_1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%gradphi_1)

            call sm%compute_mode_grad_SR(m1, mode_1, sm%gSR_1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%gSR_1)

            do m2 = -l2,  l2
                write(*,*)'   ', m2
                ! Displacement: store the radial component before rotating
                call sm%compute_mode_displacement(m2, mode_2, sm%disp2)
                s_rad_2(:,:,:,:) = sm%disp2(1,:,:,:,:)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%disp2)

                ! Grad disp and strain deviator
                call sm%compute_mode_gradS(m2, mode_2, sm%gradS_2)
                call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_2)
                call sm%compute_strain_deviator(sm%gradS_2, sm%devE_2)

                ! Get 2nd mode's gradient of the gravitational potential 
                call sm%compute_mode_gradphi(m2, mode_2, sm%gradphi_2)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%gradphi_2)

                call sm%compute_mode_grad_SR(m2, mode_2, sm%gSR_2)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%gSR_2)

                sum = SPLINE_iZERO
                ! Here we are ignoring perturbations to the boundaries
                do ispec = 1, sm%nspec 
                    do i = 1, sm%ngllx
                        do j = 1, sm%nglly
                            do k = 1, sm%ngllz

                                ! Divergences of \nabla\bfs^*_1 and \nabla\bfs_2
                                ! for Kappa kernel 
                                kap_div1 = SPLINE_iZERO
                                kap_div2 = SPLINE_iZERO

                                ! Contraction of strain deviators for mu kernel
                                mucont = SPLINE_iZERO

                                ! Rho kernel contraction 
                                rho_cont = SPLINE_iZERO
                                ups_cont = SPLINE_iZERO

                                do pp = 1, 3
                                    ! Kappa kernel stuff
                                    kap_div1 = kap_div1 + conjg(sm%gradS_1(pp,pp,i,j,k,ispec))
                                    kap_div2 = kap_div2 + sm%gradS_2(pp,pp,i,j,k,ispec)
                            
                                    do qq = 1,3 
                                        ! Mu kernel stuff
                                        mucont = mucont + & 
                                                 (conjg(sm%devE_1(pp,qq,i,j,k,ispec)) * & 
                                                        sm%devE_2(pp,qq,i,j,k,ispec)) 
                                    enddo 

                                    ! Rho kernel stuff
                                    rho_cont = rho_cont + &
                                               conjg(sm%disp1(pp,i,j,k,ispec))*sm%gradphi_2(pp,i,j,k,ispec) &
                                             + conjg(sm%gradphi_1(pp,i,j,k,ispec))*sm%disp2(pp,i,j,k,ispec)
                                
                                    ups_cont = ups_cont + (conjg(sm%disp1(pp,i,j,k,ispec)) * sm%gSR_2(pp,i,j,k,ispec)) + & 
                                                          (conjg(sm%disp2(pp,i,j,k,ispec)) * sm%gSR_1(pp,i,j,k,ispec)) 
                                enddo
                                
                                ups_cont = half*(ups_cont                   - &
                                                 conjg(s_rad_1(i,j,k,ispec)) * kap_div2  - &
                                                 s_rad_2(i,j,k,ispec) * kap_div1)        - & 
                                                two*conjg(s_rad_1(i,j,k,ispec))*s_rad_2(i,j,k,ispec)/sm%rstore(i,j,k,ispec)

                                if(sm%rstore(i,j,k,ispec).eq.zero)ups_cont=SPLINE_ZERO

                                ! Add 4 pi G Rho term for rho knl: 
                                rid = sm%rad_id(i,j,k,ispec)

                                rho_cont = rho_cont  + sm%gmag_at_r(rid)*ups_cont + &    ! +  
                                          (eight * rho_spl(rid)         * &
                                           conjg(s_rad_1(i,j,k,ispec)) * &
                                           s_rad_2(i,j,k,ispec)) 
                                

                                sum = sum + ((delta_kap(i,j,k,ispec)*kap_div1*kap_div2)  + & 
                                             (delta_mu(i,j,k,ispec)*mucont*SPLINE_TWO)   + &
                                             (delta_rho(i,j,k,ispec)*rho_cont)             &
                                            ) * sm%wglljac(i,j,k,ispec)


                            if(sum.ne.sum)then 
                                write(*,*)rho_cont
                                stop 
                            endif 



                            enddo 
                        enddo
                    enddo

                enddo

                Viso(m1+l1+1, m2+l2+1) = Viso(m1+l1+1, m2+l2+1) + sum
            enddo ! m2
    enddo ! m1

end subroutine compute_Viso_matrix






subroutine save_Viso_matrix(l1, l2, fname)
    use params, only: Viso
    implicit none 
    include "constants.h"
    character(len=*) :: fname
    integer :: l1, l2
    integer :: col, row

    open(1,file=trim(fname))
    ! Write the real matrix 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')real(Viso(row,col))
            else 
                write(1,'(E15.6)', advance='yes')real(Viso(row,col))
            endif
        enddo 
    enddo 
    do row =1, 2*l1 + 1
        do col = 1, 2*l2 + 1
            if (col .lt. 2*l2+1)then 
            write(1,'(E15.6)', advance='no')aimag(Viso(row,col))
            else 
                write(1,'(E15.6)', advance='yes')aimag(Viso(row,col))
            endif
        enddo 
    enddo 
    close(1)

end subroutine save_Viso_matrix