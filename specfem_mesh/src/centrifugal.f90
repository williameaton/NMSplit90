
subroutine compute_grad_centrifugal(SM, gpsi, ggpsi)
    ! We only need the gradient of the centrifugal potential (2.115)
    ! ∇Ψ = Ω x (Ω x r) where r is the position vector 
    ! Since the rotation axis is aligned with the vertical this should
    ! be pretty straight-forward: 
    ! Should be - w^2 (x xhat + y yhat )
    use specfem_mesh, only: SetMesh
    use params, only: verbose, myrank, Z_AXIS_EARTH_ROTATION

    implicit none 
    include "constants.h"

    type(SetMesh)          :: SM
    real(kind=CUSTOM_REAL) :: gpsi(3, SM%ngllx, SM%nglly, SM%ngllz, SM%nspec), result
    real(kind=CUSTOM_REAL) :: ggpsi(3, 3)

    ! Local: 
    real(kind=CUSTOM_REAL) :: omvec(3), tmp1vec(3), resvec(3), posvec(3), tmpgrad(3), minom2
    integer :: i,j,k,ispec, m , pp, qq


    if(verbose.ge.0.and.myrank.eq.0)then 
        write(*,*)'Computing centrifugal gradient...'
        write(*,*)
    endif 

    minom2 = -OMEGA * OMEGA

    if(Z_AXIS_EARTH_ROTATION)then 
        gpsi(1,:,:,:,:) = minom2 * sm%xstore(:,:,:,:)
        gpsi(2,:,:,:,:) = minom2 * sm%ystore(:,:,:,:)
        gpsi(3,:,:,:,:) = zero

        ! The gradient of this is just -omega^2 multiplied by the Identity matrix
        ! except the last diagonal element (z) is 0
        ! and should not be spatially variable 
        ggpsi = zero 
        ggpsi(1, 1) = minom2
        ggpsi(2, 2) = minom2

    else 
            
        ! Probably quicker to compute at the GLL level than global and THEN
        ! map to local level? 
        omvec    = zero 
        omvec(3) = OMEGA

        do ispec = 1, sm%nspec
            do i = 1, sm%ngllx
                do j = 1, sm%nglly
                    do k = 1, sm%ngllz
                        ! Position vetor
                        posvec(1) = sm%xstore(i,j,k,ispec)
                        posvec(2) = sm%ystore(i,j,k,ispec)
                        posvec(3) = sm%zstore(i,j,k,ispec)
                        ! First cross product 
                        tmp1vec(1) = omvec(2) * posvec(3) - omvec(3) * posvec(2)
                        tmp1vec(2) = omvec(3) * posvec(1) - omvec(1) * posvec(3)
                        tmp1vec(3) = omvec(1) * posvec(2) - omvec(2) * posvec(1)
                        ! Second cross product
                        gpsi(1,i,j,k,ispec) = omvec(2) * tmp1vec(3) - omvec(3) * tmp1vec(2)
                        gpsi(2,i,j,k,ispec) = omvec(3) * tmp1vec(1) - omvec(1) * tmp1vec(3)
                        gpsi(3,i,j,k,ispec) = omvec(1) * tmp1vec(2) - omvec(2) * tmp1vec(1)
                    enddo
                enddo
            enddo
        enddo


        ! Comptuing second gradient of psi: 
        do ispec = 1, sm%nspec
            do i = 1, sm%ngllx
                do j = 1, sm%nglly
                    do k = 1, sm%ngllz
                        do pp = 1, 3        ! p and q are the individual grad grad psi 
                            do qq = 1, 3    ! elements (result)

                                tmpgrad = 0 
                                do m = 1, sm%ngllx
                                    tmpgrad(1) = tmpgrad(1) +  gpsi(pp,m,j,k,ispec) * sm%dgll(m, i)
                                    tmpgrad(2) = tmpgrad(2) +  gpsi(pp,i,m,k,ispec) * sm%dgll(m, j)
                                    tmpgrad(3) = tmpgrad(3) +  gpsi(pp,i,j,m,ispec) * sm%dgll(m, k)
                                enddo 
                                result =  tmpgrad(1) * sm%jacinv(1,qq,i,j,k,ispec) + &  ! d xi /d qq
                                          tmpgrad(2) * sm%jacinv(2,qq,i,j,k,ispec) + &  ! d eta /d qq
                                          tmpgrad(3) * sm%jacinv(3,qq,i,j,k,ispec) 
                                ! CURRENTLY NOT RETURNING FOR THIS CASE but hsould work
                            enddo 
                        enddo 
                    enddo 
                enddo 
            enddo 
        enddo

    endif

end subroutine compute_grad_centrifugal





subroutine compute_Vcentrifugal(sm, interp, n1, t1, l1, n2, t2, l2, store)
    use params, only: Vcen, datadir, rho_spl
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
    integer               :: i, j, k, ispec, m, l_loop, pp, qq
    logical               :: self_coupling
    complex(SPLINE_REAL)  :: sum, sum2, cont, sum, tracegrad1, tracegrad2
    type(Mode)            :: mode_1, mode_2
    real(SPLINE_REAL)     :: posvec(3), rhol
     
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

    ! Displacement
    allocate(sm%disp1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%disp2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

    ! Grad displacement
    allocate(sm%gradS_1(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    allocate(sm%gradS_2(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )



    ! Compute diagonal only
    if (l2 .gt. l1) then 
        l_loop = l1 
    else 
        l_loop = l2
    endif
    

    do m = -l_loop, l_loop

            write(*,*)"m ", m

            ! Get 1st displacement 
            call sm%compute_mode_displacement(m, mode_1, sm%disp1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%disp1)
            if(store)then
                call sm%save_mode_disp_binary(n1, t1, l1, m, 1)
            endif

            call sm%compute_mode_gradS(m, mode_1, sm%gradS_1)
            call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_1)



            ! Get 2nd displacement
            if(self_coupling)then 
                sm%disp2(:,:,:,:,:) = sm%disp1(:,:,:,:,:)
                sm%gradS_2(:,:,:,:,:,:) = sm%gradS_1(:,:,:,:,:,:)

            else
                call sm%compute_mode_displacement(m, mode_2, sm%disp2)
                call sm%rotate_complex_vector_rtp_to_xyz(sm%disp2)

                if(store)then
                    call sm%save_mode_disp_binary(n2, t2, l2, m, 2)
                endif

                ! Compute gradient of of S 
                call sm%compute_mode_gradS(m, mode_2, sm%gradS_2)
                call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_2)

            endif

            sum = SPLINE_iZERO
            
            do ispec = 1, sm%nspec 
                do i = 1, sm%ngllx
                    do j = 1, sm%nglly
                        do k = 1, sm%ngllz

                            ! Position vetor (r_x, r_y, r_z)
                            posvec(1) = sm%xstore(i,j,k,ispec)
                            posvec(2) = sm%ystore(i,j,k,ispec)
                            posvec(3) = sm%zstore(i,j,k,ispec)
                            
                            rhol      = rho_spl(sm%rad_id(i,j,k,ispec))

                            sum2 =  - ( conjg(sm%disp1(1,i,j,k,ispec))*sm%disp2(1,i,j,k,ispec) + & 
                                        conjg(sm%disp1(2,i,j,k,ispec))*sm%disp2(2,i,j,k,ispec) )
                            
      

                            tracegrad1 = conjg(sm%gradS_1(1,1,i,j,k,ispec)) + conjg(sm%gradS_1(2,2,i,j,k,ispec)) + conjg(sm%gradS_1(3,3,i,j,k,ispec))
                            tracegrad2 = sm%gradS_2(1,1,i,j,k,ispec) + sm%gradS_2(2,2,i,j,k,ispec) + sm%gradS_2(3,3,i,j,k,ispec)

                            ! contraction 
                            cont      = SPLINE_ZERO  
                            ! Only contract 1-2 because grad psi (z dir) is 0 
                            do pp = 1,2
                                do qq = 1,3
                                   cont = cont - posvec(pp)* &     
                                                 (  conjg(sm%disp1(qq,i,j,k,ispec)) * sm%gradS_2(pp, qq, i,j,k,ispec)  & 
                                                   + sm%disp2(qq,i,j,k,ispec) * conjg(sm%gradS_1(pp, qq, i,j,k,ispec)) & 
                                                   - tracegrad2 * conjg(sm%disp1(pp, i,j,k,ispec))                       & 
                                                   - tracegrad1 * sm%disp2(pp, i,j,k,ispec)                       & 
                                                 )
                                enddo
                            enddo  
                        
                            ! This is the integrand
                            sum = sum + rhol*(half*cont + sum2)*OMEGA*OMEGA*sm%wglljac(i,j,k,ispec) 

                        enddo 
                    enddo
                enddo
            enddo


            Vcen(m+l1+1, m+l2+1) = Vcen(m+l1+1, m+l2+1) + sum

    enddo ! m

end 