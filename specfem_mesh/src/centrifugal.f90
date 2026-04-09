


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


    ! Load modes and interpolates eigenfunctions:  
    mode_1 =  get_mode(n1, t1, l1, mineos_ptr)
    call interp%interpolate_mode_eigenfunctions(mode_1)

    if(.not.self_coupling)then 
        write(*,*)"Only for self coupling currently"
        stop  
    endif

    ! Displacement
    allocate(sm%disp1(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    !allocate(sm%disp2(3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )

    ! Grad displacement
    allocate(sm%gradS_1(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )
    !allocate(sm%gradS_2(3,3, sm%ngllx, sm%nglly, sm%ngllz, sm%nspec) )


    ! Compute diagonal only
    if (l2 .gt. l1) then 
        l_loop = l1 
    else 
        l_loop = l2
    endif
    

    do m = -l_loop, l_loop

            ! Get 1st displacement 
            call sm%compute_mode_displacement(m, mode_1, sm%disp1)
            call sm%rotate_complex_vector_rtp_to_xyz(sm%disp1)
            if(store)then
                call sm%save_mode_disp_binary(n1, t1, l1, m, 1)
            endif

            call sm%compute_mode_gradS(m, mode_1, sm%gradS_1)
            call sm%rotate_complex_matrix_rtp_to_xyz(sm%gradS_1)

            write(*,*)"done loading "

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

                            sum2 =  - ( conjg(sm%disp1(1,i,j,k,ispec))*sm%disp1(1,i,j,k,ispec) + & 
                                        conjg(sm%disp1(2,i,j,k,ispec))*sm%disp1(2,i,j,k,ispec) )
      
                            tracegrad1 = conjg(sm%gradS_1(1,1,i,j,k,ispec)) + conjg(sm%gradS_1(2,2,i,j,k,ispec)) + conjg(sm%gradS_1(3,3,i,j,k,ispec))

                            tracegrad2 = sm%gradS_1(1,1,i,j,k,ispec) + sm%gradS_1(2,2,i,j,k,ispec) + sm%gradS_1(3,3,i,j,k,ispec)


                            ! contraction 
                            cont      = SPLINE_ZERO  
                            ! Only contract 1-2 because grad psi (z dir) is 0 
                            do pp = 1,2
                                do qq = 1,3
                                   cont = cont - posvec(pp)* &     
                                                 (  conjg(sm%disp1(qq,i,j,k,ispec)) * sm%gradS_1(qq, pp, i,j,k,ispec)  & 
                                                   + sm%disp1(qq,i,j,k,ispec) * conjg(sm%gradS_1(qq, pp, i,j,k,ispec)) )
                                enddo
                                cont = cont + posvec(pp)* ( tracegrad2 * conjg(sm%disp1(pp, i,j,k,ispec))  & 
                                                          + tracegrad1 * sm%disp1(pp, i,j,k,ispec)   )

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