
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


