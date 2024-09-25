module mesh_utils
    
    use math, only: cosp, atan2p, acosp

    implicit none 
    include "constants.h"


    interface map_local_global
        module procedure map_local_global_double_precision
        module procedure map_local_global_real_4
        module procedure map_local_global_complex_4
    end interface map_local_global


    interface map_complex_vector
        module procedure map_complex_vector_4
        module procedure map_complex_vector_8
    end interface map_complex_vector

    contains 
    


    real(kind=SPLINE_REAL) function delta_spline(i,j)
    ! Delta function with same precision as eigenfunctions
    integer :: i,j
        if (i.eq.j)then 
            delta_spline = ONE
            return  
        else
            delta_spline = ZERO
            return  
        endif 
    end function delta_spline


    integer function delta_int(i,j)
    integer :: i,j
        if (i.eq.j)then 
            delta_int = 1
            return  
        else
            delta_int = 0
            return  
        endif 
    end function delta_int



    subroutine map_local_global_real_4(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        real(kind=4) :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
        integer :: direction

        ! Local variables: 
        integer :: ispec, i, j, k
    
        ! Ensure we have ibool ready to use
        call check_ibool_is_defined()

        if (direction.eq.0)then 
            ! Map local to global
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            glob(ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                        enddo 
                    enddo 
                enddo 
            enddo
        elseif(direction.eq.1)then 
            ! Map global to local
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            loc(i,j,k,ispec) = glob(ibool(i,j,k,ispec))
                        enddo 
                    enddo 
                enddo 
            enddo
        else
            write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
            stop
        endif 

        return 
    end subroutine map_local_global_real_4

    subroutine map_local_global_double_precision(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        real(kind=8) :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
        integer :: direction

        ! Local variables: 
        integer :: ispec, i, j, k
    
        ! Ensure we have ibool ready to use
        call check_ibool_is_defined()

        if (direction.eq.0)then 
            ! Map local to global
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            glob(ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                        enddo 
                    enddo 
                enddo 
            enddo
        elseif(direction.eq.1)then 
            ! Map global to local
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            loc(i,j,k,ispec) = glob(ibool(i,j,k,ispec))
                        enddo 
                    enddo 
                enddo 
            enddo
        else
            write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
            stop
        endif 

        return 
    end subroutine map_local_global_double_precision



    subroutine map_local_global_complex_4(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        complex(kind=4) :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
        integer :: direction

        ! Local variables: 
        integer :: ispec, i, j, k
    
        ! Ensure we have ibool ready to use
        call check_ibool_is_defined()

        if (direction.eq.0)then 
            ! Map local to global
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            glob(ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                        enddo 
                    enddo 
                enddo 
            enddo
        elseif(direction.eq.1)then 
            ! Map global to local
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            loc(i,j,k,ispec) = glob(ibool(i,j,k,ispec))
                        enddo 
                    enddo 
                enddo 
            enddo
        else
            write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
            stop
        endif 
        return 
    end subroutine map_local_global_complex_4





    subroutine map_local_global_complex_8(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        complex(kind=8) :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
        integer :: direction

        ! Local variables: 
        integer :: ispec, i, j, k
    
        ! Ensure we have ibool ready to use
        call check_ibool_is_defined()

        if (direction.eq.0)then 
            ! Map local to global
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            glob(ibool(i,j,k,ispec)) = loc(i,j,k,ispec)
                        enddo 
                    enddo 
                enddo 
            enddo
        elseif(direction.eq.1)then 
            ! Map global to local
            do ispec = 1, nspec 
                do i = 1, ngllx 
                    do j = 1, nglly
                        do k = 1, ngllz
                            loc(i,j,k,ispec) = glob(ibool(i,j,k,ispec))
                        enddo 
                    enddo 
                enddo 
            enddo
        else
            write(*,*)'ERROR: direction can only be 1 or 0 but has value ', direction
            stop
        endif 
        return 
    end subroutine map_local_global_complex_8


    subroutine map_complex_vector_4(dim1, loc, glob, direction)
        ! Wrapper of map_local_global_complex that maps vector
        ! arrays from local <--> global 
        ! direction: 0   local  --> global 
        !            1   global --> local 
        use params, only: ngllx, nglly, ngllz, nspec, nglob
        implicit none 
        
        include "precision.h"

        ! IO variables:
        integer :: dim1, direction
        complex(kind=4) :: loc(dim1, ngllx, nglly, ngllz, nspec)
        complex(kind=4) :: glob(dim1, nglob)

        ! Local: 
        complex(kind=4) :: loc_tmp(ngllx, nglly, ngllz, nspec)
        complex(kind=4) :: glob_tmp(nglob)
        integer                   :: i

        if (direction.eq.0)then
            ! local -> global
            do i = 1, dim1
                loc_tmp(:, :, :, :) = loc(i, :, :, :, :)
                call map_local_global_complex_4(loc_tmp, glob_tmp, direction)
                glob(i,:) = glob_tmp
            enddo
        elseif(direction.eq.1)then
            ! global -> local
            do i = 1, dim1
                glob_tmp(:) = glob(i, :)
                call map_local_global_complex_4(loc_tmp, glob_tmp, direction)
                loc(i, :, :, :, :) = loc_tmp(:, :, :, :)
            enddo
        endif  
    end subroutine map_complex_vector_4



    subroutine map_complex_vector_8(dim1, loc, glob, direction)
        ! Wrapper of map_local_global_complex that maps vector
        ! arrays from local <--> global 
        ! direction: 0   local  --> global 
        !            1   global --> local 
        use params, only: ngllx, nglly, ngllz, nspec, nglob
        implicit none 
        
        include "precision.h"

        ! IO variables:
        integer :: dim1, direction
        complex(kind=8) :: loc(dim1, ngllx, nglly, ngllz, nspec)
        complex(kind=8) :: glob(dim1, nglob)

        ! Local: 
        complex(kind=8) :: loc_tmp(ngllx, nglly, ngllz, nspec)
        complex(kind=8) :: glob_tmp(nglob)
        integer                   :: i

        if (direction.eq.0)then
            ! local -> global
            do i = 1, dim1
                loc_tmp(:, :, :, :) = loc(i, :, :, :, :)
                call map_local_global_complex_8(loc_tmp, glob_tmp, direction)
                glob(i,:) = glob_tmp
            enddo
        elseif(direction.eq.1)then
            ! global -> local
            do i = 1, dim1
                glob_tmp(:) = glob(i, :)
                call map_local_global_complex_8(loc_tmp, glob_tmp, direction)
                loc(i, :, :, :, :) = loc_tmp(:, :, :, :)
            enddo
        endif  
    end subroutine map_complex_vector_8








    subroutine compute_rtp_from_xyz(iproc, save)
        ! Converts the stored xyz to r theta phi
        use allocation_module, only: allocate_if_unallocated
        use params, only: xstore, ystore, zstore, rstore, thetastore, phistore, & 
                        ngllx, nglly, ngllz, nspec     
        implicit none
        include "constants.h" 
        integer :: iproc
        logical :: save

        integer ispec, i ,j , k

        call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, rstore)
        call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, thetastore)
        call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, phistore)


        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        rstore(i,j,k,ispec) = (xstore(i,j,k,ispec)**TWO + &
                                               ystore(i,j,k,ispec)**TWO + & 
                                               zstore(i,j,k,ispec)**TWO)**HALF
                                            
                        ! 0 <= phi <= 2pi
                        phistore(i,j,k,ispec) = atan2(ystore(i,j,k,ispec),&
                                                      xstore(i,j,k,ispec))   
                        if(phistore(i,j,k,ispec) .lt. ZERO)& 
                             phistore(i,j,k,ispec)  = TWO_PI + phistore(i,j,k,ispec) 

                        ! 0 <= theta <= pi                        
                        thetastore(i,j,k,ispec) = PI_OVER_TWO -  atan2(zstore(i,j,k,ispec),&
                                                                     (xstore(i,j,k,ispec)**TWO +&
                                                                      ystore(i,j,k,ispec)**TWO)**HALF) 

                        !if (rstore(i,j,k,ispec) .eq. zero)then 
                        !    thetastore(i,j,k,ispec) = zero
                        !    phistore(i,j,k,ispec) = zero
                        !else 
                        !    thetastore(i,j,k,ispec) = acosp(zstore(i,j,k,ispec)/rstore(i,j,k,ispec))
                        !    phistore(i,j,k,ispec)   = atan2p(ystore(i,j,k,ispec),xstore(i,j,k,ispec))
                        !    if (phistore(i,j,k,ispec)  .lt. zero) phistore(i,j,k,ispec)  = phistore(i,j,k,ispec)  + TWO_PI
                        !endif 

                    enddo 
                enddo 
            enddo
        enddo 


        if(save)call save_elem_rtp(iproc)

    end subroutine compute_rtp_from_xyz



    subroutine compute_rotation_matrix()
        ! Computes rotation matrix R that converts a vector in r, theta, phi
        ! to cartesian coordinates
        ! R holds the normalised unit vectors r, theta, phi in cartesian coords
        ! as columns e.g. 
        !     ( r_x  θ_x  φ_x )
        ! R = ( r_y  θ_y  φ_y ) 
        !     ( r_z  θ_z  φ_z )

        use params, only :  x_glob, y_glob, z_glob, nglob, Rmat, verbose
        use allocation_module, only: allocate_if_unallocated
        use math, only: sinp, cosp, sqrtp
        implicit none 
        include "constants.h"

        real(kind=CUSTOM_REAL) :: x, y, z, r, theta, phi, ct, cp, st, sp, norm
        integer :: iglob, p

        !TODO:  check allocations of x, y, z glob

        if(verbose.ge.2) write(*,'(/,a)')'• Computing rotation matrix'


        if(.not.allocated(x_glob))then 
            write(*,*)'ERROR: X_glob not allocated. Stop.'
            stop
        endif 
        if(.not.allocated(y_glob))then 
            write(*,*)'ERROR: Y_glob not allocated. Stop.'
            stop
        endif 
        if(.not.allocated(z_glob))then 
            write(*,*)'ERROR: Z_glob not allocated. Stop.'
            stop
        endif 
        

        
        call allocate_if_unallocated(3, 3, nglob, Rmat)

        do iglob = 1, nglob
            ! Get x, y, z coordinates: 
            x = real(x_glob(iglob), kind=CUSTOM_REAL)
            y = real(y_glob(iglob), kind=CUSTOM_REAL)
            z = real(z_glob(iglob), kind=CUSTOM_REAL)

            r = sqrtp(x ** 2 + y ** 2 + z ** 2)

            if (r.eq.zero)then 
                ! For now we will wont rotate it if the central GLL point
                Rmat(:, :, iglob) = zero
                do p = 1, 3
                    Rmat(p, p, iglob) = one
                enddo 
            else
                ! Note that theta and phi here are not necessarily the same
                ! as in DT98 i.e. theta here will be between -pi/2 and pi/2
                theta = acosp(z/r)
                phi   = atan2p(y,x)
                if (phi .lt. zero) phi = phi + TWO_PI

                ct = real(cosp(theta), kind=CUSTOM_REAL)
                cp = real(cosp(phi),   kind=CUSTOM_REAL)
                st = real(sinp(theta), kind=CUSTOM_REAL)
                sp = real(sinp(phi),   kind=CUSTOM_REAL)

                ! Radial vector
                !norm =((st*cp)**TWO + (st*sp)**TWO + ct**TWO)**half
                Rmat(1, 1, iglob) = st*cp !/ norm
                Rmat(2, 1, iglob) = st*sp !/ norm
                Rmat(3, 1, iglob) = ct    !/ norm

                ! Theta vector
                !norm =((ct*cp)**TWO + (ct*sp)**TWO + st**TWO)**half
                Rmat(1, 2, iglob) = ct*cp !/ norm
                Rmat(2, 2, iglob) = ct*sp !/ norm
                Rmat(3, 2, iglob) = -st   !/ norm

                ! Phi vector
                !norm =(sp*sp + cp*cp)**half
                Rmat(1, 3, iglob) = -sp !/ norm
                Rmat(2, 3, iglob) =  cp !/ norm
                Rmat(3, 3, iglob) = zero
            endif
        enddo 

    end subroutine compute_rotation_matrix


    subroutine setup_global_coordinate_arrays(iproc, save)
        use params, only: nglob, x_glob, y_glob, z_glob, xstore, ystore, zstore
        use allocation_module, only: allocate_if_unallocated
        implicit none 
        integer :: iproc
        logical :: save 

        call allocate_if_unallocated(nglob, x_glob)
        call allocate_if_unallocated(nglob, y_glob)
        call allocate_if_unallocated(nglob, z_glob)
        call map_local_global(xstore, x_glob, 0)
        call map_local_global(ystore, y_glob, 0)
        call map_local_global(zstore, z_glob, 0)
    
        if(save)call save_global_xyz(iproc)


    end subroutine setup_global_coordinate_arrays




    subroutine rotate_complex_vector_rtp_to_xyz(vector)
        use params, only: Rmat, nspec, ngllx, nglly, ngllz, ibool
        implicit none  

        complex(kind=SPLINE_REAL) :: vector(3, ngllx, nglly, ngllz, nspec)

        integer :: i,j,k,ispec


        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        vector(:, i,j,k,ispec) = matmul(real(Rmat(:,:,ibool(i,j,k,ispec)), kind=SPLINE_REAL), & 
                                                        vector(:,i,j,k,ispec))
                    enddo
                enddo
            enddo 
        enddo

    end subroutine rotate_complex_vector_rtp_to_xyz




    subroutine rotate_complex_sym_matrix_rtp_to_xyz(symmat)
        ! Assumes the matrix is symmetric 3 x 3 in voigt ordering for matrix M
        ! 1  =  M_rr     4  =  M_tp
        ! 2  =  M_tt     5  =  M_rp
        ! 3  =  M_pp     6  =  M_rt

        use params, only: Rmat, nspec, ngllx, nglly, ngllz, ibool
        implicit none 

        complex(kind=SPLINE_REAL) :: symmat(6, ngllx, nglly, ngllz, nspec), M(3,3)
        real(kind=SPLINE_REAL) ::  R(3,3)

        integer :: i, j, k, ispec

        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz

                        R = real(Rmat(:,:,ibool(i,j,k,ispec)), kind=SPLINE_REAL)

                        ! Build back to 3 x 3
                        M(1,1) = symmat(1, i, j, k, ispec)
                        M(2,2) = symmat(2, i, j, k, ispec)
                        M(3,3) = symmat(3, i, j, k, ispec)

                        M(2,3) = symmat(4, i, j, k, ispec)
                        M(3,2) = symmat(4, i, j, k, ispec)

                        M(1,3) = symmat(5, i, j, k, ispec)
                        M(3,1) = symmat(5, i, j, k, ispec)

                        M(1,2) = symmat(6, i, j, k, ispec)
                        M(2,1) = symmat(6, i, j, k, ispec)

                        M = matmul(matmul(R, M), transpose(R))

                        ! Store back in symmat
                        symmat(1, i, j, k, ispec) = M(1,1) ! xx
                        symmat(2, i, j, k, ispec) = M(2,2) ! yy
                        symmat(3, i, j, k, ispec) = M(3,3) ! zz
                        symmat(4, i, j, k, ispec) = M(2,3) ! yz
                        symmat(5, i, j, k, ispec) = M(1,3) ! xz
                        symmat(6, i, j, k, ispec) = M(1,2) ! xy

                        

                    enddo
                enddo
            enddo 
        enddo


    end subroutine rotate_complex_sym_matrix_rtp_to_xyz



    subroutine cleanup_for_mode()
        use params
        use allocation_module
        implicit none 
    
        call deallocate_if_allocated(xstore)
        call deallocate_if_allocated(ystore)
        call deallocate_if_allocated(zstore)
        call deallocate_if_allocated(x_glob)
        call deallocate_if_allocated(y_glob)
        call deallocate_if_allocated(z_glob)
        call deallocate_if_allocated(ibool)
        call deallocate_if_allocated(strain1)
        call deallocate_if_allocated(globalstrain)

        call deallocate_if_allocated(disp1)
        call deallocate_if_allocated(disp2)
        call deallocate_if_allocated(rad_id)
        call deallocate_if_allocated(unique_r)
        call deallocate_if_allocated(rstore)
        call deallocate_if_allocated(thetastore)
        call deallocate_if_allocated(phistore)
        call deallocate_if_allocated(u_spl)
        call deallocate_if_allocated(udot_spl)
        call deallocate_if_allocated(v_spl)
        call deallocate_if_allocated(vdot_spl)
        call deallocate_if_allocated(rho_spl)
        call deallocate_if_allocated(interp_id_r)
        call deallocate_if_allocated(xx)
        call deallocate_if_allocated(zz)

        call deallocate_if_allocated(detjac)
        call deallocate_if_allocated(jacinv)
        call deallocate_if_allocated(jac)


        call deallocate_if_allocated(Rmat)



    end subroutine cleanup_for_mode
    
    
    subroutine check_ibool_is_defined()
        ! Checks that ibool is allocated and not just zero
        ! TODO: call ibool loading if not allocated? 
        use params, only: ibool,nglob
        implicit none 
    
        if(.not. allocated(ibool))then 
            write(*,*)'ERROR: ibool is not allocated but is about to be used.'
            stop
        else 
            if (nglob.ne.maxval(ibool))then
                write(*,*)'ERROR: nglob is not equal to the maximum ibool value'
                write(*,*)'nglob    : ', nglob
                write(*,*)'max ibool: ', maxval(ibool)
                stop
            endif
        endif 
    end subroutine check_ibool_is_defined
    
    
    
    subroutine load_ibool(iproc, region)
        use params
        use allocation_module, only: allocate_if_unallocated
        implicit none 
        integer :: iproc, region
        character(len=250) :: varname 
    
        call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, ibool)
        varname = 'ibool'
        call read_integer_proc_variable(iproc, region, ibool, varname)
    
        call check_ibool_is_defined()
    end subroutine
    

    subroutine compute_jacobian(iproc, save_jac)
        use params, only: jac, detjac, nspec, ngllx, nglly, ngllz, & 
                          xstore, ystore, zstore, dgll,verbose, jacinv
        use allocation_module, only: allocate_if_unallocated
        use math, only: mat_inv
        implicit none 
        include "constants.h"

        ! Save: 
        logical :: save_jac
        integer :: iproc 
        ! Local variables: 
        integer :: i, j, s, t, n, p, ispec
        real(kind=CUSTOM_REAL) :: val, jl(3,3), tmp(3,3)
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), pointer :: cc


        if(verbose.ge.2)then
            write(*,'(/,a)')'• Setting up jacobian'
        endif 

        call allocate_if_unallocated(3, 3, ngllx, nglly, ngllz, nspec, jac)
        call allocate_if_unallocated(3, 3, ngllx, nglly, ngllz, nspec, jacinv)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, detjac)

        detjac = zero 
        jac    = zero 

        do ispec = 1, nspec
            do s = 1, ngllx
                do t = 1, nglly
                    do n = 1, ngllz
                        ! Jij = del x_i/del Xi_j

                        ! Loop over x,y,z
                        do i = 1,3 
                            ! cc becomes x, y, or z 
                            if     (i .eq. 1) then 
                                cc => xstore
                            elseif (i .eq. 2) then 
                                cc => ystore
                            else
                                cc => zstore
                            endif 
                            
                            ! Loop over xi, eta, zeta
                            do j = 1, 3
                                val = zero

                                do p = 1, ngllx
                                    if (j.eq.1) then     ! xi 
                                        val = real(val + cc(p, t, n, ispec) * dgll(p, s), kind=CUSTOM_REAL)
                                    elseif (j.eq.2) then ! eta 
                                        val = real(val + cc(s, p, n, ispec) * dgll(p, t), kind=CUSTOM_REAL)
                                    else                 ! zeta
                                        val = real(val + cc(s, t, p, ispec) * dgll(p, n), kind=CUSTOM_REAL)
                                    endif 
                                enddo                                 
                                jac(i,j,s,t,n,ispec) = val 
                            enddo !j 
                        enddo! i 

                        ! Compute determinant for this GLL point 
                        jl(:,:) = jac(:,:,s,t,n,ispec)

                        detjac(s,t,n,ispec) = jl(1,1) * ( jl(2,2)*jl(3,3) - jl(2,3)*jl(3,2))  - & 
                                              jl(1,2) * ( jl(2,1)*jl(3,3) - jl(2,3)*jl(3,1))  + & 
                                              jl(1,3) * ( jl(2,1)*jl(3,2) - jl(2,2)*jl(3,1)) 
                        
                        ! Compute the inverse jacobian: 
                        ! J^{-1}ij = del Xi_i/del x_j
                        jacinv(:,:,s,t,n,ispec) = mat_inv(jl)

                    enddo ! n 
                enddo !t
            enddo ! s
        enddo ! ispec

        if(verbose.ge.2)write(*,'(a)')'  --> done'


        if(save_jac)then 
            if(verbose.ge.2)write(*,'(a)')'  --> saving jacobian to binary'
            call save_jacobian(iproc)
        endif 
    end subroutine compute_jacobian
    

    
    subroutine read_proc_coordinates(iproc, region)
        use params, only: ngllx, nglly, ngllz, nspec, & 
                          xstore, ystore, zstore, nglob, &
                          datadir, verbose, IIN
        use allocation_module, only: allocate_if_unallocated
        implicit none 
        
        ! IO variables: 
        integer, intent(in) :: iproc
        integer, intent(in) :: region
    
        ! Local variables
        integer :: ier
        character(len=550) :: binname
        
        double precision, allocatable :: xstore_dp(:,:,:,:)
        double precision, allocatable :: ystore_dp(:,:,:,:)
        double precision, allocatable :: zstore_dp(:,:,:,:)
    
        ! File name prefix: 
        write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',iproc,'_'//'reg',region,'_'


        if(verbose.gt.1)then
            write(*,'(/,a,/)')'• Reading files from '//trim(datadir)
            write(*,'(a,i1)')'  -- region      : ', region
            write(*,'(a,i0.6,/)')'  -- processor id: ', iproc
        endif
    
        ! Read processor info to get ngll and nspec
        open(unit=IIN,file=trim(binname)//'info.bin', &
        status='unknown',form='unformatted',action='read',iostat=ier)
        if (ier.ne.0)then 
            write(*,'(a, i0.6)')'Couldnt read info file for proc ', iproc
            stop
        endif 
        read(IIN)nglob
        read(IIN)nspec
        read(IIN)ngllx
        read(IIN)nglly
        read(IIN)ngllz
        close(IIN)
            
        if(verbose.ge.3)then
            write(*,'(a)')'  -- Info: '
            write(*,*)'    --> nglob: ', nglob
            write(*,*)'    --> nspec: ', nspec
            write(*,'(a,i1)')'     --> ngllx: ', ngllx
            write(*,'(a,i1)')'     --> nglly: ', nglly
            write(*,'(a,i1)')'     --> ngllz: ', ngllz
        endif 
    
        ! Allocate mesh arrays:        
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, xstore_dp)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, xstore)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, ystore_dp)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, ystore)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, zstore_dp)
        call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, zstore)
    
    
        ! Open the x coordinate and load: 
        open(unit=IIN,file=trim(binname)//'xstore.bin', &
        status='unknown',form='unformatted',action='read',iostat=ier)
        if (ier.ne.0)then 
            write(*,'(a, i0.6)')'Couldnt read xstore file for proc ', iproc
            stop
        endif 
        read(IIN)xstore_dp
        close(IIN)

        ! Open the y coordinate and load: 
        open(unit=IIN,file=trim(binname)//'ystore.bin', &
        status='unknown',form='unformatted',action='read',iostat=ier)
        if (ier.ne.0)then 
            write(*,'(a, i0.6)')'Couldnt read ystore file for proc ', iproc
            stop
        endif 
        read(IIN)ystore_dp
        close(IIN)
        
        ! Open the z coordinate and load: 
        open(unit=IIN,file=trim(binname)//'zstore.bin', &
        status='unknown',form='unformatted',action='read',iostat=ier)
        if (ier.ne.0)then 
            write(*,'(a, i0.6)')'Couldnt read zstore file for proc ', iproc
            stop
        endif 
        read(IIN)zstore_dp
        close(IIN)


      ! Cast from DP to CUSTOM_REAL (could also be DP)
        xstore(:,:,:,:) = real(xstore_dp(:,:,:,:), kind=CUSTOM_REAL)
        ystore(:,:,:,:) = real(ystore_dp(:,:,:,:), kind=CUSTOM_REAL)
        zstore(:,:,:,:) = real(zstore_dp(:,:,:,:), kind=CUSTOM_REAL)


        if(verbose.ge.3)then
            write(*,'(/,a)')'  -- X coordinates:'
            write(*,*)'     --> min. value: ', minval(xstore)
            write(*,*)'     --> max. value: ', maxval(xstore)

            write(*,'(/,a)')'  -- Y coordinates:'
            write(*,*)'     --> min. value: ', minval(ystore)
            write(*,*)'     --> max. value: ', maxval(ystore)
    
            write(*,'(/,a)')'  -- Z coordinates:'
            write(*,*)'     --> min. value: ', minval(zstore)
            write(*,*)'     --> max. value: ', maxval(zstore)
        endif 

        
        deallocate(xstore_dp)
        deallocate(ystore_dp)
        deallocate(zstore_dp)
    
        return 
    end subroutine read_proc_coordinates
    
    
    
    
    
    subroutine read_proc_variable(iproc, region, variable, varname)
        use params, only: ngllx, nglly, ngllz, nspec, datadir, verbose
    
        implicit none 
        include "precision.h"
    
        ! IO variables: 
        character(len=250) :: varname 
        integer            :: iproc
        integer            :: region
        real(kind=CUSTOM_REAL)   :: variable(ngllx, nglly, ngllz, nspec)
    
        ! Local variables
        integer ::  ier
        character(len=250) :: binname
    
    
        ! File name prefix: 
        write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',iproc,'_'//'reg',region,'_'//trim(varname)//'.bin'
        
        if(verbose.ge.3)then
            write(*,'(/,/,a)')'• Reading variable called '//trim(varname)
            write(*,'(/,a,i1)')'  -- data type : CUSTOM_REAL of length ', CUSTOM_REAL 
            write(*,'(a)')'  -- file name : '//trim(binname)
        endif 
    
        ! Open the variable file and load: 
        open(unit=iproc,file=trim(binname), &
        status='unknown',form='unformatted',action='read',iostat=ier)
        if (ier.ne.0)then 
            write(*,'(a, i0.6)')'Couldnt read "'//trim(varname)//'" file for proc ', iproc
            stop
        endif 
        read(iproc)variable
    
        if(verbose.ge.2)then
            write(*,'(a)')'  -- '//trim(varname)//' :'
            write(*,*)'     --> min. value: ', minval(variable)
            write(*,*)'     --> max. value: ', maxval(variable)
        endif 
        
        close(iproc)

        return 
    
    end subroutine read_proc_variable
    
    
    
    subroutine read_integer_proc_variable(iproc, region, variable, varname)
        use params, only: ngllx, nglly, ngllz, nspec,datadir, verbose, IIN
    
        implicit none 
        include "precision.h"
    
        ! IO variables: 
        character(len=250) :: varname 
        integer            :: iproc
        integer            :: region
        integer            :: variable(ngllx, nglly, ngllz, nspec)
    
        ! Local variables
        integer :: ier
        character(len=250) :: binname
    
    
        ! File name prefix: 
        write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',iproc,'_'//'reg',region,'_'//trim(varname)//'.bin'
        if(verbose.ge.3)then
            write(*,'(/,/,a)')'• Reading variable called '//trim(varname)
            write(*,'(a,i1)')'  -- data type : integer'
            write(*,'(a)')'  -- file name : '//trim(binname)
        endif 
    
        ! Open the variable file and load: 
        open(unit=IIN,file=trim(binname), &
        status='unknown',form='unformatted',action='read',iostat=ier)
        if (ier.ne.0)then 
            write(*,'(a, i0.6)')'Couldnt read "'//trim(varname)//'" file for proc ', iproc
            stop
        endif 
        read(IIN)variable
    
        if(verbose.ge.2)then
            write(*,'(a)')'  -- '//trim(varname)//' :'
            write(*,*)'     --> min. value: ', minval(variable)
            write(*,*)'     --> max. value: ', maxval(variable)
        endif 
    
        close(IIN)
        
        return 
    
    end subroutine read_integer_proc_variable
    







end module mesh_utils


