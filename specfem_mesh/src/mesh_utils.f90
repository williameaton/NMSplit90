module mesh_utils
    
    implicit none 
    include "constants.h"



    interface delta
        module procedure delta_real
    end interface delta

    contains 
    


    real(kind=CUSTOM_REAL) function delta_real(i,j)
    integer :: i,j
        if (i.eq.j)then 
            delta_real = ONE
            return  
        else
            delta_real = ZERO
            return  
        endif 
    end function delta_real

    
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



    subroutine map_local_global_custom_real(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        real(kind=CUSTOM_REAL) :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
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
    end subroutine map_local_global_custom_real

    subroutine map_local_global_double_precision(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        double precision :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
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








    subroutine map_local_global_complex(loc, glob, direction)
        ! mappings between local and global variables
        ! direction: 0   local  --> global 
        !            1   global --> local 

        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! I/O variables: 
        complex(kind=CUSTOM_REAL) :: loc(ngllx,nglly,ngllz,nspec), glob(nglob)
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
    end subroutine map_local_global_complex




    subroutine map_complex_vector(dim1, loc, glob, direction)
        ! Wrapper of map_local_global_complex that maps vector
        ! arrays from local <--> global 
        ! direction: 0   local  --> global 
        !            1   global --> local 
        
        use params, only: ngllx, nglly, ngllz, nspec, nglob, ibool
        implicit none 
        
        include "precision.h"

        ! IO variables:
        integer :: dim1, direction
        complex(kind=CUSTOM_REAL) :: loc(dim1, ngllx, nglly, ngllz, nspec)
        complex(kind=CUSTOM_REAL) :: glob(dim1, nglob)

        ! Local: 
        complex(kind=CUSTOM_REAL) :: loc_tmp(ngllx, nglly, ngllz, nspec)
        complex(kind=CUSTOM_REAL) :: glob_tmp(nglob)
        integer                   :: i


        if (direction.eq.0)then
            ! local -> global
            do i = 1, dim1
                loc_tmp(:, :, :, :) = loc(i, :, :, :, :)
                call map_local_global_complex(loc_tmp, glob_tmp, direction)
                glob(i,:) = glob_tmp
            enddo
        elseif(direction.eq.1)then
            ! global -> local
            do i = 1, dim1
                glob_tmp(:) = glob(i, :)
                call map_local_global_complex(loc_tmp, glob_tmp, direction)
                loc(i, :, :, :, :) = loc_tmp(:, :, :, :)
            enddo
        endif  
    
    end subroutine map_complex_vector










    subroutine compute_rtp_from_xyz()
        ! Converts the stored xyz to r theta phi
        use params, only: xstore, ystore, zstore, rstore, thetastore, phistore, & 
                        ngllx, nglly, ngllz, nspec     
        implicit none
        include "constants.h" 

        integer ispec, i ,j , k

        allocate(rstore(ngllx,nglly,ngllz,nspec))
        allocate(thetastore(ngllx,nglly,ngllz,nspec))
        allocate(phistore(ngllx,nglly,ngllz,nspec))


        do ispec = 1, nspec
            do i = 1, ngllx
                do j = 1, nglly
                    do k = 1, ngllz
                        rstore(i,j,k,ispec) = (xstore(i,j,k,ispec)**TWO + ystore(i,j,k,ispec)**TWO + zstore(i,j,k,ispec)**TWO)**HALF
                                            
                        ! 0 <= phi <= 2pi
                        phistore(i,j,k,ispec) = atan2(ystore(i,j,k,ispec), xstore(i,j,k,ispec))   
                        if(phistore(i,j,k,ispec) .lt. ZERO) phistore(i,j,k,ispec)  = TWO_PI + phistore(i,j,k,ispec) 

                        ! 0 <= theta <= pi                        
                        thetastore(i,j,k,ispec) = PI_OVER_TWO -  atan2(zstore(i,j,k,ispec), (xstore(i,j,k,ispec)**TWO + ystore(i,j,k,ispec)**TWO)**HALF) 

                    enddo 
                enddo 
            enddo
        enddo 

    end subroutine compute_rtp_from_xyz



    subroutine compute_rotation_matrix()
        ! Computes rotation matrix R that converts a vector in r, theta, phi
        ! to cartesian coordinates
        ! R holds the normalised unit vectors r, theta, phi in cartesian coords
        ! as columns e.g. 
        !     ( r_x  θ_x  φ_x )
        ! R = ( r_y  θ_y  φ_y ) 
        !     ( r_z  θ_z  φ_z )

        use params, only :  x_glob, y_glob, z_glob, nglob, Rmat
        use allocation_module, only: allocate_if_unallocated
        implicit none 
        include "constants.h"

        double precision :: x, y, z, r, theta, phi, ct, cp, st, sp, norm
        integer :: iglob


        !TODO:  check allocations of x, y, z glob


        call allocate_if_unallocated(3, 3, nglob, Rmat)

        do iglob = 1, nglob

            ! Get x, y, z coordinates: 
            x = x_glob(iglob) * TWO
            y = y_glob(iglob) * TWO
            z = z_glob(iglob) * TWO

            r = sqrt(x ** 2 + y ** 2 + z ** 2)

            ! Note that theta and phi here are not necessarily the same
            ! as in DT98 i.e. theta here will be between -pi/2 and pi/2
            theta = dacos(z/r)
            phi   = datan2(y,x)

            ct = dcos(theta);  cp = dcos(phi)
            st = dsin(theta);  sp = dsin(phi)

            ! Radial vector
            norm =((st*cp)**TWO + (st*sp)**TWO + ct**TWO)**half
            Rmat(1, 1, iglob) = st*cp / norm
            Rmat(2, 1, iglob) = st*sp / norm
            Rmat(3, 1, iglob) = ct    / norm

            ! Theta vector
            norm =((ct*cp)**TWO + (ct*sp)**TWO + st**TWO)**half
            Rmat(1, 2, iglob) = ct*cp / norm
            Rmat(2, 2, iglob) = ct*sp / norm
            Rmat(3, 2, iglob) = -st   / norm

            ! Phi vector
            norm =(sp*sp + cp*cp)**half
            Rmat(1, 3, iglob) = -sp / norm
            Rmat(2, 3, iglob) =  cp / norm
            Rmat(3, 3, iglob) = 0.0d0
        enddo 

    end subroutine compute_rotation_matrix

end module mesh_utils

