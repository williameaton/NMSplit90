subroutine map_local_global_custom_real(loc, glob, direction)
    ! mappings between local and global variables
    ! direction: 0   local  --> global 
    !            1   global --> local 
    use params, only: ngllx, nglly, ngllz, nspec, nglob, CUSTOM_REAL, ibool

    implicit none 

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



