program split_mesh
    use params, only: nglob, nspec, ngllx, nglly, ngllz, rstore, nspec, & 
                      ibool, x_glob, y_glob, z_glob, xstore, ystore, zstore, & 
                      thetastore, datadir, verbose
    use allocation_module, only: allocate_if_unallocated, deallocate_if_allocated
    use mesh_utils, only: compute_rtp_from_xyz, map_local_global_double_precision,& 
                          load_ibool, read_proc_coordinates
    implicit none 
    include "constants.h"
    
    integer :: iproc, i, j, k, ispec, iii, ib_orig, ib_ctr, & 
               this_nodes_id, i_newspec, nprocs, region,    &
               nspec_new, IOUT          
    character(len=250) :: varname, proc_fmtname
    double precision :: eproc, tolerance
    real(kind=CUSTOM_REAL) :: rr

    real(kind=CUSTOM_REAL), allocatable :: meshtag(:,:), tag_global(:), tag(:)
    double precision, allocatable :: elmtsum(:), xcoord_new(:,:,:,:), & 
                                     ycoord_new(:,:,:,:), zcoord_new(:,:,:,:)
    integer, allocatable :: ibool_map(:), ibool_new(:,:,:,:)


    ! Setup parameters: 
    region = 3      ! Inner core
    nprocs = 6
    
    tolerance = 1.0e-9 ! Tolerance for matching elements 

    ! Loop through all processors
    do iproc = 0, nprocs -1 

        ! Get mesh for this processor
        call read_proc_coordinates(iproc, region)
        call load_ibool(iproc, region)
        call compute_rtp_from_xyz()

        if(iproc.eq.0)then 
            ! if on first processor
            ! get the element sum for the initial node
            ! hopefully this is unique for each element in the mesh
            ! to within the tolerance defined
            allocate(meshtag(nprocs, nspec))
            allocate(elmtsum(nspec))
            allocate(tag_global(nglob))

            meshtag = zero
            elmtsum = zero
            ! Loop through the elements 
            do ispec = 1, nspec 
                do i = 1, ngllx
                    do j = 1, nglly 
                        do k = 1, ngllz 
                            elmtsum(ispec) = elmtsum(ispec) + & 
                                             xstore(i,j,k,ispec) + ystore(i,j,k,ispec) + zstore(i,j,k,ispec)& 
                                           + rstore(i,j,k,ispec) + thetastore(i,j,k,ispec)/PI
                        enddo 
                    enddo 
                enddo 
            enddo 
        else
            ! After we have done proc 0 
            ! For the other processors: 

            ! Get mesh for this processor
            call read_proc_coordinates(iproc, region)
            call load_ibool(iproc, region)
            call compute_rtp_from_xyz()

            ! Now we have the original elmtsum we can see if it is present in another processor
            ! I.e. does the sum in an element on this proc match a sum on the main proc (0)
            do ispec = 1, nspec 

                ! Compute element sum for element on this processor
                eproc = zero
                do i = 1, ngllx
                    do j = 1, nglly 
                        do k = 1, ngllz 
                            eproc = eproc + & 
                            xstore(i,j,k,ispec) + ystore(i,j,k,ispec) + zstore(i,j,k,ispec) & 
                            + rstore(i,j,k,ispec) + thetastore(i,j,k,ispec)/PI
                        enddo 
                    enddo 
                enddo 

                ! Now we need to search for this element sum in the original array 
                ! Loop over each element in main proc (0) and see if the sum
                ! matches the one on this processor
                do iii = 1, nspec
                    if ( abs(eproc-elmtsum(iii)) .lt. tolerance) then 
                        ! Elements are shared on these procs -- inner cube
                        ! Tag this element on both
                        meshtag(1,iii) = ONE
                        meshtag(iproc+1,ispec) = ONE
                    endif 
                enddo
            enddo 
        endif 


        ! For each processor other than the main proc, we now have tags 
        ! of every element that is shared in the inner cube
        if (iproc.gt.0)then
            
            ! Ibool map: 
            !  has dimensions of the original ibool array and stores the
            !  new ibool value in the same location that ibool stores
            !  the old ibool id 
            call allocate_if_unallocated(nglob, ibool_map)
            call allocate_if_unallocated(nspec, tag)
            ibool_map = 0 
            tag = 0

            ! Tags for this processor
            tag(:) = meshtag(iproc+1, :)

            ! The 1's in tag are the elements we DONT want so subtract this many
            ! from the original nspec
            nspec_new = nspec - int(sum(tag))

            if (verbose.ge.2)then
                write(*,*)'Original nspec  ', nspec
                write(*,*)'Shared elements ', int(sum(tag))
                write(*,*)'Updated nspec   ', nspec_new
            endif

            ! allocate the new arrays 
            call allocate_if_unallocated(ngllx, nglly, ngllz, nspec_new,  ibool_new)
            call allocate_if_unallocated(ngllx, nglly, ngllz, nspec_new,  xcoord_new)
            call allocate_if_unallocated(ngllx, nglly, ngllz, nspec_new,  ycoord_new)
            call allocate_if_unallocated(ngllx, nglly, ngllz, nspec_new,  zcoord_new)
            ibool_new = 0

            ! This will keep track of the new ibool IDs
            ib_ctr = 1 
            ! This will keep track of the new ispec ID 
            i_newspec = 1

            ! Loop over all elements in this processor's mesh
            do ispec = 1, nspec
                if(tag(ispec).eq.0)then 
                    ! Element not shared with main proc: keep it! 
                    do i = 1, ngllx
                        do j = 1, nglly
                            do k = 1, ngllz
                                ! get original ibool 
                                ib_orig = ibool(i,j,k,ispec)

                                ! Set the ibool map or retrieve id if already
                                ! present
                                if(ibool_map(ib_orig).eq.0)then
                                    ! ibool map stores the mapping from old 
                                    ! to new ibool
                                    ! if it is equal to zero then this global
                                    ! node has not been given its new ibool yet
                                    ! -- if it is non zero then we want to use
                                    ! the new id it already has 
                                    ibool_map(ib_orig) = ib_ctr
                                    this_nodes_id      = ib_ctr
                                    ib_ctr             = ib_ctr + 1
                                else
                                    ! Use preexisting ID
                                    this_nodes_id = ibool_map(ib_orig)
                                endif

                                ibool_new(i, j, k, i_newspec) = this_nodes_id

                                ! Copy over the relevant coordinates
                                xcoord_new(i,j,k,i_newspec) = xstore(i,j,k,ispec)
                                ycoord_new(i,j,k,i_newspec) = ystore(i,j,k,ispec)
                                zcoord_new(i,j,k,i_newspec) = zstore(i,j,k,ispec)
                            enddo 
                        enddo 
                    enddo 
                    ! Update the new ispec ID
                    i_newspec = i_newspec + 1
                endif 
            enddo !nspec

            IOUT = 1

            write(*,*)iproc 

            
            write(proc_fmtname,'(a,i0.6,a,i1,a)')trim(datadir)//'/sliced/'//'/proc',iproc,'_'//'reg',region,'_'

            open(unit=IOUT,file=trim(proc_fmtname)//'xstore.bin', &
            status='unknown',form='unformatted',action='write')
            write(IOUT) xcoord_new
            close(IOUT)

            open(unit=IOUT,file=trim(proc_fmtname)//'ystore.bin', &
            status='unknown',form='unformatted',action='write')
            write(IOUT) ycoord_new
            close(IOUT)

            open(unit=IOUT,file=trim(proc_fmtname)//'zstore.bin', &
            status='unknown',form='unformatted',action='write')
            write(IOUT) zcoord_new
            close(IOUT)
            
            open(unit=IOUT,file=trim(proc_fmtname)//'ibool.bin', &
            status='unknown',form='unformatted',action='write')
            write(IOUT) ibool_new
            close(IOUT)

            open(unit=IOUT,file=trim(proc_fmtname)//'info.bin', &
            status='unknown',form='unformatted',action='write')
            write(IOUT) maxval(ibool_new) ! new nglob
            write(IOUT) nspec_new
            write(IOUT) NGLLX
            write(IOUT) NGLLY
            write(IOUT) NGLLZ
            close(IOUT)

            ! Clean up for this processor just in case each proc has
            ! different nspec_new such that the new arrays are of
            ! different dimensions: 
            call deallocate_if_allocated(xcoord_new)
            call deallocate_if_allocated(ycoord_new)
            call deallocate_if_allocated(zcoord_new)
            call deallocate_if_allocated(ibool_new)


        endif 
    enddo
    
end program split_mesh