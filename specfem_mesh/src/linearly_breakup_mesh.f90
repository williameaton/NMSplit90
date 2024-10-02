! This program is designed as a 'one-time run' for any NEX
! Split_mesh originally breaksup the SPECFEM mesh so that only unique
! Elements exist for one region (e.g. IC)
! The functionality of this program is to recast that mesh into N linear
! sets so that each can run on one process. The aim is to have each set
! contain an equal number of elements, so that the time for each set 
! should be ~ equal. 
! Currently this will just break up the sets into sequential elements of
! a certain length - there may be a better way to do this, e.g., by radii


program linearly_breakup_mesh 
use params, only: nspec, ngllx, nglly, ngllz, xstore, ystore, zstore, & 
                  ibool, datadir
use mesh_utils, only: read_proc_coordinates, cleanup_for_mode
implicit none 

integer :: region, nprocs_before, nsets, iproc, ntot_elem,   & 
           nspec_per_set, rem_el_in_set, rem_in_proc, & 
           nelms_to_copy, set_id_b, set_id_ff,        & 
           set_id_end, proc_id_ff, proc_id_end, iset, & 
           ib_orig, maxnglob, ib_proc, min_ib_proc,   & 
           max_ib_proc, min_ib, max_ib, nibool, iiprc,&
           ispec, i,j,k,iib, count, start_reg, end_reg, & 
           iregion, procs_in_set, this_nodes_id, hhh, & 
           max_val_to_start, procctr
double precision, allocatable :: xcoord_new(:,:,:,:), & 
                                 ycoord_new(:,:,:,:), & 
                                 zcoord_new(:,:,:,:)
integer, allocatable :: ibool_new(:,:,:,:), ib_store(:,:,:,:), & 
                        ib_store_proc(:,:,:,:), ibool_map(:,:),& 
                        ib_ctr(:), ibool_new_pproc(:,:,:,:,:)
logical :: get_new_proc, filled_set
character(len=400) :: outname
character(len=30) :: fmtstr

logical :: calc_ibool

! Setup parameters: 
calc_ibool       = .true.
start_reg        = 1
end_reg          = 3     ! Inner core
nprocs_before    = 6     ! Current setup 
nsets            = 8     ! new setup 

ntot_elem = 0 

! Count the total number of elements
do iregion = start_reg, end_reg
    do iproc = 0, nprocs_before - 1 
        call read_proc_coordinates(iproc, iregion)  
        ntot_elem = ntot_elem + nspec
        call cleanup_for_mode()

        write(*,*)'Nspec = ', nspec
    enddo 
enddo 

write(*,*)'Total elements: ', ntot_elem



! TODO: Could be more flexible to cases without absolutley perfect division
if(mod(ntot_elem,nsets).ne.0)then 
    write(*,*)'Error: total number of elements is not divisible by nsets'
    stop 
endif 

nspec_per_set = ntot_elem/nsets
write(*,*)' Total nspec      : ', ntot_elem
write(*,*)' Number of set    : ', nsets
write(*,*)' Elements per set : ', nspec_per_set


maxnglob = ngllx*nglly*ngllz*nspec_per_set

! Allocate arrays of correct size for a single set: 
allocate(xcoord_new(ngllx, nglly, ngllz, nspec_per_set)) 
allocate(ycoord_new(ngllx, nglly, ngllz, nspec_per_set)) 
allocate(zcoord_new(ngllx, nglly, ngllz, nspec_per_set)) 
allocate(ibool_new(ngllx, nglly, ngllz, nspec_per_set)) 
allocate(ib_store(ngllx, nglly, ngllz,   nspec_per_set)) 
allocate(ib_store_proc(ngllx, nglly, ngllz,nspec_per_set)) 



region = start_reg

! Will set remaining_in_processor to the nspec of the prod
! Will set processor_id_fill_from to 1 
iproc = 0
procctr = 0 
call load_new_proc(iproc, region, rem_in_proc, proc_id_ff)
iset  = 0 
rem_el_in_set = nspec_per_set   ! Remaining elements in a set
set_id_ff     = 1               ! ID to fill from in set



do while (iset.lt.nsets) ! is this the correct finish? 

    write(*,*)' Remaining elements in set    :', rem_el_in_set
    write(*,*)' Remaining elements in proc   :', rem_in_proc
    

    if (rem_el_in_set.gt.rem_in_proc) then 
        ! We need more than just this processor - copy over rest of 
        ! this proc 
        
        nelms_to_copy =  rem_in_proc
        get_new_proc  = .true. 

        write(*,*)' More elements needed to fill set than in proc'
        write(*,*)'  -- Set nelms_to_copy to number remaining in PROC = ', rem_in_proc 
        write(*,*)'  -- Set flag to get new proc'
    else  
        ! Can fill set with this processor
        nelms_to_copy =  rem_el_in_set
        get_new_proc  = .false. 

        write(*,*)' Enough elements in proc to fill set'
        write(*,*)'  -- Set nelms_to_copy to number remaining in SET = ', rem_el_in_set 
        write(*,*)'  -- Set flag to keep same proc'
    endif 


    ! Copy over 'nelms_to_copy' pieces of data into the set 
    set_id_end   = set_id_ff  + nelms_to_copy - 1
    proc_id_end  = proc_id_ff + nelms_to_copy - 1 


    write(*,*)' set_id_end  = ', set_id_end
    write(*,*)' proc_id_end = ', proc_id_end

    write(*,*)' copying over the coordinates...'
    xcoord_new(:,:,:,set_id_ff:set_id_end) = xstore(:,:,:,proc_id_ff:proc_id_end)
    ycoord_new(:,:,:,set_id_ff:set_id_end) = ystore(:,:,:,proc_id_ff:proc_id_end)
    zcoord_new(:,:,:,set_id_ff:set_id_end) = zstore(:,:,:,proc_id_ff:proc_id_end)

    ! Store the ibool values AND the proc it was on 
    ib_store(:,:,:,set_id_ff:set_id_end)       = ibool(:,:,:,proc_id_ff:proc_id_end)
    ib_store_proc(:,:,:,set_id_ff:set_id_end)  = procctr


    set_id_ff  = set_id_end  + 1 
    proc_id_ff = proc_id_end + 1 

    rem_el_in_set = rem_el_in_set - nelms_to_copy
    rem_in_proc   = rem_in_proc   - nelms_to_copy

    write(*,*)' New ID to fill set from  = ', set_id_ff



    if(set_id_end.eq.nspec_per_set)then 
        write(*,*)' Set is now filled...'
        ! Compute ibool for this set & output

        ! Determine how many different processors were used to fill set
        procs_in_set  = maxval(ib_store_proc) - minval(ib_store_proc) + 1
        
        ! Reduce iprocs first index to 1 
        ib_store_proc = ib_store_proc - minval(ib_store_proc) + 1

        ! Reduce ibstore so it starts at 1: 
        ib_store = ib_store - minval(ib_store) + 1

        allocate(ibool_map(maxval(ib_store), procs_in_set)) 
        ibool_map = 0

        allocate(ibool_new_pproc(ngllx, nglly, ngllz, nspec_per_set, procs_in_set))
        ibool_new_pproc = 0 

        allocate(ib_ctr(procs_in_set))
        ib_ctr(:) = 1


        
        ! Loop over each GLL in the set: 
        do ispec = 1, nspec_per_set
            do i = 1, ngllx
                do j = 1, nglly   
                    do k = 1, ngllz 
                        ib_orig = ib_store(i,j,k,ispec)
                        ib_proc = ib_store_proc(i,j,k,ispec)

                        if(ibool_map(ib_orig, ib_proc).eq.0)then 
                            ! New ibool value to have 
                            ibool_map(ib_orig, ib_proc) = ib_ctr(ib_proc)
                            this_nodes_id               = ib_ctr(ib_proc)
                            ib_ctr(ib_proc)             = ib_ctr(ib_proc) + 1
                        else 
                            ! Use preexisting ID
                            this_nodes_id = ibool_map(ib_orig, ib_proc)
                        endif

                        ! Store for now
                        ibool_new_pproc(i, j, k, ispec, ib_proc) = this_nodes_id
                    enddo 
                enddo
            enddo 
        enddo

        !At this point we have the ibools for each of the processors that contribute
        !to this, in some continuous order each starting at 1. The next step is 
        !to ensure continuity between the processes
        !if(minval(ibool_new_pproc(:,:,:,:,1)).ne.1)then 
        !    write(*,*)'Error: minval(ibool_new_pproc(:,:,:,:,1).ne.1'
        !    write(*,*)"min value is ", minval(ibool_new_pproc(:,:,:,:,1))
        !    stop 
        !endif
        if(procs_in_set.gt.1)then 
            do hhh = 2, procs_in_set
                max_val_to_start = maxval(ibool_new_pproc(:,:,:,:,hhh-1))
                do ispec = 1, nspec_per_set
                    do i = 1, ngllx
                        do j = 1, nglly   
                            do k = 1, ngllz 
                                if (ibool_new_pproc(i,j,k,ispec,hhh).ne.0)then  
                                    ibool_new_pproc(i,j,k,ispec,hhh) =    & 
                                        ibool_new_pproc(i,j,k,ispec,hhh)+ & 
                                        max_val_to_start
                                endif 
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        endif

        ! Compile into the new ibool 
        do ispec = 1, nspec_per_set
            do i = 1, ngllx
                do j = 1, nglly   
                    do k = 1, ngllz 
                        ibool_new(i,j,k,ispec) = ibool_new_pproc(i,j,k,ispec,ib_store_proc(i,j,k,ispec))
                    enddo 
                enddo 
            enddo 
        enddo 

        deallocate(ibool_map)
        deallocate(ib_ctr)
        deallocate(ibool_new_pproc)

        ! min_ib_proc = minval(ib_store_proc)
        ! max_ib_proc = maxval(ib_store_proc)
        ! min_ib = minval(ib_store)
        ! max_ib = maxval(ib_store)

        ! if(calc_ibool)then 
        !     write(*,*,advance='no')' Reformatting ibool...'
        !     nibool = 1
        !     do iiprc = min_ib_proc, max_ib_proc
        !         do iib = min_ib, max_ib
        !             count = 0
        !             do ispec = 1, nspec_per_set
        !                 do i = 1, ngllx
        !                     do j = 1, nglly
        !                         do k = 1, ngllz
        !                             if(ib_store(i,j,k,ispec).eq.iib .and. ib_store_proc(i,j,k,ispec).eq.iiprc)then
        !                                 ibool_new(i,j,k,ispec) = nibool 
        !                                 count = count + 1
        !                             endif
        !                         enddo 
        !                     enddo 
        !                 enddo
        !             enddo 
        !             if (count.gt.0)  nibool = nibool + 1
        !         enddo !iib
        !     enddo !iiproc
        !     write(*,*,advance='yes')' done'
        ! endif 


        ! Write out to disc
        if(nsets.lt.10)then 
            fmtstr = '(a,i1,a,i0.6,a,i1,a)'
        elseif(nsets.ge.10 .and. nsets.lt.100)then
            fmtstr = '(a,i2,a,i0.6,a,i1,a)'
        elseif(nsets.ge.100 .and. nsets.lt.1000)then
            fmtstr = '(a,i3,a,i0.6,a,i1,a)'
        elseif(nsets.ge.1000 .and. nsets.lt.10000)then
            fmtstr = '(a,i4,a,i0.6,a,i1,a)'
        elseif(nsets.ge.10000 .and. nsets.lt.100000)then
            fmtstr = '(a,i5,a,i0.6,a,i1,a)'
        endif 


        ! If cross multiple regions then use region 0
        ! Else just use the region read in 
        if(end_reg-start_reg.ne.0)then
            iregion = 0 
        else
            iregion = region 
        endif 

        write(outname,trim(fmtstr))trim(datadir)//'/linear/sets',nsets,'/proc',iset,'_'//'reg',iregion,'_'

        open(1,file=trim(outname)//'xstore.bin', form='UNFORMATTED')
        write(1)xcoord_new
        close(1)

        open(2,file=trim(outname)//'ystore.bin', form='UNFORMATTED')
        write(2)ycoord_new
        close(2)

        open(3,file=trim(outname)//'zstore.bin', form='UNFORMATTED')
        write(3)zcoord_new
        close(3)

        if(calc_ibool)then
            open(4,file=trim(outname)//'ibool.bin', form='UNFORMATTED')
            write(4)ibool_new
            close(4)
        endif


        open(5,file=trim(outname)//'info.bin', form='UNFORMATTED')
        write(5)maxval(ibool_new)
        write(5) nspec_per_set
        write(5) NGLLX
        write(5) NGLLY
        write(5) NGLLZ
        close(5)


        ! Reset the set parameters : 
        write(*,*)'  --> resets set_id_ff to 1'

        iset = iset + 1
        rem_el_in_set = nspec_per_set   ! Remaining elements in a set
        set_id_ff     = 1               ! ID to fill from in set

        write(*,*)'  --> sets rem_el_in_set to ', nspec_per_set
        write(*,*)'  --> iset is now  ', iset
    endif 


    if(get_new_proc)then 
        call cleanup_for_mode()
        if (iproc.eq.nprocs_before-1)then 
            ! Need a new region:
            iproc = 0
            region = region + 1
        else 
            iproc = iproc + 1
        endif 
        procctr = procctr + 1
        call load_new_proc(iproc, region, rem_in_proc, proc_id_ff)
        get_new_proc = .false.

        write(*,*)' Loaded new proc...', iproc
        write(*,*)' For region...', region


    endif 

    
    write(*,*)
    write(*,*)

enddo 







end program linearly_breakup_mesh



subroutine load_new_proc(iproc, region, rem_in_proc, proc_id_ff)
    use params, only: nspec, ibool
    use mesh_utils, only: load_ibool, read_proc_coordinates
    implicit none 
    integer :: iproc, region, rem_in_proc, proc_id_ff

    call read_proc_coordinates(iproc, region)
    call load_ibool(iproc, region)

    write(*,*)'min ibool = ', minval(ibool)
    write(*,*)'max ibool = ', maxval(ibool)
    
    rem_in_proc = nspec
    proc_id_ff  = 1
end subroutine load_new_proc