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
use params, only: datadir
use specfem_mesh, only: SetMesh, create_setmesh
implicit none 

integer :: region, nprocs_before, nsets, iproc, ntot_elem,   & 
           nspec_per_set, rem_el_in_set, rem_in_proc, & 
           nelms_to_copy, set_id_b, set_id_ff,        & 
           set_id_end, proc_id_ff, proc_id_end, iset, & 
           ib_orig, maxnglob, ib_proc, min_ib_proc,   & 
           max_ib_proc, min_ib, max_ib, nibool, iiprc,&
            ispec, i,j,k,iib, count, ntot_top, ntot_bottom, & 
            itopstart, topidx, ielemtop, ibl, itop, ielembottom, bottomidx, & 
            ibottom
double precision, allocatable :: xcoord_new(:,:,:,:), & 
                                 ycoord_new(:,:,:,:), & 
                                 zcoord_new(:,:,:,:)
integer, allocatable :: ibool_new(:,:,:,:), ib_store(:,:,:,:), & 
                        ib_store_proc(:,:,:,:), ib_map(:), ib_top(:),& 
                        ib_bottom(:), ib_bottomproc(:), ib_topproc(:), new_ib_top(:), new_ib_bottom(:)
logical :: get_new_proc, filled_set
character(len=400) :: outname
character(len=30) :: fmtstr
type(SetMesh) :: sm 

integer :: ibsctr, igllx, iglly, isetib, iblspec

logical, parameter :: include_boundaries = .false.

! Setup parameters: 
region        = 3      ! CM
nprocs_before = 6      ! Current setup 
nsets         = 4      ! new setup 

ntot_elem   = 0 ! Total elements
ntot_top    = 0 ! Top boundary
ntot_bottom = 0 ! Bottom boundary

! Count the total number of elements
do iproc = 0, nprocs_before - 1 
    sm = create_setmesh(iproc, region)
    call sm%read_proc_coordinates()  
    ntot_elem    = ntot_elem   + sm%nspec
    ntot_top     = ntot_top    + sm%nspec2D_top     ! Should be same on all procs originally
    ntot_bottom  = ntot_bottom + sm%nspec2D_bottom
    call sm%cleanup()
enddo 


! TODO: Could be more flexible to cases without absolutley perfect division
if(mod(ntot_elem,nsets).ne.0)then 
    write(*,*)'Error: total number of elements is not divisible by nsets'
    write(*,*)'Total number of elements: ', ntot_elem
    write(*,*)'Requested sets: ',nsets
    do i = 1, nsets-1
        if (mod(ntot_elem,i).eq.0)then 
            write(*,*)'Nets = : ',i
        endif 
    enddo 
    stop 
endif 

nspec_per_set = ntot_elem/nsets
write(*,*)' Total nspec              : ', ntot_elem
write(*,*)' Total on top boundary    : ', ntot_top
write(*,*)' Total on bottom boundary : ', ntot_bottom
write(*,*)' Number of set            : ', nsets
write(*,*)' Elements per set         : ', nspec_per_set


! Allocate arrays of correct size for a single set: 
allocate(xcoord_new(sm%ngllx, sm%nglly, sm%ngllz, nspec_per_set)) 
allocate(ycoord_new(sm%ngllx, sm%nglly, sm%ngllz, nspec_per_set)) 
allocate(zcoord_new(sm%ngllx, sm%nglly, sm%ngllz, nspec_per_set)) 
allocate(ibool_new(sm%ngllx, sm%nglly, sm%ngllz, nspec_per_set)) 
allocate(ib_store(sm%ngllx, sm%nglly, sm%ngllz,   nspec_per_set)) 
allocate(ib_store_proc(sm%ngllx, sm%nglly, sm%ngllz,nspec_per_set)) 

! There can now be an unequal number of boundary elements in each set
! Because of this, we will play it safe because we know there cant be more
! Elements with a boundary than elements in the set
if(include_boundaries)then 
    allocate(ib_top(nspec_per_set)) 
    allocate(ib_topproc(nspec_per_set)) 
    allocate(ib_bottom(nspec_per_set)) 
    allocate(ib_bottomproc(nspec_per_set)) 
    ib_top    = -1 
    ib_bottom = -1 
endif


! Will set remaining_in_processor to the nspec of the proc
! Will set processor_id_fill_from to 1 
iproc = 0
call load_new_proc(sm, iproc, region, rem_in_proc, proc_id_ff)
iset  = 0 
rem_el_in_set = nspec_per_set   ! Remaining elements in a set
set_id_ff     = 1               ! ID to fill from in set
topidx = 0


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
    xcoord_new(:,:,:,set_id_ff:set_id_end) = sm%xstore(:,:,:,proc_id_ff:proc_id_end)
    ycoord_new(:,:,:,set_id_ff:set_id_end) = sm%ystore(:,:,:,proc_id_ff:proc_id_end)
    zcoord_new(:,:,:,set_id_ff:set_id_end) = sm%zstore(:,:,:,proc_id_ff:proc_id_end)

    ! Store the ibool values AND the proc it was on 
    ib_store(:,:,:,set_id_ff:set_id_end)       = sm%ibool(:,:,:,proc_id_ff:proc_id_end)
    ib_store_proc(:,:,:,set_id_ff:set_id_end)  = iproc


    ! Get the min and max ibool values for this set: 
    !minsetib = minval(ib_store(:,:,:,set_id_ff:set_id_end))
    !maxsetib = maxval(ib_store(:,:,:,set_id_ff:set_id_end))

    ! Copy over the boundary ibools: 
    ! Loop through the range of ibool values that are permissible based on the elements that
    ! are added 
    do ielemtop = 1, sm%nspec2D_top
        iblspec = sm%ibelm_top(ielemtop) ! ispec of surface element
        if(iblspec.ge.proc_id_ff .and.iblspec.le.proc_id_end)then 
            topidx = topidx + 1
            ib_top(topidx)      = iblspec-proc_id_ff + set_id_ff
            !ib_topproc(topidx)  = iproc
        endif 
    enddo


    do ielembottom = 1, sm%nspec2D_bottom
        iblspec = sm%ibelm_bottom(ielembottom) ! ispec of surface element
        if(iblspec.ge.proc_id_ff .and.iblspec.le.proc_id_end)then 
            bottomidx = bottomidx + 1
            ib_bottom(bottomidx)   = iblspec-proc_id_ff + set_id_ff
        endif 
    enddo


    ! Update where to start filling the new set from/where to start getting data from on proc
    set_id_ff  = set_id_end  + 1 
    proc_id_ff = proc_id_end + 1 

    rem_el_in_set = rem_el_in_set - nelms_to_copy
    rem_in_proc   = rem_in_proc   - nelms_to_copy

    write(*,*)' New ID to fill set from  = ', set_id_ff



    if(set_id_end.eq.nspec_per_set)then 

        write(*,*)' Set is now filled...'
        ! Compute ibool for this set & output
        min_ib_proc = minval(ib_store_proc) ! First processor ID that ibool is taken from
        max_ib_proc = maxval(ib_store_proc) ! Last  processor ID that ibool is taken from
        min_ib = minval(ib_store)   ! Minimum ibool value
        max_ib = maxval(ib_store)   ! Maximum ibool value


        !write(*,*)'Topidx is now: ', topidx
        !write(*,*)'Bottomidx is now: ', bottomidx
        !write(*,*)'allocating with size: ', topidx, bottomidx
        !if (topidx.gt.0)allocate(new_ib_top(topidx))
        !if(bottomidx.gt.0)allocate(new_ib_bottom(bottomidx))


        write(*,*,advance='no')' Reformatting ibool...'
        nibool = 1

        ! There is the possibility that elements originally on 2 different
        ! processors have been combined into the same set and that these 
        ! elements have the same ibool value even though they repreent different
        ! gll points. I.e. ibool=1000 from proc0 and proc1 would be different gll points
        ! hence we are not just looking for where it equals the ibool but also where the 
        ! proc was originally the same.
        do iiprc = min_ib_proc, max_ib_proc
            do iib = min_ib, max_ib

                ! Loop over every gll in the mesh
                ! for any gll points that are equal to the procID (iiprc) AND ibool (iib) value then we set their new ibool values
                ! to be the same/shared with value nibool 
                ! Any time that we loop through and find that we have added to the new ibool, we then update nibool 
                ! so that a NEW, unique ibool number is generated for the next case
                count = 0
                do ispec = 1, nspec_per_set
                    do i = 1, sm%ngllx
                        do j = 1, sm%nglly
                            do k = 1, sm%ngllz
                                if(ib_store(i,j,k,ispec).eq.iib .and. ib_store_proc(i,j,k,ispec).eq.iiprc)then
                                    ibool_new(i,j,k,ispec) = nibool 
                                    count = count + 1
                                endif
                            enddo 
                        enddo 
                    enddo
                enddo 
                                    
                if (count.gt.0)nibool = nibool + 1
                    ! ! We know that in this case iiprc and iib relates to gll points
                    ! ! i.e. we know that 'iib' and 'iiprc' maps to nibool
                    ! if(topidx.gt.0)then
                    !     do itop = 1, topidx
                    !         if(ib_top(itop).eq.iib .and. ib_topproc(itop).eq.iiprc)then
                    !             new_ib_top(itop) = nibool
                    !         endif 
                    !     enddo 
                    ! endif 
                    ! if(bottomidx.gt.0)then 
                    !     do ibottom = 1, bottomidx
                    !         if(ib_bottom(ibottom).eq.iib .and. ib_bottomproc(ibottom).eq.iiprc)then
                    !             new_ib_bottom(ibottom) = nibool
                    !         endif 
                    !     enddo 
                    ! endif

            enddo !iib
        enddo !iiproc
        write(*,*)' done'



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

        write(outname,trim(fmtstr))trim(datadir)//'/linear/sets',nsets,'/proc',iset,'_'//'reg',region,'_'

        open(1,file=trim(outname)//'xstore.bin', form='UNFORMATTED')
        write(1)xcoord_new
        close(1)

        open(2,file=trim(outname)//'ystore.bin', form='UNFORMATTED')
        write(2)ycoord_new
        close(2)

        open(3,file=trim(outname)//'zstore.bin', form='UNFORMATTED')
        write(3)zcoord_new
        close(3)

        open(4,file=trim(outname)//'ibool.bin', form='UNFORMATTED')
        write(4)ibool_new
        close(4)


        open(5,file=trim(outname)//'info.bin', form='UNFORMATTED')
        write(5)maxval(ibool_new)
        write(5) nspec_per_set
        write(5) sm%NGLLX
        write(5) sm%NGLLY
        write(5) sm%NGLLZ
        close(5)

        ! Top boundary
        !print *, "topidx =", topidx
        !print *, "size of new_ib_top =", size(new_ib_top)
       ! Top boundary
        !print *, "bottomidx =", bottomidx
        !print *, "size of bottomidx =", size(new_ib_bottom)


        ! DO not use id of 6 -- this causes a malloc error
        open(8,file=trim(outname)//'ibool_top.bin', form='UNFORMATTED')
        write(8) topidx ! number of nspectop on this set
        if(topidx.gt.0)then
            write(8) ib_top(1:topidx)
            write(*,*)'topset: ', topidx
            write(*,*)'topset vals: ', ib_top(1:topidx)
        endif 
        close(8)
 

        ! Bottom boundary
        open(7,file=trim(outname)//'ibool_bottom.bin', form='UNFORMATTED')
        write(7) bottomidx ! number of nspec bottom on this set
        if(bottomidx.gt.0)then
            write(7) ib_bottom(1:bottomidx)
            write(*,*)'bottom vals: ', ib_bottom(1:bottomidx)
        endif 
        close(7)


        ! Reset the set parameters : 
        write(*,*)'  --> resets set_id_ff to 1'

        iset = iset + 1
        rem_el_in_set = nspec_per_set   ! Remaining elements in a set
        set_id_ff     = 1               ! ID to fill from in set


        ! Reset boundary ibools: 
        write(*,*)' deallocating'

        !if(allocated(new_ib_top))deallocate(new_ib_top)
        !if(allocated(new_ib_bottom))deallocate(new_ib_bottom)
        write(*,*)' done'

        ib_top    = -1 
        ib_bottom = -1 
        topidx    = 0
        bottomidx = 0


        !write(*,*)'  --> sets rem_el_in_set to ', nspec_per_set
        !write(*,*)'  --> iset is now  ', iset
    endif 


    if(get_new_proc)then 
        write(*,*)' Loading new proc...'
        call sm%cleanup()
        iproc = iproc + 1
        call load_new_proc(sm, iproc, region, rem_in_proc, proc_id_ff)
        get_new_proc = .false.
    endif 

    
    write(*,*)

enddo 


end program linearly_breakup_mesh



subroutine load_new_proc(sm, iproc, region, rem_in_proc, proc_id_ff)
    use specfem_mesh, only: SetMesh, create_setmesh
    implicit none 
    integer :: iproc, region, rem_in_proc, proc_id_ff
    type(SetMesh) :: sm 
    sm = create_setmesh(iproc, region)
    call sm%read_proc_coordinates()
    call sm%load_ibool()
    !call sm%load_original_boundaries()
    
    rem_in_proc = sm%nspec
    proc_id_ff  = 1
end subroutine load_new_proc