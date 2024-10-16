! Breaks mesh in a slow way by comparing each element. Note that the goal
! here is to be flexible to meshes produced by SPECFEM with a high number 
! of NPROC_XI such that they have many overlapping points 

program break_mesh 

    use specfem_mesh, only: SetMesh, create_setmesh
    use params, only: nprocs, datadir
    implicit none 
    include "constants.h"
    integer       :: iproc,unq_in_proc, i, j, k, ispec, region, nmaxel, ngll, n_unq_el, & 
                     minspec, iuq_spec,el_ctr, IOUT, nproc_out, iunq_ib, n_unq_ibools, nglob_max, iunq, IPROG
    type(SetMesh) :: sm

    double precision :: precision, eid
    double precision, allocatable :: unq_x(:,:,:,:), unq_y(:,:,:,:), unq_z(:,:,:,:), unq_gll_coord(:,:),  unq_eids(:)
    integer, allocatable          :: unq_ib(:,:,:,:), gll_ib(:)
    logical :: test_unique_elements, match_element, escape_element, escape_element_search, ib_match, escape_ibool_loop, match
    character(len=400) :: proc_fmtname
    precision = 1e-8
    ngll   = 5
    region = 3

    IPROG = 4

    open(IPROG,file='slurm_output.txt', form='formatted')


    nmaxel = 0 
    do iproc = 0, nprocs -1
        sm = create_setmesh(iproc, region)
        call sm%read_proc_coordinates()
        nmaxel = nmaxel + sm%nspec
        call sm%cleanup()
    end do

    write(IPROG,*)'Max. number of elements:', nmaxel

    allocate(unq_ib(ngll,ngll,ngll,nmaxel))
    allocate(unq_x(ngll,ngll,ngll,nmaxel))
    allocate(unq_y(ngll,ngll,ngll,nmaxel))
    allocate(unq_z(ngll,ngll,ngll,nmaxel))
    allocate(unq_eids(nmaxel))

    nglob_max = ngll*ngll*ngll*nmaxel

    allocate(unq_gll_coord(3, nglob_max))
    allocate(gll_ib(nglob_max))

    n_unq_el = 0 

    do iproc = 0, nprocs -1
        unq_in_proc = 0 

        ! Loop over elements
        sm = create_setmesh(iproc, region)
        call sm%read_proc_coordinates()
        call sm%compute_rtp_from_xyz(.false.)

        do ispec = 1, sm%nspec 
            match_element = .false.

            ! Compute element id: 
            eid = zero 
            do i = 1, ngll
                do j = 1, ngll
                    do k = 1, ngll
                        eid = eid + &
                                sm%xstore(i,j,k,ispec)/four + sm%ystore(i,j,k,ispec) + sm%zstore(i,j,k,ispec) & 
                                + sm%rstore(i,j,k,ispec) + sm%thetastore(i,j,k,ispec)/PI
                    enddo 
                enddo 
            enddo 

            ! Loop over the number of unique elements so far: 
            if(n_unq_el.gt.0)then
                do iuq_spec = 1, n_unq_el
                    if ( abs(eid - unq_eids(iuq_spec)).lt.precision )then 
                        match_element         = .true.
                    endif 
                enddo !iuq_spec
            endif 

            ! If matching another element then we simply want to ignore this extra
            ! copy of it - otherwise we need to work out which of the GLL points on 
            ! this new element already have an ibool and which dont 
            if(match_element.ne..true.)then 
                n_unq_el              = n_unq_el + 1
                unq_x(:,:,:,n_unq_el) = sm%xstore(:,:,:,ispec)
                unq_y(:,:,:,n_unq_el) = sm%ystore(:,:,:,ispec)
                unq_z(:,:,:,n_unq_el) = sm%zstore(:,:,:,ispec)
                unq_eids(n_unq_el)    = eid
                unq_in_proc = unq_in_proc + 1
            endif ! it not matching another element 

            !write(*,*)'eid for this element', eid
            !write(*,*)'unq_eids(n_unq_el)', unq_eids(1:n_unq_el)


        enddo !ispec
        write(IPROG,*)'done proc: ', iproc, 'with unique els: ', unq_in_proc
        call sm%cleanup()
    end do ! iproc


    flush(IPROG)


    write(*,*)'N unique elements    : ', n_unq_el
    n_unq_ibools = 0 
    ! We should now have a set of unqiue elements
    ! Next we want to get their ibool 
    do ispec = 1, n_unq_el
        do i = 1, ngll 
            do j = 1, ngll 
                do k = 1, ngll 

                    ! for each GLL point - loop through the current unique ibools
                    match = .false.
                    if(n_unq_ibools.gt.0)then 
                        do iunq = 1, n_unq_ibools
                            if( abs(unq_x(i,j,k,ispec)-unq_gll_coord(1,iunq)).gt.precision .or. & 
                                abs(unq_y(i,j,k,ispec)-unq_gll_coord(2,iunq)).gt.precision .or. & 
                                abs(unq_z(i,j,k,ispec)-unq_gll_coord(3,iunq)).gt.precision) then 
                                    ! Doesnt match this 
                            else 
                                ! Matches this coordinate
                                unq_ib(i,j,k,ispec) = gll_ib(iunq)
                                match = .true.
                            endif

                            if(match)exit
                        enddo !iunq
                    endif

                    if(.not.match)then 
                        n_unq_ibools = n_unq_ibools + 1
                        unq_gll_coord(1,n_unq_ibools) = unq_x(i,j,k,ispec)
                        unq_gll_coord(2,n_unq_ibools) = unq_y(i,j,k,ispec)
                        unq_gll_coord(3,n_unq_ibools) = unq_z(i,j,k,ispec)
                        gll_ib(n_unq_ibools) = n_unq_ibools
                        unq_ib(i,j,k,ispec)  = n_unq_ibools
                    endif 

                enddo 
            enddo 
        enddo 

        if(mod(ispec,1000).eq.0)then
             write(IPROG,*)'ispec done: ', ispec 
             flush(IPROG)
        endif 
    enddo

    write(IPROG,*)'n_unq_ibools = ', n_unq_ibools 
    write(IPROG,*)'min of unq_ib = ', minval(unq_ib) 
    write(IPROG,*)'max of unq_ib = ', maxval(unq_ib)
    
    write(IPROG,*)'min of unq_ib in range = ', minval(unq_ib(:,:,:,1:n_unq_el)) 
    write(IPROG,*)'max of unq_ib in range = ', maxval(unq_ib(:,:,:,1:n_unq_el))



    ! Write it out: 
    iproc=0
    IOUT=1

    write(proc_fmtname,'(a,i0.6,a,i1,a)')trim(datadir)//'/sliced/'//'/proc',iproc,'_'//'reg',region,'_'

    open(unit=IOUT,file=trim(proc_fmtname)//'xstore.bin', &
    status='unknown',form='unformatted',action='write')
    write(IOUT) unq_x(:,:,:,1:n_unq_el)
    close(IOUT)

    open(unit=IOUT,file=trim(proc_fmtname)//'ystore.bin', &
    status='unknown',form='unformatted',action='write')
    write(IOUT) unq_y(:,:,:,1:n_unq_el)
    close(IOUT)

    open(unit=IOUT,file=trim(proc_fmtname)//'zstore.bin', &
    status='unknown',form='unformatted',action='write')
    write(IOUT) unq_z(:,:,:,1:n_unq_el)
    close(IOUT)
    
    open(unit=IOUT,file=trim(proc_fmtname)//'ibool.bin', &
    status='unknown',form='unformatted',action='write')
    write(IOUT) unq_ib(:,:,:,1:n_unq_el)
    close(IOUT)

    open(unit=IOUT,file=trim(proc_fmtname)//'info.bin', &
    status='unknown',form='unformatted',action='write')
    write(IOUT) maxval(unq_ib) ! new nglob
    write(IOUT) n_unq_el
    write(IOUT) ngll
    write(IOUT) ngll
    write(IOUT) ngll
    close(IOUT)


    ! Determine the possible number of procs: 
    do i = 1, n_unq_el
        if (mod(n_unq_el,i).eq.0)then
            write(*,*)i
        endif
    enddo 



end program break_mesh