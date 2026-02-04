


subroutine load_ACLNF_from_files(directory, nlen, suffix, myrank)
    use params, only: Arad, Crad, Lrad, Nrad, Frad, cluster_size
    use allocation_module, only: deallocate_if_allocated
    implicit none 
    ! IO variables
    character(len=*) :: directory
    integer ::  nlen, myrank

    character(len=5) :: suffix
    ! Local variables
    integer :: ios, i, icl
    character(len=300) :: fname
    character(len=5), parameter :: ACLNF='ACLNF'

    write(*,*)'Loading ACLNF from files...'

    call deallocate_if_allocated(Arad)
    allocate(Arad(nlen))
    call deallocate_if_allocated(Crad)
    allocate(Crad(nlen))
    call deallocate_if_allocated(Lrad)
    allocate(Lrad(nlen))
    call deallocate_if_allocated(Nrad)
    allocate(Nrad(nlen))
    call deallocate_if_allocated(Frad)
    allocate(Frad(nlen))


    if(cluster_size.eq.0)cluster_size = 1

    ! Avoid reading issues with MPI? 
    do icl = 0, cluster_size-1 
        if(myrank.eq.icl)then 
                
            do i = 1, 5
                write(fname,'(a)')trim(directory)//'/'//ACLNF(i:i)//suffix
                open(i, file=trim(fname), status = 'old', iostat = ios)
                if(ios.ne.0)then
                    write(*,*)'Error opening ', trim(fname)
                    stop 
                endif
            enddo 

            do i = 1, nlen
                read(1,*)Arad(i)
                read(2,*)Crad(i)
                read(3,*)Lrad(i)
                read(4,*)Nrad(i)
                read(5,*)Frad(i)
            enddo 
        endif 

    enddo 
    ! Reset 
    cluster_size = 0

end subroutine load_ACLNF_from_files