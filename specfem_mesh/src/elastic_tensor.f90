


subroutine load_ACLNF_from_files(directory, iproc, nlen)
    use params, only: Arad, Crad, Lrad, Nrad, Frad
    use allocation_module, only: deallocate_if_allocated
    implicit none 
    ! IO variables
    character(len=*) :: directory
    integer :: iproc, nlen
    
    ! Local variables
    integer :: ios, i
    character(len=300) :: fname
    character(len=4) :: proc_str
    character(len=5), parameter :: ACLNF='ACLNF'

    call buffer_int(proc_str, iproc)

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

    ! Read in the A values 
    do i = 1, 5
        write(fname,'(a)')trim(directory)//'/'//ACLNF(i:i)//'_'//trim(proc_str)
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



end subroutine load_ACLNF_from_files