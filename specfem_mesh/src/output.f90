subroutine format_mode_string(n, t, l, outstr)

    implicit none 

    integer   :: n 
    integer   :: l
    character :: t
    character(len=*) :: outstr

    character(len=2), parameter :: fmti1 = "i1"
    character(len=2), parameter :: fmti2 = "i2"
    character(len=2), parameter :: fmti3 = "i3"
    character(len=2) :: fmt_n, fmt_l
    character(len=9) :: fmtstr


    if(n.ge.0 .and. n.lt.10)then 
        fmt_n = fmti1
    elseif (n.ge.10 .and. n.lt.100)then
        fmt_n = fmti2
    elseif (n.ge.100 .and. n.lt.1000)then 
        fmt_n = fmti3
    else 
    write(*,*)'Error, limit of n is 999.'; stop
    endif

    if(l.ge.0 .and. l.lt.10)then 
        fmt_l = fmti1
    elseif (l.ge.10 .and. l.lt.100)then
        fmt_l = fmti2
    elseif (l.ge.100 .and. l.lt.1000)then 
        fmt_l = fmti3
    else 
    write(*,*)'Error, limit of l is 999.'; stop
    endif

    fmtstr = '('//fmt_n//',a,'//fmt_l//')'

    write(outstr, fmtstr)n, t, l 
end subroutine



subroutine buffer_int(str, myint)
    implicit none 

    character(len=*):: str
    integer :: myint

    if (myint.ge.0 .and. myint.lt.10)then
        write(str,'(i1)')myint
    elseif(myint.ge.10 .and. myint.lt.100)then 
        write(str,'(i2)')myint
    elseif(myint.ge.100 .and. myint.lt.1000)then 
        write(str,'(i3)')myint    
    elseif(myint.ge.1000 .and. myint.lt.10000)then 
        write(str,'(i5)')myint  


    elseif (myint.gt.-10 .and. myint.lt.0)then
            write(str,'(i2)')myint
    elseif (myint.gt.-100 .and. myint.lt.10)then
        write(str,'(i3)')myint
    elseif (myint.gt.-1000 .and. myint.lt.100)then
        write(str,'(i4)')myint
    endif

end subroutine buffer_int



subroutine create_mode_binary_fname(n, t, l, m, iproc, fname)
    ! Creates string used in reading/writing mode displacement binaries
    ! format e.g. 0S4_3_proc5.bin 
    implicit none 
    ! IO variables 
    integer :: n, l, m 
    character :: t 
    integer :: iproc
    character(len=*) :: fname

    ! Local variables:
    character(len = 5) proc_str, m_str
    character(len = 20) mode_label

    call buffer_int(proc_str, iproc)
    call buffer_int(m_str, m)
    call format_mode_string(n, t, l, mode_label)

    write(fname, '(a)') trim(mode_label)//'_'//trim(m_str)//'_proc'//trim(proc_str)//'.bin'
end subroutine create_mode_binary_fname


subroutine create_rhospline_fname(iproc, fname)
    ! Creates string used in reading/writing density spline binaries
    ! format e.g. rhospline_4.bin
    implicit none 
    
    ! IO variables 
    integer :: iproc
    character(len=*) :: fname

    ! Local variables:
    character(len = 5) proc_str, m_str
    character(len = 20) mode_label
    call buffer_int(proc_str, iproc)

    write(fname,'(a)')'rhospline_'// trim(proc_str)// '.bin'
end subroutine create_rhospline_fname


subroutine save_rhospline_binary(rvals, spline, length, iproc)
    use params, only: datadir
    implicit none
    include "constants.h"

    integer :: length, iproc
    real(kind=CUSTOM_REAL):: rvals(length)
    real(kind=SPLINE_REAL):: spline(length)

    character(len = 5) :: iproc_str
    character(len=250) :: out_name

    call create_rhospline_fname(iproc, out_name)
    
    open(1, file=trim(datadir)//'/store/rho/'//trim(out_name), form='unformatted')
    write(1)length
    write(1)rvals
    write(1)spline
    close(1)
end subroutine save_rhospline_binary
