program extract_mode
    use params, only: NL 

    implicit none
    include "constants.h"

    integer :: i, j, l, n , lentrim, lenl, lenn
    character(len=1)   :: char, type
    character(len=10)  :: mode
    character(len=250) :: arg
    character(len=4), parameter :: fmti1 = "(i1)"
    character(len=4), parameter :: fmti2 = "(i2)"
    character(len=4), parameter :: fmti3 = "(i3)"
    character(len=4) :: fmt_n, fmt_l
    real(SPLINE_REAL) :: wcom, qmod
    real(SPLINE_REAL), allocatable :: u(:), v(:), du(:), dv(:)
    
    logical :: found



    ! Read mineos model 
    call process_mineos_model()

    allocate(u(NL), v(NL), du(NL), dv(NL))

    i = 1
    do
      call get_command_argument(i, arg)
      if (len_trim(arg) == 0) exit
      mode=trim(arg)
      
      lentrim = len_trim(mode) 

      found = .false. 
      ! Process the mode string 
      do j = 1, lentrim
        char = mode(j:j)
        if (char.eq.'T' .or. char.eq.'S' .or. char.eq.'C') then 
            type = char
            found = .true.

            lenn = j-1 
            lenl = lentrim - j

            SELECT CASE (lenl)
            CASE (1)
                fmt_l = fmti1
            CASE (2)
                fmt_l = fmti2
            CASE (3)
                fmt_l = fmti3
            CASE DEFAULT
               write(*,*)'Error: could not process mode ', trim(mode), ' should be format nSl or nTl'
               write(*,*)'Error, limit of l is 999.'; stop
            END SELECT

            SELECT CASE (lenn)
            CASE (1)
                fmt_n = fmti1
            CASE (2)
                fmt_n = fmti2
            CASE (3)
                fmt_n = fmti3
            CASE DEFAULT
               write(*,*)'Error: could not process mode ', trim(mode), ' should be format nSl or nTl'
               write(*,*)'Error, limit of n is 999.'; stop
            END SELECT
            
            read(mode(1:lenn), fmt_n) n
            read(mode(j+1:lentrim),fmt_l) l
            write(*,*)'n is', n
            write(*,*)'l is', l
        endif
      enddo 

      if(.not.found)then
        write(*,*)'Error: could not process mode ', trim(mode), ' should be format nSl or nTl'
      endif 
    

      ! Now we can load the mode: 
      call get_mode(type,n,l,wcom,qmod, u, du, v, dv, .true.)


      
      i = i+1
    enddo

    write(*,'(a, i2, a)')'Got ', i-1, ' modes'



end program extract_mode