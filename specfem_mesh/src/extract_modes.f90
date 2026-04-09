! Simple script to retrieve the MINEOS eigenfunctions in ASCII format 
program extract_mode

    use mineos_model, only: mineos, mineos_ptr
    use modes, only: Mode, get_mode
    
    implicit none
    include "constants.h"

    integer, dimension(33), parameter :: modeNs = (/7, 21,  7, 13, 15, 20, 25, 27, 21, 16, 3, 8,13,6,8,21,22,9,13,13,15,18,18,23,23,2,3,5,3,9,11,11,16/)
    integer, dimension(33), parameter :: modeLs = (/4,  7,  5,  6,  4,  1,  2,  2,  8,  7, 1, 1, 1,3,5, 6, 1,2, 2, 3, 3, 3, 4, 4, 5,3,2,3,8,3, 4, 5, 5/)
    integer :: imode, thisn, thisl, ir 
    character(len=1) :: thist
    type(Mode)            :: mode_1
    character(len=4) :: nstr, lstr
    character(len=400) :: out_name

    ! Read mineos model 
    call mineos%process_mineos_model(.false.)
    mineos_ptr => mineos


    do imode = 1, 33
        thisn = modeNs(imode)
        thist = 'S'
        thisl = modeLs(imode)
        write(*,*)thisn, "S", thisl

        mode_1 = get_mode(thisn, thist, thisl, mineos_ptr)

        ! For each mode save to output/eigenfunctions
        call buffer_int(nstr, thisn)
        call buffer_int(lstr, thisl)

        out_name = 'output/eigenfunctions/'//trim(nstr)//trim(thist)//trim(lstr)//".txt"


        open(1,file=trim(out_name), form='formatted')
        do ir = 1, mineos%NR
            write(1,'(5E15.7)') mineos%rad_mineos(ir), mode_1%u(ir), mode_1%v(ir), mode_1%du(ir), mode_1%dv(ir)
        enddo 
        
    enddo 

end program extract_mode