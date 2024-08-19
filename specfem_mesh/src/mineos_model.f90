subroutine process_mineos_model()

    use params, only: rad_mineos, NR, NL, ddir, model_fname, RA, radius, &
                      verbose, IC_ID, CMB_ID

    implicit none
    include "constants.h"

    ! Local variables: 
    character(len=200) :: model_file
    integer :: iomod, ios, intjunk, i 
    real(kind=CUSTOM_REAL) :: realjunk, rho, vpv, vph, vsv, vsh, vsv_prev

    iomod = 1

    ! read model file to get radius array (normalised 0 to 1)
    model_file = trim(ddir)//trim(model_fname)
    
    if(verbose.ge.2) write(*,'(/,a)')'- Reading Model file ....'

    open(iomod, file = model_file, status = 'old', iostat = ios)
    ! Error check
    if (ios .ne. 0) stop 'Error reading model file'

    NR = -1
    do while (ios.eq.0)
        read(iomod, *, iostat=ios) intjunk, realjunk
        NR = NR + 1
    enddo 
    close(1)


    ! Setup model params based on this first read
    NL = NR
    RA = realjunk
    allocate(radius(NR), rad_mineos(NR))

    ! Warning if not radius 6371 km 
    if (RA.ne.SCALE_R .and. verbose.ge.1)then 
        write(*,*)'Warning: Radius is not 6371 km but that is being used for general non-dimensionalisation scaling'
        write(*,*)'         Maximum value may of non-dimensional rad_mineos may not be 1'
    endif 



    ! Second read of the file: 
    open(iomod, file = model_file, status = 'old', iostat = ios)

    vsv_prev = 99.99d0
    do i = 1 , NR
        read(iomod,*,iostat=ios)  intjunk, radius(i), rho, vpv, vph, vsv, vsh

        rad_mineos(i) = radius(i) / SCALE_R
        
        ! Detect ICB 
        if(vsv_prev.gt.ZERO .and. vsv.eq.ZERO)then
            IC_ID = i - 1 
        endif
        ! Detect ICB 
        if(vsv_prev.eq.ZERO .and. vsv.gt.ZERO)then
            CMB_ID = i - 1 
        endif

        vsv_prev = vsv
    enddo
    close(iomod)



    ! Find the discontinuities
    call find_disc

    if(verbose.ge.2) write(*,'(a,/)')'--> finished reading model file.'
end subroutine process_mineos_model




subroutine find_disc()
    ! Modified version of find_disc in get_mode
    use params, only: NR, disc, rdisc, ndisc, radius
    implicit none
    include "precision.h"
    
    real(8) :: eps, radius_last
    integer :: nlayer,jdisc,iomod,i
    
    ! Temporary discontinuity arrays
    integer :: disctmp(NR)
    real(kind=CUSTOM_REAL) :: rdisctmp(NR)

    radius_last = -1.0
    jdisc = 0
    iomod = 11
 
 
    do i = 1, NR
       if (abs(radius(i)-radius_last) .lt. 1.0d-6) then
          jdisc = jdisc + 1
          rdisctmp(jdisc) = radius(i)
          disctmp(jdisc) = i - 1
       endif
       radius_last = radius(i)
    enddo
    close(iomod)
    
    ! Allocate arrays of the length ndisc
    ndisc = jdisc
    allocate(rdisc(ndisc), disc(ndisc))
    rdisc(1:ndisc) = rdisctmp(1:ndisc)
    disc(1:ndisc) = disctmp(1:ndisc)

    
end subroutine find_disc
 