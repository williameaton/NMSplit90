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


subroutine save_mineos_model()
    use params, only: NR, IC_ID, CMB_ID, rad_mineos, radius, rho_mineos, & 
    vp, disc, rdisc, ndisc, datadir
    implicit none 
    integer :: iproc

    ! Radius and knot info
    open(1, file=trim(datadir)//'/store/mineos_model/radial_data', form='unformatted')
    write(1)NR
    write(1)IC_ID
    write(1)CMB_ID
    write(1)rad_mineos
    write(1)radius
    close(1)

    ! Density
    open(1, file=trim(datadir)//'/store/mineos_model/rho_mineos', form='unformatted')
    write(1)rho_mineos
    close(1)

    ! Vp 
    open(1, file=trim(datadir)//'/store/mineos_model/vp_mineos', form='unformatted')
    write(1)vp
    close(1)

    ! Discontinuity info
    open(1, file=trim(datadir)//'/store/mineos_model/disc', form='unformatted')
    write(1)ndisc
    write(1)rdisc
    write(1)disc
    close(1)

end subroutine save_mineos_model


subroutine load_mineos_radial_info()
    use params, only: NR, IC_ID, CMB_ID, rad_mineos, radius, datadir
    implicit none 

    open(1, file=trim(datadir)//'/store/mineos_model/radial_data', form='unformatted')
    read(1)NR
    read(1)IC_ID
    read(1)CMB_ID

    allocate(radius(NR), rad_mineos(NR))

    read(1)rad_mineos
    read(1)radius
    close(1)

end subroutine load_mineos_radial_info


subroutine load_mineos_density()
    use params, only: vp, datadir, NR
    implicit none 

    allocate(vp(NR))

    open(1, file=trim(datadir)//'/store/mineos_model/rho_mineos', form='unformatted')
    read(1)vp
    close(1)
end subroutine load_mineos_density


subroutine load_mineos_vp()
    use params, only: rho_mineos, datadir, NR
    implicit none 

    allocate(rho_mineos(NR))

    open(1, file=trim(datadir)//'/store/mineos_model/vp_mineos', form='unformatted')
    read(1)rho_mineos
    close(1)
end subroutine load_mineos_vp



subroutine load_mineos_discontinuities()
    use params, only: ndisc, rdisc, disc, datadir
    implicit none 
    open(1, file=trim(datadir)//'/store/mineos_model/disc', form='unformatted')
    read(1)ndisc
    read(1)rdisc
    read(1)disc
    close(1)
end subroutine load_mineos_discontinuities




subroutine save_mode_strain_binary(n, type, l, m, strain, iproc)
    use params, only: ngllx, nglly, ngllz, nspec, datadir
    implicit none 
    include "constants.h"

    integer :: n, l, m, iproc
    character :: type 
    complex(kind=SPLINE_REAL) :: strain(6, ngllx, nglly, ngllz, nspec)

    character(len = 250) binname
    
    call create_mode_binary_fname(n, type, l, m, iproc, binname)


    open(1, file=trim(datadir)//'/store/strain/'//trim(binname), form='unformatted')
    write(1)strain
    close(1)
end subroutine save_mode_strain_binary


subroutine load_mode_strain_binary(n, type, l, m, strain, iproc)
    
    use params, only: ngllx, nglly, ngllz, nspec, datadir
    implicit none 
    include "constants.h"

    integer :: n, l, m, iproc, ios
    character :: type 
    complex(kind=SPLINE_REAL) :: strain(6, ngllx, nglly, ngllz, nspec)

    character(len = 250) binname
    
    call create_mode_binary_fname(n, type, l, m, iproc, binname)

    open(1, file=trim(datadir)//'/store/strain/'//trim(binname), form='unformatted', iostat=ios)
    if(ios.ne.0)then
        write(*,*)'Could not open mode', trim(binname)
        stop
    endif 
    read(1)strain
    close(1)

end subroutine load_mode_strain_binary


subroutine save_mode_disp_binary(n, type, l, m, disp, iproc)
    
    use params, only: ngllx, nglly, ngllz, nspec, datadir
    implicit none 
    include "constants.h"

    integer :: n, l, m, iproc
    character :: type 
    complex(kind=SPLINE_REAL) :: disp(3, ngllx, nglly, ngllz, nspec)

    character(len = 250) binname
    
    call create_mode_binary_fname(n, type, l, m, iproc, binname)


    open(1, file=trim(datadir)//'/store/disp/'//trim(binname), form='unformatted')
    write(1)disp
    close(1)

end subroutine


subroutine load_mode_disp_binary(n, type, l, m, disp, iproc)
    
    use params, only: ngllx, nglly, ngllz, nspec, datadir
    implicit none 
    include "constants.h"

    integer :: n, l, m, iproc
    character :: type 
    complex(kind=SPLINE_REAL) :: disp(3, ngllx, nglly, ngllz, nspec)

    character(len = 250) binname
    
    call create_mode_binary_fname(n, type, l, m, iproc, binname)

    open(1, file=trim(datadir)//'/store/disp/'//trim(binname), form='unformatted')
    read(1)disp
    close(1)

end subroutine