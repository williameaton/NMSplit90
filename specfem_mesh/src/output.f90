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
    use params, only: NR, IC_ID, CMB_ID, rad_mineos, radius, datadir, IIN
    implicit none 

    open(IIN, file=trim(datadir)//'/store/mineos_model/radial_data', form='unformatted')
    read(IIN)NR
    read(IIN)IC_ID
    read(IIN)CMB_ID
    allocate(radius(NR), rad_mineos(NR))
    read(IIN)rad_mineos
    read(IIN)radius
    close(IIN)
end subroutine load_mineos_radial_info


subroutine load_mineos_radial_info_MPI()
    use params, only: NR, IC_ID, CMB_ID, rad_mineos, radius, myrank, & 
                      datadir, MPI_CUSTOM_REAL
    implicit none 
#ifdef WITH_MPI
    include 'mpif.h'
    
    integer :: ierr

    if(myrank.eq.0)then 
        call load_mineos_radial_info()
    endif 

    call MPI_Bcast(NR,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(IC_ID,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(CMB_ID, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if(myrank.ne.0) allocate(radius(NR), rad_mineos(NR))

    call MPI_Bcast(radius, NR, MPI_CUSTOM_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(rad_mineos, NR, MPI_CUSTOM_REAL, 0, MPI_COMM_WORLD, ierr)
#endif 
end subroutine load_mineos_radial_info_MPI





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



subroutine load_jacobian(iproc)
    use params, only: datadir, jacinv, detjac, jac, ngllx, nglly, & 
                      ngllz, nspec, IIN
    use allocation_module, only: allocate_if_unallocated
    implicit none 
    integer :: iproc 
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    call allocate_if_unallocated(3, 3, ngllx, nglly, ngllz, nspec, jac)
    call allocate_if_unallocated(3, 3, ngllx, nglly, ngllz, nspec, jacinv)
    call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, detjac)

    ! Density
    open(IIN, file=trim(datadir)//'/store/jacobian/jacdata_'//proc_str, form='unformatted')
    read(IIN)jac
    read(IIN)detjac
    read(IIN)jacinv
    close(IIN)

end subroutine load_jacobian


subroutine save_jacobian(iproc)
    use params, only: datadir, jacinv, detjac, jac, IOUT
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    open(IOUT, file=trim(datadir)//'/store/jacobian/jacdata_'//proc_str, form='unformatted')
    write(IOUT)jac
    write(IOUT)detjac
    write(IOUT)jacinv
    close(IOUT)


end subroutine save_jacobian

subroutine save_wglljac(iproc)
    use params, only: datadir, wglljac, IOUT
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    open(IOUT, file=trim(datadir)//'/store/jacobian/wglljac_'//proc_str, form='unformatted')
    write(IOUT)wglljac
    close(IOUT)
end subroutine save_wglljac


subroutine load_wglljac(iproc)
    use params, only: datadir, ngllx, nglly, ngllz, nspec, wglljac, IIN
    use allocation_module, only: deallocate_if_allocated
    implicit none 
    integer :: iproc 
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)
    call deallocate_if_allocated(wglljac)
    allocate(wglljac(ngllx,nglly,ngllz,nspec))

    ! Density
    open(IIN, file=trim(datadir)//'/store/jacobian/wglljac_'//proc_str, form='unformatted')
    read(IIN)wglljac
    close(IIN)

end subroutine load_wglljac


subroutine save_global_xyz(iproc)
    use params, only: datadir, x_glob, y_glob, z_glob, IOUT
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    open(IOUT, file=trim(datadir)//'/store/global_xyz/coords_'//proc_str, form='unformatted')
    write(IOUT)x_glob
    write(IOUT)y_glob
    write(IOUT)z_glob
    close(IOUT)
end subroutine save_global_xyz


subroutine load_global_xyz(iproc)
    use params, only: datadir, nglob, x_glob, y_glob, z_glob, IIN
    use allocation_module, only: allocate_if_unallocated
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    call allocate_if_unallocated(nglob, x_glob)
    call allocate_if_unallocated(nglob, y_glob)
    call allocate_if_unallocated(nglob, z_glob)

    open(IIN, file=trim(datadir)//'/store/global_xyz/coords_'//proc_str, form='unformatted')
    read(IIN)x_glob
    read(IIN)y_glob
    read(IIN)z_glob
    close(IIN)
end subroutine load_global_xyz





subroutine save_elem_rtp(iproc)
    use params, only: datadir, rstore, thetastore, phistore, IOUT
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    open(IOUT, file=trim(datadir)//'/store/elemental_rtp/elem_coords_'//proc_str, form='unformatted')
    write(IOUT)rstore
    write(IOUT)thetastore
    write(IOUT)phistore
    close(IOUT)
end subroutine save_elem_rtp


subroutine load_elem_rtp(iproc)
    use params, only: datadir, rstore, thetastore, phistore,& 
                      ngllx, nglly, ngllz, nspec, IIN
    use allocation_module, only: allocate_if_unallocated
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, rstore)
    call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, thetastore)
    call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, phistore)

    open(IIN, file=trim(datadir)//'/store/elemental_rtp/elem_coords_'//proc_str, form='unformatted')
    read(IIN)rstore
    read(IIN)thetastore
    read(IIN)phistore
    close(IIN)
end subroutine load_elem_rtp

subroutine save_get_mesh_radii_results(iproc)
    use params, only: datadir, unique_r, rad_id, interp_id_r, & 
                      n_unique_rad, IOUT
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call buffer_int(proc_str, iproc)

    open(IOUT, file=trim(datadir)//'/store/mesh_radii_data/data_'//proc_str, form='unformatted')
    write(IOUT)n_unique_rad
    write(IOUT)rad_id
    write(IOUT)unique_r
    write(IOUT)interp_id_r
    close(IOUT)
end subroutine save_get_mesh_radii_results



subroutine load_get_mesh_radii_results(iproc)
    use params, only: datadir, unique_r, rad_id, interp_id_r, & 
                      ngllx, nglly, ngllz, nspec, n_unique_rad, IIN 
    use allocation_module, only: deallocate_if_allocated
    implicit none 
    integer :: iproc
    character(len = 5) proc_str

    call deallocate_if_allocated(unique_r)
    call deallocate_if_allocated(rad_id)
    call deallocate_if_allocated(interp_id_r)
    allocate(rad_id(ngllx, nglly, ngllz, nspec))

    call buffer_int(proc_str, iproc)

    open(IIN, file=trim(datadir)//'/store/mesh_radii_data/data_'//proc_str, form='unformatted')
    read(IIN)n_unique_rad
    read(IIN)rad_id

    allocate(unique_r(n_unique_rad)) 
    allocate(interp_id_r(n_unique_rad))

    read(IIN)unique_r
    read(IIN)interp_id_r
    close(IIN)
end subroutine load_get_mesh_radii_results