program read_specfem_mesh
use params 
use allocation_module

implicit none 

integer :: iproc, i                  
integer :: nprocs                ! Number of processors used by mesher 
integer :: region                ! Region code
character(len=250) :: varname

integer :: n, l, m 
! Setup parameters: 
region = 3      ! Inner core
iproc  = 1
nprocs = 6

n = 10
l = 5
m = 2


! Read mineos model 
call process_mineos_model()

iproc = 0

! Read the mesh info and coordinates
call read_proc_coordinates(iproc, region)

! Load density variable: 
!allocate(rho(ngllx,nglly,ngllz,nspec))
!varname = 'rho'
!call read_proc_variable(iproc, region, rho, varname)

! Load ibool variable: 
call load_ibool(iproc, region)




! Get unique mesh radii that are present
call get_mesh_radii()

! We will need the r theta phi coordinates: 
call compute_rtp_from_xyz()

! Compute strain tensor for mode
!allocate(strain1(6,ngllx, nglly, ngllz, nspec), & 
!        globalstrain(6, nglob))
!call compute_gll_mode_strain('S', n, l, m, strain1)
!call map_complex_vector(6, strain1, globalstrain, 0)

! needed for ensight geo file
call allocate_if_unallocated(nglob, x_glob)
call allocate_if_unallocated(nglob, y_glob)
call allocate_if_unallocated(nglob, z_glob)
call map_local_global_double_precision(xstore, x_glob, 0)
call map_local_global_double_precision(ystore, y_glob, 0)
call map_local_global_double_precision(zstore, z_glob, 0)

call compute_rotation_matrix()

call allocate_if_unallocated(nglob, theta_glob)
call allocate_if_unallocated(nglob, phi_glob)
call map_local_global_double_precision(thetastore, theta_glob, 0)
call map_local_global_double_precision(phistore,   phi_glob,   0)


call create_ensight_file_prefix(iproc, region)
call create_proc_case_file()
call create_proc_geo_file(iproc, region, 1)
!call write_complex_symtensor_to_ensight(globalstrain, 'strain', 1)

! Compute global mode displacement 
call allocate_if_unallocated(3, nglob, globaldisp)
call allocate_if_unallocated(3, ngllx, nglly, ngllz, nspec, disp1)

call compute_global_mode_displacement('S', n, l, m, disp1)
call map_complex_vector(3, disp1, globaldisp, 0)
call write_complex_vector_to_ensight(globaldisp, 'disp', 1)


call cleanup_for_mode()


end program



subroutine cleanup_for_mode()
    use params
    use allocation_module
    implicit none 

    call deallocate_if_allocated(zstore)
    call deallocate_if_allocated(ibool)
    call deallocate_if_allocated(strain1)
    call deallocate_if_allocated(rad_id)
    call deallocate_if_allocated(unique_r)
    call deallocate_if_allocated(rstore)
    call deallocate_if_allocated(thetastore)
    call deallocate_if_allocated(phistore)
    call deallocate_if_allocated(u_spl)
    call deallocate_if_allocated(udot_spl)
    call deallocate_if_allocated(v_spl)
    call deallocate_if_allocated(vdot_spl)
    call deallocate_if_allocated(interp_id_r)
    call deallocate_if_allocated(globalstrain)
    call deallocate_if_allocated(xx)
    call deallocate_if_allocated(zz)
    call deallocate_if_allocated(x_glob)
    call deallocate_if_allocated(y_glob)
    call deallocate_if_allocated(z_glob)



end subroutine cleanup_for_mode




















subroutine check_ibool_is_defined()
    ! Checks that ibool is allocated and not just zero
    ! TODO: call ibool loading if not allocated? 
    use params, only: ibool,nglob
    implicit none 

    if(.not. allocated(ibool))then 
        write(*,*)'ERROR: ibool is not allocated but is about to be used.'
        stop
    else 
        if (nglob.ne.maxval(ibool))then
            write(*,*)'ERROR: nglob is not equal to the maximum ibool value'
            write(*,*)'nglob    : ', nglob
            write(*,*)'max ibool: ', maxval(ibool)
            stop
        endif
    endif 

end subroutine check_ibool_is_defined



subroutine load_ibool(iproc, region)
    use params
    use allocation_module, only: allocate_if_unallocated
    implicit none 
    integer :: iproc, region
    character(len=250) :: varname 

    call allocate_if_unallocated(ngllx,nglly,ngllz,nspec, ibool)
    varname = 'ibool'
    call read_integer_proc_variable(iproc, region, ibool, varname)

    call check_ibool_is_defined()
end subroutine





subroutine read_proc_coordinates(iproc, region)
    use params, only: ngllx, nglly, ngllz, nspec, & 
                      xstore, ystore, zstore, nglob, &
                      datadir 
    use allocation_module, only: allocate_if_unallocated
    implicit none 
    
    ! IO variables: 
    integer            :: iproc
    integer            :: region

    ! Local variables
    integer :: IIN, ier
    character(len=250) :: binname

    IIN = 1

    ! File name prefix: 
    write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',iproc,'_'//'reg',region,'_'

    write(*,'(/,a,/)')'• Reading files from '//trim(datadir)
    write(*,'(a,i1)')'  -- region      : ', region
    write(*,'(a,i0.6,/)')'  -- processor id: ', iproc


    ! Read processor info to get ngll and nspec
    open(unit=IIN,file=trim(binname)//'info.bin', &
    status='unknown',form='unformatted',action='read',iostat=ier)
    if (ier.ne.0)then 
        write(*,'(a, i0.6)')'Couldnt read info file for proc ', iproc
        stop
    endif 
    read(IIN)nglob
    read(IIN)nspec
    read(IIN)ngllx
    read(IIN)nglly
    read(IIN)ngllz
        
    write(*,'(a)')'  -- Info: '
    write(*,*)'    --> nglob: ', nglob
    write(*,*)'    --> nspec: ', nspec
    write(*,'(a,i1)')'     --> ngllx: ', ngllx
    write(*,'(a,i1)')'     --> nglly: ', nglly
    write(*,'(a,i1)')'     --> ngllz: ', ngllz


    ! Allocate mesh arrays: 
    call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, xstore)
    call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, ystore)
    call allocate_if_unallocated(ngllx, nglly, ngllz, nspec, zstore)


    ! Open the x coordinate and load: 
    open(unit=IIN,file=trim(binname)//'xstore.bin', &
    status='unknown',form='unformatted',action='read',iostat=ier)
    if (ier.ne.0)then 
        write(*,'(a, i0.6)')'Couldnt read xstore file for proc ', iproc
        stop
    endif 
    read(IIN)xstore

    write(*,'(/,a)')'  -- X coordinates:'
    write(*,*)'     --> min. value: ', minval(xstore)
    write(*,*)'     --> max. value: ', maxval(xstore)



    ! Open the y coordinate and load: 
    open(unit=IIN,file=trim(binname)//'ystore.bin', &
    status='unknown',form='unformatted',action='read',iostat=ier)
    if (ier.ne.0)then 
        write(*,'(a, i0.6)')'Couldnt read ystore file for proc ', iproc
        stop
    endif 
    read(IIN)ystore

    write(*,'(/,a)')'  -- Y coordinates:'
    write(*,*)'     --> min. value: ', minval(ystore)
    write(*,*)'     --> max. value: ', maxval(ystore)


    ! Open the z coordinate and load: 
    open(unit=IIN,file=trim(binname)//'zstore.bin', &
    status='unknown',form='unformatted',action='read',iostat=ier)
    if (ier.ne.0)then 
        write(*,'(a, i0.6)')'Couldnt read zstore file for proc ', iproc
        stop
    endif 
    read(IIN)zstore


    write(*,'(/,a)')'  -- Z coordinates:'
    write(*,*)'     --> min. value: ', minval(zstore)
    write(*,*)'     --> max. value: ', maxval(zstore)


    return 

end subroutine read_proc_coordinates





subroutine read_proc_variable(iproc, region, variable, varname)
    use params, only: ngllx, nglly, ngllz, nspec, datadir

    implicit none 
    include "precision.h"

    ! IO variables: 
    character(len=250) :: varname 
    integer            :: iproc
    integer            :: region
    real(kind=CUSTOM_REAL)   :: variable(ngllx, nglly, ngllz, nspec)

    ! Local variables
    integer :: IIN, ier
    character(len=250) :: binname

    IIN = 1

    ! File name prefix: 
    write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',iproc,'_'//'reg',region,'_'//trim(varname)//'.bin'
    write(*,'(/,/,a)')'• Reading variable called '//trim(varname)
    write(*,'(/,a,i1)')'  -- data type : CUSTOM_REAL of length ', CUSTOM_REAL 
    write(*,'(a)')'  -- file name : '//trim(binname)

    ! Open the variable file and load: 
    open(unit=IIN,file=trim(binname), &
    status='unknown',form='unformatted',action='read',iostat=ier)
    if (ier.ne.0)then 
        write(*,'(a, i0.6)')'Couldnt read "'//trim(varname)//'" file for proc ', iproc
        stop
    endif 
    read(IIN)variable

    write(*,'(a)')'  -- '//trim(varname)//' :'
    write(*,*)'     --> min. value: ', minval(variable)
    write(*,*)'     --> max. value: ', maxval(variable)

    return 

end subroutine read_proc_variable



subroutine read_integer_proc_variable(iproc, region, variable, varname)
    use params, only: ngllx, nglly, ngllz, nspec,datadir

    implicit none 
    include "precision.h"

    ! IO variables: 
    character(len=250) :: varname 
    integer            :: iproc
    integer            :: region
    integer            :: variable(ngllx, nglly, ngllz, nspec)

    ! Local variables
    integer :: IIN, ier
    character(len=250) :: binname

    IIN = 1

    ! File name prefix: 
    write(binname,'(a,i0.6,a,i1,a)')trim(datadir)//'/proc',iproc,'_'//'reg',region,'_'//trim(varname)//'.bin'
    write(*,'(/,/,a)')'• Reading variable called '//trim(varname)
    write(*,'(a,i1)')'  -- data type : integer'
    write(*,'(a)')'  -- file name : '//trim(binname)

    ! Open the variable file and load: 
    open(unit=IIN,file=trim(binname), &
    status='unknown',form='unformatted',action='read',iostat=ier)
    if (ier.ne.0)then 
        write(*,'(a, i0.6)')'Couldnt read "'//trim(varname)//'" file for proc ', iproc
        stop
    endif 
    read(IIN)variable

    write(*,'(a)')'  -- '//trim(varname)//' :'
    write(*,*)'     --> min. value: ', minval(variable)
    write(*,*)'     --> max. value: ', maxval(variable)

    return 

end subroutine read_integer_proc_variable
