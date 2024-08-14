program read_specfem_mesh
use params 

implicit none 
include "precision.h"


integer :: iproc                  
integer :: nprocs                ! Number of processors used by mesher 
integer :: region                ! Region code
character(len=250) :: varname 

! Setup parameters: 
region = 3      ! Inner core
iproc  = 3
nprocs = 54


! Read the mesh info and coordinates
call read_proc_coordinates(iproc, region)

! Load density variable: 
allocate(rho(ngllx,nglly,ngllz,nspec))
varname = 'rho'
call read_proc_variable(iproc, region, rho, varname)


! Load ibool variable: 
call load_ibool(iproc, region)



! Example mapping of local to global
allocate(globalrho(nglob))
call map_local_global_custom_real(rho, globalrho, 0)



end program





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
    implicit none 
    integer :: iproc, region
    character(len=250) :: varname 

    allocate(ibool(ngllx,nglly,ngllz,nspec))
    varname = 'ibool'
    call read_integer_proc_variable(iproc, region, ibool, varname)

    call check_ibool_is_defined()
end subroutine





subroutine read_proc_coordinates(iproc, region)
    use params, only: ngllx, nglly, ngllz, nspec, & 
                      xstore, ystore, zstore, nglob, &
                      datadir 

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
    allocate(xstore(ngllx,nglly,ngllz,nspec))
    allocate(ystore(ngllx,nglly,ngllz,nspec))
    allocate(zstore(ngllx,nglly,ngllz,nspec))


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
