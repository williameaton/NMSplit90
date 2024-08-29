
subroutine setup_ensight_for_proc(iproc, region, part)
    
    implicit none 
    integer :: iproc, region, part

    call create_ensight_file_prefix(iproc, region)
    call create_proc_case_file()
    call create_proc_geo_file(part)
    return 
end subroutine setup_ensight_for_proc




subroutine create_proc_geo_file(part)
    use params, only: GEOUNIT, nglob, x_glob, y_glob, z_glob, nspec, & 
                      ngllx, nglly, ngllz, en_fname, intfmt, realfmt, en_dir, ibool
    implicit none 
    
    ! IO variables
    integer :: part

    ! Local variables
    character(len=78) :: buffer
    integer :: i, nhex, ispec, j, k


    open(unit=GEOUNIT,file=trim(en_dir)//trim(en_fname)//'.geo', status='unknown',form='formatted',action='write')

    buffer = 'WE description 1'
    write(GEOUNIT, '(a)')buffer
    buffer = 'WE description 2'
    write(GEOUNIT, '(a)')buffer
    buffer = 'node id off'
    write(GEOUNIT, '(a)')buffer
    buffer = 'element id off'
    write(GEOUNIT, '(a)')buffer
    buffer = 'part'
    write(GEOUNIT, '(a)')buffer
    write(GEOUNIT, intfmt) part 

    ! COORDINATES
    buffer = 'coord_desc_line_we'   ! coordinate description line
    write(GEOUNIT, '(a)')buffer
    buffer = 'coordinates'          ! line with 'coordinates'
    write(GEOUNIT, '(a)')buffer
    write(GEOUNIT, intfmt) nglob    ! number of coordinates
   
    do i = 1, nglob
        write(GEOUNIT, realfmt)x_glob(i)
    enddo 
    do i = 1, nglob
        write(GEOUNIT, realfmt)y_glob(i)
    enddo 
    do i = 1, nglob
        write(GEOUNIT, realfmt)z_glob(i)
    enddo 


    ! ELEMENTS 
    ! Using Hexa8 where each element is subdivided into 8-node
    ! hexes. There will therefore be (ngllx-1)(nglly-1)(ngllz-1) hexes
    ! in each spectral element 
    nhex = nspec*(ngllx-1)*(nglly-1)*(ngllz-1)

    buffer = 'hexa8'   
    write(GEOUNIT, '(a)')buffer ! Element type
    write(GEOUNIT, intfmt) nhex ! Number of elements
    
    do ispec = 1, nspec
        do i = 1, ngllx-1
            do j = 1, nglly-1
                do k = 1, ngllz-1
                    write(GEOUNIT,'(8i10)') ibool(i,j,k,ispec),       & 
                                            ibool(i+1,j,k,ispec),     &
                                            ibool(i+1,j+1,k,ispec),   &
                                            ibool(i,j+1,k,ispec),     &
                                            ibool(i,j,k+1,ispec),     & 
                                            ibool(i+1,j,k+1,ispec),   &
                                            ibool(i+1,j+1,k+1,ispec), &
                                            ibool(i,j+1,k+1,ispec)
                enddo 
            enddo 
        enddo 
    enddo

    close(GEOUNIT)


end subroutine create_proc_geo_file



subroutine create_proc_case_file()
    use params, only: CASEUNIT, en_fname, intfmt, realfmt, en_dir
    implicit none 
    
    ! IO variables

    open(unit=CASEUNIT,file=trim(en_dir)//trim(en_fname)//'.case', & 
          status='unknown',form='formatted',action='write')


    write(CASEUNIT,'(a)')'FORMAT'
    write(CASEUNIT,'(a,/)')'type:  ensight gold'

    write(CASEUNIT,'(a)')'GEOMETRY'
    write(CASEUNIT,'(a,a/)')'model:    ',trim(en_fname)//'.geo'

    write(CASEUNIT,'(a)')'VARIABLE'

    close(CASEUNIT)

end subroutine create_proc_case_file



subroutine create_ensight_file_prefix(iproc, region)
    
    use params 
    implicit none 

    integer :: iproc, region
    character(len=250) :: format_string

    ! Create geometry file 
    if (iproc .lt. 10) then
        format_string = "(a,I1,a,I1)"
    elseif (iproc.ge.10 .and. iproc.lt.100)then
        format_string = "(a,I1,a,I2)"
    elseif (iproc.ge.100 .and. iproc.lt.1000)then
        format_string = "(a,I1,a,I3)"
    endif
    write(en_fname,format_string)'reg',region,'proc', iproc

    if(verbose.ge.2)write(*,'(/,a)')'Ensight prefix: '//trim(en_fname)

end subroutine




subroutine write_complex_symtensor_to_ensight(symten, suffix, part)
    ! Write a complex symmetric tensor to ensight
    ! Note that ensight gold supports complex scalars and vectors
    ! as well as symmetric tensors, but NOT complex tensors! 
    ! Here the choice is to output the real and the imaginary 
    ! parts of the tensor separately
    use params, only: nglob, en_dir, en_fname, TENSORSYMOUT_R, & 
                      TENSORSYMOUT_I, intfmt, realfmt, CASEUNIT

    implicit none 
    include "constants.h"

    ! IO variables: 
    integer :: part
    complex(kind=SPLINE_REAL) :: symten(6, nglob)
    character(len=*) :: suffix

    ! Local 
    integer :: i, component
    character(len=79)   :: buffer
    character(len=250)  :: fname 


    ! Case file: 
    open(unit=CASEUNIT,file=trim(en_dir)//trim(en_fname)//'.case', & 
          status='old', form='formatted', action='write', position='append')

    ! Real variable file 
    fname = trim(en_dir)//trim(en_fname)//'.'//trim(suffix)//'_real'
    open(unit=TENSORSYMOUT_R, file=trim(fname), & 
         status='unknown',form='formatted',action='write')

    write(CASEUNIT,'(a/)')'tensor symm per node: '//suffix//'_real '//trim(en_fname)//'.'//trim(suffix)//'_real'

    ! Imaginary variable file 
    fname = trim(en_dir)//trim(en_fname)//'.'//trim(suffix)//'_imag'
    open(unit=TENSORSYMOUT_I,file=trim(fname), & 
         status='unknown',form='formatted',action='write')
    
    write(CASEUNIT,'(a/)')'tensor symm per node: '//suffix//'_imag '//trim(en_fname)//'.'//trim(suffix)//'_imag'

    ! We no longer need the case file 
    close(CASEUNIT)


    ! Write the bits we need 
    buffer = 'SymTensor Real for '//trim(suffix)
    write(TENSORSYMOUT_R, '(a)')buffer
    buffer = 'SymTensor Imag for '//trim(suffix)
    write(TENSORSYMOUT_I, '(a)')buffer

    ! Add the part lines
    buffer = 'part'
    write(TENSORSYMOUT_R, '(a)')buffer
    write(TENSORSYMOUT_I, '(a)')buffer

    write(TENSORSYMOUT_R, intfmt)part
    write(TENSORSYMOUT_I, intfmt)part

    ! Add coordinate tag line
    buffer = 'coordinates'
    write(TENSORSYMOUT_R, '(a)')buffer
    write(TENSORSYMOUT_I, '(a)')buffer


    do component = 1, 6
        do i = 1, nglob
            write(TENSORSYMOUT_R, realfmt)  real(symten(component, i))
            write(TENSORSYMOUT_I, realfmt) aimag(symten(component, i))
        enddo 
    enddo 

    close(TENSORSYMOUT_R)
    close(TENSORSYMOUT_I)

end subroutine




subroutine write_complex_vector_to_ensight(comvec, suffix, part)
    ! Write a complex vector to ensight
    use params, only: nglob, en_dir, en_fname, VECOUT_R, VECOUT_I, & 
                      intfmt, realfmt, CASEUNIT

    implicit none 
    include "constants.h"

    ! IO variables: 
    integer :: part
    complex(kind=SPLINE_REAL) :: comvec(3, nglob)
    character(len=*) :: suffix

    ! Local 
    integer :: i, component
    character(len=79)   :: buffer
    character(len=250)  :: r_fname, i_fname 

    ! Case file: 
    open(unit=CASEUNIT,file=trim(en_dir)//trim(en_fname)//'.case', & 
          status='old', form='formatted', action='write', position='append')

    ! Real variable file 
    r_fname = trim(en_fname)//'.'//trim(suffix)//'_real'
    open(unit=VECOUT_R, file=trim(en_dir)//trim(r_fname), & 
         status='unknown',form='formatted', action='write')
    
    ! Imaginary variable file 
    i_fname = trim(en_fname)//'.'//trim(suffix)//'_imag'
    open(unit=VECOUT_I,file=trim(en_dir)//trim(i_fname), & 
         status='unknown',form='formatted', action='write')
    
    write(CASEUNIT,'(a/)')'complex vector per node: '//suffix//'  '//trim(r_fname)//'  '//trim(i_fname)//'  UNDEFINED'

    ! We no longer need the case file 
    close(CASEUNIT)


    ! Write the bits we need 
    buffer = 'Complex vec real for '//trim(suffix)
    write(VECOUT_R, '(a)')buffer
    buffer = 'Complex vec imag for '//trim(suffix)
    write(VECOUT_I, '(a)')buffer

    ! Add the part lines
    buffer = 'part'
    write(VECOUT_R, '(a)')buffer
    write(VECOUT_I, '(a)')buffer

    write(VECOUT_R, intfmt)part
    write(VECOUT_I, intfmt)part

    ! Add coordinate tag line
    buffer = 'coordinates'
    write(VECOUT_R, '(a)')buffer
    write(VECOUT_I, '(a)')buffer

    do component = 1, 3
        do i = 1, nglob
            write(VECOUT_R, realfmt)  real(comvec(component, i))
            write(VECOUT_I, realfmt) aimag(comvec(component, i))
        enddo 
    enddo 

    close(VECOUT_R)
    close(VECOUT_I)

end subroutine write_complex_vector_to_ensight





subroutine write_real_scalar_to_ensight(realscal, suffix, part)
    ! Write a real scalar array to ensight
    use params, only: nglob, en_dir, en_fname, REALSCALOUT, & 
                      intfmt, realfmt, CASEUNIT

    implicit none 
    include "constants.h"

    ! IO variables: 
    integer :: part
    real(kind=CUSTOM_REAL) :: realscal(nglob)
    character(len=*) :: suffix

    ! Local 
    integer :: i
    character(len=79)   :: buffer
    character(len=250)  :: fname


    ! Variable file 
    fname = trim(en_fname)//'.'//trim(suffix)
    open(unit=REALSCALOUT, file=trim(en_dir)//trim(fname), & 
         status='unknown',form='formatted', action='write')

    ! Case file: 
    open(unit=CASEUNIT,file=trim(en_dir)//trim(en_fname)//'.case', & 
          status='old', form='formatted', action='write', position='append')
    write(CASEUNIT,'(a/)')'scalar per node: '//suffix//'  '//trim(fname)
    close(CASEUNIT)

    ! Write the bits we need 
    buffer = 'Scalar for '//trim(suffix)
    write(REALSCALOUT, '(a)')buffer

    ! Add the part lines
    buffer = 'part'
    write(REALSCALOUT, '(a)')buffer
    write(REALSCALOUT, intfmt)part

    ! Add coordinate tag line
    buffer = 'coordinates'
    write(REALSCALOUT, '(a)')buffer

    do i = 1, nglob
        write(REALSCALOUT, realfmt)realscal(i)
    enddo 

    close(REALSCALOUT)

end subroutine write_real_scalar_to_ensight