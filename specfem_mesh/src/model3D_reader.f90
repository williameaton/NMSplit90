module model3d
    use params, only: verbose, myrank
    use m_KdTree, only: KdTree, KdTreeSearch
    use dArgDynamicArray_Class, only: dArgDynamicArray
    use specfem_mesh, only: SetMesh


    implicit none
    include "constants.h"


    type :: M3D


        integer :: npts                                     ! Global points in model
        integer :: nspat                                    ! Number of spatially-varying params
        integer, allocatable :: idspats(:)                  !   --> IDs for variable type
        integer :: nconst                                   ! Number of spatially-constant params
        integer, allocatable :: idconsts(:)                 !   --> IDs for variable type
        real(kind=SPLINE_REAL), allocatable :: valconsts(:) !   --> the constant value 

        logical :: exists_spat
        logical :: exists_const


        real(kind=SPLINE_REAL), allocatable :: xcoord(:)
        real(kind=SPLINE_REAL), allocatable :: ycoord(:)
        real(kind=SPLINE_REAL), allocatable :: zcoord(:)
        real(kind=SPLINE_REAL), allocatable :: valspats(:,:)


        type(KdTree) :: kdtree

        character(len=400) :: filename

        ! Ok this is a bit insane right now but should allow for later
        ! incorporation for up to to 50 variables (inc 21 elastic parameters if
        ! fully anisotropic) required without rearranging of ID orders
        character(len=10) :: varIDcodes(50) = &
        [character(len=10) :: "A", "C", "L", "N",  "F", & ! 1  - 5
                             "κ", "μ", "ρ", "Vp", "Vs", & ! 6  - 10
                             "Vb", "",  "",  "",   "",  & ! 11 - 15
                             "",   "",  "",  "",   "",  & ! 16 - 20
                             "η1", "η2",  "",  "", "",  & ! 21 - 25
                             "",   "",  "",  "",   "",  & ! 26 - 30
                             "",   "",  "",  "",   "",  & ! 31 - 35
                             "",   "",  "",  "",   "",  & ! 36 - 40
                             "",   "",  "",  "",   "",  & ! 41 - 45
                             "",   "",  "",  "",   ""]


        contains
            procedure :: read_model_from_file
            procedure :: create_KDtree
            procedure :: project_to_gll
    end type  M3D


    contains 




    subroutine read_model_from_file(self)
        ! Reads 3D model points from text file. The model text file has 
        ! the following format: 
        ! Line  1           : Comment line denoted with a hashtag # max length 300 
        ! Line  2           : npts (integer)
        ! Line  3           : nspat
        ! Line  4           : idspat(1)
        !    ...  
        ! Line  [4 + nspat] : idspat(nspat)
        ! Line  [5 + nspat] : nconstants
        ! Line  [6 + nspat] : idcons(1)  valcons(1)
        !    ... 
        ! Line  [6 + nspat + nconstants] : idcons(nconstants)  valcons(nconstants)
        ! Following then, the lines are a list of coordinates and the values 
        ! in the order of the variables listed in the idspat's
        ! px(1) py(1) pz(1)  v1(1)  v2(1) ... vm(1)
        ! px(2) py(2) pz(2)  v1(2)  v2(2) ... vm(2)
        ! ...
        ! px(npts) py(npts) pz(npts)  v1(npts)  v2(npts) ... vm(npts) 
        implicit none
        include "constants.h"
        class(M3D) :: self

        character(len=350) :: trash
        integer :: i, ios
        character(len=20) :: fmtstr

        real(kind=SPLINE_REAL) :: t1, t2


        if (verbose.ge.2)then
            if(myrank.eq.0)write(*,*)' - Loading 3D model'
        endif 


        ! Open the model file: 
        open(1, file=trim(self%filename), form='formatted', iostat=ios)
        if (ios.ne.0) then
           write(*, *) "Error: Unable to open file ", trim(self%filename)
           stop
        endif 


        ! First line should be a comment line: 
        read(1, *)trash  

        ! Number of model grid points: 
        read(1, *)self%npts

        ! Determine number of spatially-varying model params 
        ! and then load their IDs (variable types)
        read(1, *)self%nspat
        ! Check if the number of spatially varying is >= 0
        ! if so then all fine. if less than 0 then throw error
        !
        if (self%nspat.gt.0) then 
            allocate(self%idspats(self%nspat))
            do i = 1, self%nspat
                read(1,*)self%idspats(i)
            enddo 
            self%exists_spat = .true.

            allocate(self%valspats(self%npts, self%nspat))

        elseif(self%nspat.eq.0) then 
            ! No spatially-varying variables: 
            self%exists_spat = .false.
        else 
            write(*,*)'Error reading 3D model. Nspat should be >= 0 but was "', self%nspat,' " ' 
        endif 



        ! Same for variables that are constant but also load their values: 
        read(1, *)self%nconst
        if(self%nconst.gt.0)then 
            allocate(self%idconsts(self%nconst))
            allocate(self%valconsts(self%nconst))
            do i = 1, self%nconst
                read(1,*)self%idconsts(i), self%valconsts(i)
            enddo 
            self%exists_const = .true.
        elseif(self%nconst.eq.0)then
            self%exists_const = .false.
        else 
            write(*,*)'Error reading 3D model. Nconst should be >= 0 but was "', self%nconst,' " ' 
        endif 


        if (verbose.gt.2)then
            if(myrank.eq.0)then
                write(*,*)
                write(*,*)'3D Model parameters: '
                write(*,*)' -- Number of points  ', self%npts
                write(*,'(a, i3)')' -- Number of constant variables: ', self%nconst
                do i = 1, self%nconst
                    write(*,'(a, i2, a, f12.6)')'     --> Value: '// trim(self%varIDcodes(self%idconsts(i)))//" (", self%idconsts(i) ,") --",  self%valconsts(i)
                enddo 
                write(*,'(a, i2)')' -- Number of spatially-variable variables:  ', self%nspat
                write(*,*)'     --> IDs:  '
                do i = 1, self%nspat
                    write(*,'(a)') '              '//trim(self%varIDcodes(self%idspats(i)))
                enddo

                write(*,*)' -- Reading spatially variable content ... '

            endif 
        endif 




        ! Allocate memory for the x, y, z and spatially varying materials: 
        allocate(self%xcoord(self%npts))
        allocate(self%ycoord(self%npts))
        allocate(self%zcoord(self%npts))

        ! Read in each of the coordinates and it associated values:

        ! Check at least one variable is meant to be read -- could make this a warning: 
        ! Use this loop to distinguish how the variable/coord part is read
        if(.not.self%exists_spat)then
            ! ERROR check here: 
            if(.not.self%exists_const)then
                write(*,*)'Error in reading 3D model: there are no variables for either spatially-varying or constant. Stop.'
                stop 
            endif 
            ! Read with no spatially varying variables: 
            do i = 1, self%npts
                read(1, *)self%xcoord(i), self%ycoord(i), self%zcoord(i)
            enddo 
        else 
            ! Read with variables! 
            do i = 1, self%npts
                read(1, *)self%xcoord(i), self%ycoord(i), self%zcoord(i), self%valspats(i, :)
            enddo 
        endif


    end subroutine read_model_from_file



    subroutine create_KDtree(self)
        implicit none
        include "constants.h"
        class(M3D) :: self

        self%kdtree = KdTree(self%xcoord, self%ycoord, self%zcoord) 

        if(verbose.ge.2)write(*,*)'Created KD tree for 3D Model...'
    end subroutine create_KDtree



    subroutine project_to_gll(self, sm, globvar, id)
        implicit none
        include "constants.h"
        class(M3D) :: self
        type(KdTreeSearch)     :: search
        type(dArgDynamicArray) :: da
        type(SetMesh)          :: sm 
        real(kind=CUSTOM_REAL) :: globvar(sm%nglob)
        integer :: id , point, i, cluster_size, ierr
        
        globvar = zero
        
        ! For each global coordinate we want to find its nearest neighbour 
        write(*,*)'Projecting 3D model to GLL...'

        do i = 1, sm%nglob
            da = search%kNearest(self%kdtree, self%xcoord, self%ycoord, self%zcoord, & 
                                 xQuery = sm%x_glob(i), &
                                 yQuery = sm%y_glob(i), &
                                 zQuery = sm%z_glob(i), & 
                                 k = 1)
                                 
            ! Index of the closest voronoi point 
            point =  da%i%values(1)
            globvar(i) = self%valspats(point,id)

            if(verbose.ge.3)then 
                if (MOD(i,(sm%nglob/20)).eq.0)then 
                    write(*,*)' Progress: ', (1 +int(100.0d0*real(i)/real(sm%nglob))), '%'
                endif
            endif 
        enddo 

        if(verbose.ge.3)then 
            write(*,*)'Projected to GLL: '
            write(*,*)'     3D Model variable ID    -   ', id
            write(*,*)'     min. value:    -   ', minval(globvar)
            write(*,*)'     max. value:    -   ', maxval(globvar)
            write(*,*)
        endif


    end subroutine project_to_gll





end module 
