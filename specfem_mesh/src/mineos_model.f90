module mineos_model

    implicit none 
    include "constants.h"

    type ::  MineosModel
        integer :: NR
        integer :: IC_ID
        integer :: CMB_ID
        integer :: ndisc

        real(kind=CUSTOM_REAL) :: RA, min_r, max_r
        character(len=250)     :: model_file

        ! Radial discontinuities
        !   ndisc: number of discontinuities (plus 0 and surface)
        !   rdisc: radius of discontinuities (plus 0 and surface)
        !   disc: knot ID of discontinuity (-ve side)
        integer, allocatable                :: disc(:)  
        real(kind=CUSTOM_REAL), allocatable :: rdisc_norm(:)           
        real(kind=CUSTOM_REAL), allocatable :: rdisc(:)  

        real(kind=CUSTOM_REAL), allocatable :: radius(:)
        real(kind=CUSTOM_REAL), allocatable :: rad_mineos(:)
        real(kind=CUSTOM_REAL), allocatable :: rho_mineos(:)
        real(kind=CUSTOM_REAL), allocatable :: vp_mineos(:)

        contains 
            procedure :: process_mineos_model
            procedure :: load_mineos_radial_info_MPI
            procedure :: load_mineos_radial_info
            procedure :: find_disc
            procedure :: save_mineos_model
            procedure :: allocate_radii
    end type MineosModel

    ! Globally available
    type(MineosModel), target  :: mineos
    type(MineosModel), pointer :: mineos_ptr 

    contains 

    subroutine process_mineos_model(self, save_model)
        use params, only: ddir, model_fname, verbose, IIN
        implicit none

        ! IO :
        class(MineosModel) :: self
        logical :: save_model

        ! Local variables: 
        integer                :: iomod, ios, intjunk, i 
        real(kind=CUSTOM_REAL) :: realjunk, rho,  vph, vsv, vsh, vsv_prev

        iomod = IIN

        ! read model file to get radius array (normalised 0 to 1)
        self%model_file = trim(ddir)//trim(model_fname)
        
        if(verbose.ge.2) write(*,'(/,a)')'- Reading Model file ....', self%model_file

        open(iomod, file = trim(self%model_file), status = 'old', iostat = ios)
        ! Error check
        if (ios .ne. 0) stop 'Error reading model file'

        ! Read once to get basic parameters
        self%NR = -1
        do while (ios.eq.0)
            read(iomod, *, iostat=ios) intjunk, realjunk
            self%NR = self%NR + 1
        enddo 
        close(iomod)
        
        ! This is the final radius:
        self%RA = realjunk

        allocate(self%radius(self%NR))
        allocate(self%rad_mineos(self%NR))
        allocate(self%rho_mineos(self%NR))
        allocate(self%vp_mineos(self%NR))

        ! Warning if not radius 6371 km 
        if (self%RA.ne.SCALE_R .and. verbose.ge.1)then 
            write(*,*)'Warning: Radius is not 6371 km but that is being used for general non-dimensionalisation scaling'
            write(*,*)'         Maximum value may of non-dimensional rad_mineos may not be 1'
        endif 

        ! Second read of the file: 
        open(iomod, file = self%model_file, status = 'old', iostat = ios)

        vsv_prev = 99.99_CUSTOM_REAL
        do i = 1 , self%NR
            read(iomod,*,iostat=ios)  intjunk, self%radius(i), rho, self%vp_mineos(i), vph, vsv, vsh

            self%rho_mineos(i) = rho / RHOAV

            ! Detect ICB 
            if(vsv_prev.gt.ZERO .and. vsv.eq.ZERO)then
                self%IC_ID = i - 1 
            endif
            ! Detect ICB 
            if(vsv_prev.eq.ZERO .and. vsv.gt.ZERO)then
                self%CMB_ID = i - 1 
            endif

            vsv_prev = vsv
        enddo
        close(iomod)

        ! Non dimensionalised mineos radius
        self%rad_mineos = self%radius / SCALE_R

        ! Non dimensionalise the Vp: 
        self%vp_mineos = self%vp_mineos * (SCALE_T/SCALE_R)

        ! Find the discontinuities
        call self%find_disc()

        if(verbose.ge.2) write(*,'(a,/)')'--> finished reading model file.'

        ! Saving model if needed
        if(save_model)then 
            call self%save_mineos_model()
        endif 


        self%min_r = minval(self%rad_mineos)
        self%max_r = maxval(self%rad_mineos)
        if(verbose.ge.4)then 
            write(*,*)'-- Model details: '
            write(*,*)'    - NR                : ', self%NR
            write(*,*)'    - Minimum rad_mineos: ', self%min_r
            write(*,*)'    - Maximum rad_mineos: ', self%max_r
        endif

    end subroutine process_mineos_model



    subroutine find_disc(self)
        ! Modified version of find_disc in get_mode
        ! Note that this is updated so that the first 'discontinuity' 
        ! is now the centre (r=0) and the last is the surface (r=R)
        implicit none
        include "precision.h"
        
        class(MineosModel) :: self

        real(8) :: radius_last
        integer :: jdisc, iomod, i
        
        ! Temporary discontinuity arrays
        integer :: disctmp(self%NR)
        real(kind=CUSTOM_REAL) :: rdisctmp(self%NR)

        radius_last = -1.0
        jdisc = 0
        iomod = 11
    
        do i = 1, self%NR
        if (abs(self%radius(i)-radius_last) .lt. 1.0d-6) then
            jdisc = jdisc + 1
            rdisctmp(jdisc) = self%radius(i)
            disctmp(jdisc) = i - 1
        endif
        radius_last = self%radius(i)
        enddo
        close(iomod)
        
        ! Allocate arrays of the length ndisc
        ! + 2 for the centre and the surface
        self%ndisc = jdisc + 2
        allocate(self%rdisc(self%ndisc), self%disc(self%ndisc))
        self%rdisc(2:self%ndisc-1) = rdisctmp(1:self%ndisc-2)
        self%disc(2:self%ndisc-1)  = disctmp(1:self%ndisc-2)

        ! centre
        self%disc(1)      = 0
        self%rdisc(1)     = self%radius(1) 
        
        ! surface
        self%disc(self%ndisc)  = self%NR
        self%rdisc(self%ndisc) = self%radius(self%NR)   
        self%rdisc_norm = self%rdisc / self%radius(self%NR)   

    end subroutine find_disc


    subroutine allocate_radii(self)
        implicit none 
        class(MineosModel) :: self 
        allocate(self%radius(self%NR))
        allocate(self%rad_mineos(self%NR))
    end subroutine allocate_radii


    subroutine load_mineos_radial_info(self)
        use params, only:  datadir, IIN
        implicit none 
        class(MineosModel) :: self 

        open(IIN, file=trim(datadir)//'/store/mineos_model/radial_data', form='unformatted')
        read(IIN)self%NR
        read(IIN)self%IC_ID
        read(IIN)self%CMB_ID

        call self%allocate_radii()

        read(IIN)self%rad_mineos
        read(IIN)self%radius
        close(IIN)
    end subroutine load_mineos_radial_info



    subroutine load_mineos_radial_info_MPI(self)
        use params, only: myrank, datadir, MPI_CUSTOM_REAL
        implicit none 
        class(MineosModel) :: self 

    #ifdef WITH_MPI
        include 'mpif.h'
        
        integer :: ierr
    
        if(myrank.eq.0)then 
            call self%load_mineos_radial_info()
        endif 
    
        call MPI_Bcast(self%NR,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(self%IC_ID,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(self%CMB_ID, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
        if(myrank.ne.0) allocate(self%radius(self%NR), self%rad_mineos(self%NR))
    
        call MPI_Bcast(self%radius, self%NR, MPI_CUSTOM_REAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(self%rad_mineos, self%NR, MPI_CUSTOM_REAL, 0, MPI_COMM_WORLD, ierr)
    #endif 
    end subroutine load_mineos_radial_info_MPI
    
    
    
    
    
    subroutine load_mineos_density(self)
        use params, only: datadir
        implicit none 
        class(MineosModel) :: self
    
        allocate(self%vp_mineos(self%NR))
    
        open(1, file=trim(datadir)//'/store/mineos_model/rho_mineos', form='unformatted')
        read(1)self%rho_mineos
        close(1)
    end subroutine load_mineos_density
    
    
    subroutine load_mineos_vp(self)
        use params, only: datadir
        implicit none 
        class(MineosModel) :: self
    
        allocate(self%rho_mineos(self%NR))
    
        open(1, file=trim(datadir)//'/store/mineos_model/vp_mineos', form='unformatted')
        read(1)self%vp_mineos
        close(1)
    end subroutine load_mineos_vp
    
    subroutine load_mineos_discontinuities(self)
        use params, only: datadir
        implicit none 
        class(MineosModel) :: self

        open(1, file=trim(datadir)//'/store/mineos_model/disc', form='unformatted')
        read(1)self%disc
        read(1)self%rdisc
        read(1)self%rdisc_norm
        read(1)self%disc
        close(1)
    end subroutine load_mineos_discontinuities
    
    subroutine save_mineos_model(self)
        use params, only: datadir
        implicit none 
        class(MineosModel) :: self

        ! Radius and knot info
        open(1, file=trim(datadir)//'/store/mineos_model/radial_data', form='unformatted')
        write(1)self%NR
        write(1)self%IC_ID
        write(1)self%CMB_ID
        write(1)self%rad_mineos
        write(1)self%radius
        close(1)

        ! Density
        open(1, file=trim(datadir)//'/store/mineos_model/rho_mineos', form='unformatted')
        write(1)self%rho_mineos
        close(1)

        ! Vp 
        open(1, file=trim(datadir)//'/store/mineos_model/vp_mineos', form='unformatted')
        write(1)self%vp_mineos
        close(1)

        ! Discontinuity info
        open(1, file=trim(datadir)//'/store/mineos_model/disc', form='unformatted')
        write(1)self%ndisc
        write(1)self%rdisc
        write(1)self%rdisc_norm
        write(1)self%disc
        close(1)

    end subroutine save_mineos_model

    





end module mineos_model