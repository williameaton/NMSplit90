module modes
    use mineos_model, only: MineosModel
    use params, only: verbose
    implicit none
    include "constants.h"


    type :: Mode

        integer                    :: n
        integer                    :: l
        integer                    :: tl1
        integer                    :: len
        integer                    :: spl_len
        character(len=1)           :: t 
        real(kind=4)               :: wcom
        real(kind=4)               :: qmod
        real(kind=SPLINE_REAL)     :: lf         ! float version of l
        real(kind=SPLINE_REAL)     :: kf         ! float of k = √(l(l+1))
        type(MineosModel), pointer :: mineos

        real(kind=SPLINE_REAL), allocatable :: u(:)
        real(kind=SPLINE_REAL), allocatable :: v(:)
        real(kind=SPLINE_REAL), allocatable :: w(:)
        real(kind=SPLINE_REAL), allocatable :: du(:)
        real(kind=SPLINE_REAL), allocatable :: dv(:)
        real(kind=SPLINE_REAL), allocatable :: dw(:)

        real(kind=SPLINE_REAL), allocatable :: u_spl(:)
        real(kind=SPLINE_REAL), allocatable :: v_spl(:)
        real(kind=SPLINE_REAL), allocatable :: w_spl(:)
        real(kind=SPLINE_REAL), allocatable :: du_spl(:)
        real(kind=SPLINE_REAL), allocatable :: dv_spl(:)
        real(kind=SPLINE_REAL), allocatable :: dw_spl(:)
    
        ! Auxillary spline variables for strain 
        real(kind=SPLINE_REAL), allocatable :: aux_x(:), aux_z(:)
   
        contains
            procedure :: get_mineos_mode
            procedure :: write_spline
    end type Mode

    contains 


    subroutine get_mineos_mode(self, save_mode, out_dir)
        ! Loads a mode from MINEOS
        ! Usage:
        !    *  input  * 
        !  type --- can be 't','s',  upper case is o.k. (char*1)
        !  nord --- order number                        (integer)
        !  l --- degree number                          (integer)
        ! 
        !    *  output *
        !  wcom --- normalized anglar frequency         (real*8)
        !  qmod --- normalized Q factor                 (real*8)
        !  NR --- number of layers of output array     (integer)
        !  u(NR),du(NR),v(NR),dv(NR) ---
        !     in case of toroidal and radial modes, only u,du
        !     is valid, for spheroidal modes, all are valid.
        use params, only: ddir, model_fname, verbose
        implicit none
        include "constants.h"
        class(Mode) :: self


        logical                    :: save_mode
        character(len=*), optional :: out_dir
    
        character(len=200) :: catalogue, bin_file, eigstring
        integer            :: ntype,nvec,i,j
        character(len=1)   :: type1,type2,char
        integer            :: ieigtxt,iomod,iocat,iobin,nrec,ios,nn,ll
        real(8)            :: wwmhz,ttcom,ggcom,qqmod,eps, raylquo
        real(4)            :: buf(6*self%len),  cg4, phsv
        integer(4)         :: n4, l4, n4old, l4old
    
        eps = 1.0d-4
    
        ! File IDs 
        iomod   = 9
        iocat   = 10
        iobin   = 11
        ieigtxt = 12
    
        ! set up correct catalogue and bin file for mode reading
        if (self%t == 'T' .or. self%t == 't') then
            ntype = 1
            catalogue = trim(ddir)//'prem_ani_att_T'
            nvec = 2 * self%len
            type1 = 'T'; type2 = 't';
        
            bin_file = trim(ddir)//'prem_ani_att_T.bin'
        
        else if (self%t == 'S' .or. self%t == 's') then
            type1 = 'S'; type2 = 's';
    
            if (self%l == 0) then
                ! Radial modes
                ntype = 2
                nvec = 2 * self%len
                catalogue = trim(ddir)//'prem_ani_att_R'
                bin_file  = trim(ddir)//'prem_ani_att_R.bin'
            else
                ! Spheroidal not radial
                ntype = 3
                catalogue = trim(ddir)//'prem_ani_att_S'
                    bin_file = trim(ddir)//'prem_ani_att_S.bin'
                    nvec = 6 * self%len 
            endif ! l=0

        elseif (self%t == 'C' .or. self%t == 'c') then
            type1 = 'C'; type2 = 'c';
            ! Inner core toroidal modes
            ntype = 4
            catalogue = trim(ddir)//'prem_ani_att_C'
            nvec = 2 * self%len        
            bin_file = trim(ddir)//'prem_ani_att_C.bin'
        else
            write(*,*)'Error in get_mineos_mode: mode type must be S, T or C but was '//self%t
        endif ! T or C or S 
    
        ! get the record number of the desired mode from catalogue file
        if(verbose.ge.3) write(*,'(a,/,2x,a)')'- Reading Catalogue file: ', catalogue
        nrec = 0
        open(iocat,file=catalogue,status='old',iostat=ios)
    
        do while (ios .eq. 0) 
        nrec = nrec + 1
        read(iocat,'(i5,1x,a1,i5,6g16.7)',iostat=ios) nn,char,ll,phsv,wwmhz,ttcom,ggcom,qqmod,raylquo
    
            ! Checks that the type of mode in the 'catalogue file' is correctly S or T 
            if (nrec == 1 .and. (char /= type1 .and. char /= type2)) then
                write(*,*)'Incorrect mode catalogue: ', nn, ll, ';', char, ';', type1, ';', type2
                stop 
            endif
    
            ! Exits when nn is self%n and l is ll 
            if (nn == self%n .and. ll == self%l) exit
        enddo
        close(iocat)
    
        if (ios .gt. 0) stop 'Error reading 1'
        if (ios .lt. 0) stop 'Mode not found in the catalogue' !when would this ever trigger? 
        
        if(verbose.ge.3) write(*,'(a,1x,i6)')'Found mode at ID', nrec
    
    
        ! reading binary file for mode eigenfunction
        !print *,' Reading Binary file ....'
        open(iobin,file=bin_file,status='old',iostat=ios,form='unformatted')
        if (ios .ne. 0) stop 'Error opening catalogue file'
        !read(iobin,rec = nrec,iostat = ios)junk1,nnn,lll,wcom,qmod,av,ah,(buf(i),i=1,nvec),junk2
    
        do j = 1, nrec
        read(iobin,iostat = ios)n4,l4, self%wcom, self%qmod, cg4,(buf(i),i=1,nvec)
    
        ! Check if we have gone past the data in our reads
        if (n4old.eq.n4.and.l4old.eq.l4)then
            write(*,*)'Mode not in catalogue. Looped past end of file binary read.'
            stop
        endif 
        n4old = n4
        l4old = l4
        enddo 
        close(iobin)
    
        
        if(verbose.ge.3)then
            write(*,*)
            write(*,*)'Angular Freq in rad :', self%wcom
            write(*,*)'Frequency in mHz    :', wwmhz
            write(*,*)'Period in seconds   :', 2*PI/self%wcom
            write(*,*)'Q value             :', self%qmod
            write(*,*)'Group velocity      :', cg4
            write(*,*)
        endif
    
        
        if(save_mode)then
            ! Save eigenfunctions to text file in column format 
            write(eigstring,'(i0,a,i0,a)')n4,type1,l4,'.txt'
            open(ieigtxt,file=trim(trim(out_dir)//eigstring), iostat=ios)
            write(ieigtxt,*)n4,  type1, l4
            do i =1, self%len
               if(type1=='S')then 
                ! Radius, U, U', V, V', P, P'
                  write(ieigtxt,'(7e18.8)')self%mineos%rad_mineos(i), buf(i), buf(i + self%len), buf(i + 2*self%len), & 
                                           buf(i + 3*self%len), buf(i + 4*self%len), buf(i + 5*self%len)
               elseif (type1=='T' .or. type1=='C')then 
                ! Radius, W, W' ?
                  write(ieigtxt,'(7e15.8)' )self%mineos%rad_mineos(i), buf(i), buf(i + self%len)
               endif 
            enddo 
            close(ieigtxt)
        endif
    
        ! Compare the freq in mHz (wwmhz) to ang freq read from binaries
        ! as well as Q values as a sanity check
        if (abs(wwmhz*2*PI/1000 - self%wcom) > eps  .or. abs(self%qmod - qqmod) > eps) then      
        stop 'Error: frequencies or Q values do not match'
        endif 
        
        ! Format buffer into arrays 
        if(ntype.eq.1 .or. ntype.eq.4)then
            ! Toroidal only
            self%w(1  : self%len) = real(buf(1 : self%len), kind=SPLINE_REAL)
            self%dw(1 : self%len) = real(buf(self%len + 1 : 2 * self%len), kind=SPLINE_REAL)
        else
            ! Spheroidal 
            self%u(1  : self%len) = real(buf(1 : self%len), kind=SPLINE_REAL)
            self%du(1 : self%len) = real(buf(self%len + 1 : 2 * self%len), kind=SPLINE_REAL)
            if (ntype == 2) then
                ! Radial
                self%v(1 : self%len)  = SPLINE_ZERO
                self%dv(1 : self%len) = SPLINE_ZERO
            else if (ntype == 3) then
                ! Spheroidal
                self%v(1 : self%len)  = real(buf(2 * self%len + 1 : 3 * self%len), kind=SPLINE_REAL)
                self%dv(1 : self%len) = real(buf(3 * self%len + 1 : 4 * self%len), kind=SPLINE_REAL)
            endif
        endif

    
    end subroutine get_mineos_mode



    subroutine write_spline(self, radius, name)
        ! Output the spline values: 
        ! Save eigenfunctions to text file in column format 
        use params, only: verbose
        implicit none 

        class(Mode) :: self 
        real(kind=CUSTOM_REAL) :: radius(self%spl_len) 
        character(len=*), optional :: name

        ! Local :
        integer           :: i
        character(len=30) :: eigstring
        character(len=2)  :: nfmt, lfmt
        character(len=13) :: fmtstring


        if (len_trim(name) .eq. 0)then
            if(self%n.ge.0 .and. self%n.lt.10)then
                nfmt = 'i1'
            elseif(self%n.ge.10 .and. self%n.lt.100)then
                nfmt = 'i2'
            else
                write(*,*)'Format not set for n = ', self%n
            endif
            if(self%l.ge.0 .and. self%l.lt.10)then
                lfmt = 'i1'
            elseif(self%l.ge.10 .and. self%l.lt.100)then
                lfmt = 'i2'
            else
                write(*,*)'Format not set for l = ', self%l
            endif
            write(fmtstring, '(a)') '(a,' // nfmt // ',a,' // lfmt // ',a)'
            write(eigstring,fmtstring) 'spline_', self%n, self%t, self%l, '.txt'
        else 
            write(eigstring,'(a)')name
        endif 

        if (verbose.ge.3)then
        write(*,'(/,a)')'Saving spline to '//trim(eigstring)
        endif 
        

        open(1, file=trim('./spline/'//eigstring))
        do i =1, self%spl_len
            if(self%t=='S')then 
                write(1,'(f12.6,f12.6,f12.6,f12.6,f12.6)')radius(i), self%u_spl(i), self%du_spl(i), self%v_spl(i), self%dv_spl(i)
            elseif (self%t=='T')then 
                write(1,'(f12.6,f12.6,f12.6)')radius(i), self%w_spl(i), self%dw_spl(i)
            endif 
        enddo 
        close(1)
    end subroutine write_spline


    function get_mode(n, t, l, mineos, out_dir) result(m)

        implicit none 
        ! in 
        integer                    :: n, l 
        character                  :: t
        type(MineosModel), pointer :: mineos
        character(len=*), optional :: out_dir
        
        ! Out
        type(Mode) :: m 


        if(verbose.ge.2) write(*,*)'Getting mode '

        m%n      = n
        m%t      = t
        m%l      = l
        m%tl1    = 2*l + 1
        m%mineos => mineos
        m%len    = m%mineos%NR 

        ! Float versions 
        ! DT98 below D.1: k = sqrt(l(l+1))
        m%lf  = real(m%l, kind=SPLINE_REAL)
        m%kf  = sqrt(m%lf*(m%lf+SPLINE_ONE))
        

        ! Toroidal mantle or core mode:
        if(t.eq.'T' .or. t.eq.'t' .or. t.eq.'C' .or. t.eq.'c')then 
            allocate(m%w(m%len))
            allocate(m%dw(m%len))
        else 
            allocate(m%u(m%len))
            allocate(m%v(m%len))
            allocate(m%du(m%len))
            allocate(m%dv(m%len))
        endif       
        
        ! If outdir then will save: 
        if (present(out_dir)) then
            write(*,*)'Saving mode to: ', trim(out_dir)
            call m%get_mineos_mode(.true., out_dir)
        else
            call m%get_mineos_mode(.false., "")  
        endif

    end function get_mode


end module 