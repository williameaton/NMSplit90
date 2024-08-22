    subroutine get_mode(type,nord,l,wcom,qmod,u,du,v,dv, save_mode)
        !
        ! Usage:
        !    *  input  * 
        !  type --- can be 't','s',  upper case is o.k. (char*1)
        !  nord --- order number                        (integer)
        !  l --- degree number                          (integer)
        ! 
        !    *  output *
        !  wcom --- normalized anglar frequency         (real*8)
        !  qmod --- normalized Q factor                 (real*8)
        !  nnl --- number of layers of output array     (integer)
        !  u(nnl),du(nnl),v(nnl),dv(nnl) ---
        !     in case of toroidal and radial modes, only u,du
        !     is valid, for spheroidal modes, all are valid.
        !
        use params, only: NL, NR,ddir,model_fname, rad_mineos, verbose
        implicit none
        include "constants.h"

        character(len=1)   :: type
        integer            :: nord,l
        real(4)            :: u(NL),du(NL)
        real(4), optional  :: v(NL),dv(NL)
        logical :: save_mode
    
        character(len=200) :: catalogue,bin_file, eigstring
        integer            :: ntype,nvec,i,j
        integer            :: nnl 
        character(len=1)   :: type1,type2,char
        integer            :: ieigtxt,iomod,iocat,iobin,nrec,ios,nn,ll
        real(8)            :: wwmhz,ttcom,ggcom,qqmod,eps, raylquo
        real(4)            :: buf(6*NL),  wcom, qmod,  cg4, phsv
    
        integer(4)         :: n4, l4, n4old, l4old
    
        nnl = NL
        eps = 1.0d-4
    
        ! File IDs 
        iomod   = 9
        iocat   = 10
        iobin   = 11
        ieigtxt = 12
    
    
        ! set up correct catalogue and bin file for mode reading
        if (type == 'T' .or. type == 't') then
        ntype = 1
        catalogue = trim(ddir)//'prem_ani_att_T'
        nvec = 2 * NL
        type1 = 'T';type2 = 't';
    
        bin_file = trim(ddir)//'prem_ani_att_T.bin'
    
        else if (type == 'S' .or. type == 's') then
        type1 = 'S';type2 = 's';
    
        if (l == 0) then
            ! Radial modes
            ntype = 2
            nvec = 2 * NL
            catalogue = trim(ddir)//'prem_ani_att_R'
            bin_file  = trim(ddir)//'prem_ani_att_R.bin'
        else
            ! Spheroidal not radial
            ntype = 3
            catalogue = trim(ddir)//'prem_ani_att_S'
            if (NR == NL) then
                bin_file = trim(ddir)//'prem_ani_att_S.bin'
                nvec = 6 * NL 
            else
                bin_file = trim(ddir)//'prem_ani_att_S.bin'
                nvec = 4 * NL
            endif
        endif ! l=0
        endif ! T or S 
    
    
    
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
    
            ! Exits when nn is nord and l is ll 
            if (nn == nord .and. ll == l) exit
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
        read(iobin,iostat = ios)n4,l4,wcom,qmod, cg4,(buf(i),i=1,nvec)
    
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
            write(*,*)'Angular Freq in rad :', wcom
            write(*,*)'Frequency in mHz    :', wwmhz
            write(*,*)'Period in seconds   :', 2*PI/wcom
            write(*,*)'Q value             :', qmod
            write(*,*)'Group velocity      :', cg4
            write(*,*)
        endif
    
        
        if(save_mode)then
            ! Save eigenfunctions to text file in column format 
            write(eigstring,'(i0,a,i0,a)')n4,type1,l4,'.txt'
            open(ieigtxt,file=trim(eigstring), iostat=ios)
            write(ieigtxt,*)n4,  type1, l4
            do i =1,NR
               if(type1=='S')then 
                ! Radius, U, U', V, V', P, P'
                  write(ieigtxt,*)rad_mineos(i), buf(i), buf(i + NR), buf(i + 2*NR), buf(i + 3*NR), buf(i + 4*NR), buf(i + 5*NR)
               elseif (type1=='T')then 
                ! Radius, W, W' ?
                  write(ieigtxt,*)rad_mineos(i), buf(i), buf(i + NR)
               endif 
            enddo 
            close(ieigtxt)
        endif
    
    
        ! Compare the freq in mHz (wwmhz) to ang freq read from binaries
        ! as well as Q values as a sanity check
        if (abs(wwmhz*2*PI/1000 -wcom) > eps  .or. abs(qmod - qqmod) > eps) then      
        stop 'Error: frequencies or Q values do not match'
        endif 
        
        ! Format buffer into arrays 
        ! Toroidal only
        u(1 : NL)  = buf(1 : NL)
        du(1 : NL) = buf (NL + 1 : 2 * NL)
    
        if (ntype == 2) then
            ! Radial
            v(1 : NL) = 0.
           dv(1 : NL) = 0.
        else if (ntype == 3) then
            ! Spheroidal
            v(1 : NL) = buf (2 * NL + 1 : 3 * NL)
           dv(1 : NL) = buf (3 * NL + 1 : 4 * NL)
        endif
    
    end subroutine get_mode