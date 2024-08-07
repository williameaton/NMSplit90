
! TODO: 
! Parses value nnl but then sets it inside the subroutine? Could 
! just use the NL global value? 

! ADd in a safety mechanism that it doesnt keep looping/reading
! past the last mode - becasue it will just use the last mode's content



subroutine get_mode(type,nord,l,wcom,qmod,nnl,r,u,du,v,dv)
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
   !  r(nnl) --- normalized radius of each layer   (rest are all real*8)
   !  u(nnl),du(nnl),v(nnl),dv(nnl) ---
   !     in case of toroidal and radial modes, only u,du
   !     is valid, for spheroidal modes, all are valid.
   !
   implicit none
   include "modes.h"

   character(len=1)   :: type
   integer            :: nord,l,nnl
   real               :: radius(NR),r(NL),u(NL),du(NL)
   real, optional     :: v(NL),dv(NL)

   character(len=200) ::  model_file,catalogue,bin_file, ddir, eigstring
   integer            ::  ntype,nvec,reclen,i,j
   character(len=1)   :: type1,type2,char
   integer            :: ieigtxt,iomod,iocat,iobin,nrec,ios,nn,ll,nnn,lll,junk1,junk2
   real(8)            :: av,ah,bb,wwmhz,ttcom,ggcom,qqmod,eps, raylquo
   real(4)            :: buf(6*NL),  wcom, qmod,av4,ah4, rn4, vnorm4, anorm4,  cg4, phsv

   integer(4)         :: n4, l4, n4old, l4old


   nnl = NL
   eps = 1.0d-4

   ! File IDs 
   iomod = 9
   iocat = 10
   iobin = 11
   ieigtxt = 12

   ddir = '/Users/eaton/Documents/Software/NMSplit90/databases/prem_ani_att_database/'


   ! read model file to get radius array (normalised 0 to 1)
   model_file = trim(ddir)//'model'

   write(*,'(a,/)')' Reading Model file ....'
   open(iomod, file = model_file, status = 'old',iostat = ios)
   ! Error check
   if (ios .ne. 0) stop 'Error reading model file'

   do i = 1 , NR
      if (ios .eq. 0)  read(iomod,*,iostat=ios)  junk1, radius(i)
      j = i - NR + NL 
      if ( j > 0 ) r(j) = radius(i) / RA
   enddo
   close(iomod)

   write(*,'(a,/)')' Finished reading model file.'


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

   write(*,*)'Made it here!'


   ! get the record number of the desired mode from catalogue file
   write(*,*) ' Reading Catalogue file .... ', catalogue
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
   
   write(*,*) ' Record Number of ',nord,char,l,'is:  ',nrec
   write(*,'(a,2g16.7)')  '  Catalogue file:', bb, qqmod
   

   !reclen = 12 * 4+  nvec*4 
   reclen = (5*4) +  nvec*4 


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

   write(*,*) 
   write(*,*)'n', n4
   write(*,*)'l', l4
   write(*,*)'Angular Freq in mHz', wcom, 'so Period (s) is ', 2*PI/wcom
   write(*,*)'Q (qmod)', qmod
   write(*,*)'Group velocity(cg4)', cg4
   write(*,*)'wwmhz', wwmhz

   ! Save eigenfunctions to text file in column format 
   write(eigstring,'(i0,a,i0,a)')n4,type1,l4,'.txt'
   write(*,*)' EIGstring: ', trim(eigstring)
   open(ieigtxt,file=trim(eigstring), iostat=ios)
   write(ieigtxt,*)n4,  type1, l4
   do i =1,NR
      if(type1=='S')then 
         write(ieigtxt,*)r(i), buf(i), buf(i + NR), buf(i + 2*NR), buf(i + 3*NR), buf(i + 4*NR), buf(i + 5*NR)
      elseif (type1=='T')then 
         write(ieigtxt,*)r(i), buf(i), buf(i + NR), buf(i + 2*NR), buf(i + 3*NR)
      endif 
   enddo 
   close(ieigtxt)

   

   !if (junk1 /= junk2 .or. junk1 /= reclen) then
   !   print *,junk1,junk2,reclen
   !   stop 'Incorrect reading'
   !endif
   !write(*,300) nnn,char,lll,wcom,qmod
   !300 format('  Binary Record:',i4,a2,i4,'  --   wcom:',g16.7, &
   !           '   Q : ', g16.7)
   ! Double checking the read values are consistent from the text and binary files



   ! Compare the freq in mHz (wwmhz) to ang freq read from binaries
   ! as well as Q values as a sanity check
   if (abs(wwmhz*2*PI/1000 -wcom) > eps  .or. abs(qmod - qqmod) > eps) then      
      stop 'Error reading 2'
   endif 
   
   ! Format buffer into arrays 
   u(1 : NL) = buf(1 : NL)
   du(1 : NL) = buf (NL + 1 : 2 * NL)
   
   if (ntype == 2) then
      v(1 : NL) = 0.
      dv(1 : NL) = 0.
   else if (ntype == 3) then
      v(1 : NL) = buf (2 * NL + 1 : 3 * NL)
      dv(1 : NL) = buf (3 * NL + 1 : 4 * NL)
   endif

end subroutine get_mode
   


