function check_norm(type,nord,l,wcom,u,v)

  implicit none
  include 'modes.h'
  include 'model.h'

  real(4) :: check_norm
  character(len=*) type
  integer :: nord,l
  real(4) :: wcom, u(NL), v(NL), kernel(NL), norm, rho(NL), angunorm
  
  character(len=150) model_file
  integer :: iomod,junk1,i,j,ios, ndisc, disc(NR)
  real(8) :: radius(NR),  r(NL),  rdisc(NR), rho_all(NR) 
  

   ! reading model file to get the radius and the density
   iomod = 10
   model_file = trim(ddir)//trim(model_fname)
   open(iomod, file = model_file, status = 'old',iostat = ios)
   do i = 1 , NR
      if (ios .eq. 0)  read(iomod,*,iostat=ios)  junk1, radius(i),rho_all(i)
      j = i - NR + NL
      if ( j > 0 ) then 
         r(j) = radius(i) / RA
         rho(j) = real(rho_all(i) / RHOAV, 4)
      endif
   enddo
   close(iomod)


   ! set up kernel array  
   ! See MINEOS manual page 22-23
   kernel = 0.0
   do i = 1, NR
      if (type .eq. 'S' .or. type .eq. 's') then
         kernel(i) = rho(i) * (u(i)*u(i)+v(i)*v(i)) * r(i) * r(i)
      else if (type .eq. 'T' .or. type .eq. 't') then
         ! in this case u stores the W eigenfunctions
         kernel(i) = rho(i) * u(i) * u(i) * r(i) * r(i)
      else
         stop 'Incorrect mode type - should be S or T'
      endif
   enddo


   ! Save eigenfunctions to text file in column format 
   open(75,file='kernel.txt', iostat=ios)
   do i =1,NR
         write(75,*)kernel(i)
   enddo 
   write(*,'(a,/)')'written kernel to kernel.txt'
   close(75)


   ! Find the discontinuties: 
   ! outputs are  disc: the knot IDs of the disconinuities (- side)
   !             rdisc: radius of discontinuities
   !             ndisc: number of discontinuities
   call find_disc(disc,rdisc,ndisc)

   ! Integrate the kernel the entire radius with knowledge of the 
   ! discontinuities
   call intgrl_disc_trapz(norm,NR,r,disc,ndisc,1,NR,kernel)

   ! Normalise - see MINEOS manual 
   angunorm =  real(dsqrt (PI * GRAV * RHOAV),4)
   norm = norm * wcom * wcom / (angunorm * angunorm)
   check_norm = norm

end function check_norm
  

