program test_modes
  ! Program checks a mode and its normalisation 

  implicit none

  integer :: n,l,nl
  real(4)    :: wcom, qmod, r(1000), rho(1000), u(1000), du(1000), v(1000), dv(1000)
  character(len=1) :: type
  integer :: i
  real(4)    :: norm, check_norm

  write(*,*)'Enter n and l'
  read(*,*) n,l
 
  write(*,*)'Enter Type (S or T)'
  read(*,'(a)') type

  write(*,'(a,i4,3x,a,i4)')'Searching for: ', n,type,l
  

  call get_mode(type,n,l,wcom,qmod,nl,r,u,du,v,dv)

  ! Write U, U', V, V' to file 
  open(11,file='mode.data')
  do i = 1, nl
     write(11,*) r(i),u(i),du(i),v(i),dv(i)
  enddo
  close(11)

  write(*,*)'Written to mode.data'

  ! Compute norm and print
  norm = check_norm(type,n,l,wcom,u,v)
  write(*,*) 'norm = ', norm

end program test_modes
