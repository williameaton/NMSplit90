subroutine decipher_cluster(cluster,nmodes,ndim,n,type,l,om,alpha,om0,alpha0)
  
  include "coupling.h"

  integer :: nmodes, ndim, n(MMAX), type(MMAX), l(MMAX)
  real    :: om(MMAX), alpha(MMAX), om0, alpha0
  character(len=80) cluster

  integer :: lcluster,m1,m2,lmode,ind,i
  integer ::lnblnk,idigit
  real    :: q, gv, av, ahs, aht
  real    :: r(NR),u(NR),du(NR),v(NR),dv(NR),w(NR),dw(NR),p(NR),dp(NR)
  character(len=10) mode
  character(len=1) ctype

  lcluster=lnblnk(cluster)
  nmodes=0
  m1=1
  m2=m1+2
  om0=0.0
  alpha0=0.0
  ndim=0
  do while(m2.le.lcluster)
    nmodes=nmodes+1
    do while((cluster(m2:m2).ne.'-').and.(m2.le.lcluster))
        m2=m2+1
    enddo
    m2=m2-1
    mode=cluster(m1:m2)
    lmode=lnblnk(mode)
    if(index(mode(1:lmode),'S').gt.0) then
      type(nmodes)=1
      ind=index(mode(1:lmode),'S')
    else
      type(nmodes)=0
      ind=index(mode(1:lmode),'T')
    endif
    if(ind.eq.0) then
        write(6,"('Unknown mode type')")
        call exit(1)
    endif
    n(nmodes)=0
    do i=1,ind-1
        n(nmodes)=n(nmodes)+idigit(mode(i:i))*10**(ind-1-i)
    enddo
    l(nmodes)=0
    do i=ind+1,lmode
        l(nmodes)=l(nmodes)+idigit(mode(i:i))*10**(lmode-i)
    enddo
    ! use new prem and catalogue in the datalib
    ! call readmode(n(nmodes),type(nmodes),l(nmodes),om(nmodes), &
    !              q,gv,av,ahs,aht,r,u,du,v,dv,w,dw,p,dp)

    if(type(nmodes) .eq. 0) then
      ctype = 't'
    else
      ctype = 's'
    endif

    call get_mode(n(nmodes),ctype,l(nmodes),om(nmodes), &
                  q,av,ahs,aht,r,u,du,v,dv,w,dw,p,dp)
    call latest_fq(n(nmodes),type(nmodes),l(nmodes),om(nmodes),q)
    if(l(nmodes).gt.LMAX) then
      write(*,*)'angular degree of mode exceeds LMAX = ', LMAX
      call exit(1)
    endif
    alpha(nmodes)=0.5*om(nmodes)/q
    om0=om0+om(nmodes)
    alpha0=alpha0+alpha(nmodes)
    ndim=ndim+(2*l(nmodes)+1)
    m1=m2+2
    m2=m1+2
  enddo
  if(nmodes.gt.MMAX) then
    write(*,*)'number of modes in cluster exceeds M', MMAX
    call exit(1)
  endif
  om0=om0/nmodes
  alpha0=alpha0/nmodes
  return
end subroutine decipher_cluster

! ______________________________________________________________________

integer function idigit(ch)

  character(len=1) ch
  if(ch.eq.'0') then
    idigit=0
    return 
  else if(ch.eq.'1') then
    idigit=1
    return 
  else if(ch.eq.'2') then
    idigit=2
    return 
  else if(ch.eq.'3') then
    idigit=3
    return 
  else if(ch.eq.'4') then
    idigit=4
    return 
  else if(ch.eq.'5') then
    idigit=5
    return 
  else if(ch.eq.'6') then
    idigit=6
    return 
  else if(ch.eq.'7') then
    idigit=7
    return 
  else if(ch.eq.'8') then
    idigit=8
    return 
  else if(ch.eq.'9') then
    idigit=9
    return 
  else
    write(*,*) 'not an integer'
    call exit(1)
  endif
  return
  end function

! ______________________________________________________________________

  subroutine latest_fq(n,type,l,om,q)
!
!       use the latest center frequency measurements of Masters
!       and quality factors predicted by QL6 of Durek & Ekstrom (1996)
!
  include "coupling.h"

  integer :: n, type, l
  real    :: om, q
  integer :: ios, nn, ll
  real    :: f, qq
  logical ::  notfound
!

  if(type.eq.0) then ! toroidal modes
    open(3,file='modecat_dir/toroidal_fq',iostat=ios)
    rewind(3)
    notfound=.true.
    do while(ios.eq.0.and.notfound)
      read(3,*,end=10) nn,ll,f,qq
      if(nn.eq.n.and.ll.eq.l) then
        om=PI2*f*1.0E-03
        q=qq
        notfound=.false.
      endif
    enddo
10   continue
    close(3)
  else ! spheroidal modes
    open(4,file='modecat_dir/spheroidal_fq',iostat=ios)
    rewind(4)
    notfound=.true.
    do while(ios.eq.0.and.notfound)
      read(4,*,end=20) nn,ll,f,qq
      if(nn.eq.n.and.ll.eq.l) then
        om=PI2*f*1.0E-03
        q=qq
        notfound=.false.
      endif
    enddo
20   continue
    close(4)
  endif
return
end subroutine
