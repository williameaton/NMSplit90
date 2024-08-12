subroutine get_mode(nord,type,l,omega,q,aver,ahors,ahort,r,u,du,v,dv,w,dw,p,dp)
    implicit none

    include "coupling.h"

    character(len=1) type
    integer :: nord,l
    real :: omega,q,aver,ahors,ahort
    real :: r(NL),u(NL),du(NL),v(NL),dv(NL),w(NL),dw(NL),p(NL),dp(NL)

    real :: wcom,qmod,radius(NR)
    character(len=200)  model_file,catalogue,bin_file
    integer ::  ntype,nvec,reclen,i,j
    character(len=1) type1,type2,char
    integer :: iomod,iocat,iobin,nrec,ios,nn,ll,nnn,lll,junk1,junk2
    real :: av,ah,bb,wwmhz,ttcom,ggcom,qqmod,eps
    real :: k,oms
    real :: buf(6*NL)

    eps = 1.0d-3

    iomod = 9
    iocat = 10
    iobin = 11

    ! read model radius
    model_file = '/opt/Modes/nmcat/aniprem808.dk'
    !  print *,' Reading Model file ....'
    open(iomod, file = model_file, status = 'old',iostat = ios)
    do i = 1 , NR
       if (ios .eq. 0)  read(iomod,*,iostat=ios)  junk1, radius(i)
       j = i - NR + NL
       if ( j > 0 ) r(j) = sngl(radius(i) / RA)
    enddo
    close(iomod)
    if (ios .ne. 0) stop 'Error reading model file'

    ! set up correct catalogue and bin file for mode reading
    if (type .eq. 'T' .or. type .eq. 't') then
       ntype = 1
       catalogue = '/opt/Modes/nmcat/toroidal.aniprem808_40s'
       if (NR .eq. NL) then
          bin_file = '/opt/Modes/nmcat/toroidal.aniprem808_40s_all.bin'
       else
          bin_file = '/opt/Modes/nmcat/toroidal.aniprem808_40s.bin'
       endif
       nvec = 2 * NL
       type1 = 'T';type2 = 't';
       else if (type .eq. 'S' .or. type .eq. 's') then
          if (l .eq. 0) then
             ntype = 2
             catalogue = '/opt/Modes/nmcat/radial.aniprem808_40s'
             if (NR .eq. NL) then
                bin_file = '/opt/Modes/nmcat/radial.aniprem808_40s_all.bin'
             else
                bin_file = '/opt/Modes/nmcat/radial.aniprem808_40s.bin'
             endif
             nvec = 2 * NL
          else
             ntype = 3
             catalogue = '/opt/Modes/nmcat/spheroidal.aniprem808_40s'
             if (NR .eq. NL) then
                bin_file = '/opt/Modes/nmcat/spheroidal.aniprem808_40s_all.bin'
                nvec = 6 * NL
             else
                bin_file = '/opt/Modes/nmcat/spheroidal.aniprem808_40s.bin'
                nvec = 4 * NL
             endif
          endif
          type1 = 'S';type2 = 's';
       endif
       reclen = 2 * 4 + 4 * 8 + nvec * 4

       ! get the record number of the desired mode from catalogue file
       !print *, ' Reading Catalogue file ....'
       nrec = 0
       open(iocat,file=catalogue,status='old',iostat=ios)
        do while (ios .eq. 0)
          nrec = nrec + 1
          read(iocat,'(i5,a2,i5,5g16.7)',iostat=ios) nn,char,ll,bb,wwmhz,ttcom,ggcom,qqmod
             if (nrec .eq. 1 .and. (char .ne. type1 .and. char .ne. type2)) then
                print *,nn,ll,';',char,';',type1,';',type2
                stop 'Incorrect mode catalogue'
             endif
             if (nn .eq. nord .and. ll .eq. l) exit
        enddo

        close(iocat)
        if (ios .gt. 0) stop 'Error reading 1 '
        if (ios .lt. 0) stop 'Mode not found in the catalogue'

        open(iobin,file=bin_file,status='old',iostat=ios,form='unformatted',& 
!   ----- modify to fit the system in reading binary file
!          open(iobin,file=bin_file,status='old',iostat=ios,form='binary',
           access='direct',recl=reclen+2*4)

        if (ios .ne. 0) stop 'Error opening catalogue file'
        read(iobin,rec = nrec,iostat = ios)&
            junk1,nnn,lll,wcom,qmod,av,ah,(buf(i),i=1,nvec),junk2
        close(iobin)
        if (junk1 .ne. junk2 .or. junk1 .ne. reclen) then
           print *,junk1,junk2,reclen
           stop 'Incorrect reading'
        endif
        if (abs(wcom-bb) > eps  .or. abs(qmod - qqmod) > eps) then
          print*, nrec,wcom,bb,qmod,qqmod
          stop 'Error reading 2'
        endif

        k = dsqrt(dble(l*(l+1)))
        oms =  wcom / sqrt(PI*GRAV*RHOAV)
        do i = 1,NL
          if (ntype .eq. 1) then ! toroidal modes
            w(i) = sngl(oms * buf(i) / k)
            dw(i) = sngl(oms * buf (i + NL) / k)
            u(i) = 0.
            du(i) = 0.
            v(i) = 0.
            dv(i) = 0.
            p(i)= 0.
            dp(i)= 0.
          else if (ntype .eq. 2) then ! radial modes
            w(i) = 0.
            dw(i) = 0.
            u(i) = sngl(oms * buf(i))
            du(i) = sngl(oms * buf (i + NL))
            v(i) = 0.
            dv(i) = 0.
          else if (ntype .eq. 3) then ! spheroidal modes
            w(i) = 0.
            dw(i) = 0.
            u(i) = sngl(oms * buf(i))
            du(i) = sngl(oms * buf (i + NL))
            v(i) = sngl(oms * buf (i + 2 * NL) / k)
            dv(i) = sngl(oms * buf (i + 3 * NL) / k)
            p(i)= sngl(oms * buf (i + 4*NL) / k )
            dp(i)= sngl(oms * buf (i + 5*NL) / k)
          endif
        enddo

        omega = sngl(wcom)
        q = sngl(qmod)

        if(ntype .eq. 1) then
          ahort = sngl(ah * oms)
          ahors = 0.
          aver = 0.
        else if(ntype .eq. 2) then
          ahort = 0.
          ahors = 0.
          aver = sngl(av * oms)
        else if(ntype .eq. 3) then
          ahort = 0.
          ahors = sngl(ah * oms)
          aver = sngl(av * oms)
        endif

      return
      end

subroutine intgrl(sum,r,nir,ner,f,s1,s2,s3)
    ! Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
    ! radii values as in model PREM_an222

    include "coupling.h"

    integer :: i,j

    real :: r(NR), f(NR), yprime(NR), kdis(28), s1(NR), s2(NR), s3(NR)
    !dimension r(NR),f(NR),yprime(NR),kdis(28),s1(NR),s2(NR),s3(NR)
!-----    use new prem layer model (808 layer) and catalogue in the datalib
!        data kdis/55,109,114,175,180,184,193,202,209,213,216,219,16*0/,ndis/12/,n/NR/

  ! Original definitions: 
  !      data kdis/165,330,355,630,650,660,695,740,770,785,795,805,16*0/,ndis/12/,n/NR/
    integer, parameter ::  ndis = 12
    integer, parameter ::  n    = NR

    kdis(1)  = 165
    kdis(2)  = 330
    kdis(3)  = 355
    kdis(4)  = 630
    kdis(5)  = 650
    kdis(6)  = 660
    kdis(7)  = 695
    kdis(8)  = 740
    kdis(9)  = 770
    kdis(10) = 785
    kdis(11) = 795
    kdis(12) = 805
    kdis(13:28) = 0

    third = 1.0/3.0
    fifth = 1.0/5.0
    sixth = 1.0/6.0

    call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)

    nir1 = nir + 1
    sum = 0.0
    do i=nir1,ner
            j = i-1
            rji = r(i) - r(j)
            
        ! ORIGINAL by JT: 
        ! sum=sum+r(j)*r(j)*rji*(f(j)+rji*(.50*s1(j)+rji*(third*s2(j)+rji*
        !       .250*s3(j))))+2.0*r(j)*rji*rji*(.50*f(j)+rji*(third*s1(j)+rji*
        !      (.250*s2(j)+rji*fifth*s3(j))))+rji*rji*rji*(third*f(j)+rji*
        !         (.250*s1(j)+rji*(fifth*s2(j)+rji*sixth*s3(j))))
        sum = sum + r(j) * r(j) * rji * (f(j) + rji*(0.50*s1(j) + rji*(third*s2(j) + rji*0.25*s3(j)) )) &
                    + 2.0 * r(j) * rji * rji * (0.50*f(j) + rji*(third*s1(j) + rji*(0.250*s2(j) + rji*fifth*s3(j)) )) &
                    + rji * rji * rji * (third*f(j) + rji*(.250*s1(j) + rji*(fifth*s2(j) + rji*sixth*s3(j)) ))
        enddo
        
    return
    end



subroutine deriv(y,yprime,n,r,ndis,kdis,s1,s2,s3)
    !     slightly altered dp version of subroutine rspln.
    real :: y(n),yprime(n),r(n),f(3,1000),kdis(28),yy(3)
    real :: s1(n),s2(n),s3(n)

    ! NOTE: equivalence is still legal but maybe change? 
    equivalence (yy(1),y0)


    yy(1:3) = 0.0 

    ndp=ndis+1
    do 3 nd=1,ndp
    if(nd.eq.1) go to 4
    if(nd.eq.ndp) go to 5
    j1=kdis(nd-1)+1
    j2=kdis(nd)-2
    go to 6
  4 j1=1
    j2=kdis(1)-2
    go to 6
  5 j1=kdis(ndis)+1
    j2=n-2
  6 if((j2+1-j1).gt.0) go to 11
    j2=j2+2
    y0=(y(j2)-y(j1))/(r(j2)-r(j1))
  s1(j1)=yy(1)
  s1(j2)=yy(1)
  s2(j1)=yy(2)
  s2(j2)=yy(2)
  s3(j1)=yy(3)
  s3(j2)=yy(3)
    go to 3
 11 a0=0.0
    if(j1.eq.1) go to 7
    h=r(j1+1)-r(j1)
    h2=r(j1+2)-r(j1)
    y0=h*h2*(h2-h)
    h=h*h
    h2=h2*h2
    b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
    go to 8
  7 b0=0.0
  8 b1=b0
    if(j2 .gt. 1000)write(0,'("error:deriv:j2= ",i5)')j2
    do 1 i=j1,j2
    h=r(i+1)-r(i)
    y0=y(i+1)-y(i)
    h2=h*h
    ha=h-a0
    h2a=h-2.0*a0
    h3a=2.0*h-3.0*a0
    h2b=h2*b0
    s1(i)=h2/ha
    s2(i)=-ha/(h2a*h2)
    s3(i)=-h*h2a/h3a
    f(1,i)=(y0-h*b0)/(h*ha)
    f(2,i)=(h2b-y0*(2.0*h-a0))/(h*h2*h2a)
    f(3,i)=-(h2b-3.0*y0*ha)/(h*h3a)
    a0=s3(i)
  1 b0=f(3,i)
    i=j2+1
    h=r(i+1)-r(i)
    y0=y(i+1)-y(i)
    h2=h*h
    ha=h-a0
    h2a=h*ha
    h2b=h2*b0-y0*(2.*h-a0)
    s1(i)=h2/ha
    f(1,i)=(y0-h*b0)/h2a
    ha=r(j2)-r(i+1)
    y0=-h*ha*(ha+h)
    ha=ha*ha
    y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
    s3(i)=(y0*h2a+h2b)/(h*h2*(h-2.0*a0))
!     the following statements were expanded to prevent register overflow on unix
    s13=s1(i)*s3(i)
    s2(i)=f(1,i)-s13
    do 2 j=j1,j2
    k=i-1
    s32=s3(k)*s2(i)
    s1(i)=f(3,k)-s32
    s21=s2(k)*s1(i)
    s3(k)=f(2,k)-s21
    s13=s1(k)*s3(k)
    s2(k)=f(1,k)-s13
  2 i=k
    s1(i)=b1
    j2=j2+2
  s1(j2)=yy(1)
  s2(j2)=yy(2)
  s3(j2)=yy(3)
  3 continue
    do 20 i=1,n
 20 yprime(i)=s1(i)
    return
end


 real function thrj(j1,j2,j3,m1,m2,m3)
    real ::  y(100)      

!     evaluation of Wigner 3-j coefficients
    thrj=0.
    
    if(j1+j2-j3.lt.0.or.j2+j3-j1.lt.0.or.j3+j1-j2.lt.0) return
    if(j1-iabs(m1).lt.0) return
    if(j2-iabs(m2).lt.0) return
    if(j3-iabs(m3).lt.0) return
    if(m1+m2+m3.ne.0) return
    
c-----
!     use symmetries to make j3 largest of j1,j2,j3

    jc=max(j1,j2,j3)
    
    if(jc.eq.j3) then
    ja=j1
    jb=j2
    jc=j3
    ma=m1
    mb=m2
    mc=m3
    elseif(jc.eq.j2) then
    ja=j3
    jb=j1
    jc=j2
    ma=m3
    mb=m1
    mc=m2
    else
    ja=j2
    jb=j3
    jc=j1
    ma=m2
    mb=m3
    mc=m1
    endif
    
    lm2=-jb
    if(ja+mc-jb.lt.0) lm2=-mc-ja
    lm1=-mc-lm2
    m=lm1+jb+mc+1
    if(ja-jb-mc.lt.0) m=lm1+ja+1
    
    y(1)=0.
    y(2)=1.
    ss=1.
    if(m.eq.1) goto 20
    
    do 10 n=2,m
    mx=lm1-n+2
    my=lm2+n-2
    alpha=sqrt(real((ja-mx+1)*(ja+mx)*(jb+my+1)*(jb-my)))
    beta=real(jc*(jc+1)-ja*(ja+1)-jb*(jb+1)-2*mx*my)
    gamma=sqrt(real((ja+mx+1)*(ja-mx)*(jb-my+1)*(jb+my)))
    y(n+1)=(beta*y(n)-gamma*y(n-1))/alpha
    ss=ss+y(n+1)*y(n+1)
10    continue

20    n=lm1-ma+2
    thrj=y(n)*((-1)**(-ja+jb+mc))/sqrt(ss*real(2*jc+1))
    
    return
    end




  subroutine rotylm(alpha,beta,gamma,s,mp,t,bmp)
!       This subroutine calculates bmp = Dmpm as defined by
!	equation (4.1.12) of Edmonds (1960).
!       alpha, beta and gamma are the Euler angles (radians).
!	s = angular degree l, t = azimuthal order m, and mp = m'.
!	The relation between Y_l^m in the geographical (i.e. not
!	rotated) reference frame and Y'_l^m' in the equatorial (i.e.
!	rotated) reference frame is
c
!	    Y_l^m (geographic) = Sum_m' D_m'm Y'_l^m' (equatorial)
c
!	After the rotation, the source-receiver path is along the
!	equator (theta=pi/2) with the source at phi=0 and receiver
!	at phi=DELTA.
!       The real and imaginary parts of Dmpm are returned in bmp(2).
  common/pjaco/pjac0,pjac1
  dimension bmp(2)
  integer*4 s,t,facsym,tsave

  i1=mp-t
  i2=mp+t
!       For negative values of i1 or i2, must use symmetry properties
!       of Edmonds, 1960, eqn's (4.2.5) and (4.2.6).
  facsym=1
  mpsave=mp
  tsave=t
  if(i1.ge.0.and.i2.ge.0) go to 10
  if(i1.lt.0.and.i2.lt.0) go to 20
  if(i1.lt.0.and.i2.ge.0) go to 30
  mpp=mp
  mp=-t
  t=-mpp
  i1=mp-t
  i2=mp+t
  facsym=1
  go to 10
30      mpp=mp
  mp=t
  t=mpp
  i1=mp-t
  i2=mp+t
  facsym=(-1)**i1
  go to 10
20      mp=-mp
  t=-t
  i1=mp-t
  i2=mp+t
  facsym=(-1)**i1
10      i3=s-mp
  i4=s+mp
  i5=s+t
  i6=s-t
  cosb=cos(beta)
  call fact(i4,i5,a1)
  call fact(i3,i6,a2)
  anorm=a1+a2
!	First determine n=0 and n=1 values of pjac.
!	Then recurrence relations will get actual value.
  n0=0
  i1s=i1
  i2s=i2
  xs=cosb 
  call pnab(n0,i1s,i2s,anorm,xs,pjac)
  if(i3.eq.0) go to 11
  n1=1
  i1s=i1
  i2s=i2
  xs=cosb 
  call pnab(n1,i1s,i2s,anorm,xs,pjac)
  if(i3.eq.1) go to 11
  call pnab(i3,i1,i2,anorm,cosb,pjac)
!	Guard against underflow.
11	continue 
  cb2=cos(beta/2.)
  sb2=sin(beta/2.)
  if(cb2.eq.0.0.or.sb2.eq.0.0) go to 98
  dmag1=abs(float(mp+t))*alog(abs(cos(beta/2.))) 
  dmag2=abs(float(mp-t))*alog(abs(sin(beta/2.))) 
  if(dmag1.lt.-30.0.or.dmag2.lt.-30.0) go to 101
98	dmpm=(cos(beta/2.))**(mp+t)*(sin(beta/2))**(mp-t)*pjac
  arg=float(mpsave)*gamma+float(tsave)*alpha
  bmp(1)=cos(arg)*dmpm*float(facsym)
  bmp(2)=sin(arg)*dmpm*float(facsym)
  go to 100
99      write(6,*)'t or mp > s'
  stop
101	bmp(1)=0.
  bmp(2)=0.
100     mp=mpsave
  t=tsave
  return
  end





  subroutine pnab(n,nalpha,nbeta,anorm,x,pjac)
!       Returns Jacobi polynomial P(n,[alpha,beta])(x) -- Edmonds, 1960, p. 58.
!	multiplied by [(s+mp)!(s-mp)!/(s+t)!(s-t)!]**0.5
  common/pjaco/pjac0,pjac1
  if(n.gt.1) go to 12	
  flog2=float(n)*alog(2.)
  ix=1
  if(x.ge.0) go to 3
!       Use symmetry relation (4.1.20)
  ix=(-1)**n
  nbeta1=nbeta
  nbeta=nalpha
  nalpha=nbeta1
  x=-x
3       sum=0.
  np1=n+1
  do 10 i=1,np1
  nu=i-1
  npa=n+nalpha
  npb=n+nbeta
  nmnu=n-nu
  npamnu=n+nalpha-nu
  nupb=nu+nbeta
  call fact(npa,npamnu,a1)
  n1=1
  call fact(n1,nu,a2)
  call fact(npb,nupb,a3)
  call fact(n1,nmnu,a4)
  opx=1.+x
  omx=1.-x
  if(opx.eq.0.0) opx=1.e-8
  if(omx.eq.0.0) omx=1.e-8
  dmag1=float(nu)*(alog(opx))
  dmag2=float(n-nu)*(alog(omx))
  xlg=0.5*anorm+a1+a2+a3+a4+dmag1+dmag2-flog2
  if(xlg.lt.-20.0) go to 10
  sum=sum+exp(xlg)*float((-1)**(n-nu))
10      continue
  pjac=sum*float(ix)
!	  write(6,*)'pnab: n,nalpha,nbeta,x,pjac=',n,nalpha,nbeta,x,pjac
  if(n.eq.0) pjac0=pjac
  if(n.eq.1) pjac1=pjac 
  return
12	continue
!	Use recurrance relations.
  apb=float(nalpha+nbeta)
  a2b2=float(nalpha*nalpha-nbeta*nbeta)  
  fn0=0.
  fn1=1.
  n2=1
13	n2=n2+1
  fn2=float(n2)
  a1n=2.*fn2*(fn1+apb+1.)*(2.*fn1+apb)
  a2n=(2.*fn1+apb+1.)*(a2b2)
  a3n=(2.*fn1+apb)*(2.*fn1+apb+1.)*(2.*fn1+apb+2.)
  a4n=2.*(fn1+float(nalpha))*(fn1+float(nbeta))*(2.*fn1+apb+2.)
  pjac=((a2n+a3n*x)/a1n)*pjac1 - (a4n/a1n)*pjac0
  if(n2.eq.n) return  
  fn0=fn1
  fn1=fn2
  pjac0=pjac1
  pjac1=pjac
  go to 13
  end




  subroutine fact(n1,n2,a)
!       Finds log[n1(factorial)/n2(factorial)], returned in a.
  m1=n1
  m2=n2
  if(n1.eq.0) n1=1
  if(n2.eq.0) n2=1
  if(n2.gt.n1) go to 10
  a=-alog(float(n2))
  m=n1-n2+1
  do 5 i=1,m
5       a=a+alog(float(n2+i-1))
  go to 20
10      b=alog(float(n1))
  m=n2-n1+1
  do 15 i=1,m
15      b=b-alog(float(n1+i-1))
  a=b
20      n1=m1
  n2=m2
  return
  end





    subroutine lgndr(l,c,s,x,dx)
c
!    computes legendre function x(l,m,theta)
!    theta=colatitude,c=cos(theta),s=sin(theta),l=angular order,
!    sin(theta) restricted so that sin(theta).ge.1.e-7
!    x(1) contains m=0, x(2) contains m=1, x(k+1) contains m=k
!    m=azimuthal(longitudenal) order 0.le.m.le.l
!    dx=dx/dtheta
!    subroutine originally came from Physics Dept. Princeton through
!    Peter Davis
!    modified to run stably on the Perkin-Elmer by Jeffrey Park 11/23/85
c
!     implicit real*8(a-h,o-z)
    implicit integer*4(i-n)
    dimension x(2),dx(2)
    data tol/1.e-05/,rfpi/.282094791773880/,root3/1.73205080756890/
    data boeing/.707106781186550/
    if(s.ge.1.0-tol) s=1.0-tol
    lsave=l
    if(l.lt.0) l=-1-l
    if(l.gt.0) go to 1
    x(1)=rfpi
    dx(1)=0.0
    l=lsave
    return
  1 if(l.ne.1) go to 2
    c1=root3*rfpi
    c2=boeing*c1
    x(1)=c1*c
    x(2)=-c2*s
    dx(1)=-c1*s
    dx(2)=-c2*c
    l=lsave
    return
  2 sos=s
    if(s.lt.tol) s=tol
    cot=c/s
    ct=2.0*c
    ss=s*s
    lp1=l+1
    g3=0.0
    g2=1.0
    f3=0.0
!  evaluate m=l value, sans (sin(theta))**l
    do 100 i=1,l
100 g2=g2*(1.0-1.0/(2.0*i))
    g2=rfpi*sqrt((2*l+1)*g2)
    f2=l*cot*g2
    x(lp1)=g2
    dx(lp1)=f2
    w=0.0
    v=1.0
    y=2.0*l
    z=y+1.0
    d=sqrt(v*y)
    t=0.0
    mp1=l
    m=l-1
!  these recursions are similar to ordinary m-recursions, but since we
!  have taken the s**m factor out of the xlm's, the recursion has the powers
!  of sin(theta) instead
  3 g1=-(ct*mp1*g2+ss*t*g3)/d
    f1=(mp1*(2.0*s*g2-ct*f2)-t*ss*(f3+cot*g3))/d-cot*g1
    x(mp1)=g1
    dx(mp1)=f1
    if(m.eq.0) go to 4
    mp1=m
    m=m-1
    v=v+1.0
    y=y-1.0
    t=d
    d=sqrt(v*y)
    g3=g2
    g2=g1
    f3=f2
    f2=f1
     go to 3
  4 maxsin=-72.0/log10(s)
!  maxsin is the max exponent of sin(theta) without underflow
    lpsafe=min0(lp1,maxsin)
    stom=1.0
    fac=sign(1.0,(l/2)*2-l+.50)
!  multiply xlm by sin**m
    do 5 m=1,lpsafe
    x(m)=fac*x(m)*stom
    dx(m)=fac*dx(m)*stom
    stom=stom*s
  5 continue
!  set any remaining xlm to zero
    if(maxsin.le.l) then
    mmm=maxsin+1
    do 200 m=mmm,lp1
    x(m)=0.0
200 dx(m)=0.0
    endif
    s=sos
    l=lsave
    return
    end








    real function splbasis(x,x0,dx,np,i)
c
    nx=np-1
    xdiff=x-x0
    interval=1+int((x-x0)/dx)
    if(abs(xdiff-float(nx)*dx).lt.0.0000001) interval=nx
    xd=x-x0-float(interval-1)*dx
    h=1./dx
    hsq=1./dx**2
    hcu=1./dx**3
    xdsq=xd**2
    xdcu=xd**3
c
c---- return the value of the i-th basis element
c
    value=0.
    if(i.eq.0) then
      if(interval.eq.1) then
        value=0.25*hcu*xdcu-1.5*h*xd+1.5
      else if(interval.eq.2) then
        value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
      else
        value=0.
      endif
    else if(i.eq.1) then
      if(interval.eq.1) then
        value=-0.5*hcu*xdcu+1.5*h*xd
      else if(interval.eq.2) then
        value=0.75*hcu*xdcu-1.5*hsq*xdsq+1.
      else if(interval.eq.3) then
        value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
      else
       value=0.
      endif
    else if(i.gt.1.and.i.lt.nx-1) then
      if(interval.eq.i-1) then
        value=0.25*hcu*xdcu
      else if(interval.eq.i) then
        value=-0.75*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
      else if(interval.eq.i+1) then
        value=0.75*hcu*xdcu-1.5*hsq*xdsq+1.
      else if(interval.eq.i+2) then
        value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
      else
        value=0.
      endif
    else if(i.eq.nx-1) then
      if(interval.eq.nx-2) then
        value=0.25*hcu*xdcu
      else if(interval.eq.nx-1) then
        value=-0.75*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
      else if(interval.eq.nx) then
        value=0.5*hcu*xdcu-1.5*hsq*xdsq+1.
      else
        value=0.
      endif
    else if(i.eq.nx) then
      if(interval.eq.nx-1) then
        value=0.25*hcu*xdcu
      else if(interval.eq.nx) then
        value=-0.25*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
      else
        value=0.
      endif
    endif
    splbasis=value
    return
    end

