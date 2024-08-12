!
!    source calculates s = M : E^* where M is the moment tensor and
!    E^* is the complex conjugate of the strain tensor.
!    Input are the overtone number n, the type (1 for S, 0 for T),
!    the angular order l, the latitude and longitude of the source in 
!    degrees and the depth in km, and the six elements of the moment tensor
!    such that moment_tensor(1)=M_rr, moment_tensor(2)=M_tt,
!    moment_tensor(3)=M_pp, moment_tensor(4)=M_rt, moment_tensor(5)=M_rp,
!    moment_tensor(6)=M_tp, 
!    the moment tensor is expressed with equivalent str,dip,rake
!    Output is s, such that s(1) corresponds to m=0, s(2) to m=1
!    and s(l+1) to m=l.
!
subroutine kernel(n,type,l,lat,lon,depth,Mo,str,dip,rake,s_kmo,s_kstr,s_kdip,s_krake)
    
    ! NOTE I HAVE ADDED THIS BUT SOME ITEMS I DONT YET KNOW WHAT THEIR TYPE IS 
    ! SEE WHAT IS PARSED IN TO kernel 
    implicit none 

    include "coupling.h"

    integer :: n,type,l
    real    :: lat,lon,depth,moment_tensor(6),dmo(6),dstr(6),ddip(6),drake(6)
    complex :: s(1),kmo(1),kstr(1),kdip(1),krake(1)

    integer :: m,i,j
    real    :: om,q,gv,av,ahs,aht
    real    :: r(NR),u(NR),du(NR),v(NR),dv(NR),w(NR),dw(NR),p(NR),dp(NR)
    real    :: bl,rs
    real    :: us,dus,vs,dvs,ws,dws,xs,zs
    real    :: theta,phi,sint,cost,cosect,cott
    real    :: x(LMAX+1),dx(LMAX+1)
    complex ::  expphi,expmmphi,e(6)
    character(1) ::  ctype


!----- use new prem layer model and catalogues in the datalib
    if(type .eq. 0) then
    ctype = 't'
    else
    ctype = 's'
    endif

    ! NOTE MAYBE SOME ISSUES WITH REAL(4) vs 8 now 
    call get_mode(n,ctype,l,om,q,av,ahs,aht,r,u,du,v,dv,w,dw,p,dp)

    om=om/sqrt(PI*GRAV*RHOAV)
    bl=l*(l+1.0)
!
    rs=1.0-depth*1000.0/RA ! non-dimensionalized source radius
    i=NR

    do while(r(i).gt.rs) ! find model radius r(i) just below source radius rs
        i=i-1
    enddo

    if(r(i).eq.rs) then ! rs coincides with a model radius
        us=u(i)
        dus=du(i)
        vs=v(i)
        dvs=dv(i)
        ws=w(i)
        dws=dw(i)
        xs=dvs+(us-vs)/rs
        zs=dws-ws/rs
    else ! interpolate to find us, vs, etc. at the source radius
        call interpolate(r(i),u(i),du(i),r(i+1),u(i+1),du(i+1),rs,us,dus)
        call interpolate(r(i),v(i),dv(i),r(i+1),v(i+1),dv(i+1),rs,vs,dvs)
        call interpolate(r(i),w(i),dw(i),r(i+1),w(i+1),dw(i+1),rs,ws,dws)
        xs=dvs+(us-vs)/rs
        zs=dws-ws/rs
    endif

    theta  = PI/2.0 - lat*PI/180.0
    phi    = lon*PI/180.0
    sint   = sin(theta)
    cost   = cos(theta)
    cosect = 1.0/sint
    cott   = cost/sint

    call lgndr(l,cost,sint,x,dx)
    expphi=cexp(cmplx(0.0,-phi))
    expmmphi=cmplx(1.0,0.0)

    dip  = dip*PI/180.0
    str  = str*PI/180.0
    rake = rake*PI/180.0

!---  compute kernel for moment

    dmo(1) = sin(2.0*dip)*sin(rake)
    dmo(2) = sin(dip)*cos(rake)*sin(2.0*str)-sin(2.0*dip)*sin(rake)*(cos(str)**2.0)
    dmo(3) = -1.0*(sin(dip)*cos(rake)*sin(2.0*str)+sin(2.0*dip)*sin(rake)*(sin(str)**2.0))
    dmo(4) = cos(dip)*cos(rake)*sin(str) - cos(2.0*dip)*sin(rake)*cos(str)
    dmo(5) = -1.0*(cos(dip)*cos(rake)*cos(str)+cos(2.0*dip)*sin(rake)*sin(str))
    dmo(6) = -1.0*(sin(dip)*cos(rake)*cos(2.0*str)+0.5*sin(2.0*dip)*sin(rake)*sin(2.0*str))

!---  compute kernel for rake

    drake(1) = sin(2.0*dip)*cos(rake)
    drake(2) = -1.0*(sin(dip)*sin(rake)*sin(2.0*str)+sin(2.0*dip)*cos(rake)*cos(str)**2.0)
    drake(3) = sin(dip)*sin(rake)*sin(2.0*str)-sin(2.0*dip)*cos(rake)*(sin(str)**2.0)
    drake(4) = -1.0*(cos(dip)*sin(rake)*sin(str)+cos(2.0*dip)*cos(rake)*cos(str))
    drake(5) = cos(dip)*sin(rake)*cos(str)-cos(2.0*dip)*cos(rake)*sin(str)
    drake(6) = sin(dip)*sin(rake)*cos(2.0*str)-0.5*sin(2.0*dip)*cos(rake)*sin(2.0*str)

!---- compute kernel for dip

    ddip(1) = 2.0*cos(2.0*dip)*sin(rake)
    ddip(2) = cos(dip)*sin(rake)*sin(2.0*str) + sin(2.0*dip)*cos(rake)*(cos(str)**2.0)
    ddip(3) = -1.0*(cos(dip)*cos(rake)*sin(2.0*str)+2.0*cos(2.0*dip)*(sin(rake)*sin(str)**2.0))
    ddip(4) = -1.0*(sin(dip)*cos(rake)*sin(str)-2.0*sin(2.0*dip)*sin(rake)*cos(str))
    ddip(5) = sin(dip)*cos(rake)*cos(str)+2.0*sin(2.0*dip)*sin(rake)*sin(str)
    ddip(6) = -1.0*(cos(dip)*cos(rake)*cos(2.0*str)+cos(2.0*dip)*sin(rake)*sin(2.0*str))

!---  compute kernel for str

    dstr(1) = 0.0
    dstr(2) = 2.0*sin(dip)*cos(rake)*cos(2.0*str)+sin(2.0*dip)*sin(rake)*sin(2.0*str)
    dstr(3) = -1.0*(2.0*sin(dip)*cos(rake)*cos(2.0*str)+sin(2.0*dip)*sin(rake)*sin(2.0*str))
    dstr(4) = cos(dip)*cos(rake)*cos(str)+cos(2.0*dip)*sin(rake)*sin(str)
    dstr(5) = cos(dip)*cos(rake)*sin(str)-cos(2.0*dip)*sin(rake)*cos(str)
    dstr(6) = 2.0*sin(dip)*cos(rake)*sin(2.0*str)-sin(2.0*dip)*sin(rake)*cos(2.0*str)

    do m=0,l
        i=m+1
        e(1)=cmplx(dus*x(i),0.0)
        e(2)=cmplx((us*x(i)-vs*(cott*dx(i)-((m*cosect)**2.0-bl)*x(i)))/rs,0.0) + & 
             cmplx(0.0,-m*ws*cosect*(dx(i)-cott*x(i))/rs)
        e(3)=cmplx((us*x(i)+vs*(cott*dx(i)-((m*cosect)**2.0)*x(i)))/rs,0.0) - &
             cmplx(0.0,-m*ws*cosect*(dx(i)-cott*x(i))/rs)
        e(4)=cmplx(xs*dx(i),0.0)+cmplx(0.0,-m*zs*cosect*x(i))
        e(5)=cmplx(0.0,-m*xs*cosect*x(i))-cmplx(zs*dx(i),0.0)
        e(6)=cmplx(0.0,-2.0*m*vs*cosect*(dx(i)-cott*x(i))/rs) + &
             cmplx(ws*(2.0*cott*dx(i)-(2.0*((m*cosect)**2.0)-bl)*x(i))/rs,0.0)

        s_kmo(i)   = cmplx(0.0,0.0)
        s_krake(i) = cmplx(0.0,0.0)
        s_kdip(i)  = cmplx(0.0,0.0)
        s_str(i)   = cmplx(0.0,0.0)

        do j=1,6
            k_mo(i)=kmo(i)+cmplx(dmo(j),0.0)*e(j)
            k_rake(i)=krake(i)+cmplx(drake(j),0.0)*e(j)
            k_dip(i)=kdip(i)+cmplx(ddip(j),0.0)*e(j)
            k_str(i)=kstr(i)+cmplx(dstr(j),0.0)*e(j)
        enddo

        k_mo(i)   = k_mo(i)*expmmphi
        k_rake(i) = Mo*k_rake(i)*expmmphi
        k_dip(i)  = Mo*k_dip(i)*expmmphi
        k_str(i)  = Mo*k_str(i)*expmmphi
        expmmphi  = expmmphi*expphi
    enddo

    return
end subroutine kernel




subroutine interpolate(r1,f1,df1,r2,f2,df2,rs,fs,dfs)
    !     interpolate finds the value of a function fs and its derivative dfs
    !     at the radius rs by interpolating given values of the function and its
    !     derivative just below the source at radius r1, given by f1 and df1,
    !     and values just above the source ar radius r2, given by f2 and df2.
    !     Note that we should have r2 > rs > r1.
    implicit none 

    real :: r1,f1,df1,r2,f2,df2,rs,fs,dfs

    real :: h,dr

    h=rs-r1
    dr=r2-r1
    if((h.lt.0.0).or.(dr.lt.0.0)) then
        write(6,"('wrong input in interpolate')")
        call exit(1)
    endif

    fs=f1+h*df1+0.5*h*h*(df2-df1)/dr
    dfs=df1+h*(df2-df1)/dr
    return
end subroutine interpolate
