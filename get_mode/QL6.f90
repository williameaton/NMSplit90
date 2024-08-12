subroutine q_ql6(n,type,l,q)
    ! calculate the quality factor q of a mode for radial Q model QL6
    ! of Durek & Ekstrom (1996)
    include "coupling.h"
    
    integer :: n,type,l
    real    :: q

    character(1) :: ctype
    integer :: i
    real :: omega,om,gv,av,ahs,aht
    real :: r(NR),u(NR),du(NR),v(NR),dv(NR),w(NR),dw(NR),p(NR),dp(NR)
    real :: bl,rin,f,x,z,kernel(NR)
    real :: kappa,mu,k_kappa,k_mu
    real :: rho(NR),epsilon(NR),eta(NR),g(NR)
    real :: c01(NR),c02(NR),c03(NR),c04(NR),c05(NR),Qkappa(NR),Qmu(NR)
    real :: s1(NR),s2(NR),s3(NR)

    ! Common named list block (prem) with variables 
    common/prem/rho,epsilon,eta,g,c01,c02,c03,c04,c05,Qkappa,Qmu
  
    ! Determine mode type and get the mode 
    if(type .eq. 0) then
    ctype = 't'
    else
    ctype = 's'
    endif
    call get_mode(n,ctype,l,omega,q,av,ahs,aht,r,u,du,v,dv,w,dw,p,dp)

    ! Non dimensionalise the frequency
    om=omega/sqrt(PI*GRAV*RHOAV)

    ! k in DT98
    bl=l*(l+1.0)

    do i=2,NR
        rin=1.0/r(i)
        f=(2.0*u(i)-bl*v(i))*rin
        x=dv(i)+(u(i)-v(i))*rin
        z=dw(i)-w(i)*rin
        k_kappa=(du(i)+f)*(du(i)+f)
        k_mu=(1.0/3.0)*(2.0*du(i)-f)*(2.0*du(i)-f)+bl*(x*x+z*z) + &
             bl*(bl-2.0)*(v(i)*v(i)+w(i)*w(i))*rin*rin   
        kappa=(1.0+(2.0/(PI*Qkappa(i)))*log(omega/PI2))* &
              (4.0*c03(i)+c01(i)-4.0*c04(i))/9.0
        if(Qmu(i).ne.0.0) then
            mu=(1.0+(2.0/(PI*Qmu(i)))*log(omega/PI2))* &
                (c03(i)+3.0*c02(i)+c01(i)+2.0*c04(i)-6.0*c05(i))/15.0
        else
            mu=0.0
        endif

        if(Qmu(i).ne.0.0) then
            kernel(i)=kappa*k_kappa/Qkappa(i)+mu*k_mu/Qmu(i)
        else
            kernel(i)=kappa*k_kappa/Qkappa(i)
        endif
    enddo

    kernel(1)=kernel(2)
    
    call intgrl(q,r,1,NR,kernel,s1,s2,s3)
    q=(om*om)/q
    return
end
