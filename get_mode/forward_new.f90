program forward_new
  
  include "coupling.h"

  integer :: smax,model_type,ntclusters,nevents,source
  integer :: ncomponents(EMAX),nstations,components(EMAX,3),nsources
  integer :: nclusters(EMAX,3,NWIN),ntime_windows(EMAX,3)
  integer :: clusterp !,nobs_cluster(CMAX)
  integer :: dummy,lnblnk,idisc(ND),npar,no
  integer :: nps,yr,jda,ho,mi,yr1,jda1,ho1,mi1
  integer :: s20rts_type !,smax_s20rts,latest_type,smax_latest
  real :: stele,stbur,timemin,f_max,sec,sec1,spsec,theta,phi
  real :: model0(NMMAX),model(NMMAX),dmodel(NMMAX)
  real :: t_cmt(SSMAX),hdurs(SSMAX),elats(SSMAX),elons(SSMAX),depths(SSMAX)
  real :: moment_tensors(6,SSMAX),dt0(SSMAX)
  real :: rlats(NSMAX),rlons(NSMAX)
  real :: stazi(NSMAX),stdip(NSMAX),spsecs(NSMAX),ntfs(NSMAX)
  real :: time_windows(EMAX,3,NWIN,2),freq_windows(EMAX,3,NWIN,CMAX,2)
  real :: noise_windows(EMAX,3,NWIN,CMAX,2)
  real :: mantle_basis(NR,0:NK),radius(NR)
  character(80) :: stns(NSMAX),chns(NSMAX),netw
  character(80) :: events(EMAX)
  character(80) :: clusters(EMAX,3,NWIN,CMAX)
  character(80) :: tclusters(CMAX),output1,output2
  logical :: splitting,even,threeD,latest
  integer :: ios,datasource,lstr,mo,da,julian_day,ntshift,ntshift1,ntshift2
  real :: mb,ms
  character(len=24) ::reg,eventname
  character(len=128)::  string

  integer :: i,j,p,qq,event,component,cluster,tcluster,time_window,station
  integer :: nmodes,n(CMAX),type(CMAX),l(CMAX)
  integer :: ndim
  integer :: k,kp
  integer :: ntf,nts,nsf,nss,npad,log2npad,ipadmin,info
  integer :: itime_window_l,itime_window_r,ifreq_window_l,ifreq_window_r
  integer :: LDA,LDVL,LDVR,LWORK,IPIV(MMAX*(2*LMAX+1))
  integer :: lstn,lcha
  real :: om0,alpha0,om(CMAX),alpha(CMAX)
  real :: rlat,rlon,nu(3),t0,dtf,dom_max
  real :: tsf,tss,df
  real :: seismogram(NTFMAX)
  real :: taper_l,taper_r
  character(len=80):: stn,chn,clust
  character(len=80):: station_file,cmt_file,filen
  character(len=1) :: answer
  complex :: WORK(2*MMAX*(2*LMAX+1))
  real :: RWORK(2*MMAX*(2*LMAX+1))
  complex :: dom(MMAX*(2*LMAX+1))
  complex :: dom0(MMAX*(2*LMAX+1))
  complex :: C(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1)),temp(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: QR(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: QL(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: Qs(CMAX,MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: QRs(CMAX,MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: QLs(CMAX,MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: doms(CMAX,MMAX*(2*LMAX+1))
  complex :: Cp(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: Cm(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: dC(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: dCp(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: dCps(CMAX,MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: Q(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: Qp(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: Qm(MMAX*(2*LMAX+1),MMAX*(2*LMAX+1))
  complex :: slow(NTSMAX)
  complex :: r(MMAX*(2*LMAX+1))
  complex :: s(MMAX*(2*LMAX+1))
  complex :: rp(MMAX*(2*LMAX+1))
  complex :: sp(MMAX*(2*LMAX+1))
  logical :: accurate,add,sks12,bs12,new_cluster,s20rts

!---------------------------------------------------------------------
! forward parameter setup
! including what mantle model will be used, the largest angular degree
!---------------------------------------------------------------------
  write(*,*)
  write(*,*)"type of inversion"
  write(*,*)
  write(*,*)"    S: 1"
  write(*,*)"S & B: 2:"

  read(*,*)model_type
  if((model_type.lt.1).or.(model_type.gt.2)) then
    write(6,"('wrong input value for inversion type')")
    call exit(1)
  endif


  write(6,"(' ')")
      write(6,"('maximum angular degree (even, between 2 and 20):')")
      read *, smax
      if((smax.gt.20).or.(smax.lt.2)) then
        write(6,"('wrong input value for maximum angular degree')")
        call exit(1)
      endif

      write(6,"(' ')")
      write(6,"('even degree model?')")
      read(*,'(a)') answer

      if(answer.eq.'y') then
        even=.true.
      endif
      if((model_type*(smax+1)*(smax+1)*(smax+2)/2).gt.NMMAX) then
         write(6,"('length of model vector exceeds ',i6)") NMMAX
         call exit(1)
      endif
     
      if(answer.eq.'n') then 
        even=.false.
      endif
      if((model_type*(smax+1)*(smax+1)*(smax+1)).gt.NMMAX) then
         write(6,"('length of model vector exceeds ',i6)") NMMAX
         call exit(1)
      endif

      write(6,"(' ')")
      write(6,"('use s20rts mantle model as starting model?')")
      read *,answer
      if(answer.eq.'y') then
        s20rts=.false.
        sks12=.true.
        bs12=.false.
        write(6,"(' ')")
      else
        s20rts=.false.
        sks12=.false.
        bs12=.false.
      endif

      write(6,"(' ')")
      write(6,"('use 3-D model or just 1D prem ?')")
      read *,answer
      if(answer.eq.'y') then
       threeD=.true.
      else
       threeD=.false.
      endif
      write(6,"(' ')")

      write(6,"(' ')")
      write(6,"('consider rotation/ellipticity or not')")
      read *,answer
      if(answer.eq.'y') then
       splitting=.true.
      else
       splitting=.false.
      endif
      write(6,"(' ')")
!---------------------------------------------------------------------------
! read the receiver information
!---------------------------------------------------------------------------
       station_file = 'RECORDHEADERS'
       open(25,file = station_file,status='old')
       read(25,*) nstations

       do station = 1,nstations
        read(25,"(a8,1x,a5,1x,a8,1x,f8.4,1x,f9.4,1x,f6.1,1x,f6.1,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4,1x,i3,1x,i2,1x,i2,1x,f6.3)") & 
                stns(station),netw,chns(station),rlats(station),rlons(station),stele,stbur,&
                stazi(station),stdip(station),spsecs(station),nps,yr1,jda1,ho1,mi1,sec1
       enddo

!---------------------------------------------------------------------------
! read the source information, which could contain multiple point sources
! for the same earthquake, with different hdur and t_cmt
!---------------------------------------------------------------------------
        nevents = 1
        cmt_file='CMTSOLUTION'
        open(8,file = cmt_file)
        read(8,*) nsources
        do source = 1,nsources
         print *,'source: ',source,'   line 1'
         read(8,101)&
              datasource,yr,mo,da,ho,mi,sec,elats(source),elons(source),depths(source)&
             ,mb,ms,reg
         print *,&
              datasource,yr,mo,da,ho,mi,sec,elats(source),elons(source),depths(source)&
              ,mb,ms,reg
         jda=julian_day(yr,mo,da)
         do i=1,12
         print *,'source: ',source,'   line 2'
         read(8,"(a)") string
         lstr = lnblnk(string)
         print *,string(1:lstr)
         if(string(1:10).eq.'event name') then
         else if(string(1:10).eq.'time shift') then
         print *,'source: ',source,'   line 3'
          read(string(12:lstr),*) t_cmt(source)
         else if(string(1:13).eq.'half duration') then
         print *,'source: ',source,'   line 4'
          read(string(15:lstr),*) hdurs(source)
         else if(string(1:8).eq.'latitude') then
         print *,'source: ',source,'   line 5'
          read(string(10:lstr),*) elats(source)
         else if(string(1:9).eq.'longitude') then
         print *,'source: ',source,'   line 6'
          read(string(11:lstr),*) elons(source)
         else if(string(1:5).eq.'depth') then
         print *,'source: ',source,'   line 7'
          read(string(7:lstr),*) depths(source)
         else if(string(1:3).eq.'Mrr') then
         print *,'source: ',source,'   line 8'
          read(string(5:lstr),*) moment_tensors(1,source)
         else if(string(1:3).eq.'Mtt') then
         print *,'source: ',source,'   line 9'
          read(string(5:lstr),*) moment_tensors(2,source)
         else if(string(1:3).eq.'Mpp') then
         print *,'source: ',source,'   line 10'
          read(string(5:lstr),*) moment_tensors(3,source)
         else if(string(1:3).eq.'Mrt') then
         print *,'source: ',source,'   line 11'
          read(string(5:lstr),*) moment_tensors(4,source)
         else if(string(1:3).eq.'Mrp') then
         print *,'source: ',source,'   line 12'
          read(string(5:lstr),*) moment_tensors(5,source)
         else if(string(1:3).eq.'Mtp') then
         print *,'source: ',source,'   line 13'
          read(string(5:lstr),*) moment_tensors(6,source)
         endif
         enddo
         print *,'source: ',source,'   line 14'
!
!       scale the moment-tensor
!
         do i=1,6
          moment_tensors(i,source)=moment_tensors(i,source)*1.0E-30
         enddo

!      
!       centroid time for each point source
!
        sec = sec+t_cmt(source)-hdurs(source)
!        dt0(source) = sec - sec1
        dt0(source) = t_cmt(source)
        print*,'t_start=',dt0(source)

       enddo
       close(8)

! ----harvard format ----
100     format(a4,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1,f4.1,f4.1,1x,a)
! --- chen's  format ----
101     format(a3,i5,i3,i3,i3,i3,f6.2,f7.3,f7.3,f6.1,f4.1,f4.1,1x,a)
!---------------------------------------------------------------------
! read setup file
!---------------------------------------------------------------------
        ntclusters=0  ! total number of independent clusters
        filen = 'setup'
        do event = 1,nevents
        open(3,file=filen)
        rewind(3)

      
        read(3,"(i1)") ncomponents(event)

        if(ncomponents(event).lt.1.or.ncomponents(event).gt.3) then
          write(6,"('wrong number of components (1--3)')")
          call exit(1)
        endif
        do component=1,ncomponents(event)

          read(3,"(i1)") components(event,component)
          read(3,"(i1)") ntime_windows(event,component)
 

          if(ntime_windows(event,component).gt.NWIN) then
            write(6,*)'number of time windows exceeds',  NWIN
            call exit(1)
          endif
          write(6,"(' ')")
          write(6,"('     time windows:')")
          write(6,"(' ')")

          do time_window=1,ntime_windows(event,component)

            read(3,"(f9.1,1x,f9.1)") taper_l,taper_r

            timemin=2.0*taper_r
            time_windows(event,component,time_window,1)=taper_l
            time_windows(event,component,time_window,2)=taper_r

            read(3,"(i3)") nclusters(event,component,time_window)

            write(6,"('          ',f5.1,' --',f5.1,' hours (',i2,' clusters)')")&
                      taper_l,taper_r,nclusters(event,component,time_window)
            write(6,"(' ')")
            write(6,"('                      cluster   f1    f2')")
            write(6,"(' ')")
            f_max=0.0 ! keep track of highest frequency window in clusters

            do cluster=1,nclusters(event,component,time_window)
              read(3,"(a45,1x,4f6.3)") clust,&
                     freq_windows(event,component,time_window,cluster,1),  &
                     freq_windows(event,component,time_window,cluster,2),  &
                     noise_windows(event,component,time_window,cluster,1), & 
                     noise_windows(event,component,time_window,cluster,2)
              i=1
              do while(clust(i:i).eq.' ')
                i=i+1
              enddo
              clusters(event,component,time_window,cluster)=clust(i:lnblnk(clust))
              write(6,"(a45,1x,2f6.3)")&
               clusters(event,component,time_window,cluster)&
                 (1:lnblnk(clusters(event,component,time_window,cluster))),&
                     freq_windows(event,component,time_window,cluster,1),&
                     freq_windows(event,component,time_window,cluster,2)
              if(freq_windows(event,component,time_window,cluster,2).gt.f_max)then
                f_max=freq_windows(event,component,time_window,cluster,2)
              endif 
              new_cluster=.true.
              clusterp=0

              do while(new_cluster.and.(clusterp.lt.ntclusters))
                clusterp=clusterp+1
                if(tclusters(clusterp).eq.clusters(event,component,time_window,cluster))then 
                    new_cluster=.false.
                endif
              enddo

              if(new_cluster) then
                ntclusters=ntclusters+1
                tclusters(ntclusters)=&
                     clusters(event,component,time_window,cluster)&
                      (1:lnblnk(clusters(event,component,time_window,cluster)))
              endif
            enddo

            write(6,"(' ')")
            if(components(event,component).eq.1) then
              write(6,"('        vertical component records for this window:')")
            elseif(components(event,component).eq.2) then
              write(6,"('    longitudinal component records for this window:')")
            elseif(components(event,component).eq.3) then
              write(6,"('      transverse component records for this window:')")
            elseif(components(event,component).eq.4) then
              write(6,"('             N-S component records for this window:')")
            elseif(components(event,component).eq.5) then
              write(6,"('             E-W component records for this window:')")
            endif
            write(6,"(' ')")
           enddo
          dummy=unlink(filen)
          enddo
         close(3)
        enddo

!---------------------------------------------------------------------
!  read mantle model vector
!---------------------------------------------------------------------
!----- choose 3-D or 1-D mantle model
!      choose initial 3-D model or latest 3-D model

      s20rts_type = 1
!      smax = 12

      if (threeD) then
           call zero_mantle_model(mantle_basis,idisc,radius)
!           call s20rts_mantle_model(smax,s20rts_type)
!           call latest_mantle_model(smax,model_type)
           call  mantle_model(smax,model_type,even,bs12)
           call initialize_model_vector(model0,model,dmodel,smax,npar,model_tpye,even)
           call decipher_model_vector(model,smax,model_type,even)
      else
           call zero_mantle_model(mantle_basis,idisc,radius)
           call initialize_model_vector(model0,model,dmodel,smax,npar,model_tpye,even)
           print*,npar
           call decipher_model_vector(model,smax,model_type,even)
      endif

      accurate=.true.
      dom_max=0.0

!---------------------------------------------------------------------
! compute splitting matrix C, eigenvalue dom and 
! left,right eigenvector QL, QR and renormalized Q
!---------------------------------------------------------------------
      do tcluster=1,ntclusters

        call decipher_cluster(tclusters(tcluster),nmodes,ndim,n,type,l, &
                             om,alpha,om0,alpha0)
        call initialize_splitting(tcluster,nmodes,n,type,l,om,om0,smax, &
                                 mantle_basis,idisc)
        call build_splitting_matrix(tcluster,nmodes,n,type,l,om,alpha,  &
                                   om0,alpha0,smax,C,Q,ndim)

        no = MMAX*(2*LMAX+1)
        LDVL = no
        LDVR = no
        LDA = no
        LWORK = no*2

! ----- compute the righthand eigenvector QR

       call cgeev('V','V',no,C,no,dom,QL,no,QR,no,WORK,2*no,RWORK,info)


!------ compute the lefthand eigenvector QL, conj(QR)^-1
!        call cgetrf(no,no,QR,no,IPIV,info)
!        do i = 1,no
!          do j = 1,no
!           QL(i,j) = QR(i,j)
!         enddo
!       enddo
!       call cgetri(no,QL,no,IPIV,WORK,2*no,info)
        do i = 1,no
         do j = 1,no
          if(i.ne.j) then
            QL(i,j) = QR(j,i)
          else
            QL(i,j) = QR(i,j)
         endif
         enddo
        enddo

!---------------------------------------------------------------------------
        do i=1,ndim
          doms(tcluster,i)=dom(i)
          if(real(dom(i)).gt.dom_max) dom_max=real(dom(i))
          do j=1,ndim
            QLs(tcluster,i,j)=QL(i,j)
            QRs(tcluster,i,j)=QR(i,j)
            Qs(tcluster,i,j)=Q(i,j)
          enddo
        enddo

      print*,'cluster no:',tcluster

      enddo
!---------------------------------------------------------------------------
! compute synthetics for each event/component/time window/station
!---------------------------------------------------------------------------

        do event=1,nevents
          do component=1,ncomponents(event)
            do time_window=1,ntime_windows(event,component)
              do station=1,nstations

                 rlat=rlats(station)
                 rlon=rlons(station)
                 stn=stns(station)
                 chn=chns(station)
               write(6,*) 'compute station ',stn,' component ',chn

               theta=(90.0+stdip(station))*PI2/360.0
               phi=stazi(station)*PI2/360.0
!              vertical component of seismometer
               nu(1)=-cos(theta)
!              N-S component of seismometer
               nu(2)=sin(theta)*cos(phi)
!              E-W component of seismometer
               nu(3)=-sin(theta)*sin(phi)


!-----------   multiple sources ------------
               do source=1,nsources

                t0=dt0(source)
                spsec=spsecs(station)
                dtf=1.0/spsec
                ntshift = int(t0/dtf)
          
                ntf=int(time_windows(event,component,time_window,2)*3600.0/dtf)

!                if (ntf.ne.nps) then
                if(ntf.gt.nps) then
                  print*,ntf,nps
                  call exit(1)
                endif

                call get_start(ntf,nts,nsf,nss,t0,dtf,tsf,tss,dom_max,accurate)

                print*,'source=',source
                do cluster=1,nclusters(event,component,time_window)
                    if(cluster.eq.1.and.source.eq.1) then
                      add=.false.
                    else
                      add=.true.
                    endif
                    call find_tcluster(clusters(event,component,time_window,cluster),&
                                      tclusters,tcluster)
                    call decipher_cluster(tclusters(tcluster),nmodes,ndim,n,type,l,&
                                         om,alpha,om0,alpha0)
                    do i=1,ndim
                      dom(i)=doms(tcluster,i)
                      dom0(i) = cmplx(0.0,0.0)

                      do j=1,ndim
                        QL(i,j)=QLs(tcluster,i,j)
                        QR(i,j)=QRs(tcluster,i,j)
                        Q(i,j)=Qs(tcluster,i,j)
                      enddo
                    enddo

!                    do i=1,ndim
!                     do j=1,ndim
!                      temp(i,j) = temp(i,j)+ QL(i,j)*QR(j,i)
!                      print*,'i=',i,' j=',j,' QL*QR=', temp(i,j)
!                     enddo
!                    enddo
!---------------------------------------------------------------------------------
! compute synthetics from source vector, receiver vector and splitting matrix
!---------------------------------------------------------------------------------
                    call build_s(nmodes,n,type,l,elats(source),elons(source),&
                                depths(source),moment_tensors(1,source),s)
                    call build_sp(ndim,QL,Q,s,sp)
                    call build_r(nmodes,n,type,l,rlat,rlon,nu,r)
                    call build_rp(ndim,QR,Q,r,rp)
                    if (splitting) then
                     call slow_series(ndim,nsf,tsf,ntf,dtf,nts,dom,rp,sp,slow)
                    else
                     call slow_series(ndim,nsf,tsf,ntf,dtf,nts,dom0,r,s,slow)
                    endif
                    call time_series(nsf,tsf,ntf,dtf,nts,om0,alpha0,slow,&
                                    hdurs(source),seismogram,add,ntshift)
                enddo
               enddo
!------------------------------------------------------------------------
!      output the file name
!------------------------------------------------------------------------
                output1=' '
                lstn=lnblnk(stn)
                output1(1:lstn)=stn(1:lstn)
                lnew=lstn
                lcha=lnblnk(chn)
                if(lcha.gt.0) then
                 output1(lnew+1:lnew+lcha+1)='.'//chn(1:lcha)
                endif

                open(21,file=output1,iostat=ios,status='unknown')

                rewind(21)
                do i=1,ntf
                  write(21,*) i*dtf,seismogram(i)/SCALE_FAC
                 enddo
                close(21)
              enddo
            enddo
          enddo
        enddo

end program forward_new




subroutine find_tcluster(cluster,tclusters,tcluster)
      include "coupling.h"

      integer tcluster
      character(80) cluster,tclusters(CMAX)

      tcluster=1
      do while(cluster.ne.tclusters(tcluster))
        tcluster=tcluster+1
      enddo
      return
end subroutine find_tcluster




integer function julian_day(yr,mo,da)
  integer yr,mo,da

  integer mon(12)
  integer lpyr
  data mon/0,31,59,90,120,151,181,212,243,273,304,334/
  julian_day=da+mon(mo)
  if(mo.gt.2) julian_day=julian_day+lpyr(yr)
  return
  end
  integer function lpyr(yr)
  integer yr
!
!---- returns 1 if yr is a leap year
!
  lpyr=0
  if(mod(yr,400).eq.0) then
    lpyr=1
  else if(mod(yr,4).eq.0) then
    lpyr=1
    if(mod(yr,100).eq.0) then
      lpyr=0
    endif
  endif
  return
end function