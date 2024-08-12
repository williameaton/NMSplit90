subroutine get_cmt(cmt_file,yr,jda,ho,mi,sec,t_cmt,hdur, &
                   elat,elon,depth,moment_tensor,event)
    integer :: event,nlen,yr,jda,ho,mi,nevent
    real    :: sec,t_cmt,hdur,elat,elon,depth
    real    :: moment_tensor(6)
    character(len=80) :: cmt_file,temp

    integer :: i,ios,datasource,lstr,mo,da,julian_day
    real :: mb,ms
    character*24 reg
    character*128 string
    integer lnblnk
c
c---- first read hypocenter info
c
     open(1,file=cmt_file,iostat=ios)
     if(ios.ne.0) then
       write(6,"('error opening CMT file ',a80)") cmt_file
       stop
     endif
c      rewind(1)
     read(1,*) nevent

c      if(event.eq.1) then
c      else
c       nlen = (event-1)*13
c       do i=1,nlen
c        read(1,"(a)") temp
c       enddo
c      endif

     read(1,"(a4,i5,i3,i3,i3,i3,f6.2,f9.4,
    #                    f10.4,f6.1,f4.1,f4.1,1x,a)") 
    #         datasource,yr,mo,da,ho,mi,sec,elat,elon,depth,mb,ms,reg
     jda=julian_day(yr,mo,da)
     ios=0
     do while(ios.eq.0) 
     read(1,"(a)",iostat=ios) string
     lstr=lnblnk(string)
     if(string(1:10).eq.'event name') then
     else if(string(1:10).eq.'time shift') then
       read(string(12:lstr),*) t_cmt
     else if(string(1:13).eq.'half duration') then
       read(string(15:lstr),*) hdur
     else if(string(1:8).eq.'latitude') then
       read(string(10:lstr),*) elat
     else if(string(1:9).eq.'longitude') then
       read(string(11:lstr),*) elon
     else if(string(1:5).eq.'depth') then
       read(string(7:lstr),*) depth
     else if(string(1:3).eq.'Mrr') then
       read(string(5:lstr),*) moment_tensor(1)
     else if(string(1:3).eq.'Mtt') then
       read(string(5:lstr),*) moment_tensor(2)
     else if(string(1:3).eq.'Mpp') then
       read(string(5:lstr),*) moment_tensor(3)
     else if(string(1:3).eq.'Mrt') then
       read(string(5:lstr),*) moment_tensor(4)
     else if(string(1:3).eq.'Mrp') then
       read(string(5:lstr),*) moment_tensor(5)
     else if(string(1:3).eq.'Mtp') then
       read(string(5:lstr),*) moment_tensor(6)
     endif
     enddo
     close(1)
c
c       scale the moment-tensor
c
     do i=1,6
       moment_tensor(i)=moment_tensor(i)*1.0E-30
     enddo
     return
end subroutine get_cmt





integer function julian_day(yr,mo,da)
integer yr,mo,da
c
integer mon(12)
integer lpyr
data mon/0,31,59,90,120,151,181,212,243,273,304,334/
c
julian_day=da+mon(mo)
if(mo.gt.2) julian_day=julian_day+lpyr(yr)
return
end

integer function lpyr(yr)
integer yr
c
c---- returns 1 if yr is a leap year
c
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
end function julian_day
