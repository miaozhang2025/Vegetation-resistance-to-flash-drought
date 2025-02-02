module fd_def
contains
subroutine flashdrought(n,a,et,b,drt)
!======================================================================
! References:
! 1. Yuan, X., L. Wang, P. Wu, P. Ji, J. Sheffield, and M. Zhang, 2019: 
! Anthropogenic shift towards higher risk of flash drought over China. 
! Nature Communications, 10, 4661, 
! https://doi.org/10.1038/s41467-019-12692-7
! 2. Yuan, X., Y. Wang, P. Ji, P. Wu, J. Sheffield, J. Otkin, 2023: A 
! global transition to flash droughts under climate change. Science, 
! https://doi.org/10.1126/science.abn6301
!======================================================================
! n, number of pentads for a given year
! a(n), pentad soil moisture percentile [0-100] 
! b(1), number of flash droughts; 
! b(2), mean duration of flash droughts
! b(3), mean severity of flash droughts 
! b(4), mean speed of flash droughts
! drt, 1-drought onset, 2-drought recovery, 0-no drought
!======================================================================
implicit none
real,parameter::thresh=40.   ! drought start threshold [percentile]
real,parameter::thresh1=20.  ! drought end threshold [percentile]
real,parameter::speed=5.     ! drought speed threshold [percentile]
integer,parameter::td=4      ! drought duration threshold [pentads]
integer n,i,j,k
real a(n),et(n) !input
real b(4,46)    !output
real tmin,et0
integer drt(n)
integer si,cnt,flag1,cnt20
integer flag   ! flag=1, begin drought identification
integer flag2  ! flag2=1, onset period; flag2=0, recovery period
real odur(n)   ! total duration of each flash drought
real osev(n)   ! total severity of each flash drought
real ospd(n)   ! onset speed of each flash drought

cnt=0; flag=0; si=0
odur=0; osev=0.; ospd=0.
tmin=99.
flag1=0; flag2=0
et0 = 0.
drt=0

do j=2,n
  if(a(j).gt.0)then
   if(a(j)<thresh)then
      if(flag==0)then                ! drought start
         flag1=1
         if(j>1)then
            if(a(j-1)<thresh) flag1=0
         endif
         if(flag1==1)then
            flag=1
            flag2=1                  ! onset period
            cnt=cnt+1
            si=j                     ! si is drought start time
            osev(cnt)=thresh-a(j)
            odur(cnt)=1
            ospd(cnt)=thresh-a(j)
            drt(j)=1
            tmin=a(j)
            et0 = et(j)
            cnt20 = 0
            if(a(j)<20)   cnt20 = 1
         endif
      else
       if(et(j).lt.0)then
         print*,"error"
         exit
       end if
        if(a(j)<20)then
          cnt20 = cnt20 + 1
        end if 
         if(odur(cnt).ge.8)then
             flag = 0
             if(odur(cnt)<td.or.cnt20.lt.3.or.et0>50.000*odur(cnt).or.tmin>thresh1)then
                osev(cnt)=0.
                odur(cnt)=0
                ospd(cnt)=0.
              drt(si:j-1)=0
                      et0=0.
                      cnt=cnt-1
            end if
            tmin = 99.
         else 
         if(flag2==1)then            ! onset period
            if(thresh-a(j)>=(j-si+1)*speed)then
               if(tmin<=thresh1)then        ! enter thresh1 (e.g., <20%)
                  if(a(j)>=a(j-1))then      ! stop onset period
                     flag2=0
                     if(a(j)<=thresh1)then  ! enter recovery period
                        osev(cnt)=osev(cnt)+(thresh-a(j))
                        odur(cnt)=odur(cnt)+1
                        et0      = (et0 + et(j))
                        drt(j)=2
                        if(a(j)<tmin) tmin=a(j)
                     else                   ! no recovery period
                        flag=0
                        tmin=99.
                        if(odur(cnt)<td.or.cnt20.lt.3.or.et0>50.000*odur(cnt))then
                           osev(cnt)=0.
                           odur(cnt)=0
                           ospd(cnt)=0.
                           drt(si:j-1)=0
                           et0      =0.
                           cnt=cnt-1
                        endif
                     endif
                  else                      ! continue onset period
                     osev(cnt)=osev(cnt)+(thresh-a(j))
                     odur(cnt)=odur(cnt)+1
                     ospd(cnt)=(thresh-a(j))
                     drt(j)=1
                     et0      = (et0 + et(j))
                     if(a(j)<tmin) tmin=a(j)
                  endif
               else                         ! before entering thresh1, continue onset period
                   osev(cnt)=osev(cnt)+(thresh-a(j))
                   odur(cnt)=odur(cnt)+1
                   ospd(cnt)=(thresh-a(j))/real(j-si+1)
                   drt(j)=1
                   et0      = (et0 + et(j))
                   if(a(j)<tmin) tmin=a(j)
               endif
            else               ! speed does not meet      
               flag2=0         ! so enter recovery period
               if(tmin<=thresh1)then
                  if(a(j)<=thresh1)then
                     osev(cnt) = osev(cnt)+(thresh-a(j))
                     odur(cnt) = odur(cnt)+1
                     drt(j)    = 2
                     et0       = (et0 + et(j))!/real(odur(cnt))
                  else            !drought dismiss
                     flag=0
                     tmin=99.
                     if(odur(cnt)<td.or.cnt20.lt.3.or.et0>50.000*odur(cnt))then
                        osev(cnt)=0.
                        odur(cnt)=0
                        ospd(cnt)=0.
                        drt(si:j-1)=0
                        cnt=cnt-1
                     endif
                  endif
               else            ! tmin>thresh1, does not meet drought criterion 
                  drt(si:j-1)=0
                  osev(cnt)=0.
                  odur(cnt)=0
                  ospd(cnt)=0.
                  cnt=cnt-1
                  flag=0
                  tmin=99.
               end if
            endif
         else  ! flag2=0, recovery period
            if(a(j)<=thresh1)then
               osev(cnt)=osev(cnt)+(thresh-a(j))
               odur(cnt)=odur(cnt)+1
               drt(j)=2
               et0      = (et0 + et(j))
            else                  ! drought dismiss
               if(odur(cnt)<td.or.cnt20.lt.3.or.et0>50.000*odur(cnt))then
                   osev(cnt)=0.
                   odur(cnt)=0
                   ospd(cnt)=0.
                   drt(si:j-1)=0
                   cnt=cnt-1
               endif
               flag=0
               tmin=99.
            end  if
         end if  
        end if  !dur<=8
      end if  !yes or not in drought period 
      if(j==n .and. flag==1)then
         if(odur(cnt)<td .or. tmin>thresh1.or.cnt20.lt.3.or.et0>50.000*odur(cnt))then
             osev(cnt)=0.
             odur(cnt)=0
             ospd(cnt)=0.
             drt(si:j)=0
             cnt=cnt-1
         end if
      end if
   else                           ! a(j)>=thresh
      if(flag==1)then             ! drought dismiss
         flag=0
         if(odur(cnt)<td .or. tmin>thresh1.or.cnt20.lt.3.or.et0>50.000*odur(cnt))then
            drt(si:j)=0
            osev(cnt)=0.
            odur(cnt)=0
            ospd(cnt)=0.
            cnt=cnt-1
         endif
         tmin=99.
      endif
   end if
  end if  
end do

if(cnt>0)then
   b(1,1)=real(cnt)
   b(2,1:cnt)=real(odur(1:cnt))!/real(cnt)
   b(3,1:cnt)=real(osev(1:cnt))!/real(odur(1:cnt)))/real(cnt)
   b(4,1:cnt)=real(ospd(1:cnt))!/real(cnt)
else
   b=0
endif
end subroutine flashdrought
end module fd_def
