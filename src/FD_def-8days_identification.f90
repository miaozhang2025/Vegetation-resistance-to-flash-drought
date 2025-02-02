program drt_identification
use netcdf
use fd_def , only : flashdrought
use sd_def , only : slowdrought
implicit none
integer,parameter::time=46,nyr=2022-2001+1,wlon=1440,wlat=720,dims=3
character*255 ext
integer ncid,varid(7),lon_dimid,lat_dimid,cha_dimid,year_dimid,time_dimid,count_dimid,ncid1,varid1
integer :: start(4),count(4),dimids(dims)
integer y,t,i,j,d,cont,bb,y2,indexs(wlon,wlat,time,nyr),yy,yy2,l,t1,t2
character*4 cyr
character*1 cl
real,allocatable::data(:,:,:,:,:),data2(:,:,:,:),data3(:,:,:,:),et(:,:,:,:),et1d(:,:,:),sm1d(:,:,:)
real a,aa,dep
real::sm(wlon,wlat,1,time,nyr),sm2(wlon,wlat,time,nyr),smp(wlon,wlat,time,nyr),smp2(wlon,wlat,time,nyr),b_fd(wlon,wlat,4,time,nyr+1),b_sd(wlon,wlat,4,time,nyr+1),tvar(wlon,wlat,20,nyr),lon(wlon),lat(wlat)
real::b1(wlon,wlat,4,time,nyr),b2(4,time),dd(47),tt(47),fre_fd(wlon,wlat),fre_sd(wlon,wlat)
integer::drt_fd(wlon,wlat,time*nyr),drt_sd(wlon,wlat,time*nyr),fre_fd1,fre_fd2,drt0(47)
integer::mask(wlon,wlat),sos(wlon,wlat),eos(wlon,wlat),maxid(wlon,wlat)

    call check( nf90_open('/data/workspace/mzhang/analysis/FD_impact/global_mask_025deg.nc',NF90_NOWRITE,ncid1))
    call check( nf90_inq_varid(ncid1,'mask',varid1))
    call check( nf90_get_var(ncid1,varid1,mask))
    call check( nf90_close(ncid1) )
    
    call check( nf90_open("/data/workspace/mzhang/analysis/FD_impact/Global_growing_season_025deg.nc",NF90_NOWRITE,ncid1))
    call check( nf90_inq_varid(ncid1,'sos',varid1))
    call check( nf90_get_var(ncid1,varid1,sos))
    call check( nf90_inq_varid(ncid1,'eos',varid1))
    call check( nf90_get_var(ncid1,varid1,eos))
    call check( nf90_inq_varid(ncid1,'maxid',varid1))
    call check( nf90_get_var(ncid1,varid1,maxid))
    call check( nf90_close(ncid1) )
ext = "/data/workspace/mzhang/analysis/FD_impact/FD_charac/"

allocate(data2(wlon,wlat,time,nyr))
allocate(et(wlon,wlat,time,nyr))
allocate(et1d(wlon,wlat,time*nyr))
allocate(sm1d(wlon,wlat,time*nyr))

    call check( nf90_open(trim(ext)//'GLDASv2.1_ET_percentile_8day_2001_2022.nc',NF90_NOWRITE,ncid1))
    call check( nf90_inq_varid(ncid1,'ETpercentile',varid1))
    start = (/1,1,1,1/); count = (/wlon,wlat,46,nyr/)
    call check( nf90_get_var(ncid1,varid1,et,start,count))
    call check( nf90_close(ncid1) )
    
    call check( nf90_open(trim(ext)//'SMpercentile_GLDASv2.1_SMroot_8day_2001_2022.nc',NF90_NOWRITE,ncid1))
    call check( nf90_inq_varid(ncid1,'SMpercentile',varid1))
    start = (/1,1,1,1/); count = (/wlon,wlat,46,nyr/)
    call check( nf90_get_var(ncid1,varid1,data2(:,:,:,:),start,count))
    call check( nf90_close(ncid1) )

where(data2.gt.100.or.data2.lt.0) data2 = -999.
smp = data2
cont = 1
do y=1,nyr
  do t=1,time
   sm1d(:,:,cont) = smp(:,:,t,y) 
   et1d(:,:,cont) = et(:,:,t,y) 
   cont = cont + 1
  end do
end do

ext = "/data/workspace/mzhang/analysis/FD_impact/FD_charac/"
drt_sd = 0
drt_fd   = 0

call check( nf90_create(trim(ext)//'FD_SD_GLDASv2.1_SM1m_8day_2001_2022_SM3_20%_ET_50%.nc',or(NF90_64BIT_OFFSET,nf90_clobber),ncid))
call check( nf90_def_dim( ncid,'lon',wlon,lon_dimid ) )
call check( nf90_def_dim( ncid,'lat',wlat,lat_dimid ) )
call check( nf90_def_dim( ncid,'time',time*nyr,cha_dimid ) )

  dimids = (/lon_dimid,lat_dimid,cha_dimid/)
call check( nf90_def_var(ncid, "SMpercentile", NF90_REAL, dimids,varid(1)) )
call check( nf90_put_att(ncid, varid(1), 'time','pentad') )
call check( nf90_put_att(ncid, varid(1), '_FillValue',-999.))
call check( nf90_def_var(ncid, "flashdrt", NF90_INT, dimids,varid(2)) )
call check( nf90_put_att(ncid, varid(2), 'time','pentad') )
call check( nf90_put_att(ncid, varid(2), 'long_name','flash drought state 1 True/0 False'))
call check( nf90_def_var(ncid, "slowdrt", NF90_INT, dimids,varid(3)) )
call check( nf90_put_att(ncid, varid(3), 'units','-') )
call check( nf90_put_att(ncid, varid(3), 'long_name','slow drought state 1 True/0 False'))
call check( nf90_put_att(ncid, nf90_global, 'description',&
      '8day flash drought during 2001-2022') )
call check( nf90_def_var(ncid, "frequency_fd", NF90_REAL, (/lon_dimid,lat_dimid/),varid(4)) )
call check( nf90_put_att(ncid, varid(4), '_FillValue',-999.))
call check( nf90_def_var(ncid, "frequency_sd", NF90_REAL,(/lon_dimid,lat_dimid/),varid(5)) )
call check( nf90_put_att(ncid, varid(5), '_FillValue',-999.))
call check( nf90_enddef(ncid) )

  do i=1,wlon
    do j=1,wlat
      if(mask(i,j).gt.0.1.and.sos(i,j).gt.-0.1.and.sos(i,j).ne.eos(i,j))then
        if(maxid(i,j).le.33.and.maxid(i,j).ge.11.or.maxid(i,j).lt.-0.1)then !north hemisphere & whole year = growing season
          do y=1,nyr
           if(sos(i,j).gt.0)then
             t1 = sos(i,j) + 46*(y-1) !t-1
             t2 = eos(i,j)+1 + 46*(y-1)
           else if(sos(i,j).eq.0)then
             t1 = max(1,sos(i,j) + 46*(y-1))
             t2 = eos(i,j)+1 + 46*(y-1)
           end if
             dd(1:t2-t1+1) = sm1d(i,j,t1:t2)
             tt(1:t2-t1+1) = et1d(i,j,t1:t2)
            call flashdrought(t2-t1+1,dd(1:t2-t1+1),tt(1:t2-t1+1),b_fd(i,j,:,:,y),drt0(1:t2-t1+1))
            drt_fd(i,j,t1+1:t2) = drt0(2:t2-t1+1)
            call slowdrought(t2-t1+1,dd(1:t2-t1+1),tt(1:t2-t1+1),b_sd(i,j,:,:,y),drt0(1:t2-t1+1))
            drt_sd(i,j,t1+1:t2) = drt0(2:t2-t1+1)
          end do!year
        else ! southern hemisphere
          do y=1,nyr+1
            if(y==1)then
              t1 = 1
              t2 = sos(i,j)+1
            else if(y==nyr+1)then
              t1 = eos(i,j) + 46*(y-2)
              t2 = nyr*46 
            else
              t1 = eos(i,j) + 46*(y-2)
              t2 = sos(i,j) +1 + 46*(y-1)
            end if
             dd(1:t2-t1+1) = sm1d(i,j,t1:t2)
             tt(1:t2-t1+1) = et1d(i,j,t1:t2)
            call flashdrought(t2-t1+1,dd(1:t2-t1+1),tt(1:t2-t1+1),b_fd(i,j,:,:,y),drt0(1:t2-t1+1))
            drt_fd(i,j,t1:t2) = drt0(1:t2-t1+1)
            call slowdrought(t2-t1+1,dd(1:t2-t1+1),tt(1:t2-t1+1),b_sd(i,j,:,:,y),drt0(1:t2-t1+1))
            drt_sd(i,j,t1:t2) = drt0(1:t2-t1+1)
          end do  
        end if
          fre_fd(i,j) = 0.
          fre_sd(i,j) = 0.
        do y=1,nyr+1
           fre_fd(i,j) = fre_fd(i,j) + b_fd(i,j,1,1,y)
           fre_sd(i,j) = fre_sd(i,j) + b_sd(i,j,1,1,y)
        end do
      else 
        drt_sd(i,j,:) = 0
        drt_fd(i,j,:) = 0
        fre_fd(i,j)   = -999.
        fre_sd(i,j)   = -999.
      end if !mask 
     end do !lat
   end do !lon

     call check( nf90_put_var(ncid, varid(1),sm1d) )
     call check( nf90_put_var(ncid, varid(2),drt_fd))
     call check( nf90_put_var(ncid, varid(3),drt_sd))
     call check( nf90_put_var(ncid, varid(4),fre_fd))
     call check( nf90_put_var(ncid, varid(5),fre_sd))
     call check( nf90_close(ncid) )
deallocate(et)
deallocate(sm1d)
deallocate(et1d)
deallocate(data2)
print*,"done"
contains
        subroutine check(status)
        integer, intent ( in) :: status
        if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop "Stopped"
        end if
        end subroutine check

end


