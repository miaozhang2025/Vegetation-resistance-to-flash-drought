begin

 f = addfile("input_dataset.nc","r")
 gpp_percent = f->SIFpercentile
 smp = f->SMpercentile
 drt_fd = f->flashdrt
 drt_sd = f->slowdrt

 glon = 40
 glat  = 40 
;******************************************************
;       Data processing
;******************************************************
   drt = drt_fd + drt_sd
   drt = where(ismissing(drt),0,drt)
   drt_fd = where(ismissing(drt_fd),0,drt_fd)
   drt_sd = where(ismissing(drt_sd),0,drt_sd)
   gpp2 = reshape(gpp_percent,(/22*46,glat,glon/))

;******************************************************
;calculation of GPP resistance to flash drought/slow drought 
;******************************************************

resist =  new((/22*15,glat,glon/),"float")
resist!0 = "frequency"
resist!1 = "lat"
resist!2 = "lon"
onsettime = resist
response_time = resist
response_intensity = resist
drt_type  = resist

  do j=0,glat-1
   print(j)
    do i=0,glon-1
      if(.not.ismissing(mask1(j,i)).and.mask1(j,i).gt.0.2.and..not.ismissing(sos(j,i)).and.sos(j,i).ne.eos(j,i))then
         cont    =  0
         drt0    = drt(:,j,i)
         drt_fd0 = drt_fd(:,j,i)
         drt_sd0 = drt_sd(:,j,i)
         smp0    = smp(:,j,i)
         gpp0    = gpp2(:,j,i)         
        do t=1,46*22-1
          fd_flag   = 0 
          onset_bg  = -999
          onset_ed  = -999
          recover_bg= -999
          recover_ed= -999
          drt_type_flag = 0
          if(drt0(t).eq.1.and.drt0(t-1).eq.0)then 
            if(drt_fd0(t).eq.1)then
              drt_type_flag = 1
            else if(drt_sd0(t).eq.1)
              drt_type_flag = 2
            end if
            end if
             fd_flag  = 1
             onset_bg = t
            do tt=onset_bg,46*22-1
               if(drt0(tt).eq.1)then
                 onset_ed = tt
               else
                 break
               end if
            end do
            do tt=onset_ed+1,46*22-1
               if(drt0(onset_ed+1).eq.2)then
                 recover_bg = onset_ed+1 
               end if
               if(drt0(tt).eq.2)then
                 recover_ed = tt
               else 
                 break
               end if
            end do
          end if
          
         if(fd_flag.eq.1)then
            if(recover_ed.gt.0)then
              fd_ed = recover_ed
            else
              fd_ed = onset_ed
            end if
            if(fd_ed-onset_bg.gt.7)then
             fd_ed = onset_bg + 7 
            end if
               fd_ed = min((/fd_ed,46*22-1/))
            gpp_pre  = gpp0(onset_bg-1)
            if(.not.ismissing(gpp_pre))then
               if(gpp_pre.gt.40)then
                  onsettime(cont,j,i) = onset_bg
                  gpp_min = min(gpp0(onset_bg-1:fd_ed))
                  min_id  = dim_minind(gpp0(onset_bg:fd_ed),0)
                  response_time(cont,j,i)         = min_id+1
                  response_intensity(cont,j,i)   = (gpp_pre-gpp_min)
                 if(min_id.gt.1)then
                   resist(cont,j,i)    = regCoef(ispan(onset_bg-1,onset_bg-1+min_id,1),gpp0(onset_bg-1:onset_bg-1+min_id))
                 else if(min_id.eq.1)
                   resist(cont,j,i)    = (gpp_min-gpp0(onset_bg-1))/2.
                 end if
                 end if 
                  drt_type(cont,j,i)  = drt_type_flag     
                  cont = cont + 1
               end if 
            end if
          end if

        end do
      end if
  end do 
end do
resist2 = new((/22,15,glat,glon/),"float")
resist2!0 = "year"
resist2!1 = "frequency"
resist2!2 = "lat"
resist2!3 = "lon"
drt_type2 = resist2 
response_time2 = resist2
 response_intensity2 = resist2

do j=0,glat-1
do i=0,glon-1
do y=0,22-1
  t1 = y*46
  t2 = y*46+45
  t0 = onsettime(:,j,i)
  id0 = ind(t0.ge.t1.and.t0.le.t2)
  num0= num(t0.ge.t1.and.t0.le.t2)
  if(.not.ismissing(num0).and.num0.gt.1)then
    resist2(y,0:num0-1,j,i)       = resist(id0,j,i) 
    response_time2(y,0:num0-1,j,i)       = response_time(id0,j,i) 
    response_intensity2(y,0:num0-1,j,i) = response_intensity(id0,j,i) 
    drt_type2(y,0:num0-1,j,i)     = drt_type(id0,j,i)
  else if (.not.ismissing(num0).and.num0.eq.1)then
    resist2(y,0,j,i)       = resist(id0,j,i)
    response_time2(y,0,j,i)       = response_time(id0,j,i) 
    response_intensity2(y,0,j,i) = response_intensity(id0,j,i) 
    drt_type2(y,0,j,i)     = drt_type(id0,j,i)
  end if
  end if
  delete(id0)
  delete(num0)
end do
end do
end do

;******************************************************
;   Output file
;******************************************************
outfile = "Resistance2Drt_2001_2022.nc" 
system("rm -f "+ outfile)
fn = addfile(outfile,"c")
fn->resist   = resist2
fn->responsetime = response_time2
fn->response_intensity = response_time2
fn->drt_type = (/drt_type2/)
print("done")
end 
