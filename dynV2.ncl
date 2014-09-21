 load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 begin 
pss=addfile("../../../chenjh/DATA/pres.ps.sea.nc","r")
      ps= (pss->pres(2,:,:)*10 + 359650.0)/100.0 ;;;; 
;      und!0 ="time"
      ps!0 ="lat"
      ps!1 ="lon"
      ps&lat = pss->lat
      ps&lon = pss->lon
      ps@units="hPa"
;
 ;;; first draw the mean condition of wind, temperature, vapor
  umf=addfile("../../../chenjh/DATA/Muliyearmean/uwnd.ymon.nc","r")
      und= (umf->uwnd(:,:,:,:)*0.01 + 187.65) ;;;; uwnd(t,z,y,x)
      und!0 ="time"
      und!1 ="level"
      und!2 ="lat"
      und!3 ="lon"
      und&time = umf->time
      und&level = umf->level
      und&lat = umf->lat
      und&lon = umf->lon
      und@units="m s~S5~-1"
      level=umf->level 
;      
    vmf=addfile("../../../chenjh/DATA/Muliyearmean/vwnd.ymon.nc","r")
      vnd= (vmf->vwnd(:,:,:,:)*0.01 + 187.65)
      vnd!0 ="time"
      vnd!1 ="level"
      vnd!2 ="lat"
      vnd!3 ="lon"
      vnd&time = vmf->time
      vnd&level = vmf->level
      vnd&lat = vmf->lat
      vnd&lon = vmf->lon
      vnd@units="m s~S5~-1" 
;
  tmf=addfile("../../../chenjh/DATA/Muliyearmean/air.ymon.nc","r")
    tm=(tmf->air(:,:,:,:)*0.01 + 465.15) 
      tm!0 ="time"
      tm!1 ="level"     
      tm!2 ="lat"
      tm!3 ="lon"
      tm&time = tmf->time
      tm&level = tmf->level
      tm&lat = tmf->lat
      tm&lon = tmf->lon
;      tm@units="K"
;
  hmf=addfile("../../../chenjh/DATA/Muliyearmean/rhum.ymon.nc","r")
    hm=(hmf->rhum(:,:,:,:)*0.01 + 302.65) 
      hm!0 ="time"
      hm!1 ="level"     
      hm!2 ="lat"
      hm!3 ="lon"
      hm&time = hmf->time
      hm&level = hmf->level
      hm&lat = hmf->lat
      hm&lon = hmf->lon
;      hm@units="%"
dimth = dimsizes(umf->uwnd)
     nt=dimth(0)    ; time
     nl=dimth(1)    ; level
     nyy=dimth(2)    ; lat
     nxx=dimth(3)    ; lon  
      divs   = new ( (/nl,nyy,nxx /), float)  ; divergence
      vorts  = new ( (/nl,nyy,nxx /), float)  ; vorticity  (relative)
;
      div   = new ( (/nt,nl,nyy,nxx /), float)  ; divergence
      vort  = new ( (/nt,nl,nyy,nxx /), float)  ; vorticity  (relative)
      div!0 ="time"
      div!1 ="level"     
      div!2 ="lat"
      div!3 ="lon"
      div&time = umf->time
      div&level = umf->level
      div&lat = umf->lat
      div&lon = umf->lon
      div@units=" s~S5~-1"
      vort!0 ="time"
      vort!1 ="level"     
      vort!2 ="lat"
      vort!3 ="lon"
      vort&time = umf->time
      vort&level = umf->level
      vort&lat = umf->lat
      vort&lon = umf->lon
      vort@units=" s~S5~-1"
      uv2vrdvf (und,vnd,vort,div)
;      print(vort)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      divs!0 ="time"
      divs!0 ="level"     
      divs!1 ="lat"
      divs!2 ="lon"
;      divs&time = umf->time
      divs&level = umf->level
      divs&lat = umf->lat
      divs&lon = umf->lon
      divs@units=" s~S5~-1"
;      vorts!0 ="time"
      vorts!0 ="level"     
      vorts!1 ="lat"
      vorts!2 ="lon"
;      vorts&time = umf->time
      vorts&level = umf->level
      vorts&lat = umf->lat
      vorts&lon = umf->lon
      vorts@units=" s~S5~-1"
; 
       unds   = new ( (/nl,nyy,nxx /), float)  ;
       vnds  = new ( (/nl,nyy,nxx /), float)  ; 
       tms   = new ( (/nl,nyy,nxx /), float)  ;
       hms  = new ( (/nl,nyy,nxx /), float)  ;
      unds!0 ="level"     
      unds!1 ="lat"
      unds!2 ="lon"
      unds&level = umf->level
      unds&lat = umf->lat
      unds&lon = umf->lon
      unds@units="m s~S5~-1"
;     
      vnds!0 ="level"     
      vnds!1 ="lat"
      vnds!2 ="lon"
      vnds&level = umf->level
      vnds&lat = umf->lat
      vnds&lon = umf->lon
      vnds@units="m s~S5~-1"
; 
      tms!0 ="level"     
      tms!1 ="lat"
      tms!2 ="lon"
      tms&level = umf->level
      tms&lat = umf->lat
      tms&lon = umf->lon
;      tms@units="K"
;     
      hms!0 ="level"     
      hms!1 ="lat"
      hms!2 ="lon"
      hms&level = umf->level
      hms&lat = umf->lat
      hms&lon = umf->lon
;      hms@units="%"
;
; the average period
its=4 ;;may
ite=8  ;; Sep.
       vorts=0.0
       divs=0.0
       unds=0.0
       vnds=0.0
       tms=0.0
       hms=0.0
      do itt=its,ite,1
       vorts(:,:,:)=vorts(:,:,:)+vort(itt,:,:,:)
       divs(:,:,:)=divs(:,:,:)+div(itt,:,:,:)
       unds(:,:,:)=unds(:,:,:)+und(itt,:,:,:)
       vnds(:,:,:)=vnds(:,:,:)+vnd(itt,:,:,:)
       tms(:,:,:)=tms(:,:,:)+tm(itt,:,:,:)
       hms(:,:,:)=hms(:,:,:)+hm(itt,:,:,:)      
      end do
      vorts(:,:,:)=vorts(:,:,:)/(ite-its+1)
       divs(:,:,:)=divs(:,:,:)/(ite-its+1)
       unds(:,:,:)=unds(:,:,:)/(ite-its+1)
       vnds(:,:,:)=vnds(:,:,:)/(ite-its+1)
       tms(:,:,:)=tms(:,:,:)/(ite-its+1)
       hms(:,:,:)=hms(:,:,:)/(ite-its+1)
;       print(unds)
; We generate plots, but what kind do we prefer?
;  type = "x11"
   type = "pdf"
;  type = "ps"
; type = "ncgm"
;----------------------------------------------------------
;- TITLE OF THE WINDOW
do lv=0,11,1
    lev=sprintf("%3.0f",level(lv))

;    if (6 .gt. lv ) then
;      vorts=mask(vorts,(ps .gt. level(lv)), True)
;      divs=mask(divs,(ps .gt. level(lv)), True)
;      unds=mask(unds,(ps .gt. level(lv)), True)
;      vnds=mask(vnds,(ps .gt. level(lv)), True)
;      hms=mask(hms,(ps .gt. level(lv)), True)
;      tms=mask(tms,(ps .gt. level(lv)), True)
;    end if
     divs@description = "Divergence "+lev+"hPa"
     vorts@description = "Vorticity "+lev+"hPa"
;     tms@description = "Temperature "+lev+"hPa"
;     hms@description = "Relative Humidity "+lev+"hPa"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;     
  titleOW="../Pics/Div"+lev
  wks = gsn_open_wks(type,titleOW)
  res =True
  FirstTime = True
  ntimes=1
  ;;;;;; set for draw ----------------------
;  wks = gsn_open_wks ("ps", "mask")
  gsn_define_colormap(wks,"GreenMagenta16")

; Set up resource list for plot options.
  opts =res   
  opts@gsnPaperOrientation = "auto"
  opts@gsnMaximize         = True   ;                     
  opts@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  opts@tmYLLabelFontHeightF = 0.02
  opts@lbLabelFontHeightF = 0.02
  opts@gsnStringFont = 25
  opts@cnLinesOn  = False
; opts@cnLineColor  = "Black" ;colors(0)    ;blue color for gpm line
  opts@gsnFrame = False
  opts@cnFillOn             = True
  opts@gsnSpreadColors  = False
  opts@lbLabelAutoStride  = True
  opts@cnLevelSelectionMode ="AutomaticLevels" ; "ExplicitLevels"  ;"
;  print(n)
;  opts@cnLevels                = (/-1e-6,-8e-7,-6e-7,-4e-7,-3e-7,-2e-7,-1e-7, 0 ,1e-7,2e-7,3e-7,4e-7,6e-7,8e-7,1e-6/)
  opts@cnLinesOn                   = False  ; Turn off contour lines
;  opts@cnFillColors                = (/3,4,5,6,7,12,13,14,15,16/)
  opts@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
                                               ; last.
  opts@gsnDraw = False                        ; do not draw the plot
  opts@gsnFrame = False  
;;;;  map
      LAT0=umf->lat
      LON0=umf->lon
      ds1= dimsizes(umf->lat)
      ds2= dimsizes(umf->lon)
      opts@mpLimitMode ="Corners"
      opts@mpLeftCornerLatF = LAT0({5})
      opts@mpLeftCornerLonF = LON0({60})
      opts@mpRightCornerLatF = LAT0({60})
      opts@mpRightCornerLonF = LON0({150})
;for map
      opts@mpFillOn =False
      opts@mpGeophysicalLineColor      = "Black"
      opts@mpGridLineColor             = "Black"
      opts@mpLimbLineColor             = "Black"
      opts@mpNationalLineColor         = "Black"
      opts@mpPerimLineColor            = "Black"
      opts@mpUSStateLineColor          = "Black"
      opts@mpOutlineSpecifiers =(/"China"/)  
      plot = gsn_csm_contour_map(wks,divs(lv,:,:), opts) 
      draw(plot)
      frame (wks)
;- TITLE OF THE WINDOW
  titleOW="../Pics/Vort"+lev
  wks = gsn_open_wks(type,titleOW)
  res =True
  FirstTime = True
  ntimes=1
  ;;;;;; set for draw ----------------------
;  wks = gsn_open_wks ("ps", "mask")
  gsn_define_colormap(wks, "GreenMagenta16")

; Set up resource list for plot options.
  opts =res   
  opts@gsnPaperOrientation = "auto"
  opts@gsnMaximize         = True   ;                     
  opts@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  opts@tmYLLabelFontHeightF = 0.02
  opts@lbLabelFontHeightF = 0.02
  opts@gsnStringFont = 25
  opts@cnLinesOn  = False
; opts@cnLineColor  = "Black" ;colors(0)    ;blue color for gpm line
  opts@gsnFrame = False
  opts@cnFillOn             = True
  opts@gsnSpreadColors  = False
  opts@lbLabelAutoStride  = True
  opts@cnLevelSelectionMode = "AutomaticLevels" ;"ExplicitLevels"  ;
;  print(n)
;   opts@cnLevels                  = (/-1e-6,-8e-7,-6e-7,-4e-7,-3e-7,-2e-7,-1e-7, 0 ,1e-7,2e-7,3e-7,4e-7,6e-7,8e-7,1e-6/)
   opts@cnLinesOn                   = False  ; Turn off contour lines
;   opts@cnFillColors                = (/4,5,6,8,10,11,12,14,16,17/)
   opts@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
                                               ; last.
  opts@gsnDraw = False                        ; do not draw the plot
  opts@gsnFrame = False  
;;;;  map
      LAT0=umf->lat
      LON0=umf->lon
      ds1= dimsizes(umf->lat)
      ds2= dimsizes(umf->lon)
      opts@mpLimitMode ="Corners"
      opts@mpLeftCornerLatF = LAT0({5})
      opts@mpLeftCornerLonF = LON0({60})
      opts@mpRightCornerLatF = LAT0({60})
      opts@mpRightCornerLonF = LON0({150})
;for map
      opts@mpFillOn =False
      opts@mpGeophysicalLineColor      = "Black"
      opts@mpGridLineColor             = "Black"
      opts@mpLimbLineColor             = "Black"
      opts@mpNationalLineColor         = "Black"
      opts@mpPerimLineColor            = "Black"
      opts@mpUSStateLineColor          = "Black"
      opts@mpOutlineSpecifiers =(/"China"/)  
      plot = gsn_csm_contour_map(wks,vorts(lv,:,:), opts) 
      draw(plot)
      frame (wks)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;- TITLE OF THE WINDOW
  titleOW="../Pics/T+humdity"+lev
  wks = gsn_open_wks(type,titleOW)
  res =True
  FirstTime = True
  ntimes=1
  ;;;;;; set for draw ----------------------
;  wks = gsn_open_wks ("ps", "overlay")
;  gsn_define_colormap(wks, "GMT_haxby")

;  
  vcres           = True
  vcres@gsnDraw   = False
  vcres@gsnFrame  = False
  stres           = vcres
  cnres           = vcres
  mpres           = vcres
; Set up some vector resources.
;
  vcres@vcLevelSelectionMode     = "ManualLevels"
  vcres@vcMinLevelValF           = 10.0
  vcres@vcMaxLevelValF           =  100.0
  vcres@vcLevelSpacingF          = 10.0    ; 20.
;  vcres@gsnSpreadColors = True
;  vcres@gsnSpreadColorEnd = 23
   vcres@vcLevelPalette = "percent_11lev" ;"CBR_wet"
;
; Vector lengths and color
;
  vcres@vcFillArrowsOn           = True
  vcres@vcLineArrowThicknessF    =  2.5

  vcres@vcMinFracLengthF         = 0.33
  vcres@vcMinMagnitudeF          = 0.001
  vcres@vcMonoFillArrowFillColor = False
  vcres@vcMonoLineArrowColor     = False

  vcres@vcRefLengthF             = 0.050
  vcres@vcRefMagnitudeF          = 20.0
  vcres@vcRefAnnoOrthogonalPosF  = -0.12
  vcres@vcRefAnnoParallelPosF    =  0.997
  vcres@vcRefAnnoFontHeightF     = 0.015

  vcres@lbTitleString            = "Relative Humidity"
  vcres@lbTitleOffsetF           = -0.25
  vcres@lbTitleFontHeightF       = 0.015
  vcres@lbLabelFontHeightF       = 0.015
  vcres@lbLabelAutoStride        = True
  vcres@pmLabelBarOrthogonalPosF = .06           ; move whole thing down
  vcres@lbTitlePosition     = "Bottom"
;
; Make sure vectors are drawn in "predraw" phase.
;
  vcres@vcVectorDrawOrder        = "Predraw"
;
; Set up some contour resources.
;
; Turn on contour fill, and turn other things off.
;
  cnres@cnFillOn             = True
  cnres@cnLinesOn            = False
  cnres@cnInfoLabelOn        = False
  cnres@lbOrientation        = "Vertical"
  cnres@lbTitleString        = "Temperature"
  cnres@lbTitlePosition      = "Left"
  cnres@lbBoxLinesOn         = False
  cnres@lbTitleFontHeightF   = 0.020
  cnres@lbLabelFontHeightF   = 0.020
  cnres@lbLabelAutoStride = True 
;  cnres@lbAutoManage = True
   cnres@pmLabelBarOrthogonalPosF =-0.01   
  cnres@lbBoxLinesOn          = False        ; turn off box between colors

;
; Define contour levels and their colors.
;
  cnres@cnLevelSelectionMode = "ManualLevels"
  cnres@cnLevelSpacingF      = 5.0
  if(lv .eq. 0) then
  cnres@cnMinLevelValF       = 285
  cnres@cnMaxLevelValF       = 315
  end if
  if(lv.gt.0 .and. lv .lt.3)then
  cnres@cnMinLevelValF       = 280
  cnres@cnMaxLevelValF       = 310
  end if
  if(lv.gt.2 .and. lv .lt.5)then
  cnres@cnMinLevelValF       = 245
  cnres@cnMaxLevelValF       = 275
  cnres@cnLevelSpacingF      = 3.0
  end if
  if(lv.gt.4.and. lv .lt.8)then
  cnres@cnMinLevelValF       = 210
  cnres@cnMaxLevelValF       = 240
  cnres@cnLevelSpacingF      = 3.0
  end if
  if(lv.gt.7)then
  cnres@cnMinLevelValF       = 200
  cnres@cnMaxLevelValF       = 230
  cnres@cnLevelSpacingF      = 3.0
  end if
;  
  cnres@cnMinLevelValF       = 210
  cnres@cnMaxLevelValF       = 315
  cnres@cnLevelSpacingF      = 1.0
  
; use wks color map
;  cnres@gsnSpreadColors = True
;  cnres@gsnSpreadColorStart = 4
;  cnres@gsnSpreadColorEnd = 26
;========================================
; use defined color map
;  cnres@cnFillPalette = "sunshine_9lev"
   cnres@cnFillPalette = "BlAqGrYeOrRe"
;
; Make sure contours are drawn in "predraw" phase.
;
  cnres@cnFillDrawOrder      = "Predraw"

;
; Control appearance of map.
;
;;;;  map
      LAT0=umf->lat
      LON0=umf->lon
      ds1= dimsizes(umf->lat)
      ds2= dimsizes(umf->lon)
      mpres@mpLimitMode ="Corners"
      mpres@mpLeftCornerLatF = LAT0({5})
      mpres@mpLeftCornerLonF = LON0({60})
      mpres@mpRightCornerLatF = LAT0({60})
      mpres@mpRightCornerLonF = LON0({150})
      ;for map
      mpres@mpFillOn =False
      mpres@mpGeophysicalLineColor      = "Black"
      mpres@mpGridLineColor             = "Black"
      mpres@mpLimbLineColor             = "Black"
      mpres@mpNationalLineColor         = "Black"
      mpres@mpPerimLineColor            = "Black"
      mpres@mpUSStateLineColor          = "Black"
      mpres@mpOutlineSpecifiers =(/"China"/)

  mpres@mpGridAndLimbDrawOrder = "Predraw"
  mpres@tiMainString           = "Wind Field and Humidity"  ;title
;
; Create, but don't draw, a vector, streamline, contour
; and map plot. Stride the vector and contour data to
; thin it a little.
;
  vcid = gsn_csm_vector_scalar(wks,unds(lv,:,:),vnds(lv,:,:), \
                                   hms(lv,:,:),vcres)
;  stid = gsn_csm_streamline(wks,u500(0,:,:),v500(0,:,:),stres)
  cnid = gsn_csm_contour(wks,tms(lv,:,:),cnres)
  mpid = gsn_csm_map(wks,mpres)

;
; Overlay contour, streamline, and vector plots on the map plot.
;
  overlay(mpid,cnid)
;  overlay(mpid,stid)
  overlay(mpid,vcid)
;  draw(mpid)
;  frame(wks)

  maximize_output(wks,True)

   end do  ;;;; lv loop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; draw the events   
edate1=asciiread("ETP_EventsDate_cloudsat.txt",-1,"integer")
ednn=dimsizes(edate1)
edate=new((/ednn/4,4/),integer)
ii=0
;print(ednn)

do i=0,ednn-1,4
;print(i)
edate(ii,0)=edate1(i)
;print(edate(ii,0))
edate(ii,1)=edate1(i+1)
edate(ii,2)=edate1(i+2)
edate(ii,3)=edate1(i+3)
ii=ii+1
end do
;;;;;/1405,4/,  /:,0/ year, /:,1/ month , /:,2/ day /:,3/ rainfall (float actually)

su=new((/17,73,144/),float)
sv=new((/17,73,144/),float)
st=new((/17,73,144/),float)
sh=new((/17,73,144/),float)
sdiv=new((/17,73,144/),float)
svor=new((/17,73,144/),float)
;;;
      su!0 ="level"
      su!1 ="lat"
      su!2 ="lon"
      su&level = umf->level
      su&lat = umf->lat
      su&lon = umf->lon
      su@units="m s~S5~-1"
;;;
      sv!0 ="level"
      sv!1 ="lat"
      sv!2 ="lon"
      sv&level = vmf->level
      sv&lat = vmf->lat
      sv&lon = vmf->lon
      sv@units="m s~S5~-1"
;;;
      st!0 ="level"
      st!1 ="lat"
      st!2 ="lon"
      st&level = tmf->level
      st&lat = tmf->lat
      st&lon = tmf->lon
;      st@units="K"
;;;
      sh!0 ="level"
      sh!1 ="lat"
      sh!2 ="lon"
      sh&level = hmf->level
      sh&lat = hmf->lat
      sh&lon = hmf->lon
;      sh@units="%"
;;;
      sdiv!0 ="level"
      sdiv!1 ="lat"
      sdiv!2 ="lon"
      sdiv&level = umf->level
      sdiv&lat = umf->lat
      sdiv&lon = umf->lon
      sdiv@units=" s~S5~-1"
;;;
      svor!0 ="level"
      svor!1 ="lat"
      svor!2 ="lon"
      svor&level = umf->level
      svor&lat =umf->lat
      svor&lon = umf->lon
      svor@units=" s~S5~-1"
 su=0.0      
 sv=0.0
 st=0.0
 sh=0.0
 sdiv=0.0
 svor=0.0
 isum=0
do i = 0, ednn/4-1
  iopen=0
  print(i)
   if (i.eq.0) then
   year=sprintf("%4.0f",edate(i,0))
;  print(year)
;   print(edate(i,0))
   iopen=1
   end if
   if (i.gt.0 .and. (edate(i,0).ne.edate(i-1,0))) then
   year=sprintf("%4.0f",edate(i,0))
   iopen=1
   end if
   if(iopen.eq.1) then
;;;;; close the old files
    if(i.gt.0) then
    delete(tempu)
    delete(tempv)
    delete(tempt)
    delete(temph)
    delete(utm)
    delete(vtm)
    delete(ttm)
    delete(htm)
    delete(divtm)
    delete(vortm)
    end if
;;; open the new files   
   tempu=addfile("../../../chenjh/DATA/Dailymean/uwnd."+year+".nc","r")
   tempv=addfile("../../../chenjh/DATA/Dailymean/vwnd."+year+".nc","r")
   tempt=addfile("../../../chenjh/DATA/Dailymean/air."+year+".nc","r")
   temph=addfile("../../../chenjh/DATA/Dailymean/rhum."+year+".nc","r")
   end if
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; digest data from files
      utm= (tempu->uwnd(:,:,:,:)*0.01 + 187.65) ;;;; uwnd(t,z,y,x)
      utm!0 ="time"
      utm!1 ="level"
      utm!2 ="lat"
      utm!3 ="lon"
      utm&time = tempu->time
      utm&level = tempu->level
      utm&lat = tempu->lat
      utm&lon = tempu->lon
      utm@units="m s~S5~-1"
;      
      vtm= (tempv->vwnd(:,:,:,:)*0.01 + 187.65)
      vtm!0 ="time"
      vtm!1 ="level"
      vtm!2 ="lat"
      vtm!3 ="lon"
      vtm&time = tempv->time
      vtm&level = tempv->level
      vtm&lat = tempv->lat
      vtm&lon = tempv->lon
      vtm@units="m s~S5~-1"     
;
      htm= (temph->rhum(:,:,:,:)*0.01 + 302.65)
      htm!0 ="time"
      htm!1 ="level"
      htm!2 ="lat"
      htm!3 ="lon"
      htm&time = temph->time
      htm&level = temph->level
      htm&lat = temph->lat
      htm&lon = temph->lon
;      htm@units="%" 
;
      ttm= (tempt->air(:,:,:,:)*0.01 + 465.15)
      ttm!0 ="time"
      ttm!1 ="level"
      ttm!2 ="lat"
      ttm!3 ="lon"
      ttm&time = tempt->time
      ttm&level = tempt->level
      ttm&lat = tempt->lat
      ttm&lon = tempt->lon
;      ttm@units="K" 
      nnn=dimsizes(ttm) 
;
      divtm   = new ( (/nnn(0),nnn(1),nnn(2),nnn(3)/), float)  ; divergence
      vortm  = new ( (/nnn(0),nnn(1),nnn(2),nnn(3)/), float)  ; vorticity  (relative)
      uv2vrdvf (utm,vtm,vortm,divtm)
      divtm!0 ="time"
      divtm!1 ="level"
      divtm!2 ="lat"
      divtm!3 ="lon"
      divtm&time = tempu->time
      divtm&level = tempu->level
      divtm&lat = tempu->lat
      divtm&lon = tempu->lon
      divtm@units=" s~S5~-1" 
;
      vortm!0 ="time"
      vortm!1 ="level"
      vortm!2 ="lat"
      vortm!3 ="lon"
      vortm&time = tempu->time
      vortm&level = tempu->level
      vortm&lat = tempu->lat
      vortm&lon = tempu->lon
      vortm@units=" s~S5~-1"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
     dys=new(12,integer)
     ny=edate(i,0)
     do im=0,11,1
     dys(im)=31
     end do
     dys(2)=28
     dys(3)=30
     dys(5)=30
     dys(8)=30
     dys(10)=30
     if ( ( mod(ny,4) .eq. 0) .and. ( mod(ny,400) .ne. 0) ) then
     dys(1)=29
     end if
     if ( mod(ny,400) .eq. 0 ) then
     dys(1)=29
     end if
     idex=0
     do j= 0, edate(i,1)-2  ;;;;;the current is not a complete month and counts in the day
        idex=dys(j)+idex
     end do 
     idex=idex+edate(i,2)-1
;     do id=0,3
       su(:,:,:)=su(:,:,:)+utm(idex,:,:,:)
       sv(:,:,:)=sv(:,:,:)+vtm(idex,:,:,:)
       st(:,:,:)=st(:,:,:)+ttm(idex,:,:,:)
       sh(:,:,:)=sh(:,:,:)+htm(idex,:,:,:)
       sdiv(:,:,:)=sdiv(:,:,:)+divtm(idex,:,:,:)
       svor(:,:,:)=svor(:,:,:)+vortm(idex,:,:,:)
       isum=isum+1
;     end do


  end do  ; i loop 
       su(:,:,:)=su(:,:,:)/isum
       sv(:,:,:)=sv(:,:,:)/isum
       st(:,:,:)=st(:,:,:)/isum
       sh(:,:,:)=sh(:,:,:)/isum
       sdiv(:,:,:)=sdiv(:,:,:)/isum
       svor(:,:,:)=svor(:,:,:)/isum
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot setting    
; We generate plots, but what kind do we prefer?
;  type = "x11"
   type = "pdf"
;  type = "ps"
; type = "ncgm"
;----------------------------------------------------------
;- TITLE OF THE WINDOW
do lv=0,11,1
    lev=sprintf("%3.0f",level(lv))

;    if (6 .gt. lv ) then
;      svor=mask(svor,(ps .gt. level(lv)), True)
;      sdiv=mask(sdiv,(ps .gt. level(lv)), True)
;      su=mask(su,(ps .gt. level(lv)), True)
;      sv=mask(sv,(ps .gt. level(lv)), True)
;      sh=mask(sh,(ps .gt. level(lv)), True)
;      st=mask(st,(ps .gt. level(lv)), True)
;    end if
     divs@description = "Divergence "+lev+"hPa"
     vorts@description = "Vorticity "+lev+"hPa"
;     tms@description = "Temperature "+lev+"hPa"
;     hms@description = "Relative Humidity "+lev+"hPa"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;     
  titleOW="../Pics/Div_events"+lev
  wks = gsn_open_wks(type,titleOW)
  res =True
  FirstTime = True
  ntimes=1
  ;;;;;; set for draw ----------------------
;  wks = gsn_open_wks ("pdf", "mask")
  gsn_define_colormap(wks,"GreenMagenta16")

; Set up resource list for plot options.
  opts =res   
  opts@gsnPaperOrientation = "auto"
  opts@gsnMaximize         = True   ;                     
  opts@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  opts@tmYLLabelFontHeightF = 0.02
  opts@lbLabelFontHeightF = 0.02
  opts@gsnStringFont = 25
  opts@cnLinesOn  = False
; opts@cnLineColor  = "Black" ;colors(0)    ;blue color for gpm line
  opts@gsnFrame = False
  opts@cnFillOn             = True
  opts@gsnSpreadColors  = False
  opts@lbLabelAutoStride  = True
  opts@cnLevelSelectionMode = "AutomaticLevels" ;"ExplicitLevels"  ;
;  print(n)
;  opts@cnLevels                = (/-1e-6,-8e-7,-6e-7,-4e-7,-3e-7,-2e-7,-1e-7, 0 ,1e-7,2e-7,3e-7,4e-7,6e-7,8e-7,1e-6/)
  opts@cnLinesOn                   = False  ; Turn off contour lines
;  opts@cnFillColors                = (/3,4,5,6,7,12,13,14,15,16/)
  opts@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
                                               ; last.
  opts@gsnDraw = False                        ; do not draw the plot
  opts@gsnFrame = False  
;;;;  map
      LAT0=umf->lat
      LON0=umf->lon
      ds1= dimsizes(umf->lat)
      ds2= dimsizes(umf->lon)
      opts@mpLimitMode ="Corners"
      opts@mpLeftCornerLatF = LAT0({5})
      opts@mpLeftCornerLonF = LON0({60})
      opts@mpRightCornerLatF = LAT0({60})
      opts@mpRightCornerLonF = LON0({150})
;for map
      opts@mpFillOn =False
      opts@mpGeophysicalLineColor      = "Black"
      opts@mpGridLineColor             = "Black"
      opts@mpLimbLineColor             = "Black"
      opts@mpNationalLineColor         = "Black"
      opts@mpPerimLineColor            = "Black"
      opts@mpUSStateLineColor          = "Black"
      opts@mpOutlineSpecifiers =(/"China"/)  
      plot = gsn_csm_contour_map(wks,sdiv(lv,:,:), opts) 
      draw(plot)
      frame (wks)
;- TITLE OF THE WINDOW
  titleOW="../Pic/Vort_events"+lev
  wks = gsn_open_wks(type,titleOW)
  res =True
  FirstTime = True
  ntimes=1
  ;;;;;; set for draw ----------------------
;  wks = gsn_open_wks ("pdf", "mask")
;  gsn_define_colormap(wks, "GreenMagenta16")

; Set up resource list for plot options.
  opts =res   
  opts@gsnPaperOrientation = "auto"
  opts@gsnMaximize         = True   ;                     
  opts@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  opts@tmYLLabelFontHeightF = 0.02
  opts@lbLabelFontHeightF = 0.02
  opts@gsnStringFont = 25
  opts@cnLinesOn  = False
; opts@cnLineColor  = "Black" ;colors(0)    ;blue color for gpm line
  opts@gsnFrame = False
  opts@cnFillOn             = True
  opts@gsnSpreadColors  = False
  opts@lbLabelAutoStride  = True
  opts@cnLevelSelectionMode = "AutomaticLevels"  ;"ExplicitLevels"  ;
;  print(n)
;   opts@cnLevels                  = (/-1e-6,-8e-7,-6e-7,-4e-7,-3e-7,-2e-7,-1e-7, 0 ,1e-7,2e-7,3e-7,4e-7,6e-7,8e-7,1e-6/)
   opts@cnLinesOn                   = False  ; Turn off contour lines
;   opts@cnFillColors                = (/4,5,6,8,10,11,12,14,16,17/)
   opts@cnFillDrawOrder             = "PreDraw"  ; Make sure map fill happens
                                               ; last.
  opts@gsnDraw = False                        ; do not draw the plot
  opts@gsnFrame = False  
;;;;  map
      LAT0=umf->lat
      LON0=umf->lon
      ds1= dimsizes(umf->lat)
      ds2= dimsizes(umf->lon)
      opts@mpLimitMode ="Corners"
      opts@mpLeftCornerLatF = LAT0({5})
      opts@mpLeftCornerLonF = LON0({60})
      opts@mpRightCornerLatF = LAT0({60})
      opts@mpRightCornerLonF = LON0({150})
;for map
      opts@mpFillOn =False
      opts@mpGeophysicalLineColor      = "Black"
      opts@mpGridLineColor             = "Black"
      opts@mpLimbLineColor             = "Black"
      opts@mpNationalLineColor         = "Black"
      opts@mpPerimLineColor            = "Black"
      opts@mpUSStateLineColor          = "Black"
      opts@mpOutlineSpecifiers =(/"China"/)  
      plot = gsn_csm_contour_map(wks,svor(lv,:,:), opts) 
      draw(plot)
      frame (wks)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;- TITLE OF THE WINDOW
  titleOW="../Pics/T+humdity_events"+lev
  wks = gsn_open_wks(type,titleOW)
  res =True
  FirstTime = True
  ntimes=1
  ;;;;;; set for draw ----------------------
;  wks = gsn_open_wks ("pdf", "overlay")
;  gsn_define_colormap(wks, "BlueDarkRed18")

;  
  vcres           = True
  vcres@gsnDraw   = False
  vcres@gsnFrame  = False
  stres           = vcres
  cnres           = vcres
  mpres           = vcres
; Set up some vector resources.
;
  vcres@vcLevelSelectionMode     = "ManualLevels"
  vcres@vcMinLevelValF           = 10.0
  vcres@vcMaxLevelValF           =  100.0
  vcres@vcLevelSpacingF          = 10.0    ; 20.
;  vcres@gsnSpreadColors = True
; vcres@gsnSpreadColorEnd = 17
;
; Vector lengths and color
;
  vcres@vcLevelPalette = "percent_11lev"
  vcres@vcFillArrowsOn           = True
  vcres@vcLineArrowThicknessF    =  2.0

  vcres@vcMinFracLengthF         = 0.33
  vcres@vcMinMagnitudeF          = 0.001
  vcres@vcMonoFillArrowFillColor = False
  vcres@vcMonoLineArrowColor     = False

  vcres@vcRefLengthF             = 0.045
  vcres@vcRefMagnitudeF          = 20.0
  vcres@vcRefAnnoOrthogonalPosF  = -0.12
  vcres@vcRefAnnoParallelPosF    =  0.997
  vcres@vcRefAnnoFontHeightF     = 0.015

  vcres@lbTitleString            = "Relative Humidity"
  vcres@lbTitleOffsetF           = -0.25
  vcres@lbTitleFontHeightF       = 0.015
  vcres@lbLabelFontHeightF       = 0.015
  vcres@lbLabelAutoStride        = True
  vcres@pmLabelBarOrthogonalPosF = .06           ; move whole thing down
  vcres@lbTitlePosition     = "Bottom"
;
; Make sure vectors are drawn in "predraw" phase.
;
  vcres@vcVectorDrawOrder        = "Predraw"
;
; Set up some contour resources.
;
; Turn on contour fill, and turn other things off.
;
  cnres@cnFillOn             = True
  cnres@cnLinesOn            = False
  cnres@cnInfoLabelOn        = False
  cnres@lbOrientation        = "Vertical"
  cnres@lbTitleString        = "Temperature"
  cnres@lbTitlePosition      = "Left"
  cnres@lbBoxLinesOn         = False
  cnres@lbTitleFontHeightF   = 0.020
  cnres@lbLabelFontHeightF   = 0.020
  cnres@lbLabelAutoStride = True 
;  cnres@lbAutoManage = True
cnres@pmLabelBarOrthogonalPosF =-0.01
  cnres@lbBoxLinesOn          = False        ; turn off box between colors

;
; Define contour levels and their colors.
;
  cnres@cnLevelSelectionMode = "ManualLevels"
   cnres@cnLevelSpacingF      = 5.0
  if(lv .eq. 0) then
  cnres@cnMinLevelValF       = 285
  cnres@cnMaxLevelValF       = 315
  end if
  if(lv.gt.0 .and. lv .lt.3)then
  cnres@cnMinLevelValF       = 280
  cnres@cnMaxLevelValF       = 310
  end if
  if(lv.gt.2 .and. lv .lt.5)then
  cnres@cnMinLevelValF       = 245
  cnres@cnMaxLevelValF       = 275
  cnres@cnLevelSpacingF      = 3.0
  end if
  if(lv.gt.4.and. lv .lt.8)then
  cnres@cnMinLevelValF       = 210
  cnres@cnMaxLevelValF       = 240
  cnres@cnLevelSpacingF      = 3.0
  end if
  if(lv.gt.7)then
  cnres@cnMinLevelValF       = 200
  cnres@cnMaxLevelValF       = 230
  cnres@cnLevelSpacingF      = 3.0
  end if
;
  cnres@cnMinLevelValF       = 210
  cnres@cnMaxLevelValF       = 315
  cnres@cnLevelSpacingF      = 1.0
;  cnres@cnFillPalette = "sunshine_9lev"
  cnres@cnFillPalette = "BlAqGrYeOrRe"
;  cnres@gsnSpreadColors = True
;  cnres@gsnSpreadColorStart = 18
;  cnres@gsnSpreadColorEnd = 42

;
; Make sure contours are drawn in "predraw" phase.
;
  cnres@cnFillDrawOrder      = "Predraw"

;
; Control appearance of map.
;
;;;;  map
      LAT0=umf->lat
      LON0=umf->lon
      ds1= dimsizes(umf->lat)
      ds2= dimsizes(umf->lon)
      mpres@mpLimitMode ="Corners"
      mpres@mpLeftCornerLatF = LAT0({5})
      mpres@mpLeftCornerLonF = LON0({60})
      mpres@mpRightCornerLatF = LAT0({60})
      mpres@mpRightCornerLonF = LON0({150})
      ;for map
      mpres@mpFillOn =False
      mpres@mpGeophysicalLineColor      = "Black"
      mpres@mpGridLineColor             = "Black"
      mpres@mpLimbLineColor             = "Black"
      mpres@mpNationalLineColor         = "Black"
      mpres@mpPerimLineColor            = "Black"
      mpres@mpUSStateLineColor          = "Black"
      mpres@mpOutlineSpecifiers =(/"China"/)

  mpres@mpGridAndLimbDrawOrder = "Predraw"
  mpres@tiMainString           ="Wind Field and Humidity"; title
;
; Create, but don't draw, a vector, streamline, contour
; and map plot. Stride the vector and contour data to
; thin it a little.
;
  vcid = gsn_csm_vector_scalar(wks,su(lv,:,:),sv(lv,:,:), \
                                   sh(lv,:,:),vcres)
;  stid = gsn_csm_streamline(wks,u500(0,:,:),v500(0,:,:),stres)
  cnid = gsn_csm_contour(wks,st(lv,:,:),cnres)
  mpid = gsn_csm_map(wks,mpres)

;
; Overlay contour, streamline, and vector plots on the map plot.
;
  overlay(mpid,cnid)
;  overlay(mpid,stid)
  overlay(mpid,vcid)
;  draw(mpid)
;  frame(wks)

  maximize_output(wks,True)

   end do  ;;;; lv loop

 end




