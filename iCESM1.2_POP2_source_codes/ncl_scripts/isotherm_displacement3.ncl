

;*************************************************
; NCL Graphics: iso_1.ncl
;
; Concepts illustrated: 
;   - Using "int2p_n_Wrap" to interpolate to user specified levels
;   - Drawing contour plot over a map
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;************************************************
; read POP netCDF file
;************************************************

 pop_pi_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
   pop_lgm_tamasd = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")
    pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")
    ;;pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/sel_buo_test_cmake_1750_1799_timmean.nc", "r")


   pop_pi_n2_idemix = addfile ("/mnt/d/lgm_data/work4_work5_n2_pi_idemix_isotope_continued_5000_5099_timmean.nc", "r")
   pop_lgm_n2_idemix = addfile ("/mnt/d/lgm_data/work4_lgm_5849_5899_timmean.nc", "r")


   pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
     pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")
  ;pop_lgm = addfile ("/mnt/d/lgm_data/surf_buo_dinezio_continued.pop.h.1810_1859_timmean.nc", "r")

    pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")


   pop_pi_tamas_temp_salt = addfile ("/mnt/d/lgm_data2/everything_except_moc_piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
   pop_lgm_tamas_temp_salt = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")
   pop_pi_temp_salt = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
   pop_lgm_temp_salt = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")


   pop_lgm_dinezio_temp_salt_wo = addfile ("/mnt/d/lgm_data2/everything_except_moc_LGM21ka_frUsp.p875.001.pop.h.5400_5449_timmean.nc", "r") ;surf_buo_LGM21ka_frUsp.p875.001.pop.h.5450_5499.deflated_timmean.nc

   pop_lgm_dinezio_temp_salt_w = addfile ("/mnt/d/lgm_data/select_dinezio_continued_1810_1859_timmean.nc", "r")


; read temperature (temp(time, depth, lat, lon))
temp1 = pop_pi_tamas_temp_salt->TEMP(:,:,:,:)
temp2 = pop_pi_temp_salt->TEMP(:,:,:,:)
temp3 = pop_lgm_tamas_temp_salt->TEMP(:,:,:,:)
temp4 = pop_lgm_temp_salt->TEMP(:,:,:,:)
temp5 = pop_lgm_dinezio_temp_salt_wo->TEMP(:,:,:,:)
temp6 = pop_lgm_dinezio_temp_salt_w->TEMP(:,:,:,:)


region_pi  = pop_pi_region->REGION_MASK
region_lgm  = pop_lgm_region->REGION_MASK

region_atlantic_pi= new((/384,320/),float)
region_atlantic_lgm= new((/384,320/),float)

region_atlantic_pi        = where(region_pi.eq.6,1,9.969209968386869e+36)

region_atlantic_lgm        = where(region_lgm.eq.6,1,9.969209968386869e+36)



pi_temp_tamas_na = mask(temp1,region_atlantic_pi,1)

pi_temp_tamas_mean    = dim_avg_n(pi_temp_tamas_na(0,:,:,:),2)

pi_temp_tamas_mean2    = dim_avg_n(temp1(0,:,:,:),2)


pi_temp_na = mask(temp2,region_atlantic_pi,1)

pi_temp_mean    = dim_avg_n(pi_temp_na(0,:,:,:),2)

pi_temp_mean2    = dim_avg_n(temp2(0,:,:,:),2)


lgm_temp_tamas_na = mask(temp3,region_atlantic_lgm,1)

lgm_temp_tamas_mean    = dim_avg_n(lgm_temp_tamas_na(0,:,:,:),2)

lgm_temp_tamas_mean2    = dim_avg_n(temp3(0,:,:,:),2)


lgm_temp_na = mask(temp4,region_atlantic_lgm,1)

lgm_temp_mean    = dim_avg_n(lgm_temp_na(0,:,:,:),2)

lgm_temp_mean2    = dim_avg_n(temp4(0,:,:,:),2)



lgm_temp_tamas_nad = mask(temp5,region_atlantic_lgm,1)

lgm_temp_tamas_meand    = dim_avg_n(lgm_temp_tamas_nad(0,:,:,:),2)

lgm_temp_tamas_meand2    = dim_avg_n(temp5(0,:,:,:),2)


lgm_temp_nad = mask(temp6,region_atlantic_lgm,1)

lgm_temp_meand    = dim_avg_n(lgm_temp_nad(0,:,:,:),2)

lgm_temp_meand2    = dim_avg_n(temp6(0,:,:,:),2)



sub1 = pi_temp_mean - pi_temp_tamas_mean
sub2 = lgm_temp_mean - lgm_temp_tamas_mean
sub3 = lgm_temp_meand - lgm_temp_tamas_meand


sub4 = pi_temp_mean2 - pi_temp_tamas_mean2
sub5 = lgm_temp_mean2 - lgm_temp_tamas_mean2
sub6 = lgm_temp_meand2 - lgm_temp_tamas_meand2



depth =  pop_pi->z_t          ; or whatever your depth variable name is

  lon   = pop_pi->TLONG                   ; read in original lat/lon data
  lat   = pop_pi->TLAT

; dimsizes(temp1) should return: (nt, nz, ny, nx)
;print("TEMP dims: " + dimsizes(temp1))

Tiso = 10.0   ; 10°C isotherm

;************************************************
; Interpolate to specific [constant] depths
; The returned vertical coordinate is Z_T
;************************************************
 ; zlev = (/ 735, 12580/)    ; cm
 
  zlev = new((/60/), typeof(depth), -999) 
  nt = dimsizes(depth) 

 do i = 0, nt-1
  zlev(i) = depth(i)    ; cm
 end do  

  print(zlev)

  zlev@units = depth@units
  zlev!0 = "zlev"
  iso1 = int2p_n_Wrap(depth,temp1, zlev, 0, 1)

  printVarSummary(iso1) ; [time | 1] x [zlev | 2] x [nlat | 384] x [nlon | 320]  
  printMinMax(iso1, True) 

  iso2 = int2p_n_Wrap(depth,temp2, zlev, 0, 1)
  iso3 = int2p_n_Wrap(depth,temp3, zlev, 0, 1)
  iso4 = int2p_n_Wrap(depth,temp4, zlev, 0, 1)
  iso5 = int2p_n_Wrap(depth,temp5, zlev, 0, 1)
  iso6 = int2p_n_Wrap(depth,temp6, zlev, 0, 1)


;region_pi  = pop_pi_region->REGION_MASK
;region_lgm  = pop_lgm_region->REGION_MASK

;region_atlantic_pi= new((/384,320/),float)
;region_atlantic_lgm= new((/384,320/),float)

;region_atlantic_pi        = where(region_pi.eq.6,1,9.969209968386869e+36)

;region_atlantic_lgm        = where(region_lgm.eq.6,1,9.969209968386869e+36)



pi_tamas_na = mask(iso1,region_atlantic_pi,1)

printVarSummary(pi_tamas_na)

pi_tamas_mean    = dim_avg_n(pi_tamas_na(0,:,:,:),2)

pi_tamas_mean2    = dim_avg_n(iso1(0,:,:,:),2)

printVarSummary(pi_tamas_mean)


pi_na = mask(iso2,region_atlantic_pi,1)
pi_mean    = dim_avg_n(pi_na(0,:,:,:),2) 

pi_mean2    = dim_avg_n(iso2(0,:,:,:),2)


lgm_tamas_na = mask(iso3,region_atlantic_lgm,1)
lgm_tamas_mean    = dim_avg_n(lgm_tamas_na(0,:,:,:),2)

lgm_tamas_mean2    = dim_avg_n(iso3(0,:,:,:),2)


lgm_na = mask(iso4,region_atlantic_lgm,1)
lgm_mean    = dim_avg_n(lgm_na(0,:,:,:),2)

lgm_mean2    = dim_avg_n(iso4(0,:,:,:),2)


lgm_tamas_nad = mask(iso5,region_atlantic_lgm,1)
lgm_tamas_meand    = dim_avg_n(lgm_tamas_nad(0,:,:,:),2)

lgm_tamas_meand2    = dim_avg_n(iso5(0,:,:,:),2)

lgm_nad = mask(iso6,region_atlantic_lgm,1)
lgm_meand    = dim_avg_n(lgm_nad(0,:,:,:),2)

lgm_meand2    = dim_avg_n(iso6(0,:,:,:),2)




;************************************************
; Interpolate to specific [constant] TEMP levels
; The default returned vertical coordinate is Z_T but change to 'tlev'
;************************************************
 ; tlev    = (/ 13.5 , 18/)    ; same units [here, C] as TEMP
 ; tlev@units = x@units
;  tlev!0     = "tlev"
;
;  depth   = conform(x, z_t, 1)
;  copy_VarCoords(x, depth)
;  printVarSummary(depth)
  
;  isot    = int2p_n_Wrap(x,depth, tlev, 0, 1)

 ;printVarSummary(isot) ; [time | 1] x [tlev | 2] x [nlat | 384] x [nlon | 320]  
 ;printMinMax(isot, True) 

;************************************************ 
; create plot
;************************************************
  plot = new( 3, "graphic")

  wks = gsn_open_wks("pdf","iso")                 ; send graphics to PNG file

  res                      = True                 ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "temp_diff_18lev"   ;"BlAqGrYeOrReVi200"  ; set color map
  res@cnLinesOn            = False                ; turn off contour lines
  ;res@mpCenterLonF         = 210                  ; center map at 210
  ;res@lbLabelAutoStride    = True  
  ;res@lbOrientation        = "vertical"
  ;res@gsnAddCyclic         = True                 

res@mpLandFillColor = "grey56"
res@cnMissingValFillColor = "grey56"

res@lbLabelBarOn        = True
res@cnInfoLabelOn       = False
res@description = "  "

pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False

  
res@vpWidthF           = 0.15
res@vpHeightF          = 0.10
;res@lbTitleString = "Salt [g/kg]"
res@lbTitleFontHeightF= 0.004
;res@lbTitleOffsetF = 0.01

; res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res@cnMinLevelValF            = -1.             ; min level
;    res@cnMaxLevelValF            =  1.             ; max level
;    res@cnLevelSpacingF           =  0.1            ; contour interval


res@tmXBTickSpacingF = 30.0 ; North Atlantic
;re1@tmXBTickSpacingF = 10.0
res@tmYLTickSpacingF = 30.0

res@gsnYAxisIrregular2Linear = True
res@trYReverse           = True       ; reverses y-axis
res@tmYLLabelsOn  = True
res@cnMissingValFillPattern = "SolidFill"

res@tmXBTickSpacingF     = 15.0
res@vpWidthF           = 0.15
res@vpHeightF          = 0.1


res@gsnYAxisIrregular2Linear = True
res@trYReverse           = True       ; reverses y-axis

    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            = -2.0             ; min level
    res@cnMaxLevelValF            =  2.0             ; max level
    res@cnLevelSpacingF           =  0.2            ; contour interval


res@tmXBTickSpacingF = 30.0 ; North Atlantic
;res@tmXBTickSpacingF = 10.0
res@tmYLTickSpacingF = 30.0


  resP = True
  resP@gsnMaximize         = True                 ; make as large as possible   
  ;resP@gsnPaperOrientation = "portrait"           ; force portrait
  resP@gsnMaximize        = False
  resP@gsnFrame         = False    
 
  ;isox@lat2d = b->TLAT         ; used in graphics
  ;isox@lon2d = b->TLONG



;  do kl=0,dimsizes(zlev)-1
;     res@gsnCenterString   = "lev="+zlev(kl)+" cm"
;     plot(kl) = gsn_csm_contour_map(wks,isox(0,kl,:,:),res)  
;  end do
 ; resP@gsnPanelMainString  = "TEMP at User Specified Depths "
 ; gsn_panel(wks, plot, (/2,1/), resP)
 ; delete(isox)

 ; isot@lat2d = b->TLAT         ; used in graphics
 ; isot@lon2d = b->TLONG

;  do kt=0,dimsizes(tlev)-1
;     res@gsnCenterString   = "tlev="+tlev(kt)+" C"
;     plot(kt) = gsn_csm_contour_map(wks,isot(0,kt,:,:),res)   
;  end do


nlt = 9
nlt2 = 224
nlt3 = 224


latz = pop_pi->TLAT(:,0)

lat_a    = pop_pi->TLAT(:,nlt)
lat_atlantic        = lat_a
 
printVarSummary(lat_atlantic)

z_t           = pop_pi->z_t
z_t           = z_t*1.E-2       ; convert from cm to m

z_t@units     = "m"             ; fix units attribute
z_t@tiYAxisPosition = "top"
z_t@tiYAxisFontHeightF = 0.018
z_t@long_name = "Depth (m)"; overwrite long name




sub1!0    = "z_t"
sub1&z_t  =  z_t
sub1!1    = "lat_t"
sub1&lat_t=  lat_atlantic



sub2!0    = "z_t"
sub2&z_t  =  z_t
sub2!1    = "lat_t"
sub2&lat_t=  lat_atlantic



sub3!0    = "z_t"
sub3&z_t  =  z_t
sub3!1    = "lat_t"
sub3&lat_t=  lat_atlantic


sub4!0    = "z_t"
sub4&z_t  =  z_t
sub4!1    = "lat_t"
sub4&lat_t=  latz


sub5!0    = "z_t"
sub5&z_t  =  z_t
sub5!1    = "lat_t"
sub5&lat_t=  latz


sub6!0    = "z_t"
sub6&z_t  =  z_t
sub6!1    = "lat_t"
sub6&lat_t=  latz


pi_tamas_mean2!0    = "z_t"
pi_tamas_mean2&z_t  =  z_t
pi_tamas_mean2!1    = "lat_t"
pi_tamas_mean2&lat_t=  latz
;pi_tamas_mean&lat_t=  lat(:,0)

pi_mean2!0    = "z_t"
pi_mean2&z_t  =  z_t
pi_mean2!1    = "lat_t"
pi_mean2&lat_t=  latz

lgm_tamas_mean2!0    = "z_t"
lgm_tamas_mean2&z_t  =  z_t
lgm_tamas_mean2!1    = "lat_t"
lgm_tamas_mean2&lat_t=  latz


lgm_mean2!0    = "z_t"
lgm_mean2&z_t  =  z_t
lgm_mean2!1    = "lat_t"
lgm_mean2&lat_t=  latz


lgm_tamas_meand2!0    = "z_t"
lgm_tamas_meand2&z_t  =  z_t
lgm_tamas_meand2!1    = "lat_t"
lgm_tamas_meand2&lat_t=  latz


lgm_meand2!0    = "z_t"
lgm_meand2&z_t  =  z_t
lgm_meand2!1    = "lat_t"
lgm_meand2&lat_t=  latz




pi_tamas_mean!0    = "z_t"
pi_tamas_mean&z_t  =  z_t
pi_tamas_mean!1    = "lat_t"
pi_tamas_mean&lat_t=  lat_atlantic
;pi_tamas_mean&lat_t=  lat(:,0)

pi_mean!0    = "z_t"
pi_mean&z_t  =  z_t
pi_mean!1    = "lat_t"
pi_mean&lat_t=  lat_atlantic


lgm_tamas_mean!0    = "z_t"
lgm_tamas_mean&z_t  =  z_t
lgm_tamas_mean!1    = "lat_t"
lgm_tamas_mean&lat_t=  lat_atlantic


lgm_mean!0    = "z_t"
lgm_mean&z_t  =  z_t
lgm_mean!1    = "lat_t"
lgm_mean&lat_t=  lat_atlantic


lgm_tamas_meand!0    = "z_t"
lgm_tamas_meand&z_t  =  z_t
lgm_tamas_meand!1    = "lat_t"
lgm_tamas_meand&lat_t=  lat_atlantic


lgm_meand!0    = "z_t"
lgm_meand&z_t  =  z_t
lgm_meand!1    = "lat_t"
lgm_meand&lat_t=  lat_atlantic

pi_temp_mean!0    = "z_t"
pi_temp_mean&z_t  =  z_t
pi_temp_mean!1    = "lat_t"
pi_temp_mean&lat_t=  lat_atlantic

lgm_temp_mean!0    = "z_t"
lgm_temp_mean&z_t  =  z_t
lgm_temp_mean!1    = "lat_t"
lgm_temp_mean&lat_t=  lat_atlantic

lgm_temp_meand!0    = "z_t"
lgm_temp_meand&z_t  =  z_t
lgm_temp_meand!1    = "lat_t"
lgm_temp_meand&lat_t=  lat_atlantic


  plot(0) = gsn_csm_contour(wks,sub4,res)
  plot(1) = gsn_csm_contour(wks,sub5,res)
  plot(2) = gsn_csm_contour(wks,sub6,res)



  overlay_res                          = True
  overlay_res@cnInfoLabelOn       = False
  overlay_res@description = "  "
  overlay_res@gsnDraw              = False                ; don't draw yet
  overlay_res@gsnFrame             = False                ; don't advance frame yet
  overlay_res@cnLineDrawOrder         = "PostDraw"
  overlay_res@cnLinesOn               = True
  overlay_res@cnMonoLineColor  = False 
  overlay_res@cnLineColors             = span_color_rgba ("BlAqGrYeOrReVi200",8)

;  overlay_res@cnLineColors             = "WhViBlGrYeOrRe"
 ; overlay_res@cnLineColor             = "black"
  overlay_res@cnLineThicknessF        = 1.5
  ;overlay_res@cnLevelSelectionMode    = "ExplicitLevels"
  ;overlay_res@cnLevels                = (/5/)     
  ;overlay_res@cnLevels                = (/1/)
  ;overlay_res@sfXArray                = lon
  ;overlay_res@sfYArray                = lat
  ;overlay_res@cnLineLabelPlacementMode = "ConstantSpacing"
  overlay_res@tiMainString  = " "             ; add title
  overlay_res@cnLineDashPattern = 2
  overlay_res@gsnYAxisIrregular2Linear = True
  overlay_res@trYReverse           = True       ; reverses y-axis
  overlay_res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    overlay_res@cnMinLevelValF            = 0.             ; min level
    overlay_res@cnMaxLevelValF            = 25.             ; max level
    overlay_res@cnLevelSpacingF           =  5            ; contour interval
 


 overlay_res2 =  overlay_res
 ;overlay_res2@cnLineColor             = "green"
 delete(overlay_res2@cnLineDashPattern)
 overlay_res2@cnLineDashPattern = 0


  ax1 = gsn_csm_contour(wks,pi_tamas_mean2,overlay_res)
  ax2 = gsn_csm_contour(wks,pi_mean2,overlay_res2)
  ax3 = gsn_csm_contour(wks,lgm_tamas_mean2,overlay_res)
  ax4 = gsn_csm_contour(wks,lgm_mean2,overlay_res2)
  ax5 = gsn_csm_contour(wks,lgm_tamas_meand2,overlay_res)
  ax6 = gsn_csm_contour(wks,lgm_meand2,overlay_res2)

  print(min(lgm_tamas_mean2))
  print(min(lgm_mean2))

  ;  print(pi_tamas_mean)

    overlay(plot(0), ax1)
    overlay(plot(0), ax2)
  
    overlay(plot(1), ax3)
    overlay(plot(1), ax4)

    overlay(plot(2), ax5)
    overlay(plot(2), ax6)
  
   

  ;resP@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  gsn_panel(wks, plot, (/3,1/), resP)

  ;draw(plot)                     ; note we are drawing the first one!
  frame(wks)

end

