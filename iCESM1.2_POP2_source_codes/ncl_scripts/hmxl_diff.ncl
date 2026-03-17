

;----------------------------------------------------------------------
; addfiles_1.ncl
;
; Concepts illustrated:
;   - Reading a variable off several NetCDF files using 'addfiles' and 'cat' mode
;   - To view the plot xdg-open xy.png &
;----------------------------------------------------------------------
; The files read by this example are called pottmp.YYYY.nc, where
; 'YYYY' represents a four-digit year like '1980'.
;----------------------------------------------------------------------



;load "/home/hbkheves/panel_two_sets.ncl"


begin


  ;pop_pi = addfile ("/mnt/d/lgm_data/zhxml_summer_pi_idemix_4900_4654_4703_timmean.nc", "r")

   pop_pi = addfile ("/mnt/d/lgm_data2/hxml_pi_idemix_isotope_continued_5220_5269_winter_timmean.nc", "r")

   ;pop_lgm = addfile ("/mnt/d/lgm_data2/hxml_lgm_idemix_isotope_continued_5920_5969_summer_timmean.nc", "r")


   ;pop_pi_tamas = addfile ("/mnt/d/lgm_data/montly_meanof50years_hmxl_piControl.p875.001.pop.h.4200_4249_summer_timmean.nc", "r")

   ;pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/zhxml_test_cmake_1830_1879_summer_timmean.nc", "r")

   
    pop_lgm = addfile ("/mnt/d/lgm_data2/hxml_lgm_idemix_isotope_continued_5920_5969_winter_timmean.nc", "r")

    pop_pi_tamas = addfile ("/mnt/d/lgm_data/montly_meanof50years_hmxl_piControl.p875.001.pop.h.4200_4249_winter_timmean.nc", "r")

   pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/zhxml_test_cmake_1830_1879_winter_timmean.nc", "r")   



   wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"
wks  = gsn_open_wks(wks_type,"hmxl_diff_winter_lastnew")

;  cmap = read_colormap_file("cmocean_balance")       ; read color data
   cmap = read_colormap_file("MPL_RdYlBu")
   ncmap = dimsizes(cmap(:,0))

  res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = cmap(::-1,:)         ; set color map
;  res@cnFillPalette        = read_colormap_file("MPL_rainbow")
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels

res@vpWidthF           = 0.15
res@vpHeightF          = 0.20
res@tmXBLabelsOn  = True
res@tmXBLabelFontHeightF = 0.003
res@tmYLLabelFontHeightF = 0.003


pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(2,graphic)
;plots = new(2,graphic)

res@mpLandFillColor = "lightgrey"

res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "
; res@gsnLeftString =  "Vertical Diffusivity"
; res@gsnRightString           = "(m^2/s)"
 res@gsnCenterString = "  "
res@tiMainString  = " "             ; add title
  ;res@gsnMaximize         = True                 ; make as large as possible
res@tmXBLabelFontHeightF = 0.0025
res@tmYLLabelFontHeightF = 0.0025
res@lbLabelBarOn        = False
res@cnMissingValFillColor = "grey56"
res@mpLandFillColor  = "grey56"



;res@cnMissingValFillPattern = "SolidFill"


;res@cnMissingValFillColor = "lightgrey"
;res@cnFillPalette        = read_colormap_file("MPL_rainbow")
res@vpWidthF           = 0.1
res@vpHeightF          = 0.18

res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True
res@tmXBTickSpacingF     = 45.0

res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/ -0.75, -0.5, -0.25, -0.05, -0.0075, -0.0005, -0.0025, -0.0005, 0 , 0.0005, 0.0025, 0.005, 0.0075, 0.05, 0.25, 0.50, 0.75/)
;res@cnLevels    = (/ 0 , 0.000005, 0.000025, 0.00005, 0.000075, 0.0005, 0.0025, 0.005, 0.0075, 0.05, 0.25, 0.50, 0.75, 1/)
res@cnLevels    = ispan(-150,150,30)
;res@cnLevels    = ispan(-20,20,2)


res1 = res

res1@tmYLLabelsOn  = False
res1@tmXBLabelsOn  = True
;res1@tmXBTickSpacingF     = 30.0


res2 = res1

res2@tmYLLabelsOn  = True
res2@tmXBLabelsOn  = True
;res2@tmXBTickSpacingF     = 45.0

res3 =res2
res3@tmYLLabelsOn  = False
res3@tmXBLabelsOn  = True
;res3@tmXBTickSpacingF     = 45.0


  lon1   = pop_pi->TLONG                   ; read in original lat/lon data
  lat1   = pop_pi->TLAT

  ; ================================================;
; when a coordinate array contains one more element than the corresponding
; dimension in the data array, the cell boundary is automatically turned on,
; which causes discrete rasterization to use the coordinate values as
; divisions between data cells, rather than to treat each coordinate value
; as a cell center.
; ================================================;
  dims      = dimsizes(lon1)         ; get dimension sizes
  dims(1)   = dims(1) + 1               ; increase longitude by 1
  nlat      = dims(0)
  nlon      = dims(1)

  lon      = new(dims,typeof(lon1)) ; create new arrays
  lat      = new(dims,typeof(lat1))

  lat(:,1:nlon-1)  = lat1        ; fill in arrays
  lon(:,1:nlon-1)  = lon1
  lon(:,0)         = lon1(:,nlon-2)
  lat(:,0)         = lat1(:,nlon-2)




  pi_hmxl_tamas     = pop_pi_tamas->HMXL         
 
  pi_hmxl_tamas_m2 = 0.01 * pi_hmxl_tamas

  pi_hmxl_tamas_mean = dim_avg_n(pi_hmxl_tamas_m2,0)


  pi_hmxl_tamas_mean@lat2d = lat         ; used in graphics
  pi_hmxl_tamas_mean@lon2d = lon


  lgm_hmxl_tamas     = pop_lgm_tamas->HMXL       

  lgm_hmxl_tamas_m2 = 0.01 * lgm_hmxl_tamas

  lgm_hmxl_tamas_mean = dim_avg_n(lgm_hmxl_tamas_m2,0)  


   print(max(lgm_hmxl_tamas_mean))
    print(min(lgm_hmxl_tamas_mean))  

 

  pi_hmxl     = pop_pi->HMXL

  pi_hmxl_m2 = 0.01 * pi_hmxl

  pi_hmxl_mean = dim_avg_n(pi_hmxl_m2,0)


  pi_hmxl_mean@lat2d = lat    ; used in graphics
  pi_hmxl_mean@lon2d = lon


  lgm_hmxl     = pop_lgm->HMXL

  lgm_hmxl_m2 = 0.01 * lgm_hmxl

  lgm_hmxl_mean = dim_avg_n(lgm_hmxl_m2,0)

  
  lgm_hmxl_mean@lat2d = lat       ; used in graphics
  lgm_hmxl_mean@lon2d = lon


   sub1 = pi_hmxl_mean - pi_hmxl_tamas_mean
   sub2 = lgm_hmxl_mean - lgm_hmxl_tamas_mean


   sub1@lat2d = lat
   sub1@lon2d = lon

   sub2@lat2d = lat
   sub2@lon2d = lon
  
;   print(min(pi_bflux_depth))
;   print(max(pi_bflux_depth)) 
;   print(min(pi_bflux2_depth))
;   print(max(pi_bflux2_depth))   


;   print(min(lgm_bflux_depth))
;   print(max(lgm_bflux_depth))
;   print(min(lgm_bflux2_depth))
;   print(max(lgm_bflux2_depth))

;;  res@gsnPolar             = "NH"               ; specify the hemisphere
;;  res1@gsnPolar             = "NH"
 ;; res2@gsnPolar             = "SH"
;;  res3@gsnPolar             = "SH"
;  res@mpMinLatF = 60.
;  res@mpMaxLatF = 90.
;  res1@mpMinLatF = 60.
;  res1@mpMaxLatF = 90.
  
 ; res2@mpMinLatF = -60.
 ; res2@mpMaxLatF = -90.
 ; res3@mpMinLatF = -60.
 ; res3@mpMaxLatF = -90.
  

;res@trGridType = "TriangularMesh"     
;  res@gsnPolarLabelFontHeightF = 0.0035
; res1@gsnPolarLabelFontHeightF = 0.0035
;res2@gsnPolarLabelFontHeightF = 0.0035
;res3@gsnPolarLabelFontHeightF = 0.0035 

;    plot(0) = gsn_csm_contour_map_polar(wks,pi_hmxl_mean,res)
;    plot(1) = gsn_csm_contour_map_polar(wks,lgm_hmxl_mean,res1)
;    plot(2) = gsn_csm_contour_map_polar(wks,sub1,res2)
;    plot(3) = gsn_csm_contour_map_polar(wks,sub2,res3)


    plot(0) = gsn_csm_contour_map(wks,sub1(1:nlat-1,:),res)
    plot(1) = gsn_csm_contour_map(wks,sub2(1:nlat-1,:),res1)




  resP = True
  resP@gsnMaximize         = True                 ; make as large as possible   
;  resP@gsnPaperOrientation = "portrait"           ; force portrait
;  resP@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  resP@gsnMaximize        = True
  resP@gsnFrame         = False               ; save panel until both ready
;  resP@gsnPaperOrientation = "landscape"
  resP@gsnPanelLabelBar   = True
  ;resP@pmLabelBarWidthF   = 0.07
  resP@pmLabelBarHeightF = 0.056
   resP@pmLabelBarWidthF   = 0.65
  ; resP@pmLabelBarHeightF = 0.33    
   resP@pmLabelBarOrthogonalPosF = -0.047
  ; resP@pmLabelBarParallelPosF = -0.08
   resP@lbTitleExtentF = 0.01
;resP@lbOrientation        = "Vertical"
;resP@lbLabelPosition      = "Right"           ; label position

resP@lbTitleString = "Difference in Mixed Layer Depth [m]"
resP@lbTitleFontHeightF = 0.0121
;resP@lbTitlePosition     = "Top"
;resP@lbTitleDirection     = "Across"                          ; letter angle
;;resP@lbTitleAngleF        = 270.                               ; title angle
;;resP@lbTitleDirection = "Down"
;;resP@lbTitleOffsetF = 0.025
;resP@lbTitleFont   = "helvetica-bold"

;resP@gsnPanelYWhiteSpacePercent = 5
resP@gsnMaximize         = True                 ; make as large as possible
resP@lbLabelFontHeightF= 0.010
;resP@gsnPanelFigureStringsFontHeightF = 0.011
;resP@gsnPanelFigureStrings= (/"PI-w/oIDEMIX", "LGM-w/oIDEMIX", "PI-w/IDEMIX", "LGM-w/IDEMIX"/)



;resP@gsnPanelTop = 0.85
;resP@gsnPanelBottom = 0.65

;gsn_panel(wks, plot, (/1,2/), resP)


resP1 = True
  resP1@gsnMaximize         = True                 ; make as large as possible
  resP1@gsnPaperOrientation = "portrait"           ; force portrait
;  resP1@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  resP1@gsnMaximize        = True
  resP1@gsnFrame         = False               ; save panel until both ready
;  resP1@gsnPaperOrientation = "landscape"
  resP1@gsnPanelLabelBar   = True
  ;resP1@pmLabelBarWidthF   = 0.07
  ;resP1@pmLabelBarHeightF = 0.23
  ; resP1@pmLabelBarWidthF   = 0.07
  ; resP1@pmLabelBarHeightF = 0.33
   resP1@pmLabelBarOrthogonalPosF = 0.015
   ;resP1@pmLabelBarParallelPosF = -0.01
   resP1@lbTitleExtentF = 0.01
resP1@lbOrientation        = "Vertical"
resP1@lbLabelPosition      = "Right"           ; label position

resP1@lbTitleString = " [m]"
resP1@lbTitleFontHeightF = 0.010
resP1@lbTitlePosition     = "Right"
resP1@lbTitleDirection     = "Across"                          ; letter angle
resP1@lbTitleAngleF        = 270.                               ; title angle
;resP1@lbTitleDirection = "Down"
resP1@lbTitleOffsetF = 0.03
;resP1@gsnPanelYWhiteSpacePercent = 5
resP1@gsnMaximize         = True                 ; make as large as possible
resP1@lbLabelFontHeightF= 0.009
resP1@gsnPanelFigureStringsFontHeightF = 0.009

;resP1@lbTitleString = "Depth of 17degC isotherm [m]"


;resP1@PanelTop = 0.55
;resP1@gsnPanelBottom = 0.45

;gsn_panel(wks, plots, (/1,2/), resP1)

;  pres1 = True        ; Resources specific to first set of paneled plots
;  pres2 = True        ; Resources specific to second set of paneled plots
;  pres  = True        ; For setting the area of paneled plots

;  pres1@pmLabelBarOrthogonalPosF = -0.02    ; Move closer to plots
;  pres2@pmLabelBarOrthogonalPosF = -0.02

; pres@gsnPanelTop    = 1.0  ; set these if desired
; pres@gsnPanelBottom = 0.0  ; must be between 0 and 1
; pres@gsnPanelLeft   = 0.0
; pres@gsnPanelRight  = 1.0

;  panel_two_sets(wks,plot,plots,(/1,2/),(/2,1/),pres1,pres2,pres)


txres2 = True
txres2@txFontHeightF = 0.014
txres2@txFont   = "helvetica-bold"

gsn_text_ndc(wks,"Mixed Layer Depth (w/IDEMIX - w/oIDEMIX)",0.24,0.98,txres2)
gsn_text_ndc(wks,"m",0.82,0.98,txres2)


gsn_text_ndc(wks,"PI (w/IDEMIX - w/oIDEMIX)",0.26,0.69,txres2)
gsn_text_ndc(wks,"LGM (w/IDEMIX - w/oIDEMIX)",0.77,0.69,txres2)



gsn_panel(wks, plot, (/1,2/), resP)

;gsn_panel(wks, plot, (/1,2/), resP)

;gsn_panel(wks, plots, (/1,2/), resP1)


end
