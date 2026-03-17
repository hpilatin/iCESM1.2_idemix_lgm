

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



   ;pop_pi_heat_tamas = addfile ("/mnt/d/lgm_data2/b.ie12.Bi1850C5BPRDCISO.f19_g16.piControl.p875.001_N_HEAT_4149-4199_deflated_timmean.nc", "r")
   ;pop_pi_salt_tamas = addfile ("/mnt/d/lgm_data2/b.ie12.Bi1850C5BPRDCISO.f19_g16.piControl.p875.001_N_SALT_4149-4199_deflated_timmean.nc", "r")

   pop_lgm_heat_tamas = addfile ("/mnt/d/lgm_data2/b.ie12.Bi1850C5BPRDCISO.f19_g16.LGM21ka_frUsp.p875.001_N_HEAT_5449-5499_deflated.nc", "r")
   pop_lgm_salt_tamas = addfile ("/mnt/d/lgm_data2/b.ie12.Bi1850C5BPRDCISO.f19_g16.LGM21ka_frUsp.p875.001_N_SALT_5449-5499_deflated.nc", "r")

   ;! pop_lgm_heat_tamas = addfile ("/mnt/d/lgm_data2/select_test_cmake_1670_1719_timmean.nc", "r")
   ;! pop_lgm_salt_tamas = addfile ("/mnt/d/lgm_data2/select_test_cmake_1670_1719_timmean.nc", "r")
   
    ;pop_lgm_heat_tamas = addfile ("/mnt/d/lgm_data2/sel_buo_test_cmake_1750_1799_timmean.nc", "r")
    ;pop_lgm_salt_tamas = addfile ("/mnt/d/lgm_data2/sel_buo_test_cmake_1750_1799_timmean.nc", "r")

    pop_pi_heat_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")
    pop_pi_salt_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")





   ;pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
  ; pop_pi = addfile ("/mnt/d/lgm_data/temp_salt_n2_pi_4949_4999_timmean.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data/temp_salt_n2_lgm_5849_5899_timmean.nc", "r")
  

   pop_pi = addfile ("/mnt/d/lgm_data/temp_salt_n2_lgm_5849_5899_timmean.nc", "r")
   pop_lgm = addfile ("/mnt/d/lgm_data/select_dinezio_continued_1810_1859_timmean.nc", "r")

   pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")


  wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"
   
 
   wks  = gsn_open_wks(wks_type,"heat_salt_flux_lat_pi_lgm_thesis_testcmake_lastnew")

  

;  cmap = read_colormap_file("cmocean_balance")       ; read color data
   cmap = read_colormap_file("MPL_RdYlBu")
   ncmap = dimsizes(cmap(:,0))





res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = False               ; turn on color fill
  res@cnFillPalette        = cmap(::-1,:)         ; set color map
;  res@cnFillPalette        = read_colormap_file("MPL_rainbow")
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels

;res@vpWidthF           = 0.15
;res@vpHeightF          = 0.20
res@tmXBLabelsOn  = True

pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plots = new(6,graphic)



res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "
; res@gsnLeftString =  "Vertical Diffusivity"
; res@gsnRightString           = "(m^2/s)"
 res@gsnCenterString = "  "
res@tiMainString  = " "             ; add title
;res@gsnMaximize         = True                 ; make as large as possible   
;res@tmXBLabelFontHeightF = 0.0025
;res@tmYLLabelFontHeightF = 0.0025
res@lbLabelBarOn        = False

res@cnMissingValFillPattern = "SolidFill"


res@cnMissingValFillColor = "lightgrey"
;res@cnFillPalette        = read_colormap_file("MPL_rainbow")
res@vpWidthF           = 0.13
res@vpHeightF          = 0.18
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/-13,-12,-11,-10,-9,-8,-7,-6/)
res@gsnLeftString =  ""
res@gsnRightString =  ""



iw1        = pop_pi_heat_tamas->N_HEAT(:,0,0,:)
iw1_serie=dim_avg_n(iw1,0)




iw1_atl        = pop_pi_heat_tamas->N_HEAT(:,1,1,:) ;Atlantic

iw1_serie_atl=dim_avg_n(iw1_atl,0)



iw1_pac = iw1 - iw1_atl ;Pacific
iw1_serie_pac=dim_avg_n(iw1_pac,0)

  

iw2        = pop_lgm_heat_tamas->N_HEAT(:,0,0,:) ;global

iw2_serie=dim_avg_n(iw2,0)



iw2_atl        = pop_lgm_heat_tamas->N_HEAT(:,1,1,:) ;Atlantic

iw2_serie_atl=dim_avg_n(iw2_atl,0)


iw2_pac = iw2 - iw2_atl ;Pacific
iw2_serie_pac=dim_avg_n(iw2_pac,0)



iw_pi        = pop_pi->N_HEAT(:,0,0,:) ;global

iw_pi_serie=dim_avg_n(iw_pi,0)




iw_pi_atl        = pop_pi->N_HEAT(:,1,1,:) ;Atlantic

iw_pi_serie_atl=dim_avg_n(iw_pi_atl,0)



iw_pi_pac = iw_pi - iw_pi_atl ;Pacific

iw_pi_serie_pac=dim_avg_n(iw_pi_pac,0)




 region  = pop_pi_region->REGION_MASK




iw_lgm        = pop_lgm->N_HEAT(:,0,0,:) ;global


iw_lgm_serie=dim_avg_n(iw_lgm,0)


iw_lgm_atl        = pop_lgm->N_HEAT(:,1,1,:) ;Atlantic


iw_lgm_serie_atl=dim_avg_n(iw_lgm_atl,0)


iw_lgm_pac = iw_lgm - iw_lgm_atl ;Pacific

iw_lgm_serie_pac=dim_avg_n(iw_lgm_pac,0)




lat_pi = pop_pi->TLAT(285:365,8)
lat_lgm = pop_lgm->TLAT(285:365,8)





iw1_salt        = pop_pi_salt_tamas->N_SALT(:,0,0,:)

iw1_serie_salt=dim_avg_n(iw1_salt,0)*1.0e-9 ; unit conversion to m3/s




iw1_salt_atl        = pop_pi_salt_tamas->N_SALT(:,1,1,:)

iw1_serie_salt_atl=dim_avg_n(iw1_salt_atl,0)*1.0e-9 ; unit conversion to m3/s



iw1_salt_pac        = iw1_salt - iw1_salt_atl
iw1_serie_salt_pac=dim_avg_n(iw1_salt_pac,0)*1.0e-9 ; unit conversion to m3/s


iw2_salt        = pop_lgm_salt_tamas->N_SALT(:,0,0,:)

iw2_serie_salt=dim_avg_n(iw2_salt,0)*1.0e-9 ; unit conversion to m3/s



iw2_salt_atl        = pop_lgm_salt_tamas->N_SALT(:,1,1,:)

iw2_serie_salt_atl=dim_avg_n(iw2_salt_atl,0)*1.0e-9 ; unit conversion to m3/s



iw2_salt_pac        = iw2_salt - iw2_salt_atl
iw2_serie_salt_pac=dim_avg_n(iw2_salt_pac,0)*1.0e-9 ; unit conversion to m3/s




iw_pi_salt        = pop_pi->N_SALT(:,0,0,:)

iw_pi_serie_salt=dim_avg_n(iw_pi_salt,0)*1.0e-9 ; unit conversion to m3/s



iw_pi_salt_atl        = pop_pi->N_SALT(:,1,1,:)

iw_pi_serie_salt_atl=dim_avg_n(iw_pi_salt_atl,0)*1.0e-9 ; unit conversion to m3/s



iw_pi_salt_pac        = iw_pi_salt - iw_pi_salt_atl
iw_pi_serie_salt_pac=dim_avg_n(iw_pi_salt_pac,0)*1.0e-9 ; unit conversion to m3/s



iw_lgm_salt        = pop_lgm->N_SALT(:,0,0,:)

iw_lgm_serie_salt=dim_avg_n(iw_lgm_salt,0)*1.0e-9 ; unit conversion to m3/s



iw_lgm_salt_atl        = pop_lgm->N_SALT(:,1,1,:)

iw_lgm_serie_salt_atl=dim_avg_n(iw_lgm_salt_atl,0)*1.0e-9 ; unit conversion to m3/s


iw_lgm_salt_pac        = iw_lgm_salt - iw_lgm_salt_atl
iw_lgm_serie_salt_pac=dim_avg_n(iw_lgm_salt_pac,0)*1.0e-9 ; unit conversion to m3/s



serie1 =  new( (/4,395/),float)
serie1(0,:)= iw1_serie
serie1(1,:)= iw2_serie
serie1(2,:)= iw_pi_serie
serie1(3,:)= iw_lgm_serie


serie2 =  new( (/4,395/),float)
serie2(0,:)= iw1_serie_salt
serie2(1,:)= iw2_serie_salt
serie2(2,:)= iw_pi_serie_salt
serie2(3,:)= iw_lgm_serie_salt


;serie1 =  new( (/2,395/),float)
;serie1(0,:)= iw_pi_serie
;serie1(1,:)= iw_lgm_serie


;serie2 =  new( (/2,395/),float)
;serie2(0,:)= iw_pi_serie_salt
;serie2(1,:)= iw_lgm_serie_salt


serie3 =  new( (/4,395/),float)
;serie3(2,:)= iw_pi_serie_atl
;serie3(3,:)= iw_lgm_serie_atl
;serie3(0,:)= iw1_serie_atl
;serie3(1,:)= iw2_serie_atl

serie3(0,:)= iw1_serie_atl
serie3(1,:)= iw2_serie_atl
serie3(2,:)= iw_pi_serie_atl
serie3(3,:)= iw_lgm_serie_atl


serie4 =  new( (/4,395/),float)
serie4(0,:)= iw1_serie_salt_atl
serie4(1,:)= iw2_serie_salt_atl
serie4(2,:)= iw_pi_serie_salt_atl
serie4(3,:)= iw_lgm_serie_salt_atl


serie5 =  new( (/4,395/),float)
serie5(0,:)= iw1_serie_pac
serie5(1,:)= iw2_serie_pac
serie5(2,:)= iw_pi_serie_pac
serie5(3,:)= iw_lgm_serie_pac

serie6 =  new( (/4,395/),float)
serie6(0,:)= iw1_serie_salt_pac
serie6(1,:)= iw2_serie_salt_pac
serie6(2,:)= iw_pi_serie_salt_pac
serie6(3,:)= iw_lgm_serie_salt_pac



   resL = True
  resL@gsnDraw              = False                ; don't draw yet
  resL@gsnFrame             = False                ; don't advance frame yet
  resL@tiXAxisFontHeightF   = 0.022
   resL@tiYAxisFontHeightF   = 0.022
   resL@tmXBLabelFontHeightF  = 0.018
   resL@tmYLLabelFontHeightF  = 0.018
;   resL@tmYLLabelFontH@gsnDraw              = False                ; don't draw yet
    resL@tiYAxisString = " Northward Heat Transport" +"  "+"[PW]"  ; axis string
   resL@tiXAxisString = "Latitude"  ; axis string
;   resL@gsnLeftString = "1-2 km"
;  resL@xyMonoDashPattern = True
  resL@gsnLeftStringParallelPosF = 0.02

  resL@txFontHeightF      = 0.015
  resL@pmLegendDisplayMode    = "Always"            ; turn on legend
  resL@lgPerimOn              = False               ; turn off box around
  resL@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
;  resL@xyLineColors           = (/"darkorange1","forestgreen","gray40","dodgerblue2"/)
;   resL@xyLineColors           = (/"coral1","dodgerblue","firebrick3","blue4"/)

  resL@xyLineColors           = (/"dodgerblue","limegreen","blue4","green4"/)
 resL@lgPerimOn              = False                 ; no box around

;;

  resL@lgJustification        = "BottomRight"
  resL@pmLegendOrthogonalPosF = -1.48     ; Move into plot
  resL@pmLegendParallelPosF   = 1.1       ; Move to right
  resL@pmLegendWidthF         = 0.3       ; Decrease width
  resL@pmLegendHeightF        = 0.12       ; Decrease height
  resL@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL@lgLabelFontHeightF     = 0.019     ; Change the font size
; labels for the legend
 ;resL@xyExplicitLegendLabels = (/"PI - with IDEMIX", "LGM - with IDEMIX", "PI - with IDEMIX", "LGM - with IDEMIX", "PI - with IDEMIX", "LGM - with IDEMIX" /)


 resL@tmYLMode        = "Explicit"                     ; Force labels where we want them
; resL@tmYLValues      = (/-2,-1.5,-1,-0.5,0,0.5,1,1.5,2/)
 resL@xyDashPatterns   = (/1,1,0,0/)   


  resL@xyExplicitLegendLabels = (/"W.-w/oIDEMIX","D.-w/oIDEMIX","W.-w/IDEMIX","D.-w/IDEMIX"/)
;  resL@xyExplicitLegendLabels = (/"PI", "LGM"/) 
 resL@tiMainString    = "Global Ocean"   ; add titles
  


  resL@lgLabelFontThicknessF = 1.2
  resL@lgLineLabelFontThicknessF = 1.3

resL@lgLineThicknesses = 1.3
resL@gsnXRefLine = 1000
resL@gsnXRefLineDashPattern = 1
resL@gsnXRefLineColor = "black"
resL@gsnXRefLineThicknessF = 2.0

 resL@tmYLLabels = (/-2,-1.5,-1,-0.5,0,0.5,1,1.5,2/)



resL1 = True
  resL1@gsnDraw              = False                ; don't draw yet
  resL1@gsnFrame             = False                ; don't advance frame yet
  resL1@tiXAxisFontHeightF   = 0.022
   resL1@tiYAxisFontHeightF   = 0.022
   resL1@tmXBLabelFontHeightF  = 0.018
   resL1@tmYLLabelFontHeightF  = 0.018
   resL1@tiXAxisString = "Latitude"  ; axis string
;  resL1@xyMonoDashPattern = True
  resL1@gsnLeftStringParallelPosF = 0.02
  resL1@txFontHeightF      = 0.015
  resL1@pmLegendDisplayMode    = "NoCreate"            ; turn on legend
  resL1@lgPerimOn              = False               ; turn off box around
  resL1@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL1@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
;  resL1@xyLineColors           =  (/"darkorange1","forestgreen","gray40","dodgerblue2"/)

  ;resL1@xyLineColors           = (/"coral1","dodgerblue","firebrick3","blue4"/)
 
 resL1@xyLineColors           = (/"dodgerblue","limegreen","blue4","green4"/)
resL1@lgPerimOn              = False                 ; no box around
resL1@lgJustification        = "BottomRight"
  resL1@pmLegendOrthogonalPosF = -1.48     ; Move into plot
  resL1@pmLegendParallelPosF   = 1.1       ; Move to right
  resL1@pmLegendWidthF         = 0.3       ; Decrease width
  resL1@pmLegendHeightF        = 0.12       ; Decrease height
  resL1@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL1@lgLabelFontHeightF     = 0.019     ; Change the font size
 resL1@xyDashPatterns   = (/1,1,0,0/)
  resL1@xyExplicitLegendLabels = (/"PI-w/oIDEMIX","LGM-w/oIDEMIX","PI", "LGM"/)
   resL1@lgLabelFontThicknessF = 1.2
  resL1@lgLineLabelFontThicknessF = 1.2
resL1@lgLineThicknesses = 1.3
resL1@gsnXRefLine = 1000
resL1@gsnXRefLineDashPattern = 1
resL1@gsnXRefLineColor = "black"
resL1@gsnXRefLineThicknessF = 2.0


resL1@tiMainString    = "Atlantic Ocean"   ; add titles
resL1@tiYAxisString = ""
resL1@tmYLMode        = "Explicit"                     ; Force labels where we want them
;resL1@tmYLValues      = (/-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2/)

resL1@tmYLLabels = (/-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2/)




resL2 = True

resL2@gsnDraw              = False                ; don't draw yet
  resL2@gsnFrame             = False                ; don't advance frame yet
  resL2@tiXAxisFontHeightF   = 0.022
   resL2@tiYAxisFontHeightF   = 0.022
   resL2@tmXBLabelFontHeightF  = 0.018
   resL2@tmYLLabelFontHeightF  = 0.018
   resL2@tiXAxisString = "Latitude"  ; axis string
;  resL2@xyMonoDashPattern = True
  resL2@gsnLeftStringParallelPosF = 0.02
  resL2@txFontHeightF      = 0.015
  resL2@pmLegendDisplayMode    = "NoCreate"            ; turn on legend
  resL2@lgPerimOn              = False               ; turn off box around
  resL2@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL2@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
;  resL2@xyLineColors           =  (/"darkorange1","forestgreen","gray40","dodgerblue2"/)
;  resL2@xyLineColors           =  (/"coral1","dodgerblue","firebrick3","blue4"/)

 resL2@xyLineColors           = (/"dodgerblue","limegreen","blue4","green4"/)
 resL2@lgPerimOn              = False                 ; no box around
resL2@lgJustification        = "BottomRight"
  resL2@pmLegendOrthogonalPosF = -1.48     ; Move into plot
  resL2@pmLegendParallelPosF   = 1.1       ; Move to right
  resL2@pmLegendWidthF         = 0.3       ; Decrease width
  resL2@pmLegendHeightF        = 0.12       ; Decrease height
  resL2@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL2@lgLabelFontHeightF     = 0.019     ; Change the font size
 resL2@xyDashPatterns   = (/1,1,0,0/)
  resL2@xyExplicitLegendLabels = (/"PI-w/oIDEMIX","LGM-w/oIDEMIX","PI", "LGM"/)
   resL2@lgLabelFontThicknessF = 1.2
  resL2@lgLineLabelFontThicknessF = 1.2
resL2@lgLineThicknesses = 1.3
resL2@gsnXRefLine = 1000
resL2@gsnXRefLineDashPattern = 1
resL2@gsnXRefLineColor = "black"
resL2@gsnXRefLineThicknessF = 2.0

resL2@tiMainString    = "Indo-Pacific Ocean"   ; add titles
resL2@tiYAxisString = ""


resL2@tmYLMode        = "Explicit"
resL2@tmYLLabels = (/-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2/)


resL3 = True

resL3@gsnDraw              = False                ; don't draw yet
  resL3@gsnFrame             = False                ; don't advance frame yet
  resL3@tiXAxisFontHeightF   = 0.022
   resL3@tiYAxisFontHeightF   = 0.022
   resL3@tmXBLabelFontHeightF  = 0.018
   resL3@tmYLLabelFontHeightF  = 0.018
   resL3@tiXAxisString = "Latitude"  ; axis string
;  resL3@xyMonoDashPattern = True
  resL3@gsnLeftStringParallelPosF = 0.02
  resL3@txFontHeightF      = 0.015
  resL3@pmLegendDisplayMode    = "NoCreate"            ; turn on legend
  resL3@lgPerimOn              = False               ; turn off box around
  resL3@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL3@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
;  resL3@xyLineColors           =  (/"darkorange1","forestgreen","gray40","dodgerblue2"/)
; resL3@xyLineColors           =  ()
  resL3@xyLineColors           = (/"dodgerblue","limegreen","blue4","green4"/)
 resL3@lgPerimOn              = False                 ; no box around
resL3@lgJustification        = "BottomRight"
  resL3@pmLegendOrthogonalPosF = -1.48     ; Move into plot
  resL3@pmLegendParallelPosF   = 1.1       ; Move to right
  resL3@pmLegendWidthF         = 0.3       ; Decrease width
  resL3@pmLegendHeightF        = 0.12       ; Decrease height
  resL3@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL3@lgLabelFontHeightF     = 0.019     ; Change the font size
 resL3@xyDashPatterns   = (/1,1,0,0/)
  resL3@xyExplicitLegendLabels = (/"PI-w/oIDEMIX","LGM-w/oIDEMIX","PI", "LGM"/)
   resL3@lgLabelFontThicknessF = 1.2
  resL3@lgLineLabelFontThicknessF = 1.2
resL3@lgLineThicknesses = 1.3
resL3@gsnXRefLine = 1000
resL3@gsnXRefLineDashPattern = 1
resL3@gsnXRefLineColor = "black"
resL3@gsnXRefLineThicknessF = 2.0


resL3@tiYAxisString = " Northward Salt Transport" +"  "+"[PW]"  ; axis string
resL3@tiMainString    = ""   ; add titles
resL3@tiXAxisString = "Latitude"  ; axis string
resL3@tmYLMode        = "Explicit"                     ; Force labels where we want them
;resL3@tmYLValues      = ispan(-35,35,5)

resL3@tmYLLabels = ispan(30,-30,10)


resL4 = True

  resL4@gsnDraw              = False                ; don't draw yet
  resL4@gsnFrame             = False                ; don't advance frame yet
  resL4@tiXAxisFontHeightF   = 0.022
   resL4@tiYAxisFontHeightF   = 0.022
   resL4@tmXBLabelFontHeightF  = 0.018
   resL4@tmYLLabelFontHeightF  = 0.018
   resL4@tiXAxisString = "Latitude"  ; axis string
;  resL4@xyMonoDashPattern = True
  resL4@gsnLeftStringParallelPosF = 0.02
  resL4@txFontHeightF      = 0.015
  resL4@pmLegendDisplayMode    = "NoCreate"            ; turn on legend
  resL4@lgPerimOn              = False               ; turn off box around
  resL4@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL4@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
   resL4@xyLineColors           = (/"dodgerblue","limegreen","blue4","green4"/)
;  resL4@xyLineColors           =  
 resL4@lgPerimOn              = False                 ; no box around
resL4@lgJustification        = "BottomRight"
  resL4@pmLegendOrthogonalPosF = -1.48     ; Move into plot
  resL4@pmLegendParallelPosF   = 1.1       ; Move to right
  resL4@pmLegendWidthF         = 0.3       ; Decrease width
  resL4@pmLegendHeightF        = 0.12       ; Decrease height
  resL4@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL4@lgLabelFontHeightF     = 0.019     ; Change the font size
 resL4@xyDashPatterns   = (/1,1,0,0/)
  resL4@xyExplicitLegendLabels = (/"PI-w/oIDEMIX","LGM-w/oIDEMIX","PI", "LGM"/)
   resL4@lgLabelFontThicknessF = 1.2
  resL4@lgLineLabelFontThicknessF = 1.2
resL4@lgLineThicknesses = 1.3
resL4@gsnXRefLine = 1000
resL4@gsnXRefLineDashPattern = 1
resL4@gsnXRefLineColor = "black"
resL4@gsnXRefLineThicknessF = 2.0






resL4@tiYAxisString = ""  ; axis string
resL4@tiMainString    = ""   ; add titles
resL4@tmYLMode        = "Explicit"                     ; Force labels where we want them
;resL4@tmYLValues      = ispan(-45,15,5)

resL4@tmYLLabels = ispan(-50,10,10)

resL5 = True

 resL5@gsnDraw              = False                ; don't draw yet
  resL5@gsnFrame             = False                ; don't advance frame yet
  resL5@tiXAxisFontHeightF   = 0.022
   resL5@tiYAxisFontHeightF   = 0.022
   resL5@tmXBLabelFontHeightF  = 0.018
   resL5@tmYLLabelFontHeightF  = 0.018
   resL5@tiXAxisString = "Latitude"  ; axis string
;  resL5@xyMonoDashPattern = True
  resL5@gsnLeftStringParallelPosF = 0.02
  resL5@txFontHeightF      = 0.015
  resL5@pmLegendDisplayMode    = "NoCreate"            ; turn on legend
  resL5@lgPerimOn              = False               ; turn off box around
  resL5@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL5@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
   resL5@xyLineColors           = (/"dodgerblue","limegreen","blue4","green4"/)
;  resL5@xyLineColors           =  (/"coral1","dodgerblue","firebrick3","blue4"/)
 resL5@lgPerimOn              = False                 ; no box around
resL5@lgJustification        = "BottomRight"
  resL5@pmLegendOrthogonalPosF = -1.48     ; Move into plot
  resL5@pmLegendParallelPosF   = 1.1       ; Move to right
  resL5@pmLegendWidthF         = 0.3       ; Decrease width
  resL5@pmLegendHeightF        = 0.12       ; Decrease height
  resL5@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL5@lgLabelFontHeightF     = 0.019     ; Change the font size
 resL5@xyDashPatterns   = (/1,1,0,0/)
  resL5@xyExplicitLegendLabels = (/"PI-w/oIDEMIX","LGM-w/oIDEMIX","PI", "LGM"/)
   resL5@lgLabelFontThicknessF = 1.2
  resL5@lgLineLabelFontThicknessF = 1.2
resL5@lgLineThicknesses = 1.3
resL5@gsnXRefLine = 1000
resL5@gsnXRefLineDashPattern = 1
resL5@gsnXRefLineColor = "black"
resL5@gsnXRefLineThicknessF = 2.0

resL5@tiYAxisString = ""  ; axis string
resL5@tiMainString    = ""   ; add titles
resL5@tmYLMode        = "Explicit"
resL5@tmYLLabels = ispan(-10,60,10)




   plots(0) = gsn_csm_xy(wks,iw_pi&lat_aux_grid ,serie1,resL)
   plots(1) = gsn_csm_xy(wks,iw_pi_atl&lat_aux_grid,serie3,resL1)
   plots(2) = gsn_csm_xy(wks,iw_pi_atl&lat_aux_grid,serie5,resL2) 
   plots(3) = gsn_csm_xy(wks,iw_pi_salt&lat_aux_grid,serie2,resL3)
   plots(4) = gsn_csm_xy(wks,iw_pi_salt_atl&lat_aux_grid,serie4,resL4)
   plots(5) = gsn_csm_xy(wks,iw_pi_salt_atl&lat_aux_grid,serie6,resL5)   


printVarSummary(iw1_serie_atl)
print(max(iw1_serie_atl))
print(min(iw1_serie_atl))


printVarSummary(iw2_serie_atl)
print(max(iw2_serie_atl))
print(min(iw2_serie_atl))



printVarSummary(iw_pi_serie_atl)
print(max(iw_pi_serie_atl))
print(min(iw_pi_serie_atl))


printVarSummary(iw_lgm_serie_atl)
print(max(iw_lgm_serie_atl))
print(min(iw_lgm_serie_atl))


printVarSummary(iw_pi_serie_salt_atl)
print(max(iw1_serie_salt_atl))
print(min(iw1_serie_salt_atl))


printVarSummary(iw2_serie_salt_atl)
print(max(iw2_serie_salt_atl))
print(min(iw2_serie_salt_atl))



printVarSummary(iw_pi_serie_salt_atl)
print(max(iw_pi_serie_salt_atl))
print(min(iw_pi_serie_salt_atl))


printVarSummary(iw_lgm_serie_salt_atl)
print(max(iw_lgm_serie_salt_atl))
print(min(iw_lgm_serie_salt_atl))








resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
;resP@gsnPaperOrientation = "vertical"
;resP@gsnPanelLabelBar   = True
;;resP@pmLabelBarWidthF   = 0.07
;;resP@pmLabelBarHeightF = 0.23
;resP@pmLabelBarWidthF   = 0.07
;resP@pmLabelBarHeightF = 0.33

;;resP@pmLabelBarOrthogonalPosF = -0.005
;resP@pmLabelBarParallelPosF = 0.05
;resP@lbTitleExtentF = 0.18


;;resP@lbLabelOffsetF       = 0.01
;resP@lbTitleOffsetF = 0.02


;resP@lbOrientation        = "Vertical"
;resP@lbLabelPosition      = "Right"           ; label position
;resP@lbTitleOffsetF = 0.02


;resP@lbOrientation        = "Vertical"
;resP@lbLabelPosition      = "Right"           ; label position


;resP@lbTitleString = "log10[m^2/s]"
;resP@lbTitleFontHeightF = 0.010
;resP@lbTitlePosition	 = "Right"
;resP@lbTitleDirection     = "Across"                          ; letter angle
;resP@lbTitleAngleF        = 270.                               ; title angle
;resP@lbTitleDirection = "Down"
;resP@lbTitleOffsetF = 0.15

resP@gsnMaximize         = True                 ; make as large as possible   					
;resP@lbLabelFontHeightF= 0.011 


;;resP@gsnPanelXF = (/0.08,0.48,0.08,0.48/) 

;;resP@gsnPanelYF = (/0.75,0.75,0.55,0.55/)

;;resP@pmLabelBarOrthogonalPosF = -0.5
;resP@gsnPanelFigureStrings= (/"PI - w/oIDEMIX", "LGM - w/oIDEMIX", "PI - w/IDEMIX", "LGM - w/IDEMIX"/) 
;resP@gsnPanelFigureStringsFontHeightF = 0.010
;resP@cnLevelSelectionMode = "ExplicitLevels"
;;resP@cnLevels    = (/ -7, -6, -5, -4, -3, -2 /)   ; set levels
	
;;resP@gsnPanelBottom = 0.9
;;resP@gsnPanelTop = 0.92
;;resP@gsnPanelLeft = 0.02
;resP@gsnLeftString =  "Tidal Energy Dissipation"
;resP@gsnRightString           = "(W/m^2)"  

;;resP@gsnPanelXF = (/0.08,0.48,0.08,0.48/) 
;;resP@gsnPanelYF = (/0.75,0.75,0.55,0.55/)

	
;;resP@gsnPanelBottom = 0.9
;;resP@gsnPanelTop = 0.92
;;resP@gsnPanelLeft = 0.02
;;resP@gsnLeftString =  "Tidal Energy Dissipation"
;;resP@gsnRightString           = "(W/m^2)"  


;;txres3 = True
;;txres3@txFontHeightF = 0.011

;;gsn_text_ndc(wks,"Section near 60 degN",0.62,0.78,txres3)



txres4 = True
txres4@txFontHeightF = 0.011

;gsn_text_ndc(wks,"Near 30 degW in Atlantic Ocean",0.26,0.65,txres4)


txres5 = True
txres5@txFontHeightF = 0.011
  txres5@txAngleF      = 90.     ; rotate 90 degrees
;;  txres5@txJust        = "TopLeft"

;gsn_text_ndc(wks,"Depth (m)",0.02,0.53,txres5)

;gsn_text_ndc(wks,"Depth (m)",0.01,0.54,txres5)



gsn_panel(wks,plots,(/2,3/),resP)

 
; now we draw the whole thing and advance the frame.


  frame(wks)
  end

