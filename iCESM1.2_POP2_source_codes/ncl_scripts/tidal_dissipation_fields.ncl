

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


 pop_pi = addfile ("/mnt/d/tidal_energy_gx1v6_20090205_miss.nc", "r")

 pop_lgm_dinezio = addfile ("/mnt/d/diss_binary_tonetcdf2_miss.nc", "r")

 pop_lgm_wilmes = addfile ("/mnt/d/lgm_data/wave_heves08_bilinear_miss.nc", "r")

 pi_area = addfile ("/mnt/d/lgm_data/pi_tarea_m2.nc", "r")

 lgm_area = addfile ("/mnt/d/lgm_data/lgm_tarea_m2.nc", "r")


 pop_pi2 = addfile ("/mnt/d/tidal_energy_gx1v6_20090205.nc", "r")

 pop_lgm_dinezio2 = addfile ("/mnt/d/diss_binary_tonetcdf2.nc", "r")

 pop_lgm_wilmes2 = addfile ("/mnt/d/lgm_data/wave_heves08_bilinear.nc", "r")


wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"tidal_dissipation_field")


   ;cmap = read_colormap_file("temp_19lev")
   ;ncmap = dimsizes(cmap(:,0))
   
  cmap2 = read_colormap_file("MPL_RdYlBu")
  ncmap = dimsizes(cmap2(:,0))

res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
;  res@cnFillPalette        = "BlWhRe"         ; set color map
;  res@cnFillPalette        = read_colormap_file("MPL_rai")
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels

res@vpWidthF           = 0.15
res@vpHeightF          = 0.20
res@tmXBLabelsOn  = True
res@tmXBLabelFontHeightF = 0.006
res@tmYLLabelFontHeightF = 0.006
res@cnFillPalette        = cmap2(::-1,:)
res@mpLandFillColor = "grey56"
res@cnMissingValFillColor = "grey56"
res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "

pltres = True
mpres = True
pltres@PanelPlot  = True
res@FieldTitle = ""
res@UnitLabel  = ""
res@Footer    = False
plot = new(3,graphic)
;plots = new(2,graphic)
res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True
res@tmXBTickSpacingF     = 30.0
res@vpWidthF           = 0.3
res@vpHeightF          = 0.15
res@tmXBLabelFontHeightF = 0.007
res@tmYLLabelFontHeightF = 0.007

;res@gsnYAxisIrregular2Linear = False
;res@trYReverse           = False
res@lbLabelFontHeightF= 0.008
;res@pmLabelBarParallelPosF = -0.01
;res@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
res@lbTitleFontHeightF= 0.009
res@lbTitleOffsetF = 0.01

res@tmXBTickSpacingF = 60.0
res@tmYLTickSpacingF = 30.0
res@cnMissingValFillPattern = "SolidFill"
res@cnMissingValFillColor = "grey56"

   res@cnLevelSelectionMode    = "ManualLevels"  ; manually specify contour levels
   res@cnMinLevelValF            = -5             ; min level
   res@cnMaxLevelValF            =  -1             ; max level
   res@cnLevelSpacingF           =  0.5            ; contour interval


 lon1   = pop_pi->TLONG                   ; read in original lat/lon data
  lat1   = pop_pi->TLAT


res2 = res
;res2@gsnYAxisIrregular2Linear = True
;res2@trYReverse           = True
res2@tmYLLabelsOn  = False

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


pi         = pop_pi->wave_dissipation(1:nlat-1,:)
lgm_dinezio         = pop_lgm_dinezio->wave_dissipation(1:nlat-1,:)
lgm_wilmes         = pop_lgm_wilmes->wave_dissipation(1:nlat-1,:)


pi_area_m2 = pi_area->TAREA

lgm_area_m2 = lgm_area->TAREA


pi2 = log10(pi)
lgm_dinezio2 = log10(lgm_dinezio)
lgm_wilmes2 = log10(lgm_wilmes)

pi2@lat2d = lat
pi2@lon2d = lon

lgm_dinezio2@lat2d = lat
lgm_dinezio2@lon2d = lon

lgm_wilmes2@lat2d = lat
lgm_wilmes2@lon2d = lon



pi3         = pop_pi->wave_dissipation
lgm_dinezio3         = pop_lgm_dinezio->wave_dissipation
lgm_wilmes3         = pop_lgm_wilmes->wave_dissipation

  pi_sum = sum(pi3*pi_area_m2)*1d-12
  lgm_dinezio_sum = sum(lgm_dinezio3*lgm_area_m2)*1d-12
  lgm_wilmes_sum = sum(lgm_wilmes3*lgm_area_m2)*1d-12

  pi_avg = log10(pi_sum)
  lgm_dinezio_avg = log10(lgm_dinezio_sum)
  lgm_wilmes_avg = log10( lgm_wilmes_sum)

 ;pi_avg = sum(pi3)
 ;lgm_dinezio_avg = sum(lgm_dinezio3)
 ;lgm_wilmes_avg = sum(lgm_wilmes3)




print("pi_sum")
print(pi_sum)
print(pi_avg)

print("lgm_wilmes_sum")
print(lgm_wilmes_sum)
print(lgm_wilmes_avg)

print("lgm_dinezio_sum")
print(lgm_dinezio_sum)
print(lgm_dinezio_avg)


plot(0) = gsn_csm_contour_map(wks,pi2,res)
plot(1) = gsn_csm_contour_map(wks,lgm_dinezio2,res2)
plot(2) = gsn_csm_contour_map(wks,lgm_wilmes2,res2)


resP = True
;  resP@gsnMaximize         = True                 ; make as large as possible
;  resP@gsnPaperOrientation = "portrait"           ; force portrait
;  resP@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  resP@gsnMaximize        = False
  resP@gsnFrame         = False               ; save panel until both ready
  resP@gsnPaperOrientation = "landscape"
  resP@gsnPanelLabelBar   = True
;  resP@pmLabelBarWidthF   = 0.60
;  resP@pmLabelBarHeightF = 0.043
   resP@pmLabelBarWidthF   = 0.55
   resP@pmLabelBarHeightF = 0.041
   resP@pmLabelBarOrthogonalPosF = -0.042
   ;resP@pmLabelBarParallelPosF = -0.01
;   resP@lbTitleExtentF = 0.01
;resP@lbOrientation        = "Vertical"
;resP@lbLabelPosition      = "Right"           ; label position



resP@lbTitleString = "Tidal Energy Dissipation log10[W/m^2]"
;resP@lbTitleString = "Salt Flux [g/kg]"
;resP@lbTitleString = "Heat Flux [10^-6 kg/m2/s]"
;resP@lbTitleString = "Freshwater Flux [10^-6 kg/m2/s]"
;resP@gsnPaperHeight = 80
;resP@gsnPaperWidth = 80

resP@lbTitleFontHeightF = 0.010
;resP@lbTitlePosition     = "Top"
;resP@lbTitleDirection     = "Across"                          ; letter angle
;resP@lbTitleAngleF        = 270.                               ; title angle
;resP@lbTitleDirection = "Down"
;resP@lbTitleOffsetF = -0.003
;!resP@lbTitleFont   = "helvetica-bold"

;resP@gsnPanelYWhiteSpacePercent = 1
;resP@gsnPanelXWhiteSpacePercent = 1
resP@lbLabelFontHeightF= 0.009
resP@gsnPanelFigureStringsFontHeightF = 0.009
resP@gsnPanelFigureStrings= (/"(a)", "(b)", "(c)"/)

resP@gsnPanelTop = 0.95
resP@gsnPanelBottom = 0.61

;resP@gsnPanelFigureStringsFontHeightF = 0.012
resP@gsnPanelYWhiteSpacePercent = 0.3



;resP@gsnPanelXF = (/0.02,0.55/)
;resP@gsnPanelYF = (/0.95,0.95/)


;resP@lbLabelOffsetF = -0.2

resP1 = True

resP1@lbTitleString = "Difference in Salinity [g/kg]"
;resP1@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
resP1@gsnPanelTop = 0.64
resP1@gsnPanelBottom = 0.19 
;  resP1@gsnMaximize         = True                 ; make as large as possible
;  resP1@gsnPaperOrientation = "portrait"           ; force portrait
;  resP1@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  resP1@gsnMaximize        = False
  resP1@gsnFrame         = False               ; save panel until both ready
  resP1@gsnPaperOrientation = "landscape"
  resP1@gsnPanelLabelBar   = True
  resP1@pmLabelBarWidthF   = 0.60
  resP1@pmLabelBarHeightF = 0.043
;  resP1@pmLabelBarWidthF   = 0.07
;   resP1@pmLabelBarHeightF = 0.3
   resP1@pmLabelBarOrthogonalPosF = -0.04
   ;resP1@pmLabelBarParallelPosF = -0.01
;   resP1@lbTitleExtentF = 0.01
;resP1@lbOrientation        = "Vertical"
;resP1@lbLabelPosition      = "Right"           ; label position
;resP1@lbLabelOffsetF = -0.3


;resP1@gsnPaperHeightF = 80
;resP1@gsnPaperWidthF = 80

resP1@lbTitleFontHeightF = 0.011
;resP1@lbTitlePosition     = "Top"
;resP1@lbTitleDirection     = "Across"                          ; letter angle
;resP1@lbTitleAngleF        = 270.                               ; title angle
;resP1@lbTitleDirection = "Down"
;resP1@lbTitleOffsetF = -0.003
;!resP1@lbTitleFont   = "helvetica-bold"

resP1@gsnPanelYWhiteSpacePercent = 3
;resP1@gsnPanelXWhiteSpacePercent = 1

resP1@lbLabelFontHeightF= 0.010
resP1@gsnPanelFigureStrings= (/"(c)", "(d)"/)
resP1@gsnPanelFigureStringsFontHeightF = 0.011

;resP1@gsnPanelXF = (/0.02,0.55/)
;resP1@gsnPanelYF = (/0.5,0.5/)


gsn_panel(wks, plot, (/1,3/), resP)
;gsn_panel(wks, plots, (/1,2/), resP1)


txres2 = True
txres2@txFontHeightF = 0.012
txres2@txFont   = "helvetica-bold"

;gsn_text_ndc(wks,"PI (w/IDEMIX-w/oIDEMIX)",0.27,0.96,txres2)
;gsn_text_ndc(wks,"LGM (w/IDEMIX-w/oIDEMIX)",0.77,0.96,txres2)

;gsn_text_ndc(wks,"LGM - PI (w/oIDEMIX)",0.27,0.99,txres2)
;gsn_text_ndc(wks,"LGM - PI (w/IDEMIX)",0.77,0.99,txres2)

gsn_text_ndc(wks,"PI",0.18,0.91,txres2)
gsn_text_ndc(wks,"(LGM - DiNezio)",0.51,0.91,txres2)
gsn_text_ndc(wks,"(LGM - Wilmes)",0.84,0.91,txres2)

;  psres = True                                                               
;  maximize_output(wks,psres)  ; calls draw and frame for you 


 
frame(wks)
end

