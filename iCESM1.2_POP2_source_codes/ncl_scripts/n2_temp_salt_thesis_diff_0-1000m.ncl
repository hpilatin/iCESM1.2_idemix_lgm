


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


   ;pop_pi_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
   ;;pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")
   ; pop_pi_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")
    ;;pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/sel_buo_test_cmake_1750_1799_timmean.nc", "r")   


;;   pop_pi_tamas_n2 = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
;;   pop_lgm_tamas_n2 = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")


   pop_pi_n2_idemix = addfile ("/mnt/d/lgm_data/work4_work5_n2_pi_idemix_isotope_continued_5000_5099_timmean.nc", "r")
   pop_lgm_n2_idemix = addfile ("/mnt/d/lgm_data/work4_lgm_5849_5899_timmean.nc", "r")

  ;; pop_pi = addfile ("/mnt/d/lgm_data/temp_salt_n2_pi_4949_4999_timmean.nc", "r")
  ;; pop_lgm = addfile ("/mnt/d/lgm_data/temp_salt_n2_lgm_5849_5899_timmean.nc", "r")
 

   ;pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
  ;pop_pi = addfile ("/mnt/d/lgm_data/alldata_PI_ciso_j5_idemix_mvio_3920_3969_timmean.nc", "r")
  ;  pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_isotope_continued_4949_4999_timmean.nc", "r")
  ;  pop_pi = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")
  ;pop_lgm = addfile ("/mnt/d/lgm_data/surf_buo_dinezio_continued.pop.h.1810_1859_timmean.nc", "r")


  ; pop_pi_tamas_n2 = addfile ("/mnt/d/lgm_data2/monthly_pi_tamas_4150_4199_short_timmean.nc", "r")    
  ;; pop_lgm_tamas_n2 = addfile ("/mnt/d/lgm_data2/monthly_lgm_tamas_5450_5499_short_timmean.nc", "r")

   ;  pop_lgm_tamas_n2 = addfile ("/mnt/d/lgm_data2/select_test_cmake_1670_1719_timmean.nc", "r")


  
  

      
  pop_pi_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
   pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")
   pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")
   pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")

  ;pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_LGM21ka_frUsp.p875.001.pop.h.5400_5449_timmean.nc", "r") ;surf_buo_LGM21ka_frUsp.p875.001.pop.h.5450_5499.deflated_timmean.nc
   
   ;pop_lgm = addfile ("/mnt/d/lgm_data/select_dinezio_continued_1810_1859_timmean.nc", "r")

   ;pop_pi = addfile ("/mnt/d/lgm_data/sel_buo_pi_idemix_isotope_continued_5290_5299_timmean.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data/n2_temp_salt_dinezio_continued.pop.h.1850_1859_timmean.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data/sel_buo_idemix_lgm_idemix_isotope_continued_5900_5909_timmean.nc", "r")



    pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")


 ;!  pop_pi_tamas_temp_salt = addfile ("/mnt/d/lgm_data2/buo_piControl.p875.001.pop.h.4190_4199_timmean.nc", "r")
 ;!  pop_lgm_tamas_temp_salt = addfile ("/mnt/d/lgm_data2/buo_test_cmake_1860_1869_timmean.nc", "r")
 ;!  pop_pi_temp_salt = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4694_4703_timmean.nc", "r") 
 ;;pop_pi_temp_salt = addfile ("/mnt/d/lgm_data/n2_heat_salt_PI_ciso_j5_idemix_mvio_3960_3969_timmean.nc", "r")  
;;pop_pi_temp_salt = addfile ("/mnt/d/lgm_data/sel_buo_pi_idemix_isotope_continued_5290_5299_timmean.nc", "r")
 ;!  pop_lgm_temp_salt = addfile ("/mnt/d/lgm_data/sel_buo_idemix_lgm_idemix_isotope_continued_5900_5909_timmean.nc", "r")


   pop_lgm_tamas_temp_salt = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")
       
  pop_pi_tamas_temp_salt = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")

  pop_lgm_temp_salt = addfile ("/mnt/d/lgm_data/select_dinezio_continued_1810_1859_timmean.nc", "r")          
  
  pop_pi_temp_salt = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")

wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"n2_temp_salt_thesis_diff_0-5000m_noQ_lastnew")


  ; wks  = gsn_open_wks(wks_type,"n2_temp_salt_thesis_diff_column_pimvio_noQ")

;   cmap = read_colormap_file("cmocean_delta")
    cmap = read_colormap_file("cmocean_haline")
   ncmap = dimsizes(cmap(:,0))

cmap2 = read_colormap_file("MPL_RdYlBu")

res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = cmap2(::-1,:)         ; set color map
;  res@cnFillPalette        = read_colormap_file("MPL_rai")
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels


res@vpWidthF           = 0.15
res@vpHeightF          = 0.20
res@tmXBLabelsOn  = True
res@tmXBLabelFontHeightF = 0.005
res@tmYLLabelFontHeightF = 0.005
;res@cnFillPalette        = cmap
res@mpLandFillColor = "grey56"
res@cnMissingValFillColor = "grey56"
res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "

pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(6,graphic)
;plots = new(2,graphic)
;plots2 = new(2,graphic)
res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True
;res@tmXBTickSpacingF     = 15.0
res@vpWidthF           = 0.115
res@vpHeightF          = 0.08
res@tmXBLabelFontHeightF = 0.004

res@tmYLLabelFontHeightF = 0.0038


;res@gsnYAxisIrregular2Linear = False
;res@trYReverse           = False
res@lbLabelFontHeightF= 0.0035
;res@pmLabelBarParallelPosF = -0.01
;res@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
res@lbTitleString = "[10^-7 1/s^2]"
res@lbTitleFontHeightF= 0.004
;res@lbTitleOffsetF = 0.01

res@tmXBTickSpacingF = 30.0 ; North Atlantic
;res@tmXBTickSpacingF = 10.0
res@tmYLTickSpacingF = 30.0
res@cnMissingValFillPattern = "SolidFill"
res@cnMissingValFillColor = "grey56"

;res@pmLabelBarOrthogonalPosF = .03             ; position label bar

res@pmLabelBarOrthogonalPosF = .04

    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            =  -2             ; min level
    res@cnMaxLevelValF            =  2             ; max level
    res@cnLevelSpacingF           =  0.25            ; contour interval

; res@cnLevelSelectionMode = "ManualLevels"   ; set explicit contour levels
; res@cnMinLevelValF            =  -0.6             ; min level
; res@cnMaxLevelValF            =  0.6            ; max level
; res@cnLevelSpacingF           =  0.1            ; contour interval 


;res@cnLevels    = (/-1,-0.98,-0.96,-0.94,-0.92,-0.9,-0.88,-0.86,-0.84,-0.82,-0.8,-0.78,-0.76,-0.74,-0.72,-0.7,-0.68,-0.66,-0.64,-0.62,-0.6,-0.58,-0.56,-0.54,-0.52,-0.5,-0.48,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1/)   ; set levels

;!res@cnLevels    = (/-0.6,-0.58,-0.56,-0.54,-0.52,-0.5,-0.48,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6/)   ; set levels ; North Atlantic

;res@cnLevels    = (/-0.6,-0.58,-0.56,-0.54,-0.52,-0.5,-0.48,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1/)   ; set levels

;res@cnLevels    = (/-3,-2.8,-2.6,-2.4,-2.2,-2,-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3/)  

;!res@cnLevels    = (/-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3/) ; Southern Ocean


;res@cnLevels    = (/-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3/)   ; set levels ; Difference

;!!res@cnLevels    = (/-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2/)   ; set levels ; Difference original scale


;res@cnLevels    = (/-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2/)   ; set levels ; Difference original scale

;res@cnLevels    = (/-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2/)


res@gsnYAxisIrregular2Linear = True
res@trYReverse           = True       ; reverses y-axis

; res@lbOrientation        = "vertical"  ; vertical label bar


res1                      = True                 ; plot mods desired
res1@mpLandFillColor = "grey56"
res1@cnMissingValFillColor = "grey56"

  res1@gsnDraw              = False                ; don't draw yet
  res1@gsnFrame             = False                ; don't advance frame yet
  res1@cnFillOn             = True                 ; turn on color fill
   cmap2 = read_colormap_file("MPL_RdYlBu")
  res1@cnFillPalette        = cmap2(::-1,:)         ; set color map
;  res1@cnFillPalette        = read_colormap_file("MPL_RdYlBu")
  res1@cnLinesOn            = False                ; turn off contour lines
  res1@cnLineLabelsOn       = False                ; turn off contour line labels
res1@vpWidthF           = 0.2
res1@vpHeightF          = 0.12
res1@tmXBLabelsOn  = True
res1@tmXBLabelFontHeightF = 0.005
res1@tmYLLabelFontHeightF = 0.0048

;res1@mpLandFillColor = "grey"

res1@lbLabelBarOn        = False
res1@cnInfoLabelOn       = False
res1@description = "  "

 res1@gsnCenterString = "  "
res1@tiMainString  = " "             ; add title

res1@cnMissingValFillPattern = "SolidFill"
res1@tmYLLabelsOn  = True
res1@tmXBLabelsOn  = True
res1@tmXBLabelFontHeightF = 0.005
res1@tmYLLabelFontHeightF = 0.0045


res1@vpWidthF           = 0.15
res1@vpHeightF          = 0.10
res1@lbTitleString = "Salt [g/kg]"
res1@lbTitleFontHeightF= 0.004
;res1@lbTitleOffsetF = 0.01

; res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res1@cnMinLevelValF            = -1.             ; min level
;    res1@cnMaxLevelValF            =  1.             ; max level
;    res1@cnLevelSpacingF           =  0.1            ; contour interval


res1@tmXBTickSpacingF = 30.0 ; North Atlantic
;res1@tmXBTickSpacingF = 10.0
res1@tmYLTickSpacingF = 30.0

res1@gsnYAxisIrregular2Linear = True
res1@trYReverse           = True       ; reverses y-axis
res1@tmYLLabelsOn  = False
res1@cnMissingValFillPattern = "SolidFill"

res1@tmXBTickSpacingF     = 15.0
res1@vpWidthF           = 0.15
res1@vpHeightF          = 0.1

; res1@lbOrientation        = "vertical"  ; vertical label bar

res1@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
res1@cnLevels    = (/-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5/)   ; set levels ; Difference

res2 = True
res2@gsnDraw              = False                ; don't draw yet
res2@gsnFrame             = False                ; don't advance frame yet

res2@tmXBTickSpacingF     = 30.0 ; North Atlantic
;res2@tmXBTickSpacingF = 10.0 ; SO

res2@vpWidthF           = 0.5
res2@vpHeightF          = 0.1
res2@tmXBLabelFontHeightF = 0.005
res2@tmYLLabelFontHeightF = 0.004

res2@tmXBLabelsOn  = True
;res2@tmXBLabelFontHeightF = 0.005
;res2@tmYLLabelFontHeightF = 0.005
res2@cnFillPalette        = cmap
res2@mpLandFillColor = "grey56"
res2@cnMissingValFillColor = "grey56"
res2@lbLabelBarOn        = False
res2@cnInfoLabelOn       = False

res2@cnLinesOn            = False                ; turn off contour lines
res2@cnLinesOn            = False                ; turn off contour lines
res2@cnLineLabelsOn       = False                ; turn off contour line labels
res2@tmXBLabelsOn  = True
res2@tmYLabelsOn  = False
res2@lbLabelFontHeightF= 0.008
res2@mpLandFillColor = "grey56"
res2@cnMissingValFillColor = "grey56"
res2@cnFillPalette        = cmap
res2@mpLandFillColor = "grey56"
res2@cnMissingValFillColor = "grey56"
res2@lbLabelBarOn        = False
res2@cnInfoLabelOn       = False
res2@cnLinesOn            = False                ; turn off contour lines
res2@cnLineLabelsOn       = False                ; turn off contour line labels
res2@tmXBLabelsOn  = True
res2@tmYLabelsOn  = False
res2@lbLabelFontHeightF= 0.008
res2@mpLandFillColor = "grey56"
res2@cnMissingValFillColor = "grey56"
res2@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
res2@lbTitleFontHeightF= 0.004
res2@lbTitleOffsetF = 0.01
res2@cnFillOn             = True
res2@cnFillPalette        = cmap(::-1,:)
;res2@tmXBTickSpacingF = 60.0
;res2@tmYLTickSpacingF = 30.0
; res2@cnMinLevelValF            = -3             ; min level
; res2@cnMaxLevelValF            =  3             ; max level
; res2@cnLevelSpacingF           =  0.2            ; contour interval
res2@tmXBLabelsOn  = True

res2@gsnYAxisIrregular2Linear = True
res2@trYReverse           = True       ; reverses y-axis
res2@lbLabelFontHeightF= 0.003

; res2@lbOrientation        = "vertical"  ; vertical label bar

res4 = True
res4@gsnDraw              = False                ; don't draw yet
res4@gsnFrame             = False                ; don't advance frame yet
res4@tmXBTickSpacingF     = 5.0
res4@vpWidthF           = 0.115
res4@vpHeightF          = 0.08
res4@tmXBLabelFontHeightF = 0.004
res4@tmYLLabelFontHeightF = 0.0038
res4@tmXBLabelsOn  = True

res4@mpLandFillColor = "grey56"
res4@cnMissingValFillColor = "grey56"
res4@lbLabelBarOn        = False
res4@cnInfoLabelOn       = False
res4@cnLinesOn            = False                ; turn off contour lines
res4@cnLineLabelsOn       = False                ; turn off contour line labels
res4@tmXBLabelsOn  = True
res4@tmYLabelsOn  = False
res4@lbTitleString = "[g/kg]"
res4@lbTitleFontHeightF= 0.004
;res4@lbTitleOffsetF = 0.01
res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
res4@mpLandFillColor = "grey56"
res4@cnMissingValFillColor = "grey56"
res4@cnFillOn             = True
res4@cnFillPalette        = cmap2(::-1,:)

    res4@cnMinLevelValF            = -0.3             ; min level
    res4@cnMaxLevelValF            =  0.3             ; max level
    res4@cnLevelSpacingF           =  0.04            ; contour interval


;res4@cnFillPalette        = cmap(::-1,:)

res4@tmXBTickSpacingF = 30.0 ; North Atlantic
;res4@tmXBTickSpacingF = 10.0  ; SO
res4@tmYLTickSpacingF = 30.0

res4@mpLandFillColor = "grey56"
res4@cnMissingValFillColor = "grey56"
res4@tmXBLabelsOn  = True
res4@tmXBLabelsOn  = True
res4@tmYLLabelsOn  = False
res4@gsnYAxisIrregular2Linear = True
res4@trYReverse           = True       ; reverses y-axis

res4@lbLabelFontHeightF= 0.0035
;res4@pmLabelBarOrthogonalPosF = .03             ; position label bar

res4@pmLabelBarOrthogonalPosF = .04
; res4@lbOrientation        = "vertical"  ; vertical label bar

;res4@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
;res4@cnLevels    = (/-0.5,-0.48,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5/)   ; set levels ; Difference

;res4@cnLevels    = (/-0.3,-0.26,-0.22,-0.18,-0.14,-0.1,-0.06,-0.02,0,0.02,0.06,0.1,0.14,0.18,0.22,0.26,0.3/)   ; set levels ; Difference


res3 = res4
res3@lbTitleString = "[degC]"

;res3@lbTitleString = "In-Situ Density [kg/m^3]"
delete(res3@cnLevelSelectionMode)
delete(res3@cnLevels)
;res3@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
;res3@cnLevels    = (/-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5/)   ; set levels ; Difference

    res3@cnLevelSelectionMode      = "ManualLevels"
    res3@cnMinLevelValF            = -2             ; min level
    res3@cnMaxLevelValF            =  2             ; max level
    res3@cnLevelSpacingF           =  0.2            ; contour interval

  lon1   = pop_pi->TLONG                   ; read in original lat/lon data
  lat1   = pop_pi->TLAT


region_pi  = pop_pi_region->REGION_MASK
region_lgm  = pop_lgm_region->REGION_MASK

region_atlantic_pi= new((/384,320/),float)
region_atlantic_lgm= new((/384,320/),float)

;;;;############################ Assumptions########################

;;; beta= 0.801, alfa= 781*10^-7 (1/K) (these coefficients are taken from Atmosphere Ocean Dynamics (Adrian E. Gill) Table A3.1

;;; Cp= 4000 J/kg*K, sea water density (swd)= 1022.4 kg/m^3
 
;;;################################################################


;### ATLANTIC OCEAN;###




region_atlantic_pi= new((/384,320/),float)
region_atlantic_lgm= new((/384,320/),float)
region_po_pi= new((/384,320/),float)
region_po_lgm= new((/384,320/),float)
region_so_pi= new((/384,320/),float)
region_so_lgm= new((/384,320/),float)



region_atlantic_pi        = where(region_pi.eq.6,1,9.969209968386869e+36)

region_atlantic_lgm        = where(region_lgm.eq.6,1,9.969209968386869e+36)

;### SO OCEAN;###

region_so_pi        = where(region_pi.ne.1,9.969209968386869e+36,1)

region_so_lgm        = where(region_lgm.ne.1,9.969209968386869e+36,1)



;### PO OCEAN;###

region_po_pi        =  where(region_pi.ne.2,9.969209968386869e+36,1)

region_po_lgm        =  where(region_lgm.ne.2,9.969209968386869e+36,1)


pi_temp_tamas         = pop_pi_tamas_temp_salt->TEMP(0,:,:,:) 
pi_salt_tamas         = pop_pi_tamas_temp_salt->SALT(0,:,:,:) 

pi_stability_tamas         = pop_pi_tamas->Q(:,:,:,:)
;pi_pd_tamas         = pop_pi_tamas->PD(:,:,:,:)*1000 ; from g/cm3 to kg/m3
;pi_sta_tamas_mean    = abs(dim_avg_n(pi_stability_tamas,0)*1d5) ;kg/m4
pi_sta_tamas_mean    = dim_avg_n(pi_stability_tamas,0)*1d5 ;kg/m4
;pi_pd_tamas_mean    = dim_avg_n(pi_pd_tamas,0)
;pi_n2_tamas    = doubletofloat((9.80665*pi_sta_tamas_mean*10000)/pi_pd_tamas_mean)

;pi_n2_tamas    = (9.80665*pi_sta_tamas_mean*10000)/pi_pd_tamas_mean

;pi_n2_tamas         = pop_pi_n2_idemix->WORK5(0,:,:,:)


pi_rho_tamas         = pop_pi_tamas->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

z    = pop_pi_tamas->z_t

;--- Dimensions
  nlev = dimsizes(z)
  nlat = dimsizes(pi_rho_tamas(1,:,0))
  nlon = dimsizes(pi_rho_tamas(1,0,:))

depth           = pop_pi->z_t
depth_m2          = depth*1.E-2       ; convert from cm to m

;pi_rho_tamas_depth = new((/60,384,320/),float)
;do i=0,59
;pi_rho_tamas_depth(i,:,:) = pi_rho_tamas(i,:,:)/depth(i)
;end do

water_ref_density = 1025 ;kg/m3

z2          = z*1.E-2


  ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

  ;pi_rho_depth_tamas = new((/nlev-1, nlat, nlon/), typeof(pi_rho_tamas), "float")

   pi_rho_depth_tamas = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z2(k+1) - z2(k)
    pi_rho_depth_tamas(k,:,:) = (pi_rho_tamas(k+1,:,:) - pi_rho_tamas(k,:,:)) / dz
  end do


;water_ref_density = 1025 kg/m3

pi_n2_tamas    = -(((9.80665*pi_rho_depth_tamas)/water_ref_density)/0.0000001)


print(min(pi_n2_tamas))
print(max(pi_n2_tamas))






lgm_temp_tamas         = pop_lgm_tamas_temp_salt->TEMP(0,:,:,:) 
lgm_salt_tamas         = pop_lgm_tamas_temp_salt->SALT(0,:,:,:) 
;lgm_stability_tamas         = pop_lgm_tamas->Q(:,:,:,:)
;lgm_pd_tamas         = pop_lgm_tamas->PD(:,:,:,:)*1000 ; from g/cm3 to kg/m3
;lgm_sta_tamas_mean    = abs(dim_avg_n(lgm_stability_tamas,0)*1d5) ;kg/m4
;lgm_sta_tamas_mean    = dim_avg_n(lgm_stability_tamas,0)*1d5 ;kg/m4
;lgm_pd_tamas_mean    = dim_avg_n(lgm_pd_tamas,0)
;lgm_n2_tamas    = doubletofloat((9.80665*lgm_sta_tamas_mean*10000)/pi_pd_tamas_mean)
;lgm_n2_tamas    = (-1*9.80665*lgm_sta_tamas_mean*10000)/pi_pd_tamas_mean

;lgm_n2_tamas         = pop_lgm->WORK5(0,:,:,:)


lgm_rho_tamas         = pop_lgm_tamas->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3


  ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

;  lgm_rho_tamas_depth = new((/nlev-1, nlat, nlon/), typeof(lgm_rho_tamas), "float")

   lgm_rho_tamas_depth = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z2(k+1) - z2(k)
    lgm_rho_tamas_depth(k,:,:) = (lgm_rho_tamas(k+1,:,:) - lgm_rho_tamas(k,:,:)) / dz
  end do

;water_ref_density = 1025 kg/m3

lgm_n2_tamas    = -(((9.80665*lgm_rho_tamas_depth)/water_ref_density)/0.0000001)



pi_temp         = pop_pi_temp_salt->TEMP(0,:,:,:) 
pi_salt         = pop_pi_temp_salt->SALT(0,:,:,:) 

pi_stability         = pop_pi->Q(:,:,:,:)
pi_pd         = pop_pi->PD(:,:,:,:)*1000 ; from g/cm3 to kg/m3
;pi_sta_mean    = abs(dim_avg_n(pi_stability,0)*1d5) ;kg/m4
pi_sta_mean    = dim_avg_n(pi_stability,0)*1d5 ;kg/m4
pi_pd_mean    = dim_avg_n(pi_pd,0)
;pi_n2    = doubletofloat((9.80665*pi_sta_mean*10000)/pi_pd_mean)
;pi_n2    = (9.80665*pi_sta_mean*10000)/pi_pd_mean

;pi_n2         = pop_pi_n2_idemix->WORK4(0,:,:,:)

pi_rho         = pop_pi->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

;pi_rho_depth = new((/60,384,320/),float)
;do i=0,59
;pi_rho_depth(i,:,:) = pi_rho(i,:,:)/depth(i)
;end do

;water_ref_density = 1000 kg/m3


  ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

  ;pi_rho_depth = new((/nlev-1, nlat, nlon/), typeof(pi_rho), "float")

   pi_rho_depth = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z2(k+1) - z2(k)
    pi_rho_depth(k,:,:) = (pi_rho(k+1,:,:) - pi_rho(k,:,:)) / dz
  end do


;water_ref_density = 1025 kg/m3

pi_n2    = -(((9.80665*pi_rho_depth)/water_ref_density)/0.0000001)


lgm_temp       = pop_lgm_temp_salt->TEMP(0,:,:,:) 
lgm_salt         = pop_lgm_temp_salt->SALT(0,:,:,:) 

lgm_stability         = pop_lgm->Q(:,:,:,:)
;lgm_pd         = pop_lgm->PD(:,:,:,:)*1000 ; from g/cm3 to kg/m3
;lgm_sta_mean    = abs(dim_avg_n(lgm_stability,0)*1d5) ;kg/m4
lgm_sta_mean    = dim_avg_n(lgm_stability,0)*1d5 ;kg/m4
;lgm_pd_mean    = dim_avg_n(lgm_pd,0)
;lgm_n2    = doubletofloat((9.80665*lgm_sta_mean*10000)/lgm_pd_mean)
;lgm_n2    = (9.80665*lgm_sta_mean*10000)/lgm_pd_mean


lgm_rho         = pop_lgm->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

  ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

  ;lgm_rho_depth = new((/nlev-1, nlat, nlon/), typeof(lgm_rho), "float")

   lgm_rho_depth = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z2(k+1) - z2(k)
    lgm_rho_depth(k,:,:) = (lgm_rho(k+1,:,:) - lgm_rho(k,:,:)) / dz
  end do


;water_ref_density = 1025 kg/m3

lgm_n2    = -(((9.80665*lgm_rho_depth)/water_ref_density)/0.0000001)


printVarSummary(pi_n2_tamas)
print(min(pi_n2_tamas))
print(max(pi_n2_tamas))

printVarSummary(lgm_n2_tamas)
print(min(lgm_n2_tamas))
print(max(lgm_n2_tamas))


printVarSummary(pi_n2)
print(min(pi_n2))
print(max(pi_n2))

printVarSummary(lgm_n2)
print(min(lgm_n2))
print(max(lgm_n2))


pi_n2_v2         = pop_pi_n2_idemix->WORK4(0,:,:,:)
printVarSummary(pi_n2_v2)

print(min(pi_n2_v2))

print(max(pi_n2_v2))

lgm_n2_v2         = pop_lgm_n2_idemix->WORK4(0,:,:,:)
printVarSummary(lgm_n2_v2)
print(min(lgm_n2_v2))
print(max(lgm_n2_v2))



pi_temp_tamas_na = mask(pi_temp_tamas,region_atlantic_pi,1)

lgm_temp_tamas_na = mask(lgm_temp_tamas,region_atlantic_lgm,1)

pi_temp_na = mask(pi_temp,region_atlantic_pi,1)

lgm_temp_na = mask(lgm_temp,region_atlantic_lgm,1)


pi_temp_tamas_so = mask(pi_temp_tamas,region_so_pi,1)

lgm_temp_tamas_so = mask(lgm_temp_tamas,region_so_lgm,1)

pi_temp_so = mask(pi_temp,region_so_pi,1)

lgm_temp_so = mask(lgm_temp,region_so_lgm,1)


pi_salt_tamas_na = mask(pi_salt_tamas,region_atlantic_pi,1)

lgm_salt_tamas_na = mask(lgm_salt_tamas,region_atlantic_lgm,1)

pi_salt_na = mask(pi_salt,region_atlantic_pi,1)

lgm_salt_na = mask(lgm_salt,region_atlantic_lgm,1)


pi_salt_tamas_so = mask(pi_salt_tamas,region_so_pi,1)

lgm_salt_tamas_so = mask(lgm_salt_tamas,region_so_lgm,1)

pi_salt_so = mask(pi_salt,region_so_pi,1)

lgm_salt_so = mask(lgm_salt,region_so_lgm,1)


pi_n2_tamas_na = mask(pi_n2_tamas,region_atlantic_pi,1)

lgm_n2_tamas_na = mask(lgm_n2_tamas,region_atlantic_lgm,1)

pi_n2_na = mask(pi_n2,region_atlantic_pi,1)

lgm_n2_na = mask(lgm_n2,region_atlantic_lgm,1)

printVarSummary(lgm_n2_na)
print(min(lgm_n2_na))
print(max(lgm_n2_na))

pi_n2_tamas_so = mask(pi_n2_tamas,region_so_pi,1)

lgm_n2_tamas_so = mask(lgm_n2_tamas,region_so_lgm,1)

pi_n2_so = mask(pi_n2,region_so_pi,1)

lgm_n2_so = mask(lgm_n2,region_so_lgm,1)


pi_rho_tamas_na = mask(pi_rho_tamas,region_atlantic_pi,1)

lgm_rho_tamas_na = mask(lgm_rho_tamas,region_atlantic_lgm,1)

pi_rho_na = mask(pi_rho,region_atlantic_pi,1)

lgm_rho_na = mask(lgm_rho,region_atlantic_lgm,1)



pi_temp_tamas_mean_gb    = dim_avg_n(pi_temp_tamas,2)
pi_salt_tamas_mean_gb    = dim_avg_n(pi_salt_tamas,2)
pi_n2_tamas_mean_gb    = dim_avg_n(pi_n2_tamas,2)

lgm_temp_tamas_mean_gb    = dim_avg_n(lgm_temp_tamas,2)
lgm_salt_tamas_mean_gb    = dim_avg_n(lgm_salt_tamas,2)
lgm_n2_tamas_mean_gb    = dim_avg_n(lgm_n2_tamas,2)

pi_temp_mean_gb    = dim_avg_n(pi_temp,2)
pi_salt_mean_gb    = dim_avg_n(pi_salt,2)
pi_n2_mean_gb    = dim_avg_n(pi_n2,2)

lgm_temp_mean_gb    = dim_avg_n(lgm_temp,2)
lgm_salt_mean_gb    = dim_avg_n(lgm_salt,2)
lgm_n2_mean_gb    = dim_avg_n(lgm_n2,2)



pi_rho_tamas_mean    = dim_avg_n(pi_rho_tamas_na,2)
lgm_rho_tamas_mean    = dim_avg_n(lgm_rho_tamas_na,2)
pi_rho_mean    = dim_avg_n(pi_rho_na,2)
lgm_rho_mean    = dim_avg_n(lgm_rho_na,2)


pi_temp_tamas_mean    = dim_avg_n(pi_temp_tamas_na,2) 
pi_salt_tamas_mean    = dim_avg_n(pi_salt_tamas_na,2) 
pi_n2_tamas_mean    = dim_avg_n(pi_n2_tamas_na,2) 

lgm_temp_tamas_mean    = dim_avg_n(lgm_temp_tamas_na,2) 
lgm_salt_tamas_mean    = dim_avg_n(lgm_salt_tamas_na,2) 
lgm_n2_tamas_mean    = dim_avg_n(lgm_n2_tamas_na,2)

pi_temp_mean    = dim_avg_n(pi_temp_na,2) 
pi_salt_mean    = dim_avg_n(pi_salt_na,2) 
pi_n2_mean    = dim_avg_n(pi_n2_na,2) 

lgm_temp_mean    = dim_avg_n(lgm_temp_na,2) 
lgm_salt_mean    = dim_avg_n(lgm_salt_na,2) 
lgm_n2_mean    = dim_avg_n(lgm_n2_na,2) 


pi_temp_tamas_mean_so    = dim_avg_n(pi_temp_tamas_so,2)
pi_salt_tamas_mean_so    = dim_avg_n(pi_salt_tamas_so,2)
pi_n2_tamas_mean_so    = dim_avg_n(pi_n2_tamas_so,2)

lgm_temp_tamas_mean_so    = dim_avg_n(lgm_temp_tamas_so,2)
lgm_salt_tamas_mean_so    = dim_avg_n(lgm_salt_tamas_so,2)
lgm_n2_tamas_mean_so    = dim_avg_n(lgm_n2_tamas_so,2)

pi_temp_mean_so    = dim_avg_n(pi_temp_so,2)
pi_salt_mean_so    = dim_avg_n(pi_salt_so,2)
pi_n2_mean_so    = dim_avg_n(pi_n2_so,2)

lgm_temp_mean_so    = dim_avg_n(lgm_temp_so,2)
lgm_salt_mean_so    = dim_avg_n(lgm_salt_so,2)
lgm_n2_mean_so    = dim_avg_n(lgm_n2_so,2)


pi_temp_mean_so2 = pi_temp_mean_so(:,0:110)
lgm_temp_mean_so2 = lgm_temp_mean_so(:,0:110)
pi_temp_tamas_mean_so2 = pi_temp_tamas_mean_so(:,0:110)
lgm_temp_tamas_mean_so2 = lgm_temp_tamas_mean_so(:,0:110)

pi_salt_mean_so2 = pi_salt_mean_so(:,0:110)
lgm_salt_mean_so2 = lgm_salt_mean_so(:,0:110)
pi_salt_tamas_mean_so2 = pi_salt_tamas_mean_so(:,0:110)
lgm_salt_tamas_mean_so2 = lgm_salt_tamas_mean_so(:,0:110)

pi_n2_mean_so2 = pi_n2_mean_so(:,0:110)
lgm_n2_mean_so2 = lgm_n2_mean_so(:,0:110)
pi_n2_tamas_mean_so2 = pi_n2_tamas_mean_so(:,0:110)
lgm_n2_tamas_mean_so2 = lgm_n2_tamas_mean_so(:,0:110)



;sub1 = (pi_n2_mean(0:41,:)) - (pi_n2_tamas_mean(0:41,:)) ;0-1000 m 
sub1 = (pi_n2_mean(0:53,:)) - (pi_n2_tamas_mean(0:53,:)) ;0-3000 m

;sub1 = abs(pi_n2_mean(:,:)) - abs(pi_n2_tamas_mean(:,:))

;sub1 = (pi_n2_mean_gb(:,:)) - (pi_n2_tamas_mean_gb(:,:)) ;global

;sub1 = abs(pi_n2_mean_so2(0:41,:)) - abs(pi_n2_tamas_mean_so2(0:41,:))
printVarSummary(sub1)
print(avg(sub1))
print(max(sub1))
print(min(sub1))


;!sub2 = (pi_rho_mean(0:41,:)) - (pi_rho_tamas_mean(0:41,:))

;sub2 = (pi_temp_mean(0:41,:)) - (pi_temp_tamas_mean(0:41,:)) ;0-1000 m 

sub2 = (pi_temp_mean(0:53,:)) - (pi_temp_tamas_mean(0:53,:))

;sub2 = abs(pi_temp_mean(:,:)) - abs(pi_temp_tamas_mean(:,:))

;sub2 = (pi_temp_mean_gb(:,:)) - (pi_temp_tamas_mean_gb(:,:)) ;global

;sub2 = abs(pi_temp_mean_so2(0:41,:)) - abs(pi_temp_tamas_mean_so2(0:41,:))

printVarSummary(sub2)
print(avg(sub2))
print(max(sub2))
print(min(sub2))

;sub3 = (pi_salt_mean(0:41,:)) - (pi_salt_tamas_mean(0:41,:)) ;0-1000 m

sub3 = (pi_salt_mean(0:53,:)) - (pi_salt_tamas_mean(0:53,:))

;sub3 = (pi_salt_mean(:,:)) - (pi_salt_tamas_mean(:,:))

;sub3 = (pi_salt_mean_gb(:,:)) - (pi_salt_tamas_mean_gb(:,:)) ; global


;sub3 = abs(pi_salt_mean_so2(0:41,:)) - abs(pi_salt_tamas_mean_so2(0:41,:))

printVarSummary(sub3)
print(avg(sub3))
print(max(sub3))
print(min(sub3))


;sub4 = (lgm_n2_mean(0:41,:)) - (lgm_n2_tamas_mean(0:41,:)) ; 0-1000 m

sub4 = (lgm_n2_mean(0:53,:)) - (lgm_n2_tamas_mean(0:53,:))

;sub4 = abs(lgm_n2_mean(:,:)) - abs(lgm_n2_tamas_mean(:,:))

;sub4 = (lgm_n2_mean_gb(:,:)) - (lgm_n2_tamas_mean_gb(:,:)) ; global

;sub4 = abs(lgm_n2_mean_so2(0:41,:)) - abs(lgm_n2_tamas_mean_so2(0:41,:))

printVarSummary(sub4)
print(avg(sub4))
print(max(sub4))
print(min(sub4))


;!sub5 = (lgm_rho_mean(0:41,:)) - (lgm_rho_tamas_mean(0:41,:)) ; 0-1000 m

;sub5 = (lgm_temp_mean(0:41,:)) - (lgm_temp_tamas_mean(0:41,:)) ; 0-1000 m

sub5 = (lgm_temp_mean(0:53,:)) - (lgm_temp_tamas_mean(0:53,:))

;sub5 = abs(lgm_temp_mean(:,:)) - abs(lgm_temp_tamas_mean(:,:))

;sub5 = (lgm_temp_mean_gb(:,:)) - (lgm_temp_tamas_mean_gb(:,:))  ; global

;sub5 = abs(lgm_temp_mean_so2(0:41,:)) - abs(lgm_temp_tamas_mean_so2(0:41,:))
printVarSummary(sub5)
print(avg(sub5))
print(max(sub5))
print(min(sub5))

;sub6 = (lgm_salt_mean(0:41,:)) - (lgm_salt_tamas_mean(0:41,:)) ; 0-1000 m

sub6 = (lgm_salt_mean(0:53,:)) - (lgm_salt_tamas_mean(0:53,:))

;sub6 = (lgm_salt_mean(:,:)) - (lgm_salt_tamas_mean(:,:))

;sub6 = (lgm_salt_mean_gb(:,:)) - (lgm_salt_tamas_mean_gb(:,:))

;sub6 = abs(lgm_salt_mean_so2(0:41,:)) - abs(lgm_salt_tamas_mean_so2(0:41,:))
printVarSummary(sub6)
print(avg(sub6))
print(max(sub6))
print(min(sub6))

nlt = 9
nlt2 = 224
nlt3 = 224

lat_a    = pop_lgm->TLAT(:,nlt)
lat_atlantic        = lat_a
lat_p    = pop_lgm->TLAT(:,nlt2)
lat_pacific        = lat_p

lat_s    = pop_lgm->TLAT(0:110,0)
lat_so        = lat_s

lat_g    = pop_lgm->TLAT(:,0)
lat_gb        = lat_g



z_t           = pop_lgm->z_t
z_t           = z_t*1.E-2       ; convert from cm to m
;z_t@long_name = "Depth (T grid)"; overwrite long name
z_t@long_name = ""; overwrite long name
z_t@units     = "m"             ; fix units attribute
z_t@tiYAxisPosition = "top"
z_t@tiYAxisFontHeightF = 0.018

;pi_temp_tamas_mean!0    = "z_t"
;pi_temp_tamas_mean&z_t  =  z_t
;pi_temp_tamas_mean!1    = "lat_t"
;pi_temp_tamas_mean&lat_t=  lat_atlantic

;pi_temp_tamas_mean_so2!0    = "z_t"
;pi_temp_tamas_mean_so2&z_t  =  z_t
;pi_temp_tamas_mean_so2!1    = "lat_t"
;pi_temp_tamas_mean_so2&lat_t=  lat_so

;pi_temp_tamas_mean_gb!0    = "z_t"
;pi_temp_tamas_mean_gb&z_t  =  z_t
;pi_temp_tamas_mean_gb!1    = "lat_t"
;pi_temp_tamas_mean_gb&lat_t=  lat_gb

sub2!0    = "z_t"
;sub2&z_t  =  z_t
;sub2&z_t  =  z_t(0:41)
sub2&z_t  =  z_t(0:53)
sub2!1    = "lat_t"
sub2&lat_t=  lat_atlantic

;sub2&lat_t=  lat_so


;#####################################

;pi_salt_tamas_mean!0    = "z_t"
;pi_salt_tamas_mean&z_t  =  z_t
;pi_salt_tamas_mean!1    = "lat_t"
;pi_salt_tamas_mean&lat_t=  lat_atlantic

;pi_salt_tamas_mean_so2!0    = "z_t"
;pi_salt_tamas_mean_so2&z_t  =  z_t
;pi_salt_tamas_mean_so2!1    = "lat_t"
;pi_salt_tamas_mean_so2&lat_t=  lat_so

;pi_salt_tamas_mean_gb!0    = "z_t"
;pi_salt_tamas_mean_gb&z_t  =  z_t
;pi_salt_tamas_mean_gb!1    = "lat_t"
;pi_salt_tamas_mean_gb&lat_t=  lat_gb

sub3!0    = "z_t"
;sub3&z_t  =  z_t
;sub3&z_t  =  z_t(0:41)
sub3&z_t  =  z_t(0:53)
sub3!1    = "lat_t"
sub3&lat_t=  lat_atlantic

;sub3&lat_t=  lat_so

;#####################################

;pi_n2_tamas_mean!0    = "z_t"
;pi_n2_tamas_mean&z_t  =  z_t
;pi_n2_tamas_mean!1    = "lat_t"
;pi_n2_tamas_mean&lat_t=  lat_atlantic

;pi_n2_tamas_mean_so2!0    = "z_t"
;pi_n2_tamas_mean_so2&z_t  =  z_t
;pi_n2_tamas_mean_so2!1    = "lat_t"
;pi_n2_tamas_mean_so2&lat_t=  lat_so

;pi_n2_tamas_mean_gb!0    = "z_t"
;pi_n2_tamas_mean_gb&z_t  =  z_t
;pi_n2_tamas_mean_gb!1    = "lat_t"
;pi_n2_tamas_mean_gb&lat_t=  lat_gb


sub1!0    = "z_t"
;sub1&z_t  =  z_t
;sub1&z_t  =  z_t(0:41)
sub1&z_t  =  z_t(0:53)
sub1!1    = "lat_t"
sub1&lat_t=  lat_atlantic
;sub1&lat_t=  lat_so

;#####################################

;lgm_temp_tamas_mean!0    = "z_t"
;lgm_temp_tamas_mean&z_t  =  z_t
;lgm_temp_tamas_mean!1    = "lat_t"
;lgm_temp_tamas_mean&lat_t=  lat_atlantic

;lgm_temp_tamas_mean_so2!0    = "z_t"
;lgm_temp_tamas_mean_so2&z_t  =  z_t
;lgm_temp_tamas_mean_so2!1    = "lat_t"
;lgm_temp_tamas_mean_so2&lat_t=  lat_so

;lgm_temp_tamas_mean_gb!0    = "z_t"
;lgm_temp_tamas_mean_gb&z_t  =  z_t
;lgm_temp_tamas_mean_gb!1    = "lat_t"
;lgm_temp_tamas_mean_gb&lat_t=  lat_gb

sub5!0    = "z_t"
;sub5&z_t  =  z_t
;sub5&z_t  =  z_t(0:41)
sub5&z_t  =  z_t(0:53)
sub5!1    = "lat_t"
sub5&lat_t=  lat_atlantic

;sub5&lat_t=  lat_so

;#####################################

;lgm_salt_tamas_mean!0    = "z_t"
;lgm_salt_tamas_mean&z_t  =  z_t
;lgm_salt_tamas_mean!1    = "lat_t"
;lgm_salt_tamas_mean&lat_t=  lat_atlantic

;lgm_salt_tamas_mean_so2!0    = "z_t"
;lgm_salt_tamas_mean_so2&z_t  =  z_t
;lgm_salt_tamas_mean_so2!1    = "lat_t"
;lgm_salt_tamas_mean_so2&lat_t=  lat_so

;lgm_salt_tamas_mean_gb!0    = "z_t"
;lgm_salt_tamas_mean_gb&z_t  =  z_t
;lgm_salt_tamas_mean_gb!1    = "lat_t"
;lgm_salt_tamas_mean_gb&lat_t=  lat_gb


sub6!0    = "z_t"
;sub6&z_t  =  z_t
;sub6&z_t  =  z_t(0:41)
sub6&z_t  =  z_t(0:53)
sub6!1    = "lat_t"
sub6&lat_t=  lat_atlantic

;sub6&lat_t=  lat_so

;#####################################

;lgm_n2_tamas_mean!0    = "z_t"
;lgm_n2_tamas_mean&z_t  =  z_t
;lgm_n2_tamas_mean!1    = "lat_t"
;lgm_n2_tamas_mean&lat_t=  lat_atlantic

;lgm_n2_tamas_mean_so2!0    = "z_t"
;lgm_n2_tamas_mean_so2&z_t  =  z_t
;lgm_n2_tamas_mean_so2!1    = "lat_t"
;lgm_n2_tamas_mean_so2&lat_t=  lat_so

;lgm_n2_tamas_mean_gb!0    = "z_t"
;lgm_n2_tamas_mean_gb&z_t  =  z_t
;lgm_n2_tamas_mean_gb!1    = "lat_t"
;lgm_n2_tamas_mean_gb&lat_t=  lat_gb


sub4!0    = "z_t"
;sub4&z_t  =  z_t(0:58)
;sub4&z_t  =  z_t(0:41)
sub4&z_t  =  z_t(0:53)
sub4!1    = "lat_t"
sub4&lat_t=  lat_atlantic

;sub4&lat_t=  lat_so


;#####################################

;pi_temp_mean!0    = "z_t"
;pi_temp_mean&z_t  =  z_t
;pi_temp_mean!1    = "lat_t"
;pi_temp_mean&lat_t=  lat_atlantic

;pi_temp_mean_so2!0    = "z_t"
;pi_temp_mean_so2&z_t  =  z_t
;pi_temp_mean_so2!1    = "lat_t"
;pi_temp_mean_so2&lat_t=  lat_so

;pi_temp_mean_gb!0    = "z_t"
;pi_temp_mean_gb&z_t  =  z_t
;pi_temp_mean_gb!1    = "lat_t"
;pi_temp_mean_gb&lat_t=  lat_gb

;#####################################

;pi_salt_mean!0    = "z_t"
;pi_salt_mean&z_t  =  z_t
;pi_salt_mean!1    = "lat_t"
;pi_salt_mean&lat_t=  lat_atlantic

;pi_salt_mean_so2!0    = "z_t"
;pi_salt_mean_so2&z_t  =  z_t
;pi_salt_mean_so2!1    = "lat_t"
;pi_salt_mean_so2&lat_t=  lat_so

;pi_salt_mean_gb!0    = "z_t"
;pi_salt_mean_gb&z_t  =  z_t
;pi_salt_mean_gb!1    = "lat_t"
;pi_salt_mean_gb&lat_t=  lat_gb

;#####################################

;pi_n2_mean!0    = "z_t"
;pi_n2_mean&z_t  =  z_t
;pi_n2_mean!1    = "lat_t"
;pi_n2_mean&lat_t=  lat_atlantic

;pi_n2_mean_so2!0    = "z_t"
;pi_n2_mean_so2&z_t  =  z_t
;pi_n2_mean_so2!1    = "lat_t"
;pi_n2_mean_so2&lat_t=  lat_so

;pi_n2_mean_gb!0    = "z_t"
;pi_n2_mean_gb&z_t  =  z_t
;pi_n2_mean_gb!1    = "lat_t"
;pi_n2_mean_gb&lat_t=  lat_gb


;#####################################

;lgm_temp_mean!0    = "z_t"
;lgm_temp_mean&z_t  =  z_t
;lgm_temp_mean!1    = "lat_t"
;lgm_temp_mean&lat_t=  lat_atlantic

;lgm_temp_mean_so2!0    = "z_t"
;lgm_temp_mean_so2&z_t  =  z_t
;lgm_temp_mean_so2!1    = "lat_t"
;lgm_temp_mean_so2&lat_t=  lat_so

;lgm_temp_mean_gb!0    = "z_t"
;lgm_temp_mean_gb&z_t  =  z_t
;lgm_temp_mean_gb!1    = "lat_t"
;lgm_temp_mean_gb&lat_t=  lat_gb

;#####################################

;lgm_salt_mean!0    = "z_t"
;lgm_salt_mean&z_t  =  z_t
;lgm_salt_mean!1    = "lat_t"
;lgm_salt_mean&lat_t=  lat_atlantic

;lgm_salt_mean_so2!0    = "z_t"
;lgm_salt_mean_so2&z_t  =  z_t
;lgm_salt_mean_so2!1    = "lat_t"
;lgm_salt_mean_so2&lat_t=  lat_so

;lgm_salt_mean_gb!0    = "z_t"
;lgm_salt_mean_gb&z_t  =  z_t
;lgm_salt_mean_gb!1    = "lat_t"
;lgm_salt_mean_gb&lat_t=  lat_gb

;#####################################

;lgm_n2_mean!0    = "z_t"
;lgm_n2_mean&z_t  =  z_t
;lgm_n2_mean!1    = "lat_t"
;lgm_n2_mean&lat_t=  lat_atlantic

;lgm_n2_mean_so2!0    = "z_t"
;lgm_n2_mean_so2&z_t  =  z_t
;lgm_n2_mean_so2!1    = "lat_t"
;lgm_n2_mean_so2&lat_t=  lat_so

;lgm_n2_mean_gb!0    = "z_t"
;lgm_n2_mean_gb&z_t  =  z_t
;lgm_n2_mean_gb!1    = "lat_t"
;lgm_n2_mean_gb&lat_t=  lat_gb

;    nrow = 4
;    ncol = 3                           
;    plot = new(nrow*ncol , "graphic")



res_v2 = res
res_v2@lbLabelBarOn        = True
;res_v2@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;res_v2@cnMinLevelValF            = 0.01             ; min level
;res_v2@cnMaxLevelValF            =  10.             ; max level
;res_v2@cnLevelSpacingF           =  0.02            ; contour interval


res3_v2 = res3
res3_v2@lbLabelBarOn        = True
res4_v2 = res4
res4_v2@lbLabelBarOn        = True


;plot(0) = gsn_csm_contour(wks,pi_n2_tamas_mean,res)
;plot(2) = gsn_csm_contour(wks,lgm_n2_tamas_mean,res)
;plot(1) = gsn_csm_contour(wks,pi_n2_mean,res)
;plot(3) = gsn_csm_contour(wks,lgm_n2_mean,res_v2)


;plots(0) = gsn_csm_contour(wks,pi_temp_tamas_mean,res3)
;plots(2) = gsn_csm_contour(wks,lgm_temp_tamas_mean,res3)
;plots(1) = gsn_csm_contour(wks,pi_temp_mean,res3)
;plots(3) = gsn_csm_contour(wks,lgm_temp_mean,res3_v2)


;plots2(0) = gsn_csm_contour(wks,pi_salt_tamas_mean,res4)
;plots2(2) = gsn_csm_contour(wks,lgm_salt_tamas_mean,res4)
;plots2(1) = gsn_csm_contour(wks,pi_salt_mean,res4)
;plots2(3) = gsn_csm_contour(wks,lgm_salt_mean,res4_v2)




;plot(0) = gsn_csm_contour(wks,pi_n2_tamas_mean_so2,res)
;plot(2) = gsn_csm_contour(wks,lgm_n2_tamas_mean_so2,res)
;plot(1) = gsn_csm_contour(wks,pi_n2_mean_so2,res)
;plot(3) = gsn_csm_contour(wks,lgm_n2_mean_so2,res_v2)


;plots(0) = gsn_csm_contour(wks,pi_temp_tamas_mean_so2,res3)
;plots(2) = gsn_csm_contour(wks,lgm_temp_tamas_mean_so2,res3)
;plots(1) = gsn_csm_contour(wks,pi_temp_mean_so2,res3)
;plots(3) = gsn_csm_contour(wks,lgm_temp_mean_so2,res3_v2)


;plots2(0) = gsn_csm_contour(wks,pi_salt_tamas_mean_so2,res4)
;plots2(2) = gsn_csm_contour(wks,lgm_salt_tamas_mean_so2,res4)
;plots2(1) = gsn_csm_contour(wks,pi_salt_mean_so2,res4)
;plots2(3) = gsn_csm_contour(wks,lgm_salt_mean_so2,res4_v2)




;plot(0) = gsn_csm_contour(wks,pi_n2_tamas_mean_gb,res)
;plot(2) = gsn_csm_contour(wks,lgm_n2_tamas_mean_gb,res)
;plot(1) = gsn_csm_contour(wks,pi_n2_mean_gb,res)
;plot(3) = gsn_csm_contour(wks,lgm_n2_mean_gb,res_v2)


;plots(0) = gsn_csm_contour(wks,pi_temp_tamas_mean_gb,res3)
;plots(2) = gsn_csm_contour(wks,lgm_temp_tamas_mean_gb,res3)
;plots(1) = gsn_csm_contour(wks,pi_temp_mean_gb,res3)
;plots(3) = gsn_csm_contour(wks,lgm_temp_mean_gb,res3_v2)


;plots2(0) = gsn_csm_contour(wks,pi_salt_tamas_mean_gb,res4)
;plots2(2) = gsn_csm_contour(wks,lgm_salt_tamas_mean_gb,res4)
;plots2(1) = gsn_csm_contour(wks,pi_salt_mean_gb,res4)
;plots2(3) = gsn_csm_contour(wks,lgm_salt_mean_gb,res4_v2)



;plot(0) = gsn_csm_contour(wks,sub1,res)
;plot(1) = gsn_csm_contour(wks,sub4,res_v2)


;plots(0) = gsn_csm_contour(wks,sub2,res3)
;plots(1) = gsn_csm_contour(wks,sub5,res3_v2)


;plots2(0) = gsn_csm_contour(wks,sub3,res4)
;plots2(1) = gsn_csm_contour(wks,sub6,res4_v2)


plot(0) = gsn_csm_contour(wks,sub1,res)
plot(3) = gsn_csm_contour(wks,sub4,res_v2) 

plot(1) = gsn_csm_contour(wks,sub2,res3)
plot(4) = gsn_csm_contour(wks,sub5,res3_v2) 

plot(2) = gsn_csm_contour(wks,sub3,res4)
plot(5) = gsn_csm_contour(wks,sub6,res4_v2) 


printVarSummary(pi_n2_na) 
print(min(pi_n2_na))
print(max(pi_n2_na))

printVarSummary(lgm_n2_na)
print(min(lgm_n2_na))
print(max(lgm_n2_na))


print(depth)

resP = True

  resP@gsnMaximize        = False
  resP@gsnFrame         = False               ; save panel until both ready
  resP@gsnPaperOrientation = "landscape"
  resP@gsnPanelLabelBar   = False
  resP@pmLabelBarWidthF   = 0.60
  resP@pmLabelBarHeightF = 0.043

  resP@pmLabelBarOrthogonalPosF = -0.042
  resP@lbTitleString = "10^-9 [1/s^2]"

  resP@lbTitleFontHeightF = 0.011
;resP@lbTitlePosition     = "Top"
;resP@lbTitleDirection     = "Across"                          ; letter angle
;resP@lbTitleAngleF        = 270.                               ; title angle
;resP@lbTitleDirection = "Down"
;resP@lbTitleOffsetF = -0.003
;!resP@lbTitleFont   = "helvetica-bold"
resP@lbLabelFontHeightF= 0.010
;!resP@gsnPanelFigureStringsFontHeightF = 0.011
resP@gsnPanelFigureStringsFontHeightF = 0.010
resP@gsnPanelFigureStrings= (/"(a)", "(b)", "(c)", "(d)","(e)", "(f)" /)

;resP@gsnPanelTop = 0.95
;resP@gsnPanelBottom = 0.61
;resP@gsnPanelFigureStringsFontHeightF = 0.012

resP@gsnPanelYWhiteSpacePercent = 3
;!resP@gsnPanelXWhiteSpacePercent = 2
;!resP@gsnPanelXF = (/0.1,0.1,0.1,0.1/)
;!resP@gsnPanelYF = (/0.94,0.73,0.52,0.31/)

;resP@gsnPanelXF = (/0.1,0.1/)
;resP@gsnPanelYF = (/0.94,0.73/)




resP1 = True


resP1@lbTitleString = "Temperature [degC]"

;resP1@lbTitleString = "In-Situ Density [kg/m^3]"

;resP1@gsnPanelTop = 0.64
;resP1@gsnPanelBottom = 0.19


  resP1@gsnMaximize        = False
  resP1@gsnFrame         = False               ; save panel until both ready
  resP1@gsnPaperOrientation = "landscape"
  resP1@gsnPanelLabelBar   = False
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



;!resP1@gsnPanelYWhiteSpacePercent = 3
;!resP1@gsnPanelXWhiteSpacePercent = 2

resP1@lbLabelFontHeightF= 0.010
resP1@gsnPanelFigureStrings= (/"(c)", "(d)"/)
resP1@gsnPanelFigureStringsFontHeightF = 0.011

;resP1@gsnPanelXF = (/0.4,0.4,0.4,0.4/)
;resP1@gsnPanelYF = (/0.94,0.73,0.52,0.31/)


;resP1@gsnPanelXF = (/0.4,0.4/)
;resP1@gsnPanelYF = (/0.94,0.73/)

resP2 = True

resP2@lbTitleString = "Salt Flux [g/kg]"

;resP2@gsnPanelTop = 0.64
;resP2@gsnPanelBottom = 0.19


  resP2@gsnMaximize        = False
  resP2@gsnFrame         = False               ; save panel until both ready
  resP2@gsnPaperOrientation = "landscape"
  resP2@gsnPanelLabelBar   = False
  resP2@pmLabelBarWidthF   = 0.60
  resP2@pmLabelBarHeightF = 0.043
;  resP2@pmLabelBarWidthF   = 0.07
;   resP2@pmLabelBarHeightF = 0.3
   resP2@pmLabelBarOrthogonalPosF = -0.04
   ;resP2@pmLabelBarParallelPosF = -0.01
;   resP2@lbTitleExtentF = 0.01
;resP2@lbOrientation        = "Vertical"
;resP2@lbLabelPosition      = "Right"           ; label position
;resP2@lbLabelOffsetF = -0.3
;resP2@gsnPaperHeightF = 80
;resP2@gsnPaperWidthF = 80

resP2@lbTitleFontHeightF = 0.011
;resP2@lbTitlePosition     = "Top"
;resP2@lbTitleDirection     = "Across"                          ; letter angle
;resP2@lbTitleAngleF        = 270.                               ; title angle
;resP2@lbTitleDirection = "Down"
;resP2@lbTitleOffsetF = -0.003
;!resP2@lbTitleFont   = "helvetica-bold"


;!resP2@gsnPanelYWhiteSpacePercent = 3
;!resP2@gsnPanelXWhiteSpacePercent = 2

resP2@lbLabelFontHeightF= 0.010

resP2@gsnPanelFigureStrings= (/"(e)", "(f)"/)
resP2@gsnPanelFigureStringsFontHeightF = 0.011

;resP2@gsnPanelXF = (/0.7,0.7,0.7,0.7/)
;resP2@gsnPanelYF = (/0.94,0.73,0.52,0.31/)


;resP2@gsnPanelXF = (/0.7,0.7/)
;resP2@gsnPanelYF = (/0.94,0.73/)

;gsn_panel(wks, plot, (/1,2/), resP)
;gsn_panel(wks, plots, (/1,2/), resP1)
;gsn_panel(wks, plots2, (/1,2/), resP2)

gsn_panel(wks, plot, (/2,3/), resP)

;gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots

;txres2 = True
;txres2@txFontHeightF = 0.013
;txres2@txFont   = "helvetica-bold"

;gsn_text_ndc(wks,"N2",0.25,0.98,txres2)
;gsn_text_ndc(wks,"Temperature",0.55,0.98,txres2)
;gsn_text_ndc(wks,"Salt",0.85,0.98,txres2)


txres2 = True
txres2@txFontHeightF = 0.012
txres2@txFont   = "helvetica-bold"

gsn_text_ndc(wks,"N^2",0.19,0.76,txres2)
gsn_text_ndc(wks,"Temperature",0.52,0.76,txres2)
gsn_text_ndc(wks,"Salt",0.86,0.76,txres2)

frame(wks)
end


