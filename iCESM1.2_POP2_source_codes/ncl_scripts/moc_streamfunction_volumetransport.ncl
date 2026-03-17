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


  ;---Open files using addfiles




 ;;pop_pi_tamas = addfile ("/mnt/d/lgm_data2/moc_pi_tamas_last100_timmean.nc", "r")



    pop_pi_tamas = addfile ("/mnt/c/Users/pilat/Downloads/pi_indensitycoord2_new_noidemix.nc","r")


    pop_lgm_tamas = addfile ("/mnt/c/Users/pilat/Downloads/lgm_indensitycoord2_new_noidemixwilmes.nc","r")

    pop_lgm_tamas2 = addfile ("/mnt/c/Users/pilat/Downloads/lgm_indensitycoord2_new_noidemixdinezio.nc","r")   
  

   pop_lgm = addfile ("/mnt/c/Users/pilat/Downloads/lgm_indensitycoord2_new_idemixwilmes.nc", "r")  

    pop_lgm2 = addfile ("/mnt/c/Users/pilat/Downloads/lgm_indensitycoord2_new_idemixdinezio.nc", "r")  
  
  ;;   pop_lgm = addfile ("/mnt/d/lgm_data/moc_lgm_year1000.1900_1999_timmean.nc", "r")
  
    pop_pi = addfile ("/mnt/c/Users/pilat/Downloads/pi_indensitycoord2_new_idemix.nc", "r")   

   

    pop_nc_example = addfile ("/mnt/c/Users/pilat/Downloads/kvmix_pi_tamas_4190_4199_timmean.nc", "r")


  
   wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"

  
   wks  = gsn_open_wks(wks_type,"moc_streamfunction_volumetransport_6panel_new")
 

   cmap = read_colormap_file("MPL_RdYlGn")
   ncmap = dimsizes(cmap(:,0))
 

   wkWidth = 1100
   wkHeight = 680 

res = True
pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(6,graphic)
;plots = new(2,graphic)

res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "
res@gsnPaperMargin = 0.2
res@gsnPaperWidth = 10.0

nltp    = 224               ; latitude subscript
ntm    = 0                  ; time subscript
nlt    = 9               ; latitude subscript
lat    = pop_nc_example->TLAT(nlt,:)
latavg = avg(lat)
lon_o  = pop_nc_example->TLONG(nlt,:)
nlon   = dimsizes (lon_o)
lon_t        = lon_o          ; trick to conserve meta data
lon_t@units  = "degrees_east"
lon_t( 0:283) = lon_o(36:319)
lon_t(284:319) = lon_o( 0:35)
ii           = ind(lon_o.GT. 320.0)
lon_t(ii)    = 360.-lon_o(ii)
lonavg = lon_t(nlt-1)


lon_oh  = pop_nc_example->TLONG(:,:)
nlonh   = dimsizes (lon_oh)
lon_th        = lon_oh          ; trick to conserve meta data
lon_th@units  = "degrees_east"
lon_th( :,0:283) = lon_oh(:,36:319)
lon_th(:,284:319) = lon_oh(:, 0:35)
lon_th3 = lon_th(:,nlt)

lat_t        = lat          ; trick to conserve meta data
lat_t@units  = "degrees_north"
iih           = ind(lon_th3.GT. 180.0)
lon_th3(iih)    = lon_th3(iih)-360.
lonavg2 = avg(lon_th3)

print(lon_th(:,53))


;moc         = pop_pi_tamas->MOC

;printVarSummary(moc)

;mc1         = pop_pi_tamas->MOC(:,1,0,:,115:375)

;mc1         = pop_pi_tamas->MOC(:,1,0,:,105:360)



sigma1         = pop_pi_tamas->sigma
lat_aux_grid1 = pop_pi_tamas->lat_aux_grid


print(sigma1)

amoc1         = pop_pi_tamas->MOC(0,:,:)


;mc1_mask = mask(mc1,(mc1.gt.0.000001.or.mc1.lt.-0.00001),True)
;mc1_mask = mask(mc1,(mc1.gt.0.000085.or.mc1.lt.-0.00001),True)
;printVarSummary(mc1_mask)

 rho_bins = fspan(30, 38, 86)


lat1p    = pop_nc_example->TLAT(:,nltp)
lat_t1p        = lat1p
;lat_aux_grid6p          = mc1_mask&lat_aux_grid
z_1tp           =  rho_bins
z_1tp@long_name = "sigma"; overwrite long name
z_1tp@units     = "kg/m^3"             ; fix units attribute

;amoc1 = dim_avg_n(mc1_mask,0)
printVarSummary(amoc1)

amoc1!0    = "sigma"
amoc1&sigma  =  rho_bins 
amoc1!1    = "lat_aux_grid"
amoc1&lat_aux_grid= lat_aux_grid1



;mc2         = pop_lgm_tamas->MOC(:,1,0,:,115:375)

;mc2         = pop_lgm_tamas->MOC(:,1,0,:,105:360)

;mc2         = pop_lgm_tamas->MOC(:,0,0,:,:)

;mc2_mask = mask(mc2,(mc2.gt.0.000003.or.mc2.lt.-0.00003),True)
;printVarSummary(mc2_mask)


lat2p    = pop_nc_example->TLAT(:,nltp)
lat_t2p        = lat2p
;lat_aux_grid2p          = mc2_mask&lat_aux_grid
;z_2tp           = pop_pi->moc_z
;z_2tp           = z_2tp*1.E-2       ; convert from cm to m
;z_2tp@long_name = ""; overwrite long name
;z_2tp@units     = "m"             ; fix units attribute

;amoc2 = dim_avg_n(mc2_mask,0)
;printVarSummary(amoc2)

amoc2         = pop_lgm_tamas->MOC(0,:,:)

amoc2!0    = "sigma"
amoc2&sigma  =  rho_bins
amoc2!1    = "lat_aux_grid"
amoc2&lat_aux_grid= lat_aux_grid1


ntmp    = 0                  ; time subscript
nltp    = 224               ; latitude subscript
latp    = pop_nc_example->TLAT(nltp,:)
latavgp = avg(latp)
lon_op  = pop_nc_example->TLONG(nltp,:)
nlonp   = dimsizes (lon_op)
lon_tp        = lon_op          ; trick to conserve meta data
lon_tp@units  = "degrees_east"
lon_tp( 0:283) = lon_op(36:319)
lon_tp(284:319) = lon_op( 0:35)
iip           = ind(lon_op.GT. 320.0)
lon_tp(iip)    = 360.-lon_op(iip)
lonavgp = lon_tp(nltp-1)


lon_ohp  = pop_nc_example->TLONG(:,:)
nlonhp   = dimsizes (lon_ohp)
lon_thp        = lon_ohp         ; trick to conserve meta data
lon_thp@units  = "degrees_east"
lon_thp( :,0:283) = lon_ohp(:,36:319)
lon_thp(:,284:319) = lon_ohp(:, 0:35)
lon_th3p = lon_thp(:,nltp)





;mc3          = pop_pi[:]->MOC(:,0,0,:,:)

;moc3 = dim_avg_n(mc3,0)
;lat2p    = pop_lgm10[0]->TLAT(:,nltp)
;lat_t2p        = lat2p
;lat_aux_gridp          = mc3&lat_aux_grid
;z_tp           = pop_lgm10[0]->moc_z
;z_tp           = z_tp*1.E-2       ; convert from cm to m
;z_tp@long_name = ""; overwrite long name
;z_tp@units     = "m"             ; fix units attribute
;moc3!0    = "moc_z"
;moc3&moc_z  =  z_t
;moc3!1    = "lat_aux_grid"
;moc3&lat_aux_grid=  lat_aux_gridp


;mc4         = pop_pi[:]->MOC(:,1,0,:,115:375)

;mc44 = pop_pi->

;mc4         = pop_pi->MOC(:,1,0,:,105:360)

;mc4         = pop_pi[:]->MOC(:,0,0,:,:)

;amoc1 = dim_avg_n(mc4,0)

;mc4_mask = mask(mc4,(mc4.gt.0.00001.or.mc4.lt.-0.00001),True)
;mc4_mask = mask(mc4,(mc4.gt.0.000085.or.mc4.lt.-0.00001),True)

;mc4_mask = mask(mc4,(mc4.gt.0.000085.or.mc4.lt.-0.0001),True)
;printVarSummary(mc4_mask)
;mc4_mask = mask(mc4,(mc4.gt.0.000085.or.mc4.lt.-0.00001),True)


;z_tp4           = pop_pi[0]->moc_z
;z_tp4           = z_tp4*1.E-2 


;mc4_mask!1    = "moc_z"
;mc4_mask&moc_z  =  z_tp4
;mc4_mask!2    = "lat_aux_grid"

;amoc4 = dim_avg_n(mc4_mask,0)
;printVarSummary(amoc1)

;lat4p    = pop_nc_example->TLAT(:,nltp)
;lat_t4p        = lat4p
;lat_aux_grid4p          = mc4&lat_aux_grid
;z_4tp           = pop_pi->moc_z
;z_4tp           = z_4tp*1.E-2       ; convert from cm to m
;z_4tp@long_name = ""; overwrite long name
;z_4tp@units     = "m"             ; fix units attribute


amoc4         = pop_pi->MOC(0,:,:)


amoc4!0    = "sigma"
amoc4&sigma  =  rho_bins
amoc4!1    = "lat_aux_grid"
amoc4&lat_aux_grid= lat_aux_grid1


;mc5         = pop_tm->MOC(:,1,0,:,:)

;amoc2 = dim_avg_n(mc5,0)
;lat5p    = pop_lgm5[0]->TLAT(:,nltp)
;lat_t5p        = lat5p
;lat_aux_grid5p          = mc5&lat_aux_grid
;z_5tp           = pop_lgm10[0]->moc_z
;z_5tp           = z_5tp*1.E-2       ; convert from cm to m
;z_5tp@long_name = ""; overwrite long name
;z_5tp@units     = "m"             ; fix units attribute
;amoc2!0    = "moc_z"
;amoc2&moc_z  =  z_t
;amoc2!1    = "lat_aux_grid"
;amoc2&lat_aux_grid=  lat_aux_grid5p




;mc6         = pop_lgm5->MOC(49:99,1,0,:,115:375)

;mc6         = pop_lgm->MOC(:,1,0,:,105:360)

;mc6         = pop_lgm5->MOC(49:99,0,0,:,:)

;mc6_mask = mask(mc6,(mc6.gt.0.000085.or.mc6.lt.-0.00001),True)

;mc6_mask = mask(mc6,(mc6.gt.0.000003.or.mc6.lt.-0.00003),True)

;mc6_mask = mask(mc6,(mc6.gt.0.000085.or.mc6.lt.-0.0001),True)
;printVarSummary(mc6_mask)



;lat6p    = pop_nc_example->TLAT(:,nltp)
;lat_t6p        = lat6p

;lat_aux_grid6p          = mc6_mask&lat_aux_grid
;z_6tp           = pop_pi->moc_z
;z_6tp           = z_6tp*1.E-2       ; convert from cm to m
;z_6tp@long_name = ""; overwrite long name
;z_6tp@units     = "m"             ; fix units attribute

;amoc6 = dim_avg_n(mc6_mask,0)
;printVarSummary(amoc6)

amoc6         = pop_lgm->MOC(0,:,:)

amoc6!0    = "sigma"
amoc6&sigma  =  rho_bins
amoc6!1    = "lat_aux_grid"
amoc6&lat_aux_grid= lat_aux_grid1



amoc8         = pop_lgm_tamas2->MOC(0,:,:)

amoc8!0    = "sigma"
amoc8&sigma  =  rho_bins
amoc8!1    = "lat_aux_grid"
amoc8&lat_aux_grid= lat_aux_grid1


amoc10         = pop_lgm2->MOC(0,:,:)

amoc10!0    = "sigma"
amoc10&sigma  =  rho_bins
amoc10!1    = "lat_aux_grid"
amoc10&lat_aux_grid= lat_aux_grid1


print(amoc10&lat_aux_grid)

lat_aux         = pop_lgm2->lat_aux_grid



amoc1_30 = amoc1(:,91)
amoc1_50 = amoc1(:,54)

amoc2_30 = amoc2(:,91)
amoc2_50 = amoc2(:,54)

amoc4_30 = amoc4(:,91)
amoc4_50 = amoc4(:,54)

amoc6_30 = amoc6(:,91)
amoc6_50 = amoc6(:,54)

amoc8_30 = amoc8(:,91)
amoc8_50 = amoc8(:,54)

amoc10_30 = amoc10(:,91)
amoc10_50 = amoc10(:,54)

print(avg(amoc4))

res                      = True       ; plot mods desired
res@gsnDraw              = False                ; don't draw yet
res@gsnFrame             = False                ; don't advance frame yet
;res@lbLabelBarOn        = True
;res@cnFillOn             = True       ; turn on color fill
;res@cnFillPalette        = "cmocean_balance" ; set color map
;res@cnFillPalette        = cmap(::-1,:)
;res@cnFillPalette        = "cmp_b2r"
;res@cnLinesOn            = False      ; turn off contour lines
;res@cnLinesOn            = True      ; turn off contour lines
;res@cnLineLabelsOn       = True      ; turns off contour line labels
;res@cnMinLevelValF       = -2
;res@cnMaxLevelValF       =  2
;res@cnLevelSpacingF      =  0.4
;res@cnLabelMasking        = True   
;res@cnLineLabelDensityF = 0.8
;res@cnLineLabelAngleF = 0.0
;res@cnLineLabelPlacementMode = "computed"
;res@cnLineLabelPlacementMode = "constant"
;res@cnLineLabelFontColor = "black"
;res@cnLineLabelFontThicknessF = 1.3
;res@cnLineLabelBackgroundColor = "white"
;res@cnLineLabelFontHeightF = 0.016
;res@gsnContourZeroLineThicknessF = 1.3       ; eliminates zero contour
;res@gsnContourNegLineDashPattern = 1       ; sets negative contours to dash pattern 1
;res@lbOrientation        = "horizontal" ; vertical label bar
;res@gsnYAxisIrregular2Linear = True 
;res@trYReverse           = True       ; reverses y-axis
;res@cnMissingValFillColor = "gray"
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/ -8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24 /)   ; set levels

;res@cnLevelSelectionMode = "ManualLevels" ; manual contour levels
;res@cnMinLevelValF = -22. ; set min contour level
;res@cnMaxLevelValF = 22. ; set max contour level
;res@cnLevelSpacingF = 2 ; set contour spacing


res@gsnLeftString =  ""
res@gsnRightString = ""
;res@tiMainString         = "cross-section (approx): "+sprintf("%5.1f", lonavg)+"W"
;res@tiMainFontHeightF = 0.023
;res@gsnLeftString = "AMOC" 
res@description = "  "
;res@gsnRightString = "[Sv]"
;res@vpXF = 0.4
;res@vpYF = 0.6
res@tiYAxisString = ""
res@tiXAxisString = ""
;res@vpWidthF = 0.6
;res@vpHeightF = 0.6

;res@tmXBTickSpacingF = 20.0
;res@cnMissingValFillColor = "grey52"

;res@vpWidthF           = 0.68 ; longer graphs
;res@vpHeightF          = 0.3


;res2 = res
;res2@tmXBTickSpacingF = 20.0
;res2@cnMissingValFillColor = "grey52"
;res2@tmYLLabelsOn = False
;res2@tmXBBorderOn = False



res1                      = True       ; plot mods desired
res1@gsnDraw              = False                ; don't draw yet
res1@gsnFrame             = False                ; don't advance frame yet
res1@lbLabelBarOn        = True
res1@cnFillOn             = True       ; turn on color fill
;res1@cnFillPalette        = cmap(::-1,:)
;res1@cnFillPalette        = "MPL_RdYlGn" ; set color map
res1@cnFillPalette        = "cmp_b2r"
res1@cnLinesOn            = True      ; turn off contour lines
res1@cnLineLabelsOn       = True      ; turns off contour line labels
res1@cnMinLevelValF       = -16
res1@cnMaxLevelValF       =  16
res1@cnLevelSpacingF      =  2 
res1@cnLineLabelAngleF = 0.0
;res1@cnLineLabelPlacementMode = "computed"
res1@cnLineLabelPlacementMode = "constant"
;res1@cnLineLabelPlacementMode = "randomized"
res1@cnLineLabelFontColor = "white"
res1@cnLineLabelFontThicknessF = 1.3
res1@cnLineLabelBackgroundColor = "white"
res1@cnLabelMasking        = True   
res1@cnLineLabelFontHeightF = 0.0165
res1@cnLineLabelDensityF = 1.0
res1@cnLineLabelInterval = 3
res1@gsnContourZeroLineThicknessF = 1.3 	    ; eliminates zero contour
res1@gsnContourNegLineDashPattern = 1       ; sets negative contours to dash pattern 1
res1@cnInfoLabelOn        = False      ; turns off contour info label
res1@lbOrientation        = "horizontal" ; vertical label bar
res1@gsnYAxisIrregular2Linear = True
res1@trYReverse           = True       ; reverses y-axis
res1@cnMissingValFillColor = "gray"
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/ -8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24 /)   ; set levels
res1@description = "  "
res1@gsnLeftString =  ""
res1@gsnRightString= ""
;res1@tiMainString         = "cross-section (approx): "+sprintf("%5.1f", lonavg)+"W"
;res1@tiMainFontHeightF = 0.023
;res1@gsnLeftString = "AMOC"

;res1@gsnRightString = "[Sv]"
;res1@vpXF = 0.4
;res1@vpYF = 0.6
res1@tiYAxisString = ""
res1@tiXAxisString = ""

;res1@vpWidthF = 0.6
;res1@vpHeightF = 0.6

;;res1@cnLevelSelectionMode = "ManualLevels" ; manual contour levels
;res1@cnMinLevelValF = -16. ; set min contour level
;res1@cnMaxLevelValF = 16. ; set max contour level
;res1@cnLevelSpacingF = 2 ; set contour spacing

;res1@cnLevelSpacingF = 2.0 

res1@tmXBTickSpacingF = 20.0
res1@cnMissingValFillColor = "grey52"


res1@pmLabelBarWidthF   = 0.09
res1@pmLabelBarHeightF = 0.35
;res1@pmLabelBarParallelPosF = 0.08
res1@lbTitleExtentF = 0.18
res1@lbTitleOffsetF = 0.02
res1@lbOrientation        = "Vertical"
res1@lbLabelPosition      = "Right"           ; label position
;res1@lbTitleString = "Sv"
res1@lbTitleFontHeightF = 0.010
res1@lbTitlePosition     = "Right"
res1@lbTitleDirection     = "Across"                          ; letter angle
res1@lbTitleAngleF        = 270.                               ; title angle
res1@lbTitleOffsetF = 0.15
res1@lbLabelFontHeightF= 0.017
res1@tmYLLabelsOn = False

res1@vpWidthF           = 0.6 ;longer graphs
res1@vpHeightF          = 0.3
res1@pmLabelBarOrthogonalPosF = 0.01
;res1@cnLineLabelDensityF = 2.0
;res1@cnLineDashSegLenF  = 0.4
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/ -16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16 /)   ; set levels
res1@cnLineColor = "gray22"
res1@cnLineLabelFontColor = "black"
res1@tmXBLabelFontHeightF =  0.0165
res1@tmYLLabelFontHeightF =  0.0165


   res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   res1@cnMinLevelValF            = -30             ; min level
   res1@cnMaxLevelValF            =  30             ; max level
   res1@cnLevelSpacingF           =  3            ; contour interval


res2 = res1
res2@tmXBTickSpacingF = 20.0
res2@cnMissingValFillColor = "grey52"
res2@tmYLLabelsOn = True
;res2@tmXBBorderOn = False
res2@lbLabelFontHeightF= 0.017

res2@lbLabelBarOn        = False
res2@pmLabelBarWidthF   = 0.12
res2@pmLabelBarHeightF = 0.22
res2@pmLabelBarParallelPosF = 0.1
res2@lbTitleExtentF = 0.18
res2@lbTitleOffsetF = 0.02
res2@lbOrientation        = "Vertical"
res2@lbLabelPosition      = "Right"           ; label position
;res2@lbTitleString = "Sv"
res2@lbTitleFontHeightF = 0.010
res2@lbTitlePosition     = "Right"
res2@lbTitleDirection     = "Across"                          ; letter angle
res2@lbTitleAngleF        = 270.                               ; title angle
res2@lbTitleOffsetF = 0.15
res2@lbLabelFontHeightF= 0.0165
res2@pmLabelBarOrthogonalPosF = 0.01           ; move whole thing up
;res2@pmLabelBarParallelPosF = 0.1

res2@vpWidthF           = 0.6 ;longer graphs
res2@vpHeightF          = 0.3
;res2@cnLineLabelDensityF = 2.0
;res2@cnLineDashSegLenF  = 0.4
res2@tmXBLabelFontHeightF =  0.0165
res2@tmYLLabelFontHeightF =  0.0165





;res5@tmXBLabelFontHeightF =  0.022
;res5@tmYLLabelFontHeightF =  0.022

res5 = res
res5                      = True       ; plot mods desired
res5@gsnDraw              = False                ; don't draw yet
res5@gsnFrame             = False                ; don't advance frame yet
res5@lbLabelBarOn        = True
res5@cnFillOn             = True       ; turn on color fill

res5@cnFillPalette        = "cmocean_balance"
res5@pmLabelBarWidthF   = 0.09
res5@pmLabelBarHeightF = 0.33
res5@lbLabelFontHeightF= 0.015
res5@pmLabelBarOrthogonalPosF = 0.01
res5@tmXBLabelFontHeightF =  0.0165
res5@tmYLLabelFontHeightF =  0.0165
res5@tmXBTickSpacingF = 20.00
res5@cnMissingValFillColor = "grey52"
res5@tmYLLabelsOn = False
;res5@tmXBBorderOn = False
res5@lbLabelBarOn        = True
;res5@cnMinLevelValF       = -1.5
;res5@cnMaxLevelValF       =  1.5
;res5@cnLevelSpacingF      =  0.3


;!!!res5@cnLevelSelectionMode = "ExplicitLevels"
;res5@cnLevels    = (/ -1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5 /)   ; set levels
;!!!res5@cnLevels    = (/ -2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5,1.8,2.1 /)   ; set levels


   res5@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   res5@cnMinLevelValF            = -3             ; min level
   res5@cnMaxLevelValF            =  3             ; max level
   res5@cnLevelSpacingF           =  0.2            ; contour interval

res5@pmLabelBarHeightF = 0.35
;res5@pmLabelBarParallelPosF = 0.1
res5@lbTitleExtentF = 0.18
res5@lbTitleOffsetF = 0.02
res5@lbOrientation        = "Vertical"
res5@lbLabelPosition      = "Right"           ; label position
;res5@lbTitleString = "Sv"
res5@lbTitleFontHeightF = 0.010
res5@lbTitlePosition     = "Right"
res5@lbTitleDirection     = "Across"                          ; letter angle
res5@lbTitleAngleF        = 270.                               ; title angle
res5@lbTitleOffsetF = 0.15
res5@lbLabelFontHeightF= 0.017
res5@pmLabelBarOrthogonalPosF = 0.01           ; move whole thing up
;;res5@pmLabelBarParallelPosF = 0.1

;res5@cnLineLabelPlacementMode = "computed"
res5@cnLineLabelPlacementMode = "constant"
;;res5@cnLineLabelPlacementMode = "randomized"
;res5@cnLineColor = "gray22"
;res5@cnLineLabelFontColor = "black"
;res5@cnLineLabelFontThicknessF = 1.3
res5@cnLineLabelBackgroundColor = "white"
;res5@cnLabelMasking        = True
;res5@cnLineLabelFontHeightF = 0.016
;;res5@cnLineLabelInterval = 2

;res5@cnLineLabelFormat    = "0@;*.2f"

res5@vpWidthF           = 0.6 ;longer graphs
res5@vpHeightF          = 0.3
;;res5@cnLevelSpacingF = 2.0
;res5@cnLineLabelDensityF = 2.0

res5@cnLineLabelsOn             = False  
res5@cnLinesOn             = False  


res4 = res5
res4@tmXBTickSpacingF = 20.0
res4@cnMissingValFillColor = "grey52"
res4@tmYLLabelsOn = True
;res4@tmXBBorderOn = False
res4@lbLabelBarOn        = False
;res4@lbLabelBarOn        = True
;res4@cnMinLevelValF       = -1.5
;res4@cnMaxLevelValF       =  1.5
;res4@cnLevelSpacingF      =  0.3


;!!!res4@cnLevelSelectionMode = "ExplicitLevels"
;res4@cnLevels    = (/ -1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5 /)   ; set levels
;!!!res4@cnLevels    = (/ -2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5,1.8,2.1 /)   ; set levels



res4@pmLabelBarHeightF = 0.33
;res4@pmLabelBarParallelPosF = 0.08
res4@lbTitleExtentF = 0.18
res4@lbTitleOffsetF = 0.02
res4@lbOrientation        = "Vertical"
res4@lbLabelPosition      = "Right"           ; label position
;res4@lbTitleString = "Sv"
res4@lbTitleFontHeightF = 0.010
res4@lbTitlePosition     = "Right"
res4@lbTitleDirection     = "Across"                          ; letter angle
res4@lbTitleAngleF        = 270.                               ; title angle
res4@lbTitleOffsetF = 0.15
res4@lbLabelFontHeightF= 0.017
;res4@pmLabelBarOrthogonalPosF = 0.01           ; move whole thing up


;res4@cnLineLabelDensityF = 2.0

res4@cnLineLabelFormat    = "0@;*.2f"



   res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   res4@cnMinLevelValF            = -3             ; min level
   res4@cnMaxLevelValF            =  3             ; max level
   res4@cnLevelSpacingF           =  0.2            ; contour interval



;mm = pop_pi[:]->MOC(:,1,0,:,115:375)
;m4 = dim_avg_n_Wrap(mm,(/1,2/))
;printVarSummary(m4)


;data=new( (/1,156/),float)   ; Method 1
;data(0,0:155)= m4


res4@vpWidthF           = 0.6 ;longer graphs
res4@vpHeightF          = 0.3


print(max(amoc1))
print(min(amoc1))

print(max(amoc4))
print(min(amoc4))


print(max(amoc2))
print(min(amoc2))


print(max(amoc6))
print(min(amoc6))




amoc1_v2 = abs(amoc1)
amoc2_v2 = abs(amoc2)
amoc4_v2 = abs(amoc4)
amoc6_v2 = abs(amoc6)



diff1 = amoc4_v2 - amoc1_v2

diff2 = amoc6_v2 - amoc2_v2

printVarSummary(diff1)
printVarSummary(diff2)

print(max(diff1))
print(min(diff1))

print(max(diff2))
print(min(diff2))




diff1!0    = "sigma"
diff1&sigma  = rho_bins 
diff1!1    = "lat_aux_grid"
printVarSummary(diff1)
;diff1&lat_aux_grid= pop_pi->lat_aux_grid(105:360) ; AMOC
diff1&lat_aux_grid= pop_pi->lat_aux_grid  ; global MOC


diff2!0    = "sigma"
diff2&sigma  =  rho_bins
diff2!1    = "lat_aux_grid"
;diff2&lat_aux_grid= pop_pi_tamas->lat_aux_grid(105:360) ; AMOC
diff2&lat_aux_grid= pop_pi_tamas->lat_aux_grid ; global MOC


;!plot(0) = gsn_csm_contour(wks,diff1,res)     ; create plot
;!plot(1) = gsn_csm_contour(wks,diff2,res2)     ; create plot
;plot(0) = gsn_csm_contour(wks,amoc4,res2)     ; create plot
;!plot(1) = gsn_csm_contour(wks,amoc2,res1)     ; create plot
;plot(1) = gsn_csm_contour(wks,amoc6,res1)     ; create plot
;plot(2) = gsn_csm_contour(wks,diff1,res4)     ; create plot
;plot(3) = gsn_csm_contour(wks,diff2,res5)     ; create plot
;!plot(3) = gsn_csm_contour(wks,diff2,res4)     ; create plot

;plot(0) = gsn_csm_contour(wks,amoc1,res2)     ; create plot
;plot(1) = gsn_csm_contour(wks,amoc2,res1)     ; create plot



amoc_sigma = amoc6&sigma



 data1      = new((/2,86/),double)

 printVarSummary(data1)
 printVarSummary(amoc1_30)

 data1(0,:) = amoc1_30(:)
 data1(1,:) = amoc1_50(:)
 
 data2      = new((/2,86/),double)

 data2(0,:) = amoc2_30(:)
 data2(1,:) = amoc2_50(:)

 data4      = new((/2,86/),double)

 printVarSummary(data4)
 printVarSummary(amoc4_30)

 data4(0,:) = amoc4_30(:)
 data4(1,:) = amoc4_50(:)

 data6      = new((/2,86/),double)

 printVarSummary(data6)
 printVarSummary(amoc6_30)

 data6(0,:) = amoc6_30(:)
 data6(1,:) = amoc6_50(:)

 data8      = new((/2,86/),double)

 printVarSummary(data8)
 printVarSummary(amoc8_30)

 data8(0,:) = amoc8_30(:)
 data8(1,:) = amoc8_50(:)

 data10     = new((/2,86/),double)

 printVarSummary(data10)
 printVarSummary(amoc10_30)

 data10(0,:) = amoc10_30(:)
 data10(1,:) = amoc10_50(:)


;#######################################################################

  lat30_ind1 = ind_nearest_coord(-30.,amoc1&lat_aux_grid, 0)
  lat50_ind1 = ind_nearest_coord(-50.,amoc1&lat_aux_grid, 0)

  psi30_v1 = amoc1(:,lat30_ind1)
  psi50_v1 = amoc1(:,lat50_ind1)

   print(psi30_v1)
   print(psi50_v1)

lat30_ind2 = ind_nearest_coord(-30.,amoc2&lat_aux_grid, 0)
  lat50_ind2 = ind_nearest_coord(-50.,amoc2&lat_aux_grid, 0)

  psi30_v2 = amoc2(:,lat30_ind2)
  psi50_v2 = amoc2(:,lat50_ind2)


lat30_ind4 = ind_nearest_coord(-30.,amoc4&lat_aux_grid, 0)
  lat50_ind4 = ind_nearest_coord(-50.,amoc4&lat_aux_grid, 0)

  psi30_v4 = amoc4(:,lat30_ind4)
  psi50_v4 = amoc4(:,lat50_ind4)


lat30_ind6 = ind_nearest_coord(-30.,amoc6&lat_aux_grid, 0)
  lat50_ind6 = ind_nearest_coord(-50.,amoc6&lat_aux_grid, 0)

  psi30_v6 = amoc6(:,lat30_ind6)
  psi50_v6 = amoc6(:,lat50_ind6)

lat30_ind8 = ind_nearest_coord(-30.,amoc8&lat_aux_grid, 0)
  lat50_ind8 = ind_nearest_coord(-50.,amoc8&lat_aux_grid, 0)

  psi30_v8 = amoc8(:,lat30_ind8)
  psi50_v8 = amoc8(:,lat50_ind8)

lat30_ind10 = ind_nearest_coord(-30.,amoc10&lat_aux_grid, 0)
  lat50_ind10 = ind_nearest_coord(-50.,amoc10&lat_aux_grid, 0)

  psi30_v10 = amoc10(:,lat30_ind10)
  psi50_v10 = amoc10(:,lat50_ind10)



;-----------------------------
; 3. Define 6 isopycnal bins
;   (EDIT TO MATCH YOUR CASE)
;-----------------------------
  ;sigma_edges = (/ 36.5, 36.8, 37.0, 37.2, 37.4, 37.6, 37.8 /)
  
  ;sigma_edges = (/ 32, 32.5, 33.8, 34.1, 34.4, 34.7, 35, 35.3, 35.6, 35.9, 36.2, 36.5, 36.8, 37.0, 37.2, 37.4, 37.6, 37.8 /)
  
   ;sigma_edges = (/ 30, 31.5, 32.5, 33.2, 34.1, 34.7, 35, 36.2, 36.8, 37.0, 37.8 /)


;sigma_edges = (/ 30, 31.2, 32.1, 33.6, 36, 37.5, 38 /)

   ;sigma_edges = (/ 32, 34.3, 36.6, 36.9, 37.1, 37.3, 37.5, 37.7, 37.9  /)

   sigma_edges = (/ 30, 31, 32, 33, 34, 35, 36, 37, 38 , 39/)



   nbins = dimsizes(sigma_edges) - 1


  T1 = new((/2,nbins/), typeof(psi30_v1))   ; dimension (latitude, bin)
  T2 = new((/2,nbins/), typeof(psi30_v2))   ; dimension (latitude, bin)
T4 = new((/2,nbins/), typeof(psi30_v4))   ; dimension (latitude, bin)
T6 = new((/2,nbins/), typeof(psi30_v6))   ; dimension (latitude, bin)
T8 = new((/2,nbins/), typeof(psi30_v8))   ; dimension (latitude, bin)
T10 = new((/2,nbins/), typeof(psi30_v10))   ; dimension (latitude, bin)


sigma_mid = new(nbins, typeof(amoc6&sigma))

  do i=0,nbins-1
    sigma_mid(i) = 0.5*(sigma_edges(i) + sigma_edges(i+1))
  end do



;-----------------------------
; 4. Compute transports
;-----------------------------
  do i = 0, nbins-1

    ind1 = ind_nearest_coord(sigma_edges(i),   amoc1&sigma, 0)
    ind2 = ind_nearest_coord(sigma_edges(i+1), amoc1&sigma, 0)

    print(ind1)
    print(ind2)

    T1(0,i) = psi30_v1(ind2) - psi30_v1(ind1)
    T1(1,i) = psi50_v1(ind2) - psi50_v1(ind1)

  end do


  do i = 0, nbins-1

    ind1 = ind_nearest_coord(sigma_edges(i),   amoc2&sigma, 0)
    ind2 = ind_nearest_coord(sigma_edges(i+1), amoc2&sigma, 0)

    print(ind1)
    print(ind2)

    T2(0,i) = psi30_v2(ind2) - psi30_v2(ind1)
    T2(1,i) = psi50_v2(ind2) - psi50_v2(ind1)

  end do


  do i = 0, nbins-1

    ind1 = ind_nearest_coord(sigma_edges(i),   amoc4&sigma, 0)
    ind2 = ind_nearest_coord(sigma_edges(i+1), amoc4&sigma, 0)

    print(ind1)
    print(ind2)

    T4(0,i) = psi30_v4(ind2) - psi30_v4(ind1)
    T4(1,i) = psi50_v4(ind2) - psi50_v4(ind1)

  end do


  do i = 0, nbins-1

    ind1 = ind_nearest_coord(sigma_edges(i),   amoc6&sigma, 0)
    ind2 = ind_nearest_coord(sigma_edges(i+1), amoc6&sigma, 0)

    print(ind1)
    print(ind2)

    T6(0,i) = psi30_v6(ind2) - psi30_v6(ind1)
    T6(1,i) = psi50_v6(ind2) - psi50_v6(ind1)

  end do


  do i = 0, nbins-1

    ind1 = ind_nearest_coord(sigma_edges(i),   amoc8&sigma, 0)
    ind2 = ind_nearest_coord(sigma_edges(i+1), amoc8&sigma, 0)

    print(ind1)
    print(ind2)

    T8(0,i) = psi30_v8(ind2) - psi30_v8(ind1)
    T8(1,i) = psi50_v8(ind2) - psi50_v8(ind1)

  end do


  do i = 0, nbins-1

    ind1 = ind_nearest_coord(sigma_edges(i),   amoc10&sigma, 0)
    ind2 = ind_nearest_coord(sigma_edges(i+1), amoc10&sigma, 0)

    print(ind1)
    print(ind2)

    T10(0,i) = psi30_v10(ind2) - psi30_v10(ind1)
    T10(1,i) = psi50_v10(ind2) - psi50_v10(ind1)

  end do






;#######################################################################



 resx = res
 resx@vpHeightF = 0.43        ; Changes the aspect ratio
 resx@vpWidthF  = 0.65
 resx@gsnMaximize = True
 resx@tiYAxisString = "" ; y-axis label      
 resx@tiXAxisString = ""
; resx@tiMainString  = "Darwin Southern Oscillation Index" ; title
 resx@gsnXRefLine           = 0.              ; reference line   
 resx@gsnXYBarChart         = True            ; create bar chart 
 ; resx@gsnRightXRefLineColor = "red"            ; right of ref line fill red
 ;resx@gsnLeftXRefLineColor  = "blue"          ; left of ref line fill blue
 resx@gsnXYBarChartColors = (/"red","blue"/)

;**********************************************************************
; now, change the outline to be the color of the fill. This is useful when
; you have many data points b/c the outline will normally obscure the colors
;**********************************************************************
 ;resx@xyLineColors = (/"red","blue"/)         ; colors for the above/blow

resx@gsnYAxisIrregular2Linear = True 
resx@trYReverse           = True       ; reverses y-axis

 resx2 = resx
 ;resx2@tmYLLabelsOn = False
;resx2@tmXBBorderOn = False
resx2@gsnXYBarChart         = True 
resx2@gsnXRefLine           = 0.              ; reference line   
 resx2@gsnXYBarChart         = True            ; create bar chart 
 ;resx2@gsnRightXRefLineColor = (/"red","orange","black"/)           ; right of ref line fill red

 resx2@gsnXYBarChartColors = (/"red","blue"/)

; res2@gsnYRefLine           = 0.              ; reference line   
 res2@gsnXYBarChart         = True            ; create bar chart 
; res2@gsnAboveYRefLineColor = "red"           ; above ref line fill red
; res2@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue


resx2@gsnYAxisIrregular2Linear = True 
resx2@trYReverse           = True       ; reverses y-axis

; res1@gsnYRefLine           = 0.              ; reference line   
 res1@gsnXYBarChart         = True            ; create bar chart 
; res1@gsnAboveYRefLineColor = "red"           ; above ref line fill red
; res1@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue


resx@trXMinF   = -20.0        ; min value on y-axis
resx@trXMaxF   =  30.0        ; max value on y-axis

resx2@trXMinF   = -20.0        ; min value on y-axis
resx2@trXMaxF   =  30.0        ; max value on y-axis


 resx@tmYLMode   = "Explicit"
  resx@tmYLValues = sigma_mid
  resx@tmYLLabels = sprintf("%.2f", sigma_mid)

  resx@tmXMajorGrid                = True          ; implement x grid 
  resx@tmXMajorGridThicknessF      = 0.5           ; 2.0 is default
  resx@tmXMajorGridLineDashPattern = 0             ; select short dash lines
  
  resx@tmYMajorGrid                = True          ; implement y grid
  resx@tmYMajorGridThicknessF      = 0.5           ; 2.0 is default
  resx@tmYMajorGridLineDashPattern = 0             ; select short dash lines

;-------------------------------------
; Make bar thickness match density bin
;--------------------------------------
  bin_height = sigma_edges(1) - sigma_edges(0)
  resx@gsnXYBarChartBarWidth = bin_height * 0.9

; Optional: emphasize zero line
  resx2@tmXBMinorOn = True


 resx2@tmYLMode   = "Explicit"
  resx2@tmYLValues = sigma_mid
  resx2@tmYLLabels = sprintf("%.2f", sigma_mid)



  resx2@tmXMajorGrid                = True          ; implement x grid 
  resx2@tmXMajorGridThicknessF      = 0.5           ; 2.0 is default
  resx2@tmXMajorGridLineDashPattern = 0             ; select short dash lines
  
  resx2@tmYMajorGrid                = True          ; implement y grid
  resx2@tmYMajorGridThicknessF      = 0.5           ; 2.0 is default
  resx2@tmYMajorGridLineDashPattern = 0             ; select short dash lines

;--------------------------------------
; Make bar thickness match density bin
;--------------------------------------
  bin_height = sigma_edges(1) - sigma_edges(0)
  resx2@gsnXYBarChartBarWidth = bin_height * 0.9

; Optional: emphasize zero line
  resx2@tmXBMinorOn = True





 xxx = ispan(1,nbins,1)

sigma_xxx = new(nbins, typeof(amoc6&sigma))

do i=0,nbins-1
    sigma_xxx(i) = sigma_edges(i)
  end do

  x2 = (/  31.5, 32.5, 33.2, 34.1, 34.7, 35, 36.2, 36.8, 37.0, 37.8 /)

plot(0) = gsn_csm_xy(wks,T1,sigma_xxx,resx)     ; create plot
plot(2) = gsn_csm_xy(wks,T2,sigma_xxx,resx)     ; create plot
plot(1) = gsn_csm_xy(wks,T4,sigma_xxx,resx)     ; create plot
plot(3) = gsn_csm_xy(wks,T6,sigma_xxx,resx2)     ; create plot
plot(4) = gsn_csm_xy(wks,T8,sigma_xxx,resx)     ; create plot
plot(5) = gsn_csm_xy(wks,T10,sigma_xxx,resx2)     ; create plot


printVarSummary(diff1)
printVarSummary(diff2)

print(max(diff1))
print(min(diff1))

print(max(diff2))
print(min(diff2))



nboxes = dimsizes(resx@gsnXYBarChartColors)

   labels1 = (/"30S"/)
   labels2 = (/"50S"/)


   lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.05     ; labelbar width
  lbres@vpHeightF          = 0.012      ; labelbar height
  ;lbres@lbBoxMajorExtentF  = 0.75          ; puts space between color boxes
  ;lbres@lbFillColors       = resx@gsnXYBarChartColors ; labelbar colors
  lbres@lbFillColors       = "red"
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.0076         ; font height. default is small
  ;lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
   
  gsn_labelbar_ndc(wks,1,labels1,0.435,0.835,lbres)
   
lbres1 = lbres 
lbres1@lbFillColors       = "blue"

gsn_labelbar_ndc(wks,1,labels2,0.435,0.815,lbres1)  
  

resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
;resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = False
resP@pmLabelBarWidthF   = 0.07
resP@pmLabelBarHeightF = 0.33
;resP@lbOrientation   = "vertical"
;resP@pmLabelBarOrthogonalPosF = -0.006
resP@lbTitleString = " Sverdrups"
resP@lbTitleFontHeightF = 0.012
resP@gsnMaximize         = True                 ; make as large as possible



resP@pmLabelBarParallelPosF = 0.05
resP@lbTitleExtentF = 0.18
resP@lbTitleOffsetF = 0.02
resP@lbLabelFontHeightF= 0.011
resP@lbOrientation        = "Vertical"
resP@lbLabelPosition      = "Right"           ; label position





;resP@lbTitleString = "log10[m^2/s]"
resP@lbTitleFontHeightF = 0.010
resP@lbTitlePosition     = "Right"
resP@lbTitleDirection     = "Across"                          ; letter angle
resP@lbTitleAngleF        = 270.                               ; title angle
resP@lbTitleOffsetF = 0.15





;resP@gsnPanelFigureStringsFontHeightF = 0.012
resP@gsnPanelYWhiteSpacePercent = 2
;resP@gsnPanelXWhiteSpacePercent = 0.4



;resP@gsnPanelFigureStrings= (/"PI", "LGM", "PI", "LGM", "PI", "LGM"/)
resP@gsnPanelFigureStrings= (/"a) PI","b) PI","c) LGM-Wilmes","d) LGM-Wilmes", "e) LGM-Dinezio","f) LGM-Dinezio"/)
resP@gsnPanelFigureStringsFontHeightF = 0.0085
;;resP@cnLevelSelectionMode = "ExplicitLevels"
;;resP@cnLevels    = (/ -7, -6, -5, -4, -3, -2 /)   ; set levels
;;resP@gsnPanelTop = 0.2
;resP@gsnPanelYF = (/0.7,0.7/)
;resP@gsnPanelXF = (/0.25,0.55/)

;resP@gsnPanelBottom = 0.45
resP@gsnPanelTop = 0.85
resP@gsnPanelBottom = 0.15
;!resP@gsnPanelLeft = 0.1
resP@gsnLeftString =  " Tamas"
resP@gsnRightString           = "(Sv)"
resP@tiXAxisString = ""

resPs                 = False
resPs@gsnMaximize        = True
resPs@gsnFrame         = False               ; save panel until both ready
resPs@gsnPaperOrientation = "landscape"
resPs@gsnPanelLabelBar   = True
resPs@pmLabelBarWidthF   = 0.38
resPs@pmLabelBarHeightF = 0.04
resPs@lbOrientation   = "horizontal"
resPs@pmLabelBarOrthogonalPosF = -0.01
resPs@lbTitleString = " Sverdrups"
resPs@lbTitleFontHeightF = 0.012
resPs@gsnMaximize         = True                 ; make as large as possible

resPs@lbLabelFontHeightF= 0.009
resPs@gsnPanelFigureStringsFontHeightF = 0.012
;resPs@gsnPanelFigureStringsFontHeightF = 0.014
;resPs@gsnPanelYWhiteSpacePercent = 0.5
resPs@gsnPanelFigureStrings= (/"PI - with IDEMIX", "LGM - with IDEMIX", "AMOC - LGM-j10"/)
resPs@gsnPanelFigureStringsFontHeightF = 0.012
;resPs@cnLevelSelectionMode = "ExplicitLevels"
;resPs@cnLevels    = (/ -7, -6, -5, -4, -3, -2 /)   ; set levels
;resPs@gsnPanelXWhiteSpacePercent = 0.4

;;resPs@gsnPanelBottom = 0.1
resPs@gsnPanelTop = 0.45
;resPs@gsnPanelLeft = 0.1
resPs@gsnLeftString =  " Wilmes"
resPs@gsnRightString           = "(Sv)"
;resPs@gsnPanelRight = 0.65
resPs@tiYAxisString = ""
;resPs@gsnPanelYF = (/0.4,0.4/)
;resPs@gsnPanelXF = (/0.25,0.55/)


txres1 = True
txres1@txFontHeightF = 0.011
txres1@txFont   = "helvetica-bold"

;gsn_text_ndc(wks,"MOC",0.12,0.69,txres1) ; to plot normal data
;gsn_text_ndc(wks,"AMOC, w/IDEMIX - w/oIDEMIX",0.24,0.41,txres1)
;gsn_text_ndc(wks,"AMOC, w/oIDEMIX",0.19,0.95,txres1)
;gsn_text_ndc(wks,"AMOC, w/IDEMIX",0.19,0.68,txres1)

;gsn_text_ndc(wks,"w/IDEMIX - w/oIDEMIX",0.24,0.41,txres1)
gsn_text_ndc(wks,"w/oIDEMIX",0.35,0.87,txres1)
gsn_text_ndc(wks,"w/IDEMIX",0.7,0.87,txres1)


;gsn_text_ndc(wks,"AMOC, w/IDEMIX - w/oIDEMIX",0.21,0.43,txres1)
;gsn_text_ndc(wks,"AMOC, w/IDEMIX",0.15,0.83,txres1)
;gsn_text_ndc(wks,"AMOC, w/oIDEMIX",0.15,0.67,txres1)


txres4 = True
txres4@txFontHeightF = 0.009

gsn_text_ndc(wks,"[Sv]",0.35,0.15,txres4)
gsn_text_ndc(wks,"[Sv]",0.7,0.15,txres4)
;gsn_text_ndc(wks,"[Sv]",0.95,0.29,txres4)

 
;  txres5 = True

;  txid3 = gsn_create_text(wks, "Depth (m)", txres5)


txres5 = True
txres5@txFontHeightF = 0.010
txres5@txAngleF      = 90.     ; rotate 90 degrees
 ;gsn_text_ndc(wks,":F8:s_2000",0.135,0.51,txres5)
gsn_text_ndc(wks," Potential Density - 1000 (kg/m^3)",0.135,0.51,txres5)
;gsn_text_ndc(wks,"Depth (m)",0.04,0.54,txres5)
;gsn_text_ndc(wks,"Depth (m)",0.032,0.31,txres5)





;;resP@gsnPanelYF = (/0.88,0.88,0.61,0.61,0.34,0.34/)
;resP@gsnPanelXF = (/0.07,0.51,0.07,0.51,0.07,0.51/)

;;resP@gsnPanelXF = (/0.12,0.51,0.12,0.51,0.12,0.51/)

 pres  = True

;resP@pmLabelBarOrthogonalPosF = -0.02    ; Move closer to plots
;resPs@pmLabelBarOrthogonalPosF = -0.02

  gsn_panel(wks,plot,(/3,2/),resP)
;  gsn_panel(wks,plot,(/3,2/),resP)




; pres@gsnPanelTop    = 0.91  ; set these if desired
;; pres@gsnPanelBottom = 0.1  ; must be between 0 and 1
;; pres@gsnPanelLeft   = 0.65
;; pres@gsnPanelRight  = 1.0

;pres@gsnPanelRowSpec = True


;  panel_two_sets(wks,plot,plots,(/1,2/),(/1,2/),resP,resPs,pres)


  frame(wks)

end
