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



function get_color_rgba1(value, levels, cmap)

  local nlev
 
 begin
    nlev = dimsizes(levels)

    ; below first level
    if (value.lt.levels(0)) then
      return(cmap(0,:))
    end if

    ; above last level
    if (value.gt.levels(nlev-1)) then
      return(cmap(nlev-1,:))
    end if

    ; find interval
    do i = 0, nlev-2
      if (value.ge.levels(i) .and. value.lt.levels(i+1)) then
        return(cmap(i+1,:))   ; RETURN RGBA vector
      end if
    end do

    return(cmap(nlev-1,:))
  
end


begin




   pop_pi_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_piControl.p875.001.pop.h.4100_4199_timmean.nc","r")


    pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_LGM21ka_frUsp.p875.001.pop.h.5400_5499_timmean.nc","r")

   
    pop_pi = addfile ("/mnt/d/lgm_data/sselect_pi_idemix_4900_4604_4703_timmean.nc", "r")   

    pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5800_5899_timmean.nc", "r")
   

    pop_nc_example = addfile ("/mnt/d/lgm_data/kvmix_pi_tamas_4190_4199_timmean.nc", "r")

    pi_atlantic_mask = addfile ("/mnt/d/lgm_data/pi_atl_mask_miss.nc", "r")

    lgm_atlantic_mask = addfile ("/mnt/d/lgm_data/lgm_atl_mask_miss.nc", "r")

    atlantic_mask = addfile ("/mnt/c/Users/pilat/Downloads/d13c_atlantic_mask.nc", "r")    


   wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"

;   wks_type@wkOrientation = "landscape"
  
   wks  = gsn_open_wks(wks_type,"d13c_proxy_pi_lgm_last100years_thesis_30degN_noscale_v2_-3and3")
 


  ; cmap = read_colormap_file("MPL_RdYlGn")
   
   

   ;ncmap = dimsizes(cmap(:,0))
 


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


;mmol/m^3

nlt    = 9 

;mc1         = pop_pi_tamas->CO3(0,:,:,nlt)

;mc2         = pop_lgm_tamas->CO3(0,:,:,nlt)

;mc3         = pop_pi->CO3(0,:,:,nlt)

;mc4         = pop_lgm->CO3(0,:,:,nlt)


;mc1         = pop_pi_tamas->AOU(0,:,:,nlt)

;mc2         = pop_lgm_tamas->AOU(0,:,:,nlt)

;mc3         = pop_pi->AOU(0,:,:,nlt)

;mc4         = pop_lgm->AOU(0,:,:,nlt)


mc1         = pop_pi_tamas->CISO_DIC_d13C(0,:,:,:)

mc2         = pop_lgm_tamas->CISO_DIC_d13C(0,:,:,:)

mc3         = pop_pi->CISO_DIC_d13C(0,:,:,:)

mc4         = pop_lgm->CISO_DIC_d13C(0,:,:,:)


;mc1         = pop_pi_tamas->PO4(0,:,:,nlt)

;mc2         = pop_lgm_tamas->PO4(0,:,:,nlt)

;mc3         = pop_pi->PO4(0,:,:,nlt)

;mc4         = pop_lgm->PO4(0,:,:,nlt)


 atl_mask = pi_atlantic_mask->mask

printVarSummary(mc1)
printVarSummary(atl_mask)
 
mc1_mask = mask(mc1,atl_mask,1)
mc2_mask = mask(mc2,atl_mask,1)
mc3_mask = mask(mc3,atl_mask,1)
mc4_mask = mask(mc4,atl_mask,1)



lat_a    = pop_lgm->TLAT(:,nlt)
lat_atlantic        = lat_a



z_t           = pop_lgm->z_t
z_t           = z_t*1.E-2       ; convert from cm to m
;z_t@long_name = "Depth (T grid)"; overwrite long name
z_t@long_name = ""; overwrite long name
z_t@units     = "m"             ; fix units attribute
z_t@tiYAxisPosition = "top"
z_t@tiYAxisFontHeightF = 0.018

mc1_mask!0    = "z_t"
mc1_mask&z_t  =  z_t
mc1_mask!1    = "lat_t"
mc1_mask&lat_t=  lat_atlantic


mc2_mask!0    = "z_t"
mc2_mask&z_t  =  z_t
mc2_mask!1    = "lat_t"
mc2_mask&lat_t=  lat_atlantic



mc3_mask!0    = "z_t"
mc3_mask&z_t  =  z_t
mc3_mask!1    = "lat_t"
mc3_mask&lat_t=  lat_atlantic



mc4_mask!0    = "z_t"
mc4_mask&z_t  =  z_t
mc4_mask!1    = "lat_t"
mc4_mask&lat_t=  lat_atlantic



; proxy_data = asciiread("/mnt/c/Users/pilat/Downloads/ts01_excel2.txt",-1,"string")
 
 proxy_data = asciiread("/mnt/c/Users/pilat/Downloads/ts01_excel2.txt",(/401,6/),"float")

 printVarSummary(proxy_data)


 lat1 = proxy_data(:,0)
 lon1 = proxy_data(:,1)
 depth1 = proxy_data(:,2) ; m
 hol1 = proxy_data(:,3)
 lgm1 = proxy_data(:,4) 
 prx1 = proxy_data(:,5)

 ;print(proxy_data(:,5)) 

  print(prx1)
 
  print(lat1)


ii2 =  ind(lat1.ge.-40 .and. lat1.le.60)

 ;if (lat1.ge.-40 .and. lat1.le.70) then
 ;     lat1_v2 = ind(lat1)
 ;end if

;ii =  ind(lat1.ge.-40)
;ii2 = ind(ii.le.70)




lat_sel   = lat1(ii2) 
depth_sel   = depth1(ii2)
proxy_sel   = prx1(ii2)
hol_sel = hol1(ii2) 
lgm_sel = lgm1(ii2) 

print(lat_sel)
print(depth_sel)


;idx = qsort(depth_sel)

;lat_sorted   = lat_sel(idx)
;lon_sorted   = lon_sel(idx)
;depth_sorted = depth_sel(idx)    ; shallow → deep
;d13_sorted   = d13_sel(idx)


;proxy_data@lon1d = lon1
;proxy_data@lat1d = lat1

 mkres               = True
  mkres@gsMarkerIndex = 16 
  mkres@gsMarkerSizeF = 5.0 
  mkres@gsMarkerColor  = "black"



;lines = asciiread("/mnt/c/Users/pilat/Downloads/ts01_excel2.txt",-1,"string")
;
; Use "str_get_field" to indicate which fields to read in.
; Each field is separated by an arbitrary number of spaces.
;
;  stano       = tofloat(str_get_field(lines(1:),5," "))
 ; stano@lon1d = tofloat(str_get_field(lines(1:),1," "))  ; Attach as lat1d, lon1d
 ; stano@lat1d = tofloat(str_get_field(lines(1:),0," "))  ; for plotting later

;printVarSummary(stano)
;print(stano)




res                      = True       ; plot mods desired
res@gsnDraw              = False                ; don't draw yet
res@gsnFrame             = False                ; don't advance frame yet
res@lbLabelBarOn        = True
res@cnFillOn             = True       ; turn on color fill
;res@cnFillPalette        = "cmocean_balance" ; set color map
;res@cnFillPalette        = cmap(::-1,:)


;res@cnFillPalette        = cmap
res@cnFillPalette        = "cmp_b2r"
res@cnLinesOn            = False      ; turn off contour lines
res@cnLinesOn            = True      ; turn off contour lines
res@cnLineLabelsOn       = True      ; turns off contour line labels
res@cnMinLevelValF       = -2
res@cnMaxLevelValF       =  2
res@cnLevelSpacingF      =  0.4
res@cnLabelMasking        = True   
res@cnLineLabelDensityF = 0.8
res@cnLineLabelAngleF = 0.0
;res@cnLineLabelPlacementMode = "computed"
;res@cnLineLabelPlacementMode = "constant"
res@cnLineLabelFontColor = "black"
res@cnLineLabelFontThicknessF = 1.3
res@cnLineLabelBackgroundColor = "white"
res@cnLineLabelFontHeightF = 0.016
res@gsnContourZeroLineThicknessF = 1.3       ; eliminates zero contour
res@gsnContourNegLineDashPattern = 1       ; sets negative contours to dash pattern 1
res@lbOrientation        = "horizontal" ; vertical label bar
res@gsnYAxisIrregular2Linear = True 
res@trYReverse           = True       ; reverses y-axis
res@cnMissingValFillColor = "gray"
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/ -8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24 /)   ; set levels

res@cnLevelSelectionMode = "ManualLevels" ; manual contour levels
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

res@tmXBTickSpacingF = 20.0
res@cnMissingValFillColor = "grey52"

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
;res1@cnFillPalette        = cmap

 
res1@cnFillPalette        = "cmp_b2r"
res1@cnLinesOn            = True      ; turn off contour lines
res1@cnLineLabelsOn       = False      ; turns off contour line labels
;res1@cnMinLevelValF       = -16
;res1@cnMaxLevelValF       =  16
;res1@cnLevelSpacingF      =  2 
res1@cnLineLabelAngleF = 0.0
;res1@cnLineLabelPlacementMode = "computed"
;res1@cnLineLabelPlacementMode = "constant"
;res1@cnLineLabelPlacementMode = "randomized"
res1@cnLineLabelFontColor = "grey"
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

   ;res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   ;res1@cnMinLevelValF            =   40             ; min level
   ;res1@cnMaxLevelValF            =  200             ; max level
   ;res1@cnLevelSpacingF           =  20            ; contour interval


;res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;   res1@cnMinLevelValF            =   40             ; min level
;   res1@cnMaxLevelValF            =  300             ; max level
;   res1@cnLevelSpacingF           =  40            ; contour interval



   res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   res1@cnMinLevelValF            =   -2             ; min level
   res1@cnMaxLevelValF            =  2             ; max level
  res1@cnLevelSpacingF           =  0.2            ; contour interval

  ; res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  ; res1@cnMinLevelValF            =  0             ; min level
  ; res1@cnMaxLevelValF            =  120             ; max level
  ; res1@cnLevelSpacingF           =  20            ; contour interval

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


  ; res2@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  ; res2@cnMinLevelValF            =  0             ; min level
  ; res2@cnMaxLevelValF            =  120             ; max level
  ; res2@cnLevelSpacingF           =  20            ; contour interval


;res5@tmXBLabelFontHeightF =  0.022
;res5@tmYLLabelFontHeightF =  0.022

res5 = res
res5                      = True       ; plot mods desired
res5@gsnDraw              = False                ; don't draw yet
res5@gsnFrame             = False                ; don't advance frame yet
res5@lbLabelBarOn        = True
res5@cnFillOn             = True       ; turn on color fill

;res5@cnFillPalette        = "cmocean_balance"

;res5@cnFillPalette        = cmap
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


res5@cnLineLabelsOn       = False

;!!!res5@cnLevelSelectionMode = "ExplicitLevels"
;res5@cnLevels    = (/ -1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5 /)   ; set levels
;!!!res5@cnLevels    = (/ -2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5,1.8,2.1 /)   ; set levels


   ;res5@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   ;res5@cnMinLevelValF            = -50             ; min level
   ;res5@cnMaxLevelValF            =  50             ; max level
   ;res5@cnLevelSpacingF           =  10            ; contour interval

 ;res5@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
 ;  res5@cnMinLevelValF            = -80             ; min level
 ;  res5@cnMaxLevelValF            =  80             ; max level
 ;  res5@cnLevelSpacingF           =  20            ; contour interval

  
  res5@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   res5@cnMinLevelValF            = -3             ; min level
   res5@cnMaxLevelValF            =  3             ; max level
  res5@cnLevelSpacingF           =  0.5            ; contour interval


  ; res5@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  ; res5@cnMinLevelValF            = -1             ; min level
  ; res5@cnMaxLevelValF            =  1             ; max level
  ; res5@cnLevelSpacingF           =  0.2            ; contour interval

;   res5@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;   res5@cnMinLevelValF            = -2             ; min level
;   res5@cnMaxLevelValF            =  2             ; max level
;   res5@cnLevelSpacingF           =  0.5            ; contour interval
  

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
;res5@cnLineLabelPlacementMode = "constant"
;res5@cnLineLabelPlacementMode = "randomized"
res5@cnLineColor = "gray22"
res5@cnLineLabelFontColor = "black"
res5@cnLineLabelFontThicknessF = 1.3
res5@cnLineLabelBackgroundColor = "white"
res5@cnLabelMasking        = True
res5@cnLineLabelFontHeightF = 0.016
;res5@cnLineLabelInterval = 2

res5@cnLineLabelFormat    = "0@;*.2f"

res5@vpWidthF           = 0.6 ;longer graphs
res5@vpHeightF          = 0.3
;res5@cnLevelSpacingF = 2.0
res5@cnLineLabelDensityF = 2.0



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


res4@cnLineLabelDensityF = 2.0

res4@cnLineLabelFormat    = "0@;*.2f"



  ; res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  ; res4@cnMinLevelValF            = -50             ; min level
  ; res4@cnMaxLevelValF            =  50             ; max level
  ; res4@cnLevelSpacingF           =  10            ; contour interval

   ; res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   ; res4@cnMinLevelValF            = -80             ; min level
   ; res4@cnMaxLevelValF            =  80             ; max level
   ; res4@cnLevelSpacingF           =  20            ; contour interval

   res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   res4@cnMinLevelValF            = -3             ; min level
   res4@cnMaxLevelValF            =  3             ; max level
   res4@cnLevelSpacingF           =  0.5            ; contour interval



;res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;   res4@cnMinLevelValF            = -1             ; min level
;   res4@cnMaxLevelValF            =  1             ; max level
;  res4@cnLevelSpacingF           =  0.2            ; contour interval



;   res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;   res4@cnMinLevelValF            = -2             ; min level
;   res4@cnMaxLevelValF            =  2             ; max level
;   res4@cnLevelSpacingF           =  0.5            ; contour interval


res4@vpWidthF           = 0.6 ;longer graphs
res4@vpHeightF          = 0.3


print(max(mc1_mask))
print(min(mc1_mask))

print(max(mc2_mask))
print(min(mc2_mask))


print(max(mc3_mask))
print(min(mc3_mask))


print(max(mc4_mask))
print(min(mc4_mask))




;amoc1_v2 = abs(amoc1)
;amoc2_v2 = abs(amoc2)
;amoc4_v2 = abs(amoc4)
;amoc6_v2 = abs(amoc6)



diff1 = mc2_mask - mc1_mask

diff2 = mc4_mask - mc3_mask

;printVarSummary(diff1)
printVarSummary(diff2)

print(max(diff1))
print(min(diff1))

print(max(diff2))
print(min(diff2))


diff1!0    = "z_t"
diff1&z_t  = z_t  
diff1!1    = "lat_t"
diff1&lat_t= lat_atlantic


diff2!0    = "z_t"
diff2&z_t  = z_t  
diff2!1    = "lat_t"
diff2&lat_t= lat_atlantic

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

plot(0) = gsn_csm_contour(wks,mc1_mask(:,:,nlt),res2)     ; create plot
plot(1) = gsn_csm_contour(wks,mc3_mask(:,:,nlt),res1)     ; create plot
plot(2) = gsn_csm_contour(wks,mc2_mask(:,:,nlt),res2)     ; create plot
plot(3) = gsn_csm_contour(wks,mc4_mask(:,:,nlt),res1)     ; create plot
plot(4) = gsn_csm_contour(wks,diff1(:,:,nlt),res4)     ; create plot
plot(5) = gsn_csm_contour(wks,diff2(:,:,nlt),res5)     ; create plot


printVarSummary(diff1)
printVarSummary(diff2)

print(max(diff1))
print(min(diff1))

print(max(diff2))
print(min(diff2))



pres = True
pres@gsMarkerIndex = 16
pres@gsMarkerSizeF = 5.0


cmap = read_colormap_file("cmp_b2r")   ; shape (ncolors,4)


npts = dimsizes(proxy_sel)

ncmap = dimsizes(cmap(:,0))

marker_rgba = new((/npts,4/), float)


dmin = min(proxy_sel)
dmax = max(proxy_sel)


mnmx = nice_mnmxintvl(-3, 3, ncmap, False)

;mnmx = nice_mnmxintvl(-1, 1, ncmap, False)

;mnmx = nice_mnmxintvl(-2, 2, ncmap, False)

print(mnmx(0))
print(mnmx(1))
print(mnmx(2))

minlev = mnmx(0)
maxlev = mnmx(1)
dlev   = mnmx(2)

nlev = toint( 1 + (maxlev - minlev) / dlev )

levels = fspan(minlev, maxlev, nlev)
print(levels)



renk = span_color_rgba (cmap,ncmap)
do i = 0, npts-1
  marker_rgba(i,:) = get_color_rgba1(proxy_sel(i), levels, renk)
end do

;print(marker_rgba)

printVarSummary(marker_rgba)

  dum1 = new(npts,graphic)
  dum2 = new(npts,graphic)  

do i = 0, npts-1
 pres@gsMarkerColor = marker_rgba(i,:)
  dum1(i) = gsn_add_polymarker(wks, plot(4), lat_sel(i), depth_sel(i), pres)
end do


do i = 0, npts-1
 pres@gsMarkerColor = marker_rgba(i,:)
  dum2(i) = gsn_add_polymarker(wks, plot(5), lat_sel(i), depth_sel(i), pres)
end do


print(proxy_sel)


;########### pi proxie ###########

pres1 = True
pres1@gsMarkerIndex = 16
pres1@gsMarkerSizeF = 5.0


cmap1 = read_colormap_file("cmp_b2r")   ; shape (ncolors,4)


npts1 = dimsizes(hol_sel)

ncmap1 = dimsizes(cmap1(:,0))

marker_rgba1 = new((/npts1,4/), float)


dmin1 = min(hol_sel)
dmax1 = max(hol_sel)


mnmx1 = nice_mnmxintvl(-2, 2, ncmap, False)

print(mnmx1(0))
print(mnmx1(1))
print(mnmx1(2))

minlev1 = mnmx1(0)
maxlev1 = mnmx1(1)
dlev1   = mnmx1(2)

nlev1 = toint( 1 + (maxlev1 - minlev1) / dlev1 )

levels1 = fspan(minlev1, maxlev1, nlev1)
;print(levels1)



renk1 = span_color_rgba (cmap,ncmap)
do i = 0, npts1-1
  marker_rgba1(i,:) = get_color_rgba1(hol_sel(i), levels1, renk1)
end do

;print(marker_rgba)

printVarSummary(marker_rgba1)

  dum11 = new(npts1,graphic)
  dum21 = new(npts1,graphic)  

do i = 0, npts1-1
 pres1@gsMarkerColor = marker_rgba1(i,:)
  dum11(i) = gsn_add_polymarker(wks, plot(0), lat_sel(i), depth_sel(i), pres1)
end do


do i = 0, npts1-1
 pres1@gsMarkerColor = marker_rgba1(i,:)
  dum21(i) = gsn_add_polymarker(wks, plot(1), lat_sel(i), depth_sel(i), pres1)
end do




;########### lgm proxie ###########

pres2 = True
pres2@gsMarkerIndex = 16
pres2@gsMarkerSizeF = 5.0


cmap2 = read_colormap_file("cmp_b2r")   ; shape (ncolors,4)


npts2 = dimsizes(lgm_sel)

ncmap2 = dimsizes(cmap2(:,0))

marker_rgba2 = new((/npts2,4/), float)


dmin2 = min(lgm_sel)
dmax2 = max(lgm_sel)


mnmx2 = nice_mnmxintvl(-2, 2, ncmap, False)

print(mnmx2(0))
print(mnmx2(1))
print(mnmx2(2))

minlev2 = mnmx2(0)
maxlev2 = mnmx2(1)
dlev2   = mnmx2(2)

nlev2 = toint( 1 + (maxlev2 - minlev2) / dlev2 )

levels2 = fspan(minlev2, maxlev2, nlev2)
;print(levels2)



renk2 = span_color_rgba (cmap,ncmap)
do i = 0, npts2-1
  marker_rgba2(i,:) = get_color_rgba1(lgm_sel(i), levels2, renk2)
end do

;print(marker_rgba2)

printVarSummary(marker_rgba2)

  dum12 = new(npts2,graphic)
  dum22 = new(npts2,graphic)  

do i = 0, npts2-1
 pres2@gsMarkerColor = marker_rgba2(i,:)
  dum12(i) = gsn_add_polymarker(wks, plot(2), lat_sel(i), depth_sel(i), pres2)
end do


do i = 0, npts2-1
 pres2@gsMarkerColor = marker_rgba2(i,:)
  dum22(i) = gsn_add_polymarker(wks, plot(3), lat_sel(i), depth_sel(i), pres2)
end do




resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = False
resP@pmLabelBarWidthF   = 0.07
resP@pmLabelBarHeightF = 0.33
resP@lbOrientation   = "vertical"
;resP@pmLabelBarOrthogonalPosF = -0.006
resP@lbTitleString = "[mmol/m^3]"
resP@lbTitleFontHeightF = 0.012
resP@gsnMaximize         = True                 ; make as large as possible



resP@pmLabelBarParallelPosF = 0.05
resP@lbTitleExtentF = 0.18
resP@lbTitleOffsetF = 0.04
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
;resP@gsnPanelYWhiteSpacePercent = 0.2
;resP@gsnPanelXWhiteSpacePercent = 0.4


resP@gsnPanelFigureStrings= (/"a) PI", "d) PI", "b) LGM", "e) LGM", "c) LGM - PI", "f) LGM - PI"/)
;resP@gsnPanelFigureStrings= (/"W", "D", "W", "D", "W", "D"/)
resP@gsnPanelFigureStringsFontHeightF = 0.011
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
;resPs@gsnPanelFigureStrings= (/"PI - with IDEMIX", "LGM - with IDEMIX", "AMOC - LGM-j10"/)
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
txres1@txFontHeightF = 0.013
txres1@txFont   = "helvetica-bold"

;gsn_text_ndc(wks,"MOC",0.12,0.69,txres1) ; to plot normal data
;gsn_text_ndc(wks,"AMOC, w/IDEMIX - w/oIDEMIX",0.24,0.41,txres1)
gsn_text_ndc(wks,"w/oIDEMIX",0.29,0.98,txres1)
gsn_text_ndc(wks,"w/IDEMIX",0.69,0.98,txres1)

;Carbonate Ion Concentration
;Apparent O2 Utilization

;gsn_text_ndc(wks,"AMOC, w/IDEMIX - w/oIDEMIX",0.21,0.43,txres1)
;gsn_text_ndc(wks,"AMOC, w/IDEMIX",0.15,0.83,txres1)
;gsn_text_ndc(wks,"AMOC, w/oIDEMIX",0.15,0.67,txres1)


txres4 = True
txres4@txFontHeightF = 0.0115
txres4@txAngleF      = 270.     ; rotate 90 degrees



;gsn_text_ndc(wks,"[mmol/m^3]",0.95,0.83,txres4)
;gsn_text_ndc(wks,"[mmol/m^3]",0.95,0.56,txres4)
;gsn_text_ndc(wks,"[mmol/m^3]",0.95,0.29,txres4)

gsn_text_ndc(wks,"[permil]",0.95,0.83,txres4)
gsn_text_ndc(wks,"[permil]",0.95,0.56,txres4)
gsn_text_ndc(wks,"[permil]",0.95,0.29,txres4)
 
;  txres5 = True

;  txid3 = gsn_create_text(wks, "Depth (m)", txres5)


txres5 = True
txres5@txFontHeightF = 0.012
txres5@txAngleF      = 90.     ; rotate 90 degrees
gsn_text_ndc(wks,"Depth (m)",0.04,0.54,txres5)
;gsn_text_ndc(wks,"Depth (m)",0.032,0.31,txres5)





resP@gsnPanelYF = (/0.92,0.92,0.65,0.65,0.38,0.38/)
;resP@gsnPanelXF = (/0.07,0.51,0.07,0.51,0.07,0.51/)

resP@gsnPanelXF = (/0.12,0.51,0.12,0.51,0.12,0.51/)

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


;draw(plot)
  frame(wks)   ; Now advance the frame.

end
