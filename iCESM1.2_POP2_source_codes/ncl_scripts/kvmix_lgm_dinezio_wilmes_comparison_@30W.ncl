
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



    pop = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4694_4703_timmean.nc", "r")

;    pop_pi_tamas = addfile ("/mnt/d/lgm_data/kvmix_pi_tamas_4190_4199_timmean.nc", "r")

;    pop_pi_tamas = addfile ("/mnt/d/lgm_data/kvmix_lgm_tamas_5490_5499_timmean.nc", "r")

   pop_lgm_tamas = addfile ("/mnt/d/lgm_data/kvmix_lgm_tamas_5490_5499_timmean.nc", "r")


     pop_pi_tamas = addfile ("/mnt/d/lgm_data2/kvmix_test_cmake_1790_1799_timmean.nc", "r")

     pop_lgm = addfile ("/mnt/d/lgm_data2/dinezio_continued_kvmix_1850_1859.nc", "r")  


;    pop_pi_tamas = addfile ("/mnt/d/lgm_data/se 
                                                                                                                                                    ;    pop_lgm_tamas = addfile ("/mnt/d/lgm_data/select_lgm_year1000_idemix_isotope_continued_2229_2338.nc", "r")

   ;pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4694_4703_timmean.nc", "r")

   ;pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5800_5899.nc", "r")

   pop_pi = addfile ("/mnt/d/lgm_data/kvmix_lgm_idemix_isotope_continued_5900_5909_timmean.nc", "r")


     


   wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"
   
    wks  = gsn_open_wks(wks_type,"kvmix_pi_lgm_thesis_@30W_lastnew_dinezio") 

;   wks  = gsn_open_wks(wks_type,"kvmix_lgm_dinezio_wilmes_comparison_@30W_v4")
  
;   wks  = gsn_open_wks("png","kvmix_pi_lgm_idemix_latest_2023_fr")

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

pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
;plot = new(2,graphic)
plots = new(4,graphic)

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

res@cnMissingValFillPattern = "SolidFill"


res@cnMissingValFillColor = "grey56"
;res@cnFillPalette        = read_colormap_file("MPL_rainbow")
res@vpWidthF           = 0.13
res@vpHeightF          = 0.18
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/-13,-12,-11,-10,-9,-8,-7,-6/)
;res@cnLevels    = (/-13,-12.5,-12,-11,-11.5,-10,-10.5,-9,-9.5,-8,-8.5,-7,-7.5,-6,-6.5,-5/)   ; set levels
res@gsnLeftString =  ""
res@gsnRightString =  ""



res1 = True
res1@vpWidthF           = 0.5
res1@vpHeightF          = 0.23
  res1@gsnDraw              = False                ; don't draw yet
  res1@gsnFrame             = False                ; don't advance frame
res1@Footer = "  "
res1@tmXBLabelFontHeightF = 0.011
res1@tmYLLabelFontHeightF = 0.011
res1@cnMissingValFillPattern = "SolidFill"
res1@lbLabelBarOn        = False
res1@gsnYAxisIrregular2Linear = True
res1@trYReverse           = True       ; reverses y-axis
res1@cnMissingValFillColor = "grey56"
res1@cnInfoLabelOn   = False
res1@tmXBLabelsOn  = True
res1@cnLinesOn            = False                ; turn off contour lines
res1@cnLineLabelsOn       = False                ; turn off contour line
res1@cnFillOn             = True                 ; turn on color fill
;res1@cnFillPalette        = read_colormap_file("MPL_rainbow")         ; set color map
res1@cnFillPalette        = cmap(::-1,:)
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/-13,-12,-11,-10,-9,-8,-7,-6/)
;res1@tiYAxisString = "Depth (m)"
;res1@cnLevels    = (/-13,-12.5,-12,-11,-11.5,-10,-10.5,-9,-9.5,-8,-8.5,-7,-7.5,-6,-6.5,-5/)  

;res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;res1@cnMinLevelValF            = -10             ; min level
;res1@cnMaxLevelValF            = -4             ; max level
;res1@cnLevelSpacingF           =  0.5            ; contour interval



res2 = True
res2@vpWidthF           = 0.35
res2@vpHeightF          = 0.25
  res2@gsnDraw              = False                ; don't draw yet
  res2@gsnFrame             = False                ; don't advance frame yet
res2@lbLabelBarOn        = False

;res2@tmXBLabelFontHeightF = 0.007
;res2@tmYLLabelFontHeightF = 0.007
res2@cnMissingValFillPattern = "SolidFill"

res2@cnMissingValFillColor = "grey56"
res2@cnLinesOn            = False                ; turn off contour lines
res2@cnLineLabelsOn       = False                ; turn off contour line
res2@cnFillOn             = True                 ; turn on color fill
;res2@cnFillPalette        = cmap(:ncmap-1,:)         ; set color map
res2@cnFillPalette        = cmap(:ncmap-1,:)
;res2@cnLevelSelectionMode = "ExplicitLevels"
;res2@cnLevels    = (/ -10.5,-10,-9.5,-9,-8.5,-8,-7.5,-7,-6.5,-6 /)   ; set levels


res3 = res1
res3@lbLabelBarOn        = True
  
  iw1        = pop_pi_tamas->KVMIX
;mean_pi_tamas = iw1(0,43,:,:)
;mean_pi_tamas_m2 = mean_pi_tamas * 0.0001
;printVarSummary(mean_pi_tamas_m2)
;mean_pi_tamas_log = log10(mean_pi_tamas_m2)
;printVarSummary(mean_pi_tamas_log)



;iw1        = pop_pi_tamas->KVMIX(80:99,:,:,:) 
 
;;mean_temp2 = dim_avg_n_Wrap(iw1, (/0,1/))
;mean_pi_tamas = iw1(:,43,:,:)
;mean_temp1 = dim_avg_n(mean_pi_tamas, 0)
;mean_pi_tamas_m2 = mean_temp1 * 0.0001
;printVarSummary(mean_pi_tamas_m2)
;mean_pi_tamas_log = log10(mean_pi_tamas_m2)
;printVarSummary(mean_pi_tamas_log)


  iw2        = pop_lgm_tamas->KVMIX

mean_lgm_tamas = iw2(0,43,:,:)
;;mean_temp2 = dim_avg_n_Wrap(iw1, (/0,1/))

mean_lgm_tamas_m2 = mean_lgm_tamas * 0.0001
printVarSummary(mean_lgm_tamas_m2)
mean_lgm_tamas_log = log10(mean_lgm_tamas_m2)
printVarSummary(mean_lgm_tamas_log)



;iw2        = pop_lgm_tamas->KVMIX(80:99,:,:,:)
;mean_lgm_tamas = iw2(:,43,:,:)
;mean_temp2 = dim_avg_n(mean_lgm_tamas, 0)
;mean_lgm_tamas_m2 = mean_temp2 * 0.0001
;printVarSummary(mean_lgm_tamas_m2)
;mean_lgm_tamas_log = log10(mean_lgm_tamas_m2)
;printVarSummary(mean_lgm_tamas_log)



;!  iw_pi        = pop_pi->KVMIX(80:99,:,:,:)

iw_pi        = pop_pi->KVMIX(:,:,:,:)

mean_pi = iw_pi(:,43,:,:)
;mean_temp2_int = dim_avg_n_Wrap(iw_int, (/0,1/))
mean_pi_time = dim_avg_n(iw_pi, 0)
mean_pi_m2 = mean_pi_time * 0.0001
print(max(mean_pi_m2))
print(min(mean_pi_m2))
printVarSummary(mean_pi_time)
mean_pi_log = log10(mean_pi_m2)





  iw_lgm        = pop_lgm->KVMIX(:,:,:,:)


mean_lgm = iw_lgm(:,43,:,:)
;mean_temp2_int = dim_avg_n_Wrap(iw_int, (/0,1/))
mean_lgm_time = dim_avg_n(iw_lgm, 0)
mean_lgm_m2 = mean_lgm_time * 0.0001
print(max(mean_lgm_m2))
print(min(mean_lgm_m2))
printVarSummary(mean_lgm_time)
mean_lgm_log = log10(mean_lgm_m2)




ntm    = 0                  ; time subscript
nlt    = 9               ; latitude subscript
;nlt = 10

ntm2   = 0                  ; time subscript
nlt2    = 224               ; latitude subscript

nlt3    = 305 
nlt4 = 312  ; along 50N


lat    = pop_lgm->TLAT(nlt3,:)
latavg = avg(lat)
print(latavg)

lon_o  = pop_lgm->TLONG(nlt3,:)
print(lon_o)


nlon   = dimsizes (lon_o)
lon_t        = lon_o          ; trick to conserve meta data
lon_t@units  = "degrees_east"
lon_t( 0:283) = lon_o(36:319)
lon_t(284:319) = lon_o( 0:35)
ii           = ind(lon_o.GT. 320.0)
lon_t(ii)    = 360.-lon_o(ii)
lonavg = lon_t(nlt-1)

print(lon_t)

print("heves")
printVarSummary(iw1)


kvmix1         = (dim_avg_n(iw1, 0))*0.0001
;printVarSummary(kvmix1)
kvmix12 = kvmix1(:,:,nlt)
;kvmix12 = kvmix1(:,nlt4,:)
;kvmix12 = log(kvmix1)
;kvmix12_atlantic_log = kvmix12(:,:,nlt3)
kvmix12_atlantic_log = log(kvmix12)


lat_a    = pop_lgm->TLAT(:,nlt)
lat_atlantic        = lat_a
;lat_a    = pop_lgm->TLAT(nlt4,:)
;lat_atlantic        = lat_a
lat_p    = pop_lgm->TLAT(:,nlt2)
lat_pacific        = lat_p

kvmix2         = (dim_avg_n(iw2, 0))*0.0001
kvmix22 = kvmix2(:,:,nlt)
;kvmix22 = kvmix2(:,nlt4,:)
;kvmix22 = log(kvmix2)
kvmix22_atlantic_log = log(kvmix22)

kvmix3         = (dim_avg_n(iw_pi, 0))*0.0001
;printVarSummary(kvmix3)
kvmix32 = kvmix3(:,:,nlt)
;kvmix32 = kvmix3(:,nlt4,:)
;kvmix32 = log(kvmix3)
;kvmix32_atlantic_log = kvmix32(:,:,nlt)
kvmix32_atlantic_log = log(kvmix32)


kvmix4         = (dim_avg_n(iw_lgm, 0))*0.0001
kvmix42 = kvmix4(:,:,nlt)
;kvmix42 = kvmix4(:,nlt4,:)
;kvmix42 = log(kvmix4)
kvmix42_atlantic_log = log(kvmix42)





lat_aw    = pop_lgm->TLAT(:,nlt)
lat_atlanticw        = lat_aw
lat_pw    = pop_lgm->TLAT(:,nlt2)
lat_pacificw        = lat_pw


printVarSummary(kvmix12_atlantic_log)



z_t           = pop->z_t
z_t           = z_t*1.E-2       ; convert from cm to m
;z_t@long_name = "Depth (T grid)"; overwrite long name
z_t@long_name = ""; overwrite long name
z_t@units     = "m"             ; fix units attribute
z_t@tiYAxisPosition = "top"
z_t@tiYAxisFontHeightF = 0.018
kvmix12_atlantic_log!0    = "z_t"
kvmix12_atlantic_log&z_t  =  z_t
kvmix12_atlantic_log!1    = "lat_t"
kvmix12_atlantic_log&lat_t=  lat_atlantic
kvmix22_atlantic_log!0    = "z_t"
kvmix22_atlantic_log&z_t  =  z_t
kvmix22_atlantic_log!1    = "lat_t"
kvmix22_atlantic_log&lat_t=  lat_atlantic
kvmix32_atlantic_log!0    = "z_t"
kvmix32_atlantic_log&z_t  =  z_t
kvmix32_atlantic_log!1    = "lat_t"
kvmix32_atlantic_log&lat_t=  lat_atlantic
kvmix42_atlantic_log!0    = "z_t"
kvmix42_atlantic_log&z_t  =  z_t
kvmix42_atlantic_log!1    = "lat_t"
kvmix42_atlantic_log&lat_t=  lat_atlantic




;mean_temp2_log@lat2d   = ref->TLAT
;mean_temp2_log@lon2d   = ref->TLONG


;iw_temp2_log@lat2d   = ref->TLAT
;iw_temp2_log@lon2d   = ref->TLONG


;=================================================;
; Create plot


;  plot(0) = gsn_csm_contour_map(wks, iw1_log,res) 

  
;  plot(0) = gsn_csm_contour_map(wks, mean_temp2_log,res)

;  plot(1) = gsn_csm_contour_map(wks, iw_temp2_log,res)

;   plots(0) = gsn_csm_contour(wks, kvmix12_atlantic_log,res1)

;   plots(1) = gsn_csm_contour(wks, kvmix22_atlantic_log,res1)
 
  
   plots(0) = gsn_csm_contour(wks, kvmix12_atlantic_log,res1)

   plots(1) = gsn_csm_contour(wks, kvmix22_atlantic_log,res1)
   
   plots(2) = gsn_csm_contour(wks, kvmix32_atlantic_log,res1)

   plots(3) = gsn_csm_contour(wks, kvmix42_atlantic_log,res1)




print(max(kvmix12_atlantic_log))
print(min(kvmix12_atlantic_log))
print(max(kvmix22_atlantic_log))
print(min(kvmix22_atlantic_log))

 
resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = True
; resP@pmLabelBarWidthF   = 0.60
;resP@pmLabelBarHeightF = 0.048
resP@lbOrientation   = "Horizontal"   
resP@pmLabelBarOrthogonalPosF = -0.025
resP@lbTitleString = " Vertical Diffusivity {log10[K(m^2/s)]}"
resP@lbTitleFontHeightF = 0.011
resP@gsnMaximize         = True                 ; make as large as possible   					
resP@lbLabelFontHeightF= 0.010 
resP@gsnPanelFigureStringsFontHeightF = 0.007
resP@gsnPanelYWhiteSpacePercent = 3
resP@gsnPanelXWhiteSpacePercent = 1


resP@gsnPanelFigureStrings= (/"w/oIDEMIX", "w/oIDEMIX", "w/IDEMIX", "w/IDEMIX"/) 
;resP@gsnPanelFigureStrings= (/"yr1000", "yr1000", "4900", "5900"/)
resP@gsnPanelFigureStringsFontHeightF = 0.011
resP@cnLevelSelectionMode = "ExplicitLevels"
;resP@cnLevels    = (/ -7, -6, -5, -4, -3, -2 /)   ; set levels
;resP@cnLevels    = (/ -13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2 /)

	
;;resP@gsnPanelBottom = 0.9
;resP@gsnPanelTop = 0.92
resP@gsnPanelLeft = 0.02
resP@gsnLeftString =  "Tidal Energy Dissipation"
resP@gsnRightString           = "(W/m^2)"  

;; resP@gsnPanelTop       = 0.92
;;  resP@gsnPanelBottom    = 0.45
;resP@gsnPanelTop       = 0.92
;resP@gsnPanelBottom    = 0.34

;  resP@gsnPanelLabelBar  = True
  resP@pmLabelBarWidthF  = 0.6

txres1 = True
txres1@txFontHeightF = 0.011

;gsn_text_ndc(wks,"near 45 degW",0.12,0.78,txres1)



;txres2 = True
;txres2@txFontHeightF = 0.011

;gsn_text_ndc(wks,"Near 30.4 degW in Atlantic Ocean",0.16,0.45,txres2)


txres22 = True
txres22@txFontHeightF = 0.014
txres22@txFont   = "helvetica-bold"
;gsn_text_ndc(wks,"PI",0.28,0.8,txres22)
;gsn_text_ndc(wks,"LGM",0.78,0.8,txres22)
gsn_text_ndc(wks,"(Di-Nezio et al., 2016)",0.78,0.8,txres22)
gsn_text_ndc(wks,"(Wilmes et al., 2021)",0.28,0.8,txres22)

;txres3 = True
;txres3@txFontHeightF = 0.011

;gsn_text_ndc(wks,"near 45 degW",0.62,0.78,txres3)



;txres4 = True
;txres4@txFontHeightF = 0.011

;gsn_text_ndc(wks,"Near 30.4 degW in Atlantic Ocean",0.66,0.45,txres4)


txres5 = True
txres5@txFontHeightF = 0.012
  txres5@txAngleF      = 90.     ; rotate 90 degrees
;  txres5@txJust        = "TopLeft"

gsn_text_ndc(wks,"Depth (m)",0.02,0.55,txres5)



;pres  = True

;resP@pmLabelBarOrthogonalPosF = -0.02    ; Move closer to plots
;resPs@pmLabelBarOrthogonalPosF = -0.02

;; pres@gsnPanelTop    = 1.0  ; set these if desired
;; pres@gsnPanelBottom = 0.0  ; must be between 0 and 1
;; pres@gsnPanelLeft   = 0.0
;; pres@gsnPanelRight  = 1.0

;pres@gsnPanelRowSpec = False


;  panel_two_sets(wks,plot,plots,(/1,2/),(/1,2/),resP,resPs,pres)
 




;  gsn_panel(wks,plot,(/1,2/),resP)
;  gsn_panel(wks,plots,(/1,2/),resPs)

   gsn_panel(wks,plots,(/2,2/),resP)

 
; now we draw the whole thing and advance the frame.


  frame(wks)
  end










