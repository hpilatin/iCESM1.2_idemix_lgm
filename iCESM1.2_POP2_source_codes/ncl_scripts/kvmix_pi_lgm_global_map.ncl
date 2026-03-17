
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



;---Get list of files from years 5016 to 5026


  ; pop_pi_tamas = addfile ("/home/hbkheves/from_Tamas/piControl_rho_kvmix_q.p875.001.pop.h_4190_4199.deflated.nc", "r")
   
   ;pop_lgm_tamas = addfile ("/home/hbkheves/from_Tamas/lgm_rho_kvmix_q.p875.001.pop.h_5490_5499.deflated.nc", "r")
 
  ;---Get list of files from years 2850 to 3005 (PI)

 ;  yyyy      = ispan(3550,3559,1) ; Last 10 years
 ;  all_files = "/scratch/projects/hbk00090/b.ie12.Bi1850C5BPRDCISO.f19_g16_PI_ciso_j5_idemix/ocn/hist/b.ie12.Bi1850C5BPRDCISO.f19_g16_PI_ciso_j5_idemix.pop.h." + yyyy + "_deflated"+".nc"
  
 ; pop_pi = addfiles (all_files, "r")
  
  ;  yyyy2      = ispan(2060,2069,1) ; Last 10 years
  ; all_files2 = "/scratch/projects/hbk00090/b.ie12.Bi1850C5BPRDCISO.f19_g16_LGM21ka_ciso_idemix_j5_year1000/ocn/hist/b.ie12.Bi1850C5BPRDCISO.f19_g16_LGM21ka_ciso_idemix_j5_year1000.pop.h." + yyyy2 + "_deflated" +".nc"
   
  
  ;pop_lgm = addfiles (all_files2, "r")
      
;  pop_pi_tamas = addfile ("/mnt/d/lgm_data/kvmix_pi_tamas_4190_4199_timmean.nc", "r")                                                                                                                                                                                                                                                                      
pop_lgm_tamas = addfile ("/mnt/d/lgm_data/kvmix_lgm_tamas_5490_5499_timmean.nc", "r")                                                                                                                                                          
pop_pi_tamas = addfile ("/mnt/d/lgm_data2/kvmix_test_cmake_1790_1799_timmean.nc", "r")
      
;pop_pi = addfile ("/mnt/d/lgm_data/kvmix_pi_idemix_isotope_continued_5310_5319_timmean.nc", "r")

pop_pi = addfile ("/mnt/d/lgm_data/kvmix_lgm_idemix_isotope_continued_5900_5909_timmean.nc", "r")

pop_lgm = addfile ("/mnt/d/lgm_data/kvmix_dinezio_continued.pop.h.1850_1859_timmean.nc", "r")


    wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"



   wks  = gsn_open_wks(wks_type,"kvmix_pi_lgm_global_map_thesis_test_cmake_surface_lastnew_dinezio")

;  cmap = read_colormap_file("NCV_jet")       ; read color data
;   cmap = read_colormap_file("cmocean_balance")
    cmap = read_colormap_file("MPL_RdYlBu") 
  ncmap = dimsizes(cmap(:,0))


  res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
;  res@cnFillPalette        = cmap(:ncmap-10,:)         ; set color map
  res@cnFillPalette        = cmap(::-1,:)
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels

res@vpWidthF           = 0.14
res@vpHeightF          = 0.15
res@tmXBLabelsOn  = True


pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(2,graphic)
plots = new(2,graphic)

res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "
; res@gsnLeftString =  "Vertical Diffusivity"
; res@gsnRightString           = "(m^2/s)"
 res@gsnCenterString = "  "
res@tiMainString  = " "             ; add title
  ;res@gsnMaximize         = True                 ; make as large as possible   
res@tmXBLabelFontHeightF = 0.012
res@tmYLLabelFontHeightF = 0.012
res@lbLabelBarOn        = False

res@cnMissingValFillPattern = "SolidFill"
res@cnMissingValFillColor = "grey56"
 res@mpLandFillColor  = "grey56"


;;res@vpWidthF           = 0.15
;;res@vpHeightF          = 0.18
res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels    = (/-7,-6.5,-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2/)   ; set levels
;res@cnLevels    = (/-6.8,-6.4,-6,-5.6,-5.2,-4.8,-4.4,-4,-3.6,-3.2,-2.8,-2.4/)
res@cnLevels    = (/-6,-5.8,-5.6,-5.4,-5.2,-5,-4.8,-4.6,-4.4,-4.2,-4,-3.8,-3.6,-3.4,-3.2,-3,-2.8,-2.6,-2.4,-2.2,-2/)


;res@cnLevels    = (/-0.1,-0.01,0,0.01,0.1/)
res@gsnLeftString =  ""
res@gsnRightString =  ""

res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True

res@tmXBTickSpacingF = 60.0
res@tmYLTickSpacingF = 30.0

res@vpWidthF          = 0.48 ;longer graphs
res@vpHeightF          = 0.25





res1 = True
;res1@vpWidthF           = 0.15
;res1@vpHeightF          = 0.18
res1@vpWidthF           = 0.11
res1@vpHeightF          = 0.15
 
 res1@gsnDraw              = False                ; don't draw yet
  res1@gsnFrame             = False                ; don't advance frame yet

res1@tmXBLabelFontHeightF = 0.0020
res1@tmYLLabelFontHeightF = 0.0020
res1@cnMissingValFillPattern = "SolidFill"
res1@lbLabelBarOn        = False

res1@cnMissingValFillColor = "grey56"
 res1@mpLandFillColor  = "grey56"


res1@tmXBLabelsOn  = True
res1@cnLinesOn            = False                ; turn off contour lines
res1@cnLineLabelsOn       = False                ; turn off contour line
res1@cnFillOn             = True                 ; turn on color fill
;cmap2                     = read_colormap_file("MPL_BuPu")         ; set color map
;res1@cnFillPalette        = cmap2(::-1,:)
;res1@cnFillPalette        = cmap2(:ncmap-1,:)

;;cmap2                     = read_colormap_file("MPL_Purples")         ; set color map
;cmap2                     = read_colormap_file("BlueYellowRed")
;res1@cnFillPalette        = cmap2
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/-11,-10.5,-10,-9.5,-9,-8.5,-8/)   ; set levels
;res1@cnLevels    = (/-7,-6,-5,-4,-3,-2/)   ; set levels



res2 = True
res2@vpWidthF           = 0.35
res2@vpHeightF          = 0.24
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
;res2@cnLevels    = (/-7,-6,-5,-4,-3,-2/)   ; set levels

res3 = res

res3@tmYLLabelsOn  = False

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

;************************************************************
; now elimininate the 0th element in the latitudinal direction for the
; data. This means that lon and lat have one more element along both
; dimensions of the data array. This will automatically cause the
; discrete raster fill to treat X and Y coordianate arrays as cell bounds
;************************************************************


  lon1b   = pop_lgm->TLONG                   ; read in original lat/lon data
  lat1b   = pop_lgm->TLAT


; ================================================;
; when a coordinate array contains one more element than the corresponding
; dimension in the data array, the cell boundary is automatically turned on,
; which causes discrete rasterization to use the coordinate values as
; divisions between data cells, rather than to treat each coordinate value
; as a cell center.
; ================================================;
  dimsb      = dimsizes(lon1b)         ; get dimension sizes
  dimsb(1)   = dimsb(1) + 1               ; increase longitude by 1
  nlatb      = dimsb(0)
  nlonb      = dimsb(1)

  lonb      = new(dimsb,typeof(lon1b)) ; create new arrays
  latb      = new(dimsb,typeof(lat1b))

  latb(:,1:nlon-1)  = lat1b        ; fill in arrays
  lonb(:,1:nlon-1)  = lon1b
  lonb(:,0)         = lon1b(:,nlon-2)
  latb(:,0)         = lat1b(:,nlon-2)





  iw1        = pop_pi_tamas->KVMIX(:,:,1:nlat-1,:)

  

;temp1 = iw1(:,42,:,:)
;temp1 = iw1(:,39,:,:)
temp1 = iw1(:,0,:,:)
mean_temp1 = dim_avg_n(temp1, 0)

;mean_temp1 = dim_avg_n(temp1, 0)

mean_temp1_m2 = mean_temp1 * 0.0001

mean_temp1_log = log10(mean_temp1_m2)




  iw2        = pop_lgm_tamas->KVMIX(:,:,1:nlatb-1,:)


;temp2 = iw2(:,42,:,:)
;temp2 = iw2(:,39,:,:)
temp2 = iw2(:,0,:,:)
mean_temp2 = dim_avg_n(temp2, 0)

;mean_temp2 = dim_avg_n(temp2, 0)

mean_temp2_m2 = mean_temp2 * 0.0001

mean_temp2_log = log10(mean_temp2_m2)



iw3        = pop_pi->KVMIX(:,:,1:nlat-1,:)



;temp3 = iw3(50:99,42,:,:)
;temp3 = iw3(:,39,:,:)
temp3 = iw3(:,0,:,:)
mean_temp3 = dim_avg_n(temp3, 0)

;mean_temp3 = dim_avg_n(temp3, 0)

mean_temp3_m2 = mean_temp3 * 0.0001

mean_temp3_log = log10(mean_temp3_m2)




  iw4        = pop_lgm->KVMIX(:,:,1:nlatb-1,:)


;temp4 = iw4(50:99,42,:,:)
;temp4 = iw4(:,39,:,:)
temp4 = iw4(:,0,:,:)
mean_temp4 = dim_avg_n(temp4, 0)

;mean_temp4 = dim_avg_n(temp4, 0)

mean_temp4_m2 = mean_temp4 * 0.0001

mean_temp4_log = log10(mean_temp4_m2)



mean_temp1_log@lat2d   = lat
mean_temp1_log@lon2d   = lon


mean_temp3_log@lat2d   = lat
mean_temp3_log@lon2d   = lon


mean_temp2_log@lat2d   = latb
mean_temp2_log@lon2d   = lonb


mean_temp4_log@lat2d   = latb
mean_temp4_log@lon2d   = lonb




;=================================================;
; Create plot


res4 = res

res4@tmYLLabelsOn  = False
res4@tmXBLabelsOn  = True
res4@cnLevelSelectionMode = "ExplicitLevels"
;res4@cnLevels    = (/-7,-6.5,-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2/)   ; set levels
;res4@cnLevels    = (/-6.8,-6.4,-6,-5.6,-5.2,-4.8,-4.4,-4,-3.6,-3.2,-2.8,-2.4/)
res4@cnLevels    = (/-6,-5.8,-5.6,-5.4,-5.2,-5,-4.8,-4.6,-4.4,-4.2,-4,-3.8,-3.6,-3.4,-3.2,-3,-2.8,-2.6,-2.4,-2.2,-2/)
 res4@mpLandFillColor  = "grey56"


  plot(0) = gsn_csm_contour_map(wks, mean_temp1_log,res)

  plot(1) = gsn_csm_contour_map(wks, mean_temp2_log,res4)

  plots(0) = gsn_csm_contour_map(wks, mean_temp3_log,res)

  plots(1) = gsn_csm_contour_map(wks, mean_temp4_log,res4)
  
resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = False
; resP@pmLabelBarWidthF   = 0.60
; resP@pmLabelBarHeightF = 0.043
resP@lbOrientation   = "Horizontal"   
resP@pmLabelBarOrthogonalPosF = -0.03
resP@lbTitleString = "Vertical Diffusivity {log10[K(m^2/s)]}"
resP@lbTitleFontHeightF = 0.013
resP@gsnMaximize         = True                 ; make as large as possible   					
resP@lbLabelFontHeightF= 0.012 
resP@gsnPanelFigureStringsFontHeightF = 0.010
resP@gsnPanelYWhiteSpacePercent = 0.3

;resP@gsnPanelXF = (/-0.2,0.53/)


resP@gsnPanelFigureStrings= (/"w/oIDEMIX", "w/oIDEMIX", "w/IDEMIX", "w/IDEMIX"/) 
resP@gsnPanelFigureStringsFontHeightF = 0.012
;resP@cnLevelSelectionMode = "ExplicitLevels"
;resP@cnLevels    = (/ -7, -6, -5, -4, -3, -2 /)   ; set levels


	
;;resP@gsnPanelBottom = 0.9
;resP@gsnPanelTop = 0.92
;;resP@gsnPanelLeft = 0.1
;resP@gsnLeftString =  "Ve"
;resP@gsnRightString           = "(m^2/s)"  

 resP@gsnPanelTop       = 0.92
  resP@gsnPanelBottom    = 0.45
  resP@pmLabelBarWidthF  = 0.63

resPs                 = True
resPs@gsnMaximize        = True
resPs@gsnFrame         = False               ; save panel until both ready
resPs@gsnPaperOrientation = "landscape"
resPs@gsnPanelLabelBar   = True
;resPs@pmLabelBarWidthF   = 0.50
;resPs@pmLabelBarHeightF = 0.043
resPs@lbOrientation   = "Horizontal"
resPs@pmLabelBarOrthogonalPosF = -0.04
resPs@lbTitleString = "Vertical Diffusivity {log10[K(m^2/s)]}"
resPs@lbTitleFontHeightF = 0.013
resPs@gsnMaximize         = True                 ; make as large as possible

resPs@lbLabelFontHeightF= 0.012
resPs@gsnPanelFigureStringsFontHeightF = 0.012
resPs@gsnPanelYWhiteSpacePercent = 0.3
resPs@gsnPanelFigureStrings= (/"with IDEMIX", "with IDEMIX"/)
;;resPs@cnLevelSelectionMode = "ExplicitLevels"
;;resPs@cnLevels    = (/ -7, -6, -5, -4, -3, -2 /)   ; set levels


;resPs@gsnPanelBottom = 0.2
;;resPs@gsnPanelTop = 0.95
;resPs@gsnPanelRight = 0.6
;resPs@gsnLeftString =  "Tidal Dissipation Rate - Wilmes"
;resPs@gsnRightString           = "(W/kg)"

;resPs@gsnPanelTop       = 0.64
resPs@gsnPanelTop       = 0.75
resPs@pmLabelBarWidthF  = 0.6

txres1 = True
txres1@txFontHeightF = 0.014
txres1@txFont   = "helvetica-bold"
;gsn_text_ndc(wks,"w/oIDEMIX",0.15,0.85,txres1)
;gsn_text_ndc(wks,"PI",0.26,0.85,txres1)
;gsn_text_ndc(wks,"LGM",0.77,0.85,txres1)

gsn_text_ndc(wks,"(Di-Nezio et al., 2016)",0.77,0.85,txres1)
gsn_text_ndc(wks,"(Wilmes et al., 2021)",0.26,0.85,txres1)

txres2 = True
txres2@txFontHeightF = 0.013

;gsn_text_ndc(wks,"w/IDEMIX",0.15,0.50,txres2)


txres3 = True
txres3@txFontHeightF = 0.011

;gsn_text_ndc(wks,"(W/m^2)",0.94,0.84,txres3)



txres4 = True
txres4@txFontHeightF = 0.011

;gsn_text_ndc(wks,"(W/kg)",0.94,0.47,txres4)





;pres  = True

;resP@pmLabelBarOrthogonalPosF = -0.02    ; Move closer to plots
;resPs@pmLabelBarOrthogonalPosF = -0.02

;; pres@gsnPanelTop    = 1.0  ; set these if desired
;; pres@gsnPanelBottom = 0.0  ; must be between 0 and 1
;; pres@gsnPanelLeft   = 0.0
;; pres@gsnPanelRight  = 1.0

;pres@gsnPanelRowSpec = False


;  panel_two_sets(wks,plot,plots,(/1,2/),(/1,2/),resP,resPs,pres)
 




  gsn_panel(wks,plot,(/1,2/),resP)
  gsn_panel(wks,plots,(/1,2/),resPs)


 
; now we draw the whole thing and advance the frame.


  frame(wks)
  end










