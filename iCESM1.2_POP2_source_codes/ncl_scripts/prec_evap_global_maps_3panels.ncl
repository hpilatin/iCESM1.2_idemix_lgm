

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



     pop_pi_tamas = addfile ("/mnt/d/lgm_data2/buo_piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
    pop_lgm_tamasd = addfile ("/mnt/d/lgm_data2/buo_LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")

    pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/sel_buo_test_cmake_1750_1799_timmean.nc", "r")

    ;pop_pi = addfile ("/mnt/d/lgm_data/alldata_PI_ciso_j5_idemix_mvio_3920_3969_timmean.nc", "r")
     
    ;pop_pi = addfile ("/mnt/d/lgm_data/buo_pi_idemix_isotope_continued_4949_4999_timmean.nc", "r")

     pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
    pop_lgm = addfile ("/mnt/d/lgm_data/buo_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")
     pop_lgmd = addfile ("/mnt/d/lgm_data/surf_buo_dinezio_continued.pop.h.1810_1859_timmean.nc", "r")

wks_type = "pdf"


   wks  = gsn_open_wks(wks_type,"prec_evap_global_maps_pi_lgm_new2_3panels")


   cmap = read_colormap_file("cmocean_delta")
   ncmap = dimsizes(cmap(:,0))


res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = cmap(::-1,:)         ; set color map
;  res@cnFillPalette        = read_colormap_file("MPL_rai")
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels



;res@gsnPaperHeight = 14
;res@gsnPaperWidth = 9.5

;res@wkPaperHeightF = 12
;res@wkPaperWidthF = 20

   


res@vpWidthF           = 0.15
res@vpHeightF          = 0.20
res@tmXBLabelsOn  = True
res@tmXBLabelFontHeightF = 0.005
res@tmYLLabelFontHeightF = 0.005
res@cnFillPalette        = cmap
res@mpLandFillColor = "grey56"
res@cnMissingValFillColor = "grey56"
res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "

pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(3,graphic)
plots = new(3,graphic)
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
res@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
res@lbTitleFontHeightF= 0.009
res@lbTitleOffsetF = 0.01

res@tmXBTickSpacingF = 60.0
res@tmYLTickSpacingF = 30.0
res@cnMissingValFillPattern = "SolidFill"
res@cnMissingValFillColor = "grey56"


;    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res@cnMinLevelValF            =  20             ; min level
;    res@cnMaxLevelValF            =  200             ; max level
;    res@cnLevelSpacingF           =  20            ; contour interval



    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            = -16             ; min level
    res@cnMaxLevelValF            =  16             ; max level
    res@cnLevelSpacingF           =  2            ; contour interval


    ;res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    ;res@cnMinLevelValF            = -2             ; min level
    ;res@cnMaxLevelValF            =  2             ; max level
    ;res@cnLevelSpacingF           =  0.2            ; contour interval

  
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet


res1                      = True                 ; plot mods desired
res1@mpLandFillColor = "grey56"
res1@cnMissingValFillColor = "grey56"
 
  res1@gsnDraw              = False                ; don't draw yet
  res1@gsnFrame             = False                ; don't advance frame yet
  res1@cnFillOn             = True                 ; turn on color fill
   cmap2 = read_colormap_file("MPL_RdYlBu")
;   res1@cnFillPalette        = cmap(::-1,:)
  res1@cnFillPalette        = cmap2(::-1,:)         ; set color map
;!  res1@cnFillPalette        = read_colormap_file("MPL_RdYlBu")
  res1@cnLinesOn            = False                ; turn off contour lines
  res1@cnLineLabelsOn       = False                ; turn off contour line labels

res1@vpWidthF           = 0.35
res1@vpHeightF          = 0.20
res1@tmXBLabelsOn  = True
res1@tmXBLabelFontHeightF = 0.003
res1@tmYLLabelFontHeightF = 0.003

;res1@mpLandFillColor = "grey"

res1@lbLabelBarOn        = False
res1@cnInfoLabelOn       = False
res1@description = "  "
; res1@gsnLeftString =  "Vertical Diffusivity"
; res1@gsnRightString           = "(m^2/s)"
 res1@gsnCenterString = "  "
res1@tiMainString  = " "             ; add title
  ;res1@gsnMaximize         = True                 ; make as large as possible
;res1@tmXBLabelFontHeightF = 0.0025
;res1@tmYLLabelFontHeightF = 0.0025
;res1@lbLabelBarOn        = False

res1@cnMissingValFillPattern = "SolidFill"


;res1@cnMissingValFillColor = "grey52"
;res1@cnFillPalette        = read_colormap_file("MPL_rainbow")
res1@vpWidthF           = 0.13
res1@vpHeightF          = 0.18

res1@tmYLLabelsOn  = True
res1@tmXBLabelsOn  = True
res1@tmXBTickSpacingF     = 20.0
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/ -0.75, -0.5, -0.25, -0.05, -0.0075, -0.0005, -0.0025, -0.0005, 0 , 0.0005, 0.0025, 0.005, 0.0075, 0.05, 0.25, 0.50, 0.75/)
;res1@cnLevels    = (/ 0 , 0.000005, 0.000025, 0.00005, 0.000075, 0.0005, 0.0025, 0.005, 0.0075, 0.05, 0.25, 0.50, 0.75, 1/)
;res1@cnLevels    = ispan(-12,12,2)


res1@vpWidthF           = 0.35
res1@vpHeightF          = 0.20
res1@tmXBLabelFontHeightF = 0.007
res1@tmYLLabelFontHeightF = 0.007



;res1@tmXBLabelsOn  = True
res1@tmXBTickSpacingF     = 30.0
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/ -0.75, -0.5, -0.25, -0.05, -0.0075, -0.0005, -0.0025, -0.0005, 0 , 0.0005, 0.0025, 0.005, 0.0075, 0.05, 0.25, 0.50, 0.75/)
;res1@cnLevels    = (/ 0 , 0.000005, 0.000025, 0.00005, 0.000075, 0.0005, 0.0025, 0.005, 0.0075, 0.05, 0.25, 0.50, 0.75, 1/)
;res1@cnLevels    = ispan(-12,12,2)


res1@vpWidthF           = 0.3
res1@vpHeightF          = 0.15
;res1@tmXBLabelFontHeightF = 0.014
;res1@tmYLLabelFontHeightF = 0.014

;res1@gsnYAxisIrregular2Linear = False
;res1@trYReverse           = False
;res1@cnLineLabelFormat    = "0@;*.4f"
;res1@cnLineLabelFormat    = "*+^sg"
;res1@cnLineLabelFormat    =  "%6.1f"
;res1@cnConstFLabelFormat = "0@;*.8f"
;res1@cnConstFLabelFormat =  "%2f"
res1@lbLabelFontHeightF= 0.008

;res1@pmLabelBarParallelPosF = -0.01
res1@lbTitleString = "Salt [g/kg]"
res1@lbTitleFontHeightF= 0.009
res1@lbTitleOffsetF = 0.01


; res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res1@cnMinLevelValF            = -1.2             ; min level
;    res1@cnMaxLevelValF            =  1.2             ; max level
;    res1@cnLevelSpacingF           =  0.1            ; contour interval



;    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res1@cnMinLevelValF            = -2.             ; min level
;    res1@cnMaxLevelValF            =  2.             ; max level
;    res1@cnLevelSpacingF           =  0.2            ; contour interval

;    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res1@cnMinLevelValF            = -3.             ; min level
;    res1@cnMaxLevelValF            =  3.             ; max level
;    res1@cnLevelSpacingF           =  0.25            ; contour interval


    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res1@cnMinLevelValF            = -16.             ; min level
    res1@cnMaxLevelValF            =  16.             ; max level
    res1@cnLevelSpacingF           =  2            ; contour interval


;    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res1@cnMinLevelValF            =  20.             ; min level
;    res1@cnMaxLevelValF            =  200.             ; max level
;    res1@cnLevelSpacingF           =  20            ; contour interval



res1@tmXBTickSpacingF = 60.0
res1@tmYLTickSpacingF = 30.0


res1@cnMissingValFillPattern = "SolidFill"
;res1@cnMissingValFillColor = "lightgrey"


res2 = True
res2@tmXBTickSpacingF     = 30.0
res2@vpWidthF           = 0.3
res2@vpHeightF          = 0.15
res2@tmXBLabelFontHeightF = 0.007
res2@tmYLLabelFontHeightF = 0.007

res2@tmXBLabelsOn  = True
;res2@tmXBLabelFontHeightF = 0.005
;res2@tmYLLabelFontHeightF = 0.005
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
;res2@pmLabelBarParallelPosF = -0.01
;res2@lbTitleString = "Wind Stress [m/s]"
res2@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
res2@lbTitleFontHeightF= 0.009
res2@lbTitleOffsetF = 0.01
res2@cnFillOn             = True
;res2@cnFillPalette        = cmap(::-1,:)
res2@tmXBTickSpacingF = 60.0
res2@tmYLTickSpacingF = 30.0
; res2@cnMinLevelValF            = -3             ; min level
; res2@cnMaxLevelValF            =  3             ; max level
; res2@cnLevelSpacingF           =  0.2            ; contour interval

;res2@tmYLLabelsOn  = False
res2@tmXBLabelsOn  = True
;res2@tmXBLabelFontHeightF = 0.0010
;res2@tmYLLabelFontHeightF = 0.0010
res2@tmYLLabelsOn  = False
res2@tmXBLabelsOn  = True

;    res2@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res2@cnMinLevelValF            = -2.             ; min level
;    res2@cnMaxLevelValF            =  2.             ; max level
;    res2@cnLevelSpacingF           =  0.2            ; contour interval


  res2@gsnDraw              = False                ; don't draw yet
  res2@gsnFrame             = False                ; don't advance frame yet

res5 = res
res5@tmYLLabelsOn  = False

res4 = True
res4@tmXBTickSpacingF     = 30.0
res4@vpWidthF           = 0.3
res4@vpHeightF          = 0.15
res4@tmXBLabelFontHeightF = 0.007
res4@tmYLLabelFontHeightF = 0.007

res4@tmXBLabelsOn  = True
;res4@tmXBLabelFontHeightF = 0.005
;res4@tmYLLabelFontHeightF = 0.005
;res4@cnFillPalette        = cmap
res4@mpLandFillColor = "grey56"
res4@cnMissingValFillColor = "grey56"
res4@lbLabelBarOn        = False
res4@cnInfoLabelOn       = False                                                                                                                                                            
res4@cnLinesOn            = False                ; turn off contour lines
res4@cnLineLabelsOn       = False                ; turn off contour line labels
res4@tmXBLabelsOn  = True
res4@tmYLabelsOn  = False
;res4@pmLabelBarParallelPosF = -0.01
res4@lbTitleString = "Salt [g/kg]"
res4@lbTitleFontHeightF= 0.009
res4@lbTitleOffsetF = 0.01
;res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res4@cnMinLevelValF            = -1.2             ; min level
;    res4@cnMaxLevelValF            =  1.2             ; max level
;    res4@cnLevelSpacingF           =  0.1            ; contour interval
res4@mpLandFillColor = "grey56"
res4@cnMissingValFillColor = "grey56"
res4@cnFillOn             = True
res4@cnFillPalette        = cmap2(::-1,:)
;res4@cnFillPalette        = cmap(::-1,:)


;    res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels 
;    res4@cnMinLevelValF            = -2.             ; min level
;    res4@cnMaxLevelValF            =  2.             ; max level
;    res4@cnLevelSpacingF           =  0.2            ; contour interval


;    res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res4@cnMinLevelValF            = -3.             ; min level
;    res4@cnMaxLevelValF            =  3.             ; max level
;    res4@cnLevelSpacingF           =  0.25            ; contour interval


    res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res4@cnMinLevelValF            = -16.             ; min level
    res4@cnMaxLevelValF            =  16.             ; max level
    res4@cnLevelSpacingF           =  2            ; contour interval

;    res4@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res4@cnMinLevelValF            =  20.             ; min level
;    res4@cnMaxLevelValF            =  200.             ; max level
;    res4@cnLevelSpacingF           =  20            ; contour interval



;res4@cnFillPalette        = cmap(::-1,:)

res4@tmXBTickSpacingF = 60.0
res4@tmYLTickSpacingF = 30.0

res4@mpLandFillColor = "grey56"
res4@cnMissingValFillColor = "grey56"

;res4@tmYLLabelsOn  = False
;res4@tmXBLabelFontHeightF = 0.0010
;res4@tmYLLabelFontHeightF = 0.0010
;res4@tmYLLabelsOn  = False
res4@tmXBLabelsOn  = True


  res4@gsnDraw              = False                ; don't draw yet
  res4@gsnFrame             = False                ; don't advance frame yet


;res@gsnPaperOrientation  = "landscape"
;res1@gsnPaperOrientation  = "landscape"
res4@tmXBLabelsOn  = True
res4@tmYLLabelsOn  = False

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



;; Take the mean only for the Southern Ocean;;;
;;;########## Data (PI-without IDEMIX) ############;;;

;pi_stress_tamas         = cam_pi_tamas->TAUX(:,:,:)  ; Zonal surface wind stress (N/m2)
;pi_shf_tamas        = pop_pi_tamas->SHF(:,0:145,:)  ; Total surface heat flux (W/m2)
;pi_sfwf_tamas         = pop_pi_tamas->SFWF(:,0:145,:)  ; Freshwater Flux (kg/m2/s)
;pi_salt_tamas         = pop_pi_tamas->SALT(:,0,0:145,:)
;pi_prec_tamas         = pop_pi_tamas->PREC_F(:,0:145,:)
;pi_evap_tamas         = pop_pi_tamas->EVAP_F(:,0:145,:)

;pi_pd_tamas         = pop_pi_tamas->PD(:,0,0:145,:) ; surface density (g/cm3)

pi_shf_tamas        = pop_pi_tamas->SHF(:,1:nlat-1,:)  ; Total surface heat flux (W/m2)
pi_sfwf_tamas         = pop_pi_tamas->SFWF(:,1:nlat-1,:)  ; Freshwater Flux (kg/m2/s)
pi_salt_tamas         = pop_pi_tamas->SALT(:,0,1:nlat-1,:)
pi_prec_tamas         = pop_pi_tamas->PREC_F(:,1:nlat-1,:)
pi_evap_tamas         = pop_pi_tamas->EVAP_F(:,1:nlat-1,:)

pi_pd_tamas         = pop_pi_tamas->PD(:,0,1:nlat-1,:) ; surface density (g/cm3)

pi_salt2_tamas      = dim_avg_n(pi_salt_tamas,0)




;;;########## Means ############;;;
;pi_stress_mean_tamas    = dim_avg_n(pi_stress_tamas,0)
pi_shf_mean_tamas    = dim_avg_n(pi_shf_tamas,0)
pi_sfwf_mean_tamas    = dim_avg_n(pi_sfwf_tamas,0)
pi_salt_mean_tamas    = dim_avg_n(pi_salt_tamas,0)
pi_prec_mean_tamas    = dim_avg_n(pi_prec_tamas,0)
pi_evap_mean_tamas    = dim_avg_n(pi_evap_tamas,0)
pi_pd_mean_tamas    = dim_avg_n(pi_pd_tamas,0)*1000 ; conversion to kg/m3

;pi_stress_mean2_tamas    = dim_avg_n(pi_stress_mean_tamas,1)
pi_shf_mean2_tamas    = dim_avg_n(pi_shf_mean_tamas,1)
pi_sfwf_mean2_tamas    = dim_avg_n(pi_sfwf_mean_tamas,1)
pi_salt_mean2_tamas    = dim_avg_n(pi_salt_mean_tamas,1)
pi_prec_mean2_tamas    = dim_avg_n(pi_prec_mean_tamas,1)
pi_evap_mean2_tamas   = dim_avg_n(pi_evap_mean_tamas,1)
pi_pd_mean2_tamas    = dim_avg_n(pi_pd_mean_tamas,1)

;;;########## Data (LGM-without IDEMIX) ############;;;

;lgm_stress_tamas         = cam_lgm_tamas->TAUX(:,:,:)  ; Zonal surface wind stress (N/m2)
;lgm_shf_tamas         = pop_lgm_tamas->SHF(:,0:145,:)  ; Total surface heat flux (W/m2)
;lgm_sfwf_tamas         = pop_lgm_tamas->SFWF(:,0:145,:)  ; Freshwater Flux (kg/m2/s)
;lgm_salt_tamas         = pop_lgm_tamas->SALT(:,0,0:145,:)
;lgm_prec_tamas         = pop_lgm_tamas->PREC_F(:,0:145,:)
;lgm_evap_tamas         = pop_lgm_tamas->EVAP_F(:,0:145,:)

;lgm_pd_tamas         = pop_lgm_tamas->PD(:,0,0:145,:) ; surface density (g/cm3)


lgm_shf_tamas         = pop_lgm_tamas->SHF(:,1:nlat-1,:)  ; Total surface heat flux (W/m2)
lgm_sfwf_tamas         = pop_lgm_tamas->SFWF(:,1:nlat-1,:)  ; Freshwater Flux (kg/m2/s)
lgm_salt_tamas         = pop_lgm_tamas->SALT(:,0,1:nlat-1,:)
lgm_prec_tamas         = pop_lgm_tamas->PREC_F(:,1:nlat-1,:)
lgm_evap_tamas         = pop_lgm_tamas->EVAP_F(:,1:nlat-1,:)

lgm_pd_tamas         = pop_lgm_tamas->PD(:,0,1:nlat-1,:) ; surface density (g/cm3)

lgm_salt2_tamas      = dim_avg_n(lgm_salt_tamas,0)


print(max(lgm_salt2_tamas))
print(min(lgm_salt2_tamas))

;;;########## Means ############;;;
;lgm_stress_mean_tamas    = dim_avg_n(lgm_stress_tamas,0)
lgm_shf_mean_tamas   = dim_avg_n(lgm_shf_tamas,0)
lgm_sfwf_mean_tamas    = dim_avg_n(lgm_sfwf_tamas,0)
lgm_salt_mean_tamas    = dim_avg_n(lgm_salt_tamas,0)
lgm_prec_mean_tamas    = dim_avg_n(lgm_prec_tamas,0)
lgm_evap_mean_tamas    = dim_avg_n(lgm_evap_tamas,0)
lgm_pd_mean_tamas    = dim_avg_n(lgm_pd_tamas,0)*1000 ; conversion to kg/m3

;lgm_stress_mean2_tamas    = dim_avg_n(lgm_stress_mean_tamas,1)
lgm_shf_mean2_tamas    = dim_avg_n(lgm_shf_mean_tamas,1)
lgm_sfwf_mean2_tamas    = dim_avg_n(lgm_sfwf_mean_tamas,1)
lgm_salt_mean2_tamas    = dim_avg_n(lgm_salt_mean_tamas,1)
lgm_prec_mean2_tamas    = dim_avg_n(lgm_prec_mean_tamas,1)
lgm_evap_mean2_tamas    = dim_avg_n(lgm_evap_mean_tamas,1)
lgm_pd_mean2_tamas    = dim_avg_n(lgm_pd_mean_tamas,1)


;;;########## Data (PI-with IDEMIX) ############;;;

;pi_stress         = cam_pi->TAUX(:,:,:)  ; Zonal surface wind stress (N/m2)
;pi_shf         = pop_pi[:]->SHF(:,0:145,:)  ; Total surface heat flux (W/m2)
;pi_sfwf         = pop_pi[:]->SFWF(:,0:145,:)  ; Freshwater Flux (kg/m2/s)
;pi_salt         = pop_pi[:]->SALT(:,0,0:145,:)
;pi_prec         = pop_pi[:]->PREC_F(:,0:145,:)
;pi_evap         = pop_pi[:]->EVAP_F(:,0:145,:)
;pi_pd         = pop_pi[:]->PD(:,0,0:145,:) ; surface density (g/cm3)

pi_shf         = pop_pi->SHF(:,1:nlat-1,:)  ; Total surface heat flux (W/m2)
pi_sfwf         = pop_pi->SFWF(:,1:nlat-1,:)  ; Freshwater Flux (kg/m2/s)
pi_salt         = pop_pi->SALT(:,0,1:nlat-1,:)
pi_prec         = pop_pi->PREC_F(:,1:nlat-1,:)
pi_evap         = pop_pi->EVAP_F(:,1:nlat-1,:)

pi_pd         = pop_pi->PD(:,0,1:nlat-1,:) ; surface density (g/cm3)

pi_salt2      = dim_avg_n(pi_salt,0)

;;;########## Means ############;;;
;pi_stress_mean    = dim_avg_n(pi_stress,0)
pi_shf_mean    = dim_avg_n(pi_shf,0)
pi_sfwf_mean    = dim_avg_n(pi_sfwf,0)
pi_salt_mean    = dim_avg_n(pi_salt,0)
pi_prec_mean    = dim_avg_n(pi_prec,0)
pi_evap_mean    = dim_avg_n(pi_evap,0)
pi_pd_mean    = dim_avg_n(pi_pd,0)*1000 ; conversion to kg/m3

;pi_stress_mean2    = dim_avg_n(pi_stress_mean,1)
pi_shf_mean2    = dim_avg_n(pi_shf_mean,1)
pi_sfwf_mean2    = dim_avg_n(pi_sfwf_mean,1)
pi_salt_mean2    = dim_avg_n(pi_salt_mean,1)
pi_prec_mean2    = dim_avg_n(pi_prec_mean,1)
pi_evap_mean2    = dim_avg_n(pi_evap_mean,1)
pi_pd_mean2    = dim_avg_n(pi_pd_mean,1)

;;;########## Data (LGM-with IDEMIX) ############;;;

;lgm_stress         = cam_lgm->TAUX(:,:,:)  ; Zonal surface wind stress (N/m2)
;lgm_shf         = pop_lgm[:]->SHF(:,0:145,:)  ; Total surface heat flux (W/m2)
;lgm_sfwf         = pop_lgm[:]->SFWF(:,0:145,:)  ; Freshwater Flux (kg/m2/s)
;lgm_salt         = pop_lgm[:]->SALT(:,0,0:145,:)
;lgm_prec         = pop_lgm[:]->PREC_F(:,0:145,:)
;lgm_evap         = pop_lgm[:]->EVAP_F(:,0:145,:)

;lgm_pd         = pop_lgm[:]->PD(:,0,0:145,:) ; surface density (g/cm3)

lgm_shf         = pop_lgm->SHF(:,1:nlat-1,:)  ; Total surface heat flux (W/m2)
lgm_sfwf         = pop_lgm->SFWF(:,1:nlat-1,:)  ; Freshwater Flux (kg/m2/s)
lgm_salt         = pop_lgm->SALT(:,0,1:nlat-1,:) ; surface density (g/cm3)
lgm_prec         = pop_lgm->PREC_F(:,1:nlat-1,:)
lgm_evap         = pop_lgm->EVAP_F(:,1:nlat-1,:)

lgm_pd         = pop_lgm->PD(:,0,1:nlat-1,:) ; surface density (g/cm3)

lgm_salt2      = dim_avg_n(lgm_salt,0)

print(max(lgm_salt2))
print(min(lgm_salt2))


;;;########## Means ############;;;
;lgm_stress_mean    = dim_avg_n(lgm_stress,0)
lgm_shf_mean    = dim_avg_n(lgm_shf,0)
lgm_sfwf_mean    = dim_avg_n(lgm_sfwf,0)
lgm_salt_mean    = dim_avg_n(lgm_salt,0)
lgm_prec_mean    = dim_avg_n(lgm_prec,0)
lgm_evap_mean    = dim_avg_n(lgm_evap,0)
lgm_pd_mean    = dim_avg_n(lgm_pd,0)*1000 ; conversion to kg/m3

;lgm_stress_mean2    = dim_avg_n(lgm_stress_mean,1)
lgm_shf_mean2    = dim_avg_n(lgm_shf_mean,1)
lgm_sfwf_mean2    = dim_avg_n(lgm_sfwf_mean,1)
lgm_salt_mean2    = dim_avg_n(lgm_salt_mean,1)
lgm_prec_mean2    = dim_avg_n(lgm_prec_mean,1)
lgm_evap_mean2    = dim_avg_n(lgm_evap_mean,1)
lgm_pd_mean2    = dim_avg_n(lgm_pd_mean,1)

;printVarSummary(lgm_stress)



lgm_shf_tamasd         = pop_lgm_tamasd->SHF(:,1:nlat-1,:)  ; Total surface heat flux (W/m2)
lgm_sfwf_tamasd         = pop_lgm_tamasd->SFWF(:,1:nlat-1,:)  ; Freshwater Flux (kg/m2/s)
lgm_salt_tamasd         = pop_lgm_tamasd->SALT(:,0,1:nlat-1,:)
lgm_prec_tamasd         = pop_lgm_tamasd->PREC_F(:,1:nlat-1,:)
lgm_evap_tamasd         = pop_lgm_tamasd->EVAP_F(:,1:nlat-1,:)

lgm_pd_tamasd         = pop_lgm_tamasd->PD(:,0,1:nlat-1,:) ; surface density (g/cm3)

lgm_salt2_tamasd      = dim_avg_n(lgm_salt_tamasd,0)


print(max(lgm_salt2_tamasd))
print(min(lgm_salt2_tamasd))

;;;########## Means ############;;;
;lgm_stress_mean_tamas    = dim_avg_n(lgm_stress_tamas,0)
lgm_shf_mean_tamasd   = dim_avg_n(lgm_shf_tamasd,0)
lgm_sfwf_mean_tamasd    = dim_avg_n(lgm_sfwf_tamasd,0)
lgm_salt_mean_tamasd    = dim_avg_n(lgm_salt_tamasd,0)
lgm_prec_mean_tamasd    = dim_avg_n(lgm_prec_tamasd,0)
lgm_evap_mean_tamasd    = dim_avg_n(lgm_evap_tamasd,0)
lgm_pd_mean_tamasd   = dim_avg_n(lgm_pd_tamasd,0)*1000 ; conversion to kg/m3

;lgm_stress_mean2_tamas    = dim_avg_n(lgm_stress_mean_tamas,1)
lgm_shf_mean2_tamasd    = dim_avg_n(lgm_shf_mean_tamasd,1)
lgm_sfwf_mean2_tamasd    = dim_avg_n(lgm_sfwf_mean_tamasd,1)
lgm_salt_mean2_tamasd    = dim_avg_n(lgm_salt_mean_tamasd,1)
lgm_prec_mean2_tamasd    = dim_avg_n(lgm_prec_mean_tamasd,1)
lgm_evap_mean2_tamasd    = dim_avg_n(lgm_evap_mean_tamasd,1)
lgm_pd_mean2_tamasd    = dim_avg_n(lgm_pd_mean_tamasd,1)



lgm_shfd         = pop_lgmd->SHF(:,1:nlat-1,:)  ; Total surface heat flux (W/m2)
lgm_sfwfd         = pop_lgmd->SFWF(:,1:nlat-1,:)  ; Freshwater Flux (kg/m2/s)
lgm_saltd         = pop_lgmd->SALT(:,0,1:nlat-1,:) ; surface density (g/cm3)
lgm_precd         = pop_lgmd->PREC_F(:,1:nlat-1,:)
lgm_evapd         = pop_lgmd->EVAP_F(:,1:nlat-1,:)

lgm_pdd         = pop_lgmd->PD(:,0,1:nlat-1,:) ; surface density (g/cm3)

lgm_salt2d      = dim_avg_n(lgm_saltd,0)

print(max(lgm_salt2d))
print(min(lgm_salt2d))

lgm_shf_meand    = dim_avg_n(lgm_shfd,0)
lgm_sfwf_meand    = dim_avg_n(lgm_sfwfd,0)
lgm_salt_meand    = dim_avg_n(lgm_saltd,0)
lgm_prec_meand    = dim_avg_n(lgm_precd,0)
lgm_evap_meand    = dim_avg_n(lgm_evapd,0)
lgm_pd_meand    = dim_avg_n(lgm_pdd,0)*1000 ; conversion to kg/m3

;lgm_stress_mean2    = dim_avg_n(lgm_stress_mean,1)
lgm_shf_mean2d    = dim_avg_n(lgm_shf_meand,1)
lgm_sfwf_mean2d    = dim_avg_n(lgm_sfwf_meand,1)
lgm_salt_mean2d    = dim_avg_n(lgm_salt_meand,1)
lgm_prec_mean2d    = dim_avg_n(lgm_prec_meand,1)
lgm_evap_mean2d    = dim_avg_n(lgm_evap_meand,1)
lgm_pd_mean2d    = dim_avg_n(lgm_pd_meand,1)

;;;;############################ Assumptions########################

;;; beta= 0.801, alfa= 781*10^-7 (1/K) (these coefficients are taken from Atmosphere Ocean Dynamics (Adrian E. Gill) Table A3.1

;;; Cp= 4000 J/kg*K, sea water density (swd)= 1022.4 kg/m^3 , ref density water = 1000 kg/m3

;;;################################################################

pd_ref_water = 1000 ; kg/m3

;;;; without IDEMIX ;;;;

pi_buoyancy_heat_flux_tamas = -(781*1d-7*pi_shf_mean_tamas)/4000  ; (alpha * surface heat flux) / Cp  ; (kg/m2/s)

printVarSummary(pi_buoyancy_heat_flux_tamas)
print(max(pi_buoyancy_heat_flux_tamas))
print(min(pi_buoyancy_heat_flux_tamas))

pi_buoyancy_freshwater_flux_tamas = pi_sfwf_mean_tamas*7.8*(1d-4)*avg(pi_salt_mean_tamas)*(pi_pd_mean_tamas/pd_ref_water) ; (kg/m2/s)


printVarSummary(pi_buoyancy_freshwater_flux_tamas)
print(max(pi_buoyancy_freshwater_flux_tamas))
print(min(pi_buoyancy_freshwater_flux_tamas))
print(avg(pi_buoyancy_freshwater_flux_tamas))

pi_buoyancy_flux_tamas = (pi_buoyancy_heat_flux_tamas+pi_buoyancy_freshwater_flux_tamas)*1d6

printVarSummary(pi_buoyancy_flux_tamas)
print(max(pi_buoyancy_flux_tamas))
print(min(pi_buoyancy_flux_tamas))
;;;;

lgm_buoyancy_heat_flux_tamas = -(781*1d-7*lgm_shf_mean_tamas)/4000  ; (alpha * surface heat flux) / Cp  ; (kg/m2/s)

printVarSummary(lgm_buoyancy_heat_flux_tamas)
print(max(lgm_buoyancy_heat_flux_tamas))
print(min(lgm_buoyancy_heat_flux_tamas))

lgm_buoyancy_freshwater_flux_tamas = lgm_sfwf_mean_tamas*7.8*(1d-4)*avg(lgm_salt_mean_tamas)*(lgm_pd_mean_tamas/pd_ref_water) ; (kg/m2/s)

printVarSummary(lgm_buoyancy_freshwater_flux_tamas)
print(max(lgm_buoyancy_freshwater_flux_tamas))
print(min(lgm_buoyancy_freshwater_flux_tamas))
print(avg(lgm_buoyancy_freshwater_flux_tamas))

lgm_buoyancy_flux_tamas = (lgm_buoyancy_heat_flux_tamas+lgm_buoyancy_freshwater_flux_tamas)*1d6

printVarSummary(lgm_buoyancy_flux_tamas)
print(max(lgm_buoyancy_flux_tamas))
print(min(lgm_buoyancy_flux_tamas))


;;;; with IDEMIX ;;;;
pi_buoyancy_heat_flux = -(781*1d-7*pi_shf_mean)/4000  ; (alpha * surface heat flux) / Cp  ; (kg/m2/s)

printVarSummary(pi_buoyancy_heat_flux)
print(max(pi_buoyancy_heat_flux))
print(min(pi_buoyancy_heat_flux))

pi_buoyancy_freshwater_flux = pi_sfwf_mean*7.8*(1d-4)*avg(pi_salt_mean)*(pi_pd_mean/pd_ref_water) ; (kg/m2/s)




printVarSummary(pi_buoyancy_freshwater_flux)
print(max(pi_buoyancy_freshwater_flux))
print(min(pi_buoyancy_freshwater_flux))
print(avg(pi_buoyancy_freshwater_flux))

pi_buoyancy_flux = (pi_buoyancy_heat_flux+pi_buoyancy_freshwater_flux)*1d6

printVarSummary(pi_buoyancy_flux)
print(max(pi_buoyancy_flux))
print(min(pi_buoyancy_flux))
;;;;

lgm_buoyancy_heat_flux = -(781*1d-7*lgm_shf_mean)/4000  ; (alpha * surface heat flux) / Cp  ; (kg/m2/s)

printVarSummary(lgm_buoyancy_heat_flux)
print(max(lgm_buoyancy_heat_flux))
print(min(lgm_buoyancy_heat_flux))

lgm_buoyancy_freshwater_flux = lgm_sfwf_mean*7.8*(1d-4)*avg(lgm_salt_mean)*(lgm_pd_mean/pd_ref_water) ; (kg/m2/s)

printVarSummary(lgm_buoyancy_freshwater_flux)
print(max(lgm_buoyancy_freshwater_flux))
print(min(lgm_buoyancy_freshwater_flux))
print(avg(lgm_buoyancy_freshwater_flux))

lgm_buoyancy_flux = (lgm_buoyancy_heat_flux+lgm_buoyancy_freshwater_flux)*1d6

printVarSummary(lgm_buoyancy_flux)
print(max(lgm_buoyancy_flux))
print(min(lgm_buoyancy_flux))


pi_buoyancy_flux_tamas2 = abs(doubletofloat(pi_buoyancy_flux_tamas))
lgm_buoyancy_flux_tamas2 = abs(doubletofloat(lgm_buoyancy_flux_tamas))

pi_buoyancy_flux2 = abs(doubletofloat(pi_buoyancy_flux))
lgm_buoyancy_flux2 = abs(doubletofloat(lgm_buoyancy_flux))


;;sub1 = flt2dble(abs(pi_stress_mean) - abs(pi_stress_mean_tamas))    
;sub1 = abs(pi_salt2) - abs(pi_salt2_tamas) ;correct
;sub1 = abs(pi_buoyancy_heat_flux*1d6) - abs(pi_buoyancy_heat_flux_tamas*1d6)

;sub1 = abs(pi_buoyancy_freshwater_flux*1d6) - abs(pi_buoyancy_freshwater_flux_tamas*1d6)
;sub2 = abs(pi_buoyancy_flux) - abs(pi_buoyancy_flux_tamas) ;correct
;;sub2 = abs(pi_buoyancy_flux) - abs(pi_buoyancy_flux_tamas)

diff1 = ((pi_prec_mean_tamas*1d6) - (pi_evap_mean_tamas*1d6))
diff2 = ((pi_prec_mean*1d6) - (pi_evap_mean*1d6))



diff3 = ((lgm_prec_mean_tamas*1d6) - (lgm_evap_mean_tamas*1d6))
diff4 = ((lgm_prec_mean*1d6) - (lgm_evap_mean*1d6))


;diff5 = abs((pi_prec_mean*1d6)-(pi_prec_mean_tamas*1d6)) - abs((pi_evap_mean*1d6)-(pi_evap_mean_tamas*1d6))
;diff6 = abs((lgm_prec_mean*1d6)-(lgm_prec_mean_tamas*1d6)) - abs((lgm_evap_mean*1d6)-(lgm_evap_mean_tamas*1d6))


diff5 = diff2 - diff1
diff6 = diff4 - diff3


diff7 = ((lgm_prec_mean_tamasd*1d6) - (lgm_evap_mean_tamasd*1d6))
diff8 = ((lgm_prec_meand*1d6) - (lgm_evap_meand*1d6))

diff9 = diff8 - diff7


;sub1 = (-pi_evap_mean*1d6) - (-pi_evap_mean_tamas*1d6)
;;sub1 = (pi_prec_mean*1d6) - (pi_prec_mean_tamas*1d6)
;;sub1 = abs(lgm_evap_mean_tamas*1d6) - abs(pi_evap_mean_tamas*1d6)
;;sub1 = abs(lgm_prec_mean_tamas*1d6) - abs(pi_prec_mean_tamas*1d6)
;;sub2 = diff5
;sub2 = (-pi_prec_mean*1d6) - (-pi_prec_mean_tamas*1d6)

;sub1 = (pi_prec_mean_tamas*1d6) - (pi_evap_mean_tamas*1d6)
;sub2 = (lgm_prec_mean_tamas*1d6) - (lgm_evap_mean_tamas*1d6)

;sub1 = diff5
;sub2 = diff5


sub1 = (pi_evap_mean*1d6) - (pi_evap_mean_tamas*1d6)
sub2 = (pi_prec_mean*1d6) - (pi_prec_mean_tamas*1d6)

printVarSummary(diff1)
print(max(diff1))
print(min(diff1))

printVarSummary(diff2)
print(max(diff2))
print(min(diff2))


printVarSummary(diff3)
print(max(diff3))
print(min(diff3))

printVarSummary(diff4)
print(max(diff4))
print(min(diff4))



printVarSummary(pi_buoyancy_freshwater_flux)
print(max(pi_buoyancy_freshwater_flux*1d6))
print(min(pi_buoyancy_freshwater_flux*1d6))
print(max(pi_buoyancy_freshwater_flux_tamas*1d6))
print(min(pi_buoyancy_freshwater_flux_tamas*1d6))


;sub1 = abs(lgm_stress_mean_tamas) - abs(pi_stress_mean_tamas) ;correct
;sub2 = lgm_buoyancy_flux_tamas2 - pi_buoyancy_flux_tamas2 ;correct


;;sub3 = flt2dble(abs(lgm_stress_mean) - abs(lgm_stress_mean_tamas))    
;sub3 = abs(lgm_salt2) - abs(lgm_salt2_tamas) ;correct
;sub3 = abs(lgm_buoyancy_heat_flux*1d6) - abs(lgm_buoyancy_heat_flux_tamas*1d6)

;sub3 = abs(lgm_buoyancy_freshwater_flux*1d6) - abs(lgm_buoyancy_freshwater_flux_tamas*1d6)
;sub4 = lgm_buoyancy_flux2 - lgm_buoyancy_flux_tamas2     ;correct
;sub4 = abs(lgm_buoyancy_flux) - abs(lgm_buoyancy_flux_tamas)



;sub3 = (-lgm_evap_mean*1d6) - (-lgm_evap_mean_tamas*1d6)
;;sub3 = (lgm_prec_mean*1d6) - (lgm_prec_mean_tamas*1d6)
;;sub3 = abs(lgm_evap_mean*1d6) - abs(pi_evap_mean*1d6)
;;sub3 = abs(lgm_prec_mean*1d6) - abs(pi_prec_mean*1d6)
;;sub4 = diff6
;sub4 = (-lgm_prec_mean*1d6) - (-lgm_prec_mean_tamas*1d6)

;sub3 = (pi_prec_mean*1d6) - (pi_evap_mean*1d6)
;sub4 = (lgm_prec_mean*1d6) - (lgm_evap_mean*1d6)

;sub3 = diff6
;sub4 = diff6

sub3 = (lgm_evap_mean*1d6) - (lgm_evap_mean_tamas*1d6)
sub4 = (lgm_prec_mean*1d6) - (lgm_prec_mean_tamas*1d6)

printVarSummary(lgm_buoyancy_freshwater_flux)
print(max(lgm_buoyancy_freshwater_flux*1d6))
print(min(lgm_buoyancy_freshwater_flux*1d6))
print(max(lgm_buoyancy_freshwater_flux_tamas*1d6))
print(min(lgm_buoyancy_freshwater_flux_tamas*1d6))



;sub3 = abs(lgm_stress_mean) - abs(pi_stress_mean) ;correct
;sub4 = lgm_buoyancy_flux2 - pi_buoyancy_flux2 ;correct

printVarSummary(sub3)
print(max(sub3))
print(min(sub3))



sub5 = (lgm_evap_meand*1d6) - (lgm_evap_mean_tamasd*1d6)
sub6 = (lgm_prec_meand*1d6) - (lgm_prec_mean_tamasd*1d6)


;pi_buoyancy_flux_tamas2@lat2d = lat
;pi_buoyancy_flux_tamas2@lon2d = lon

;lgm_buoyancy_flux_tamas2@lat2d = lat
;lgm_buoyancy_flux_tamas2@lon2d = lon



;pi_buoyancy_flux2@lat2d = lat
;pi_buoyancy_flux2@lon2d = lon

;lgm_buoyancy_flux2@lat2d = lat
;lgm_buoyancy_flux2@lon2d = lon



;sub1@lat2d = cam_pi->lat
;sub1@lon2d =cam_pi->lon

;sub3@lat2d = cam_pi->lat
;sub3@lon2d = cam_pi->lon


;sub2@lat2d = doubletofloat(lat)
;sub2@lon2d = doubletofloat(lon)

;sub4@lat2d = doubletofloat(lat)
;sub4@lon2d = doubletofloat(lon)


sub1@lat2d = lat
sub1@lon2d = lon

sub3@lat2d = lat
sub3@lon2d = lon


sub2@lat2d = lat
sub2@lon2d = lon

sub4@lat2d = lat
sub4@lon2d = lon

sub5@lat2d = lat
sub5@lon2d = lon

sub6@lat2d = lat
sub6@lon2d = lon


printVarSummary(sub1)
printVarSummary(sub2)
printVarSummary(sub3)
printVarSummary(sub4)




;res@gsnPaperOrientation = "landspace"
;res1@gsnPaperOrientation = "landspace"
;res2@gsnPaperOrientation = "landspace"
;res4@gsnPaperOrientation = "landspace"

;res@wkPaperHeightF = 12
;res@wkPaperWidthF = 20

;res1@wkPaperHeightF = 12
;res1@wkPaperWidthF = 20
;res2@wkPaperHeightF = 12
;res2@wkPaperWidthF = 20

;res4@wkPaperHeightF = 12
;res4@wkPaperWidthF = 20


res@pmLabelBarOrthogonalPosF = 0.2
res1@pmLabelBarOrthogonalPosF = 0.2
res2@pmLabelBarOrthogonalPosF = 0.2
res4@pmLabelBarOrthogonalPosF = 0.2


plots(0) = gsn_csm_contour_map(wks,sub1,res)
plots(1) = gsn_csm_contour_map(wks,sub3,res5)
plots(2) = gsn_csm_contour_map(wks,sub5,res5)
plot(0) = gsn_csm_contour_map(wks,sub2,res1)
plot(1) = gsn_csm_contour_map(wks,sub4,res4)
plot(2) = gsn_csm_contour_map(wks,sub6,res4)



;lat = pop_pi[0]->TLAT(0:145,0)




print(max(sub2))
print(min(sub2))

print(max(sub4))
print(min(sub4))



printVarSummary(sub2)

resP = True
;  resP@gsnMaximize         = True                 ; make as large as possible
;  resP@gsnPaperOrientation = "portrait"           ; force portrait
;  resP@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  resP@gsnMaximize        = False
  resP@gsnFrame         = False               ; save panel until both ready
  resP@gsnPaperOrientation = "landscape"
  resP@gsnPanelLabelBar   = True
  resP@pmLabelBarWidthF   = 0.60
  ;!resP@pmLabelBarHeightF = 0.043

resP@pmLabelBarHeightF = 0.04

;   resP@pmLabelBarWidthF   = 0.07
;   resP@pmLabelBarHeightF = 0.33
 ;!resP@pmLabelBarOrthogonalPosF = -0.042
   ;resP@pmLabelBarParallelPosF = -0.01
;   resP@lbTitleExtentF = 0.01
;resP@lbOrientation        = "Vertical"
;resP@lbLabelPosition      = "Right"           ; label position


;resP@lbTitleString = "P - E [10^-6 kg/m2/s]"
;resP@lbTitleString = "Precipitation Flux [g/kg]"

resP@lbTitleString = "Precipitation Flux [10^-6 kg/m2/s]"
;resP@lbTitleString = " Evaporation Flux [10^-6 kg/m2/s]"
;resP@lbTitleString = "Freshwater Flux [10^-6 kg/m2/s]"
;resP@gsnPaperHeight = 80
;resP@gsnPaperWidth = 80

resP@lbTitleFontHeightF = 0.0095
resP@lbLabelFontHeightF= 0.0085

resP@pmLabelBarOrthogonalPosF = -0.032

;!resP@lbTitleFontHeightF = 0.011
;resP@lbTitlePosition     = "Top"
;resP@lbTitleDirection     = "Across"                          ; letter angle
;resP@lbTitleAngleF        = 270.                               ; title angle
;resP@lbTitleDirection = "Down"
;resP@lbTitleOffsetF = -0.003
;!resP@lbTitleFont   = "helvetica-bold"

;resP@gsnPanelYWhiteSpacePercent = 1
;resP@gsnPanelXWhiteSpacePercent = 1
;!resP@lbLabelFontHeightF= 0.010
resP@gsnPanelFigureStringsFontHeightF = 0.0085
;resP@gsnPanelFigureStrings= (/"(a)", "(b)"/)
resP@gsnPanelFigureStrings= (/"PI", "LGM-W", "LGM-D"/)

resP@gsnPanelTop = 0.95
resP@gsnPanelBottom = 0.61

;resP@gsnPanelFigureStringsFontHeightF = 0.012
;!resP@gsnPanelYWhiteSpacePercent = 0.3

resP@gsnPanelXWhiteSpacePercent = 0.1

;resP@gsnPanelXF = (/0.02,0.55/)
;resP@gsnPanelYF = (/0.95,0.95/)


;resP@lbLabelOffsetF = -0.2

resP1 = True


resP1@lbTitleString = " Evaporation Flux [10^-6 kg/m2/s]"
;resP1@lbTitleString = "Precipitation Flux [10^-6 kg/m2/s]"
;resP1@lbTitleString = "P - E [10^-6 kg/m2/s]"
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
resP1@gsnPanelFigureStrings= (/"(a)", "(b)", "(c)"/)
;resP1@gsnPanelFigureStrings= (/"(w/IDEMIX-w/oIDEMIX)", "(w/IDEMIX-w/oIDEMIX)"/)
resP1@gsnPanelFigureStringsFontHeightF = 0.011

;resP1@gsnPanelXF = (/0.02,0.55/)
;resP1@gsnPanelYF = (/0.5,0.5/)


gsn_panel(wks, plot, (/1,3/), resP)
gsn_panel(wks, plots, (/1,3/), resP1)


txres2 = True
txres2@txFontHeightF = 0.011
txres2@txFont   = "helvetica-bold"



;gsn_text_ndc(wks,"(w/IDEMIX-w/oIDEMIX)",0.12,0.91,txres2)


;gsn_text_ndc(wks,"PI (w/IDEMIX-w/oIDEMIX)",0.27,0.96,txres2)
;gsn_text_ndc(wks,"LGM (w/IDEMIX-w/oIDEMIX)",0.77,0.96,txres2)

;gsn_text_ndc(wks,"LGM - PI (w/oIDEMIX)",0.27,0.99,txres2)
;gsn_text_ndc(wks,"LGM - PI (w/IDEMIX)",0.77,0.99,txres2)

;  psres = True                                                               
;  maximize_output(wks,psres)  ; calls draw and frame for you 

;frame(wks)
end

