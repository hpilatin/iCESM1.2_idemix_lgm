

begin


;   pop_pi_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
   pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/pre_evap_cam_test_cmake_1838_1887_timmean.nc", "r")
;   pop_pi = addfile ("/mnt/d/lgm_data/buo_pi_idemix_isotope_continued_4900_4999.nc", "r")
   pop_lgm = addfile ("/mnt/d/lgm_data/pre_evap_cam_lgm_idemix_isotope_continued_5944_5993_timmean.nc", "r")

  wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"prec_evap_cam_lgm_idemix")


;   cmap = read_colormap_file("cmocean_delta")
 
  cmap2 = read_colormap_file("MPL_RdYlBu")
  ncmap = dimsizes(cmap2(:,0))


res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = cmap2(::-1,:)        ; set color map
;  res@cnFillPalette        = read_colormap_file("MPL_rai")
  res@cnLinesOn            = False                ; turn off contour lines
  res@cnLineLabelsOn       = False                ; turn off contour line labels


pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(2,graphic)
res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True
res@tmXBTickSpacingF     = 30.0
res@vpWidthF           = 0.3
res@vpHeightF          = 0.15
res@tmXBLabelFontHeightF = 0.007
res@tmYLLabelFontHeightF = 0.007
res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
res@cnMinLevelValF            = -6.             ; min level
res@cnMaxLevelValF            =  6.             ; max level
res@cnLevelSpacingF           =  0.5            ; contour interval


lgm_prec_conv_tamas        = pop_lgm_tamas->PRECC(0,:,:)  ;Convective precipitation rate (liq + ice), m/s
                                                                                                                                                  lgm_prec_large_tamas        = pop_lgm_tamas->PRECL(0,:,:)  ; Large-scale (stable) precipitation rate (liq + ice), m/s


lgm_prec_total_tamas = lgm_prec_conv_tamas + lgm_prec_large_tamas


lgm_latent_tamas            = pop_lgm_tamas->LHFLX(0,:,:)  ; Surface latent heat flux, W/m2
                                                                         
lgm_evap_tamas = lgm_latent_tamas / 2.5e9  ; LHFLX in W/m^2 → evap in m/s

lgm_diff_tamas = lgm_prec_total_tamas - lgm_evap_tamas ; m/s                               
lgm_diff_unit_tamas = (lgm_diff_tamas*1000)*86400     ; from m/s to mm/day

printVarSummary(lgm_diff_unit_tamas)

print(max(lgm_diff_unit_tamas))

print(min(lgm_diff_unit_tamas))




lgm_prec_conv        = pop_lgm->PRECC(0,:,:)  ;Convective precipitation rate (liq + ice), m/s

lgm_prec_large        = pop_lgm->PRECL(0,:,:)  ; Large-scale (stable) precipitation rate (liq + ice), m/s


lat1 = pop_lgm->lat
lon1 = pop_lgm->lon


lgm_prec_total = lgm_prec_conv + lgm_prec_large 


lgm_latent            = pop_lgm->LHFLX(0,:,:)  ; Surface latent heat flux, W/m2


lgm_evap = lgm_latent / 2.5e9  ; LHFLX in W/m^2 → evap in m/s

lgm_diff = lgm_prec_total - lgm_evap ; m/s

lgm_diff_unit = (lgm_diff*1000)*86400     ; from m/s to mm/day

printVarSummary(lgm_diff_unit)

print(max(lgm_diff_unit))

print(min(lgm_diff_unit))

;lgm_diff_unit@lat2d = lat1
;lgm_diff_unit@lon2d = lon1

plot(0) = gsn_csm_contour_map(wks,lgm_diff_unit,res)
plot(1) = gsn_csm_contour_map(wks,lgm_diff_unit_tamas,res)


resP = True
gsn_panel(wks, plot, (/1,2/), resP)

end

