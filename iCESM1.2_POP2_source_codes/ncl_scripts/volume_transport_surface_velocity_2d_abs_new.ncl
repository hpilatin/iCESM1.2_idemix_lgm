


begin

  pop_pi_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
  pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/everything_except_moc_test_cmake.pop.h.1830_1879_timmean.nc", "r")

  pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")

  pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")

  pop_lgm_tamasd = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")

  pop_lgmd = addfile ("/mnt/d/lgm_data/select_dinezio_continued_1810_1859_timmean.nc", "r")

   pop_lgm_tamasd_rho = addfile ("/mnt/d/lgm_data/select_dinezio_continued_1810_1859_timmean.nc", "r")



   pop_ref = addfile ("/mnt/d/b.ie12.Bi1850C5BPRDCISO.f19_g16_LGM21ka_ciso_idemix_j5.pop.h.5374.nc", "r")


      pop_pi_angle  = addfile ("/mnt/d/lgm_data/angle_pi.nc", "r")
   pop_lgm_angle  = addfile ("/mnt/d/lgm_data/angle_lgm.nc", "r")


wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"volume_transport_surface_velocity_6panel_2d_diff_new")


   cmap = read_colormap_file("temp_19lev")
   ncmap = dimsizes(cmap(:,0))

  cmap2 = read_colormap_file("MPL_RdYlBu")

res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
 
    cmap2 = read_colormap_file("MPL_RdYlBu")
  ;res@cnFillPalette        = "BlWhRe"         ; set color map
  res@cnFillPalette        = cmap2(::-1,:)
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
res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True
res@tmXBTickSpacingF     = 15.0
res@vpWidthF           = 0.3
res@vpHeightF          = 0.15
res@tmXBLabelFontHeightF = 0.006
res@tmYLLabelFontHeightF = 0.006


res@lbLabelFontHeightF= 0.008
;res@pmLabelBarParallelPosF = -0.01
res@lbTitleString = "Volume Transport [m3/s]"
res@lbTitleFontHeightF= 0.009
res@lbTitleOffsetF = 0.01

res@tmXBTickSpacingF = 30.0
res@tmYLTickSpacingF = 15.0
res@cnMissingValFillPattern = "SolidFill"
res@cnMissingValFillColor = "grey56"

;  res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
;    res@cnMinLevelValF            = 50             ; min level
;   res@cnMaxLevelValF            =  600             ; max level
;   res@cnLevelSpacingF           =  25            ; contour interval


  res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels 
  res@cnMinLevelValF            = -150             ; min level          
  res@cnMaxLevelValF            =  150             ; max level                                                     
 res@cnLevelSpacingF           =  25            ; contour interval

 
  res@mpMaxLatF             = 90.                  ; zoom in on map
  res@mpMinLatF             =  30.
  res@mpMinLonF             = -60.
  res@mpMaxLonF             = 30.



rho_pi_tamas          = pop_pi_tamas->RHO(0,:,:,:)
rho_pi          = pop_pi->RHO(0,:,:,:)
rho_lgm_tamas          = pop_lgm_tamas->RHO(0,:,:,:)
rho_lgm          = pop_lgm->RHO(0,:,:,:)
rho_lgm_tamasd          = pop_lgm_tamasd_rho->RHO(0,:,:,:)
rho_lgmd          = pop_lgmd->RHO(0,:,:,:)


uvel_pi_tamas          = pop_pi_tamas->UVEL(0,:,:,:)
vvel_pi_tamas          = pop_pi_tamas->VVEL(0,:,:,:)

uvel_pi          = pop_pi->UVEL(0,:,:,:)
vvel_pi          = pop_pi->VVEL(0,:,:,:)

uvel_lgm_tamas          = pop_lgm_tamas->UVEL(0,:,:,:)
vvel_lgm_tamas          = pop_lgm_tamas->VVEL(0,:,:,:)

uvel_lgm          = pop_lgm->UVEL(0,:,:,:)
vvel_lgm          = pop_lgm->VVEL(0,:,:,:)


uvel_lgm_tamasd          = pop_lgm_tamasd->UVEL(0,:,:,:)
vvel_lgm_tamasd          = pop_lgm_tamasd->VVEL(0,:,:,:)

uvel_lgmd          = pop_lgmd->UVEL(0,:,:,:)
vvel_lgmd          = pop_lgmd->VVEL(0,:,:,:)


rot_pi         = pop_pi_angle->ANGLE
  rot_lgm         = pop_lgm_angle->ANGLE

  urot_pi        = uvel_pi(0,:,:)*cos(rot_pi) - vvel_pi(0,:,:)*sin(rot_pi)  ; perform standard rotation
  vrot_pi        = uvel_pi(0,:,:)*sin(rot_pi) + vvel_pi(0,:,:)*cos(rot_pi)

  urot_lgm        = uvel_lgm(0,:,:)*cos(rot_lgm) - vvel_lgm(0,:,:)*sin(rot_lgm)  ; perform standard rotation
  vrot_lgm        = uvel_lgm(0,:,:)*sin(rot_lgm) + vvel_lgm(0,:,:)*cos(rot_lgm)

   urot_pi_tamas        = uvel_pi_tamas(0,:,:)*cos(rot_pi) - vvel_pi_tamas(0,:,:)*sin(rot_pi)  ; perform standard rotation
  vrot_pi_tamas        = uvel_pi_tamas(0,:,:)*sin(rot_pi) + vvel_pi_tamas(0,:,:)*cos(rot_pi)

  urot_lgm_tamas    = uvel_lgm_tamas(0,:,:)*cos(rot_lgm) - vvel_lgm_tamas(0,:,:)*sin(rot_lgm)  ; perform standard rotation
  vrot_lgm_tamas    = uvel_lgm_tamas(0,:,:)*sin(rot_lgm) + vvel_lgm_tamas(0,:,:)*cos(rot_lgm)

  
  urot_lgm_tamasd    = uvel_lgm_tamasd(0,:,:)*cos(rot_lgm) - vvel_lgm_tamasd(0,:,:)*sin(rot_lgm)  ; perform standard rotation
  vrot_lgm_tamasd    = uvel_lgm_tamasd(0,:,:)*sin(rot_lgm) + vvel_lgm_tamasd(0,:,:)*cos(rot_lgm)  

  urot_lgmd        = uvel_lgmd(0,:,:)*cos(rot_lgm) - vvel_lgmd(0,:,:)*sin(rot_lgm)  ; perform standard rotation
  vrot_lgmd        = uvel_lgmd(0,:,:)*sin(rot_lgm) + vvel_lgmd(0,:,:)*cos(rot_lgm)


  tarea = pop_ref->TAREA          ; horizontal cell area [lat, lon]
  dz = pop_ref->dz                ; vertical thickness [depth]


rho_pi_tamas2 = rho_pi_tamas*1000 ; from g/cm3 to kg/m3
rho_pi2 = rho_pi*1000 ; from g/cm3 to kg/m3

rho_lgm_tamas2 = rho_lgm_tamas*1000 ; from g/cm3 to kg/m3
rho_lgm2 = rho_lgm*1000 ; from g/cm3 to kg/m3

rho_lgm_tamasd2 = rho_lgm_tamasd*1000 ; from g/cm3 to kg/m3
rho_lgmd2 = rho_lgmd*1000 ; from g/cm3 to kg/m3

dim_sizes = dimsizes(uvel_lgm)

printVarSummary(dim_sizes)


 tarea = pop_ref->TAREA          ; horizontal cell area [lat, lon]
  dz = pop_ref->dz                ; vertical thickness [depth]

lgm_vel = sqrt((uvel_lgm*uvel_lgm*1d-4)+(vvel_lgm*vvel_lgm*1d-4))
pi_vel = sqrt((uvel_pi*uvel_pi*1d-4)+(vvel_pi*vvel_pi*1d-4))
lgmd_vel = sqrt((uvel_lgmd*uvel_lgmd*1d-4)+(vvel_lgmd*vvel_lgmd*1d-4))


lgm_vel_tamas = sqrt((uvel_lgm_tamas*uvel_lgm_tamas*1d-4)+(vvel_lgm_tamas*vvel_lgm_tamas*1d-4))
pi_vel_tamas = sqrt((uvel_pi_tamas*uvel_pi_tamas*1d-4)+(vvel_pi_tamas*vvel_pi_tamas*1d-4))
lgm_vel_tamasd = sqrt((uvel_lgm_tamasd*uvel_lgm_tamasd*1d-4)+(vvel_lgm_tamasd*vvel_lgm_tamasd*1d-4))

  lon   = pop_pi->TLONG                   ; read in original lat/lon data
  lat   = pop_pi->TLAT


; ================================================;
; when a coordinate array contains one more element than the corresponding
; dimension in the data array, the cell boundary is automatically turned on,
; which causes discrete rasterization to use the coordinate values as
; divisions between data cells, rather than to treat each coordinate value
; as a cell center.
; ================================================;
  dims      = dimsizes(lon)         ; get dimension sizes
  dims(1)   = dims(1) + 1               ; increase longitude by 1
  nlat      = dims(0)
  nlon      = dims(1)

  lon1      = new(dims,typeof(lon)) ; create new arrays
  lat1      = new(dims,typeof(lat))

  lat1(:,1:nlon-1)  = lat        ; fill in arrays
  lon1(:,1:nlon-1)  = lon
  lon1(:,0)         = lon(:,nlon-2)
  lat1(:,0)         = lat(:,nlon-2)



tarea_m2 = tarea*1d-4 ;m2

scalar_total_transport_pi = pi_vel(0,:,:) * tarea_m2 * 1d-6 ;m3/s to Sv
scalar_total_transport_lgm = lgm_vel(0,:,:) * tarea_m2 * 1d-6  ;m3/s to Sv
scalar_total_transport_lgmd = lgmd_vel(0,:,:) * tarea_m2 * 1d-6 ;m3/s to Sv
scalar_total_transport_pi_tamas = pi_vel_tamas(0,:,:) * tarea_m2 * 1d-6 ;m3/s to Sv
scalar_total_transport_lgm_tamas = lgm_vel_tamas(0,:,:) * tarea_m2 * 1d-6 ;m3/s to Sv
scalar_total_transport_lgm_tamasd = lgm_vel_tamasd(0,:,:) * tarea_m2 * 1d-6 ;m3/s to Sv


scalar_total_transport_pi@lon2d  = lon1
scalar_total_transport_pi@lat2d  = lat1

scalar_total_transport_lgm@lon2d  = lon1
scalar_total_transport_lgm@lat2d  = lat1

scalar_total_transport_lgmd@lon2d  = lon1
scalar_total_transport_lgmd@lat2d  = lat1


scalar_total_transport_pi_tamas@lon2d  = lon1
scalar_total_transport_pi_tamas@lat2d  = lat1

scalar_total_transport_lgm_tamas@lon2d  = lon1
scalar_total_transport_lgm_tamas@lat2d  = lat1

scalar_total_transport_lgm_tamasd@lon2d  = lon1
scalar_total_transport_lgm_tamasd@lat2d  = lat1





res2 = res
res2@tmYLLabelsOn  = False

;plots0 = gsn_csm_contour_map(wks,scalar_total_transport_pi_tamas(1:nlat-1,:),res)
;plots2 = gsn_csm_contour_map(wks,scalar_total_transport_lgm_tamas(1:nlat-1,:),res)
;plots4 = gsn_csm_contour_map(wks,scalar_total_transport_lgm_tamasd(1:nlat-1,:),res)
;plots1 = gsn_csm_contour_map(wks,scalar_total_transport_pi(1:nlat-1,:),res2)
;plots3 = gsn_csm_contour_map(wks,scalar_total_transport_lgm(1:nlat-1,:),res2)
;plots5 = gsn_csm_contour_map(wks,scalar_total_transport_lgmd(1:nlat-1,:),res2)

diff1 = abs(scalar_total_transport_pi) - abs(scalar_total_transport_pi_tamas)
diff2 = abs(scalar_total_transport_lgm) - abs(scalar_total_transport_lgm_tamas)
diff3 = abs(scalar_total_transport_lgmd) - abs(scalar_total_transport_lgm_tamasd)




diff1@lon2d  = lon1
diff1@lat2d  = lat1
diff2@lon2d  = lon1
diff2@lat2d  = lat1
diff3@lon2d  = lon1
diff3@lat2d  = lat1


plots0 = gsn_csm_contour_map(wks,diff1(1:nlat-1,:),res)
plots2 = gsn_csm_contour_map(wks,diff2(1:nlat-1,:),res)
plots4 = gsn_csm_contour_map(wks,diff3(1:nlat-1,:),res)
plots1 = gsn_csm_contour_map(wks,diff1(1:nlat-1,:),res2)
plots3 = gsn_csm_contour_map(wks,diff2(1:nlat-1,:),res2)
plots5 = gsn_csm_contour_map(wks,diff3(1:nlat-1,:),res2)


 ;;-- Vector plot resources
  vecres                  = True
  vecres@gsnDraw          = False
  vecres@gsnFrame         = False
  ;vecres@gsnAddCyclic     = True
  vecres@vcGlyphStyle          = "CurlyVector"        ; turn on curly vectors
  vecres@vcRefMagnitudeF  = 10.0              ; define vector ref mag
  vecres@vcRefLengthF     = 0.04           ; define length of vec ref
  vecres@gsnRightString   = " "             ; turn off right string
  vecres@gsnLeftString    = " "             ; turn off left string
  vecres@tiXAxisString    = " "             ; turn off axis label
  vecres@vcLineArrowColor        = "black"           ; change vector color
  vecres@vcLineArrowThicknessF   = 2.0               ; change vector thickne

  vecres@vcRefAnnoOn = False
;vecres@vcRefAnnoOrthogonalPosF = 0.2    ; move ref vector into plot

 ; vecres@vcRefMagnitudeF  = 0.2   ; reference vector length
 ; vecres@vcRefLengthF     = 0.04  ; reference arrow size
  vecres@vcMinDistanceF   = 0.025  ; minimum distance between vectors
 ; vecres@vcGlyphStyle     = "Arrow"
 ; vecres@vcLineArrowColor = "black"

    vecres@mpMaxLatF             = 90.                  ; zoom in on map  
    vecres@mpMinLatF             =  30.  
    vecres@mpMinLonF             = -60.
    vecres@mpMaxLonF             = 30.


  urot_pi_tamas@lon2d  = lon
  urot_pi_tamas@lat2d  = lat
  vrot_pi_tamas@lon2d  = lon
  vrot_pi_tamas@lat2d  = lat

  urot_lgm_tamas@lon2d  = lon
  urot_lgm_tamas@lat2d  = lat
  vrot_lgm_tamas@lon2d  = lon
  vrot_lgm_tamas@lat2d  = lat

  urot_pi@lon2d  = lon
  urot_pi@lat2d  = lat
  vrot_pi@lon2d  = lon
  vrot_pi@lat2d  = lat

  urot_lgm@lon2d  = lon
  urot_lgm@lat2d  = lat
  vrot_lgm@lon2d  = lon
  vrot_lgm@lat2d  = lat

  urot_lgmd@lon2d  = lon
  urot_lgmd@lat2d  = lat
  vrot_lgmd@lon2d  = lon
  vrot_lgmd@lat2d  = lat


  urot_lgm_tamasd@lon2d  = lon
  urot_lgm_tamasd@lat2d  = lat
  vrot_lgm_tamasd@lon2d  = lon
  vrot_lgm_tamasd@lat2d  = lat


  vectors1 = gsn_csm_vector(wks, urot_pi_tamas,  vrot_pi_tamas, vecres)
  vectors2 = gsn_csm_vector(wks, urot_lgm_tamas,  vrot_lgm_tamas, vecres)
  vectors3 = gsn_csm_vector(wks, urot_lgm_tamasd,  vrot_lgm_tamasd, vecres)  
  vectors4 = gsn_csm_vector(wks, urot_pi,  vrot_pi, vecres)
  vectors5 = gsn_csm_vector(wks, urot_lgm,  vrot_lgm, vecres)
  vectors6 = gsn_csm_vector(wks, urot_lgmd,  vrot_lgmd, vecres)






 overlay(plots0, vectors1)
 overlay(plots1, vectors2)
 overlay(plots2, vectors3)
 overlay(plots3, vectors4)
 overlay(plots4, vectors5) 
 overlay(plots5, vectors6)




  plot(0) = plots0
  plot(1) = plots1
  plot(2) = plots2
  plot(3) = plots3
  plot(4) = plots4
  plot(5) = plots5



resP = True
;  resP@gsnMaximize         = True                 ; make as large as possible
;  resP@gsnPaperOrientation = "portrait"           ; force portrait
;  resP@gsnPanelMainString  = "DEPTH of Specified TEMP Levels"
  resP@gsnMaximize        = False
  resP@gsnFrame         = False               ; save panel until both ready
  resP@gsnPaperOrientation = "landscape"
  resP@gsnPanelLabelBar   = True
  resP@pmLabelBarWidthF   = 0.60
  resP@pmLabelBarHeightF = 0.043
;   resP@pmLabelBarWidthF   = 0.07
;   resP@pmLabelBarHeightF = 0.33
   resP@pmLabelBarOrthogonalPosF = -0.042
   ;resP@pmLabelBarParallelPosF = -0.01
;   resP@lbTitleExtentF = 0.01
;resP@lbOrientation        = "Vertical"
;resP@lbLabelPosition      = "Right"   


;resP@lbTitleString = "Difference in Potential Temperature [degC]"
resP@lbTitleFontHeightF = 0.011

resP@lbLabelFontHeightF= 0.010
resP@gsnPanelFigureStringsFontHeightF = 0.011
resP@gsnPanelFigureStrings= (/"(a)", "(b)" , "(c)", "(d)", "(e)", "(f)"/)

;resP@gsnPanelTop = 0.95
;resP@gsnPanelBottom = 0.61



resP@gsnPanelYWhiteSpacePercent = 3
                                                     

gsn_panel(wks, plot, (/3,2/), resP)


txres2 = True
txres2@txFontHeightF = 0.014
txres2@txFont   = "helvetica-bold"

;gsn_text_ndc(wks,"PI (w/IDEMIX-w/oIDEMIX)",0.27,0.96,txres2)
;gsn_text_ndc(wks,"LGM (w/IDEMIX-w/oIDEMIX)",0.77,0.96,txres2)

;gsn_text_ndc(wks,"LGM - PI (w/oIDEMIX)",0.27,0.99,txres2)
;gsn_text_ndc(wks,"LGM - PI (w/IDEMIX)",0.77,0.99,txres2)


;gsn_text_ndc(wks,"(Di-Nezio et al., 2016)",0.77,0.96,txres2)
;gsn_text_ndc(wks,"(Wilmes et al., 2021)",0.27,0.96,txres2)


frame(wks)
end

