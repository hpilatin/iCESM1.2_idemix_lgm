
begin


   ;pop_lgm_tamas = addfile ("/mnt/d/lgm_data/buo_test_cmake.pop.h.1550_1599.nc", "r")

   ;pop_pi = addfile ("/mnt/d/lgm_data/buo_pi_idemix_isotope_continued_4900_4999_timmean.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data/buo_lgm_idemix_isotope_continued_5800_5899_timmean.nc", "r"

 ;pop_pi_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")
pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")


   pop_pi_tamas = addfile ("/mnt/d/lgm_data/streamfunctionpiControl.p875.001.pop.h.4100_4199.deflated2_timmean.nc", "r")
   


   pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")


pi_tarea = addfile ("/mnt/d/lgm_data/pi_tarea.nc", "r")
pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")

pi_kmt = addfile ("/mnt/d/lgm_data/kmt_pi.nc", "r")
pi_grid = addfile ("/mnt/d/lgm_data/pi_grid_spacing_and_thickness.nc", "r")


wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"deneme5")

   
   cmap = read_colormap_file("temp_diff_18lev")
   ;cmap = read_colormap_file("cmocean_delta")
   ncmap = dimsizes(cmap(:,0))

 res                      = True                 ; plot mods desired
  res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = cmap(:,:)
  ;res@cnFillPalette        = cmap(::-1,:)         ; set color map
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
plot = new(4,graphic)
;plots = new(4,graphic)                                                        

 ; Read variables
  VVEL = pop_pi_tamas->VVEL(0,:,:,:)      ; [z,lat,lon], take first timestep
  DXT  = pi_grid->DXT(:,:)           ; [lat,lon], zonal grid spacing (m)
  dz  = pi_grid->dz        ; [z], layer thickness (m)

VVEL2 = pop_pi_tamas->VVEL(:,:,:,:)      ; [z,lat,lon], take first timestep



  print(max(VVEL))
  print(min(VVEL))

  z_t           = pop_pi_tamas->z_t
z_t           = z_t*1.E-2       ; convert from cm to m
;z_t@long_name = "Depth (T grid)"; overwrite long name
z_t@long_name = ""; overwrite long name
z_t@units     = "m"             ; fix units attribute
z_t@tiYAxisPosition = "top"
z_t@tiYAxisFontHeightF = 0.018
;kvmix12_atlantic_log!0    = "z_t"
;kvmix12_atlantic_log&z_t  =  z_t
;kvmix12_atlantic_log!1    = "lat_t"
;kvmix12_atlantic_log&lat_t=  lat_atlantic

  ; Create 3D version of dz (broadcast)
  nz   = dimsizes(VVEL(:,0,0))
  nlat = dimsizes(VVEL(0,:,0))
  nlon = dimsizes(VVEL(0,0,:))

  dz3 = new((/nz,nlat,nlon/), typeof(dz))
  do k = 0, nz-1
    dz3(k,:,:) = dz(k)
  end do

  ; Compute cross-sectional area (m²)
  area = DXT            ; [lat,lon]
  area3 = new((/nz,nlat,nlon/), typeof(area))
  do k = 0, nz-1
    area3(k,:,:) = area
  end do

  ; Compute meridional volume transport (m³/s)
  vmo_calc = VVEL * area3 * dz3 * 1d-6; cm3/s to m3/s

   

  ; Convert to Sverdrups
  vmo_calc = vmo_calc / 1e6
  vmo_calc@units = "Sv"

  printVarSummary(vmo_calc)
  print(max(vmo_calc))
  print(min(vmo_calc))


  ; Zonal integration (sum over lon)
  vmo_zonal = dim_sum_n_Wrap(vmo_calc, 2)    ; (z,lat)
  
 ;vmo_zonal = dim_avg_n_Wrap(vmo_calc, 2)    ; (z,lat)

; --- assume vmo_zonal has dims (z,lat)
nz = dimsizes(vmo_zonal(:,0))
nlat = dimsizes(vmo_zonal(0,:))

psi = new((/nz,nlat/), typeof(vmo_zonal))
psi = 0.0

if (nz .eq. 1) then
  ; single level case: psi is just that level
  psi(0,:) = vmo_zonal(0,:)
else
  ; set bottom level (index nz-1)
  psi(nz-1,:) = vmo_zonal(nz-1,:)

  ; positive-stride loop: iterate ip = 1..(nz-1), compute k = nz-1-ip
  do ip = 1, nz-1
    k = nz-1 - ip
    ; psi(k,:) = psi(k+1,:) + vmo_zonal(k,:)
    psi(k,:) = psi(k+1,:) + vmo_zonal(k,:)
  end do
end if

; convert to Sv if vmo_zonal was in m3/s
;psi = psi / 1e6
psi@units = "Sv"

  ; Vertical cumulative sum from bottom to surface
;  psi = new(dimsizes(vmo_zonal), typeof(vmo_zonal))
;  nz = dimsizes(vmo_zonal(:,0))
;  psi(nz-1,:) = vmo_zonal(nz-1,:)
;  do k = nz-2, 0, -1
;    psi(k,:) = psi(k+1,:) + vmo_zonal(k,:)
;  end do

;  psi@units = "Sv"
  printVarSummary(psi)
  print(avg(psi))
  printMinMax(psi, True)

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
;res1@cnFillPalette        = cmap(::-1,:)
res1@cnFillPalette        = cmap(:,:)
;res1@cnLevelSelectionMode = "ExplicitLevels"
;res1@cnLevels    = (/-13,-12,-11,-10,-9,-8,-7,-6/)
;res1@tiYAxisString = "Depth (m)"
;res1@cnLevels    = (/-13,-12.5,-12,-11,-11.5,-10,-10.5,-9,-9.5,-8,-8.5,-7,-7.5,-6,-6.5,-5/)  

res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
res1@cnMinLevelValF            = -20             ; min level
res1@cnMaxLevelValF            = 20             ; max level
res1@cnLevelSpacingF           =  5            ; contour interval

ntm    = 0                  ; time subscript
nlt    = 9               ; latitude subscript
;nlt = 10

ntm2   = 0                  ; time subscript
nlt2    = 224               ; latitude subscript

nlt3    = 305 
nlt4 = 312  ; along 50N

lat_a    = pop_lgm->TLAT(:,nlt)
lat_atlantic        = lat_a
;lat_a    = pop_lgm->TLAT(nlt4,:)
;lat_atlantic        = lat_a
lat_p    = pop_lgm->TLAT(:,nlt2)
lat_pacific        = lat_p


latz = pop_pi_tamas->TLAT(:,0)

psi!0    = "z_t"
psi&z_t  =  z_t
psi!1    = "lat_t"
psi&lat_t=  lat_atlantic

  
plot(0) = gsn_csm_contour(wks, psi,res1)
plot(1) = gsn_csm_contour(wks, psi,res1)
plot(2) = gsn_csm_contour(wks, psi,res1)
plot(3) = gsn_csm_contour(wks, psi,res1)


 PD_mean   = pop_pi_tamas->PD(0,:,:,:)  
  PD2   = pop_pi_tamas->PD(:,:,:,:)    

PD2_v2 = (PD2*1000)-1000

rho1 = PD_mean*1000        ; convert to kg/m^3
sigma1 = rho1-1000           ; convert to sigma0 (potential density anomaly)


printMinMax(sigma1, 0)




depth          = pop_pi_tamas->z_t

Tiso = 27.7        ; 27.7°C isotherm

;************************************************
; Interpolate to specific [constant] depths
; The returned vertical coordinate is Z_T
;************************************************
 ; zlev = (/ 735, 12580/)    ; cm
 
  zlev = new((/60/), typeof(depth), -999) 
  nt = dimsizes(depth) 

 do i = 0, nt-1
  zlev(i) = depth(i)    ; cm
 end do  

  print(zlev)

  zlev@units = depth@units
  zlev!0 = "zlev"
  iso1 = int2p_n_Wrap(depth,PD2_v2, zlev, 0, 1)

  printVarSummary(iso1) ; [time | 1] x [zlev | 2] x [nlat | 384] x [nlon | 320]  
  printMinMax(iso1, True) 

   pi_tamas_mean2    = dim_avg_n(iso1(0,:,:,:),2)




lat_a    = pop_pi->TLAT(:,nlt)
lat_atlantic        = lat_a

z_t12           = pop_pi_tamas->z_t
z_t12           = z_t12*1.E-2       ; convert from cm to m

z_t12@units     = "m"             ; fix units attribute
z_t12@tiYAxisPosition = "top"
z_t12@tiYAxisFontHeightF = 0.018
z_t12@long_name = "Depth (m)"; overwrite long name


pi_tamas_mean2!0    = "z_t"
pi_tamas_mean2&z_t  =  z_t
pi_tamas_mean2!1    = "lat_t"
pi_tamas_mean2&lat_t=  latz



overlay_res                          = True
  overlay_res@cnInfoLabelOn       = False
  overlay_res@description = "  "
  overlay_res@gsnDraw              = False                ; don't draw yet
  overlay_res@gsnFrame             = False                ; don't advance frame yet
  overlay_res@cnLineDrawOrder         = "PostDraw"
  overlay_res@cnLinesOn               = True
  overlay_res@cnMonoLineColor  = True 
  ;!overlay_res@cnLineColors             = span_color_rgba ("BlAqGrYeOrReVi200",8)

;  overlay_res@cnLineColors             = "WhViBlGrYeOrRe"
 overlay_res@cnLineColor             = "black"
  overlay_res@cnLineThicknessF        = 1.5
  overlay_res@tiMainString  = " "             ; add title
  overlay_res@cnLineDashPattern = 0
  overlay_res@gsnYAxisIrregular2Linear = True
  overlay_res@trYReverse           = True       ; reverses y-axis
  ;overlay_res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
   ; overlay_res@cnMinLevelValF            =27.7             ; min level
   ; overlay_res@cnMaxLevelValF            = 27.7             ; max level
   ; overlay_res@cnLevelSpacingF           =  0            ; contour interval
 
  overlay_res@cnLevelSelectionMode      = "ExplicitLevels"
  overlay_res@cnLevels = (/27.7/)
  overlay_res@cnLineLabelsOn = False

;pi_tamas_mean2!0    = "z_t"
;pi_tamas_mean2&z_t  =  z_t
;pi_tamas_mean2!1    = "lat_t"
;pi_tamas_mean2&lat_t=  lat_atlantic



  ax1 = gsn_csm_contour(wks,pi_tamas_mean2,overlay_res)
  ax2 = gsn_csm_contour(wks,pi_tamas_mean2,overlay_res)
  ax3 = gsn_csm_contour(wks,pi_tamas_mean2,overlay_res)
  ax4 = gsn_csm_contour(wks,pi_tamas_mean2,overlay_res)

    overlay(plot(0), ax1)
    overlay(plot(1), ax2)
    overlay(plot(2), ax3)
    overlay(plot(3), ax4)



resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = True
; resP@pmLabelBarWidthF   = 0.60
;resP@pmLabelBarHeightF = 0.048
resP@lbOrientation   = "Horizontal"   
resP@pmLabelBarOrthogonalPosF = -0.025
resP@lbTitleString = " Averaged Streamfunctions of Overturning [Sv]"
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

txres5 = True
txres5@txFontHeightF = 0.012
  txres5@txAngleF      = 90.     ; rotate 90 degrees
;  txres5@txJust        = "TopLeft"

gsn_text_ndc(wks,"Depth (m)",0.02,0.55,txres5)

gsn_panel(wks,plot,(/2,2/),resP)

 
; now we draw the whole thing and advance the frame.


  frame(wks)
  end


