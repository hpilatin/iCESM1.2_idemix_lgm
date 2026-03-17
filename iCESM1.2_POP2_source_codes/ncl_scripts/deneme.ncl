
begin

  
   ;pop_lgm_tamas = addfile ("/mnt/d/lgm_data/buo_test_cmake.pop.h.1550_1599.nc", "r")
   
   ;pop_pi = addfile ("/mnt/d/lgm_data/buo_pi_idemix_isotope_continued_4900_4999_timmean.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data/buo_lgm_idemix_isotope_continued_5800_5899_timmean.nc", "r"

 pop_pi_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")
pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")


   pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")

wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"deneme")


   cmap = read_colormap_file("cmocean_delta")
   ncmap = dimsizes(cmap(:,0))


pltres = True
mpres = True
pltres@PanelPlot  = True

plot = new(4,graphic)


;==============================================================
; Compute MOC streamfunction Ψ(y,z) from POP2 output variables
;==============================================================


  ;--- Read required variables
  vvel1  = pop_pi_tamas->WVEL(0,:,:,:) ; [z_t, nlat, nlon]  or [time, z_t, nlat, nlon]
  vvel2  = pop_lgm_tamas->WVEL(0,:,:,:)
  vvel3  = pop_pi->WVEL(0,:,:,:)
  vvel4  = pop_lgm-WVEL(0,:,:,:)


  tlat = pop_pi_tamas->TLAT               ; [nlat, nlon]
  tlong= pop_pi_tamas->TLONG              ; [nlat, nlon]
  z_t2  = pop_pi_tamas->z_t                ; vertical grid (cm or m)
  mask1 = pop_pi_region->REGION_MASK        ; region mask (to isolate basins)
  mask2 = pop_lgm_region->REGION_MASK        ; region mask (to isolate basins)



  z_t = z_t2*0.01 ; from cm to m

  ;--- Convert depth to meters if needed 
  ;if (max(z_t) > 2000) then 
  depth =z_t / 100.0 ;else
  ;  depth = z_t
  ;end if
  ;depth@units = "m"

  ;--- Select Indo-Pacific basin
  ; In POP2, region masks often: 1=Atlantic, 2=Pacific, 3=Indian, 6=Global
  ; Combine Pacific (2) + Indian (3) for Indo-Pacific
  ;indo_mask1 = where(mask1.eq.2 .or. mask1.eq.3, 1, 0)
  ;indo_mask2 = where(mask2.eq.2 .or. mask2.eq.3, 1, 0)




  ;--- If multiple time steps, choose one or take mean
  if (dimsizes(dimsizes(vvel1)) .eq. 4) then
    vvel1 = dim_avg_Wrap(vvel1(time|:,z_t|:, nlat|:, nlon|:))
  end if

   if (dimsizes(dimsizes(vvel2)) .eq. 4) then
    vvel2 = dim_avg_Wrap(vvel2(time|:,z_t|:, nlat|:, nlon|:))
  end if

 if (dimsizes(dimsizes(vvel3)) .eq. 4) then
    vvel3 = dim_avg_Wrap(vvel3(time|:,z_t|:, nlat|:, nlon|:))
  end if

 if (dimsizes(dimsizes(vvel4)) .eq. 4) then
   vvel4  = dim_avg_Wrap(vvel4(time|:,z_t|:, nlat|:, nlon|:))
  end if


region_po_pi= new((/384,320/),float)
region_po_lgm= new((/384,320/),float)

region_po_pi        = where(mask1.eq.2,1,9.969209968386869e+36)
region_po_lgm        = where(mask2.eq.2,1,9.969209968386869e+36)

pi_vvel_tamas_po = mask(vvel1,region_po_pi,1)
lgm_vvel_tamas_po = mask(vvel2,region_po_lgm,1)
pi_vvel_po = mask(vvel3,region_po_pi,1)
lgm_vvel_po = mask(vvel4,region_po_lgm,1)



vvel1_indo    = dim_avg_n_Wrap(pi_vvel_tamas_po,(/1,2/))
vvel2_indo    = dim_avg_n_Wrap(lgm_vvel_tamas_po,(/1,2/))
vvel3_indo    = dim_avg_n_Wrap(pi_vvel_po,(/1,2/))
vvel4_indo    = dim_avg_n_Wrap(lgm_vvel_po,(/1,2/))


  ;--- Apply mask
 ; vvel1_indo = vvel1 * indo_mask1
 ; vvel2_indo = vvel2 * indo_mask2
 ; vvel3_indo = vvel3 * indo_mask1
 ; vvel4_indo = vvel4 * indo_mask2

printVarSummary(vvel1_indo)
printVarSummary(vvel4_indo)

  ;--- Zonal sum (integrate along longitude)
  ; dim syntax: assume vmo(z, y, x)
  vvel1_zonal = dim_sum_n_Wrap(vvel1_indo, 2)   ; integrate over x
  vvel2_zonal = dim_sum_n_Wrap(vvel2_indo, 2)   ; integrate over x
  vvel3_zonal = dim_sum_n_Wrap(vvel3_indo, 2)   ; integrate over x
  vvel4_zonal = dim_sum_n_Wrap(vvel4_indo, 2)   ; integrate over x


  ;--- Now integrate vertically from bottom to top
  nz = dimsizes(depth)
  nlat = dimsizes(vvel1_zonal(1,:))
  psi1 = new((/nz, nlat/), typeof(vvel1_zonal))
  psi1 = 0.0

  psi2 = new((/nz, nlat/), typeof(vvel2_zonal))
  psi2 = 0.0  

  psi3 = new((/nz, nlat/), typeof(vvel3_zonal))
  psi3 = 0.0


  psi4 = new((/nz, nlat/), typeof(vvel4_zonal))
  psi4 = 0.0


  do k = 0, nz-1
    if (k.eq.0) then
      psi1(k,:) = vvel1_zonal(k,:)
    else
      psi1(k,:) = psi1(k-1,:) + vvel1_zonal(k,:)
    end if
  end do

  do k = 0, nz-1
    if (k.eq.0) then
      psi2(k,:) = vvel2_zonal(k,:)
    else
      psi2(k,:) = psi2(k-1,:) + vvel2_zonal(k,:)
    end if
  end do


  do k = 0, nz-1
    if (k.eq.0) then
      psi3(k,:) = vvel3_zonal(k,:)
    else
      psi3(k,:) = psi3(k-1,:) + vvel3_zonal(k,:)
    end if
  end do


  do k = 0, nz-1
    if (k.eq.0) then
      psi4(k,:) = vvel4_zonal(k,:)
    else
      psi4(k,:) = psi4(k-1,:) + vvel4_zonal(k,:)
    end if
  end do

 ;--- Convert to Sverdrups
  psi1 = psi1 / 1e6
  psi1@long_name = "Meridional Overturning Streamfunction"
  psi1@units = "Sv"

  psi2 = psi2 / 1e6
  psi2@long_name = "Meridional Overturning Streamfunction"
  psi2@units = "Sv"


  psi3 = psi3 / 1e6
  psi3@long_name = "Meridional Overturning Streamfunction"
  psi3@units = "Sv"


  psi4 = psi4 / 1e6
  psi4@long_name = "Meridional Overturning Streamfunction"
  psi4@units = "Sv"

  ;--- Plot section (latitude–depth)
  ;wks = gsn_open_wks("pdf","MOC_IndoPacific_computed")

  res                     = True
  res@gsnMaximize         = True
  res@cnFillOn            = True
  res@cnFillPalette       = "MPL_BrBG"
  res@cnLinesOn           = False
  res@cnLevelSpacingF     = 2
  res@tiMainString        = "Computed Indo-Pacific MOC"
  res@tiXAxisString       = "Latitude (°)"
  res@tiYAxisString       = "Depth (m)"
  res@trYReverse          = True
  res@sfXArray            = tlat(:,0)
  res@sfYArray            = depth
  res@lbLabelBarOn        = True
  res@gsnRightString      = "Ψ [Sv]"

 
plot(0) = gsn_csm_contour(wks,psi1,res)     ; create plot
plot(1) = gsn_csm_contour(wks,psi2,res)     ; create plot
plot(2) = gsn_csm_contour(wks,psi3,res)     ; create plot
plot(3) = gsn_csm_contour(wks,psi4,res)     ; create plot


resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = False
resP@pmLabelBarWidthF   = 0.07
resP@pmLabelBarHeightF = 0.33
resP@lbOrientation   = "vertical"
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





resP@lbTitleString = "log10[m^2/s]"
resP@lbTitleFontHeightF = 0.010
resP@lbTitlePosition     = "Right"
resP@lbTitleDirection     = "Across"                          ; letter angle
resP@lbTitleAngleF        = 270.                               ; title angle
resP@lbTitleOffsetF = 0.15



;resP@gsnPanelFigureStringsFontHeightF = 0.012
;resP@gsnPanelYWhiteSpacePercent = 0.2
;resP@gsnPanelXWhiteSpacePercent = 0.4


;resP@gsnPanelFigureStrings= (/"PI - w/oIDEMIX", "LGM - w/oIDEMIX", "PI - w/IDEMIX", "LGM - w/IDEMIX"/)
resP@gsnPanelFigureStrings= (/"W", "D", "W", "D", "W", "D"/)
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

;resP@gsnPanelYF = (/0.92,0.92,0.65,0.65,0.38,0.38/)
;!resP@gsnPanelXF = (/0.07,0.51,0.07,0.51,0.07,0.51/)

;resP@gsnPanelXF = (/0.12,0.51,0.12,0.51,0.12,0.51/)


  gsn_panel(wks,plot,(/2,2/),resP)


 ; draw(plot)
 frame(wks)
end

