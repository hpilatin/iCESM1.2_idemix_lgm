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




; pop_pi_tamas = addfile ("/home/hbkheves/from_Tamas/check_o2_moc_pi_tamas_restart_files/b.ie12.Bi1850C5BPRDCISO.f19_g16.piControl.p875.001_MOC_mean_2850-2899.nc", "r")


;   pop_lgm_tamas = addfile ("/home/hbkheves/from_Tamas/check_o2_moc_pi_tamas_restart_files/230913_LGM_means_0950-1000/b.ie12.Bi1850C5BPRDCISO.f19_g16.LGM21ka_frUsp.p875.001_MOC_mean_0950-1000.nc","r")


;---Open files using addfiles
;   pop_lgm5 = addfile ("/scratch/usr/hbkheves/work/CESM/archive/b.ie12.Bi1850C5BPRDCISO.f19_g16_LGM21ka_ciso_idemix_j5_year1000/ocn/hist/lgm_1700_1799_moc.nc", "r")

 
   pop_pi_tamas = addfile ("/mnt/d/lgm_data/moc_pi_idemix_3900_4604_4703.nc", "r")
   
    ;pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/moc_lgm_idemix_isotope_continued_4900_5899.nc", "r")
  
   pop_lgm_tamas = addfile ("/mnt/d/lgm_data/moc_lgm_idemix_isotope_continued_5800_5899.nc", "r")
   

   pop_pi = addfile ("/mnt/d/lgm_data/moc_pi_idemix_3900_4604_4703.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data2/moc_lgm_idemix_isotope_continued_4900_5899.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data2/moc_lgm_idemix_isotope_continued_5944_5993.nc", "r")
  
   
   ;pop_lgm = addfile ("/mnt/d/lgm_data/moc_lgm_year1000.1000_2000.nc", "r")

   pop_lgm = addfile ("/mnt/d/lgm_data/moc_lgm_idemix_isotope_continued_5800_5899.nc", "r")

   pop_lgm_coord = addfile ("/mnt/d/lgm_data2/select_test_cmake_1670_1719_timmean.nc", "r")



   pop_pi_new = addfile ("/mnt/d/lgm_data/moc_pi_idemix_3900_4604_4703.nc", "r")

 
   pi_tamas = addfile ("/mnt/d/lgm_data2/b.ie12.Bi1850C5BPRDCISO.f19_g16.piControl.p875.001_MOC_mean_4150-4199.nc", "r")
   lgm_tamas = addfile ("/mnt/d/lgm_data2/hmoc_test_cmake.pop.h.1830_1879.nc", "r")



pop_pi_tamas_no = addfile ("/mnt/d/lgm_data2/moc_pi_tamas_last100_timmean.nc", "r")                                      
pop_lgm_wilmes_no = addfile ("/mnt/d/lgm_data2/moc_test_cmake.pop.h.1780_1879_timmean.nc","r")
pop_pi_ide = addfile ("/mnt/d/lgm_data/moc_pi_idemix_3900_4604_4703_timmean.nc", "r")
pop_lgm_wilmes = addfile ("/mnt/d/lgm_data/moc_lgm_idemix_isotope_continued_5800_5899_timmean.nc", "r")
pop_nc_example1 = addfile ("/mnt/d/lgm_data/kvmix_pi_tamas_4190_4199_timmean.nc", "r")

  pop_lgm_dinezio_no = addfile ("/mnt/d/lgm_data2/moc_lgm_tamas_last100_timmean.nc","r")
  pop_lgm_dinezio = addfile ("/mnt/d/lgm_data2/moc_dinezio_continued.pop.h.1760_1859_timmean.nc", "r") 





 
   wks_type = "pdf"
  ; wks_type@wkPaperSize = "A4"


   wks  = gsn_open_wks(wks_type,"nadw_aabw_pi_lgm_timeseries_thesis_lastnew")


 ;wks  = gsn_open_wks("png","moc_pi_lgm_2023_bathymetry")


   cmap = read_colormap_file("MPL_RdYlGn")
   ncmap = dimsizes(cmap(:,0))
 

res = True
pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False

;plot = new(4,graphic)
plots = new(3,graphic)

res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "


nltp    = 224               ; latitude subscript
ntm    = 0                  ; time subscript
nlt    = 9               ; latitude subscript
lat    = pop_lgm_coord->TLAT(nlt,:)
latavg = avg(lat)
lon_o  = pop_lgm_coord->TLONG(nlt,:)
nlon   = dimsizes (lon_o)
lon_t        = lon_o          ; trick to conserve meta data
lon_t@units  = "degrees_east"
lon_t( 0:283) = lon_o(36:319)
lon_t(284:319) = lon_o( 0:35)
ii           = ind(lon_o.GT. 320.0)
lon_t(ii)    = 360.-lon_o(ii)
lonavg = lon_t(nlt-1)


lon_oh  = pop_lgm_coord->TLONG(:,:)
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

;print(lon_th(:,53))



mc1_new         = pop_pi_new->MOC(90:99,1,0,33:60,:)
amoc_pi_v2_26n_new = mc1_new(:,:,274) ; near 26N

ntime_pi1_new = dimsizes(amoc_pi_v2_26n_new(:,0))
mc_pi_time_26n_new = new(ntime_pi1_new,float)


do i = 0,ntime_pi1_new-1
mc_pi_time_26n_new(i)= max(amoc_pi_v2_26n_new(i,:))
end do

;amoc_pi_26n_mean = avg(mc_pi_time_26n_new)

amoc_pi_26n_mean = avg(amoc_pi_v2_26n_new)


amoc_pi_v2_32s_new = mc1_new(:,:,89) ; near 32S

amoc_pi_32s_mean = avg(amoc_pi_v2_32s_new)



mc_pi_tamas = pi_tamas->MOC(:,1,0,33:60,:)

amoc1_pi_v2_26n_tamas = mc_pi_tamas(:,:,274) ; near 26N
ntime_pi1_tamas = dimsizes(amoc1_pi_v2_26n_tamas(:,0))
mc1_time_26n_tamas = new(ntime_pi1_tamas,float)


do i = 0,ntime_pi1_tamas-1
mc1_time_26n_tamas(i) = max(amoc1_pi_v2_26n_tamas(i,:))
end do



amoc1_pi_v2_20s_tamas = mc_pi_tamas(:,:,115)
mc1_time_20s_tamas = new(ntime_pi1_tamas,float)

do i = 0,ntime_pi1_tamas-1
mc1_time_20s_tamas(i) = min(amoc1_pi_v2_20s_tamas(i,:))
end do






mc_lgm_tamas = lgm_tamas->MOC(:,1,0,33:60,:)

amoc6_lgm_v2_26n_tamas = mc_lgm_tamas(:,:,274) ; near 26N
ntime_lgm1_tamas = dimsizes(amoc6_lgm_v2_26n_tamas(:,0))
mc6_time_26n_tamas = new(ntime_lgm1_tamas,float)


do i = 0,ntime_lgm1_tamas-1
mc6_time_26n_tamas(i) = max(amoc6_lgm_v2_26n_tamas(i,:))
end do


amoc6_lgm_v2_20s_tamas = mc_lgm_tamas(:,:,115) ; near 20S
mc6_time_20s_tamas = new(ntime_lgm1_tamas,float)


do i = 0,ntime_lgm1_tamas-1
mc6_time_20s_tamas(i) = min(amoc6_lgm_v2_20s_tamas(i,:))
end do




mc1         = pop_pi_tamas->MOC(:,1,0,33:60,:)


;mc1_mask = mask(mc1,(mc1.gt.0.000001.or.mc1.lt.-0.00001),True)
;printVarSummary(mc1_mask)


;lat1p    = pop_pi[0]->TLAT(:,nltp)
;lat_t1p        = lat1p
;lat_aux_grid6p          = mc1_mask&lat_aux_grid
;z_1tp           = pop_pi[0]->moc_z
;z_1tp           = z_1tp*1.E-2       ; convert from cm to m
;z_1tp@long_name = ""; overwrite long name
;z_1tp@units     = "m"             ; fix units attribute

;amoc1 = dim_avg_n(mc1_mask,0)
;printVarSummary(amoc1)

;amoc1!0    = "moc_z"
;amoc1&moc_z  =  z_1tp
;amoc1!1    = "lat_aux_grid"
;amoc1&lat_aux_grid= mc1&lat_aux_grid






amoc1_pi_v2_26n = mc1(:,:,274) ; near 26N
ntime_pi1 = dimsizes(amoc1_pi_v2_26n(:,0))
mc1_time_26n = new(ntime_pi1,float)


do i = 0,ntime_pi1-1
mc1_time_26n(i) = max(amoc1_pi_v2_26n(i,:))
end do
;print(mc1_time_26n)


;mc1_time_26n_cycle = new(50,float)
;do i = 0,49
;a = avg(mc1_time_26n((20*i):(20*i+18)))
;mc1_time_26n_cycle(i) = a
;end do



amoc1_pi_v2_20s = mc1(:,:,115)
mc1_time_20s = new(ntime_pi1,float)

do i = 0,ntime_pi1-1
mc1_time_20s(i) = min(amoc1_pi_v2_20s(i,:))
end do
;print(mc1_time_20s)



amoc1_pi_v2_32s = mc1(:,:,89)
mc1_time_32s = new(ntime_pi1,float)

do i = 0,ntime_pi1-1
mc1_time_32s(i) = max(amoc1_pi_v2_32s(i,:))
end do
;print(mc1_time_32s)

printVarSummary(mc1_time_32s)

;mc1_time_32s_cycle = new(50,float)
;do i = 0,49
;a = avg(mc1_time_32s((20*i):(20*i+18)))
;mc1_time_32s_cycle(i) = a
;end do

;printVarSummary(mc1_time_32s_cycle)

mc2         = pop_lgm_tamas->MOC(:,1,0,33:60,:)


;mc2_mask = mask(mc2,(mc2.gt.0.000003.or.mc2.lt.-0.00003),True)
;printVarSummary(mc2_mask)


;lat2p    = pop_pi[0]->TLAT(:,nltp)
;lat_t2p        = lat2p
;;lat_aux_grid2p          = mc2_mask&lat_aux_grid
;z_2tp           = pop_pi[0]->moc_z
;z_2tp           = z_2tp*1.E-2       ; convert from cm to m
;z_2tp@long_name = ""; overwrite long name
;z_2tp@units     = "m"             ; fix units attribute

;amoc2 = dim_avg_n(mc2_mask,0)
;printVarSummary(amoc2)

;amoc2!0    = "moc_z"
;amoc2&moc_z  =  z_2tp
;amoc2!1    = "lat_aux_grid"
;amoc2&lat_aux_grid= mc2&lat_aux_grid




amoc1_lgm_v2_26n = mc2(:,:,274) ; near 26N
ntime_lgm1 = dimsizes(amoc1_lgm_v2_26n(:,0))
mc2_time_26n = new(ntime_lgm1,float)

do i = 0,ntime_lgm1-1
mc2_time_26n(i) = max(amoc1_lgm_v2_26n(i,:))
end do
;print(mc2_time_26n)



amoc1_lgm_v2_20s = mc2(:,:,115)
mc2_time_20s = new(ntime_lgm1,float)

do i = 0,ntime_lgm1-1
mc2_time_20s(i) = min(amoc1_lgm_v2_20s(i,:))
end do
;print(mc2_time_20s)



ntmp    = 0                  ; time subscript
nltp    = 224               ; latitude subscript
latp    = pop_lgm_coord->TLAT(nltp,:)
latavgp = avg(latp)
lon_op  = pop_lgm_coord->TLONG(nltp,:)
nlonp   = dimsizes (lon_op)
lon_tp        = lon_op          ; trick to conserve meta data
lon_tp@units  = "degrees_east"
lon_tp( 0:283) = lon_op(36:319)
lon_tp(284:319) = lon_op( 0:35)
iip           = ind(lon_op.GT. 320.0)
lon_tp(iip)    = 360.-lon_op(iip)
lonavgp = lon_tp(nltp-1)


lon_ohp  = pop_lgm_coord->TLONG(:,:)
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


mc4         = pop_pi->MOC(:,1,0,33:60,:) ;below 500m

;amoc1 = dim_avg_n(mc4,0)
mc4_mask = mask(mc4,(mc4.gt.0.00001.or.mc4.lt.-0.00001),True)
printVarSummary(mc4_mask)


;z_tp4           = pop_pi[0]->moc_z
;z_tp4           = z_tp4*1.E-2 


;mc4_mask!1    = "moc_z"
;mc4_mask&moc_z  =  z_tp4
;mc4_mask!2    = "lat_aux_grid"

amoc4 = dim_avg_n(mc4_mask,0)
;printVarSummary(amoc1)


;lat4p    = pop_pi[0]->TLAT(:,nltp)
;lat_t4p        = lat4p
;lat_aux_grid4p          = mc4&lat_aux_grid
;z_4tp           = pop_pi[0]->moc_z
;z_4tp           = z_4tp*1.E-2       ; convert from cm to m
;z_4tp@long_name = ""; overwrite long name
;z_4tp@units     = "m"             ; fix units attribute
;amoc4!0    = "moc_z"
;amoc4&moc_z  =  z_4tp
;amoc4!1    = "lat_aux_grid"
;amoc4&lat_aux_grid= mc4&lat_aux_grid



;amoc_pi_v2_25n = mc4(:,:,272) ;near 25N

amoc_pi_v2_26n = mc4(:,:,274) ; near 26N
ntime_pi = dimsizes(amoc_pi_v2_26n(:,0))
mc4_time_26n = new(ntime_pi,float)

do i = 0,ntime_pi-1
mc4_time_26n(i) = max(amoc_pi_v2_26n(i,:))
end do
;print(mc4_time_26n)






amoc_pi_v2_40n = mc4(:,:,301)
mc4_time_40n = new(ntime_pi,float)

do i = 0,ntime_pi-1
mc4_time_40n(i) = max(amoc_pi_v2_40n(i,:))
end do
;print(mc4_time_40n)


amoc_pi_v2_20s = mc4(:,:,115)
mc4_time_20s = new(ntime_pi,float)

do i = 0,ntime_pi-1
mc4_time_20s(i) = min(amoc_pi_v2_20s(i,:))
end do
;print(mc4_time_20s)



amoc_pi_v2_32s = mc4(:,:,89)
mc4_time_32s = new(ntime_pi,float)

do i = 0,ntime_pi-1
mc4_time_32s(i) = min(amoc_pi_v2_32s(i,:))
end do
;print(mc4_time_32s)



;mc4_time_32s_cycle = new(50,float)
;do i = 0,49
;a = avg(mc4_time_32s((20*i):(20*i+18)))
;mc4_time_32s_cycle(i) = a
;end do


;print(max(mc4_time_32s_cycle))
;print(min(mc4_time_32s_cycle))

;print(mc4_time_32s_cycle(0))
;print(mc4_time_32s_cycle(49))

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




mc6         = pop_lgm->MOC(:,1,0,33:60,:) ; below 500m
;mc6_mask = mask(mc6,(mc6.gt.0.000085.or.mc6.lt.-0.00001),True)
mc6_mask = mask(mc6,(mc6.gt.0.000003.or.mc6.lt.-0.00003),True)
printVarSummary(mc6_mask)


;lat6p    = pop_pi[0]->TLAT(:,nltp)
;lat_t6p        = lat6p
;;lat_aux_grid6p          = mc6_mask&lat_aux_grid
;z_6tp           = pop_pi[0]->moc_z
;z_6tp           = z_6tp*1.E-2       ; convert from cm to m
;z_6tp@long_name = ""; overwrite long name
;z_6tp@units     = "m"             ; fix units attribute

amoc6 = dim_avg_n(mc6_mask,0)
printVarSummary(amoc6)

;amoc6!0    = "moc_z"
;amoc6&moc_z  =  z_6tp
;amoc6!1    = "lat_aux_grid"
;amoc6&lat_aux_grid= mc6&lat_aux_grid


;amoc_lgm_v2_25n = mc6(:,:,272) ; near 25N

amoc_lgm_v2_26n = mc6(:,:,274) ; near 26N
ntime_lgm = dimsizes(amoc_lgm_v2_26n(:,0))
mc6_time_26n = new(ntime_lgm,float)

do i = 0,ntime_lgm-1
mc6_time_26n(i) = max(amoc_lgm_v2_26n(i,:))
end do

print(max(mc6_time_26n))
print(min(mc6_time_26n))


; when the entire data is ready, activate the following section

;mc6_time_26n_cycle = new(50,float)
;do i = 0,49
;a = avg(mc6_time_26n((20*i):(20*i+18)))
;mc6_time_26n_cycle(i) = a
;end do

;print(max(mc6_time_26n_cycle))
;print(min(mc6_time_26n_cycle))


;mc1_new_lgm         = pop_lgm->MOC(990:999,1,0,33:60,:)
mc1_new_lgm         = pop_lgm->MOC(:,1,0,33:60,:)
amoc_lgm_v2_26n_new = mc1_new_lgm(:,:,274) ; near 26N

ntime_lgm1_new = dimsizes(amoc_lgm_v2_26n_new(:,0))
mc_lgm_time_26n_new = new(ntime_lgm1_new,float)


do i = 0,ntime_lgm1_new-1
mc_lgm_time_26n_new(i)= max(amoc_lgm_v2_26n_new(i,:))
end do

;amoc_lgm_26n_mean = avg(mc_lgm_time_26n_new)

amoc_lgm_26n_mean = avg(amoc_lgm_v2_26n_new)




amoc_lgm_v2_40n = mc6(:,:,301)
mc6_time_40n = new(ntime_lgm,float)

do i = 0,ntime_lgm-1
mc6_time_40n(i) = max(amoc_lgm_v2_40n(i,:))
end do
;print(mc6_time_40n)

mc6_time_40n_cycle = new(50,float)
do i = 0,49
mc6_time_40n_cycle(i) = avg(amoc_lgm_v2_40n(20i:20i+19,:))
end do



amoc_lgm_v2_20s = mc6(:,:,115)
mc6_time_20s = new(ntime_lgm,float)

do i = 0,ntime_lgm-1
mc6_time_20s(i) = min(amoc_lgm_v2_20s(i,:))
end do
;print(mc6_time_20s)



amoc_lgm_v2_32s = mc6(:,:,89)
mc6_time_32s = new(ntime_lgm,float)

do i = 0,ntime_lgm-1
mc6_time_32s(i) = min(amoc_lgm_v2_32s(i,:))
end do
;print(mc6_time_32s)

printVarSummary(mc6_time_32s)

;mc6_time_32s_cycle = new(50,float)
;do i = 0,49
;a = avg(mc6_time_32s((20*i):(20*i+18)))
;mc6_time_32s_cycle(i) = a
;end do


amoc_lgm_v2_32s_new = mc1_new_lgm(:,:,89) ; near 32S

amoc_lgm_32s_mean = avg(amoc_lgm_v2_32s_new)

;printVarSummary(mc6_time_32s_cycle)















;mc4_time = dim_avg_n_Wrap(mc4,(/1,2/))
;mc6_time = dim_avg_n_Wrap(mc6,(/1,2/))


print(avg(amoc4))

res                      = True       ; plot mods desired
res@gsnDraw              = False                ; don't draw yet
res@gsnFrame             = False                ; don't advance frame yet
res@lbLabelBarOn        = False
res@cnFillOn             = True       ; turn on color fill
;res@cnFillPalette        = "cmocean_balance" ; set color map
;res@cnFillPalette        = cmap(::-1,:)
res@cnFillPalette        = "cmp_b2r"
res@cnLinesOn            = False      ; turn off contour lines
res@cnLinesOn            = True      ; turn off contour lines
res@cnLineLabelsOn       = True      ; turns off contour line labels
;res@cnMinLevelValF       = -4
;res@cnMaxLevelValF       =  20
res@cnLevelSpacingF      =  2
res@cnLabelMasking        = True   
res@cnLineLabelDensityF = 2.0
res@cnLineLabelAngleF = 0.0
;res@cnLineLabelPlacementMode = "computed"
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
res@cnMinLevelValF = -22. ; set min contour level
res@cnMaxLevelValF = 22. ; set max contour level
res@cnLevelSpacingF = 2 ; set contour spacing


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



res1                      = True       ; plot mods desired
res1@gsnDraw              = False                ; don't draw yet
res1@gsnFrame             = False                ; don't advance frame yet
res1@lbLabelBarOn        = False
res1@cnFillOn             = True       ; turn on color fill
;res1@cnFillPalette        = cmap(::-1,:)
;res1@cnFillPalette        = "MPL_RdYlGn" ; set color map
res1@cnFillPalette        = "cmp_b2r"
res1@cnLinesOn            = True      ; turn off contour lines
res1@cnLineLabelsOn       = True      ; turns off contour line labels
;res1@cnMinLevelValF       = -4
;res1@cnMaxLevelValF       =  20
res1@cnLevelSpacingF      =  2 
res1@cnLineLabelAngleF = 0.0
;res1@cnLineLabelPlacementMode = "computed"
res1@cnLineLabelBackgroundColor = "white"
res1@cnLabelMasking        = True   
res1@cnLineLabelFontHeightF = 0.016
res1@cnLineLabelDensityF = 2.0
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

res1@cnLevelSelectionMode = "ManualLevels" ; manual contour levels
res1@cnMinLevelValF = -22. ; set min contour level
res1@cnMaxLevelValF = 22. ; set max contour level
res1@cnLevelSpacingF = 2 ; set contour spacing



;time = ispan (1,1859,1)

;time = ispan (1,1000,1)
;time2 = ispan (1,1000,20)

time = ispan (1,100,1)

;global_serie=new( (/1,60/),float)   ; Method 1
;global_serie(0,0:59)= avg_time_global


;region_serie=new( (/1,60/),float)   ; Method 1
;region_serie(0,0:59)= avg_time_region

printVarSummary(mc4_time_26n)
printVarSummary(mc6_time_26n)

;combined_serie=new( (/2,1859/),float)   ; Method 1
;combined_serie(0,0:999) = mc1_time_26n
;combined_serie(0,1000:1709)= mc4_time_26n
;combined_serie(1,0:999)= mc2_time_26n
;combined_serie(1,1000:1858)= mc6_time_26n



;printvarSummary(mc1_time_26n_cycle)
;print(max(mc1_time_26n_cycle))

;printvarSummary(mc6_time_26n_cycle)
;print(max(mc6_time_26n_cycle))

;combined_serie=new( (/2,100/),float)   ; Method 1
;combined_serie(0,:)= mc1_time_26n_cycle
;combined_serie(1,:)= mc6_time_26n_cycle


;combined_serie2=new( (/2,100/),float)   ; Method 1
;combined_serie2(0,:)= mc6_time_40n_cycle
;combined_serie2(1,:)= mc6_time_40n_cycle

printVarSummary(mc1_time_26n)
printVarSummary(mc6_time_26n)

combined_serie=new( (/2,100/),float)   ; Method 1
combined_serie(0,:)= mc1_time_26n
combined_serie(1,:)= mc6_time_26n

print(avg(mc1_time_26n(50:99)))
print(avg(mc6_time_26n(50:99)))

print(avg(mc1_time_26n_tamas))
print(avg(mc6_time_26n_tamas))

print(avg(mc1_time_20s(50:99)))
print(avg(mc6_time_20s(50:99)))


combined_serie2=new( (/2,100/),float)   ; Method 1
combined_serie2(0,:)= mc1_time_20s
combined_serie2(1,:)= mc6_time_20s

print(avg(mc1_time_20s_tamas))
print(avg(mc6_time_20s_tamas))


;combined_serie3=new( (/2,1859/),float)   ; Method 1
;combined_serie3(0,0:999) = mc1_time_20s
;combined_serie3(0,1000:1709)= mc4_time_20s
;combined_serie3(1,0:999)= mc2_time_20s
;combined_serie3(1,1000:1858)= mc6_time_20s



;combined_serie4=new( (/2,100/),float)   ; Method 1
;combined_serie4(0,:)= mc4_time_32s_cycle
;combined_serie4(1,:)= mc6_time_32s_cycle

combined_serie4=new( (/2,100/),float)   ; Method 1
combined_serie4(0,:)= mc4_time_32s
combined_serie4(1,:)= mc6_time_32s


;print(max(mc6_time_26n_cycle))
;print(min(mc6_time_26n_cycle))


printVarSummary(amoc_pi_26n_mean)
print(amoc_pi_26n_mean)

printVarSummary(amoc_lgm_26n_mean)
print(amoc_lgm_26n_mean)


printVarSummary(amoc_pi_32s_mean)
print(amoc_pi_32s_mean)

printVarSummary(amoc_lgm_32s_mean)
print(amoc_lgm_32s_mean)


latts = pop_pi_new->lat_aux_grid
print(latts(274))
print(latts(89))
print(latts(115))
print(latts(73))
print(latts(291))

depthh = pop_pi_new->moc_z


print(depthh)


depthh2 = depthh(33:60)

pi_val1 = pop_pi_tamas_no->MOC(0,1,0,33:60,:)
print("pi_no_idemix: 26n,32s,20s,40s,35n")
print(max(pi_val1(:,274)))
print(min(pi_val1(:,89)))
print(min(pi_val1(:,115)))
print(min(pi_val1(:,73)))
print(max(pi_val1(:,291)))


lgm_val1 = pop_lgm_wilmes_no->MOC(0,1,0,33:60,:)
print("lgm_wilmes_no_idemix: 26n,32s,20s,40s,35n")
print(max(lgm_val1(:,274)))
print(min(lgm_val1(:,89)))
print(min(lgm_val1(:,115)))
print(min(lgm_val1(:,73)))
print(max(lgm_val1(:,291)))

pi_val2 = pop_pi_ide->MOC(0,1,0,33:60,:)
print("pi_idemix: 26n,32s,20s,40s,35n")
print(max(pi_val2(:,274)))
print(min(pi_val2(:,89)))
print(min(pi_val2(:,115)))
print(min(pi_val2(:,73)))
print(max(pi_val2(:,291)))


lgm_val2 = pop_lgm_wilmes->MOC(0,1,0,33:60,:)
print("lgm_wilmes_idemix: 26n,32s,20s,40s,35n")
print(max(lgm_val2(:,274)))
print(min(lgm_val2(:,89)))
print(min(lgm_val2(:,115)))
print(min(lgm_val2(:,73)))
print(max(lgm_val2(:,291)))


lgm_val3 = pop_lgm_dinezio_no->MOC(0,1,0,33:60,:)
print("lgm_dinezio_no_idemix: 26n,32s,20s,40s,35n")
print(max(lgm_val3(:,274)))
print(min(lgm_val3(:,89)))
print(min(lgm_val3(:,115)))
print(min(lgm_val3(:,73)))
print(max(lgm_val3(:,291)))


lgm_val4 = pop_lgm_dinezio->MOC(0,1,0,33:60,:)
print("lgm_dinezio_idemix: 26n,32s,20s,40s,35n")
print(max(lgm_val4(:,274)))
print(min(lgm_val4(:,89)))
print(min(lgm_val4(:,115)))
print(min(lgm_val4(:,73)))
print(max(lgm_val4(:,291)))


pi_val1_avg = dim_avg_n(pi_val1(0:15,:),0) ; mean for 500 - 2500 m
printVarSummary(pi_val1_avg)

print(max(pi_val1_avg))
pi_val1_avg_max = maxind(pi_val1_avg)          ; index of the maximum value
pi_lat_max1 = latts(pi_val1_avg_max)          ; corresponding latitude

print(pi_val1_avg_max)
print("pi_no_idemix, Maximum occurs at latitude: " + pi_lat_max1)
;print(pi_val1_avg)


pi_val1_avg_coor = pi_val1(:,pi_val1_avg_max)
print(pi_val1_avg_coor)
pi_depth1 = maxind(pi_val1_avg_coor)
;print(depth1)
pi_depth1_max = depthh2(pi_depth1)

print("pi_no_idemix, Maximum occurs at that latitude at depth (m): " + pi_depth1_max*0.01)




pi_val2_avg = dim_avg_n(pi_val2(0:15,:),0) ; mean for 500 - 2500 m
printVarSummary(pi_val2_avg)

print(max(pi_val2_avg))
pi_val2_avg_max = maxind(pi_val2_avg)          ; index of the maximum value
pi_lat_max2 = latts(pi_val2_avg_max)          ; corresponding latitude

print(pi_val2_avg_max)
print("pi_idemix, Maximum occurs at latitude: " + pi_lat_max2)
;print(pi_val2_avg)


pi_val2_avg_coor = pi_val2(:,pi_val2_avg_max)
;print(pi_val2_avg_coor)
pi_depth2 = maxind(pi_val2_avg_coor)
;print(depth2)
pi_depth2_max = depthh2(pi_depth2)

print("pi_idemix, Maximum occurs at that latitude at depth (m): " + pi_depth2_max*0.01)


lgm_val1_avg = dim_avg_n(lgm_val1(0:15,:),0) ; mean for 500 - 2500 m
printVarSummary(lgm_val1_avg)

print(max(lgm_val1_avg))
lgm_val1_avg_max = maxind(lgm_val1_avg)          ; index of the maximum value
lat_max1 = latts(lgm_val1_avg_max)          ; corresponding latitude

print(lgm_val1_avg_max)
print("lgm_wilmes_no_idemix, Maximum occurs at latitude: " + lat_max1)
;print(lgm_val1_avg)


depthh3 = depthh2(0:15)
;print(depthh3)

lgm_val1_avg_coor = lgm_val1(:,lgm_val1_avg_max)
;print(lgm_val1_avg_coor)
depth1 = maxind(lgm_val1_avg_coor)
print(depth1)
depth1_max = depthh2(depth1)

print("lgm_wilmes_no_idemix, Maximum occurs at that latitude at depth (m): " + depth1_max*0.01)




lgm_val2_avg = dim_avg_n(lgm_val2(0:15,:),0) ; mean for 500 - 2500 m
printVarSummary(lgm_val2_avg)

print(max(lgm_val2_avg))
lgm_val2_avg_max = maxind(lgm_val2_avg)          ; index of the maximum value
lat_max2 = latts(lgm_val2_avg_max)          ; corresponding latitude

print(lgm_val2_avg_max)
print("lgm_wilmes_idemix, Maximum occurs at latitude: " + lat_max2)
;print(lgm_val2_avg)


lgm_val2_avg_coor = lgm_val2(:,lgm_val2_avg_max)
;print(lgm_val2_avg_coor)
depth2 = maxind(lgm_val2_avg_coor)
;print(depth2)
depth2_max = depthh3(depth2)

print("lgm_wilmes_idemix, Maximum occurs at that latitude at depth (m): " + depth2_max*0.01)



lgm_val3_avg = dim_avg_n(lgm_val3(0:15,:),0) ; mean for 500 - 2500 m
printVarSummary(lgm_val3_avg)

print(max(lgm_val3_avg))
lgm_val3_avg_max = maxind(lgm_val3_avg)          ; index of the maximum value
lat_max3 = latts(lgm_val3_avg_max)          ; corresponding latitude

print(lgm_val3_avg_max)
print("lgm_dinezio_no_idemix, Maximum occurs at latitude: " + lat_max3)
;print(lgm_val3_avg)



lgm_val3_avg_coor = lgm_val3(:,lgm_val3_avg_max)
;print(lgm_val3_avg_coor)
depth3 = maxind(lgm_val3_avg_coor)
;print(depth3)
depth3_max = depthh2(depth3)

print("lgm_dinezio_no_idemix, Maximum occurs at that latitude at depth (m): " + depth3_max*0.01)


print(max(lgm_val3(:,lgm_val3_avg_max)))
print(max(lgm_val3(depth3,:)))
print(max(lgm_val3))


lgm_val4_avg = dim_avg_n(lgm_val4(0:15,:),0) ; mean for 500 - 2500 m
printVarSummary(lgm_val4_avg)

print(max(lgm_val4_avg))
lgm_val4_avg_max = maxind(lgm_val4_avg)          ; index of the maximum value
lat_max4 = latts(lgm_val4_avg_max)          ; corresponding latitude

print(lgm_val4_avg_max)
print("lgm_dinezio_idemix, Maximum occurs at latitude: " + lat_max4)
;print(lgm_val4_avg)


lgm_val4_avg_coor = lgm_val4(:,lgm_val4_avg_max)
;print(lgm_val4_avg_coor)
depth4 = maxind(lgm_val4_avg_coor)
;print(depth4)
depth4_max = depthh2(depth4)

print("lgm_dinezio_idemix, Maximum occurs at that latitude at depth (m): " + depth4_max*0.01)

resL = True
   resL@tiXAxisFontHeightF   = 0.019
   resL@tiYAxisFontHeightF   = 0.019
   resL@tmXBLabelFontHeightF  = 0.017
   resL@tmYLLabelFontHeightF  = 0.017
  resL@tiYAxisString = "Dissolved Oxygen" +"  "+"(mmol/m3)"  ; axis string
  resL@tiXAxisString = "Time (yrs)"  ; axis string
;  resL@gsnLeftString = "1-2 km"
  resL@xyMonoDashPattern = True
  resL@gsnLeftStringParallelPosF = 0.01

    resL@gsnDraw              = False                ; don't draw yet
   resL@gsnFrame             = False                ; don't advance frame yet
   resL@tiXAxisFontHeightF   = 0.020
   resL@tiYAxisFontHeightF   = 0.020
   resL@tmXBLabelFontHeightF  = 0.017
   resL@tmYLLabelFontHeightF  = 0.017
  resL@tiYAxisString = " NADW" +"  "+"[Sv]"  ; axis string
  resL@tiXAxisString = "Simulation year"  ; axis string
;  resL@gsnLeftString = "1-2 km"
  resL@xyMonoDashPattern = True
  resL@gsnLeftStringParallelPosF = 0.03
   resL@tiXAxisOffsetYF = -0.019
   resL@tiYAxisOffsetXF = -0.013
  

   labels = ispan(0,999,20)
   print(labels)
 ;   resL@tmXUseBottom  = False      ; Keep top axis independent of bottom.
 ;  resL@tmXTLabelsOn  = True       ; have tick mark labels
 ;  resL@tmXTOn        = True       ; have tick marks
 ;  resL@tmXBMode  = "Explicit"
 ;  resL@tmXBValues    = labels 
 ; resL@tmXBLabels = labels

  resL@txFontHeightF      = 0.015
  resL@pmLegendDisplayMode    = "Always"            ; turn on legend
  resL@lgPerimOn              = False               ; turn off box around
  resL@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL@xyLineThicknesses      = (/1.5,1.5/)     ; line thickness
  resL@xyLineColors           = (/"red","blue"/)
 resL@lgPerimOn              = False                 ; no box around

;;

  resL@lgJustification        = "TopRight"; 
  resL@pmLegendOrthogonalPosF = -1.42     ; Move into plot
  resL@pmLegendParallelPosF   = 1.1       ; Move to right
  resL@pmLegendWidthF         = 0.3       ; Decrease width
  resL@pmLegendHeightF        = 0.12       ; Decrease height
  resL@lgBoxMinorExtentF      = 0.10       ; Shorten the legend lines
  resL@lgLabelFontHeightF     = 0.017     ; Change the font size
  resL@lgLabelFontThicknessF = 1.2
  resL@lgLineLabelFontThicknessF = 1.2

  resL@tmYLMode  = "Manual"
  resL@trYMinF = 14
  resL@trYMaxF = 22
  
;  resL@tmXBMode   = "Explicit"			; explicitly set the Y-axis labels
;  resL@tmXBValues = (/0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960/)
 ;  resL@tmXBTickSpacingF = 200
;  resL@tmYLLabels = "" + resL@tmYLValues 
  resL@tmXBMode  = "Manual"

  resL@trXMinF = 0
  resL@trXMaxF = 100
  resL@tmXBTickSpacingF = 5

;!  resL@trXMinF = 0
;!  resL@trXMaxF = 980
;!resL@tmXBTickSpacingF = 200


; labels for the legend
; resL@xyExplicitLegendLabels = (/"Global"/)


resL@lgLineThicknesses = 1.3 

resL@xyExplicitLegendLabels = (/"PI - w/IDEMIX", "LGM - w/IDEMIX"/)

resL@gsnXRefLine = 1000

resL@gsnXRefLineDashPattern = 1

resL@gsnXRefLineColor = "black"

resL@gsnXRefLineThicknessF = 2.0

;plot = gsn_csm_xy(wks,time,combined_serie,resL)

;print(max(amoc1))
;print(min(amoc1))

;print(max(amoc4))
;print(min(amoc4))


;print(max(amoc2))
;print(min(amoc2))


;print(max(amoc6))
;print(min(amoc6))



txres1 = True
txres1@txFontHeightF = 0.017

;gsn_text_ndc(wks,"AMOC",0.26,0.86,txres1)



txres4 = True
txres4@txFontHeightF = 0.017

;gsn_text_ndc(wks,"[Sv]",0.79,0.86,txres4)

 

txres5 = True
txres5@txFontHeightF = 0.013
txres5@txAngleF      = 0.     ; rotate 0 degrees
gsn_text_ndc(wks,"@26.5N",0.27,0.98,txres5)
gsn_text_ndc(wks,"@32S",0.76,0.98,txres5)
;gsn_text_ndc(wks,"w/o IDEMIX",0.19,0.6,txres5)
;gsn_text_ndc(wks,"w/ IDEMIX",0.32,0.60,txres5)


;
; Add another text string.
;
;  txres3@txBackgroundFillColor = "Orange"
;  txres3@txFont                = "times-bold"

;  txres4 = True
;  txid2 = gsn_create_text(wks, "Depth (m)", txres4)

; amres = True

;  amres@amParallelPosF   = -0.5    ; This is the left edge of the plot.
;  amres@amOrthogonalPosF =  0.6    ; This is the bottom edge of the plot.
;  amres@amJust           = "Left"
;  annoid2 = gsn_add_annotation(plot(0), txid2, amres)


;  txres5 = True
;  txid3 = gsn_create_text(wks, "Depth (m)", txres5)

;  amres2 = True

;  amres2@amParallelPosF   = -0.5    ; This is the left edge of the plot.
;  amres2@amOrthogonalPosF =  0.4    ; This is the bottom edge of the plot.
;  amres2@amJust           = "Left"
;  annoid2 = gsn_add_annotation(plots(0), txid3, amres2)

;  txres6 = True

;  txres6@txAngleF      = 90.     ; rotate 90 degrees
;  txres6@txJust        = "TopLeft"
;   txres6@txHeightF = 0.011
;  gsn_text_ndc(wks,"Depth (m)",0.17,0.33,txres6)

;  txres7 = True

;   txres7@txAngleF      = 90.     ; rotate 90 degrees
;  txres7@txJust        = "TopLeft"
;   txres7@txFontHeightF = 0.011
;  gsn_text_ndc(wks,"Depth (m)",0.17,0.62,txres7)
  
;   resL@trYMinF = -8
;   resL@trYMaxF =  2
   resL@tiXAxisOffsetYF = -0.019
   resL@tiYAxisOffsetXF = -0.021
 resL@tiYAxisString = " NADW" +"  "+"[Sv]"  ; axis string
resL@xyDashPatterns         = (/1,1/)


;  frame(wks)

resL2 = resL
delete(resL2@trYMinF)
delete(resL2@trYMaxF)
resL2@tmYLMode  = "Manual"
;resL2@tmYLValues = (/12,14,16,18,20,22,24/)
   resL2@trYMinF = -4.6
   resL2@trYMaxF = -1.6   
   resL2@tiXAxisOffsetYF = -0.019
   resL2@tiYAxisOffsetXF = -0.017
 ;resL2@tiYAxisString = " NADW" +"  "+"[Sv]"  ; axis string
  resL2@tiYAxisString = " AABW" +"  "+"[Sv]"  ; axis string
  resL2@xyDashPatterns         = (/1,1/)



plot(0) = gsn_csm_xy(wks,time,combined_serie,resL)
plot(1) = gsn_csm_xy(wks,time,combined_serie2,resL2)
plot(2) = gsn_csm_xy(wks,time,combined_serie4,resL2)
plot(3) = gsn_csm_xy(wks,time,combined_serie4,resL2)



resP                 = True
resP@gsnMaximize        = True
resP@gsnFrame         = False               ; save panel until both ready
resP@gsnPaperOrientation = "landscape"
resP@gsnPanelLabelBar   = False

  gsn_panel(wks,plot,(/2,2/),resP)

end
