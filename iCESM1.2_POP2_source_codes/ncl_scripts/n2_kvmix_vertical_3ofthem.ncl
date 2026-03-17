

begin


   pop_pi_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
;   pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")
    pop_lgm_tamas = addfile ("/mnt/d/lgm_data2/select_test_cmake_1670_1719_timmean.nc", "r")

;   pop_pi_tamas_n2 = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
;   pop_lgm_tamas_n2 = addfile ("/mnt/d/lgm_data2/montly_parameters.LGM21ka_frUsp.p875.001.pop.h.5450_5499_timmean.nc", "r")


   pop_pi_n2_idemix = addfile ("/mnt/d/lgm_data/work4_work5_n2_pi_idemix_isotope_continued_5000_5099_timmean.nc", "r")
   pop_lgm_n2_idemix = addfile ("/mnt/d/lgm_data/work4_lgm_5849_5899_timmean.nc", "r")

    
    pop_pi = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4654_4703_timmean.nc", "r")
    pop_lgm = addfile ("/mnt/d/lgm_data/select_lgm_idemix_isotope_continued_5849_5899_timmean.nc", "r")

   
   pop_pi_tamas_n2 = addfile ("/mnt/d/lgm_data2/monthly_pi_tamas_4150_4199_short_timmean.nc", "r")
   ;pop_lgm_tamas_n2 = addfile ("/mnt/d/lgm_data2/monthly_lgm_tamas_5450_5499_short_timmean.nc", "r")
    pop_lgm_tamas_n2 = addfile ("/mnt/d/lgm_data2/select_test_cmake_1670_1719_timmean.nc", "r")

    pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")




  pop_pi_tamas_kvmix = addfile ("/mnt/d/lgm_data/kvmix_pi_tamas_4190_4199_timmean.nc", "r")                                                                                                                                                                                                                   
pop_lgm_tamas_kvmixd = addfile ("/mnt/d/lgm_data/kvmix_lgm_tamas_5490_5499_timmean.nc", "r")                                                              
pop_lgm_tamas_kvmix = addfile ("/mnt/d/lgm_data2/kvmix_test_cmake_1790_1799_timmean.nc", "r")

pop_pi_kvmix = addfile ("/mnt/d/lgm_data/select_pi_idemix_4900_4694_4703_timmean.nc", "r")

pop_lgm_kvmix = addfile ("/mnt/d/lgm_data/kvmix_lgm_idemix_isotope_continued_5900_5909_timmean.nc", "r")

pop_lgm_kvmixd = addfile ("/mnt/d/lgm_data/kvmix_dinezio_continued.pop.h.1850_1859_timmean.nc", "r")


   vol_pi = addfile("/mnt/d/lgm_data/grid_pi_volume.nc","r")


   vol_lgm = addfile("/mnt/d/lgm_data/grid_lgm_volume.nc","r")

  pi_mask = addfile ("/mnt/d/lgm_data/mask_salt.nc", "r")

  lgm_mask = addfile ("/mnt/d/lgm_data/mask_salt_lgm.nc","r")



   wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"n2_kvmix_vertical_profile_lastnew_3ofthem")



;   cmap = read_colormap_file("cmocean_delta")
    cmap = read_colormap_file("cmocean_haline")
   ncmap = dimsizes(cmap(:,0))

resL3 = True
pltres = True
mpres = True
pltres@PanelPlot  = True
resL3@Footer    = False
plot = new(6,graphic)


pi_rho_tamas         = pop_pi_tamas_n2->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

z    = pop_pi_tamas->z_t

;--- Dimensions
  nlev = dimsizes(z)
  nlat = dimsizes(pi_rho_tamas(1,:,0))
  nlon = dimsizes(pi_rho_tamas(1,0,:))


depth           = pop_pi->z_t(0:58)
depth_m2          = depth*1.E-2       ; convert from cm to m

 ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

  ;pi_rho_depth_tamas = new((/nlev-1, nlat, nlon/), typeof(pi_rho_tamas), "float")

   pi_rho_depth_tamas = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z(k+1) - z(k)
    pi_rho_depth_tamas(k,:,:) = (pi_rho_tamas(k+1,:,:) - pi_rho_tamas(k,:,:)) / dz
  end do


water_ref_density = 1025 ;kg/m3


pi_n2_tamas   = -1*(((9.80665*pi_rho_depth_tamas)/water_ref_density)/0.000001)


iw1        = pop_pi_tamas_kvmix->KVMIX(0,:,:,:)

mean_temp1 = dim_avg_n_Wrap(iw1, (/1,2/))

mean_temp1_m2 = mean_temp1 * 0.0001

mean_temp1_log = log10(mean_temp1_m2)





lgm_rho_tamas         = pop_lgm_tamas_n2->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

;water_ref_density = 1000 kg/m3

 ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

  ;lgm_rho_depth_tamas = new((/nlev-1, nlat, nlon/), typeof(lgm_rho_tamas), "float")

   lgm_rho_depth_tamas = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z(k+1) - z(k)
    lgm_rho_depth_tamas(k,:,:) = (lgm_rho_tamas(k+1,:,:) - lgm_rho_tamas(k,:,:)) / dz
  end do


;water_ref_density = 1025 kg/m3

lgm_n2_tamas    = -1*(((9.80665*lgm_rho_depth_tamas)/water_ref_density)/0.000001)




iw2        = pop_lgm_tamas_kvmix->KVMIX(0,:,:,:)

mean_temp2 = dim_avg_n_Wrap(iw2, (/1,2/))

mean_temp2_m2 = mean_temp2 * 0.0001

mean_temp2_log = log10(mean_temp2_m2)


pi_rho         = pop_pi->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

;water_ref_density = 1000 kg/m3



   pi_rho_depth = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z(k+1) - z(k)
    pi_rho_depth(k,:,:) = (pi_rho(k+1,:,:) - pi_rho(k,:,:)) / dz
  end do


;water_ref_density = 1025 kg/m3

pi_n2    = -1*(((9.80665*pi_rho_depth)/water_ref_density)/0.000001)



iw3        = pop_pi_kvmix->KVMIX(0,:,:,:)

mean_temp3 = dim_avg_n_Wrap(iw3, (/1,2/))

mean_temp3_m2 = mean_temp3 * 0.0001

mean_temp3_log = log10(mean_temp3_m2)



lgm_rho         = pop_lgm->RHO(0,:,:,:)*1000 ; from g/cm3 to kg/m3

  ;--- Create output array for drho/dz, which is lgm_rho_tamas_depth here

  ;lgm_rho_depth = new((/nlev-1, nlat, nlon/), typeof(lgm_rho), "float")

   lgm_rho_depth = new((/nlev-1, nlat, nlon/), "float")

  ;--- Compute drho/dz
  do k = 0, nlev-2
    dz = z(k+1) - z(k)
    lgm_rho_depth(k,:,:) = (lgm_rho(k+1,:,:) - lgm_rho(k,:,:)) / dz
  end do


;water_ref_density = 1025 kg/m3

lgm_n2    = -1*(((9.80665*lgm_rho_depth)/water_ref_density)/0.000001)



iw4        = pop_lgm_kvmix->KVMIX(0,:,:,:)

mean_temp4 = dim_avg_n_Wrap(iw4, (/1,2/))

;;mean_temp4_m2 = mean_temp4 * 0.0001

;mean_temp4_log = log10(mean_temp4_m2)



region_pi  = pop_pi_region->REGION_MASK
region_lgm  = pop_lgm_region->REGION_MASK
region_atlantic_pi= new((/384,320/),float)
region_atlantic_lgm= new((/384,320/),float) 
region_po_pi= new((/384,320/),float)
region_po_lgm= new((/384,320/),float)
region_so_pi= new((/384,320/),float)
region_so_lgm= new((/384,320/),float)


;### ATLANTIC OCEAN;###
region_atlantic_pi        = where(region_pi.eq.6,1,9.969209968386869e+36)
region_atlantic_lgm        = where(region_lgm.eq.6,1,9.969209968386869e+36)

;### SO OCEAN;###
region_so_pi        = where(region_pi.ne.1,9.969209968386869e+36,1)
region_so_lgm        = where(region_lgm.ne.1,9.969209968386869e+36,1)


;### PO OCEAN;###
region_po_pi        = where(region_pi.ne.2,9.969209968386869e+36,1)
region_po_lgm        = where(region_lgm.ne.2,9.969209968386869e+36,1)




vol_pi_ocean = vol_pi->dz(:,:,:)


vol_lgm_ocean = vol_lgm->dz(:,:,:)


mask_pi = pi_mask->mask

mask_lgm = lgm_mask->mask


volume_pi = new((/60,384,320/),float)

volume_pi = vol_pi_ocean(:,:,:) * mask_pi(0,:,:,:) * (10^-6) ; ocean mask and unit conversion from cm3 to m3


volume_lgm = new((/60,384,320/),float)

volume_lgm = vol_lgm_ocean(:,:,:) * mask_pi(0,:,:,:) * (10^-6) ; ocean mask and unit conversion from cm3 to m3



vol_pi_v2 = dim_sum_n(volume_pi,0)

print(vol_pi_v2)


vol_lgm_v2 = dim_sum_n(volume_lgm,0)


iw1_vol = iw1*volume_pi
iw2_vol = iw2*volume_lgm
iw3_vol = iw3*volume_pi
iw4_vol = iw4*volume_lgm

vol_pi_v2_na = mask(volume_pi,region_atlantic_pi,1)
vol_lgm_v2_na = mask(volume_lgm,region_atlantic_lgm,1)
vol_pi_v2_po = mask(volume_pi,region_po_pi,1)
vol_lgm_v2_po = mask(volume_lgm,region_po_lgm,1)

vol_pi_v2_so = mask(volume_pi,region_so_pi,1)
vol_lgm_v2_so = mask(volume_lgm,region_so_lgm,1)


pi_n2_tamas_na = mask(pi_n2_tamas,region_atlantic_pi,1)
lgm_n2_tamas_na = mask(lgm_n2_tamas,region_atlantic_lgm,1)
pi_n2_na = mask(pi_n2,region_atlantic_pi,1)
lgm_n2_na = mask(lgm_n2,region_atlantic_lgm,1)

pi_n2_tamas_mean_na    = dim_avg_n_Wrap(pi_n2_tamas_na,(/1,2/))
lgm_n2_tamas_mean_na    = dim_avg_n_Wrap(lgm_n2_tamas_na,(/1,2/))
pi_n2_mean_na    = dim_avg_n_Wrap(pi_n2_na,(/1,2/))
lgm_n2_mean_na    = dim_avg_n_Wrap(lgm_n2_na,(/1,2/))



printVarSummary(lgm_n2_mean_na)
print(max(lgm_n2_mean_na))
print(min(lgm_n2_mean_na))

pi_n2_tamas_mean_gb    = dim_avg_n_Wrap(pi_n2_tamas,(/1,2/))
lgm_n2_tamas_mean_gb    = dim_avg_n_Wrap(lgm_n2_tamas,(/1,2/))
pi_n2_mean_gb    = dim_avg_n_Wrap(pi_n2,(/1,2/))
lgm_n2_mean_gb    = dim_avg_n_Wrap(lgm_n2,(/1,2/))




pi_n2_tamas_so = mask(pi_n2_tamas,region_so_pi,1)
lgm_n2_tamas_so = mask(lgm_n2_tamas,region_so_lgm,1)
pi_n2_so = mask(pi_n2,region_so_pi,1)
lgm_n2_so = mask(lgm_n2,region_so_lgm,1)

pi_n2_tamas_mean_so    = dim_avg_n_Wrap(pi_n2_tamas_so,(/1,2/))
lgm_n2_tamas_mean_so    = dim_avg_n_Wrap(lgm_n2_tamas_so,(/1,2/))
pi_n2_mean_so    = dim_avg_n_Wrap(pi_n2_so,(/1,2/))
lgm_n2_mean_so    = dim_avg_n_Wrap(lgm_n2_so,(/1,2/))




pi_kvmix_tamas_na = mask(iw1,region_atlantic_pi,1)
lgm_kvmix_tamas_na = mask(iw2,region_atlantic_lgm,1)
pi_kvmix_na = mask(iw3,region_atlantic_pi,1)
lgm_kvmix_na = mask(iw4,region_atlantic_lgm,1)

pi_kvmix_tamas_mean_na    = dim_avg_n_Wrap(pi_kvmix_tamas_na,(/1,2/))
lgm_kvmix_tamas_mean_na    = dim_avg_n_Wrap(lgm_kvmix_tamas_na,(/1,2/))
pi_kvmix_mean_na    = dim_avg_n_Wrap(pi_kvmix_na,(/1,2/))
lgm_kvmix_mean_na    = dim_avg_n_Wrap(lgm_kvmix_na,(/1,2/))

pi_kvmix_tamas_po = mask(iw1,region_po_pi,1)
lgm_kvmix_tamas_po = mask(iw2,region_po_lgm,1)
pi_kvmix_po = mask(iw3,region_po_pi,1)
lgm_kvmix_po = mask(iw4,region_po_lgm,1)


pi_kvmix_tamas_mean_po    = dim_avg_n_Wrap(pi_kvmix_tamas_po,(/1,2/))
lgm_kvmix_tamas_mean_po    = dim_avg_n_Wrap(lgm_kvmix_tamas_po,(/1,2/))
pi_kvmix_mean_po    = dim_avg_n_Wrap(pi_kvmix_po,(/1,2/))
lgm_kvmix_mean_po    = dim_avg_n_Wrap(lgm_kvmix_po,(/1,2/))


pi_kvmix_tamas_mean_gb    = dim_avg_n_Wrap(iw1,(/1,2/))
lgm_kvmix_tamas_mean_gb    = dim_avg_n_Wrap(iw2,(/1,2/))
pi_kvmix_mean_gb    = dim_avg_n_Wrap(iw3,(/1,2/))
lgm_kvmix_mean_gb    = dim_avg_n_Wrap(iw4,(/1,2/))



pi_kvmix_tamas_so = mask(iw1,region_so_pi,1)
lgm_kvmix_tamas_so = mask(iw2,region_so_lgm,1)
pi_kvmix_so = mask(iw3,region_so_pi,1)
lgm_kvmix_so = mask(iw4,region_so_lgm,1)

pi_kvmix_tamas_mean_so    = dim_avg_n_Wrap(pi_kvmix_tamas_so,(/1,2/))
lgm_kvmix_tamas_mean_so    = dim_avg_n_Wrap(lgm_kvmix_tamas_so,(/1,2/))
pi_kvmix_mean_so    = dim_avg_n_Wrap(pi_kvmix_so,(/1,2/))
lgm_kvmix_mean_so    = dim_avg_n_Wrap(lgm_kvmix_so,(/1,2/))






pi_kvmix_tamas_mean_gb_m2 = pi_kvmix_tamas_mean_gb * 0.0001

pi_kvmix_tamas_mean_gb_m2_log = log10(pi_kvmix_tamas_mean_gb_m2)


lgm_kvmix_tamas_mean_gb_m2 = lgm_kvmix_tamas_mean_gb * 0.0001

lgm_kvmix_tamas_mean_gb_m2_log = log10(lgm_kvmix_tamas_mean_gb_m2)



pi_kvmix_mean_gb_m2 = pi_kvmix_mean_gb * 0.0001

pi_kvmix_mean_gb_m2_log = log10(pi_kvmix_mean_gb_m2)


lgm_kvmix_mean_gb_m2 = lgm_kvmix_mean_gb * 0.0001

lgm_kvmix_mean_gb_m2_log = log10(lgm_kvmix_mean_gb_m2)




pi_kvmix_tamas_mean_na_m2 = pi_kvmix_tamas_mean_na * 0.0001

pi_kvmix_tamas_mean_na_m2_log = log10(pi_kvmix_tamas_mean_na_m2)


lgm_kvmix_tamas_mean_na_m2 = lgm_kvmix_tamas_mean_na * 0.0001

lgm_kvmix_tamas_mean_na_m2_log = log10(lgm_kvmix_tamas_mean_na_m2)


pi_kvmix_mean_na_m2 = pi_kvmix_mean_na * 0.0001

pi_kvmix_mean_na_m2_log = log10(pi_kvmix_mean_na_m2)


lgm_kvmix_mean_na_m2 = lgm_kvmix_mean_na * 0.0001

lgm_kvmix_mean_na_m2_log = log10(lgm_kvmix_mean_na_m2)




pi_kvmix_tamas_mean_so_m2 = pi_kvmix_tamas_mean_so * 0.0001

pi_kvmix_tamas_mean_so_m2_log = log10(pi_kvmix_tamas_mean_so_m2)


lgm_kvmix_tamas_mean_so_m2 = lgm_kvmix_tamas_mean_so * 0.0001

lgm_kvmix_tamas_mean_so_m2_log = log10(lgm_kvmix_tamas_mean_so_m2)


pi_kvmix_mean_so_m2 = pi_kvmix_mean_so * 0.0001

pi_kvmix_mean_so_m2_log = log10(pi_kvmix_mean_so_m2)


lgm_kvmix_mean_so_m2 = lgm_kvmix_mean_so * 0.0001

lgm_kvmix_mean_so_m2_log = log10(lgm_kvmix_mean_so_m2)


vol_pi_mean_gb    = dim_avg_n_Wrap(volume_pi,(/1,2/))
vol_lgm_mean_gb    = dim_avg_n_Wrap(volume_lgm,(/1,2/))


vol_pi_mean_na    = dim_avg_n_Wrap(vol_pi_v2_na,(/1,2/))
vol_lgm_mean_na    = dim_avg_n_Wrap(vol_lgm_v2_na,(/1,2/))


vol_pi_mean_po    = dim_avg_n_Wrap(vol_pi_v2_po,(/1,2/))
vol_lgm_mean_po    = dim_avg_n_Wrap(vol_lgm_v2_po,(/1,2/))

vol_pi_mean_so    = dim_avg_n_Wrap(vol_pi_v2_so,(/1,2/))
vol_lgm_mean_so    = dim_avg_n_Wrap(vol_lgm_v2_so,(/1,2/))


d1_na1 = (avg(pi_kvmix_tamas_mean_na(19:40))*0.0001)
d1_na2 = (avg(pi_kvmix_tamas_mean_na(40:46))*0.0001)
d1_na3 = (avg(pi_kvmix_tamas_mean_na(46:54))*0.0001)
d1_na4 = (avg(pi_kvmix_tamas_mean_na)*0.0001)

d1_po1 = (avg(pi_kvmix_tamas_mean_po(19:40))*0.0001)
d1_po2 = (avg(pi_kvmix_tamas_mean_po(40:46))*0.0001)
d1_po3 = (avg(pi_kvmix_tamas_mean_po(46:54))*0.0001)
d1_po4 = (avg(pi_kvmix_tamas_mean_po)*0.0001)

d1_so1 = (avg(pi_kvmix_tamas_mean_so(19:40))*0.0001)
d1_so2 = (avg(pi_kvmix_tamas_mean_so(40:46))*0.0001)
d1_so3 = (avg(pi_kvmix_tamas_mean_so(46:54))*0.0001)
d1_so4 = (avg(pi_kvmix_tamas_mean_so)*0.0001)

d1_gb1 = (avg(pi_kvmix_tamas_mean_gb(19:40))*0.0001)
d1_gb2 = (avg(pi_kvmix_tamas_mean_gb(40:46))*0.0001)
d1_gb3 = (avg(pi_kvmix_tamas_mean_gb(46:54))*0.0001)
d1_gb4 = (avg(pi_kvmix_tamas_mean_gb)*0.0001)

printVarSummary(d1_na1)
print(log10(d1_na1))
print(log10(d1_na2))
print(log10(d1_na3))
print(log10(d1_po1))
print(log10(d1_po2))
print(log10(d1_po3))
print(log10(d1_so1))
print(log10(d1_so2))
print(log10(d1_so3))
print(log10(d1_gb1))
print(log10(d1_gb2))
print(log10(d1_gb3))


d2_na1 = (avg(lgm_kvmix_tamas_mean_na(19:40))*0.0001)
d2_na2 = (avg(lgm_kvmix_tamas_mean_na(40:46))*0.0001)
d2_na3 = (avg(lgm_kvmix_tamas_mean_na(46:54))*0.0001)
d2_na4 = (avg(lgm_kvmix_tamas_mean_na)*0.0001)

d2_po1 = (avg(lgm_kvmix_tamas_mean_po(19:40))*0.0001)
d2_po2 = (avg(lgm_kvmix_tamas_mean_po(40:46))*0.0001)
d2_po3 = (avg(lgm_kvmix_tamas_mean_po(46:54))*0.0001)
d2_po4 = (avg(lgm_kvmix_tamas_mean_po)*0.0001)

d2_so1 = (avg(lgm_kvmix_tamas_mean_so(19:40))*0.0001)
d2_so2 = (avg(lgm_kvmix_tamas_mean_so(40:46))*0.0001)
d2_so3 = (avg(lgm_kvmix_tamas_mean_so(46:54))*0.0001)
d2_so4 = (avg(lgm_kvmix_tamas_mean_so)*0.0001)

d2_gb1 = (avg(lgm_kvmix_tamas_mean_gb(19:40))*0.0001)
d2_gb2 = (avg(lgm_kvmix_tamas_mean_gb(40:46))*0.0001)
d2_gb3 = (avg(lgm_kvmix_tamas_mean_gb(46:54))*0.0001)
d2_gb4 = (avg(lgm_kvmix_tamas_mean_gb)*0.0001)

printVarSummary(d2_na1)
print(log10(d2_na1))
print(log10(d2_na2))
print(log10(d2_na3))
print(log10(d2_po1))
print(log10(d2_po2))
print(log10(d2_po3))
print(log10(d2_so1))
print(log10(d2_so2))
print(log10(d2_so3))
print(log10(d2_gb1))
print(log10(d2_gb2))
print(log10(d2_gb3))


d3_na1 = (avg(pi_kvmix_mean_na(19:40))*0.0001)
d3_na2 = (avg(pi_kvmix_mean_na(40:46))*0.0001)
d3_na3 = (avg(pi_kvmix_mean_na(46:54))*0.0001)
d3_na4 = (avg(pi_kvmix_mean_na)*0.0001)

d3_po1 = (avg(pi_kvmix_mean_po(19:40))*0.0001)
d3_po2 = (avg(pi_kvmix_mean_po(40:46))*0.0001)
d3_po3 = (avg(pi_kvmix_mean_po(46:54))*0.0001)
d3_po4 = (avg(pi_kvmix_mean_po)*0.0001)

d3_so1 = (avg(pi_kvmix_mean_so(19:40))*0.0001)
d3_so2 = (avg(pi_kvmix_mean_so(40:46))*0.0001)
d3_so3 = (avg(pi_kvmix_mean_so(46:54))*0.0001)
d3_so4 = (avg(pi_kvmix_mean_so)*0.0001)

d3_gb1 = (avg(pi_kvmix_mean_gb(19:40))*0.0001)
d3_gb2 = (avg(pi_kvmix_mean_gb(40:46))*0.0001)
d3_gb3 = (avg(pi_kvmix_mean_gb(46:54))*0.0001)
d3_gb4 = (avg(pi_kvmix_mean_gb)*0.0001)

printVarSummary(d3_na1)
print(log10(d3_na1))
print(log10(d3_na2))
print(log10(d3_na3))
print(log10(d3_po1))
print(log10(d3_po2))
print(log10(d3_po3))
print(log10(d3_so1))
print(log10(d3_so2))
print(log10(d3_so3))
print(log10(d3_gb1))
print(log10(d3_gb2))
print(log10(d3_gb3))

d4_na1 = (avg(lgm_kvmix_mean_na(19:40))*0.0001)
d4_na2 = (avg(lgm_kvmix_mean_na(40:46))*0.0001)
d4_na3 = (avg(lgm_kvmix_mean_na(46:54))*0.0001)
d4_na4 = (avg(lgm_kvmix_mean_na)*0.0001)

d4_po1 = (avg(lgm_kvmix_mean_po(19:40))*0.0001)
d4_po2 = (avg(lgm_kvmix_mean_po(40:46))*0.0001)
d4_po3 = (avg(lgm_kvmix_mean_po(46:54))*0.0001)
d4_po4 = (avg(lgm_kvmix_mean_po)*0.0001)

d4_so1 = (avg(lgm_kvmix_mean_so(19:40))*0.0001)
d4_so2 = (avg(lgm_kvmix_mean_so(40:46))*0.0001)
d4_so3 = (avg(lgm_kvmix_mean_so(46:54))*0.0001)
d4_so4 = (avg(lgm_kvmix_mean_so)*0.0001)

d4_gb1 = (avg(lgm_kvmix_mean_gb(19:40))*0.0001)
d4_gb2 = (avg(lgm_kvmix_mean_gb(40:46))*0.0001)
d4_gb3 = (avg(lgm_kvmix_mean_gb(46:54))*0.0001)
d4_gb4 = (avg(lgm_kvmix_mean_gb)*0.0001)

printVarSummary(d4_na1)
print(log10(d4_na1))
print(log10(d4_na2))
print(log10(d4_na3))
print(log10(d4_po1))
print(log10(d4_po2))
print(log10(d4_po3))
print(log10(d4_so1))
print(log10(d4_so2))
print(log10(d4_so3))
print(log10(d4_gb1))
print(log10(d4_gb2))
print(log10(d4_gb3))


iw5        = pop_lgm_tamas_kvmixd->KVMIX(0,:,:,:)
iw6        = pop_lgm_kvmixd->KVMIX(0,:,:,:)

lgm_kvmix_tamasd_na = mask(iw5,region_atlantic_lgm,1)

lgm_kvmixd_na = mask(iw6,region_atlantic_lgm,1)

lgm_kvmix_tamasd_mean_na    = dim_avg_n_Wrap(lgm_kvmix_tamasd_na,(/1,2/))

lgm_kvmixd_mean_na    = dim_avg_n_Wrap(lgm_kvmixd_na,(/1,2/))


lgm_kvmix_tamasd_po = mask(iw5,region_po_lgm,1)

lgm_kvmixd_po = mask(iw6,region_po_lgm,1)

lgm_kvmix_tamasd_mean_po  = dim_avg_n_Wrap(lgm_kvmix_tamasd_po,(/1,2/))

lgm_kvmixd_mean_po    = dim_avg_n_Wrap(lgm_kvmixd_po,(/1,2/))


lgm_kvmix_tamasd_so = mask(iw5,region_so_lgm,1)

lgm_kvmixd_so = mask(iw6,region_so_lgm,1)

lgm_kvmix_tamasd_mean_so = dim_avg_n_Wrap(lgm_kvmix_tamasd_so,(/1,2/))

lgm_kvmixd_mean_so = dim_avg_n_Wrap(lgm_kvmixd_so,(/1,2/))


lgm_kvmix_tamasd_mean_gb    = dim_avg_n_Wrap(iw5,(/1,2/))
lgm_kvmixd_mean_gb    = dim_avg_n_Wrap(iw6,(/1,2/))




d5_na1 = (avg(lgm_kvmix_tamasd_mean_na(19:40))*0.0001)
d5_na2 = (avg(lgm_kvmix_tamasd_mean_na(40:46))*0.0001)
d5_na3 = (avg(lgm_kvmix_tamasd_mean_na(46:54))*0.0001)
d5_na4 = (avg(lgm_kvmix_tamasd_mean_na)*0.0001)

d5_po1 = (avg(lgm_kvmix_tamasd_mean_po(19:40))*0.0001)
d5_po2 = (avg(lgm_kvmix_tamasd_mean_po(40:46))*0.0001)
d5_po3 = (avg(lgm_kvmix_tamasd_mean_po(46:54))*0.0001)
d5_po4 = (avg(lgm_kvmix_tamasd_mean_po)*0.0001)

d5_so1 = (avg(lgm_kvmix_tamasd_mean_so(19:40))*0.0001)
d5_so2 = (avg(lgm_kvmix_tamasd_mean_so(40:46))*0.0001)
d5_so3 = (avg(lgm_kvmix_tamasd_mean_so(46:54))*0.0001)
d5_so4 = (avg(lgm_kvmix_tamasd_mean_so)*0.0001)


d5_gb1 = (avg(lgm_kvmix_tamasd_mean_gb(19:40))*0.0001)
d5_gb2 = (avg(lgm_kvmix_tamasd_mean_gb(40:46))*0.0001)
d5_gb3 = (avg(lgm_kvmix_tamasd_mean_gb(46:54))*0.0001)
d5_gb4 = (avg(lgm_kvmix_tamasd_mean_gb)*0.0001)



d6_na1 = (avg(lgm_kvmixd_mean_na(19:40))*0.0001)
d6_na2 = (avg(lgm_kvmixd_mean_na(40:46))*0.0001)
d6_na3 = (avg(lgm_kvmixd_mean_na(46:54))*0.0001)
d6_na4 = (avg(lgm_kvmixd_mean_na)*0.0001)

d6_po1 = (avg(lgm_kvmixd_mean_po(19:40))*0.0001)
d6_po2 = (avg(lgm_kvmixd_mean_po(40:46))*0.0001)
d6_po3 = (avg(lgm_kvmixd_mean_po(46:54))*0.0001)
d6_po4 = (avg(lgm_kvmixd_mean_po)*0.0001)

d6_so1 = (avg(lgm_kvmixd_mean_so(19:40))*0.0001)
d6_so2 = (avg(lgm_kvmixd_mean_so(40:46))*0.0001)
d6_so3 = (avg(lgm_kvmixd_mean_so(46:54))*0.0001)
d6_so4 = (avg(lgm_kvmixd_mean_so)*0.0001)

d6_gb1 = (avg(lgm_kvmixd_mean_gb(19:40))*0.0001)
d6_gb2 = (avg(lgm_kvmixd_mean_gb(40:46))*0.0001)
d6_gb3 = (avg(lgm_kvmixd_mean_gb(46:54))*0.0001)
d6_gb4 = (avg(lgm_kvmixd_mean_gb)*0.0001)




;p1 = ((abs(d3_na1) - abs(d1_na1))*100)/(abs(d1_na1))
;p2 = ((abs(d3_na2) - abs(d1_na2))*100)/(abs(d1_na2))
;p3 = ((abs(d3_na3) - abs(d1_na3))*100)/(abs(d1_na3))

;p4 = ((abs(d3_po1) - abs(d1_po1))*100)/(abs(d1_po1))
;p5 = ((abs(d3_po2) - abs(d1_po2))*100)/(abs(d1_po2))
;p6 = ((abs(d3_po3) - abs(d1_po3))*100)/(abs(d1_po3))

;p7 = ((abs(d3_so1) - abs(d1_so1))*100)/(abs(d1_so1))
;p8 = ((abs(d3_so2) - abs(d1_so2))*100)/(abs(d1_so2))
;p9 = ((abs(d3_so3) - abs(d1_so3))*100)/(abs(d1_so3))

;p10 = ((abs(d3_gb1) - abs(d1_gb1))*100)/(abs(d1_gb1))
;p11 = ((abs(d3_gb2) - abs(d1_gb2))*100)/(abs(d1_gb2))
;p12 = ((abs(d3_gb3) - abs(d1_gb3))*100)/(abs(d1_gb3))

;p13 = ((abs(d4_na1) - abs(d2_na1))*100)/(abs(d2_na1))
;p14 = ((abs(d4_na2) - abs(d2_na2))*100)/(abs(d2_na2))
;p15 = ((abs(d4_na3) - abs(d2_na3))*100)/(abs(d2_na3))

;p16 = ((abs(d4_po1) - abs(d2_po1))*100)/(abs(d2_po1))
;p17 = ((abs(d4_po2) - abs(d2_po2))*100)/(abs(d2_po2))
;p18 = ((abs(d4_po3) - abs(d2_po3))*100)/(abs(d2_po3))

;p19 = ((abs(d4_so1) - abs(d2_so1))*100)/(abs(d2_so1))
;p20 = ((abs(d4_so2) - abs(d2_so2))*100)/(abs(d2_so2))
;p21 = ((abs(d4_so3) - abs(d2_so3))*100)/(abs(d2_so3))

;p22 = ((abs(d4_gb1) - abs(d2_gb1))*100)/(abs(d2_gb1))
;p23 = ((abs(d4_gb2) - abs(d2_gb2))*100)/(abs(d2_gb2))
;p24 = ((abs(d4_gb3) - abs(d2_gb3))*100)/(abs(d2_gb3))



p1 = (((d3_na1) - (d1_na1))*100)/((d1_na1))
p2 = (((d3_na2) - (d1_na2))*100)/((d1_na2))
p3 = (((d3_na3) - (d1_na3))*100)/((d1_na3))

p4 = (((d3_po1) - (d1_po1))*100)/((d1_po1))
p5 = (((d3_po2) - (d1_po2))*100)/((d1_po2))
p6 = (((d3_po3) - (d1_po3))*100)/((d1_po3))


p7 = (((d3_so1) - (d1_so1))*100)/((d1_so1))
p8 = (((d3_so2) - (d1_so2))*100)/((d1_so2))
p9 = (((d3_so3) - (d1_so3))*100)/((d1_so3))


p10 = (((d3_gb1) - (d1_gb1))*100)/((d1_gb1))
p11 = (((d3_gb2) - (d1_gb2))*100)/((d1_gb2))
p12 = (((d3_gb3) - (d1_gb3))*100)/((d1_gb3))


p13 = (((d4_na1) - (d2_na1))*100)/((d2_na1))
p14 = (((d4_na2) - (d2_na2))*100)/((d2_na2))
p15 = (((d4_na3) - (d2_na3))*100)/((d2_na3))

p16 = (((d4_po1) - (d2_po1))*100)/((d2_po1))
p17 = (((d4_po2) - (d2_po2))*100)/((d2_po2))
p18 = (((d4_po3) - (d2_po3))*100)/((d2_po3))


p19 = (((d4_so1) - (d2_so1))*100)/((d2_so1))
p20 = (((d4_so2) - (d2_so2))*100)/((d2_so2))
p21 = (((d4_so3) - (d2_so3))*100)/((d2_so3))


p22 = (((d4_gb1) - (d2_gb1))*100)/((d2_gb1))
p23 = (((d4_gb2) - (d2_gb2))*100)/((d2_gb2))
p24 = (((d4_gb3) - (d2_gb3))*100)/((d2_gb3))


p25 = (((d6_na1) - (d5_na1))*100)/((d5_na1))
p26 = (((d6_na2) - (d5_na2))*100)/((d5_na2))
p27 = (((d6_na3) - (d5_na3))*100)/((d5_na3))

p28 = (((d6_po1) - (d5_po1))*100)/((d5_po1))
p29 = (((d6_po2) - (d5_po2))*100)/((d5_po2))
p30 = (((d6_po3) - (d5_po3))*100)/((d5_po3))

p31 = (((d6_so1) - (d5_so1))*100)/((d5_so1))
p32 = (((d6_so2) - (d5_so2))*100)/((d5_so2))
p33 = (((d6_so3) - (d5_so3))*100)/((d5_so3))

p34 = (((d6_gb1) - (d5_gb1))*100)/((d5_gb1))
p35 = (((d6_gb2) - (d5_gb2))*100)/((d5_gb2))
p36 = (((d6_gb3) - (d5_gb3))*100)/((d5_gb3))

p36 = (((d3_na4) - (d1_na4))*100)/((d1_na4))
p37 = (((d3_po4) - (d1_po4))*100)/((d1_po4))
p38 = (((d3_so4) - (d1_so4))*100)/((d1_so4))
p39 = (((d3_gb4) - (d1_gb4))*100)/((d1_gb4))


p40 = (((d4_na4) - (d2_na4))*100)/((d2_na4))
p41 = (((d4_po4) - (d2_po4))*100)/((d2_po4))
p42 = (((d4_so4) - (d2_so4))*100)/((d2_so4))
p43 = (((d4_gb4) - (d2_gb4))*100)/((d2_gb4))

p44 = (((d6_na4) - (d5_na4))*100)/((d5_na4))
p45 = (((d6_po4) - (d5_po4))*100)/((d5_po4))
p46 = (((d6_so4) - (d5_so4))*100)/((d5_so4))
p47 = (((d6_gb4) - (d5_gb4))*100)/((d5_gb4))



print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
print(p12)
print(p13)
print(p14)
print(p15)
print(p16)
print(p17)
print(p18)
print(p19)
print(p20)
print(p21)
print(p22)
print(p23)
print(p24)
print(p25)
print(p26)
print(p27)
print(p28)
print(p29)
print(p30)
print(p31)
print(p32)
print(p33)
print(p34)
print(p35)
print(p36)
print(p37)
print(p38)
print(p39)
print(p40)
print(p41)
print(p42)
print(p43)
print(p44)
print(p45)
print(p46)
print(p47)


print(d3_po1)
print(d1_po1)

print(d3_po2)
print(d1_po2)

resL3@gsnDraw              = False                ; don't draw yet
  resL3@gsnFrame             = False                ; don't advance frame yet
  resL3@tiXAxisFontHeightF   = 0.023
   resL3@tiYAxisFontHeightF   = 0.023
   resL3@tmXBLabelFontHeightF  = 0.021
   resL3@tmYLLabelFontHeightF  = 0.021
   resL3@tiXAxisString = "Latitude"  ; axis string
;  resL3@xyMonoDashPattern = True
  resL3@gsnLeftStringParallelPosF = 0.02
  ;resL3@txFontHeightF      = 0.015
  resL3@pmLegendDisplayMode    = "Always"            ; turn on legend
  resL3@lgPerimOn              = True               ; turn off box around
  resL3@xyMarkLineModes        = (/"Lines"/)         ; line style
  resL3@xyLineThicknesses      = (/2.,2.,2.,2./)     ; line thickness
;  resL3@xyLineColors           =  (/"darkorange1","forestgreen","gray40","dodgerblue2"/)
 resL3@xyLineColors           =  (/"firebrick3","blue4","coral1","dodgerblue"/)
 resL3@lgPerimOn              = False                 ; no box around
resL3@lgJustification        = "BottomRight"
  resL3@pmLegendOrthogonalPosF = 0     ; Move into plot
  resL3@pmLegendParallelPosF   = 1.0       ; Move to right
  resL3@pmLegendWidthF         = 0.3       ; Decrease width
  resL3@pmLegendHeightF        = 0.12       ; Decrease height
  resL3@lgBoxMinorExtentF      = 0.12       ; Shorten the legend lines
  resL3@pmLegendFontHeightF     = 0.024
  resL3@lgLabelFontHeightF     = 0.0215     ; Change the font size
 resL3@xyDashPatterns   = (/1,1,0,0/)
  resL3@xyExplicitLegendLabels = (/"PI-w/oIDEMIX","LGM-w/oIDEMIX","PI-w/IDEMIX", "LGM-w/IDEMIX"/)
   resL3@lgLabelFontThicknessF = 1.2
  resL3@lgLineLabelFontThicknessF = 1.4
resL3@lgLineThicknesses = 1.3
;resL3@gsnXRefLine = 1000
resL3@gsnXRefLineDashPattern = 1
resL3@gsnXRefLineColor = "black"
;resL3@gsnXRefLineThicknessF = 2.0
resL3@tiXAxisString = "N2" +"  "+"10^-6 [1/s^2]"  ; axis string
resL3@tiYAxisString = "Depth" +"  "+"[m]"  ; axis string
resL3@tiMainString    = ""   ; add titles

resL3@tiXAxisString = ""  ; axis string
resL3@tmYLMode        = "Explicit"                     ; Force labels where we want them
;resL3@tmYLValues      = ispan(-2,2,0.2)
;resL3@tmYLLabels = ispan(30,-30,10)

;  resL3@trXMaxF           = 0                    ; axis max
;  resL3@trXMinF           = 2                    ; axis min
  resL3@tmXTLabelConstantSpacingF = 0.5
;  resL3@xyLineColors      = "red"                    ; line color
; resL3@tmXBTickSpacing = 0.5
   resL3@trYReverse        = True

  resL3@lgTitleConstantSpacingF = 0.2



resL3_v2 = resL3
;resL3@tiXAxisString = "N2" +"  "+"[1/s^2]"  ; axis string



resL4 = resL3
resL4@tiXAxisString = "Vertical Diffusivitiy" +"  "+"[m^2/s]"  ; axis string
delete(resL4@pmLegendOrthogonalPosF)
delete(resL4@pmLegendParallelPosF)
resL4@pmLegendOrthogonalPosF = 0.05     ; Move into plot
resL4@pmLegendParallelPosF   = 1.0       ; Move to right


delete(resL4@trXMinF)
delete(resL4@trXMaxF)
;delete(resL4@xyExplicitLegendLabels)
;delete(resL4@tiYAxisString)



resL4@trXMaxF           = -2                    ; axis max
resL4@trXMinF           = -6                    ; axis min



resL4_v2 = resL4
delete(resL4_v2@tiYAxisString)


 data1      = new((/4,59/),double)

 printVarSummary(data1)
 printVarSummary(pi_n2_tamas_mean_gb)

 data1(0,:) = pi_n2_tamas_mean_gb(:)
 data1(1,:) = lgm_n2_tamas_mean_gb(:)
 data1(2,:) = pi_n2_mean_gb(:)
 data1(3,:) = lgm_n2_mean_gb(:)


 data2      = new((/4,59/),double)

 printVarSummary(data2)
 printVarSummary(pi_n2_tamas_mean_na)

 data2(0,:) = pi_n2_tamas_mean_na(:)
 data2(1,:) = lgm_n2_tamas_mean_na(:)
 data2(2,:) = pi_n2_mean_na(:)
 data2(3,:) = lgm_n2_mean_na(:)




data3      = new((/4,59/),double)

 printVarSummary(data3)
 printVarSummary(pi_n2_tamas_mean_na)

 data3(0,:) = pi_n2_tamas_mean_so(:)
 data3(1,:) = lgm_n2_tamas_mean_so(:)
 data3(2,:) = pi_n2_mean_so(:)
 data3(3,:) = lgm_n2_mean_so(:)



data4      = new((/4,60/),double)

 printVarSummary(data4)
 printVarSummary(pi_kvmix_tamas_mean_gb_m2_log)

 data4(0,:) = pi_kvmix_tamas_mean_gb_m2_log
 data4(1,:) = lgm_kvmix_tamas_mean_gb_m2_log
 data4(2,:) = pi_kvmix_mean_gb_m2_log
 data4(3,:) = lgm_kvmix_mean_gb_m2_log


data5      = new((/4,60/),double)

 printVarSummary(data5)
 printVarSummary(pi_n2_tamas_mean_na)

 data5(0,:) = pi_kvmix_tamas_mean_na_m2_log
 data5(1,:) = lgm_kvmix_tamas_mean_na_m2_log
 data5(2,:) = pi_kvmix_mean_na_m2_log
 data5(3,:) = lgm_kvmix_mean_na_m2_log



data6      = new((/4,60/),double)

 printVarSummary(data6)
 printVarSummary(pi_n2_tamas_mean_na)

 data6(0,:) = pi_kvmix_tamas_mean_so_m2_log
 data6(1,:) = lgm_kvmix_tamas_mean_so_m2_log
 data6(2,:) = pi_kvmix_mean_so_m2_log
 data6(3,:) = lgm_kvmix_mean_so_m2_log


depth           = pop_pi->z_t(0:58)
depth_m2          = depth*1.E-2       ; convert from cm to m

depth2           = pop_pi->z_t(0:59)
depth2_m2          = depth2*1.E-2       ; convert from cm to m


 plot(3)  = gsn_csm_xy(wks,data1,depth_m2,resL3) ; create plot
 plot(4) = gsn_csm_xy(wks,data2,depth_m2,resL3_v2)
 plot(5)  = gsn_csm_xy(wks,data3,depth_m2,resL3_v2) ; create plot
 plot(0) = gsn_csm_xy(wks,data4,depth2_m2,resL4)
 plot(1)  = gsn_csm_xy(wks,data5,depth2_m2,resL4_v2) ; create plot
 plot(2) = gsn_csm_xy(wks,data6,depth2_m2,resL4_v2)


 resP = True
  resP@gsnMaximize        = False
  resP@gsnFrame         = False               ; save panel until both ready
  gsn_panel(wks, plot, (/2,3/), resP)


txres2 = True
txres2@txFontHeightF = 0.011
txres2@txFont   = "helvetica-bold"

gsn_text_ndc(wks,"Global Ocean",0.19,0.9,txres2)
gsn_text_ndc(wks,"Atlantic Ocean",0.52,0.9,txres2)
gsn_text_ndc(wks,"Southern Ocean",0.86,0.9,txres2)

frame(wks)

end

