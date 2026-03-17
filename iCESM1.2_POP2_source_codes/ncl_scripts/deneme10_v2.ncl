

;========================================================================
; moc_density_remap.ncl
; Compute MOC in density space (Psi(lat, sigma)) and remap to depth.
; Assumptions:
; - VVEL(time,z,lat,lon) in m/s
; - dz(z) in m (1D)
; - DXT(lat,lon) in m (zonal grid spacing)
; - PD2 or PD(time,z,lat,lon) in g/cm3 for potential density referenced to 2000 dbar (or PD)
; - We compute sigma2 = PD2*1000 - 1000 (kg/m3 - 1000)
; - Time-mean over last Nyears (Nyears default 100)
;========================================================================

begin


   ;pop_lgm_tamas = addfile ("/mnt/d/lgm_data/buo_test_cmake.pop.h.1550_1599.nc", "r")

   ;pop_pi = addfile ("/mnt/d/lgm_data/buo_pi_idemix_isotope_continued_4900_4999_timmean.nc", "r")
   ;pop_lgm = addfile ("/mnt/d/lgm_data/buo_lgm_idemix_isotope_continued_5800_5899_timmean.nc", "r"

 ;pop_pi_tamas = addfile ("/mnt/d/lgm_data2/montly_parameters.piControl.p875.001.pop.h.4150_4199_timmean.nc", "r")
pop_pi = addfile ("/mnt/d/lgm_data/streamfunction_pi_idemix_4900_4604_4703_timmean.nc", "r")
pop_lgm = addfile ("/mnt/d/lgm_data/streamfunction_lgm_idemix_isotope_continued.5800_5899_timmean.nc", "r")
pop_lgm_tamas = addfile ("/mnt/d/lgm_data/streamfunction_test_cmake.1780_1879_timmean.nc", "r")


   pop_pi_tamas = addfile ("/mnt/d/lgm_data/streamfunctionpiControl.p875.001.pop.h.4100_4199.deflated2_timmean.nc", "r")
   


   pop_pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")
   pop_lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")


pi_tarea = addfile ("/mnt/d/lgm_data/pi_tarea.nc", "r")
pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")

pi_kmt = addfile ("/mnt/d/lgm_data/kmt_pi.nc", "r")
pi_grid = addfile ("/mnt/d/lgm_data/pi_grid_spacing_and_thickness.nc", "r")

lgm_tarea = addfile ("/mnt/d/lgm_data/lgm_tarea.nc", "r")
lgm_region = addfile ("/mnt/d/lgm_data/region_mask_lgm.nc", "r")

lgm_kmt = addfile ("/mnt/d/lgm_data/kmt_lgm.nc", "r")
lgm_grid = addfile ("/mnt/d/lgm_data/lgm_grid_spacing_and_thickness.nc", "r")

pi_atlantic_mask = addfile ("/mnt/d/lgm_data/pi_atl_mask_miss.nc", "r")

lgm_atlantic_mask = addfile ("/mnt/d/lgm_data/lgm_atl_mask_miss.nc", "r")


pi_indop_mask = addfile ("/mnt/d/lgm_data/pi_indo_pacific_masknew.nc", "r")

lgm_indop_mask = addfile ("/mnt/d/lgm_data/lgm_indo_pacific_masknew.nc", "r")


wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"deneme10_v2")

   
   cmap = read_colormap_file("temp_diff_18lev")
   ;cmap = read_colormap_file("cmocean_delta")
   ncmap = dimsizes(cmap(:,0))


 ;file      = "pop2_output.nc"
  ;Vname     = "VVEL"      ; try VVEL or VVEL2
  ;dzname    = "dz"        ; vertical thickness (1D)
  ;dxname    = "DXT"       ; zonal grid spacing (2D)
  ;PDnames   = (/"PD2","PD"/)  ; try PD2 first (sigma2), else PD
  Nyears    = 100
  nbins     = 301         ; number of sigma bins
  ;sigma_min = 25.0        ; kg/m3-1000
  ;sigma_max = 28.0
  ;----------------------------------------------


  z_t1           = pop_pi_tamas->z_t

  vvel1 = pop_pi_tamas->VVEL(0,:,:,:)          ; [z, y, x], m/s
  dxt1  = pi_grid->DXT(:,:)               ; cm
  dz1   = pi_grid->dz(:)                  ; cm
  kmt1  = pi_kmt->KMT(:,:)               ; ocean depth mask

 ; Convert cm → m
  dxt1 = dxt1 * 0.01
  dz1  = dz1  * 0.01
  vvel1 = vvel1 * 0.01

  z_t2           = pop_lgm_tamas->z_t

  vvel2 = pop_lgm_tamas->VVEL(0,:,:,:)          ; [z, y, x], m/s
  dxt2  = lgm_grid->DXT(:,:)               ; cm
  dz2   = lgm_grid->dz(:)                  ; cm
  kmt2  = lgm_kmt->KMT(:,:)               ; ocean depth mask

 ; Convert cm → m
  dxt2 = dxt2 * 0.01
  dz2  = dz2  * 0.01
  vvel2 = vvel2 * 0.01

  
 vvel3 = pop_pi->VVEL(0,:,:,:)          ; [z, y, x], m/s
 vvel3 = vvel3 * 0.01 


 vvel4 = pop_lgm->VVEL(0,:,:,:)          ; [z, y, x], m/s
 vvel4 = vvel4 * 0.01

  nlat1 = dimsizes(vvel1(0,:,0))
  nlon1 = dimsizes(vvel1(0,0,:))
  nlev1 = dimsizes(vvel1(:,0,0))

  nlat2 = dimsizes(vvel2(0,:,0))
  nlon2 = dimsizes(vvel2(0,0,:))
  nlev2 = dimsizes(vvel2(:,0,0))


   pi_atlantic_mask1 = pi_atlantic_mask->mask
    lgm_atlantic_mask1 = lgm_atlantic_mask->mask
    pi_indop_mask1 = pi_indop_mask->mask
    lgm_indop_mask1 = lgm_indop_mask->mask

 nlatcp = dimsizes(pi_atlantic_mask1(:,0))
 nloncp = dimsizes(pi_atlantic_mask1(0,:))

 nlatcl = dimsizes(lgm_atlantic_mask1(:,0))
 nloncl = dimsizes(lgm_atlantic_mask1(0,:))

vvel1v = new((/nlev1, nlat1, nlon1/), float)
vvelv = 0.0

vvel2v = new((/nlev2, nlat2, nlon2/), float)
vvel2v = 0.0

vvel3v = new((/nlev1, nlat1, nlon1/), float)
vvel3v = 0.0

vvel4v = new((/nlev2, nlat2, nlon2/), float)
vvel4v = 0.0

vvel1vv = vvel1v
vvel2vv = vvel2v
vvel3vv = vvel3v
vvel4vv = vvel4v

printVarSummary(vvel1v)
printVarSummary(vvel1)


region_atlantic_pi= new((/384,320/),float)
region_atlantic_lgm= new((/384,320/),float)

region_indop_pi= new((/384,320/),float)
region_indop_lgm= new((/384,320/),float)


;### ATLANTIC OCEAN;###

region_atlantic_pi        = where(pi_atlantic_mask1.eq.1,0,1)
region_atlantic_lgm        = where(lgm_atlantic_mask1.eq.1,0,1)

region_indop_pi        = where(pi_indop_mask1.eq.1,0,1)
region_indop_lgm        = where(lgm_indop_mask1.eq.1,0,1)


 ;do k = 0, nlev1
  vvel1v = mask(vvel1,region_atlantic_pi,0)
  vvel2v = mask(vvel2,region_atlantic_lgm,0)
  vvel3v = mask(vvel3,region_atlantic_pi,0)
  vvel4v = mask(vvel4,region_atlantic_lgm,0)
 ;end do

printVarSummary(vvel1v)
print(max(vvel1v))
print(min(vvel1v))

 ;do k = 0, nlev1
  vvel1vv = mask(vvel1,region_indop_pi,0)
  vvel2vv = mask(vvel2,region_indop_lgm,0)
  vvel3vv = mask(vvel3,region_indop_pi,0)
  vvel4vv = mask(vvel4,region_indop_lgm,0)
 ;end do

printVarSummary(vvel1vv)


  ;nlat1 = dimsizes(vvel1(0,:,0))
  ;nlon1 = dimsizes(vvel1(0,0,:))
  ;nlev1 = dimsizes(vvel1(:,0,0))

  ;nlat2 = dimsizes(vvel2(0,:,0))
  ;nlon2 = dimsizes(vvel2(0,0,:))
  ;nlev2 = dimsizes(vvel2(:,0,0))


; Create an array for volume transport (m3/s)
  transport1 = new((/nlev1, nlat1, nlon1/), float)
  transport1 = 0.0

  ; Calculate transport per grid cell
  do k = 0, nlev1-1
    do j = 0, nlat1-1
      do i = 0, nlon1-1
        if (k .lt. kmt1(j,i)) then
          transport1(k,j,i) = tofloat(vvel1v(k,j,i) * dxt1(j,i) * dz1(k))
        else
          transport1(k,j,i) = 0.0
        end if
      end do
    end do
  end do

 ; Sum over depth and longitude for meridional transport (e.g., across latitude)

  
   printVarSummary(transport1)
   print(max(transport1))
   print(min(transport1))

  
  transport2 = new((/nlev2, nlat2, nlon2/), float)
  transport2 = 0.0

  ; Calculate transport per grid cell
  do k = 0, nlev2-1
    do j = 0, nlat2-1
      do i = 0, nlon2-1
        if (k .lt. kmt2(j,i)) then
          transport2(k,j,i) = tofloat(vvel2v(k,j,i) * dxt2(j,i) * dz2(k))
        else
          transport2(k,j,i) = 0.0
        end if
      end do
    end do
  end do

 ; Sum over depth and longitude for meridional transport (e.g., across latitude)

  
   printVarSummary(transport2)
   print(max(transport2))
   print(min(transport2))

 transport3 = new((/nlev1, nlat1, nlon1/), float)
  transport3 = 0.0

  ; Calculate transport per grid cell
  do k = 0, nlev1-1
    do j = 0, nlat1-1
      do i = 0, nlon1-1
        if (k .lt. kmt1(j,i)) then
          transport3(k,j,i) = tofloat(vvel3v(k,j,i) * dxt1(j,i) * dz1(k))
        else
          transport3(k,j,i) = 0.0
        end if
      end do
    end do
  end do

 ; Sum over depth and longitude for meridional transport (e.g., across latitude)

  
   printVarSummary(transport3)
   print(max(transport3))
   print(min(transport3))

 transport4 = new((/nlev2, nlat2, nlon2/), float)
  transport4 = 0.0

  ; Calculate transport per grid cell
  do k = 0, nlev2-1
    do j = 0, nlat2-1
      do i = 0, nlon2-1
        if (k .lt. kmt2(j,i)) then
          transport4(k,j,i) = tofloat(vvel4v(k,j,i) * dxt2(j,i) * dz2(k))
        else
          transport4(k,j,i) = 0.0
        end if
      end do
    end do
  end do

 ; Sum over depth and longitude for meridional transport (e.g., across latitude)

  
   printVarSummary(transport4)
   print(max(transport4))
   print(min(transport4))


 


  ; trans_y = new((/nlat/), float)

  ; do j = 0, nlat-1
  ;  trans_y(j) = sum(transport(:,j,:))
  ; end do
  
  ; Units: m3/s

  ; Optional: Convert to Sverdrup (1 Sv = 1e6 m3/s)
 
  ; trans_y_Sv = trans_y / 1e6

  ;printVarSummary(trans_y_Sv)
  ;print("Meridional volume transport (Sv):")
  ;print(trans_y_Sv)
  ;print(max(trans_y_Sv))
  ;print(min(trans_y_Sv))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PART 2 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  z_t1 = z_t1 * 0.01
  z_t2 = z_t2 * 0.01

  rho1   = pop_pi_tamas->PD(0,:,:,:)
  rho2   = pop_lgm_tamas->PD(0,:,:,:)
  rho3   = pop_pi->PD(0,:,:,:)
  rho4   = pop_lgm->PD(0,:,:,:)
   
  ;-----------------------------------------------------------------
  ; Define density bins (adjust to your data range)
  ;-----------------------------------------------------------------
  nbins = 200
  rhomin1 = min(rho1)
  rhomax1 = max(rho1)
  drho1 = (rhomax1 - rhomin1) / nbins
  rho_bins1 = fspan(rhomin1, rhomax1, nbins)

  rhomin2 = min(rho2)
  rhomax2 = max(rho2)
  drho2 = (rhomax1 - rhomin1) / nbins
  rho_bins2 = fspan(rhomin2, rhomax2, nbins)


  rhomin3 = min(rho3)
  rhomax3 = max(rho3)
  drho3 = (rhomax3 - rhomin3) / nbins
  rho_bins3 = fspan(rhomin3, rhomax3, nbins)


  rhomin4 = min(rho4)
  rhomax4 = max(rho4)
  drho4 = (rhomax4 - rhomin4) / nbins
  rho_bins4 = fspan(rhomin4, rhomax4, nbins)


 ;-----------------------------------------------------------------
  ; Initialize transport in density bins
  ;-----------------------------------------------------------------
  trans_rho1 = new((/nbins, nlat1/), float)
  trans_rho1 = 0.0

  trans_rho2 = new((/nbins, nlat2/), float)
  trans_rho2 = 0.0

  trans_rho3 = new((/nbins, nlat1/), float)
  trans_rho3 = 0.0

  trans_rho4 = new((/nbins, nlat2/), float)
  trans_rho4 = 0.0
 


 ;-----------------------------------------------------------------
  ; Loop through grid and bin by density
  ;-----------------------------------------------------------------
  do j = 0, nlat1-1
    do i = 0, nlon1-1
      do k = 0, nlev1-1
        if (k .lt. kmt1(j,i)) then
          rho_val1 = rho1(k,j,i)
          bin_index = toint((rho_val1 - rhomin1) / drho1)
          if (bin_index.ge.0 .and. bin_index.lt.nbins) then
            vol = vvel1(k,j,i) * dxt1(j,i) * dz1(k)    ; m3/s per grid cell   
            trans_rho1(bin_index,j) = trans_rho1(bin_index,j) + tofloat(vol)
          end if
        end if
      end do
    end do
  end do

  printVarSummary(trans_rho1)
  print(max(trans_rho1))
  print(min(trans_rho1))

  do j = 0, nlat2-1
    do i = 0, nlon2-1
      do k = 0, nlev2-1
        if (k .lt. kmt2(j,i)) then
          rho_val2 = rho2(k,j,i)
          bin_index = toint((rho_val2 - rhomin2) / drho2)
          if (bin_index.ge.0 .and. bin_index.lt.nbins) then
            vol = vvel2(k,j,i) * dxt2(j,i) * dz2(k)    ; m3/s per grid cell   
            trans_rho2(bin_index,j) = trans_rho2(bin_index,j) + tofloat(vol)
          end if
        end if
      end do
    end do
  end do

  printVarSummary(trans_rho2)
  print(max(trans_rho2))
  print(min(trans_rho2))


do j = 0, nlat1-1
    do i = 0, nlon1-1
      do k = 0, nlev1-1
        if (k .lt. kmt1(j,i)) then
          rho_val3 = rho3(k,j,i)
          bin_index = toint((rho_val3 - rhomin3) / drho1)
          if (bin_index.ge.0 .and. bin_index.lt.nbins) then
            vol = vvel3(k,j,i) * dxt1(j,i) * dz1(k)    ; m3/s per grid cell   
            trans_rho3(bin_index,j) = trans_rho3(bin_index,j) + tofloat(vol)
          end if
        end if
      end do
    end do
  end do

  printVarSummary(trans_rho3)
  print(max(trans_rho3))
  print(min(trans_rho3))

  do j = 0, nlat2-1
    do i = 0, nlon2-1
      do k = 0, nlev2-1
        if (k .lt. kmt2(j,i)) then
          rho_val4 = rho4(k,j,i)
          bin_index = toint((rho_val4 - rhomin4) / drho2)
          if (bin_index.ge.0 .and. bin_index.lt.nbins) then
            vol = vvel4(k,j,i) * dxt2(j,i) * dz2(k)    ; m3/s per grid cell   
            trans_rho4(bin_index,j) = trans_rho4(bin_index,j) + tofloat(vol)
          end if
        end if
      end do
    end do
  end do

  printVarSummary(trans_rho4)
  print(max(trans_rho4))
  print(min(trans_rho4))




;-----------------------------------------------------------------
  ; Compute cumulative streamfunction in density space
  ;-----------------------------------------------------------------
  ; Integrate from high to low density (deep → surface)
  psi_rho1 = new((/nbins, nlat1/), float)
  psi_rho1 = 0.0
  psi_rho21 = new((/nbins, nlat1/), float)
  psi_rho21 = 0.0


   do j = 0, nlat1-1
    psi_rho21(:,j) = dim_cumsum(-trans_rho1(::-1,j), 0) / 1e6  ; Sv
    psi_rho1(:,j) = psi_rho21(::-1,j)
  end do


  psi_rho2 = new((/nbins, nlat2/), float)
  psi_rho2 = 0.0
  psi_rho22 = new((/nbins, nlat2/), float)
  psi_rho22 = 0.0


   do j = 0, nlat2-1
    psi_rho22(:,j) = dim_cumsum(-trans_rho2(::-1,j), 0) / 1e6  ; Sv
    psi_rho2(:,j) = psi_rho22(::-1,j)
  end do


  psi_rho3 = new((/nbins, nlat1/), float)
  psi_rho3 = 0.0
  psi_rho23 = new((/nbins, nlat1/), float)
  psi_rho23 = 0.0


   do j = 0, nlat1-1
    psi_rho23(:,j) = dim_cumsum(-trans_rho3(::-1,j), 0) / 1e6  ; Sv
    psi_rho3(:,j) = psi_rho23(::-1,j)
  end do


  psi_rho4 = new((/nbins, nlat2/), float)
  psi_rho4 = 0.0
  psi_rho24 = new((/nbins, nlat2/), float)
  psi_rho24 = 0.0


   do j = 0, nlat2-1
    psi_rho24(:,j) = dim_cumsum(-trans_rho4(::-1,j), 0) / 1e6  ; Sv
    psi_rho4(:,j) = psi_rho24(::-1,j)
  end do




;  do j = 0, nlat-1
;    psi_rho(:,j) = dim_cumsum(-trans_rho(::-1,j), 0) / 1e6  ; Sv
;  end do

  psi_rho1@long_name = "Density streamfunction"
  psi_rho1@units = "Sv"

  printVarSummary(psi_rho1)
  print(max(psi_rho1))
  print(min(psi_rho1)) 

  psi_rho2@long_name = "Density streamfunction"
  psi_rho2@units = "Sv"

  printVarSummary(psi_rho2)
  print(max(psi_rho2))
  print(min(psi_rho2)) 
 
 
  psi_rho3@long_name = "Density streamfunction"
  psi_rho3@units = "Sv"

  printVarSummary(psi_rho3)
  print(max(psi_rho3))
  print(min(psi_rho3)) 

  psi_rho4@long_name = "Density streamfunction"
  psi_rho4@units = "Sv"

  printVarSummary(psi_rho4)
  print(max(psi_rho4))
  print(min(psi_rho4)) 
 
    lat1 = pop_pi_tamas->TLAT(:,0)


 ;-----------------------------------------------------------------
  ; Convert density bins to mean depths
  ;-----------------------------------------------------------------

  ;  pi_atlantic_mask1 = pi_atlantic_mask->mask
  ;  lgm_atlantic_mask1 = lgm_atlantic_mask->mask
  ;  pi_indop_mask1 = pi_indop_mask->mask
  ;  lgm_indop_mask1 = lgm_indop_mask->mask


;   psi_depth1_atl = mask(psi_rho1,pi_atlantic_mask1,0)  ; same streamfunction, but now plotted vs depth
 

   rho_bins6 = fspan(34, 36, nbins)


  psi_depth1_atl = psi_rho1
  psi_depth1_atl!0 = "density"
   psi_depth1_atl!1 = "lat"
 ;  psi_depth1_atl&density = (rho_bins1*1000)-1000
  psi_depth1_atl&density = rho_bins6 
 psi_depth1_atl&lat = lat1


  printVarSummary(psi_depth1_atl)
  print(max(psi_depth1_atl))
  print(min(psi_depth1_atl))

   

    psi_depth1_ip = psi_rho1  ; same streamfunction, but now plotted vs depth
   psi_depth1_ip!0 = "density"
   psi_depth1_ip!1 = "lat"
   ;psi_depth1_ip&density = (rho_bins1*1000)-1000
   psi_depth1_ip&density = rho_bins6
   psi_depth1_ip&lat = lat1


  printVarSummary(psi_depth1_ip)
  print(max(psi_depth1_ip))
  print(min(psi_depth1_ip))



   psi_depth2_atl = psi_rho2   ; same streamfunction, but now plotted vs depth
   psi_depth2_atl!0 = "density"
   psi_depth2_atl!1 = "lat"
   ;psi_depth2_atl&density = (rho_bins2*1000)-1000
   psi_depth2_atl&density = rho_bins6
   psi_depth2_atl&lat = lat1


  printVarSummary(psi_depth2_atl)
  print(max(psi_depth2_atl))
  print(min(psi_depth2_atl))


   psi_depth2_ip = psi_rho2  ; same streamfunction, but now plotted vs depth
   psi_depth2_ip!0 = "density"
   psi_depth2_ip!1 = "lat"
   ;psi_depth2_ip&density = (rho_bins2*1000)-1000
    psi_depth2_ip&density = rho_bins6
   psi_depth2_ip&lat = lat1


  printVarSummary(psi_depth2_ip)
  print(max(psi_depth2_ip))
  print(min(psi_depth2_ip))


   psi_depth3_atl = psi_rho3   ; same streamfunction, but now plotted vs depth
   psi_depth3_atl!0 = "density"
   psi_depth3_atl!1 = "lat"
   ;psi_depth3_atl&density = (rho_bins3*1000)-1000
   psi_depth3_atl&density = rho_bins6
    psi_depth3_atl&lat = lat1


  printVarSummary(psi_depth3_atl)
  print(max(psi_depth3_atl))
  print(min(psi_depth3_atl))


 psi_depth3_ip = psi_rho3  ; same streamfunction, but now plotted vs depth
   psi_depth3_ip!0 = "density"
   psi_depth3_ip!1 = "lat"
   ;psi_depth3_ip&density = (rho_bins3*1000)-1000
  psi_depth3_ip&density = rho_bins6
   psi_depth3_ip&lat = lat1


  printVarSummary(psi_depth1_ip)
  print(max(psi_depth1_ip))
  print(min(psi_depth1_ip))



   psi_depth4_atl =  psi_rho4   ; same streamfunction, but now plotted vs depth
   psi_depth4_atl!0 = "density"
   psi_depth4_atl!1 = "lat"
   ;psi_depth4_atl&density = (rho_bins4*1000)-1000
   psi_depth4_atl&density = rho_bins6
   psi_depth4_atl&lat = lat1


  printVarSummary(psi_depth4_atl)
  print(max(psi_depth4_atl))
  print(min(psi_depth4_atl))

   psi_depth4_ip = psi_rho4  ; same streamfunction, but now plotted vs depth
   psi_depth4_ip!0 = "density"
   psi_depth4_ip!1 = "lat"
   ;psi_depth4_ip&density = (rho_bins4*1000)-1000
   psi_depth4_ip&density = rho_bins6
   psi_depth4_ip&lat = lat1


  printVarSummary(psi_depth4_ip)
  print(max(psi_depth4_ip))
  print(min(psi_depth4_ip))

  printVarSummary(rho_bins4)
  print(max(rho_bins4-1000))
  print(min(rho_bins4-1000))


;-----------------------------------------------------------------
  ; Optional: plot streamfunction in depth coordinates
  ;-----------------------------------------------------------------


  res = True
    res@gsnDraw              = False                ; don't draw yet
  res@gsnFrame             = False                ; don't advance frame yet

  res@gsnMaximize = True
  res@cnFillOn = True             ; filled contours
  res@cnLinesOn = False           ; remove contour lines
  res@trYReverse = True           ; depth increasing downward
  ;res@tiMainString = "Meridional Overturning Streamfunction (Sv)"
  res@tiXAxisString = ""
  ;res@tiYAxisString = "Potential Density - 1000 (kg/m^3)"

  ; Optional: set contour levels manually
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -30
  res@cnMaxLevelValF = 30
  res@cnLevelSpacingF = 3

 ; res@vpWidthF           = 0.15
;res@vpHeightF          = 0.20

res@tmXBLabelsOn  = True
res@tmXBLabelFontHeightF = 0.005
res@tmYLLabelFontHeightF = 0.005
;res@cnFillPalette        = cmap
res@mpLandFillColor = "grey56"
res@cnMissingValFillColor = "grey56"
res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "
res@tmXBLabelFontHeightF = 0.004
res@tmYLLabelFontHeightF = 0.0038


pltres = True
mpres = True
pltres@PanelPlot  = True
res@Footer    = False
plot = new(4,graphic)
plots = new(2,graphic)
;plots2 = new(2,graphic)
res@tmYLLabelsOn  = True
res@tmXBLabelsOn  = True
;res@tmXBTickSpacingF     = 15.0
res@vpWidthF           = 0.115
res@vpHeightF          = 0.08
res@tmXBLabelFontHeightF = 0.004

res@tmYLLabelFontHeightF = 0.0038


res@lbLabelBarOn        = False
res@cnInfoLabelOn       = False
res@description = "  "

;res@gsnYAxisIrregular2Linear = False
;res@trYReverse           = False
res@lbLabelFontHeightF= 0.0035
;res@pmLabelBarParallelPosF = -0.01
;res@lbTitleString = "Buoyancy Flux [10^-6 kg/m2/s]"
res@lbTitleString = "[Sv]"
res@lbTitleFontHeightF= 0.004
;res@lbTitleOffsetF = 0.01

res@tmXBTickSpacingF = 30.0 ; North Atlantic
;res@tmXBTickSpacingF = 10.0
res@tmYLTickSpacingF = 30.0
res@cnMissingValFillPattern = "SolidFill"
res@cnMissingValFillColor = "grey56"

;res@pmLabelBarOrthogonalPosF = .03             ; position label bar

res@pmLabelBarOrthogonalPosF = .04
 res@cnLineLabelsOn  = False


res1 = res
res1@tmYLLabelsOn  = False
res1@tmXBLabelsOn  = True
  res1@tiXAxisString = ""
  res1@tiYAxisString = ""

res2 = res
;res2@tiXAxisString = "Latitude (degC)"
;res2@tiYAxisString = "Potential Density - 1000 (kg/m^3)"

res3 = res
;res3@tiXAxisString = "Latitude (degC)"
;res3@tiYAxisString = ""
res3@tmYLLabelsOn  = False
res3@tmXBLabelsOn  = True

sub1 = abs(psi_depth3_atl) - abs(psi_depth1_atl) 
sub2 = abs(psi_depth4_atl) - abs(psi_depth2_atl)


res5 = res
  res5@cnLevelSelectionMode = "ManualLevels"
  res5@cnMinLevelValF = -5
  res5@cnMaxLevelValF = 5
  res5@cnLevelSpacingF = 1

   
   sub1!0 = "density"
   sub1!1 = "lat"
   sub1&density = rho_bins6
   sub1&lat = lat1

   sub2!0 = "density"
   sub2!1 = "lat"
   sub2&density = rho_bins6
   sub2&lat = lat1




  ; Plot
  plot(0) = gsn_csm_contour(wks, psi_depth1_atl, res)
  plot(1) = gsn_csm_contour(wks, psi_depth2_atl, res1) 
  plot(2) = gsn_csm_contour(wks, psi_depth3_atl, res2)
  plot(3) = gsn_csm_contour(wks, psi_depth4_atl, res3)

 ;  plot(0) = gsn_csm_contour(wks, psi_depth1_ip, res)
 ;  plot(1) = gsn_csm_contour(wks, psi_depth2_ip, res)
 ;  plot(2) = gsn_csm_contour(wks, psi_depth3_ip, res)
 ;  plot(3) = gsn_csm_contour(wks, psi_depth4_ip, res)

   plots(0) = gsn_csm_contour(wks, sub1, res5)
   plots(1) = gsn_csm_contour(wks, sub2, res5)

resP = True
resP@gsnMaximize        = False
  resP@gsnFrame         = False               ; save panel until both ready
  resP@gsnPaperOrientation = "landscape"
  resP@gsnPanelLabelBar   = True
  resP@pmLabelBarWidthF   = 0.60
  resP@pmLabelBarHeightF = 0.043

;  resP@pmLabelBarOrthogonalPosF = -0.042
  resP@lbTitleString = "[Sv]"

  resP@lbTitleFontHeightF = 0.013
resP@gsnPanelFigureStringsFontHeightF = 0.012
resP@gsnPanelFigureStrings= (/"(a)", "(b)", "(c)", "(d)"/)

resP@gsnPanelTop = 0.95
resP@gsnPanelBottom = 0.61
;resP@gsnPanelFigureStringsFontHeightF = 0.012

resP@gsnPanelYWhiteSpacePercent = 3


gsn_panel(wks, plot, (/2,2/), resP)

resPs = True

resPs@gsnMaximize        = False
  resPs@gsnFrame         = False               ; save panel until both ready
  resPs@gsnPaperOrientation = "landscape"
  resPs@gsnPanelLabelBar   = True
  resPs@pmLabelBarWidthF   = 0.60
  resPs@pmLabelBarHeightF = 0.043

 resPs@lbTitleString = "[Sv]"

  resPs@lbTitleFontHeightF = 0.013
resPs@gsnPanelFigureStringsFontHeightF = 0.012

resPs@gsnPanelTop = 0.59
resPs@gsnPanelBottom = 0.33
resPs@gsnPanelFigureStrings= (/"(e)", "(f)"/)

resPs@gsnPanelYWhiteSpacePercent = 3

gsn_panel(wks, plots, (/1,2/), resPs)


txres2 = True
txres2@txFontHeightF = 0.012
txres2@txFont   = "helvetica-bold"

;gsn_text_ndc(wks,"N2",0.19,0.76,txres2)
;gsn_text_ndc(wks,"Temperature",0.52,0.76,txres2)
;gsn_text_ndc(wks,"Salt",0.86,0.76,txres2)

frame(wks)


end
