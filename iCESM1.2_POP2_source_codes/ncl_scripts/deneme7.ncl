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


wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"deneme7")

   
   cmap = read_colormap_file("temp_diff_18lev")
   ;cmap = read_colormap_file("cmocean_delta")
   ncmap = dimsizes(cmap(:,0))

;--- Read variables
  vvel = pop_pi_tamas->VVEL(0,:,:,:)        ; [time|optional, z_t, nlat, nlon]  (cm/s)
  rho  = pop_pi_tamas->PD(0,:,:,:)          ; potential density (g/cm3)
  lat  = pop_pi_tamas->TLAT        ; (nlat, nlon) degrees_north
  zt   = pop_pi_tamas->z_t         ; (z) cm
  dz   = pi_grid->dz          ; (z) cm
  dxt  = pi_grid->DXT         ; (nlon) cm

;--- Select time average if multiple timesteps
;  if (dimsizes(dimsizes(vvel)) .eq. 4) then
;    vvel = dim_avg_n(vvel,0)    ; average over time
;    rho  = dim_avg_n(rho,0)
;  end if

;--- Convert units
  vvel = vvel * 1.0e-2          ; cm/s → m/s
  dz   = dz   * 1.0e-2          ; cm → m
  dxt  = dxt  * 1.0e-2          ; cm → m
  zt   = zt   * 1.0e-2          ; cm → m
  rho  = (rho * 1000.0) - 1000          ; g/cm3 → kg/m3

;--- Define region 20°S–40°N
  lat1 = -20.
  lat2 =  40.
  lat_mask = where((lat(:,0).ge.lat1).and.(lat(:,0).le.lat2),1,0)
  jmask = ind(lat_mask.eq.1)
  nlats = dimsizes(jmask)

  print(jmask)

;--- Define vertical and zonal sizes
  nz   = dimsizes(zt)
  nlon = dimsizes(dxt)

;--- Density bins
  nbins = 301
  sigma_min = min(rho) - 0.01 
  sigma_max = max(rho) + 0.01 
  ;sigma_centers = fspan(sigma_min, sigma_max, nbins)

  sigma_edges = sigma_min + (sigma_max - sigma_min) * (ispan(0, nbins, 1) / nbins)  ; nbins+1 edges
  sigma_centers = 0.5 * (sigma_edges(0:nbins-1) + sigma_edges(1:nbins))


  printVarSummary(sigma_centers)
  print(sigma_centers)
  print(max(sigma_centers))
  print(min(sigma_centers))


  mins = min(sigma_centers)
   maxs = max(sigma_centers)

;--- Preallocate arrays
  vtrans_zonal = new((/nz, nlats/), float, 0.0)
  psi_sigma    = new((/nbins, nlats/), float, 0.0)

  printVarSummary(dxt)

;--- Compute zonal transport (integrate over longitude)

  print("Computing zonal volume transport...")
  print("error")
 ;--- vtrans_zonal: zonal integral of VVEL

 do j = 0, nlats-1
  jj = jmask(j)
  do k = 0, nz-1
    v_slice   = reshape(vvel(k,jj,:), (/dimsizes(vvel(k,jj,:))/))
    dxt_slice = reshape(dxt(jj,:), (/dimsizes(dxt(jj,:))/))

    v_slice = where(ismissing(v_slice), 0.0, v_slice)
    dxt_slice = where(ismissing(dxt_slice), 0.0, dxt_slice)

    vtrans_zonal(k,j) = tofloat(sum(v_slice * dxt_slice))   ; scalar → 2D array element
  end do
end do


;--- Conform dz to vtrans_zonal for later multiplication
  dz_2d = conform(vtrans_zonal, dz, 0)

;--- Compute zonal-mean potential density
  print("Computing zonal mean potential density...")


  rho_zonal = dim_avg_n_Wrap(rho(:,jmask,:),2)    ; (z, lat)

  printVarSummary(rho_zonal)
  print("rho_zonal")
  print(avg(rho_zonal))

  printVarSummary(sigma_centers)
  print("sigma_centers")
  print(avg(sigma_centers))


  print("vtrans_zonal")
  print(avg(vtrans_zonal))


;--- Bin transports into density space
;  print("Mapping transports into density bins...")
;  do jj = 0, nlats-1
;    do k = 0, nz-1
;      sigma_layer = tofloat(rho_zonal(k,jj))
;      if (.not. ismissing(sigma_layer)) then
;        if (sigma_layer.ge.min(sigma_centers) .and. sigma_layer.le.max(sigma_centers)) then
;          ib = dim_gte_n(sigma_centers, sigma_layer, 0)
;          if (.not. ismissing(ib)) then
;            psi_sigma(ib:,jj) = psi_sigma(ib:,jj) + vtrans_zonal(k,jj) * dz_2d(k,jj)
;          end if
;        end if
;      end if
;    end do
;  end do


psi_sigma  = new((/nbins, nlats/), float, 0.0)


do j = 0, nlats-1
  do k = 0, nz-1
    sigma_layer = rho_zonal(k,j)    ; zonal mean potential density at this depth/latitude
    if (.not. ismissing(sigma_layer)) then
      if (sigma_layer .ge. mins .and. sigma_layer .le. maxs) then
;     print("error")

ib = -1
do bb = 0, nbins-1
  if (sigma_centers(bb) .ge. sigma_layer) then
    ib = bb
    exit
  end if
end do

;   print("error")

; If sigma_layer exceeds max, assign to top bin

if (ib .eq. -1) then
  ib = nbins-1
end if

 psi_sigma = 0.0
; Add transport
do ii = ib, nbins-1
  psi_sigma(ii, j) = psi_sigma(ii, j) + (vtrans_zonal(k,j) * dz(k))
end do

print("error")

      end if
    end if
  end do
end do


printVarSummary(psi_sigma)
print(max(psi_sigma))
print(min(psi_sigma))


;do j = 0, nlats-1
;  do k = 0, nz-1
;    sigma_layer = rho_zonal(k,j)       ; already σ₂ space
;    if (.not. ismissing(sigma_layer)) then
;      if (sigma_layer.ge.sigma_min .and. sigma_layer.le.sigma_max) then
;       ib = ind(sigma_centers.ge.sigma_layer)
;        if (dimsizes(ib).gt.0) then
;          do ii = 0, dimsizes(ib)-1
;          psi_sigma(ib(ii), j) = psi_sigma(ib(ii), j) + vtrans_zonal(k,j) * dz(k)
;          end do
;        end if
;      end if
;    end if
;  end do
;end do


;--- Loop over latitude and depth


;do j = 0, nlats-1
;  do k = 0, nz-1
;    sigma_layer = rho_zonal(k,j)    ; zonal mean potential density at this depth/latitude
;    if (.not. ismissing(sigma_layer)) then
;      if (sigma_layer.ge.min(sigma_centers) .and. sigma_layer.le.max(sigma_centers)) then
;
;       print("error")
;
;        ;--- Manual search for first density bin >= sigma_layer
;
;        ib = -1                            ; source of the problem starts from here
;        do bb = 0, nbins-1
;          if (sigma_centers(bb) .ge. sigma_layer) then
;            ib = bb
;            exit
;          end if
;        end do


        ;--- Add transport to all bins >= sigma_layer

;        if (ib .ge. 0) then
;          psi_sigma(ib:, j) = psi_sigma(ib:, j) + vtrans_zonal(k,j) * dz(k)
;        end if
;
;      end if
;    end if
;  end do
;end do

print("Heves")
print(avg(vtrans_zonal))

;--- Convert to Sverdrups
  psi_sigma = psi_sigma * 1.0e-6      ; m3/s → Sv

;--- Plotting
  print("Plotting MOC(σ₂)...")
;  wks = gsn_open_wks("png","MOC_sigma2")

  res                    = True
  res@gsnMaximize        = True
  res@cnFillOn           = True
  res@cnLinesOn          = False
  res@cnLevelSelectionMode = "ManualLevels"
 ; res@cnMinLevelValF     = -20.
 ; res@cnMaxLevelValF     =  20.
 ; res@cnLevelSpacingF    =   2.
  res@gsnLeftString      = "MOC (σ₂ space)"
  res@gsnRightString     = "Sv"
  res@tiYAxisString      = "σ₂ (kg m~S~-3~N~)"
  res@tiXAxisString      = "Latitude (°N)"
  res@trYReverse         = True

  ;res@sfXArray           = lat(jmask)
  ;res@sfYArray           = sigma_centers
  
  
  
   
  ;lat_subset = lat(where(lat_mask, -1))
  
 ; mask1 = (lat .ge. -20) .and. (lat .le. 40)
 ; lat_indices = ind(mask1)
 ; lat_subset = lat(lat_indices)
 
lat1D = lat(:,0)
mask1  = (lat1D .ge. -20) .and. (lat1D .le. 40)
lat_indices = ind(mask1)
lat_subset  = lat1D(lat_indices)

  res@sfXArray = lat_subset
  res@sfYArray = tofloat(sigma_centers)

  plot = gsn_csm_contour(wks, psi_sigma, res)

  draw(plot)
  frame(wks)

  print("Done.")
end
