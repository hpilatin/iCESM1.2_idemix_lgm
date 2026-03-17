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

   wks  = gsn_open_wks(wks_type,"deneme8")

   
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
  rho  = (rho * 1000.0)           ; g/cm3 → kg/m3

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

   nlat = dimsizes(vvel(0,:,0))
  nlon = dimsizes(vvel(0,0,:))

;--- Density bins
  nbins = 301
;  sigma_min = min(rho) - 1000 
;  sigma_max = max(rho) - 1000 
;  sigma_centers = fspan(sigma_min, sigma_max, nbins)

;  printVarSummary(sigma_centers)
;  print(sigma_centers)


;--- Preallocate arrays
  vtrans_zonal = new((/nz, nlats/), float, 0.0)
  trans_bin   = new((/nbins, nlats/), float, 0.0)

Psi_sigma  = new((/nbins, nlats/), float, 0.0)


  printVarSummary(dxt)

;--- Compute zonal transport (integrate over longitude)

  print("Computing zonal volume transport...")

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


     sigma2_3d = rho(:,jmask,:) - 1000.0   ; [z,lat,lon] = dim_avg_n_Wrap(rho(:,jmask,:),2) - 1000   ; (z, lat)
     
     rho_zonal = dim_avg_n_Wrap(rho(:,jmask,:),2) - 1000


;  printVarSummary(sigma_centers)
;  print("sigma_centers")
;  print(avg(sigma_centers))


  print("vtrans_zonal")
  print(avg(vtrans_zonal))

  ;nbins = 301
  sigma_min = min(rho) - 1000 
  sigma_max = max(rho) - 1000 


   


;  printVarSummary(sigma_centers)
;  print(sigma_centers)

  sigma_edges = sigma_min + (sigma_max - sigma_min) * (ispan(0, nbins, 1) / nbins)  ; nbins+1 edges
  sigma_centers = 0.5 * (sigma_edges(0:nbins-1) + sigma_edges(1:nbins))
  nb = nbins

   ;---------------- Loop over grid cells and accumulate transport into sigma bins ----------------
  ; We go through all (k,j,i) and add vmo3(k,j,i) to bin where sigma2_3d(k,j,i) falls
   ; Now trans_bin(b,j) is the zonally-summed transport (m3/s) of grid cells at latitude j with sigma in bin b.


nlats = dimsizes(jmask)     ; already defined earlier

do j = 0, nlat-1
  do k = 0, nz-1
    s = rho_zonal(k,j)   ; or rho_zonal(k,j) if you already have zonal mean
    if (.not. ismissing(s)) then
      if (s.ge.sigma_min .and. s.lt.sigma_max) then
        b = toint( floor( (s - sigma_min) / ((sigma_max - sigma_min) / nb) ) )
        if (b.lt.0) then
          b = 0
        end if
        if (b.ge.nb) then
          b = nb-1
        end if
        trans_bin(b,j) = trans_bin(b,j) + (vtrans_zonal(k,j) * dz(k))
      end if
    end if
  end do
end do




  ;---------------- Compute cumulative overturning in density space ----------------
  ; Convention: integrate from dense -> light (sigma large -> small)
  ; Psi_density(lat, sigma_k) = sum_{b >= k} trans_bin(b,lat)  (i.e. cumulative sum from heavy to light)
  Psi_sigma = new((/nb,nlat/), 0.0)
  do j=0,nlat-1
    cum = 0.0
    do b = nb-1, 0, -1   ; start from heaviest bin
      cum = cum + trans_bin(b,j)
      Psi_sigma(b,j) = cum
    end do
  end do

; Convert to Sv for plotting
  Psi_sigma_Sv = Psi_sigma / 1e6
  Psi_sigma_Sv@long_name = "MOC streamfunction in density space (Sv)"
  ; sigma_centers array available

  ;---------------- Remap Psi_sigma (lat,sigma) to depth levels ----------------
  ; Approach: use zonal-mean sigma2 profile (sigma2_zonal[z,lat]); for each (z,lat) find sigma at that depth
  ; and then sample Psi_sigma at that sigma (linear interpolation across sigma bins).
  sigma2_zonal = dim_avg_n(sigma2_3d, 2)   ; [z,lat]

  ; Create output Psi_remap[z,lat]
 

 Psi_remap = new((/nz,nlat/), -9999.0)


  do j = 0, nlat-1
    ; build vector sigma_profile = sigma2_zonal(:,j) and corresponding psi_sigma vs sigma for this latitude
    sigma_profile = sigma2_zonal(:,j)    ; [z]
    ; for interpolation we need Psi_sigma as a function of sigma; get arrays for this latitude
    psi_vs_sigma = Psi_sigma_Sv(:,j)      ; [nb]

    ; For each depth level k, determine sigma_at_depth and then interpolate psi_vs_sigma to that sigma
    do k = 0, nz-1
      s_at_z = sigma_profile(k)
      if (ismissing(s_at_z)) then
        Psi_remap(k,j) = 0.0
      else
        ; if s_at_z outside bin range, set to edge value
        
        if (s_at_z.le.sigma_centers(0)) then
          Psi_remap(k,j) = psi_vs_sigma(0)
        else if (s_at_z.ge.sigma_centers(nb-1)) then
          Psi_remap(k,j) = psi_vs_sigma(nb-1)
        else
        
         ; find bin b such that sigma_centers(b) <= s_at_z < sigma_centers(b+1)
          ; linear search (nb usually a few hundred, so OK)
        
          found = False
          do b = 0, nb-2
            if (s_at_z.ge.sigma_centers(b) .and. s_at_z.lt.sigma_centers(b+1)) then
              
           ; linear interpolate psi between centers b and b+1
           
             frac = (s_at_z - sigma_centers(b)) / (sigma_centers(b+1) - sigma_centers(b) + 1e-12)
              Psi_remap(k,j) = psi_vs_sigma(b) + frac * (psi_vs_sigma(b+1) - psi_vs_sigma(b))
              found = True
              exit
            end if
          end do
          if (.not. found) then
            Psi_remap(k,j) = psi_vs_sigma(0)
          end if
        end if
      end if
    end if
   end do
  end do

  ; Psi_remap is in Sv (because psi_vs_sigma was in Sv)
  Psi_remap@long_name = "Psi remapped to depth (Sv)"

  ;---------------- Plotting (two panels): Psi_sigma (density space) and Psi_remap (depth) ----------------


  ; plot 1: density-space (sigma on y-axis) use sigma_centers vertical
  res1 = True
  res1@gsnMaximize = True
  res1@cnFillOn = True
  res1@trYReverse = False         ; sigma increases downward? choose what you prefer
  res1@sfXArray = lat
  res1@sfYArray = sigma_centers
  res1@tiMainString = "Psi(lat, sigma2) (Sv)"
  plot1 = gsn_csm_contour(wks, Psi_sigma_Sv, res1)

  ; plot 2: remapped to depth (depth on y)
  res2 = True
  res2@cnFillOn = True
  res2@trYReverse = True
  res2@sfXArray = lat
  res2@sfYArray = zt
  res2@tiMainString = "Psi remapped to depth (Sv)"
  plot2 = gsn_csm_contour(wks, Psi_remap, res2)

  ; Panel layout
  panelres = True
  panelres@gsnPanelXWhiteSpacePercent = 5
  panelres@gsnPanelYWhiteSpacePercent = 7
  gsn_panel(wks, (/plot1, plot2/), (/2,1/), panelres)

  print("Done. Psi_sigma_Sv (nb,lat) and Psi_remap (z,lat) available in this NCL session.")
end
