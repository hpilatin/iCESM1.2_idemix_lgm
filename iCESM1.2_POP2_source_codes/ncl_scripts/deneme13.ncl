

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

   wks  = gsn_open_wks(wks_type,"deneme13")

   
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


  z_t           = pop_pi_tamas->z_t

  vvel = pop_pi_tamas->VVEL(0,:,:,:)          ; [z, y, x], m/s
  dxt  = pi_grid->DXT(:,:)               ; cm
  dz   = pi_grid->dz(:)                  ; cm
  kmt  = pi_kmt->KMT(:,:)               ; ocean depth mask

 ; Convert cm → m
  dxt = dxt * 0.01
  dz  = dz  * 0.01
  vvel = vvel * 0.01


  nlat = dimsizes(vvel(0,:,0))
  nlon = dimsizes(vvel(0,0,:))
  nlev = dimsizes(vvel(:,0,0))

; Create an array for volume transport (m3/s)
  transport = new((/nlev, nlat, nlon/), float)
  transport = 0.0

  ; Calculate transport per grid cell
  do k = 0, nlev-1
    do j = 0, nlat-1
      do i = 0, nlon-1
        if (k .lt. kmt(j,i)) then
          transport(k,j,i) = tofloat(vvel(k,j,i) * dxt(j,i) * dz(k))
        else
          transport(k,j,i) = 0.0
        end if
      end do
    end do
  end do

 ; Sum over depth and longitude for meridional transport (e.g., across latitude)

  
   printVarSummary(transport)
   print(max(transport))
   print(min(transport))

  
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


  z_t = z_t * 0.01

  rho   = pop_pi_tamas->PD(0,:,:,:)

  rho = rho * 1000 ; g/cm3 to kg/m3  

  ;-----------------------------------------------------------------
  ; Define density bins (adjust to your data range)
  ;-----------------------------------------------------------------
  nbins = 200
  rhomin = min(rho)
  rhomax = max(rho)
  drho = (rhomax - rhomin) / nbins
  rho_bins = fspan(rhomin, rhomax, nbins)

 ;-----------------------------------------------------------------
  ; Initialize transport in density bins
  ;-----------------------------------------------------------------
  trans_rho = new((/nbins, nlat/), float)
  trans_rho = 0.0


 ;-----------------------------------------------------------------
  ; Loop through grid and bin by density
  ;-----------------------------------------------------------------
  do j = 0, nlat-1
    do i = 0, nlon-1
      do k = 0, nlev-1
        if (k .lt. kmt(j,i)) then
          rho_val = rho(k,j,i)
          if (.not. ismissing(rho_val)) then
          bin_index = toint((rho_val - rhomin) / drho)
          if (bin_index.ge.0 .and. bin_index.lt.nbins) then
            vol = vvel(k,j,i) * dxt(j,i) * dz(k)    ; m3/s per grid cell   
            trans_rho(bin_index,j) = trans_rho(bin_index,j) + tofloat(vol)
          end if
         end if
        end if
      end do
    end do
  end do


  printVarSummary(trans_rho)
  print(max(trans_rho))
  print(min(trans_rho))
 ; print(trans_rho)

;-----------------------------------------------------------------
  ; Compute cumulative streamfunction in density space
  ;-----------------------------------------------------------------
  ; Integrate from high to low density (deep → surface)
  psi_rho = new((/nbins, nlat/), float)
  psi_rho = 0.0
  psi_rho2 = new((/nbins, nlat/), float)
  psi_rho2 = 0.0

; Normally

   do j = 0, nlat-1
    psi_rho2(:,j) = dim_cumsum(-trans_rho(::-1,j), 0) / 1e6  ; Sv
    psi_rho(:,j) = psi_rho2(::-1,j)
  end do



;  do j = 0, nlat-1
;    psi_rho(:,j) = dim_cumsum(-trans_rho(::-1,j), 0) / 1e6  ; Sv
;  end do

  psi_rho@long_name = "Density streamfunction"
  psi_rho@units = "Sv"



  printVarSummary(psi_rho)
  print(max(psi_rho))
  print(min(psi_rho)) 

  printVarSummary(psi_rho2)
  print(max(psi_rho2))
  print(min(psi_rho2))
 

 ;-----------------------------------------------------------------
  ; Convert density bins to mean depths
  ;-----------------------------------------------------------------
 ; mean_depth = new((/nbins/), float)
 ; mean_depth = 0.0

;nbins1 = dimsizes(rho_bins)
;mean_depth = new((/nbins1/), float)

print(z_t)




print(drho)
print(rho_bins)
print(rhomin)
print(rhomax)


;z3d = conform(rho, z_t, 0)   ; broadcast depth to 3D
;mean_depth = new((/nbins/), float)

;do n = 0, nbins-1
;  if (n.lt.nbins-1) then
;    weights = where(rho.ge.rho_bins(n) .and. rho.lt.(rho_bins(n)+drho), 1, 0)
;  else
;    weights = where(rho.ge.rho_bins(n), 1, 0)
 ; end if

 ; total_weight = sum(weights)   ; sum over all points

 ; if (total_weight.gt.0) then
 ;   mean_depth(n) = sum((z3d*weights)) / total_weight
;  else
 ;   if (n.gt.0) then
 ;     mean_depth(n) = mean_depth(n-1)
  ;  else
 ;     mean_depth(n) = 0.0
;    end if
 ; end if
;end do


printMinMax(rho_bins,0)

printMinMax(rho, 0)

print(rho_bins(0))

print(rho_bins(199))


 miss = -9999.0

  ;==============================
  ;--- Initialize output array ---
  ;==============================
  psi_rho_z = new((/nlev,nlat,nlon/), typeof(psi_rho))
  psi_rho_z = miss   ; fill with missing

  ;==============================
  ;--- Interpolation loop ---
  ;==============================
  do j = 0, nlat-1
    do i = 0, nlon-1
      do k = 0, nlev-1
        rho_val = rho(k,j,i)

        if (.not. ismissing(rho_val)) then
          found = False

          do b = 0, nbins-2
            denom = rho_bins(b+1) - rho_bins(b)
            if (denom.ne.0. .and. .not.ismissing(psi_rho(b,j)) .and. .not.ismissing(psi_rho(b+1,j))) then
              if (rho_val .ge. rho_bins(b) .and. rho_val .le. rho_bins(b+1)) then
                tfrac = (rho_val - rho_bins(b)) / denom
                psi_rho_z(k,j,i) = (psi_rho(b,j)*(1-tfrac)) + (psi_rho(b+1,j)*tfrac)
                found = True
                exit
              end if
            end if
          end do

          ; Extrapolate if outside bins
          if (.not. found) then
            if (rho_val < rho_bins(0) .and. .not. ismissing(psi_rho(0,j))) then
              psi_rho_z(k,j,i) = psi_rho(0,j)
            else
              if (rho_val > rho_bins(nbins-1) .and. .not. ismissing(psi_rho(nbins-1,j))) then
              psi_rho_z(k,j,i) = psi_rho(nbins-1,j)
            else
              psi_rho_z(k,j,i) = miss
            end if
          end if
         end if
        end if
      end do
    end do
  end do

  ;==============================
  ;--- Optional: mask small values ---
  ;==============================
  threshold = 1e-3
  psi_plot = psi_rho_z
  mask1 = abs(psi_rho_z) .lt. threshold
  psi_plot(mask1) = miss

  ;==============================
  ;--- Plotting latitude-depth section ---
  ;==============================
  lat_index = nlat/2  ; example: plot middle latitude
  psi_section = psi_plot(:,lat_index,:)  ; dimensions: depth x lon

 ; wks = gsn_open_wks("png","psi_depth_plot")  ; output PNG
  res = True

  res@gsnLeftString  = "Density Streamfunction (m^3/s)"
  res@gsnRightString = "Depth (m)"
  res@trYReverse     = True     ; depth increasing downward
  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels = (/ -20, -10, -5, 0, 5, 10, 20 /)  ; example contour levels
  res@gsnMaximize    = True

  plot = gsn_csm_contour(wks, psi_section, res)

end

