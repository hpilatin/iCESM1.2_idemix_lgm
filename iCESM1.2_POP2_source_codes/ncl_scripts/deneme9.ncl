

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

   wks  = gsn_open_wks(wks_type,"deneme9")

   
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


;  z_t           = pop_pi_tamas->z_t

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

  
   trans_y = new((/nlat/), float)

   do j = 0, nlat-1
    trans_y(j) = sum(transport(:,j,:))
   end do
  
  ; Units: m3/s

  ; Optional: Convert to Sverdrup (1 Sv = 1e6 m3/s)
  trans_y_Sv = trans_y / 1e6

  printVarSummary(trans_y_Sv)
  print("Meridional volume transport (Sv):")
  print(trans_y_Sv)
  print(max(trans_y_Sv))
  print(min(trans_y_Sv))
end


