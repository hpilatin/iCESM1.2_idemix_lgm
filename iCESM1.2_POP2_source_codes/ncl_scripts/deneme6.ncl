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

   wks  = gsn_open_wks(wks_type,"deneme6")

   
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

 ; f = addfile(file,"r")

 ; ; -- find VVEL
  ;if (.not. isvar(f,Vname)) then
   ; if (isvar(f,"VVEL2")) then
    ;  Vname = "VVEL2"
   ; else
   ;   print("ERROR: no VVEL or VVEL2 found in file"); exit
  ;  end if
 ; end if
 ; V = f->$Vname$    ; [time,z,lat,lon] or maybe [time,lat,lon,z]
 ; printVarSummary(V)

 ; ; -- dz 1D
 ; if (.not. isvar(f,dzname)) then
  ;  print("ERROR: dz not found. Check variable name"); exit
 ; end if
 ; dz = f->$dzname$    ; [z]

 ; ; -- DXT 2D
 ; if (.not. isvar(f,dxname)) then
  ;  print("ERROR: DXT not found. Check variable name"); exit
 ; end if
 ; DXT = f->$dxname$   ; [lat,lon] (or [j,i] dims)
  
 ; ; -- depth coordinates / z_t (optional)
 ; if (isvar(f,"z_t")) then
 ;   zt = f->z_t
 ; else if (isvar(f,"zt")) then
 ;   zt = f->zt
 ; else
 ;   ; approximate depth at layer midpoints from dz
  ;  nz = dimsizes(dz)
 ;   zt = new(nz, typeof(dz))
  ;  zt(0) = dz(0)/2.0
  ;  do k=1,nz-1
  ;    zt(k) = zt(k-1) + 0.5*(dz(k-1)+dz(k))
  ;  end do
 ; end if


  V = pop_pi_tamas->VVEL(:,:,:,:)      ; [z,lat,lon], take first timestep
  DXT  = pi_grid->DXT(:,:)           ; [lat,lon], zonal grid spacing (m)
  dz  = pi_grid->dz        ; [z], layer thickness (m)

 ;VVEL2 = pop_pi_tamas->VVEL(:,:,:,:)      ; [z,lat,lon], take first timestep



  ;print(max(VVEL))
  ;print(min(VVEL))

  z_t           = pop_pi_tamas->z_t

    lat2d = pop_pi_tamas->TLAT
    lat = lat2d(:,0)
 
  PD   = pop_pi_tamas->PD(:,:,:,:)

 ;; -- lat array (pick common names)
  ;if (isvar(f,"YC")) then
  ;  lat2d = f->YC
  ;  lat = lat2d(:,0)
 ; else if (isvar(f,"TLAT")) then
  ;  lat2d = f->TLAT
  ;  lat = lat2d(:,0)
  ;else if (isvar(f,"lat")) then
  ;  lat = f->lat
  ;else
   ; print("WARNING: latitude not found, using index positions")
   ; lat = ispan(0,dimsizes(DXT(:,0))-1,1)
  ;end if

  ;; -- PD variable detection
 ;pd_found = False
  ;do i=0,dimsizes(PDnames)-1
   ; if (isvar(f,PDnames(i))) then
     ; PDname = PDnames(i)
     ; pd_found = True
    ;  break
   ; end if
  ;end do
  ;if (.not. pd_found) then
   ; print("ERROR: neither PD2 nor PD found."); exit
  ;end if
  ;PD = f->$PDname$   ; [time,z,lat,lon]
 ; print("Using density variable: " + PDname)

  ;---------------- get sizes and time selection ----------------
  Vdims = dimsizes(V)
  ; assume V has time as first dim
  nt = dimsizes(V(:,0,0,0))
  nz = dimsizes(V(0,:,0,0))
  nlat = dimsizes(V(0,0,:,0))
  nlon = dimsizes(V(0,0,0,:))

  ; time indices for last Nyears (assume V timestep = 1 per year; if different, user must adjust)
  ;tstart = max(0, nt - Nyears)
  ;tend   = nt - 1

  ;---------------- Time mean of fields ----------------
  ;V_mean  = dim_avg_n(V(tstart:tend,:,:,:), 0)    ; [z,lat,lon]
  ;PD_mean = dim_avg_n(PD(tstart:tend,:,:,:), 0)   ; [z,lat,lon]

   V_mean  = V(0,:,:,:)
   PD_mean  = PD(0,:,:,:)

  ;---------------- Convert PD to sigma2 (kg/m3 - 1000) ----------------
  ; PD is often in g/cm3 -> convert to kg/m3
  ; For PD2 (potential density referenced to 2000 dbar) this yields sigma2.
  


   sigma2_3d = (PD_mean * 1000.0) - 1000.0   ; [z,lat,lon]

   printVarSummary(sigma2_3d)
   print(max(sigma2_3d))
   print(min(sigma2_3d))

  ; handle missing
  ;---------------- Build 3D dz and DXT arrays ----------------
  ; broadcast dz (1D) to 3D
  dz3 = new((/nz,nlat,nlon/), typeof(dz))
  do k=0,nz-1
    dz3(k,:,:) = dz(k)
  end do

  area3 = new((/nz,nlat,nlon/), typeof(DXT))
  do k=0,nz-1
    area3(k,:,:) = DXT(:,:)
  end do

  ;---------------- Compute cell-wise transport (m^3/s) ----------------
  ; V_mean assumed in m/s
  
  vmo3 = V_mean * dz3 * area3 * 1d-6    ; [z,lat,lon] in m^3/s

  print("vmo3")
  printVarSummary(vmo3)
  print(min(vmo3))
  print(max(vmo3))

  ; We'll work in m^3/s for binning; convert to Sv only for plotting if desired

  ;---------------- Zonal integration per cell for simpler binning ----------------
  ; But for density binning we want to consider each grid cell separately, so keep vmo3.

  ;---------------- Build sigma bin edges & centers ----------------
 
   sigma_min = min(sigma2_3d) - 0.01
   sigma_max = max(sigma2_3d) + 0.01

  nb = nbins

  sigma_edges = sigma_min + (sigma_max - sigma_min) * (ispan(0, nb, 1) / nb)  ; nbins+1 edges
  sigma_centers = 0.5 * (sigma_edges(0:nb-1) + sigma_edges(1:nb))
  
   print(sigma_max)
   print(sigma_min)  
 
   print(max(sigma2_3d))
   print(min(sigma2_3d))

  ;---------------- Initialize arrays: transport histogram per (lat, sigma bin) ----------------
  ; We'll accumulate the zonally-integrated transport into sigma bins for each latitude and vertical level equivalent

  trans_bin = new((/nb,nlat/),float, 0.0)   ; will be m^3/s total per sigma bin and latitude
 
  printVarSummary(trans_bin)
  
  vmo4 = tofloat(vmo3)


  ;---------------- Loop over grid cells and accumulate transport into sigma bins ----------------
  ; We go through all (k,j,i) and add vmo3(k,j,i) to bin where sigma2_3d(k,j,i) falls
  

;b = toint( floor( (s - sigma_min) / ((sigma_max - sigma_min) / nb) ) )


do j = 0, nlat-1
  do i = 0, nlon-1
    do k = 0, nz-1
      s = sigma2_3d(k,j,i)
      if (.not. ismissing(s)) then
        if (s.ge.sigma_min .and. s.le.sigma_max) then
            b = toint( floor( (s - sigma_min) / ((sigma_max - sigma_min) / nb) ) )
          if (b.lt.0) then
            b = 0
          end if
          if (b.ge.nb) then
            b = nb-1
          end if          
          trans_bin(b,j) = trans_bin(b,j) + vmo4(k,j,i)
        end if
      end if
    end do
  end do
end do



  printVarSummary(trans_bin)
  print(min(trans_bin))
  print(max(trans_bin))

  ; Now trans_bin(b,j) is the zonally-summed transport (m3/s) of grid cells at latitude j with sigma in bin b.

  ;---------------- Compute cumulative overturning in density space ----------------
  ; Convention: integrate from dense -> light (sigma large -> small)
  ; Psi_density(lat, sigma_k) = sum_{b >= k} trans_bin(b,lat)  (i.e. cumulative sum from heavy to light)
  
Psi_sigma = new((/nb,nlat/), float, 0.0)
 ; do j=0,nlat-1
 ;   cum = 0.0
 ;   do b = nb-1, 0, -1   ; start from heaviest bin
 ;     cum = cum + trans_bin(b,j)
 ;     Psi_sigma(b,j) = cum
 ;   end do
 ; end do



do j=0,nlat-1
 cum = 0.0
 do b = 0, nb-1
  bb = nb-1 - b        ; reverse index
  cum = cum + trans_bin(bb,j)
  Psi_sigma(bb,j) = cum
 end do
end do

  ; Convert to Sv for plotting
  Psi_sigma_Sv = Psi_sigma / 1e6
  Psi_sigma_Sv@long_name = "MOC streamfunction in density space (Sv)"
  

    printVarSummary(Psi_sigma_Sv)
  print(min(Psi_sigma_Sv))
  print(max(Psi_sigma_Sv))


; sigma_centers array available

  ;---------------- Remap Psi_sigma (lat,sigma) to depth levels ----------------
  ; Approach: use zonal-mean sigma2 profile (sigma2_zonal[z,lat]); for each (z,lat) find sigma at that depth
  ; and then sample Psi_sigma at that sigma (linear interpolation across sigma bins).
  
sigma2_zonal = dim_avg_n(sigma2_3d, 2)   ; [z,lat]

  ; Create output Psi_remap[z,lat]

  Psi_remap = new((/nz,nlat/), float, -9999.0)

  ;rho_zonal = dim_avg_n_Wrap(rho(:,jmask,:),2)   ; (z, lat) array

;--- Loop over latitude
do j = 0, nlat-1
  sigma_profile = sigma2_zonal(:,j)    ; [nz] density profile at latitude j
  psi_vs_sigma  = Psi_sigma(:,j)    ; [nb] cumulative MOC in sigma-space

  ;--- Loop over depth
  do k = 0, nz-1
    s_at_z = sigma_profile(k)

    if (ismissing(s_at_z)) then
      Psi_remap(k,j) = 0.0
    else 
      if (s_at_z .le. sigma_centers(0)) then
      Psi_remap(k,j) = psi_vs_sigma(0)
    else 
     if (s_at_z .ge. sigma_centers(nb-1)) then
      Psi_remap(k,j) = psi_vs_sigma(nb-1)
    else
      idx = ind(sigma_centers .le. s_at_z)
      b = max(idx)   ; last bin where sigma_centers .le. s_at_z 

      frac = (s_at_z - sigma_centers(b)) / (sigma_centers(b+1) - sigma_centers(b) + 1e-12)
      Psi_remap(k,j) = psi_vs_sigma(b) + frac * (psi_vs_sigma(b+1) - psi_vs_sigma(b))
    end if
   end if
  end if
 end do
end do





;     found = False
;      do b = 0, nb-2
;        if (s_at_z .ge. sigma_centers(b) .and. s_at_z .lt. sigma_centers(b+1)) then
;          frac = (s_at_z - sigma_centers(b)) / (sigma_centers(b+1) - sigma_centers(b) + 1e-12)
;          Psi_remap(k,j) = psi_vs_sigma(b) + frac * (psi_vs_sigma(b+1) - psi_vs_sigma(b))
;          found = True
;          exit
;        end if
;      end do
;      print("check")
      
;      if (.not. found) then
;        Psi_remap(k,j) = psi_vs_sigma(0)  ; fallback
;      end if
;    end if
;   end if
;  end if 
; end do
;end do

  
 ; if (s_at_z.le.sigma_centers(0)) then
 ; Psi_remap(k,j) = psi_vs_sigma(0)
;else
;  if (s_at_z.ge.sigma_centers(nb-1)) then
;    Psi_remap(k,j) = psi_vs_sigma(nb-1)
;  else
;    ; find bin b such that sigma_centers(b) <= s_at_z < sigma_centers(b+1)
;    found = False
;    do b = 0, nb-2
;      if (s_at_z.ge.sigma_centers(b) .and. s_at_z.lt.sigma_centers(b+1)) then
;        frac = (s_at_z - sigma_centers(b)) / (sigma_centers(b+1) - sigma_centers(b) + 1e-12)
;        Psi_remap(k,j) = psi_vs_sigma(b) + frac * (psi_vs_sigma(b+1) - psi_vs_sigma(b))
;        found = True
;        exit
;      end if
;    end do
;    if (.not. found) then
;      Psi_remap(k,j) = psi_vs_sigma(0)
;    end if
;  end if
;end if






  ; Psi_remap is in Sv (because psi_vs_sigma was in Sv)
  Psi_remap@long_name = "Psi remapped to depth (Sv)"

  ;---------------- Plotting (two panels): Psi_sigma (density space) and Psi_remap (depth) ----------------
  ;wks = gsn_open_wks("pdf","Psi_density_and_remap")

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
  res2@sfYArray = nz
  res2@tiMainString = "Psi remapped to depth (Sv)"

  plot2 = gsn_csm_contour(wks, Psi_remap, res2)

  ; Panel layout
  panelres = True
  panelres@gsnPanelXWhiteSpacePercent = 5
  panelres@gsnPanelYWhiteSpacePercent = 7

  print("Done")

  gsn_panel(wks, (/plot1, plot2/), (/2,1/), panelres)

  print("Done. Psi_sigma_Sv (nb,lat) and Psi_remap (z,lat) available in this NCL session.")
end
