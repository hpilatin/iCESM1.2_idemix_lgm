

begin
  f = addfile("pop2_output.nc","r")

  ;---- VARIABLES (edit names if different) ----
  if (isvar(f,"VVEL")) then
    V = f->VVEL
  else if (isvar(f,"VVEL2")) then
    V = f->VVEL2
  else
    print("ERROR: no VVEL found"); exit
  end if

  ; depth (dz) or level depths z_t
  if (isvar(f,"dz")) then
    dz = f->dz      ; 1D [z]
  else if (isvar(f,"dzt")) then
    dz = f->dzt
  else if (isvar(f,"z_t")) then
    zt = f->z_t
    nz = dimsizes(zt)
    dz = new(nz, typeof(zt))
    dz(0:nz-2) = zt(1:nz-1) - zt(0:nz-2)
    dz(nz-1) = dz(nz-2)
  else
    print("ERROR: no dz/z_t found"); exit
  end if

  ; zonal cell width
  if (isvar(f,"e1t")) then
    e1t = f->e1t
  else if (isvar(f,"DXT")) then
    e1t = f->DXT
  else
    print("ERROR: no e1t/DXT found"); exit
  end if

  ; choose first time-step if time present
  if (dimsizes(V) .ge. 4) then
    V1 = V(0,:,:,:)
  else
    V1 = V
  end if

  ; reorder V1 to (z,lat,lon) if necessary — adapt if your file ordering differs
  ; Here we assume V1 is (z,lat,lon)
  Vz = V1

  ; broadcast dz and e1t to 3D arrays
  nz = dimsizes(Vz( :,0,0 ))
  nlat = dimsizes(Vz(0, :,0))
  nlon = dimsizes(Vz(0,0,:))
  dz3 = new((/nz,nlat,nlon/), typeof(dz))
  do k=0,nz-1
    dz3(k,:,:) = dz(k)
  end do
  e13 = new((/nz,nlat,nlon/), typeof(e1t))
  do k=0,nz-1
    e13(k,:,:) = e1t(:,:)
  end do

  ; compute cell transport (m3/s)
  vmo_calc = Vz * e13 * dz3    ; dims (z,lat,lon)

  ; apply basin mask: prefer REGION_MASK if available
  if (isvar(f,"REGION_MASK")) then
    mask = f->REGION_MASK
    ipmask = where((mask.eq.2) .or. (mask.eq.3), 1, 0)  ; Indo-Pacific choice
    ; broadcast mask
    ip3 = new((/nz,nlat,nlon/), typeof(ipmask))
    do k=0,nz-1
      ip3(k,:,:) = ipmask(:,:)
    end do
    vmo_calc = vmo_calc * ip3
  else
    ; optionally use lon bounds (not shown)
  end if

  ; zonal integration (sum over lon)
  vmo_zonal = dim_sum_n_Wrap(vmo_calc, 2)    ; => (z,lat)

  ; vertical cumulative integral from bottom up
  ; determine vertical ordering: assume index 0 = surface -> integrate from bottom (nz-1) to top (0)
  psi = new((/nz,nlat/), typeof(vmo_zonal))
  psi(nz-1,:) = vmo_zonal(nz-1,:)
  do k = nz-2,0,-1
    psi(k,:) = psi(k+1,:) + vmo_zonal(k,:)
  end do

  ; convert to Sv
  psi = psi / 1e6
  psi@units = "Sv"

  print("Psi (Sv) min/max:"); printMinMax(psi)

  ; save or plot psi as needed...
end

