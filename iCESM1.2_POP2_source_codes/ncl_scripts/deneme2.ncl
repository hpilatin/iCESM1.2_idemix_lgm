

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


pi_tarea = addfile ("/mnt/d/lgm_data/pi_tarea.nc", "r")
pi_region = addfile ("/mnt/d/lgm_data/region_mask_pi.nc", "r")

pi_kmt = addfile ("/mnt/d/lgm_data/kmt_pi.nc", "r")

wks_type = "pdf"

   wks  = gsn_open_wks(wks_type,"deneme")


   cmap = read_colormap_file("cmocean_delta")
   ncmap = dimsizes(cmap(:,0))


pltres = True
mpres = True
pltres@PanelPlot  = True

plot = new(4,graphic)



;vvel1  = pop_pi_tamas->WVEL(0,:,:,:)


  w_e    = pop_pi_tamas->WVEL		; Eulerian-mean vertical velocity
;  w_i    = pop_pi_tamas->WISOP		; Eddy-induced (bolus) vertical velocity
;  w_sm   = pop_pi_tamas->WSUBM		; Eddy-induced (bolus) vertical velocity
  tarea  = pi_tarea->TAREA             ; [60] x [384] x [320]  => double
  rmask  = pi_region->REGION_MASK       ; [384] x [320] => int
  kmt    = pi_kmt->KMT
  tlat   = pop_pi_tamas->TLAT              ; [384] x [320] => float

  dims   = dimsizes(tarea)
  ny     = dims(0)
  nx     = dims(1)
  km     = max(kmt)

k3d     = conform_dims((/km,ny,nx/),ispan(0,km-1,1),(/0/))
  kmt3d   = conform_dims((/km,ny,nx/),kmt  ,(/1,2/))
  tarea3d = conform_dims((/km,ny,nx/),tarea,(/1,2/))
  ocean   = k3d.le.kmt3d

  WORK1   = tofloat(where(ocean,w_e(0,:,:,:)*tarea3d,0.0))
  ;WORK2   = tofloat(where(ocean,w_i(0,:,:,:)*tarea3d,0.0))
  ;WORK3   = tofloat(where(ocean,w_sm(0,:,:,:)*tarea3d,0.0))
  WMSG    = WORK1@_FillValue

printVarSummary(WORK1)
print(max(WORK1))
print(min(WORK1))


  delete(w_e)
  ;delete(w_i)
  ;delete(w_sm)
  delete(k3d)
  delete(kmt3d)
  delete(tarea3d)
  delete(ocean)

end
