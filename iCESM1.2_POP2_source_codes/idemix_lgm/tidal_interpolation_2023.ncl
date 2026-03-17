


;export NCL_POP_REMAP=/mnt/c/Users/pilat/Downloads
;export NCARG_ROOT=/home/pilatin/ncl-6.5.0
;export PATH=$NCARG_ROOT/bin:$PATH
;ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/panel_two_sets.ncl"
;load "/home/hbkheves/ESMF_regrid_template.ncl"


;---Names of various input/output files
 src_file = "/home/hbkheves/diss_dir_js_1_8_rtop_21kyrBP_i6g_-t_-I6.0_7673000_10Dec_remapbil_-180180.nc"
;     src_file = "/home/hbkheves/diss_dir_js_1_8_rtop_21kyrBP_i6g_-t_-I6.0_7673000_10Dec_remapbil_-180180_ctomiss.nc"


;    dst_file = "/home/hbkheves/map_1x1d_to_gx1v6_bilin_da_100716.nc"
    wgt_file = "/home/hbkheves/wilmes_to_POP_gx1v6.nc"
	
    interp_method = "conserve"     ; interpolation method
     

    var_name = "diss"

;---Open files containing source HOMME and destination POP grids
    sfile = addfile(src_file,"r")
;    dfile = addfile(dst_file,"r")
    var_name     = "diss"
;---Get variable to regrid and the source lat/lon grid
    var      = sfile->diss          ; time (1) x ncol
    var_in   = sfile->$var_name$(:,:)

    

;     var_mask = var                                                            ; where origmoc is the one written by POP
;  var_mask = mask(var,(var.gt.0.0.or.var.lt.0.0),True)
;    var_mask = mask(var,(var.gt.0.0),True)


src_lat  = sfile->lat                 ; ncol
src_lon  = sfile->lon                 ; ncol



;  DO_REGRID = True

;  if(DO_REGRID) then


nc_file = "/scratch/usr/hbkheves/211004_LGM_setup/solFe_scenario4_current_gx1v6_AlbaniLGMratio.nc"
nc_file2 = "/scratch/usr/hbkheves/LGM_example_history_file/b.ie12.Bi1850C5BPRDCISO.f19_g16.LGM21ka_frUsp.p875.001.pop.h.1999-12.nc"


    ;---Open files containing source HOMME and destination POP grids
    ncfile = addfile(nc_file,"r")
    ncfile2 = addfile(nc_file2,"r")

present_file = "/home/hbkheves/tidal_energy_gx1v6_20090205.nc"
prfile = addfile(present_file,"r")


dst_lat  = ncfile->TLAT                 ; ncol
dst_lon  = ncfile->TLONG                 ; ncol




dimnames = getvardims(dst_lat)
dims  = dimsizes(dst_lat)
nlat = dims(0)
nlon = dims(1)



smask = sfile->mz
dmask = ncfile->REGION_MASK



print(dimnames)
	
;---Set up regridding options
    Opt                = True
    remap_method = "conserve" 

;---"bilinear" is the default. "patch" and "conserve" are other options.
    Opt@InterpMethod     = interp_method
    Opt@WgtFileName      = wgt_file

    Opt@SrcGridLat       = src_lat
    Opt@SrcGridLon       = src_lon
    Opt@SrcInputFileName = src_file

; The POP grid was obtained from an existing weight file.


;Opt@SrcMask2D  	      = smask
Opt@DstMask2D        = dmask
Opt@DstGridLat       = dst_lat	
Opt@DstGridLon       = dst_lon
Opt@DstGridType    = "curvilinear"
Opt@ForceOverwrite   = True
Opt@Debug            = True
;Opt@PrintTimings     = True
;Opt@lat2d = ncfile->TLAT
;Opt@lon2d = ncfile->TLONG

;Opt@SrcLLCorner = (/ -90.d,   -180.d/)
;Opt@SrcURCorner = (/ 90.d,   179.d/)




   var_regrid = ESMF_regrid(var_in,Opt)     ; Do the regridding
    printVarSummary(var_regrid)

;    dims_regrid = dimsizes(var_regrid)	
;    rank_regrid = dimsizes(dims_regrid)
;    var_regrid!(rank_regrid-1) = "Y"
;    delete(var_regrid@lon2d)
;    delete(var_regrid@lat2d)
    var_regrid!0 = dimnames(0) 
    var_regrid!1 = dimnames(1) 
    var_regrid@coordinates = "TLONG TLAT"
    delete(var_regrid@lon2d)
    delete(var_regrid@lat2d)
;    var_regrid@coordinates = "TLONG TLAT"

    printVarSummary(var_regrid)
  
   var_regrid(199:347,35:42)= 0.0
   var_regrid(0:25,35:42)= 0.0

   

;----Write out regridded file
; https://www.ncl.ucar.edu/Applications/method_2.shtml
    ; output file name

    f_out = "ESMF_regrid_19_script_hlrn_heves07_mask_v2.nc"
    system("/bin/rm -f " + f_out)    ; remove if exists
    fout = addfile(f_out, "c")


;    ; set up netCDF file
    setfileoption(fout, "DefineMode", True)


 
;    Y = ncfile->Y
;    X = ncfile->X
;    TLAT = ncfile->TLAT
;    TLONG = ncfile->TLONG
;    TAREA = ncfile->TAREA



;     var_regrid2(:,0:283) =   var_regrid(:,36:319)
;     var_regrid2(:,284:319) =   var_regrid(:,0:35)

    
;    wave_dissipation = var_regrid(:,:)



    ; global attributes
    fAtt = True
    fAtt@creation_date = systemfunc("date")
    fAtt@weight_file = wgt_file
    fAtt@src_file = src_file
    fAtt@author = "Heves Pilatin"
    fAtt@contact = "hpilatin@marum.de"
    fileattdef(fout, fAtt)


    ; set up dimensions
    dimNames = (/dimnames(0),dimnames(1)/)
    dimSizes = (/nlat, nlon/)
    dimUnlim = (/False, False/)
    filedimdef(fout,dimNames,dimSizes,dimUnlim)




;   var_regrid(199:347,35:42)= 9.969209968386869e+36
;   var_regrid(0:25,35:42)= 9.969209968386869e+36


;   var_regrid(199:347,35:42)= 0.0
;   var_regrid(0:25,35:42)= 0.0

;   wave_dissipation = var_regrid(:, :)

;   print(wave_dissipation(199:347,35:42))

;    wave_dissipation = var_regrid(:, :)
   

     wave_dissipation = new((/384,320/),double)
     assignFillValue(wave_dissipation,var_regrid)


    wave_dissipation!0 = dimnames(0)
    wave_dissipation!1 = dimnames(1)
    wave_dissipation@coordinates = "TLONG TLAT"
    printVarSummary(wave_dissipation)
 

;    printVarSummary(var_regrid)


;  print(wave_dissipation(199:347,35:42))
; print(wave_dissipation(319:347,37:42))

;wave_dissipation@coordinates = "TLONG TLAT"

;wave_dissipation@lat2d = ncfile->TLAT
;wave_dissipation@lon2d = ncfile->TLONG

;wave_dissipation&Y =  ncfile->Y
;wave_dissipation&X=  ncfile->X




    printVarSummary(wave_dissipation)

    ; create variables

   TLAT = ncfile->TLAT
    TLONG = ncfile->TLONG
    TAREA = ncfile->TAREA
     X = ncfile->X
     Y = ncfile->Y




    filevardef(fout,"wave_dissipation",typeof(var_regrid),getvardims(var_regrid))
     filevardef(fout,"TLAT",typeof(TLAT),getvardims(TLAT))
     filevardef(fout,"TLONG",typeof(TLONG),getvardims(TLONG))
    filevardef(fout,"TAREA",typeof(TAREA),getvardims(TAREA))
    filevardef(fout,"X",typeof(X),getvardims(X))
    filevardef(fout,"Y",typeof(Y),getvardims(Y))


   

;   filevardef(fout,"TLAT",typeof(dst_lat),getvardims(dst_lat))
;   filevardef(fout,"TLONG",typeof(dst_lon),getvardims(dst_lon))
;   filevardef(fout,"wave_dissipation",typeof(var_regrid),getvardims(var_regrid))
   



    ; copy attributes from variables
   ; filevarattdef(fout, "wave_dissipation",var_in)
    remapAtt=0
    remapAtt@remap_method = remap_method
   
   filevarattdef(fout, "wave_dissipation",var_regrid)
   filevarattdef(fout,"TLAT",dst_lat)
    filevarattdef(fout,"TLONG",dst_lon)
    filevarattdef(fout,"TAREA",TAREA)
    filevarattdef(fout,"X",X)
    filevarattdef(fout,"Y",Y)


;    filevarattdef(fout,"TLAT",dst_lat)
;    filevarattdef(fout,"TLONG",dst_lon)
;    filevarattdef(fout, "wave_dissipation",var_regrid)
    
 




    setfileoption(fout,"DefineMode",False)



   ; add to netCDF file


   
    fout->wave_dissipation = (/var_regrid/)
    fout->TLAT = (/dst_lat/)
    fout->TLONG = (/dst_lon/)
    fout->TAREA = (/TAREA/)
    fout->X = (/X/)
    fout->Y = (/Y/)


;     delete(fout)

;   end if
 
