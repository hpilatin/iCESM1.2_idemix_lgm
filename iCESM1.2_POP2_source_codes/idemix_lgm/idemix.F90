!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module idemix

!BOP
! !MODULE: idemix
!
! !DESCRIPTION:
! This module computes the vertical diffusivity as specified in
! the Internal Wave Dissipation, Energy and MIXing scheme first
! described in 
!
!    Olbers, D. and Eden, C., 2013: A Global Model for the Diapycnal
!      Diffusivity Induced by Internal Gravity Waves. J. Phys.
!      Oceanogr., 43, 1759-1779.
!
! Written: Soeren B. Nielsen, University of Copenhagen, December 2015,
! following Prof. Carsten Eden and supervised by prof. Markus Johum.

! !REVISION HISTORY:
! SVN:$Id$

! !USES

   use kinds_mod
   use domain_size
   use domain
   use blocks
   use io
   use io_types
   use constants
   use exit_mod
   use grid
   use communicate
   use global_reductions
   use broadcast
   use tavg
   use time_management
   use tidal_mixing, only: TIDAL_ENERGY_FLUX

   implicit none
   private
   save

! !PUBLIC FUNCTION MEMBERS:

   public :: init_idemix, kvmix_idemix1

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public ::  &
      lidemix,             &! namelist variable; if true, idemix is on
      enable_ide_eke,      &! if true, eke forcing is on
      enable_idemix_hor_diffusion !If true, hor diffusion is on

   real (r8), dimension(:,:,:,:), allocatable, public :: &
      KVIDEMIX,                   &
      E_iw,                       &
      E_iw2,                      &
      iw_diss_plot,               &
      forc_tidal,                 &   
      SIGMA_TOPO_MASK,            &
!      idfdiss_plot,               &
!      idczer_plot,                &
      NMASK
   
!   real (r8), dimension(:,:,:), allocatable, public :: &
!      KVIDEMIX,                   &
!      E_iw

   real (r8), dimension(:,:,:), allocatable, public :: &
      betapl,                     &! Betaplane parameter
      ide_test_niw                ! FAKE FORCING, CONSTANT
  

   real (r8), public ::           &! NOT SURE IF THESE SHOULD EVEN BE HERE 
      idmu0,                      &! namelist variable; dissipation parameter
      idjstar,                    &! namelist variable; GM spectrum bandwidth
      idgamma,                    &! namelist variable; GM spectrum shape
      idtauv,                     &! namelist variable; vertical time scale
      idtauh,                     &! namelist variable; horizontal time scale
      kvidemix_max,               &! namelist variable; maximum kappa from idemix
      kvidemix_min,               &! namelist variable; minimum kappa from idemix
      ide_mix_efficiency,         &! namelist variable; mixing efficiency in idemix
      ide_min_buoy,               &! namelist variable; minimum buoyancy freq
      ide_const_niw,              &! namelist variable; constant surface forcing value
      min_gamma_eg,               &! gamme for EKE
      c_eke_eg                     ! tuning parameter, see Eden Greatbatch 2008 pp 234

   integer (int_kind) ::          &
      i, j, k, iblock,            &
      tavg_E_iw,                  &! id for E_iw
      tavg_iw_diss,             &! id for dissipation
      tavg_tidal,               &
!      tavg_idfdiss,             &
!      tavg_idczer,              &
      tavg_KVIDEMIX

   character (char_len) :: &
      string

  
!EOP

!***********************************************************************

   contains

!***********************************************************************

!BOP
! !IROUTINE: init_idemix
! !INTERFACE:

 subroutine init_idemix

! !DESCRIPTION:
!  Initializes parameters for idemix VERSION 1
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     input namelist variables (for other public namelist variables, see above)
!
!-----------------------------------------------------------------------

!   real (r8) ::                   &! NOT SURE IF THIS IS THE RIGHT PLACE EITHER
!      idmu0,                      &! IDEMIX dissipation parameter 
!      idjstar,                    &! IDEMIX spectral bandwidth
!      idgamma,                    &! IDEMIX GM Spectrum shape constant
!      idtauv,                     &! IDEMIX vertical symmetrisation timescale
!      idtauh,                     &! IDEMIX horizontal symm. time scale
!      kvidemix_max,               &! IDEMIX max kappa
!      ide_mix_efficiency,         &! IDEMIX mixing efficiency
!      ide_min_buoy,               &! namelist variable; minimum buoyancy freq
!      ide_const_niw                ! namelist variable; constant surface forcing value


! FOLLOWING LINE IS WRITTEN TO MAKE IT EASY TO PUT VARIABLES IN NAMELIST

!   namelist /idemix_nml/ lidemix,idmu0,idjstar,idgamma,idtauv,idtauh


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     set defaults for idemix parameters, then read them from namelist
!
!-----------------------------------------------------------------------

   lidemix                      = .true.          ! Enabled if module is used
   enable_ide_eke               = .false.         ! EKE disabled  (Heves)
   enable_idemix_hor_diffusion  = .true.          ! Horizontal diffusion enabled (Heves)
   idmu0                        = 4.0/3.0_r8      ! From McComas and Muller 1981
   idjstar                      = 10.0_r8         ! Olbers and Eden 2013, [3;10]
   idgamma                      = 1.57_r8         ! From O&E 2013 appendix, a=0.1
   idtauv                       = 1.0*86400.0_r8  ! 1 day in seconds
   idtauh                       = 10.0*86400.0_r8 ! 10 days, page 1769
   kvidemix_max                 = 100.0_r8          ! max vertical mixing
   kvidemix_min                 = 0.001_r8         ! min vertical mixing, set to molecular lvl  
   ide_mix_efficiency           = 0.2_r8          ! Idemix mixing efficiency
   ide_min_buoy                 = 1e-12_r8        ! Minimum N^2 for idemix
   ide_const_niw                = 0.0_r8          ! Constant forcing surface
   min_gamma_eg                 = 200_r8
   c_eke_eg                     = 0.0_r8          ! EKE dissipation tuning parameter

!-----------------------------------------------------------------------
!
!     read namelist input and broadcast variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     exit if idemix mixing is not enabled
!
!-----------------------------------------------------------------------

   if (.not. lidemix) return

   allocate( E_iw(nx_block,ny_block,km,nblocks_clinic), &
             KVIDEMIX(nx_block,ny_block,km,nblocks_clinic), &
             ide_test_niw(nx_block,ny_block,nblocks_clinic) ) 
   E_iw = c0
   KVIDEMIX = c0
   ide_test_niw = ide_const_niw/c1000

   allocate (SIGMA_TOPO_MASK(nx_block,ny_block,km,nblocks_clinic))
   allocate (NMASK(nx_block,ny_block,km,nblocks_clinic))
   allocate(betapl(nx_block,ny_block,nblocks_clinic))

 
 

     do iblock=1,nblocks_clinic

       do k=1,km !
         where ( k < KMT(:,:,iblock) ) 
           SIGMA_TOPO_MASK(:,:,k,iblock) = c1
         elsewhere
           SIGMA_TOPO_MASK(:,:,k,iblock) = c0
         endwhere
       enddo ! k

       do k=1,km-1 
         do j=2,ny_block-1
           do i=2,nx_block-1 
             if ( k < KMT(i,j,iblock) ) then
               if ( k == KMT(i-1,j+1,iblock)  .or.  &
                    k == KMT(i  ,j+1,iblock)  .or.  &
                    k == KMT(i+1,j+1,iblock)  .or.  &
                    k == KMT(i-1,j  ,iblock)  .or.  &
                    k == KMT(i+1,j  ,iblock)  .or.  &
                    k == KMT(i-1,j-1,iblock)  .or.  &
                    k == KMT(i  ,j-1,iblock)  .or.  &
                    k == KMT(i+1,j-1,iblock) )      &
                 SIGMA_TOPO_MASK(i,j,k,iblock) = c0  
             endif 
           enddo ! i
         enddo ! j
       enddo ! k



!       do k=1,km !
!         where ( k < KMT(:,:,bid) )
!           SIGMA_TOPO_MASK(:,:,k,bid) = c1
!         elsewhere
!           SIGMA_TOPO_MASK(:,:,k,bid) = c0
!         endwhere
!       enddo ! k

!       do k=1,km-1
!         do j=2,ny_block-1
!           do i=2,nx_block-1
!             if ( k < KMT(i,j,bid) ) then
!               if ( k == KMT(i-1,j+1,bid)  .or.  &
!                    k == KMT(i  ,j+1,bid)  .or.  &
!                    k == KMT(i+1,j+1,bid)  .or.  &
!                    k == KMT(i-1,j  ,bid)  .or.  &
!                    k == KMT(i+1,j  ,bid)  .or.  &
!                    k == KMT(i-1,j-1,bid)  .or.  &
!                    k == KMT(i  ,j-1,bid)  .or.  &
!                    k == KMT(i+1,j-1,bid) )      &
!                 SIGMA_TOPO_MASK(i,j,k,bid) = c0  
!             endif
!           enddo ! i
!         enddo ! j
!       enddo ! k




!       NMASK(:,:,:,bid)=c0

! NMASK is a counter that calculates how many neighboring grid points are masks
! This is because where there are more than two neighboring masks, instabilities
! occur when coupling to biogeochemestry

!       do k=1,km-1
!         do j=2,ny_block-1
!           do i=2,nx_block-1
!             if ( k < KMT(i,j,bid) ) then
!               if ( k >= KMT(i-1,j,bid)) then
!                   NMASK(i,j,k,bid)=NMASK(i,j,k,bid)+c1
!               endif
!               if ( k >= KMT(i+1,j,bid)) then
!                   NMASK(i,j,k,bid)=NMASK(i,j,k,bid)+c1
!               endif
!               if ( k >= KMT(i,j-1,bid)) then
!                   NMASK(i,j,k,bid)=NMASK(i,j,k,bid)+c1
!               endif
!               if ( k >= KMT(i,j+1,bid)) then
!                   NMASK(i,j,k,bid)=NMASK(i,j,k,bid)+c1
!               endif
!             endif
!           enddo ! i
!         enddo ! j
!       enddo ! k
!     do k=1,km
!        where (k>=KMT(:,:,bid))
!           NMASK(:,:,k,bid)=c0
!        endwhere
!     enddo



!     call ugrid_to_tgrid(betapl(:,:,bid),ULAT(:,:,bid),bid)

!     betapl(:,:,bid) = c2*omega*cos(betapl(:,:,bid))/radius







       NMASK(:,:,:,iblock)=c0

! NMASK is a counter that calculates how many neighboring grid points are masks
! This is because where there are more than two neighboring masks, instabilities
! occur when coupling to biogeochemestry
!       do k=1,km-1  
!         do j=2,ny_block-1
!           do i=2,nx_block-1
!             if ( k < KMT(i,j,iblock) ) then
!               if ( k >= KMT(i-1,j,iblock)) then
!                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
!               endif
!               if ( k >= KMT(i+1,j,iblock)) then
!                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
!               endif
!               if ( k >= KMT(i,j-1,iblock)) then
!                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
!               endif
!               if ( k >= KMT(i,j+1,iblock)) then
!                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
!               endif
!             endif
!           enddo ! i
!         enddo ! j
!       enddo ! k
!     do k=1,km 
!        where (k>=KMT(:,:,iblock))
!           NMASK(:,:,k,iblock)=c0
!        endwhere
!     enddo


        do k=1,km-1
         do j=2,ny_block-1
           do i=2,nx_block-1
             if ( k < KMT(i,j,iblock) ) then
               if ( k >= KMT(i-1,j+1,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i,j+1,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i+1,j+1,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i-1,j ,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i+1,j ,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i-1,j-1,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i ,j-1,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
               if ( k >= KMT(i+1,j-1,iblock)) then
                   NMASK(i,j,k,iblock)=NMASK(i,j,k,iblock)+c1
               endif
             endif
          enddo ! i
        enddo ! j
      enddo ! k
     do k=1,km
        where (k>=KMT(:,:,iblock))
           NMASK(:,:,k,iblock)=c0
        endwhere
     enddo





     call ugrid_to_tgrid(betapl(:,:,iblock),ULAT(:,:,iblock),iblock)
      
     betapl(:,:,iblock) = c2*omega*cos(betapl(:,:,iblock))/radius

     enddo ! iblock

   string = 'Energy in IW field'
   call define_tavg_field(tavg_E_iw,'E_iw',3,             &
                          long_name=trim(string),             &
                          units='meters^2/s^2',                 &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   string = 'Dissipation Rate'
   call define_tavg_field(tavg_iw_diss,'iw_diss_plot',3,             &
                          long_name=trim(string),             &
                          units='centimeters^2/s^3',                 &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string = 'KVIDEMIX'
   call define_tavg_field(tavg_KVIDEMIX,'KVIDEMIX',3,             &
                          long_name=trim(string),             &
                          units='centimeters^2/s',                 &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

    string = 'Energy in IW field'
   call define_tavg_field(tavg_E_iw,'E_iw2',3,             &
                          long_name=trim(string),             &
                          units='meters^2/s^2',                 &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )

   string ='Tidal Energy Flux'
   call define_tavg_field(tavg_tidal,'forc_tidal',3,             &
                          long_name=trim(string),             &
                          units='Watts/m^2',                 &
                          grid_loc='3113',                    &
                          coordinates  ='TLONG TLAT z_w_bot time' )   


!   string = 'idczer'
!   call define_tavg_field(tavg_idczer,'idczer_plot',3,             &
!                          long_name=trim(string),             &
!                          units='cm^2/s^3',                 &
!                          grid_loc='3113',                    &
!                          coordinates  ='TLONG TLAT z_w_bot time' )


!   string = 'idfdiss'
!   call define_tavg_field(tavg_idfdiss,'idfdiss_plot',3,             &
!                          long_name=trim(string),             &
!                          units='s/cm^2',                 &
!                          grid_loc='3113',                    &
!                          coordinates  ='TLONG TLAT z_w_bot time' )




!-----------------------------------------------------------------------
!EO

 end subroutine init_idemix

!***********************************************************************


! !IROUTINE: kvmix_idemix1
! !INTERFACE:

 subroutine kvmix_idemix1( WORK3, VISC, E_iw, KVIDEMIX, this_block)

! !DESCRIPTION:
!  Calculate IDEMIX vertical mixing parameters
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(in) :: &
     VISC                ! Working array for RI_LOC
!     WORK3              ! Working array, N^2

   real (r8), dimension(nx_block,ny_block,0:km+1), intent(in) :: &
     WORK3              ! Working array, N^2


!   real (r8), dimension(nx_block,ny_block,km,nblocks_clinic), intent(in) :: &
!      WORK4                ! Buoyany frequency (1/s)



! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nblocks_clinic), intent(out) :: &
      KVIDEMIX          ! kappa from idemix


   type (block), intent(in) :: &
      this_block        ! block information for current block


   real (r8), dimension(nx_block,ny_block,km,nblocks_clinic), intent(inout) :: &
      E_iw

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     input namelist variables (for other public namelist variables, see above)
!
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::            &
      k,                            &! vertical level index
      i,                            &! longitudinal index
      j,                            &! latitudinal index
      kbot,                         &! integer of bottom lvl
      kboundl,                      &! integer of bl lvl
      bid,                          &! local block index
      iblock

   real (r8) ::                     &
      effcor,                       &! 
      gfunc,                        &!
      hfunc,                        &!
      gparm1,                       &!
      gparm2,                       &!
      seawdens,                     &! 
      fxa,                          &!
      intN                           !

   real (r8), dimension(0:km) ::    &
      a_tri,                        &! sub-diagonal
      b_tri,                        &! diagonal
      c_tri,                        &! super diagonal
      d_tri,                        &! 
      delta                          ! vertical operator

   real (r8), dimension(0:km+1) ::  &
      zgrid                          ! used to calculate N^2

   real (r8), dimension(nx_block,ny_block,km) :: &
      idfdiss,                      &! Fdiss/(E^2), Carsten's alpha
      iw_diss,                      &! Storing dissipation
      idczer,                       &! c0, from group velocity
      maxE_iw,                      &! To ensure a positive energy flux
      newE_iw,                      &! this timesteps E_iw
      forc_eke,                     &! Forcing from EKE
      v0,                           &! V0 parameter from eq. 27, maybe not important yet
      kap_idemix,                   &! for calculating kappa
      SIGMA,                        &! Eddy velocity scale
      flux_east,                    &! Flux if horizontal diffusion is enabled
      flux_north                     ! Flux if horizontal diffusion is enabled

   real (r8), dimension(nx_block,ny_block,nblocks_clinic) :: &
      forc_kbot                     ! Bottom forcing from tidal
       

   real (r8), dimension(nx_block,ny_block,km,nblocks_clinic) :: &
      iw_diss_plot,                  &
      forc_tidal

   real (r8), dimension(nx_block,ny_block,km,nblocks_clinic) :: &
       E_iw2

   real (r8), dimension(nx_block,ny_block) :: &
      forc_kbl,                     &! Boundary layer forcing
      cstar,                        &! idemix cstar
      WORK1,                        &! working array for EKE
      C_ROSSBY,                     &! integrated N
      L_ROSSBY                       ! eddy length scale
    


!-----------------------------------------------------------------------
!
!  set constants
!
!-----------------------------------------------------------------------
   gparm1 = 0.9_r8
   gparm2 = 4.3_r8
   seawdens = 1026.0_r8
   


!-----------------------------------------------------------------------
!
!  check if idemix is active
!
!-----------------------------------------------------------------------

   if( lidemix) then

!     write(stdout,*) ' Idemix count 1'

      bid = this_block%local_id

! Calculate forcing from NIWs

    

     forc_kbl = ide_test_niw(:,:,bid) 
     
! Calculate forcing from tides 


   
! do j=1,ny_block
!   do i=1,nx_block
!     do k=1,km ! k=1,km
!       where (KMT(:,:,bid) < k)  ! HEVES
!         TIDAL_ENERGY_FLUX(:,:,bid) = c0
!       endwhere
!     enddo
 !  enddo
 !enddo

!   where (KMT(:,:,bid) < k) ! HEVES
!     TIDAL_ENERGY_FLUX(:,:,bid) = c0
!   endwhere

 

     forc_kbot(:,:,bid) = TIDAL_ENERGY_FLUX(:,:,bid) / (c1000 * 1026.0_r8) ! conversion from g/cm^3 to m^3/s^3 
  
  
    do k=1,42
       where (KMT(:,:,bid) < k) ! HEVES
         forc_kbot(:,:,bid) = c0
       endwhere
    enddo



!     forc_tidal(:,:,bid) = TIDAL_ENERGY_FLUX(:,:,bid) / c1000

!-----------------------------------------------------------------------
!
!  Calculate relevant IDEMIX parameters
!
!-----------------------------------------------------------------------
    
  
     do j=1,ny_block
       do i=1,nx_block

         intN = c0 ! The value of N integrated over the water column
         do k=1,km ! without partial bottom cells  ! HEVES
           if (k < KMT(i,j,bid) ) then ! HEVES 
             intN = intN + ((max(c0, WORK3(i,j,k)))**p5) * dzw(k) * 0.01_r8 ! WORK3 (N^2) is in 1/s, dzw is in cm ! HEVES 
           endif
         enddo ! k

         cstar(i,j) = max(1e-2_r8, (intN/(pi*idjstar))) ! 

         do k=1,km ! HEVES 
           effcor = ((max(c0, WORK3(i,j,k)))**p5)/(eps + abs(FCORT(i,j,bid)))  ! FCORT is 1/s ! HEVES
           gfunc = c2/pi/(c1-(c2/pi)*asin(c1/max(c3,effcor)))*gparm1*max(c3,effcor)**(-c2/c3) * &
                    (c1-exp(-(max(c3,effcor))/gparm2)) ! g(x) is in m/s^2, but here it is converted in cm^2/s

           idczer(i,j,k) = max(c0, idgamma * cstar(i,j) * gfunc) ! idczer or c0 is in m^2/s^3, but here it is already in cm^2/s^3  ! HEVES 

           hfunc =  (c2/pi)/(c1-(c2/pi)*asin(c1/(effcor+eps))) * (effcor-c1)/(effcor+c1) ! h(x) is in m/s^2, but here it is converted in cm^2/s

          ! v0(i,j,k) = max(c0, idgamma * cstar(i,j) * hfunc) ! v0 is in m^2/s^3, but here it is already in cm^2/s^3 ! HEVES 

           idfdiss(i,j,k)  =  max(1d-4,idmu0 * acosh(max(c1,effcor)) * abs(FCORT(i,j,bid))/cstar(i,j)**2) ! idfdiss or alpha_c is in s/m^2, but, here, cstar is in cm/s. Therefore, idfdiss is also in s/cm^2 ! HEVES 
         enddo ! k
       enddo ! i
     enddo ! j
   
     if (enable_idemix_hor_diffusion) then
    ! check for stability criterium, lateral diffusion is explicit
    !  tau_h v0^2 *dt/dx^2 <= 0.5  ->   v0  <  sqrt( 0.5*dx^2/(dt tau_h)  )
       do j=1,ny_block
         do i=1,nx_block
           fxa = c2/c10*min( DXT(i,j,bid), DYT(i,j,bid) )**c2/ max(c1,dtt*idtauh )
           v0(i,j,:) = min( sqrt(fxa), v0(i,j,:) )
         enddo
       enddo
     endif


   
 
  
     do k=1,km ! k=1,km

       where ( k > KMT(:,:,bid))
         idczer(:,:,k) = c0
         idfdiss(:,:,k) = c0 
         v0(:,:,k) = c0
       endwhere
     enddo ! k

     

!  do k=1,km
!    idfdiss_plot(:,:,k,bid) = idfdiss(:,:,k)
!    idczer_plot(:,:,k,bid)  = idczer(:,:,k)
!  enddo



   dzw = dzw * 0.01_r8

!-----------------------------------------------------------------------
!
! Now calculate E in the vertical
!
!-----------------------------------------------------------------------

! TRIDIAGONAL IMPLICIT SOLVER FOR INTERNAL WAVE ENERGY

  
     maxE_iw(:,:,:) = max(c0, E_iw(:,:,:,bid))

     NewE_iw = c0

     do i=1,nx_block
       do j=1,ny_block
         a_tri = c0
         b_tri = c0
         c_tri = c0
         d_tri = c0
         kboundl = 1  !  (HEVES)  ! Set to KBL(i,j) if true boundary layer forcing is added
         kbot = KMT(i,j,bid)
         if (kbot > 0) then ! (kbot > 0) 
           do k = 2,kbot 
             delta(k) = dtt * idtauv/dzw(k-1) * p5 * ( idczer(i,j,k) + idczer(i,j,k-1) )  
           enddo ! k
           delta(1) = c0  
           do k=2,kbot    
             a_tri(k) = - delta(k)*idczer(i,j,k-1)/dzw(k)
           enddo ! k
           a_tri(1)=c0
           do k=2,kbot-1  ! k =2,kbot-1
             b_tri(k) = 1 + delta(k) * idczer(i,j,k)/dzw(k) + delta(k+1) * idczer(i,j,k)/dzw(k) + &
                          dtt * idfdiss(i,j,k)*maxE_iw(i,j,k)

           enddo

           b_tri(1) = 1 + delta(2)/(dzw(1))*idczer(i,j,1) + &
                            dtt * idfdiss(i,j,1)*maxE_iw(i,j,1)

           b_tri(kbot) = 1 + delta(kbot)/dzw(kbot)*idczer(i,j,kbot) + &
                           dtt * idfdiss(i,j,kbot)*maxE_iw(i,j,kbot)
           do k = 2,kbot-1  ! k = 2,kbot-1
             c_tri(k) = - delta(k+1)/dzw(k) * idczer(i,j,k+1)

           enddo
           c_tri(kbot) = c0
           c_tri(1) = -delta(2)/(dzw(1))*idczer(i,j,2) 
           d_tri(1:kbot) = maxE_iw(i,j,1:kbot) + dtt*forc_eke(i,j,1:kbot)
           d_tri(kboundl) = d_tri(kboundl) + dtt*forc_kbl(i,j)/dzw(kboundl)
           d_tri(kbot) = d_tri(kbot) + dtt*forc_kbot(i,j,bid)/dzw(kbot)
           
           call solve_tridiag(a_tri(1:kbot),b_tri(1:kbot),c_tri(1:kbot),d_tri(1:kbot),newE_iw(i,j,1:kbot),kbot)



         endif ! kbot > 0
       enddo ! j
     enddo ! i

 ! Store IW dissipation
 
!   do k=2,km ! k=1,km
!    where ( k < KMT(:,:,bid) ) ! k<= KMT
!      where (REGION_MASK(:,:,bid) < c0)    ! KMT is the deepest point in a grid cell (Heves)
!          E_iw2(:,:,k,bid) = 1d-7
!        elsewhere
          E_iw2(:,:,:,bid) = newE_iw
!         endwhere
!     elsewhere
!           E_iw2(:,:,k,bid) = c0
!     endwhere
!   enddo   

  
    iw_diss(:,:,:) = idfdiss * maxE_iw * newE_iw * 1d4 

 

    dzw = dzw * 1d2
    
   
! When if (enable_idemix_hor_diffusion)
! and add tendency due to lateral diffusion are deactivated, E_iw2 also gives 0
! as the min.


 if (enable_idemix_hor_diffusion) then
 !---------------------------------------------------------------------------------
 ! add tendency due to lateral diffusion
 !---------------------------------------------------------------------------------
 flux_east(:,:,:)=c0
 flux_north(:,:,:)=c0

 
  do k=1,km
   do j=1,ny_block
    do i=1,nx_block-1
      if (k<=KMT(i,j,bid) .and. k<=KMT(i+1,j,bid)) then
        flux_east(i,j,k)=idtauh*p5*(v0(i+1,j,k)+v0(i,j,k)) * &
            (v0(i+1,j,k)*E_iw(i+1,j,k,bid)-v0(i,j,k)*E_iw(i,j,k,bid))/(DXT(i,j,bid))
      else
        flux_east(i,j,k)=c0
      endif
    enddo
   enddo
   do j=1,ny_block
    do i=1,nx_block-1
      if (k<=KMT(i,j,bid) .and. k<=KMT(i,j+1,bid)) then
        flux_north(i,j,k)= idtauh*p5*(v0(i,j+1,k)+v0(i,j,k)) * &

            (v0(i,j+1,k)*E_iw(i,j+1,k,bid)-v0(i,j,k)*E_iw(i,j,k,bid))/DYT(i,j,bid)
      else
        flux_north(i,j,k) = c0
      endif
    enddo
   enddo
  enddo

 do j=2,ny_block
   do i=2,nx_block
    newE_iw(i,j,:)=newE_iw(i,j,:) + dtt* &
                                  (( flux_east(i,j,:) - flux_east(i-1,j,:))/(DXT(i,j,bid)) &
                                 +(flux_north(i,j,:) -flux_north(i,j-1,:))/(DYT(i,j,bid)) )
   enddo
  enddo
 endif


     KVIDEMIX(:,:,:,:) = c0
   
     do k=1,km
       where ( k <= KMT(:,:,bid) )
         where (REGION_MASK(:,:,bid) < c0)    ! KMT is the deepest point in a grid cell (Heves)
           KVIDEMIX(:,:,k,bid) = kvidemix_min ! for the marginal seas(Heves)
         elsewhere
           kap_idemix(:,:,k) = ide_mix_efficiency/(c1 + ide_mix_efficiency) * iw_diss(:,:,k)
           kap_idemix(:,:,k) = max(kvidemix_min, min(kvidemix_max, kap_idemix(:,:,k)/(max(ide_min_buoy, WORK3(:,:,k)))) )
           KVIDEMIX(:,:,k,bid) = min(kvidemix_max, kap_idemix(:,:,k) )
           E_iw(:,:,k,bid) = newE_iw(:,:,k)
         endwhere
       elsewhere
        KVIDEMIX(:,:,k,bid) = c0
        E_iw(:,:,k,bid) = c0
      endwhere 

! CHECK numerically unstable grid cells, i.e. grid cells with more than 2
! neighbor landmask cells
       where (NMASK(:,:,k,bid) > c2)
          KVIDEMIX(:,:,k,bid)=0.1_r8
       endwhere    
       where (WORK3(:,:,k) < c0)
          KVIDEMIX(:,:,k,bid) = kvidemix_min
       endwhere

       iw_diss_plot(:,:,k,bid) = iw_diss(:,:,k)
       forc_tidal(:,:,k,bid) = forc_kbot(:,:,bid)


       call accumulate_tavg_field(KVIDEMIX(:,:,k,bid),tavg_KVIDEMIX,bid,k)

       call accumulate_tavg_field(E_iw(:,:,k,bid),tavg_E_iw,bid,k)

       call accumulate_tavg_field(E_iw2(:,:,k,bid),tavg_E_iw,bid,k)
       
       call accumulate_tavg_field(iw_diss_plot(:,:,k,bid),tavg_iw_diss,bid,k)

       call accumulate_tavg_field(forc_tidal(:,:,k,bid),tavg_tidal,bid,k)

!       call accumulate_tavg_field(idfdiss_plot(:,:,k,bid),tavg_idfdiss,bid,k)

!       call accumulate_tavg_field(idczer_plot(:,:,k,bid),tavg_idczer,bid,k)

     enddo ! k
   
   endif !lidemix




!-----------------------------------------------------------------------
!EOC

 end subroutine kvmix_idemix1

!***********************************************************************



subroutine solve_tridiag(a,b,c,d,f,n)
      implicit none
!---------------------------------------------------------------------------------
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
!---------------------------------------------------------------------------------
        integer,intent(in) :: n
        real*8,dimension(n),intent(in) :: a,b,c,d
        real*8,dimension(n),intent(out) :: f
        real*8,dimension(n) :: q
        real*8 :: p
        integer ii

    q(1)=-c(1)/b(1)
    f(1)=d(1)/b(1)

    do ii = 2,n
      p=c1/(b(ii)+a(ii)*q(ii-1))
      q(ii) = -c(ii)*p
      f(ii) = (d(ii)-a(ii)*f(ii-1))*p
    enddo
 
    do ii = 1,n-1
      f(n-ii) = f(n-ii)+q(n-ii)*f(n-ii+1)
    enddo

 
end subroutine solve_tridiag






end module idemix
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
