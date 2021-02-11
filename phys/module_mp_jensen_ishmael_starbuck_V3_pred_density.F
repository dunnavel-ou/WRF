module module_mp_jensen_ishmael

  !------------------------------------------------------------------------------------------------------!
  ! This is the initial release of the Ice-Spheroids Habit Model with Aspect-ratio Evolution (ISHMAEL).  !
  ! This code calculates the moisture tendencies of water vapor, cloud water, rainwater, and three ice   !
  ! species which are planar-nucleated particles, columnar-nucleated particles, and aggregates. Planar-  !
  ! and columnar-nucleated ice particles are parameterized as spheroids (two growth dimensions are       !
  ! tracked), and ice particle shape and density evolution are predicted for vapor growth, sublimation,  !
  ! riming, and melting. The model can produce various vapor-grown habits (which currently depend on the !
  ! environmental temperature and supersaturation), and the model predicts the change in ice particle    !
  ! shape and density during riming. Ice particle fall speeds are diagnosed based on ice particle        !
  ! properties. Aggregates are currently parameterized as spheroids with a fixed shape and density. Ice  !
  ! is initiated as planar or columnar based on temperature. After nucleation, the shapes of particles   !
  ! are not constrained, i.e., planar-nucleated particles can attain columnar shapes if the growth       !
  ! conditions allow it.                                                                                 !
  !                                                                                                      !
  ! Predicting two ice particle dimensions requires two prognostic variables. Thus mass, number, and two !
  ! size-related mixing ratio are predicted for each ice species.  That makes 12 prognostic variables    !
  ! for ice. Just throw more processors at it. Rain mass and number and cloud mass are also predicted    !
  ! See Jensen et al. 2017 (JAS), Jensen et al. 2018 (MWR), and Jensen et al. 2018 (JAS) for details.    !
  !
  !------------------------------------------------------------------------------------------------------!
  ! WRFV4.2.1 Update: This version of ISHMAEL represents ice aggregates as reduced-density tri-axial     !
  !                 ellipsoids following Dunnavan et al. 2019 and Dunnavan 2020. Aggregate mass and      !
  !                 number are predicted and volume-weighted aggregate density is fixed. Using ellipsoids!
  !                 rather than spheres/fractals therefore provides multiplicative correction factors    !
  !                 for aggregate shapes based off Multi-Angle Snowflake Camera (MASC) data. Currently,  !
  !                 shape and density are fixed but future updates will predict shape and density.       !
  !                                                                                                      !
  !                 The tenative name for this weird aggregation scheme is:                              !
  !                 Stochastic Triple Axis Redistribution By Using Collection Kernels                    !
  !                                                  or                                                  !
  !                                          STARBUCK (Melville 1851)                                    !
  !                                                                                                      !
  !                                                                                                      !
  ! Edwin Lee Dunnavan (ELD), Cooperative Institute for Mesoscale Meteorological Studies (CIMMS),        !
  ! dunnavel@ou.edu                                                                                      !
  ! Date: 05/20/19                                                                                       !
  !                                                                                                      !
  !------------------------------------------------------------------------------------------------------!
  !                                                                                                      !
  !                                                                                                      !
  ! Author:        Anders A. Jensen, NCAR-RAL, ajensen@ucar.edu, 303-497-8484                            ! 
  ! Last modified: 07 Sep 2020 (ISHMAEL version 0.2)                                                     !
  !------------------------------------------------------------------------------------------------------!

  use module_wrf_error   

  implicit none !.. You are welcome

  !.. Constants
  real, private, parameter ::     &
       PI          = 3.14159265,  &  !.. Pi = 3.14159265... until forever
       G_HOME      = 9.8,         &  !.. Gravity (m s^-1)
       RD          = 287.15,      &  !.. Dry air gas constant (J kg^-1 K^-1)
       RV          = 461.5,       &  !.. Moist air gas constant (J kg^-1 K^-1)
       CP          = 1004.0,      &  !.. Heat capacity Air (J kg^-1 K^-1)
       CPW         = 4218.0,      &  !.. Heat capacity Water (J kg^-1 K^-1)
       RCP         = 0.285856574, &  !.. RD / CP
       RHOW        = 1000.0,      &  !.. Water density (kg m^-3)
       T0          = 273.15,      &  !.. STP temperature (K)
       P0          = 101325.0,    &  !.. STP pressure (Pa)
       RHOI        = 920.0,       &  !.. Bulk ice density (kg m^-3)
       R0          = 1.27494,     &  !.. Air density at 1000 mb, T=T0 (kg m^-3)
       NU          = 4.0,         &  !.. Shape parameter for free ice PSDs
       NU_AGG      = 1.0,         &  !.. Shape parameter for aggregate PSD
       ao          = 0.1e-6,      &  !.. min a-axis scaling factor
       zet         = 1.88,        &  !.. Particle area power-law exponent from Mitchell 1996
       sig_heyms   = 0.2285,      &  !.. Particle area power-law coefficient from Mitchell 1996
       kc          = 1.0,         &  !.. Fallspeed parameter from Dunnavan 2020
       a0          = 1.7e-3,      &  !.. Mitchell 2005 fall speed drag parameter
       b0          = 0.8,         &  !.. Mitchell 2005 fall speed drag parameter
       Co          = 0.6,         &  !.. Mitchell 2005 fall speed drag parameter
       delta_0     = 5.83,        &  !.. Mitchell 2005 fall speed drag parameter
       phinuc      = 1.0,         &
       Ec          = 1.0,         &
       betmon      = 0.5,         &
       bcp         = 16.98,       & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       ccp         = 1.023,       & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       dcp         = 3981.,       & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       fcp         = 4.e-6,       & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       E0p         = 3.e4,        & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       T1p         = -15.,        & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       gamp        = 1.3,         & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       kamp        = 1.5,         & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       chipa       = 1.,          & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       chipb       = 3.,          & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       omegp       = 0.05,        & ! Phillips et al. 2015 Snow Aggregation Sticking Efficiency
       b1p         = 0.05,        & ! Phillips et al. 20015 Snow Aggregation breakup efficiency
       b2p         = 0.05,        & ! Phillips et al. 20015 Snow Aggregation breakup efficiency
       b3p         = 0.05

!.. Define the ice species
!.. Do not change these without changing/knowing how nucleation works
  integer, parameter :: CAT  = 3 !.. Number of ice categories (species)
  integer, parameter :: IICE1 = 1 !.. Ice initiated as planar
  integer, parameter :: IICE2 = 2 !.. Ice initiated as columnar
  integer, parameter :: IICE3 = 3 !.. Aggregates

 !.. AGG Stuff
  real, private :: alph_ar_agg, bet_ar_agg, alphv_agg,chi_ba,chi_ca
  real, private :: abatog, pbmin, pcmin, fs_agg, alph_agg,Sa0_agg, Sa1_agg
  real, private :: sigma, pbeta_area, abamax, abamin,pbeta_chi
  real, private :: pba, pcb, effc, rhobar_ellip, phiagg
  real, private :: betx_agg, betm_sph, betm_agg, abet, bbet

  real, private :: C1, C2

  real, private :: aba,bba,bcb,bet_ba,bet_ca,phiba_phica_bar,pb2

  real, dimension(3), private :: ki, Cnk, pbeta_colx, pbeta_coly


  integer, private, parameter ::  &
      agg_switch = 1,             & ! 0 = Spheres/fractals, 1 = Ellipsoid dist., 2 = Oblate spheroid?
      dvt_switch = 1                ! 0 = Wisner Approx., 1 = Seifert Approx.


  !.. The lookup table arrays are itab for ice-cloud water collection 
  !.. and itabr for ice-rain collection. The lookup tables are read in 
  !.. the subroutine jensen_ishmael_init which occurs once.
  real, private, allocatable, dimension(:,:,:,:,:) :: itab, itabr

  !.. Inherent growth ratio data is stored in the igrdata array
  !.. The gamma_tab is a tabulation of the gamma function needed
  !.. because ice particle shape evolves. 
  real, private, allocatable, dimension(:) :: igrdata, gamma_tab

  !.. Aggregation lookup tables arrays
  real, private, allocatable, dimension(:,:,:) :: coltab, coltabn

  !.. Added for double precision
  integer, parameter :: LUT_KIND_R4 = selected_real_kind(6) !.. 4 byte real
  real(kind=LUT_KIND_R4), private, allocatable, dimension(:,:,:,:,:) :: itab_o, itabr_o
  real(kind=LUT_KIND_R4), private, allocatable, dimension(:) :: gamma_tab_o


contains
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. This subroutine currently reads in the ice-cloud collection,
  !.. ice-rain collecion, gamma calculations, and inherent growth ratio 
  !.. data lookup tables and calculates the aggregate and ice-ice aggregation parameters.

  subroutine jensen_ishmael_init
    implicit none

    integer :: i1, i2, i3, i4, kcol



  !.. Read in ice-cloud water collection lookup table
    if(.not.allocated(itab_o)) then 
       allocate(itab_o(51,51,51,11,2))

       open(20, FILE='ishmael-qi-qc.bin', FORM='unformatted')
       do i1 = 1,51
          do i2 = 1,51
             do i3 = 1,51
                do i4 = 1,11
                   read(20) itab_o(i1,i2,i3,i4,1), &
                        itab_o(i1,i2,i3,i4,2)
                enddo
             enddo
          enddo
       enddo
       close(20)
       print*, 'Jensen_ISHMAEL lookup table 1 of 3'
    endif

  !.. Read in ice-rain collection lookup table
    if(.not.allocated(itabr_o)) then
       allocate(itabr_o(51,51,51,11,6))

       open(30, FILE='ishmael-qi-qr.bin', FORM='unformatted')
       do i1 = 1,51
          do i2 = 1,51
             do i3 = 1,51
                do i4 = 1,11
                   read(30) itabr_o(i1,i2,i3,i4,1), &
                        itabr_o(i1,i2,i3,i4,2), &
                        itabr_o(i1,i2,i3,i4,3), &
                        itabr_o(i1,i2,i3,i4,4), &
                        itabr_o(i1,i2,i3,i4,5), &
                        itabr_o(i1,i2,i3,i4,6)
                enddo
             enddo
          enddo
       enddo
       close(30)
       print*, 'Jensen_ISHMAEL lookup table 2 of 3'
    endif

  !.. Read in gamma function lookup table for speedy code
    if(.not.allocated(gamma_tab_o)) then
       allocate(gamma_tab_o(505001))

       open(40, FILE='ishmael-gamma-tab.bin', FORM='unformatted')
       do i1 = 1,505001
          read(40) gamma_tab_o(i1)
       enddo
       close(40)
       print*, 'Jensen_ISHMAEL lookup table 3 of 3'
    endif

  !.. Read in inhernent growth ratio data (-1C to -60C)
  !.. This currently assumes columns below -20C (plates are used later in the code)
  !.. From Chen and Lamb 1994
    if(.not.allocated(igrdata)) then
       allocate(igrdata(60))

       igrdata = &
            (/ 0.910547, 0.81807, 0.6874, 0.60127, 1.59767, 2.32423, 2.08818,  &
            1.61921, 1.15865, 0.863071, 0.617586, 0.453917, 0.351975, 0.28794, &
            0.269298, 0.28794, 0.333623, 0.418883, 0.56992, 0.796458, 1.14325, &
            1.64103, 1.90138, 1.82653, 1.61921, 1.47436, 1.32463, 1.25556,     &
            1.22239, 1.206, 1.11522, 1.10751, 1.10738, 1.11484, 1.12234,       &
            1.12221, 1.14529, 1.16884, 1.20104, 1.22573, 1.25094, 1.27666,     &
            1.31183, 1.3388, 1.35704, 1.37553, 1.38479, 1.39411, 1.40349,      &
            1.41294, 1.42245, 1.43202, 1.44166, 1.45137, 1.46114, 1.47097,     &
            1.48087, 1.50105, 1.50087, 1.51098 /)
    endif

  !.. Check that tables are allocated (Thanks to JYH and Shaofeng)
    if(.not.allocated(itab)) allocate(itab(51,51,51,11,2))
    if(.not.allocated(itabr)) allocate(itabr(51,51,51,11,6))
    if(.not.allocated(gamma_tab)) allocate(gamma_tab(505001))

    itab = itab_o; itabr = itabr_o; gamma_tab = gamma_tab_o

   !.. Here, pre-calculate all the cool aggregate stuff. ELD.
   ! Here are some assumed aggregate property constants
   alphv_agg = (4./3.)*PI
   betm_sph = 2.0
   rhobar_ellip = 15.
   C2 = (delta_0**2)/4.
   C1 = 1/(C2*sqrt(Co))


   ! These are the mean aspect ratios where: pba = mean(b/a), and pcb = mean(c/b)
   ! abatog is a dispersion factor where higher values represent more concentrated
   ! ellipsoid shapes. sigma is the orientation std. (cf. Ryzhkov et al. 2011)
   pba = 0.616
   pcb = 0.6795
   !abatog = 0.001
   abatog = 0.2922
    !pba = 0.95
    !pcb = 0.6
   !abatog = 1.0

   sigma = 20.

   pbmin = pba/(1.-pba)
   pcmin = pba*pcb/(1.-pcb)
   abamin = max(pbmin,pcmin)
   abamax = 20.

   aba = abatog*abamax+(1.-abatog)
   bba = ((1./pba)-1.)*aba
   bcb = ((1./pcb)-1.)*aba/pba

   abet = 1.0
   bbet = 3280.*sigma**(-2)-1.055

    ! Put in bet_ba and bet_ca calculations. sigma.
   ! if (agg_switch.eq.1) then
   ! fs_agg = 0.9379*beta_moms(0.3574,0.1636,aba,bba,bcb)
   ! elseif (agg_switch.eq.2) then
   ! phiagg = pcb
   ! fs_agg = sqrt(1.-phiagg**2.)/asin(sqrt(1-phiagg**2.))
   ! elseif (agg_switch.eq.0) then
   ! fs_agg = 1.
   ! endif


    if (agg_switch.eq.1) then
    fs_agg = 0.9379*beta_moms(0.3574,0.1636,aba,bba,bcb)
    Sa0_agg = 0.4977*beta_moms(1.386,0.09691,aba,bba,bcb)+&
             0.4965*beta_moms(0.,1.,aba,bba,bcb)
    Sa1_agg =0.4977*beta_moms(1.386+1.,0.09691+1.,aba,bba,bcb)+&
            0.4965*beta_moms(1.,2.,aba,bba,bcb)
    elseif (agg_switch.eq.2) then
    phiagg = pcb
    fs_agg = sqrt(1.-phiagg**2.)/asin(sqrt(1-phiagg**2.))
    elseif (agg_switch.eq.0) then
    fs_agg = 1.
    Sa0_agg = 1.
    Sa1_agg = 1.
    endif

    alph_ar_agg = 2.**(zet)*10.**(2.*zet-4.)/PI*sig_heyms
    bet_ar_agg = zet-2.

    ! Calculate all the beta distribution calculations
    ! Average ellipsoid aspect ratio
    phiba_phica_bar = beta_moms(1.,1.,aba,bba,bcb)

    ! second shape moment
    pb2 = beta_moms(2.,2.,aba,bba,bcb)

    ! Need to change
    ! betm_agg = 1.016*betm_sph-0.1955 ! ellipsoidal beta_m fit which yields same Zagg/qagg**2 as spheres (Dunnavan 2020)
    ! NEW POWER LAW FIT THAT ENSURES THAT WE USE THE SAME TWO MASS MOMENTS BETWEEN SPHERES AND ELLIPSOIDS
    betm_agg = 0.8591*(log((gamma(1.+2.*betm_sph)/(gamma(1.+betm_sph)**2))*&
    (phiba_phica_bar**2/pb2))+0.533)

    ! Ellipsoid projected area power-law approx.
    bet_ba = 1.0
    bet_ca = 0.0065*sigma

    ! Fall speed power-law parameters
    betx_agg = betm_agg - kc*bet_ar_agg
    !alphx_agg = 8.*alphm_agg*alph_ar_agg**(-kc)*G_HOME*rhodum/(PI*mu**2)

    chi_ba = 1.+kc*bet_ba
    chi_ca = 1.+kc*bet_ca

    ! Ellipsoid Bulk Area calculation
    pbeta_area = beta_moms(bet_ba,bet_ca,aba,bba,bcb)

    ! Best number shape calculation
    pbeta_chi = beta_moms(chi_ba,chi_ca,aba,bba,bcb)

    ! k index
    ki = (/0., 1., 2./)

    ! Binomial coefficients (2 choose k)
    Cnk = (/1., 2., 1./)


    ! Aggregate Collection
    do kcol=1,3 ! collection kernel expansion
    pbeta_colx(kcol) = beta_moms(bet_ba*ki(kcol)/2.,bet_ca*ki(kcol)/2.,aba,bba,bcb)
    pbeta_coly(kcol) = beta_moms(bet_ba*(1.-ki(kcol)/2.),bet_ca*(1.-ki(kcol)/2.),aba,bba,bcb)

    enddo



  end subroutine jensen_ishmael_init

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. The main subroutine.

  subroutine  MP_JENSEN_ISHMAEL(ITIMESTEP, DT_IN, P, DZ, &  
       TH, QV, QC, QR, NR, QI1, NI1, AI1, CI1,           &
       QI2, NI2, AI2, CI2,                               &
       QI3, NI3, AI3, CI3,                               &
       IDS,IDE, JDS,JDE, KDS,KDE,                        &
       IMS,IME, JMS,JME, KMS,KME,                        &
       ITS,ITE, JTS,JTE, KTS,KTE,                        &
       RAINNC, RAINNCV, SNOWNC, SNOWNCV,                 &
       diag_effc3d, diag_effi3d, diag_dbz3d,             &
       diag_vmi3d_1, diag_di3d_1, diag_rhopo3d_1, diag_phii3d_1, &
       diag_vmi3d_2, diag_di3d_2, diag_rhopo3d_2, diag_phii3d_2, &
       diag_vmi3d_3, diag_di3d_3, diag_rhopo3d_3, diag_phii3d_3, &
       diag_itype_1,diag_itype_2,diag_itype_3       &
       )

  !.. ITIMESTEP - time step counter (integer)
  !.. DT_IN - model time step (s)
  !.. P  - air pressure (Pa)
  !.. DZ - difference in height over interface (m)
  !.. TH - potential temperature (K)
  !.. QV - water vapor mixing ratio (kg kg^-1)
  !.. QC - cloud water mixing ratio (kg kg^-1)
  !.. QR - rain water mixing ratio (kg kg^-1)
  !.. NR - rain number concentration (kg^-1)
  !.. QI1 - planar-nucleated mixing ratio (kg kg^-1)
  !.. NI1 - planar-nucleated number mixing ratio (kg^-1)
  !.. AI1 - planar-nucleated a^2c mixing ratio (m^3 kg^-1)
  !.. CI1 - planar-nucleated c^2a mixing ratio (m^3 kg^-1)
  !.. QI2 - columnar-nucleated mixing ratio (kg kg^-1)
  !.. NI2 - columnar-nucleated number mixing ratio (kg^-1)
  !.. AI2 - columnar-nucleated a^2c mixing ratio (m^3 kg^-1)
  !.. CI2 - columnar-nucleated c^2a mixing ratio (m^3 kg^-1)
  !.. QI3 - aggregate mixing ratio (kg kg^-1)
  !.. NI3 - aggregate number mixing ratio (kg^-1)
  !.. NMON3 - aggregate monomer mixing ratio (kg^-1)
  !.. AI3 - aggregate a^2c mixing ratio (m^3 kg^-1)
  !.. CI3 - aggregate c^2a mixing ratio (m^3 kg^-1)
  !.. IDS,IDE - I domain start/end locations
  !.. JDS,JDE - J domain start/end locations 
  !.. KDS,KDE - K domain start/end locations
  !.. IMS,IME - I memory start/end locations
  !.. JMS,JME - J memory start/end locations 
  !.. KMS,KME - K memory start/end locations
  !.. ITS,ITE - I tile start/end locations
  !.. JTS,JTE - J tile start/end locations
  !.. KTS,KTE - K tile start/end locations
  !.. RAINNC - total accumulated precipitation (mm)
  !.. RAINNCV - total accumulated precipitation in a time step (mm)
  !.. SNOWNC - total ice accumulated precipitation (mm)
  !.. SNOWNCV - total ice accumulated precipitation in a time step (mm)
  !.. DIAG_EFFC3D - Effective cloud water radius for coupling to radiation (m)
  !.. DIAG_EFFI3D - Effective ice radius for coupling to radiation (m)
  !.. DIAG_DBZ3D - 10-cm reflectivity (dBZ)
  !.. DIAG_VMI3D_1 - planar-nucleated mass-weighted fall speeds (m s^-1)
  !.. DIAG_VMI3D_2 - columnar-nucleated mass-weighted fall speeds (m s^-1)
  !.. DIAG_VMI3D_3 - aggregate mass-weighted fall speeds (m s^-1)
  !.. DIAG_DI3D_1 - planar-nucleated mass-weighted maximum diameter (m)
  !.. DIAG_DI3D_2 - columnar-nucleated mass-weighted maximum diameter (m)
  !.. DIAG_DI3D_3 - aggregate mass-weighted maximum diameter (m)
  !.. DIAG_RHOPO3D_1 - planar-nucleated mass-weighted effective density (kg m^-3)
  !.. DIAG_RHOPO3D_2 - columnar-nucleated mass-weighted effective density (kg m^-3)
  !.. DIAG_RHOPO3D_3 - aggregate mass-weighted effective density (kg m^-3)
  !.. DIAG_PHII3D_1 - planar-nucleated number-weighted aspect ratio
  !.. DIAG_PHII3D_2 - columnar-nucleated number-weighted aspect ratio
  !.. DIAG_PHII3D_3 - aggregate number-weighted aspect ratio
  !.. DIAG_NMON3D - aggregate number-weighted # of monomers
  !.. DIAG_ITYPE_1 - currently not used
  !.. DIAG_ITYPE_2 - currently not used
  !.. DIAG_ITYPE_3 - currently not used

    implicit none

  !.. Grid variables
    integer, intent(in)    ::           &
         ids, ide, jds, jde, kds, kde , &
         ims, ime, jms, jme, kms, kme , &
         its, ite, jts, jte, kts, kte

  !.. Time and timestep
    integer, intent(in) :: ITIMESTEP
    real, intent(in) :: dt_in

  !.. 3D variables
    real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: p, dz

    real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) ::     &
         th, qv, qc, qr, nr, qi1, ni1, ai1, ci1, qi2, ni2, ai2, ci2, &
         qi3, ni3, ai3, ci3,                                         &
         diag_effc3d, diag_effi3d, diag_dbz3d,                       &
         diag_vmi3d_1, diag_di3d_1, diag_rhopo3d_1, diag_phii3d_1,   &
         diag_vmi3d_2, diag_di3d_2, diag_rhopo3d_2, diag_phii3d_2,   &
         diag_vmi3d_3, diag_di3d_3, diag_rhopo3d_3, diag_phii3d_3,   &
         diag_itype_1, diag_itype_2, diag_itype_3


    real, dimension(ims:ime, kms:kme, jms:jme) :: diag_nmon3d

  !.. 2D variables
    real, dimension(ims:ime, jms:jme), intent(inout) ::              &
         rainnc, rainncv, snownc, snowncv

    integer :: I,K,J


  !.. Use tile locations for our variables
  !.. Make everthing 1D for our loop
  !.. These variables go in and out of newmicro
    real, dimension(kts:kte) ::                    &  
         th1d, qv1d, qc1d, qr1d, nr1d, p1d, dz1d,  &
         effc1d, effi1d, dbz1d, nmon1d

    real, dimension(cat, kts:kte) ::               &  
         qi1d, ni1d, ai1d, ci1d,                   &
         vmi1d, di1d, rhopo1d, phii1d, icetype1d

  !.. Precipitation
    real :: qrpre1d
    real, dimension(cat) :: qipre1d

    ! Aggregates output for fallspeed variables
!    REAL :: sigvt1D(KTE), bm_agg1D(KTE)
!    REAL :: zetba1(KTE), zetca1(KTE)

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Accumulated precipitation in one timestep (zero out before microphysics)
    qrpre1d      =0.
    qipre1d(IICE1)=0.
    qipre1d(IICE2)=0.
    qipre1d(IICE3)=0.
    
  !.. Loop over the grid and send each vertical column into the microphysics routine
    do j=jts,jte             !.. j loop (north-south)         
       do i=its,ite          !.. i loop (horizontal)

  !.. Get the 1D variables from WRF and zero out process rates
          do k=kts,kte       !.. k loop (vertical)

  !.. Thermo,rain,cloud,vapor,etc
             TH1D(k)         =  TH(i,k,j)
             QV1D(k)         =  QV(i,k,j)
             QC1D(k)         =  QC(i,k,j)
             QR1D(k)         =  QR(i,k,j)
             NR1D(k)         =  NR(i,k,j)
             P1D(k)          =   P(i,k,j)
             DZ1D(k)         =  DZ(i,k,j)
             EFFC1D(k)       =  0.
             EFFI1D(k)       =  0.
             DBZ1D(k)        =  -35.
  !.. Planar-nucleated
             QI1D(IICE1,k)    = QI1(i,k,j)
             NI1D(IICE1,k)    = NI1(i,k,j)
             AI1D(IICE1,k)    = AI1(i,k,j)
             CI1D(IICE1,k)    = CI1(i,k,j)
             VMI1D(IICE1,k)   = 0.
             DI1D(IICE1,k)    = 0.
             RHOPO1D(IICE1,k) = 0.
             PHII1D(IICE1,k)  = 0.
             ICETYPE1D(IICE1,k) =  0.
  !.. Columnar-nucleated
             QI1D(IICE2,k)    = QI2(i,k,j)
             NI1D(IICE2,k)    = NI2(i,k,j)
             AI1D(IICE2,k)    = AI2(i,k,j)
             CI1D(IICE2,k)    = CI2(i,k,j)
             VMI1D(IICE2,k)   = 0.
             DI1D(IICE2,k)    = 0.
             RHOPO1D(IICE2,k) = 0.
             PHII1D(IICE2,k)  = 0.
             ICETYPE1D(IICE2,k) =  0.
  !.. Aggregates
             QI1D(IICE3,k)    = QI3(i,k,j)
             NI1D(IICE3,k)    = NI3(i,k,j)
             AI1D(IICE3,k)    = AI3(i,k,j)
             CI1D(IICE3,k)    = 0.
             NMON1D(k)        = CI3(i,k,j)
             VMI1D(IICE3,k)   = 0.
             DI1D(IICE3,k)    = 0.
             RHOPO1D(IICE3,k) = 0.
             PHII1D(IICE3,k)  = 0.
             ICETYPE1D(IICE3,k) =  0.
          enddo !.. (vertical k-loop)

  !.. Call me Ishmael.
          call me_ishmael(itimestep, kts, kte, i, j, dt_in, p1d, dz1d,   &
               th1d, qv1d, qc1d, qr1d, nr1d, qi1d, ni1d, nmon1d, ai1d, ci1d,&
               effc1d, effi1d, dbz1d, icetype1d,  &
               vmi1d, di1d, rhopo1d, phii1d, qrpre1d, qipre1d)

  !.. Send the 1D columns from newmicro back to WRF
          do k=kts,kte !.. k loop (vertical)
             TH(i,k,j)     =         TH1D(k)
             QV(i,k,j)     =         QV1D(k)       
             QC(i,k,j)     =         QC1D(k)
             QR(i,k,j)     =         QR1D(k)
             NR(i,k,j)     =         NR1D(k) 
             DIAG_EFFC3D(i,k,j)  =   max((min(EFFC1D(k),50.e-6)),2.51e-6)
             DIAG_EFFI3D(i,k,j)  =   max((min(EFFI1D(k),999.e-6)),5.e-6)
             DIAG_DBZ3D(i,k,j)   =   max(dbz1d(k),-35.)
  !.. Planar-nucleated
             QI1(i,k,j)    =    QI1D(IICE1,k)
             NI1(i,k,j)    =    NI1D(IICE1,k)
             AI1(i,k,j)    =    AI1D(IICE1,k)
             CI1(i,k,j)    =    CI1D(IICE1,k)
             DIAG_VMI3D_1(i,k,j)  = vmi1d(IICE1,k)
             DIAG_DI3D_1(i,k,j)   = di1d(IICE1,k)
             DIAG_RHOPO3D_1(i,k,j)= rhopo1d(IICE1,k)
             DIAG_PHII3D_1(i,k,j) = phii1d(IICE1,k)
             DIAG_ITYPE_1(i,k,j)  = icetype1d(IICE1,k)
  !.. Columnar-nucleated
             QI2(i,k,j)    =    QI1D(IICE2,k)
             NI2(i,k,j)    =    NI1D(IICE2,k)
             AI2(i,k,j)    =    AI1D(IICE2,k)
             CI2(i,k,j)    =    CI1D(IICE2,k)
             DIAG_VMI3D_2(i,k,j)  = vmi1d(IICE2,k)
             DIAG_DI3D_2(i,k,j)   = di1d(IICE2,k)
             DIAG_RHOPO3D_2(i,k,j)= rhopo1d(IICE2,k)
             DIAG_PHII3D_2(i,k,j) = phii1d(IICE2,k)
             DIAG_ITYPE_2(i,k,j)  = icetype1d(IICE2,k)
  !.. Aggregates
             QI3(i,k,j)    =    QI1D(IICE3,k)
             NI3(i,k,j)    =    NI1D(IICE3,k)
             !NMON3(i,k,j)  =    NMON1D(k)
             AI3(i,k,j)    =    AI1D(IICE3,k)
             !CI3(i,k,j)    =    CI1D(IICE3,k)
             CI3(i,k,j)    = NMON1D(k) ! replaced with monomer number
             DIAG_VMI3D_3(i,k,j)  = vmi1d(IICE3,k)
             DIAG_DI3D_3(i,k,j)   = di1d(IICE3,k)
             DIAG_RHOPO3D_3(i,k,j)= rhopo1d(IICE3,k)
             !DIAG_PHII3D_3(i,k,j) = phii1d(IICE3,k)
             if (QI3(i,k,j).gt.1e-8.and.NI3(i,k,j).gt.1.25e-7) then
             DIAG_PHII3D_3(i,k,j) = max(NMON1D(k),1.25e-7)/max(NI1D(IICE3,k),1.25e-7) ! number-weighted # of monomers for aggs
             DIAG_PHII3D_3(i,k,j) = max(DIAG_PHII3D_3(i,k,j),2.)
             else
             DIAG_PHII3D_3(i,k,j) = 2.0
             endif

             !DIAG_NMON3D(i,k,j) = NMON1D(k)/NI1D(IICE3,k) ! number-weighted # of monomers for aggs
             !DIAG_NMON3D(i,k,j) = 0.
             DIAG_ITYPE_3(i,k,j)  = icetype1d(IICE3,k)
          enddo

  !.. Accumulated precipitation and precipitation rates
  !.. Frozen precipitation
          SNOWNC(i,j) = SNOWNC(i,j) + qipre1d(IICE1) + qipre1d(IICE2) + qipre1d(IICE3)
          SNOWNCV(i,j) = qipre1d(IICE1) + qipre1d(IICE2) + qipre1d(IICE3)

  !.. Total precipitation         
          RAINNC(i,j) = RAINNC(i,j) + qrpre1d + qipre1d(IICE1) + qipre1d(IICE2) + qipre1d(IICE3)
          RAINNCV(i,j) = qrpre1d + qipre1d(IICE1) + qipre1d(IICE2) + qipre1d(IICE3)

       enddo ! i loop (horizontal)
    enddo    ! j loop (north-south)         

  end subroutine MP_JENSEN_ISHMAEL

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. The microphysics
  subroutine me_ishmael(it, kts, kte, i, j, dt, pres_e, dzmic,          &
       theta, qv, qc, qr, nr, qi, ni, nmon, ai, ci,                     &
       effc1d, effi1d, dbz1d, icetype1d, vmi1d, di1d, rhopo1d, phii1d,  &
       qrpre1d, qipre1d)

    implicit none
    
  !.. User logicals (all should be TRUE except for testing)
    logical, parameter :: AGGREGATION = .TRUE. ! Aggregation
    logical, parameter :: RIMING      = .TRUE. ! Riming
    logical, parameter :: RAIN_ICE    = .TRUE. ! Rain-ice collection (for testing)
    logical, parameter :: FREEZE_QC   = .TRUE. ! Freeze cloud drops (Homogeneous)
    logical, parameter :: SPLINTERS   = .TRUE. ! Rime Splinters
    logical, parameter :: WET_GROW    = .TRUE. ! Wet growth

  !.. Input/Output variables from WRF
    integer, intent(in) :: it, i, j  !.. Integer timestep, model i and j
    integer, intent(in) :: kts, kte  !.. Tile k-start and k-end
    !integer, intent(in) :: CAT, IICE1, IICE2, IICE3 !.. Number of categories, planar and columnar index
    real, intent(in) :: dt !.. Timestep


    real, dimension(kts:kte), intent(in) :: pres_e, dzmic !.. Pressure and dz
    real, dimension(kts:kte), intent(inout) :: theta !.. Theta (K) 
    real, dimension(kts:kte), intent(inout) :: qv    !.. Vapor mixing ratio (kg kg^-1)
    real, dimension(kts:kte), intent(inout) :: qc    !.. Cloud water mixing ratio (kg kg^-1)
    real, dimension(kts:kte), intent(inout) :: qr    !.. Rain mixing ratio (kg kg^-1)
    real, dimension(kts:kte), intent(inout) :: nr    !.. Rain number mixing ratio (# kg^-1)
    real, dimension(kts:kte), intent(inout) :: nmon  !.. Aggregate monomer number mixing ratio (# kg^-1)

    real, dimension(cat, kts:kte), intent(inout) :: qi !.. Ice mass mixing ratio (kg kg^-1)
    real, dimension(cat, kts:kte), intent(inout) :: ni !.. Ice number mixing ratio (# kg^-1)
    real, dimension(cat, kts:kte), intent(inout) :: ai !.. Ice a^2c mixing ratio (m^3 kg^-1)
    real, dimension(cat, kts:kte), intent(inout) :: ci !.. Ice c^2a mixing ratio (m^3 kg^-1)

  !.. Ice diagnostics
    real, dimension(kts:kte), intent(inout) :: effc1d, effi1d, dbz1d
    real, dimension(cat, kts:kte), intent(inout) :: rhopo1d 
    real, dimension(cat, kts:kte), intent(inout) :: phii1d
    real, dimension(cat, kts:kte), intent(inout) :: vmi1d  
    real, dimension(cat, kts:kte), intent(inout) :: di1d 
    real, dimension(cat, kts:kte), intent(inout) :: icetype1d
    real :: dbzr,dbzsum

    real :: qrpre1d                 !.. Rain precipitation accumulation (mm)
    real, dimension(cat) :: qipre1d !.. Ice liquid-equivalent precipitation accumulation (mm)

  !.. The Parameters                                            
!    real, parameter :: NU = 4.0          !.. Ice distribution parameter (keep at 4)

    real, parameter :: QSMALL = 1.e-12   !.. Smallest ice mass (kg kg^-1)
    real, parameter :: QNSMALL= 1.25e-7  !.. Smallest ice number (# kg^-1)
    real, parameter :: QASMALL= 1.e-24   !.. Smallest ice volume (m^-3 kg^-1)
 !   real, parameter :: ao = 0.1e-6       !.. alphstr=ao^(1-deltastr)
    real, parameter :: AR = 149.1        !.. 'a' parameter in rain fallspeed-size relation (m^(1-br)/s)
    real, parameter :: BR = 0.5          !.. 'b' parameter in rain fallspeed-size relation   
    real, parameter :: F1R = 0.78        !.. Ventilation coefficient for rain    
    real, parameter :: F2R = 0.308       !.. Ventilation coefficient for rain  
    real, parameter :: LAMMAXR = 1./20.E-6   !.. Max allowed value for rain slope parameter (m^-1)     
    real, parameter :: LAMMINR = 1./2800.E-6 !.. Min allowed value for rain slope parameter (m^-1)

  !.. Local/Dummy variables
    logical :: vgflag  !.. Vapor growth flag to determine if vapor growth
                       !.. is producing subsaturated conditions

    logical, dimension(cat) :: has_ice    !.. true if ice exists for an ice species
    logical, dimension(cat) :: dry_growth !.. Wet growth/dry growth riming flag
    integer :: k, cc, ccvar     !.. k is the k-loop index, cc is the ice species index
    integer :: gi         !.. Index for the Gamma(NU) lookup
    integer :: gi2, gi3   !.. Index for the Gamma(NU) lookup
    integer :: numice     !.. Number of ice species currently with ice
    integer :: current_index  !.. Used to determine number of ice species wi th mass > QSMALL
    real, allocatable, dimension(:) :: icearray !.. Array with size of number
                                                !.. of ice species with mass > QSMALL
    real :: wetg          !.. wet growth check
    real :: gamma_arg     !.. Gamma argument: Gamma(gamma_arg)
    real, dimension(cat) :: ani, cni, rni !.. Characteristic a-,c-,r-axis sizes
    real, dimension(cat) :: aniold, cniold, dsold
    real :: alphstr !.. alpha_* = ao^(1-delta_*)
    real :: alphv   !.. alphv   = 4/3*PI*alpha_*
    real :: betam   !.. betam   = 2 + delta_*
    real :: anf, cnf, rnf, iwcf     !.. a-,c-,r-axis and iwc after vapor growth
    real :: anfr, cnfr, rnfr, iwcfr !.. a-,c-,r-axis and iwc after riming
    real :: phibr, phifr            !.. aspect ratio before and after riming
    real :: capgam  !.. distribution-weighted capacitance (m)
    real :: i_dt    !.. inverse timestep 
    real :: i_cp    !.. inverse heat capacity of air
    real :: sink, source, ratio !.. Terms used to prevent over-depletion
    real :: rrr                 !.. rain characteristic radius

  !.. Fall speed parameters
    real :: phiivt, bl, al, aa, ba, qe, xn, bx, xm, f_c1, f_c2, bm, am, nre
                
  !.. Local variables for riming
    real :: vir, vif, vfr, rhobarrime, rimedr, rimedrr, rimec1, rhobarwg

  !.. Local variables for ice nucleation, multiplication
    real :: anuc, nibnuc(CAT), amass, cmass, fmult, nucfrac
    real :: a_demott, b_demott, c_demott, d_demott
    real :: curnum, ratel, ratekg, inrate

  !.. Local variables for vapor growth
    real :: vi, iwci, nim3dum, prdsum, qcrimesum, qcrimesumr

  !.. Local variables for melting
    real :: mfact, nimelt

  !.. Local variables for evaporation
    real :: epsr, dumqv, dumqvs, arn
    
  !.. Some dummy variables
    real :: dum, dum1, dumt, dumgam, dumgam2
    real :: tmpsum, tmpsumr, tmpmelt

  !.. Local variables for Radiation
    real :: lrsig, sig, lwc, ncm3dum, r_n, phirad, radw1, radw2, radw3, totaliwc

    real, dimension(cat) :: alph_area, bet_area, alph_m, alph_v, bet_m, betac, alphc
    real, dimension(cat) :: Viadum, Vicdum

  !.. Local variables for Aggregation
    real :: an1, an2, an3, phiagg1, phiagg2, phiagg3, alphm_agg, alphx_agg, alphvt_agg
    real :: vtn, Nmonb

    real, dimension(6) :: dvtn, dvtm, dvtc


  !.. Fall speeds (some not used) and effectice ice radius variables
    real, dimension(kts:kte) :: vtrm, vtrn, vtrmc, vtrnc, nc
    real, dimension(cat, kts:kte) :: vtn1,vtrni1, effi, vtrmi1, vtrzi1, betvt1
    real, dimension(cat, kts:kte) :: zetba1, zetca1


  !.. Lookup table values 
    real ::    rrri, rrni, rqci, rdsi, rrho, proc(cat,2), procr(cat,6)
    integer :: irri, irni, iqci, idsi, irho, iti
    real :: qi_qc_nrm(cat), qi_qc_nrd(cat)
    real :: qi_qr_nrm(cat), qi_qr_nrd(cat), qi_qr_nrn(cat)

  !.. Sedimentation variables
    integer :: nn, nstep
    real, dimension(kts:kte) :: fluxqr, fluxnr, fluxnmon
    real, dimension(cat, kts:kte) :: fluxqi, fluxni, fluxai, fluxci
    real, dimension(cat) :: falltndqi, falltndni, falltndai, falltndci
    real :: falltndqr, falltndnr, falltndmon, maxfall
    logical :: sedi
    logical :: domicro

  !.. The rest of the variables
    real, dimension(kts:kte) :: rhoair, i_rhoair !.. Air density (kg m^-3) and inverse air density
    real :: qcrimefrac(cat)   !.. Qc_Rime/(Qc_Rime + Qr_Rime)
    real :: gammnu     !.. Gamma(NU)
    real :: i_gammnu   !.. 1/Gamma(NU)
    real :: fourthirdspi !.. 4/3 * PI
    real :: mnuccd  !.. Deposition nucleation rate (kg kg^-1 s^-1)
    real :: nnuccd  !.. Deposition nucleation rate (# kg^-1)                               
    real :: pcc     !.. Cloud mass growth rate (kg kg^-1 s^-1)                                            
    real :: qagg(cat)    !.. Aggregation mass tendency (kg kg^-1)
    real :: nagg(cat)    !.. Aggregation number tendency (# kg^-1)
    real :: vagg(cat)    !.. Aggregation volume tendency (m^3 kg^-1)
    real :: vcagg(cat)   !.. Aggregation volume tendency (m^3 kg^-1)
    real :: prd(cat)     !.. Vapor growth rate (kg kg^-1 s^-1)                       
    real :: nrd(cat)     !.. Vapor growth number rate sublimation (# kg^-1 s^-1)                       
    real :: ard(cat)     !.. A-axis vapor growth rate mixing ratio (m kg^-1 s^-1) 
    real :: crd(cat)     !.. C-axis vapor growth rate mixing ratio (m kg^-1 s^-1) 
    real :: prdr(cat)    !.. Riming growth rate (kg kg^-1 s^-1)                               
    real :: ardr(cat)    !.. A-axis riming growth rate mixing ratio (m kg^-1 s^-1)
    real :: crdr(cat)    !.. C-axis riming growth rate mixing ratio (m kg^-1 s^-1)
    real :: qmlt(cat)    !.. Mass melting rate (kg kg^-1 s^-1)
    real :: nmlt(cat)    !.. Number melting rate (# kg^-1 s^-1)
    real :: amlt(cat)    !.. A-axis melting rate (m kg^-1 s^-1)
    real :: cmlt(cat)    !.. C-axis melting rate (m kg^-1 s^-1)
    real :: rimesum(cat)  !.. Ice-cloud rime mass rate (kg m^-3 s^-1)
    real :: rimesumr(cat)  !.. Ice-rain rime mass rate (kg m^-3 s^-1)
    real :: rimetotal(cat) !.. Total rime mass rate (kg m^-3 s^-1)
    real :: dQRfzri(cat)   !.. Change in rain mass from ice-rain (kg kg^-1 s^-1)
    real :: dQIfzri(cat)   !.. Change in ice mass from ice-rain (kg kg^-1 s^-1)
    real :: dNfzri(cat)    !.. Change in ice number from ice-rain (# kg^-1 s^-1)
    real :: dQImltri(cat)  !.. Change in ice mass from ice-rain melt (kg kg^-1 s^-1)
    real :: dNmltri(cat)   !.. Change in ice number from ice-rain melt (# kg^-1 s^-1)
    real :: numrateri(cat) !.. Number rate from ice-rain (value from lookkup table)
    real :: rainrateri(cat) !.. Rain mass rate from ice-rain (value from lookkup table)
    real :: icerateri(cat)  !.. Ice mass rate from ice-rain (value from lookkup table)
    real :: totalnuc(cat)   !.. Total mass nucleated to ice (kg kg^-1 s^-1)
    real :: totalnucn(cat)  !.. Total number nucleated to ice (# kg^-1 s^-1)
    real :: nim   !.. Number of homogeneously frozen cloud drops (# kg^-1 s^-1)
    real :: mim   !.. Mass of homogeneously frozen cloud drops (kg kg^-1 s^-1)
    real :: nimr    !.. Number of homogeneously frozen rain drops (# kg^-1 s^-1)
    real :: mimr    !.. Mass of homogeneously frozen rain drops (kg kg^-1 s^-1)
    real :: nmon3 ! Monomer transfer rate from aggregation (kg kg^-1 s^-1)
    real :: deltastr(cat)  !.. Delta_* (Average inherent growth ratio)
    real :: rhobar(cat)    !.. Ice effective density (kg m^-3)
    real :: rhodepout(cat) !.. Ice effective density (kg m^-3)
    real :: dsdepout(cat) !.. Ice effective density (kg m^-3)
    real :: rhorimeout(cat) !.. Ice effective density (kg m^-3)
    real :: gdentotal     !.. Total rime density from ice-cloud and ice-rain (kg m^-3)
    real :: gdenavg(cat)  !.. Rime density from ice-cloud (kg m^-3)
    real :: gdenavgr(cat) !.. Rime density from ice-rain (kg m^-3)
    real :: fv(cat)       !.. Ventilation coefficient for diffusion
    real :: fh(cat)       !.. Ventilation coefficient for heat
    real :: temp, i_temp  !.. Air temperature (K) and inverse
    real :: qvs     !.. Liquid saturation vapor mixing ratio (kg kg^-1)
    real :: qvi     !.. Ice saturation vapor mixing ratio (kg kg^-1) 
    real :: sup     !.. Liquid supersaturation ratio 
    real :: sui     !.. Ice supersaturation ratio      
    real :: svpi    !.. Ice saturation vapor pressure
    real :: svpl    !.. Liquid saturation vapor pressure
    real :: qs0     !.. Saturation vapor mixing ratio at 0C (kg kg^-1) 
    real :: xxls    !.. Latent heat of sublimation (J kg^-1)
    real :: xxlv    !.. Latent heat of vaporization (J kg^-1) 
    real :: xxlf    !.. Latent heat of fusion (J kg^-1)
    real :: ab      !.. Latent heat correction parameter for growth   
    real :: abi     !.. Latent heat correction parameter for growth (ice)  
    real :: mu      !.. Dynamic visocity of air (kg m^-1 s^-1)
    real :: dv      !.. Diffusivity of water vapor in air (m^2 s^-1)
    real :: kt      !.. Thermal conductivity (W m^-1 K^-1)
    real :: nsch    !.. Schmidt number          
    real :: npr     !.. Prandlt number              
    real :: lamr    !.. Rain slope parameter n(D)=n0rr*exp(-lamr*D)dD  
    real :: n0rr    !.. Rain distribution parameter         
    real :: prc     !.. Autoconversion cloud water to rain (kg kg^-1 s^-1)
    real :: nprc    !.. Autoconversion cloud number loss (# kg^-1 s^-1)
    real :: nprc1   !.. Autoconversion rain number gain (# kg^-1 s^-1)
    real :: pre     !.. Rain evpaoration rate (kg kg^-1 s^-1)    
    real :: npre    !.. Rain evaporation number  (# kg^-1 s^-1)    
    real :: pra     !.. Rain accretion of cloud (kg kg^-1 s^-1)
    real :: npra    !.. Rain accretion of cloud number (# kg^-1 s^-1)
    real :: nragg   !.. Rain self-collection (# kg^-1 s^-1)         
    real :: igr     !.. Inherent growth ratio
    real :: qmult(cat)  !.. Ice multiplication rate (kg kg^-1 s^-1)
    real :: nmult(cat)  !.. Ice number multiplication rate (# kg^-1 s^-1)
    real :: mbiggr !.. Biggs freezing (kg kg^-1 s^-1)
    real :: nbiggr !.. Biggs number freezing (# kg^-1 s^-1)


  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. NU=4 is constant so gammnu=6
    gammnu   = gamma(NU)
    i_gammnu = 1./gammnu !.. Inverse gamma
    i_dt     = 1./dt     !.. Inverse timestep
    i_cp     = 1./CP     !.. Inverse heat capacity
    fourthirdspi = 4./3.*PI !.. Four divided by three times pi

  !.. No need to loop through microphysics and sedimentation unless
  !.. condensate exists
    sedi = .false.

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
    do k=kts,kte !.. Grid K loop             

  !.. Check the column for condensate or supersaturation wrt liquid.
  !.. If it is there, then do the microphysics
       domicro= .false.
       temp   = theta(k)/(100000./pres_e(k))**(RCP)
       svpl   = polysvp(temp,0)
       qvs    = 0.622*svpl/(pres_e(k)-svpl)
       sup    = qv(k)/qvs-1.

       svpi  = polysvp(temp,1)
       qvi   = 0.622*svpi/(pres_e(k)-svpi)
       if(temp.gt.T0) qvi=qvs
       sui   = qv(k)/qvi-1.

       if(qc(k).gt.QSMALL.or.qr(k).gt.QSMALL.or.qi(IICE1,k).gt.QSMALL.or.&
            qi(IICE2,k).gt.QSMALL.or.qi(IICE3,k).gt.QSMALL.or.sup.ge.0.) then
          domicro= .true.
       endif
       
  !..  If condensate or supersaturation wrt liquid, do microphysics
       if(domicro) then
          
          if(allocated(icearray)) deallocate(icearray)
          
          numice = 0     !.. Number of species with ice
          do cc = 1, cat !.. Loop over all ice species
             has_ice(cc) = .false.
             if(qi(cc,k).gt.QSMALL) then !.. Check for ice
                
  !.. Get the number of species with ice
                has_ice(cc) = .true.
                numice = numice + 1
             endif
          enddo

  !.. Determine the number of ice species to loop over
  !.. For vapor growth/riming, etc
          if(numice.gt.0.and..not.allocated(icearray)) then 
             allocate(icearray(numice))
             current_index = 1
             do cc = 1, cat
                if(has_ice(cc)) then
                   icearray(current_index) = cc
                   current_index = current_index + 1
                endif
             enddo
          endif

  !.. Initialize variables for all ice species and process rates                               
          do cc = 1, cat !.. Ice variables
             dry_growth(cc)=.true.
             qcrimefrac(cc)=0.
             rimetotal(cc)=0.
             nmult(cc)=0.
             qmult(cc)=0.
             ani(cc)=0.
             cni(cc)=0.
             rni(cc)=0.
             qagg(cc)=0.
             nagg(cc)=0.
             vagg(cc)= 0.
             vcagg(cc)=0.
             prd(cc) =0.
             nrd(cc) =0.
             ard(cc) =0.
             crd(cc) =0.
             prdr(cc)=0.
             ardr(cc)=0.
             crdr(cc)=0.
             qmlt(cc)=0.
             nmlt(cc)=0.
             amlt(cc)=0.
             cmlt(cc)=0.
             deltastr(cc)=1.
             dsold(cc)=1.
             dsdepout(cc)=1.
             rhobar(cc)=RHOI
             rhodepout(cc)=RHOI
             rhorimeout(cc)=RHOI
             if(cc.eq.IICE3) then
             rhobar(cc)=rhobar_ellip
             endif
             fv(cc)=1.
             fh(cc)=1.
             rimesum(cc)=0.
             rimesumr(cc)=0.
             numrateri(cc)=0.
             rainrateri(cc)=0.
             icerateri(cc)=0.
             dQRfzri(cc)=0.
             dQIfzri(cc)=0.
             dNfzri(cc)=0.
             dQImltri(cc)=0.
             dNmltri(cc)=0.
             gdenavg(cc)=400.
             gdenavgr(cc)=400.
             totalnuc(cc)=0.
             totalnucn(cc)=0.
             nibnuc(cc)=0.
             aniold(cc)=0.
             cniold(cc)=0.
             proc(cc,1)=0.
             proc(cc,2)=0.
             qi_qc_nrm(cc)=0.
             qi_qc_nrd(cc)=0.
             qi_qr_nrm(cc)=0.
             qi_qr_nrn(cc)=0.
             qi_qr_nrd(cc)=0.
             procr(cc,1)=0.
             procr(cc,2)=0.
             procr(cc,3)=0.
             procr(cc,4)=0.
             procr(cc,5)=0.
             procr(cc,6)=0.
             fluxqi(cc,k)=0.
             fluxni(cc,k)=0.
             fluxai(cc,k)=0.
             fluxci(cc,k)=0.
             falltndqi(cc)=0.
             falltndni(cc)=0.
             falltndai(cc)=0.
             falltndci(cc)=0.
             di1d(cc,k)=0.
             vmi1d(cc,k)=0.
             vtn1(cc,k)=0.
             vtrmi1(cc,k)=0.
             vtrni1(cc,k)=0.
             vtrzi1(cc,k)=0.
             betvt1(cc,k)=1.
             zetba1(cc,k)=0.
             zetca1(cc,k)=0.
             rhopo1d(cc,k)=0.
             phii1d(cc,k)=0.
             qipre1d(cc)=0.
             effi(cc,k)=0.
             icetype1d(cc,k)=0.
          enddo
          dbz1d(k)=-35.
          qrpre1d=0.
          fluxqr(k)=0.
          fluxnr(k)=0.
          fluxnmon(k)=0.
          falltndqr=0.; falltndnr=0.
          falltndmon=0.
          vtrm(k)=0.; vtrn(k)=0.
          effc1d(k) = 0.
          effi1d(k) = 0.
          nprc =0.; npre=0.; npra=0.; pra=0.; prc=0.
          nprc1=0.; nragg=0.
          nim =0.; mim =0.
          nimr=0.; mimr=0.
          mbiggr=0.; nbiggr=0.
          mnuccd=0.; nnuccd=0.; pcc=0.; fmult=0.
          gdentotal=400.
          rimedr=0.; rimedrr=0.
          dbzr=0.; dbzsum=0.
          nmon3 = 0.
          Nmonb = 2.


            ! Output for aggregate fallspeed stuff
         ! alphvt_agg1D(k) = 0.
         ! betvt_agg1D(k) = 1.0
         ! zetba1(k) = 1.0
         ! zetca1(k) = 1.0

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Temperature and temperature-dependent variables
          i_temp= 1./temp
          rhoair(k)= pres_e(k)/(RD*temp)
          i_rhoair(k) = 1./rhoair(k)
          svpi  = polysvp(temp,1)
          qvi   = 0.622*svpi/(pres_e(k)-svpi)
          if(temp.gt.T0) qvi=qvs
          qs0   = 0.622*polysvp(T0,0)/(pres_e(k)-polysvp(T0,0))
          sui   = qv(k)/qvi-1.
          xxls  = 3.15e6-2370.*temp+0.3337e6
          xxlv  = 3.1484E6-2370.*temp         
          xxlf  = xxls - xxlv
          ab    = 1.+(xxlv*qvs/(RV*temp**2))*xxlv*i_cp 
          abi   = 1.+(xxls*qvi/(RV*temp**2))*xxls*i_cp
          mu    = 1.496E-6*temp**1.5/(temp+120.)            
          arn   = ar*(R0*i_rhoair(k))**0.5
          dv    = 8.794E-5*temp**1.81/pres_e(k)
          kt    = 2.3823e-2 +7.1177e-5*(temp-T0)
          nsch  = mu*i_rhoair(k)/dv
          npr   = mu*i_rhoair(k)/kt

  !.. Inherent growth ratio data from Lamb and Scott 1972 and Chen and Lamb 1994
  !.. Assumes planar particles at temperatures below -20C (see Bailey and Hallett)
  !.. Function of temperature only (thanks to JHY and GHB for moving this fuction to where it belongs)        
          igr = get_igr(igrdata, temp)
          if((temp-T0).lt.-20.) then
             igr=0.7
          endif
                
  !.. Assume constant cloud drop number concentration (200 cm^-3)
          nc(k) = 200.e6*i_rhoair(k)
       
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Ice Microphysics

          if(numice.gt.0) then     !.. If one or more ice species have ice
             tmpsum = 0.
             tmpsumr= 0.
             do ccvar = 1, numice  
                cc = icearray(ccvar)
             
  !.. Variable check
                ni(cc,k)=max(ni(cc,k),QNSMALL)
                ai(cc,k)=max(ai(cc,k),QASMALL)
                ci(cc,k)=max(ci(cc,k),QASMALL)
                if (cc.ne.IICE3) then
                ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
                cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)

                if(ai(cc,k).lt.1.e-12.or.ci(cc,k).lt.1.e-12) then
                ai(cc,k)=min(ai(cc,k),ci(cc,k))
                   ci(cc,k)=ai(cc,k)

                   ani(cc)=(ai(cc,k)/ni(cc,k))**0.3333333333
                   cni(cc)=(ci(cc,k)/ni(cc,k))**0.3333333333
                endif

                if(ani(cc).lt.2.e-6.or.cni(cc).lt.2.e-6) then
                   ani(cc)=2.e-6
                   cni(cc)=2.e-6
                endif

                ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)
                ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)

                deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))

  !.. Check to keep ice properties in bounds
                call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                     ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                     alphstr, alphv, betam)

                dsold(cc) = deltastr(cc)

                elseif(cc.eq.IICE3) then
                 call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))
                alphm_agg = (qi(cc,k)/ni(cc,k))*(1./phiba_phica_bar)*ani(cc)**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)

                endif

  !.. Set incoming volume limit after advection
  !.. This is the smallest bulk volume for which we care about shape              
  !.. If smaller than this limit, assume spherical

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Ice-cloud riming rates from lookup table (currently excluding riming of aggregates)
                if(RIMING.and.qc(k).gt.1.e-7) then
                
  !.. The lookup table is in log space.  Do not change the following 4 lines without
  !.. changing the offline-built lookup table
                   rrni = 13.498*log10(0.5e6*rni(cc))
                   rqci = 8.776*log10(1.0e7*(exp(qc(k))-1.))
                   rdsi = 50.*(deltastr(cc) - 0.5)
                   rrho = 7.888*log10(0.02*rhobar(cc))
                   irni = int(rrni)
                   iqci = int(rqci)
                   idsi = int(rdsi)
                   irho = int(rrho)
                
  !.. Small limit (stay in the bounds of the lookup table)
                   rrni = max(rrni,1.)
                   rqci = max(rqci,1.)
                   rdsi = max(rdsi,1.)
                   rrho = max(rrho,1.)
                   
                   irni = max(irni,1)
                   iqci = max(iqci,1)
                   idsi = max(idsi,1)
                   irho = max(irho,1)
                
  !.. Large limit (stay in the bounds of the lookup table)
                   rrni = min(rrni,real(size(itab,1))-1.)
                   rqci = min(rqci,real(size(itab,2))-1.)
                   rdsi = min(rdsi,real(size(itab,3))-1.)
                   rrho = min(rrho,real(size(itab,4))-1.)
                   
                   irni = min(irni,size(itab,1)-1)
                   iqci = min(iqci,size(itab,2)-1)
                   idsi = min(idsi,size(itab,3)-1)
                   irho = min(irho,size(itab,4)-1)
                   
                   do iti = 1, 2
                      call access_lookup_table(itab,irni,iqci,idsi,irho,iti, &
                           rdsi,rrho,rqci,rrni,proc(cc,iti))
                   enddo
                   
  !.. proc(cc,1) is the normalized riming rate
  !.. proc(cc,2) is related to the normalized rime density
  !.. rimesum is the riming rate in kg m^-3 s^-1
                   rimesum(cc)  = max((proc(cc,1)*ni(cc,k)*nc(k)*rhoair(k)**2),0.)
                   qi_qc_nrm(cc)= proc(cc,1)
                   qi_qc_nrd(cc)= proc(cc,2)

  !.. Limit riming rate if too small
                   if((rimesum(cc)*i_rhoair(k)).lt.QSMALL.or.cc.eq.IICE3) then
                      rimesum(cc)  =0.
                      qi_qc_nrm(cc)=0.
                      qi_qc_nrd(cc)=0.
                   endif
                endif   !.. End ice-cloud riming
             
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Ice-rain riming rates from lookup table
                if(RAIN_ICE.and.qr(k).gt.QSMALL) then
                
                   nr(k) = max(nr(k),QNSMALL)
  !.. Get rain size distribution properties
                   lamr = (PI*RHOW*nr(k)/qr(k))**0.333333333
                   n0rr = nr(k)*lamr
                   if(lamr.LT.LAMMINR) then
                      lamr = LAMMINR
                      n0rr = lamr**4*qr(k)/(PI*RHOW)
                      nr(k) = n0rr/lamr
                   elseif(lamr.gt.LAMMAXR) then
                      lamr = LAMMAXR
                      n0rr = lamr**4*qr(k)/(PI*RHOW)
                      nr(k) = n0rr/lamr
                   endif
                   rrr = 0.5 * (1./lamr)
                   
  !.. The lookup table is in log space.  Do not change the following 4 lines without
  !.. changing the offline-built lookup table
                   rrni = 13.498*log10(0.5e6*rni(cc))
                   rrri= 23.273*log10(1.e5*rrr)
                   rdsi = 50.*(deltastr(cc) - 0.5)
                   rrho = 7.888*log10(0.02*rhobar(cc))
                   
                   irni = int(rrni)
                   irri = int(rrri)
                   idsi = int(rdsi)
                   irho = int(rrho)
                   
  !.. Small lookup table limit
                   rrni = max(rrni,1.)
                   rrri = max(rrri,1.)
                   rdsi = max(rdsi,1.)
                   rrho = max(rrho,1.)
                   
                   irni = max(irni,1)
                   irri = max(irri,1)
                   idsi = max(idsi,1)
                   irho = max(irho,1)
                   
  !.. Large lookup table limit
                   rrni = min(rrni,real(size(itabr,1))-1.)
                   rrri = min(rrri,real(size(itabr,2))-1.)
                   rdsi = min(rdsi,real(size(itabr,3))-1.)
                   rrho = min(rrho,real(size(itabr,4))-1.)
                
                   irni = min(irni,size(itabr,1)-1)
                   irri = min(irri,size(itabr,2)-1)
                   idsi = min(idsi,size(itabr,3)-1)
                   irho = min(irho,size(itabr,4)-1)
                   
                   do iti = 1, 6
                      call access_lookup_table(itabr,irni,irri,idsi,irho,iti, &
                           rdsi,rrho,rrri,rrni,procr(cc,iti))
                   enddo
                   
  !.. procr(cc,4) is the normalized ice number tendency from ice-rain
  !.. procr(cc,5) is the normalized rain mass tendency from ice-rain
  !.. procr(cc,6) is the normalized ice mass tendency from ice-rain
                   numrateri(cc) = max((procr(cc,4)*ni(cc,k)*nr(k)*rhoair(k)),0.)
                   rainrateri(cc)= max((procr(cc,5)*ni(cc,k)*nr(k)*rhoair(k)),0.)
                   icerateri(cc) = max((procr(cc,6)*ni(cc,k)*nr(k)*rhoair(k)),0.)
                   
  !.. Limit ice-rain rates if to small
                   if(rainrateri(cc).lt.QSMALL.or.icerateri(cc).lt.QSMALL.or.cc.eq.IICE3) then
                      rainrateri(cc)=0.
                      icerateri(cc) =0.
                      numrateri(cc) =0.
                   endif
                   
                   if(temp.le.T0) then
  !.. For temperatures less than melting point ice-rain freezes rain  
                      if(qr(k).gt.0.1e-3.and.qi(cc,k).gt.0.1e-3) then
                         dQRfzri(cc) = rainrateri(cc)
                         dQIfzri(cc) = icerateri(cc)
                         dNfzri(cc)  = numrateri(cc)
                      endif
                   else
  !.. For temperatures greater than melting point ice-rain melts ice  
                      dQImltri(cc) = icerateri(cc)
                      dNmltri(cc)  = numrateri(cc)
                   endif
                   
  !.. procr(cc,1) is the normalized riming rate
  !.. procr(cc,2) is related to the normalized rime density
  !.. procr(cc,3) is the normalized loss in rain number from riming
  !.. rimesumr is in kg m^-3 s^-1

  !.. Ice-rain leads to riming when the fall speed of the collected rain is slower than
  !.. the fall speed of the ice collector.  Otherwise, ice-rain freezes rain
                   rimesumr(cc) = max((procr(cc,1)*ni(cc,k)*nr(k)*rhoair(k)**2),0.)
                   qi_qr_nrm(cc)= procr(cc,1)
                   qi_qr_nrd(cc)= procr(cc,2)
                   qi_qr_nrn(cc)= procr(cc,3)

  !.. Limit ice-rain riming if small
                   if((rimesumr(cc)*i_rhoair(k)).lt.QSMALL.or.cc.eq.IICE3) then
                      rimesumr(cc) =0.
                      qi_qr_nrm(cc)=0.
                      qi_qr_nrd(cc)=0.
                      qi_qr_nrn(cc)=0.
                   endif
                
                endif !.. End ice-rain collection
                
  !.. Do not over-deplete cloud water or rain from ice-cloud and ice-rain
                tmpsum = tmpsum + rimesum(cc)
                tmpsumr= tmpsumr+ rimesumr(cc)
                
                if(ccvar.eq.numice) then !.. Check at end

  !.. Do not let ice-cloud collection overdeplete cloud water
                   sink = tmpsum*i_rhoair(k)
                   source = qc(k)*i_dt 
                   if(sink.gt.source.and.qc(k).gt.QSMALL) then
                      ratio = source / sink
                      rimesum  = rimesum*ratio
                      qi_qc_nrm= qi_qc_nrm*ratio
                      qi_qc_nrd= qi_qc_nrd*ratio
                   endif
                   
  !.. Do not let ice-rain collection overdeplete rainwater
                   sink = tmpsum*i_rhoair(k)
                   source = qr(k)*i_dt 
                   if(sink.gt.source.and.qr(k).gt.QSMALL) then
                      ratio = source / sink
                      rimesumr  = rimesumr*ratio
                      qi_qr_nrm = qi_qr_nrm*ratio
                      qi_qr_nrd = qi_qr_nrd*ratio
                      qi_qr_nrn = qi_qr_nrn*ratio
                   endif
                endif
                
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Inital volume, iwc, #/m^3 for vapor growth
               if (cc.ne.IICE3) then
                gamma_arg = NU+2.+deltastr(cc)
                gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)

                vi = fourthirdspi*rni(cc)**3*gamma_tab(gi)*i_gammnu     
                iwci = ni(cc,k)*rhobar(cc)*vi*rhoair(k)
                nim3dum=ni(cc,k)*rhoair(k)
                
  !.. Total riming rate (the heating goes into vapor growth)
                rimetotal(cc) = rimesum(cc)+rimesumr(cc)
                alphstr=ao**(1.-deltastr(cc))
  !.. Distribution-averaged capacitance (see Harrington et al. 2013)
                capgam = capacitance_gamma(ani(cc), deltastr(cc), alphstr, i_gammnu)
               elseif (cc.eq.IICE3) then
                nim3dum = ni(cc,k)*rhoair(k)
                iwci = qi(cc,k)*rhoair(k)
                call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))
                capgam = NU_AGG*ani(cc)*fs_agg
                !ani(cc) = ((qi(cc,k)/ni(cc,k))*gamma(NU_AGG)/(alphv_agg*phiba_phica_bar*rhobar_ellip*gamma(NU_AGG+3.)))**(1./3.)
                !ani(cc) = min(max(ani(cc),1.e-8),1.e-2)
                alphm_agg = (qi(cc,k)/ni(cc,k))*(1./phiba_phica_bar)*ani(cc)**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)

               endif


  !.. This subroutine first calculates the bulk fall speeds, then the ventilation, 
  !.. and finally the vapor growth (or sublimation) rate. (See Harrington et al. 2013)
                call vaporgrow(cc, dt, ani(cc), cni(cc), rni(cc), igr, nim3dum, temp, rimetotal(cc),&
                     pres_e(k), alphstr, sui, sup, qvs, qvi, mu, iwci, rhoair(k), qi(cc,k),    &
                     dv, kt, nsch, npr, gammnu, i_gammnu, fourthirdspi, svpi, xxls, xxlv, xxlf,&
                     capgam, vtn1(cc,k), vtrni1(cc,k), vtrmi1(cc,k), vtrzi1(cc,k), betvt1(cc,k), &
                     anf, cnf, rnf, iwcf, fv(cc),fh(cc), rhobar(cc), deltastr(cc),rhodepout(cc),&
                     dsdepout(cc),alphm_agg,zetba1(cc,k),zetca1(cc,k))

  !.. Vapor growth mass growth and axes rates
                prd(cc)=(iwcf-iwci)*i_rhoair(k)*i_dt
                prd(cc)=max(prd(cc),-qi(cc,k)*i_dt)
                if (cc.ne.IICE3) then
                ard(cc)=(2.*(anf-ani(cc))*cni(cc)+(cnf-cni(cc))*ani(cc))*ani(cc)*ni(cc,k)*i_dt
                crd(cc)=(2.*(cnf-cni(cc))*ani(cc)+(anf-ani(cc))*cni(cc))*cni(cc)*ni(cc,k)*i_dt
                elseif(cc.eq.IICE3) then
                ard(cc) = 0.
                crd(cc) = 0.
                endif

  !.. Sublimation number loss
                if(prd(cc).lt.0..and.qi(cc,k).gt.QSMALL) then
                   nrd(cc)=prd(cc)*ni(cc,k)/qi(cc,k)
                endif
                
  !.. Limit vapor growth if values small
                if(abs(prd(cc)*dt).lt.(QSMALL*0.01)) then
                   prd(cc)=0.
                   nrd(cc)=0.
                   ard(cc)=0.
                   crd(cc)=0.
                endif

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Now calculated the mass growth and axes growth rates and density for riming

  !.. Check for wet-growth riming conditions
  !.. See equation 12.13 in Lamb and Verlinde (2011)
                if(WET_GROW.and.cc.ne.IICE3) then
                   rimetotal(cc)=rimesum(cc)+rimesumr(cc)
                   call wet_growth_check(temp, rhoair(k), xxlv, xxlf, qv(k), dv, kt, qs0,  &
                        fv(cc), fh(cc), rimetotal(cc), rni(cc), ni(cc,k), dry_growth(cc))


                endif
                
  !.. If above the melting point, assume high-density riming
                if(temp.gt.T0) then
                   dry_growth(cc) = .false.
                endif
                
  !.. Initial rime volume and r-axis length
                vir = vi
                rnfr= rni(cc)
             
  !.. Calculate the ice-cloud rime density
  !.. This is parameterized from Macklin (1962) and is described in
  !.. Jensen et al. 2017 (JAS) 
                if(RIMING.and.qc(k).gt.QSMALL.and.qi_qc_nrm(cc).gt.0..and.cc.ne.IICE3) then
                   rimec1 = 0.0066 !.. Value at -10C
                   if((temp-T0).lt.-30.) then
                      rimec1 = 0.0036
                   elseif((temp-T0).lt.-20.and.(temp-T0).ge.-30.) then
                      dum = (abs((temp-T0))-20.) / 10.
                      rimec1 = dum*(0.0036) + (1.-dum)*(0.004)
                   elseif((temp-T0).lt.-15.and.(temp-T0).ge.-20.) then
                      dum = (abs((temp-T0))-15.) / 5.
                      rimec1 = dum*(0.004) + (1.-dum)*(0.005)
                   elseif((temp-T0).lt.-10.and.(temp-T0).ge.-15.) then
                      dum = (abs((temp-T0))-10.) / 5.
                      rimec1 = dum*(0.005) + (1.-dum)*(0.0066)
                   elseif((temp-T0).lt.-5.and.(temp-T0).ge.-10.) then
                      dum = (abs((temp-T0))-5.) / 5.
                      rimec1 = dum*(0.0066) + (1.-dum)*(0.012)
                   elseif((temp-T0).ge.-5.) then
                      rimec1 = 0.012
                   endif
                   gdenavg(cc) = (1000.*(0.8*tanh(rimec1*qi_qc_nrd(cc)/qi_qc_nrm(cc))+0.1))
                   
                   if((temp-T0).gt.-5..and.(temp-T0).le.0.) then
                      dum = (abs((temp-T0))-0.) / 5.
                      gdenavg(cc) = dum*gdenavg(cc) + (1.-dum)*900.
                   endif
                   
                   if((temp-T0).gt.0..or..NOT.dry_growth(cc)) then
                      gdenavg(cc)=900.
                   endif
  !.. Keep rime density between 50 and 900 kg m^-3
                   gdenavg(cc) = max(gdenavg(cc),50.)
                   gdenavg(cc) = min(gdenavg(cc),900.)
                   
                   gamma_arg = NU+2.+deltastr(cc)
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   
  !.. Change in r-axis from riming (based on the rime density)
                   rimedr = max(((qi_qc_nrm(cc)/gdenavg(cc))/((gamma_tab(gi)*i_gammnu)* &
                        4.*PI*rni(cc)**2)),0.)
                   rnfr   = rnfr + max((rimedr*nc(k)*rhoair(k)),0.)*dt
                   
                endif !.. End ice-cloud rime density calculation
                
  !.. Calculate the ice-rain rime density
  !.. This is parameterized from Macklin (1962)
                if(RAIN_ICE.and.qr(k).gt.QSMALL.and.qi_qr_nrm(cc).gt.0..and.cc.ne.IICE3) then
                   rimec1 = 0.0066 !.. Value at -10C
                   if((temp-T0).lt.-30.) then
                      rimec1 = 0.0036
                   elseif((temp-T0).lt.-20.and.(temp-T0).ge.-30.) then
                      dum = (abs((temp-T0))-20.) / 10.
                      rimec1 = dum*(0.0036) + (1.-dum)*(0.004)
                   elseif((temp-T0).lt.-15.and.(temp-T0).ge.-20.) then
                      dum = (abs((temp-T0))-15.) / 5.
                      rimec1 = dum*(0.004) + (1.-dum)*(0.005)
                   elseif((temp-T0).lt.-10.and.(temp-T0).ge.-15.) then
                      dum = (abs((temp-T0))-10.) / 5.
                      rimec1 = dum*(0.005) + (1.-dum)*(0.0066)
                   elseif((temp-T0).lt.-5.and.(temp-T0).ge.-10.) then
                      dum = (abs((temp-T0))-5.) / 5.
                      rimec1 = dum*(0.0066) + (1.-dum)*(0.012)
                   elseif((temp-T0).ge.-5.) then
                      rimec1 = 0.012
                   endif
                   gdenavgr(cc) = (1000.*(0.8*tanh(rimec1*qi_qr_nrd(cc)/qi_qr_nrm(cc))+0.1))
                   
                   if((temp-T0).gt.-5..and.(temp-T0).le.0.) then
                      dum = (abs((temp-T0))-0.) / 5.
                      gdenavgr(cc) = dum*gdenavgr(cc) + (1.-dum)*900.
                   endif
                   
                   if((temp-T0).gt.0..or..NOT.dry_growth(cc)) then
                      gdenavgr(cc)=900.
                   endif
  !.. Keep rime density between 50 and 900 kg m^03
                   gdenavgr(cc) = max(gdenavgr(cc),50.)
                   gdenavgr(cc) = min(gdenavgr(cc),900.)
                   
                   gamma_arg = NU+2.+deltastr(cc)
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   
  !.. Change in r-axs from riming (based on rime density)
                   rimedrr = max(((qi_qr_nrm(cc)/gdenavgr(cc))/((gamma_tab(gi)*i_gammnu)* &
                        4.*PI*rni(cc)**2)),0.)
                   rnfr    = rnfr + max((rimedrr*nr(k)*rhoair(k)),0.)*dt
                   
                endif !.. End ice-rain rime density calculation
             
  !.. Caculate the volume after riming
                if (cc.ne.IICE3) then
                gamma_arg = NU+2.+deltastr(cc)
                gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                vfr = fourthirdspi*rnfr**3*gamma_tab(gi)*i_gammnu
                vfr = max(vfr,vi)
                rnfr= max(rnfr,rni(cc))
                
  !.. Determine the fraction of rime mass from cloud versus small rain drops
  !.. Update the ice density and iwc
                if(rimetotal(cc).gt.0.and.vfr.gt.vi) then
                   qcrimefrac(cc) = rimesum(cc)/rimetotal(cc)
                   qcrimefrac(cc) = max(qcrimefrac(cc),0.)
                   qcrimefrac(cc) = min(qcrimefrac(cc),1.)
                   gdentotal = qcrimefrac(cc)*gdenavg(cc) + (1.-qcrimefrac(cc))*gdenavgr(cc)
                   rhobarrime = rhobar(cc)
                   rhorimeout(cc) = (rhobar(cc)*(vi/vfr))+(gdentotal*(1.-(vi/vfr)))
                   rhorimeout(cc) = min(rhorimeout(cc),RHOI)
                   iwcfr = rhorimeout(cc) * vfr * nim3dum
                   rhorimeout(cc) = gdentotal
                   rhorimeout(cc) = min(rhorimeout(cc),RHOI)
                else
                   iwcfr= iwci
                   rnfr = rni(cc)
                endif
             
  !.. For dry growth the drops freeze on contact and one axis (the axis perpendicular to the
  !.. fall direction) grows
                if(dry_growth(cc).and.cc.ne.IICE3) then
                   cnfr =cni(cc)
                   anfr =ani(cc)
                   phibr= 1.
                   if(rimetotal(cc).gt.0.and.vfr.gt.vi) then
                      gamma_arg = NU-1.+deltastr(cc)
                      gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                      phibr = cni(cc)/ani(cc)*gamma_tab(gi)*i_gammnu                   
                      if(phibr.eq.1.) then !.. Spherical ice
                         phifr= 1.
                         anfr = rnfr
                         cnfr = rnfr
                      elseif(phibr.gt.1.25) then !.. Prolate ice
                         phifr = phibr*(((rnfr/rni(cc))**3)**(-0.5))
                         anfr = (ani(cc)/rni(cc)**(1.5))*rnfr**(1.5)
                      elseif(phibr.lt.0.8) then !.. Oblate ice
                         phifr = phibr*(rnfr/rni(cc))**3
                         gamma_arg = NU-1.+deltastr(cc)
                         gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         cnfr = phifr*anfr*gammnu/gamma_tab(gi)
                      else
                         phifr = phibr !.. Constant aspect ratio if quasi-spherical
                         gamma_arg = NU-1.+deltastr(cc)
                         gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         gamma_arg = NU+2.+deltastr(cc)
                         gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         anfr = ((vfr*gamma_tab(gi))/(fourthirdspi*phifr*gamma_tab(gi2)))**0.3333333333
                         cnfr = phifr*anfr*gammnu/gamma_tab(gi)
                      endif

  !.. For long time steps, riming can cause oblate ice to become prolate or vise versa
  !.. This prevents that from occuring
                      if(phibr.le.1..and.phifr.gt.1.) then
                         phifr = 0.99
                         gamma_arg = NU-1.+deltastr(cc)
                         gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         gamma_arg = NU+2.+deltastr(cc)
                         gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         anfr = ((vfr*gamma_tab(gi))/(fourthirdspi*phifr*gamma_tab(gi2)))**0.3333333333
                         cnfr = phifr*anfr*gammnu/gamma_tab(gi)
                      endif
                      
                      if(phibr.ge.1..and.phifr.lt.1.) then
                         phifr = 1.01
                         gamma_arg = NU-1.+deltastr(cc)
                         gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         gamma_arg = NU+2.+deltastr(cc)
                         gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                         anfr = ((vfr*gamma_tab(gi))/(fourthirdspi*phifr*gamma_tab(gi2)))**0.3333333333
                         cnfr = phifr*anfr*gammnu/gamma_tab(gi)
                      endif
                   endif !.. End axis evolution from riming

  !.. Don't let the axes shrink from riming
                   cnfr=max(cnfr,cni(cc))
                   anfr=max(anfr,ani(cc))

  !.. Riming mass growth and axes rates
                   if (cc.ne.IICE3) then
                   prdr(cc)=(iwcfr-iwci)*i_rhoair(k)*i_dt
                   elseif(cc.eq.IICE3) then
                   prdr(cc) = 0.

                   endif

                   ardr(cc)=(2.*(anfr-ani(cc))*cni(cc)+(cnfr-cni(cc))*ani(cc))*ani(cc)*ni(cc,k)*i_dt
                   crdr(cc)=(2.*(cnfr-cni(cc))*ani(cc)+(anfr-ani(cc))*cni(cc))*cni(cc)*ni(cc,k)*i_dt

  !.. These can never be negative
                   prdr(cc) = max(prdr(cc),0.)
                   ardr(cc) = max(ardr(cc),0.)
                   crdr(cc) = max(crdr(cc),0.)
                   
  !.. If no rime mass, no axis growth
                   if(prdr(cc).eq.0.0) then
                      ardr(cc) = 0.
                      crdr(cc) = 0.
                   endif

                else !.. Wet Growth

  !.. For wet growth, mass is added but soaks into the particle and thus
  !.. the density increases but the size does not
                   if (cc.ne.IICE3) then
                   prdr(cc)=(iwcfr-iwci)*i_rhoair(k)*i_dt
                   elseif (cc.eq.IICE3) then
                   prdr(cc) = 0.
                    endif
                   ardr(cc)=0.
                   crdr(cc)=0.
                   
                endif !.. End dry versus wet growth

                endif !.. End non-aggregages riming

                if(temp.gt.T0) then

  !.. Melting: ice shapes become more spherical
                   if (cc.ne.IICE3) then
                   qmlt(cc)=2.*PI*(kt*fh(cc)*(T0-temp)+rhoair(cc)*xxlv*dv*fv(cc)*(qs0-qv(k)))/   & 
                        xxlf*(ni(cc,k)*NU*max(ani(cc),cni(cc))) - &
                        (CPW/xxlf*(temp-T0)*(rimetotal(cc)/rhoair(k) + dQImltri(cc)))
                   else

                    qmlt(cc)=2.*PI*(kt*fh(cc)*(T0-temp)+rhoair(cc)*xxlv*dv*fv(cc)*(qs0-qv(k)))/   &
                        xxlf*(ni(cc,k)*NU_AGG*ani(cc))

                   endif
                   
                   qmlt(cc)=min(qmlt(cc),0.)
                   qmlt(cc)=max(qmlt(cc),(-qi(cc,k)*i_dt))

  !.. Don't let very small number mixing ratios (large sizes) cause ice precipitation

                    if(qmlt(cc).lt.0.) then

                     if(cc.lt.IICE3) then
                      if(ai(cc,k).lt.1.e-12.or.ci(cc,k).lt.1.e-12) then
                         qmlt(cc)=-qi(cc,k)*i_dt
                      endif

                     elseif(cc.eq.IICE3) then

                      if(qi(cc,k).lt.QSMALL) then

                        qmlt(cc)=-qi(cc,k)*i_dt

                      endif

                     endif
                   endif

                   nmlt(cc) = max((-ni(cc,k)*i_dt),(ni(cc,k)*qmlt(cc)/qi(cc,k)-dNmltri(cc)))

                   if (cc.ne.IICE3) then
                   gamma_arg = NU+2.+deltastr(cc)
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   tmpmelt=fourthirdspi*alphstr*gamma_tab(gi)*i_gammnu
                   
                   amlt(cc) = (1./tmpmelt)*(1./rhobar(cc))*qmlt(cc) + ai(cc,k)*nmlt(cc)/ni(cc,k)
                   cmlt(cc) = amlt(cc)*cni(cc)/ani(cc)*(1.+(2.*deltastr(cc)))/(2.+deltastr(cc))
                   endif

                endif
             enddo    !.. End the loop over ice speices
          endif       !.. End if ice
          
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Nucleation: this occurs outside of the ice species loop
          
  !.. Freeze evertying at temperatures below -35C 
  !.. That is how it happens in nature, right?
          if(FREEZE_QC.and.qc(k).gt.QSMALL.and.temp.lt.(T0-35.)) then
             mim = qc(k)*i_dt
             nim = nc(k)*i_dt
          endif
          
          if(FREEZE_QC.and.qr(k).gt.QSMALL.and.temp.lt.(T0-35.)) then
             mimr = qr(k)*i_dt
             nimr = nr(k)*i_dt
          endif

  !.. Bigg (1953) freezing below -4C
          if(qr(k).gt.QSMALL.and.temp.lt.(T0-4.)) then
             nr(k) = max(nr(k),QNSMALL)
             lamr = (PI*RHOW*nr(k)/qr(k))**0.333333333
             n0rr = nr(k)*lamr
             if(lamr.LT.LAMMINR) then
                lamr = LAMMINR
                n0rr = lamr**4*qr(k)/(PI*RHOW)
                nr(k) = n0rr/lamr
             elseif(lamr.gt.LAMMAXR) then
                lamr = LAMMAXR
                n0rr = lamr**4*qr(k)/(PI*RHOW)
                nr(k) = n0rr/lamr
             endif
             
             mbiggr = 20.*PI*PI*RHOW*100.*nr(k)*(exp(0.66*(T0-temp))-1.)/lamr**3/lamr**3          
             nbiggr = PI*100.*nr(k)*(exp(0.66*(T0-temp))-1.)/lamr**3
             mbiggr = min(mbiggr,qr(k)*i_dt)
             nbiggr = min(nbiggr,nr(k)*i_dt)
          endif

  !.. DeMott (2010) PNAS, Heterogeneous nuclation
          if(temp.lt.T0.and.sup.ge.0.) then

             a_demott = 0.0000594
             b_demott = 3.33
             c_demott = 0.0264
             d_demott = 0.0033
    
  !.. The 0.03 is the number of large aerosol (>0.5 microns)
  !.. See Chagnon and Junge 1962
             inrate = a_demott * (273.16 - temp)**b_demott * &
                  0.03**((c_demott*(273.16 - temp)) + d_demott)
             inrate = inrate*1000.0 !.. #/L to #/m^3
             inrate = inrate*i_rhoair(k) !.. # kg^-1
             
             if((((ni(IICE1,k)+ni(IICE2,k)+inrate)*rhoair(k))/1000.).le.10000.) then
                mnuccd=(inrate)*(fourthirdspi*RHOI*(2.e-6)**3)*i_dt
                nnuccd=inrate*i_dt
             else
                curnum = (ni(IICE1,k)+ni(IICE2,k))*rhoair(k)/1000.
                ratel = max(0.,(10000.-curnum))
                ratekg = ratel*1000.*i_rhoair(k)
                mnuccd=(ratekg)*(fourthirdspi*RHOI*(2.e-6)**3)*i_dt
                nnuccd=ratekg*i_dt
             endif
          endif

  !.. Rime splinter (Hallett and Mossop 1974, Nature)
          if(SPLINTERS) then
             if(temp.lt.270.16.and.temp.gt.265.16) then
                if (temp.GT.270.16) then
                   fmult = 0.
                elseif(temp.le.270.16.and.temp.gt.268.16) then
                   fmult = (270.16-temp)/2.
                elseif(temp.ge.265.16.and.temp.le.268.16) then
                   fmult = (temp-265.16)/3.
                elseif(temp.lt.265.16) then
                   fmult = 0.
                endif
             endif

             do cc = 1, cat
                if(prdr(cc).gt.0.) then

                   nmult(cc) = 35.E4*prdr(cc)*fmult*1000.
                   qmult(cc) = nmult(cc)*fourthirdspi*RHOI*(5.e-6)**3
                   qmult(cc) = min(qmult(cc),prdr(cc))
                   prdr(cc)  = prdr(cc)-qmult(cc)
                endif
             enddo
          endif




  !.. Aggregation (currently only for T < T0)
        if(AGGREGATION) then
          if(temp.le.T0) then

             do cc = 1, cat
                ni(cc,k) = max(ni(cc,k),QNSMALL)
                ai(cc,k) = max(ai(cc,k),QASMALL)
                ci(cc,k) = max(ci(cc,k),QASMALL)

             enddo

                
  !.. Get characteristic diameters for aggregation
             an1      = ((ai(IICE1,k)**2)/(ci(IICE1,k)*ni(IICE1,k)))**0.333333333333
             an2      = ((ai(IICE2,k)**2)/(ci(IICE2,k)*ni(IICE2,k)))**0.333333333333

             !an3      = 2.*((ai(IICE3,k)**2)/(ci(IICE3,k)*ni(IICE3,k)))**0.333333333333
             an3 = ((qi(IICE3,k)/ni(IICE3,k))*gamma(NU_AGG)/(alphv_agg*phiba_phica_bar*rhobar_ellip*gamma(NU_AGG+3.)))**(1./3.)
             an3 = min(max(an3,1.e-8),1.e-2)
             alphm_agg = (qi(IICE3,k)/ni(IICE3,k))*(1./phiba_phica_bar)*an3**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)


             an1 = MIN(an1,0.5*1.e-2)
             an2 = MIN(an2,0.5*1.e-2)
             !an3 = MIN(an3,0.5*1.e-2)
                
             an1 = MAX(an1,0.5*1.e-6)
             an2 = MAX(an2,0.5*1.e-6)
             !an3 = MAX(an3,0.5*1.e-6)
                
  !.. Get ice-one and ice-two aspect ratios 
  !.. Spherical paritcles aggregate at a slower rate than less spherical ones
             gamma_arg = NU-1.+deltastr(IICE1)
             gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
             phiagg1=ci(IICE1,k)/ai(IICE1,k)*gamma_tab(gi)*i_gammnu
                
             gamma_arg = NU-1.+deltastr(IICE2)
             gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
             phiagg2=ci(IICE2,k)/ai(IICE2,k)*gamma_tab(gi)*i_gammnu
                
             phiagg1=MIN(phiagg1,100.)
             phiagg1=MAX(phiagg1,0.01)
                
             phiagg2=MIN(phiagg2,100.)
             phiagg2=MAX(phiagg2,0.01)


    ! projected area-dimensional relationships for free ice categories

            do cc=1,cat

            if (cc.ne.IICE3) then

            alphstr = ao**(1.-deltastr(cc))

            if (deltastr(cc).gt.1.0) then
                alph_area(cc) = PI*alphstr
                bet_area(cc) = 1.+deltastr(cc)
            else
                alph_area(cc) = PI
                bet_area(cc) = 2.0
            endif

    ! Mass-dimensional relationships
            alph_m(cc) = fourthirdspi*alphstr*rhobar(cc)
            alph_v(cc) = fourthirdspi*alphstr
            bet_m(cc) = 2.+deltastr(cc)
            alphc(cc) = fourthirdspi*alphstr**2
            betac(cc) = 1.+2.*deltastr(cc)



            elseif (cc.eq.IICE3) then

            alph_m(cc) = alphm_agg
            bet_m(cc) = betm_agg
            alph_area(cc) = PI
            bet_area(cc) = 2.0


            endif

            enddo


            Viadum(1) = alph_v(1)*an1**bet_m(1)*ni(1,k)*gamma(NU+bet_m(1))/gammnu
            Viadum(2) = alph_v(2)*an2**bet_m(2)*ni(2,k)*gamma(NU+bet_m(2))/gammnu

            Vicdum(1) = alphc(1)*an1**betac(1)*ni(1,k)*gamma(NU+betac(1))/gammnu
            Vicdum(2) = alphc(2)*an2**betac(2)*ni(2,k)*gamma(NU+betac(2))/gammnu




            !---------------------------------------------------!
            ! So total combinations for aggregation are:
            !---------------------------------------------------!
            !    Collected + Collectee = New Category
            !---------------------------------------------------!
            !    Plates + Plates   = Agg                    (1)
            !    Columns + Columns = Agg                    (2)
            !    Plates + Columns  = Agg                    (3)
            !    Columns + Plates  = Agg (mass/volume only) (4)
            !    Plates + Agg      = Agg                    (5)
            !    Columns + Agg     = Agg                    (6)
            !    Agg + Agg         = Agg (number only)      (7)
            !---------------------------------------------------!
            !---------------------------------------------------!
            ! total dvt terms
            !---------------------------------------------------!
            ! Number: 6 combos
            !---------------------------------------------------!
            ! dvtn: Plates + Plates (1)
            ! dvtn: Columns + Columns (2)
            ! dvtn: Plates + Columns (3)
            ! dvtn: Plates + Aggs (4)
            ! dvtn: Columns + Aggs (5)
            ! dvtn: Aggs + Aggs (6)
            !---------------------------------------------------!
            ! Mass: 6 combos
            ! dvtm: Plates + Plates (1)
            ! dvtm: Columns + Columns (2)
            ! dvtm: Plates + Columns (3)
            ! dvtm: Columns + Plates (4)
            ! dvtm: Plates + Aggs (5)
            ! dvtm: Columns + Aggs (6)
            !---------------------------------------------------!
            ! Volume (same as mass):
            !---------------------------------------------------!
            ! Volume*aspect ratio: 6 combos
            ! dvtc: Plates + Plates (1)
            ! dvtc: Columns + Columns (2)
            ! dvtc: Plates + Columns (3)
            ! dvtc: Columns + Plates (4)
            ! dvtc: Plates + Aggs (5)
            ! dvtc: Columns + Aggs (6)
            !---------------------------------------------------!
            !---------------------------------------------------!

            if(dvt_switch.eq.0) then

            ! Calculate fallspeed difference term for each specie combo
            ! Mass
            dvtm(1) = delta_vtm(vtrmi1(IICE1,k),vtrmi1(IICE1,k)) ! P-P
            dvtm(2) = delta_vtm(vtrmi1(IICE2,k),vtrmi1(IICE2,k)) ! C-C
            dvtm(3) = delta_vtm(vtrmi1(IICE1,k),vtrmi1(IICE2,k)) ! P-C and C-P
            dvtm(4) = delta_vtm(vtrmi1(IICE1,k),vtrmi1(IICE3,k)) ! P-A and A-P
            dvtm(5) = delta_vtm(vtrmi1(IICE2,k),vtrmi1(IICE3,k)) ! C-A and A-C

            ! Volume*Aspect Ratio
            dvtc(1) = dvtm(1)
            dvtc(2) = dvtm(2)
            dvtc(3) = dvtm(3)
            dvtc(4) = dvtm(4)
            dvtc(5) = dvtm(5)

            ! Number
            dvtn(1) = delta_vtm(vtrni1(IICE1,k),vtrni1(IICE1,k)) ! P-P
            dvtn(2) = delta_vtm(vtrni1(IICE2,k),vtrni1(IICE2,k)) ! C-C
            dvtn(3) = delta_vtm(vtrni1(IICE1,k),vtrni1(IICE2,k)) ! P-C and C-P
            dvtn(4) = delta_vtm(vtrni1(IICE1,k),vtrni1(IICE3,k)) ! P-A and A-P
            dvtn(5) = delta_vtm(vtrni1(IICE2,k),vtrni1(IICE3,k)) ! C-A and A-C
            dvtn(6) = delta_vtm(vtrni1(IICE3,k),vtrni1(IICE3,k)) ! A-A

            elseif(dvt_switch.eq.1) then
            ! Calculate fallspeed difference term for each specie combo
            dvtm(1) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),bet_m(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE1,k),bet_area(IICE1),betvt1(IICE1,k),&
                      zetba1(IICE1,k),zetca1(IICE1,k),0) !P-P

            dvtm(2) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),bet_m(IICE2),&
                      betvt1(IICE2,k),vtn1(IICE2,k),bet_area(IICE2),betvt1(IICE2,k),&
                      zetba1(IICE2,k),zetca1(IICE2,k),0) ! C-C

            dvtm(3) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),bet_m(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE2,k),bet_area(IICE2),betvt1(IICE2,k),&
                      zetba1(IICE2,k),zetca1(IICE2,k),0) !P-C

            dvtm(4) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),bet_m(IICE2),&
                      betvt1(IICE2,k),vtn1(IICE1,k),bet_area(IICE1),betvt1(IICE1,k),&
                      zetba1(IICE1,k),zetca1(IICE1,k),0) !C-P

            dvtm(5) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),bet_m(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) !P-A

            dvtm(6) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),bet_m(IICE2),&
                      betvt1(IICE2,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) !C-A


            dvtc(1) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),betac(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE1,k),bet_area(IICE1),betvt1(IICE1,k),&
                      zetba1(IICE1,k),zetca1(IICE1,k),0) !P-P

            dvtc(2) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),betac(IICE2),&
                      betvt1(IICE2,k),vtn1(IICE2,k),bet_area(IICE2),betvt1(IICE2,k),&
                      zetba1(IICE2,k),zetca1(IICE2,k),0) !C-C

            dvtc(3) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),betac(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE2,k),bet_area(IICE2),betvt1(IICE2,k),&
                      zetba1(IICE2,k),zetca1(IICE2,k),0) !P-C

            dvtc(4) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),betac(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) !C-P

            dvtc(5) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),betac(IICE1),&
                      betvt1(IICE1,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) !P-A

            dvtc(6) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),betac(IICE2),&
                      betvt1(IICE2,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) !C-A



            dvtn(1) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),0.,&
                      betvt1(IICE1,k),vtn1(IICE1,k),bet_area(IICE1),betvt1(IICE1,k),&
                      zetba1(IICE1,k),zetca1(IICE1,k),0) ! P-P

            dvtn(2) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),0.,&
                      betvt1(IICE2,k),vtn1(IICE2,k),bet_area(IICE2),betvt1(IICE2,k),&
                      zetba1(IICE2,k),zetca1(IICE2,k),0) !C-C

            dvtn(3) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),0.,&
                      betvt1(IICE1,k),vtn1(IICE2,k),bet_area(IICE2),betvt1(IICE2,k),&
                      zetba1(IICE2,k),zetca1(IICE2,k),0) ! P-C

            dvtn(4) = delta_vtm_seifert(vtn1(IICE1,k),bet_area(IICE1),0.,&
                      betvt1(IICE1,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) ! P-A

            dvtn(5) = delta_vtm_seifert(vtn1(IICE2,k),bet_area(IICE2),0.,&
                      betvt1(IICE2,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),1) ! C-A

            dvtn(6) = delta_vtm_seifert(vtn1(IICE3,k),bet_area(IICE3),0.,&
                      betvt1(IICE3,k),vtn1(IICE3,k),bet_area(IICE3),betvt1(IICE3,k),&
                      zetba1(IICE3,k),zetca1(IICE3,k),2) ! A-A




            endif


                
  !.. Ice-1 and ice-2 mass (volume), number, and aspect ratio-weighted volume (Vc)
  !.. loss to aggregates occurs inside of aggregation
  !.. After aggregation, ice-1 and ice-2 bulk distribution properties are adjusted 
  !.. according the new mass, volumes, and number mixing ratios (shape_diag).

             if(qi(IICE1,k).gt.1.e-8.or.qi(IICE2,k).gt.1.e-8.or.qi(IICE3,k).gt.1.e-8) then

                call starbuck(rhoair(k),temp,phiagg1,phiagg2,rhobar(IICE1),rhobar(IICE2),qi(IICE1,k),Viadum(IICE1),Vicdum(IICE1), &
                              ni(IICE1,k),an1,alph_area(IICE1),bet_area(IICE1),alph_v(IICE1),alph_m(IICE1),bet_m(IICE1),alphc(IICE1),&
                              betac(IICE1),qi(IICE2,k),Viadum(IICE2),Vicdum(IICE2),ni(IICE2,k),an2,alph_area(IICE2),bet_area(IICE2),&
                              alph_v(IICE2),alph_m(IICE2),bet_m(IICE2),alphc(IICE2),betac(IICE2),qi(IICE3,k),Viadum(IICE3),&
                              Vicdum(IICE3),ni(IICE3,k),an3,alph_area(IICE3),bet_area(IICE3),dvtn,dvtm,dvtc,qagg(IICE1),&
                              qagg(IICE2),qagg(IICE3),nagg(IICE1),nagg(IICE2),nagg(IICE3),nmon3,Vagg(IICE1),Vagg(IICE2),Vagg(IICE3),&
                              Vcagg(IICE1),Vcagg(IICE2),Vcagg(IICE3))



             endif  !.. If ice mass mixing ratios are high enough to aggregate
          endif
        endif !.. If aggregation
          
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Warm Microphysics

  !.. Get raindrop size distribution parameters
          if (qr(k).gt.QSMALL) then
             nr(k) = max(nr(k),QNSMALL)
             lamr = (PI*RHOW*nr(k)/qr(k))**0.333333333
             n0rr = nr(k)*lamr
             if(lamr.LT.LAMMINR) then
                lamr = LAMMINR
                n0rr = lamr**4*qr(k)/(PI*RHOW)
                nr(k) = n0rr/lamr
             elseif(lamr.gt.LAMMAXR) then
                lamr = LAMMAXR
                n0rr = lamr**4*qr(k)/(PI*RHOW)
                nr(k) = n0rr/lamr
             endif
          endif
          
  !.. note: nprc1 is change in Nr,                                                   
  !.. nprc is change in Nc                         
  !.. NOTE: NPRC IS CALCULATED BUT NOT CURRENTLY USED SINCE NC IS NOT PROGNOSED 
  !.. Auto conversion of cloud water to rain (Khairoutdinov and Kogan 2000)
          if(qc(k).ge.1.E-6) then
             PRC=1350.*qc(k)**2.47*  &
                  (nc(k)/1.e6*rhoair(k))**(-1.79)
             NPRC1 = PRC/(4./3.*PI*RHOW*(25.E-6)**3)
             NPRC  = PRC/(qc(k)/nc(k))
             NPRC  = min(NPRC,nc(k)*i_dt)
             NPRC1 = min(NPRC1,NPRC)
          endif
       
  !.. Accretion of cloud water by rain Khairoutdinov and Kogan 2000)
  !.. NOTE: NPRA CURRENTLY NOT USED SINCE NC IS NOT PROGNOSED
          if(qr(k).ge.1.e-8.and.qc(k).ge.1.e-8) then
             DUM=(qc(k)*qr(k))
             PRA = 67.*(DUM)**1.15
             NPRA = PRA/(qc(k)/nc(k))                                    
          endif
          
  !.. Self-collection (Beheng 1994), Breakup (Verlinde and Cotton 1993)
          if(qr(k).ge.1.e-8) then
             dum1=300.e-6
             if(1./lamr.lt.dum1) then
                dum=1.
             elseif(1./lamr.ge.dum1) then
                dum=2.-exp(2300.*(1./lamr-dum1))
             endif
             nragg = -5.78*dum*nr(k)*qr(k)*rhoair(k)
          endif
          
  !.. Evaporation of rain (Rutledge and Hobbs 1983)
          if(qr(k).gt.QSMALL) then
             epsr = 2.*PI*n0rr*rhoair(k)*dv* &
                  (F1R/(LAMR*LAMR)+ &
                  F2R*(arn*rhoair(k)/mu)**0.5* &
                  nsch**0.333333333*(gamma(5./2.+BR/2.))/ &
                  (lamr**(5./2.+BR/2.)))
          else
             epsr = 0.
          endif
          
          if(qv(k).lt.qvs) then
             pre = epsr*(qv(k)-qvs)/ab
             pre = MIN(PRE,0.)
             dumt = temp+xxlv*i_cp*pre*dt
             dumqv = qv(k)-pre*dt
             dumqvs = 0.622*polysvp(dumt,0)/ &
                  (pres_e(k)-polysvp(dumt,0))
             if(pre.lt.0..and.dumqv/dumqvs.gt.1.) then
                pre=(qv(k)-qvs)*i_dt/ab
             endif
          else
             pre = 0.
          endif
          
  !.. Rain number mixing ratio loss from evaporation
          if(pre.lt.0.) then
             dum = pre*dt/qr(k)
             dum = max(-1.,dum)
             npre = dum*nr(k)/dt
          endif
       
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Updates after process-rate calculations
  !.. For ice, updates to mass, number, volume, and volume
  !.. times aspect ratio occur first for vapor growth, 
  !.. riming, and melting (shape is predicted to evolve).
  !.. Next updates the bulk distribution shapes are done for 
  !.. nucleation and aggregation. For nucleation, ice particle shape
  !.. itself does not change but the distribution shape becomes 
  !.. more spherical as nucleated particles are added.   

  !.. Do not over-deplete
  
  !.. Cloud water
  !.. First sum the total loss from riming
          qcrimesum=0.
          do cc = 1, cat
             qcrimesum = qcrimesum + prdr(cc)*qcrimefrac(cc)
          enddo

          sink  =(prc+pra+qcrimesum+mim)*dt
          source=qc(k)
          if(sink.gt.source.and.qc(k).gt.QSMALL) then
             ratio = source/sink
             prc = prc*ratio
             pra = pra*ratio
             mim = mim*ratio
             do cc = 1, cat
                prdr(cc) = prdr(cc)*ratio
                ardr(cc) = ardr(cc)*ratio
                crdr(cc) = crdr(cc)*ratio
             enddo
          endif

  !.. Rain
  !.. First sum the total loss from riming
          qcrimesumr=0.
          do cc = 1, cat
             qcrimesumr = qcrimesumr + prdr(cc)*(1.-qcrimefrac(cc))
          enddo

          sink = (-pre+qcrimesumr+mimr+mbiggr)*dt
          source = qr(k) + (pra+prc)*dt
          do cc = 1, cat
             source = source + (-qmlt(cc)*dt)
          enddo
          do cc = 1, cat
             sink = sink + (dQRfzri(cc))*dt
          enddo
          if(sink.gt.source.and.qr(k).gt.QSMALL) then
             ratio= source/sink
             do cc = 1, cat
                prdr(cc) = prdr(cc)*ratio
                ardr(cc) = ardr(cc)*ratio
                crdr(cc) = crdr(cc)*ratio
                qi_qr_nrn(cc)=qi_qr_nrn(cc)*ratio
                dQRfzri(cc) = dQRfzri(cc)*ratio
                dQIfzri(cc) = dQIfzri(cc)*ratio
                dNfzri(cc)  = dNfzri(cc)*ratio
             enddo
             mimr       = mimr*ratio
             nimr       = nimr*ratio
             mbiggr     =mbiggr*ratio
             nbiggr     =nbiggr*ratio
             pre        = pre*ratio
          endif

  !.. Ice 
          if(igr.le.1.) then !.. Nucleation initiation to ice 1
  !.. Ice 1 
             sink  =(-qmlt(IICE1)-qagg(IICE1))*dt
             source=qi(IICE1,k) + (prdr(IICE1))*dt
             if(prd(IICE1).gt.0.) then
                source = source + prd(IICE1)*dt
             else
                sink = sink + (-prd(IICE1))*dt
             endif

  !.. Added because of ice initiation to ice 1
             source = source + (mnuccd+mim+mimr+mbiggr+dQIfzri(IICE2)+dQIfzri(IICE3))*dt
             do cc = 1, cat
                source = source + (qmult(cc)+dQRfzri(cc))*dt
             enddo

             if(sink.gt.source.and.qi(IICE1,k).gt.QSMALL) then
                ratio= source/sink
                qmlt(IICE1) = qmlt(IICE1)*ratio
                qagg(IICE1) = qagg(IICE1)*ratio
                nagg(IICE1) = nagg(IICE1)*ratio
                if(prd(IICE1).lt.0.) then
                   prd(IICE1) = prd(IICE1)*ratio
                endif
             endif

  !.. Ice 2
             sink  =(dQIfzri(IICE2)-qmlt(IICE2)-qagg(IICE2))*dt
             source=qi(IICE2,k) + (prdr(IICE2))*dt
             if(prd(IICE2).gt.0.) then
                source = source + prd(IICE2)*dt
             else
                sink = sink + (-prd(IICE2))*dt
             endif

             if(sink.gt.source.and.qi(IICE2,k).gt.QSMALL) then
                ratio= source/sink
                dQIfzri(IICE2) = dQIfzri(IICE2)*ratio
                dQRfzri(IICE2) = dQRfzri(IICE2)*ratio
                dNfzri(IICE2)  = dNfzri(IICE2)*ratio
                qmlt(IICE2) = qmlt(IICE2)*ratio
                qagg(IICE2) = qagg(IICE2)*ratio
                nagg(IICE2) = nagg(IICE2)*ratio
                if(prd(IICE2).lt.0.) then
                   prd(IICE2) = prd(IICE2)*ratio
                endif
             endif
  !.. Ice 3
             sink  =(dQIfzri(IICE3)-qmlt(IICE3))*dt
             source=qi(IICE3,k) + (prdr(IICE3)+qagg(IICE3))*dt
             if(prd(IICE3).gt.0.) then
                source = source + prd(IICE3)*dt
             else
                sink = sink + (-prd(IICE3))*dt
             endif

             if(sink.gt.source.and.qi(IICE3,k).gt.QSMALL) then
                ratio= source/sink
                dQIfzri(IICE3) = dQIfzri(IICE3)*ratio
                dQRfzri(IICE3) = dQRfzri(IICE3)*ratio
                dNfzri(IICE3)  = dNfzri(IICE3)*ratio
                qmlt(IICE3) = qmlt(IICE3)*ratio
                if(prd(IICE3).lt.0.) then
                   prd(IICE3) = prd(IICE3)*ratio
                endif
             endif
          else !.. nucleation to ice 2
  !.. Ice 2 
             sink  =(-qmlt(IICE2)-qagg(IICE2))*dt
             source=qi(IICE2,k) + (prdr(IICE2))*dt
             if(prd(IICE2).gt.0.) then
                source = source + prd(IICE2)*dt
             else
                sink = sink + (-prd(IICE2))*dt
             endif

  !.. Added because of ice initiation to ice 2
             source = source + (mnuccd+mim+mimr+mbiggr+dQIfzri(IICE1)+dQIfzri(IICE3))*dt
             do cc = 1, cat
                source = source + (qmult(cc)+dQRfzri(cc))*dt
             enddo

             if(sink.gt.source.and.qi(IICE2,k).gt.QSMALL) then
                ratio= source/sink
                qmlt(IICE2) = qmlt(IICE2)*ratio
                qagg(IICE2) = qagg(IICE2)*ratio
                nagg(IICE2) = nagg(IICE2)*ratio
                if(prd(IICE2).lt.0.) then
                   prd(IICE2) = prd(IICE2)*ratio
                endif
             endif

  !.. Ice 1
             sink  =(dQIfzri(IICE1)-qmlt(IICE1)-qagg(IICE1))*dt
             source=qi(IICE1,k) + (prdr(IICE1))*dt
             if(prd(IICE1).gt.0.) then
                source = source + prd(IICE1)*dt
             else
                sink = sink + (-prd(IICE1))*dt
             endif

             if(sink.gt.source.and.qi(IICE1,k).gt.QSMALL) then
                ratio= source/sink
                dQIfzri(IICE1) = dQIfzri(IICE1)*ratio
                dQRfzri(IICE1) = dQRfzri(IICE1)*ratio
                dNfzri(IICE1)  = dNfzri(IICE1)*ratio
                qmlt(IICE1) = qmlt(IICE1)*ratio
                qagg(IICE1) = qagg(IICE1)*ratio
                nagg(IICE1) = nagg(IICE1)*ratio
                if(prd(IICE1).lt.0.) then
                   prd(IICE1) = prd(IICE1)*ratio
                endif
             endif
  !.. Ice 3
             sink  =(dQIfzri(IICE3)-qmlt(IICE3))*dt
             source=qi(IICE3,k) + (prdr(IICE3)+qagg(IICE3))*dt
             if(prd(IICE3).gt.0.) then
                source = source + prd(IICE3)*dt
             else
                sink = sink + (-prd(IICE3))*dt
             endif

             if(sink.gt.source.and.qi(IICE3,k).gt.QSMALL) then
                ratio= source/sink
                dQIfzri(IICE3) = dQIfzri(IICE3)*ratio
                dQRfzri(IICE3) = dQRfzri(IICE3)*ratio
                dNfzri(IICE3)  = dNfzri(IICE3)*ratio
                qmlt(IICE3) = qmlt(IICE3)*ratio
                if(prd(IICE3).lt.0.) then
                   prd(IICE3) = prd(IICE3)*ratio
                endif
             endif
          endif

  !.. Do not let vapor growth of ice create subsaturated conditions          
          prdsum=0.
          do cc = 1, cat
             prdsum = prdsum + prd(cc)
          enddo
          vgflag = .false.
          sink = prdsum
          source=0.99*sui*qvi/abi*i_dt
          if(sink.gt.source.and.sui.gt.0.) then
             vgflag = .true.
             ratio = source/sink
             do cc = 1, cat
                prd(cc) = prd(cc)*ratio
             enddo
          endif

  !.. Update to ice
          do cc = 1, cat

  !.. Vapor growth update
             if (cc.ne.IICE3) then
             if((qi(cc,k)+(prd(cc))*dt).gt.QSMALL.and.qi(cc,k).gt.QSMALL) then
                rhobar(cc) = rhobar(cc)*(qi(cc,k)/(qi(cc,k)+(prd(cc))*dt)) + &
                     rhodepout(cc)*(1.-(qi(cc,k)/(qi(cc,k)+(prd(cc))*dt)))
             else
                rhobar(cc) = RHOI
             endif

            endif


            if (cc.ne.IICE3) then

             qi(cc,k)=qi(cc,k)+(prd(cc))*dt
             ni(cc,k)=ni(cc,k)+(nrd(cc))*dt
             ai(cc,k)=ai(cc,k)+(ard(cc))*dt
             ci(cc,k)=ci(cc,k)+(crd(cc))*dt

            elseif(cc.eq.IICE3) then

             Nmonb = max(nmon(k)/ni(cc,k),2.) ! Temp store of avg. mon number
             qi(cc,k)=qi(cc,k)+(prd(cc))*dt
             ni(cc,k)=ni(cc,k)+(nrd(cc))*dt

             nmon(k)=ni(cc,k)*Nmonb ! Re-diagnose number of monomers keeping Nmon_avg const.

             endif

             
             ni(cc,k) = max(ni(cc,k),QNSMALL)
             nmon(k) = max(nmon(k),2.*QNSMALL)
             ai(cc,k) = max(ai(cc,k),QASMALL)
             ci(cc,k) = max(ci(cc,k),QASMALL)



             
  !.. If vapor growth of ice is going to cause
  !.. subsaturated condtions, reduce vapor growth
  !.. rates (done above) and rediagnose the axes
             if(vgflag.and.qi(cc,k).gt.QSMALL) then
                if (cc.ne.IICE3) then
                alphstr=ao**(1.-deltastr(cc))
                alphv=fourthirdspi*alphstr
                betam=2.+deltastr(cc)
                gamma_arg = NU+2.+deltastr(cc)
                gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                ani(cc)=((qi(cc,k)*gammnu)/(rhobar(cc)*ni(cc,k)*alphv* &
                     gamma_tab(gi)))**(1./betam)
                cni(cc)=ao**(1.-deltastr(cc))*ani(cc)**deltastr(cc)
                ai(cc,k)=ani(cc)**2*cni(cc)*ni(cc,k)
                ci(cc,k)=cni(cc)**2*ani(cc)*ni(cc,k)

                else

                call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))
                alphm_agg = (qi(cc,k)/ni(cc,k))*(1./phiba_phica_bar)*ani(cc)**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)


                endif

             endif
          enddo

  !.. Update to ice
          do cc = 1, cat
  !.. Riming update

             if (cc.ne.IICE3) then

             if((qi(cc,k)+(prdr(cc))*dt).gt.QSMALL) then
                if(dry_growth(cc)) then
                   rhobar(cc) = rhobar(cc)*(qi(cc,k)/(qi(cc,k)+(prdr(cc))*dt)) + &
                        rhorimeout(cc)*(1.-(qi(cc,k)/(qi(cc,k)+(prdr(cc))*dt)))
                else
                   rhobar(cc) = rhobar(cc)*(qi(cc,k)/(qi(cc,k)+(prdr(cc))*dt)) + &
                        RHOW*(1.-(qi(cc,k)/(qi(cc,k)+(prdr(cc))*dt)))
                   rhobar(cc) = max(rhobar(cc),RHOI)
                endif
             else
                rhobar(cc) = RHOI
             endif

             qi(cc,k)=qi(cc,k)+(prdr(cc))*dt
             ai(cc,k)=ai(cc,k)+(ardr(cc))*dt
             ci(cc,k)=ci(cc,k)+(crdr(cc))*dt

             endif

          enddo

  !.. Keep everything in bounds after vapor growth and riming
          do cc = 1, cat
             if(qi(cc,k).gt.QSMALL.and.cc.ne.IICE3) then
                ni(cc,k)=max(ni(cc,k),QNSMALL)
                ai(cc,k)=max(ai(cc,k),QASMALL)
                ci(cc,k)=max(ci(cc,k),QASMALL)
                ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
                cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)

                ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)                                          
                ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)

                deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))
                   
                call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                     ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                     alphstr, alphv, betam)

                dsold(cc) = deltastr(cc)

            elseif (cc.eq.IICE3) then

                ni(cc,k)=max(ni(cc,k),QNSMALL)
                nmon(k)=max(nmon(k),2.*QNSMALL)

                call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))

                alphm_agg = (qi(cc,k)/ni(cc,k))*(1./phiba_phica_bar)*ani(cc)**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)

             endif
          enddo

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Update to ice after melting
          if(temp.gt.T0) then
             do cc = 1, cat

                if (cc.ne.IICE3) then
                if(qi(cc,k).gt.QSMALL) then
                   rhobar(cc) = rhobar(cc)*((qi(cc,k)+(qmlt(cc))*dt)/qi(cc,k)) + &
                        RHOW*((-qmlt(cc)*dt)/qi(cc,k))
                   rhobar(cc) = min(rhobar(cc),RHOI)
                else
                   rhobar(cc) = RHOI
                endif
                endif



                qi(cc,k)=qi(cc,k)+qmlt(cc)*dt
                ai(cc,k)=ai(cc,k)+amlt(cc)*dt
                ci(cc,k)=ci(cc,k)+cmlt(cc)*dt


                if (cc.ne.IICE3) then
                ni(cc,k)=ni(cc,k)+nmlt(cc)*dt
                elseif (cc.eq.IICE3) then
                Nmonb = max(nmon(k)/ni(cc,k),2.)
                ni(cc,k)=ni(cc,k)+nmlt(cc)*dt
                nmon(k)=Nmonb*ni(cc,k)
                endif


                
  !.. Small masses of leftover ice can become rain
                if(qi(cc,k).le.QSMALL) then
                   qr(k) = qr(k) + qi(cc,k)
                   qi(cc,k)=0.
                endif
                
  !.. Keep everything in bounds after melting
               if(qi(cc,k).gt.QSMALL.and.cc.ne.IICE3) then
                   ni(cc,k)=max(ni(cc,k),QNSMALL)
                   ai(cc,k)=max(ai(cc,k),QASMALL)
                   ci(cc,k)=max(ci(cc,k),QASMALL)
                   ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
                   cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
                   
                   ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)                                   
                   ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)

                   deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))

  !.. Melting should only make the average distribution more spherical
                   if(dsold(cc).le.1.) then
                      if(deltastr(cc).lt.dsold(cc)) then
                         deltastr(cc)=dsold(cc)
                      endif
                      if(deltastr(cc).gt.1.) then
                         deltastr(cc)=1.
                      endif
                   else
                      if(deltastr(cc).gt.dsold(cc)) then
                         deltastr(cc)=dsold(cc)
                      endif
                      if(deltastr(cc).lt.1.) then
                         deltastr(cc)=1.
                      endif
                   endif
                   cni(cc)=ao**(1.-deltastr(cc))*ani(cc)**deltastr(cc)
                   ai(cc,k)=ani(cc)**2*cni(cc)*ni(cc,k)
                   ci(cc,k)=cni(cc)**2*ani(cc)*ni(cc,k)

                   call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                        ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                        alphstr, alphv, betam)

                   dsold(cc) = deltastr(cc)

              elseif(qi(cc,k).gt.QSMALL.and.cc.eq.IICE3) then
                call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))
                alphm_agg = (qi(cc,k)/ni(cc,k))*(1./phiba_phica_bar)*ani(cc)**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)
                endif  !.. End if ice exists after melting
             enddo     !.. Loop over all ices species
          endif        !.. End temperature above melting
       
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Update to vapor
          qv(k)=qv(k)+(-pre-mnuccd)*dt
          do cc = 1, cat
             qv(k)=qv(k)+(-prd(cc))*dt
          enddo

  !.. Update to cloud
          qc(k)=qc(k)+(-prc-pra-mim)*dt
          do cc = 1, cat
             qc(k) = qc(k)-(prdr(cc)*qcrimefrac(cc)*dt)
          enddo

  !.. Update to rain
          qr(k)=qr(k)+(pre+prc+pra-mimr-mbiggr)*dt
          nr(k)=nr(k)+(nprc1+nragg+npre-nimr-nbiggr)*dt
          do cc = 1, cat
             qr(k)=qr(k)+((-prdr(cc)*(1.-qcrimefrac(cc)))-qmlt(cc)-dQRfzri(cc))*dt
             nr(k)=max((nr(k)+((-qi_qr_nrn(cc)*ni(cc,k)*nr(k)*rhoair(k))-nmlt(cc)-dNfzri(cc))*dt),0.)
          enddo

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Store the number before nucleation for later use
          do cc = 1, cat
             nibnuc(cc)=ni(cc,k)
          enddo

          if(igr.le.1.) then
  !.. Nucleation to IICE1 (planar-nuclated)
             totalnuc(IICE1) =(mnuccd+mim+mimr+mbiggr+dQIfzri(IICE2)+dQIfzri(IICE3))*dt
             totalnucn(IICE1)=(nnuccd+nim+nimr+nbiggr+dNfzri(IICE2)+dNfzri(IICE3))*dt

             do cc = 1, cat
                totalnuc(IICE1) =totalnuc(IICE1)  + (qmult(cc)+dQRfzri(cc))*dt
                totalnucn(IICE1)=totalnucn(IICE1) + (nmult(cc))*dt
             enddo

             qi(IICE1,k)=qi(IICE1,k)+totalnuc(IICE1)
             ni(IICE1,k)=ni(IICE1,k)+totalnucn(IICE1)

  !.. Loss of IICE2 from IICE2-rain producing IICE1
             qi(IICE2,k)=qi(IICE2,k)-(dQIfzri(IICE2))*dt
             ni(IICE2,k)=ni(IICE2,k)-(dNfzri(IICE2))*dt

  !.. Loss of IICE3 from IICE3-rain producing IICE1
             qi(IICE3,k)=qi(IICE3,k)-(dQIfzri(IICE3))*dt
             Nmonb = max(nmon(k)/ni(IICE3,k),2.)
             ni(IICE3,k)=ni(IICE3,k)-(dNfzri(IICE3))*dt
             nmon(k) = ni(IICE3,k)*Nmonb
          else

  !.. Nucleation to IICE2 (columnar-nucleated)
             totalnuc(IICE2) =(mnuccd+mim+mimr+mbiggr+dQIfzri(IICE1)+dQIfzri(IICE3))*dt
             totalnucn(IICE2)=(nnuccd+nim+nimr+nbiggr+dNfzri(IICE1)+dNfzri(IICE3))*dt

             do cc = 1, cat
                totalnuc(IICE2) =totalnuc(IICE2)  + (qmult(cc)+dQRfzri(cc))*dt
                totalnucn(IICE2)=totalnucn(IICE2) + (nmult(cc))*dt
             enddo

             qi(IICE2,k)=qi(IICE2,k)+totalnuc(IICE2)
             ni(IICE2,k)=ni(IICE2,k)+totalnucn(IICE2)

  !.. Loss of IICE1 from IICE1-rain producing IICE2
             qi(IICE1,k)=qi(IICE1,k)-(dQIfzri(IICE1))*dt
             ni(IICE1,k)=ni(IICE1,k)-(dNfzri(IICE1))*dt

  !.. Loss of IICE3 from IICE3-rain producing IICE2
             qi(IICE3,k)=qi(IICE3,k)-(dQIfzri(IICE3))*dt
             Nmonb = max(nmon(k)/ni(IICE3,k),2.)
             ni(IICE3,k)=ni(IICE3,k)-(dNfzri(IICE3))*dt
             nmon(k) = ni(IICE3,k)*Nmonb

          endif

  !.. Limit ice number concentrations to 1000/L
          do cc = 1, cat
             if (cc.ne.IICE3) then
             ni(cc,k)=min(ni(cc,k),(1000.*1000.*i_rhoair(k)))

             elseif (cc.eq.IICE3) then
             Nmonb = max(nmon(k)/ni(cc,k),2.)
             ni(cc,k)=min(ni(cc,k),(1000.*1000.*i_rhoair(k)))
             nmon(k) = ni(cc,k)*Nmonb

             endif

          enddo

  !.. Final check on rain and fall speed calculation
          if (qr(k).gt.QSMALL) then
             nr(k) = max(nr(k),QNSMALL)
             lamr = (PI*RHOW*nr(k)/qr(k))**0.333333333
             n0rr = nr(k)*lamr
             if(lamr.LT.LAMMINR) then
                lamr = LAMMINR
                n0rr = lamr**4*qr(k)/(PI*RHOW)
                nr(k) = n0rr/lamr
             elseif(lamr.gt.LAMMAXR) then
                lamr = LAMMAXR
                n0rr = lamr**4*qr(k)/(PI*RHOW)
                nr(k) = n0rr/lamr
             endif
             vtrn(k) = ARN*(gamma(1.+BR))/LAMR**BR
             vtrm(k) = ARN*(gamma(4.+BR))/6./LAMR**BR
             vtrn(k) = min(vtrn(k),9.1)
             vtrm(k) = min(vtrm(k),9.1)

             dbzr = rhoair(k)*nr(k)/lamr**3/lamr**3*720.
             dbzr = 1.e18*dbzr

          else
  !.. If only small mass of rain, push to vapor
             qv(k)=qv(k)+qr(k)
             theta(k)=theta(k)+theta(k)*i_temp*(qr(k)*xxlv)*i_cp
             qr(k)   =0.
             nr(k)   =0.
             vtrm(k) =0.
             vtrn(k) =0.
          endif

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
          do cc = 1, cat
  !.. Check for ice after nucleation and then update density and bulk distribution sizes
             if(qi(cc,k).ge.QSMALL) then
                ni(cc,k) =max(ni(cc,k),QNSMALL)
                nucfrac = totalnuc(cc)/qi(cc,k)
                nucfrac = min(nucfrac, 1.0)
                nucfrac = max(nucfrac, 0.0)

  !.. Do not nucleate to aggregates
                if(cc.eq.IICE3) then
                   nucfrac=0.
                endif

  !.. PUT IN NEW METHODOLOGY FOR NUCLEATION AND AGGREGATION

                
  !.. Update density from nucleation
                rhobar(cc) = RHOI*nucfrac + (1.-nucfrac)*rhobar(cc)
                
  !.. Keep everything in bounds
                ai(cc,k)  =max(ai(cc,k),QASMALL)
                ci(cc,k)  =max(ci(cc,k),QASMALL)
                nibnuc(cc)=max(nibnuc(cc),QNSMALL)

                ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*nibnuc(cc)))**0.333333333333),2.e-6)
                cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*nibnuc(cc)))**0.333333333333),2.e-6)
                
                ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)                                      
                ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)

                aniold(cc)  = ani(cc)
                cniold(cc)  =ao**(1.-dsold(cc))*aniold(cc)**dsold(cc)

  !.. Get the size distribution properties of nucleated particles
                if(nucfrac.gt.1.e-5.and.totalnucn(cc).gt.0.) then
                   gamma_arg = NU+3.
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   gamma_arg = NU+4.
                   gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   anuc=(((totalnuc(cc)*gammnu)/(RHOI*totalnucn(cc)*fourthirdspi* &
                        gamma_tab(gi)))**0.3333333)*gamma_tab(gi2)/gamma_tab(gi)
                   anuc = MAX(anuc,2.e-6)
                   anuc = MIN(anuc,3.e-3)

                    ! Shape_diag takes in the new volumes, mass, and number and calculates deltastar and ani (cni)

                    ! Calculate new volumes
                    ! cnuc =
                    ! anuc = cnuc/phinuc
                    ! Vanuc = fourthirdspi*anuc**3*phinuc*totalnucn(cc)

                    !  call shape_diag(NU,gammnu,ni(cc,k),qi(cc,k),Viadum,Vicdum,ao,deltastr(cc),alphstr,alphv,betam,&
                    !  ai(cc,k),ci(cc,k),ani,cni)
                   
                   gamma_arg = NU+2.+deltastr(cc)
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   gamma_arg = NU+3.+deltastr(cc)
                   gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   gamma_arg = NU+2.+(2.*deltastr(cc))
                   gi3=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)

                   amass=ani(cc)*gamma_tab(gi2)/gamma_tab(gi)
                   cmass=cni(cc)*gamma_tab(gi3)/gamma_tab(gi)
                   
  !.. Mass-weighted size of old distribution and nucleated distribution
                   amass = amass*(1.-nucfrac) + anuc*nucfrac
                   cmass = cmass*(1.-nucfrac) + anuc*nucfrac
                   
  !.. New ani and cni after nucleation
                   ani(cc) = max((amass*gamma_tab(gi)/gamma_tab(gi2)),2.e-6)
                   cni(cc) = max((cmass*gamma_tab(gi)/gamma_tab(gi3)),2.e-6)
     
  !.. New delta* (shape) after nucleation
                   if(ani(cc).gt.(1.1*ao).and.cni(cc).gt.(1.1*ao)) then
                      deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))
                   else
                      deltastr(cc)=1.
                   endif

  !.. Nucleation should only make the average distribution more spherical
                   if(dsold(cc).le.1.) then
                      if(deltastr(cc).lt.dsold(cc)) then
                         deltastr(cc)=dsold(cc)
                      endif
                      if(deltastr(cc).gt.1.) then
                         deltastr(cc)=1.
                      endif
                   else
                      if(deltastr(cc).gt.dsold(cc)) then
                         deltastr(cc)=dsold(cc)
                      endif
                      if(deltastr(cc).lt.1.) then
                         deltastr(cc)=1.
                      endif
                   endif






                endif    !.. End nucleation update distribution properties

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. This needs to be done because of transfer from e.g.
  !.. IICE1 to IICE2 during nucleation.
  !.. Assume constant density/delta*, loss of mass and number

  !..With the new delta* and density after nucleation re-diagnose a and c           
  !..Use the old density and deltastr for loss values e.g. q1 = q1 - (transfered to q2)

                if(cc.ne.IICE3) then

                alphstr=ao**(1.-deltastr(cc))
                alphv=fourthirdspi*alphstr
                betam=2.+deltastr(cc)
                gamma_arg = NU+2.+deltastr(cc)
                gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                ani(cc)=((qi(cc,k)*gammnu)/(rhobar(cc)*ni(cc,k)*alphv* &
                     gamma_tab(gi)))**(1./betam)
                cni(cc)=ao**(1.-deltastr(cc))*ani(cc)**deltastr(cc)

  !.. volume and volume*aspect ratio process rates
                ai(cc,k) = ai(cc,k) + cniold(cc)*nibnuc(cc)*2.*aniold(cc)*(ani(cc)-aniold(cc)) + &
                     nibnuc(cc)*aniold(cc)**2*(cni(cc)-cniold(cc)) + &
                     aniold(cc)**2*cniold(cc)*(ni(cc,k)-nibnuc(cc))

                ci(cc,k) = ci(cc,k) + aniold(cc)*nibnuc(cc)*2.*cniold(cc)*(cni(cc)-cniold(cc)) + &
                     nibnuc(cc)*cniold(cc)**2*(ani(cc)-aniold(cc)) + &
                     cniold(cc)**2*aniold(cc)*(ni(cc,k)-nibnuc(cc))

                endif

  !.. Keep everything in bounds
                ni(cc,k)=max(ni(cc,k),QNSMALL)
                ai(cc,k)=max(ai(cc,k),QASMALL)
                ci(cc,k)=max(ci(cc,k),QASMALL)

                if (cc.ne.IICE3) then
                if(ani(cc).gt.(1.1*ao).and.cni(cc).gt.(1.1*ao)) then
                   deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao)) 
                else
                   deltastr(cc)=1.
                endif
                
                call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                     ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                     alphstr, alphv, betam)

                dsold(cc) = deltastr(cc)

                elseif (cc.eq.IICE3) then


                call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))


                endif

  !.. Limit the fall speeds here to 25 m s^-1
                vtrzi1(cc,k)=min(vtrzi1(cc,k),25.)
                vtrmi1(cc,k)=min(vtrmi1(cc,k),25.)
                vtrni1(cc,k)=min(vtrni1(cc,k),25.)
             else
             
  !.. If there is no ice even after nucleation, 
  !.. add remainder to vapor and update temperature
                qv(k)=qv(k)+qi(cc,k)
                theta(k)=theta(k)+theta(k)*i_temp*(qi(cc,k)*xxls)*i_cp
                qi(cc,k)=0.
                ai(cc,k)=0.
                ni(cc,k)=0.
                if (cc.eq.IICE3) then
                nmon(k)=0.
                endif
                ci(cc,k)=0.
                vtrzi1(cc,k)=0.
                vtrmi1(cc,k)=0.
                vtrni1(cc,k)=0.
             endif  !.. End check for ice
          enddo     !.. Loop over all ice species

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
          if(temp.le.T0.and.AGGREGATION) then
  !.. Aggregation occurs similarly to nulceation
  !.. We assume that the loss of IICE1 and IICE2 to aggregates does not
  !.. change the shape and density of the remaining particles.

             ! Transfer
             nmon(k) = nmon(k) + nmon3*dt
             nmon(k) = max(nmon(k),2.*QNSMALL)

             do cc = 1, cat

                if (cc.ne.IICE3) then
                if(ani(cc).gt.(1.1*ao).and.cni(cc).gt.(1.1*ao)) then
                deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))
                else
                deltastr(cc)=1.
                endif

                call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                alphstr, alphv, betam)

                dsold(cc) = deltastr(cc)

                alph_v(cc) = fourthirdspi*alphstr
                bet_m(cc) = 2.+deltastr(cc)
                alphc(cc) = fourthirdspi*alphstr**2
                betac(cc) = 1.+2.*deltastr(cc)


                Viadum(cc) = alph_v(cc)*ani(cc)**betam*ni(cc,k)*gamma(NU+betam)/gammnu
                Vicdum(cc) = alphc(cc)*ani(cc)**betac(cc)*ni(cc,k)*gamma(NU+betac(cc))/gammnu

                endif

                nibnuc(cc) = max(ni(cc,k),QNSMALL)

                qi(cc,k) = qi(cc,k) + qagg(cc)*dt
                ni(cc,k) = ni(cc,k) + nagg(cc)*dt
                Viadum(cc) = Viadum(cc)+vagg(cc)*dt
                Vicdum(cc) = Vicdum(cc)+vcagg(cc)*dt

                if (cc.ne.IICE3.and.ni(cc,k).gt.QNSMALL.and.qi(cc,k).gt.QSMALL) then



                !if (ani(cc).gt.2.e-3) then
                !print*, 'ds_before=',deltastr(cc)
                !print*, 'ni_before=',ni(cc,k)
                !print*, 'qi_before=',qi(cc,k)
                !print*, 'ani_before=',ani(cc)
                !print*, 'rhobar_before=',rhobar(cc)
                !print*, 'Viadum_before=',Viadum(cc)
                !print*, 'Vicdum_before=',Vicdum(cc)
                !endif

                ! New diagnose of shape (ai, ci, ani, cni)
                call shape_diag(gammnu,ni(cc,k),qi(cc,k),Viadum(cc),Vicdum(cc),ao,deltastr(cc),&
                                alphstr,alphv,betam,ai(cc,k),ci(cc,k),ani(cc),cni(cc))

                !call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                !ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
               ! alphstr, alphv, betam)


                !alphstr=ao**(1.-deltastr(cc))
                !alphv = fourthirdspi*alphstr
                !betam = 2.+deltastr(cc)

                ! Hold deltastar constant but diagnose new ani
                !ani(cc) = ((Viadum(cc)**2/(fourthirdspi*Vicdum(cc)*ni(cc,k)))*&
                !(gamma(nu+2.*deltastr(cc)+1.)*gammnu/gamma(nu+2.+deltastr(cc))**2))**(1./3.)


                !rhobar(cc) = (qi(cc,k)*gammnu/(alphv*ani(cc)**betam*ni(cc,k)*gamma(nu+betam)))


                !if (ani(cc).gt.2.e-3) then
                !print*, 'ds_after=',deltastr(cc)
                !print*, 'ani_after=',ani(cc)
                !print*, 'rhobar_after=',rhobar(cc)
                !print*, 'Viadum_after=',Viadum(cc)
                !print*, 'Vicdum_after=',Vicdum(cc)
                !stop
                !endif

                endif

               ! ai(cc,k) = ai(cc,k) + vagg(cc)*dt
               ! ci(cc,k) = ci(cc,k) + vcagg(cc)*dt
                if(qi(cc,k).gt.QSMALL) then
                   
  !.. Aggregate (unmelted) density and aspect ratio are assumed
  !.. Force the shape of aggregates
                   if(cc.eq.IICE3) then
                      rhobar(cc) = rhobar_ellip



                      call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))


                   else

  !.. Variable check on ice-1 and ice-2 after aggregation

                      !ni(cc,k)=max(ni(cc,k),QNSMALL)
                      !ai(cc,k)=max(ai(cc,k),QASMALL)
                      !ci(cc,k)=max(ci(cc,k),QASMALL)
                      !ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*nibnuc(cc)))**0.333333333333),2.e-6)
                      !cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*nibnuc(cc)))**0.333333333333),2.e-6)

                      !aniold(cc)  = ani(cc)
                      !cniold(cc)  = cni(cc)
                      !deltastr(cc)= dsold(cc)
                
  !.. Update IICE1 and IICE2 sizes after aggregation from
  !.. updated mass and number mixing ratios (transfer to aggregates) and an
  !.. assumed constant deltastar
                      !alphstr=ao**(1.-deltastr(cc))
                      !gamma_arg = NU+2.+deltastr(cc)
                      !gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                      !ani(cc) = ((qi(cc,k))/(ni(cc,k)*rhobar(cc)*fourthirdspi*alphstr*gamma_tab(gi)*i_gammnu))**&
                      !     (1./(2.+deltastr(cc)))
                      !ani(cc)=max(ani(cc),2.e-6)
                      !cni(cc)=alphstr*ani(cc)**deltastr(cc)

  !.. volume and volume*aspect ratio updates
                      !ai(cc,k) = ai(cc,k) + cniold(cc)*nibnuc(cc)*2.*aniold(cc)*(ani(cc)-aniold(cc)) + &
                      !     nibnuc(cc)*aniold(cc)**2*(cni(cc)-cniold(cc)) + &
                      !     aniold(cc)**2*cniold(cc)*(ni(cc,k)-nibnuc(cc))
                      
                      !ci(cc,k) = ci(cc,k) + aniold(cc)*nibnuc(cc)*2.*cniold(cc)*(cni(cc)-cniold(cc)) + &
                      !     nibnuc(cc)*cniold(cc)**2*(ani(cc)-aniold(cc)) + &
                      !     cniold(cc)**2*aniold(cc)*(ni(cc,k)-nibnuc(cc))
                      
  !.. Keep everything in bounds
                      ni(cc,k)=max(ni(cc,k),QNSMALL)
                      ai(cc,k)=max(ai(cc,k),QASMALL)
                      ci(cc,k)=max(ci(cc,k),QASMALL)

                      !if(ani(cc).gt.(1.1*ao).and.cni(cc).gt.(1.1*ao)) then
                     !    deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))
                     ! else
                     !    deltastr(cc)=1.
                     ! endif


                      !.. Aggregation should only make the average distribution more spherical
                      if(dsold(cc).le.1.) then
                      if(deltastr(cc).lt.dsold(cc)) then
                         deltastr(cc)=dsold(cc)
                      endif
                      if(deltastr(cc).gt.1.) then
                         deltastr(cc)=1.
                      endif
                      else
                      if(deltastr(cc).gt.dsold(cc)) then
                         deltastr(cc)=dsold(cc)
                      endif
                      if(deltastr(cc).lt.1.) then
                         deltastr(cc)=1.
                      endif
                      endif

                      if (cc.ne.IICE3) then
                      call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                           ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                           alphstr, alphv, betam)


                      endif

                   endif
                endif  !.. If ice exist
             enddo     !.. Loop over all species



          endif        !.. If temperature < 273.15 K and Aggregation turned on.
          
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Update the thermodynamics
          theta(k)=theta(k)+theta(k)*i_temp*((mnuccd)*xxls*i_cp*dt)
          do cc = 1, cat
             theta(k)=theta(k)+theta(k)*i_temp*((prd(cc))*xxls*i_cp*dt)
             theta(k)=theta(k)+theta(k)*i_temp*((qmlt(cc))*xxlf*i_cp*dt)
          enddo

          theta(k)=theta(k)+theta(k)*i_temp*(pre*xxlv)*i_cp*dt
          
          if(temp.le.T0) then
             theta(k)=theta(k)+theta(k)*i_temp*(mim+mimr+mbiggr)*xxlf*i_cp*dt

             do cc = 1, cat
                theta(k)=theta(k)+theta(k)*i_temp*(prdr(cc)+dQRfzri(cc))*xxlf*i_cp*dt
             enddo
          endif

          
          temp  = theta(k)/(100000./pres_e(k))**(RCP)
          i_temp= 1./temp
          qvs   = 0.622*polysvp(temp,0)/(pres_e(k)-polysvp(temp,0))
          xxlv  = 3.1484E6-2370.*temp         
          rhoair(k)= pres_e(k)/(RD*temp)
          i_rhoair(k) = 1./rhoair(k)

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Saturation adjustment for qc
          pcc = (qv(k)-qvs) / (1.0 + xxlv**2*qvs/(CP*rv*temp**2))*i_dt
          if (pcc*dt+qc(k).lt.0.) then
             pcc=-qc(k)/dt
          end if
          theta(k)=theta(k) + theta(k)*i_temp*pcc*xxlv*i_cp*dt
          qv(k) = qv(k) - pcc*dt
          qc(k) = qc(k) + pcc*dt
          
          if(qc(k).le.QSMALL) then
             qv(k)=qv(k)+qc(k)
             theta(k)=theta(k)+theta(k)*i_temp*(qc(k)*xxlv)*i_cp
             qc(k)=0.
             nc(k)=0.
          endif
       
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Cloud water effective radius
          lrsig = 1.3
          sig = log(lrsig)         !.. Distribution standard deviation                                      
          lwc = qc(k)*rhoair(k)    !.. liquid water content
          ncm3dum = nc(k)*rhoair(k)
          
          if(ncm3dum .gt. 0.) then
             r_n = ((lwc)/(fourthirdspi*ncm3dum*RHOW*exp(4.5*(sig**2))))**(0.33333333333)
          else
             r_n = 0.
          endif

          effc1d(k) = r_n * (5.*EXP((sig**2)/2.))
          effc1d(k) = MIN(effc1d(k),50.e-6)
          effc1d(k) = MAX(effc1d(k),2.51e-6)
          
          ! For now, assume cloud droplet fallspeed is zero
          vtrmc(k) = 0.0
          vtrnc(k) = 0.0

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Final check on ice before sedimentation
          do cc = 1, cat

             if(qi(cc,k).gt.QSMALL) then
                ni(cc,k)=max(ni(cc,k),QNSMALL)
                if(cc.eq.IICE3) then
                nmon(k)=max(nmon(k),2.*QNSMALL)
                endif
                ai(cc,k)=max(ai(cc,k),QASMALL)
                ci(cc,k)=max(ci(cc,k),QASMALL)

                if (cc.ne.IICE3) then
                ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
                cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
                
                ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)                                              
                ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)

                deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao)) 
                
                call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
                     ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
                     alphstr, alphv, betam)

                endif
                
  !.. Final check on aggregates
                if(cc.eq.IICE3) then
                      rhobar(cc) = rhobar_ellip

                      call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))

                   
                endif ! Check on aggregates



  !.. Fall speed moved to here
  !.. See Harrington et al. 2013 and Mitchell and Heymsfield 2005
                fv(cc) = 1.0
                fh(cc) = 1.0

            if (cc.ne.IICE3) then
                
                gamma_arg = NU-1.+deltastr(cc)
                gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                phiivt = cni(cc)/ani(cc)*gamma_tab(gi)*i_gammnu
                alphstr=ao**(1.-deltastr(cc))
                   
  !.. Determine  coefficeients for the Best number for fall speeds
                if(phiivt.lt.1.0)then
                   bl = 1.
                   al = 2.
                   aa = PI
                   ba = 2.
                   qe = (1.-phiivt)*(rhobar(cc)/RHOI) + phiivt
                else if(phiivt .gt. 1.0)then
                   al = 2.0
                   bl = 1.
                   aa = PI*alphstr
                   ba = deltastr(cc) + 1.
                  qe = 1.0
                else if(phiivt .eq. 1.0)then
                   bl = 1.
                   al = 2.
                   aa = PI
                   ba = 2.
                   qe = 1.0
                endif
                qe = min(qe, 1.0)
                
  !.. Fall speed and ventilation (Best number formulation)
                xn =  2./rhobar(cc)*(rhobar(cc)-rhoair(k))*G_HOME*rhoair(k)/mu**2 * &
                     (fourthirdspi*rhobar(cc))*alphstr*al**2/(aa*qe) * qe**(3./4.)
                bx = deltastr(cc)+2.+2.*bl-ba
                
  !.. Number-average Best Number
                xm = xn*ani(cc)**bx * (gamma(NU+bx))*i_gammnu
                
  !.. The following fall speed coefficients are from 
  !.. Mitchell and Heymsfiled 2005
                f_c1 = 4.0 / (5.83 * 5.83 * SQRT(0.6))
                f_c2 = (5.83 * 5.83) / 4.0
                bm = ((f_c1 * SQRT(xm)) / &
                     (2.0 * (SQRT(1.0 + f_c1 * SQRT(xm)) - 1.0) * &
                     (SQRT(1.0 + f_c1 * SQRT(xm))))) - &
                     (1.0e-5 * xm) / &
                     (f_c2 * (SQRT(1.0 + f_c1 * SQRT(xm)) - 1.0)**2.0)
                am = ((f_c2 * (SQRT(1.0 + f_c1 * SQRT(xm)) - 1.0)**2.0) - &
                     (1.0e-5 * xm)) / (xm**bm)
                
                if(xm.gt.1.e8)then
                   am = 1.0865
                   bm = 0.499
                endif
                
  !.. Reynolds Number
                Nre = am*xm**bm
                
  !.. Number-averaged ice fall speed
                vtrni1(cc,k) = min((mu/rhoair(k)*0.5 * am*(xn)**bm*ani(cc)**(bx*bm-1.) * &
                     (gamma(NU+bx*bm-1.))*i_gammnu),25.)
                
  !.. Mass-averaged ice fall speed
                vtrmi1(cc,k) = min((mu/rhoair(k)*0.5 * am*(xn)**bm*ani(cc)**(bx*bm-1.) * &
                     (gamma(NU+bx*bm-1.+2.+deltastr(cc)))/(gamma(NU+2.+deltastr(cc)))),25.)

                elseif (cc.eq.IICE3) then ! Aggregate fallspeed calc for sedimentation

                alphm_agg = (qi(cc,k)/ni(cc,k))*(1./phiba_phica_bar)*ani(cc)**(-betm_agg)*gamma(NU_AGG)/gamma(NU_AGG+betm_agg)

                   ! Coefficients for projected area from Mitchell 1996
                alphx_agg = 8.*alphm_agg*alph_ar_agg**(-kc)*G_HOME*rhoair(k)/(PI*mu**2)

                !Average Best Number
                xm = alphx_agg*ani(cc)**betx_agg * gamma(NU_AGG+betx_agg)/gamma(NU_AGG) * &
                pbeta_chi

                bm = (C1*sqrt(xm)/&
                (2.*(sqrt(1+C1*sqrt(xm)) - 1.)*sqrt(1.+C1*sqrt(xm)))) - &
                ((a0*b0*xm**b0)/(C2*(sqrt(1.+C1*sqrt(xm))-1.)**2))

                am = &
                (C2*((sqrt(1.+C1*sqrt(xm))-1.)**2) - a0*xm**b0)/(xm**bm)

                zetba1(cc,k) = chi_ba*bm-0.5*bet_ba
                zetca1(cc,k) = chi_ca*bm-0.5*bet_ca

                !Average Reynolds Number
                Nre = am*xm**bm

                alphvt_agg = mu/rhoair(k)*0.5*am*(alphx_agg)**bm
                betvt1(cc,k) = betx_agg*bm-1.

                vtn = alphvt_agg*ani(cc)**betvt1(cc,k)


            !    print*, 'an=',ani(cc)
            !    print*, 'alphx_agg=',alphx_agg
            !    print*, 'alphm_agg=',alphm_agg
            !    print*, 'C1=',C1
            !    print*, 'C2=',C2
            !    print*, 'bm=',bm
            !    print*, 'xm=',xm
            !    print*, 'pbeta_chi=',pbeta_chi
            !    print*, 'chi_ba=',chi_ba
            !    print*, 'chi_ca=',chi_ca
            !    print*, 'zetaba=',zetba1(cc,k)
            !    print*, 'zetaca=',zetca1(cc,k)


                vtrni1(cc,k) = vtn*gamma(NU_AGG+betvt1(cc,k))/gamma(NU_AGG) * &
                beta_moms(zetba1(cc,k),zetca1(cc,k),aba,bba,bcb)

                vtrmi1(cc,k) = vtn*gamma(NU_AGG+betvt1(cc,k)+betm_agg)/gamma(NU_AGG+betm_agg) * &
                beta_moms(zetba1(cc,k)+1.,zetca1(cc,k)+1.,aba,bba,bcb)/phiba_phica_bar


                endif
                
  !.. Remove size sorting when melting
                if(temp.gt.T0.and.qmlt(cc).lt.0.) then
                   vtrni1(cc,k) = vtrmi1(cc,k)
                endif

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Determine the ice radiative properties 
  !.. For simlicity this is 
  !.. Planar ice:   integral of a^2c / integral of a^2
  !.. Columnar ice: integral of a^2c / integral of a*c

                if (cc.ne.IICE3) then

                gamma_arg = NU-1.+deltastr(cc)
                gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                phirad = cni(cc)/ani(cc)*gamma_tab(gi)*i_gammnu                   
                if(phirad.le.1.) then 
                   gamma_arg = NU+2.+deltastr(cc)
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   gamma_arg = NU+2.
                   gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   effi(cc,k) = cni(cc)*gamma_tab(gi)/gamma_tab(gi2)
                else
                   gamma_arg = NU+2.+deltastr(cc)
                   gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   gamma_arg = NU+1.+deltastr(cc)
                   gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
                   effi(cc,k) = ani(cc)*gamma_tab(gi)/gamma_tab(gi2)
                endif
                
  !.. For advection use volume and volume times aspect ratio mixing ratios
                ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)
                ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)

                endif

                ai(cc,k)=max(ai(cc,k),QASMALL)
                ci(cc,k)=max(ci(cc,k),QASMALL)
                ni(cc,k)=max(ni(cc,k),QNSMALL)
                if(cc.eq.IICE3) then
                nmon(k)=max(nmon(k),2.*QNSMALL)
                endif

             else
  !.. If not enough ice convert to vapor
                qv(k)=qv(k)+qi(cc,k)
                theta(k)=theta(k)+theta(k)*i_temp*(qi(cc,k)*xxls)*i_cp
                qi(cc,k)=0.
                ai(cc,k)=0.
                ni(cc,k)=0.
                if(cc.eq.IICE3) then
                nmon(k)=0.
                endif
                ci(cc,k)=0.
                vtrzi1(cc,k)=0.
                vtrmi1(cc,k)=0.
                vtrni1(cc,k)=0.
             endif !.. If ice exists

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Output diagnostics
             gamma_arg = NU-1.+deltastr(cc)
             gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
             
  !.. For aspect ratio diagnostic set lower size limit
             if(ai(cc,k).gt.1.e-12.and.ci(cc,k).gt.1.e-12.and. &
                  ani(cc).gt.2.e-6.and.cni(cc).gt.2.e-6.and.qi(cc,k).gt.1.e-9) then
                phii1d(cc,k) = cni(cc)/ani(cc)*gamma_tab(gi)*i_gammnu
                rhopo1d(cc,k) = rhobar(cc)
             endif

             if (cc.eq.IICE3.and.qi(cc,k).gt.1.e-9.and.ani(cc).gt.2.e-6) then

             rhopo1d(cc,k) = rhobar(cc)
             di1d(cc,k) = 2.*NU_AGG*ani(cc) ! Currently number-weighted.
             vmi1d(cc,k) = vtrmi1(cc,k)

             endif

            ! LEE NOTE: probably can reduce these calculations using gamma function properties.
             
             gamma_arg = NU+3.+deltastr(cc)
             gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
             gamma_arg = NU+2.+deltastr(cc)
             gi2=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
             gamma_arg = NU+2.+(2.*deltastr(cc))
             gi3=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)

             if(cc.ne.IICE3.and.ai(cc,k).gt.1.e-12.and.ci(cc,k).gt.1.e-12.and. &
                  ani(cc).gt.2.e-6.and.cni(cc).gt.2.e-6.and.qi(cc,k).gt.1.e-9) then
                di1d(cc,k)  = MAX((2.*ani(cc)*gamma_tab(gi)/gamma_tab(gi2)), &
                     (2.*cni(cc)*gamma_tab(gi3)/gamma_tab(gi2)))
                vmi1d(cc,k) = vtrmi1(cc,k)
             endif

             if(ni(cc,k).gt.0..and.ani(cc).gt.0.) then
                if (cc.ne.IICE3) then
                vtrzi1(cc,k) = 0.176/0.93*(6./PI)**2*((fourthirdspi*rhobar(cc)* &
                     (ao)**(1.-deltastr(cc)))/(2.**(2.+deltastr(cc))))**2 * &
                     rhoair(k)/900./900.*ni(cc,k)*(2.*ani(cc))**(2.*(2.+deltastr(cc))) * &
                     gamma(NU+2.*(2.+deltastr(cc)))*i_gammnu
                elseif (cc.eq.IICE3) then

                vtrzi1(cc,k) = 0.176/0.93*(6./PI)**2*((fourthirdspi*rhobar(cc)/&
                (2.**3)))**2 * &
                rhoair(k)/900./900.*ni(cc,k)*(2.*ani(cc))**(6) * &
                gamma(NU_AGG+4.)/gamma(NU_AGG)*pb2

                endif
                vtrzi1(cc,k) = 1.e18*vtrzi1(cc,k)
             else
                vtrzi1(cc,k) = 0.
             endif
          enddo    !.. Loop over all ice species
          
          dbzsum = dbzr + vtrzi1(IICE1,k) + vtrzi1(IICE2,k) + vtrzi1(IICE3,k)
          if(dbzsum.ne.0.) then
             dbz1d(k) = 10.*Log10(dbzsum)
             dbz1d(k) = MAX(-35.,dbz1d(k))
          endif

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Ice effective radius is a mass-weighting of all ice species
          totaliwc = 0.
          do cc = 1, cat
             totaliwc = totaliwc + qi(cc,k)
          enddo
          if(totaliwc.gt.0.) then
             radw1 = MAX((MIN((qi(IICE1,k)/totaliwc),1.)),0.)
             radw2 = MAX((MIN((qi(IICE2,k)/totaliwc),1.)),0.)
             radw3 = MAX((MIN((qi(IICE3,k)/totaliwc),1.)),0.)
             effi1d(k) = MAX((MIN((effi(IICE1,k)*radw1 + effi(IICE2,k)*radw2 + &
                  effi(IICE3,k)*radw3),999.e-6)),5.e-6)
          else
             effi1d(k)=0.
          endif
          
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Check for rain or ice to do sedimentation
          if(qr(k).gt.QSMALL) then
             sedi=.true.
          endif
          do cc = 1, cat
             if(qi(cc,k).gt.QSMALL) then
                sedi=.true.
             endif
          enddo
          
  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. else on domicro
       else 
  
  !.. If we don't do all of the above
  !.. Push to vapor and get density for seimentation in the column

          xxls  = 3.15e6-2370.*temp+0.3337e6
          xxlv  = 3.1484E6-2370.*temp         
          
          i_temp= 1./temp
          rhoair(k)= pres_e(k)/(RD*temp)
          i_rhoair(k) = 1./rhoair(k)
          
          qv(k)=qv(k)+qc(k)+qr(k)
          theta(k)=theta(k)+theta(k)*i_temp*((qc(k)+qr(k))*xxlv)*i_cp
          do cc = 1, cat
             qv(k)=qv(k)+qi(cc,k)
             theta(k)=theta(k)+theta(k)*i_temp*((qi(cc,k))*xxls)*i_cp
             qi(cc,k)=0.
             ni(cc,k)=0.
             ai(cc,k)=0.
             ci(cc,k)=0.
             vtrmi1(cc,k)=0.
             vtrni1(cc,k)=0.
             vtrzi1(cc,k)=0.
             vmi1d(cc,k)=0.
             di1d(cc,k)=0.
             phii1d(cc,k)=0.
             rhopo1d(cc,k)=0.
             qipre1d(cc)  =0.
             icetype1d(cc,k)=0.
          enddo

          effc1d(k)=0.
          effi1d(k)=0.
          dbz1d(k)=-35.
          qc(k)=0.
          qr(k)=0.
          nr(k)=0.
          vtrm(k)=0.
          vtrn(k)=0.
          qrpre1d=0.
          nmon(k)=0.

       endif  !.. If domicro
    enddo     !.. do k = kts, kte

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Sedimentation (Reisner et al. 1998)
  !.. Fallout terms are calculated on split time steps to ensure numerical stability (Courant < 1)

    nstep = 1
    maxfall=0.
    
    if(sedi) then
       do k = kte,kts,-1
          
  !.. Calculate number of split timesteps
          maxfall = MAX(vtrmi1(IICE1,k),vtrni1(IICE1,k),vtrmi1(IICE2,k),vtrni1(IICE2,k), &
               vtrmi1(IICE3,k),vtrni1(IICE3,k),vtrm(k),vtrn(k))
          nstep = MAX(INT(maxfall*dt/dzmic(k)+1.),nstep)
          
          do cc = 1, cat
             qi(cc,k) = qi(cc,k)*rhoair(k)
             ni(cc,k) = ni(cc,k)*rhoair(k)
             ai(cc,k) = ai(cc,k)*rhoair(k)
             ci(cc,k) = ci(cc,k)*rhoair(k)
          enddo

          nmon(k) = nmon(k)*rhoair(k)
          qr(k) = qr(k)*rhoair(k)
          nr(k) = nr(k)*rhoair(k)



       enddo
       
       do nn = 1, nstep
          do k = kts, kte
             
             do cc = 1, cat
                fluxqi(cc,k) = vtrmi1(cc,k)*qi(cc,k)
                fluxni(cc,k) = vtrni1(cc,k)*ni(cc,k)
                fluxai(cc,k) = vtrmi1(cc,k)*ai(cc,k)
                fluxci(cc,k) = vtrmi1(cc,k)*ci(cc,k)
                
             enddo
             fluxnmon(k) = vtrni1(IICE3,k)*nmon(k)
             fluxqr(k) = vtrm(k)*qr(k)
             fluxnr(k) = vtrn(k)*nr(k)
          enddo
          
          do cc = 1, cat
             falltndqi(cc) = fluxqi(cc,KTE)/dzmic(KTE)
             falltndni(cc) = fluxni(cc,KTE)/dzmic(KTE)
             falltndai(cc) = fluxai(cc,KTE)/dzmic(KTE)
             falltndci(cc) = fluxci(cc,KTE)/dzmic(KTE)
          enddo
          falltndmon = fluxnmon(KTE)/dzmic(KTE)
          falltndqr = fluxqr(KTE)/dzmic(KTE)
          falltndnr = fluxnr(KTE)/dzmic(KTE)
          
          do cc = 1, cat
             qi(cc,KTE) = qi(cc,KTE)-falltndqi(cc)*dt/nstep
             ni(cc,KTE) = ni(cc,KTE)-falltndni(cc)*dt/nstep
             ai(cc,KTE) = ai(cc,KTE)-falltndai(cc)*dt/nstep
             ci(cc,KTE) = ci(cc,KTE)-falltndci(cc)*dt/nstep
          enddo
          nmon(KTE) = nmon(KTE)-falltndmon*dt/nstep
          qr(KTE) = qr(KTE)-falltndqr*dt/nstep
          nr(KTE) = nr(KTE)-falltndnr*dt/nstep
          
          
          do k = kte-1,kts,-1
             do cc = 1, cat
                falltndqi(cc) = (fluxqi(cc,k+1)-fluxqi(cc,k))/dzmic(k)
                falltndni(cc) = (fluxni(cc,k+1)-fluxni(cc,k))/dzmic(k)
                falltndai(cc) = (fluxai(cc,k+1)-fluxai(cc,k))/dzmic(k)
                falltndci(cc) = (fluxci(cc,k+1)-fluxci(cc,k))/dzmic(k)
              enddo

             falltndmon = (fluxnmon(k+1)-fluxnmon(k))/dzmic(k)
             falltndqr = (fluxqr(k+1)-fluxqr(k))/dzmic(k)
             falltndnr = (fluxnr(k+1)-fluxnr(k))/dzmic(k)
             
             do cc = 1, cat
                qi(cc,k) = qi(cc,k)+falltndqi(cc)*dt/nstep
                ni(cc,k) = ni(cc,k)+falltndni(cc)*dt/nstep
                ai(cc,k) = ai(cc,k)+falltndai(cc)*dt/nstep
                ci(cc,k) = ci(cc,k)+falltndci(cc)*dt/nstep

  !.. GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
  !.. FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
  !.. OF LIQUID WATER CANCELS THIS FACTOR OF 1000
                
  !.. Ice precipitation rate
                if(k.eq.kts.and.qi(cc,k).gt.1.e-9) then
                   qipre1d(cc) = qipre1d(cc)+fluxqi(cc,k)*dt/nstep
                endif
                
             enddo
             nmon(k) = nmon(k)+falltndmon*dt/nstep
             qr(k) = qr(k)+falltndqr*dt/nstep
             nr(k) = nr(k)+falltndnr*dt/nstep
             
  !.. Rain precipitation rate
             if(k.eq.kts.and.qr(k).gt.1.e-9) then
                qrpre1d = qrpre1d + fluxqr(k)*dt/nstep
             endif
          enddo
       enddo
          
       do k = kts, kte
          do cc = 1, cat
             qi(cc,k) = qi(cc,k)*i_rhoair(k)
             ni(cc,k) = ni(cc,k)*i_rhoair(k)
             ai(cc,k) = ai(cc,k)*i_rhoair(k)
             ci(cc,k) = ci(cc,k)*i_rhoair(k)
          enddo
          nmon(k) = nmon(k)*i_rhoair(k)
          qr(k) = qr(k)*i_rhoair(k)
          nr(k) = nr(k)*i_rhoair(k)


       enddo


    endif

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Check after sedimentation
    do k = kts, kte
       temp  = theta(k)/(100000./pres_e(k))**(RCP)
       i_temp= 1./temp
       xxls  = 3.15e6-2370.*temp+0.3337e6
       xxlv  = 3.1484E6-2370.*temp         
       
       do cc = 1, CAT
          if(qi(cc,k).gt.QSMALL) then

  !.. Keep variables in bounds
             ni(cc,k)=max(ni(cc,k),QNSMALL)
             ai(cc,k)=max(ai(cc,k),QASMALL)
             ci(cc,k)=max(ci(cc,k),QASMALL)
             if (cc.ne.IICE3) then
             ani(cc) =max((((ai(cc,k)**2)/(ci(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)
             cni(cc) =max((((ci(cc,k)**2)/(ai(cc,k)*ni(cc,k)))**0.333333333333),2.e-6)

  !.. Set incoming limit after sedimentation
  !.. smallest values for which we care about shape              
             if(ai(cc,k).lt.1.e-12.or.ci(cc,k).lt.1.e-12) then
                ai(cc,k)=min(ai(cc,k),ci(cc,k))
                ci(cc,k)=ai(cc,k)
                ani(cc)=(ai(cc,k)/ni(cc,k))**0.3333333333
                cni(cc)=(ci(cc,k)/ni(cc,k))**0.3333333333
             endif
             
             if(ani(cc).lt.2.e-6.or.cni(cc).lt.2.e-6) then
                ani(cc)=2.e-6
                cni(cc)=2.e-6
             endif

             ci(cc,k) = cni(cc)**2*ani(cc)*ni(cc,k)                                           
             ai(cc,k) = ani(cc)**2*cni(cc)*ni(cc,k)
             
             deltastr(cc)=(log(cni(cc))-log(ao))/(log(ani(cc))-log(ao))

            call var_check(fourthirdspi, gammnu, qi(cc,k), deltastr(cc),       &
            ani(cc), cni(cc), rni(cc), rhobar(cc), ni(cc,k), ai(cc,k), ci(cc,k),  &
            alphstr, alphv, betam)
             
  !.. Final check on aggregates
             elseif(cc.eq.IICE3) then
                nmon(k)=max(nmon(k),2.*QNSMALL)
                if(temp.le.T0) then
                   rhobar(cc) = rhobar_ellip
                endif

                deltastr(cc) = 1.0
                
  !.. Final check on aggregate size (implicit breakup)

              call var_check_agg(fourthirdspi,qi(cc,k),ni(cc,k),ani(cc),rhobar(cc))



             endif !.. Check on aggregates



          else
             qv(k)=qv(k)+qi(cc,k)
             theta(k)=theta(k)+theta(k)*i_temp*(qi(cc,k)*xxls)*i_cp
             qi(cc,k)=0.
             ni(cc,k)=0.
             if (cc.eq.IICE3) then
             nmon(k)=0.
             endif
             ai(cc,k)=0.
             ci(cc,k)=0.
          endif
       enddo
       
       if(qr(k).gt.QSMALL) then
          nr(k) = max(nr(k),QNSMALL)
          lamr = (PI*RHOW*nr(k)/qr(k))**0.333333333
          n0rr = nr(k)*lamr
          if(lamr.LT.LAMMINR) then
             lamr = LAMMINR
             n0rr = lamr**4*qr(k)/(PI*RHOW)
             nr(k) = n0rr/lamr
          elseif(lamr.gt.LAMMAXR) then
             lamr = LAMMAXR
             n0rr = lamr**4*qr(k)/(PI*RHOW)
             nr(k) = n0rr/lamr
          endif
       else
          qv(k)=qv(k)+qr(k)
          theta(k)=theta(k)+theta(k)*i_temp*(qr(k)*xxlv)*i_cp
          qr(k)   =0.
          nr(k)   =0.
       endif

    enddo

  end subroutine me_ishmael

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. This subroutine keeps the ice properties in check

  subroutine var_check(fourthirdspi, gammnu, qidum, dsdum, ani,  &
       cni, rni, rbdum, nidum, aidum, cidum, alphstr, alphv, betam)
    
    implicit none
    
    REAL, INTENT(IN) ::  fourthirdspi, gammnu, qidum
    REAL :: voltmp, maxsize, gamma_arg
    INTEGER :: gi
    REAL, INTENT(INOUT) :: dsdum, ani, cni, rbdum, nidum, aidum, cidum
    REAL, INTENT(OUT) :: rni, alphstr, alphv, betam
    
  !.. Deltastr check
    if(dsdum.lt.0.55) then 
       dsdum=0.55
       ani=(cni/(ao**(1.-dsdum)))**(1./dsdum)
       aidum=ani**2*cni*nidum
       cidum=cni**2*ani*nidum
    else if (dsdum.gt.1.3) then
       dsdum=1.3
       cni=ao**(1.-dsdum)*ani**dsdum
       aidum=ani**2*cni*nidum
       cidum=cni**2*ani*nidum
    endif

    alphstr=ao**(1.-dsdum)
    alphv = fourthirdspi*alphstr
    betam = 2.0 + dsdum
    gamma_arg = NU+2.+dsdum
    gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)

  !.. Ice density check
  !.. Keep ice density between 50 and RHOI=920 kg m^-3
    if(ani.gt.2.e-6) then
       rbdum = qidum*gammnu/(nidum*alphv* &
            ani**betam*gamma_tab(gi))
    else
       rbdum = RHOI
    endif

    if(rbdum.gt.RHOI) then
       rbdum=RHOI
       ani=((qidum*gammnu)/(rbdum*nidum*alphv* & 
            gamma_tab(gi)))**(1./betam)
       cni=ao**(1.-dsdum)*ani**dsdum
       aidum=ani**2*cni*nidum
       cidum=cni**2*ani*nidum
    elseif(rbdum.lt.50.) then
       rbdum=50.
       ani=((qidum*gammnu)/(rbdum*nidum*alphv* &
            gamma_tab(gi)))**(1./betam)
       cni=ao**(1.-dsdum)*ani**dsdum
       aidum=ani**2*cni*nidum
       cidum=cni**2*ani*nidum
    endif
    
  !.. Small ice limit
  !.. Keep rni > 2 microns
    rni= (qidum*3./(nidum*rbdum*4.*PI* &
         (gamma_tab(gi)/gammnu)))**0.333333333333 
    if(rni.lt.2.e-6) then
       rni=2.e-6
       nidum=3.*qidum*gammnu/(4.*pi*rbdum*rni**3* &
            (gamma_tab(gi)))
       ani=((qidum*gammnu)/(rbdum*nidum*alphv* &
              gamma_tab(gi)))**(1./betam)
       cni=ao**(1.-dsdum)*ani**dsdum
       aidum=ani**2*cni*nidum
       cidum=cni**2*ani*nidum
    endif

  !.. Large ice limit
  !.. This is a number weighted diameter of 8 mm and 
  !.. a spherical mass weighted diameter of 14 mm 
    maxsize = max(ani,cni)
    if (maxsize.gt.1.e-3) then
       if(ani.ge.cni) then
          ani = 1.e-3
          nidum=qidum*gammnu/(fourthirdspi*rbdum*ao**(1.-dsdum)*ani**(2.+dsdum)* &
               (gamma_tab(gi)))
          cni=ao**(1.-dsdum)*ani**dsdum
          aidum=ani**2*cni*nidum
          cidum=cni**2*ani*nidum
       else
          cni = 1.e-3
          ani = (cni/(ao**(1.-dsdum)))**(1./dsdum)
          nidum=qidum*gammnu/(fourthirdspi*rbdum*ao**(1.-dsdum)*ani**(2.+dsdum)* &
               (gamma_tab(gi)))
          aidum=ani**2*cni*nidum
          cidum=cni**2*ani*nidum
       endif
       rni= (qidum/(nidum*rbdum*fourthirdspi* &
            (gamma_tab(gi)/gammnu)))**0.333333333333 
    end if
    
    return
  end subroutine var_check


!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
  subroutine var_check_agg(fourthirdspi,qidum,nidum,anagg,rhobar)
  !Make sure ice variables are in line

  implicit none

  REAL, INTENT(INOUT) :: anagg, nidum, qidum
  REAL, INTENT(IN)    :: fourthirdspi, rhobar
  real, parameter :: anmax = 1.e-2, anmin = 1.e-6
  real, parameter :: QSMALL = 1.e-12, QNSMALL= 1.25e-7  !.. Smallest ice mass (kg kg^-1)


  if(qidum.lt.0.) then
    qidum = QSMALL
  endif


 ! diagnose anagg
    anagg = ((qidum/nidum)/(fourthirdspi*rhobar*phiba_phica_bar)*&
              gamma(NU_AGG)/(gamma(NU_AGG+3.)))**(0.33333333)


!.. Size check

  if (anagg.gt.anmax) then
   anagg = anmax
   nidum = qidum/(fourthirdspi*anagg**3. * rhobar*&
           (gamma(nu_agg+3.)/gamma(nu_agg))*phiba_phica_bar)
  elseif(anagg.lt.anmin) then
   anagg = anmin
   nidum = qidum/(fourthirdspi*anagg**3. * rhobar*&
           (gamma(nu_agg+3.)/gamma(nu_agg))*phiba_phica_bar)


  !else

 !   anagg = ((qidum/nidum)/(fourthirdspi*rhobar*phiba_phica_bar)*&
 !             gamma(NU_AGG)/(gamma(NU_AGG+3.)))**(0.33333333)

  endif

 nidum = max(nidum,QNSMALL)

  return
  end subroutine var_check_agg



  !----------------------------------------------------------------------!

  !----------------------------------------------------------------------!

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. This is the ventilation, fall speed, vapor growth routine
  !.. See Harrington et al. 2013

  subroutine vaporgrow(cc, dt, ani, cni, rni, igr, nidum, temp, rimesum, presdum,&
       alphstr, sui, sup, qvs, qvi, mu, iwci, rhodum, qidum, dv, kt,     &
       nsch, npr, gammnu, i_gammnu, fourthirdspi, svpi, xxls, xxlv, xxlf,        &
       capgam, vtn, vtbarb, vtbarbm, vtbarbz, betvt, anf, cnf, rnf, iwcf, fvdum, &
       fhdum, rbdum, dsdum, rdout, dsdumout,alphm_agg, zet_ba, zet_ca)

    implicit none
    
    real, parameter :: QASMALL= 1.e-19   ! Smallest ice size for sublimation  (squared)
    integer, intent(IN) :: cc
    REAL, INTENT(IN) :: dt, ani, cni, rni, igr, nidum, temp, presdum
    REAL, INTENT(IN) :: alphstr, alphm_agg, mu, iwci, rhodum, qidum
    REAL, INTENT(IN) :: dv, kt, nsch, npr, gammnu, i_gammnu, fourthirdspi, svpi
    REAL, INTENT(IN) :: xxls, xxlv, xxlf, capgam, rimesum, sui, sup, qvs, qvi, dsdum
    REAL :: gammnubet, phii, fs, alphanr, bl, al, aa, ba
    REAL :: xn, bx, xm, am, bm, f_c1, f_c2, del1, del2, alpha, afn, del
    REAL :: xvent, ntherm, bv1, bv2, gv, bt1, bt2, gt, vtbranch, rhodep
    REAL :: videp, vmin, betavol, vi, vf, cf, maxsui, qe, Nre, phif, gamma_arg, rbdumtmp
    INTEGER :: gi, gi2
    REAL, INTENT(OUT) :: vtn, vtbarb, vtbarbm, vtbarbz
    REAL, INTENT(OUT) :: betvt, anf, cnf, rnf, iwcf, rdout, zet_ba, zet_ca
    REAL, INTENT(INOUT) :: fvdum, fhdum, dsdumout, rbdum

    real :: alphx_agg, alphvt_agg
   ! real :: C2, C1


    fvdum = 1.0
    fhdum = 1.0

    if(cc.lt.IICE3) then

    dsdumout = dsdum
    gamma_arg = NU-1.+dsdum
    gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
    phii = cni/ani*gamma_tab(gi)*i_gammnu
    fs = capgam/rni
    alphanr = ani/rni**(3./(2.+igr))

  !.. Determine  coefficients for the Best number for fall speeds
    if(phii.lt.1.0)then
       bl = 1.
       al = 2.
       aa = PI
       ba = 2.
       qe = (1.-phii)*(rbdum/RHOI) + phii
    else if(phii .gt. 1.0)then
       al = 2.0
       bl = 1.
       aa = PI*alphstr
       ba = dsdum + 1.
       qe = 1.0
    else if(phii .eq. 1.0)then
       bl = 1.
       al = 2.
       aa = pi
       ba = 2.
       qe = 1.0
    endif
    qe = min(qe, 1.0)
    
  !.. Fall speed and ventilation (Best number formulation)
    xn =  2./rbdum*(rbdum-rhodum)*G_HOME*rhodum/mu**2 * &
         (fourthirdspi*rbdum)*alphstr*al**2/(aa*qe) * qe**(3./4.)
    bx = dsdum+2.+2.*bl-ba
    
  !.. Number-average Best Number
    xm = xn*ani**bx * (gamma(nu+bx))*i_gammnu
    
  !.. The following fall speed coefficients are from 
  !.. Mitchell and Heymsfiled 2005
    f_c1 = 4.0 / (5.83 * 5.83 * SQRT(0.6))
    f_c2 = (5.83 * 5.83) / 4.0
    bm = ((f_c1 * SQRT(xm)) / &
         (2.0 * (SQRT(1.0 + f_c1 * SQRT(xm)) - 1.0) * &
         (SQRT(1.0 + f_c1 * SQRT(xm))))) - &
         (1.0e-5 * xm) / &
         (f_c2 * (SQRT(1.0 + f_c1 * SQRT(xm)) - 1.0)**2.0)
    am = ((f_c2 * (SQRT(1.0 + f_c1 * SQRT(xm)) - 1.0)**2.0) - &
         (1.0e-5 * xm)) / (xm**bm)
    
    if(xm.gt.1.e8)then
       am = 1.0865
       bm = 0.499
    endif
    
  !.. Reynolds Number
    Nre = am*xm**bm

  ! Fall speed-size exponent
    betvt =bx*bm-1.

  !.. Characteristic fall speed
    vtn = mu/rhodum*0.5 * am*(xn)**bm*ani**(betvt)
    
  !.. Number-averaged ice fall speed
   ! vtbarb = mu/rhodum*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * &
    !     (gamma(nu+bx*bm-1.))*i_gammnu

    vtbarb = vtn *(gamma(nu+betvt))*i_gammnu


  !.. Mass-averaged ice fall speed
   ! vtbarbm = mu/rhodum*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * &
   !      (gamma(nu+bx*bm-1.+2.+dsdum))/(gamma(nu+2.+dsdum))

    vtbarbm = vtn *(gamma(nu+betvt+2.+dsdum))/(gamma(nu+2.+dsdum))

    
  !.. Reflectivity-weighted fall speed (if needed)
    !vtbarbz = mu/rhodum*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * &
   !      exp(gammln(nu+bx*bm-1.+4.+2.*dsdum))/exp(gammln(nu+4.+2.*dsdum))

    vtbarbz = vtn *exp(gammln(nu+betvt+4.+2.*dsdum))/exp(gammln(nu+4.+2.*dsdum))


  !.. Calculate Ventilation
    xvent = nsch**0.333333333*Nre**0.5
    ntherm = Nre**0.5*npr**0.333333333
    if(xvent.le.1.0)then
       bv1 = 1.0
       bv2 = 0.14
       gv = 2.
    else
       bv1 = 0.86
       bv2 = 0.28
       gv = 1.
    endif
    if(ntherm.lt.1.4)then
       bt1 = 1.0
       bt2 = 0.108
       gt = 2.0
    else
       bt1 = 0.78
       bt2 = 0.308
       gt = 1.0
    endif
    fvdum = bv1 + bv2*xvent**gv
    fhdum = bt1 + bt2*ntherm**gt
    
    if(temp.le.T0) then !.. If T < T0 do vapor growth/sublimation (otherwise assume water on surface)
       
  !.. Fall speed needed to determine when branching occurs
       vtbranch = vtbarbm
       
       if(sup.ge.0.) then
          maxsui = 1.
       elseif(sui.ge.0..and.qvi.lt.qvs) then
          maxsui = ((sui+1.)*qvi)/(qvs-qvi) - qvi/(qvs-qvi)
          maxsui = min(maxsui,1.)
          maxsui = max(maxsui,0.)
       else
          maxsui = 0.
       endif

  !.. Vapor growth density 
       if(igr.le.1.) then !.. Planar
          if(vtbranch.gt.0.) then
             if(ani.gt.SQRT((dv*PI*2.*cni)/(vtbranch*nu))) then
                rhodep = (RHOI*igr)*maxsui + RHOI*(1.-maxsui)
             else
                rhodep = RHOI
             endif
          else
             rhodep = RHOI
          endif
       else !.. Columnar
          rhodep = (RHOI/igr)*maxsui + RHOI*(1.-maxsui)
       endif
  !.. high limit on rhodep for vapor growth = 700 kg m^-3 
  !.. Cotton et al. 2012 (QJRMS)
       rhodep = min(rhodep,700.)
       
  !.. Vapor growth solution which includes heating from rime mass
  !.. See J.P. Chen's thesis
       alpha = (dv*fvdum*svpi*xxls)/(RV*kt*fhdum*temp)
       if(nidum.gt.0.0.and.rimesum.gt.0.0) then
          del1 = (xxlf*(rimesum/nidum)/(4.*PI*kt*fhdum*capgam)) * &
               ((temp + alpha*((xxls/(RV*temp))-1.))**(-1.0))
       else
          del1 = 0.0
       endif
       del2 = sui * &
            (((temp/alpha) + ((xxls/(RV*temp))-1.))**(-1.0))
       del = del1 + del2
       afn = ((dv*fvdum*polysvp(temp,1))/(RV*temp)) * &
            (sui - del*((xxls/(RV*temp))-1.))
 
  !.. During sublimation using polynomial removal of density
       if(afn.lt.0.0) then   
          rhodep = rbdum     
          videp = rni**3
          vmin  = (10.e-6)**3
          if(Vmin.lt.videp)then
             betavol = log(RHOI/rbdum)*1./(log(Vmin/videp))
             rhodep = rbdum*(1.+betavol)
          else
             rhodep = rbdum
          endif
       endif
       rhodep=max(rhodep,50.)
       rhodep=min(rhodep,RHOI)
       
       gamma_arg = NU+2.+dsdum
       gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
       gammnubet = gamma_tab(gi)

  !..  characteristic r-axis and a-axis after growth timestep
       rnf = (max((rni**2 + 2.*afn*fs/rhodep*gammnu/gammnubet*dt),QASMALL))**(0.5) 
       anf = alphanr*rnf**(3./(2.+igr)) 
       
  !.. Do not sublimation change the shape of ice from 
  !.. prolate to oblate or vice versa
       phif = phii*(rnf**3/rni**3)**((igr-1.)/(igr+2.))
       if(sui.lt.0.0.or.afn.lt.0.0) then
          if(phii.gt.1.0.and.phif.lt.1.0) then
             phif = phii
             alphanr = ani/rni
             anf = alphanr*rnf
          endif
          if(phii.lt.1.0.and.phif.gt.1.0) then
             phif = phii
             alphanr = ani/rni
             anf = alphanr*rnf
          endif

  !.. Do not let sublimation create extreme shapes
          if(phii.gt.1.) then
             if(phif.gt.phii) then
                phif = phii
                alphanr = ani/rni
                anf = alphanr*rnf
             endif
          else
             if(phif.lt.phii) then
                phif = phii
                alphanr = ani/rni
                anf = alphanr*rnf
             endif
          endif
       endif
       
       vi = fourthirdspi*rni**3*gamma_tab(gi)*i_gammnu 
       vf = fourthirdspi*rnf**3*gamma_tab(gi)*i_gammnu
       rdout = rhodep
       rbdumtmp = rbdum*(vi/vf) + rhodep*(1.-vi/vf)
       rbdumtmp = min(rbdumtmp,RHOI)
       iwcf = nidum*rbdumtmp*vf
       
  !.. Update delta* from vapor growth/sublimation
       if(igr.ne.1.0)then

          if(anf.gt.(1.1*ao)) then
             dsdumout = (3.*log(rnf)-2.*log(anf)-log(ao))/ &
                  (log(anf)-log(ao))
          else
             dsdumout=1.
          endif
       endif

  !.. Do not let particles sublimate to sizes that are too small
       if(afn.lt.0..and.rnf.lt.1.e-6) then
          rbdum  =RHOI
          rdout  =RHOI
          dsdumout   =1.
          phif    =1.
          alphanr = ani/rni
          anf = alphanr*rnf
       endif
       
       if(afn.lt.0..and.anf.lt.1.e-6) then
          rbdum  =RHOI
          rdout  =RHOI
          dsdumout   =1.
          phif    =1.
          alphanr = ani/rni
          anf = alphanr*rnf
       endif
      
       
  !.. Sublimation check
       if(afn.lt.0..and.dsdumout.le.0.) then
          dsdumout=1.
          anf  =rnf
       endif
       
  !.. C-axis after vapor growth
       gamma_arg = NU-1.+dsdumout
       gi=MIN((MAX((NINT((gamma_arg*100000.)-355000+1)),1)),505001)
       cnf = phif*anf*gammnu/gamma_tab(gi)
       
    else  !.. T > T0
       
       anf =ani
       cnf =cni
       rnf =rni
       iwcf=iwci
    endif


    else

      if (qidum.ge.1.e-12) then


        ! Coefficients for projected area from Mitchell 1996
        alphx_agg = 8.*alphm_agg*alph_ar_agg**(-kc)*G_HOME*rhodum/(PI*mu**2)

        !Average Best Number
        xm = alphx_agg*ani**betx_agg * gamma(NU_AGG+betx_agg)/gamma(NU_AGG) * &
        pbeta_chi

        bm = (C1*sqrt(xm)/&
        (2.*(sqrt(1+C1*sqrt(xm)) - 1.)*sqrt(1.+C1*sqrt(xm)))) - &
        ((a0*b0*xm**b0)/(C2*(sqrt(1.+C1*sqrt(xm))-1.)**2))

        am = &
        (C2*((sqrt(1.+C1*sqrt(xm))-1.)**2) - a0*xm**b0)/(xm**bm)

        zet_ba = chi_ba*bm-0.5*bet_ba
        zet_ca = chi_ca*bm-0.5*bet_ca

        !Average Reynolds Number
        Nre = am*xm**bm

        alphvt_agg = mu/rhodum*0.5*am*(alphx_agg)**bm
        betvt = betx_agg*bm-1.

        vtn = alphvt_agg*ani**betvt

        vtbarb = vtn*gamma(NU_AGG+betvt)/gamma(NU_AGG) * &
        beta_moms(zet_ba,zet_ca,aba,bba,bcb)

        vtbarbm = vtn*gamma(NU_AGG+betvt+betm_agg)/gamma(NU_AGG+betm_agg) * &
        beta_moms(zet_ba+1.,zet_ca+1.,aba,bba,bcb)/phiba_phica_bar

    ! Relative fall speed dispersion (Coefficient of Variation; function of size)
   ! sig_vt = sqrt((gamma(NU_AGG+2.*betvt_agg)/gamma(NU_AGG+betvt_agg)**2)*&
    !               ((beta_moms(2.*zet_ba,2.*zet_ca,aba,bba,bcb)/&
    !                beta_moms(zet_ba,zet_ca,aba,bba,bcb)**2))-1)

        !sig_vt = sqrt((beta_moms(2.*zet_ba,2.*zet_ca,aba,bba,bcb)/&
         !           (beta_moms(zet_ba,zet_ca,aba,bba,bcb)**2))-1.)


  !  if (bm_agg.lt.0.5) then

   ! print*, 'bm_agg=',bm_agg
   ! print*, 'zetba=',zet_ba,'zetca=',zet_ca
   ! print*, 'numterm=',beta_moms(2.*zet_ba,2.*zet_ca,aba,bba,bcb)
   ! print*, 'demterm=',beta_moms(zet_ba,zet_ca,aba,bba,bcb)**2
   ! print*, 'sig_vt=',sig_vt
   ! stop
   ! endif


 if (vtbarb/=vtbarb.or.vtbarbm/=vtbarbm) then
  print*, 'vt nan'
  print*, 'cc=',cc
  print*, 'qidum=',qidum
  print*, 'nidum=',nidum
  print*, 'alphm=',alphm_agg
  print*, 'betx=',betx_agg
  print*, 'alphx=',alphx_agg
  print*, 'bm=',bm
  print*, 'am=',am
  print*,'rhoa=',rhodum
  print*, 'rhoagg=',rbdum
  print*,'mu=',mu
  print*,'an=',ani
  print*, 'xm=',xm
  print*, 'NRe=',Nre
  print*, 'alphvt=',alphvt_agg
  print*, 'betvt=',betvt
  print*, 'vtbarb=',vtbarb
  print*, 'vtbarbm=',vtbarbm
  print*, 'AR=',alph_ar_agg*ani**bet_ar_agg
  print*, 'qagg_mom=',qidum/nidum
  print*, 'qagg_int=',(1/rhodum)*alphm_agg*ani**betm_agg * gamma(NU_AGG+betm_agg)/gamma(NU_AGG) * phiba_phica_bar
  print*, 'phi_const=', beta_moms(zet_ba,zet_ca,aba,bba,bcb)
  stop
 endif



     ! Ventilation
     xvent = nsch**0.333333333*Nre**0.5
     ntherm = Nre**0.5*npr**0.333333333
     if(xvent.le.1.0)then
     bv1 = 1.0
     bv2 = 0.14
     gv = 2.
     else
     bv1 = 0.86
     bv2 = 0.28
     gv = 1.
     endif
     if(ntherm.lt.1.4)then
     bt1 = 1.0
     bt2 = 0.108
     gt = 2.0
     else
     bt1 = 0.78
     bt2 = 0.308
     gt = 1.0
    endif
    fvdum = bv1 + bv2*xvent**gv
    fhdum = bt1 + bt2*ntherm**gt


   if(temp.le.T0) then

!.. Vapor growth solution which includes heating from rime mass
!.. See J.P. Chen's thesis
    alpha = (dv*fvdum*polysvp(temp,1)*xxls)/(RV*kt*fhdum*temp)
    if(nidum.gt.0.0.and.rimesum.gt.0.0) then

    del1 = (xxlf*(rimesum/nidum)/(4.*PI*kt*fhdum*capgam)) * &
    ((temp + alpha*((xxls/(RV*temp))-1.))**(-1.0))
    else
    del1 = 0.0
    endif
    del1 =0.0
    del2 = sui * &
    (((temp/alpha) + ((xxls/(RV*temp))-1.))**(-1.0))
    del = del1 + del2
    afn = ((dv*fvdum*polysvp(temp,1))/(RV*temp)) * &
    (sui - del*((xxls/(RV*temp))-1.))

    ! NEW AGGREGATE MASS FROM VAPOR
    iwcf = iwci+4.*PI*nidum*afn*capgam*dt

    else

     iwcf=iwci
    !vtbarb=vtbarbm !TEST TO REMOVE SIZE SORTING FOR MELTING ICE

    endif




    endif



    endif

  end subroutine vaporgrow

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. This routine is used to calculate the distribution-averaged 
  !.. capacitance (see Harrington et al. 2013)

  real function capacitance_gamma(ani,dsdum,alphstr,i_gammnu)
    
    implicit none
    
    REAL, INTENT(IN) :: ani, dsdum, alphstr, i_gammnu
    REAL :: a1, a2, b1, b2, c1, c2, d1, d2, gammad1, gammad2
    
  !.. Oblate Spheroid
    if (dsdum.le.1.0) then
       a1 = 0.6369427      
       a2 = 0.57*a1
       b1 = 0.0
       b2 = 0.95   
       c1 = a1*alphstr**b1
       c2 = a2*alphstr**b2
       d1 = b1*(dsdum - 1.0) + 1.0
       d2 = b2*(dsdum - 1.0) + 1.0
  !.. Prolate Spheroid
    elseif (dsdum.gt.1.) then
       a1 = 0.5714285
       a2 = 0.75*a1
       b1 = -1.0
       b2 = -0.18 
       c1 = a1*alphstr**(b1+1.0)
       c2 = a2*alphstr**(b2+1.0)
       d1 = b1*(dsdum - 1.0) + dsdum
       d2 = b2*(dsdum - 1.0) + dsdum
    endif
      
    if(dsdum.eq.1..and.nu.eq.1.) then
       gammad1 = (gamma(nu+1.0))
       capacitance_gamma = ani*gammad1*i_gammnu
    elseif(dsdum.le.1.)then
       gammad1 = (gamma(nu+d1))
       gammad2 = (gamma(nu+d2))
       capacitance_gamma = c1*ani**d1 * gammad1*i_gammnu + & 
            c2*ani**d2 * gammad2*i_gammnu
    elseif(dsdum.gt.1.) then
       gammad1 = (gamma(nu+d1))
       gammad2 = (gamma(nu+d2))
       capacitance_gamma = c1*ani**d1 * gammad1*i_gammnu  + &
            c2*ani**d2 * gammad2*i_gammnu
    endif
    
    return
  end function capacitance_gamma

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  SUBROUTINE wet_growth_check(temp, rhodum, xxlv, xxlf, qvdum, dv, kt, qs0, &
       fvdum, fhdum, rimedum, rni, nidum, dgflag)
    
    IMPLICIT NONE

  !.. Check for wet growth conditions
  !.. See Lamb and Verlinde (2011)
    REAL, INTENT(IN) :: temp, rhodum, xxlv, xxlf, qvdum, dv, kt
    REAL, INTENT(IN) :: fvdum, fhdum, rni, nidum
    REAL, INTENT(IN) :: qs0, rimedum
    REAL :: wetg, dum
    LOGICAL, INTENT(INOUT) :: dgflag
    
    dum=nidum*NU*2.*rni
    wetg=2.*PI*dum*(kt*fhdum*(T0-temp)+rhodum*xxlv*dv*fvdum*(qs0-qvdum))/(xxlf + (CPW*(temp-T0)))
    if(rimedum/rhodum.gt.wetg) then
       dgflag = .false.
    endif
    
  end subroutine wet_growth_check

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  real function polysvp(T,TYPE)
    
  !--------------------------------------------------------------
  !.. Taken from 'module_mp_morr_two_moment.F' (WRFV3.4)

  !..  COMPUTE SATURATION VAPOR PRESSURE

  !..  POLYSVP RETURNED IN UNITS OF PA.
  !..  T IS INPUT IN UNITS OF K.
  !..  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

  !.. REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992,
  !.. TABLE 4 (RIGHT-HAND COLUMN)
  !--------------------------------------------------------------

    IMPLICIT NONE
    
    REAL DUM
    REAL T
    INTEGER TYPE
  !.. ice
    real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
    data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
         6.11147274, 0.503160820, 0.188439774e-1, &
         0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
         0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
    
  !.. liquid
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    
  !.. V1.7
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
         6.11239921, 0.443987641, 0.142986287e-1, &
         0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
         0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
    real dt
    
  !.. ICE
    IF (TYPE.EQ.1) THEN

       dt = max(-80.,t-273.16)
       polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
       polysvp = polysvp*100.
       
    END IF

  !.. LIQUID
    IF (TYPE.EQ.0) THEN
       
       dt = max(-80.,t-273.16)
       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp = polysvp*100.
         
    END IF
      
  end function polysvp

 !-------------------------------------------------------------------------
 real function beta(x,y)
! Beta function
 implicit none

 real :: x, y,xmax
 !double precision :: xd, yd

 xmax = 50.

 !beta = gamma_tab(MIN((MAX((NINT((x*100000.)-355000+1)),1)),505001))*&
 !       gamma_tab(MIN((MAX((NINT((y*100000.)-355000+1)),1)),505001))/&
 !       gamma_tab(MIN((MAX((NINT(((x+y)*100000.)-355000+1)),1)),505001))

 beta = exp(gammln(x)+gammln(y)-gammln(x+y))

 !beta = gamma(x)*gamma(y)/gamma(x+y)

 if(beta/=beta.or.beta.gt.1.e20) then
 print*, 'NaN or Inf expression in beta'
 print*,'x=',x,'y=',y
 stop
 endif

 return
 end function beta



!-------------------------------------------------------------------------
 real function beta_moms(m,n,aba,bba,bcb)
 ! Equation 8 from Dunnavan et al. 2019 JAS
 implicit none

 !real :: aba, bba, bcb
 real :: m, n, aba, bba, bcb

 !m = dble(m)
 !n = dble(n)
 !aba = dble(aba)
 !bba = dble(bba)
 !bcb = dble(bcb)

 ! agg_switch toggles between spheres/fractals (0) and ellipsoids (1)

 if (m.eq.0..and.n.eq.0.) then

 beta_moms = 1.0

 else

 if (agg_switch.eq.1) then

 beta_moms = beta(m+n+aba,bba)*beta(n+aba+bba,bcb)/(beta(aba,bba)*beta(aba+bba,bcb))

 elseif (agg_switch.eq.2) then

 beta_moms = phiagg**n

 else

 beta_moms = 1.0

 endif
endif

 return
 end function beta_moms

!------------------------------------------------------------------------

!------------------------------------------------------------------------
 real function delta_vtm(vtmx,vtmy)
 ! Calculates fallspeed difference approximation following Milbrandt and Yau (2005b)
 ! (which follows Wisner et al. 1972 and Murakami 1990)
 implicit none

 real :: r, vtmx, vtmy

 r = 0.04
 delta_vtm = sqrt((vtmx-vtmy)**2. + r*vtmx*vtmy )

 return
 end function delta_vtm


!------------------------------------------------------------------------
 real function delta_vtm_seifert(vtnx,betax,betmx,betvx,vtny,betay,betvy,vba,vca,&
                                S_switch)
 ! Calculates fallspeed difference approximation (which follows Seifert et al. 2006/2014)
 implicit none

 real :: vtnx, vtny, dS, betvx, betmx, betvy, betax, betay, vba, vca

 integer :: S_switch
 integer :: kcol

 real, dimension (1:3) :: Cnk, ki


 ! Binomial coefficients (2 choose k)
 Cnk = (/1., 2., 1./)
 ! k index
 ki = (/0., 1., 2./)

 dS = 0.

 if(S_switch.eq.0) then
    ! P-P Aggregation

    do kcol=1,3
        dS = dS + Cnk(kcol)*(-1.)**ki(kcol)*gamma(NU+(2.-ki(kcol))*betvx+betax+betmx)*gamma(NU+ki(kcol)*betvx+betax+betmx)
    enddo

    dS = vtnx/gamma(NU+betax+betmx)*sqrt(dS)


 elseif(S_switch.eq.1) then
    ! P-A Aggregation

    do kcol=1,3
        dS = dS + Cnk(kcol)*(-1.)**(ki(kcol))*vtny**(2.-ki(kcol))*vtnx**(ki(kcol))*gamma(NU+ki(kcol)*betvx+betax+betmx)*&
            gamma(NU_AGG+(2.-ki(kcol))*betvy+betay)*beta_moms((2.-ki(kcol))*vba+bet_ba,(2.-ki(kcol))*&
                  vca+bet_ca,aba,bba,bcb)
    enddo

    dS =sqrt(dS/(gamma(NU+betax+betmx)*pbeta_area*gamma(NU_AGG+betay)))


 elseif(S_switch.eq.2) then


    do kcol=1,3

    dS = dS + Cnk(kcol)*(-1.)**ki(kcol)*gamma(NU_AGG+(2.-ki(kcol))*betvy+betay)*gamma(NU_AGG+ki(kcol)*betvy+betay)*&
              beta_moms(ki(kcol)*vba+bet_ba,ki(kcol)*vca+bet_ca,aba,bba,bcb)*&
              beta_moms((2.-ki(kcol))*vba+bet_ba,(2.-ki(kcol))*vca+bet_ca,aba,bba,bcb)

    enddo

    dS = vtny/(gamma(NU_AGG+betay)*pbeta_area)*sqrt(dS)


 endif

! A-A Aggregation

 delta_vtm_seifert = dS


 return
 end function delta_vtm_seifert

!--------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------!

  subroutine starbuck(rhoairdum,tempdum,phidum1,phidum2,rhodum1,rhodum2,&
                      qdum1,vadum1,vcdum1,ndum1,andum1,alphadum1,betadum1,alphvdum1,alphmdum1,betmdum1,alphcdum1,&
                      betcdum1,qdum2,vadum2,vcdum2,ndum2,andum2,alphadum2,betadum2,alphvdum2,alphmdum2,betmdum2, &
                      alphcdum2,betcdum2,qdum3,vadum3,vcdum3,ndum3,andum3,alphadum3,betadum3,dvtn_dum,&
                      dvtm_dum,dvtc_dum,qagg1,qagg2,qagg3,nagg1,nagg2,nagg3,nmon3,Vagg1,Vagg2,&
                      Vagg3,Vcagg1,Vcagg2,Vcagg3)

       ! This is the new aggregation method. First, I'm going to follow the
       ! approach of Milbrandt and Yau (2005) where I'll use the Wisner approx.
       ! for that nasty fall speed term. This will let me analytically calculate
       ! the transfer rates, thereby bypassing the aggregation LU. Future updates
       ! will incorporate the full hydrodynamic kernel.
       ! ELD

      implicit none

      real, intent(IN) :: qdum1,ndum1,qdum2,ndum2,qdum3,ndum3,andum1,andum2,andum3
      real, intent(IN) :: vadum1, vcdum1, vadum2, vcdum2, vadum3, vcdum3
      real, intent(IN) :: alphadum1, betadum1, alphadum2, betadum2
      real, intent(IN) :: alphmdum1, betmdum1,alphmdum2, betmdum2
      real, intent(IN) :: alphvdum1, alphvdum2
      real, intent(IN) :: alphcdum1, betcdum1,alphcdum2, betcdum2
      real, intent(IN) :: alphadum3, betadum3
      real, intent(INOUT) :: qagg1,qagg2,qagg3,nagg1,nagg2,nagg3, nmon3
      real, intent(INOUT) :: Vagg1,Vagg2,Vagg3,Vcagg1,Vcagg2,Vcagg3
      real, intent(IN) :: rhoairdum, tempdum
      real, intent(IN) :: rhodum1, rhodum2, phidum1, phidum2

      real, intent(IN), dimension(1:6) :: dvtn_dum,dvtm_dum,dvtc_dum


      ! New mass transfer rates
      real :: pprxfer, pcrxfer, parxfer
      real :: cprxfer, ccrxfer, carxfer
      real :: aprxfer, acrxfer

      ! New Volume*aspect ratio transfer rates
      real :: ppcxfer, pccxfer, pacxfer
      real :: cpcxfer, cccxfer, cacxfer
      real :: apcxfer, accxfer

      ! New number transfer rates
      real :: ppenxfer, pcenxfer, paenxfer
      real :: ccenxfer, caenxfer
      real :: apenxfer, acenxfer, aaenxfer

      ! New rho/phi dependent efficiencies based on Jensen et al. 2018
      real :: effnew_pp, effnew_cc, effnew_pc
      real :: rhoeffmax, phieffmax, rhoeff, phieff

      real :: sourcep, sinkp, sourcec, sinkc, ratioagg

       real, parameter :: QNSMALL= 1.25e-7   ! Smallest ice number
       real :: eff,ncrossloss,rcrossloss,nselfloss, &
            rselfloss,ncrossgain,rcrossgain,nselfgain,rselfgain

       integer :: kcol ! index for collection kernel expansion
       integer, parameter :: eff_switch = 0 ! Switch for Eagg
       real, parameter :: effc = 0.05


      ! Calculate efficiencies

      if (eff_switch.eq.0) then
       ! Set agg efficiency to constant
            eff = effc
      else
    ! Connolly 2012 aggregation parameterization
            eff = min(0.2,10.**(0.035*tempdum-0.7))

      endif

 !---------------------------------------------------------!
 !---------------------------------------------------------!

      ! New efficiencies
      !.. Plates + Plates = Aggregates
       rhoeffmax = rhodum1
       phieffmax = min(phidum1,1./phidum1)

       if(rhoeffmax.le.400.) then
          rhoeff = 1.
       else
          rhoeff = -0.001923*rhoeffmax + 1.76916
       endif
       rhoeff = max(rhoeff,0.)
       rhoeff = min(rhoeff,1.)

  !.. Based on Connolly et al 2012
  !.. See Jensen et al. 2018b
       if(phieffmax.le.0.03) then
          phieff = 1.
       elseif(phieffmax.gt.0.03.and.phieffmax.lt.0.85) then
          phieff = 0.0001*((phieffmax+0.07)**(-4.))
       else
          phieff = 0.
       endif

       if(phieff.le.0.001) then
          phieff = 0.
       endif

       phieff = max(phieff,0.0)
       phieff = min(phieff,1.)
       effnew_pp = phieff * rhoeff

!---------------------------------------------------------!

      ! New efficiencies
      !.. Columns + Columns = Aggregates
       rhoeffmax = rhodum2
       phieffmax = min(phidum2,1./phidum2)

       if(rhoeffmax.le.400.) then
          rhoeff = 1.
       else
          rhoeff = -0.001923*rhoeffmax + 1.76916
       endif
       rhoeff = max(rhoeff,0.)
       rhoeff = min(rhoeff,1.)

  !.. Based on Connolly et al 2012
  !.. See Jensen et al. 2018b
       if(phieffmax.le.0.03) then
          phieff = 1.
       elseif(phieffmax.gt.0.03.and.phieffmax.lt.0.85) then
          phieff = 0.0001*((phieffmax+0.07)**(-4.))
       else
          phieff = 0.
       endif

       if(phieff.le.0.001) then
          phieff = 0.
       endif

       phieff = max(phieff,0.0)
       phieff = min(phieff,1.)
       effnew_cc = phieff * rhoeff

 !--------------------------------------------------------!

      ! New efficiencies
      !.. Plates + Columns = Aggregates
       rhoeffmax = max(rhodum1,rhodum2)
       phieffmax = max(min(phidum1,1./phidum1),min(phidum2,1./phidum2))
       if(rhoeffmax.le.400.) then
          rhoeff = 1.
       else
          rhoeff = -0.001923*rhoeffmax + 1.76916
       endif
       rhoeff = max(rhoeff,0.)
       rhoeff = min(rhoeff,1.)

  !.. Based on Connolly et al 2012
  !.. See Jensen et al. 2018b
       if(phieffmax.le.0.03) then
          phieff = 1.
       elseif(phieffmax.gt.0.03.and.phieffmax.lt.0.85) then
          phieff = 0.0001*((phieffmax+0.07)**(-4.))
       else
          phieff = 0.
       endif

       if(phieff.le.0.001) then
          phieff = 0.
       endif

       phieff = max(phieff,0.00)
       phieff = min(phieff,1.)
       effnew_pc = phieff * rhoeff


 !--------------------------------------------------------!
 !--------------------------------------------------------!


   ! ELD: Apparently these multipliers are correct (according to JYH)
   !      We're not entirely sure though but we can check this out numerically.
   !      RAMS apparently has different coefficients. ELD.

        nselfloss = 1.0
        rselfloss = 1.0
        nselfgain = 1./2.
        rselfgain = 1.0

        ncrossloss = 1.0
        rcrossloss = 1.0
        ncrossgain = 1./2.
        rcrossgain = 1.0

! Fill up transfer variables with zeros
        ! Mass
        pprxfer  = 0.0
        ccrxfer  = 0.0
        pcrxfer  = 0.0
        cprxfer  = 0.0
        parxfer  = 0.0
        carxfer  = 0.0

        ! Number
        ppenxfer = 0.0
        ccenxfer = 0.0
        pcenxfer = 0.0
        paenxfer = 0.0
        caenxfer = 0.0
        aaenxfer = 0.0

        ! Volume*aspect ratio
        ppcxfer = 0.0
        cccxfer = 0.0
        pccxfer = 0.0
        cpcxfer = 0.0
        pacxfer = 0.0
        cacxfer = 0.0

!----------------------------------------------!
! So total combinations for aggregation are:
!----------------------------------------------!
!    Collected + Collectee = New Category
!----------------------------------------------!
!    Plates + Plates   = Agg               (1)
!    Columns + Columns = Agg               (2)
!    Plates + Columns  = Agg               (3)
!    Columns + Plates  = Agg               (4)
!    Plates + Agg      = Agg               (5)
!    Columns + Agg     = Agg               (6)
!    Agg + Agg         = Agg (number only) (7)
!----------------------------------------------!


! LEE NOTE: Make sure to square out the dvtn/dvtm/dvtc etc. stuff

        do kcol=1,3 ! collection kernel expansion

            ! Plates + Plates = Agg
            ! Mass, Volume*Aspect Ratio, Number
            pprxfer = pprxfer + &
            Cnk(kcol)*gamma(nu+betmdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum1*ki(kcol))

            ppcxfer = ppcxfer + &
            Cnk(kcol)*gamma(nu+betcdum1+2.-ki(kcol))*gamma(nu+ki(kcol))

            ppenxfer = ppenxfer + &
            Cnk(kcol)*gamma(nu+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum1*ki(kcol))

            ! Columns + Columns = Agg
            ! Mass, Volume*Aspect Ratio, Number

            ccrxfer = ccrxfer + &
            Cnk(kcol)*gamma(nu+betmdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum2*ki(kcol))

            cccxfer = cccxfer + &
            Cnk(kcol)*gamma(nu+betcdum2+2.-ki(kcol))*gamma(nu+ki(kcol))

            ccenxfer = ccenxfer + &
            Cnk(kcol)*gamma(nu+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum2*ki(kcol))

            ! Plates + Columns = Agg
            ! Mass, Volume*Aspect Ratio, Number

            pcrxfer = pcrxfer + &
            Cnk(kcol)*alphadum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      alphadum2**(0.5*betadum2*ki(kcol))*&
                      andum1**(betmdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      andum2**(0.5*betadum2*ki(kcol))*&
                      gamma(nu+betmdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum2*ki(kcol))

            pccxfer = pccxfer + &
            Cnk(kcol)*alphadum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      alphadum2**(0.5*betadum2*ki(kcol))*&
                      andum1**(betcdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      andum2**(0.5*betadum2*ki(kcol))*&
                      gamma(nu+betcdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum2*ki(kcol))

            pcenxfer = pcenxfer + &
            Cnk(kcol)*alphadum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      alphadum2**(0.5*betadum2*ki(kcol))*&
                      andum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      andum2**(0.5*betadum2*ki(kcol))*&
                      gamma(nu+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum2*ki(kcol))


            ! Columns + Plates = Agg
            ! Mass, Volume*Aspect Ratio, Number

           cprxfer = cprxfer + &
            Cnk(kcol)*alphadum2**(0.5*betadum2*(2.-ki(kcol)))*&
                      alphadum1**(0.5*betadum1*ki(kcol))*&
                      andum2**(betmdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      andum1**(0.5*betadum1*ki(kcol))*&
                      gamma(nu+betmdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum1*ki(kcol))


            cpcxfer = cpcxfer + &
            Cnk(kcol)*alphadum2**(0.5*betadum2*(2.-ki(kcol)))*&
                      alphadum1**(0.5*betadum1*ki(kcol))*&
                      andum2**(betcdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      andum1**(0.5*betadum1*ki(kcol))*&
                      gamma(nu+betcdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu+0.5*betadum1*ki(kcol))

            ! Plates + Agg = Agg
            ! Mass, Volume*Aspect Ratio, Number

            parxfer = parxfer +   &
            Cnk(kcol)*alphadum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      alphadum3**(0.5*betadum3*ki(kcol))*&
                      andum1**(betmdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      andum3**(0.5*betadum3*ki(kcol))*&
                      gamma(nu+betmdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)

            pacxfer = pacxfer +   &
            Cnk(kcol)*alphadum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      alphadum3**(0.5*betadum3*ki(kcol))*&
                      andum1**(betcdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      andum3**(0.5*betadum3*ki(kcol))*&
                      gamma(nu+betcdum1+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)


            paenxfer = paenxfer + &
            Cnk(kcol)*alphadum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      alphadum3**(0.5*betadum3*ki(kcol))*&
                      andum1**(0.5*betadum1*(2.-ki(kcol)))*&
                      andum3**(0.5*betadum3*ki(kcol))*&
                      gamma(nu+0.5*betadum1*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)


            ! Columns + Agg = Agg
            ! Mass, Volume*Aspect Ratio, Number


            carxfer = carxfer + &
            Cnk(kcol)*alphadum2**(0.5*betadum2*(2.-ki(kcol)))*&
                      alphadum3**(0.5*betadum3*ki(kcol))*&
                      andum2**(betmdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      andum3**(0.5*betadum3*ki(kcol))*&
                      gamma(nu+betmdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)


            cacxfer = cacxfer + &
            Cnk(kcol)*alphadum2**(0.5*betadum2*(2.-ki(kcol)))*&
                      alphadum3**(0.5*betadum3*ki(kcol))*&
                      andum2**(betcdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      andum3**(0.5*betadum3*ki(kcol))*&
                      gamma(nu+betcdum2+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)


            caenxfer = caenxfer + &
            Cnk(kcol)*alphadum2**(0.5*betadum2*(2.-ki(kcol)))*&
                      alphadum3**(0.5*betadum3*ki(kcol))*&
                      andum2**(0.5*betadum2*(2.-ki(kcol)))*&
                      andum3**(0.5*betadum3*ki(kcol))*&
                      gamma(nu+0.5*betadum2*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)

            ! Agg + Agg = Agg
            ! Number

            aaenxfer = aaenxfer + &
            Cnk(kcol)*gamma(nu_agg+0.5*betadum3*(2.-ki(kcol)))*&
                      gamma(nu_agg+0.5*betadum3*ki(kcol))*&
                      pbeta_colx(kcol)*pbeta_coly(kcol)


        enddo


!---------------------------------------------------------------------------!
        ! Self Collection
        ! Mass
        pprxfer = ndum1*ndum1*(Eff*effnew_pp*dvtm_dum(1)/(gamma(nu)*gamma(nu)))*&
                   alphmdum1*alphadum1*andum1**(betmdum1+betadum1)*pprxfer

        ccrxfer = ndum2*ndum2*(Eff*effnew_cc*dvtm_dum(2)/(gamma(nu)*gamma(nu)))*&
                   alphmdum2*alphadum2*andum2**(betmdum2+betadum2)*ccrxfer

        ! Volume*Aspect Ratio
        ppcxfer = ndum1*ndum1*(Eff*effnew_pp*dvtc_dum(1)/(gamma(nu)*gamma(nu)))*&
                  alphcdum1*alphadum1*andum1**(betcdum1+betadum1)*ppcxfer

        cccxfer = ndum2*ndum2*(Eff*effnew_cc*dvtc_dum(2)/(gamma(nu)*gamma(nu)))*&
                  alphcdum2*alphadum2*andum2**(betcdum2+betadum2)*cccxfer


        ! Number
        ppenxfer = ndum1*ndum1*(Eff*effnew_pp*dvtn_dum(1)/(gamma(nu)*gamma(nu)))*&
                   alphadum1*andum1**(betadum1)*ppenxfer

        ccenxfer = ndum2*ndum2*(Eff*effnew_cc*dvtn_dum(2)/(gamma(nu)*gamma(nu)))*&
                   alphadum2*andum2**(betadum2)*ccenxfer


        aaenxfer = ndum3*ndum3*(Eff*dvtn_dum(6)/(gamma(nu_agg)*gamma(nu_agg)))*&
                   alphadum3*andum3**(betadum3)*aaenxfer





!---------------------------------------------------------------------------!

!---------------------------------------------------------------------------!
       ! Cross Collection
       ! Mass
       pcrxfer =  ndum1*ndum2*(Eff*effnew_pc*dvtm_dum(3)/(gamma(nu)*gamma(nu)))*alphmdum1*pcrxfer

       cprxfer =  ndum2*ndum1*(Eff*effnew_pc*dvtm_dum(4)/(gamma(nu)*gamma(nu)))*alphmdum2*cprxfer

       parxfer =  ndum1*ndum3*(Eff*effnew_pp*dvtm_dum(5)/(gamma(nu)*gamma(nu_agg)))*alphmdum1*parxfer

       carxfer =  ndum2*ndum3*(Eff*effnew_cc*dvtm_dum(6)/(gamma(nu)*gamma(nu_agg)))*alphmdum2*carxfer

       ! Volume*Aspect Ratio

       pccxfer =  ndum1*ndum2*(Eff*effnew_pc*dvtc_dum(3)/(gamma(nu)*gamma(nu)))*alphcdum1*pccxfer

       cpcxfer =  ndum2*ndum1*(Eff*effnew_pc*dvtc_dum(4)/(gamma(nu)*gamma(nu)))*alphcdum2*cpcxfer

       pacxfer =  ndum1*ndum3*(Eff*effnew_pp*dvtc_dum(5)/(gamma(nu)*gamma(nu_agg)))*alphcdum1*pacxfer

       cacxfer =  ndum2*ndum3*(Eff*effnew_cc*dvtc_dum(6)/(gamma(nu)*gamma(nu_agg)))*alphcdum2*cacxfer

       ! Number

       pcenxfer = ndum1*ndum2*(Eff*effnew_pc*dvtn_dum(3)/(gamma(nu)*gamma(nu)))*pcenxfer

       paenxfer = ndum1*ndum3*(Eff*effnew_pp*dvtn_dum(4)/(gamma(nu)*gamma(nu_agg)))*paenxfer

       caenxfer = ndum2*ndum3*(Eff*effnew_cc*dvtn_dum(5))/(gamma(nu)*gamma(nu_agg))*caenxfer

!---------------------------------------------------------------------------!


  !.. Do not over-deplete mass from aggregation
          sinkp = rcrossloss*pcrxfer + &
                  rcrossloss*parxfer + &
                  rselfloss*pprxfer
          sourcep = qdum1

          sinkc = rcrossloss*cprxfer + &
                  rcrossloss*carxfer + &
                  rselfloss*ccrxfer
          sourcec = qdum2

          if(sinkp.gt.sourcep.and.sourcep.gt.1.e-8) then
             ratioagg= sourcep/sinkp
             pcrxfer = pcrxfer*ratioagg
             parxfer = parxfer*ratioagg
             pprxfer = pprxfer*ratioagg
          endif

          if(sinkc.gt.sourcec.and.sourcec.gt.1.e-8) then
             ratioagg= sourcec/sinkc
             cprxfer = cprxfer*ratioagg
             carxfer = carxfer*ratioagg
             ccrxfer = ccrxfer*ratioagg
          endif



  !.. Do not over-deplete Volume*Aspect Ratio from aggregation
        sinkp = rcrossloss*pccxfer + &
        rcrossloss*pacxfer + &
        rselfloss*ppcxfer
        sourcep = vcdum1

        if(sinkp.gt.sourcep.and.sourcep.gt.1.e-8) then
        ratioagg= sourcep/sinkp
        pccxfer = pccxfer*ratioagg
        pacxfer = pacxfer*ratioagg
        ppcxfer = ppcxfer*ratioagg
        endif



 !.. (NEW) qagg is the mass tendency, nagg is the number tendency, vagg is the volume tendency, vcagg is the phi*V tendency
          qagg1 = -rselfloss*pprxfer  - rcrossloss*pcrxfer  - rcrossloss*parxfer

  !.. Do not over-deplete mass from aggregation
          sinkp = (alphvdum1/alphmdum1)*qagg1
          sourcep = Vadum1

          if(sinkp.gt.sourcep.and.sourcep.gt.1.e-8) then
             ratioagg= sourcep/sinkp
             pcrxfer = pcrxfer*ratioagg
             parxfer = parxfer*ratioagg
             pprxfer = pprxfer*ratioagg
          endif

          vagg1 = (alphvdum1/alphmdum1)*&
                  (-rselfloss*pprxfer  - rcrossloss*pcrxfer  - rcrossloss*parxfer)

          vcagg1 = -rselfloss*ppcxfer  - rcrossloss*pccxfer  - rcrossloss*pacxfer
          nagg1 = -nselfloss*ppenxfer - ncrossloss*pcenxfer - ncrossloss*paenxfer

          qagg2 = -rselfloss*ccrxfer  - rcrossloss*cprxfer  - rcrossloss*carxfer

          sinkc = (alphvdum2/alphmdum2)*qagg2
          sourcec = Vadum2

         if(sinkc.gt.sourcec.and.sourcec.gt.1.e-8) then
             ratioagg= sourcec/sinkc
             cprxfer = cprxfer*ratioagg
             carxfer = carxfer*ratioagg
             ccrxfer = ccrxfer*ratioagg
         endif

          vagg2 = (alphvdum2/alphmdum2)*&
                  (-rselfloss*ccrxfer  - rcrossloss*cprxfer  - rcrossloss*carxfer)

          vcagg2 = -rselfloss*cccxfer  - rcrossloss*cpcxfer  - rcrossloss*cacxfer
          nagg2 = -nselfloss*ccenxfer - ncrossloss*pcenxfer - ncrossloss*caenxfer

          qagg3 = -qagg1-qagg2
          vagg3 = -vagg1-vagg2
          vcagg3 = -vcagg1-vcagg2
         !nagg3 = -nagg1-nagg2-nselfgain*aaenxfer !! NOTE: make sure to fix 1/2 coeff. (and rhoairdum) for gains in other versions.

          nagg3 =nselfgain*(ppenxfer+ccenxfer+aaenxfer) + &
                 ncrossgain*(pcenxfer+paenxfer+caenxfer)

          nmon3 = -nagg1-nagg2

         ! if (abs(nagg3).gt.0.) then
         ! print*, 'nagg=',nagg3
         ! print*, 'nmon3=',nmon3
         ! stop
         ! endif

          qagg1 = min(qagg1,0.)
          qagg2 = min(qagg2,0.)
          qagg3 = max(qagg3,0.)

          vagg1 = min(vagg1,0.)
          vagg2 = min(vagg2,0.)
          vagg3 = max(vagg3,0.)

          vcagg1 = min(vcagg1,0.)
          vcagg2 = min(vcagg2,0.)
          vcagg3 = max(vcagg3,0.)

          ! Scale by air density for numbers?
          !nagg1 = nagg1/rhoairdum
          !nagg2 = nagg2/rhoairdum
          !nagg3 = nagg3/rhoairdum

          !nmon3 = nmon3/rhoairdum



!          nagg1 = min(nagg1,0.)
!          nagg2 = min(nagg2,0.)
!          nagg3 = max(nagg3,0.)



         ! if (qdum3.gt.1e-6) then
         ! print*, '-----------------'
         ! print*, 'Mass before:',qdum1,qdum2,qdum3
         ! print*, 'Number before:',ndum1,ndum2,ndum3
         ! print*, 'Volume before:',vadum1,vadum2,vadum3
         ! print*, 'Volume-C before:',vcdum1,vcdum2,vcdum3
         ! print*, 'Mass transfer rates: Plates:', 'pprxfer=',pprxfer
         ! print*, 'pcrxfer=',pcrxfer,'parxfer=',parxfer
         ! print*, 'Number transfer rates: Plates:'
         ! print*, 'ppenxfer=',ppenxfer,'pcenxfer=',pcenxfer
         ! print*, 'Mass transfer rates: Columns:','ccrxfer=',ccrxfer
         ! print*, 'cprxfer=',cprxfer, 'carxfer=',carxfer
         ! print*, 'Number transfer rates: Columns:'
         ! print*, 'ccenxfer=',ccenxfer,'caenxfer=',caenxfer
         ! print*, 'Aggregate transfer:','aaenxfer=',aaenxfer
         ! print*, 'total mass transfer:', qagg1+qagg2, qagg3
         ! print*, 'new agg mass:', qdum1+qagg1,qdum2+qagg2,qdum3+qagg3
         ! print*, 'new agg num:', ndum1+nagg1,ndum2+nagg2,ndum3+nagg3
         ! print*, 'fallspeed terms:',dvtm_dum, dvtn_dum, dvtc_dum
         ! stop

         ! endif

  end subroutine starbuck




!----------------------------------------------------------------------------
! Subroutine shape_diag
!
! ELD October 2019
!
! This subroutine uses new mass and volume mixing ratios to appropriately
! diagnose an, ds, and rhobar after aggregation. This is necessary because
! aggregation decreases both mass and volume, thus implicitly changing rhobar.
! However, this diagnosis is tricky because the integrated volume mixing
! ratios yield gamma function factors that incorporate deltastar. Therefore,
! deltastar is diagnosed by taking the natural logarithm of the integrated
! volumes and approximating loggamma factors using second order polynomials.
! an is diagnosed from the new deltastar and the integrated volume (Va).
! rhobar is more simply given as qi/Va.
!-----------------------------------------------------------------------------

 subroutine shape_diag(gammnu,ni,qi,Via,Vic,ao,ds,astar,alphv,betam,ai,ci,ani,cni)

 real, parameter :: qasmall=1.e-20

 real, intent(INOUT) :: ani, cni,ai,ci,astar,alphv,betam
 real, intent(IN) :: gammnu, ao, ni, qi, Via, Vic
 real, intent(INOUT) :: ds
 real :: p1, p2, p3, Vi0, dsold, aniold, cniold,aiold,ciold

 ! Solve for Vi0

 Vi0 = (4./3.)*PI*ao**3.*ni

 dsold = ds
 aniold = ani
 cniold =cni
 aiold = ai
 ciold = ci

 ! First, get appropriate polynomial coefficients
 if(nu.eq.1.) then
    p1 = 1.278
    p2 = -0.5793
    p3 = -0.6989
 elseif(nu.eq.4.) then

    p1 = 0.6888
    p2 = -0.536
    p3 = -0.1528

 else
 print*, 'input nu needs to be either 1 or 4 currently'
 endif

 p2 = p2+log(Via**2/(Vi0*Vic))
 p3 = p3+log(Vi0*Via/(Vic**2))


 !p2 = p2+2.*log(Via)-log(Vi0)-log(Vic)
 !p3 = p3+log(Vi0)+log(Via)-2.*log(Vic)

 !print*, 'p2=',p2
 !print*, 'p3=',p3
 !print*, 'Via=',log(Via)
 !print*, 'Vic=',log(Vic)
 !print*, 'ao=',ao
 !print*, 'nu=',nu

 ! Now, use quadratic formula to get new ds

 ds = (-p2-sqrt(p2*p2 - 4.*p1*p3))/(2.*p1)

 if (ds.lt.0..or.ds.gt.5.) then

  ds =(-p2+sqrt(p2*p2 - 4.*p1*p3))/(2.*p1)

 endif



 !print*,'ds=',ds

 ds = min(max(ds,0.55),1.5)

 !print*, 'ds_diag=', ds


 !print*, 'Viadum=',via
 !print*, 'ai_old=',ai
 !print*, 'ani_old=',ani

 betam = 2.+ds
 astar = ao**(1.-ds)
 alphv = (4./3.)*PI*astar
 ani = ao*((Via/Vi0)*(gammnu/gamma(nu+betam)))**(1./betam)
 cni =  astar*ani**ds
 ai = max(ani**2*cni*ni,qasmall) ! Volume-a (as described by an^2*cn)
 ci = max(ani*cni**2*ni,qasmall) ! Volume-c (as described by cn^2*an)



 !if(dsold.lt.0.99) then
 ! print*, 'dsold=',dsold
 ! print*, 'dsnew=',ds
 ! print*, 'aniold=',aniold
 ! print*, 'ani=',ani
 ! print*, 'cniold=',cniold
 ! print*, 'cni=',cni
 ! print*, 'aiold=',aiold
 ! print*, 'ai=',ai
 ! print*, 'ciold=',ciold
 ! print*, 'ci=',ci

  !print*, p2,p3
 ! print*, Via, Vic
 ! stop
 !endif

 !ai = max(NU*ni*ani,qasmall)
 !ci = max(NU*ni*cni,qasmall)




 end subroutine shape_diag
!------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Access and interpolate lookup table values (from Hugh Morrison)
  subroutine access_lookup_table(itabdum,dumjj,dumii,dumi,dumk,index, &
       dum1,dum2,dum4,dum5,proc)
    
    implicit none
    
    real :: itabdum(:,:,:,:,:), dum1,dum2,dum4,dum5,proc
    real :: dproc1,dproc2,iproc1,gproc1,tmp1,tmp2
    integer :: dumjj,dumii,dumi,dumk,index
    
  !.. get value at current density index
  !.. first interpolate for current rimed fraction index
    dproc1 =itabdum(dumjj,dumii,dumi,dumk,index)+(dum1-real(dumi))* &
         (itabdum(dumjj,dumii,dumi+1,dumk,index)-itabdum(dumjj,dumii,dumi,dumk,index))
    
    dproc2 =itabdum(dumjj,dumii,dumi,dumk+1,index)+(dum1-real(dumi))* &
         (itabdum(dumjj,dumii,dumi+1,dumk+1,index)-itabdum(dumjj,dumii,dumi,dumk+1,index))
    
    iproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)
    
  !.. linearly interpolate to get process rates for rimed fraction index + 1
    dproc1 =itabdum(dumjj,dumii+1,dumi,dumk,index)+(dum1-real(dumi))* &
         (itabdum(dumjj,dumii+1,dumi+1,dumk,index)- &
         itabdum(dumjj,dumii+1,dumi,dumk,index))
    
    dproc2 =itabdum(dumjj,dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))* &
         (itabdum(dumjj,dumii+1,dumi+1,dumk+1,index)- &
         itabdum(dumjj,dumii+1,dumi,dumk+1,index))
    
    gproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)
    
    tmp1=iproc1+(dum4-real(dumii))*(gproc1-iproc1)
    
  !.. get value at density index + 1
  !.. first interpolate for current rimed fraction index
    dproc1 =itabdum(dumjj+1,dumii,dumi,dumk,index)+(dum1-real(dumi))* &
         (itabdum(dumjj+1,dumii,dumi+1,dumk,index)- &
         itabdum(dumjj+1,dumii,dumi,dumk,index))
    
    dproc2 =itabdum(dumjj+1,dumii,dumi,dumk+1,index)+ &
         (dum1-real(dumi))*(itabdum(dumjj+1,dumii,dumi+1,dumk+1,index)- &
         itabdum(dumjj+1,dumii,dumi,dumk+1,index))
    
    iproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)
    
  !.. linearly interpolate to get process rates for rimed fraction index + 1
    dproc1 =itabdum(dumjj+1,dumii+1,dumi,dumk,index)+(dum1-real(dumi))* &
         (itabdum(dumjj+1,dumii+1,dumi+1,dumk,index)- &
         itabdum(dumjj+1,dumii+1,dumi,dumk,index))
    
    dproc2 =itabdum(dumjj+1,dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))* &
         (itabdum(dumjj+1,dumii+1,dumi+1,dumk+1,index)- &
         itabdum(dumjj+1,dumii+1,dumi,dumk+1,index))
    
    gproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)
    
    tmp2=iproc1+(dum4-real(dumii))*(gproc1-iproc1)
    
  !.. get final process rate
    proc=tmp1+(dum5-real(dumjj))*(tmp2-tmp1)
       
  end subroutine access_lookup_table

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Computes ln(gamma(xx)) using Lanczos                                
  real function gammln(xx)
    
    implicit none
    
    REAL :: xx
    REAL*8 :: cof(6), stp, x, y, tmp, ser 
    INTEGER :: j
    DATA cof, stp/76.18009172947146d0, &
         -86.50532032941677d0,         &
         24.01409824083091d0,          &
         -1.231739572450155d0,         &
         0.1208650973866179d-2,        &
         -0.5395239384953d-5,          &
         2.5066282746310005d0/
    x = xx
    y = x
    tmp = x + 5.5d0
    tmp = (x + 0.5d0) * log(tmp) - tmp
    ser = 1.000000000190015d0
    do j = 1, 6
       y = y + 1.0d0
       ser = ser + cof(j) / y
    enddo
    
    gammln = tmp + log(stp * ser / x)

    return
  end function gammln

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  REAL FUNCTION GAMMA(X)
    !----------------------------------------------------------------------
    !
    ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
    !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
    !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
    !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
    !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
    !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
    !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
    !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
    !   MACHINE-DEPENDENT CONSTANTS.
    !
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
    !
    ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
    ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
    ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
    !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
    !                  GAMMA(XBIG) = BETA**MAXEXP
    ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
    !          APPROXIMATELY BETA**MAXEXP
    ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    !          1.0+EPS .GT. 1.0
    ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
    !          1/XMININ IS MACHINE REPRESENTABLE
    !
    !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
    !
    !                            BETA       MAXEXP        XBIG
    !
    ! CRAY-1         (S.P.)        2         8191        966.961
    ! CYBER 180/855
    !   UNDER NOS    (S.P.)        2         1070        177.803
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (S.P.)        2          128        35.040
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (D.P.)        2         1024        171.624
    ! IBM 3033       (D.P.)       16           63        57.574
    ! VAX D-FORMAT   (D.P.)        2          127        34.844
    ! VAX G-FORMAT   (D.P.)        2         1023        171.489
    !
    !                            XINF         EPS        XMININ
    !
    ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
    ! CYBER 180/855
    !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
    ! IEEE (IBM/XT,
    !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
    ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
    ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
    ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! ERROR RETURNS
    !
    !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
    !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
    !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
    !
    !
    !  INTRINSIC FUNCTIONS REQUIRED ARE:
    !
    !     INT, DBLE, EXP, LOG, REAL, SIN
    !
    !
    ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
    !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
    !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
    !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
    !
    !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
    !              SONS, NEW YORK, 1968.
    !
    !  LATEST MODIFICATION: OCTOBER 12, 1989
    !
    !  AUTHORS: W. J. CODY AND L. STOLTZ
    !           APPLIED MATHEMATICS DIVISION
    !           ARGONNE NATIONAL LABORATORY
    !           ARGONNE, IL 60439
    !
    !----------------------------------------------------------------------
    implicit none
    INTEGER I,N
    LOGICAL PARITY
    REAL                                                          &
         CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
         TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
    REAL, DIMENSION(7) :: C
    REAL, DIMENSION(8) :: P
    REAL, DIMENSION(8) :: Q
    REAL, PARAMETER :: xxx = 0.9189385332046727417803297
    !----------------------------------------------------------------------
    !  MATHEMATICAL CONSTANTS
    !----------------------------------------------------------------------
    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/
    
    
    !----------------------------------------------------------------------
    !  MACHINE DEPENDENT PARAMETERS
    !----------------------------------------------------------------------
    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
    !----------------------------------------------------------------------
    !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
    !     APPROXIMATION OVER (1,2).
    !----------------------------------------------------------------------
    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
         -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
         8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
         -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
         -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
         2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
         -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
    !----------------------------------------------------------------------
    !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
    !----------------------------------------------------------------------
    DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
         -5.952379913043012E-04,7.93650793500350248E-04,	           &
         -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
         5.7083835261E-03/
    !----------------------------------------------------------------------
    !  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
    !----------------------------------------------------------------------
    CONV(I) = REAL(I)
    PARITY=.FALSE.
    FACT=ONE
    N=0
    Y=X
    IF(Y.LE.ZERO)THEN
       !----------------------------------------------------------------------
       !  ARGUMENT IS NEGATIVE
       !----------------------------------------------------------------------
       Y=-X
       Y1=AINT(Y)
       RES=Y-Y1
       IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
       ELSE
          RES=XINF
          GOTO 900
       ENDIF
    ENDIF
    !----------------------------------------------------------------------
    !  ARGUMENT IS POSITIVE
    !----------------------------------------------------------------------
    IF(Y.LT.EPS)THEN
       !----------------------------------------------------------------------
       !  ARGUMENT .LT. EPS
       !----------------------------------------------------------------------
       IF(Y.GE.XMININ)THEN
          RES=ONE/Y
       ELSE
          RES=XINF
          GOTO 900
       ENDIF
    ELSEIF(Y.LT.TWELVE)THEN
       Y1=Y
       IF(Y.LT.ONE)THEN
          !----------------------------------------------------------------------
          !  0.0 .LT. ARGUMENT .LT. 1.0
          !----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
       ELSE
          !----------------------------------------------------------------------
          !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
          !----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
       ENDIF
       !----------------------------------------------------------------------
       !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
       !----------------------------------------------------------------------
       XNUM=ZERO
       XDEN=ONE
       DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
       END DO
       RES=XNUM/XDEN+ONE
       IF(Y1.LT.Y)THEN
          !----------------------------------------------------------------------
          !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
          !----------------------------------------------------------------------
          RES=RES/Y1
       ELSEIF(Y1.GT.Y)THEN
          !----------------------------------------------------------------------
          !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
          !----------------------------------------------------------------------
          DO I=1,N
             RES=RES*Y
             Y=Y+ONE
          END DO
       ENDIF
    ELSE
       !----------------------------------------------------------------------
       !  EVALUATE FOR ARGUMENT .GE. 12.0,
       !----------------------------------------------------------------------
       IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
             SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+xxx
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
       ELSE
          RES=XINF
          GOTO 900
       ENDIF
    ENDIF
    !----------------------------------------------------------------------
    !  FINAL ADJUSTMENTS AND RETURN
    !----------------------------------------------------------------------
    IF(PARITY)RES=-RES
    IF(FACT.NE.ONE)RES=FACT/RES
900 GAMMA=RES
    RETURN
    ! ---------- LAST LINE OF GAMMA ----------
  END FUNCTION GAMMA

  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!
  !.. Get the inherent growth ratio data as a function of temperature
  real function get_igr(igrdatadum, temp)
    
    implicit none
    
    real :: igrdatadum(:), dum, igr1, igr2, temp
    
  !.. Inherent growth ratio (IGR) from data
  !.. See Lamb and Scott (1972) and Chen and Lamb (1994)
    if((temp-T0).ge.-59..and.(temp-T0).le.-1.) then
       dum = (abs(real(int(temp-T0))) + 1.) - abs(temp-T0)
       igr1 = igrdatadum(max((int(temp-T0)*(-1)),1))
       igr2 = igrdatadum(min(((int(temp-T0)*(-1))+1),60))
       get_igr = dum*igr1 + (1.-dum)*igr2
    elseif((temp-T0).gt.-1..and.(temp-T0).le.0.) then
       dum = 1.0 - abs(temp-T0)
       igr1 = 1.0
       igr2 = igrdatadum(1)
       get_igr = dum*igr1 + (1.-dum)*igr2
    elseif((temp-T0).ge.-60..and.(temp-T0).lt.-59.) then
       dum = 60.0 - abs(temp-T0)
       igr1 = igrdatadum(59)
       igr2 = igrdatadum(60)
       get_igr = dum*igr1 + (1.-dum)*igr2
    elseif((temp-T0).lt.-60.) then
       get_igr = igrdatadum(60)
    else
       get_igr = 1.
    endif
    
  end function get_igr



  !--------------------------------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------------!

end module module_mp_jensen_ishmael

  !.. Full stop

