module statevarsmod

use parametersmod, only : sp,dp,npft,ncvar,climbuf,ns,nl

implicit none

public :: initstatevars

!---------------------------------------------------------------
!module shared variables

!We put all of the state variables in a few large structures

type gridcell_state_vars

  !daily gridcell meteorological state variables

  real(sp), dimension(365) :: dtemp       !daily temperature (deg C)
  real(sp), dimension(365) :: dprec       !daily precipitation (mm)
  real(sp), dimension(365) :: dsunp       !daily sunshine (% full sunshine)
  real(sp), dimension(365) :: dpet        !daily potential evapotranspiration (mm)
  real(sp), dimension(365) :: dmelt       !daily snowmelt (mm)
  real(sp), dimension(365) :: snowpack    !daily snow water equivalent (mm)
  real(sp), dimension(365) :: dtemp_soil  !daily soil temperature (deg C)
  real(sp), dimension(365) :: dw1         !daily soil layer 1 water content
  real(sp), dimension(365) :: dwind       !daily mean wind speed (m s-1)
  real(sp), dimension(365) :: dlight      !daily lightning flashes

  !Gridcell state variables
  
  real(sp) :: grid_area   !gridcell area (m2)
  real(sp) :: slopeangle  !mean gridcell slope angle (radians)
  
  real(sp) :: mtemp_min20 !20-year average minimum monthly temperature (deg C)
  real(sp) :: gdd20       !20-year average growing degree days
  
  real(sp) :: lu_turnover !fraction of the gridcell that is cycled through each year by anthropogenic land use (1/years)

  !state variables describing the gridcell radiation

  real(sp) :: Pjj     !precipitation equitability index
  real(sp) :: Ratm    !relative atmospheric pressure 1=sea level

  real(dp) :: mbar    !daytime mean optical air mass (unitless, 1 at equatorial noon)
  real(dp) :: mo      !air mass at cosine zenith angle maximum
  real(dp) :: mc      !air mass at cosine zenith angle medium
  real(dp) :: ml      !air mass at cosine zenith angle bottom quarter range point

  real(dp) :: direct  !direct beam downwelling shortwave (kJ m-2 d-1)
  real(dp) :: diffuse !diffuse downwelling shortwave (kJ m-2 d-1)

  real(dp) :: lw_rad  !longwave radiation (kJ m-2 d-1)
  real(dp) :: lvap    !Latent heat of vaporization of water (temperature dependent) (kJ kg-1)
  real(dp) :: gamma   !psychrometer constant (Pa K-1)
  real(dp) :: ss      !rate of increase of saturated vapor pressure with temperature (desdT) (Pa K-1)

  !other gridcell state variables
  
  real(sp) :: PD        !human population density
  real(sp) :: cropfrac  !crop fraction

  !decadal gridcell state variables

  real(sp), dimension(climbuf) :: gdd_buf        !buffer to store 'climbuf' years of GDD totals
  real(sp), dimension(climbuf) :: mtemp_min_buf  !buffer to store 'climbuf' years of coldest month temperatures

  real(sp), dimension(12) :: temp0  !temperature of the previous year
  
  !anthropogenic wood usage pools
  
  real(sp) :: wood_fast
  real(sp) :: wood_slow
  real(sp) :: crop_harvest

  !anthropogenic wood decay flux
  
  real(sp), dimension(3) :: prod_flux  !product flux(3) = respiration of harvest

end type gridcell_state_vars

!--------------------------------------------

type state_variables

  !scalars
  
  real(sp) :: afire_frac    !fraction of gridcell burnt this year
  real(sp) :: k_fast_ave    !running average k_fast for subroutine littersom
  real(sp) :: k_slow_ave    !running average k_slow for subroutine littersom
  real(sp) :: litterC_bg
  real(sp) :: litterC_fast
  real(sp) :: litterC_slow
  real(sp) :: livebiomass   !sum of living plant carbon across all PFTs
  real(sp) :: snow0         !Dec 31 snow pack
  real(sp) :: coverfrac     !land use cover fraction
  real(sp) :: tilecarbon    !total carbon by land use tile (living + littler + SOM)
  real(sp) :: soilerosion   !total annual soilerosion by tile (g m-2)
  real(sp) :: NI            !Nesterov index (because this is cumulative, we have to track it separately for each tile)
  real(sp) :: erflux        !soil erosion carbon flux (g m-2)
  real(sp) :: forager_pd    !forager potential population density (persons 100km-2)

  !soil layer

  real(sp), dimension(2) :: zpos   !depth to soil layer midpoint
  real(sp), dimension(2) :: dz     !soil layer thickness (m)
  real(sp), dimension(2) :: sand   !mass %
  real(sp), dimension(2) :: silt   !mass %
  real(sp), dimension(2) :: clay   !mass %
  real(sp), dimension(2) :: OM     !mass %
  real(sp), dimension(2) :: OrgM   !organic matter (g m-2)
  real(sp), dimension(2) :: bulk   !bulk density (g cm-3)
  real(sp), dimension(2) :: Tsat   !saturation (volume fraction)
  real(sp), dimension(2) :: T33    !33 kPa (volume fraction)
  real(sp), dimension(2) :: T1500  !1500 kPa (volume fraction)
  real(sp), dimension(2) :: whc    !water holding capacity T33 - T1500 (volume fraction)
  real(sp), dimension(2) :: Ksat   !saturated conductivity (mm hr-1)
  real(sp), dimension(2) :: w      !instantaneous soil water content

  !npft  

  logical,  dimension(npft) :: present     !whether PFT present in gridcell
  logical,  dimension(npft) :: leafon
  integer,  dimension(npft) :: leafondays
  integer,  dimension(npft) :: leafoffdays
  real(sp), dimension(npft) :: crownarea    !crown area (m2)
  real(sp), dimension(npft) :: dwscal365    !daily water scalar for day 365
  real(sp), dimension(npft) :: fpc_grid     !gridcell foliar projective cover (FPC)
  real(sp), dimension(npft) :: fpc_inc      !increment (if +ve) in FPC since last year
  real(sp), dimension(npft) :: height       !tree height (m)
  real(sp), dimension(npft) :: lai_ind      !individual leaf area index
  real(sp), dimension(npft) :: nind         !gridcell individual density (indiv/m2)
  real(sp), dimension(npft) :: plant_carbon
  real(sp), dimension(12,npft) :: mlai	     !monthly mean LAI, per pft

  !ncvar
  
  real(sp), dimension(ncvar) :: acflux_conv       !C flux to atmosphere due to anthropogenic deforestation (gC/m2)
  real(sp), dimension(ncvar) :: acflux_estab      !annual biomass increment due to establishment (gC/m2)
  real(sp), dimension(ncvar) :: acflux_fire       !C flux to atmosphere due to fire (gC/m2)
  real(sp), dimension(ncvar) :: arh               !annual heterotrophic respiration (gC/m2)
  real(sp), dimension(ncvar) :: cpool_fast        !fast-decomposing soil C pool (gC/m2)
  real(sp), dimension(ncvar) :: cpool_slow        !slow-decomposing soil C pool (gC/m2)
  real(sp), dimension(ncvar) :: grid_npp          !gridcell total npp (gC/m2)
  real(sp), dimension(ncvar) :: litter_decom_ave  !running average litter_decom for subroutine littersom
  
  !npft,ncvar
  
  real(sp), dimension(npft,ncvar) :: anpp           !annual gridcell NPP (gC/m2)
  real(sp), dimension(npft,ncvar) :: lm_ind         !individual leaf mass (gC)
  real(sp), dimension(npft,ncvar) :: sm_ind         !individual sapwood mass (gC)
  real(sp), dimension(npft,ncvar) :: hm_ind         !individual heartwood mass (gC)
  real(sp), dimension(npft,ncvar) :: rm_ind         !individual fine root mass (gC)
  real(sp), dimension(npft,ncvar) :: litter_ag_fast !gridcell above-ground litter (gC/m2)
  real(sp), dimension(npft,ncvar) :: litter_ag_slow !gridcell above-ground litter (gC/m2)
  real(sp), dimension(npft,ncvar) :: litter_bg      !gridcell below-ground litter (gC/m2)

end type state_variables

!----------

type output_variables
  real(sp) :: livebiomass
  real(sp) :: litterC_fast
  real(sp) :: litterC_slow
  real(sp) :: litterC_bg
  real(sp) :: soilC_fast
  real(sp) :: soilC_slow
  real(sp) :: NPP
  real(sp) :: NEP
  real(sp) :: NBP
  real(sp) :: AET
  real(sp) :: PET
  real(sp) :: aprec
  real(sp) :: soilerosion
end type output_variables

!--------------------------------------------
!type declarations

type(gridcell_state_vars), save, allocatable, target, dimension(:)   :: gsv
type(state_variables),     save, allocatable, target, dimension(:,:) :: sv
type(output_variables),    save, allocatable, target, dimension(:)   :: ov

real(sp), dimension(ncvar) :: co2

!--------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine initstatevars(bufsize,ntiles)

use iovariablesmod, only : lu_turn_yrs

implicit none

integer, intent(in) :: bufsize  !the number of gridcells that will be run by each worker
integer, intent(in) :: ntiles   !the number of (sub-gridcell) tiles that will be run on each grid

integer :: i
integer :: j

!---------

allocate(gsv(bufsize))
allocate(sv(ntiles,bufsize))
allocate(ov(bufsize))

!---------

if (lu_turn_yrs > 0.) then
  gsv%lu_turnover = 1. /  lu_turn_yrs
else
  gsv%lu_turnover = 0.
end if

gsv%gdd20       = 0.
gsv%mtemp_min20 = 0.

gsv%wood_fast    = 0.
gsv%wood_slow    = 0.
gsv%prod_flux(1) = 0.
gsv%prod_flux(2) = 0.

forall (i=1:climbuf)
  gsv%gdd_buf(i) = -9999.
  gsv%mtemp_min_buf(i) = -9999.
end forall

sv%w(1) = 1.
sv%w(2) = 1.
sv%snow0 = 0.

sv%k_fast_ave  = 0.
sv%k_slow_ave  = 0.

forall (i=1:npft)
  sv%dwscal365(i)    = 1.0
  sv%crownarea(i)    = 0.
  sv%fpc_grid(i)     = 0.
  sv%fpc_inc(i)      = 0.
  sv%height(i)       = 2.0
  sv%lai_ind(i)      = 0.
  sv%nind(i)         = 0.
  sv%plant_carbon(i) = 0.
  
  forall (j=1:ncvar)
    sv%hm_ind(i,j) = 0.
    sv%lm_ind(i,j) = 0.
    sv%sm_ind(i,j) = 0.
    sv%rm_ind(i,j) = 0.
    
    sv%litter_ag_fast(i,j) = 0.
    sv%litter_ag_slow(i,j) = 0.
    sv%litter_bg(i,j) = 0.

    sv%anpp(i,j) = 0.

  end forall
end forall

forall (j=1:ncvar)
  sv%grid_npp(j)     = 0.
  sv%acflux_estab(j) = 0.
  sv%acflux_fire(j)  = 0.
  sv%acflux_conv(j)  = 0.
  sv%arh(j)          = 0.

  sv%cpool_fast(j)       = 0.   
  sv%cpool_slow(j)       = 0.
  sv%litter_decom_ave(j) = 0.
end forall

sv%livebiomass  = 0.
sv%litterC_fast = 0.
sv%litterC_slow = 0.
sv%litterC_bg   = 0.

sv%forager_pd   = 0.

sv%coverfrac      = 0.   !initialize all categories to zero
sv(1,:)%coverfrac = 1.   !set natural vegetation to 100% cover

end subroutine initstatevars

!-------------------------------------------------------------------------------

end module statevarsmod
