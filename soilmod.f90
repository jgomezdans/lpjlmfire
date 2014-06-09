module newsoilmod

use parametersmod, only : sp,long

implicit none

public :: initializesoil
public :: updatesoil
public :: rusle

contains

!---------------------------------------------------------------------------

subroutine initializesoil(ncells)
!FLAG this subroutine is currently not used
use statevars,   only : gsv,sv
use iovariables, only : cellindex,soil
use lpjparams,   only : ntiles
use ptfmod,      only : fbulk,calctheta,fKsat

implicit none

integer(long), intent(in) :: ncells

!values from input dataset (n-layers)

real(sp), pointer, dimension(:) :: sand_in
real(sp), pointer, dimension(:) :: clay_in

!calculated values (scalars)

real(sp), pointer :: zpos
real(sp), pointer :: dz
real(sp), pointer :: sand
real(sp), pointer :: silt
real(sp), pointer :: clay
real(sp), pointer :: OM
real(sp), pointer :: bulk
real(sp), pointer :: Tsat
real(sp), pointer :: T33
real(sp), pointer :: T1500
real(sp), pointer :: whc
real(sp), pointer :: Ksat
real(sp), pointer :: slopeangle

!local variables

integer :: i
integer :: j
integer :: x
integer :: y
integer :: l
integer :: nl

real(sp) :: blk0
real(sp) :: orgC

integer :: soiltype

real(sp), dimension(2) :: sand2
real(sp), dimension(2) :: clay2

integer :: it

!-----------------------

nl = size(sv(1,1)%dz)

soiltype = 1  !set to typical for now

do j = 1,ncells

  x = cellindex(j,1)
  y = cellindex(j,2)

  !get mineral values from original input soil dataset

  sand_in => soil(x,y)%sand
  clay_in => soil(x,y)%clay
   
  !aggregate these values to two layers
  
  !no soil
  
  if (.not.(any(sand_in >= 0.))) cycle

  !topsoil
    
  sand2(1) = sum(sand_in(1:2),mask=sand_in(1:2) >= 0.) / count(sand_in(1:2) >= 0.)
  clay2(1) = sum(clay_in(1:2),mask=clay_in(1:2) >= 0.) / count(clay_in(1:2) >= 0.)

  !no subsoil, for now just set values equal to topsoil

  if (any(sand_in(3:5) >= 0.)) then
    sand2(2) = sum(sand_in(3:5),mask=sand_in(3:5) >= 0.) / count(sand_in(3:5) >= 0.)
    clay2(2) = sum(clay_in(3:5),mask=clay_in(3:5) >= 0.) / count(clay_in(3:5) >= 0.)
  else
    sand2(2) = sand2(1)
    clay2(2) = clay2(1)
  end if
  
  !transfer slope angle to statevars
  
  gsv(j)%slopeangle = soil(x,y)%slopeangle
  
  do i = 1,ntiles
    
    !loop over layers
    
    do l = 1,nl
      
      zpos  => sv(i,j)%zpos(l)
      dz    => sv(i,j)%dz(l)
      sand  => sv(i,j)%sand(l)
      silt  => sv(i,j)%silt(l)
      clay  => sv(i,j)%clay(l)
      OM    => sv(i,j)%OM(l)
      bulk  => sv(i,j)%bulk(l)
      Tsat  => sv(i,j)%Tsat(l)
      T33   => sv(i,j)%T33(l)
      T1500 => sv(i,j)%T1500(l)
      whc   => sv(i,j)%whc(l)
      Ksat  => sv(i,j)%Ksat(l)
      
      if (l == 1) then
        zpos = 0.15
        dz   = 0.3
      else
        zpos = 1.65
        dz   = 2.7
      end if
      
      sand = sand2(l) * 100.
      clay = clay2(l) * 100.
      silt = 100. - (sand + clay)
      OM   = 0.
            
      !initialize soil physical properties theta(x) and Ksat for mineral soil with no organic matter
      !because bulk still depends weakly on T1500, iterate to a stable solution for bulk

      T1500 = 0.1  !initial guess value for wilting point

      blk0 = fbulk(0.,T1500*100,clay,zpos,silt)  !bulk density of mineral soil (g cm-3)
      
      it = 1

      do
         
        if (it > 50) then
          exit
        else
          it = it + 1
        end if

        call calctheta(sand/100.,clay/100.,OM/100.,blk0,Tsat,T33,T1500,soiltype)

        !recalculate bulk

        bulk = fbulk(orgC,T1500*100.,clay,zpos,silt)  !units (g cm-3)

        if (abs(bulk - blk0) < 0.001) exit

        blk0 = bulk

      end do

      !with the final value for bulk density, recalculate porosity

      call calctheta(sand/100.,clay/100.,OM/100.,bulk,Tsat,T33,T1500)
      
      !calculate layer-integrated WHC
      
      whc = 1000. * dz * (T33 - T1500)  !mm

      !calculate saturated conductivity

      Ksat = fKsat(Tsat,T33,T1500)

    end do  !layer
    
  end do    !tile
end do      !grid

end subroutine initializesoil

!---------------------------------------------------------------------------

subroutine updatesoil(i,j)

use statevars, only : sv
use ptfmod,    only : fbulk,calctheta,fKsat,ombd,omcf

implicit none

!arguments

integer, intent(in) :: i
integer, intent(in) :: j

!parameter

real(sp) :: orgz = 1. / (ombd * 100.**3)  !inverse of density of organic matter (m3 g-1) multiply by OM g m-2 to get dz

!pointers

real(sp), pointer :: cpool_fast
real(sp), pointer :: cpool_slow

real(sp), pointer :: sand
real(sp), pointer :: silt
real(sp), pointer :: clay
real(sp), pointer :: OM    !(mass %)
real(sp), pointer :: bulk
real(sp), pointer :: Tsat
real(sp), pointer :: T33
real(sp), pointer :: T1500
real(sp), pointer :: whc
real(sp), pointer :: Ksat

real(sp), pointer, dimension(:) :: dz
real(sp), pointer, dimension(:) :: zpos
real(sp), pointer, dimension(:) :: OrgM  !(g m-2)

!local variables

integer :: l
integer :: nl = size(sv(1,1)%dz)
integer :: soiltype

real(sp) :: Csoil     !(g m-2)
real(sp) :: orgC      !(mass %)
real(sp) :: soilmass  !(g m-3)
real(sp) :: blk0      !(g m-3)
real(sp) :: dOM       !change in organic matter (g m-2)
real(sp) :: dzOM      !interlayer transport of SOM (g m-2)
real(sp) :: dzx       !excess change in top layer thickness

integer :: it

!-----------------------------------------------------------
!point pointers

dz   => sv(i,j)%dz
zpos => sv(i,j)%zpos
OrgM => sv(i,j)%OrgM
OM   => sv(i,j)%OM(1)

!--removed-- all SOM comes into top layer and is transported down column based on fixed thickness for top layer (<=30 cm)
!all fast SOM goes into top layer, all slow SOM goes into lower layer (try to rebalance SOM distribution a bit better)

cpool_fast => sv(i,j)%cpool_fast(1)
cpool_slow => sv(i,j)%cpool_slow(1)

do l = 1,nl

  bulk => sv(i,j)%bulk(l)

  !topsoil mass before addition of new carbon

  soilmass = bulk * 1.e6 * dz(l)     !(g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g

  !mass of new organic matter
  
  if (l == 1) then
    Csoil = cpool_fast
  else
    Csoil = cpool_slow  !total soil C (g m-2)
  end if

  !write(0,*)
  !write(0,*)'old bulk',bulk
  !write(0,*)'cpool total',Csoil

  dOM = Csoil * omcf - OrgM(l)  !total soil change in organic matter (g m-2)

  !since the formulation with all SOM going to top layer seems to be a bit unbalanced to the top layer,
  !send 10% of the new OM directly to the lower layer, e.g., for deep fine root turnover (could do this based on root vert dist).

  OrgM(l) = OrgM(l) + dOM

  !write(0,*)'dOM',dOM

  !mass weighted average bulk density of combined new carbon with previous soil

  bulk = (soilmass * bulk + dOM * ombd) / (soilmass + dOM)  !(g cm-3)

  !write(0,*)'new bulk',bulk

  !update layer thickness

  !dz(1) = (soilmass + dOM) / (1.e6 * bulk)  !g m-3 / g m-2 = m
  dz(1) = dz(1) + dOM / (1.e6 * ombd)  !g m-2 / g m-3 = m

  !write(0,*)'new dz1',dz(1)

  !update soil mass

  soilmass = bulk * 1.e6 * dz(1)     !(g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g
  
end do

!excess layer thickness

dzx   = dz(1) - 0.3
dz(1) = 0.3
dz(2) = max(dz(2) + dzx,0.1)  !FLAG fixed for now so that the minumum subsoil depth is 10 cm so we avoid negative soil depth

zpos(1) = 0.15
zpos(2) = dz(1) + 0.5 * dz(2)

!mass transfer of OM from top layer to bottom

OM = 100. * OrgM(1) / soilmass  !mass %

dzOM = dzx * bulk * 1.e6 * OM * 0.01 !mass (g m-2)

OrgM(1) = OrgM(1) - dzOM
OrgM(2) = OrgM(2) + dzOM

!write(0,*)'dz  ',dz,dzx
!write(0,*)'OrgM',OrgM,dzOM

do l = 1,nl

  !point pointers

  sand  => sv(i,j)%sand(l)
  silt  => sv(i,j)%silt(l)
  clay  => sv(i,j)%clay(l)
  OM    => sv(i,j)%OM(l)
  bulk  => sv(i,j)%bulk(l)
  Tsat  => sv(i,j)%Tsat(l)
  T33   => sv(i,j)%T33(l)
  T1500 => sv(i,j)%T1500(l)
  whc   => sv(i,j)%whc(l)
  Ksat  => sv(i,j)%Ksat(l)

  Csoil = OrgM(l) / omcf  !layer
  
  !write(0,*)'lyr ',l,' tile',i
  !write(0,*)'dz  ',dz(l)
  !write(0,*)'Csol',Csoil

  !because bulk density depends strongly on organic matter content and
  !weakly on wilting point water content, we iterate to a stable solution for bulk density
  
  it = 1

  do
    
    if (it > 50) then
      exit
    else
      it = it + 1
    end if

    !convert SOM to mass fraction and calculate the difference in OM content

    soilmass = bulk * 1.e6 * dz(l)     !(g) (m-2)  note: g cm-3 * (m*100=cm) 100cm * 100cm = g

    orgC = max(100. * Csoil / soilmass,0.)  !(mass %)

    !recalculate bulk

    blk0 = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  !units (g cm-3)

    !calculate wilting point, field capacity, and saturation, needs input in fractions not percent

    OM = orgC * omcf             !(mass %)

    if (OM >= 30.) soiltype = 3  !humic soil

    call calctheta(sand/100.,clay/100.,OM/100.,blk0,Tsat,T33,T1500,soiltype)

    !recalculate bulk

    bulk = fbulk(orgC,T1500*100.,clay,zpos(l),silt)  !units (g cm-3)

    if (abs(bulk - blk0) < 0.001) exit

    blk0 = bulk

  end do

  !with the final value for bulk density, recalculate porosity

  call calctheta(sand/100.,clay/100.,OM/100.,bulk,Tsat,T33,T1500)

  !update layer-integrated WHC
      
  whc = 1000. * dz(l) * (T33 - T1500)  !mm

  !calculate saturated conductivity

  Ksat = fKsat(Tsat,T33,T1500)

  !write(0,*)'sand',sand
  !write(0,*)'silt',silt
  !write(0,*)'clay',clay
  !write(0,*)'OM  ',OM
  !write(0,*)'bulk',bulk

  !write(0,*)'Tsat',Tsat
  !write(0,*)'T33 ',T33
  !write(0,*)'Twp ',T1500
  !write(0,*)'Ksat',Ksat
  
end do  !layers

!read(*,*)

end subroutine updatesoil

!---------------------------------------------------------------------------

subroutine rusle(i,j)

use statevars, only : gsv,sv,ov
use ptfmod,    only : omcf

implicit none

!arguments

integer, intent(in) :: i  !tile
integer, intent(in) :: j  !gridcell index

!pointers

real(sp), pointer :: beta  !slope angle (radians)
real(sp), pointer :: sand  !(mass %)
real(sp), pointer :: silt  !(mass %)
real(sp), pointer :: clay  !(mass %)
real(sp), pointer :: OM    !(mass %)
real(sp), pointer :: Pann  !total annual precipitation (mm)
real(sp), pointer :: bulk

real(sp), pointer, dimension(:) :: dz
real(sp), pointer, dimension(:) :: zpos
real(sp), pointer, dimension(:) :: fpc_grid  !(fraction)

real(sp), pointer :: cpool_fast
real(sp), pointer :: cpool_slow

real(sp), pointer :: erosion !g m-2 = A * 1e6 (g ton-1) / 100**2 (m2 ha-1) = 100
real(sp), pointer :: erflux  !g m-2 Carbon

!parameters

real(sp), parameter :: dclay = log10(sqrt(1.e-7))  !=-3.5  factors from Torri et al, 1997, 1998
real(sp), parameter :: dsilt = log10(sqrt(1.e-4))  !=-2.0
real(sp), parameter :: dsand = log10(sqrt(1.e-1))  !=-0.5

!local variables

real(sp) :: fsand   !(mass fraction)
real(sp) :: fsilt   !(mass fraction)
real(sp) :: fclay   !(mass fraction)
real(sp) :: vegfrac !vegetation cover fraction

real(sp) :: A  !annual potential soil erosion (ton ha-1 yr-1)
real(sp) :: R  !rainfall erosivity factor (MJ mm ha-1 yr-1)
real(sp) :: L  !slope length factor
real(sp) :: S  !slope steepness factor
real(sp) :: K  !soil erodibility factor (ton ha h ha-1 MJ-1 mm-1)
real(sp) :: C  !land cover and management factor
real(sp) :: P  !conservation practice factor

real(sp) :: F       !intermediate calculation variable
real(sp) :: m       !slope length exponent
real(sp) :: Dg      !intermediate calculation variable
real(sp) :: lambda  !slope length (m)

real(sp) :: dzx     !excess change in top layer thickness
real(sp) :: dze     !thickness of eroded soil

real(sp) :: erodedC

!write(0,*)'rusle'

!-------------
!point pointers

Pann => ov(j)%aprec

beta     => gsv(j)%slopeangle
fpc_grid => sv(i,j)%fpc_grid
erosion  => sv(i,j)%soilerosion
erflux   => sv(i,j)%erflux

dz   => sv(i,j)%dz
zpos => sv(i,j)%zpos
bulk => sv(i,j)%bulk(1)

!soil erosion will always happen from the top layer of the soil, so only specify layer 1 here

sand => sv(i,j)%sand(1)
silt => sv(i,j)%silt(1)
clay => sv(i,j)%clay(1)
OM   => sv(i,j)%OM(1)

fsand = 0.01 * sand
fsilt = 0.01 * silt
fclay = 0.01 * clay

!-------------
!slope length and steepness

F = (sin(beta) / 0.0896) / (3. * sin(beta)**0.8 + 0.56)  !eqn. 3

m = F / (1. + F)  !eqn. 3

lambda = 100.  !fixed for now (m)

L = (lambda / 22.13)**m  !eqn. 2

if (tan(beta) < 0.09) then
  S = 10.8 * sin(beta) + 0.03  !eqn. 4
else
  S = 16.8 * sin(beta) - 0.50
end if

!soil erodibility

Dg = dsand * fsand + dsilt * fsilt + dclay * fclay  !eqn. 6

K = 0.0293 * (0.65 - Dg + 0.24 * Dg**2) * exp(-0.0021 * (OM / fclay) - 0.00037 * (OM / fclay)**2 - 4.02 * fclay + 1.72 * fclay**2)  !eqn. 5

!rainfall erosivity

if (Pann > 850.) then
  R = 587.8 - 1.219 * Pann + 0.004105 * Pann**2  !eqn. 8
else
  R = 0.0483 * Pann**1.610
end if

!land use and conservation practice (table 1)

if (i == 2) then !crop or pasture tile (assuming that plowing or animals have a greater impact on soil than natural veg.)

  C = 0.5
  P = 0.5

else

  vegfrac = sum(fpc_grid)
  
  C = 0.35 - 0.349 * vegfrac  !ranges from 0.35 with bare ground to 0.001 with full veg cover
  P = 1.0

end if

!final calculation

A = R * L * S * K * C * P  !eqn. 1

!convert A to g m-2

erosion = A * 100.

!calculate thickness of soil removed

dze = erosion / (1.e6 * bulk)  !g m-2 / g m-3 = m

dz(1) = dz(1) - dze

!repartition soil layers

dzx   = dz(1) - 0.3
dz(1) = 0.3
dz(2) = max(dz(2) + dzx,0.1)  !FLAG fixed for now so that the minumum subsoil depth is 10 cm (just in case so much would erode)

zpos(1) = 0.15
zpos(2) = dz(1) + 0.5 * dz(2)

!oxidation of SOM (20%) following Yang et al 2003

cpool_fast => sv(i,j)%cpool_fast(1)
cpool_slow => sv(i,j)%cpool_slow(1)

erodedC = erosion * OM * 0.01 / omcf  !eroded mass * om frac = g m-2 C

!20% of the soil carbon that eroded is oxidized, this is the remaining fraction assumed redeposited on the same landscape

erflux = 0.2 * erosion  !carbon flux to atmosphere from erosion

cpool_fast = max(cpool_fast - 0.5 * erflux,0.)
cpool_slow = max(cpool_slow - 0.5 * erflux,0.)

!----------------------

!write(0,*)'R ',R
!write(0,*)'L ',L
!write(0,*)'S ',S
!write(0,*)'Dg',Dg
!write(0,*)'K ',K
!write(0,*)'C ',C
!write(0,*)'P ',P
!write(0,*)'A ',A
!write(0,*)'erosion ',erosion
!write(0,*)'dz  ',dz

end subroutine rusle

!---------------------------------------------------------------------------

end module newsoilmod

