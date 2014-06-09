c     SUBROUTINE PFTPARAMETERS
c     Assignment of PFT-specific parameters and bioclimatic limits
c     Definition of initial sapling and grass mass structure
c     Calculation of maximum crown area for woody PFTs

! Vegetation is represented by a combination of the following plant functional
! types (PFTs)

!       1. tropical broadleaved evergreen tree
!       2. tropical broadleaved raingreen tree
!       3. temperate needleleaved evergreen tree
!       4. temperate broadleaved evergreen tree
!       5. temperate broadleaved summergreen tree
!       6. boreal needleleaved evergreen tree
!       7. boreal summergreen tree
!       8. C3 perennial grass
!       9. C4 perennial grass

      subroutine pftparameters(pftpar,sla,tree,evergreen,
     *  summergreen,raingreen,needle,boreal,lm_sapl,sm_sapl,hm_sapl,
     *  rm_sapl,latosa,allom1,allom2,allom3,
     *  wooddens,reinickerp,co2)

      implicit none

c     PARAMETERS:
      integer npft,npftpar,ncvar
        parameter (npft=9,npftpar=32,ncvar=3)
      real pi
        parameter (pi=3.14159265)

c     ARGUMENTS:
      real pftpar(1:npft,1:npftpar),sla(1:npft),co2(1:3)
      logical tree(1:npft),evergreen(1:npft)
      logical summergreen(1:npft),raingreen(1:npft),needle(1:npft)
      logical boreal(1:npft)
      real lm_sapl(1:npft,1:ncvar),sm_sapl(1:npft,1:ncvar)
      real rm_sapl(1:npft,1:ncvar),hm_sapl(1:npft,1:ncvar)
      real latosa,wooddens,reinickerp
      real allom1,allom2,allom3

c     LOCAL VARIABLES:
      integer n,pft
      real table(1:npft,1:npftpar)
      real lai_sapl       !sapling or initial grass LAI
      real x
      real lmtorm         !non-waterstressed leafmass to rootmass ratio
      real stemdiam       !sapling stem diameter
      real height_sapl    !sapling height

c-----------------------------------------------------------------------------
 
c     PFT PARAMETERS
 
c      1  fraction of roots in upper soil layer
c      2  plants with C4 (1) or C3 (0) photosynthetic pathway
c      3  water scalar value at which leaves shed by drought deciduous PFT
c      4  canopy conductance component (gmin, mm/s) not associated with
c         photosynthesis (Haxeltine & Prentice 1996, Table 4)
c      5  maintenance respiration coefficient
c      6  flammability threshold
c      7  maximum foliar N content (mg/g)
c         (Haxeltine & Prentice 1996a, Fig 4)
c      8  fire resistance index
c      9  leaf turnover period (years)
c     10  leaf longevity (years)
c     11  sapwood turnover period (sapwood converted to heartwood) (years)
c     12  root turnover period (years)
c     13  leaf C:N mass ratio
c     14  sapwood C:N mass ratio
c     15  root C:N mass ratio
c     16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
c     17  phenology type: evergreen (1), summergreen (2), raingreen (3),
c         any type (4) 
c     18  leaf to root ratio under non-water stressed conditions
c     19  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
c     20  tree maximum crown area (m2)
c     21  sapling (or grass on initialisation) LAI
c     22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
c     23  boreal pft (1), non-boreal pft (0)     
c     24  low temperature limit for CO2 uptake
c     25  lower range of temperature optimum for photosynthesis
c     26  upper range of temperature optimum for photosynthesis
c     27  high temperature limit for CO2 unptake
  
c     BIOCLIMATIC LIMITS
 
c     28 minimum coldest monthly mean temperature
c     29 maximum coldest monthly mean temperature
c     30 minimum growing degree days (at or above 5 deg C)
c     31 upper limit of temperature of the warmest month 
c     32 lower limit of growth efficiency (g/m2)
c

      data ((table(pft,n),n=1,8),pft=1,npft) /

c     ---------------------------------------------------------------------
c          1      2      3      4     5      6      7      8          PFT  !original LPJ values for respcoef, replaced with Wania values
c     ---------------------------------------------------------------------

     *  0.85,   0.0,  0.00,   0.5,  0.1,  0.15, 100.0,  0.12,       !  1  0.20,
     *  0.70,   0.0,  0.35,   0.5,  0.1,  0.15, 100.0,  0.50,       !  2  0.20,
     *  0.70,   0.0,  0.00,   0.3,  1.0,  0.15, 100.0,  0.12,       !  3  1.20,
     *  0.70,   0.0,  0.00,   0.5,  1.0,  0.15, 100.0,  0.50,       !  4  1.20,
     *  0.80,   0.0,  0.00,   0.5,  1.0,  0.15, 120.0,  0.12,       !  5  1.20,
     *  0.90,   0.0,  0.00,   0.3,  1.2,  0.15, 100.0,  0.12,       !  6  1.20,
     *  0.90,   0.0,  0.00,   0.3,  1.2,  0.15, 100.0,  0.12,       !  7  1.20,
     *  0.90,   0.0,  0.35,   0.5,  0.7,  0.15, 100.0,  1.00,       !  8  1.20,
     *  0.90,   1.0,  0.35,   0.5,  0.2,  0.15, 100.0,  1.00/       !  9  1.20,

      data ((table(pft,n),n=9,17),pft=1,npft) /

c     ---------------------------------------------------------------------
c          9     10     11     12     13     14     15     16     17   PFT
c     ---------------------------------------------------------------------

     *   2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  1
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, !  2
     *   2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, !  3
     *   1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  4
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, !  5
     *   2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   2.0,   1.0, !  6
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, !  7
     *   1.0,  1.00,   1.0,   2.0,  29.0,   0.0,  29.0,   3.0,   4.0, !  8
     *   1.0,  1.00,   1.0,   2.0,  29.0,   0.0,  29.0,   3.0,   4.0/ !  9 

      data ((table(pft,n),n=18,23),pft=1,npft) /

c     ---------------------------------------------------
c           18      19     20      21     22     23     PFT
c     ---------------------------------------------------

     *     1.0,  1000.0,  30.0,  4.000,  3.0,    0.0,    !  1
     *     1.0,  1000.0,  30.0,  4.000,  3.0,    0.0,    !  2
     *     1.0,  1000.0,  30.0,  4.000,  3.0,    0.0,    !  3
     *     1.0,  1000.0,  30.0,  4.000,  3.0,    0.0,    !  4
     *     1.0,   200.0,  30.0,  4.000,  3.0,    0.0,    !  5
     *     1.0,  1000.0,  30.0,  4.000,  3.0,    1.0,    !  6
     *     1.0,   200.0,  30.0,  4.000,  3.0,    1.0,    !  7
     *    0.75,   100.0,   0.0,  0.001,  3.0,    1.0,    !  8
     *    0.75,   100.0,   0.0,  0.001,  3.0,    0.0/    !  9  

      data ((table(pft,n),n=24,27),pft=1,npft) /
c     -------------------------------------
c          24     25     26      27    PFT
c     -------------------------------------
     *    2.0,   25.0,  30.0,   55.0, ! 1
     *    2.0,   25.0,  30.0,   55.0, ! 2
     *   -4.0,   20.0,  30.0,   42.0, ! 3
     *   -4.0,   20.0,  30.0,   42.0, ! 4
     *   -4.0,   20.0,  25.0,   38.0, ! 5
     *   -4.0,   15.0,  25.0,   38.0, ! 6
     *   -4.0,   15.0,  25.0,   38.0, ! 7
     *   -4.0,   10.0,  30.0,   45.0, ! 8
     *    6.0,   20.0,  45.0,   55.0/ ! 9     

      data ((table(pft,n),n=28,npftpar),pft=1,npft) /
      
c     --------------------------------------------------------
c          28       29       30       31       32       PFT
c     --------------------------------------------------------
     *    15.5,  1000.0,    0.0,    1000.0,    0.0,  !  1
     *    15.5,  1000.0,    0.0,    1000.0,    0.0,  !  2
     *    -2.0,    22.0,  900.0,    1000.0,    0.0,  !  3 
     *     3.0,    18.8, 1200.0,    1000.0,    0.0,  !  4
     *   -17.0,    15.5, 1200.0,    1000.0,    0.0,  !  5
     *   -32.5,    -2.0,  600.0,      23.0,    0.0,  !  6
     * -1000.0,    -2.0,  350.0,      23.0,    0.0,  !  7
     * -1000.0,    15.5,    0.0,    1000.0,    0.0,  !  8
     *    15.5,  1000.0,    0.0,    1000.0,    0.0/  !  9
c----------------------------------------------------------------------------

      do pft=1,npft

c       Transfer parameter values to array pftpar

        do n=1,npftpar
          pftpar(pft,n)=table(pft,n)
        enddo
        

c       Assign leaf and phenology logicals

        if (pftpar(pft,16).le.2.0) then
          tree(pft)=.true.
          if (pftpar(pft,16).eq.2.0) then
            needle(pft)=.true.
          else
            needle(pft)=.false.
          endif
        else
          tree(pft)=.false.
          needle(pft)=.false.
        endif

        if (pftpar(pft,17).eq.1.0) then
          evergreen(pft)=.true.
          summergreen(pft)=.false.
          raingreen(pft)=.false.
        elseif (pftpar(pft,17).eq.2.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.true.
          raingreen(pft)=.false.
        elseif (pftpar(pft,17).eq.3.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.false.
          raingreen(pft)=.true.
        else
          evergreen(pft)=.true.
          summergreen(pft)=.true.
          raingreen(pft)=.true.
        endif
        
        if (pftpar(pft,23).eq.1.0) then
          boreal(pft)=.true.
        else
          boreal(pft)=.false.
        endif            

c       Calculate specific leaf area (SLA) for each PFT from leaf longevity
c       Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
c       Equation based on Reich et al 1997, Fig 1f:

c       SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

c       SLA in m2/gC, leaf_longevity in years

        sla(pft)=2.e-4*exp(6.15-0.46*log(pftpar(pft,10)*12.))

c       Define initial mass structure

        lai_sapl=pftpar(pft,21)

        if (tree(pft)) then  !woody PFTs

c         Calculate leafmass for a sapling individual
c          (1) lai = leafmass * sla / (crown area)
c          (2) (leaf area) = latosa * (sapwood xs area)
c                 (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c          (3) (crown area) = allom1 * (stem diameter) ** reinickerp
c                 (Reinickes theory)
c         From (1),
c          (4) leafmass = lai * (crown area) / sla
c         From (1) & (3),
c          (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
c         From (2),
c          (6) leafmass = latosa * (sapwood xs area) / sla
c          (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
c         From (6) and (7),
c          (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
c         From (8),
c          (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
c         (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
c         Define x,
c         (11) x = [ (sapwood diameter)+(heartwood diameter) ] / 
c                  (sapwood diameter)
c         From (10) & (11),
c         (12) (stem diameter) = x * (sapwood diameter)
c         From (5), (9) & (12),
c         (13) leafmass = lai * allom1 * x**reinickerp * 
c                       (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
c         From (13),
c         (14) leafmass = [ lai * allom1 * x**reinickerp *
c                (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(1-1/reinickerp)

          x=pftpar(pft,22)

          lm_sapl(pft,1)=(lai_sapl*allom1*x**reinickerp*
     *      (4.0*sla(pft)/pi/latosa)**(reinickerp*0.5) / sla(pft))**
     *      (1.0-1.0/reinickerp)  !eqn 14
          lm_sapl(pft,2)=17.8-co2(2)   ! initial 13C value from llyod & farquhar, 1994
          lm_sapl(pft,3)=0.

c         Calculate sapling stem diameter
c         From (9) & (12),
c         (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

          stemdiam=x*(4.0*lm_sapl(pft,1)*sla(pft)/pi/latosa)**0.5  !Eqn 15

c         Calculate sapling height
c         (16) height = allom2 * (stem diameter)**allom3 (source?)

          height_sapl=allom2*stemdiam**allom3   !Eqn 16

c         Calculate sapling sapwood mass
c         (17) (sapwood volume) = height * (sapwood xs area)
c         (18) (sapwood xs area) = leafmass * sla / latosa
c         From (17) & (18),

c     (19) (sapwood volume) = height * leafmass * sla / latosa
c         (20) (sapwood mass) = (wood density) * (sapwood volume)
c         From (19) & (20),
c         (21) (sapwood mass) = (wood density) * height * leafmass * sla /
c                latosa

          sm_sapl(pft,1)=wooddens*height_sapl*lm_sapl(pft,1)*sla(pft)/
     *      latosa   !Eqn 21
          sm_sapl(pft,2)=lm_sapl(pft,2)   ! 13C value in permille
          sm_sapl(pft,3)=lm_sapl(pft,3)

c         Calculate sapling heartwood mass
c         From (11),
c         (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(pft,1)=(x-1.0)*sm_sapl(pft,1)  !Eqn 22
          hm_sapl(pft,2)=sm_sapl(pft,2)   ! 13C value in permille
          hm_sapl(pft,3)=sm_sapl(pft,3)

        else ! grass PFT

          lm_sapl(pft,1)=lai_sapl/sla(pft)

c         Set initial 13C values for saplings, grass

          if (pftpar(pft,2).eq.1) then   !C4 plants
            lm_sapl(pft,2)=3.6-co2(2)  !from lloyd & farquhar,1994           
            lm_sapl(pft,3)=0.
          else                           !C3 plpants
            lm_sapl(pft,2)=17.8-co2(2)  !from lloyd & farquhar,1994          
            lm_sapl(pft,3)=0.
          endif

          sm_sapl(pft,2)=0.0             ! no sapwood and hartwood
          hm_sapl(pft,2)=0.0             ! for grass PFT
          sm_sapl(pft,3)=0.0             ! no sapwood and hartwood
          hm_sapl(pft,3)=0.0             ! for grass PFT

        endif

c       Calculate sapling or initial grass rootmass
c       (23) lmtorm = (leafmass) / (rootmass)

        lmtorm=pftpar(pft,18) 
        rm_sapl(pft,1)=(1.0/lmtorm)*lm_sapl(pft,1)  !From Eqn 23
        rm_sapl(pft,2)=lm_sapl(pft,2)       ! 13C value in permille

      !write(0,*)pft,1000.*stemdiam,height_sapl,lai_sapl,lm_sapl(pft,1)

      enddo ! pft loop

      !write(0,*)'done pftinitassign'
      !read(*,*)

      return
      end
