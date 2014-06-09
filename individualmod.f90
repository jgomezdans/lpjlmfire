module individualmod

use parametersmod, only : sp

implicit none

public :: allomind

type sizeind
  real(sp) :: massf
  real(sp) :: height
  real(sp) :: lcrown
  real(sp) :: stemd
  real(sp) :: barkt
end type sizeind

!-------------------------------------------------------------------

contains

subroutine allomind(pft,height,ind) !hclass,cheight,clength,stemd,barkt)

use parametersmod, only : sp,npft,nhclass

implicit none

!arguments

integer,  intent(in) :: pft
real(sp), intent(in) :: height

type(sizeind), dimension(:), intent(out) :: ind


!parameters

real(sp), parameter :: ot = 1. / 3.      !one third

real(sp), parameter :: hsapl = 1.

real(sp), parameter, dimension(npft) :: CLf    = [     ot, 0.1,        ot,     ot,     ot,     ot,     ot, -99., -99.  ]  !crown length
real(sp), parameter, dimension(npft) :: par1   = [ 0.0301, 0.1085, 0.0670, 0.0451, 0.0347, 0.0292, 0.0347, -99., -99.  ]  !bark thickness 1
real(sp), parameter, dimension(npft) :: par2   = [ 0.0281, 0.2120, 0.5590, 0.1412, 0.1086, 0.1086, 0.1086, -99., -99.  ]  !bark thickness 2

!slope and intercept are the result of a comparison between LPJ-modeled maximum PFT-height and empirically observed maximum PFT-height

real(sp), parameter, dimension(npft) :: hs = [  1.75,  1.90,  2.57,  1.23,  1.68,  1.23,  1.32, -99., -99. ]  !height slope (hmax_r-hsap_r)/(hmax_m-hsap_m)
real(sp), parameter, dimension(npft) :: hi = [ -3.75, -4.20, -6.20, -2.20, -3.53, -2.20, -2.46, -99., -99. ]  !height intercept (hsap_r - hs * hsap_m)

real(sp), parameter, dimension(npft) :: ds = [  2.52,  3.40,  4.60,  6.86,  3.40,  7.95,  3.35, -99., -99. ]  !diameter slope (cm diameter gained per m height)
real(sp), parameter, dimension(npft) :: di = [ -0.78, -2.11, -3.90, -7.30, -2.11, -8.92, -2.03, -99., -99. ]  !diameter intercept

real(sp), parameter, dimension(nhclass) :: hsd    = [ -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0 ]

!local variables

integer :: i

real(sp) :: hbar
real(sp) :: hmax
real(sp) :: hsl

!------
!initialize height class biomass fraction (based on skewed normal distribution)

ind%massf = [ 0.0185, 0.0557, 0.1312, 0.2419, 0.3339, 0.1974, 0.0154 ]

!calculate mean of the height class distribution based on the average individual height

hbar = hs(pft) * height + hi(pft)

!maximum individual height can be 25% taller than the current mean height of the height class distribution

hmax = 1.25 * hbar  

!calculate difference

hsl = hmax - hbar !well, that just equals 0.25 * hbar

!calculate height of each height class

ind%height = hsd * hsl + hbar

ind%lcrown = ind%height * CLf(pft)

!calculate fraction of total biomass in each height class
!check if one of the lower height classes is shorter than a sapling, if so, roll the distribution up to taller height classes


do i = 1,6 ! we have 7 height classes, so go up to 6 to roll up to class 7 max and not step over the array boundary (M.P. 12.07.2013)
  if (ind(i)%height < hsapl) then
    
    ind(i+1)%massf = ind(i+1)%massf + ind(i)%massf
    
   ! write(0,*) i, ind(i)%height, hsapl, ind(i)%massf, ind(i+1)%massf
    
    ind(i)%massf = 0.
    
  end if
end do


!calculate stem diameter

ind%stemd = max(1.,ds(pft) * ind%height + di(pft)) !constrained to be not less than 1 cm

!calculate bark thickness

ind%barkt = par1(pft) * ind%stemd + par2(pft)                              !bark thickness (cm) eqn. 21

end subroutine allomind

end module individualmod
