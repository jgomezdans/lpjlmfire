module initsoilmod

implicit none

public :: initsoil

contains

!-----------------------------------------------------------------------------------------------

subroutine initsoil()

use netcdf
use typesizes

use errormod,       only : netcdf_err,ncstat
use iovariablesmod, only : srtx,cntx,srty,cnty,inputlonlen,inputlatlen,lonvect,latvect,bounds,soilfile,soilfid,soil
use parametersmod,  only : sp

implicit none


real(sp), allocatable, dimension(:) :: depth
real(sp), allocatable, dimension(:,:,:) :: sand
real(sp), allocatable, dimension(:,:,:) :: clay

integer :: i,j
integer :: varid
integer :: dimid

integer :: layers

integer, dimension(1) :: pos
integer, dimension(2) :: xpos,ypos


!-------------------------
!generate file names. Regardless of the path name the input files must always have these names.

write(0,'(a,a)')'using soil file: ',trim(soilfile)

!-------------------------
!open soil file

ncstat = nf90_open(soilfile,nf90_nowrite,soilfid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_inq_dimid(soilfid,'lon',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(soilfid,dimid,len=inputlonlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_dimid(soilfid,'lat',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(soilfid,dimid,len=inputlatlen)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

allocate(lonvect(inputlonlen))
allocate(latvect(inputlatlen))

ncstat = nf90_inq_varid(soilfid,'lon',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(soilfid,varid,lonvect)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)
  
ncstat = nf90_inq_varid(soilfid,'lat',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)  

ncstat = nf90_get_var(soilfid,varid,latvect)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!-------------------------
!calculate the number and indices of the pixels to be calculated

!gridres(1) = 60. * (lonvect(2) - lonvect(1))
!gridres(2) = 60. * (maxval(latvect) - minval(latvect)) / inputlatlen
!lboundlon = minval(lonvect)
!uboundlat = maxval(latvect)

pos = minloc(abs(lonvect - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lonvect - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(latvect - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(latvect - bounds(4)))
ypos(2) = pos(1)

srtx = minval(xpos)
srty = minval(ypos)

if (bounds(1) == bounds(2) .and. bounds(3) == bounds(4)) then  !special case, run just one nearest gridcell even if coords were ambiguous
  
  cntx = 1
  cnty = 1
  
else

  if (lonvect(srtx) < bounds(1)) srtx = srtx + 1
  cntx = 1 + abs(maxval(xpos) - srtx)

  if (latvect(srty) < bounds(3)) srty = srty + 1
  cnty = 1 + abs(maxval(ypos) - srty)

end if

!--------------------------
!get soil properties from file

write(0,*)'reading soil data'

ncstat = nf90_inq_dimid(soilfid,'layer',dimid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inquire_dimension(soilfid,dimid,len=layers)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)  
  
allocate(depth(layers))

allocate(sand(cntx,cnty,layers))
allocate(clay(cntx,cnty,layers))
!allocate(sorg(cntx,cnty,layers))
!allocate(whc (cntx,cnty,layers))

allocate(soil(cntx,cnty))

ncstat = nf90_inq_varid(soilfid,'zpos',varid)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_get_var(soilfid,varid,depth)
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(soilfid,'sand',varid)
ncstat = nf90_get_var(soilfid,varid,sand,start=[srtx,srty,1],count=[cntx,cnty,layers])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

ncstat = nf90_inq_varid(soilfid,'clay',varid)
ncstat = nf90_get_var(soilfid,varid,clay,start=[srtx,srty,1],count=[cntx,cnty,layers])
if (ncstat /= nf90_noerr) call netcdf_err(ncstat)

!------------

do j=1,cnty
  do i=1,cntx

!    if (whc(i,j,1) < 0.) then
!
!      soil(i,j)%whc = whc(i,j,:)
!      soil(i,j)%cond = -1.
!
   !else
      !----------------------------------------------------------------------------------------
      !LGM soil file has only values for 2 layers, so redistribute accordingly on 5 layers 
      
      if(layers == 2) then 
        
        soil(i,j)%sand(3:5) = sand(i,j,2)
        soil(i,j)%sand(1:2) = sand(i,j,1)
        soil(i,j)%clay(3:5) = clay(i,j,2)
        soil(i,j)%clay(1:2) = clay(i,j,1) 
             
      else if(layers == 5) then
        
        soil(i,j)%sand = sand(i,j,:)
        soil(i,j)%clay = clay(i,j,:)
        
      end if  
       
      !----------------------------------------------------------------------------------------
      

      soil(i,j)%whc = 0 ! this value is calculated later in soilmod.f90 !  max(whc(i,j,:),0.)          !whc in cm m-1
      
      soil(i,j)%cond(1) = fKsat(soil(i,j)%sand(1),0.10) !cond in mm s-1
      soil(i,j)%cond(2) = fKsat(soil(i,j)%sand(2),0.30)
      soil(i,j)%cond(3) = fKsat(soil(i,j)%sand(3),0.50)
      soil(i,j)%cond(4) = fKsat(soil(i,j)%sand(4),0.70)
      soil(i,j)%cond(5) = fKsat(soil(i,j)%sand(5),0.90)

!    end if

  end do
end do

do i = 1,5
  soil%sand(i) = soil%sand(i) * 0.01  !convert to fraction
  soil%clay(i) = soil%clay(i) * 0.01  !convert to fraction
  soil%orgm(i) = 0. !    sorg(:,:,i) * 1.724 !convert organic carbon to organic matter
  
  if(layers == 2) then
    
    soil%zpos(3) = depth(2)
    soil%zpos(4) = depth(2)
    soil%zpos(5) = depth(2)
    soil%zpos(1) = depth(1)
    soil%zpos(2) = depth(1)
    
  else if(layers == 5) then 
     
    soil%zpos(i) = depth(i)
    
  end if
end do

!-------------------------

deallocate(sand)
deallocate(clay)
!deallocate(sorg)
!deallocate(whc)

write(0,*) 'Done reading soil data'

!-------------------------

contains

!----------------------------------------------------------------------
function fKsat(sand,depth)

!function to calculate saturated conductivity (mm s-1)
!related to pore geometry and particle shape also packing and structure
!was from  Vereecken et al. 1990 (removed in favor of CLM scheme)

implicit none

real(sp), parameter :: a =  0.0070556
real(sp), parameter :: b = -0.8840000
real(sp), parameter :: c =  0.0153000
real(sp), parameter :: zs = 0.5000000  !length scale for depth decrease (m)

real(sp) :: fKsat

real(sp), intent(in) :: sand   !soil sand concent (percent)
real(sp), intent(in) :: depth  !soil depth (m)

!---

!function from CLM3.0
 
fKsat = a * (10.d0**(b + c * sand)) * exp(-depth / zs)

 
end function fKsat

!----------------------------------------------------------------------

end subroutine initsoil

end module initsoilmod
