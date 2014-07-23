!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2012.
! However, it has been created by John Thuburn.
!-------------------------------------------------------------------------------
!>  @brief   Routines for coordinate transformations.
!!
!!  @details Contains routines for conversion of e.g. lat-long to Cartesian XYZ.
!-------------------------------------------------------------------------------
module coord_algorithms_mod
use constants_mod, only: dp
implicit none
private

!Public subroutines
public :: ll2xyz
public :: llr2xyz
public :: xyz2ll
public :: starea2
public :: spdist

!--------------------------------------------------------------------------------
! Contained functions / subroutines
!--------------------------------------------------------------------------------
contains

!--------------------------------------------------------------------------------
!>  @brief  Convert longitude and latitude to Cartesian coordinates on the unit
!!          sphere.
!!
!!  @param[in]   long  The longitude to convert.
!!  @param[in]   lat   The latitude  to convert.
!!  @param[out]  x     The x coordinate.
!!  @param[out]  y     The y coordinate.
!!  @param[out]  z     The z coordinate.
!--------------------------------------------------------------------------------

subroutine ll2xyz(long,lat,x,y,z)
  implicit none
 
  !Arguments
  real(kind=dp), intent(in)  :: long,lat
  real(kind=dp), intent(out) :: x,y,z 
 
  !Internal variables
  real(kind=dp) :: cln, sln, clt, slt
 
  sln=sin(long)
  cln=cos(long)
  slt=sin(lat)
  clt=cos(lat)
 
  x=cln*clt
  y=sln*clt
  z=slt
 
  return
end subroutine ll2xyz

!--------------------------------------------------------------------------------
!>  @brief  Convert longitude and latitude to Cartesian coordinates on 
!!          a sphere with some specified radius.
!!
!!  @param[in]   long    The longitude to convert.
!!  @param[in]   lat     The latitude  to convert.
!!  @param[in]   radius  The radius of the sphere.
!!  @param[out]  x       The x coordinate.
!!  @param[out]  y       The y coordinate.
!!  @param[out]  z       The z coordinate.
!--------------------------------------------------------------------------------

subroutine llr2xyz(long,lat,radius,x,y,z)
  implicit none
 
  !Arguments
  real(kind=dp), intent(in)  :: long,lat,radius
  real(kind=dp), intent(out) :: x,y,z 
 
  !Internal variables
  real(kind=dp) :: cln, sln, clt, slt
 
  sln = sin(long)
  cln = cos(long)
  slt = sin(lat)
  clt = cos(lat)
 
  x = radius * cln * clt
  y = radius * sln * clt
  z = radius * slt
 
  return
end subroutine llr2xyz

!--------------------------------------------------------------------------------
!>  @brief  Convert Cartesian coordinates to longitude and latitude.
!!
!!  @param[in]   x     The x coordinate to convert.
!!  @param[in]   y     The y coordinate to convert.
!!  @param[in]   z     The z coordinate to convert.
!!  @param[out]  long  The longitude.
!!  @param[out]  lat   The latitude.
!--------------------------------------------------------------------------------

subroutine xyz2ll(x,y,z,long,lat)
  use constants_mod, only: pi
  implicit none

  !Arguments
  real(kind=dp), intent(in)  :: x,y,z
  real(kind=dp), intent(out) :: long, lat

  !Internal variables
  real(kind=dp) :: tln, tlt, r

  if (x == 0.0_dp) then
    if (y >= 0.0_dp) then
      long = 0.5_dp*pi
    else
      long = 1.5_dp*pi
    end if
  else
    tln=y/x
    long=atan(tln)
    if (x < 0.0_dp) then
      long=long+pi
    endif
    if (long < 0.0_dp) then
      long=long+2.0_dp*pi
    endif
  end if

  r=sqrt(x*x+y*y)
  if (r == 0.0_dp) then
    if (z > 0.0_dp) then
      lat=0.5_dp*pi
    else
      lat=-0.5_dp*pi
    end if
  else
    tlt=z/r
    lat=atan(tlt)
  end if

  return
end subroutine xyz2ll

!-------------------------------------------------------------------------------
!>  @brief  Calculate the area of a spherical triangle.
!!
!!  @details Calculate the area of the spherical triangle whose corners have
!!           Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2).
!!           The formula below is more robust to roundoff error than the better
!!           known sum of angle - PI formula
!!
!!  @param[in]   x0    Coordinate of a triangle corner.
!!  @param[in]   y0    Coordinate of a triangle corner.
!!  @param[in]   z0    Coordinate of a triangle corner.
!!  @param[in]   x1    Coordinate of a triangle corner.
!!  @param[in]   y1    Coordinate of a triangle corner.
!!  @param[in]   z1    Coordinate of a triangle corner.
!!  @param[in]   x2    Coordinate of a triangle corner.
!!  @param[in]   y2    Coordinate of a triangle corner.
!!  @param[in]   z2    Coordinate of a triangle corner.
!!  @param[out]  area  Area of the spherical triangle.
!-------------------------------------------------------------------------------

subroutine starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area)
  implicit none

  !Arguments
  real(kind=dp), intent(in)  :: x0, y0, z0   
  real(kind=dp), intent(in)  :: x1, y1, z1
  real(kind=dp), intent(in)  :: x2, y2, z2
  real(kind=dp), intent(out) :: area

  !Internal variables
  real(kind=dp) :: d0,d1,d2,s,t0,t1,t2,t3

  !Distances between pairs of points
  call spdist(x0,y0,z0,x1,y1,z1,d2)
  call spdist(x1,y1,z1,x2,y2,z2,d0)
  call spdist(x2,y2,z2,x0,y0,z0,d1)

  !Half perimeter
  s=0.5_dp*(d0+d1+d2)

  !Tangents
  t0 = tan(0.5_dp*(s-d0))
  t1 = tan(0.5_dp*(s-d1))
  t2 = tan(0.5_dp*(s-d2))
  t3 = tan(0.5_dp*s)

  !Area
  area = 4.0_dp*atan(sqrt(t0*t1*t2*t3))

  return
end subroutine starea2

!-------------------------------------------------------------------------------
!> @brief  Calculate the spherical distance between to points.
!!
!! @details  Calculate the spherical distance s between two points with 
!!           Cartesian coordinates (x1,y1,z1), (x2,y2,z2) on the unit sphere.
!!
!! @param[in]  x1  First Cartesian coordinate.
!! @param[in]  x2  First Cartesian coordinate.
!! @param[in]  x3  First Cartesian coordinate.
!! @param[in]  y1  Second Cartesian coordinate.
!! @param[in]  y2  Second Cartesian coordinate.
!! @param[in]  y3  Second Cartesian coordinate.
!! @param[out] s   Spherical distance between the points.
!-------------------------------------------------------------------------------

subroutine spdist(x1,y1,z1,x2,y2,z2,s)
  implicit none

  !Arguments
  real(kind=dp), intent(in)  :: x1, y1, z1, x2, y2, z2
  real(kind=dp), intent(out) :: s

  !Internal variables
  real(kind=dp) :: dx, dy, dz
  real(kind=dp) :: ad

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1

  ad = sqrt(dx*dx + dy*dy + dz*dz)
  s = 2.0_dp*asin(0.5_dp*ad)

  return
end subroutine spdist

end module coord_algorithms_mod

