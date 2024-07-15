!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the surface area of a cell projected to mean sea level
!>        i.e. ignoring the orographic effect on the area
module sci_calc_da_msl_proj_kernel_mod

use argument_mod,          only: arg_type,                               &
                                 GH_FIELD, GH_SCALAR, GH_READ, GH_WRITE, &
                                 CELL_COLUMN, ANY_DISCONTINUOUS_SPACE_1, &
                                 GH_REAL, GH_INTEGER
use fs_continuity_mod,     only: W2
use constants_mod,         only: r_def, i_def
use kernel_mod,            only: kernel_type
use reference_element_mod, only: B

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_da_msl_proj_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                       &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,   W2),                         &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),  &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                               &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                               &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                            &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                             &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_da_msl_proj_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: calc_da_msl_proj_code

contains

!> @details Using the input field dA_at_w2, we take the area of the model top
!>          level (which is always flat) and scale this by
!>          (radius_surface/radius_top)**2 to obtain the area of a cell at
!>          mean-sea-level in the absence of orography.
!> @param[in]     nlayers     Number of layers
!> @param[in,out] dA_msl_proj dA projected to mean sea level
!> @param[in]     dA_at_w2    dA at the w2 points of the mesh
!> @param[in]     radius      radius at mean sea level
!> @param[in]     domain top  height of domain top above radius
!> @param[in]     ndf_2d      Number of degrees of freedom per cell for 2d
!> @param[in]     undf_2d     Number of unique degrees of freedom for 2d
!> @param[in]     map_2d      Dofmap for the cell at the base of the column for 2d
!> @param[in]     ndf_w2  Number of degrees of freedom per cell for w2
!> @param[in]     undf_w2 Number of unique degrees of freedom for w2
!> @param[in]     map_w2  Dofmap for the cell at the base of the column for w2
subroutine calc_da_msl_proj_code(nlayers, dA_at_w2, dA_msl_proj, &
                                 radius, domain_top,             &
                                 geometry, geometry_spherical,   &
                                 ndf_w2, undf_w2, map_w2,        &
                                 ndf_2d, undf_2d, map_2d)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_2d, undf_w2, undf_2d

  integer(kind=i_def), dimension(ndf_2d), intent(in) :: map_2d
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_2d), intent(inout)  :: dA_msl_proj
  real(kind=r_def), dimension(undf_w2), intent(in)     :: dA_at_w2

  real(kind=r_def)    :: radius, domain_top
  integer(kind=i_def) :: geometry, geometry_spherical

  if (geometry == geometry_spherical) then
    dA_msl_proj(map_2d(1)) = dA_at_w2(map_w2(B) + nlayers) *           &
                             ( radius / (radius+domain_top) )**2_i_def
  else
    dA_msl_proj(map_2d(1)) = dA_at_w2(map_w2(B) + nlayers)
  end if

end subroutine calc_da_msl_proj_code

end module sci_calc_da_msl_proj_kernel_mod
