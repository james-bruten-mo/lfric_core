!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief The unweighted restriction operation for multidata scalar fields.
!> @details Restrict a scalar multidata field on a fine grid to a coarse grid multidata field. No
!!          weighting is used for the averaging.
!!          This kernel only works for the lowest-order W3 and Wtheta spaces
!!          and for column first multidata fields.

module sci_restrict_multidata_scalar_kernel_mod

use constants_mod,           only: i_def, r_double, r_single
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   GH_SCALAR, GH_INTEGER,     &
                                   GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: restrict_multidata_scalar_kernel_type
   private
   type(arg_type) :: meta_args(2) = (/                                   &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                              mesh_arg=GH_COARSE),       &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, &
                                              mesh_arg=GH_FINE  )        &
        /)
  integer :: operates_on = CELL_COLUMN
end type restrict_multidata_scalar_kernel_type

public :: restrict_multidata_scalar_kernel_code

  ! Generic interface for real32 and real64 types
  interface restrict_multidata_scalar_kernel_code
    module procedure  &
      restrict_multidata_scalar_code_r_single, &
      restrict_multidata_scalar_code_r_double
  end interface

contains

  !> @brief The unweighted restriction operation for scalar multidata fields.
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid field to write to
  !> @param[in]     fine_field               The fine grid field to restrict
  !> @param[in]     ndata                    Number fields in multidata field
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid

  ! R_SINGLE PRECISION
  ! ==================
  subroutine restrict_multidata_scalar_code_r_single(               &
                                           nlayers,                 &
                                           cell_map,                &
                                           ncell_fine_per_coarse_x, &
                                           ncell_fine_per_coarse_y, &
                                           ncell_fine,              &
                                           coarse_field,            &
                                           fine_field,              &
                                           ndata,                   &
                                           undf_coarse,             &
                                           map_coarse,              &
                                           ndf,                     &
                                           undf_fine,               &
                                           map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_single), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_single), intent(in)    :: fine_field(undf_fine)
    integer(kind=i_def), intent(in)    :: ndata

    integer(kind=i_def) :: k, x_idx, y_idx, k_start, i, top_df
    real(kind=r_single) :: denom
    integer(kind=i_def), parameter :: df = 1 ! Lowest order function space

    denom = 1.0_r_single/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, kind=r_single)

    ! Assume lowest order W3 or Wtheta space
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do i = 1, ndata
          k_start = (i-1)*(top_df+1)+1
          do k = 0, top_df
            coarse_field(map_coarse(df) +k+k_start-1) =  coarse_field(map_coarse(df) +k+k_start-1) + &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+k+k_start-1)*denom
          end do
        end do
      end do
    end do

  end subroutine restrict_multidata_scalar_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  subroutine restrict_multidata_scalar_code_r_double(               &
                                           nlayers,                 &
                                           cell_map,                &
                                           ncell_fine_per_coarse_x, &
                                           ncell_fine_per_coarse_y, &
                                           ncell_fine,              &
                                           coarse_field,            &
                                           fine_field,              &
                                           ndata,                   &
                                           undf_coarse,             &
                                           map_coarse,              &
                                           ndf,                     &
                                           undf_fine,               &
                                           map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_double), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_double), intent(in)    :: fine_field(undf_fine)
    integer(kind=i_def), intent(in)    :: ndata

    integer(kind=i_def) :: k, x_idx, y_idx, k_start, i, top_df
    real(kind=r_double) :: denom
    integer(kind=i_def), parameter :: df = 1 ! Lowest order function space

    denom = 1.0_r_double/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, kind=r_double)

    ! Assume lowest order W3 or Wtheta space
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do i = 1, ndata
          k_start = (i-1)*(top_df+1)+1
          do k = 0, top_df
            coarse_field(map_coarse(df) +k+k_start-1) =  coarse_field(map_coarse(df) +k+k_start-1) + &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+k+k_start-1)*denom
          end do
        end do
      end do
    end do

  end subroutine restrict_multidata_scalar_code_r_double


end module sci_restrict_multidata_scalar_kernel_mod
