!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets a given dof of a field to the input value.
!>        This only works for the lowest-order W1 fields
!>
module sci_set_w1h_dofs_kernel_mod

  use argument_mod,      only: arg_type,            &
                               GH_FIELD, GH_SCALAR, &
                               GH_REAL, GH_INTEGER, &
                               GH_WRITE, GH_READ,   &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type
  use fs_continuity_mod, only: W1

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: set_w1h_dofs_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                        &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W1),           &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                 &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_w1h_dofs_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: set_w1h_dofs_code

contains

!> @details For the input dof number, set the field at this dof to the real
!!          value provided.
!! @param[in] nlayers       The number of layers
!! @param[in,out] field     A field on any functions space
!! @param[in] value         The value to assign to the given dof
!! @param[in] ndf           The number of degrees of freedom per cell for any space
!! @param[in] undf          The number of unique degrees of freedom for any space
!! @param[in] map           Array holding the dofmap for the cell at the
!!                          base of the column for any space
subroutine set_w1h_dofs_code(nlayers,                &
                             field,                  &
                             value,                  &
                             ndf, undf, map          &
                             )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf, undf

  real(kind=r_def), dimension(undf), intent(inout) :: field
  real(kind=r_def), intent(in)                     :: value
  integer(kind=i_def), dimension(ndf),  intent(in) :: map

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 0, nlayers-1
    do df =1,4
      field(map(df) + k ) = value
      field(map(df+8) + k ) = value
    end do
  end do

end subroutine set_w1h_dofs_code

end module sci_set_w1h_dofs_kernel_mod
