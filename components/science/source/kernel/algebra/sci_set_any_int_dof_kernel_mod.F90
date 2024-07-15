
!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets given local DoFs of an integer field to the input value. For
!!        example, set the W DoFs for an integer W2 field to a particular value.
!>
module sci_set_any_int_dof_kernel_mod

  use argument_mod,      only: arg_type,            &
                               GH_FIELD, GH_SCALAR, &
                               GH_REAL, GH_INTEGER, &
                               GH_WRITE, GH_READ,   &
                               ANY_SPACE_1, CELL_COLUMN
  use constants_mod,     only: i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: set_any_int_dof_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                          &
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_SPACE_1), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),               &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_any_int_dof_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: set_any_int_dof_code

contains

!> @details For the input local DoF number, set the integer field at these DoFs
!!          to the integer value provided. For example, set the W DoFs for an
!!          integer W2 field to a particular value.
!! @param[in] nlayers       The number of layers
!! @param[in,out] field     A field on any functions space
!! @param[in] dof_to_update The DoF of the field to update
!! @param[in] value         The value to assign to the given DoF
!! @param[in] ndf           The number of degrees of freedom per cell for any space
!! @param[in] undf          The number of unique degrees of freedom for any space
!! @param[in] map           Array holding the dofmap for the cell at the
!!                          base of the column for any space
subroutine set_any_int_dof_code(nlayers,                &
                                field,                  &
                                dof_to_update,          &
                                value,                  &
                                ndf, undf, map          &
                                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf, undf

  integer(kind=i_def), dimension(undf), intent(inout) :: field
  integer(kind=i_def),                  intent(in)    :: dof_to_update
  integer(kind=i_def),                  intent(in)    :: value
  integer(kind=i_def), dimension(ndf),  intent(in)    :: map

  ! Internal variables
  integer(kind=i_def) :: k

  do k = 0, nlayers-1

     field(map(dof_to_update) + k) = value

  end do

end subroutine set_any_int_dof_code

end module sci_set_any_int_dof_kernel_mod
