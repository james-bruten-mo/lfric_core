!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the dot product of a vector field defined by three components in a
!>        scalar space (with itself)
!>
!> \f$ f . f \f$
!>
module sci_dot_product_of_components_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_WRITE, GH_READ, &
                                GH_REAL, CELL_COLUMN,        &
                                ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only : r_def, i_def
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: dot_product_of_components_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                      &
        arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  &
        arg_type(GH_FIELD*3,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1) &
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: dot_product_of_components_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: dot_product_of_components_code

contains

!> @brief The kernel computes the dot product of a vector field defined by scalar components
!>        with itself
!! @param[in] nlayers   Number of layers
!! @param[in,out] dot_prod Dot product
!! @param[in] field1    Field component 1
!! @param[in] field2    Field component 2
!! @param[in] field3    Field component 3
!! @param[in] ndf1    Number of degrees of freedom per cell
!! @param[in] undf1   Number of unique degrees of freedom
!! @param[in] map1    Dofmap for the cell at the base of the column
subroutine dot_product_of_components_code(nlayers,         &
                                          dot_prod,        &
                                          field1,          &
                                          field2,          &
                                          field3,          &
                                          ndf1,            &
                                          undf1,           &
                                          map1)


  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf1
  integer(kind=i_def), intent(in) :: undf1

  integer(kind=i_def), dimension(ndf1), intent(in) :: map1

  real(kind=r_def), dimension(undf1),  intent(inout) :: dot_prod
  real(kind=r_def), dimension(undf1),  intent(in)    :: field1
  real(kind=r_def), dimension(undf1),  intent(in)    :: field2
  real(kind=r_def), dimension(undf1),  intent(in)    :: field3

  ! Internal variables
  integer(kind=i_def) :: df, k

  real(kind=r_def), dimension(3) :: vector

  do k = 0, nlayers-1
    do df = 1, ndf1
      ! Combine components into a vector
      vector(1) = field1(map1(df)+k)
      vector(2) = field2(map1(df)+k)
      vector(3) = field3(map1(df)+k)

      ! Compute dot product
      dot_prod(map1(df)+k) =  dot_product(vector, vector)
    end do
  end do

end subroutine dot_product_of_components_code

end module sci_dot_product_of_components_kernel_mod
