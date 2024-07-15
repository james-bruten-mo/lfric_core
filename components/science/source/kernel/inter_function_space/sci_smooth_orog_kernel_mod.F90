!-----------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------
module sci_smooth_orog_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_WRITE,         &
                                ANY_DISCONTINUOUS_SPACE_1, &
                                CELL_COLUMN, STENCIL, REGION
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: smooth_orog_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: smooth_orog_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: smooth_orog_code

contains

subroutine smooth_orog_code( nlayers,                      &
                             orog_out, orog_in,            &
                             stencil_w3_size, stencil_map_w3,          &
                             ndf_w3, undf_w3, map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: stencil_w3_size

  integer(kind=i_def), dimension(ndf_w3,stencil_w3_size), intent(in)  :: stencil_map_w3
  integer(kind=i_def), dimension(ndf_w3),                 intent(in)  :: map_w3

  real(kind=r_def), dimension(undf_w3),  intent(inout) :: orog_out
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: orog_in

  ! If the full stencil isn't available, we must be at the domain edge.
  if ( stencil_w3_size < 9_i_def ) then

    orog_out( map_w3(1) ) = orog_in( map_w3(1) )

    return
  end if

  ! Gaussian blur (1-2-1) filter

  orog_out( map_w3(1) ) = ( 4.0_r_def * orog_in( stencil_map_w3(1,1) ) + &
                            2.0_r_def * orog_in( stencil_map_w3(1,2) ) + &
                            2.0_r_def * orog_in( stencil_map_w3(1,4) ) + &
                            2.0_r_def * orog_in( stencil_map_w3(1,6) ) + &
                            2.0_r_def * orog_in( stencil_map_w3(1,8) ) + &
                            1.0_r_def * orog_in( stencil_map_w3(1,3) ) + &
                            1.0_r_def * orog_in( stencil_map_w3(1,5) ) + &
                            1.0_r_def * orog_in( stencil_map_w3(1,7) ) + &
                            1.0_r_def * orog_in( stencil_map_w3(1,9) )  &
                          ) / 16.0_r_def

end subroutine smooth_orog_code

end module sci_smooth_orog_kernel_mod
