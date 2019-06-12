!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Steps the gravity wave miniapp through one timestep

!> @details Handles the stepping (for a single timestep) of the
!>          gravity wave app

module step_gravity_wave_mod

  use field_mod,                      only : field_type
  use field_collection_mod,           only : field_collection_type
  use gravity_wave_alg_mod,           only : gravity_wave_alg_step

  implicit none

  private
  public step_gravity_wave

  contains

  !> @brief Steps the gravity wave miniapp through one timestep
  !> @param [inout] prognostic_fields A collection of all the prognostic fields
  subroutine step_gravity_wave(prognostic_fields)

    implicit none
    type( field_collection_type ), intent(inout) :: prognostic_fields

    type( field_type), pointer :: wind => null()
    type( field_type), pointer :: pressure => null()
    type( field_type), pointer :: buoyancy => null()

    wind => prognostic_fields%get_field('wind')
    pressure => prognostic_fields%get_field('pressure')
    buoyancy => prognostic_fields%get_field('buoyancy')

    call gravity_wave_alg_step(wind, pressure, buoyancy)

  end subroutine step_gravity_wave

end module step_gravity_wave_mod
