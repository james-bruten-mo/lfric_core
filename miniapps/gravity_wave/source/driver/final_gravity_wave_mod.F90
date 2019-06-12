!-----------------------------------------------------------------------------
! Copyright (c) 2019,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Finalise the gravity wave simulation

!> @details Generates checksum output and finalises the algorithm
!>          (deallocates objects set up for the solver api)

module final_gravity_wave_mod

  use field_mod,                      only : field_type
  use field_collection_mod,           only : field_collection_type
  use checksum_alg_mod,               only : checksum_alg
  use gravity_wave_alg_mod,           only : gravity_wave_alg_final

  implicit none

  contains

  !> @brief Finalise the gravity wave simulation
  !> @param [inout] prognostics A collection of all the prognostic fields.
  !> @param [in] program_name An identifier given to the model begin run
  subroutine final_gravity_wave( prognostic_fields, program_name )

    implicit none

    type(field_collection_type ), intent(inout) :: prognostic_fields
    character(*), intent(in) :: program_name

    type( field_type), pointer :: wind => null()
    type( field_type), pointer :: pressure => null()
    type( field_type), pointer :: buoyancy => null()

    wind => prognostic_fields%get_field('wind')
    pressure => prognostic_fields%get_field('pressure')
    buoyancy => prognostic_fields%get_field('buoyancy')

    ! Write checksums to file
    call checksum_alg( program_name, wind, 'wind', buoyancy, 'buoyancy', pressure, 'pressure')

    call gravity_wave_alg_final()

  end subroutine final_gravity_wave

end module final_gravity_wave_mod
