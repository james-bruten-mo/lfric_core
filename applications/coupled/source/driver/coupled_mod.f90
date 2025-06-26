!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Definition of what configuration the coupled app requires.
!>
module coupled_mod

  implicit none

  private

  character(*), public, parameter ::                        &
      coupled_required_namelists(5) =   [ 'base_mesh     ', &
                                          'extrusion     ', &
                                          'finite_element', &
                                          'partitioning  ', &
                                          'planet        ']

end module coupled_mod
