!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page transport Transport miniapp
!> Program file for running transport miniapp. Subroutine calls include initialise_transport(),
!> run_transport() and finalise_transport().
program transport

  use transport_driver_mod, only: initialise_transport, &
                                  run_transport,        &
                                  finalise_transport

  implicit none

  call initialise_transport()

  call run_transport()

  call finalise_transport()

end program transport
