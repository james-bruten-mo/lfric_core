!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp multires_coupling program

!> @brief Main program used for multires_coupling miniapp

!> @details Calls init, run and finalise routines from a driver module

program multires_coupling

  use multires_coupling_driver_mod, only : initialise, run, finalise

  implicit none

  call initialise()

  call run()

  call finalise()

end program multires_coupling
