!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
module kernel_interface_mod
use constants_mod, only: dp
use argument_mod,  only: argument_type
implicit none
private

!-------------------------------------------------------------------------------
! Public parameters
!-------------------------------------------------------------------------------

integer, parameter, public :: CELLS = 0

!-------------------------------------------------------------------------------
! Public types.
!-------------------------------------------------------------------------------

type :: kernel_interface_type
  private
  integer                          :: iterates_over
  type(argument_type), allocatable :: meta_arguments(:)
contains
  procedure :: validate
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(kernel_type) function constructor(iterates_over, meta_arguments) result(self)
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  implicit none

  !Arguments
  integer, intent(in) :: iterates_over
  type(argument_type) :: meta_arguments(:)

  self%iterates_over  = iterates_over

  !Guard against incomplete F2003 compiler support with an explicit allocate.
  allocate(self%meta_arguments(1:size(meta_arguments))
  self%meta_arguments = meta_arguments

  return
end function constructor

subroutine validate(self)
  !-----------------------------------------------------------------------------
  ! Placeholder for a potential validator routine. e.g. check that kernel
  ! arguments specified by a specific kernel developer are compatible with each
  ! other.
  !-----------------------------------------------------------------------------
  implicit none

  class(kernel_interface_type) :: self

  return
end subroutine validate

end module kernel_interface_mod
