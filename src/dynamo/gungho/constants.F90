!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines various constants.

!> @details Various physical and geometrical constants are defined in this module.
!! Their values are also set here.
module constants_mod
implicit none

!Working precision
integer,       parameter :: dp=8                    !< working precision

!Numerical constants
real(kind=dp), parameter :: pi=4.0_dp*atan(1.0_dp)  !< pi value
real(kind=dp), parameter :: eps=3.0E-15_dp          !< relative precision

! Physical constants
real(kind=dp) :: earth_radius = 6371229.0_dp

integer,       parameter :: max_iter = 10 ! maximum iteration number for solver

end module constants_mod

