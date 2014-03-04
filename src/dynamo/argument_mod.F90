!-----------------------------------------------------------------------------
! (c) The copyright relating to this information/data is jointly owned by
! the Crown, Met Office and NERC 2013.
! The contribution of STFC in creating this information/data is acknowledged.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! DESCRIPTION
!   An arg has a function space (where dofs live on cell) what the stencil is,
!   this is not Grad Phi, but which facets are touched.  could be simple, for
!   example the FE cell integral stencil.  intent: what happens to the data
!   members
!-----------------------------------------------------------------------------
module argument_mod
use iso_c_binding
use lfric
implicit none
private

enum, bind(c) 
  !The following value is valid for any arg.
  enumerator :: GH_READ

  !The following values are only valid for fields.
  enumerator :: GH_WRITE, GH_READWRITE, GH_INC

  !The following values are only valid for globals.
  enumerator :: GH_MIN, GH_MAX, GH_SUM
end enum

  !args(fs,stencil,arg_intent) ! this need defining
type, public :: argument_type
  public
  integer(kind(GH_READ)) :: arg_intent
  integer                :: element
  integer(kind(FE))      :: stencil=0
end type 

!-------------------------------------------------------------------------------
! Expose public type parameters
!-------------------------------------------------------------------------------

! Types to enable declarations of elements.
integer, public, parameter :: CG(3)   = [1,2,3]
integer, public, parameter :: DG(0:3) = [0,1,2,3]
integer, public, parameter :: R=0

public :: GH_READ, GH_WRITE, GH_READWRITE, GH_INC
public :: GH_SUM, GH_MIN, GH_MAX 

end module argument_mod
