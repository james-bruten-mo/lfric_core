!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief The argument type to hold kernel metadata required by the psy layer.

!> @details Metadata for the kernels. For each field passed to a kernel 
!> the psy layer needs to know how this field is to be accessed.
!> read, write etc, to which function space it belongs and what
!> stencil it operates over. These are the three integers of
!> the arg_type and the values are then one of the parameters
!> defined in this module. 
!> field metadata also has three logicals controlling whether the psy layer
!> needs to pass the basis function, the differential basis function,
!> and the guassian quadrature type.
!> Another metadatum which describes the kernel, not the fields
!> is what the kernel will iterate over. Usually cells, sometimes
!> all the dofs.

module argument_mod

  use function_space_mod, only : W0, W1, W2, W3
  implicit none

! access descriptors
  integer, public, parameter :: GH_READ  = 1 
  integer, public, parameter :: GH_WRITE = 2
  integer, public, parameter :: GH_RW    = 3
  integer, public, parameter :: GH_INC   = 4
  integer, public, parameter :: GH_SUM   = 5
  integer, public, parameter :: GH_MIN   = 6
  integer, public, parameter :: GH_MAX   = 7

  integer, public, parameter :: ANY_SPACE = 0    

! stencil label
  integer, public, parameter :: FE = 1 

! kernel iterator
  integer, public, parameter :: CELLS     = 1
  integer, public, parameter :: ALL_DOFS  = 2

  type, public :: arg_type
     integer :: arg_intent
     integer :: vspace
     integer :: stencil
     logical :: basis
     logical :: diff_basis
     logical :: nodal_coords
     logical :: gaussian_quad
  end type arg_type

end module argument_mod

