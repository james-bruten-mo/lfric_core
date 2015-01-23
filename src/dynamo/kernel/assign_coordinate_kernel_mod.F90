!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the coordinate field from the grid node values

module assign_coordinate_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def, earth_radius
use argument_mod,            only : arg_type, &          ! the type
                                    GH_READ, GH_WRITE, ANY_SPACE, FE, CELLS ! the enums       
use mesh_mod,                only : dz, l_spherical  

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: assign_coordinate_kernel_type
  private
  type(arg_type) :: meta_args(3) = [  &
       arg_type(GH_WRITE,ANY_SPACE,FE,.false.,.false.,.false.,.false.),        &
       arg_type(GH_WRITE,ANY_SPACE,FE,.false.,.false.,.false.,.false.),        &
       arg_type(GH_WRITE,ANY_SPACE,FE,.false.,.false.,.false.,.false.)         &
       ]
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: assign_coordinate_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface assign_coordinate_kernel_type
   module procedure assign_coordinate_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public assign_coordinate_code
contains

type(assign_coordinate_kernel_type) function assign_coordinate_kernel_constructor() result(self)
  return
end function assign_coordinate_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[out] chi_1 Real array the first component of the genralised physical coordinates 
!! @param[out] chi_2 Real array the second component of the genralised physical coordinates 
!! @param[out] chi_3 Real array the third component of the genralised physical coordinates 
!! @param[in] vertex_coords Real array. the cordinates of the vertices
!! @param[in] x_node Real array. the cordinates of the nodal points in the function space
subroutine assign_coordinate_code(nlayers,ndf,nverts,map,chi_1,chi_2,chi_3, &
                                  vertex_coords,chi_hat_node,chi_hat_vert)

  !Arguments
  integer, intent(in) :: nlayers, ndf, nverts
  integer, intent(in) :: map(ndf)  
  real(kind=r_def), intent(out) :: chi_1(*), chi_2(*), chi_3(*)
  real(kind=r_def), intent(in)  :: vertex_coords(3,nverts,nlayers)
  real(kind=r_def), intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)

  !Internal variables
  integer          :: k, df, dfk, vert
  
  real(kind=r_def) :: interp_weight, x, y, z, radius_correction

  radius_correction = 1.0_r_def
  
  ! compute the representation of the coordinate field
  do k = 0, nlayers-1
    do df = 1, ndf 
! compute interpolation weights
      x = 0.0_r_def
      y = 0.0_r_def
      z = 0.0_r_def
      do vert = 1,nverts
        interp_weight = (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                       *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                       *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))
      
        x = x + interp_weight*vertex_coords(1,vert,k+1)
        y = y + interp_weight*vertex_coords(2,vert,k+1)
        z = z + interp_weight*vertex_coords(3,vert,k+1)
      end do
! For spherical domains we need to project x,y,z back onto spherical shells
      if ( l_spherical ) then
        radius_correction = earth_radius + (real(k) + chi_hat_node(3,df))*dz
        radius_correction = radius_correction/sqrt(x*x + y*y + z*z)
      end if  
      dfk = map(df)+k 
      chi_1(dfk) = x*radius_correction
      chi_2(dfk) = y*radius_correction
      chi_3(dfk) = z*radius_correction
    end do
  end do
  
end subroutine assign_coordinate_code

end module assign_coordinate_kernel_mod
