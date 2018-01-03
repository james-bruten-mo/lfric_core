!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the absolute vorticity

!> @details Compute the vorticity ( =curl(u) ) using a finite difference
!>          approximation on a uniform grid
module fd_vorticity_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    ANY_SPACE_9, W1, W2,                     &
                                    GH_DIFF_BASIS, CELLS,                    &
                                    GH_QUADRATURE_XYoZ

use constants_mod,           only : r_def
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: fd_vorticity_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W1),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: fd_vorticity_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface fd_vorticity_kernel_type
   module procedure fd_vorticity_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public fd_vorticity_code
contains

type(fd_vorticity_kernel_type) function fd_vorticity_kernel_constructor() result(self)
  return
end function fd_vorticity_kernel_constructor

!> @brief Compute the projection of curl(u) into the vorticity function space
!! @param[in] nlayers Number of layers
!! @param[in] ndf_xi Number of degrees of freedom per cell for w1
!! @param[in] undf_xi Unique number of degrees of freedom for w1
!! @param[in] map_xi Dofmap for the cell at the base of the column for w1
!! @param[in] diff_basis_xi Differential of the basis functions evaluated at gaussian quadrature point
!! @param[inout] rhs Right hand side to be computed
!! @param[in] ndf_u Number of degrees of freedom per cell for the velocity field
!! @param[in] undf_u Unique number of degrees of freedom for the velocity field
!! @param[in] map_u Dofmap for the cell at the base of the column for the velocity field
!! @param[in] basis_u Basis functions evaluated at gaussian quadrature points
!! @param[in] u Velocity field
!! @param[in] ndf_chi Number of degrees of freedom per cell for the function space containing chi
!! @param[in] undf_chi Unique number of degrees of freedom for the chi arrays
!! @param[in] map_chi Dofmap for the cell at the base of the column for the function space containing chi
!! @param[in] diff_basis_chi Differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of the horizontal quadrature points
!! @param[in] wqp_v Weights of the vertical quadrature points
subroutine fd_vorticity_code(nlayers,                           &
                         xi, u, chi_1, chi_2, chi_3,            &
                         ndf_xi, undf_xi, map_xi,               &
                         ndf_u, undf_u, map_u,                  &
                         ndf_chi, undf_chi, map_chi,            &
                         diff_basis_chi,                        &    
                         nqp_h, nqp_v, wqp_h, wqp_v             &                       
                         )
                           
  use coordinate_jacobian_mod,  only: coordinate_jacobian
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_chi, ndf_u, ndf_xi, undf_chi, undf_u, undf_xi 
  integer, intent(in) :: nqp_h, nqp_v 
  integer, dimension(ndf_xi),  intent(in) :: map_xi
  integer, dimension(ndf_u),   intent(in) :: map_u
  integer, dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), dimension(undf_xi),               intent(inout) :: xi
  real(kind=r_def), dimension(undf_u),                intent(in)    :: u 
  real(kind=r_def), dimension(undf_chi),              intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  !Internal variables
  integer               :: df, k, loc
  
  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)        :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)    :: jac
  real(kind=r_def), dimension(ndf_u)  :: u_cell(ndf_u) 
  real(kind=r_def), dimension(ndf_xi) :: xi_cell(ndf_xi)
  real(kind=r_def)                    :: dx, dy, dz
  real(kind=r_def), parameter         :: m1 = -1.0_r_def
  real(kind=r_def), parameter         ::  h =  0.5_r_def

  ! m1 parameter (-1) is used to cancel out the fact that the basis vectors for 
  ! y components of the velocity (dofs 2 & 4) point in the negative y direction
  ! If m1 is ignored then the computation is the same as would be obtained in a 
  ! standard C-grid discretisation

  xi_cell(:) = 0.0_r_def

  ! Extract element arrays of chi
  k = 0
  do df = 1, ndf_chi
    loc = map_chi(df) + k
    chi_1_e(df) = chi_1( loc )
    chi_2_e(df) = chi_2( loc )
    chi_3_e(df) = chi_3( loc )
  end do
  call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                           diff_basis_chi, jac, dj)   
  dx = jac(1,1,1,1)
  dy = jac(2,2,1,1)
  dz = jac(3,3,1,1)
  ! compute the vorticity
  ! xi(1) & xi(3) & xi(9) & xi(11)  =  du/dz - dw/dx
  ! xi(2) & xi(4) & xi(10) & xi(12) = -dv/dz + dw/dy
  ! xi(5) & xi(6) & xi(7) & xi(8)   = dv/dx - du/dy
  k = 0
  do df = 1, ndf_u
    u_cell(df) = u( map_u(df) + k )
  end do    
  ! x - components
  xi_cell(2)  =   u_cell(5)*dx/(dx*dy*dy)
  xi_cell(4)  = - u_cell(5)*dx/(dx*dy*dy)
  xi_cell(10) =   u_cell(6)*dx/(dx*dy*dy) + m1*u_cell(2)*dx/(dx*dz*dz)
  xi_cell(12) = - u_cell(6)*dx/(dx*dy*dy) + m1*u_cell(4)*dx/(dx*dz*dz)

  ! y - components
  xi_cell(1)  = - u_cell(5)*dy/(dx*dy*dx)
  xi_cell(3)  =   u_cell(5)*dy/(dx*dy*dx)
  xi_cell(9)  = - u_cell(6)*dy/(dx*dy*dx) - u_cell(1)*dy/(dy*dz*dz)
  xi_cell(11) =   u_cell(6)*dy/(dx*dy*dx) - u_cell(3)*dy/(dy*dz*dz)
 
  ! z - components
  xi_cell(5) =   m1*u_cell(2)*dz/(dx*dz*dx) - u_cell(1)*dz/(dy*dz*dy)
  xi_cell(6) = - m1*u_cell(2)*dz/(dx*dz*dx) - u_cell(3)*dz/(dy*dz*dy)
  xi_cell(7) = - m1*u_cell(4)*dz/(dx*dz*dx) + u_cell(3)*dz/(dy*dz*dy)
  xi_cell(8) =   m1*u_cell(4)*dz/(dx*dz*dx) + u_cell(1)*dz/(dy*dz*dy)

  do df = 1,ndf_xi
    xi(map_xi(df) + k) = xi(map_xi(df) + k) + h*xi_cell(df)
  end do


  do k = 1,nlayers-2
    do df = 1, ndf_u
      u_cell(df) = u( map_u(df) + k )
    end do  

    ! x - components
    xi_cell(2)  =   u_cell(5)*dx/(dx*dy*dy) - m1*u_cell(2)*dy/(dy*dz*dz)
    xi_cell(4)  = - u_cell(5)*dx/(dx*dy*dy) - m1*u_cell(4)*dy/(dy*dz*dz)
    xi_cell(10) =   u_cell(6)*dx/(dx*dy*dy) + m1*u_cell(2)*dx/(dx*dz*dz)
    xi_cell(12) = - u_cell(6)*dx/(dx*dy*dy) + m1*u_cell(4)*dx/(dx*dz*dz)

    ! y - components
    xi_cell(1)  = - u_cell(5)*dy/(dx*dy*dx) + u_cell(1)*dy/(dy*dz*dz)
    xi_cell(3)  =   u_cell(5)*dy/(dx*dy*dx) + u_cell(3)*dy/(dy*dz*dz)
    xi_cell(9)  = - u_cell(6)*dy/(dx*dy*dx) - u_cell(1)*dy/(dy*dz*dz)
    xi_cell(11) =   u_cell(6)*dy/(dx*dy*dx) - u_cell(3)*dy/(dy*dz*dz)
 
    ! z - components
    xi_cell(5) =   m1*u_cell(2)*dz/(dx*dz*dx) - u_cell(1)*dz/(dy*dz*dy)
    xi_cell(6) = - m1*u_cell(2)*dz/(dx*dz*dx) - u_cell(3)*dz/(dy*dz*dy)
    xi_cell(7) = - m1*u_cell(4)*dz/(dx*dz*dx) + u_cell(3)*dz/(dy*dz*dy)
    xi_cell(8) =   m1*u_cell(4)*dz/(dx*dz*dx) + u_cell(1)*dz/(dy*dz*dy)

    do df = 1,ndf_xi
      xi(map_xi(df) + k) = xi(map_xi(df) + k) + h*xi_cell(df)
    end do
  end do

  k = nlayers-1
  do df = 1, ndf_u
    u_cell(df) = u( map_u(df) + k )
  end do  

  ! x - components
  xi_cell(2)  =   u_cell(5)*dx/(dx*dy*dy) - m1*u_cell(2)*dy/(dy*dz*dz)
  xi_cell(4)  = - u_cell(5)*dx/(dx*dy*dy) - m1*u_cell(4)*dy/(dy*dz*dz)
  xi_cell(10) =   u_cell(6)*dx/(dx*dy*dy)
  xi_cell(12) = - u_cell(6)*dx/(dx*dy*dy)

  ! y - components
  xi_cell(1)  = - u_cell(5)*dy/(dx*dy*dx) + u_cell(1)*dy/(dy*dz*dz)
  xi_cell(3)  =   u_cell(5)*dy/(dx*dy*dx) + u_cell(3)*dy/(dy*dz*dz)
  xi_cell(9)  = - u_cell(6)*dy/(dx*dy*dx)
  xi_cell(11) =   u_cell(6)*dy/(dx*dy*dx)
 
  ! z - components
  xi_cell(5) =   m1*u_cell(2)*dz/(dx*dz*dx) - u_cell(1)*dz/(dy*dz*dy)
  xi_cell(6) = - m1*u_cell(2)*dz/(dx*dz*dx) - u_cell(3)*dz/(dy*dz*dy)
  xi_cell(7) = - m1*u_cell(4)*dz/(dx*dz*dx) + u_cell(3)*dz/(dy*dz*dy)
  xi_cell(8) =   m1*u_cell(4)*dz/(dx*dz*dx) + u_cell(1)*dz/(dy*dz*dy)

  do df = 1,ndf_xi
    xi(map_xi(df) + k) = xi(map_xi(df) + k) + h*xi_cell(df)
  end do
end subroutine fd_vorticity_code

end module fd_vorticity_kernel_mod
