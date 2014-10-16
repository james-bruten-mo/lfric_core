!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Kernel which operates on w3 field. Determines the RHS of Galerkin projection

!> @detail The kernel computes the integral of rho_df * P 
!! with P_analytic over a single column


module w3_rhs_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use gaussian_quadrature_mod, only : ngp_h, ngp_v, gaussian_quadrature_type
use argument_mod,            only: arg_type, &          ! the type
                                   gh_rw, gh_read, w3, w0, fe, cells ! the enums                                       

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: w3_rhs_kernel_type
  private
  type(arg_type) :: meta_args(4) = [ &
       arg_type(gh_rw,   w3,fe,.true., .false.,.false.,.true.),  &
       arg_type(gh_read, w0,fe,.false.,.true., .false.,.false.), &
       arg_type(gh_read, w0,fe,.false.,.false.,.false.,.false.), &
       arg_type(gh_read, w0,fe,.false.,.false.,.false.,.false.)  &
       ]
  integer :: iterates_over = cells

contains
  procedure, nopass :: rhs_w3_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface w3_rhs_kernel_type
   module procedure w3_rhs_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rhs_w3_code
contains

type(w3_rhs_kernel_type) function w3_rhs_kernel_constructor() result(self)
  return
end function w3_rhs_kernel_constructor

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w3 The number of degrees of freedom per cell
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column
!! @param[in] w3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[inout] X Real array, the actual data
!! @param[inout] gq Type, gaussian quadrature rule
!! @param[in] ndf_w0 The number of degrees of freedom per cell
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column
!! @param[in] w0_diff_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[inout] chi_1 Real array, the x component of the w0 coordinate field
!! @param[inout] chi_2 Real array, the y component of the w0 coordinate field
!! @param[inout] chi_3 Real array, the z component of the w0 coordinate field

subroutine rhs_w3_code(nlayers, &
                       ndf_w3, map_w3, w3_basis, x, gq, &
                       ndf_w0, map_w0, w0_diff_basis, chi_1, chi_2, chi_3 &
                       )
                       
  use coordinate_jacobian_mod, only: coordinate_jacobian                       
                         
  ! needs to compute the integral of rho_df * P 
  ! P_analytic over a single column

  !Arguments
  integer, intent(in) :: nlayers, ndf_w3, ndf_w0
  integer, intent(in) :: map_w3(ndf_w3), map_w0(ndf_w0)
  real(kind=r_def), intent(in), dimension(1,ndf_w3,ngp_h,ngp_v) :: w3_basis 
  real(kind=r_def), intent(inout) :: x(*)
  real(kind=r_def), intent(in)    :: chi_1(*), chi_2(*), chi_3(*)
  type(gaussian_quadrature_type), intent(inout) :: gq
  
  real(kind=r_def), intent(in), dimension(3,ndf_w0,ngp_h,ngp_v) :: w0_diff_basis 

  !Internal variables
  integer               :: df, k
  integer               :: qp1, qp2
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac
  real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), pointer                    :: wgp_h(:), wgp_v(:)

  wgp_h => gq%get_wgp_h()
  wgp_v => gq%get_wgp_v()  
  ! compute the analytic R integrated over one cell
  do k = 0, nlayers-1
    do df = 1, ndf_w0
      chi_1_e(df) = chi_1( map_w0(df) + k)
      chi_2_e(df) = chi_2( map_w0(df) + k)
      chi_3_e(df) = chi_3( map_w0(df) + k)
    end do
    call coordinate_jacobian(ndf_w0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e, w0_diff_basis, jac, dj)
    do df = 1, ndf_w3
       x(map_w3(df) + k) = 0.0_r_def
       do qp2 = 1, ngp_v
          do qp1 = 1, ngp_h
             integrand = w3_basis(1,df,qp1,qp2) * real(k+1) * dj(qp1,qp2)
             x(map_w3(df) + k) = x(map_w3(df) + k) &
                               + 0.125*wgp_h(qp1)*wgp_v(qp2)*integrand
          end do
       end do
    end do
  end do
  
end subroutine rhs_w3_code

end module w3_rhs_kernel_mod
