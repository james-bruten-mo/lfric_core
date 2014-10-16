!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides access to the members of the w1_kernel class.

!> @details Accessor functions for the w1_kernel class are defined in this module.

!> @param RHS_w1_code              Code to implement the RHS for a w1 field
!> @param gaussian_quadrature      Contains result of gaussian quadrature

module w1_kernel_mod
use gaussian_quadrature_mod, only: gaussian_quadrature_type, &
                                   ngp_h, ngp_v ! parameter for how many GQ points
use argument_mod,            only: arg_type, &          ! the type
                                   gh_inc, w1, fe, cells ! the enums
use kernel_mod,              only: kernel_type
use constants_mod,           only: r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: w1_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/ &
       arg_type(gh_inc,w1,fe, .true., .false., .false., .true.) &
       /)
  integer :: iterates_over = cells

contains
  procedure, nopass :: rhs_w1_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface w1_kernel_type
   module procedure w1_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rhs_w1_code              
contains

type(w1_kernel_type) function w1_kernel_constructor() result(self)
  return
end function w1_kernel_constructor
  
subroutine rhs_w1_code(nlayers,ndf,undf,map,w1_basis,x,gq)

!> @brief This subroutine calculates the RHS of Galerkin projection on W1 space.
!! @param[in] nlayers Integer: The number of layers.
!! @param[in] ndf Integer: The number of degrees of freedom per cell.
!! @param[in] map Integer: Array holding the dofmap for the cell at the base of the column.
!! @param[in] w3_basis Real: 4-dim array holding VECTOR basis functions evaluated at quadrature points.
!! @param[inout] X Real: The array of actual data.
!! @param[inout] gq Type: Quadrature rule (here Gaussian).

  ! Needs to compute the integral of c_df * P 
  ! P_analytic over a single column
  
  !Arguments
  integer, intent(in) :: nlayers, ndf,undf
  integer, intent(in) :: map(ndf)
  real(kind=r_def), intent(in), dimension(3,ndf,ngp_h,ngp_v) :: w1_basis 
  real(kind=r_def), intent(inout) :: x(undf)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer                                  :: df, k
  integer                                  :: qp1, qp2
  real(kind=r_def), dimension(1,3)         :: basisfun_i
  real(kind=r_def), dimension(3,1)         :: constantvec
  real(kind=r_def), dimension(1,1)         :: T_1
  real(kind=r_def), pointer                :: wgp_h(:), wgp_v(:)

  wgp_h => gq%get_wgp_h()
  wgp_v => gq%get_wgp_v()

  constantvec(1,1) =  1.0_r_def;
  constantvec(2,1) =  2.0_r_def;  
  constantvec(3,1) =  4.0_r_def;  
  ! Compute the analytic R integrated over one cell
  do k = 0, nlayers-1
    do df = 1, ndf
       do qp1 = 1, ngp_h
          do qp2 = 1, ngp_v
             basisfun_i(1,1:3) = w1_basis(1:3,df,qp1,qp2)
             T_1 = MATMUL(basisfun_i,constantvec)
             x(map(df) + k) = x(map(df) + k) &
                            + 0.125_r_def*wgp_h(qp1)*wgp_v(qp2)*T_1(1,1)
          end do
       end do
    end do
  end do
end subroutine rhs_w1_code

end module w1_kernel_mod
