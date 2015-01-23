!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module calc_exner_pointwise_mod

use constants_mod, only : r_def, KAPPA, Rd, P_ZERO

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Function to compute the exner pressure from the linear equation of state at a
!! point
!! @param[in]  rho      real, the density perturbation
!! @param[in]  theta    real, the potential temperature perturbation
!! @param[in]  exner_s  real, the reference exner pressure
!! @param[in]  rho_s    real, the reference density
!! @param[in]  theta_s  real, the reference potential temperature
!! @param[out] exner    real, the exner pressure perturbation
function calc_exner_pointwise(rho, theta, exner_s, rho_s, theta_s) result(exner)

  real(kind=r_def)              :: exner
  real(kind=r_def), intent(in)  :: rho, theta, exner_s, rho_s, theta_s

! linear
  exner = KAPPA / ( 1.0_r_def - KAPPA ) * exner_s * ( rho/rho_s + theta/theta_s )
  
! nonlinear  
!   exner = ( Rd/P_ZERO * rho * theta ) ** (  ( 1.0_r_def - KAPPA ) / KAPPA )

end function calc_exner_pointwise

end module calc_exner_pointwise_mod
