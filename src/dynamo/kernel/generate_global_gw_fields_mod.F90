!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module generate_global_gw_fields_mod
!> @brief module that contains routines taken to initialise fields based upon
!> DCMIP test 31 - non-hydrostatic gravity waves
!> @detail The non-hydrostatic gravity wave test examines the response of models to short time-scale wavemotion
!> triggered by a localized perturbation. The formulation presented in this document is new,
!> but is based on previous approaches by Skamarock et al. (JAS 1994), Tomita and Satoh (FDR 2004), and
!> Jablonowski et al. (NCAR Tech Report 2008) 

use constants_mod, only: r_def, PI, GRAVITY, Cp, P_ZERO, &
                         omega, earth_radius, N_SQ, KAPPA, Rd

implicit none

contains

subroutine generate_global_gw_fields (lat, z, exner, u, theta, rho)

implicit none
        
  real(kind=r_def), intent(in)  :: lat, z ! Latitude (radians) and Height (m)
                                   
  real(kind=r_def), intent(out) :: u(3), &               ! (Zonal,Meridional,Vertical) wind (m s^-1)
                                   theta, &              ! potential Temperature (K)
                                   exner, &              ! exner pressure
                                   rho                   ! density (kg m^-3)

  real(kind=r_def), parameter :: U0        = 0.0_r_def,    &     ! Reference Velocity 
                                 T_EQUATOR = 300.0_r_def,   &     ! Temperature at Equator    
                                 ZTOP      = 10000.0_r_def        ! Model Top       
                           
  real(kind=r_def) :: bigG = (GRAVITY*GRAVITY)/(N_SQ*Cp)      ! G constant from DCMIP formulation                            
  real(kind=r_def) :: tsurf, psurf                            ! Surface temperature (k) and pressure (Pa)
  real(kind=r_def) :: temperature, pressure                   ! temperature(k) and pressure (Pa)
  real(kind=r_def) :: exp_fac
  real(kind=r_def) :: p_equator = P_ZERO

! intialise wind field
  u(1) = U0 * cos(lat)
  u(2) = 0.0_r_def
  u(3) = 0.0_r_def

! 
  exp_fac = (U0+2.0_r_def*omega*earth_radius)*(cos(2.0_r_def*lat)-1.0_r_def)

! Compute surface temperture
  tsurf = bigG + (T_EQUATOR - bigG)*exp( -(U0*N_SQ/(4.0_r_def*GRAVITY*GRAVITY))*exp_fac ) 

! Compute surface pressure
  psurf = p_equator*exp( (U0/(4.0_r_def*bigG*Rd))*exp_fac  ) * (tsurf/T_EQUATOR)**(Cp/Rd)

! Compute pressure and temperature
  pressure = psurf*( (bigG/tsurf)*exp(-N_SQ*z/GRAVITY)+1.0_r_def - (bigG/tsurf)  )**(Cp/Rd)

  temperature = bigG*(1.0_r_def - exp(N_SQ*z/GRAVITY))+ tsurf*exp(N_SQ*z/GRAVITY)

! Compute density from equation of state
  rho = pressure/(Rd*temperature)

! convert pressure to exner pressure and temperature to potential temperature
  exner = (pressure/P_ZERO)**KAPPA
  theta = temperature/exner

end subroutine generate_global_gw_fields

!=================================================================================

pure function generate_global_gw_pert(lon, lat, z) result(theta)
!> @brief Function to generate the potential temperature pertubation for 
!> the global gravity wave test
implicit none

  real(kind=r_def)              :: theta
  real(kind=r_def), intent(in)  :: lon, lat, z

  real(kind=r_def) :: sin_tmp, cos_tmp, r, shape_function

  real(kind=r_def), parameter :: LAMBDAC = 2.0_r_def*PI/3.0_r_def,     &     ! Lon of Pert Center
                                 D       = 5000.0_r_def,               &     ! Width for Pert
                                 PHIC    = 0.0_r_def,                  &     ! Lat of Pert Center
                                 DELTA_THETA = 1.0_r_def,              &     ! Max Amplitude of Pert
                                 LZ      = 10000.0_r_def                     ! Vertical half-Wavelength of Pert
 
  sin_tmp = sin(lat) * sin(PHIC)
  cos_tmp = cos(lat) * cos(PHIC)

! great circle distance  
  r  = earth_radius * acos (sin_tmp + cos_tmp*cos(lon-LAMBDAC)) 

  shape_function = (D**2)/(D**2 + r**2)

  theta = DELTA_THETA*shape_function*sin(PI*z/LZ)
end function generate_global_gw_pert

end module generate_global_gw_fields_mod

