!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PsyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using 
!> a hand-rolled Psy layer, Psy-lite. A PsyKAl-lite needs a dynamo!
!> Eventually, PsyKAl-lite will be replaced with the real Psy and Dynamo
!> will become Gung Ho.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Creates the function spaces and calls <code>set_up</code> to
!> populate them (either read in or compute) then individual calls to the
!> psy-layer with kernels as if the code has been pre-processed by Psyclone.

program dynamo

  use dynamo_algorithm_mod, only : dynamo_algorithm
  use field_mod,            only : field_data_type
  use lfric
  use log_mod,              only : log_event, LOG_LEVEL_INFO
  use set_up_mod,           only : set_up

  implicit none

  type( function_space_type )      :: v3_function_space, v2_function_space, &
                                      v1_function_space, v0_function_space
  type( field_data_type )          :: pressure_density, rhs
  type( gaussian_quadrature_type ) :: gq
  integer                          :: num_layers

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  call set_up( v0_function_space, v1_function_space, v2_function_space, &
               v3_function_space, num_layers )

  gq = gaussian_quadrature_type( )

  pressure_density = field_data_type( vector_space = v3_function_space, &
                                      gq = gq,                          &
                                      num_layers = num_layers)

  rhs = field_data_type( vector_space = v3_function_space, &
                         gq = gq,                          &
                         num_layers = num_layers )

  call dynamo_algorithm( pressure_density%new_proxy( ), rhs%new_proxy( ) )

  call rhs%print_field( 'RHS field...' )
  call pressure_density%print_field( 'LHS field...' )

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

end program dynamo
