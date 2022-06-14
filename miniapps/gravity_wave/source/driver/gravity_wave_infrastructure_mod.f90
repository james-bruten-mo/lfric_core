!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls infrastructure related information used by the model

module gravity_wave_infrastructure_mod

  use check_configuration_mod,    only : get_required_stencil_depth
  use cli_mod,                    only : get_initial_filename
  use clock_mod,                  only : clock_type
  use configuration_mod,          only : final_configuration
  use constants_mod,              only : i_def, i_native, PRECISION_REAL, r_def
  use convert_to_upper_mod,       only : convert_to_upper
  use derived_config_mod,         only : set_derived_config
  use gravity_wave_mod,           only : load_configuration
  use log_mod,                    only : log_event,          &
                                         log_set_level,      &
                                         log_scratch_space,  &
                                         initialise_logging, &
                                         finalise_logging,   &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_WARNING,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_TRACE
  use mesh_mod,                   only : mesh_type
  use mpi_mod,                    only : get_comm_size, get_comm_rank
  use io_context_mod,             only : io_context_type
  use field_mod,                  only : field_type
  use mesh_collection_mod,        only : mesh_collection, &
                                         mesh_collection_type
  use local_mesh_collection_mod,  only : local_mesh_collection, &
                                         local_mesh_collection_type
  use driver_comm_mod,            only : init_comm, final_comm
  use driver_fem_mod,             only : init_fem
  use driver_io_mod,              only : init_io, final_io, get_clock
  use driver_mesh_mod,            only : init_mesh
  use runtime_constants_mod,      only : create_runtime_constants
  use formulation_config_mod,     only : l_multigrid

  implicit none

  private
  public initialise_infrastructure, finalise_infrastructure

  character(*), public, parameter   :: xios_ctx = 'gravity_wave'

contains

  !> @brief Initialises the infrastructure used by the model
  !> @param [in]     program_name  An identifier given to the model begin run
  !> @param [in,out] mesh          The model prime mesh
  !> @param [in,out] twod_mesh     The model prime 2D mesh
  subroutine initialise_infrastructure(program_name, &
                                       mesh,         &
                                       twod_mesh)

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING
    implicit none

    character(*),      intent(in) :: program_name
    type(mesh_type),   intent(inout), pointer :: mesh
    type(mesh_type),   intent(inout), pointer :: twod_mesh

    type(field_type), target :: chi(3)
    type(field_type), target :: panel_id

    integer(i_def),   allocatable :: multigrid_mesh_ids(:)
    integer(i_def),   allocatable :: multigrid_2d_mesh_ids(:)
    type(field_type), allocatable :: chi_mg(:,:)
    type(field_type), allocatable :: panel_id_mg(:)

    integer(i_def)    :: stencil_depth
    integer(i_native) :: log_level, comm

    class(clock_type), pointer :: clock
    real(r_def)                :: dt_model

    character(:), allocatable :: filename

    ! Set up the communicator for later use
    call init_comm(program_name, comm)

    call get_initial_filename( filename )
    call load_configuration( filename )

    call initialise_logging(get_comm_rank(), get_comm_size(), program_name)

    select case (run_log_level)
    case( RUN_LOG_LEVEL_ERROR )
      log_level = LOG_LEVEL_ERROR
    case( RUN_LOG_LEVEL_WARNING )
      log_level = LOG_LEVEL_WARNING
    case( RUN_LOG_LEVEL_INFO )
      log_level = LOG_LEVEL_INFO
    case( RUN_LOG_LEVEL_DEBUG )
      log_level = LOG_LEVEL_DEBUG
    case( RUN_LOG_LEVEL_TRACE )
      log_level = LOG_LEVEL_TRACE
    end select

    call log_set_level( log_level )

    write(log_scratch_space,'(A)')                              &
        'Runtime message logging severity set to log level: '// &
        convert_to_upper(key_from_run_log_level(run_log_level))
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .false. )


    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------
    allocate( mesh_collection, &
              source=mesh_collection_type() )
    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )

    stencil_depth = get_required_stencil_depth()

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), stencil_depth, mesh,  &
                    twod_mesh             = twod_mesh,             &
                    multigrid_mesh_ids    = multigrid_mesh_ids,    &
                    multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                    use_multigrid         = l_multigrid )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id,                           &
                   multigrid_mesh_ids    = multigrid_mesh_ids,    &
                   multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                   chi_mg                = chi_mg,                &
                   panel_id_mg           = panel_id_mg,           &
                   use_multigrid         = l_multigrid )

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------
    call init_io( xios_ctx, comm, chi, panel_id )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------
    clock => get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field and limited area masks.
    call create_runtime_constants( mesh, twod_mesh,                        &
                                   chi, panel_id, dt_model,                &
                                   mg_mesh_ids    = multigrid_mesh_ids,    &
                                   mg_2D_mesh_ids = multigrid_2D_mesh_ids, &
                                   chi_mg         = chi_mg,                &
                                   panel_id_mg    = panel_id_mg )


  end subroutine initialise_infrastructure


  !> @brief Finalises infrastructure used by the model
  subroutine finalise_infrastructure()

    implicit none

    ! Finalise I/O
    call final_io()

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise communicator
    call final_comm()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_infrastructure

end module gravity_wave_infrastructure_mod
