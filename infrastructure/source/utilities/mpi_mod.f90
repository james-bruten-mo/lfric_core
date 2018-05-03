!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief Provides access to the MPI related functionality
!>
!> Provides access to global reduction functions and generation of the halo
!> exchange redistribution object. In order for that functionality to work, the
!> subroutine store_comm must first be called to store a valid MPI communicator.
!>
module mpi_mod
  use constants_mod, only : r_def, i_def, i_native, l_def, i_halo_index
  use mpi
  use yaxt,          only: xt_redist, xt_idxlist, xt_xmap, &
                           xt_idxvec_new, xt_xmap_dist_dir_new, &
                           xt_redist_p2p_off_new, &
                           xt_xmap_delete, xt_idxlist_delete
  use log_mod,       only: log_event, LOG_LEVEL_ERROR
  implicit none

  private
  public initialise_comm, finalise_comm, store_comm,&
         global_sum, global_min, global_max, &
         generate_redistribution_map, &
         get_mpi_datatype

  ! The mpi communicator
  integer(i_def), private :: comm
  ! Flag marks whether an MPI communicator has been stored
  logical(l_def), private :: comm_set = .false.

contains

  !> Initialises MPI and returns mpi_comm_world as the communicator 
  !> as well as storing it in a private variable, ready for later use.
  !> Note: Only to be used when XIOS isn't being used, Normally XIOS
  !>       initialises MPI and returns a communicator.
  !>
  !> @param out_comm The MPI communicator.
  !>
  subroutine initialise_comm(out_comm)
    implicit none
    integer(i_def), intent(out) :: out_comm
    integer(i_def) :: ierr

    call mpi_init(ierr)
    if (ierr /= mpi_success) &
          call log_event('Unable to initialise MPI', LOG_LEVEL_ERROR )
    out_comm = mpi_comm_world
    comm = out_comm
    comm_set = .true.
  end subroutine initialise_comm

  !> Finalises MPI 
  !> Note: Only to be used when XIOS isn't being used.
  !>
  subroutine finalise_comm()
    implicit none
    integer(i_def) :: ierr

    call mpi_finalize(ierr)
    if (ierr /= mpi_success) &
          call log_event('Unable to finalise MPI', LOG_LEVEL_ERROR )
  end subroutine finalise_comm

  !> Stores the MPI communicator in a private variable, ready for later use.
  !>
  !> @param in_comm The MPI communicator to be stored.
  !>
  subroutine store_comm(in_comm)
    implicit none
    integer(i_def), intent(in) :: in_comm

    comm = in_comm
    comm_set = .true.
  end subroutine store_comm


  !> Calculates the global sum of a collection of local arrays
  !>
  !> @param l_array The local arrays that are to be summed
  !> @param g_sum The calculated global sum 
  !> @param count The number of elements in the local array
  !> @param err The return code
  !>
  subroutine global_sum(l_array, g_sum, count, err)
    implicit none
    integer(i_def), intent(in) :: count
    integer(i_def), intent(out) :: err
    real(r_def), intent(in) :: l_array(count)
    real(r_def), intent(out) :: g_sum

    integer(i_def) :: i
    real(r_def) :: l_sum

    if(comm_set)then
      ! Generate local sum
      l_sum = 0.0
      do i = 1, count
        l_sum = l_sum + l_array(i)
      end do

      ! Generate global sum
      call mpi_allreduce( l_sum, g_sum, 1, get_mpi_datatype(r_def), &
                          mpi_sum, comm, err )
    else
      call log_event( &
      'Call to global_sum failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_sum


  !> Calculates the global minimum of a collection of local arrays
  !>
  !> @param l_array The local arrays from which the minimum is to be found
  !> @param g_min The calculated global minimum
  !> @param count The number of elements in the local array
  !> @param err The return code
  !>
  subroutine global_min(l_array, g_min, count, err)
    implicit none
    integer(i_def), intent(in) :: count
    integer(i_def), intent(out) :: err
    real(r_def), intent(in) :: l_array(count)
    real(r_def), intent(out) :: g_min

    integer(i_def) :: i
    real(r_def) :: l_min

    if(comm_set)then
      ! Generate local min
      l_min = l_array(1)
      do i = 2, count
        if( l_array(i) < l_min ) l_min = l_array(i)
      end do

      ! Generate global min
      call mpi_allreduce( l_min, g_min, 1, get_mpi_datatype(r_def), &
                          mpi_min, comm, err )
    else
      call log_event( &
      'Call to global_min failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_min


  !> Calculates the global maximum of a collection of local arrays
  !>
  !> @param l_array The local arrays from which the maximum is to be found
  !> @param g_max The calculated global maximum
  !> @param count The number of elements in the local array
  !> @param err The return code
  !>
  subroutine global_max(l_array, g_max, count, err)
    implicit none
    integer(i_def), intent(in) :: count
    integer(i_def), intent(out) :: err
    real(r_def), intent(in) :: l_array(count)
    real(r_def), intent(out) :: g_max

    integer(i_def) :: i
    real(r_def) :: l_max

    if(comm_set)then
      ! Generate local max
      l_max = l_array(1)
      do i = 2, count
        if( l_array(i) > l_max ) l_max = l_array(i)
      end do

      ! Generate global max
      call mpi_allreduce( l_max, g_max, 1, get_mpi_datatype(r_def), &
                          mpi_max, comm, err )
    else
      call log_event( &
      'Call to global_max failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_max


  !> Generate the halo exchange redistribution object fo be used for future
  !> halo exchanges
  !>
  !> @param src_indices The global indices of all the owned points in this
  !>                    MPI task
  !> @param tgt_indices The global indices of all the halo points in this
  !>                    MPI task
  !> @param datatype    The MPI datatype of a single element in the data to be
  !>                    exchanged
  !> @return redist     The halo exchange redistribution object
  !>
  function generate_redistribution_map(src_indices, tgt_indices, datatype) &
                                       result(redist)
    implicit none
    integer(i_halo_index), intent(in) :: src_indices(:), tgt_indices(:)
    integer(i_def),        intent(in) :: datatype
    type(xt_redist) :: redist

    type(xt_idxlist) :: src_idxlist, tgt_idxlist
    type(xt_xmap) :: xmap
    integer(i_def), pointer :: src_offsets(:)
    integer(i_def), pointer :: tgt_offsets(:)
    integer(i_def) :: i

    if(comm_set)then
      ! create decomposition descriptors
      src_idxlist = xt_idxvec_new( src_indices, size(src_indices) )
      tgt_idxlist = xt_idxvec_new( tgt_indices, size(tgt_indices) )

      ! generate exchange map
      xmap = xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, comm)

      allocate(src_offsets( size(src_indices) ))
      allocate(tgt_offsets( size(tgt_indices) ))

      src_offsets = (/(i, i = 0, size(src_indices) - 1)/)
      tgt_offsets = (/(i, i = size(src_indices) , &
                              size(src_indices) + size(tgt_indices) - 1 )/)

      redist = xt_redist_p2p_off_new(xmap, src_offsets,tgt_offsets, datatype)

      call xt_xmap_delete(xmap)
      call xt_idxlist_delete(tgt_idxlist)
      call xt_idxlist_delete(src_idxlist)
      deallocate(src_offsets)
      deallocate(tgt_offsets)
    else
      call log_event( &
      'Call to generate_redistribution_map failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end function generate_redistribution_map


  !> Returns the appropriate MPI datatype enumerator for all the Fortran
  !> kinds supported by the LFRic distributed memory code
  !>
  !> @param f_kind A Fortran kind variable
  !> @return mpi_datatype The MPI datatype enumerator associated with the
  !>                      given Fortran kind
  function get_mpi_datatype(f_kind) result(mpi_datatype)
    use, intrinsic :: iso_fortran_env, only : real64, int32
    implicit none
    integer(i_native), intent(in) :: f_kind
    integer(i_native)             :: mpi_datatype

   ! Determine MPI datatype enumerator from a Fortran kind.
   ! (To support a new Fortran kind, just add a new case clause)
    select case (f_kind)
      case (int32)
        mpi_datatype = MPI_INTEGER
      case (real64)
        mpi_datatype = MPI_DOUBLE_PRECISION
      case default
        call log_event( 'Unrecognised Fortran kind used for MPI comms', &
                         LOG_LEVEL_ERROR )
    end select

  end function get_mpi_datatype

end module mpi_mod
