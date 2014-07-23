!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
program generate_cubedsphere
use gencube_mod,    only: gencube_type 
use ugrid_file_mod, only: ugrid_file_type
use ugrid_2d_mod,   only: ugrid_2d_type
use ncdf_quad_mod,  only: ncdf_quad_type
implicit none

type(gencube_type)  :: cubed_sphere_generator 
type(ugrid_2d_type) :: ugrid_2d
!type(ugrid_2d_type) :: ugrid_2d_recover

class(ugrid_file_type), allocatable :: ugrid_file_handler

character(len=*), parameter :: filename = 'ugrid_quads_2d.nc'

integer :: istat
integer :: mesh_res
character(len=10) :: mesh_res_string

!If argument is present on the command line, use that value. Otherwise go with
!default of 2.
call get_command_argument(1, value=mesh_res_string, status=istat)
if(istat == 0) then
  read(mesh_res_string,*) mesh_res
  write(*,'(a,2x,i0)') 'Argument supplied:', mesh_res
else
  mesh_res = 2
  write(*,'(a,2x,i0)') 'No argument found. Defaulting to: ', mesh_res
end if

!Set up the file handler
allocate(ncdf_quad_type :: ugrid_file_handler)
call ugrid_2d%set_file_handler(ugrid_file_handler)

!Generate the grid
cubed_sphere_generator = gencube_type(half_num_cells_on_panel_edge=mesh_res)
call ugrid_2d%set_by_generator(cubed_sphere_generator)

!Write data to file
call ugrid_2d%write_to_file(trim(filename))

!! Demonstrate re-read.
!!Set up the file handler
!allocate(ncdf_quad_type :: ugrid_file_handler)
!call ugrid_2d_recover%set_file_handler(ugrid_file_handler)
!
!call ugrid_2d_recover%read_from_file(trim(filename))
!
!call ugrid_2d%write_coordinates()

stop
end program generate_cubedsphere


