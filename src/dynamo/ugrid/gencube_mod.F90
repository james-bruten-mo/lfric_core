!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2012.
! However, it has been created by John Thuburn.
!-------------------------------------------------------------------------------
!>  @brief    Produce a cubed-sphere grid.
!!
!!  @details  Implements the mesh generator interface to produce a cubed-sphere
!!            grid.
!!
!!  <pre>
!!  Program to generate a cubed sphere grid, including cross-
!!  reference tables of adjacent faces, edges and vertices,
!!  coordinates of faces and vertices, and lengths of edges
!!  and areas of faces.
!!
!!          .....
!!         :  3  |
!!         :     |
!!          -----
!!   -----  .....  .....  -----
!!  |  5  :|  1  :|  2  :|  4  :
!!  |     :|     :|     :|     :
!!   .....  -----  -----  .....
!!          .....
!!         |  6  :
!!         |     :
!!          -----
!!
!!  Solid lines: left and bottom edges of panel
!!  Dotted lines: right and top edges of panel
!!
!!  John Thuburn Nov 2011
!!  </pre>
!-------------------------------------------------------------------------------
module gencube_mod
use constants_mod,       only : r_def, str_def, PI
use ugrid_generator_mod, only : ugrid_generator_type
implicit none
private

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

!Number of grids in multigrid hierarchy
integer, parameter :: NGRIDS = 2

!Number of smoothing iterations. must be at least 1 for consistency of h
!operator.
integer, parameter :: NSMOOTH = 1

real(kind=r_def), parameter :: PIBY4 = PI / 4.0_r_def

!-------------------------------------------------------------------------------
!> @brief    Cubed sphere generator type (function object).
!!
!! @details  Implements abstract mesh generator type.
!-------------------------------------------------------------------------------

! Some representative variable descriptions. (Namings follow a pattern.)

! nfacex : number of faces
! neoff  : number of edges of each face
! neofv  : number of edges associated with each vertex 
! eofv   : edges corresponding to each vertex
! fnxtf  : faces next to each face
! vlong  : Longitude of each vertex

type, extends(ugrid_generator_type), public :: gencube_type
  private

  !n x n cells on each panel smallest and largest n
  integer                       :: n0
  integer                       :: nfacex
  integer                       :: nedgex
  integer                       :: nvertx
  integer         , allocatable :: neoff(:,:)
  integer         , allocatable :: neofv(:,:)
  integer         , allocatable :: nface(:)
  integer         , allocatable :: nedge(:)
  integer         , allocatable :: nvert(:)
  integer         , allocatable :: fnxtf(:,:,:)
  integer         , allocatable :: eoff(:,:,:)
  integer         , allocatable :: voff(:,:,:)
  integer         , allocatable :: fnxte(:,:,:)
  integer         , allocatable :: vofe(:,:,:)
  integer         , allocatable :: fofv(:,:,:)
  integer         , allocatable :: eofv(:,:,:)
  real(kind=r_def), allocatable :: flong(:,:)
  real(kind=r_def), allocatable :: flat(:,:)
  real(kind=r_def), allocatable :: vlong(:,:)
  real(kind=r_def), allocatable :: vlat(:,:)
  real(kind=r_def), allocatable :: farea(:,:)
  real(kind=r_def), allocatable :: ldist(:,:)
  real(kind=r_def), allocatable :: ddist(:,:)

contains
  procedure :: get_dimensions
  procedure :: get_coordinates
  procedure :: get_connectivity
  procedure :: generate
  procedure :: write_data
end type

integer, parameter :: QUADS_NUM_NODES_PER_FACE  = 4
integer, parameter :: QUADS_NUM_EDGES_PER_FACE  = 4
integer, parameter :: QUADS_NUM_NODES_PER_EDGE  = 2
integer, parameter :: QUADS_NUM_EDGE_NEIGHBOURS = 2
integer, parameter :: QUADS_NUM_FACE_NEIGHBOURS = 4 

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

interface gencube_type
  module procedure constructor
end interface gencube_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
!>  @brief   Constructor for cubed-sphere grids.
!!
!!  @param[in]  half_num_cells_on_panel_edge   Half the number of cells along 
!!                                             the edge of one panel.
!-------------------------------------------------------------------------------

type(gencube_type) function constructor(half_num_cells_on_panel_edge) result(self)
  implicit none

  !Arguments
  integer, intent(in) :: half_num_cells_on_panel_edge

  !Internal variables
  integer :: n0, nx

  !n x n cells on each panel smallest and largest n
  n0          = half_num_cells_on_panel_edge
  nx          = n0*(2**(NGRIDS-1))

  self%n0     = n0
  self%nfacex = 6*nx*nx
  self%nedgex = 2*self%nfacex
  self%nvertx = self%nfacex + 2

  !Allocate storage
  allocate (self%neoff (self%nfacex,    NGRIDS))
  allocate (self%neofv (self%nvertx,    NGRIDS))
  allocate (self%nface (                NGRIDS))
  allocate (self%nedge (                NGRIDS))
  allocate (self%nvert (                NGRIDS))
  allocate (self%fnxtf (self%nfacex, 4, NGRIDS))
  allocate (self%eoff  (self%nfacex, 4, NGRIDS))
  allocate (self%voff  (self%nfacex, 4, NGRIDS))
  allocate (self%fnxte (self%nedgex, 2, NGRIDS))
  allocate (self%vofe  (self%nedgex, 2, NGRIDS))
  allocate (self%fofv  (self%nvertx, 4, NGRIDS))
  allocate (self%eofv  (self%nvertx, 4, NGRIDS))
  allocate (self%flong (self%nfacex,    NGRIDS))
  allocate (self%flat  (self%nfacex,    NGRIDS))
  allocate (self%vlong (self%nvertx,    NGRIDS))
  allocate (self%vlat  (self%nvertx,    NGRIDS))
  allocate (self%farea (self%nfacex,    NGRIDS))
  allocate (self%ldist (self%nedgex,    NGRIDS))
  allocate (self%ddist (self%nedgex,    NGRIDS))

  return
end function constructor

!-------------------------------------------------------------------------------
!>  @brief    Gets number of nodes and other mesh dimensions.
!!
!!  @details  Routine provides the number of nodes, edges and faces; also the 
!!            number of nodes per face etc.
!!
!!  @param[in]   self                 The cubed-sphere function object.
!!  @param[out]  num_nodes            The number of nodes on the mesh.
!!  @param[out]  num_edges            The number of edges on the mesh.
!!  @param[out]  num_faces            The number of faces on the mesh.
!!  @param[out]  num_nodes_per_face   The number of nodes around each face.
!!  @param[out]  num_edges_per_face   The number of edges around each face.
!!  @param[out]  num_nodes_per_face   The number of nodes around each edge.
!-------------------------------------------------------------------------------
  
subroutine get_dimensions(self, num_nodes, num_edges, num_faces,    &
                           num_nodes_per_face,  num_edges_per_face, &
                           num_nodes_per_edge)
  implicit none

  !Arguments
  class(gencube_type), intent(in)  :: self

  integer, intent(out) :: num_nodes
  integer, intent(out) :: num_edges
  integer, intent(out) :: num_faces
  integer, intent(out) :: num_nodes_per_face
  integer, intent(out) :: num_edges_per_face
  integer, intent(out) :: num_nodes_per_edge

  num_nodes = self%nvert(NGRIDS)
  num_edges = self%nedge(NGRIDS)
  num_faces = self%nface(NGRIDS)

  num_nodes_per_face = QUADS_NUM_NODES_PER_FACE 
  num_edges_per_face = QUADS_NUM_EDGES_PER_FACE 
  num_nodes_per_edge = QUADS_NUM_NODES_PER_EDGE 

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!>  @brief   Gets the coordinates of nodes, edges and faces.
!!
!!  @details Places node, edge and face coordinates into the passed arrays.
!!
!!  @param[in]   self                The cubed-sphere function object.
!!  @param[out]  node_coordinates    The node coordinates.
!-------------------------------------------------------------------------------

subroutine get_coordinates(self, node_coordinates)
  implicit none

  class(gencube_type), intent(in) :: self
  real(kind=r_def),   intent(out) :: node_coordinates(:,:)

  !Internal variables
  integer :: iv
  real(kind=r_def) :: long, lat

  do iv = 1, self%nvert(NGRIDS)
    long    = self%vlong(iv,NGRIDS)
    lat     = self%vlat(iv,NGRIDS)
    node_coordinates(:,iv) = [long, lat]
  end do
 
  return
end subroutine get_coordinates

!-------------------------------------------------------------------------------
!>  @brief   Gets the connectivity information of nodes, edges and faces.
!!
!!  @details Places face-node, edge-node and face-face connectivity to the
!!           passed arrays. The left-most loop index covers e.g. all nodes
!!           around one face.
!!
!!  @param[in]   self                     The cubed-sphere function object.
!!  @param[out]  face_node_connectivity   Face-node connectivity.
!!  @param[out]  edge_node_connectivity   Edge-node connectivity.
!!  @param[out]  face_edge_connectivity   Face-edge connectivity.
!!  @param[out]  face_face_connectivity   Face-face connectivity.
!-------------------------------------------------------------------------------

subroutine get_connectivity(self,                            &
           face_node_connectivity, edge_node_connectivity,   &
           face_edge_connectivity, face_face_connectivity)
  implicit none

  !Arguments
  class(gencube_type), intent(in)  :: self
  integer, intent(out) :: face_node_connectivity(:,:)
  integer, intent(out) :: edge_node_connectivity(:,:) 
  integer, intent(out) :: face_edge_connectivity(:,:)
  integer, intent(out) :: face_face_connectivity(:,:)

  !Internal variables
  integer :: iface, iedge, inode
  integer :: iface1, iface2

  !Face-face connectivity
  do iface2 = 1, self%nface(NGRIDS)
    do iface1 = 1, QUADS_NUM_FACE_NEIGHBOURS
      face_face_connectivity(iface1, iface2) = self%fnxtf(iface2, iface1, NGRIDS)
    end do
  end do

  !Face-node connectivity
  do iface = 1, self%nface(NGRIDS)
    do inode = 1, QUADS_NUM_NODES_PER_FACE
      face_node_connectivity(inode, iface) = self%voff(iface, inode, NGRIDS)
    end do
  end do

  !Edge-node connectivity
  do iedge = 1, self%nedge(NGRIDS)
    do inode = 1, QUADS_NUM_NODES_PER_EDGE
      edge_node_connectivity(inode, iedge) = self%vofe(iedge, inode, NGRIDS)
    end do
  end do

  !Face-edge connectivity
  do iface = 1, self%nface(NGRIDS)
    do iedge = 1, QUADS_NUM_EDGES_PER_FACE
      face_edge_connectivity(iedge, iface) = self%eoff(iface, iedge, NGRIDS)
    end do
  end do

  return
end subroutine get_connectivity

!-------------------------------------------------------------------------------
!>  @brief     Generate a cubed-sphere mesh.
!!
!!  @details   Implements the deferred generate method to generate a
!!             cubed-sphere mesh, as configured in the constructor. Calls down
!!             to other subroutines in this module in the correct sequence.
!!
!!  @param[in,out]   self   The cubed-sphere function object.
!-------------------------------------------------------------------------------

subroutine generate(self)
  implicit none

  !Arguments
  class(gencube_type), intent(inout)  :: self

  !Internal variables
  integer :: igrid

  integer          :: n, n2
  real(kind=r_def) :: dlambda

  do igrid = 1, NGRIDS

    !Size of panels on this grid
    n       = self%n0*(2**(igrid-1))
    n2      = n*n
    dlambda = 0.5_r_def*PI/n

    self%nface(igrid) = 6*n2
    self%nedge(igrid) = 2*self%nface(igrid)

    self%nvert(igrid) = self%nface(igrid) + 2

    call part1(self, igrid, dlambda, n, n2)
    call part2(self, igrid, n, n2)
    call part3(self, igrid, n, n2)

  end do ! end of main loop over grids

  call part4(self)
  call part5(self)
  call part6(self)
  call part7(self)
  call part8(self)

  return
end subroutine generate

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part1.
!!
!!  @details   First unit of work to generate the mesh. Called by the
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!!  @param[in]      igrid    Mesh index in the multi-grid hierarchy.
!!  @param[in]      n        
!!  @param[in]      n2        
!!  @param[in]      dlambda  
!-------------------------------------------------------------------------------

subroutine part1(self, igrid, dlambda, n, n2)
  use coord_algorithms_mod, only: xyz2ll
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self
  integer,          intent(in)       :: igrid
  integer,          intent(in)       :: n
  integer,          intent(in)       :: n2
  real(kind=r_def), intent(in)       :: dlambda

  !Internal variables
  real(kind=r_def) :: lambda1, lambda2
  real(kind=r_def) :: lat, long
  real(kind=r_def) :: t1, t2
  real(kind=r_def) :: x1, y1, z1
  real(kind=r_def) :: x2, y2, z2
  real(kind=r_def) :: x3, y3, z3
  real(kind=r_def) :: x4, y4, z4
  real(kind=r_def) :: x5, y5, z5
  real(kind=r_def) :: x6, y6, z6

  integer          :: i, j, ixv, p1

  !Loop over vertices/faces of one panel
  do j = 1, n
    lambda2 = (j-1)*dlambda - PIBY4
    t2 = tan(lambda2)

    do i = 1, n
      lambda1 = (i-1)*dlambda - PIBY4
      t1 = tan(lambda1)

      !Set up coordinates of vertices
      !Panel 1 ...

      !Index of vertex
      ixv = (j-1)*n + i

      !Cartesian coordinates of vertex
      x1 = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
      y1 = x1*t1
      z1 = x1*t2

      !Lat long coordinates of vertex
      call xyz2ll(x1,y1,z1,long,lat)
      self%vlong(ixv,igrid) = long
      self%vlat(ixv,igrid) = lat

      !Panel 2 ...
      !Index of vertex
      ixv = ixv + n2

      !Cartesian coordinates of vertex
      x2 = -y1
      y2 = x1
      z2 = z1

      !Lat long coordinates of vertex
      call xyz2ll(x2,y2,z2,long,lat)
      self%vlong(ixv,igrid) = long
      self%vlat(ixv,igrid) = lat

      !Panel 3 ...
      !Index of vertex
      ixv = ixv + n2

      !Cartesian coordinates of vertex
      x3 = x2
      y3 = -z2
      z3 = y2

      !Lat long coordinates of vertex
      call xyz2ll(x3,y3,z3,long,lat)
      self%vlong(ixv,igrid) = long
      self%vlat(ixv,igrid) = lat

      !Panel 4
      !Index of vertex
      ixv = ixv + n2

      !Cartesian coordinates of vertex
      x4 = -z3
      y4 = y3
      z4 = x3

      !Lat long coordinates of vertex
      call xyz2ll(x4,y4,z4,long,lat)
      self%vlong(ixv,igrid) = long
      self%vlat(ixv,igrid)  = lat

      !Panel 5
      !Index of vertex
      ixv = ixv + n2

      !Cartesian coordinates of vertex
      x5 = -y4
      y5 = x4
      z5 = z4

      !Lat long coordinates of vertex
      call xyz2ll(x5,y5,z5,long,lat)
      self%vlong(ixv,igrid) = long
      self%vlat(ixv,igrid)  = lat

      !Panel 6
      !Index of vertex
      ixv = ixv + n2

      !Cartesian coordinates of vertex
      x6 = x5
      y6 = -z5
      z6 = y5

      !Lat long coordinates of vertex
      call xyz2ll(x6,y6,z6,long,lat)
      self%vlong(ixv,igrid) = long
      self%vlat(ixv,igrid)  = lat

      !Set up incidence tables ignoring complications at
      !panel edges
      do p1 = 1, 6
        ixv = (p1 - 1)*n2 + (j-1)*n + i
        !Edges of the face
        self%eoff(ixv,1,igrid) = 2*ixv - 1
        self%eoff(ixv,2,igrid) = 2*ixv
        self%eoff(ixv,3,igrid) = 2*ixv + 1
        self%eoff(ixv,4,igrid) = 2*ixv + 2*n

        !Vertices of the face
        self%voff(ixv,1,igrid) = ixv
        self%voff(ixv,2,igrid) = ixv + 1
        self%voff(ixv,3,igrid) = ixv + n + 1
        self%voff(ixv,4,igrid) = ixv + n

        !Faces neighboring this face
        self%fnxtf(ixv,1,igrid) = ixv - 1
        self%fnxtf(ixv,2,igrid) = ixv - n
        self%fnxtf(ixv,3,igrid) = ixv + 1
        self%fnxtf(ixv,4,igrid) = ixv + n

        !Edges incident on the vertex
        self%eofv(ixv,1,igrid) = 2*ixv - 2
        self%eofv(ixv,2,igrid) = 2*(ixv-n) - 1
        self%eofv(ixv,3,igrid) = 2*ixv 
        self%eofv(ixv,4,igrid) = 2*ixv - 1
      end do

    end do
  end do

  return
end subroutine part1

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part2.
!!
!!  @details   Second unit of work to generate the mesh. Called by the
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!!  @param[in]      igrid    Mesh index in the multi-grid hierarchy.
!!  @param[in]      n        
!!  @param[in]      n2        
!-------------------------------------------------------------------------------

subroutine part2(self, igrid, n, n2)
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self
  integer,             intent(in)    :: igrid
  integer,             intent(in)    :: n, n2

  !Internal variables
  integer :: j
  integer :: jr
  integer :: pp
  integer :: p1, p2
  integer :: ixv

  !Now sort out complications at panel edges
  do j = 1, n
    jr = n + 1 - j

    do pp = 1, 3

      !Odd numbered panels
      p1 = 2*pp - 1

      !Left edge of panel p1 joins to top edge of panel p1 - 2
      !Reverse order
      p2 = modulo(p1 + 3, 6) + 1
      ixv = (p1 - 1)*n2 + n*(j - 1) + 1
      self%fnxtf(ixv,1,igrid) = p2*n2 - n + jr
      self%eofv(ixv,1,igrid)  = 2*p2*n2 - 2*n - 1 + 2*(jr + 1)

      !Bottom edge of panel p1 joins to top edge of panel p1 - 1
      p2 = modulo(p1 + 4, 6) + 1
      ixv = (p1 - 1)*n2 + j
      self%fnxtf(ixv,2,igrid) = p2*n2 - n + j
      self%eofv(ixv,2,igrid)  = 2*p2*n2 - 2*n - 1 + 2*j

      !Right edge of panel p1 joins to left edge of panel p1 + 1
      p2 = modulo(p1, 6) + 1
      ixv = (p1 - 1)*n2 + n*j
      self%eoff(ixv,3,igrid)  = 2*(p2 - 1)*n2 + 2*(j - 1)*n + 1
      self%voff(ixv,2,igrid)  = (p2 - 1)*n2 + (j - 1)*n + 1
      self%voff(ixv,3,igrid)  = (p2 - 1)*n2 + j*n + 1
      self%fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1

      !Top edge of panel p1 joins to left edge of panel p1 + 2
      !Reverse order
      p2 = modulo(p1 + 1, 6) + 1
      ixv = p1*n2 - n + j
      self%eoff(ixv,4,igrid)  = 2*(p2 - 1)*n2 + 2*(jr - 1)*n + 1
      self%voff(ixv,3,igrid)  = (p2 - 1)*n2 + (jr - 1)*n + 1
      self%voff(ixv,4,igrid)  = (p2 - 1)*n2 + jr*n + 1
      self%fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1

      !Even numbered panels
      p1 = 2*pp

      !Left edge of panel p1 joins to right edge of panel p1 - 1
      p2 = modulo(p1 + 4, 6) + 1
      ixv = (p1 - 1)*n2 + n*(j - 1) + 1
      self%fnxtf(ixv,1,igrid) = (p2 - 1)*n2 + n*j
      self%eofv(ixv,1,igrid)  = 2*(p2 - 1)*n2 + 2*n*j

      !Bottom edge of panel p1 joins to right edge of panel p1 - 2
      !Reverse order
      p2 = modulo(p1 + 3, 6) + 1
      ixv = (p1 - 1)*n2 + j
      self%fnxtf(ixv,2,igrid) = (p2 - 1)*n2 + n*jr
      self%eofv(ixv,2,igrid)  = 2*(p2 - 1)*n2 + 2*n*(jr + 1)

      !Top edge of panel p1 joins to bottom edge of panel p1 + 1
      p2 = modulo(p1, 6) + 1
      ixv = p1*n2 - n + j
      self%eoff(ixv,4,igrid)  = 2*(p2 - 1)*n2 + 2*j
      self%voff(ixv,3,igrid)  = (p2 - 1)*n2 + j + 1
      self%voff(ixv,4,igrid)  = (p2 - 1)*n2 + j
      self%fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + j

      !Right edge of panel p1 joins to bottom edge of panel p1 + 2
      !Reverse order
      p2 = modulo(p1 + 1, 6) + 1
      ixv = (p1-1)*n2 + n*j
      self%eoff(ixv,3,igrid)  = 2*(p2 - 1)*n2 + 2*jr
      self%voff(ixv,2,igrid)  = (p2 - 1)*n2 + jr + 1
      self%voff(ixv,3,igrid)  = (p2 - 1)*n2 + jr
      self%fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + jr

    end do
  end do

  return
end subroutine part2

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part3.
!!
!!  @details   Third unit of work to generate the mesh. Called by the
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!!  @param[in]      igrid    Mesh index in the multi-grid hierarchy.
!!  @param[in]      n        
!!  @param[in]      n2        
!-------------------------------------------------------------------------------

subroutine part3(self, igrid, n, n2)
  use coord_algorithms_mod, only: xyz2ll
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self
  integer,             intent(in)    :: igrid
  integer,             intent(in)    :: n
  integer,             intent(in)    :: n2

  !Internal variables
  integer          :: pp, p1
  integer          :: ixv, iv
  real(kind=r_def) :: t1, t2

  real(kind=r_def) :: lambda1
  real(kind=r_def) :: lambda2

  real(kind=r_def) :: lat, long

  real(kind=r_def) :: x1, y1, z1

  !All faces have 4 edges and vertices
  self%neoff(1:self%nface(igrid),igrid) = 4

  !Almost all vertices have 4 edges (exceptions dealt with below)
  self%neofv(1:self%nvert(igrid),igrid) = 4

  !Vertices not correctly captured by the above
  !Corner vertices only have 3 edges
  do pp = 1, 3
    !Bottom left of odd numbered panels
    p1 = 2*pp - 1
    ixv = (p1 - 1)*n2 + 1
    !First edge needs to be deleted
    self%eofv(ixv,1,igrid) = self%eofv(ixv,2,igrid)
    self%eofv(ixv,2,igrid) = self%eofv(ixv,3,igrid)
    self%eofv(ixv,3,igrid) = self%eofv(ixv,4,igrid)
    self%eofv(ixv,4,igrid) = 0
    self%neofv(ixv,igrid)  = 3

    !Bottom left of even numbered panels
    p1 = 2*pp
    ixv = (p1 - 1)*n2 + 1

    !Second edge needs to be deleted
    self%eofv(ixv,2,igrid) = self%eofv(ixv,3,igrid)
    self%eofv(ixv,3,igrid) = self%eofv(ixv,4,igrid)
    self%eofv(ixv,4,igrid) = 0
    self%neofv(ixv,igrid)  = 3
  end do

  !Vertex 6*n2 + 1 is at top left of panels 1, 3, and 5
  iv      = 6*n2 + 1
  lambda2 = PIBY4
  t2      = tan(lambda2)
  lambda1 = - PIBY4
  t1      = tan(lambda1)

  !Cartesian coordinates of vertex
  x1 = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
  y1 = x1*t1
  z1 = x1*t2

  !Lat long coordinates of vertex
  call xyz2ll(x1,y1,z1,long,lat)
  self%vlong(iv,igrid) = long
  self%vlat(iv,igrid) = lat

  do pp = 1, 3
    p1 = 2*pp - 1
    ixv = p1*n2 - n + 1
    self%voff(ixv,4,igrid) = iv
    self%eofv(iv,pp,igrid) = 2*p1*n2 - 2*n + 1
  end do

  self%neofv(iv,igrid) = 3

  !Vertex 6*n2 + 2 is at bottom right of panels 2, 4, and 6
  iv = 6*n2 + 2
  x1 = -x1
  y1 = -y1
  z1 = -z1

  !Lat long coordinates of vertex
  call xyz2ll(x1,y1,z1,long,lat)

  self%vlong(iv,igrid) = long
  self%vlat(iv,igrid) = lat
  do pp = 1, 3
    p1 = 2*pp
    ixv = (p1 - 1)*n2 + n
    self%voff(ixv,2,igrid) = iv
    self%eofv(iv,pp,igrid) = 2*(p1 - 1)*n2 + 2*n
  end do

  self%neofv(iv,igrid) = 3

  return
end subroutine part3

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part4.
!!
!!  @details   Fourth unit of work to generate the mesh. Called by the
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!!  @param[in]      igrid    Mesh index in the multi-grid hierarchy.
!!  @param[in]      i        
!!  @param[in]      j        
!!  @param[in]      ixv        
!!  @param[in]      ie1        
!-------------------------------------------------------------------------------

subroutine part4(self)
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self

  !Internal variables
  integer :: igrid
  integer :: i, j
  integer :: ixv
  integer :: ie1
 
  !Now construct inverse tables
  !First initialize entries to zero
  self%fnxte = 0
  self%vofe  = 0
  self%fofv  = 0

  do igrid = 1, NGRIDS
    do j = 1, self%nface(igrid)
      do i = 1, 4
        ie1 = self%eoff(j,i,igrid)
        call addtab(self%fnxte(1,1,igrid),ie1,j,self%nedgex,2)
        ixv = self%voff(j,i,igrid)
        call addtab(self%fofv(1,1,igrid),ixv,j,self%nvertx,4)
      end do
    end do

    do j = 1, self%nvert(igrid)
      do i = 1, 4
        ixv = self%eofv(j,i,igrid)
        if (ixv > 0) then
          call addtab(self%vofe(1,1,igrid),ixv,j,self%nedgex,2)
        end if
      end do
    end do
  end do

  return
end subroutine part4

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part5.
!!
!!  @details   Fifth unit of work to generate the mesh. Called by the
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!-------------------------------------------------------------------------------

subroutine part5(self)
  use coord_algorithms_mod, only : xyz2ll, ll2xyz, starea2, spdist
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self

  !Internal variables
  integer :: igrid, ismooth

  integer :: if1, if2, ixv,  ie0, ie1, i
  integer :: iv1, iv2

  real(kind=r_def) :: lat, long
  real(kind=r_def) :: rmag

  real(kind=r_def) :: xc, yc, zc
  real(kind=r_def) :: x0, y0, z0
  real(kind=r_def) :: x1, y1, z1
  real(kind=r_def) :: x2, y2, z2

  real(kind=r_def) :: aface, atri

  real(kind=r_def) :: lmn
  real(kind=r_def) :: lmx
  real(kind=r_def) :: dmn
  real(kind=r_def) :: dmx
  real(kind=r_def) :: dav

  real(kind=r_def) :: s

  !Calculate geometrical quantities
  do igrid = 1, NGRIDS

    !Smoothing iterations
    do ismooth = 1, NSMOOTH

      !First locate face centres at barycentres of 
      !surrounding vertices
      do if1 = 1, self%nface(igrid)
        xc = 0.0_r_def
        yc = 0.0_r_def
        zc = 0.0_r_def
        do i = 1, 4
          ixv  = self%voff(if1,i,igrid)
          long = self%vlong(ixv,igrid)
          lat  = self%vlat(ixv,igrid)
          call ll2xyz(long,lat,x1,y1,z1)
          xc = xc + x1
          yc = yc + y1
          zc = zc + z1
        end do
        rmag = 1.0_r_def/sqrt(xc*xc + yc*yc + zc*zc)
        xc = xc*rmag
        yc = yc*rmag
        zc = zc*rmag
        call xyz2ll(xc,yc,zc,long,lat)
        self%flong(if1,igrid) = long
        self%flat(if1,igrid)  = lat
      end do

      !next relocate vertices at barycentres of 
      !surrounding face centres - needed for h operator
!      do iv1 = 1, self%nvert(igrid)
!        xc = 0.0_r_def
!        yc = 0.0_r_def
!        zc = 0.0_r_def
!        do i = 1, self%neofv(iv1,igrid)
!          if1  = self%fofv(iv1,i,igrid)
!          long = self%flong(if1,igrid)
!          lat  = self%flat(if1,igrid)
!          call ll2xyz(long,lat,x1,y1,z1)
!          xc = xc + x1
!          yc = yc + y1
!          zc = zc + z1
!        enddo
!        rmag = 1.0_r_def/sqrt(xc*xc + yc*yc + zc*zc)
!        xc = xc*rmag
!        yc = yc*rmag
!        zc = zc*rmag
!        call xyz2ll(xc,yc,zc,long,lat)
!        self%vlong(iv1,igrid) = long
!        self%vlat(iv1,igrid) = lat
!      end do

    end do

    !Tabulate areas
    do if1 = 1, self%nface(igrid)
      long = self%flong(if1,igrid)
      lat  = self%flat(if1,igrid)
      call ll2xyz(long,lat,x0,y0,z0)
      !Compute face area
      aface = 0.0_r_def
      do i = 1, 4

        ie1  = self%eoff(if1,i,igrid)
        ixv  = self%vofe(ie1,1,igrid)
        long = self%vlong(ixv,igrid)
        lat  = self%vlat(ixv,igrid)
        call ll2xyz(long,lat,x1,y1,z1)

        ixv  = self%vofe(ie1,2,igrid)
        long = self%vlong(ixv,igrid)
        lat  = self%vlat(ixv,igrid)
        call ll2xyz(long,lat,x2,y2,z2)
        call starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,atri)

        aface = aface + atri
      end do

      self%farea(if1,igrid) = aface
    end do

    !Tabulate lengths of edges and distances between face centres
    !across each edge
    lmn = 5.0_r_def
    lmx = 0.0_r_def
    dmn = 5.0_r_def
    dmx = 0.0_r_def
    dav = 0.0_r_def

    do ie0 = 1, self%nedge(igrid)

      !Vertices at ends of this edge
      iv1  = self%vofe(ie0,1,igrid)
      iv2  = self%vofe(ie0,2,igrid)
      long = self%vlong(iv1,igrid)
      lat  = self%vlat(iv1,igrid)
      call ll2xyz(long,lat,x1,y1,z1)

      long = self%vlong(iv2,igrid)
      lat  = self%vlat(iv2,igrid)
      call ll2xyz(long,lat,x2,y2,z2)

      call spdist(x1,y1,z1,x2,y2,z2,s)

      self%ldist(ie0,igrid) = s
      lmn = min(lmn,self%ldist(ie0,igrid))
      lmx = min(lmx,self%ldist(ie0,igrid))

      !faces either side of this edge
      if1 = self%fnxte(ie0,1,igrid)
      if2 = self%fnxte(ie0,2,igrid)
      long = self%flong(if1,igrid)
      lat  = self%flat(if1,igrid)
      call ll2xyz(long,lat,x1,y1,z1)
      long = self%flong(if2,igrid)
      lat  = self%flat(if2,igrid)
      call ll2xyz(long,lat,x2,y2,z2)
      call spdist(x1,y1,z1,x2,y2,z2,s)
      self%ddist(ie0,igrid) = s
      dmn = min(dmn,self%ddist(ie0,igrid))
      dmx = min(dmx,self%ddist(ie0,igrid))
      dav = dav + self%ddist(ie0,igrid)/self%nedge(igrid)
    end do
  end do

  return
end subroutine part5

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part6.
!!
!!  @details   Sixth unit of work to generate the mesh. Called by the
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!-------------------------------------------------------------------------------

subroutine part6(self)
  use coord_algorithms_mod, only : ll2xyz
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self

  !Internal variables
  integer          :: igrid
  logical          :: lfound

  integer          :: if0, if1, if2, if3
  integer          :: ix1, ix2
  integer          :: ie1
  integer          :: ixmin, ifmin

  integer          :: if21, if22
  integer          :: iv0

  real(kind=r_def) :: x0, y0, z0
  real(kind=r_def) :: x1, y1, z1
  real(kind=r_def) :: x2, y2, z2

  real(kind=r_def) :: d1x, d1y, d1z
  real(kind=r_def) :: d2x, d2y, d2z

  real(kind=r_def) :: long, lat
  
  real(kind=r_def) :: thetamin, theta

  real(kind=r_def) :: sn, cs


  !Sort FNXTF into anticlockwise order on each grid
  !and sort EOFF to correspond to FNXTF
  !Also sort fofv into anticlockwise order
  do igrid = 1, NGRIDS
    do if0 = 1, self%nface(igrid)

      !Coordinates of face if0
      long = self%flong(if0,igrid)
      lat  = self%flat (if0,igrid)
      call ll2xyz(long,lat,x0,y0,z0)

      do ix1 = 1, 2

        !Coordinates of IX1'th neighbour
        if1  = self%fnxtf(if0,ix1,igrid)
        long = self%flong(if1,igrid)
        lat  = self%flat(if1,igrid)
        call ll2xyz(long,lat,x1,y1,z1)

        d1x = x1 - x0
        d1y = y1 - y0
        d1z = z1 - z0

        !Find next neighbour (anticlockwise)
        thetamin = PI
        ixmin = 0
        ifmin = 0

        do ix2 = ix1 + 1, 4

          !Coordinates of IX2'th neighbour
          if2  = self%fnxtf(if0,ix2,igrid)
          long = self%flong(if2,igrid)
          lat  = self%flat(if2,igrid)
          CALL ll2xyz(long,lat,x2,y2,z2)

          d2x=x2 - x0
          d2y=y2 - y0
          d2z=z2 - z0

          cs = d1x*d2x + d1y*d2y + d1z*d2z
          sn = x0*(d1y*d2z - d1z*d2y) &
             + y0*(d1z*d2x - d1x*d2z) &
             + z0*(d1x*d2y - d1y*d2x)
          theta = atan2(sn,cs)

          if ((theta < thetamin) .and. (theta > 0.0_r_def)) then
            ixmin    = ix2
            ifmin    = if2
            thetamin = theta
          end if
        end do

        !The face in position IXMIN belongs in position IX1+1 so swap them
        if3 = self%fnxtf(if0,ix1+1,igrid)
        self%fnxtf(if0,ix1+1,igrid) = ifmin

        self%fnxtf(if0,ixmin,igrid) = if3

      end do

      do ix1 = 1, 4
        if1 = self%fnxtf(if0,ix1,igrid)
        ix2 = ix1 - 1
        lfound = .false.

        do while (.not. lfound)
          ix2  = ix2 + 1
          ie1  = self%eoff(if0,ix2,igrid)
          if21 = self%fnxte(ie1,1,igrid)
          if22 = self%fnxte(ie1,2,igrid)
          if ((if21 + if22) == (if0 + if1)) lfound = .true.
        end do
  !     edge ie2 corresponds to face if1
        self%eoff(if0,ix2,igrid) = self%eoff(if0,ix1,igrid)
        self%eoff(if0,ix1,igrid) = ie1
      end do

    end do

    do iv0 = 1, self%nvert(igrid)

      !coordinates of vertex iv0
      long = self%vlong(iv0,igrid)
      lat  = self%vlat(iv0,igrid)
      call ll2xyz(long,lat,x0,y0,z0)
      do ix1 = 1, self%neofv(iv0,igrid) - 2

        !coordinates of ix1'th face
        if1  = self%fofv(iv0,ix1,igrid)
        long = self%flong(if1,igrid)
        lat  = self%flat(if1,igrid)

        call ll2xyz(long,lat,x1,y1,z1)
        d1x = x1 - x0
        d1y = y1 - y0
        d1z = z1 - z0

        !find next neighbour (anticlockwise)
        thetamin = PI
        ixmin = 0
        ifmin = 0

        do ix2 = ix1 + 1, self%neofv(iv0,igrid)
          !coordinates of ix2'th neighbour
          if2  = self%fofv(iv0,ix2,igrid)
          long = self%flong(if2,igrid)
          lat  = self%flat(if2,igrid)
          call ll2xyz(long,lat,x2,y2,z2)

          d2x=x2 - x0
          d2y=y2 - y0
          d2z=z2 - z0

          cs = d1x*d2x + d1y*d2y + d1z*d2z
          sn = x0*(d1y*d2z - d1z*d2y) &
             + y0*(d1z*d2x - d1x*d2z) &
             + z0*(d1x*d2y - d1y*d2x)
          theta = atan2(sn,cs)
          if ((theta < thetamin) .and. (theta > 0.0_r_def)) then
            ixmin = ix2
            ifmin = if2
            thetamin = theta
          end if
        end do

        !the face in position ixmin belongs in position ix1+1 so swap them
        if3 = self%fofv(iv0,ix1+1,igrid)
        self%fofv(iv0,ix1+1,igrid) = ifmin
        self%fofv(iv0,ixmin,igrid) = if3

      end do
    end do

  end do

  return
end subroutine part6

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part7.
!!
!!  @details   Seventh unit of work to generate the mesh. Called by the 
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!-------------------------------------------------------------------------------

subroutine part7(self)
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self

  !Internal variables
  integer :: iv0, iv11, iv12, iv21, iv22
  integer :: if0
  integer :: ix1, ix2
  integer :: ie1, ie2

  integer :: igrid

  !Order VOFF so that the k'th vertex is between the
  !k'th and (k+1)'th edges in EOFF
  do igrid = 1, NGRIDS
    do if0 = 1, self%nface(igrid)
      do ix1 = 1, self%neoff(if0,igrid)
        ix2 = ix1 + 1
        if (ix2 > self%neoff(if0,igrid)) ix2 = 1
        ie1 = self%eoff(if0,ix1,igrid)
        ie2 = self%eoff(if0,ix2,igrid)

        !Find the common vertex of ie1 and ie2
        iv11 = self%vofe(ie1,1,igrid)
        iv12 = self%vofe(ie1,2,igrid)
        iv21 = self%vofe(ie2,1,igrid)
        iv22 = self%vofe(ie2,2,igrid)

        if ((iv11 == iv21) .or. (iv11 == iv22)) then
          iv0 = iv11
        else if ((iv12 == iv21) .or. (iv12 == iv22)) then
          iv0 = iv12
        else
          print *,'Common vertex not found'
          stop
        end if
        self%voff(if0,ix1,igrid) = iv0
      end do
    end do
  end do

  return
end subroutine part7

!-------------------------------------------------------------------------------
!>  @brief     Mesh generation part8.
!!
!!  @details   Eighth unit of work to generate the mesh. Called by the 
!!             generate method. 
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!-------------------------------------------------------------------------------

subroutine part8(self)
  use coord_algorithms_mod, only: ll2xyz
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self

  !Internal variables
  real(kind=r_def) :: sn

  real(kind=r_def) :: x0,  y0,  z0
  real(kind=r_def) :: x1,  y1,  z1
  real(kind=r_def) :: d1x, d1y, d1z
  real(kind=r_def) :: d2x, d2y, d2z
  real(kind=r_def) :: lat, long

  integer          :: igrid
  integer          :: ie0

  integer          :: if1, if2
  integer          :: iv1, iv2

  !Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
  !is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)

  do igrid = 1, NGRIDS
    do ie0 = 1, self%nedge(igrid)

      if1  = self%fnxte(ie0,1,igrid)
      if2  = self%fnxte(ie0,2,igrid)
      long = self%flong(if1,igrid)
      lat  = self%flat(if1,igrid)
      call ll2xyz(long,lat,x0,y0,z0)

      long = self%flong(if2,igrid)
      lat  = self%flat(if2,igrid)
      call ll2xyz(long,lat,x1,y1,z1)
      d1x  = x1 - x0
      d1y  = y1 - y0
      d1z  = z1 - z0
      iv1  = self%vofe(ie0,1,igrid)
      iv2  = self%vofe(ie0,2,igrid)
      long = self%vlong(iv1,igrid)
      lat  = self%vlat(iv1,igrid)

      call ll2xyz(long,lat,x0,y0,z0)
      long = self%vlong(iv2,igrid)
      lat  = self%vlat(iv2,igrid)
      call ll2xyz(long,lat,x1,y1,z1)

      d2x = x1 - x0
      d2y = y1 - y0
      d2z = z1 - z0

      sn = x0*(d1y*d2z - d1z*d2y) &
         + y0*(d1z*d2x - d1x*d2z) &
         + z0*(d1x*d2y - d1y*d2x)

      if (sn < 0.0_r_def) then
        !swap the two vertices
        self%vofe(ie0,1,igrid) = iv2
        self%vofe(ie0,2,igrid) = iv1
      end if
    end do
  end do

  return
end subroutine part8

!-------------------------------------------------------------------------------
!>  @brief     Write mesh data to rough-and-ready files.
!!
!!  @details   Write mesh data to files, largely for code testing purposes.
!!             File unit numbers are hard-coded for the moment.
!!
!!  @param[in,out]  self     The cubed-sphere function object.
!-------------------------------------------------------------------------------

subroutine write_data(self)
  use coord_algorithms_mod, only: ll2xyz
  implicit none

  !Arguments
  class(gencube_type), intent(inout) :: self

  integer          :: iv, if1, if2, j
  integer          :: ie1, iv1, iv2
  integer          :: igrid

  real(kind=r_def) :: long, lat

  real(kind=r_def) :: x1, y1, z1
  real(kind=r_def) :: x2, y2, z2

  character(len=str_def) :: ygridfile

  !Write out coordinates of edges for plotting
  open(44,file='primalgrid.dat',form='formatted')

  do j = 1, self%nedgex
    iv   = self%vofe(j,1,NGRIDS)
    long = self%vlong(iv,NGRIDS)
    lat  = self%vlat(iv,NGRIDS)
    call ll2xyz(long,lat,x1,y1,z1)

    iv   = self%vofe(j,2,NGRIDS)
    long = self%vlong(iv,NGRIDS)
    lat  = self%vlat(iv,NGRIDS)
    call ll2xyz(long,lat,x2,y2,z2)

    write(44,'(3e15.7)') x1,y1,z1
    write(44,'(3e15.7)') x2,y2,z2
  end do

  close(44)

  open(44,file='dualgrid.dat',form='formatted')
  do j = 1, self%nedgex
    if1  = self%fnxte(j,1,NGRIDS)
    long = self%flong(if1,NGRIDS)
    lat  = self%flat(if1,NGRIDS)
    call ll2xyz(long,lat,x1,y1,z1)

    if1  = self%fnxte(j,2,NGRIDS)
    long = self%flong(if1,NGRIDS)
    lat  = self%flat(if1,NGRIDS)
    call ll2xyz(long,lat,x2,y2,z2)
    write(44,'(3e15.7)') x1,y1,z1
    write(44,'(3e15.7)') x2,y2,z2
  end do

  close(44)

  !output gridmap file
  write(ygridfile,'(''gridmap_cube_'',i10.10,''.dat'')') self%nfacex
  open(22,file=trim(ygridfile),form='unformatted')
  
  ! write(22,*) 'gridmap for NGRIDS=',NGRIDS
  write(22) NGRIDS
  write(22) self%nface
  write(22) self%nedge
  write(22) self%nvert

  ! write(22,*) 'number of edges of each face - all grids'
  write(22) ((self%neoff(if1,igrid),            &
                 if1 = 1, self%nface(igrid)),   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'number of edges of each vertex - all grids'
  write(22) ((self%neofv(iv1,igrid),            &
                 iv1 = 1, self%nvert(igrid)),   &
                 igrid=1, NGRIDS)
  ! write(22,*) 'faces next to each face - all grids'
  write(22) (((self%fnxtf(if1,if2,igrid),       &
                 if1 = 1, self%nface(igrid)),   &
                 if2 = 1, 4),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'edges of each face - all grids'
  write(22) (((self%eoff(if1,ie1,igrid),        &
                 if1 = 1, self%nface(igrid)),   &
                 ie1 = 1, 4),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'vertices of each face - all grids'
  write(22) (((self%voff(if1,iv1,igrid),        &
                 if1 = 1, self%nface(igrid)),   &
                 iv1 = 1, 4),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'faces next to each edge - all grids'
  write(22) (((self%fnxte(ie1,if2,igrid),       &
                 ie1 = 1, self%nedge(igrid)),   &
                 if2 = 1, 2),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'vertices of each edge - all grids'
  write(22) (((self%vofe(ie1,iv2,igrid),        &
                 ie1 = 1, self%nedge(igrid)),   &
                 iv2 = 1, 2),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'faces around each vertex - all grids'
  write(22) (((self%fofv(iv1,if2,igrid),        &
                 iv1 = 1, self%nvert(igrid)),   &
                 if2 = 1, 4),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'edges around each vertex - all grids'
  write(22) (((self%eofv(iv1,ie1,igrid),        &
                 iv1 = 1, self%nvert(igrid)),   &
                 ie1 = 1, 4),                   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'longitudes of faces - all grids'
  write(22) ((self%flong(if1,igrid),            &
                 if1 = 1, self%nface(igrid)),   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'latitudes of faces - all grids'
  write(22) ((self%flat(if1,igrid),             &
                 if1 = 1, self%nface(igrid)),   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'longitudes of vertices - all grids'
  write(22) ((self%vlong(iv1,igrid),            &
                 iv1 = 1, self%nvert(igrid)),   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'latitudes of vertices - all grids'
  write(22) ((self%vlat(iv1,igrid),             &
                 iv1 = 1, self%nvert(igrid)),   &
                 igrid = 1, NGRIDS)
  ! write(22,*) 'areas of faces - all grids'
  write(22) ((self%farea(if1,igrid),            &
                if1 = 1, self%nface(igrid)),    &
                igrid = 1, NGRIDS)
  ! write(22,*) 'lengths of edges - all grids'
  write(22) ((self%ldist(ie1,igrid),            &
                ie1 = 1, self%nedge(igrid)),    &
                igrid = 1, NGRIDS)
  ! write(22,*) 'distance between faces across edges - all grids'
  write(22) ((self%ddist(ie1,igrid),            &
                ie1 = 1, self%nedge(igrid)),    &
                igrid = 1, NGRIDS)

  return
end subroutine write_data

!-------------------------------------------------------------------------------
!>  @brief     Add an entry to a table.
!!
!!  @param[in]      dim1 
!!  @param[in]      dim2 
!!  @param[in]      index 
!!  @param[in]      entry 
!!  @param[in,out]  tab 
!-------------------------------------------------------------------------------

subroutine addtab(tab,index,entry,dim1,dim2)
  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  implicit none

  !Arguments
  integer, intent(in)    :: dim1,dim2
  integer, intent(in)    :: index
  integer, intent(in)    :: entry
  integer, intent(inout) :: tab(dim1,dim2)

  !Internal
  integer :: i

  i=1
  do while (tab(index,i)/=0)
    if (i>dim2) then
      call log_event('**********', LOG_LEVEL_INFO)
      call log_event('TABLE FULL', LOG_LEVEL_INFO)
      call log_event('**********', LOG_LEVEL_INFO)
      write(log_scratch_space,*) index,entry,dim1,dim2
      call log_event(trim(log_scratch_space), LOG_LEVEL_INFO)
      write(log_scratch_space,*) tab(index,:)
      call log_event(trim(log_scratch_space), LOG_LEVEL_INFO)
      call log_event('Terminating execution.', LOG_LEVEL_ERROR)
      stop
    end if
    i=i+1
  end do

  tab(index,i)=entry

  return
end subroutine addtab

end module gencube_mod
