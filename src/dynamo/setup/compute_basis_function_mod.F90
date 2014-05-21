!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Computes the basis functions on quadrature points for 4 element spaces 
!-------------------------------------------------------------------------------
module compute_basis_function_mod

  use num_dof_mod
  use reference_element_mod

  use constants_mod, only: dp
  use gaussian_quadrature_mod, only: gaussian_quadrature_type, &
                                     ngp_v, ngp_h !parameter for how many GQ points
  use function_space_mod, only : function_space_type

  implicit none

contains 

  subroutine compute_basis(k,v0,v1,v2,v3, v_unique_dofs,v_dof_entity,gq)
    !-----------------------------------------------------------------------------
    ! Subroutine to compute test/trial functions on quadrature points
    !-----------------------------------------------------------------------------

    ! order of elements
    integer, intent(in) :: k
    type(function_space_type), intent(inout) :: v0, v1, v2, v3
    type(gaussian_quadrature_type), intent(in) :: gq
    integer, intent(in) :: v_unique_dofs(4,2), v_dof_entity(4,0:3)

    integer :: i, jx, jy, jz, order, idx, j1, j2, h_ctr
    integer :: j(3), j2l_edge(12,3), j2l_face(6,3), face_idx(6), edge_idx(12,2)
    integer, allocatable :: lx(:), ly(:), lz(:)
    real(kind=dp)    :: fx, fy, fz, gx, gy, gz, dfx, dfy, dfz
    real(kind=dp)    :: x1(k+2), x2(k+1)
    !real(kind=dp)    :: unit_vec_v2(nv2,3), unit_vec_v1(nv1,3)
    real(kind=dp),allocatable    :: unit_vec_v2(:,:), unit_vec_v1(:,:)


    ! Allocate to be larger than should be needed
    allocate ( lx(3*(k+2)**3) )
    allocate ( ly(3*(k+2)**3) )
    allocate ( lz(3*(k+2)**3) )

    allocate(unit_vec_v2(v_unique_dofs(3,2),3))
    allocate(unit_vec_v1(v_unique_dofs(2,2),3))

    ! positional arrays - need two, i.e quadratic and linear for RT1
    do i=1,k+2
      x1(i) = real(i-1)/real(k+1)
    end do
    do i=1,k+1
      x2(i) = real(i-1)/real(k)
    end do
    if ( k == 0 ) x2(1) = 0.5

    ! some look arrays based upon reference cube topology
    face_idx = (/ 1, k+2, k+2, 1, 1, k+2 /)

    edge_idx(:,1) = (/ 1, k+2, k+2, 1, 1, k+2, k+2, 1,   1,   k+2, k+2, 1   /)
    edge_idx(:,2) = (/ 1, 1,   1,   1, 1, 1,   k+2, k+2, k+2, k+2, k+2, k+2 /)

    j2l_face(1,:) = (/ 2, 3, 1 /)
    j2l_face(2,:) = (/ 3, 2, 1 /)
    j2l_face(3,:) = (/ 2, 3, 1 /)
    j2l_face(4,:) = (/ 3, 2, 1 /)
    j2l_face(5,:) = (/ 2, 1, 3 /)
    j2l_face(6,:) = (/ 2, 1, 3 /)

    j2l_edge(1 ,:) = (/ 1, 2, 3 /)
    j2l_edge(2 ,:) = (/ 2, 1, 3 /)
    j2l_edge(3 ,:) = (/ 1, 2, 3 /)
    j2l_edge(4 ,:) = (/ 2, 1, 3 /)
    j2l_edge(5 ,:) = (/ 2, 3, 1 /)
    j2l_edge(6 ,:) = (/ 2, 3, 1 /)
    j2l_edge(7 ,:) = (/ 2, 3, 1 /)
    j2l_edge(8 ,:) = (/ 2, 3, 1 /)
    j2l_edge(9 ,:) = (/ 1, 2, 3 /)
    j2l_edge(10,:) = (/ 2, 1, 3 /)
    j2l_edge(11,:) = (/ 1, 2, 3 /)
    j2l_edge(12,:) = (/ 2, 1, 3 /)

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of v0 fields
    !-----------------------------------------------------------------------------
    order = k+1

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    do jz=2,k+1
      do jy=2,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=2,k+1
        do j2=2,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on edges
    do i=1,nedges
      do j1=2,k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        idx = idx + 1
      end do
    end do

    ! dofs on vertices
    do i=1,nverts
    !  do j1=1,nv0_vert
      do j1=1,v_dof_entity(1,0)
        lx(idx) =  1+(k+1)*int(x_vert(i,1))
        ly(idx) =  1+(k+1)*int(x_vert(i,2))
        lz(idx) =  1+(k+1)*int(x_vert(i,3))
        idx = idx + 1
      end do
    end do
    do i=1,v_unique_dofs(1,2)
    !do i=1,nv0
       ! explicitly for quads, as ngp_h = ngp_v * ngp_v
       h_ctr = 1
       do jx=1,ngp_v
          fx = gq%poly1d(order,jx,x1(lx(i)),x1,lx(i))
          dfx = gq%poly1d_deriv(order,jx,x1(lx(i)),x1,lx(i))
          do jy=1,ngp_v
             fy = gq%poly1d(order,jy,x1(ly(i)),x1,ly(i))
             dfy = gq%poly1d_deriv(order,jy,x1(ly(i)),x1,ly(i))
             do jz=1,ngp_v
                fz = gq%poly1d(order,jz,x1(lz(i)),x1,lz(i))
                dfz = gq%poly1d_deriv(order,jz,x1(lz(i)),x1,lz(i))
                call v0%set_basis(fx*fy*fz,1,i,h_ctr,jz)
                call v0%set_diff_basis(dfx*fy*fz,1,i,h_ctr,jz)
                call v0%set_diff_basis(fx*dfy*fz,2,i,h_ctr,jz)
                call v0%set_diff_basis(fx*fy*dfz,3,i,h_ctr,jz)                
             end do
             h_ctr = h_ctr + 1 
          end do
       end do
       call v0%set_nodes(x1(lx(i)),x1(ly(i)),x1(lz(i)),i)
    end do

    !-----------------------------------------------------------------------------
    ! section for test/trial functions of v1 fields
    !-----------------------------------------------------------------------------
    order = k+1

    !do idx=1,nv1
    do idx = 1,v_unique_dofs(2,2)
      do i=1,3
        unit_vec_v1(idx,i) = 0.0
      end do
    end do

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    ! u components
    do jz=2,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v1(idx,1) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=2,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v1(idx,2) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=1,k+1
      do jy=2,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v1(idx,3) = 1.0
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=1,k+1
        do j2=2,k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_v1(idx,:) = tangent_to_edge(edge_on_face(i,1),:)
          idx = idx + 1
        end do
      end do
      do j1=2,k+1
        do j2=1,k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_v1(idx,:) = tangent_to_edge(edge_on_face(i,2),:)
          idx = idx + 1
        end do
      end do  
    end do

    ! dofs on edges
    do i=1,nedges
      do j1=1,k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        unit_vec_v1(idx,:) = tangent_to_edge(i,:)
        idx = idx + 1
      end do
    end do

    ! this needs correcting
    !do i=1,nv1  
    do i=1,v_unique_dofs(2,2)
       ! Quads only as ngp_h = ngp_v * ngp_v
       h_ctr = 1
       do jx=1,ngp_v
          fx = gq%poly1d(order,jx,x1(lx(i)),x1,lx(i))
          dfx = gq%poly1d_deriv(order,jx,x1(lx(i)),x1,lx(i))
          if (lx(i) <= order) then
             gx = gq%poly1d(order-1,jx,x2(lx(i)),x2,lx(i))
          else
             gx = 0.0
          end if
          do jy=1,ngp_v
             fy = gq%poly1d(order,jy,x1(ly(i)),x1,ly(i))
             dfy = gq%poly1d_deriv(order,jy,x1(ly(i)),x1,ly(i))
             if (ly(i) <= order) then 
                gy = gq%poly1d(order-1,jy,x2(ly(i)),x2,ly(i))
             else
                gy = 0.0
             end if
             do jz=1,ngp_v
                fz = gq%poly1d(order,jz,x1(lz(i)),x1,lz(i))
                dfz = gq%poly1d_deriv(order,jz,x1(lz(i)),x1,lz(i))
                if (lz(i) <= order) then     
                   gz = gq%poly1d(order-1,jz,x2(lz(i)),x2,lz(i))
                else
                   gz = 0.0
                end if
                
                call v1%set_basis(gx*fy*fz*unit_vec_v1(i,1),1,i,h_ctr,jz)
                call v1%set_basis(fx*gy*fz*unit_vec_v1(i,2),2,i,h_ctr,jz)
                call v1%set_basis(fx*fy*gz*unit_vec_v1(i,3),3,i,h_ctr,jz)

                call v1%set_diff_basis(                                          &
                     (fx*dfy*gz*unit_vec_v1(i,3) - fx*gy*dfz*unit_vec_v1(i,2) ), &
                     1,i,h_ctr,jz)
                call v1%set_diff_basis(                                          &
                     (gx*fy*dfz*unit_vec_v1(i,1) - dfx*fy*gz*unit_vec_v1(i,3) ), &
                     2,i,h_ctr,jz)
                call v1%set_diff_basis(                                          &
                     (dfx*gy*fz*unit_vec_v1(i,2) - gx*dfy*fz*unit_vec_v1(i,1) ), &
                     3,i,h_ctr,jz)                                          
             end do
             h_ctr = h_ctr + 1
          end do
       end do
       call v1%set_nodes(unit_vec_v2(i,1)*x2(lx(i)) + (1.0 - unit_vec_v1(i,1))*x1(lx(i)), &
            unit_vec_v1(i,2)*x2(ly(i)) + (1.0 - unit_vec_v1(i,2))*x1(ly(i)), &
            unit_vec_v1(i,3)*x2(lz(i)) + (1.0 - unit_vec_v1(i,3))*x1(lz(i)), i)
    end do


    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of v2 fields
    !-----------------------------------------------------------------------------
    order = k + 1

    !do idx=1,nv2
    do idx=1,v_unique_dofs(3,2)
      do i=1,3
        unit_vec_v2(idx,i) = 0.0
      end do
    end do

    idx = 1
    ! dofs in volume
    ! u components
    do jz=1,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v2(idx,1) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=1,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v2(idx,2) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=2,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v2(idx,3) = 1.0
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=1,k+1
        do j2=1,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_v2(idx,:) = normal_to_face(i,:)
          idx = idx + 1
        end do
      end do
    end do

    !do i=1,nv2
    do i=1,v_unique_dofs(3,2)
       ! Quads only as ngp_h = ngp_h * ngp_h
       h_ctr = 1
       do jx=1,ngp_v
          fx = gq%poly1d(order,jx,x1(lx(i)),x1,lx(i))
          dfx = gq%poly1d_deriv(order,jx,x1(lx(i)),x1,lx(i))
          if (lx(i) <= order) then
             gx = gq%poly1d(order-1,jx,x2(lx(i)),x2,lx(i))
          else
             gx = 0.0
          end if
          do jy=1,ngp_v
             fy = gq%poly1d(order,jy,x1(ly(i)),x1,ly(i))
             dfy = gq%poly1d_deriv(order,jy,x1(ly(i)),x1,ly(i))
             if (ly(i) <= order) then
                gy = gq%poly1d(order-1,jy,x2(ly(i)),x2,ly(i))
             else
                gy = 0.0
             end if
             do jz=1,ngp_v
                fz = gq%poly1d(order,jz,x1(lz(i)),x1,lz(i))
                dfz = gq%poly1d_deriv(order,jz,x1(lz(i)),x1,lz(i))
                if (lz(i) <= order) then
                   gz = gq%poly1d(order-1,jz,x2(lz(i)),x2,lz(i))
                else
                   gz = 0.0
                end if
            
                call v2%set_basis(fx*gy*gz*unit_vec_v2(i,1),1,i,h_ctr,jz)
                call v2%set_basis(gx*fy*gz*unit_vec_v2(i,2),2,i,h_ctr,jz)
                call v2%set_basis(gx*gy*fz*unit_vec_v2(i,3),3,i,h_ctr,jz)
                
                call v2%set_diff_basis( (dfx*gy*gz*unit_vec_v2(i,1) & 
                     + gx*dfy*gz*unit_vec_v2(i,2)                   &
                     + gx*gy*dfz*unit_vec_v2(i,3) ),                &
                     1,i,h_ctr,jz) 
             end do
             h_ctr = h_ctr + 1
          end do
       end do
       call v2%set_nodes(unit_vec_v2(i,1)*x1(lx(i)) + (1.0 - unit_vec_v2(i,1))*x2(lx(i)), &
            unit_vec_v2(i,2)*x1(ly(i)) + (1.0 - unit_vec_v2(i,2))*x2(ly(i)), &
            unit_vec_v2(i,3)*x1(lz(i)) + (1.0 - unit_vec_v2(i,3))*x2(lz(i)), i)
    end do

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of v3 fields
    !-----------------------------------------------------------------------------
    order = k
    ! compute indices of functions
    idx = 1
    ! dofs in volume
    do jz=1,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do

    !do i=1,nv3
    ! For Quads only as ngp_h = ngp_v * ngp_v
    do i=1,v_unique_dofs(4,2)
       h_ctr = 1
       do jx=1,ngp_v
          gx = gq%poly1d(order,jx,x2(lx(i)),x2,lx(i))
          do jy=1,ngp_v
             gy = gq%poly1d(order,jy,x2(ly(i)),x2,ly(i))
             do jz=1,ngp_v
                gz = gq%poly1d(order,jz,x2(lz(i)),x2,lz(i))
                call v3%set_basis(gx*gy*gz,1,i,h_ctr,jz)                
             end do
             h_ctr = h_ctr + 1
          end do
       end do
       call v3%set_nodes(x2(lx(i)),x2(ly(i)),x2(lz(i)),i)
    end do

    ! tidy up
    deallocate ( lx )
    deallocate ( ly )
    deallocate ( lz )

    deallocate ( unit_vec_v2, unit_vec_v1)

  end subroutine compute_basis


 

end module compute_basis_function_mod

