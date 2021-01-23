! Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

!*************************************************!
!         https://www.facebook.com/VnCFD          !
!        https://vncfdgroup.wordpress.com         !
!*************************************************!

!compiler: python3 -m numpy.f2py -c functions.f95 fluxes.f95 -m fortranlib

!*************************************************!
!              CELL-FACE TOPO                     !
! node(i,j+1)+----------------+node(i+1,j+1)      !  
!            |  cell*(i,j)    |                   !
! node(i,j)  +----------------+node(i+1,j)        !
!*************************************************!
subroutine cell_topo(node_X, cell_X, all_cell_half_size, cell_volume, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: node_X(n, m, 2)        
  real(8), intent(out) :: cell_X(n-1, m-1, 2)
  real(8), intent(out) :: all_cell_half_size(n+1, m+1, 2)
  real(8), intent(out) :: cell_volume(n-1, m-1)        
  !integer :: i,j

! cell_X(i,j) = (node(i,j)+node(i+1,j)+node(i,j+1)+node(i+1,j+1))/4
  cell_X(1:n-1,1:m-1,:) = 0.25*(node_X(1:n-1,1:m-1,:)+node_X(1:n-1,2:m,:)+node_X(2:n,1:m-1,:)+node_X(2:n,2:m,:))

! cell_volume(i,j) = abs(0.5*cross_product(node(i,j)-node(i+1,j+1), node(i,j+1)-node(i+1,j)))
  cell_volume(1:n-1,1:m-1) = 0.5*abs((node_X(1:n-1,1:m-1,1)-node_X(2:n,2:m,1))*(node_X(1:n-1,2:m,2)-node_X(2:n,1:m-1,2))&
                                  -(node_X(1:n-1,1:m-1,2)-node_X(2:n,2:m,2))*(node_X(1:n-1,2:m,1)-node_X(2:n,1:m-1,1)))

! cell_half_size(i,j)[1] = 0.25*norm2(node(i,j)-node(i+1,j) + node(i,j+1)-node(i+1,j+1))
  all_cell_half_size(1,:,:) = 0.0
  all_cell_half_size(n+1,:,:) = 0.0
  all_cell_half_size(2:n,2:m,1) = 0.25*norm2(node_X(1:n-1,1:m-1,:)-node_X(2:n,1:m-1,:)+node_X(1:n-1,2:m,:)-node_X(2:n,2:m,:), dim=3)

! cell_half_size(i,j)[2] = 0.25*norm2(node(i,j)-node(i,j+1) + node(i+1,j)-node(i+1,j+1))
  all_cell_half_size(:,1,:) = 0.0
  all_cell_half_size(:,m+1,:) = 0.0
  all_cell_half_size(2:n,2:m,2) = 0.25*norm2(node_X(1:n-1,1:m-1,:)-node_X(1:n-1,2:m,:)+node_X(2:n,1:m-1,:)-node_X(2:n,2:m,:), dim=3)

end subroutine cell_topo


!*************************************************!
!       ______+______jface______+_____            !
!             |                 |                 !
!             |iface  ^n        |--->n            !
!       ______+_______|_________+_____            !
!*************************************************!
subroutine face_topo(node_x, iface_X, iface_normal, iface_area, jface_X, jface_normal, jface_area, n, m)
  implicit none
  integer, intent(in) :: n,m
  real(8), intent(in) :: node_X(n, m, 2)
  real(8), intent(out) :: iface_X(n, m-1, 2), iface_normal(n, m-1, 2), iface_area(n, m-1)
  real(8), intent(out) :: jface_X(n-1, m, 2), jface_normal(n-1, m, 2), jface_area(n-1, m)

  iface_X(1:n,1:m-1,:) = 0.5*(node_X(1:n,1:m-1,:)+node_X(1:n,2:m,:))
  iface_area(1:n,1:m-1) = norm2(node_X(1:n,1:m-1,:)-node_X(1:n,2:m,:), dim=3)
  iface_normal(1:n,1:m-1,1) = (-node_X(1:n,1:m-1,2)+node_X(1:n,2:m,2))/iface_area(1:n,1:m-1)
  iface_normal(1:n,1:m-1,2) = ( node_X(1:n,1:m-1,1)-node_X(1:n,2:m,1))/iface_area(1:n,1:m-1)


  jface_X(1:n-1,1:m,:) = 0.5*(node_X(1:n-1,1:m,:)+node_X(2:n,1:m,:))
  jface_area(1:n-1,1:m) = norm2(node_X(1:n-1,1:m,:)-node_X(2:n,1:m,:), dim=3)
  jface_normal(1:n-1,1:m,1) = ( node_X(1:n-1,1:m,2)-node_X(2:n,1:m,2))/jface_area(1:n-1,1:m)
  jface_normal(1:n-1,1:m,2) = (-node_X(1:n-1,1:m,1)+node_X(2:n,1:m,1))/jface_area(1:n-1,1:m)

end subroutine face_topo


!*************************************************!
!               JACOBIAN MATRIX                   !
!*************************************************!
subroutine jacobian_inverse(jacob, jacob_inv)
  implicit none
  real(8), intent(in) :: jacob(2,2)
  real(8), intent(out) :: jacob_inv(2,2)
  real(8) :: det
  det = jacob(1,1)*jacob(2,2)-jacob(1,2)*jacob(2,1)
  jacob_inv(1,:) = (/jacob(2,2), -jacob(1,2)/)/det
  jacob_inv(2,:) = (/-jacob(2,1), jacob(1,1)/)/det
end subroutine jacobian_inverse


subroutine interpolation(fl, dl, fr, dr, f)
  implicit none
  real(8), intent(in) :: fl, dl, fr, dr
  real(8), intent(out) :: f
  f = (fl*dr+fr*dl)/(dl+dr)
end subroutine interpolation


subroutine extrapolation(f1, d1, f2, d2, f)
  implicit none
  real(8), intent(in) :: f1, d1, f2, d2
  real(8), intent(out) :: f
  f = f1 +(f1-f2)*d1/d2
end subroutine extrapolation


subroutine cell_jacobian(node_X, cell_J, cell_J_inv, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: node_X(n,m,2)
  real(8), intent(out) :: cell_J(n-1,m-1,2,2), cell_J_inv(n-1,m-1,2,2)
  real(8) :: jacob(2,2), dksi(2), deta(2)
  integer :: i,j

  !cell jacobian
  do i=1,n-1
    do j=1,m-1
      dksi = node_X(i+1,j,:) - node_X(i,j,:) + node_X(i+1,j+1,:) - node_X(i,j+1,:)
      deta = node_X(i,j+1,:) - node_X(i,j,:) + node_X(i+1,j+1,:) - node_X(i+1,j,:)
      jacob(:,1) = dksi/norm2(dksi)
      jacob(:,2) = deta/norm2(deta)
      cell_J_inv(i,j,:,:) = jacob
      call jacobian_inverse(jacob, cell_J(i,j,:,:))
    end do
  end do
end subroutine cell_jacobian

subroutine face_jacobian(node_X, cell_X, iface_J, jface_J, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: node_X(n,m,2), cell_X(n-1,m-1,2)
  real(8), intent(out) :: iface_J(n,m-1,2,2), jface_J(n-1,m,2,2)
  real(8) :: jacob(2,2), dksi(2), deta(2), side_X(2), dXl(2), dXr(2), dl, dr
  integer :: i,j

  !inner iface jacobian
  do i=2,n-1
    do j=1,m-1
      deta = node_X(i,j+1,:) - node_X(i,j,:)
      side_X = 0.5*(node_X(i,j+1,:) + node_X(i,j,:))
      dXl = side_X - cell_X(i-1,j,:)
      dXr = cell_X(i,j,:) - side_X
      dl = norm2(dXl)
      dr = norm2(dXr)
      dXl = dXl/dl !unit vector
      dXr = dXr/dr !unit vector
      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      dksi = (dXl*dr + dXr*dl)/(dl+dr)
      jacob(1,:) = dksi
      jacob(2,:) = deta/norm2(deta)
      call jacobian_inverse(jacob, iface_J(i,j,:,:))
    end do
  end do
  !boundary iface jacobian
  !i=1
  do j=1,m-1
    dksi = cell_X(1,j,:)-0.5*(node_X(1,j+1,:) + node_X(1,j,:))
    deta = node_X(1,j+1,:) - node_X(1,j,:)
    jacob(1,:) = dksi/norm2(dksi)
    jacob(2,:) = deta/norm2(deta)
    call jacobian_inverse(jacob, iface_J(1,j,:,:))
  end do
  !i=n
  do j=1,m-1
    dksi = 0.5*(node_X(n,j+1,:) + node_X(n,j,:)) - cell_X(n-1,j,:)
    deta = node_X(n,j+1,:) - node_X(n,j,:)
    jacob(1,:) = dksi/norm2(dksi)
    jacob(2,:) = deta/norm2(deta)
    call jacobian_inverse(jacob, iface_J(n,j,:,:))
  end do

  !inner jface jacobian
  do i=1,n-1
    do j=2,m-1
      dksi = node_X(i+1,j,:) - node_X(i,j,:)
      side_X = 0.5*(node_X(i+1,j,:) + node_X(i,j,:))
      dXl = side_X - cell_X(i,j-1,:)
      dXr = cell_X(i,j,:) - side_X
      dl = norm2(dXl)
      dr = norm2(dXr)
      dXl = dXl/dl !unit vector
      dXr = dXr/dr !unit vector
      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      deta = (dXl*dr + dXr*dl)/(dl+dr)
      jacob(1,:) = dksi/norm2(dksi)
      jacob(2,:) = deta
      call jacobian_inverse(jacob, jface_J(i,j,:,:))
    end do
  end do

  !boundary jface jacobian
  !j=1
  do i=1,n-1
    deta = cell_X(i,1,:)-0.5*(node_X(i+1,1,:) + node_X(i,1,:))
    dksi = node_X(i+1,1,:) - node_X(i,1,:)
    jacob(1,:) = dksi/norm2(dksi)
    jacob(2,:) = deta/norm2(deta)
    call jacobian_inverse(jacob, jface_J(i,1,:,:))
  end do
  !j=m
  do i=1,n-1
    deta = 0.5*(node_X(i+1,m,:) + node_X(i,m,:)) - cell_X(i,m-1,:)
    dksi = node_X(i+1,m,:) - node_X(i,m,:)
    jacob(1,:) = dksi/norm2(dksi)
    jacob(2,:) = deta/norm2(deta)
    call jacobian_inverse(jacob, jface_J(i,m,:,:))
  end do
end subroutine face_jacobian


!*************************************************!
!              CONVERT VARIABLES                  !
! P=(rho, u, v, p), U=(rho, rho*u, rho*v, rho*e)  !
!*************************************************!
subroutine P2U(cell_P, cell_U, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: cell_P(n, m, 4)        
  real(8), intent(inout) :: cell_U(n, m, 4)
        
  cell_U(:,:,1) = cell_P(:,:,1)
  cell_U(:,:,2) = cell_P(:,:,1)*cell_P(:,:,2)
  cell_U(:,:,3) = cell_P(:,:,1)*cell_P(:,:,3)
  cell_U(:,:,4) = cell_P(:,:,4)*2.5 + 0.5*(cell_P(:,:,2)**2+cell_P(:,:,3)**2)
  !1/(gamm-1.0) = 2.5

end subroutine P2U 

subroutine U2P(cell_U, all_cell_P, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: cell_U(n, m, 4)        
  real(8), intent(inout) :: all_cell_P(n+2, m+2, 4)

  all_cell_P(2:n+1,2:m+1,1) = cell_U(:,:,1)
  all_cell_P(2:n+1,2:m+1,2) = cell_U(:,:,2)/cell_U(:,:,1)
  all_cell_P(2:n+1,2:m+1,3) = cell_U(:,:,3)/cell_U(:,:,1)
  all_cell_P(2:n+1,2:m+1,4) = (cell_U(:,:,4) - 0.5*(cell_U(:,:,2)**2+cell_U(:,:,3)**2)/cell_U(:,:,1))*0.4
  !(gamm-1.0) = 0.4

end subroutine U2P

subroutine  P2F(P, normal, area, F)
  implicit none
  real(8), intent(in) :: P(4), normal(2), area
  real(8), intent(out) :: F(4)
  real(8) :: Vn
  Vn = dot_product(normal, P(2:3))
  F(1) = P(1)*Vn
  F(2) = F(1)*P(2)+P(4)*normal(1)
  F(3) = F(1)*P(3)+P(4)*normal(2)
  F(4) = F(1)*(P(4)/P(1)*3.5 + 0.5*(P(2)*P(2)+P(3)*P(3))) !rho*vn*(p/rho*(gamma/gamma-1) + 0.5*V*V)
  F = F*area
end subroutine P2F


!***************************************!
!         BOUNDARY CONDITIONS           !
!***************************************!
subroutine bc_symmetry(face_normal, icell_P, side_P, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: face_normal(n,2), icell_P(n,4)
  real(8), intent(out) :: side_P(n,4)
  real(8) :: V_in(2), Vn, V_bd(2)
  integer :: i

  do i=1,n
    V_in = icell_P(i,2:3)
    Vn = dot_product(face_normal(i,:), V_in)
    V_bd = V_in - Vn*face_normal(i,:)
    side_P(i,:) = (/ icell_P(i,1), V_bd(1), V_bd(2), icell_P(i,4) /)
  end do
end subroutine bc_symmetry

subroutine bc_inflow(face_normal, icell_P, P_freestream, side_P, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: face_normal(n,2), icell_P(n,4), P_freestream(4)
  real(8), intent(out) :: side_P(n,4)
  real(8), parameter :: k = 1.4, km1 = 0.4, km1_inv = 2.5
  real(8) :: a_out, Vn_out, a_in, Vn_in, R_plus, R_minus, a_bound, Vn_bound
  real(8) :: R_bound, rho_bound, V_bound(2), p_bound
  integer :: i

  do i=1,n
    a_out = sqrt(k*P_freestream(4)/P_freestream(1))
    Vn_out = dot_product(face_normal(i,:), P_freestream(2:3))
    if (abs(Vn_out) > a_out) then
      side_P(i,:) = P_freestream(:)
    else
      a_in = sqrt(k*icell_P(i,4)/icell_P(i,1))
      Vn_in = dot_product(face_normal(i,:), icell_P(i,2:3))
      R_plus = Vn_out + 2.0*a_out*km1_inv
      R_minus = Vn_in - 2.0*a_in*km1_inv
      a_bound = 0.25*km1*(R_plus-R_minus)
      Vn_bound = 0.5*(R_plus+R_minus)
      R_bound = P_freestream(4)/(P_freestream(1)**k)
      rho_bound = (a_bound**2/(k*R_bound))**km1_inv
      V_bound = P_freestream(2:3) + (Vn_bound-Vn_out)*face_normal(i,:)
      p_bound = R_bound*(rho_bound**k)

      side_P(i,:) = (/ rho_bound, V_bound(1), V_bound(2), p_bound /)
    end if
  end do

end subroutine bc_inflow


subroutine bc_outflow(face_normal, icell_P, pressure_exit, opt, side_P, n)
  implicit none
  integer, intent(in) :: opt, n
  real(8), intent(in) :: face_normal(n,2), icell_P(n,4), pressure_exit
  real(8), intent(out) :: side_P(n,4)
  real(8), parameter :: k = 1.4, km1 = 0.4, km1_inv = 2.5
  real(8) :: a_in, Vn_in, R_plus, a_bound, Vn_bound, rho_bound, V_bound(2), p_bound = 0.0
  integer :: i

  if (opt==1) then
    p_bound = pressure_exit
  end if

  do i=1,n
    a_in = sqrt(k*icell_P(i,4)/icell_P(i,1))
    Vn_in = dot_product(face_normal(i,:), icell_P(i,2:3))
    if (abs(Vn_in) > a_in) then
      side_P(i,:) = icell_P(i,:)
    else
      if (opt==2) then
        p_bound = 0.5*(pressure_exit+icell_P(i,4))
      end if
      rho_bound = icell_P(i,1)*(p_bound/icell_P(i,4))**(1.0/k)
      a_bound = sqrt(k*p_bound/rho_bound)
      R_plus = Vn_in + 2*a_in*km1_inv
      Vn_bound = R_plus - 2*a_bound*km1_inv
      V_bound = icell_P(i,2:3) + (Vn_bound-Vn_in)*face_normal(i,:)

      side_P(i,:) = (/ rho_bound, V_bound(1), V_bound(2), p_bound /)
    end if
  end do

end subroutine bc_outflow


subroutine bound_flux_p2f(face_normal, face_area, side_P, cell_res, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: face_normal(n, 2), face_area(n), side_P(n,4)
  real(8), intent(out) :: cell_res(n, 4)
  integer :: i

  do i=1,n
    call P2F(side_P(i,:), face_normal(i,:), face_area(i), cell_res(i,:))
  end do
end subroutine bound_flux_p2f


!***************************************!
!          GLOBAL TIME STEP             !
!***************************************!
subroutine eu_time_step(cell_P, cell_half_size, dt_min, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) ::  cell_P(n, m, 4), cell_half_size(n, m, 2)
  real(8), intent(out) :: dt_min
  real(8) :: dt, v, a
  integer :: i, j

  dt_min = 1000.0
  do i=1,n
    do j=1,m
      v = norm2(cell_P(i,j,2:3)) !air speed
      a = sqrt(1.4*cell_P(i,j,4)/cell_P(i,j,1))
      dt = min(cell_half_size(i,j,1), cell_half_size(i,j,2))/(v+a)
      dt_min = min(dt_min, dt)
    end do
  end do
  dt_min = dt_min*2

end subroutine eu_time_step

subroutine ns_time_step(cell_P, cell_half_size, dt_min, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) ::  cell_P(n, m, 4), cell_half_size(n, m, 2)
  real(8), intent(out) :: dt_min
  real(8), parameter :: Pr = 0.72, Rgas=287.052873836, c_mu=1.506196284597745e-06
  real(8) :: dt, v, a, dd, T, mu, tau_conv(2), tau_diff(2), tau_inv(2)
  integer :: i, j

  dt_min = 1000.0
  do i=1,n
    do j=1,m
      v = norm2(cell_P(i,j,2:3)) !air speed
      a = sqrt(1.4*cell_P(i,j,4)/cell_P(i,j,1))
      tau_conv = 2*cell_half_size(i,j,:)/(a+v)

      T = cell_P(i,j,4)/(Rgas*cell_P(i,j,1))
      mu = c_mu*(T**1.5)/(T+122.0)
      dd = mu/(Pr*cell_P(i,j,1))

      tau_diff = (2*cell_half_size(i,j,:))**2/dd
      tau_inv = (1.0 + sqrt(1.0+tau_diff/tau_conv))**2/tau_diff
      dt = 1.0/sum(tau_inv)

      dt_min = min(dt_min, dt)
    end do
  end do

end subroutine ns_time_step

!***************************************!
!     CALCULATION - MISCLOSURE          !
!***************************************!
subroutine misclosure_rho(cell_rho_prev, cell_rho, sum_dr, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: cell_rho_prev(n, m), cell_rho(n, m)
  real(8), intent(out) :: sum_dr
  real(8) :: dr
  integer :: i, j

  sum_dr = 1.0e-16
  do i=1,n
    do j=1,m
      if (cell_rho(i,j) == cell_rho_prev(i,j)) then
        dr = 0.0
      else
        dr = abs(cell_rho(i,j)-cell_rho_prev(i,j))/min(cell_rho(i,j),cell_rho_prev(i,j))
      end if
      
      sum_dr = sum_dr + dr
    end do
  end do
  
end subroutine misclosure_rho
