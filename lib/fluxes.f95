! Copyright (C) 2021  Nguyen Ngoc Sang, <https://github.com/SangVn>

!*************************************************!
!         https://www.facebook.com/VnCFD          !
!        https://vncfdgroup.wordpress.com         !
!*************************************************!

!compiler: python3 -m numpy.f2py -c functions.f95 fluxes.f95 -m fortranlib

!*************************************************!
!                     MINMOD                      !
!*************************************************!
subroutine minmod_kolgan(a, b, c)
  implicit none
  real(8), intent(in) :: a, b
  real(8), intent(out) :: c

  if (a==0.0 .or. b==0.0) then
    c = 0.0
  else if (a>0.0 .and. b>0.0) then
    c = min(a,b)
  else if (a<0.0 .and. b<0.0) then
    c = max(a,b)
  else
    c = 0.0
  end if

end subroutine minmod_kolgan

subroutine minmod_van_leer(a, b, c)
  implicit none
  real(8), intent(in) :: a, b
  real(8), intent(out) :: c

  if (a==0.0 .or. b==0.0) then
    c = 0.0
  else if (a>0.0 .and. b>0.0) then
    c = min(1.25*min(a,b), 0.5*(a+b))
  else if (a<0.0 .and. b<0.0) then
    c = max(1.25*max(a,b), 0.5*(a+b))
  else
    c = 0.0
  end if
end subroutine minmod_van_leer


!*************************************************!
!              RECONSTRUCTION - VECTOR            !
!*************************************************!
subroutine reconstr_vector(Pl, P, Pr, dl, d, dr, dP)
  implicit none
  real(8), intent(in) :: Pl(4), P(4), Pr(4), dl, d, dr
  real(8), intent(out) :: dP(4)
  real(8) :: dPl(4), dPr(4), dV(2), V(2), Vl(2), Vr(2), Vmod, dVmod, ddl, ddr

  ddl = d/(dl+d)
  ddr = d/(d+dr)
  dPl = (P-Pl)*ddl
  dPr = (Pr-P)*ddr
  !call minmod_kolgan(dPl(1), dPr(1), dP(1))
  !call minmod_kolgan(dPl(4), dPr(4), dP(4))
  call minmod_van_leer(dPl(1), dPr(1), dP(1))
  call minmod_van_leer(dPl(4), dPr(4), dP(4))

  dV = 0.5*(dPl(2:3)+dPr(2:3))
  V  = P(2:3)
  Vl = V-dV
  Vr = V+dV

  Vmod = norm2(V)
  !call minmod_kolgan((Vmod-norm2(Pl(2:3)))*ddl, (norm2(Pr(2:3))-Vmod)*ddr, dVmod)
  call minmod_van_leer((Vmod-norm2(Pl(2:3)))*ddl, (norm2(Pr(2:3))-Vmod)*ddr, dVmod)

  dP(2:3) = 0.5*(Vr/norm2(Vr)*(Vmod+dVmod) - Vl/norm2(Vl)*(Vmod-dVmod))
    
end subroutine reconstr_vector

subroutine reconstruction_vector(all_cell_P, all_cell_half_size, cell_dP, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: all_cell_P(n+2,m+2,4), all_cell_half_size(n+2,m+2,2)
  real(8), intent(inout) :: cell_dP(2, n,m,4)
  integer :: i, j

  ! i-direction
  do i=2,n+1
    do j=2,m+1
      call reconstr_vector(all_cell_P(i-1,j,:), all_cell_P(i,j,:), all_cell_P(i+1,j,:), &
              all_cell_half_size(i-1,j,1), all_cell_half_size(i,j,1), all_cell_half_size(i+1,j,1), cell_dP(1,i-1,j-1,:))
    end do
  end do

  ! j-direction
  do i=2,n+1
    do j=2,m+1
      call reconstr_vector(all_cell_P(i,j-1,:), all_cell_P(i,j,:), all_cell_P(i,j+1,:), &
              all_cell_half_size(i,j-1,2), all_cell_half_size(i,j,2), all_cell_half_size(i,j+1,2), cell_dP(2,i-1,j-1,:))
    end do
  end do

end subroutine reconstruction_vector


!*************************************************!
!          RECONSTRUCTION - COMPONENT             !
!*************************************************!
subroutine reconstr_component(Pl, P, Pr, dl, d, dr, J, J_inv, dP)
  implicit none
  real(8), intent(in) :: Pl(4), P(4), Pr(4), dl, d, dr, J(2,2), J_inv(2,2)
  real(8), intent(out) :: dP(4)
  real(8) :: dPl(4), dPr(4), V_tran(2), ddl, ddr
  integer :: i

  ddl = d/(dl+d)
  ddr = d/(d+dr)
  dPl = (P-Pl)*ddl
  dPr = (Pr-P)*ddr
  V_tran = matmul(J, dPl(2:3))
  dPl(2:3) = V_tran(:)
  V_tran = matmul(J, dPr(2:3))
  dPr(2:3) = V_tran(:)

  do i=1,4
    call minmod_van_leer(dPl(i), dPr(i), dP(i))
  end do
  V_tran = matmul(J_inv, dP(2:3))
  dP(2:3) = V_tran(:)
end subroutine reconstr_component

subroutine reconstruction_component(all_cell_P, all_cell_half_size, cell_J, cell_J_inv, cell_dP, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: all_cell_P(n+2,m+2,4), all_cell_half_size(n+2,m+2,2), cell_J(n,m,2,2), cell_J_inv(n,m,2,2)
  real(8), intent(inout) :: cell_dP(2, n,m,4)
  integer :: i, j

  ! i-direction
  do i=2,n+1
    do j=2,m+1
      call reconstr_component(all_cell_P(i-1,j,:), all_cell_P(i,j,:), all_cell_P(i+1,j,:), &
              all_cell_half_size(i-1,j,1), all_cell_half_size(i,j,1), all_cell_half_size(i+1,j,1), &
              cell_J(i-1,j-1,:,:), cell_J_inv(i-1,j-1,:,:), cell_dP(1,i-1,j-1,:))
    end do
  end do

  ! j-direction
  do i=2,n+1
    do j=2,m+1
      call reconstr_component(all_cell_P(i,j-1,:), all_cell_P(i,j,:), all_cell_P(i,j+1,:), &
              all_cell_half_size(i,j-1,2), all_cell_half_size(i,j,2), all_cell_half_size(i,j+1,2), &
              cell_J(i-1,j-1,:,:), cell_J_inv(i-1,j-1,:,:), cell_dP(2,i-1,j-1,:))
    end do
  end do

end subroutine reconstruction_component


!*************************************************!
!                  GRADIENT - CELL                !
!*************************************************!
subroutine gradient_cell(all_cell_P, all_cell_half_size, cell_J, cell_J_inv, cell_G, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: all_cell_P(n+2,m+2,4), all_cell_half_size(n+2,m+2,2)
  real(8), intent(in) :: cell_J(n,m,2,2), cell_J_inv(n,m,2,2)
  real(8), intent(inout) :: cell_G(2, n,m,4)
  real(8) :: dl, dr, dksil(4), dksi(4), dksir(4), detal(4), deta(4), detar(4), V_tran(2)
  integer :: i, j

  do i=2,n+1
    do j=2,m+1
      dl = all_cell_half_size(i,j,1) + all_cell_half_size(i-1,j,1)
      dr = all_cell_half_size(i+1,j,1) + all_cell_half_size(i,j,1)
      dksil = all_cell_P(i,j,:) - all_cell_P(i-1,j,:)
      dksir = all_cell_P(i+1,j,:) - all_cell_P(i,j,:)

      V_tran = matmul(cell_J(i-1,j-1,:,:), dksil(2:3))
      dksil(2:3) = V_tran(:)
      V_tran = matmul(cell_J(i-1,j-1,:,:), dksir(2:3))
      dksir(2:3) = V_tran(:)

      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      dksi = (dksil/dl*dr + dksir/dr*dl)/(dl+dr)
      V_tran = matmul(cell_J_inv(i-1,j-1,:,:), dksi(2:3))
      dksi(2:3) = V_tran(:)

      dl = all_cell_half_size(i,j,2) + all_cell_half_size(i,j-1,2)
      dr = all_cell_half_size(i,j+1,2) + all_cell_half_size(i,j,2)
      detal = all_cell_P(i,j,:) - all_cell_P(i,j-1,:)
      detar = all_cell_P(i,j+1,:) - all_cell_P(i,j,:)

      V_tran = matmul(cell_J(i-1,j-1,:,:), detal(2:3))
      detal(2:3) = V_tran(:)
      V_tran = matmul(cell_J(i-1,j-1,:,:), detar(2:3))
      detar(2:3) = V_tran(:)

      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      deta = (detal/dl*dr + detar/dr*dl)/(dl+dr)
      V_tran = matmul(cell_J_inv(i-1,j-1,:,:), deta(2:3))
      deta(2:3) = V_tran(:)

      cell_G(1,i-1,j-1,:) = dksi
      cell_G(2,i-1,j-1,:) = deta
    end do
  end do

end subroutine gradient_cell

!*************************************************!
!                  GRADIENT - FACE                !
!*************************************************!
subroutine gradient_inner_face(lcell_P, lcell_half_size, lcell_G, rcell_P, rcell_half_size, rcell_G, face_G, nn, n)
  implicit none
  integer, intent(in) :: nn, n
  real(8), intent(in) :: lcell_P(n,4), lcell_half_size(n,2), lcell_G(2,n,4)
  real(8), intent(in) :: rcell_P(n,4), rcell_half_size(n,2), rcell_G(2,n,4)
  real(8), intent(out) :: face_G(2,n,4)
  real(8) :: dl, dr, d
  integer :: i, nt=2

  if (nn==2) then
    nt = 1
  end if

  do i=1,n
    dl = lcell_half_size(i,nn)
    dr = rcell_half_size(i,nn)
    d = 1.0/(dl+dr)
    ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
    face_G(nt, i,:) = (lcell_G(nt, i,:)*dr + rcell_G(nt ,i,:)*dl)*d
    face_G(nn, i,:) = (rcell_P(i,:) - lcell_P(i,:) - 0.5*(dr-dl)*(rcell_G(nn, i,:)-lcell_G(nn, i,:)))*d
  end do
end subroutine gradient_inner_face

subroutine gradient_face(all_cell_P, all_cell_half_size, cell_G, iface_G, jface_G, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: all_cell_P(n+2,m+2,4), all_cell_half_size(n+2,m+2,2), cell_G(2, n,m,4)
  real(8), intent(inout) :: iface_G(2 ,n+1,m,4), jface_G(2 ,n,m+1,4)
  real(8) :: dl, dr, d, d1, d2
  integer :: i, j

  !inner ifaceG
  do i=2,n
    do j=1,m
      dl = all_cell_half_size(i,j+1,1)
      dr = all_cell_half_size(i+1,j+1,1)
      d = 1.0/(dl+dr)
      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      iface_G(2, i,j,:) = (cell_G(2, i-1,j,:)*dr+cell_G(2 ,i,j,:)*dl)*d
      iface_G(1, i,j,:) = (all_cell_P(i+1,j+1,:) - all_cell_P(i,j+1,:) - 0.5*(dr-dl)*(cell_G(1, i,j,:)-cell_G(1, i-1,j,:)))*d
    end do
  end do

  !bound ifaceG: i=1,n+1
  do j=1,m
    d1 = all_cell_half_size(2,j+1,1)
    d2 = all_cell_half_size(3,j+1,1)
    ! extrapolation: f = f1 +(f1-f2)*d1/d2
    iface_G(2, 1,j,:) = cell_G(2, 1,j,:) + (cell_G(2, 1,j,:) - cell_G(2 ,2,j,:))*d1/d2
    iface_G(1, 1,j,:) = (all_cell_P(2,j+1,:) - all_cell_P(1,j+1,:))/(d1 + all_cell_half_size(1,j+1,1))
  end do
  do j=1,m
    d1 = all_cell_half_size(n+1,j+1,1)
    d2 = all_cell_half_size(n,j+1,1)
    ! extrapolation: f = f1 +(f1-f2)*d1/d2
    iface_G(2, n+1,j,:) = cell_G(2, n,j,:) + (cell_G(2, n,j,:) - cell_G(2 ,n-1,j,:))*d1/d2
    iface_G(1, n+1,j,:) = (all_cell_P(n+2,j+1,:) - all_cell_P(n+1,j+1,:))/(d1 + all_cell_half_size(n+2,j+1,1))
  end do

  !inner jfaceG
  do i=1,n
    do j=2,m
      dl = all_cell_half_size(i+1,j,2)
      dr = all_cell_half_size(i+1,j+1,2)
      d = 1.0/(dl+dr)
      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      jface_G(1, i,j,:) = (cell_G(1, i,j-1,:)*dr+cell_G(1 ,i,j,:)*dl)*d
      jface_G(2, i,j,:) = (all_cell_P(i+1,j+1,:) - all_cell_P(i+1,j,:) - 0.5*(dr-dl)*(cell_G(2, i,j,:)-cell_G(2, i,j-1,:)))*d
    end do
  end do

  !bound jfaceG: j=1,m+1
  do i=1,n
    d1 = all_cell_half_size(i+1,2,2)
    d2 = all_cell_half_size(i+1,3,2)
    ! extrapolation: f = f1 +(f1-f2)*d1/d2
    jface_G(1, i,1,:) = cell_G(1, i,1,:) + (cell_G(1, i,1,:) - cell_G(1 ,i,2,:))*d1/d2
    jface_G(2, i,1,:) = (all_cell_P(i+1,2,:) - all_cell_P(i+1,1,:))/(d1 + all_cell_half_size(i+1,1,2))
  end do
  do i=1,n
    d1 = all_cell_half_size(i+1,m+1,2)
    d2 = all_cell_half_size(i+1,m,2)
    ! extrapolation: f = f1 +(f1-f2)*d1/d2
    jface_G(1, i,m+1,:) = cell_G(1, i,m,:) + (cell_G(1, i,m,:) - cell_G(1 ,i,m-1,:))*d1/d2
    jface_G(2, i,m+1,:) = (all_cell_P(i+1,m+2,:) - all_cell_P(i+1,m+1,:))/(d1 + all_cell_half_size(i+1,m+2,2))
  end do
end subroutine gradient_face


!*************************************************!
!                 CALCULATION FACE_P              !
!*************************************************!
subroutine calc_face_P(cell_P, cell_half_size, iface_P, jface_P, n,m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: cell_P(n,m,4), cell_half_size(n,m,2)
  real(8), intent(inout) :: iface_P(n+1,m,4), jface_P(n,m+1,4)
  real(8) :: dl, dr
  integer :: i, j

  !inner iface_P
  do i=2,n
    do j=1,m
      dl = cell_half_size(i-1,j,1)
      dr = cell_half_size(i,j,1)
      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      iface_P(i,j,:) = (cell_P(i-1,j,:)*dr+cell_P(i,j,:)*dl)/(dl+dr)
    end do
  end do

  !inner jface_P
  do i=1,n
    do j=2,m
      dl = cell_half_size(i,j-1,2)
      dr = cell_half_size(i,j,2)
      ! interpolation: f = (fl*dr+fr*dl)/(dl+dr)
      jface_P(i,j,:) = (cell_P(i,j-1,:)*dr+cell_P(i,j,:)*dl)/(dl+dr)
    end do
  end do

end subroutine calc_face_P


!*************************************************!
!                  FLUX - CONVECTION              !
!*************************************************!

!*************************************************!
!       --------+---jf(i,j+1)--+-----------       !
!               |              |                  !
!        if(i,j)|    c(i*j)    |if(i+1,j)         !
!               |              |                  !
!       --------+----jf(i,j)---+-----------       !
!*************************************************!

!first order - godunov reconstruction
subroutine inner_flux_reconstr_godunov(iface_normal, iface_area, jface_normal, jface_area, cell_P, cell_res, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: iface_normal(n, m-1, 2), iface_area(n, m-1)
  real(8), intent(in) :: jface_normal(n-1, m, 2), jface_area(n-1, m)
  real(8), intent(in) :: cell_P(n-1, m-1, 4)
  real(8), intent(inout) :: cell_res(n-1, m-1, 4)
  real(8) :: flux(4)
  integer :: i, j

  !inner_iface_flux
  do i=2,n-1
    do j=1,m-1
      call flux_roe(cell_P(i-1,j,:), cell_P(i,j,:), iface_normal(i,j,:), iface_area(i,j), flux)
      cell_res(i-1,j,:) = cell_res(i-1,j,:)-flux
      cell_res(i,j,:) = cell_res(i,j,:)+flux
    end do
  end do

  !inner_jface_flux
  do i=1,n-1
    do j=2,m-1
      call flux_roe(cell_P(i,j-1,:), cell_P(i,j,:), jface_normal(i,j,:), jface_area(i,j), flux)
      cell_res(i,j-1,:) = cell_res(i,j-1,:)-flux
      cell_res(i,j,:) = cell_res(i,j,:)+flux
    end do
  end do

end subroutine inner_flux_reconstr_godunov

!second order - TVD reconstruction
subroutine inner_flux_reconstr_tvd(iface_normal, iface_area, jface_normal, jface_area, cell_P, cell_dP, cell_res, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: iface_normal(n, m-1, 2), iface_area(n, m-1)
  real(8), intent(in) :: jface_normal(n-1, m, 2), jface_area(n-1, m)
  real(8), intent(in) :: cell_P(n-1, m-1, 4), cell_dP(2, n-1, m-1, 4)
  real(8), intent(inout) :: cell_res(n-1, m-1, 4)
  real(8) :: flux(4), Pl(4), Pr(4)
  integer :: i, j

  !inner_iface_flux
  do i=2,n-1
    do j=1,m-1
      Pl = cell_P(i-1,j,:)+cell_dP(1,i-1,j,:)
      Pr = cell_P(i,j,:)-cell_dP(1,i,j,:)
      call flux_roe(Pl, Pr, iface_normal(i,j,:), iface_area(i,j), flux)
      cell_res(i-1,j,:) = cell_res(i-1,j,:)-flux
      cell_res(i,j,:) = cell_res(i,j,:)+flux
    end do
  end do

  !inner_jface_flux
  do i=1,n-1
    do j=2,m-1
      Pl = cell_P(i,j-1,:)+cell_dP(2,i,j-1,:)
      Pr = cell_P(i,j,:)-cell_dP(2,i,j,:)
      call flux_roe(Pl, Pr, jface_normal(i,j,:), jface_area(i,j), flux)
      cell_res(i,j-1,:) = cell_res(i,j-1,:)-flux
      cell_res(i,j,:) = cell_res(i,j,:)+flux
    end do
  end do

end subroutine inner_flux_reconstr_tvd

!boundary flux
subroutine bound_flux_roe(face_normal, face_area, cell_Pl, cell_Pr, cell_res, n)
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: face_normal(n, 2), face_area(n)
  real(8), intent(in) :: cell_Pl(n, 4), cell_Pr(n, 4)
  real(8), intent(out) :: cell_res(n, 4)
  integer :: i

  do i=1,n
    call flux_roe(cell_Pl(i,:), cell_Pr(i,:), face_normal(i,:), face_area(i), cell_res(i,:))
  end do

end subroutine bound_flux_roe



!*************************************************!
!                  FLUX - DIFFUSION               !
!*************************************************!
subroutine viscous_strees_tensor(P, G, mu, Tau, gradT)
  implicit none
  real(8), intent(in) :: P(4), G(2,4), mu
  real(8), intent(inout) :: Tau(2,2), gradT(2)
  real(8) :: mu_m2, divV_d3, temp, Rgas=287.052873836

  mu_m2 = -mu*2.0
  divV_d3 = (G(1,2)+G(2,3))/3.0
  Tau(1,2) = -mu*(G(1,3)+G(2,2))
  Tau(2,1) = Tau(1,2)
  Tau(1,1) = mu_m2*(G(1,2) - divV_d3)
  Tau(2,2) = mu_m2*(G(2,3) - divV_d3)

  temp = 1.0/(Rgas*P(1))
  gradT(1) = (G(1,4)-G(1,1)*P(4)/P(1))*temp
  gradT(2) = (G(2,4)-G(2,1)*P(4)/P(1))*temp

end subroutine viscous_strees_tensor

subroutine p2f_diff(face_normal, face_area, face_J, face_G, face_P, flux)
  implicit none
  real(8), intent(in) :: face_normal(2), face_area, face_J(2,2), face_G(2,4), face_P(4)
  real(8), intent(inout) :: flux(4)
  real(8):: mu, Tau(2,2), gradT(2), T, q(2), heat_flux, Gxy(2,4)
  real(8), parameter :: Rgas=287.052873836, c_mu=1.506196284597745e-06, cp_dPr=1395.3959144805553

  T = face_P(4)/(Rgas*face_P(1))
  mu = c_mu*(T**1.5)/(T+122.0)
  Gxy = matmul(transpose(face_J),face_G)
  call viscous_strees_tensor(face_P, Gxy, mu, Tau, gradT)
  q = cp_dPr*mu*gradT
  heat_flux = dot_product(face_normal,q)
  flux(1) = 0.0
  flux(2) = dot_product(face_normal, tau(1,:))
  flux(3) = dot_product(face_normal, tau(2,:))
  flux(4) = flux(2)*face_P(2) + flux(3)*face_P(3) - heat_flux
  flux = flux*face_area

end subroutine p2f_diff

subroutine flux_diff(iface_normal, iface_area, iface_J, iface_G, iface_P, &
        jface_normal, jface_area, jface_J, jface_G, jface_P, cell_res, n, m)
  implicit none
  integer, intent(in) :: n, m
  real(8), intent(in) :: iface_normal(n+1,m,2), iface_area(n+1,m), iface_J(n+1,m,2,2), iface_G(2,n+1,m,4), iface_P(n+1,m,4)
  real(8), intent(in) :: jface_normal(n,m+1,2), jface_area(n,m+1), jface_J(n,m+1,2,2), jface_G(2,n,m+1,4), jface_P(n,m+1,4)
  real(8), intent(inout) :: cell_res(n,m,4)
  real(8) :: mu, Tau(2,2), gradT(2), T, q(2), heat_flux, flux(4), Gxy(2,4)
  real(8), parameter :: Rgas=287.052873836, c_mu=1.506196284597745e-06, cp_dPr=1395.3959144805553
  integer :: i, j
  !c_mu = 1.72e-5*(273+122.0)/(273.0**1.5)
  !mu = c_mu*(T**1.5)/(T+122.)
  !cp_dPr = cp/Pr

  ! iface
  do i=1,n+1
    do j=1,m
      T = iface_P(i,j,4)/(Rgas*iface_P(i,j,1))
      mu = c_mu*(T**1.5)/(T+122.0)
      Gxy = matmul(iface_J(i,j,:,:),iface_G(:,i,j,:))
      call viscous_strees_tensor(iface_P(i,j,:), Gxy, mu, Tau, gradT)
      q = cp_dPr*mu*gradT
      heat_flux = dot_product(iface_normal(i,j,:),q)
      flux(1) = 0.0
      flux(2) = dot_product(iface_normal(i,j,:), tau(1,:))
      flux(3) = dot_product(iface_normal(i,j,:), tau(2,:))
      flux(4) = flux(2)*iface_P(i,j,2) + flux(3)*iface_P(i,j,3) - heat_flux
      flux = flux*iface_area(i,j)

      if (i<n+1) then
        cell_res(i,j,:) = cell_res(i,j,:)+flux
      end if
      if (i>1) then
        cell_res(i-1,j,:) = cell_res(i-1,j,:)-flux
      end if

    end do
  end do

  ! jface
  do i=1,n
    do j=1,m+1
      T = jface_P(i,j,4)/(Rgas*jface_P(i,j,1))
      mu = c_mu*(T**1.5)/(T+122.0)
      Gxy = matmul(jface_J(i,j,:,:),jface_G(:,i,j,:))
      call viscous_strees_tensor(jface_P(i,j,:), Gxy, mu, Tau, gradT)
      q = cp_dPr*mu*gradT
      heat_flux = dot_product(jface_normal(i,j,:),q)
      flux(1) = 0.0
      flux(2) = dot_product(jface_normal(i,j,:), tau(1,:))
      flux(3) = dot_product(jface_normal(i,j,:), tau(2,:))
      flux(4) = flux(2)*jface_P(i,j,2) + flux(3)*jface_P(i,j,3) - heat_flux
      flux = flux*jface_area(i,j)

      if (j<m+1) then
        cell_res(i,j,:) = cell_res(i,j,:)+flux
      end if
      if (j>1) then
        cell_res(i,j-1,:) = cell_res(i,j-1,:)-flux
      end if

    end do
  end do
end subroutine flux_diff


!*************************************************!
!                   FLUX - ROE                    !
!   Katate Masatsuka (http://www.cfdbooks.com)    !
!   Nguyen Ngoc Sang (https://github.com/sangvn)  !
!*************************************************!
subroutine flux_roe(primL, primR, njk, area, flux)
  implicit none
  real(8), parameter :: zero = 0.0, one = 1.0, two = 2.0, half = 0.5, fifth = 1.0/ 5.0, gamma = 1.4, gamma_m1 = 0.4
 !Input
  real(8), intent(in) :: primL(4), primR(4) ! [rho, u, v, p]_{L,R}
  real(8), intent(in) :: njk(2)             ! Face normal, njk=[nx, ny]
  real(8), intent(in) :: area
 !Output
  real(8), intent(out) :: flux(4)
 !Local variables
  real(8) :: nx, ny                  ! Normal vector
  real(8) :: mx, my                  ! Tangent vector: mx*nx+my*ny = 0
  real(8) :: uL, uR, vL, vR          ! Velocity components.
  real(8) :: rhoL, rhoR, pL, pR      ! Primitive variables.
  real(8) :: unL, unR, umL, umR      ! Normal and tangent velocities
  real(8) :: aL, aR, HL, HR          ! Speeds of sound.
  real(8) :: RT,rho,u,v,H,a,un, um   ! Roe-averages
  real(8) :: drho,dun,dum,dp,LdU(4)  ! Wave strenghs
  real(8) :: ws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
  real(8) :: fL(4), fR(4), diss(4)   ! Fluxes ad dissipation term
  real(8) :: dws(4)                  ! User-specified width for entropy fix
  !real(8) :: wsn
  integer :: i, j

  nx = njk(1)
  ny = njk(2)

!Tangent vector (Do you like it? Actually, Roe flux can be implemented
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  mx = -ny
  my =  nx

!Primitive and other variables.
!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
     unL = uL*nx+vL*ny
     umL = uL*mx+vL*my
      pL = primL(4)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma_m1) + half*(uL*uL+vL*vL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
     unR = uR*nx+vR*ny
     umR = uR*mx+vR*my
      pR = primR(4)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma_m1) + half*(uR*uR+vR*vR)

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
     u = (uL+RT*uR)/(one+RT)
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT* HR)/(one+RT)
     a = sqrt( (gamma_m1)*(H-half*(u*u+v*v)) )
    un = u*nx+v*ny
    um = u*mx+v*my

!Wave Strengths
   drho = rhoR - rhoL
     dp =   pR - pL
    dun =  unR - unL
    dum =  umR - umL

  LdU(1) = (dp - rho*a*dun )/(two*a*a)
  LdU(2) = rho*dum
  LdU(3) =  drho - dp/(a*a)
  LdU(4) = (dp + rho*a*dun )/(two*a*a)

!Wave Speed
  ws(1) = abs(un-a)
  ws(2) = abs(un)
  ws(3) = abs(un)
  ws(4) = abs(un+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one
  Rv(2,1) = u - a*nx
  Rv(3,1) = v - a*ny
  Rv(4,1) = H - un*a

  Rv(1,2) = zero
  Rv(2,2) = mx
  Rv(3,2) = my
  Rv(4,2) = um

  Rv(1,3) = one
  Rv(2,3) = u
  Rv(3,3) = v
  Rv(4,3) = half*(u*u+v*v)

  Rv(1,4) = one
  Rv(2,4) = u + a*nx
  Rv(3,4) = v + a*ny
  Rv(4,4) = H + un*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*unL
  fL(2) = rhoL*unL * uL + pL*nx
  fL(3) = rhoL*unL * vL + pL*ny
  fL(4) = rhoL*unL * HL

  fR(1) = rhoR*unR
  fR(2) = rhoR*unR * uR + pR*nx
  fR(3) = rhoR*unR * vR + pR*ny
  fR(4) = rhoR*unR * HR

  flux = half * (fL + fR - diss) *area
  !wsn = half*(abs(un) + a)  !Normal max wave speed times half

end subroutine flux_roe
