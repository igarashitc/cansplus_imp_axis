module  model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains


  subroutine  model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
                         ,roi,pri,vxi,vyi,vzi,bxi,byi,bzi,phii &
                         ,x,dx,xm,y,dy,ym,z,dz,zm &
                         ,gx,gz,eta,min_dx)

!---Input & Output
  real(8),dimension(ix),      intent(out) :: x,dx
  real(8),dimension(0:ix),    intent(out) :: xm
  real(8),dimension(jx),      intent(out) :: y,dy
  real(8),dimension(0:jx),    intent(out) :: ym
  real(8),dimension(kx),      intent(out) :: z,dz
  real(8),dimension(0:kx),    intent(out) :: zm
  real(8),dimension(ix,jx,kx),intent(out) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(out) :: bx,by,bz


  real(8),dimension(ix,jx,kx),intent(out) :: eta,phi
  real(8),dimension(ix,jx,kx),intent(out) :: gx,gz
  real(8),dimension(ix,jx,kx),intent(out) :: roi,pri
  real(8),dimension(ix,jx,kx),intent(out) :: vxi,vyi,vzi
  real(8),dimension(ix,jx,kx),intent(out) :: bxi,byi,bzi
  real(8),dimension(ix,jx,kx),intent(out) :: phii
  real(8),                    intent(out) :: min_dx

  integer :: i,j,k,kzero
  integer :: ig,jg,kg
  real(8) :: ss
  real(8) :: tmp 
  real(8) :: psi0, pot0 
  real(8) :: roc,prc,vyc
  real(8) :: rod,prd,vyd,bxd,byd,bzd
  real(8) :: bbAbsMax
  real(8),dimension(ix,jx,kx) :: curx,cury,curz
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(kgx) :: zg,dzg
  real(8),dimension(0:kgx) :: zmg
  real(8),dimension(ix,jx,kx) :: gpot
  real(8),dimension(igx,jgx,kgx) :: gxg,gzg,gpotg

  real(8),dimension(ix,jx,kx) :: pr_per, flag_torus
  real(8),parameter :: x0=1.d0
  real(8) :: Lnrml,L_r

!---Step 1a.-------------------------------------------------------------|
! set global x-grid 
!$OMP PARALLEL
!$OMP DO
  do i=1,igx
     dxg(i) = dxg0
  enddo
!$OMP END DO
!$OMP END PARALLEL

  dxg(margin+1)=4.0d0*dxg0
  xmg(margin+1) = xmin + dxg(margin+1)
  do i=margin+2,igx
    xmg(i) = xmg(i-1)+dxg(i)
    if(xmg(i) >= ugrid_x .and. i < igx)then
      dxg(i+1) = min(dxg(i)*ratio_x,dxmax)
    endif
  enddo
  do i=margin,0,-1
     xmg(i) = xmg(i+1)-dxg(i+1)
  enddo

!$OMP PARALLEL DO
  do i=1,igx
     xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
  enddo
!$OMP END PARALLEL DO

!---Step 1b.-------------------------------------------------------------|
! set global y-grid

!$OMP PARALLEL DO
  do j=1,jgx
     dyg(j) = dyg0
  enddo
!$OMP END PARALLEL DO

! Y origin
  ymg(margin+1) = ymin+0.5d0*dyg(margin+1)
  do j=margin+1,jgx-1
     ymg(j+1) = ymg(j)+dyg(j+1)
  enddo
  do j=margin,0,-1
     ymg(j) = ymg(j+1)-dyg(j+1)
  enddo

!$OMP PARALLEL DO
  do j=1,jgx
     yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
  enddo
!$OMP END PARALLEL DO

!---Step 1c.-------------------------------------------------------------|
! set global z-grid
!$OMP PARALLEL DO
  do,k=1,kgx
     dzg(k) = dzg0
  enddo
!$OMP END PARALLEL DO

  kzero = int(kgx/2.0)+1
  zmg(kzero) = zmin+dzg(kzero)*(kzero-kgx/2.0)
  do k=kzero+1,kgx
     zmg(k) = zmg(k-1)+dzg(k)
     if(zmg(k) >= ugrid_z .and. k < kgx)then
       dzg(k+1) = min(dzg(k)*ratio_z,dzmax)
     endif
  enddo

  kzero = int(kgx/2.0+0.5)-1
  zmg(kzero) = zmin-dzg(kzero+1)*(kgx/2.0-kzero)
  do k=kzero-1,0,-1
     zmg(k) = zmg(k+1)-dzg(k+1)
     if(abs(zmg(k)) >= ugrid_z .and. k > 1)then
       dzg(k) = min(dzg(k+1)*ratio_z,dzmax)
     endif
  enddo

!$OMP PARALLEL DO
  do k=1,kgx
     zg(k) = 0.5d0*(zmg(k)+zmg(k-1))
  enddo
!$OMP END PARALLEL DO

!---Step 2a.-------------------------------------------------------------|
! set individual x-grid 
!$OMP PARALLEL
!$OMP DO PRIVATE(ig)
  do i=1,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
!$OMP END DO
!$OMP DO PRIVATE(ig)
  do i=0,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     xm(i) = xmg(ig)
  end do
!$OMP END DO

!---Step 2b.-------------------------------------------------------------|
! set individual y-grid 
!$OMP DO PRIVATE(jg)
  do j=1,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
!$OMP END DO
!$OMP DO PRIVATE(jg)
  do j=0,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  end do
!$OMP END DO

!---Step 2c.-------------------------------------------------------------|
! set individual z-grid 
!$OMP DO PRIVATE(kg)
  do k=1,kx
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
  enddo
!$OMP END DO
!$OMP DO PRIVATE(kg)
  do k=0,kx
     kg=mpid%mpirank_3d(3)*(kx-2*margin)+k
     zm(k) = zmg(kg)
  end do
!$OMP END DO
!$OMP END PARALLEL

! calculate min_dx
  min_dx = min(minval(dxg),minval(dzg))
!$OMP PARALLEL DO REDUCTION(min:min_dx)
  do j=1,jgx
     do i=margin+1,igx
        min_dx=min(min_dx,xg(i)*dyg(j))
     enddo
  enddo
!$OMP END PARALLEL DO
!----------------------------------------------------------------------|
! set initial model

     do k=1,kx
        do j=1,jx
           do i=1,ix
              if(y(j).lt.pi)then
                ro(i,j,k) = 1d0
                pr(i,j,k) = 1.d0
                vx(i,j,k) = 0.d0
                vy(i,j,k) = 0.d0
                vz(i,j,k) = 0.d0
                bx(i,j,k) = 0.d0
                by(i,j,k) = 0.d0
                bz(i,j,k) = 0.d0
                phi(i,j,k) = 0.d0
                eta(i,j,k) = 0.d0
              else
                ro(i,j,k) = 0.125d0
                pr(i,j,k) = 0.1d0
                vx(i,j,k) = 0.d0
                vy(i,j,k) = 0.d0
                vz(i,j,k) = 0.d0
                bx(i,j,k) = 0.d0
                by(i,j,k) = 0d0
                bz(i,j,k) = 0.d0
                phi(i,j,k) = 0.d0
                eta(i,j,k) = 0.d0
              endif
           enddo
        enddo
     enddo

end subroutine model_setup

end module model
