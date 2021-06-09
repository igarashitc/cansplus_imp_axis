module getEta

  implicit none
  private

  public :: getEta__anomalous


contains


  subroutine getEta__anomalous(ix,jx,kx,ro,bx,by,bz,x,dx,dy,dz,eta0,vc,eta,curx,cury,curz)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: eta0,vc
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(in) :: ro
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(out) :: eta,curx,cury,curz

  integer :: i,j,k
  real(8) :: cur_abs,vd,etamax,flag


! In Cylindrical coordinate
  call getcurrent__cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
                      ,curx,cury,curz)
  
  etamax=0.01d0
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j,cur_abs,vd,flag)
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           cur_abs = sqrt(+curx(i,j,k)*curx(i,j,k)+cury(i,j,k)*cury(i,j,k) &
                          +curz(i,j,k)*curz(i,j,k))
           vd = cur_abs/ro(i,j,k)

           flag = max(sign(1.0d0,ro(i,j,k)-1.0d-2),0.0d0)
           eta(i,j,k) = flag*min(etamax,eta0*(max((vd/vc-1.0d0),0.0d0)**2))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  end subroutine getEta__anomalous


  subroutine getcurrent__cyl(bx,by,bz,ix,jx,kx,x,dx,dy,dz &
                            ,curx,cury,curz)

  integer,intent(in) :: ix,jx,kx
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(ix,jx,kx),intent(in) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(out) :: curx,cury,curz

  integer :: i,j,k
  real(8) :: ddx,ddy,ddz
  real(8) :: line1,line2

! x-component
  !$OMP PARALLEL PRIVATE(i,j,ddx,ddy,ddz,line1,line2)
  !$OMP DO
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ddy = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)
           ddz = 0.5d0*dz(k-1)+dz(k)+0.5d0*dz(k+1)

           curx(i,j,k) = (bz(i,j+1,k)-bz(i,j-1,k))/(x(i)*ddy) &
                - (by(i,j,k+1)-by(i,j,k-1))/ddz
        enddo
     enddo
  enddo
  !$OMP END DO

! y-component
  !$OMP DO
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ddz = 0.5d0*dz(k-1)+dz(k)+0.5d0*dz(k+1)
           ddx = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
           
           cury(i,j,k) = (bx(i,j,k+1)-bx(i,j,k-1))/ddz &
                - (bz(i+1,j,k)-bz(i-1,j,k))/ddx
        enddo
     enddo
  enddo
  !$OMP END DO

! z-component
  !$OMP DO
  do k=2,kx-1
     do j=2,jx-1
        do i=2,ix-1
           ddx = 0.5d0*dx(i-1)+dx(i)+0.5d0*dx(i+1)
           ddy = 0.5d0*dy(j-1)+dy(j)+0.5d0*dy(j+1)

           line2 = x(i+1)
           line1 = x(i-1)
           curz(i,j,k) = (line2*by(i+1,j,k)-line1*by(i-1,j,k))/(x(i)*ddx) &
                - (bx(i,j+1,k)-bx(i,j-1,k))/(x(i)*ddy)
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  end subroutine getcurrent__cyl


end module getEta
