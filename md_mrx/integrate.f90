module integrate

  implicit none
  private

  public :: integrate__TVDRK3


contains


  subroutine integrate__TVDRK3(margin,ix,jx,kx,gm,dx,dy,dz,dt          &
                              ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cp &
                              ,eta,ccx,ccy,ccz)

  use convert
  use lr_state, only : lr_state__MP5
  use flux_calc
  use bnd
  use getEta

  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: ch,cp
  real(8),intent(in) :: dt,gm
  real(8),dimension(ix),intent(in) :: dx
  real(8),dimension(jx),intent(in) :: dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: phi
  real(8),dimension(ix,jx,kx),intent(inout) :: eta

!-- using flux
  real(8),dimension(ix,jx,kx) :: ro1,pr1,vx1,vy1,vz1
  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1
  real(8),dimension(ix,jx,kx) :: phi1
!-conserved variable
  real(8),dimension(ix,jx,kx) :: rx,ry,rz,ee
  real(8),dimension(ix,jx,kx) :: rx1,ry1,rz1,ee1
!-surface variables
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw
  real(8),dimension(ix,jx,kx) :: bx_m,by_m,bz_m,phi_m
!-Numerical flux
!x-component
  real(8),dimension(ix,jx,kx) :: frox,feex,frxx,fryx,frzx
  real(8),dimension(ix,jx,kx) :: fbyx,fbzx,fbxx,fphix
  real(8),dimension(ix,jx,kx) :: fbyxr,fbzxr,fbxxr,feexr
!y-component
  real(8),dimension(ix,jx,kx) :: froy,feey,frxy,fryy,frzy
  real(8),dimension(ix,jx,kx) :: fbyy,fbzy,fbxy,fphiy
  real(8),dimension(ix,jx,kx) :: fbyyr,fbzyr,fbxyr,feeyr
!z-component
  real(8),dimension(ix,jx,kx) :: froz,feez,frxz,fryz,frzz
  real(8),dimension(ix,jx,kx) :: fbxz,fbyz,fbzz,fphiz
  real(8),dimension(ix,jx,kx) :: fbxzr,fbyzr,fbzzr,feezr

!-other temporary variables
  integer :: mdir
  integer :: i,j,k,n
  real(8), parameter :: fac=1.D0/12.D0
  real(8) :: dtodx,dtody,dtodz,k1,k2
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

!-----Step 0.----------------------------------------------------------|
! primitive to conserve
  call convert__ptoc(ix,jx,kx,gm,ro,pr,vx,vy,vz,bx,by,bz &
                    ,rx,ry,rz,ee)

  ro1=ro
  rx1=rx
  ry1=ry
  rz1=rz
  bx1=bx
  by1=by
  bz1=bz
  ee1=ee
  phi1=phi

  do n=1,3
  
  !call getEta__anomalous(ix,jx,kx,ro,bx,by,bz,dx,dy,dz,eta0,vc,eta,curx,cury,curz)
  call getcurrent__cart(bx,by,bz,ix,jx,kx,dx,dy,dz &
                             ,curx,cury,curz)

!-----Step 1a.---------------------------------------------------------|
! Compute flux in x-direction
! set L/R state at x-direction
!
  mdir = 1

  call lr_state__MP5(mdir,ix,jx,kx,ro,pr &
       ,vx,vy,vz,bx,by,bz,phi &
       ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy,ccz)

  call flux_calc__bp(ix,jx,kx,bxw,phiw &
       ,bx_m,phi_m,ch)
  call flux_calc__glm(bx_m,phi_m,ch,fbxx,fphix,ix,jx,kx)
  
  call flux_calc__hlld(row,prw,vxw,vyw,vzw,bx_m,byw,bzw,gm,margin,ix,jx,kx &
       ,frox,feex,frxx,fryx,frzx,fbyx,fbzx)

  call flux_calc__fbres(mdir,margin,ix,jx,kx,fbyx,curz,eta,-1.0d0 &
       ,fbyxr)
  call flux_calc__fbres(mdir,margin,ix,jx,kx,fbzx,cury,eta,+1.0d0 &
       ,fbzxr)
  call flux_calc__feres(mdir,margin,ix,jx,kx,feex,curx,cury,curz,bx,by,bz,eta &
       ,feexr)
!-----Step 1b.---------------------------------------------------------|
! compute flux at y-direction
! set L/R state at y-direction
!
  mdir = 2

  call lr_state__MP5(mdir,ix,jx,kx,ro,pr &
       ,vy,vz,vx,by,bz,bx,phi &
       ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,ccx,ccy,ccz)

  call flux_calc__bp(ix,jx,kx,byw,phiw &
       ,by_m,phi_m,ch)

  call flux_calc__glm(by_m,phi_m,ch,fbyy,fphiy,ix,jx,kx)

  call flux_calc__hlld(row,prw,vyw,vzw,vxw,by_m,bzw,bxw,gm,margin,ix,jx,kx &
       ,froy,feey,fryy,frzy,frxy,fbzy,fbxy)

  call flux_calc__fbres(mdir,margin,ix,jx,kx,fbzy,curx,eta,-1.0d0 &
       ,fbzyr)
  call flux_calc__fbres(mdir,margin,ix,jx,kx,fbxy,curz,eta,+1.0d0 &
       ,fbxyr)
  call flux_calc__feres(mdir,margin,ix,jx,kx,feey,curx,cury,curz,bx,by,bz,eta &
       ,feeyr)
!-----Step 1c.---------------------------------------------------------|
! compute flux at z-direction
! set L/R state at z-direction
!
  mdir = 3
  
  call lr_state__MP5(mdir,ix,jx,kx,ro,pr &
       ,vz,vx,vy,bz,bx,by,phi &
       ,ch,gm,row,prw,vzw,vxw,vyw,bzw,bxw,byw,phiw,ccx,ccy,ccz)

  call flux_calc__bp(ix,jx,kx,bzw,phiw &
       ,bz_m,phi_m,ch)

  call flux_calc__glm(bz_m,phi_m,ch,fbzz,fphiz,ix,jx,kx)

  call flux_calc__hlld(row,prw,vzw,vxw,vyw,bz_m,bxw,byw,gm,margin,ix,jx,kx &
       ,froz,feez,frzz,frxz,fryz,fbxz,fbyz)

  call flux_calc__fbres(mdir,margin,ix,jx,kx,fbxz,cury,eta,-1.0d0 &
       ,fbxzr)
  call flux_calc__fbres(mdir,margin,ix,jx,kx,fbyz,curx,eta,+1.0d0 &
       ,fbyzr)
  call flux_calc__feres(mdir,margin,ix,jx,kx,feez,curx,cury,curz,bx,by,bz,eta &
       ,feezr)

!-----Step 2.---------------------------------------------------------|
! TVDRK substep
  k1 = fac*(-7.D0*n*n+30.D0*n-23.D0)
  k2 = fac*(+7.D0*n*n-30.D0*n+35.D0)
  do k=margin+1,kx-margin
     do j=margin+1,jx-margin
        do i=margin+1,ix-margin

           dtodx = dt/dx(i)
           dtody = dt/dy(j)
           dtodz = dt/dz(k)

           ro(i,j,k) = k1*ro1(i,j,k)+k2*(+ro(i,j,k) &
                +dtodx*(frox(i-1,j,k)-frox(i,j,k))  &
                +dtody*(froy(i,j-1,k)-froy(i,j,k))  &
                +dtodz*(froz(i,j,k-1)-froz(i,j,k)))
           ee(i,j,k) = k1*ee1(i,j,k)+k2*(+ee(i,j,k)  &
                +dtodx*(feexr(i-1,j,k)-feexr(i,j,k)) &
                +dtody*(feeyr(i,j-1,k)-feeyr(i,j,k)) &
                +dtodz*(feezr(i,j,k-1)-feezr(i,j,k)))
           rx(i,j,k) = k1*rx1(i,j,k)+k2*(+rx(i,j,k) &
                +dtodx*(frxx(i-1,j,k)-frxx(i,j,k))  &
                +dtody*(frxy(i,j-1,k)-frxy(i,j,k))  &
                +dtodz*(frxz(i,j,k-1)-frxz(i,j,k)) )
           ry(i,j,k) = k1*ry1(i,j,k)+k2*(+ry(i,j,k) &
                +dtodx*(fryx(i-1,j,k)-fryx(i,j,k))  &
                +dtody*(fryy(i,j-1,k)-fryy(i,j,k))  &
                +dtodz*(fryz(i,j,k-1)-fryz(i,j,k)) )
           rz(i,j,k) = k1*rz1(i,j,k)+k2*(+rz(i,j,k) &
                +dtodx*(frzx(i-1,j,k)-frzx(i,j,k))  &
                +dtody*(frzy(i,j-1,k)-frzy(i,j,k))  &
                +dtodz*(frzz(i,j,k-1)-frzz(i,j,k)) )
           bx(i,j,k) = k1*bx1(i,j,k)+k2*(+bx(i,j,k)  &
                +dtodx*(fbxx(i-1,j,k)-fbxx(i,j,k))   &
                +dtody*(fbxyr(i,j-1,k)-fbxyr(i,j,k)) &
                +dtodz*(fbxzr(i,j,k-1)-fbxzr(i,j,k)) )
           by(i,j,k) = k1*by1(i,j,k)+k2*(+by(i,j,k)  &
                +dtodx*(fbyxr(i-1,j,k)-fbyxr(i,j,k)) &
                +dtody*(fbyy(i,j-1,k)-fbyy(i,j,k))   &
                +dtodz*(fbyzr(i,j,k-1)-fbyzr(i,j,k)) )
           bz(i,j,k) = k1*bz1(i,j,k)+k2*(+bz(i,j,k)  &
                +dtodx*(fbzxr(i-1,j,k)-fbzxr(i,j,k)) &
                +dtody*(fbzyr(i,j-1,k)-fbzyr(i,j,k)) &
                +dtodz*(fbzz(i,j,k-1)-fbzz(i,j,k)) )
           phi(i,j,k) = k1*phi1(i,j,k)+k2*(+phi(i,j,k) &
                +dtodx*(fphix(i-1,j,k)-fphix(i,j,k))   &
                +dtody*(fphiy(i,j-1,k)-fphiy(i,j,k))   &
                +dtodz*(fphiz(i,j,k-1)-fphiz(i,j,k))   &
                )*exp(-dt*ch**2/cp**2)

        enddo
     enddo
  enddo

!-----Step 3.----------------------------------------------------------|
! conserved to primitive
!
  call convert__ctop(ix,jx,kx,gm,ro,ee,rx,ry,rz,bx,by,bz &
                    ,vx,vy,vz,pr)
  call bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  enddo

  end subroutine integrate__TVDRK3


end module integrate
