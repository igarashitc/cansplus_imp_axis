module init
  use const
  implicit none
  private

  public :: initialize

  real(8),public,dimension(ix) :: x,dx
  real(8),public,dimension(0:ix) :: xm
  real(8),public,dimension(jx) :: y,dy
  real(8),public,dimension(kx) :: z,dz
  real(8),public,dimension(5,2,ix) :: ccx
  real(8),public,dimension(5,2,jx) :: ccy
  real(8),public,dimension(5,2,kx) :: ccz

! cell center variables
  real(8),public,dimension(ix,jx,kx) :: ro,pr
  real(8),public,dimension(ix,jx,kx) :: vx,vy,vz
  real(8),public,dimension(ix,jx,kx) :: bx,by,bz,phi
  real(8),public,dimension(ix,jx,kx) :: eta

! gravitation
  real(8),public,dimension(ix,jx,kx) :: gx,gy,gz

  real(8),public :: ch,cp,min_dx

!----------------------------------------------------------------------|
!  initialize counters

  integer,public :: nd,ns,ns1,ns2
  real(8),public :: time,timep

!----------------------------------------------------------------------|
!   time control parameters
  integer,public :: mwflag,mw,nt1,nt2
  real(8),public :: dt,dtg


contains


  subroutine initialize
  use lr_state, only : reconstructionConstant, reconstructionConstant_cyl
  use model, only : model_setup
  use mpi_setup, only : mpi_setup__init, mpi_setup__init_cyl, mpid, mnull
  use bnd

  real(8),dimension(0:jx) :: ym
  real(8),dimension(0:kx) :: zm

!----------------------------------------------------------------------|
!   for MPI
  call mpi_setup__init(mpisize_x,mpisize_y,mpisize_z,pbcheck)
  call mpi_setup__init_cyl(mpisize_y)
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!   setup numerical model (grid, initial conditions, etc.)
!----------------------------------------------------------------------|
  call model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
       ,x,dx,xm,y,dy,ym,z,dz,zm,gx,gz,eta,min_dx) 

  call bnd__exec(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,x)

!----------------------------------------------------------------------|
!  cal reconstruction constant for MP5
!----------------------------------------------------------------------|
  call reconstructionConstant_cyl(margin,ix,xm,dx,ccx)
  if(mpid%l == mnull)then
    ccx(1:5,1,margin+1) = (/0.0D0,0.0D0,1.0D0,0.0D0,0.0D0/)
    ccx(1:5,2,margin  ) = (/0.0D0,0.0D0,1.0D0,0.0D0,0.0D0/)
  endif
  call reconstructionConstant(margin,jx,ym,dy,ccy)
  call reconstructionConstant(margin,kx,zm,dz,ccz)
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!  initialize counters
!----------------------------------------------------------------------|
  nd=1
  time  = 0.0d0
  timep = 0.0d0
  ns    = 0
  dt = tend

  end subroutine initialize

end module init
