module const

  implicit none
 
  integer,parameter :: margin=3
  integer,parameter :: ix = 32+2*margin,jx=128+2*margin,kx=2+2*margin
  real(8),parameter :: pi= acos(-1.0d0)            ! circular constant
  real(8),parameter :: pi2=2.0d0*pi

!--------------------------------------------------------------------------
!  MPI
  integer,parameter :: mpisize_x=4
  integer,parameter :: mpisize_y=4
  integer,parameter :: mpisize_z=1
  integer,parameter :: igx = ix*mpisize_x-2*margin*(mpisize_x-1)
  integer,parameter :: jgx = jx*mpisize_y-2*margin*(mpisize_y-1)
  integer,parameter :: kgx = kx*mpisize_z-2*margin*(mpisize_z-1)
  !TRUE if periodic boundary condition is applied. (1:x, 2:y, 3:z)
  logical,parameter :: pbcheck(3) = (/.false., .true., .false./) 

!--------------------------------------------------------------------------
!   model of accretion flow
  logical,parameter :: pol_2D = .true.     ! if ".true.", then 2D model with poloidal B-field is employed
  logical,parameter :: tor_3D = .false.     ! if ".true.", then 3D model with toroidal B-field is employed
!--------------------------------------------------------------------------
!   time control parameters

  logical :: restart = .false.     ! if ".true." then start from restart data
  character(*),parameter :: input_dir = "./data/", output_dir = "./data/"
  real(8),parameter :: tend  = 0.1d0
  real(8),parameter :: dtout = tend/1.d1
  integer,parameter :: nstop = int(2.d6)       !number of total time steps for the run

  real(8),parameter :: swtch_t=pi2*1.d2    !  cooling switching 

  real(8),parameter :: safety=0.3d0
  real(8),parameter :: dtmin=1.d-10

!--------------------------------------------------------------------------
! model parameter

  !  size 
  real(8),parameter :: xmin = 0.0d0  
  real(8),parameter :: xmax = 1.0d0
  real(8),parameter :: ymin = 0.0d0
  real(8),parameter :: ymax = pi2
  real(8),parameter :: zmin = 0.0d0  
  real(8),parameter :: zmax = 1.0d0  

  ! set uniform grid
  real(8),parameter :: dxg0 = (xmax-xmin)/real(igx-margin*2)
  real(8),parameter :: dyg0 = (ymax-ymin)/real(jgx-margin*2)
  real(8),parameter :: dzg0 = (ymax-ymin)/real(kgx-margin*2)

  ! set global non-uniform grid  
  real(8),parameter :: dxmax = 10.0d0*dxg0
  real(8),parameter :: ratio_x = 1.0d0     ! in case of uniform, ratio_x=1.0d0
  real(8),parameter :: ugrid_x = 1.0d0       !uniform grid spacing region in r
  real(8),parameter :: dzmax = 10.0d0*dxg0
  real(8),parameter :: ratio_z = 1.0d0     ! in case of uniform, ratio_z=1.0d0
  real(8),parameter :: ugrid_z = 0.5d0       !uniform grid spacing region in |z|

  real(8),parameter :: gm = 5d0/3d0

end module const 
