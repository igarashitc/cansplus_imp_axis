!======================================================================
module openfile
!======================================================================
! Name :: openfile
!         file open
!
!======================================================================
  use dac_header
  use const, only : input_dir,output_dir

  implicit none
  private

  public :: file_input, file_output, file_output_param

  integer,public :: mf_params
  integer,public :: mf_x,mf_y,mf_z

  integer,public :: mf_ro,mf_pr       ! primitive variables
  integer,public :: mf_vx,mf_vy,mf_vz !
  integer,public :: mf_bx,mf_by,mf_bz !
  integer,public :: mf_phi,mf_eta

  integer,public :: mfi_x,mfi_y,mfi_z ! coordianates

  integer,public :: mfi_ro,mfi_pr       ! primitive variables
  integer,public :: mfi_vx,mfi_vy,mfi_vz !
  integer,public :: mfi_bx,mfi_by,mfi_bz !
  integer,public :: mfi_phi,mfi_eta
  
  integer,public :: mf_gx,mf_gy,mf_gz ! graviation

  integer,public :: mf_lambda1,mf_lambda2

  character,private :: cno*4
  character,private :: cnond*4

contains

  subroutine file_input(nd,mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,ix,jx,kx)

    integer,intent(in) :: nd,mpirank,ix,jx,kx
    real(8),intent(inout),dimension(ix,jx,kx) :: ro,pr,vx,vy,vz,bx,by,bz,phi,eta
    integer :: mt,nx0,ix0,jx0,kx0

    mt=6
    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mfi_ro=70
    call dacopnr3s(mfi_ro,input_dir//cnond//'_ro_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_ro) ro
    close(mfi_ro)

    mfi_pr=71
    call dacopnr3s(mfi_pr,input_dir//cnond//'_pr_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_pr) pr
    close(mfi_pr)

    mfi_vx=72
    call dacopnr3s(mfi_vx,input_dir//cnond//'_vx_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_vx) vx
    close(mfi_vx)

    mfi_vy=73
    call dacopnr3s(mfi_vy,input_dir//cnond//'_vy_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_vy) vy
    close(mfi_vy)

    mfi_vz=74
    call dacopnr3s(mfi_vz,input_dir//cnond//'_vz_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_vz) vz
    close(mfi_vz)

    mfi_bx=75
    call dacopnr3s(mfi_bx,input_dir//cnond//'_bx_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_bx) bx
    close(mfi_bx)

    mfi_by=76
    call dacopnr3s(mfi_by,input_dir//cnond//'_by_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_by) by
    close(mfi_by)

    mfi_bz=77
    call dacopnr3s(mfi_bz,input_dir//cnond//'_bz_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_bz) bz
    close(mfi_bz)

    mfi_phi=78
    call dacopnr3s(mfi_phi,input_dir//cnond//'_phi_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_phi) phi
    close(mfi_phi)

    mfi_eta=79
    call dacopnr3s(mfi_eta,input_dir//cnond//'_eta_rank='//cno//'.dac',mt,ix0,jx0,kx0,nx0)
    read(mfi_eta) eta
    close(mfi_eta)
     
  end subroutine file_input

  subroutine file_output(nd,mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,ix,jx,kx)
    implicit none

    integer,intent(in) :: nd,mpirank,ix,jx,kx
    real(8),intent(in),dimension(ix,jx,kx) :: ro,pr,vx,vy,vz,bx,by,bz,phi,eta

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mf_ro=20
    call dacdef3s(mf_ro,output_dir//cnond//'_ro_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_ro) ro
    close(mf_ro) 

    mf_pr=21
    call dacdef3s(mf_pr,output_dir//cnond//'_pr_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_pr) pr
    close(mf_pr) 

    mf_vx=22
    call dacdef3s(mf_vx,output_dir//cnond//'_vx_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_vx) vx
    close(mf_vx)

    mf_vy=23
    call dacdef3s(mf_vy,output_dir//cnond//'_vy_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_vy) vy
    close(mf_vy)

    mf_vz=24
    call dacdef3s(mf_vz,output_dir//cnond//'_vz_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_vz) vz
    close(mf_vz)

    mf_bx=25
    call dacdef3s(mf_bx,output_dir//cnond//'_bx_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_bx) bx
    close(mf_bx)

    mf_by=26
    call dacdef3s(mf_by,output_dir//cnond//'_by_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_by) by
    close(mf_by)

    mf_bz=27
    call dacdef3s(mf_bz,output_dir//cnond//'_bz_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_bz) bz
    close(mf_bz)
    
    mf_phi=28
    call dacdef3s(mf_phi,output_dir//cnond//'_phi_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_phi) phi
    close(mf_phi)

    mf_eta=29
    call dacdef3s(mf_eta,output_dir//cnond//'_eta_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_eta) eta
    close(mf_eta)

  end subroutine file_output

  subroutine file_output_param(nd,dtout,tend,ix,jx,kx,igx,jgx,kgx,margin &
                              ,mpisize,mpirank,mpisize_x,mpisize_y,mpisize_z &
                              ,nrmlro,nrmlte,nrmlx,nrmlv,nrmlt,nrmlee,mass_bh,rg,rg_nrmlx &
                              ,RadCool,te_factor,rohalo,eta0,vc,gm,x,y,z,dx,dy,dz,gx,gz)

    integer,intent(in) :: nd,ix,jx,kx,igx,jgx,kgx,margin,mpisize,mpirank
    integer,intent(in) :: mpisize_x,mpisize_y,mpisize_z
    real(8),intent(in),dimension(ix) :: x,dx
    real(8),intent(in),dimension(jx) :: y,dy
    real(8),intent(in),dimension(kx) :: z,dz
    real(8),intent(in),dimension(ix,jx,kx) :: gx,gz
    real(8),intent(in) :: dtout,tend,nrmlro,nrmlte,nrmlx,nrmlv,nrmlt,nrmlee
    real(8),intent(in) :: mass_bh,rg,rg_nrmlx,RadCool,te_factor,rohalo,eta0,vc,gm

    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    mf_params=9
    open (mf_params,file=output_dir//'params_rank='//cno//'.txt',form='formatted')
    
    call dacputparamc(mf_params,'comment','model_machida,int_2')
    call dacputparami(mf_params,'ix',ix)
    call dacputparami(mf_params,'jx',jx)
    call dacputparami(mf_params,'kx',kx)
    call dacputparami(mf_params,'igx',igx)
    call dacputparami(mf_params,'jgx',jgx)
    call dacputparami(mf_params,'kgx',kgx)
    call dacputparami(mf_params,'margin',margin)
    call dacputparamd(mf_params,'dtout',dtout)
    call dacputparamd(mf_params,'tend',tend)
    call dacputparami(mf_params,'mpi',1)
    call dacputparami(mf_params,'mpisize',mpisize)
    call dacputparami(mf_params,'mpirank',mpirank)
    call dacputparami(mf_params,'mpix',mpisize_x)
    call dacputparami(mf_params,'mpiy',mpisize_y)
    call dacputparami(mf_params,'mpiz',mpisize_z)
    call dacputparami(mf_params,'beta',100)
    call dacputparamd(mf_params,'nrmlro',nrmlro)
    call dacputparamd(mf_params,'nrmlte',nrmlte)
    call dacputparamd(mf_params,'nrmlx',nrmlx)
    call dacputparamd(mf_params,'nrmlv',nrmlv)
    call dacputparamd(mf_params,'nrmlt',nrmlt)
    call dacputparamd(mf_params,'nrmlee',nrmlee)
    call dacputparamd(mf_params,'mass_bh',mass_bh)
    call dacputparamd(mf_params,'rg',rg)
    call dacputparamd(mf_params,'rg_nrmlx',rg_nrmlx)
    call dacputparamd(mf_params,'RadCool',RadCool)
    call dacputparamd(mf_params,'te_factor',te_factor)
    call dacputparamd(mf_params,'rohalo',rohalo)
    call dacputparamd(mf_params,'eta0',eta0)
    call dacputparamd(mf_params,'vc',vc)
    call dacputparamd(mf_params,'x(1)',x(1))
    call dacputparamd(mf_params,'y(1)',y(1))
    call dacputparamd(mf_params,'z(1)',z(1))
    call dacputparamd(mf_params,'dx(1)',dx(1))
    call dacputparamd(mf_params,'dy(1)',dy(1))
    call dacputparamd(mf_params,'dz(1)',dz(1))
    call dacputparamd(mf_params,'gm',gm)
    close(mf_params)

    mf_x=11
    call dacdef1d(mf_x,output_dir//'x_rank='//cno//'.dac',6,ix)
    write(mf_x) x
    close(mf_x)

    mf_y=12
    call dacdef1d(mf_y,output_dir//'y_rank='//cno//'.dac',6,jx)
    write(mf_y) y
    close(mf_y)

    mf_z=13
    call dacdef1d(mf_z,output_dir//'z_rank='//cno//'.dac',6,kx)
    write(mf_z) z
    close(mf_z)

    mf_gx=40
    call dacdef3s(mf_gx,output_dir//'gx_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_gx) gx
    close(mf_gx)

    close(mf_gy)

    mf_gz=42
    call dacdef3s(mf_gz,output_dir//'gz_rank='//cno//'.dac',6,ix,jx,kx)
    write(mf_gz) gz
    close(mf_gz)

  end subroutine file_output_param

end module openfile
