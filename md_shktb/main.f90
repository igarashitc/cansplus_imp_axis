!##############################################
! Model : 1D MHD shock tube problem
! Coordinate system : Cartesian grid
! Reference : Ryu et al., Astrophys. J., 1995.
!##############################################
program main

  use mpi_setup
  use openfile
  use const
  use init
  use getNewdt
  use integrate, only : integrate__TVDRK3

  implicit none

!----------------------------------------------------------------------|
!  initial setting
  call initialize
!----------------------------------------------------------------------|

!----------------------------------------------------------------------|
!  output data
  call file_output_param(nd,dtout,tend,ix,jx,kx,igx,jgx,kgx,margin &
                        ,mpid%mpisize,mpid%mpirank                 &
                        ,mpisize_x,mpisize_y,mpisize_z             &
                        ,gm,x,y,z,dx,dy,dz)
  call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                  ,ix,jx,kx)


!----------------------------------------------------------------------|
!  read-data
  if(restart)then
     open(91,file=input_dir//'readFileNumber.dat')
     read(91,*) nd
     close(91)
     call file_input(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                    ,ix,jx,kx)
     time = real(nd-1)*dtout
  endif
 
!======================================================================|
!     time integration 
!======================================================================|
  nd=nd+1
  loop: do ns=1,nstop
       
     mwflag=0

!----------------------------------------------------------------------|
!     obtain time spacing
!----------------------------------------------------------------------|
     call getNewdt__glm(margin,safety,dtmin,ix,jx,kx,gm,ro,pr &
                       ,vx,vy,vz,bx,by,bz,x,dx,y,dy,z,dz,eta  &
                       ,dt,ch,cp,min_dx)

     if(merr /= 0) exit loop
     timep = time
     time = time+dt

!---- integrate--------------------------------------------------------|
     call integrate__TVDRK3(margin,ix,jx,kx,gm,dx,dy,dz,dt    &
                           ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cp &
                           ,eta,ccx,ccy,ccz)
!----------------------------------------------------------------------|

!----- check output
! dtout
     mw=0
     nt1=int(timep/dtout)
     nt2=int(time/dtout)
     if(nt1 < nt2) mw=1
     if(mw /= 0) then
        call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                        ,ix,jx,kx)

        if(mpid%mpirank == 0)then
           write(6,913) ns,time,nd
           write(*,*) '[NORMAL] dt :: ',dt
           write(*,*) 'nd: ',nd
        endif

        open(91,file=output_dir//'readFileNumber.dat')
        write(91,*) nd
        close(91)

        nd=nd+1
        mwflag=1
     endif

     if(time > tend) exit loop

  enddo loop
!======================================================================|
! end of time integration
!======================================================================|

!  data output
  if(mwflag == 0) then
     call file_output(nd,mpid%mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta &
                     ,ix,jx,kx)
  endif

  call mpi_finalize(merr)

  write(6,915) ns,time

  if(merr == 0) then
     write(6,*) '  ### normal stop ###'
  else
     write(6,*) '  ### abnormal stop ###'
  endif
  
913 format (1x,' write    ','step=',i8,' time=',e10.3,' nd =',i3)
915 format (1x,' stop     ','step=',i8,' time=',e10.3)


end program main
