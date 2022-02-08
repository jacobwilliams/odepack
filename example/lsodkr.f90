!>
!! Demonstration program for the DLSODKR package.
!! This is the version of 27 April 2005.
!!
!! This version is in double precision.
!!
!! An ODE system is generated from the following 2-species diurnal
!! kinetics advection-diffusion PDE system in 2 space dimensions:
!!```text
!! dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
!!                 + Ri(c1,c2,t)      for i = 1,2,   where
!!   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
!!   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
!!   Kv(z) = Kv0*exp(z/5) ,
!!```
!! Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
!! vary diurnally.  The species are oxygen singlet and ozone.
!! The problem is posed on the square
!!```text
!!   0 .le. x .le. 20,    30 .le. z .le. 50   (all in km),
!!```
!! with homogeneous Neumann boundary conditions, and for time t in
!!```text
!!   0 .le. t .le. 86400 sec (1 day).
!!```
!! The PDE system is treated by central differences on a uniform
!! 10 x 10 mesh, with simple polynomial initial profiles.
!!
!! The problem is solved with DLSODKR, with the BDF/GMRES method and
!! the block-diagonal part of the Jacobian as a left preconditioner.
!! At intervals of 7200 sec (2 hrs), output includes sample values
!! of c1, c2, and c2tot = total of all c2 values.
!!
!! Roots of the function g = d(c2tot)/dt are found, i.e. the points at
!! which the total c2 (ozone) is stationary.
!!
!! Note: The preconditioner routines call LINPACK routines DGEFA, DGESL,
!! and the BLAS routines DCOPY and DSCAL.
!-----------------------------------------------------------------------
program lsodkr
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: a3 , a4 , c3 , dz , haco , hdco , om , q1 , q2 ,  &
 & q3 , q4 , vdco
integer :: mm , mx , mz
common /pcom  / q1 , q2 , q3 , q4 , a3 , a4 , om , c3 , dz ,      &
 & hdco , vdco , haco , mx , mz , mm
!
real(kind=dp) :: atol , avdim , c2tot , cx , cz , dx , t , tout ,  &
 & x , z
real(kind=dp) , save :: dkh , dkv0 , floor , halfda , pi , rtol ,  &
 & twohr , vel
integer :: iout , istate , jroot(1) , jx , jz , leniw , lenrw ,      &
 & ncfl , ncfn , neq , nfe , nge , njev , nli , nni ,     &
 & npe , nps , nsfi , nst
integer , dimension(235) :: iwork
integer , save :: jacflg , jpre , liw , lrw , mf
real(kind=dp) , dimension(4264) :: rwork
real(kind=dp) , dimension(2,10,10) :: y
external fdem , gdem , jacbd , solbd
!
data dkh/4.0d-6/ , vel/0.001d0/ , dkv0/1.0d-8/ , halfda/4.32d4/ , &
 & pi/3.1415926535898d0/ , twohr/7200.0d0/ , rtol/1.0d-5/ ,     &
 & floor/100.0d0/ , lrw/4264/ , liw/235/ , mf/22/ , jpre/1/ ,   &
 & jacflg/1/
!
! Load Common block of problem parameters.
      mx = 10
      mz = 10
      mm = mx*mz
      q1 = 1.63d-16
      q2 = 4.66d-16
      a3 = 22.62d0
      a4 = 7.601d0
      om = pi/halfda
      c3 = 3.7d16
      dx = 20.0d0/(mx-1.0d0)
      dz = 20.0d0/(mz-1.0d0)
      hdco = dkh/dx**2
      haco = vel/(2.0d0*dx)
      vdco = (1.0d0/dz**2)*dkv0
! Set other input arguments.
      atol = rtol*floor
      neq = 2*mx*mz
      iwork(1) = 8*mx*mz
      iwork(2) = neq
      iwork(3) = jpre
      iwork(4) = jacflg
      t = 0.0d0
      tout = 0.0d0
      istate = 1
! Initialize values for tout = 0 output.
      iwork(11) = 0
      iwork(14) = 0
      rwork(11) = 0.0d0
! Set initial profiles.
      do jz = 1 , mz
         z = 30.0d0 + (jz-1.0d0)*dz
         cz = (0.1d0*(z-40.0d0))**2
         cz = 1.0d0 - cz + 0.5d0*cz**2
         do jx = 1 , mx
            x = (jx-1.0d0)*dx
            cx = (0.1d0*(x-10.0d0))**2
            cx = 1.0d0 - cx + 0.5d0*cx**2
            y(1,jx,jz) = 1.0d6*cx*cz
            y(2,jx,jz) = 1.0d12*cx*cz
         enddo
      enddo
!
! Write heading, problem parameters, solution parameters.
      write (6,99001) mx , mz , mf , rtol , atol
99001 format ('Demonstration program for DLSODKR package'//             &
             &'2D diurnal kinetics-transport PDE system with 2 species'/&
             &'Spatial mesh is',i3,' by',i3/'Method flag is mf =',i3,   &
             &'   Tolerances are rtol =',d8.1,'   atol =',d8.1/         &
             &'Left preconditioner uses block-diagonal part of Jacobian'&
            & /'Root function finds stationary points of total ozone,'/ &
             &'  i.e. roots of (d/dt)(sum of c2 over all mesh points)'/)
!
! Loop over output points, call DLSODKR, print sample solution values.
      do iout = 1 , 13
         do
            call dlsodkr(fdem,[neq],y,t,tout,1,[rtol],[atol],1,istate,0,rwork,&
                       & lrw,iwork,liw,jacbd,solbd,mf,gdem,1,jroot)
            write (6,99002) t , iwork(11) , iwork(14) , rwork(11)
99002       format (/' t =',d10.3,4x,'no. steps =',i5,'   order =',i2,  &
                   &'   stepsize =',d10.3)
            call c2sum(y,mx,mz,c2tot)
            write (6,99003) y(1,1,1) , y(1,5,5) , y(1,10,10) , y(2,1,1) &
                          & , y(2,5,5) , y(2,10,10)
99003       format ('   c1 (bot.left/middle/top rt.) =',                &
                   &3d12.3/'   c2 (bot.left/middle/top rt.) =',3d12.3)
            write (6,99004) c2tot , jroot
99004       format ('   total c2 =',d15.6,'   jroot =',                 &
                   &i2' (1 = root found, 0 = no root)')
            if ( istate<0 ) then
               write (6,99005) istate
99005          format ('DLSODKR returned istate = ',i3)
               goto 100
            endif
            if ( istate==3 ) then
               istate = 2
               cycle
            endif
            tout = tout + twohr
            exit
         enddo
      enddo
!
! Print final statistics.
 100  continue
      lenrw = iwork(17)
      leniw = iwork(18)
      nst = iwork(11)
      nsfi = iwork(24)
      nfe = iwork(12)
      nge = iwork(10)
      npe = iwork(13)
      njev = iwork(25)
      nps = iwork(21)
      nni = iwork(19)
      nli = iwork(20)
      avdim = real(nli)/real(nni)
      ncfn = iwork(22)
      ncfl = iwork(23)
      write (6,99006) lenrw , leniw , nst , nsfi , nfe , nge , npe ,    &
                    & njev , nps , nni , nli , avdim , ncfn , ncfl
99006 format (//' Final statistics:'/' rwork size =',i5,5x,             &
             &' iwork size =',i4/' number of steps        =',i5,5x,     &
             &'no. fnal. iter. steps  =',i5/' number of f evals.     =',&
            & i5,5x,'number of g evals.     =',                         &
             &i5/' number of prec. evals. =',i5,5x,                     &
             &'number of Jac. evals.  =',i5/' number of prec. solves =',&
            & i5/' number of nonl. iters. =',i5,5x,                     &
             &'number of lin. iters.  =',                               &
             &i5/' average Krylov subspace dimension (nli/nni)  =',     &
             &f8.4/' number of conv. failures:  nonlinear =',i3,        &
             &'  linear =',i3)
end program lsodkr

subroutine fdem(neq,t,y,ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: a3 , a4 , c3 , dz , haco , hdco , om , q1 , q2 ,  &
 & q3 , q4 , vdco
integer :: mm , mx , mz
common /pcom  / q1 , q2 , q3 , q4 , a3 , a4 , om , c3 , dz ,      &
 & hdco , vdco , haco , mx , mz , mm
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(2,*) :: y
real(kind=dp) , dimension(2,*) :: ydot
!
real(kind=dp) :: c1 , c1dn , c1lt , c1rt , c1up , c2 , c2dn ,      &
 & c2lt , c2rt , c2up , czdn , czup , horad1 ,       &
 & horad2 , hord1 , hord2 , qq1 , qq2 , qq3 , qq4 ,  &
 & rkin1 , rkin2 , s , vertd1 , vertd2 , zdn , zup
integer :: iblok , iblok0 , idn , ileft , iright , iup , jx , jz
!
! Set diurnal rate coefficients.
      s = sin(om*t)
      if ( s>0.0d0 ) then
         q3 = exp(-a3/s)
         q4 = exp(-a4/s)
      else
         q3 = 0.0d0
         q4 = 0.0d0
      endif
! Loop over all grid points.
      do jz = 1 , mz
         zdn = 30.0d0 + (jz-1.5d0)*dz
         zup = zdn + dz
         czdn = vdco*exp(0.2d0*zdn)
         czup = vdco*exp(0.2d0*zup)
         iblok0 = (jz-1)*mx
         idn = -mx
         if ( jz==1 ) idn = mx
         iup = mx
         if ( jz==mz ) iup = -mx
         do jx = 1 , mx
            iblok = iblok0 + jx
            c1 = y(1,iblok)
            c2 = y(2,iblok)
! Set kinetic rate terms.
            qq1 = q1*c1*c3
            qq2 = q2*c1*c2
            qq3 = q3*c3
            qq4 = q4*c2
            rkin1 = -qq1 - qq2 + 2.0d0*qq3 + qq4
            rkin2 = qq1 - qq2 - qq4
! Set vertical diffusion terms.
            c1dn = y(1,iblok+idn)
            c2dn = y(2,iblok+idn)
            c1up = y(1,iblok+iup)
            c2up = y(2,iblok+iup)
            vertd1 = czup*(c1up-c1) - czdn*(c1-c1dn)
            vertd2 = czup*(c2up-c2) - czdn*(c2-c2dn)
! Set horizontal diffusion and advection terms.
            ileft = -1
            if ( jx==1 ) ileft = 1
            iright = 1
            if ( jx==mx ) iright = -1
            c1lt = y(1,iblok+ileft)
            c2lt = y(2,iblok+ileft)
            c1rt = y(1,iblok+iright)
            c2rt = y(2,iblok+iright)
            hord1 = hdco*(c1rt-2.0d0*c1+c1lt)
            hord2 = hdco*(c2rt-2.0d0*c2+c2lt)
            horad1 = haco*(c1rt-c1lt)
            horad2 = haco*(c2rt-c2lt)
! Load all terms into ydot.
            ydot(1,iblok) = vertd1 + hord1 + horad1 + rkin1
            ydot(2,iblok) = vertd2 + hord2 + horad2 + rkin2
         enddo
      enddo
end subroutine fdem

subroutine gdem(neq,t,y,ng,gout)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: a3 , a4 , c3 , dz , haco , hdco , om , q1 , q2 ,  &
 & q3 , q4 , vdco
integer :: mm , mx , mz
common /pcom  / q1 , q2 , q3 , q4 , a3 , a4 , om , c3 , dz ,      &
 & hdco , vdco , haco , mx , mz , mm
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(2,*) :: y
integer :: ng
real(kind=dp) :: gout
!
real(kind=dp) :: c1 , c2 , c2dn , c2dot , c2lt , c2rt , c2up ,     &
 & czdn , czup , horad2 , hord2 , qq1 , qq2 , qq4 ,  &
 & rkin2 , s , sum , vertd2 , zdn , zup
integer :: iblok , iblok0 , idn , ileft , iright , iup , jx , jz
!
! This routine computes the rates for c2 and adds them.
!
! Set diurnal rate coefficient q4.
      s = sin(om*t)
      if ( s>0.0d0 ) then
         q4 = exp(-a4/s)
      else
         q4 = 0.0d0
      endif
      sum = 0.0d0
! Loop over all grid points.
      do jz = 1 , mz
         zdn = 30.0d0 + (jz-1.5d0)*dz
         zup = zdn + dz
         czdn = vdco*exp(0.2d0*zdn)
         czup = vdco*exp(0.2d0*zup)
         iblok0 = (jz-1)*mx
         idn = -mx
         if ( jz==1 ) idn = mx
         iup = mx
         if ( jz==mz ) iup = -mx
         do jx = 1 , mx
            iblok = iblok0 + jx
            c1 = y(1,iblok)
            c2 = y(2,iblok)
! Set kinetic rate term for c2.
            qq1 = q1*c1*c3
            qq2 = q2*c1*c2
            qq4 = q4*c2
            rkin2 = qq1 - qq2 - qq4
! Set vertical diffusion terms for c2.
            c2dn = y(2,iblok+idn)
            c2up = y(2,iblok+iup)
            vertd2 = czup*(c2up-c2) - czdn*(c2-c2dn)
! Set horizontal diffusion and advection terms for c2.
            ileft = -1
            if ( jx==1 ) ileft = 1
            iright = 1
            if ( jx==mx ) iright = -1
            c2lt = y(2,iblok+ileft)
            c2rt = y(2,iblok+iright)
            hord2 = hdco*(c2rt-2.0d0*c2+c2lt)
            horad2 = haco*(c2rt-c2lt)
! Load all terms into c2dot and sum.
            c2dot = vertd2 + hord2 + horad2 + rkin2
            sum = sum + c2dot
         enddo
      enddo
      gout = sum
end subroutine gdem

subroutine jacbd(f,neq,t,y,ysv,rewt,f0,f1,hl0,jok,bd,ipbd,ier)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: a3 , a4 , c3 , dz , haco , hdco , om , q1 , q2 ,  &
 & q3 , q4 , vdco
integer :: mm , mx , mz
common /pcom  / q1 , q2 , q3 , q4 , a3 , a4 , om , c3 , dz ,      &
 & hdco , vdco , haco , mx , mz , mm
!
real(kind=dp) :: bd(2,2,*)
real(kind=dp) :: c1 , c2 , czdn , czup , diag , temp , zdn , zup
real(kind=dp) , dimension(neq) :: f0 , f1 , rewt , ysv
real(kind=dp) :: hl0 , t
integer :: iblok , iblok0 , jx , jz , lenbd
integer :: ier , jok , neq
integer , dimension(2,*) :: ipbd
real(kind=dp) , dimension(2,*) :: y
external f
!
      lenbd = 4*mm
! If jok = 1, copy saved block-diagonal approximate Jacobian into bd.
      if ( jok==1 ) then
         call dcopy(lenbd,bd(1,1,mm+1),1,bd,1)
         goto 100
      endif
!
! If jok = -1, compute and save diagonal Jacobian blocks
!  (using q3 and q4 values computed on last f call).
      do jz = 1 , mz
         zdn = 30.0d0 + (jz-1.5d0)*dz
         zup = zdn + dz
         czdn = vdco*exp(0.2d0*zdn)
         czup = vdco*exp(0.2d0*zup)
         diag = -(czdn+czup+2.0d0*hdco)
         iblok0 = (jz-1)*mx
         do jx = 1 , mx
            iblok = iblok0 + jx
            c1 = y(1,iblok)
            c2 = y(2,iblok)
            bd(1,1,iblok) = (-q1*c3-q2*c2) + diag
            bd(1,2,iblok) = -q2*c1 + q4
            bd(2,1,iblok) = q1*c3 - q2*c2
            bd(2,2,iblok) = (-q2*c1-q4) + diag
         enddo
      enddo
      call dcopy(lenbd,bd,1,bd(1,1,mm+1),1)
! Scale by -hl0, add identity matrix and LU-decompose blocks.
 100  continue
      temp = -hl0
      call dscal(lenbd,temp,bd,1)
      do iblok = 1 , mm
         bd(1,1,iblok) = bd(1,1,iblok) + 1.0d0
         bd(2,2,iblok) = bd(2,2,iblok) + 1.0d0
         call dgefa(bd(1,1,iblok),2,2,ipbd(1,iblok),ier)
         if ( ier/=0 ) return
      enddo
end subroutine jacbd

subroutine solbd(neq,t,y,f0,wk,hl0,bd,ipbd,v,lr,ier)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: a3 , a4 , c3 , dz , haco , hdco , om , q1 , q2 ,  &
 & q3 , q4 , vdco
integer :: mm , mx , mz
common /pcom  / q1 , q2 , q3 , q4 , a3 , a4 , om , c3 , dz ,      &
 & hdco , vdco , haco , mx , mz , mm
!
integer :: neq
real(kind=dp) :: t
real(kind=dp), dimension(neq) :: y
real(kind=dp), dimension(neq) :: f0
real(kind=dp), dimension(neq) :: wk
real(kind=dp) :: hl0
real(kind=dp), dimension(2,2,*) :: bd
integer , dimension(2,*) :: ipbd
real(kind=dp), dimension(2,*) :: v
integer :: lr
integer :: ier
!
integer :: i
!
! Solve the block-diagonal system Px = v using LU factors stored in bd
! and pivot data in ipbd, and return the solution in v.
      ier = 0
      do i = 1 , mm
         call dgesl(bd(1,1,i),2,2,ipbd(1,i),v(1,i),0)
      enddo
end subroutine solbd

subroutine c2sum(y,mx,mz,c2tot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp)     :: y(2,mx,mz)
integer           :: mx
integer           :: mz
real(kind=dp)     :: c2tot
!
integer           :: jx , jz
real(kind=dp)     :: sum
!
! Sum the c2 values.
      sum = 0.0d0
      do jz = 1 , mz
         do jx = 1 , mx
            sum = sum + y(2,jx,jz)
         enddo
      enddo
      c2tot = sum
end subroutine c2sum
