!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dstode(3f) - [M_odepack] Performs one step of an ODEPACK integration.
!!### SYNOPSIS
!!
!!    subroutine dstode(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,f,jac,pjac,slvs)
!!
!!     integer,parameter :: dp=kind(0.0d0)
!!     integer                     :: Neq(*)
!!     real(kind=dp),intent(inout) :: Y(*)
!!     integer                     :: Nyh
!!     real(kind=dp),intent(inout) :: Yh(Nyh,*)
!!     real(kind=dp),intent(inout) :: Yh1(*)
!!     real(kind=dp)               :: Ewt(*)
!!     real(kind=dp),intent(inout) :: Savf(*)
!!     real(kind=dp),intent(inout) :: Acor(*)
!!     real(kind=dp)               :: Wm(*)
!!     integer                     :: Iwm(*)
!!     external f
!!     external jac
!!     external pjac
!!     external slvs

!!### DESCRIPTION
!!
!!  DSTODE performs one step of the integration of an initial value
!!  problem for a system of ordinary differential equations.
!!
!!  Note:  DSTODE is independent of the value of the iteration method
!!  indicator dls1%MITER, when this is .ne. 0, and hence is independent
!!  of the type of chord method used, or the Jacobian structure.
!!
!!  Communication with DSTODE is done with the following variables:
!!
!!  NEQ
!!
!!  : integer array containing problem size in NEQ(1), and
!!  passed as the NEQ argument in all calls to F and JAC.
!!
!!  Y
!!
!!  : an array of length .ge. dls1%N used as the Y argument in
!!  all calls to F and JAC.
!!
!!  YH
!!
!!  : an NYH by LMAX array containing the dependent variables
!!  and their approximate scaled derivatives, where
!!  LMAX = dls1%MAXORD + 1.  YH(i,j+1) contains the approximate
!!  j-th derivative of y(i), scaled by dls1%h**j/factorial(j)
!!  (j = 0,1,...,NQ).  on entry for the first step, the first
!!  two columns of YH must be set from the initial values.
!!
!!  NYH
!!
!!  : a constant integer .ge. dls1%N, the first dimension of YH.
!!
!!  YH1
!!
!!  : a one-dimensional array occupying the same space as YH.
!!
!!  EWT
!!
!!  : an array of length dls1%N containing multiplicative weights
!!  for local error measurements.  Local errors in Y(i) are
!!  compared to 1.0/EWT(i) in various error tests.
!!
!!  SAVF
!!
!!  : an array of working storage, of length dls1%N.
!!  Also used for input of YH(*,dls1%MAXORD+2) when dls1%JSTART = -1
!!  and dls1%MAXORD .lt. the current order NQ.
!!
!!  ACOR
!!
!!  : a work array of length dls1%N, used for the accumulated
!!  corrections.  On a successful return, ACOR(i) contains
!!  the estimated one-step local error in Y(i).
!!
!!  WM,IWM
!!
!!  : real and integer work arrays associated with matrix
!!  operations in chord iteration (dls1%MITER .ne. 0).
!!
!!  PJAC
!!
!!  : name of routine to evaluate and preprocess Jacobian matrix
!!  and P = I - dls1%h*el0*JAC, if a chord method is being used.
!!
!!  SLVS
!!
!!  : name of routine to solve linear system in chord iteration.
!!
!!### GLOBAL VARIABLES
!!
!!  dls1%CCMAX
!!
!!  : maximum relative change in dls1%h*el0 before PJAC is called.
!!
!!  dls1%H
!!
!!  : the step size to be attempted on the next step.
!!  dls1%H is altered by the error control algorithm during the
!!  problem.  dls1%H can be either positive or negative, but its
!!  sign must remain constant throughout the problem.
!!
!!  dls1%HMIN
!!
!!  : the minimum absolute value of the step size dls1%h to be used.
!!
!!  dls1%HMXI
!!
!!  : inverse of the maximum absolute value of dls1%h to be used.
!!  dls1%HMXI = 0.0 is allowed and corresponds to an infinite hmax.
!!  dls1%HMIN and dls1%HMXI may be changed at any time, but will not
!!  take effect until the next change of dls1%h is considered.
!!
!!  dls1%TN
!!
!!  : the independent variable. dls1%TN is updated on each step taken.
!!
!!  dls1%JSTART
!!
!!  : an integer used for input only, with the following
!!  values and meanings:
!!            0  perform the first step.
!!        .gt.0  take a new step continuing from the last.
!!           -1  take the next step with a new value of dls1%H, dls1%MAXORD,
!!               dls1%N, dls1%METH, dls1%MITER, and/or matrix parameters.
!!           -2  take the next step with a new value of dls1%H,
!!               but with other inputs unchanged.
!!  On return, dls1%JSTART is set to 1 to facilitate continuation.
!!
!!  dls1%KFLAG
!!
!!  : a completion code with the following meanings:
!!            0  the step was succesful.
!!           -1  the requested error could not be achieved.
!!           -2  corrector convergence could not be achieved.
!!           -3  fatal error in PJAC or SLVS.
!!  A return with dls1%KFLAG = -1 or -2 means either
!!  abs(dls1%H) = dls1%HMIN or 10 consecutive failures occurred.
!!  On a return with dls1%KFLAG negative, the values of dls1%TN and
!!  the YH array are as of the beginning of the last
!!  step, and dls1%H is the last step size attempted.
!!
!!  dls1%MAXORD
!!
!!  : the maximum order of integration method to be allowed.
!!
!!  dls1%MAXCOR
!!
!!  : the maximum number of corrector iterations allowed.
!!
!!  dls1%MSBP
!!
!!  : maximum number of steps between PJAC calls (dls1%MITER .gt. 0).
!!
!!  dls1%MXNCF
!!
!!  : maximum number of convergence failures allowed.
!!
!!  dls1%METH/dls1%MITER
!!
!!  : the method flags.  See description in driver.
!!
!!  dls1%N
!!
!!  : the number of first-order differential equations.
!!
!!  The values of dls1%CCMAX, dls1%H, dls1%HMIN, dls1%HMXI, dls1%TN, dls1%JSTART, dls1%KFLAG, dls1%MAXORD,
!!  dls1%MAXCOR, dls1%MSBP, dls1%MXNCF, dls1%METH, dls1%MITER, and dls1%N are communicated via a global compound variable (DLS1)
!!
!-----------------------------------------------------------------------
! ### BEGIN PROLOGUE  DSTODE
! ### SUBSIDIARY
! ### PURPOSE  Performs one step of an ODEPACK integration.
! ### TYPE      DOUBLE PRECISION (SSTODE-S, DSTODE-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  DCFODE, DVNORM
! ### COMMON BLOCKS    DLS001
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890503  Minor cosmetic changes.  (FNF)
!     19930809  Renamed to allow single/double precision versions. (ACH)
!     20010418  Reduced size of Common block /DLS001/. (ACH)
!     20031105  Restored 'own' variables to Common block /DLS001/, to
!               enable interrupt/restart feature. (ACH)
! ### END PROLOGUE  DSTODE
! **End
!-----------------------------------------------------------------------
subroutine dstode(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,f,jac,pjac,slvs)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)

integer                     :: Neq(*)
real(kind=dp),intent(inout) :: Y(*)
integer                     :: Nyh
real(kind=dp),intent(inout) :: Yh(Nyh,*)
real(kind=dp),intent(inout) :: Yh1(*)
real(kind=dp)               :: Ewt(*)
real(kind=dp),intent(inout) :: Savf(*)
real(kind=dp),intent(inout) :: Acor(*)
real(kind=dp)               :: Wm(*)
integer                     :: Iwm(*)
external f
external jac
external pjac
external slvs

real(kind=dp) :: dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, rh, rhdn, rhsm, rhup, told
integer       :: i, i1, iredo, iret, j, jb, m, ncf, newq

dls1%kflag = 0
told = dls1%tn
ncf = 0
dls1%ierpj = 0
dls1%iersl = 0
dls1%jcur = 0
dls1%icf = 0
delp = 0.0D0
if ( dls1%jstart>0 ) goto 400
if ( dls1%jstart==-1 ) then
!-----------------------------------------------------------------------
!  The following block handles preliminaries needed when dls1%JSTART = -1.
!  IPUP is set to dls1%MITER to force a matrix update.
!  If an order increase is about to be considered (IALTH = 1),
!  IALTH is reset to 2 to postpone consideration one more step.
!  If the caller has changed dls1%METH, DCFODE is called to reset
!  the coefficients of the method.
!  If the caller has changed dls1%MAXORD to a value less than the current
!  order NQ, NQ is reduced to dls1%MAXORD, and a new dls1%H chosen accordingly.
!  If dls1%H is to be changed, YH must be rescaled.
!  If dls1%H or dls1%METH is being changed, IALTH is reset to L = NQ + 1
!  to prevent further changes in dls1%H for that many steps.
!-----------------------------------------------------------------------
   dls1%ipup = dls1%miter
   dls1%lmax = dls1%maxord + 1
   if ( dls1%ialth==1 ) dls1%ialth = 2
   if ( dls1%meth/=dls1%meo ) then
      call dcfode(dls1%meth,dls1%elco,dls1%tesco)
      dls1%meo = dls1%meth
      if ( dls1%nq<=dls1%maxord ) then
         dls1%ialth = dls1%l
         iret = 1
         goto 100
      endif
   elseif ( dls1%nq<=dls1%maxord ) then
      goto 200
   endif
   dls1%nq = dls1%maxord
   dls1%l = dls1%lmax
   do i = 1, dls1%l
      dls1%el(i) = dls1%elco(i,dls1%nq)
   enddo
   dls1%nqnyh = dls1%nq*Nyh
   dls1%rc = dls1%rc*dls1%el(1)/dls1%el0
   dls1%el0 = dls1%el(1)
   dls1%conit = 0.5D0/(dls1%nq+2)
   ddn = dvnorm(dls1%n,Savf,Ewt)/dls1%tesco(1,dls1%l)
   exdn = 1.0D0/dls1%l
   rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
   rh = min(rhdn,1.0D0)
   iredo = 3
   if ( dls1%h==dls1%hold ) then
      rh = max(rh,dls1%hmin/abs(dls1%h))
   else
      rh = min(rh,abs(dls1%h/dls1%hold))
      dls1%h = dls1%hold
   endif
   goto 300
else
   if ( dls1%jstart==-2 ) goto 200
!-----------------------------------------------------------------------
!  On the first call, the order is set to 1, and other variables are
!  initialized.  RMAX is the maximum ratio by which dls1%H can be increased
!  in a single step.  It is initially 1.E4 to compensate for the small
!  initial dls1%H, but then is normally equal to 10.  If a failure
!  occurs (in corrector convergence or error test), RMAX is set to 2
!  for the next increase.
!-----------------------------------------------------------------------
   dls1%lmax = dls1%maxord + 1
   dls1%nq = 1
   dls1%l = 2
   dls1%ialth = 2
   dls1%rmax = 10000.0D0
   dls1%rc = 0.0D0
   dls1%el0 = 1.0D0
   dls1%crate = 0.7D0
   dls1%hold = dls1%h
   dls1%meo = dls1%meth
   dls1%nslp = 0
   dls1%ipup = dls1%miter
   iret = 3
!-----------------------------------------------------------------------
!  DCFODE is called to get all the integration coefficients for the
!  current dls1%METH.  Then the EL vector and related constants are reset
!  whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
   call dcfode(dls1%meth,dls1%elco,dls1%tesco)
endif
 100  continue
do i = 1, dls1%l
   dls1%el(i) = dls1%elco(i,dls1%nq)
enddo
dls1%nqnyh = dls1%nq*Nyh
dls1%rc = dls1%rc*dls1%el(1)/dls1%el0
dls1%el0 = dls1%el(1)
dls1%conit = 0.5D0/(dls1%nq+2)
select case (iret)
case (2)
   rh = max(rh,dls1%hmin/abs(dls1%h))
   goto 300
case (3)
   goto 400
case default
endselect
!-----------------------------------------------------------------------
!  If dls1%H is being changed, the dls1%H ratio RH is checked against
!  RMAX, dls1%HMIN, and dls1%HMXI, and the YH array rescaled.  IALTH is set to
!  L = NQ + 1 to prevent a change of dls1%H for that many steps, unless
!  forced by a convergence or error test failure.
!-----------------------------------------------------------------------
 200  continue
if ( dls1%h==dls1%hold ) goto 400
rh = dls1%h/dls1%hold
dls1%h = dls1%hold
iredo = 3
 300  continue
rh = min(rh,dls1%rmax)
rh = rh/max(1.0D0,abs(dls1%h)*dls1%hmxi*rh)
r = 1.0D0
do j = 2, dls1%l
   r = r*rh
   do i = 1, dls1%n
      Yh(i,j) = Yh(i,j)*r
   enddo
enddo
dls1%h = dls1%h*rh
dls1%rc = dls1%rc*rh
dls1%ialth = dls1%l
if ( iredo==0 ) then
   dls1%rmax = 10.0D0
   goto 1200
endif
!-----------------------------------------------------------------------
!  This section computes the predicted values by effectively
!  multiplying the YH array by the Pascal Triangle matrix.
!  RC is the ratio of new to old values of the coefficient  dls1%H*EL(1).
!  When RC differs from 1 by more than dls1%CCMAX, IPUP is set to dls1%MITER
!  to force PJAC to be called, if a Jacobian is involved.
!  In any case, PJAC is called at least every dls1%MSBP steps.
!-----------------------------------------------------------------------
 400  continue
if ( abs(dls1%rc-1.0D0)>dls1%ccmax ) dls1%ipup = dls1%miter
if ( dls1%nst>=dls1%nslp+dls1%msbp ) dls1%ipup = dls1%miter
dls1%tn = dls1%tn + dls1%h
i1 = dls1%nqnyh + 1
do jb = 1, dls1%nq
   i1 = i1 - Nyh
! dir$ ivdep
   do i = i1, dls1%nqnyh
      Yh1(i) = Yh1(i) + Yh1(i+Nyh)
   enddo
enddo
!-----------------------------------------------------------------------
!  Up to dls1%MAXCOR corrector iterations are taken.  A convergence test is
!  made on the R.M.S. norm of each correction, weighted by the error
!  weight vector EWT.  The sum of the corrections is accumulated in the
!  vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 500  continue
m = 0
do i = 1, dls1%n
   Y(i) = Yh(i,1)
enddo
call f(Neq,dls1%tn,Y,Savf)
dls1%nfe = dls1%nfe + 1
if ( dls1%ipup>0 ) then
!-----------------------------------------------------------------------
!  If indicated, the matrix P = I - dls1%h*dls1%el(1)*J is reevaluated and
!  preprocessed before starting the corrector iteration.  IPUP is set
!  to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
   call pjac(Neq,Y,Yh,Nyh,Ewt,Acor,Savf,Wm,Iwm,f,jac)
   dls1%ipup = 0
   dls1%rc = 1.0D0
   dls1%nslp = dls1%nst
   dls1%crate = 0.7D0
   if ( dls1%ierpj/=0 ) goto 800
endif
do i = 1, dls1%n
   Acor(i) = 0.0D0
enddo
 600  continue
if ( dls1%miter/=0 ) then
!-----------------------------------------------------------------------
!  In the case of the chord method, compute the corrector error,
!  and solve the linear system with that as right-hand side and
!  P as coefficient matrix.
!-----------------------------------------------------------------------
   do i = 1, dls1%n
      Y(i) = dls1%h*Savf(i) - (Yh(i,2)+Acor(i))
   enddo
   call slvs(Wm,Iwm,Y,Savf)
   if ( dls1%iersl<0 ) goto 800
   if ( dls1%iersl>0 ) goto 700
   del = dvnorm(dls1%n,Y,Ewt)
   do i = 1, dls1%n
      Acor(i) = Acor(i) + Y(i)
      Y(i) = Yh(i,1) + dls1%el(1)*Acor(i)
   enddo
else
!-----------------------------------------------------------------------
!  In the case of functional iteration, update Y directly from
!  the result of the last function evaluation.
!-----------------------------------------------------------------------
   do i = 1, dls1%n
      Savf(i) = dls1%h*Savf(i) - Yh(i,2)
      Y(i) = Savf(i) - Acor(i)
   enddo
   del = dvnorm(dls1%n,Y,Ewt)
   do i = 1, dls1%n
      Y(i) = Yh(i,1) + dls1%el(1)*Savf(i)
      Acor(i) = Savf(i)
   enddo
endif
!-----------------------------------------------------------------------
!  Test for convergence.  If M.gt.0, an estimate of the convergence
!  rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
if ( m/=0 ) dls1%crate = max(0.2D0*dls1%crate,del/delp)
dcon = del*min(1.0D0,1.5D0*dls1%crate)/(dls1%tesco(2,dls1%nq)*dls1%conit)
if ( dcon<=1.0D0 ) then
!-----------------------------------------------------------------------
!  The corrector has converged.  JCUR is set to 0
!  to signal that the Jacobian involved may need updating later.
!  The local error test is made and control passes to statement 500
!  if it fails.
!-----------------------------------------------------------------------
   dls1%jcur = 0
   if ( m==0 ) dsm = del/dls1%tesco(2,dls1%nq)
   if ( m>0 ) dsm = dvnorm(dls1%n,Acor,Ewt)/dls1%tesco(2,dls1%nq)
   if ( dsm>1.0D0 ) then
!-----------------------------------------------------------------------
!  The error test failed.  dls1%KFLAG keeps track of multiple failures.
!  Restore dls1%TN and the YH array to their previous values, and prepare
!  to try the step again.  Compute the optimum step size for this or
!  one lower order.  After 2 or more failures, dls1%H is forced to decrease
!  by a factor of 0.2 or less.
!-----------------------------------------------------------------------
      dls1%kflag = dls1%kflag - 1
      dls1%tn = told
      i1 = dls1%nqnyh + 1
      do jb = 1, dls1%nq
         i1 = i1 - Nyh
! dir$ ivdep
         do i = i1, dls1%nqnyh
            Yh1(i) = Yh1(i) - Yh1(i+Nyh)
         enddo
      enddo
      dls1%rmax = 2.0D0
      if ( abs(dls1%h)<=dls1%hmin*1.00001D0 ) then
!-----------------------------------------------------------------------
!  All returns are made through this section.  dls1%H is saved in HOLD
!  to allow the caller to change dls1%H on the next step.
!-----------------------------------------------------------------------
         dls1%kflag = -1
         goto 1300
      elseif ( dls1%kflag<=-3 ) then
!-----------------------------------------------------------------------
!  Control reaches this section if 3 or more failures have occured.
!  If 10 failures have occurred, exit with dls1%KFLAG = -1.
!  It is assumed that the derivatives that have accumulated in the
!  YH array have errors of the wrong order.  Hence the first
!  derivative is recomputed, and the order is set to 1.  Then
!  dls1%H is reduced by a factor of 10, and the step is retried,
!  until it succeeds or dls1%H reaches dls1%HMIN.
!-----------------------------------------------------------------------
         if ( dls1%kflag==-10 ) then
            dls1%kflag = -1
            goto 1300
         else
            rh = 0.1D0
            rh = max(dls1%hmin/abs(dls1%h),rh)
            dls1%h = dls1%h*rh
            do i = 1, dls1%n
               Y(i) = Yh(i,1)
            enddo
            call f(Neq,dls1%tn,Y,Savf)
            dls1%nfe = dls1%nfe + 1
            do i = 1, dls1%n
               Yh(i,2) = dls1%h*Savf(i)
            enddo
            dls1%ipup = dls1%miter
            dls1%ialth = 5
            if ( dls1%nq==1 ) goto 400
            dls1%nq = 1
            dls1%l = 2
            iret = 3
            goto 100
         endif
      else
         iredo = 2
         rhup = 0.0D0
         goto 900
      endif
   else
!-----------------------------------------------------------------------
!  After a successful step, update the YH array.
!  Consider changing dls1%H if IALTH = 1.  Otherwise decrease IALTH by 1.
!  If IALTH is then 1 and NQ .lt. dls1%MAXORD, then ACOR is saved for
!  use in a possible order increase on the next step.
!  If a change in dls1%H is considered, an increase or decrease in order
!  by one is considered also.  A change in dls1%H is made only if it is by a
!  factor of at least 1.1.  If not, IALTH is set to 3 to prevent
!  testing for that many steps.
!-----------------------------------------------------------------------
      dls1%kflag = 0
      iredo = 0
      dls1%nst = dls1%nst + 1
      dls1%hu = dls1%h
      dls1%nqu = dls1%nq
      do j = 1, dls1%l
         do i = 1, dls1%n
            Yh(i,j) = Yh(i,j) + dls1%el(j)*Acor(i)
         enddo
      enddo
      dls1%ialth = dls1%ialth - 1
      if ( dls1%ialth==0 ) then
!-----------------------------------------------------------------------
!  Regardless of the success or failure of the step, factors
!  RHDN, RHSM, and RHUP are computed, by which dls1%H could be multiplied
!  at order NQ - 1, order NQ, or order NQ + 1, respectively.
!  In the case of failure, RHUP = 0.0 to avoid an order increase.
!  The largest of these is determined and the new order chosen
!  accordingly.  If the order is to be increased, we compute one
!  additional scaled derivative.
!-----------------------------------------------------------------------
         rhup = 0.0D0
         if ( dls1%l/=dls1%lmax ) then
            do i = 1, dls1%n
               Savf(i) = Acor(i) - Yh(i,dls1%lmax)
            enddo
            dup = dvnorm(dls1%n,Savf,Ewt)/dls1%tesco(3,dls1%nq)
            exup = 1.0D0/(dls1%l+1)
            rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
         endif
         goto 900
      else
         if ( dls1%ialth<=1 ) then
            if ( dls1%l/=dls1%lmax ) then
               do i = 1, dls1%n
                  Yh(i,dls1%lmax) = Acor(i)
               enddo
            endif
         endif
         goto 1200
      endif
   endif
else
   m = m + 1
   if ( m/=dls1%maxcor ) then
      if ( m<2 .or. del<=2.0D0*delp ) then
         delp = del
         call f(Neq,dls1%tn,Y,Savf)
         dls1%nfe = dls1%nfe + 1
         goto 600
      endif
   endif
endif
!-----------------------------------------------------------------------
!  The corrector iteration failed to converge.
!  If dls1%MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
!  the next try.  Otherwise the YH array is retracted to its values
!  before prediction, and dls1%H is reduced, if possible.  If dls1%H cannot be
!  reduced or dls1%MXNCF failures have occurred, exit with dls1%KFLAG = -2.
!-----------------------------------------------------------------------
 700  continue
if ( dls1%miter/=0 .and. dls1%jcur/=1 ) then
   dls1%icf = 1
   dls1%ipup = dls1%miter
   goto 500
endif
 800  continue
dls1%icf = 2
ncf = ncf + 1
dls1%rmax = 2.0D0
dls1%tn = told
i1 = dls1%nqnyh + 1
do jb = 1, dls1%nq
   i1 = i1 - Nyh
! dir$ ivdep
   do i = i1, dls1%nqnyh
      Yh1(i) = Yh1(i) - Yh1(i+Nyh)
   enddo
enddo
if ( dls1%ierpj<0 .or. dls1%iersl<0 ) then
   dls1%kflag = -3
   goto 1300
elseif ( abs(dls1%h)<=dls1%hmin*1.00001D0 ) then
   dls1%kflag = -2
   goto 1300
elseif ( ncf==dls1%mxncf ) then
   dls1%kflag = -2
   goto 1300
else
   rh = 0.25D0
   dls1%ipup = dls1%miter
   iredo = 1
   rh = max(rh,dls1%hmin/abs(dls1%h))
   goto 300
endif
 900  continue
exsm = 1.0D0/dls1%l
rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
rhdn = 0.0D0
if ( dls1%nq/=1 ) then
   ddn = dvnorm(dls1%n,Yh(1,dls1%l),Ewt)/dls1%tesco(1,dls1%nq)
   exdn = 1.0D0/dls1%nq
   rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
endif
if ( rhsm>=rhup ) then
   if ( rhsm>=rhdn ) then
      newq = dls1%nq
      rh = rhsm
      goto 1000
   endif
elseif ( rhup>rhdn ) then
   newq = dls1%l
   rh = rhup
   if ( rh<1.1D0 ) then
      dls1%ialth = 3
      goto 1200
   else
      r = dls1%el(dls1%l)/dls1%l
      do i = 1, dls1%n
         Yh(i,newq+1) = Acor(i)*r
      enddo
      goto 1100
   endif
endif
newq = dls1%nq - 1
rh = rhdn
if ( dls1%kflag<0 .and. rh>1.0D0 ) rh = 1.0D0
 1000 continue
if ( (dls1%kflag==0) .and. (rh<1.1D0) ) then
   dls1%ialth = 3
   goto 1200
else
   if ( dls1%kflag<=-2 ) rh = min(rh,0.2D0)
!-----------------------------------------------------------------------
!  If there is a change of order, reset NQ, dls1%l, and the coefficients.
!  In any case dls1%H is reset according to RH and the YH array is rescaled.
!  Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
   if ( newq==dls1%nq ) then
      rh = max(rh,dls1%hmin/abs(dls1%h))
      goto 300
   endif
endif
 1100 continue
dls1%nq = newq
dls1%l = dls1%nq + 1
iret = 2
goto 100
 1200 continue
r = 1.0D0/dls1%tesco(2,dls1%nqu)
do i = 1, dls1%n
   Acor(i) = Acor(i)*r
enddo
 1300 continue
dls1%hold = dls1%h
dls1%jstart = 1
end subroutine dstode
