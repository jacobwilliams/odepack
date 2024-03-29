!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DSTODA performs one step of the integration of an initial value
!! problem for a system of ordinary differential equations.
!!
!! Note: DSTODA is independent of the value of the iteration method
!! indicator MITER, when this is .ne. 0, and hence is independent
!! of the type of chord method used, or the Jacobian structure.
!!
!! Communication with DSTODA is done with the following variables:
!!
!! Y
!!
!! : an array of length .ge. N used as the Y argument in
!! all calls to F and JAC.
!!
!! NEQ
!!
!! : integer array containing problem size in NEQ(1), and
!! passed as the NEQ argument in all calls to F and JAC.
!!
!! YH
!!
!! : an NYH by LMAX array containing the dependent variables
!! and their approximate scaled derivatives, where
!! LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!! j-th derivative of y(i), scaled by H**j/factorial(j)
!! (j = 0,1,...,NQ).  On entry for the first step, the first
!! two columns of YH must be set from the initial values.
!!
!! NYH
!!
!! : a constant integer .ge. N, the first dimension of YH.
!!
!! YH1
!!
!! : a one-dimensional array occupying the same space as YH.
!!
!! EWT
!!
!! : an array of length N containing multiplicative weights
!! for local error measurements.  Local errors in y(i) are
!! compared to 1.0/EWT(i) in various error tests.
!!
!! SAVF
!!
!! : an array of working storage, of length N.
!!
!! ACOR
!!
!! : a work array of length N, used for the accumulated
!! corrections.  On a successful return, ACOR(i) contains
!! the estimated one-step local error in y(i).
!!
!! WM,IWM
!!
!! : real and integer work arrays associated with matrix
!! operations in chord iteration (MITER .ne. 0).
!!
!! PJAC
!!
!! : name of routine to evaluate and preprocess Jacobian matrix
!! and P = I - H*EL0*Jac, if a chord method is being used.
!! It also returns an estimate of norm(Jac) in PDNORM.
!!
!! SLVS
!!
!! : name of routine to solve linear system in chord iteration.
!!
!! CCMAX
!!
!! : maximum relative change in H*EL0 before PJAC is called.
!!
!! H
!!
!! : the step size to be attempted on the next step.
!! H is altered by the error control algorithm during the
!! problem.  H can be either positive or negative, but its
!! sign must remain constant throughout the problem.
!!
!! HMIN
!!
!! : the minimum absolute value of the step size H to be used.
!!
!! HMXI
!!
!! : inverse of the maximum absolute value of H to be used.
!! HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
!! HMIN and HMXI may be changed at any time, but will not
!! take effect until the next change of H is considered.
!!
!! TN
!!
!! : the independent variable. TN is updated on each step taken.
!!
!! JSTART
!!
!! : an integer used for input only, with the following
!! values and meanings:
!!               0  perform the first step.
!!           .gt.0  take a new step continuing from the last.
!!              -1  take the next step with a new value of H,
!!                    N, METH, MITER, and/or matrix parameters.
!!              -2  take the next step with a new value of H,
!!                    but with other inputs unchanged.
!!          On return, JSTART is set to 1 to facilitate continuation.
!!
!! KFLAG
!!
!! : a completion code with the following meanings:
!!               0  the step was succesful.
!!              -1  the requested error could not be achieved.
!!              -2  corrector convergence could not be achieved.
!!              -3  fatal error in PJAC or SLVS.
!! A return with KFLAG = -1 or -2 means either
!! ABS(H) = HMIN or 10 consecutive failures occurred.
!! On a return with KFLAG negative, the values of TN and
!! the YH array are as of the beginning of the last
!! step, and H is the last step size attempted.
!!
!! MAXORD
!!
!! : the maximum order of integration method to be allowed.
!!
!! MAXCOR
!!
!! : the maximum number of corrector iterations allowed.
!!
!! MSBP
!!
!! : maximum number of steps between PJAC calls (MITER .gt. 0).
!!
!! MXNCF
!!
!! : maximum number of convergence failures allowed.
!!
!! METH
!!
!! : current method.
!!          METH = 1 means Adams method (nonstiff)
!!          METH = 2 means BDF method (stiff)
!!          METH may be reset by DSTODA.
!!
!! MITER
!!
!! : corrector iteration method.
!! MITER = 0 means functional iteration.
!! MITER = JT .gt. 0 means a chord iteration corresponding
!! to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
!! communicated here as JTYP, but is not used in DSTODA
!! except to load MITER following a method switch.)
!! MITER may be reset by DSTODA.
!!
!! N
!!
!! : the number of first-order differential equations.
!-----------------------------------------------------------------------
subroutine dstoda(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,f,jac,pjac,slvs)
!
integer                      :: Neq(*)
real(kind=dp),intent(inout)  :: Y(*)
integer                      :: Nyh
real(kind=dp),intent(inout)  :: Yh(Nyh,*)
real(kind=dp),intent(inout)  :: Yh1(*)
real(kind=dp)                :: Ewt(*)
real(kind=dp),intent(inout)  :: Savf(*)
real(kind=dp),intent(inout)  :: Acor(*)
real(kind=dp)                :: Wm(*)
integer                      :: Iwm(*)
external f
external jac
external pjac
external slvs

real(kind=dp) :: alpha, dcon, ddn, del, delp, dm1, dm2, dsm, dup, exdn, exm1, exm2, exsm, exup, pdh, pnorm, r,      &
              & rate, rh, rh1, rh1it, rh2, rhdn, rhsm, rhup, rm, told
integer :: i, i1, iredo, iret, j, jb, lm1, lm1p1, lm2, lm2p1, m, ncf, newq, nqm1, nqm2
real(kind=dp),parameter :: sm1(12)=&
 & [0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0, 0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0]

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
!  The following block handles preliminaries needed when JSTART = -1.
!  IPUP is set to MITER to force a matrix update.
!  If an order increase is about to be considered (IALTH = 1),
!  IALTH is reset to 2 to postpone consideration one more step.
!  If the caller has changed METH, DCFODE is called to reset
!  the coefficients of the method.
!  If H is to be changed, YH must be rescaled.
!  If H or METH is being changed, IALTH is reset to L = NQ + 1
!  to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
   dls1%ipup = dls1%miter
   dls1%lmax = dls1%maxord + 1
   if ( dls1%ialth==1 ) dls1%ialth = 2
   if ( dls1%meth==dlsa%mused ) goto 200
   call dcfode(dls1%meth,dls1%elco,dls1%tesco)
   dls1%ialth = dls1%l
   iret = 1
else
   if ( dls1%jstart==-2 ) goto 200
!-----------------------------------------------------------------------
!  On the first call, the order is set to 1, and other variables are
!  initialized.  RMAX is the maximum ratio by which H can be increased
!  in a single step.  It is initially 1.E4 to compensate for the small
!  initial H, but then is normally equal to 10.  If a failure
!  occurs (in corrector convergence or error test), RMAX is set at 2
!  for the next increase.
!  DCFODE is called to get the needed coefficients for both methods.
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
   dls1%nslp = 0
   dls1%ipup = dls1%miter
   iret = 3
!  Initialize switching parameters.  METH = 1 is assumed initially. -----
   dlsa%icount = 20
   dlsa%irflag = 0
   dlsa%pdest = 0.0D0
   dlsa%pdlast = 0.0D0
   dlsa%ratio = 5.0D0
   call dcfode(2,dls1%elco,dls1%tesco)
   do i = 1, 5
      dlsa%cm2(i) = dls1%tesco(2,i)*dls1%elco(i+1,i)
   enddo
   call dcfode(1,dls1%elco,dls1%tesco)
   do i = 1, 12
      dlsa%cm1(i) = dls1%tesco(2,i)*dls1%elco(i+1,i)
   enddo
endif
!-----------------------------------------------------------------------
!  The dls1%el vector and related constants are reset
!  whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
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
!  If H is being changed, the H ratio RH is checked against
!  RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
!  L = NQ + 1 to prevent a change of H for that many steps, unless
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
!-----------------------------------------------------------------------
!  If METH = 1, also restrict the new step size by the stability region.
!  If this reduces H, set IRFLAG to 1 so that if there are roundoff
!  problems later, we can assume that is the cause of the trouble.
!-----------------------------------------------------------------------
if ( dls1%meth/=2 ) then
   dlsa%irflag = 0
   pdh = max(abs(dls1%h)*dlsa%pdlast,0.000001D0)
   if ( rh*pdh*1.00001D0>=sm1(dls1%nq) ) then
      rh = sm1(dls1%nq)/pdh
      dlsa%irflag = 1
   endif
endif
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
   goto 1300
endif
!-----------------------------------------------------------------------
!  This section computes the predicted values by effectively
!  multiplying the YH array by the Pascal triangle matrix.
!  RC is the ratio of new to old values of the coefficient  H*EL(1).
!  When RC differs from 1 by more than CCMAX, IPUP is set to MITER
!  to force PJAC to be called, if a Jacobian is involved.
!  In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
 400  continue
if ( abs(dls1%rc-1.0D0)>dls1%ccmax ) dls1%ipup = dls1%miter
if ( dls1%nst>=dls1%nslp+dls1%msbp ) dls1%ipup = dls1%miter
dls1%tn = dls1%tn + dls1%h
i1 = dls1%nqnyh + 1
do jb = 1, dls1%nq
   i1 = i1 - Nyh
! DIR$ IVDEP
   do i = i1, dls1%nqnyh
      Yh1(i) = Yh1(i) + Yh1(i+Nyh)
   enddo
enddo
pnorm = dmnorm(dls1%n,Yh1,Ewt)
!-----------------------------------------------------------------------
!  Up to MAXCOR corrector iterations are taken.  A convergence test is
!  made on the RMS-norm of each correction, weighted by the error
!  weight vector EWT.  The sum of the corrections is accumulated in the
!  vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 500  continue
m = 0
rate = 0.0D0
del = 0.0D0
do i = 1, dls1%n
   Y(i) = Yh(i,1)
enddo
call f(Neq,dls1%tn,Y,Savf)
dls1%nfe = dls1%nfe + 1
if ( dls1%ipup>0 ) then
!-----------------------------------------------------------------------
!  If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
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
   del = dmnorm(dls1%n,Y,Ewt)
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
   del = dmnorm(dls1%n,Y,Ewt)
   do i = 1, dls1%n
      Y(i) = Yh(i,1) + dls1%el(1)*Savf(i)
      Acor(i) = Savf(i)
   enddo
endif
!-----------------------------------------------------------------------
!  Test for convergence.  If M .gt. 0, an estimate of the convergence
!  rate constant is stored in CRATE, and this is used in the test.
!
!  We first check for a change of iterates that is the size of
!  roundoff error.  If this occurs, the iteration has converged, and a
!  new rate estimate is not formed.
!  In all other cases, force at least two iterations to estimate a
!  local Lipschitz constant estimate for Adams methods.
!  On convergence, form PDEST = local maximum Lipschitz constant
!  estimate.  PDLAST is the most recent nonzero estimate.
!-----------------------------------------------------------------------
if ( del<=100.0D0*pnorm*dls1%uround ) goto 900
if ( m/=0 .or. dls1%meth/=1 ) then
   if ( m/=0 ) then
      rm = 1024.0D0
      if ( del<=1024.0D0*delp ) rm = del/delp
      rate = max(rate,rm)
      dls1%crate = max(0.2D0*dls1%crate,rm)
   endif
   dcon = del*min(1.0D0,1.5D0*dls1%crate)/(dls1%tesco(2,dls1%nq)*dls1%conit)
   if ( dcon<=1.0D0 ) then
      dlsa%pdest = max(dlsa%pdest,rate/abs(dls1%h*dls1%el(1)))
      if ( dlsa%pdest/=0.0D0 ) dlsa%pdlast = dlsa%pdest
      goto 900
   endif
endif
m = m + 1
if ( m/=dls1%maxcor ) then
   if ( m<2 .or. del<=2.0D0*delp ) then
      delp = del
      call f(Neq,dls1%tn,Y,Savf)
      dls1%nfe = dls1%nfe + 1
      goto 600
   endif
endif
!-----------------------------------------------------------------------
!  The corrector iteration failed to converge.
!  If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
!  the next try.  Otherwise the YH array is retracted to its values
!  before prediction, and H is reduced, if possible.  If H cannot be
!  reduced or MXNCF failures have occurred, exit with KFLAG = -2.
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
! DIR$ IVDEP
   do i = i1, dls1%nqnyh
      Yh1(i) = Yh1(i) - Yh1(i+Nyh)
   enddo
enddo
if ( dls1%ierpj<0 .or. dls1%iersl<0 ) then
   dls1%kflag = -3
   goto 1400
elseif ( abs(dls1%h)<=dls1%hmin*1.00001D0 ) then
   dls1%kflag = -2
   goto 1400
elseif ( ncf==dls1%mxncf ) then
   dls1%kflag = -2
   goto 1400
else
   rh = 0.25D0
   dls1%ipup = dls1%miter
   iredo = 1
   rh = max(rh,dls1%hmin/abs(dls1%h))
   goto 300
endif
!-----------------------------------------------------------------------
!  The corrector has converged.  JCUR is set to 0
!  to signal that the Jacobian involved may need updating later.
!  The local error test is made and control passes to statement 500
!  if it fails.
!-----------------------------------------------------------------------
 900  continue
dls1%jcur = 0
if ( m==0 ) dsm = del/dls1%tesco(2,dls1%nq)
if ( m>0 ) dsm = dmnorm(dls1%n,Acor,Ewt)/dls1%tesco(2,dls1%nq)
if ( dsm>1.0D0 ) then
!-----------------------------------------------------------------------
!  The error test failed.  KFLAG keeps track of multiple failures.
!  Restore TN and the YH array to their previous values, and prepare
!  to try the step again.  Compute the optimum step size for this or
!  one lower order.  After 2 or more failures, H is forced to decrease
!  by a factor of 0.2 or less.
!-----------------------------------------------------------------------
   dls1%kflag = dls1%kflag - 1
   dls1%tn = told
   i1 = dls1%nqnyh + 1
   do jb = 1, dls1%nq
      i1 = i1 - Nyh
! DIR$ IVDEP
      do i = i1, dls1%nqnyh
         Yh1(i) = Yh1(i) - Yh1(i+Nyh)
      enddo
   enddo
   dls1%rmax = 2.0D0
   if ( abs(dls1%h)<=dls1%hmin*1.00001D0 ) then
!-----------------------------------------------------------------------
!  All returns are made through this section.  H is saved in HOLD
!  to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
      dls1%kflag = -1
      goto 1400
   elseif ( dls1%kflag<=-3 ) then
!-----------------------------------------------------------------------
!  Control reaches this section if 3 or more failures have occured.
!  If 10 failures have occurred, exit with KFLAG = -1.
!  It is assumed that the derivatives that have accumulated in the
!  YH array have errors of the wrong order.  Hence the first
!  derivative is recomputed, and the order is set to 1.  Then
!  H is reduced by a factor of 10, and the step is retried,
!  until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
      if ( dls1%kflag==-10 ) then
         dls1%kflag = -1
         goto 1400
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
   endif
else
!-----------------------------------------------------------------------
!  After a successful step, update the YH array.
!  Decrease ICOUNT by 1, and if it is -1, consider switching methods.
!  If a method switch is made, reset various parameters,
!  rescale the YH array, and exit.  If there is no switch,
!  consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
!  If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
!  use in a possible order increase on the next step.
!  If a change in H is considered, an increase or decrease in order
!  by one is considered also.  A change in H is made only if it is by a
!  factor of at least 1.1.  If not, IALTH is set to 3 to prevent
!  testing for that many steps.
!-----------------------------------------------------------------------
   dls1%kflag = 0
   iredo = 0
   dls1%nst = dls1%nst + 1
   dls1%hu = dls1%h
   dls1%nqu = dls1%nq
   dlsa%mused = dls1%meth
   do j = 1, dls1%l
      do i = 1, dls1%n
         Yh(i,j) = Yh(i,j) + dls1%el(j)*Acor(i)
      enddo
   enddo
   dlsa%icount = dlsa%icount - 1
   if ( dlsa%icount<0 ) then
      if ( dls1%meth==2 ) then
!-----------------------------------------------------------------------
!  We are currently using a BDF method.  Consider switching to Adams.
!  Compute the step size we could have (ideally) used on this step,
!  with the current (BDF) method, and also that for the Adams.
!  If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
!  Compare the two step sizes to decide whether to switch.
!  The step size advantage must be at least 5/RATIO = 1 to switch.
!  If the step size for Adams would be so small as to cause
!  roundoff pollution, we stay with BDF.
!-----------------------------------------------------------------------
         exsm = 1.0D0/dls1%l
         if ( dlsa%mxordn>=dls1%nq ) then
            dm1 = dsm*(dlsa%cm2(dls1%nq)/dlsa%cm1(dls1%nq))
            rh1 = 1.0D0/(1.2D0*dm1**exsm+0.0000012D0)
            nqm1 = dls1%nq
            exm1 = exsm
         else
            nqm1 = dlsa%mxordn
            lm1 = dlsa%mxordn + 1
            exm1 = 1.0D0/lm1
            lm1p1 = lm1 + 1
            dm1 = dmnorm(dls1%n,Yh(1,lm1p1),Ewt)/dlsa%cm1(dlsa%mxordn)
            rh1 = 1.0D0/(1.2D0*dm1**exm1+0.0000012D0)
         endif
         rh1it = 2.0D0*rh1
         pdh = dlsa%pdnorm*abs(dls1%h)
         if ( pdh*rh1>0.00001D0 ) rh1it = sm1(nqm1)/pdh
         rh1 = min(rh1,rh1it)
         rh2 = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
         if ( rh1*dlsa%ratio>=5.0D0*rh2 ) then
            alpha = max(0.001D0,rh1)
            dm1 = (alpha**exm1)*dm1
            if ( dm1>1000.0D0*dls1%uround*pnorm ) then
!  The switch test passed.  Reset relevant quantities for Adams. --------
               rh = rh1
               dlsa%icount = 20
               dls1%meth = 1
               dls1%miter = 0
               dlsa%pdlast = 0.0D0
               dls1%nq = nqm1
               dls1%l = dls1%nq + 1
               rh = max(rh,dls1%hmin/abs(dls1%h))
               goto 300
            endif
         endif
!-----------------------------------------------------------------------
!  We are currently using an Adams method.  Consider switching to BDF.
!  If the current order is greater than 5, assume the problem is
!  not stiff, and skip this section.
!  If the Lipschitz constant and error estimate are not polluted
!  by roundoff, go to 470 and perform the usual test.
!  Otherwise, switch to the BDF methods if the last step was
!  restricted to insure stability (dlsa%irflag = 1), and stay with Adams
!  method if not.  When switching to BDF with polluted error estimates,
!  in the absence of other information, double the step size.
!
!  When the estimates are OK, we make the usual test by computing
!  the step size we could have (ideally) used on this step,
!  with the current (Adams) method, and also that for the BDF.
!  If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
!  Compare the two step sizes to decide whether to switch.
!  The step size advantage must be at least RATIO = 5 to switch.
!-----------------------------------------------------------------------
      elseif ( dls1%nq<=5 ) then
         if ( dsm>100.0D0*pnorm*dls1%uround .and. dlsa%pdest/=0.0D0 ) then
            exsm = 1.0D0/dls1%l
            rh1 = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
            rh1it = 2.0D0*rh1
            pdh = dlsa%pdlast*abs(dls1%h)
            if ( pdh*rh1>0.00001D0 ) rh1it = sm1(dls1%nq)/pdh
            rh1 = min(rh1,rh1it)
            if ( dls1%nq<=dlsa%mxords ) then
               dm2 = dsm*(dlsa%cm1(dls1%nq)/dlsa%cm2(dls1%nq))
               rh2 = 1.0D0/(1.2D0*dm2**exsm+0.0000012D0)
               nqm2 = dls1%nq
            else
               nqm2 = dlsa%mxords
               lm2 = dlsa%mxords + 1
               exm2 = 1.0D0/lm2
               lm2p1 = lm2 + 1
               dm2 = dmnorm(dls1%n,Yh(1,lm2p1),Ewt)/dlsa%cm2(dlsa%mxords)
               rh2 = 1.0D0/(1.2D0*dm2**exm2+0.0000012D0)
            endif
            if ( rh2<dlsa%ratio*rh1 ) goto 950
         else
            if ( dlsa%irflag==0 ) goto 950
            rh2 = 2.0D0
            nqm2 = min(dls1%nq,dlsa%mxords)
         endif
!  THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
         rh = rh2
         dlsa%icount = 20
         dls1%meth = 2
         dls1%miter = dlsa%jtyp
         dlsa%pdlast = 0.0D0
         dls1%nq = nqm2
         dls1%l = dls1%nq + 1
         rh = max(rh,dls1%hmin/abs(dls1%h))
         goto 300
      endif
   endif
!
!  No method switch is being made.  Do the usual step/order selection. --
 950  continue
   dls1%ialth = dls1%ialth - 1
   if ( dls1%ialth==0 ) then
!-----------------------------------------------------------------------
!  Regardless of the success or failure of the step, factors
!  RHDN, RHSM, and RHUP are computed, by which H could be multiplied
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
         dup = dmnorm(dls1%n,Savf,Ewt)/dls1%tesco(3,dls1%nq)
         exup = 1.0D0/(dls1%l+1)
         rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
      endif
   else
      if ( dls1%ialth<=1 ) then
         if ( dls1%l/=dls1%lmax ) then
            do i = 1, dls1%n
               Yh(i,dls1%lmax) = Acor(i)
            enddo
         endif
      endif
      goto 1300
   endif
endif
exsm = 1.0D0/dls1%l
rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
rhdn = 0.0D0
if ( dls1%nq/=1 ) then
   ddn = dmnorm(dls1%n,Yh(1,dls1%l),Ewt)/dls1%tesco(1,dls1%nq)
   exdn = 1.0D0/dls1%nq
   rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
endif
!  If METH = 1, limit RH according to the stability region also. --------
if ( dls1%meth/=2 ) then
   pdh = max(abs(dls1%h)*dlsa%pdlast,0.000001D0)
   if ( dls1%l<dls1%lmax ) rhup = min(rhup,sm1(dls1%l)/pdh)
   rhsm = min(rhsm,sm1(dls1%nq)/pdh)
   if ( dls1%nq>1 ) rhdn = min(rhdn,sm1(dls1%nq-1)/pdh)
   dlsa%pdest = 0.0D0
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
      goto 1300
   else
      r = dls1%el(dls1%l)/dls1%l
      do i = 1, dls1%n
         Yh(i,newq+1) = Acor(i)*r
      enddo
      goto 1200
   endif
endif
newq = dls1%nq - 1
rh = rhdn
if ( dls1%kflag<0 .and. rh>1.0D0 ) rh = 1.0D0
!  If METH = 1 and H is restricted by stability, bypass 10 percent test.
 1000 continue
if ( dls1%meth/=2 ) then
   if ( rh*pdh*1.00001D0>=sm1(newq) ) goto 1100
endif
if ( dls1%kflag==0 .and. rh<1.1D0 ) then
   dls1%ialth = 3
   goto 1300
endif
 1100 continue
if ( dls1%kflag<=-2 ) rh = min(rh,0.2D0)
!-----------------------------------------------------------------------
!  If there is a change of order, reset NQ, L, and the coefficients.
!  In any case H is reset according to RH and the YH array is rescaled.
!  Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
if ( newq==dls1%nq ) then
   rh = max(rh,dls1%hmin/abs(dls1%h))
   goto 300
endif
 1200 continue
dls1%nq = newq
dls1%l = dls1%nq + 1
iret = 2
goto 100
 1300 continue
r = 1.0D0/dls1%tesco(2,dls1%nqu)
do i = 1, dls1%n
   Acor(i) = Acor(i)*r
enddo
 1400 continue
dls1%hold = dls1%h
dls1%jstart = 1
end subroutine dstoda
