!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DSTODI performs one step of the integration of an initial value
!! problem for a system of Ordinary Differential Equations.
!!
!! Note: DSTODI is independent of the value of the iteration method
!! indicator MITER, and hence is independent
!! of the type of chord method used, or the Jacobian structure.
!!
!! Communication with DSTODI is done with the following variables:
!!
!! NEQ
!!
!! : integer array containing problem size in NEQ(1), and
!! passed as the NEQ argument in all calls to RES, ADDA,
!! and JAC.
!!
!! Y
!!
!! : an array of length .ge. N used as the Y argument in
!! all calls to RES, JAC, and ADDA.
!!
!! NEQ
!!
!! : integer array containing problem size in NEQ(1), and
!! passed as the NEQ argument in all calls tO RES, G, ADDA,
!! and JAC.
!!
!! YH
!!
!! : an NYH by LMAX array containing the dependent variables
!! and their approximate scaled derivatives, where
!!
!! LMAX
!!
!! : MAXORD + 1.  YH(i,j+1) contains the approximate
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
!! : an array of working storage, of length N. also used for
!! input of YH(*,MAXORD+2) when JSTART = -1 and MAXORD is less
!! than the current order NQ.
!! Same as YDOTI in the driver.
!!
!! SAVR
!!
!! : an array of working storage, of length N.
!!
!! ACOR
!!
!! : a work array of length N used for the accumulated
!! corrections. On a succesful return, ACOR(i) contains
!! the estimated one-step local error in y(i).
!!
!! WM,IWM
!!
!! : real and integer work arrays associated with matrix
!! operations in chord iteration.
!!
!! PJAC
!!
!! : name of routine to evaluate and preprocess Jacobian matrix.
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
!!
!!               0  perform the first step.
!!           .gt.0  take a new step continuing from the last.
!!              -1  take the next step with a new value of H, MAXORD,
!!                  N, METH, MITER, and/or matrix parameters.
!!              -2  take the next step with a new value of H,
!!                  but with other inputs unchanged.
!!                  On return, JSTART is set to 1 to facilitate continuation.
!!
!! KFLAG
!!
!! : a completion code with the following meanings:
!!
!!               0  the step was succesful.
!!              -1  the requested error could not be achieved.
!!              -2  corrector convergence could not be achieved.
!!              -3  RES ordered immediate return.
!!              -4  error condition from RES could not be avoided.
!!              -5  fatal error in PJAC or SLVS.
!!
!! A return with KFLAG = -1, -2, or -4 means either
!!
!! ABS(H)
!!
!! : HMIN or 10 consecutive failures occurred.
!!
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
!! : maximum number of steps between PJAC calls.
!!
!! MXNCF
!!
!! : maximum number of convergence failures allowed.
!! METH/MITER = the method flags.  See description in driver.
!!
!! N
!!
!! : the number of first-order differential equations.
!-----------------------------------------------------------------------
subroutine dstodi(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Savr,Acor,Wm,Iwm,res,adda,jac,pjac,slvs)
Use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)

integer, dimension(*) :: Neq
real(kind=dp), dimension(*) :: Y
integer :: Nyh
real(kind=dp), intent(inout), dimension(Nyh,*) :: Yh
real(kind=dp), intent(inout), dimension(*) :: Yh1
real(kind=dp), dimension(*) :: Ewt
real(kind=dp), dimension(*) :: Savf
real(kind=dp), dimension(*) :: Savr
real(kind=dp), intent(inout), dimension(*) :: Acor
real(kind=dp), dimension(*) :: Wm
integer, dimension(*) :: Iwm
external res
external adda
external jac
external pjac
external slvs

real(kind=dp) :: dcon, ddn, del, delp, dsm, dup, el1h, eljh, exdn, exsm, exup, r, rh, rhdn, rhsm, rhup, told
integer :: i, i1, iredo, ires, iret, j, jb, kgo, m, ncf, newq

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
   !  If the caller has changed MAXORD to a value less than the current
   !  order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
   !  If H is to be changed, YH must be rescaled.
   !  If H or METH is being changed, IALTH is reset to L = NQ + 1
   !  to prevent further changes in H for that many steps.
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
   !  initialized.  RMAX is the maximum ratio by which H can be increased
   !  in a single step.  It is initially 1.E4 to compensate for the small
   !  initial H, but then is normally equal to 10.  If a failure
   !  occurs (in corrector convergence or error test), RMAX is set at 2
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
   !  current METH.  Then the EL vector and related constants are reset
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
!  multiplying the YH array by the Pascal triangle matrix.
!  RC is the ratio of new to old values of the coefficient  H*EL(1).
!  When RC differs from 1 by more than CCMAX, IPUP is set to MITER
!  to force PJAC to be called.
!  In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
 400  continue
if ( abs(dls1%rc-1.0D0)>dls1%ccmax ) dls1%ipup = dls1%miter
if ( dls1%nst>=dls1%nslp+dls1%msbp ) dls1%ipup = dls1%miter
dls1%tn = dls1%tn + dls1%h
i1 = dls1%nqnyh + 1
do jb = 1, dls1%nq
   i1 = i1 - Nyh
   do i = i1, dls1%nqnyh
      Yh1(i) = Yh1(i) + Yh1(i+Nyh)
   enddo
enddo
!-----------------------------------------------------------------------
!  Up to MAXCOR corrector iterations are taken.  A convergence test is
!  made on the RMS-norm of each correction, weighted by H and the
!  error weight vector EWT.  The sum of the corrections is accumulated
!  in ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 500  continue
m = 0
do i = 1, dls1%n
   Savf(i) = Yh(i,2)/dls1%h
   Y(i) = Yh(i,1)
enddo
if ( dls1%ipup>0 ) then
   !-----------------------------------------------------------------------
   !  If indicated, the matrix P = A - H*EL(1)*dr/dy is reevaluated and
   !  preprocessed before starting the corrector iteration.  IPUP is set
   !  to 0 as an indicator that this has been done.
   !-----------------------------------------------------------------------
   call pjac(Neq,Y,Yh,Nyh,Ewt,Acor,Savr,Savf,Wm,Iwm,res,jac,adda)
   dls1%ipup = 0
   dls1%rc = 1.0D0
   dls1%nslp = dls1%nst
   dls1%crate = 0.7D0
   if ( dls1%ierpj==0 ) goto 600
   if ( dls1%ierpj<0 ) goto 800
   ires = dls1%ierpj
   select case (ires)
   case (1,3)
      goto 700
   case (2)
      goto 800
   case default
   endselect
endif
!  Get residual at predicted values, if not already done in PJAC. -------
ires = 1
call res(Neq,dls1%tn,Y,Savf,Savr,ires)
dls1%nfe = dls1%nfe + 1
kgo = abs(ires)
select case (kgo)
case (2)
   goto 800
case (3)
   goto 700
case default
endselect
 600  continue
do i = 1, dls1%n
   Acor(i) = 0.0D0
enddo
do
   !-----------------------------------------------------------------------
   !  Solve the linear system with the current residual as
   !  right-hand side and P as coefficient matrix.
   !-----------------------------------------------------------------------
   call slvs(Wm,Iwm,Savr,Savf)
   if ( dls1%iersl<0 ) exit
   if ( dls1%iersl<=0 ) then
      el1h = dls1%el(1)*dls1%h
      del = dvnorm(dls1%n,Savr,Ewt)*abs(dls1%h)
      do i = 1, dls1%n
         Acor(i) = Acor(i) + Savr(i)
         Savf(i) = Acor(i) + Yh(i,2)/dls1%h
         Y(i) = Yh(i,1) + el1h*Acor(i)
      enddo
      !-----------------------------------------------------------------------
      !  Test for convergence.  If M .gt. 0, an estimate of the convergence
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
         if ( m>0 ) dsm = abs(dls1%h)*dvnorm(dls1%n,Acor,Ewt)/dls1%tesco(2,dls1%nq)
         if ( dsm>1.0D0 ) then
            !-----------------------------------------------------------------------
            !  The error test failed.  KFLAG keeps track of multiple failures.
            !  restore TN and the YH array to their previous values, and prepare
            !  to try the step again.  Compute the optimum step size for this or
            !  one lower order.  After 2 or more failures, H is forced to decrease
            !  by a factor of 0.1 or less.
            !-----------------------------------------------------------------------
            dls1%kflag = dls1%kflag - 1
            dls1%tn = told
            i1 = dls1%nqnyh + 1
            do jb = 1, dls1%nq
               i1 = i1 - Nyh
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
               goto 1300
            elseif ( dls1%kflag<=-7 ) then
               dls1%kflag = -1
               goto 1300
            else
               iredo = 2
               rhup = 0.0D0
               goto 900
            endif
         else
            !-----------------------------------------------------------------------
            !  After a successful step, update the YH array.
            !  Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
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
            do j = 1, dls1%l
               eljh = dls1%el(j)*dls1%h
               do i = 1, dls1%n
                  Yh(i,j) = Yh(i,j) + eljh*Acor(i)
               enddo
            enddo
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
                  dup = abs(dls1%h)*dvnorm(dls1%n,Savf,Ewt)/dls1%tesco(3,dls1%nq)
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
               ires = 1
               call res(Neq,dls1%tn,Y,Savf,Savr,ires)
               dls1%nfe = dls1%nfe + 1
               kgo = abs(ires)
               select case (kgo)
               case (1)
                  cycle
               case (2)
                  goto 800
               case default
               endselect
            endif
         endif
      endif
   endif
   !-----------------------------------------------------------------------
   !  The correctors failed to converge, or RES has returned abnormally.
   !  on a convergence failure, if the Jacobian is out of date, PJAC is
   !  called for the next try.  Otherwise the YH array is retracted to its
   !  values before prediction, and H is reduced, if possible.
   !  take an error exit if IRES = 2, or H cannot be reduced, or MXNCF
   !  failures have occurred, or a fatal error occurred in PJAC or SLVS.
   !-----------------------------------------------------------------------
   dls1%icf = 1
   if ( dls1%jcur==1 ) exit
   dls1%ipup = dls1%miter
   goto 500
enddo
 700  continue
dls1%icf = 2
ncf = ncf + 1
dls1%rmax = 2.0D0
 800  continue
dls1%tn = told
i1 = dls1%nqnyh + 1
do jb = 1, dls1%nq
   i1 = i1 - Nyh
   do i = i1, dls1%nqnyh
      Yh1(i) = Yh1(i) - Yh1(i+Nyh)
   enddo
enddo
if ( ires==2 ) then
   dls1%kflag = -1 - ires
elseif ( dls1%ierpj<0 .or. dls1%iersl<0 ) then
   dls1%kflag = -5
else
   if ( abs(dls1%h)>dls1%hmin*1.00001D0 ) then
      if ( ncf/=dls1%mxncf ) then
         rh = 0.25D0
         dls1%ipup = dls1%miter
         iredo = 1
         rh = max(rh,dls1%hmin/abs(dls1%h))
         goto 300
      endif
   endif
   if ( ires==3 ) then
      dls1%kflag = -1 - ires
   else
      dls1%kflag = -2
   endif
endif
goto 1300
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
      r = dls1%h*dls1%el(dls1%l)/dls1%l
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
   if ( dls1%kflag<=-2 ) rh = min(rh,0.1D0)
   !-----------------------------------------------------------------------
   !  If there is a change of order, reset NQ, L, and the coefficients.
   !  In any case H is reset according to RH and the YH array is rescaled.
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
!-----------------------------------------------------------------------------------------------------------------------------------
 1200 continue
r = dls1%h/dls1%tesco(2,dls1%nqu)
do i = 1, dls1%n
   Acor(i) = Acor(i)*r
enddo
 1300 continue
dls1%hold = dls1%h
dls1%jstart = 1
end subroutine dstodi
