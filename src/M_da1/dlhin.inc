!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### Name
!!   dlhin(3f) - [M_odepack] compute step size H0 to be attempted on
!!   the first step, when the user supplied value is absent
!!
!!### Synopsis
!!
!!        subroutine dlhin(Neq,N,T0,Y0,Ydot,f,Tout,Uround,Ewt,Itol,Atol,Y,Temp,H0,Niter,Ier)
!!
!!        integer                     :: Neq(*)
!!        integer                     :: N
!!        real(kind=dp),intent(in)    :: T0
!!        real(kind=dp)               :: Y0(*)
!!        real(kind=dp),intent(in)    :: Ydot(*)
!!        external                    :: f
!!        real(kind=dp),intent(in)    :: Tout
!!        real(kind=dp),intent(in)    :: Uround
!!        real(kind=dp)               :: Ewt(*)
!!        integer,intent(in)          :: Itol
!!        real(kind=dp),intent(in)    :: Atol(*)
!!        real(kind=dp)               :: Y(*)
!!        real(kind=dp)               :: Temp(*)
!!        real(kind=dp),intent(inout) :: H0
!!        integer, intent(out)        :: Niter
!!        integer, intent(out)        :: Ier
!!
!!### Description
!!
!! This routine computes the step size, H0, to be attempted on the
!! first step, when the user has not supplied a value for this.
!!
!! First we check that TOUT - T0 differs significantly from zero.  Then
!! an iteration is done to approximate the initial second derivative
!! and this is used to define H from WRMS-norm(H**2 * yddot / 2) = 1.
!! A bias factor of 1/2 is applied to the resulting h.
!! The sign of H0 is inferred from the initial values of TOUT and T0.
!!
!! Communication with DLHIN is done with the following variables:
!!
!! Call DLHIN(NEQ,N,T0,Y0,YDOT,F,TOUT,UROUND,EWT,ITOL,ATOL,Y,TEMP &
!! & H0,NITER,IER
!!
!! Subroutines called by DLHIN: F
!! Function routines called by DLHIN: DVNORM
!!
!!### INPUT OPTIONS
!!
!! NEQ
!!
!! : NEQ array of solver, passed to F.
!!
!! N
!!
!! : size of ODE system, input.
!!
!! T0
!!
!! : initial value of independent variable, input.
!!
!! Y0
!!
!! : vector of initial conditions, input.
!!
!! YDOT
!!
!! : vector of initial first derivatives, input.
!!
!! F
!!
!! : name of subroutine for right-hand side f(t,y), input.
!!
!! TOUT
!!
!! : first output value of independent variable
!!
!! UROUND
!!
!! : machine unit roundoff
!!
!! EWT, ITOL, ATOL
!!
!! : error weights and tolerance parameters
!! as described in the driver routine, input.
!!
!! Y, TEMP
!!
!! : work arrays of length N.
!!
!!### RETURNS
!!
!! H0
!!
!! : step size to be attempted, output.
!!
!! NITER
!!
!! : number of iterations (and of f evaluations) to compute H0,
!! output.
!!
!! IER
!!
!! : the error flag, returned with the value
!!       IER = 0  if no trouble occurred, or
!!       IER = -1 if TOUT and t0 are considered too close to proceed.
!-----------------------------------------------------------------------
subroutine dlhin(Neq,N,T0,Y0,Ydot,f,Tout,Uround,Ewt,Itol,Atol,Y,Temp,H0,Niter,Ier)

real(kind=dp),parameter :: HALF = 0.50D0 , TWO = 2.0D0 , HUN = 100.0D0 , PT1 = 0.1D0
!
integer                     :: Neq(*)
integer                     :: N
real(kind=dp),intent(in)    :: T0
real(kind=dp)               :: Y0(*)
real(kind=dp),intent(in)    :: Ydot(*)
external                    :: f
real(kind=dp),intent(in)    :: Tout
real(kind=dp),intent(in)    :: Uround
real(kind=dp)               :: Ewt(*)
integer,intent(in)          :: Itol
real(kind=dp),intent(in)    :: Atol(*)
real(kind=dp)               :: Y(*)
real(kind=dp)               :: Temp(*)
real(kind=dp),intent(inout) :: H0
integer, intent(out)        :: Niter
integer, intent(out)        :: Ier

real(kind=dp) :: afi, atoli, delyi, hg, hlb, hnew, hrat, hub, t1, tdist, tround, yddnrm
integer       :: i, iter

   Niter = 0
   tdist = abs(Tout-T0)
   tround = Uround*max(abs(T0),abs(Tout))
   if ( tdist < TWO*tround ) then
      !  Error return for TOUT - T0 too small. --------------------------------
      Ier = -1
      return
   else
      !
      !  Set a lower bound on H based on the roundoff level in T0 and TOUT. ---
      hlb = HUN*tround
      !  Set an upper bound on H based on TOUT-T0 and the initial Y and YDOT. -
      hub = PT1*tdist
      atoli = Atol(1)
      do i = 1 , N
         if ( Itol == 2 .or. Itol == 4 ) atoli = Atol(i)
         delyi = PT1*abs(Y0(i)) + atoli
         afi = abs(Ydot(i))
         if ( afi*hub > delyi ) hub = delyi/afi
      enddo
      !
      !  Set initial guess for H as geometric mean of upper and lower bounds. -
      iter = 0
      hg = sqrt(hlb*hub)
      !  If the bounds have crossed, exit with the mean value. ----------------
      if ( hub<hlb ) then
         H0 = hg
      else
         do
            !
            !  Looping point for iteration. -----------------------------------------
            !  Estimate the second derivative as a difference quotient in f. --------
            t1 = T0 + hg
            do i = 1 , N
               Y(i) = Y0(i) + hg*Ydot(i)
            enddo
            call f(Neq,t1,Y,Temp)
            do i = 1 , N
               Temp(i) = (Temp(i)-Ydot(i))/hg
            enddo
            yddnrm = dvnorm(N,Temp,Ewt)
            !  Get the corresponding new value of H. --------------------------------
            if ( yddnrm*hub*hub>TWO ) then
               hnew = sqrt(TWO/yddnrm)
            else
               hnew = sqrt(hg*hub)
            endif
            iter = iter + 1
            !-----------------------------------------------------------------------
            !  Test the stopping conditions.
            !  Stop if the new and previous H values differ by a factor of .lt. 2.
            !  Stop if four iterations have been done.  Also, stop with previous H
            !  if hnew/hg .gt. 2 after first iteration, as this probably means that
            !  the second derivative value is bad because of cancellation error.
            !-----------------------------------------------------------------------
            if ( iter>=4 ) exit
            hrat = hnew/hg
            if ( (hrat>HALF) .and. (hrat<TWO) ) exit
            if ( (iter>=2) .and. (hnew>TWO*hg) ) then
               hnew = hg
               exit
            endif
            hg = hnew
         enddo
         !
         !  Iteration done.  Apply bounds, bias factor, and sign. ----------------
         H0 = hnew*HALF
         if ( H0<hlb ) H0 = hlb
         if ( H0>hub ) H0 = hub
      endif
   endif
   H0 = sign(H0,Tout-T0)
   !  Restore Y array from Y0, then exit. ----------------------------------
   Y(1:N) = Y0(1:N)
   Niter = iter
   Ier = 0
end subroutine dlhin
