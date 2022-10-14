!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This subroutine finds the leftmost root of a set of arbitrary
!! functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
!! of odd multiplicity (i.e. changes of sign of the gi) are found.
!! Here the sign of X1 - X0 is arbitrary, but is constant for a given
!! problem, and -leftmost- means nearest to X0.
!! The values of the vector-valued function g(x) = (gi, i=1...NG)
!! are communicated through the call sequence of DROOTS.
!! The method used is the Illinois algorithm.
!!
!! #### Reference:
!! Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
!! Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
!! February 1980.
!!
!!### Description of parameters.
!!
!! NG
!!
!! : number of functions gi, or the number of components of
!! the vector valued function g(x).  Input only.
!!
!! HMIN
!!
!! : resolution parameter in X.  Input only.  When a root is
!! found, it is located only to within an error of HMIN in X.
!! Typically, HMIN should be set to something on the order of
!!               100 * UROUND * MAX(ABS(X0),ABS(X1)),
!! where UROUND is the unit roundoff of the machine.
!!
!! JFLAG
!!
!! : integer flag for input and output communication.
!!
!! On input, set JFLAG = 0 on the first call for the problem,
!! and leave it unchanged until the problem is completed.
!! (The problem is completed when JFLAG .ge. 2 on return.)
!!
!!### On output
!! JFLAG
!!
!! :JFLAG has the following values and meanings:
!!
!!          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
!!                    and call DROOTS again.
!!          JFLAG = 2 means a root has been found.  The root is
!!                    at X, and GX contains g(X).  (Actually, X is the
!!                    rightmost approximation to the root on an interval
!!                    (X0,X1) of size HMIN or less.)
!!          JFLAG = 3 means X = X1 is a root, with one or more of the gi
!!                    being zero at X1 and no sign changes in (X0,X1).
!!                    GX contains g(X) on output.
!!          JFLAG = 4 means no roots (of odd multiplicity) were
!!                    found in (X0,X1) (no sign changes).
!!
!! X0,X1
!!
!! : endpoints of the interval where roots are sought.
!! X1 and X0 are input when JFLAG = 0 (first call), and
!! must be left unchanged between calls until the problem is
!! completed.  X0 and X1 must be distinct, but X1 - X0 may be
!! of either sign.  However, the notion of -left- and -right-
!! will be used to mean nearer to X0 or X1, respectively.
!! When JFLAG .ge. 2 on return, X0 and X1 are output, and
!! are the endpoints of the relevant interval.
!!
!! G0,G1
!!
!! : arrays of length NG containing the vectors g(X0) and g(X1),
!! respectively.  When JFLAG = 0, G0 and G1 are input and
!! none of the G0(i) should be zero.
!! When JFLAG .ge. 2 on return, G0 and G1 are output.
!!
!! GX
!!
!! : array of length NG containing g(X).  GX is input
!! when JFLAG = 1, and output when JFLAG .ge. 2.
!!
!! X
!!
!! : independent variable value.  Output only.
!! When JFLAG = 1 on output, X is the point at which g(x)
!! is to be evaluated and loaded into GX.
!!
!!          When JFLAG = 2 or 3, X is the root.
!!
!!          When JFLAG = 4, X is the right endpoint of the interval, X1.
!!
!! JROOT
!!
!! : integer array of length NG.  Output only.
!! When JFLAG = 2 or 3, JROOT indicates which components
!! of g(x) have a root at X.  JROOT(i) is 1 if the i-th
!! component has a root, and JROOT(i) = 0 otherwise.
!-----------------------------------------------------------------------
subroutine droots(Ng,Hmin,Jflag,X0,X1,G0,G1,Gx,X,Jroot)
!
integer                      :: Ng
real(kind=dp),intent(in)     :: Hmin
integer,intent(inout)        :: Jflag
real(kind=dp),intent(inout)  :: X0
real(kind=dp),intent(inout)  :: X1
real(kind=dp)                :: G0(Ng)
real(kind=dp)                :: G1(Ng)
real(kind=dp)                :: Gx(Ng)
real(kind=dp),intent(out)    :: X
integer,intent(out)          :: Jroot(Ng)
!
real(kind=dp) :: fracint , fracsub , t2 , tmax
integer :: i , imxold , nxlast
logical :: sgnchg , xroot , zroot
!
real(kind=dp) , parameter :: five=5.0d0 , half=0.5d0 , tenth=0.1d0 , zero=0.0d0
!
   if ( Jflag==1 ) then
   !  Check to see in which interval g changes sign. -----------------------
      imxold = dlsr%imax
      dlsr%imax = 0
      tmax = zero
      zroot = .false.
      do i = 1 , Ng
         if ( abs(Gx(i))<=zero ) then
            zroot = .true.
         !  Neither G0(i) nor GX(i) can be zero at this point. -------------------
         elseif ( sign(1.0D0,G0(i))/=sign(1.0D0,Gx(i)) ) then
            t2 = abs(Gx(i)/(Gx(i)-G0(i)))
            if ( t2>tmax ) then
               tmax = t2
               dlsr%imax = i
            endif
         endif
      enddo
      if ( dlsr%imax>0 ) then
         sgnchg = .true.
      else
         sgnchg = .false.
         dlsr%imax = imxold
      endif
      nxlast = dlsr%last
      if ( sgnchg ) then
   !  Sign change between X0 and X2, so replace X1 with X2. ----------------
         X1 = dlsr%x2
         !X!call dcopy(Ng,Gx,1,G1,1)
         G1(1:Ng) = Gx(1:Ng)!X!
         dlsr%last = 1
         xroot = .false.
      elseif ( .not.zroot ) then
   !  No sign change between X0 and X2.  Replace X0 with X2. ---------------
         !X!call dcopy(Ng,Gx,1,G0,1)
         G0(1:Ng) = Gx(1:Ng)!X!
         X0 = dlsr%x2
         dlsr%last = 0
         xroot = .false.
      else
   !  Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
         X1 = dlsr%x2
         !X!call dcopy(Ng,Gx,1,G1,1)
         G1(1:Ng) = Gx(1:Ng)!X!
         xroot = .true.
      endif
      if ( abs(X1-X0)<=Hmin ) xroot = .true.
   else
   !  JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
      dlsr%imax = 0
      tmax = zero
      zroot = .false.
      do i = 1 , Ng
         if ( abs(G1(i))<=zero ) then
            zroot = .true.
   !  At this point, G0(i) has been checked and cannot be zero. ------------
         elseif ( sign(1.0D0,G0(i))/=sign(1.0D0,G1(i)) ) then
            t2 = abs(G1(i)/(G1(i)-G0(i)))
            if ( t2>tmax ) then
               tmax = t2
               dlsr%imax = i
            endif
         endif
      enddo
      if ( dlsr%imax>0 ) then
         sgnchg = .true.
      else
         sgnchg = .false.
      endif
      if ( .not.sgnchg ) then
   !
   !  No sign change in the interval.  Check for zero at right endpoint. ---
         if ( zroot ) then
            !
            !  Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
            X = X1
            !X!call dcopy(Ng,G1,1,Gx,1)
            Gx(1:Ng) = G1(1:Ng)!X!
            do i = 1 , Ng
               Jroot(i) = 0
               if ( abs(G1(i))<=zero ) Jroot(i) = 1
            enddo
            Jflag = 3
            return
         endif
   !
   !  No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
         !X!call dcopy(Ng,G1,1,Gx,1)
         Gx(1:Ng) = G1(1:Ng)!X!
         X = X1
         Jflag = 4
         return
      else
   !  There is a sign change.  Find the first root in the interval. --------
         xroot = .false.
         nxlast = 0
         dlsr%last = 1
      endif
   endif
   !
   !  Repeat until the first root in the interval is found.  Loop point. ---
   if ( xroot ) then
   !
   !  Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
      Jflag = 2
      X = X1
      !X!call dcopy(Ng,G1,1,Gx,1)
      Gx(1:Ng) = G1(1:Ng)!X!
      do i = 1 , Ng
         Jroot(i) = 0
         if ( abs(G1(i))>zero ) then
            if ( sign(1.0D0,G0(i))/=sign(1.0D0,G1(i)) ) Jroot(i) = 1
         else
            Jroot(i) = 1
         endif
      enddo
   else
      if ( nxlast/=dlsr%last ) then
         dlsr%alpha = 1.0D0
      elseif ( dlsr%last==0 ) then
         dlsr%alpha = 2.0D0*dlsr%alpha
      else
         dlsr%alpha = 0.5D0*dlsr%alpha
      endif
      dlsr%x2 = X1 - (X1-X0)*G1(dlsr%imax)/(G1(dlsr%imax)-dlsr%alpha*G0(dlsr%imax))
   !  If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
   !  distance that is between 0.1 and 0.5. --------------------------------
      if ( abs(dlsr%x2-X0)<half*Hmin ) then
         fracint = abs(X1-X0)/Hmin
         fracsub = tenth
         if ( fracint<=five ) fracsub = half/fracint
         dlsr%x2 = X0 + fracsub*(X1-X0)
      endif
      if ( abs(X1-dlsr%x2)<half*Hmin ) then
         fracint = abs(X1-X0)/Hmin
         fracsub = tenth
         if ( fracint<=five ) fracsub = half/fracint
         dlsr%x2 = X1 - fracsub*(X1-X0)
      endif
      Jflag = 1
      X = dlsr%x2
   !  Return to the calling routine to get a value of GX = g(X). -----------
   endif
end subroutine droots
