!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dnrm2(3f) - [M_odepack::matrix] Compute the Euclidean length (L2 norm) of a vector.
!!
!!### DESCRIPTION
!!
!!   Euclidean norm of the N-vector stored in DX with storage
!!   increment INCX.
!!
!!       If N .LE. 0, return with result = 0.
!!       If N .GE. 1, then INCX must be .GE. 1
!!
!!   Four phase method using two built-in constants that are
!!   hopefully applicable to all machines.
!!
!!         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
!!         CUTHI = minimum of  SQRT(V)      over all known machines.
!!   where
!!
!!         EPS = smallest no. such that EPS + 1. .GT. 1.
!!         U   = smallest positive no.   (underflow limit)
!!         V   = largest  no.            (overflow  limit)
!!
!!   Brief outline of algorithm.
!!
!!       Phase 1 scans zero components.
!!       move to phase 2 when a component is nonzero and .LE. CUTLO
!!       move to phase 3 when a component is .GT. CUTLO
!!       move to phase 4 when a component is .GE. CUTHI/M
!!       where M = N for X() real and M = 2*N for complex.
!!
!!### INPUT
!!
!!  N
!!  :   number of elements in input vector(s)
!!
!!  DX
!!  :   double precision vector with N elements
!!
!!  INCX
!!  :   storage spacing between elements of DX
!!
!!### OUTPUT
!!
!!  DNRM2
!!  :   double precision result (zero if N .LE. 0)
!!
!!### Values for CUTLO and CUTHI.
!!   From the environmental parameters listed in the IMSL converter
!!   document the limiting values are as follows:
!!
!!          CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
!!                        Univac and DEC at 2**(-103)
!!                        Thus CUTLO = 2**(-51) = 4.44089E-16
!!          CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
!!                        Thus CUTHI = 2**(63.5) = 1.30438E19
!!          CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
!!                        Thus CUTLO = 2**(-33.5) = 8.23181D-11
!!          CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
!!          DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
!!          DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
!!
!!### REFERENCES 
!!#### B L A S  Subprogram
!!   C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!!   Krogh, Basic linear algebra subprograms for Fortran
!!   usage, Algorithm No. 539, Transactions on Mathematical
!!   Software 5, 3 (September 1979), pp. 308-323.
!!
! ### CATEGORY  D1A3B
! ### TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
! ### KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
!               LINEAR ALGEBRA, UNITARY, VECTOR
! ### AUTHOR  Lawson, C. L., (JPL)
!             Hanson, R. J., (SNLA)
!             Kincaid, D. R., (U. of Texas)
!             Krogh, F. T., (JPL)
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791001  DATE WRITTEN
!     19890531  Changed all specific intrinsics to generic.  (WRB)
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
! 
pure function dnrm2 (n, dx, incx)

real(kind=dp)            :: dnrm2
real(kind=dp),intent(in) :: dx(*)
integer,intent(in)       :: n
integer,intent(in)       :: incx

real(kind=dp),parameter  :: zero=0.0d0, one=1.0d0
real(kind=dp),parameter  :: cutlo=8.232d-11, cuthi=1.304d19

integer                  :: next
real(kind=dp)            :: hitest, sum, xmax
integer                  :: i, j, nn

      if (n .le. 0)then
         dnrm2  = zero
         return
      endif

      assign 30 to next
      sum = zero
      nn = n * incx
      ! BEGIN MAIN LOOP
      i = 1
   20 continue
      goto next,(30, 50, 70, 110)
   30 continue
      if (abs(dx(i)) .gt. cutlo) goto 85
      assign 50 to next
      xmax = zero
   50 continue
      ! PHASE 1.  SUM IS ZERO
      if (dx(i) .eq. zero) goto 200
      if (abs(dx(i)) .gt. cutlo) goto 85
      ! PREPARE FOR PHASE 2.
      assign 70 to next
      goto 105
  105 continue
      xmax = abs(dx(i))
      sum = sum + (dx(i)/xmax)**2
      goto 200
   70 continue
      ! PHASE 2.  SUM IS SMALL.  SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
      if (abs(dx(i)) .gt. cutlo) then
         ! PREPARE FOR PHASE 3.
         sum = (sum * xmax) * xmax
         goto 85
      endif
  110 continue
      ! COMMON CODE FOR PHASES 2 AND 4.  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
      if (abs(dx(i)) .le. xmax) then
         sum = sum + (dx(i)/xmax)**2
         goto 200
      endif
      sum = one + sum * (xmax / dx(i))**2
      xmax = abs(dx(i))
      goto 200
   85 continue
      ! FOR REAL OR D.P. SET HITEST = CUTHI/N
      ! FOR COMPLEX      SET HITEST = CUTHI/(2*N)
      hitest = cuthi / n
      ! PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
      do j = i,nn,incx
         if (abs(dx(j)) .ge. hitest)then
            ! PREPARE FOR PHASE 4.
            i = j
            assign 110 to next
            sum = (sum / dx(i)) / dx(i)
            goto 105
         endif
         sum = sum + dx(j)**2
      enddo
      dnrm2 = sqrt(sum)
      return

  200 continue
      i = i + incx
      if (i .le. nn) goto 20
      ! END OF MAIN LOOP.  COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
      dnrm2 = xmax * sqrt(sum)
      return
end function dnrm2
