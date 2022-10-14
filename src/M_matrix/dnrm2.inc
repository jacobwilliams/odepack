!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dnrm2(3f) - [M_odepack::matrix] Compute the Euclidean length (L2 norm) of a vector.
!!
!!### SYNOPSIS
!!      double precision function dnrm2(n,x,incx)
!!      integer,intent(in) :: incx,n
!!      double precision,intent(in) :: x(*)
!!
!!### DESCRIPTION
!!
!!   Euclidean norm of the N-vector stored in DX with storage
!!   increment INCX, so that
!!
!!        DNRM2 := sqrt( x'*x )
!!
!!### INPUT
!!
!!  N
!!  :   number of elements in input vector(s)
!!
!!  DX
!!  :   double precision vector with N elements
!>  dimensioned to at least ( 1 + ( N - 1 )*abs( INCX ) )
!!
!!  INCX
!!  :   storage spacing between elements of DX
!!
!!      + If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!!      + If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!!      + If INCX = 0, x isn't a vector so there is no need to call
!!        this subroutine.  If you call it anyway, it will count x(1)
!!        in the vector norm N times.
!!
!!### OUTPUT
!!
!!  DNRM2
!!  :   double precision result (zero if N .LE. 0)
!-----------------------------------------------------------------------------------------------------------------------------------
! ### REFERENCES
! #### B L A S  Subprogram
!    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!    Krogh, Basic linear algebra subprograms for Fortran
!    usage, Algorithm No. 539, Transactions on Mathematical
!    Software 5, 3 (September 1979), pp. 308-323.
!
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
!  Authors:
!  ========
!
!   author Edward Anderson, Lockheed Martin
!   date August 2016
!
!   ingroup single_blas_level1
!
!  Contributors:
!  ==================
!
!  Weslley Pereira, University of Colorado Denver, USA
!
!  Further Details:
!  =====================
!
!   Anderson E. (2017)
!   Algorithm 978: Safe Scaling in the Level 1 BLAS
!   ACM Trans Math Softw 44:1--28
!   https://doi.org/10.1145/3061665
!
!   Blue, James L. (1978)
!   A Portable Fortran Program to Find the Euclidean Norm of a Vector
!   ACM Trans Math Softw 4:15--23
!   https://doi.org/10.1145/355769.355771
!
!-----------------------------------------------------------------------------------------------------------------------------------
pure function DNRM2( n, x, incx )
   real(dp) :: DNRM2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: one  = 1.0_dp
   real(dp), parameter :: maxN = huge(0.0_dp)
!  ..
!  .. Blue's ccaling constants ..
   real(dp), parameter :: tsml = real(radix(0._dp), dp)**ceiling( &
       (minexponent(0._dp) - 1) * 0.5_dp)
   real(dp), parameter :: tbig = real(radix(0._dp), dp)**floor( &
       (maxexponent(0._dp) - digits(0._dp) + 1) * 0.5_dp)
   real(dp), parameter :: ssml = real(radix(0._dp), dp)**( - floor( &
       (minexponent(0._dp) - 1) * 0.5_dp))
   real(dp), parameter :: sbig = real(radix(0._dp), dp)**( - ceiling( &
       (maxexponent(0._dp) - digits(0._dp) + 1) * 0.5_dp))
!  ..
!  .. Scalar Arguments ..
   integer,intent(in) :: incx, n
!  ..
!  .. Array Arguments ..
   real(dp),intent(in) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(dp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
   DNRM2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   DNRM2 = scl*sqrt( sumsq )
   return
end function
!-----------------------------------------------------------------------------------------------------------------------------------
