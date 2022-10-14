!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!! idamax(3f) - [M_odepack::matrix] Find the smallest index of that
!! component of a vector having the maximum magnitude.
!!
!!### SYNOPSIS
!!    function idamax(n,dx,incx)
!!    integer                  :: idamax
!!    integer,intent(in)       :: n
!!    real(kind=dp),intent(in) :: dx(*)
!!    integer , intent(in)     :: incx
!!
!!### DESCRIPTION
!!
!!  Find smallest index of maximum magnitude of double precision DX.
!!  IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!!  where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!!
!!### OPTIONS
!!
!!   N
!!   :   number of elements in input vector(s)
!!
!!   DX
!!   :   double precision vector with N elements
!!
!!   INCX
!!   :   storage spacing between elements of DX
!!
!!### RETURNS
!!
!!   IDAMAX
!!   :   smallest index (zero if N .LE. 0)
!!
!!### REFERENCES
!!#### B L A S  Subprogram
!!
!!   C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!!   Krogh, Basic linear algebra subprograms for Fortran
!!   usage, Algorithm No. 539, Transactions on Mathematical
!!   Software 5, 3 (September 1979), pp. 308-323.
!!
! ### CATEGORY  D1A2
! ### TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
! ### KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
! ### AUTHOR  Lawson, C. L., (JPL)
!            Hanson, R. J., (SNLA)
!            Kincaid, D. R., (U. of Texas)
!            Krogh, F. T., (JPL)
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791001  DATE WRITTEN
!     19890531  Changed all specific intrinsics to generic.  (WRB)
!     19890531  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19900821  Modified to correct problem with a negative increment.
!               (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
!
function idamax(N,Dx,Incx)
!
integer                  :: idamax
integer,intent(in)       :: N
real(kind=dp),intent(in) :: Dx(*)
integer,intent(in)       :: Incx
!
real(kind=dp)            :: dmax , xmag
integer                  :: i , ix
!
   idamax = 0
   if (n .le. 0) return
   idamax = 1
   if (n .eq. 1) return

   if (incx .ne. 1) then
      !
      !     Code for increments not equal to 1.
      !
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      dmax = abs(dx(ix))
      ix = ix + incx
      do i = 2,n
        xmag = abs(dx(ix))
        if (xmag .gt. dmax) then
          idamax = i
          dmax = xmag
        endif
        ix = ix + incx
      enddo
   else
      !
      !     Code for increments equal to 1.
      !
      dmax = abs(dx(1))
      do i = 2,n
        xmag = abs(dx(i))
        if (xmag .gt. dmax) then
          idamax = i
          dmax = xmag
        endif
      enddo
   endif

end function idamax
