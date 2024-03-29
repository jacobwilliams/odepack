!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dscal(3f) - [M_odepack::matrix] Multiply a vector by a constant.
!!
!!### SYNOPSIS
!!    subroutine dscal(N,Da,Dx,Incx)
!!    integer , intent(in) :: N
!!    real(kind=dp) , intent(in) :: Da
!!    real(kind=dp) , intent(inout) , dimension(*) :: Dx
!!    integer , intent(in) :: Incx
!!
!!### DESCRIPTION
!!
!!   Replace double precision DX by double precision DA*DX.
!!   For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
!!   where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!!
!!### OPTIONS
!!
!!   N
!!   :   number of elements in input vector(s)
!!
!!   DA
!!   :   double precision scale factor
!!
!!   DX
!!   :   double precision vector with N elements
!!
!!   INCX
!!   :   storage spacing between elements of DX
!!
!!### RETURNS
!!
!!   DX
!!   :   double precision result (unchanged if N.LE.0)
!!
!!### REFERENCES
!!#### B L A S  Subprogram
!!
!!   C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!!   Krogh, Basic linear algebra subprograms for Fortran
!!   usage, Algorithm No. 539, Transactions on Mathematical
!!   Software 5, 3 (September 1979), pp. 308-323.
!!
! ### CATEGORY  D1A6
! ### TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
! ### KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
! ### AUTHOR  Lawson, C. L., (JPL)
!             Hanson, R. J., (SNLA)
!             Kincaid, D. R., (U. of Texas)
!             Krogh, F. T., (JPL)
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791001  DATE WRITTEN
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19900821  Modified to correct problem with a negative increment.
!               (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
!
subroutine dscal(N,Da,Dx,Incx)
!
integer,intent(in)          :: N
real(kind=dp),intent(in)    :: Da
real(kind=dp),intent(inout) :: Dx(*)
integer,intent(in)          :: Incx
!
integer :: i , ix , m , mp1
!
if ( N<=0 ) return
if ( Incx==1 ) then
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
   m = mod(N,5)
   if ( m/=0 ) then
      do i = 1 , m
         Dx(i) = Da*Dx(i)
      enddo
      if ( N<5 ) return
   endif
   mp1 = m + 1
   do i = mp1 , N , 5
      Dx(i) = Da*Dx(i)
      Dx(i+1) = Da*Dx(i+1)
      Dx(i+2) = Da*Dx(i+2)
      Dx(i+3) = Da*Dx(i+3)
      Dx(i+4) = Da*Dx(i+4)
   enddo
else
!
!     Code for increment not equal to 1.
!
   ix = 1
   if ( Incx<0 ) ix = (-N+1)*Incx + 1
   do i = 1 , N
      Dx(ix) = Da*Dx(ix)
      ix = ix + Incx
   enddo
   return
endif
end subroutine dscal
