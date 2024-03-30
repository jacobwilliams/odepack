!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   daxpy(3f) - [M_odepack::matrix] Compute a constant times a vector plus a vector.
!!
!!### SYNOPSIS
!!    subroutine daxpy(N,Da,Dx,Incx,Dy,Incy)
!!    integer,intent(in)          :: N
!!    real(kind=dp),intent(in)    :: Da
!!    real(kind=dp),intent(in)    :: Dx(*)
!!    integer,intent(in)          :: Incx
!!    real(kind=dp),intent(inout) :: Dy(*)
!!    integer,intent(in)          :: Incy
!!
!!### DESCRIPTION
!!
!!   daxpy(3f) computes a constant times a vector plus a vector.
!!   It uses unrolled loops for increments equal to one.
!!
!!   Overwrite double precision DY with double precision DA*DX + DY.
!!   For I = 0 to N-1, replace  DY(LY+I*INCY) with
!!   DA*DX(LX+I*INCX) + DY(LY+I*INCY),
!!
!!   where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!!   defined in a similar way using INCY.
!!
!!### INPUT OPTIONS
!!
!!   N
!!
!!   :   number of elements in input vector(s)
!!
!!   DA
!!
!!   :   double precision scalar multiplier
!!
!!   DX
!!
!!   :   double precision vector with N elements
!!
!!   INCX
!!
!!   :   storage spacing between elements of DX
!!
!!   DY
!!
!!   :   double precision vector with N elements
!!
!!   INCY
!!
!!   :   storage spacing between elements of DY
!!
!!### RETURNS
!!   DY
!!
!!   :   double precision result (unchanged if N .LE. 0)
!!
!!### REFERENCES
!!####  B L A S  Subprogram
!!   C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!!   Krogh, Basic linear algebra subprograms for Fortran
!!   usage, Algorithm No. 539, Transactions on Mathematical
!!   Software 5, 3 (September 1979), pp. 308-323.
!!
! ### BEGIN PROLOGUE  DAXPY
! ### CATEGORY  D1A7
! ### TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
! ### KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
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
!     19920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)

subroutine daxpy(N,Da,Dx,Incx,Dy,Incy)
!
integer , intent(in)          :: N
real(kind=dp) , intent(in)    :: Da
real(kind=dp) , intent(in)    :: Dx(*)
integer , intent(in)          :: Incx
real(kind=dp) , intent(inout) :: Dy(*)
integer , intent(in)          :: Incy
!
integer :: i , ix , iy , m , mp1 , ns
!
if ( N<=0 .or. Da==0.0D0 ) return
if ( Incx==Incy ) then
   if ( Incx<1 ) then
   elseif ( Incx==1 ) then
!
!      Code for both increments equal to 1.
!
!      Clean-up loop so remaining vector length is a multiple of 4.
!
      m = mod(N,4)
      if ( m/=0 ) then
         do i = 1 , m
            Dy(i) = Dy(i) + Da*Dx(i)
         enddo
         if ( N<4 ) return
      endif
      mp1 = m + 1
      do i = mp1 , N , 4
         Dy(i) = Dy(i) + Da*Dx(i)
         Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
         Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
         Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
      enddo
      return
   else
!
!      Code for equal, positive, non-unit increments.
!
      ns = N*Incx
      do i = 1 , ns , Incx
         Dy(i) = Da*Dx(i) + Dy(i)
      enddo
      return
   endif
endif
!
!      Code for unequal or nonpositive increments.
!
ix = 1
iy = 1
if ( Incx<0 ) ix = (-N+1)*Incx + 1
if ( Incy<0 ) iy = (-N+1)*Incy + 1
do i = 1 , N
   Dy(iy) = Dy(iy) + Da*Dx(ix)
   ix = ix + Incx
   iy = iy + Incy
enddo
end subroutine daxpy
