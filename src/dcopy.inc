!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dcopy(3f) - [M_odepack::matrix] copy a vector
!!
!!### SYNOPSIS
!!        subroutine dcopy(N,Dx,Incx,Dy,Incy)
!!
!!        integer,intent(in)        :: N
!!        real(kind=dp),intent(in)  :: Dx(*)
!!        integer,intent(in)        :: Incx
!!        real(kind=dp),intent(out) :: Dy(*)
!!        integer,intent(in)        :: Incy
!!
!!### DESCRIPTION
!!
!!   Copy double precision DX to double precision DY.
!!   For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
!!   where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!!   defined in a similar way using INCY.
!!
!!### INPUT OPTIONS
!!
!!   N
!!
!!   :   number of elements in input vector(s)
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
!!
!!   DY
!!
!!   :   copy of vector DX (unchanged if N .LE. 0)
!!
!!### REFERENCES
!!#### B L A S  Subprogram
!!    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!!    Krogh, Basic linear algebra subprograms for Fortran
!!    usage, Algorithm No. 539, Transactions on Mathematical
!!    Software 5, 3 (September 1979), pp. 308-323.
!!
! ### BEGIN PROLOGUE  DCOPY
! ### PURPOSE  Copy a vector.
! ### CATEGORY  D1A5
! ### TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
! ### KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
! ### AUTHOR  Lawson, C. L., (JPL)
!            Hanson, R. J., (SNLA)
!            Kincaid, D. R., (U. of Texas)
!            Krogh, F. T., (JPL)
! ### DESCRIPTION
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791001  DATE WRITTEN
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
! ### END PROLOGUE  DCOPY
subroutine dcopy(N,Dx,Incx,Dy,Incy)
!
integer , intent(in) :: N
real(kind=dp) , intent(in) :: Dx(*)
integer , intent(in) :: Incx
real(kind=dp) , intent(out) :: Dy(*)
integer , intent(in) :: Incy
!
integer :: i , ix , iy , m , mp1 , ns
   !
   if ( N<=0 ) return
   if ( Incx==Incy ) then
      if ( Incx<1 ) then
      elseif ( Incx==1 ) then
   !
   !      Code for both increments equal to 1.
   !
   !      Clean-up loop so remaining vector length is a multiple of 7.
   !
         m = mod(N,7)
         if ( m/=0 ) then
            do i = 1 , m
               Dy(i) = Dx(i)
            enddo
            if ( N<7 ) return
         endif
         mp1 = m + 1
         do i = mp1 , N , 7
            Dy(i) = Dx(i)
            Dy(i+1) = Dx(i+1)
            Dy(i+2) = Dx(i+2)
            Dy(i+3) = Dx(i+3)
            Dy(i+4) = Dx(i+4)
            Dy(i+5) = Dx(i+5)
            Dy(i+6) = Dx(i+6)
         enddo
         return
      else
         !
         !      Code for equal, positive, non-unit increments.
         !
         ns = N*Incx
         do i = 1 , ns , Incx
            Dy(i) = Dx(i)
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
      Dy(iy) = Dx(ix)
      ix = ix + Incx
      iy = iy + Incy
   enddo

end subroutine dcopy
