!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   ddot(3f) - [M_odepack::matrix] Compute the inner product of two vectors.
!!
!!### DESCRIPTION
!!   Returns the dot product of double precision DX and DY.
!!   DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
!!   where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!!   defined in a similar way using INCY.
!!
!!   Description of Parameters
!!
!!### INPUT
!!   N
!!
!!   :   number of elements in input vector(s)
!!
!!   DX
!!
!!   :   double precision vector with N elements
!!   INCX
!!
!!   :   storage spacing between elements of DX
!!
!!   DY
!!   :   double precision vector with N elements
!!
!!   INCY
!!   :   storage spacing between elements of DY
!!
!!### OUTPUT
!!
!!   DDOT
!!   :   double precision dot product (zero if N .LE. 0)
!!
!!### REFERENCES
!!#### B L A S  Subprogram
!!   C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!!   Krogh, Basic linear algebra subprograms for Fortran
!!   usage, Algorithm No. 539, Transactions on Mathematical
!!   Software 5, 3 (September 1979), pp. 308-323.
!!
! ### BEGIN PROLOGUE  DDOT
! ### PURPOSE  Compute the inner product of two vectors.
! ### CATEGORY  D1A4
! ### TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
! ### KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
! ### AUTHOR  Lawson, C. L., (JPL)
!            Hanson, R. J., (SNLA)
!            Kincaid, D. R., (U. of Texas)
!            Krogh, F. T., (JPL)
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791001  DATE WRITTEN
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
! ### END PROLOGUE  DDOT

pure function ddot(N,Dx,Incx,Dy,Incy)
!
real(kind=dp)            :: ddot
integer,intent(in)       :: N
real(kind=dp),intent(in) :: Dx(*)
integer,intent(in)       :: Incx
real(kind=dp),intent(in) :: Dy(*)
integer, intent(in)      :: Incy
!
integer :: i, ix, iy, m, mp1, ns
!
   ddot = 0.0D0
   if ( N<=0 ) return
   if ( Incx==Incy ) then
      if ( Incx<1 ) then
      elseif ( Incx==1 ) then
         !
         !      Code for both increments equal to 1.
         !
         !      Clean-up loop so remaining vector length is a multiple of 5.
         !
         m = mod(N,5)
         if ( m/=0 ) then
            do i = 1 , m
               ddot = ddot + Dx(i)*Dy(i)
            enddo
            if ( N<5 ) return
         endif
         mp1 = m + 1
         do i = mp1 , N , 5
            ddot = ddot + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)*Dy(i+2) + Dx(i+3)*Dy(i+3) + Dx(i+4)*Dy(i+4)
         enddo
         return
      else
         !
         !      Code for equal, positive, non-unit increments.
         !
         ns = N*Incx
         do i = 1 , ns , Incx
            ddot = ddot + Dx(i)*Dy(i)
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
      ddot = ddot + Dx(ix)*Dy(iy)
      ix = ix + Incx
      iy = iy + Incy
   enddo
end function ddot
