!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
! ### NAME
!!   dgefa(3f) - [M_odepack::matrix] Factor a matrix using Gaussian elimination.
!!
!!      subroutine dgefa(A,Lda,N,Ipvt,Info)
!!      real(kind=dp),intent(inout) :: A(Lda,*)
!!      integer,intent(in)          :: Lda
!!      integer,intent(in)          :: N
!!      integer,intent(out)         :: Ipvt(*)
!!      integer,intent(out)         :: Info
!!
!!### DESCRIPTION
!!
!!   DGEFA factors a double precision matrix by Gaussian elimination.
!!
!!   DGEFA is usually called by DGECO, but it can be called
!!   directly with a saving in time if  RCOND  is not needed.
!!   (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!!
!!### OPTIONS
!!
!!   A
!!
!!   :   the matrix to be factored.
!!
!!   LDA
!!
!!   :   the leading dimension of the array  A .
!!
!!   N
!!
!!   :   the order of the matrix  A .
!!
!!### RETURNS
!!
!!   A
!!
!!   :   an upper triangular matrix and the multipliers
!!       which were used to obtain it.
!!       The factorization can be written  A = L*U  where
!!
!!   L
!!
!!   :   is a product of permutation and unit lower
!!       triangular matrices and  U  is upper triangular.
!!
!!   IPVT
!!
!!   :   an integer vector of pivot indices.
!!
!!   INFO
!!
!!   :
!!        = 0  normal value.
!!        = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!!             condition for this subroutine, but it does
!!             indicate that DGESL or DGEDI will divide by zero
!!             if called.  Use  RCOND  in DGECO for a reliable
!!             indication of singularity.
!!
!!### REFERENCES
!!   J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!!   Stewart, LINPACK Users' Guide, SIAM, 1979.
!!
! ### CATEGORY  D2A1
! ### TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
! ### KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!              MATRIX FACTORIZATION
! ### AUTHOR  Moler, C. B., (U. of New Mexico)
! ### ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
! ### REVISION HISTORY  (YYMMDD)
!     19780814  DATE WRITTEN
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19900326  Removed duplicate information from DESCRIPTION section.
!             (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
!
subroutine dgefa(A,Lda,N,Ipvt,Info)
!
integer,intent(in)          :: Lda
real(kind=dp),intent(inout) :: A(Lda,*)
integer,intent(in)          :: N
integer,intent(out)         :: Ipvt(*)
integer,intent(out)         :: Info
!
integer :: j , k , kp1 , l , nm1
real(kind=dp) :: t
   !
   !      GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
   !
   Info = 0
   nm1 = N - 1
   if ( nm1>=1 ) then
      do k = 1 , nm1
         kp1 = k + 1
   !
   !         FIND L = PIVOT INDEX
   !
         l = idamax(N-k+1,A(k,k),1) + k - 1
         Ipvt(k) = l
   !
   !         ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
   !
         if ( A(l,k)==0.0D0 ) then
            Info = k
         else
   !
   !            INTERCHANGE IF NECESSARY
   !
            if ( l/=k ) then
               t = A(l,k)
               A(l,k) = A(k,k)
               A(k,k) = t
            endif
   !
   !            COMPUTE MULTIPLIERS
   !
            t = -1.0D0/A(k,k)
            call dscal(N-k,t,A(k+1,k),1)
   !
   !            ROW ELIMINATION WITH COLUMN INDEXING
   !
            do j = kp1 , N
               t = A(l,j)
               if ( l/=k ) then
                  A(l,j) = A(k,j)
                  A(k,j) = t
               endif
               call daxpy(N-k,t,A(k+1,k),1,A(k+1,j),1)
            enddo
         endif
      enddo
   endif
   Ipvt(N) = N
   if ( A(N,N)==0.0D0 ) Info = N
end subroutine dgefa
