!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dgbfa(3f) - [M_odepack::matrix] Factor a band matrix using Gaussian
!!   elimination.
!!
!!    subroutine dgbfa(Abd,Lda,N,Ml,Mu,Ipvt,Info)
!!    real(kind=dp),intent(inout) :: Abd(Lda,*)
!!    integer,intent(in)    :: Lda
!!    integer,intent(in)    :: N
!!    integer,intent(in)    :: Ml
!!    integer,intent(in)    :: Mu
!!    integer,intent(inout) :: Ipvt(*)
!!    integer,intent(out)   :: Info
!!
!!### DESCRIPTION
!!
!!   DGBFA factors a double precision band matrix by elimination.
!!
!!   DGBFA is usually called by DGBCO, but it can be called
!!   directly with a saving in time if  RCOND  is not needed.
!!
!!### On Entry
!!
!!   ABD
!!
!!   :   contains the matrix in band storage.  The columns
!!   of the matrix are stored in the columns of  ABD  and
!!   the diagonals of the matrix are stored in rows
!!   ML+1 through 2*ML+MU+1 of  ABD .
!!   See the comments below for details.
!!
!!   LDA
!!
!!   :   the leading dimension of the array  ABD .
!!   LDA must be .GE. 2*ML + MU + 1 .
!!
!!   N
!!
!!   :   the order of the original matrix.
!!
!!   ML
!!
!!   :   number of diagonals below the main diagonal.
!!   0 .LE. ML .LT.  N .
!!
!!   MU
!!
!!   :   number of diagonals above the main diagonal.
!!   0 .LE. MU .LT.  N .
!!   More efficient if  ML .LE. MU .
!!
!!### On Return
!!
!!   ABD
!!
!!   :   an upper triangular matrix in band storage and
!!   the multipliers which were used to obtain it.
!!   The factorization can be written  A = L*U  where
!!   L  is a product of permutation and unit lower
!!   triangular matrices and  U  is upper triangular.
!!
!!   IPVT
!!
!!   :   an integer vector of pivot indices.
!!
!!   INFO
!!
!!        = 0  normal value.
!!        = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!!             condition for this subroutine, but it does
!!             indicate that DGBSL will divide by zero if
!!             called.  Use  RCOND  in DGBCO for a reliable
!!             indication of singularity.
!!
!!#### Band Storage
!!
!!   If  A  is a band matrix, the following program segment
!!   will set up the input.
!!
!!       ML = (band width below the diagonal)
!!       MU = (band width above the diagonal)
!!       M = ML + MU + 1
!!       DO J = 1, N
!!            I1 = MAX(1, J-MU)
!!            I2 = MIN(N, J+ML)
!!            DO I = I1, I2
!!                K = I - J + M
!!                ABD(K,J) = A(I,J)
!!            ENDDO
!!       ENDDO
!!
!!   This uses rows ML+1 through 2*ML+MU+1 of ABD.
!!   In addition, the first ML rows in ABD are used for
!!   elements generated during the triangularization.
!!   The total number of rows needed in ABD is 2*ML+MU+1.
!!   The ML+MU by ML+MU upper left triangle and the
!!   ML by ML lower right triangle are not referenced.
!!
!!### REFERENCES
!!   J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!!   Stewart, LINPACK Users' Guide, SIAM, 1979.
!!
! ### CATEGORY  D2A2
! ### TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
! ### KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
! ### AUTHOR  Moler, C. B., (U. of New Mexico)
! ### ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
! ### REVISION HISTORY  (YYMMDD)
!     19780814  DATE WRITTEN
!     19890531  Changed all specific intrinsics to generic.  (WRB)
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19900326  Removed duplicate information from DESCRIPTION section.
!               (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
!
subroutine dgbfa(Abd,Lda,N,Ml,Mu,Ipvt,Info)
!
integer , intent(in) :: Lda
real(kind=dp) , intent(inout) :: Abd(Lda,*)
integer , intent(in) :: N
integer , intent(in) :: Ml
integer , intent(in) :: Mu
integer , intent(inout) :: Ipvt(*)
integer , intent(out) :: Info
!
integer :: i , i0 , j , j0 , j1 , ju , jz , k , kp1 , l , lm , m , mm , nm1
real(kind=dp) :: t
!
m = Ml + Mu + 1
Info = 0
!
!      ZERO INITIAL FILL-IN COLUMNS
!
j0 = Mu + 2
j1 = min(N,m) - 1
if ( j1>=j0 ) then
   do jz = j0 , j1
      i0 = m + 1 - jz
      do i = i0 , Ml
         Abd(i,jz) = 0.0D0
      enddo
   enddo
endif
jz = j1
ju = 0
!
!      GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
nm1 = N - 1
if ( nm1>=1 ) then
   do k = 1 , nm1
      kp1 = k + 1
!
!         ZERO NEXT FILL-IN COLUMN
!
      jz = jz + 1
      if ( jz<=N ) then
         if ( Ml>=1 ) then
            do i = 1 , Ml
               Abd(i,jz) = 0.0D0
            enddo
         endif
      endif
!
!         FIND L = PIVOT INDEX
!
      lm = min(Ml,N-k)
      l = idamax(lm+1,Abd(m,k),1) + m - 1
      Ipvt(k) = l + k - m
!
!         ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
      if ( Abd(l,k)==0.0D0 ) then
         Info = k
      else
!
!            INTERCHANGE IF NECESSARY
!
         if ( l/=m ) then
            t = Abd(l,k)
            Abd(l,k) = Abd(m,k)
            Abd(m,k) = t
         endif
!
!            COMPUTE MULTIPLIERS
!
         t = -1.0D0/Abd(m,k)
         call dscal(lm,t,Abd(m+1,k),1)
!
!            ROW ELIMINATION WITH COLUMN INDEXING
!
         ju = min(max(ju,Mu+Ipvt(k)),N)
         mm = m
         if ( ju>=kp1 ) then
            do j = kp1 , ju
               l = l - 1
               mm = mm - 1
               t = Abd(l,j)
               if ( l/=mm ) then
                  Abd(l,j) = Abd(mm,j)
                  Abd(mm,j) = t
               endif
               call daxpy(lm,t,Abd(m+1,k),1,Abd(mm+1,j),1)
            enddo
         endif
      endif
   enddo
endif
Ipvt(N) = N
if ( Abd(m,N)==0.0D0 ) Info = N
end subroutine dgbfa
