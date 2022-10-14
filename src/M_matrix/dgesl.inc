!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dgesl(3f) - [M_odepack::matrix] Solve the real system A*X=B or TRANS(A)*X=B
!!               using the factors computed by DGECO or DGEFA.
!!### SYNOPSIS
!!    subroutine dgesl(A,Lda,N,Ipvt,B,Job)
!!    integer,intent(in)          :: Lda
!!    real(kind=dp)               :: A(Lda,*)
!!    integer,intent(in)          :: N
!!    integer,intent(in)          :: Ipvt(*)
!!    real(kind=dp),intent(inout) :: B(*)
!!    integer,intent(in)          :: Job
!!
!!### DESCRIPTION
!!
!!   DGESL solves the double precision system
!!
!!        A * X = B  or  TRANS(A) * X = B
!!
!!   using the factors computed by DGECO or DGEFA.
!!
!!#### Error Condition
!!
!!        A division by zero will occur if the input factor contains a
!!        zero on the diagonal.  Technically this indicates singularity
!!        but it is often caused by improper arguments or improper
!!        setting of LDA .  It will not occur if the subroutines are
!!        called correctly and if DGECO has set RCOND .GT. 0.0
!!        or DGEFA has set INFO .EQ. 0 .
!!
!!   To compute  INVERSE(A) * C  where  C  is a matrix
!!   with  P  columns
!!
!!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!!           IF (RCOND is too small) GO TO ...
!!           DO J = 1, P
!!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!!           ENDDO
!!
!!### ON ENTRY
!!
!!   A
!!
!!   :   the output from DGECO or DGEFA.
!!
!!   LDA
!!
!!   :   the leading dimension of the array  A .
!!
!!   N
!!
!!   :   the order of the matrix  A .
!!
!!   IPVT
!!
!!   :   the pivot vector from DGECO or DGEFA.
!!
!!   B
!!
!!   :   the right hand side vector.
!!
!!   JOB
!!   :
!!        = 0         to solve  A*X = B ,
!!        = nonzero   to solve  TRANS(A)*X = B  where
!!                    TRANS(A)  is the transpose.
!!
!!### ON RETURN
!!
!!   B
!!
!!   :   the solution vector  X .
!!
!!### REFERENCES
!!   J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!!   Stewart, LINPACK Users' Guide, SIAM, 1979.
!!
!
! ### CATEGORY  D2A1
! ### TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
! ### KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
! ### AUTHOR  Moler, C. B., (U. of New Mexico)
! ### ROUTINES CALLED  DAXPY, DDOT
! ### REVISION HISTORY  (YYMMDD)
!     19780814  DATE WRITTEN
!     19890831  Modified array declarations.  (WRB)
!     19890831  REVISION DATE from Version 3.2
!     19891214  Prologue converted to Version 4.0 format.  (BAB)
!     19900326  Removed duplicate information from DESCRIPTION section.
!               (WRB)
!     19920501  Reformatted the REFERENCES section.  (WRB)
!
subroutine dgesl(A,Lda,N,Ipvt,B,Job)
!
integer,intent(in)          :: Lda
real(kind=dp)               :: A(Lda,*)
integer,intent(in)          :: N
integer,intent(in)          :: Ipvt(*)
real(kind=dp),intent(inout) :: B(*)
integer,intent(in)          :: Job
!
integer :: k , kb , l , nm1
real(kind=dp) :: t
!
nm1 = N - 1
if ( Job==0 ) then
   !
   ! JOB = 0 , SOLVE  A * X = B
   ! FIRST SOLVE  L*Y = B
   !
   if ( nm1>=1 ) then
      do k = 1 , nm1
         l = Ipvt(k)
         t = B(l)
         if ( l/=k ) then
            B(l) = B(k)
            B(k) = t
         endif
         call daxpy(N-k,t,A(k+1,k),1,B(k+1),1)
      enddo
   endif
   !
   ! NOW SOLVE  U*X = Y
   !
   do kb = 1 , N
      k = N + 1 - kb
      B(k) = B(k)/A(k,k)
      t = -B(k)
      call daxpy(k-1,t,A(1,k),1,B(1),1)
   enddo
else
   !
   ! JOB = NONZERO, SOLVE  TRANS(A) * X = B
   ! FIRST SOLVE  TRANS(U)*Y = B
   !
   do k = 1 , N
      t = ddot(k-1,A(1,k),1,B(1),1)
      B(k) = (B(k)-t)/A(k,k)
   enddo
   !
   ! NOW SOLVE TRANS(L)*X = Y
   !
   if ( nm1>=1 ) then
      do kb = 1 , nm1
         k = N - kb
         B(k) = B(k) + ddot(N-k,A(k+1,k),1,B(k+1),1)
         l = Ipvt(k)
         if ( l==k ) cycle
         t = B(l)
         B(l) = B(k)
         B(k) = t
      enddo
   endif
endif
end subroutine dgesl
