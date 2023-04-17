!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!  dgbsl(3f) - [M_odepack::Matrix] Solve the real band system A*X=B or
!!              TRANS(A)*X=B using the factors computed by DGBCO(3f) or DGBFA().
!!### SYNOPSIS
!!         subroutine dgbsl(Abd,Lda,N,Ml,Mu,Ipvt,B,Job)
!!
!!         integer,intent(in)          :: Lda
!!         real(kind=dp)               :: Abd(Lda,*)
!!         integer,intent(in)          :: N
!!         integer,intent(in)          :: Ml
!!         integer,intent(in)          :: Mu
!!         integer,intent(in)          :: Ipvt(*)
!!         real(kind=dp),intent(inout) :: B(*)
!!         integer,intent(in)          :: Job
!!
!!### DESCRIPTION
!!
!!     DGBSL solves the double precision band system
!!
!!          A * X = B  or  TRANS(A) * X = B
!!
!!     using the factors computed by DGBCO or DGBFA.
!!
!!     To compute  INVERSE(A) * C  where  C  is a matrix
!!     with  P  columns
!!
!!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!!           IF (RCOND is too small) GO TO ...
!!           DO J = 1, P
!!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!!           ENDDO
!!
!!### Error Condition
!!
!!     A division by zero will occur if the input factor contains a
!!     zero on the diagonal.  Technically this indicates singularity
!!     but it is often caused by improper arguments or improper
!!     setting of LDA .  It will not occur if the subroutines are
!!     called correctly and if DGBCO has set RCOND .GT. 0.0
!!     or DGBFA has set INFO .EQ. 0 .
!!
!!### ON ENTRY
!!
!!   ABD
!!
!!   :   the output from DGBCO(3f) or DGBFA(3f).
!!
!!   LDA
!!
!!   :   the leading dimension of the array ABD .
!!
!!   N
!!
!!   :   the order of the original matrix.
!!
!!   ML
!!
!!   :   number of diagonals below the main diagonal.
!!
!!   MU
!!
!!   :   number of diagonals above the main diagonal.
!!
!!   IPVT
!!
!!   :   the pivot vector from DGBCO(3f) or DGBFA(3f).
!!
!!   B
!!
!!   :   the right hand side vector.
!!
!!   JOB
!!   :
!!        = 0         to solve  A*X = B ,
!!        = nonzero   to solve  TRANS(A)*X = B , where
!!                            TRANS(A)  is the transpose.
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
! ### PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
!             the factors computed by DGBCO or DGBFA.
! ### CATEGORY  D2A2
! ### TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
! ### KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
! ### AUTHOR  Moler, C. B., (U. of New Mexico)
! ### ROUTINES CALLED  DAXPY, DDOT
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
subroutine dgbsl(Abd,Lda,N,Ml,Mu,Ipvt,B,Job)
!
integer,intent(in)          :: Lda
real(kind=dp)               :: Abd(Lda,*)
integer,intent(in)          :: N
integer,intent(in)          :: Ml
integer,intent(in)          :: Mu
integer,intent(in)          :: Ipvt(*)
real(kind=dp),intent(inout) :: B(*)
integer,intent(in)          :: Job
!
integer :: k,kb,l,la,lb,lm,m,nm1
real(kind=dp) :: t
!
m = Mu + Ml + 1
nm1 = N - 1
if ( Job/=0 ) then
!
!         JOB = NONZERO, SOLVE  TRANS(A) * X = B
!         FIRST SOLVE  TRANS(U)*Y = B
!
   do k = 1 , N
      lm = min(k,m) - 1
      la = m - lm
      lb = k - lm
      t = ddot(lm,Abd(la,k),1,B(lb),1)
      B(k) = (B(k)-t)/Abd(m,k)
   enddo
!
!         NOW SOLVE TRANS(L)*X = Y
!
   if ( Ml/=0 ) then
      if ( nm1>=1 ) then
         do kb = 1 , nm1
            k = N - kb
            lm = min(Ml,N-k)
            B(k) = B(k) + ddot(lm,Abd(m+1,k),1,B(k+1),1)
            l = Ipvt(k)
            if ( l==k ) cycle
            t = B(l)
            B(l) = B(k)
            B(k) = t
         enddo
      endif
   endif
else
!
!         JOB = 0 , SOLVE  A * X = B
!         FIRST SOLVE L*Y = B
!
   if ( Ml/=0 ) then
      if ( nm1>=1 ) then
         do k = 1 , nm1
            lm = min(Ml,N-k)
            l = Ipvt(k)
            t = B(l)
            if ( l/=k ) then
               B(l) = B(k)
               B(k) = t
            endif
            call daxpy(lm,t,Abd(m+1,k),1,B(k+1),1)
         enddo
      endif
   endif
!
!         NOW SOLVE  U*X = Y
!
   do kb = 1 , N
      k = N + 1 - kb
      B(k) = B(k)/Abd(m,k)
      lm = min(k,m) - 1
      la = m - lm
      lb = k - lm
      t = -B(k)
      call daxpy(lm,t,Abd(la,k),1,B(lb),1)
   enddo
endif
end subroutine dgbsl
