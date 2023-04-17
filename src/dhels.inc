!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This is part of the LINPACK routine DGESL with changes
!! due to the fact that A is an upper Hessenberg matrix.
!!
!! DHELS solves the least squares problem
!!
!!           min (b-A*x, b-A*x)
!!
!!     using the factors computed by DHEQR.
!!
!!### On entry
!!
!!        A       DOUBLE PRECISION(LDA, N)
!!                the output from DHEQR which contains the upper
!!                triangular factor R in the QR decomposition of A.
!!
!!        LDA     INTEGER
!!                the leading dimension of the array  A .
!!
!!        N       INTEGER
!!                A is originally an (N+1) by N matrix.
!!
!!        Q       DOUBLE PRECISION(2*N)
!!                The coefficients of the N givens rotations
!!                used in the QR factorization of A.
!!
!!        B       DOUBLE PRECISION(N+1)
!!                the right hand side vector.
!!
!!### On return
!!
!!        B       the solution vector  x .
!!
!-----------------------------------------------------------------------
!    Modification of LINPACK, by Peter Brown, LLNL.
!    Written 1/13/86.  This version dated 6/20/01.
!
!    BLAS called: DAXPY
!-----------------------------------------------------------------------

subroutine dhels(A,Lda,N,Q,B)
!
integer,intent(in)          :: Lda
real(kind=dp)               :: A(Lda,*)
integer,intent(in)          :: N
real(kind=dp),intent(in)    :: Q(*)
real(kind=dp),intent(inout) :: B(*)
!
real(kind=dp) :: c , s , t , t1 , t2
integer :: iq , k , kb , kp1
!
!         Minimize (b-A*x, b-A*x)
!         First form Q*b.
!
   do k = 1 , N
      kp1 = k + 1
      iq = 2*(k-1) + 1
      c = Q(iq)
      s = Q(iq+1)
      t1 = B(k)
      t2 = B(kp1)
      B(k) = c*t1 - s*t2
      B(kp1) = s*t1 + c*t2
   enddo
   !
   !         Now solve  R*x = Q*b.
   !
   do kb = 1 , N
      k = N + 1 - kb
      B(k) = B(k)/A(k,k)
      t = -B(k)
      call daxpy(k-1,t,A(1,k),1,B(1),1)
   enddo
end subroutine dhels
