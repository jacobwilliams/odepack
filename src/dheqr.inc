!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!   This routine performs a QR decomposition of an upper
!!   Hessenberg matrix A.  There are two options available:
!!
!!          (1)  performing a fresh decomposition
!!          (2)  updating the QR factors by adding a row and a
!!               column to the matrix A.
!!
!!   DHEQR decomposes an upper Hessenberg matrix by using Givens
!!   rotations.
!!
!!### On entry
!!
!!        A       DOUBLE PRECISION(LDA, N)
!!                the matrix to be decomposed.
!!
!!        LDA     INTEGER
!!                the leading dimension of the array  A .
!!
!!        N       INTEGER
!!                A is an (N+1) by N Hessenberg matrix.
!!
!!        IJOB    INTEGER
!!                = 1     means that a fresh decomposition of the
!!                        matrix A is desired.
!!                .ge. 2  means that the current decomposition of A
!!                        will be updated by the addition of a row
!!                        and a column.
!!### On return
!!
!!        A       the upper triangular matrix R.
!!                The factorization can be written Q*A = R, where
!!                Q is a product of Givens rotations and R is upper
!!                triangular.
!!
!!        Q       DOUBLE PRECISION(2*N)
!!                the factors c and s of each Givens rotation used
!!                in decomposing A.
!!
!!        INFO    INTEGER
!!                = 0  normal value.
!!                = k  if  A(k,k) .eq. 0.0 .  This is not an error
!!                     condition for this subroutine, but it does
!!                     indicate that DHELS will divide by zero
!!                     if called.
!!
!-----------------------------------------------------------------------
!    Modification of LINPACK, by Peter Brown, LLNL.
!    Written 1/13/86.  This version dated 6/20/01.
!-----------------------------------------------------------------------
subroutine dheqr(A,Lda,N,Q,Info,Ijob)
!
integer,intent(in)          :: Lda
real(kind=dp),intent(inout) :: A(Lda,*)
integer,intent(in)          :: N
real(kind=dp),intent(inout) :: Q(*)
integer,intent(out)         :: Info
integer,intent(in)          :: Ijob
!
real(kind=dp) :: c , s , t , t1 , t2
integer :: i , iq , j , k , km1 , kp1 , nm1
!
if ( Ijob>1 ) then
!
!  The old factorization of A will be updated.  A row and a column
!  has been added to the matrix A.
!  N by N-1 is now the old size of the matrix.
!
   nm1 = N - 1
!
!  Multiply the new column by the N previous Givens rotations.
!
   do k = 1 , nm1
      i = 2*(k-1) + 1
      t1 = A(k,N)
      t2 = A(k+1,N)
      c = Q(i)
      s = Q(i+1)
      A(k,N) = c*t1 - s*t2
      A(k+1,N) = s*t1 + c*t2
   enddo
!
!  Complete update of decomposition by forming last Givens rotation,
!  and multiplying it times the column vector (A(N,N), A(N+1,N)).
!
   Info = 0
   t1 = A(N,N)
   t2 = A(N+1,N)
   if ( t2==0.0D0 ) then
      c = 1.0D0
      s = 0.0D0
   elseif ( abs(t2)<abs(t1) ) then
      t = t2/t1
      c = 1.0D0/sqrt(1.0D0+t*t)
      s = -c*t
   else
      t = t1/t2
      s = -1.0D0/sqrt(1.0D0+t*t)
      c = -s*t
   endif
   iq = 2*N - 1
   Q(iq) = c
   Q(iq+1) = s
   A(N,N) = c*t1 - s*t2
   if ( A(N,N)==0.0D0 ) Info = N
else
!
!  A new facorization is desired.
!
!      QR decomposition without pivoting
!
   Info = 0
   do k = 1 , N
      km1 = k - 1
      kp1 = k + 1
!
!            Compute kth column of R.
!            First, multiply the kth column of A by the previous
!            k-1 Givens rotations.
!
      if ( km1>=1 ) then
         do j = 1 , km1
            i = 2*(j-1) + 1
            t1 = A(j,k)
            t2 = A(j+1,k)
            c = Q(i)
            s = Q(i+1)
            A(j,k) = c*t1 - s*t2
            A(j+1,k) = s*t1 + c*t2
         enddo
      endif
!
!            Compute Givens components c and s
!
      iq = 2*km1 + 1
      t1 = A(k,k)
      t2 = A(kp1,k)
      if ( t2==0.0D0 ) then
         c = 1.0D0
         s = 0.0D0
      elseif ( abs(t2)<abs(t1) ) then
         t = t2/t1
         c = 1.0D0/sqrt(1.0D0+t*t)
         s = -c*t
      else
         t = t1/t2
         s = -1.0D0/sqrt(1.0D0+t*t)
         c = -s*t
      endif
      Q(iq) = c
      Q(iq+1) = s
      A(k,k) = c*t1 - s*t2
      if ( A(k,k)==0.0D0 ) Info = k
   enddo
   return
endif
end subroutine dheqr
