!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!   This is essentially the LINPACK routine DGESL except for changes
!!   due to the fact that A is an upper Hessenberg matrix.
!!
!!   DHESL solves the real system A * x = b
!!   using the factors computed by DHEFA.
!!
!!### On entry
!!
!!        A       DOUBLE PRECISION(LDA, N)
!!                the output from DHEFA.
!!
!!        LDA     INTEGER
!!                the leading dimension of the array  A .
!!
!!        N       INTEGER
!!                the order of the matrix  A .
!!
!!        IPVT    INTEGER(N)
!!                the pivot vector from DHEFA.
!!
!!        B       DOUBLE PRECISION(N)
!!                the right hand side vector.
!!
!!### On return
!!
!!        B       the solution vector  x .
!!
!-----------------------------------------------------------------------
!    Modification of LINPACK, by Peter Brown, LLNL.
!    Written 7/20/83.  This version dated 6/20/01.
!
!    BLAS called: DAXPY
!-----------------------------------------------------------------------
subroutine dhesl(A,Lda,N,Ipvt,B)
!
integer,intent(in)           :: Lda
real(kind=dp)                :: A(Lda,*)
integer,intent(in)           :: N
integer,intent(in)           :: Ipvt(*)
real(kind=dp),intent(inout)  :: B(*)
!
integer :: k , kb , l , nm1
real(kind=dp) :: t
!
!
nm1 = N - 1
!
!         Solve  A * x = b
!         First solve  L*y = b
!
if ( nm1>=1 ) then
   do k = 1 , nm1
      l = Ipvt(k)
      t = B(l)
      if ( l/=k ) then
         B(l) = B(k)
         B(k) = t
      endif
      B(k+1) = B(k+1) + t*A(k+1,k)
   enddo
endif
!
!         Now solve  U*x = y
!
do kb = 1 , N
   k = N + 1 - kb
   B(k) = B(k)/A(k,k)
   t = -B(k)
   call daxpy(k-1,t,A(1,k),1,B(1),1)
enddo
end subroutine dhesl
