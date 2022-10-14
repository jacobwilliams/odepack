!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine solves the linear system A * x = b using a scaled
!! preconditioned version of the Incomplete Orthogonalization Method.
!! An initial guess of x = 0 is assumed.
!!
!!### On entry
!!
!! NEQ
!!
!! : problem size, passed to F and PSOL (NEQ(1) = N).
!!
!! TN
!!
!! : current value of t.
!!
!! Y
!!
!! : array containing current dependent variable vector.
!!
!! SAVF
!!
!! : array containing current value of f(t,y).
!!
!! B
!!
!! : the right hand side of the system A*x = b.
!! B is also used as work space when computing the
!! final approximation.
!! (B is the same as V(*,MAXL+1) in the call to DSPIOM.)
!!
!! WGHT
!!
!! : array of length N containing scale factors.
!! 1/WGHT(i) are the diagonal elements of the diagonal
!! scaling matrix D.
!!
!! N
!!
!! : the order of the matrix A, and the lengths
!! of the vectors Y, SAVF, B, WGHT, and X.
!!
!! MAXL
!!
!! : the maximum allowable order of the matrix HES.
!!
!! KMP
!!
!! : the number of previous vectors the new vector VNEW
!! must be made orthogonal to.  KMP .le. MAXL.
!!
!! DELTA
!!
!! : tolerance on residuals b - A*x in weighted RMS-norm.
!!
!! HL0
!!
!! : current value of (step size h) * (coefficient l0).
!!
!! JPRE
!!
!! : preconditioner type flag.
!!
!! MNEWT
!!
!! : Newton iteration counter (.ge. 0).
!!
!! WK
!!
!! : real work array of length N used by DATV and PSOL.
!!
!! WP
!!
!! : real work array used by preconditioner PSOL.
!!
!! IWP
!!
!! : integer work array used by preconditioner PSOL.
!!
!!### On return
!!
!! X
!!
!! : the final computed approximation to the solution
!! of the system A*x = b.
!!
!! V
!!
!! : the N by (LIOM+1) array containing the LIOM
!! orthogonal vectors V(*,1) to V(*,LIOM).
!!
!! HES
!!
!! : the LU factorization of the LIOM by LIOM upper
!! Hessenberg matrix whose entries are the
!! scaled inner products of A*V(*,k) and V(*,i).
!!
!! IPVT
!!
!! : an integer array containg pivoting information.
!! It is loaded in DHEFA and used in DHESL.
!!
!! LIOM
!!
!! : the number of iterations performed, and current
!! order of the upper Hessenberg matrix HES.
!!
!! NPSL
!!
!! : the number of calls to PSOL.
!!
!! IFLAG
!!
!! : integer error flag:
!!
!!                0 means convergence in LIOM iterations, LIOM.le.MAXL.
!!                1 means the convergence test did not pass in MAXL
!!                  iterations, but the residual norm is .lt. 1,
!!                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
!!                2 means the convergence test did not pass in MAXL
!!                  iterations, residual .gt. 1, and X is undefined.
!!                3 means there was a recoverable error in PSOL
!!                  caused by the preconditioner being out of date.
!!               -1 means there was a nonrecoverable error in PSOL.
!!
!-----------------------------------------------------------------------
subroutine dspiom(Neq,Tn,Y,Savf,B,Wght,N,Maxl,Kmp,Delta,Hl0,Jpre,Mnewt,f,psol,Npsl,X,V,Hes,Ipvt,Liom,Wp,Iwp,Wk,Iflag)
!
integer , dimension(*) :: Neq
real(kind=dp) :: Tn
real(kind=dp) , dimension(*) :: Y
real(kind=dp) , dimension(*) :: Savf
real(kind=dp) , intent(inout) , dimension(*) :: B
real(kind=dp) , dimension(*) :: Wght
integer :: N
integer :: Maxl
integer :: Kmp
real(kind=dp) , intent(inout) :: Delta
real(kind=dp) :: Hl0
integer :: Jpre
integer , intent(in) :: Mnewt
external :: f
!
real(kind=dp) :: bnrm , bnrm0 , prod , rho , snormw , tem
real(kind=dp) , intent(inout) , dimension(Maxl,Maxl) :: Hes
integer :: i , ier , info , j , k , ll , lm1
integer , intent(out) :: Iflag , Liom , Npsl
integer , dimension(*) :: Ipvt , Iwp
real(kind=dp) , dimension(N,*) :: V
real(kind=dp) , dimension(*) :: Wk , Wp , X
external psol
   !
   Iflag = 0
   Liom = 0
   Npsl = 0
   !-----------------------------------------------------------------------
   !  The initial residual is the vector b.  Apply scaling to b, and test
   !  for an immediate return with X = 0 or X = b.
   !-----------------------------------------------------------------------
   do i = 1 , N
      V(i,1) = B(i)*Wght(i)
   enddo
   bnrm0 = dnrm2(N,V,1)
   bnrm = bnrm0
   IFDELTA: if ( bnrm0>Delta ) then
      !  Apply inverse of left preconditioner to vector b. --------------------
      ier = 0
      if ( Jpre/=0 .and. Jpre/=2 ) then
         call psol(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,B,1,ier)
         Npsl = 1
         if ( ier/=0 )then
            !  This block handles error returns forced by routine PSOL.
            if ( ier<0 ) Iflag = -1
            if ( ier>0 ) Iflag = 3
         return
         endif
         !  Calculate norm of scaled vector V(*,1) and normalize it. -------------
         do i = 1 , N
            V(i,1) = B(i)*Wght(i)
         enddo
         bnrm = dnrm2(N,V,1)
         Delta = Delta*(bnrm/bnrm0)
      endif
      tem = 1.0D0/bnrm
      call dscal(N,tem,V(1,1),1)
      !  Zero out the HES array. ----------------------------------------------
      do j = 1 , Maxl
         do i = 1 , Maxl
            Hes(i,j) = 0.0D0
         enddo
      enddo
      !-----------------------------------------------------------------------
      !  Main loop on LL = l to compute the vectors V(*,2) to V(*,MAXL).
      !  The running product PROD is needed for the convergence test.
      !-----------------------------------------------------------------------
      prod = 1.0D0
      do ll = 1 , Maxl
         Liom = ll
         !-----------------------------------------------------------------------
         !  Call routine DATV to compute VNEW = Abar*v(l), where Abar is
         !  the matrix A with scaling and inverse preconditioner factors applied.
         !  Call routine DORTHOG to orthogonalize the new vector vnew = V(*,l+1).
         !  Call routine DHEFA to update the factors of HES.
         !-----------------------------------------------------------------------
         call datv(Neq,Y,Savf,V(1,ll),Wght,X,f,psol,V(1,ll+1),Wk,Wp,Iwp,Hl0,Jpre,ier,Npsl)
         if ( ier/=0 ) then
            !  This block handles error returns forced by routine PSOL.
            if ( ier<0 ) Iflag = -1
            if ( ier>0 ) Iflag = 3
         return
         endif
         call dorthog(V(1,ll+1),V,Hes,N,ll,Maxl,Kmp,snormw)
         call dhefa(Hes,Maxl,ll,Ipvt,info,ll)
         lm1 = ll - 1
         if ( ll>1 .and. Ipvt(lm1)==lm1 ) prod = prod*Hes(ll,lm1)
         if ( info/=ll ) then
            !-----------------------------------------------------------------------
            !  Update RHO, the estimate of the norm of the residual b - A*x(l).
            !  test for convergence.  If passed, compute approximation x(l).
            !  If failed and l .lt. MAXL, then continue iterating.
            !-----------------------------------------------------------------------
            rho = bnrm*snormw*abs(prod/Hes(ll,ll))
            if ( rho<=Delta ) then
               call approximate()
               return
            endif
            if ( ll==Maxl ) exit
         else
            !-----------------------------------------------------------------------
            !  The last pivot in HES was found to be zero.
            !  If vnew = 0 or l = MAXL, take an error return with IFLAG = 2.
            !  otherwise, continue the iteration without a convergence test.
            !-----------------------------------------------------------------------
            if ( snormw==0.0D0 ) exit IFDELTA
            if ( ll==Maxl ) exit IFDELTA
         endif
         !  If l .lt. MAXL, store HES(l+1,l) and normalize the vector v(*,l+1).
         Hes(ll+1,ll) = snormw
         tem = 1.0D0/snormw
         call dscal(N,tem,V(1,ll+1),1)
      enddo
      !-----------------------------------------------------------------------
      !  l has reached MAXL without passing the convergence test:
      !  If RHO is not too large, compute a solution anyway and return with
      !  IFLAG = 1.  Otherwise return with IFLAG = 2.
      !-----------------------------------------------------------------------
      if ( rho<=1.0D0 ) then
         Iflag = 1
         call approximate()
         return
      elseif ( rho<=bnrm .and. Mnewt==0 ) then
         Iflag = 1
         call approximate()
         return
      endif
   elseif ( Mnewt>0 ) then
      do i = 1 , N
         X(i) = 0.0D0
      enddo
      return
   else
      !X!call dcopy(N,B,1,X,1)
      X(1:N) = B(1:N)!X!
      return
   endif IFDELTA

   Iflag = 2

contains

subroutine approximate()
!-----------------------------------------------------------------------
!  Compute the approximation x(l) to the solution.
!  Since the vector X was used as work space, and the initial guess
!  of the Newton correction is zero, X must be reset to zero.
!-----------------------------------------------------------------------
   integer :: ll
   ll = Liom
   do k = 1 , ll
      B(k) = 0.0D0
   enddo
   B(1) = bnrm
   call dhesl(Hes,Maxl,ll,Ipvt,B)
   do k = 1 , N
      X(k) = 0.0D0
   enddo
   do i = 1 , ll
      call daxpy(N,B(i),V(1,i),1,X,1)
   enddo
   do i = 1 , N
      X(i) = X(i)/Wght(i)
   enddo
   if ( Jpre<=1 ) return
   call psol(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,X,2,ier)
   Npsl = Npsl + 1
   if ( ier==0 ) return
   !-----------------------------------------------------------------------
   !  This block handles error returns forced by routine PSOL.
   !-----------------------------------------------------------------------
   if ( ier<0 ) Iflag = -1
   if ( ier>0 ) Iflag = 3

end subroutine approximate

end subroutine dspiom
