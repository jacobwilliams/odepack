!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine computes the solution to the system A*x = b using a
!! preconditioned version of the Conjugate Gradient algorithm.
!! It is assumed here that the matrix A and the preconditioner
!! matrix M are symmetric positive definite or nearly so.
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
!! R
!!
!! : the right hand side of the system A*x = b.
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
!! of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
!!
!! MAXL
!!
!! : the maximum allowable number of iterates.
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
!! : real work array used by routine DATP.
!!
!! WP
!!
!! : real work array used by preconditioner PSOL.
!!
!! IWP
!!
!! : integer work array used by preconditioner PSOL.
!!
!!#### On return
!!
!! X
!!
!! : the final computed approximation to the solution
!! of the system A*x = b.
!!
!! LPCG
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
!!                0 means convergence in LPCG iterations, LPCG .le. MAXL.
!!                1 means the convergence test did not pass in MAXL
!!                  iterations, but the residual norm is .lt. 1,
!!                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
!!                2 means the convergence test did not pass in MAXL
!!                  iterations, residual .gt. 1, and X is undefined.
!!                3 means there was a recoverable error in PSOL
!!                  caused by the preconditioner being out of date.
!!                4 means there was a zero denominator in the algorithm.
!!                  The system matrix or preconditioner matrix is not
!!                  sufficiently close to being symmetric pos. definite.
!!               -1 means there was a nonrecoverable error in PSOL.
!!
!-----------------------------------------------------------------------
subroutine dpcg(Neq,Tn,Y,Savf,R,Wght,N,Maxl,Delta,Hl0,Jpre,Mnewt,f,psol,Npsl,X,P,W,Z,Lpcg,Wp,Iwp,Wk,Iflag)
!
integer                      :: Neq(*)
real(kind=dp)                :: Tn
real(kind=dp)                :: Y(*)
real(kind=dp)                :: Savf(*)
real(kind=dp)                :: R(*)
real(kind=dp)                :: Wght(*)
integer                      :: N
integer,intent(in)           :: Maxl
real(kind=dp),intent(in)     :: Delta
real(kind=dp)                :: Hl0
integer,intent(in)           :: Jpre
integer,intent(in)           :: Mnewt
external                     :: f
external                     :: psol
integer,intent(out)          :: Npsl
real(kind=dp),intent(inout)  :: X(*)
real(kind=dp),intent(inout)  :: P(*)
real(kind=dp),intent(inout)  :: W(*)
real(kind=dp),intent(inout)  :: Z(*)
integer,intent(out)          :: Lpcg
real(kind=dp),intent(inout)  :: Wp(*)
integer                      :: Iwp(*)
real(kind=dp),intent(inout)  :: Wk(*)
integer,intent(out)          :: Iflag
!
real(kind=dp) :: alpha , beta , bnrm , ptw , rnrm , ztr , ztr0
integer :: i , ier
   !
   Iflag = 0
   Npsl = 0
   Lpcg = 0
   do i = 1 , N
      X(i) = 0.0D0
   enddo
   bnrm = dvnorm(N,R,Wght)
   !  Test for immediate return with X = 0 or X = b. -----------------------
   if ( bnrm>Delta ) then
   !
      ztr = 0.0D0
   else
      if ( Mnewt>0 ) return
      !X!call dcopy(N,R,1,X,1)
      X(1:N) = R(1:N)!X!
      return
   endif
   !  Loop point for PCG iterations. ---------------------------------------
   INFINITE: do
      Lpcg = Lpcg + 1
      !X!call dcopy(N,R,1,Z,1)
      Z(1:N) = R(1:N)!X!
      ier = 0
      if ( Jpre/=0 ) then
         call psol(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,Z,3,ier)
         Npsl = Npsl + 1
         if ( ier/=0 ) then
      !-----------------------------------------------------------------------
      !  This block handles error returns from PSOL.
      !-----------------------------------------------------------------------
            if ( ier<0 ) Iflag = -1
            if ( ier>0 ) Iflag = 3
            return
         endif
      endif
      ztr0 = ztr
      ztr = ddot(N,Z,1,R,1)
      if ( Lpcg==1 ) then
         !X!call dcopy(N,Z,1,P,1)
         P(1:N) = Z(1:N)!X!
      elseif ( ztr0==0.0D0 ) then
      !-----------------------------------------------------------------------
      !  This block handles division by zero errors.
      !-----------------------------------------------------------------------
         Iflag = 4
         return
      else
         beta = ztr/ztr0
         do i = 1 , N
            P(i) = Z(i) + beta*P(i)
         enddo
      endif
      !-----------------------------------------------------------------------
      !   Call DATP to compute A*p and return the answer in W.
      !-----------------------------------------------------------------------
      call datp(Neq,Y,Savf,P,Wght,Hl0,Wk,f,W)
      !
      ptw = ddot(N,P,1,W,1)
      if ( ptw==0.0D0 ) then
         Iflag = 4
         exit INFINITE
      else
         alpha = ztr/ptw
         call daxpy(N,alpha,P,1,X,1)
         alpha = -alpha
         call daxpy(N,alpha,W,1,R,1)
         rnrm = dvnorm(N,R,Wght)
         if ( rnrm<=Delta ) return
         if ( Lpcg<Maxl ) cycle INFINITE
         Iflag = 2
         if ( rnrm<=1.0D0 ) Iflag = 1
         if ( rnrm<=bnrm .and. Mnewt==0 ) Iflag = 1
         exit INFINITE
      endif
   enddo INFINITE
end subroutine dpcg
