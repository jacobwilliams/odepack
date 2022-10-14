!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine computes the solution to the system A*x = b using a
!! scaled preconditioned version of the Conjugate Gradient algorithm.
!!
!! It is assumed here that the scaled matrix D**-1 * A * D and the
!! scaled preconditioner D**-1 * M * D are close to being
!! symmetric positive definite.
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
!!### On return
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
!!                  the scaled matrix or scaled preconditioner is not
!!                  sufficiently close to being symmetric pos. definite.
!!               -1 means there was a nonrecoverable error in PSOL.
!!
!-----------------------------------------------------------------------
subroutine dpcgs(Neq,Tn,Y,Savf,R,Wght,N,Maxl,Delta,Hl0,Jpre,Mnewt,f,psol,Npsl,X,P,W,Z,Lpcg,Wp,Iwp,Wk,Iflag)
!
integer , dimension(*) :: Neq
real(kind=dp) :: Tn
real(kind=dp) , dimension(*) :: Y
real(kind=dp) , dimension(*) :: Savf
real(kind=dp) , dimension(*) :: R
real(kind=dp) , dimension(*) :: Wght
integer :: N
integer , intent(in) :: Maxl
real(kind=dp) , intent(in) :: Delta
real(kind=dp) :: Hl0
integer , intent(in) :: Jpre
integer , intent(in) :: Mnewt
external :: f
!
real(kind=dp) :: alpha , beta , bnrm , ptw , rnrm , ztr , ztr0
integer :: i , ier
integer , intent(out) :: Iflag , Lpcg , Npsl
integer , dimension(*) :: Iwp
real(kind=dp) , intent(inout) , dimension(*) :: P , W , Wk , Wp , X , Z
!
external psol
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
      X(1:N) = R(1:N) !X!
      X(1:N) = R(1:N)
      return
   endif
   !  Loop point for PCG iterations. ---------------------------------------
   INFINITE: do
      Lpcg = Lpcg + 1
      !X!call dcopy(N,R,1,Z,1)
      Z(1:N) = R(1:N) !X!
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
      ztr = 0.0D0
      do i = 1 , N
         ztr = ztr + Z(i)*R(i)*Wght(i)**2
      enddo
      if ( Lpcg/=1 ) then
         if ( ztr0==0.0D0 ) then
            !-----------------------------------------------------------------------
            !  This block handles division by zero errors.
            !-----------------------------------------------------------------------
            Iflag = 4
            return
         endif
         beta = ztr/ztr0
         do i = 1 , N
            P(i) = Z(i) + beta*P(i)
         enddo
      else
         !X!call dcopy(N,Z,1,P,1)
         P(1:N) = Z(1:N) !X!
      endif
      !-----------------------------------------------------------------------
      !   Call DATP to compute A*p and return the answer in W.
      !-----------------------------------------------------------------------
      call datp(Neq,Y,Savf,P,Wght,Hl0,Wk,f,W)
      !
      ptw = 0.0D0
      do i = 1 , N
         ptw = ptw + P(i)*W(i)*Wght(i)**2
      enddo
      if ( ptw==0.0D0 ) then
         !-----------------------------------------------------------------------
         !  This block handles division by zero errors.
         !-----------------------------------------------------------------------
         Iflag = 4
         return
      endif
      alpha = ztr/ptw
      call daxpy(N,alpha,P,1,X,1)
      alpha = -alpha
      call daxpy(N,alpha,W,1,R,1)
      rnrm = dvnorm(N,R,Wght)
      if ( rnrm<=Delta ) return
      if ( Lpcg>=Maxl ) exit INFINITE
   enddo INFINITE
   Iflag = 2
   if ( rnrm<=1.0D0 ) Iflag = 1
   if ( rnrm<=bnrm .and. Mnewt==0 ) Iflag = 1
end subroutine dpcgs
