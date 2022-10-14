!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine solves the linear system A * x = b using only a call
!! to the user-supplied routine PSOL (no Krylov iteration).
!! If the norm of the right-hand side vector b is smaller than DELTA,
!! the vector X returned is X = b (if MNEWT = 0) or X = 0 otherwise.
!! PSOL is called with an LR argument of 0.
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
!!
!! WGHT
!!
!! : the vector of length N containing the nonzero
!! elements of the diagonal scaling matrix.
!!
!! N
!!
!! : the order of the matrix A, and the lengths
!! of the vectors WGHT, B and X.
!!
!! DELTA
!!
!! : tolerance on residuals b - A*x in weighted RMS-norm.
!!
!! HL0
!!
!! : current value of (step size h) * (coefficient l0).
!!
!! MNEWT
!!
!! : Newton iteration counter (.ge. 0).
!!
!! WK
!!
!! : real work array used by PSOL.
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
!! NPSL
!!
!! : the number of calls to PSOL.
!!
!! IFLAG
!!
!! : integer error flag:
!!
!!                0 means no trouble occurred.
!!                3 means there was a recoverable error in PSOL
!!                  caused by the preconditioner being out of date.
!!               -1 means there was a nonrecoverable error in PSOL.
!!
!-----------------------------------------------------------------------
subroutine dusol(Neq,Tn,Y,Savf,B,Wght,N,Delta,Hl0,Mnewt,psol,Npsl,X,Wp,Iwp,Wk,Iflag)
!
integer , dimension(*) :: Neq
real(kind=dp) :: Tn
real(kind=dp) , dimension(*) :: Y
real(kind=dp) , dimension(*) :: Savf
real(kind=dp) , dimension(*) :: B
real(kind=dp) , dimension(*) :: Wght
integer :: N
real(kind=dp) , intent(in) :: Delta
real(kind=dp) :: Hl0
integer , intent(in) :: Mnewt
external :: psol
!
real(kind=dp) :: bnrm
integer :: i , ier
integer , intent(out) :: Iflag , Npsl
integer , dimension(*) :: Iwp
real(kind=dp) , dimension(*) :: Wk , Wp , X
   !
   Iflag = 0
   Npsl = 0
   !-----------------------------------------------------------------------
   !  Test for an immediate return with X = 0 or X = b.
   !-----------------------------------------------------------------------
   bnrm = dvnorm(N,B,Wght)
   if ( bnrm>Delta ) then
   !  Make call to PSOL and copy result from B to X. -----------------------
      ier = 0
      call psol(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,B,0,ier)
      Npsl = 1
      if ( ier/=0 ) then
   !-----------------------------------------------------------------------
   !  This block handles error returns forced by routine PSOL.
   !-----------------------------------------------------------------------
         if ( ier<0 ) Iflag = -1
         if ( ier>0 ) Iflag = 3
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
   endif
   !X!call dcopy(N,B,1,X,1)
   X(1:N) = B(1:N)!X!
end subroutine dusol
