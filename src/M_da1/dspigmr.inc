!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine solves the linear system A * x = b using a scaled
!! preconditioned version of the Generalized Minimal Residual method.
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
!!
!! B is also used as work space when computing
!!
!! the final approximation.
!! (B is the same as V(*,MAXL+1) in the call to DSPIGMR.)
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
!! MAXL 
!!
!! : the maximum allowable order of the matrix HES.
!!
!! MAXLP1
!!
!! : MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
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
!! : real work array used by routine DATV and PSOL.
!!
!! DL 
!!
!! : real work array used for calculation of the residual
!! norm RHO when the method is incomplete (KMP .lt. MAXL).
!! Not needed or referenced in complete case (KMP = MAXL).
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
!! LGMR 
!!
!! : the number of iterations performed and
!! the current order of the upper Hessenberg
!! matrix HES.
!!
!! NPSL 
!!
!! : the number of calls to PSOL.
!!
!! V    
!!
!! : the N by (LGMR+1) array containing the LGMR
!! orthogonal vectors V(*,1) to V(*,LGMR).
!!
!! HES  
!!
!! : the upper triangular factor of the QR decomposition
!! of the (LGMR+1) by lgmr upper Hessenberg matrix whose
!! entries are the scaled inner-products of A*V(*,i)
!! and V(*,k).
!!
!! Q    
!!
!! : real array of length 2*MAXL containing the components
!! of the Givens rotations used in the QR decomposition
!! of HES.  It is loaded in DHEQR and used in DHELS.
!!
!! IFLAG 
!!
!! : integer error flag:
!!
!!                0 means convergence in LGMR iterations, LGMR .le. MAXL.
!!                1 means the convergence test did not pass in MAXL
!!                  iterations, but the residual norm is .lt. 1,
!!                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
!!                2 means the convergence test did not pass in MAXL
!!                  iterations, residual .gt. 1, and X is undefined.
!!                3 means there was a recoverable error in PSOL
!!                  caused by the preconditioner being out of date.
!!               -1 means there was a nonrecoverable error in PSOL.
!!
!-----------------------------------------------------------------------
subroutine dspigmr(Neq,Tn,Y,Savf,B,Wght,N,Maxl,Maxlp1,Kmp,Delta,Hl0,Jpre,Mnewt,f,psol,Npsl,X,V,Hes,Q,Lgmr,Wp,Iwp,Wk,Dl,Iflag)
!
integer , dimension(*) :: Neq
real(kind=dp) :: Tn
real(kind=dp) , dimension(*) :: Y
real(kind=dp) , dimension(*) :: Savf
real(kind=dp) , intent(inout) , dimension(*) :: B
real(kind=dp) , dimension(*) :: Wght
integer :: N
integer , intent(in) :: Maxl
integer :: Maxlp1
integer :: Kmp
real(kind=dp) , intent(inout) :: Delta
real(kind=dp) :: Hl0
integer :: Jpre
integer , intent(in) :: Mnewt
external :: f
external psol
!
real(kind=dp) :: bnrm , bnrm0 , c , dlnrm , prod , rho , s , snormw , tem
real(kind=dp) , intent(inout) , dimension(*) :: Dl , Q , Wk , Wp , X
real(kind=dp) , dimension(Maxlp1,*) :: Hes
integer :: i , i2 , ier , info , ip1 , j , k , ll , llp1
integer,intent(out)  :: Iflag , Lgmr , Npsl
integer :: Iwp(*)
real(kind=dp), intent(inout), dimension(N,*) :: V
!
Iflag = 0
Lgmr = 0
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
if ( bnrm0>Delta ) then
!  Apply inverse of left preconditioner to vector b. --------------------
   ier = 0
   if ( Jpre/=0 .and. Jpre/=2 ) then
      call psol(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,B,1,ier)
      Npsl = 1
      if ( ier/=0 ) then
         !-----------------------------------------------------------------------
         !  This block handles error returns forced by routine PSOL.
         !-----------------------------------------------------------------------
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
      do i = 1 , Maxlp1
         Hes(i,j) = 0.0D0
      enddo
   enddo
!-----------------------------------------------------------------------
!  Main loop to compute the vectors V(*,2) to V(*,MAXL).
!  The running product PROD is needed for the convergence test.
!-----------------------------------------------------------------------
   prod = 1.0D0
   do ll = 1 , Maxl
      Lgmr = ll
!-----------------------------------------------------------------------
!  Call routine DATV to compute VNEW = Abar*v(ll), where Abar is
!  the matrix A with scaling and inverse preconditioner factors applied.
!  Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
!  Call routine DHEQR to update the factors of HES.
!-----------------------------------------------------------------------
      call datv(Neq,Y,Savf,V(1,ll),Wght,X,f,psol,V(1,ll+1),Wk,Wp,Iwp,Hl0,Jpre,ier,Npsl)
      if ( ier/=0 ) then
            !-----------------------------------------------------------------------
            !  This block handles error returns forced by routine PSOL.
            !-----------------------------------------------------------------------
         if ( ier<0 ) Iflag = -1
         if ( ier>0 ) Iflag = 3
         return
      endif
      call dorthog(V(1,ll+1),V,Hes,N,ll,Maxlp1,Kmp,snormw)
      Hes(ll+1,ll) = snormw
      call dheqr(Hes,Maxlp1,ll,Q,info,ll)
      if ( info==ll )then
         Iflag = 2
         return
      endif
      !-----------------------------------------------------------------------
      !  Update RHO, the estimate of the norm of the residual b - A*xl.
      !  If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
      !  necessarily orthogonal for LL .gt. KMP.  The vector DL must then
      !  be computed, and its norm used in the calculation of RHO.
      !-----------------------------------------------------------------------
      prod = prod*Q(2*ll)
      rho = abs(prod*bnrm)
      if ( (ll>Kmp) .and. (Kmp<Maxl) ) then
         if ( ll==Kmp+1 ) then
            !X!call dcopy(N,V(1,1),1,Dl,1)
	    Dl(1:N)=V(1:N,1)
            do i = 1 , Kmp
               ip1 = i + 1
               i2 = i*2
               s = Q(i2)
               c = Q(i2-1)
               do k = 1 , N
                  Dl(k) = s*Dl(k) + c*V(k,ip1)
               enddo
            enddo
         endif
         s = Q(2*ll)
         c = Q(2*ll-1)/snormw
         llp1 = ll + 1
         do k = 1 , N
            Dl(k) = s*Dl(k) + c*V(k,llp1)
         enddo
         dlnrm = dnrm2(N,Dl,1)
         rho = rho*dlnrm
      endif
      !-----------------------------------------------------------------------
      !  Test for convergence.  If passed, compute approximation xl.
      !  if failed and LL .lt. MAXL, then continue iterating.
      !-----------------------------------------------------------------------
      if ( rho<=Delta ) then
         call approximate()
         return
      endif
      if ( ll==Maxl ) exit
      !-----------------------------------------------------------------------
      !  Rescale so that the norm of V(1,LL+1) is one.
      !-----------------------------------------------------------------------
      tem = 1.0D0/snormw
      call dscal(N,tem,V(1,ll+1),1)
   enddo
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
   X(1:N) = B(1:N)
   return
endif

Iflag = 2

contains

subroutine approximate
!-----------------------------------------------------------------------
!  Compute the approximation xl to the solution.
!  Since the vector X was used as work space, and the initial guess
!  of the Newton correction is zero, X must be reset to zero.
!-----------------------------------------------------------------------
   integer :: ll
   ll = Lgmr
   llp1 = ll + 1
   do k = 1 , llp1
      B(k) = 0.0D0
   enddo
   B(1) = bnrm
   call dhels(Hes,Maxlp1,ll,Q,B)
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
   if ( ier/=0 ) then
            !-----------------------------------------------------------------------
            !  This block handles error returns forced by routine PSOL.
            !-----------------------------------------------------------------------
      if ( ier<0 ) Iflag = -1
      if ( ier>0 ) Iflag = 3
      return
   endif
end subroutine approximate

end subroutine dspigmr
