!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine orthogonalizes the vector VNEW against the previous
!! KMP vectors in the V array.  It uses a modified Gram-Schmidt
!! orthogonalization procedure with conditional reorthogonalization.
!! This is the version of 28 may 1986.
!!
!!      On entry
!!
!!         VNEW = the vector of length N containing a scaled product
!!                of the Jacobian and the vector V(*,LL).
!!
!!         V    = the N x l array containing the previous LL
!!                orthogonal vectors v(*,1) to v(*,LL).
!!
!!         HES  = an LL x LL upper Hessenberg matrix containing,
!!                in HES(i,k), k.lt.LL, scaled inner products of
!!                A*V(*,k) and V(*,i).
!!
!!        LDHES = the leading dimension of the HES array.
!!
!!         N    = the order of the matrix A, and the length of VNEW.
!!
!!         LL   = the current order of the matrix HES.
!!
!!          KMP = the number of previous vectors the new vector VNEW
!!                must be made orthogonal to (KMP .le. MAXL).
!!
!!
!!      On return
!!
!!         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
!!                where i0 = MAX(1, LL-KMP+1).
!!
!!         HES  = upper Hessenberg matrix with column LL filled in with
!!                scaled inner products of A*V(*,LL) and V(*,i).
!!
!!       SNORMW = L-2 norm of VNEW.
!!
!-----------------------------------------------------------------------
subroutine dorthog(Vnew,V,Hes,N,Ll,Ldhes,Kmp,Snormw)

real(kind=dp)               :: Vnew(*)
integer                     :: N
real(kind=dp)               :: V(N,*)
integer,intent(in)          :: Ldhes
real(kind=dp),intent(inout) :: Hes(Ldhes,*)
integer,intent(in)          :: Ll
integer,intent(in)          :: Kmp
real(kind=dp),intent(inout) :: Snormw

real(kind=dp) :: arg, sumdsq, tem, vnrm
integer       :: i, i0

   !  Get norm of unaltered VNEW for later use. ----------------------------
   vnrm = dnrm2(N,Vnew,1)
   !-----------------------------------------------------------------------
   !  Do modified Gram-Schmidt on VNEW = A*v(LL).
   !  Scaled inner products give new column of HES.
   !  Projections of earlier vectors are subtracted from VNEW.
   !-----------------------------------------------------------------------
   i0 = max(1,Ll-Kmp+1)
   do i = i0 , Ll
      Hes(i,Ll) = ddot(N,V(1,i),1,Vnew,1)
      tem = -Hes(i,Ll)
      call daxpy(N,tem,V(1,i),1,Vnew,1)
   enddo
   !-----------------------------------------------------------------------
   !  Compute SNORMW = norm of VNEW.
   !  If VNEW is small compared to its input value (in norm), then
   !  reorthogonalize VNEW to V(*,1) through V(*,LL).
   !  Correct if relative correction exceeds 1000*(unit roundoff).
   !  finally, correct SNORMW using the dot products involved.
   !-----------------------------------------------------------------------
   Snormw = dnrm2(N,Vnew,1)
   if ( vnrm+0.001D0*Snormw/=vnrm ) return
   sumdsq = 0.0D0
   do i = i0 , Ll
      tem = -ddot(N,V(1,i),1,Vnew,1)
      if ( Hes(i,Ll)+0.001D0*tem/=Hes(i,Ll) ) then
         Hes(i,Ll) = Hes(i,Ll) - tem
         call daxpy(N,tem,V(1,i),1,Vnew,1)
         sumdsq = sumdsq + tem**2
      endif
   enddo
   if ( sumdsq==0.0D0 ) return
   arg = max(0.0D0,Snormw**2-sumdsq)
   Snormw = sqrt(arg)

end subroutine dorthog
