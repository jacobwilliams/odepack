!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!###  subroutine nntc
!!### DESCRIPTION
!!   numeric solution of the transpose of a sparse nonsymmetric system
!!   of linear equations given lu-factorization (compressed pointer
!!   storage)
!!
!!
!!       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
!!       output variables.. z
!!
!!       parameters used internally..
!!```text
!! fia   - tmp   - temporary vector which gets result of solving ut y = b
!!       -           size = n.
!!```
!!
!!  internal variables..
!!    jmin, jmax - indices of the first and last positions in a row of
!!      u or l  to be used.
!!
!-----------------------------------------------------------------------
subroutine nntc(N,R,C,Il,Jl,Ijl,L,D,Iu,Ju,Iju,U,Z,B,Tmp)
!
!      real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
integer , intent(in) :: N
integer , intent(in) , dimension(*) :: R
integer , intent(in) , dimension(*) :: C
integer , intent(in) , dimension(*) :: Il
integer , intent(in) , dimension(*) :: Jl
integer , intent(in) , dimension(*) :: Ijl
real(kind=dp) , intent(in) , dimension(*) :: L
real(kind=dp) , intent(in) , dimension(*) :: D
integer , intent(in) , dimension(*) :: Iu
integer , intent(in) , dimension(*) :: Ju
integer , intent(in) , dimension(*) :: Iju
real(kind=dp) , intent(in) , dimension(*) :: U
real(kind=dp) , intent(out) , dimension(*) :: Z
real(kind=dp) , intent(in) , dimension(*) :: B
real(kind=dp) , intent(inout) , dimension(*) :: Tmp
!
integer :: i , j , jmax , jmin , k , ml , mu
real(kind=dp) :: sum , tmpk
   !
   !   ******  set tmp to reordered b  *************************************
   do k = 1 , N
      Tmp(k) = B(C(k))
   enddo
   !   ******  solve  ut y = b  by forward substitution  *******************
   do k = 1 , N
      jmin = Iu(k)
      jmax = Iu(k+1) - 1
      tmpk = -Tmp(k)
      if ( jmin<=jmax ) then
         mu = Iju(k) - jmin
         do j = jmin , jmax
            Tmp(Ju(mu+j)) = Tmp(Ju(mu+j)) + tmpk*U(j)
         enddo
      endif
   enddo
   !   ******  solve  lt x = y  by back substitution  **********************
   k = N
   do i = 1 , N
      sum = -Tmp(k)
      jmin = Il(k)
      jmax = Il(k+1) - 1
      if ( jmin<=jmax ) then
         ml = Ijl(k) - jmin
         do j = jmin , jmax
            sum = sum + L(j)*Tmp(Jl(ml+j))
         enddo
      endif
      Tmp(k) = -sum*D(k)
      Z(R(k)) = Tmp(k)
      k = k - 1
   enddo
end subroutine nntc
