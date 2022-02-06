!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!###  subroutine nnsc
!!###  numerical solution of sparse nonsymmetric system of linear
!!     equations given ldu-factorization (compressed pointer storage)
!!
!!
!!       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
!!       output variables.. z
!!
!!       parameters used internally..
!!```text
!! fia   - tmp   - temporary vector which gets result of solving  ly = b.
!!       -           size = n.
!!```
!!
!!### internal variables..
!!    jmin, jmax - indices of the first and last positions in a row of
!!      u or l  to be used.
!!
!-----------------------------------------------------------------------
subroutine nnsc(N,R,C,Il,Jl,Ijl,L,D,Iu,Ju,Iju,U,Z,B,Tmp)
!
!      real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum
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
      Tmp(k) = B(R(k))
   enddo
   !   ******  solve  ly = b  by forward substitution  *********************
   do k = 1 , N
      jmin = Il(k)
      jmax = Il(k+1) - 1
      tmpk = -D(k)*Tmp(k)
      Tmp(k) = -tmpk
      if ( jmin<=jmax ) then
         ml = Ijl(k) - jmin
         do j = jmin , jmax
            Tmp(Jl(ml+j)) = Tmp(Jl(ml+j)) + tmpk*L(j)
         enddo
      endif
   enddo
   !   ******  solve  ux = y  by back substitution  ************************
   k = N
   do i = 1 , N
      sum = -Tmp(k)
      jmin = Iu(k)
      jmax = Iu(k+1) - 1
      if ( jmin<=jmax ) then
         mu = Iju(k) - jmin
         do j = jmin , jmax
            sum = sum + U(j)*Tmp(Ju(mu+j))
         enddo
      endif
      Tmp(k) = -sum
      Z(C(k)) = -sum
      k = k - 1
   enddo
end subroutine nnsc
