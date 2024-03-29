!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!! nnfc(3f) - [M_odepack] numerical LDU-factorization of sparse nonsymmetric matrix
!!
!!### DESCRIPTION
!!   numerical ldu-factorization of sparse nonsymmetric matrix and
!!   solution of system of linear equations (compressed pointer
!!   storage)
!!```text
!!       input variables..  n, r, c, ic, ia, ja, a, b,
!!                          il, jl, ijl, lmax, iu, ju, iju, umax
!!       output variables.. z, l, d, u, flag
!!```
!!       parameters used internally..
!!```text
!! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
!! nia   - jrl       of the factorization,  jrl(k)  points to the head
!!       -           of a linked list in  jrl  of column indices j
!!       -           such j .lt. k and  l(k,j)  is nonzero.  zero
!!       -           indicates the end of the list.  irl(j)  (j.lt.k)
!!       -           points to the smallest i such that i .ge. k and
!!       -           l(i,j)  is nonzero.
!!       -           size of each = n.
!! fia   - row   - holds intermediate values in calculation of  u and l.
!!       -           size = n.
!! fia   - tmp   - holds new right-hand side  b*  for solution of the
!!       -           equation ux = b*.
!!       -           size = n.
!!```
!!
!!### internal variables..
!!
!!    jmin, jmax - indices of the first and last positions in a row to
!!      be examined.
!!
!!    sum - used in calculating  tmp.
!!
!-----------------------------------------------------------------------
subroutine nnfc(N,R,C,Ic,Ia,Ja,A,Z,B,Lmax,Il,Jl,Ijl,L,D,Umax,Iu,Ju,Iju,U,Row,Tmp,Irl,Jrl,Flag)
!
!      real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
!      real tmp(*), lki, sum, dk
integer , intent(in) :: N
integer , intent(in) , dimension(*) :: R
integer , intent(in) , dimension(*) :: C
integer , intent(in) , dimension(*) :: Ic
integer , intent(in) , dimension(*) :: Ia
integer , intent(in) , dimension(*) :: Ja
real(kind=dp) , intent(in) , dimension(*) :: A
real(kind=dp) , intent(out) , dimension(*) :: Z
real(kind=dp) , intent(in) , dimension(*) :: B
integer , intent(in) :: Lmax
integer , intent(in) , dimension(*) :: Il
integer , intent(in) , dimension(*) :: Jl
integer , intent(in) , dimension(*) :: Ijl
real(kind=dp) , intent(out) , dimension(*) :: L
real(kind=dp) , intent(out) , dimension(*) :: D
integer , intent(in) :: Umax
integer , intent(in) , dimension(*) :: Iu
integer , intent(in) , dimension(*) :: Ju
integer , intent(in) , dimension(*) :: Iju
real(kind=dp) , intent(inout) , dimension(*) :: U
real(kind=dp) , intent(inout) , dimension(*) :: Row
real(kind=dp) , intent(inout) , dimension(*) :: Tmp
integer , intent(inout) , dimension(*) :: Irl
integer , intent(inout) , dimension(*) :: Jrl
integer , intent(out) :: Flag
!
real(kind=dp) :: dk , lki , sum
integer :: i , i1 , i2 , ijlb , j , jmax , jmin , k , mu , rk
   !
   !   ******  initialize pointers and test storage  ***********************
   if ( Il(N+1)-1>Lmax ) then
   !
   !  ** error.. insufficient storage for l
      Flag = 4*N + 1
      return
   elseif ( Iu(N+1)-1>Umax ) then
   !  ** error.. insufficient storage for u
      Flag = 7*N + 1
      return
   else
      do k = 1 , N
         Irl(k) = Il(k)
         Jrl(k) = 0
      enddo
   !
   !   ******  for each row  ***********************************************
      do k = 1 , N
   !   ******  reverse jrl and zero row where kth row of l will fill in  ***
         Row(k) = 0
         i1 = 0
         if ( Jrl(k)/=0 ) then
            i = Jrl(k)
            do
               i2 = Jrl(i)
               Jrl(i) = i1
               i1 = i
               Row(i) = 0
               i = i2
               if ( i==0 ) exit
            enddo
         endif
   !   ******  set row to zero where u will fill in  ***********************
         jmin = Iju(k)
         jmax = jmin + Iu(k+1) - Iu(k) - 1
         if ( jmin<=jmax ) then
            do j = jmin , jmax
               Row(Ju(j)) = 0
            enddo
         endif
   !   ******  place kth row of a in row  **********************************
         rk = R(k)
         jmin = Ia(rk)
         jmax = Ia(rk+1) - 1
         do j = jmin , jmax
            Row(Ic(Ja(j))) = A(j)
         enddo
   !   ******  initialize sum, and link through jrl  ***********************
         sum = B(rk)
         i = i1
         if ( i/=0 ) then
            !   ******  assign the kth row of l and adjust row, sum  ****************
            do
               lki = -Row(i)
               !   ******  if l is not required, then comment out the following line  **
               L(Irl(i)) = -lki
               sum = sum + lki*Tmp(i)
               jmin = Iu(i)
               jmax = Iu(i+1) - 1
               if ( jmin<=jmax ) then
                  mu = Iju(i) - jmin
                  do j = jmin , jmax
                     Row(Ju(mu+j)) = Row(Ju(mu+j)) + lki*U(j)
                  enddo
               endif
               i = Jrl(i)
               if ( i==0 ) exit
            enddo
         endif
         !   ******  assign kth row of u and diagonal d, set tmp(k)  *************
         if ( Row(k)==0.0D0 ) then
            !  ** error.. zero pivot
            Flag = 8*N + k
         return
         else
            dk = 1.0D0/Row(k)
            D(k) = dk
            Tmp(k) = sum*dk
            if ( k==N ) cycle
            jmin = Iu(k)
            jmax = Iu(k+1) - 1
            if ( jmin<=jmax ) then
               mu = Iju(k) - jmin
               do j = jmin , jmax
                  U(j) = Row(Ju(mu+j))*dk
               enddo
            endif
            !
            !   ******  update irl and jrl, keeping jrl in decreasing order  ********
            i = i1
            if ( i/=0 )then
               do
                  Irl(i) = Irl(i) + 1
                  i1 = Jrl(i)
                  if ( Irl(i)<Il(i+1) ) then
                     ijlb = Irl(i) - Il(i) + Ijl(i)
                     j = Jl(ijlb)
                     do while ( i<=Jrl(j) )
                        j = Jrl(j)
                     enddo
                     Jrl(i) = Jrl(j)
                     Jrl(j) = i
                  endif
                  i = i1
                  if ( i==0 ) exit
               enddo
         endif
         endif
         if ( Irl(k)<Il(k+1) ) then
            j = Jl(Ijl(k))
            Jrl(k) = Jrl(j)
            Jrl(j) = k
         endif
      enddo
   !
   !   ******  solve  ux = tmp  by back substitution  **********************
      k = N
      do i = 1 , N
         sum = Tmp(k)
         jmin = Iu(k)
         jmax = Iu(k+1) - 1
         if ( jmin<=jmax ) then
            mu = Iju(k) - jmin
            do j = jmin , jmax
               sum = sum - U(j)*Tmp(Ju(mu+j))
            enddo
         endif
         Tmp(k) = sum
         Z(C(k)) = sum
         k = k - 1
      enddo
      Flag = 0
   endif
end subroutine nnfc
