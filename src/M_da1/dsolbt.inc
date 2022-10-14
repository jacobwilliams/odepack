!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! Solution of block-tridiagonal linear system.
!! Coefficient matrix must have been previously processed by DDECBT.
!! M, N, A,B,C, and IP  must not have been changed since call to DDECBT.
!! Written by A. C. Hindmarsh.
!!### Input
!!
!! M
!!
!! : order of each block.
!!
!! N
!!
!! : number of blocks in each direction of matrix.
!!
!! A,B,C
!!
!! : M by M by N arrays containing block LU decomposition
!! of coefficient matrix from DDECBT.
!!
!! IP
!!
!! : M by N integer array of pivot information from DDECBT.
!!
!! Y
!!
!! : array of length M*N containg the right-hand side vector
!! (treated as an M by N array here).
!!
!!### Output
!!
!! Y
!!
!! : solution vector, of length M*N.
!!
!! External routines required: DGESL (LINPACK) and DDOT (BLAS).
!-----------------------------------------------------------------------
subroutine dsolbt(M,N,A,B,C,Y,Ip)
!
integer                     :: M
integer,intent(in)          :: N
real(kind=dp)               :: A(M,M,N)
real(kind=dp)               :: B(M,M,N)
real(kind=dp)               :: C(M,M,N)
real(kind=dp),intent(inout) :: Y(M,N)
integer                     :: Ip(M,N)
!
real(kind=dp) :: dpr
integer :: i , k , kb , km1 , kp1 , nm1 , nm2
!
   nm1 = N - 1
   nm2 = N - 2
   ! Forward solution sweep. ----------------------------------------------
   call dgesl(A,M,M,Ip,Y,0)

   do k = 2 , nm1
      km1 = k - 1
      do i = 1 , M
         dpr = ddot(M,C(i,1,k),M,Y(1,km1),1)
         Y(i,k) = Y(i,k) - dpr
      enddo
      call dgesl(A(1,1,k),M,M,Ip(1,k),Y(1,k),0)
   enddo

   do i = 1 , M
      dpr = ddot(M,C(i,1,N),M,Y(1,nm1),1) + ddot(M,B(i,1,N),M,Y(1,nm2),1)
      Y(i,N) = Y(i,N) - dpr
   enddo

   call dgesl(A(1,1,N),M,M,Ip(1,N),Y(1,N),0)

   ! Backward solution sweep. ---------------------------------------------
   do kb = 1 , nm1
      k = N - kb
      kp1 = k + 1
      do i = 1 , M
         dpr = ddot(M,B(i,1,k),M,Y(1,kp1),1)
         Y(i,k) = Y(i,k) - dpr
      enddo
   enddo

   do i = 1 , M
      dpr = ddot(M,C(i,1,1),M,Y(1,3),1)
      Y(i,1) = Y(i,1) - dpr
   enddo

end subroutine dsolbt
