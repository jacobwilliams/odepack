!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!! ddecbt(3f) -[M_odepack] Block-tridiagonal matrix decomposition routine.
!!
!!### DESCRIPTION
!! The input matrix contains three blocks of elements in each block-row,
!! including blocks in the (1,3) and (N,N-2) block positions.
!! DDECBT uses block Gauss elimination and Subroutines DGEFA and DGESL
!! for solution of blocks.  Partial pivoting is done within
!! block-rows only.
!!
!! Note: this version uses LINPACK routines DGEFA/DGESL instead of
!! of dec/sol for solution of blocks, and it uses the BLAS routine DDOT
!! for dot product calculations.
!!
!!### INPUT
!!
!! M
!!
!! : order of each block.
!!
!! N
!!
!! : number of blocks in each direction of the matrix.
!!   N must be 4 or more.  The complete matrix has order M*N.
!!
!! A
!!
!! : M by M by N array containing diagonal blocks.
!!   A(i,j,k) contains the (i,j) element of the k-th block.
!!
!! B
!!
!! : M by M by N array containing the super-diagonal blocks
!!   (in B(*,*,k) for k = 1,...,N-1) and the block in the (N,N-2)
!!   block position (in B(*,*,N)).
!!
!! C
!!
!! : M by M by N array containing the subdiagonal blocks
!!   (in C(*,*,k) for k = 2,3,...,N) and the block in the
!!   (1,3) block position (in C(*,*,1)).
!!
!! IP
!!
!! : integer array of length M*N for working storage.
!!
!!### OUTPUT
!!
!! A,B,C
!!
!! : M by M by N arrays containing the block-LU decomposition
!!   of the input matrix.
!!
!! IP
!!
!! : M by N array of pivot information.  IP(*,k) contains
!!   information for the k-th digonal block.
!!
!! IER
!!
!! : 0  if no trouble occurred, or
!!       = -1 if the input value of M or N was illegal, or
!!       = k  if a singular matrix was found in the k-th diagonal block.
!!
!! Use DSOLBT to solve the associated linear system.
!!
!! External routines required: DGEFA and DGESL (from LINPACK) and
!! DDOT (from the BLAS, or Basic Linear Algebra package).
!-----------------------------------------------------------------------
!  Written by A. C. Hindmarsh.
!  Latest revision:  November 10, 1983 (ACH)
!  Reference:  UCID-30150
!              Solution of Block-Tridiagonal Systems of Linear
!              Algebraic Equations
!              A.C. Hindmarsh
!              February 1977
!-----------------------------------------------------------------------
subroutine ddecbt(M,N,A,B,C,Ip,Ier)
!
integer                     :: M
integer,intent(in)          :: N
real(kind=dp),intent(inout) :: A(M,M,N)
real(kind=dp),intent(inout) :: B(M,M,N)
real(kind=dp),intent(inout) :: C(M,M,N)
integer                     :: Ip(M,N)
integer, intent(inout)      :: Ier
!
real(kind=dp) :: dpp
integer :: i , j , k , km1 , nm1 , nm2
!
   if ( M<1 .or. N<4 ) then
      Ier = -1
      return
   else
      nm1 = N - 1
      nm2 = N - 2
   ! Process the first block-row. -----------------------------------------
      call dgefa(A,M,M,Ip,Ier)
      k = 1
      if ( Ier==0 ) then
         do j = 1 , M
            call dgesl(A,M,M,Ip,B(1,j,1),0)
            call dgesl(A,M,M,Ip,C(1,j,1),0)
         enddo
   ! Adjust B(*,*,2). -----------------------------------------------------
         do j = 1 , M
            do i = 1 , M
               dpp = ddot(M,C(i,1,2),M,C(1,j,1),1)
               B(i,j,2) = B(i,j,2) - dpp
            enddo
         enddo
   ! Main loop.  Process block-rows 2 to N-1. -----------------------------
         do k = 2 , nm1
            km1 = k - 1
            do j = 1 , M
               do i = 1 , M
                  dpp = ddot(M,C(i,1,k),M,B(1,j,km1),1)
                  A(i,j,k) = A(i,j,k) - dpp
               enddo
            enddo
            call dgefa(A(1,1,k),M,M,Ip(1,k),Ier)
            if ( Ier/=0 )then
               Ier = k
            return
         endif
            do j = 1 , M
               call dgesl(A(1,1,k),M,M,Ip(1,k),B(1,j,k),0)
            enddo
         enddo
   ! Process last block-row and return. -----------------------------------
         do j = 1 , M
            do i = 1 , M
               dpp = ddot(M,B(i,1,N),M,B(1,j,nm2),1)
               C(i,j,N) = C(i,j,N) - dpp
            enddo
         enddo
         do j = 1 , M
            do i = 1 , M
               dpp = ddot(M,C(i,1,N),M,B(1,j,nm1),1)
               A(i,j,N) = A(i,j,N) - dpp
            enddo
         enddo
         call dgefa(A(1,1,N),M,M,Ip(1,N),Ier)
         k = N
         if ( Ier==0 ) return
      endif
   endif
   ! Error returns. -------------------------------------------------------
   Ier = k
end subroutine ddecbt
