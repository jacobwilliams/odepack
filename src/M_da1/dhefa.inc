!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### DESCRIPTION
!!   This routine is a modification of the LINPACK routine DGEFA and
!!   performs an LU decomposition of an upper Hessenberg matrix A.
!!   There are two options available:
!!
!!          (1)  performing a fresh factorization
!!          (2)  updating the LU factors by adding a row and a
!!               column to the matrix A.
!!
!!   DHEFA factors an upper Hessenberg matrix by elimination.
!!
!!### ON ENTRY
!!
!!   A
!!
!!   : the matrix to be factored.
!!
!!   LDA
!!
!!   : the leading dimension of the array  A .
!!
!!   N
!!
!!   : the order of the matrix  A .
!!
!!   JOB
!!
!!       JOB = 1    means that a fresh factorization of the
!!                  matrix A is desired.
!!       JOB .ge. 2 means that the current factorization of A
!!                  will be updated by the addition of a row
!!                  and a column.
!!
!!### ON RETURN
!!
!!   A
!!
!!   : an upper triangular matrix and the multipliers
!!     which were used to obtain it.
!!
!!     The factorization can be written  A = L*U  where
!!     L is a product of permutation and unit lower
!!     triangular matrices and U is upper triangular.
!!
!!   IPVT
!!
!!   : an integer vector of pivot indices.
!!
!!   INFO
!!
!!            = 0  normal value.
!!            = k  if  U(k,k) .eq. 0.0 .  This is not an error
!!                 condition for this subroutine, but it does
!!                 indicate that DHESL will divide by zero if called.
!!
!-----------------------------------------------------------------------
!      Modification of LINPACK, by Peter Brown, LLNL.
!      Written 7/20/83.  This version dated 6/20/2001.
!
!      BLAS called: DAXPY, IDAMAX
!-----------------------------------------------------------------------
subroutine dhefa(A,Lda,N,Ipvt,Info,Job)
!
integer , intent(in)          :: Lda
real(kind=dp) , intent(inout) :: A(Lda,*)
integer , intent(in)          :: N
integer , intent(inout)       :: Ipvt(*)
integer , intent(out)         :: Info
integer , intent(in)          :: Job
!
integer :: j , k , km1 , kp1 , l , nm1
real(kind=dp) :: t
!
   if ( Job>1 ) then
   !
   !  The old factorization of A will be updated.  A row and a column
   !  has been added to the matrix A.
   !  N-1 is now the old order of the matrix.
   !
      nm1 = N - 1
   !
   !  Perform row interchanges on the elements of the new column, and
   !  perform elimination operations on the elements using the multipliers.
   !
      if ( nm1>1 ) then
         do k = 2 , nm1
            km1 = k - 1
            l = Ipvt(km1)
            t = A(l,N)
            if ( l/=km1 ) then
               A(l,N) = A(km1,N)
               A(km1,N) = t
            endif
            A(k,N) = A(k,N) + A(k,km1)*t
         enddo
      endif
   !
   !  Complete update of factorization by decomposing last 2x2 block.
   !
      Info = 0
   !
   !         Find L = pivot index
   !
      l = idamax(2,A(nm1,nm1),1) + nm1 - 1
      Ipvt(nm1) = l
   !
   !         Zero pivot implies this column already triangularized
   !
      if ( A(l,nm1)==0.0D0 ) then
         Info = nm1
      else
   !
   !            Interchange if necessary
   !
         if ( l/=nm1 ) then
            t = A(l,nm1)
            A(l,nm1) = A(nm1,nm1)
            A(nm1,nm1) = t
         endif
   !
   !            Compute multipliers
   !
         t = -1.0D0/A(nm1,nm1)
         A(N,nm1) = A(N,nm1)*t
   !
   !            Row elimination with column indexing
   !
         t = A(l,N)
         if ( l/=nm1 ) then
            A(l,N) = A(nm1,N)
            A(nm1,N) = t
         endif
         A(N,N) = A(N,N) + t*A(N,nm1)
      endif
      Ipvt(N) = N
      if ( A(N,N)==0.0D0 ) Info = N
   else
   !
   !  A new facorization is desired.  This is essentially the LINPACK
   !  code with the exception that we know there is only one nonzero
   !  element below the main diagonal.
   !
   !      Gaussian elimination with partial pivoting
   !
      Info = 0
      nm1 = N - 1
      if ( nm1>=1 ) then
         do k = 1 , nm1
            kp1 = k + 1
   !
   !         Find L = pivot index
   !
            l = idamax(2,A(k,k),1) + k - 1
            Ipvt(k) = l
   !
   !         Zero pivot implies this column already triangularized
   !
            if ( A(l,k)==0.0D0 ) then
               Info = k
            else
   !
   !            Interchange if necessary
   !
               if ( l/=k ) then
                  t = A(l,k)
                  A(l,k) = A(k,k)
                  A(k,k) = t
               endif
   !
   !            Compute multipliers
   !
               t = -1.0D0/A(k,k)
               A(k+1,k) = A(k+1,k)*t
   !
   !            Row elimination with column indexing
   !
               do j = kp1 , N
                  t = A(l,j)
                  if ( l/=k ) then
                     A(l,j) = A(k,j)
                     A(k,j) = t
                  endif
                  call daxpy(N-k,t,A(k+1,k),1,A(k+1,j),1)
               enddo
            endif
         enddo
      endif
      Ipvt(N) = N
      if ( A(N,N)==0.0D0 ) Info = N
      return
   endif
end subroutine dhefa
