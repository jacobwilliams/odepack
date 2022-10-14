!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This function computes the norm of a banded N by N matrix,
!! stored in the array A, that is consistent with the weighted max-norm
!! on vectors, with weights stored in the array W.
!! ML and MU are the lower and upper half-bandwidths of the matrix.
!! NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
!! In terms of the matrix elements a(i,j), the norm is given by:
!!
!!     DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
!!
!-----------------------------------------------------------------------
pure function dbnorm (n, a, nra, ml, mu, w)
integer,intent(in)       :: n
integer,intent(in)       :: nra
real(kind=dp),intent(in) :: a(nra,n)
integer,intent(in)       :: ml
integer,intent(in)       :: mu
real(kind=dp),intent(in) :: w(n)
real(kind=dp)            :: dbnorm

integer :: i, i1, jlo, jhi, j
real(kind=dp) :: an, sum

   an = 0.0d0

   do i = 1,n
      sum = 0.0d0
      i1 = i + mu + 1
      jlo = max(i-ml,1)
      jhi = min(i+mu,n)

      do j = jlo,jhi
         sum = sum + abs(a(i1-j,j))/w(j)
      enddo

      an = max(an,sum*w(i))
   enddo

   dbnorm = an

end function dbnorm
