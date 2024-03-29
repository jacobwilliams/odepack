!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This function computes the norm of a full N by N matrix,
!! stored in the array A, that is consistent with the weighted max-norm
!! on vectors, with weights stored in the array W:
!!
!!     DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
!-----------------------------------------------------------------------
pure function dfnorm (n, a, w)
integer,intent(in)       :: n
real(kind=dp),intent(in) :: a(n,n)
real(kind=dp),intent(in) :: w(n)
real(kind=dp)            :: dfnorm

real(kind=dp) an, sum
integer i, j

   an = 0.0d0
   do i = 1,n
      sum = 0.0d0
      do j = 1,n
         sum = sum + abs(a(i,j))/w(j)
      enddo
      an = max(an,sum*w(i))
   enddo
   dfnorm = an
end function dfnorm
