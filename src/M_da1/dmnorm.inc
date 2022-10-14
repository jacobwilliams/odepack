!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This function routine computes the weighted max-norm
!! of the vector of length N contained in the array V, with weights
!! contained in the array w of length N:
!!   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
!-----------------------------------------------------------------------
pure function dmnorm (n, v, w)
integer,intent(in)       :: n
real(kind=dp),intent(in) :: v(n)
real(kind=dp),intent(in) :: w(n)
real(kind=dp) :: dmnorm

integer       :: i
real(kind=dp) :: vm

   vm = 0.0d0

   do i = 1,n
      vm = max(vm,abs(v(i))*w(i))
   enddo

   dmnorm = vm
end function dmnorm
