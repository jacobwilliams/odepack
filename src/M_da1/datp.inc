!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine computes the product
!!
!!              w = (I - hl0*df/dy)*p
!!
!! This is computed by a call to F and a difference quotient.
!!
!!### On entry
!!
!!   NEQ
!!
!!   : problem size, passed to F and PSOL (NEQ(1) = N).
!!
!!   Y
!!
!!   : array containing current dependent variable vector.
!!
!!   SAVF
!!
!!   : array containing current value of f(t,y).
!!
!!   P
!!
!!   : real array of length N.
!!
!!   WGHT
!!
!!   : array of length N containing scale factors.
!!   1/WGHT(i) are the diagonal elements of the matrix D.
!!
!!   WK
!!
!!   : work array of length N.
!!
!!### On return
!!
!!   W
!!
!!   : array of length N containing desired
!!   matrix-vector product.
!!
!! In addition, this routine uses the Common variables TN, N, NFE.
!
!-----------------------------------------------------------------------
subroutine datp(Neq,Y,Savf,P,Wght,Hl0,Wk,f,W)

integer                     :: Neq(*)
real(kind=dp)               :: Y(*)
real(kind=dp),intent(in)    :: Savf(*)
real(kind=dp)               :: P(*)
real(kind=dp)               :: Wght(*)
real(kind=dp),intent(in)    :: Hl0
real(kind=dp)               :: Wk(*)
external                    :: f
real(kind=dp),intent(inout) :: W(*)

real(kind=dp) :: fac, pnrm, rpnrm
integer :: i

   pnrm = dvnorm(dls1%n,P,Wght)
   rpnrm = 1.0D0/pnrm
   !X!call dcopy(dls1%n,Y,1,W,1)
   W(1:dls1%n) = Y(1:dls1%n)!X!

   do i = 1, dls1%n
      Y(i) = W(i) + P(i)*rpnrm
   enddo

   call f(Neq,dls1%tn,Y,Wk)

   dls1%nfe = dls1%nfe + 1
   !X!call dcopy(dls1%n,W,1,Y,1)
   Y(1:dls1%n) = W(1:dls1%n)!X!
   fac = Hl0*pnrm

   do i = 1, dls1%n
      W(i) = P(i) - fac*(Wk(i)-Savf(i))
   enddo

end subroutine datp
