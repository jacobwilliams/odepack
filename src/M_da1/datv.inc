!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine computes the product
!!
!!     (D-inverse)*(P1-inverse)*(I - hl0*df/dy)*(P2-inverse)*(D*v),
!!
!! where D is a diagonal scaling matrix, and P1 and P2 are the
!! left and right preconditioning matrices, respectively.
!! v is assumed to have WRMS norm equal to 1.
!! The product is stored in z.  This is computed by a
!! difference quotient, a call to F, and two calls to PSOL.
!!
!!### On entry
!!
!! NEQ
!!
!! : problem size, passed to F and PSOL (NEQ(1) = N).
!!
!! Y
!!
!! : array containing current dependent variable vector.
!!
!! SAVF
!!
!! : array containing current value of f(t,y).
!!
!! V
!!
!! : real array of length N (can be the same array as Z).
!!
!! WGHT
!!
!! : array of length N containing scale factors.
!! 1/WGHT(i) are the diagonal elements of the matrix D.
!!
!! FTEM
!!
!! : work array of length N.
!!
!! VTEM
!!
!! : work array of length N used to store the
!! unscaled version of V.
!!
!! WP
!!
!! : real work array used by preconditioner PSOL.
!!
!! IWP
!!
!! : integer work array used by preconditioner PSOL.
!!
!! HL0
!!
!! : current value of (step size dls1%h) * (coefficient l0).
!!
!! JPRE
!!
!! : preconditioner type flag.
!!
!!### On return
!!
!! Z
!!
!! : array of length N containing desired scaled
!! matrix-vector product.
!!
!! IER
!!
!! : error flag from PSOL.
!!
!! NPSL
!!
!! : the number of calls to PSOL.
!!
!! In addition, this routine uses the Common variables TN, N, NFE.
!-----------------------------------------------------------------------
subroutine datv(Neq,Y,Savf,V,Wght,Ftem,f,psol,Z,Vtem,Wp,Iwp,Hl0,Jpre,Ier,Npsl)

integer                     :: Neq(*)
real(kind=dp)               :: Y(*)
real(kind=dp)               :: Savf(*)
real(kind=dp),intent(in)    :: V(*)
real(kind=dp),intent(in)    :: Wght(*)
real(kind=dp)               :: Ftem(*)
external                    :: f
external                    :: psol
real(kind=dp),intent(inout) :: Z(*)
real(kind=dp),intent(inout) :: Vtem(*)
real(kind=dp),intent(inout) :: Wp(*)
integer                     :: Iwp(*)
real(kind=dp)               :: Hl0
integer, intent(inout)      :: Jpre
integer, intent(inout)      :: Ier
integer, intent(inout)      :: Npsl

real(kind=dp) :: fac, rnorm, tempn
integer :: i

   !  Set VTEM = D * V.
   do i = 1, dls1%n
      Vtem(i) = V(i)/Wght(i)
   enddo

   Ier = 0

   if ( Jpre>=2 ) then

      !  JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM.
      call psol(Neq,dls1%tn,Y,Savf,Ftem,Hl0,Wp,Iwp,Vtem,2,Ier)
      Npsl = Npsl + 1

      if ( Ier/=0 ) return

      !  Calculate L-2 norm of (D-inverse) * VTEM.
      do i = 1, dls1%n
         Z(i) = Vtem(i)*Wght(i)
      enddo

      tempn = dnrm2(dls1%n,Z,1)
      rnorm = 1.0D0/tempn

      !  Save Y in Z and increment Y by VTEM/norm.
      !X!call dcopy(dls1%n,Y,1,Z,1)
      Z(1:dls1%n)=Y(1:dls1%n)

      do i = 1, dls1%n
         Y(i) = Z(i) + Vtem(i)*rnorm
      enddo

      fac = Hl0*tempn
   else

      !  JPRE = 0 or 1.  Save Y in Z and increment Y by VTEM.
      !X!call dcopy(dls1%n,Y,1,Z,1)
      Z(1:dls1%n)=Y(1:dls1%n)

      do i = 1, dls1%n
         Y(i) = Z(i) + Vtem(i)
      enddo

      fac = Hl0
   endif

   !  For all JPRE, call F with incremented Y argument, and restore Y.
   call f(Neq,dls1%tn,Y,Ftem)
   dls1%nfe = dls1%nfe + 1
   !X!call dcopy(dls1%n,Z,1,Y,1)
   Y(1:dls1%n)=Z(1:dls1%n)

   !  Set Z = (identity - hl0*Jacobian) * VTEM, using difference quotient.
   do i = 1, dls1%n
      Z(i) = Ftem(i) - Savf(i)
   enddo
   do i = 1, dls1%n
      Z(i) = Vtem(i) - fac*Z(i)
   enddo

   !  Apply inverse of left preconditioner to Z, if nontrivial.
   if ( Jpre/=0 .and. Jpre/=2 ) then
      call psol(Neq,dls1%tn,Y,Savf,Ftem,Hl0,Wp,Iwp,Z,1,Ier)
      Npsl = Npsl + 1
      if ( Ier/=0 ) return
   endif

   !  Apply D-inverse to Z and return.
   do i = 1, dls1%n
      Z(i) = Z(i)*Wght(i)
   enddo

end subroutine datv
