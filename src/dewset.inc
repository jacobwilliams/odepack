!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dewset(3f) - [M_odepack] Set error weight vector.
!!
!!### SYNOPSIS
!!        subroutine dewset(N,Itol,Rtol,Atol,Ycur,Ewt)
!!        integer,intent(in)        :: N
!!        integer,intent(in)        :: Itol
!!        real(kind=dp),intent(in)  :: Rtol(*)
!!        real(kind=dp),intent(in)  :: Atol(*)
!!        real(kind=dp),intent(in)  :: Ycur(N)
!!        real(kind=dp),intent(out) :: Ewt(N)
!!
!!### DESCRIPTION
!!  This subroutine sets the error weight vector EWT according to
!!
!!       EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!!
!!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!!  depending on the value of ITOL.
!!
! ### SUBSIDIARY
! ### PURPOSE  Set error weight vector.
! ### TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### DESCRIPTION
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890503  Minor cosmetic changes.  (FNF)
!     19930809  Renamed to allow single/double precision versions. (ACH)
!-----------------------------------------------------------------------
subroutine dewset(N,Itol,Rtol,Atol,Ycur,Ewt)
!
integer,intent(in)        :: N
integer,intent(in)        :: Itol
real(kind=dp),intent(in)  :: Rtol(*)
real(kind=dp),intent(in)  :: Atol(*)
real(kind=dp),intent(in)  :: Ycur(N)
real(kind=dp),intent(out) :: Ewt(N)
!
integer :: i
!
select case (Itol)
case (2)
   do i = 1 , N
      Ewt(i) = Rtol(1)*abs(Ycur(i)) + Atol(i)
   enddo
case (3)
   do i = 1 , N
      Ewt(i) = Rtol(i)*abs(Ycur(i)) + Atol(1)
   enddo
case (4)
   do i = 1 , N
      Ewt(i) = Rtol(i)*abs(Ycur(i)) + Atol(i)
   enddo
case default
   do i = 1 , N
      Ewt(i) = Rtol(1)*abs(Ycur(i)) + Atol(1)
   enddo
endselect
end subroutine dewset
