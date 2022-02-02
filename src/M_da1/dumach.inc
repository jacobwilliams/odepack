!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dumach(3f) - [M_odepack] Compute the unit roundoff of the machine.
!!
!!### SYNOPSIS
!!        a = dumach()
!!        double precision :: dumach
!!        double precision :: a
!!
!!### DESCRIPTION
!!   The unit roundoff is defined as the smallest positive machine
!!   number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!!   in a machine-independent manner.
!!
!!### RETURNS
!!   A
!!
!!   : the unit roundoff of the machine.
!!
!-----------------------------------------------------------------------
! ### BEGIN PROLOGUE  DUMACH
! ### PURPOSE  Compute the unit roundoff of the machine.
! ### CATEGORY  R1
! ### TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
! ### KEYWORDS  MACHINE CONSTANTS
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### REFERENCES  (NONE)
! ### ROUTINES CALLED  DUMSUM
! ### REVISION HISTORY  (YYYYMMDD)
!     19930216 DATE WRITTEN
!     19930818 Added SLATEC-format prologue.  (FNF)
!     20030707 Added DUMSUM to force normal storage of COMP.  (ACH)
!     20220122 f90+ syntax (JSU)
!     TODO: replace with EPSILON() and cached value
! ### END PROLOGUE  DUMACH
!-----------------------------------------------------------------------
function dumach()
real(kind=dp) :: dumach
real(kind=dp) :: u, comp
    u = 1.0d0
    do
        u = u * 0.5d0
        call dumsum( 1.0d0, u, comp )
        if ( comp == 1.0d0 ) then
           dumach = u * 2.0d0
           exit
        endif
    enddo
contains

subroutine dumsum(a,b,c)
!  Routine to force normal storing of A + B, for DUMACH.
real(kind=dp),intent(in)  :: a, b
real(kind=dp),intent(out) :: c
    c = a + b
end subroutine dumsum

end function dumach
