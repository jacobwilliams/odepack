!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   xsetf(3f) - [M_odepack::matrix] Reset the error print control flag.
!!
!!### SYNOPSIS
!!    subroutine xsetf (mflag)
!!    integer mflag
!!
!!### DESCRIPTION
!!   XSETF sets the error print control flag to MFLAG:
!!### OPTIONS
!!
!!   MFLAG
!!
!!   : MFLAG=1 means print all messages (the default).
!!
!!   MFLAG=0 means no printing.
!!
! ### BEGIN PROLOGUE  XSETF
! ### CATEGORY  R3A
! ### TYPE      ALL (XSETF-A)
! ### KEYWORDS  ERROR CONTROL
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### DESCRIPTION
! ### SEE ALSO  XERRWD, XERRWV
! ### REFERENCES  (NONE)
! ### ROUTINES CALLED  IXSAV
! ### REVISION HISTORY  (YYMMDD)
!     19921118  DATE WRITTEN
!     19930329  Added SLATEC format prologue. (FNF)
!     19930407  Corrected SEE ALSO section. (FNF)
!     19930922  Made user-callable, and other cosmetic changes. (FNF)
!
!  Subroutines called by XSETF.. None
!  Function routine called by XSETF.. IXSAV

subroutine xsetf (mflag)
integer mflag, junk
!!
   if (mflag .eq. 0 .or. mflag .eq. 1) junk = ixsav (2,mflag,.true.)
end subroutine xsetf
