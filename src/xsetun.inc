!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### Name
!!   xsetun(3f) - [M_odepack::matrix] Reset the logical unit number for error messages.
!!
!!### Synopsis
!!    subroutine xsetun (lun)
!!    integer,intent(in) :: lun
!!
!!### Description
!!   XSETUN sets the logical unit number for error messages to LUN.
!!
!!### Options
!!
!!   LUN
!!
!!   :   logical unit number used for error messages
!!
!-------------------------------------------------------------------------
! ### CATEGORY  R3B
! ### TYPE      ALL (XSETUN-A)
! ### KEYWORDS  ERROR CONTROL
! ### DESCRIPTION
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### SEE ALSO  XERRWD, XERRWV
! ### REFERENCES  (NONE)
! ### ROUTINES CALLED  IXSAV
! ### REVISION HISTORY  (YYMMDD)
!     19921118  DATE WRITTEN
!     19930329  Added SLATEC format prologue. (FNF)
!     19930407  Corrected SEE ALSO section. (FNF)
!     19930922  Made user-callable, and other cosmetic changes. (FNF)
!-------------------------------------------------------------------------
!  Subroutines called by XSETUN.. None
!  Function routine called by XSETUN ... IXSAV
!-------------------------------------------------------------------------

subroutine xsetun (lun)

integer,intent(in) :: lun

integer :: junk

   if (lun > 0) junk = ixsav ( 1, lun, .true. )

end subroutine xsetun
