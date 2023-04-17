!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!! dsrcom(3f) - [M_odepack] Save/restore ODEPACK COMMON blocks.
!!
!!### DESCRIPTION
!!  This routine saves or restores (depending on JOB) the contents of
!!  the internal types used to store the current state
!!  by one (or more) of the ODEPACK solvers.
!!
!!### VALUES
!!
!!  RSAV
!!
!!  : real array of length 218 or more.
!!
!!  ISAV
!!
!!  : integer array of length 37 or more.
!!
!!  JOB
!!
!!  : flag indicating to save or restore the COMMON blocks:
!!              JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
!!              JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
!!  A call with JOB = 2 presumes a prior call with JOB = 1.
!!
!-----------------------------------------------------------------------
! ### BEGIN PROLOGUE  DSRCOM
! ### SUBSIDIARY
! ### PURPOSE  Save/restore ODEPACK COMMON blocks.
! ### TYPE      DOUBLE PRECISION (SSRCOM-S, DSRCOM-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### DESCRIPTION
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  (NONE)
! ### COMMON BLOCKS    DLS001
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890503  Minor cosmetic changes.  (FNF)
!     19921116  Deleted treatment of block /EH0001/.  (ACH)
!     19930801  Reduced Common block length by 2.  (ACH)
!     19930809  Renamed to allow single/double precision versions. (ACH)
!     20010418  Reduced Common block length by 209+12. (ACH)
!     20031105  Restored 'own' variables to Common block /DLS001/, to
!               enable interrupt/restart feature. (ACH)
!     20031112  Added SAVE statement for data-loaded constants.
! ### END PROLOGUE  DSRCOM
!-----------------------------------------------------------------------
subroutine dsrcom (rsav, isav, job)
real(kind=dp)      :: Rsav(*)
integer            :: Isav(*)
integer,intent(in) :: Job

integer,parameter :: lenils=37, lenrls=218

   select case(job)

   case(1)

     rsav(1:lenrls)=return_dls1_real()
     isav(1:lenils)=return_dls1_int()

   case(2)

     call set_dls1_real(rsav(1:lenrls))
     call set_dls1_int(isav(1:lenils))

   case default

      write(*,*)'<ERROR>*dsrcom* unknown JOB value=',job
      stop 1

   end select
end subroutine dsrcom
