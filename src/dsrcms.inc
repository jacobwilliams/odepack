!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine saves or restores (depending on JOB) the contents of
!! the Common blocks DLS001, DLSS01, which are used
!! internally by one or more ODEPACK solvers.
!!
!! RSAV
!!
!! : real array of length 224 or more.
!!
!! ISAV
!!
!! : integer array of length 71 or more.
!!
!! JOB
!!
!! : flag indicating to save or restore the Common blocks:
!!
!!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
subroutine dsrcms (rsav, isav, job)
real(kind=dp),intent(inout) :: Rsav(*)
integer,intent(inout)       :: Isav(*)
integer,intent(in)          :: Job

integer,parameter ::  lenils=37, lenrls=218

   select case(job)

   case(1)

      rsav(1:lenrls) = return_dls1_real()

      rsav( 1 + lenrls) = dlss%con0
      rsav( 2 + lenrls) = dlss%conmin
      rsav( 3 + lenrls) = dlss%ccmxj
      rsav( 4 + lenrls) = dlss%psmall
      rsav( 5 + lenrls) = dlss%rbig
      rsav( 6 + lenrls) = dlss%seth

      isav(1:lenils) = return_dls1_int()

      isav( 1 + lenils) = dlss%iplost
      isav( 2 + lenils) = dlss%iesp
      isav( 3 + lenils) = dlss%istatc
      isav( 4 + lenils) = dlss%iys
      isav( 5 + lenils) = dlss%iba
      isav( 6 + lenils) = dlss%ibian
      isav( 7 + lenils) = dlss%ibjan
      isav( 8 + lenils) = dlss%ibjgp
      isav( 9 + lenils) = dlss%ipian
      isav(10 + lenils) = dlss%ipjan
      isav(11 + lenils) = dlss%ipjgp
      isav(12 + lenils) = dlss%ipigp
      isav(13 + lenils) = dlss%ipr
      isav(14 + lenils) = dlss%ipc
      isav(15 + lenils) = dlss%ipic
      isav(16 + lenils) = dlss%ipisp
      isav(17 + lenils) = dlss%iprsp
      isav(18 + lenils) = dlss%ipa
      isav(19 + lenils) = dlss%lenyh
      isav(20 + lenils) = dlss%lenyhm
      isav(21 + lenils) = dlss%lenwk
      isav(22 + lenils) = dlss%lreq
      isav(23 + lenils) = dlss%lrat
      isav(24 + lenils) = dlss%lrest
      isav(25 + lenils) = dlss%lwmin
      isav(26 + lenils) = dlss%moss
      isav(27 + lenils) = dlss%msbj
      isav(28 + lenils) = dlss%nslj
      isav(29 + lenils) = dlss%ngp
      isav(30 + lenils) = dlss%nlu
      isav(31 + lenils) = dlss%nnz
      isav(32 + lenils) = dlss%nsp
      isav(33 + lenils) = dlss%nzl
      isav(34 + lenils) = dlss%nzu

   case(2)
      call set_dls1_real(rsav(1:lenrls))

      dlss%con0    =  rsav( 1 + lenrls)
      dlss%conmin  =  rsav( 2 + lenrls)
      dlss%ccmxj   =  rsav( 3 + lenrls)
      dlss%psmall  =  rsav( 4 + lenrls)
      dlss%rbig    =  rsav( 5 + lenrls)
      dlss%seth    =  rsav( 6 + lenrls)

      call set_dls1_int(isav(1:lenils))

      dlss%iplost  =  isav( 1 + lenils)
      dlss%iesp    =  isav( 2 + lenils)
      dlss%istatc  =  isav( 3 + lenils)
      dlss%iys     =  isav( 4 + lenils)
      dlss%iba     =  isav( 5 + lenils)
      dlss%ibian   =  isav( 6 + lenils)
      dlss%ibjan   =  isav( 7 + lenils)
      dlss%ibjgp   =  isav( 8 + lenils)
      dlss%ipian   =  isav( 9 + lenils)
      dlss%ipjan   =  isav(10 + lenils)
      dlss%ipjgp   =  isav(11 + lenils)
      dlss%ipigp   =  isav(12 + lenils)
      dlss%ipr     =  isav(13 + lenils)
      dlss%ipc     =  isav(14 + lenils)
      dlss%ipic    =  isav(15 + lenils)
      dlss%ipisp   =  isav(16 + lenils)
      dlss%iprsp   =  isav(17 + lenils)
      dlss%ipa     =  isav(18 + lenils)
      dlss%lenyh   =  isav(19 + lenils)
      dlss%lenyhm  =  isav(20 + lenils)
      dlss%lenwk   =  isav(21 + lenils)
      dlss%lreq    =  isav(22 + lenils)
      dlss%lrat    =  isav(23 + lenils)
      dlss%lrest   =  isav(24 + lenils)
      dlss%lwmin   =  isav(25 + lenils)
      dlss%moss    =  isav(26 + lenils)
      dlss%msbj    =  isav(27 + lenils)
      dlss%nslj    =  isav(28 + lenils)
      dlss%ngp     =  isav(29 + lenils)
      dlss%nlu     =  isav(30 + lenils)
      dlss%nnz     =  isav(31 + lenils)
      dlss%nsp     =  isav(32 + lenils)
      dlss%nzl     =  isav(33 + lenils)
      dlss%nzu     =  isav(34 + lenils)
   case default

      write(*,*)'<ERROR>*dsrcms* unknown value for JOB=',job
      stop 1

   end select
end subroutine dsrcms
