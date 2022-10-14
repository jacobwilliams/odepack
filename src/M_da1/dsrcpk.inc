!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine saves or restores (depending on JOB) the contents of
!! the internal types used internally by the DLSODPK solver.
!!
!! RSAV
!!
!! : real array of length 222 or more.
!!
!! ISAV
!!
!! : integer array of length 50 or more.
!!
!! JOB
!!
!! : flag indicating to save or restore the Common blocks:
!!
!!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!!        A call with JOB = 2 presumes a prior call with JOB = 1.
!!
!-----------------------------------------------------------------------
subroutine dsrcpk (rsav, isav, job)

real(kind=dp)      :: Rsav(*)
integer            :: Isav(*)
integer,intent(in) :: Job

integer,parameter :: lenils=37, lenrls=218

select case(job)
case(1)

      rsav(1:lenrls) = return_dls1_real()

      rsav(lenrls+1)   =  dlpk%delt
      rsav(lenrls+2)   =  dlpk%epcon
      rsav(lenrls+3)   =  dlpk%sqrtn
      rsav(lenrls+4)   =  dlpk%rsqrtn

      isav(1:lenils) = return_dls1_int()

      isav(lenils+1)   =  dlpk%jpre
      isav(lenils+2)   =  dlpk%jacflg
      isav(lenils+3)   =  dlpk%locwp
      isav(lenils+4)   =  dlpk%lociwp
      isav(lenils+5)   =  dlpk%lsavx
      isav(lenils+6)   =  dlpk%kmp
      isav(lenils+7)   =  dlpk%maxl
      isav(lenils+8)   =  dlpk%mnewt
      isav(lenils+9)   =  dlpk%nni
      isav(lenils+10)  =  dlpk%nli
      isav(lenils+11)  =  dlpk%nps
      isav(lenils+12)  =  dlpk%ncfn
      isav(lenils+13)  =  dlpk%ncfl

case(2)

      call set_dls1_real(rsav(1:lenrls))

      dlpk%delt    =  rsav(lenrls+1)
      dlpk%epcon   =  rsav(lenrls+2)
      dlpk%sqrtn   =  rsav(lenrls+3)
      dlpk%rsqrtn  =  rsav(lenrls+4)

      call set_dls1_int(isav(1:lenils))

      dlpk%jpre    =  isav(lenils+1)
      dlpk%jacflg  =  isav(lenils+2)
      dlpk%locwp   =  isav(lenils+3)
      dlpk%lociwp  =  isav(lenils+4)
      dlpk%lsavx   =  isav(lenils+5)
      dlpk%kmp     =  isav(lenils+6)
      dlpk%maxl    =  isav(lenils+7)
      dlpk%mnewt   =  isav(lenils+8)
      dlpk%nni     =  isav(lenils+9)
      dlpk%nli     =  isav(lenils+10)
      dlpk%nps     =  isav(lenils+11)
      dlpk%ncfn    =  isav(lenils+12)
      dlpk%ncfl    =  isav(lenils+13)

case default

   write (*,*) '<ERROR>*dsrcpk* unknown value for JOB=', Job
   stop 1

end select

end subroutine dsrcpk
