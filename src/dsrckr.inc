!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine saves or restores (depending on JOB) the contents of
!! the Common blocks DLS001, DLSR01, DLPK01, which
!! are used internally by the DLSODKR solver.
!!
!! RSAV
!!
!! : real array of length 228 or more.
!!
!! ISAV
!!
!! : integer array of length 63 or more.
!!
!! JOB
!!
!! : flag indicating to save or restore the Common blocks:
!!
!!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!!
!! A call with JOB = 2 presumes a prior call with JOB = 1.
!!
!-----------------------------------------------------------------------
subroutine dsrckr(Rsav,Isav,Job)
!
integer, parameter  ::  LENRLS = 218, LENILS = 37
!
real(kind=dp),intent(inout) :: Rsav(*)
integer,intent(inout)       :: Isav(*)
integer,intent(in)          :: Job
!
integer :: ioff
!
select case (Job)
case (1)
   Rsav(1:LENRLS) = return_dls1_real()

   Rsav(LENRLS+1) = dls%stifr

   Rsav(LENRLS+1+1) = dlsr%alpha
   Rsav(LENRLS+2+1) = dlsr%x2
   Rsav(LENRLS+3+1) = dlsr%t0
   Rsav(LENRLS+4+1) = dlsr%tlast
   Rsav(LENRLS+5+1) = dlsr%toutc

   Rsav(LENRLS+5+2) = dlpk%delt
   Rsav(LENRLS+5+3) = dlpk%epcon
   Rsav(LENRLS+5+4) = dlpk%sqrtn
   Rsav(LENRLS+5+5) = dlpk%rsqrtn

   Isav(1:LENILS) = return_dls1_int()

   Isav(LENILS+1) = dls%newt
   Isav(LENILS+2) = dls%nsfi
   Isav(LENILS+3) = dls%nslj
   Isav(LENILS+4) = dls%njev
   ioff = LENILS + 2

   Isav(ioff+1) = dlsr%lg0
   Isav(ioff+2) = dlsr%lg1
   Isav(ioff+3) = dlsr%lgx
   Isav(ioff+4) = dlsr%imax
   Isav(ioff+5) = dlsr%last
   Isav(ioff+6) = dlsr%irfnd
   Isav(ioff+7) = dlsr%itaskc
   Isav(ioff+8) = dlsr%ngc
   Isav(ioff+9) = dlsr%nge

   ioff = ioff + 9

   Isav(ioff+1) = dlpk%jpre
   Isav(ioff+2) = dlpk%jacflg
   Isav(ioff+3) = dlpk%locwp
   Isav(ioff+4) = dlpk%lociwp
   Isav(ioff+5) = dlpk%lsavx
   Isav(ioff+6) = dlpk%kmp
   Isav(ioff+7) = dlpk%maxl
   Isav(ioff+8) = dlpk%mnewt
   Isav(ioff+9) = dlpk%nni
   Isav(ioff+10) = dlpk%nli
   Isav(ioff+11) = dlpk%nps
   Isav(ioff+12) = dlpk%ncfn
   Isav(ioff+13) = dlpk%ncfl

case (2)
   call set_dls1_real(Rsav(1:LENRLS))

   dls%stifr = Rsav(LENRLS+1)

   dlsr%alpha  =  Rsav(LENRLS+1+1)
   dlsr%x2     =  Rsav(LENRLS+2+1)
   dlsr%t0     =  Rsav(LENRLS+3+1)
   dlsr%tlast  =  Rsav(LENRLS+4+1)
   dlsr%toutc  =  Rsav(LENRLS+5+1)

   dlpk%delt    =  Rsav(LENRLS+5+2)
   dlpk%epcon   =  Rsav(LENRLS+5+3)
   dlpk%sqrtn   =  Rsav(LENRLS+5+4)
   dlpk%rsqrtn  =  Rsav(LENRLS+5+5)

   call set_dls1_int(Isav(1:LENILS))

   dls%newt = Isav(LENILS+1)
   dls%nsfi = Isav(LENILS+2)
   dls%nslj = Isav(LENILS+3)
   dls%njev = Isav(LENILS+4)

   ioff = LENILS + 2

   dlsr%lg0     =  Isav(ioff+1)
   dlsr%lg1     =  Isav(ioff+2)
   dlsr%lgx     =  Isav(ioff+3)
   dlsr%imax    =  Isav(ioff+4)
   dlsr%last    =  Isav(ioff+5)
   dlsr%irfnd   =  Isav(ioff+6)
   dlsr%itaskc  =  Isav(ioff+7)
   dlsr%ngc     =  Isav(ioff+8)
   dlsr%nge     =  Isav(ioff+9)

   ioff = ioff + 9

   dlpk%jpre    =  Isav(ioff+1)
   dlpk%jacflg  =  Isav(ioff+2)
   dlpk%locwp   =  Isav(ioff+3)
   dlpk%lociwp  =  Isav(ioff+4)
   dlpk%lsavx   =  Isav(ioff+5)
   dlpk%kmp     =  Isav(ioff+6)
   dlpk%maxl    =  Isav(ioff+7)
   dlpk%mnewt   =  Isav(ioff+8)
   dlpk%nni     =  Isav(ioff+9)
   dlpk%nli     =  Isav(ioff+10)
   dlpk%nps     =  Isav(ioff+11)
   dlpk%ncfn    =  Isav(ioff+12)
   dlpk%ncfl    =  Isav(ioff+13)

case default

   write (*,*) '<ERROR>*dsrckr* unknown value for JOB=', Job
   stop 1

endselect

end subroutine dsrckr
