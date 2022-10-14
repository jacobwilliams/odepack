!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine saves or restores (depending on JOB) the contents of
!! the Common blocks DLS001, type(dlsa01)::DLSA, DLSR01, which are used
!! internally by one or more ODEPACK solvers.
!!
!! RSAV
!!
!! : real array of length 245 or more.
!!
!! ISAV
!!
!! : integer array of length 55 or more.
!!
!! JOB
!!
!! : flag indicating to save or restore the Common blocks:
!!
!!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!!
!!  A call with JOB = 2 presumes a prior call with JOB = 1.
!!
!-----------------------------------------------------------------------
subroutine dsrcar(Rsav,Isav,Job)
!
integer,parameter  ::  LENRLS = 218, LENILS = 37
!
real(kind=dp),intent(inout) :: Rsav(*)
integer,intent(inout)       :: Isav(*)
integer,intent(in)          :: Job
!
integer :: ioff
!
   select case(Job)
   case(1)

      Rsav(1:LENRLS) = return_dls1_real()

      Rsav(LENRLS+1) = dlsa%tsw
      Rsav(LENRLS+2:LENRLS+13) = dlsa%cm1(1:12)
      Rsav(LENRLS+14:LENRLS+18) = dlsa%cm2(1:5)
      Rsav(LENRLS+19) = dlsa%pdest
      Rsav(LENRLS+20) = dlsa%pdlast
      Rsav(LENRLS+21) = dlsa%ratio
      Rsav(LENRLS+22) = dlsa%pdnorm

      Isav(1:LENILS) = return_dls1_int()

      Isav(LENILS+1) = dlsa%insufr
      Isav(LENILS+2) = dlsa%insufi
      Isav(LENILS+3) = dlsa%ixpr
      Isav(LENILS+4) = dlsa%icount
      Isav(LENILS+5) = dlsa%irflag
      Isav(LENILS+6) = dlsa%jtyp
      Isav(LENILS+7) = dlsa%mused
      Isav(LENILS+8) = dlsa%mxordn
      Isav(LENILS+9) = dlsa%mxords

      ioff = LENRLS + 22
      Rsav(ioff+1) = dlsr%alpha
      Rsav(ioff+2) = dlsr%x2
      Rsav(ioff+3) = dlsr%t0
      Rsav(ioff+4) = dlsr%tlast
      Rsav(ioff+5) = dlsr%toutc

      ioff = LENILS + 9
      Isav(ioff+1) = dlsr%lg0
      Isav(ioff+2) = dlsr%lg1
      Isav(ioff+3) = dlsr%lgx
      Isav(ioff+4) = dlsr%imax
      Isav(ioff+5) = dlsr%last
      Isav(ioff+6) = dlsr%irfnd
      Isav(ioff+7) = dlsr%itaskc
      Isav(ioff+8) = dlsr%ngc
      Isav(ioff+9) = dlsr%nge

   case(2)

      call set_dls1_real(Rsav(1:LENRLS))

      dlsa%tsw = Rsav(LENRLS+1)
      dlsa%cm1(1:12) = Rsav(LENRLS+2:LENRLS+13)
      dlsa%cm2(1:5) = Rsav(LENRLS+14:LENRLS+18)
      dlsa%pdest = Rsav(LENRLS+19)
      dlsa%pdlast = Rsav(LENRLS+20)
      dlsa%ratio = Rsav(LENRLS+21)
      dlsa%pdnorm = Rsav(LENRLS+22)

      dlsa%insufr = Isav(LENILS+1)
      dlsa%insufi = Isav(LENILS+2)
      dlsa%ixpr = Isav(LENILS+3)
      dlsa%icount = Isav(LENILS+4)
      dlsa%irflag = Isav(LENILS+5)
      dlsa%jtyp = Isav(LENILS+6)
      dlsa%mused = Isav(LENILS+7)
      dlsa%mxordn = Isav(LENILS+8)
      dlsa%mxords = Isav(LENILS+9)

      ioff = LENRLS + 22
      dlsr%alpha = Rsav(ioff+1)
      dlsr%x2 = Rsav(ioff+2)
      dlsr%t0 = Rsav(ioff+3)
      dlsr%tlast = Rsav(ioff+4)
      dlsr%toutc = Rsav(ioff+5)

      ioff = LENILS + 9
      dlsr%lg0 = Isav(ioff+1)
      dlsr%lg1 = Isav(ioff+2)
      dlsr%lgx = Isav(ioff+3)
      dlsr%imax = Isav(ioff+4)
      dlsr%last = Isav(ioff+5)
      dlsr%irfnd = Isav(ioff+6)
      dlsr%itaskc = Isav(ioff+7)
      dlsr%ngc = Isav(ioff+8)
      dlsr%nge = Isav(ioff+9)

      call set_dls1_int(Isav(1:LENILS))

   case default

      write(*,*)'<ERROR>*dsrcar* unknown value for JOB=',job
      stop 1

   end select

end subroutine dsrcar
