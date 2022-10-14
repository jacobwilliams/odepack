!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine saves or restores (depending on JOB) the contents of
!! the Common blocks DLS001, type(DLSA01)::DLSA, which are used
!! internally by one or more ODEPACK solvers.
!!
!! RSAV
!!
!! : real array of length 240 or more.
!!
!! ISAV
!!
!! : integer array of length 46 or more.
!!
!! JOB
!!
!! : flag indicating to save or restore the Common blocks:
!!
!!        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
!!        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
!!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
subroutine dsrcma(Rsav,Isav,Job)
!
integer,parameter  ::  LENRLS = 218 , LENILS = 37
!
real(kind=dp),intent(inout) :: Rsav(*)
integer,intent(inout)       :: Isav(*)
integer,intent(in)          :: Job
!
select case (Job)
case (1)

   Rsav(1:LENRLS) = return_dls1_real()

   Rsav(LENRLS+1)             =  dlsa%tsw
   Rsav(LENRLS+2:LENRLS+13)   =  dlsa%cm1(1:12)
   Rsav(LENRLS+14:LENRLS+18)  =  dlsa%cm2(1:5)
   Rsav(LENRLS+19)            =  dlsa%pdest
   Rsav(LENRLS+20)            =  dlsa%pdlast
   Rsav(LENRLS+21)            =  dlsa%ratio
   Rsav(LENRLS+22)            =  dlsa%pdnorm

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

case (2)

   call set_dls1_real(Rsav(1:LENRLS))

   dlsa%tsw        =  Rsav(LENRLS+1)
   dlsa%cm1(1:12)  =  Rsav(LENRLS+2:LENRLS+13)
   dlsa%cm2(1:5)   =  Rsav(LENRLS+14:LENRLS+18)
   dlsa%pdest      =  Rsav(LENRLS+19)
   dlsa%pdlast     =  Rsav(LENRLS+20)
   dlsa%ratio      =  Rsav(LENRLS+21)
   dlsa%pdnorm     =  Rsav(LENRLS+22)

   call set_dls1_int(Isav(1:LENILS))

   dlsa%insufr  =  Isav(LENILS+1)
   dlsa%insufi  =  Isav(LENILS+2)
   dlsa%ixpr    =  Isav(LENILS+3)
   dlsa%icount  =  Isav(LENILS+4)
   dlsa%irflag  =  Isav(LENILS+5)
   dlsa%jtyp    =  Isav(LENILS+6)
   dlsa%mused   =  Isav(LENILS+7)
   dlsa%mxordn  =  Isav(LENILS+8)
   dlsa%mxords  =  Isav(LENILS+9)

case default

   write (*,*) '<ERROR>*dsrcma* unknown value for JOB=' , Job
   stop 1

endselect

end subroutine dsrcma
