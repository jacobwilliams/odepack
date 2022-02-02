!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DPKSET is called by DSTODPK to interface with the user-supplied
!! routine JAC, to compute and process relevant parts of
!! the matrix P = I - H*EL(1)*J, where J is the Jacobian df/dy,
!! as need for preconditioning matrix operations later.
!!
!! In addition to variables described previously, communication
!! with DPKSET uses the following:
!!
!! Y
!!
!! : array containing predicted values on entry.
!!
!! YSV
!!
!! : array containing predicted y, to be saved (YH1 in DSTODPK).
!!
!! FTEM
!!
!! : work array of length N (ACOR in DSTODPK).
!!
!! SAVF
!!
!! : array containing f evaluated at predicted y.
!!
!! WM
!!
!! : real work space for matrices.
!! Space for preconditioning data starts at WM(LOCWP).
!!
!! IWM
!!
!! : integer work space.
!! Space for preconditioning data starts at IWM(LOCIWP).
!!
!! IERPJ
!!
!! : output error flag,  = 0 if no trouble, .gt. 0 if
!! JAC returned an error flag.
!!
!! JCUR
!!
!! : output flag = 1 to indicate that the Jacobian matrix
!! (or approximation) is now current.
!!
!! This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
!-----------------------------------------------------------------------
subroutine dpkset(Neq,Y,Ysv,Ewt,Ftem,Savf,Wm,Iwm,f,jac)

integer       :: Neq(*)
real(kind=dp) :: Y(*)
real(kind=dp) :: Ysv(*)
real(kind=dp) :: Ewt(*)
real(kind=dp) :: Ftem(*)
real(kind=dp) :: Savf(*)
real(kind=dp) :: Wm(*)
integer       :: Iwm(*)
external      :: f
external      :: jac

real(kind=dp) :: hl0
integer       :: ier

   dls1%ierpj = 0
   dls1%jcur = 1
   hl0 = dls1%el0*dls1%h
   call jac(f,Neq,dls1%tn,Y,Ysv,Ewt,Savf,Ftem,hl0,Wm(dlpk%locwp),Iwm(dlpk%lociwp),ier)
   dls1%nje = dls1%nje + 1

   if ( ier /= 0 )then
      dls1%ierpj = 1
   endif

end subroutine dpkset
