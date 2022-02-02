!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DSETPK is called by DSTOKA to interface with the user-supplied
!! routine JAC, to compute and process relevant parts of
!! the matrix P = I - H*EL(1)*J, where J is the Jacobian df/dy,
!! as need for preconditioning matrix operations later.
!!
!! In addition to variables described previously, communication
!! with DSETPK uses the following:
!!
!! Y
!!
!! : array containing predicted values on entry.
!!
!! YSV
!!
!! : array containing predicted y, to be saved (YH1 in DSTOKA).
!!
!! FTEM
!!
!! : work array of length N (ACOR in DSTOKA).
!!
!! SAVF
!!
!! : array containing f evaluated at predicted y.
!!
!! JOK
!!
!! : input flag showing whether it was judged that Jacobian matrix
!! data need not be recomputed (JOK = 1) or needs to be (JOK = -1).
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
!! : output flag to indicate whether the matrix data involved
!! is now current (JCUR = 1) or not (JCUR = 0).
!!
!! This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
!-----------------------------------------------------------------------
subroutine dsetpk(Neq,Y,Ysv,Ewt,Ftem,Savf,Jok,Wm,Iwm,f,jac)

integer       :: Neq(*)
real(kind=dp) :: Y(*)
real(kind=dp) :: Ysv(*)
real(kind=dp) :: Ewt(*)
real(kind=dp) :: Ftem(*)
real(kind=dp) :: Savf(*)
integer       :: Jok
real(kind=dp) :: Wm(*)
integer       :: Iwm(*)
external      :: f
external      :: jac

real(kind=dp) :: hl0
integer       :: ier

   dls1%ierpj = 0
   dls1%jcur = 0

   if ( Jok==-1 ) dls1%jcur = 1

   hl0 = dls1%el0*dls1%h
   call jac(f,Neq,dls1%tn,Y,Ysv,Ewt,Savf,Ftem,hl0,Jok,Wm(dlpk%locwp),Iwm(dlpk%lociwp),ier)
   dls1%nje = dls1%nje + 1

   if ( ier==0 ) return

   dls1%ierpj = 1

end subroutine dsetpk
