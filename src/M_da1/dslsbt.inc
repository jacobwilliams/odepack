!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine acts as an interface between the core integrator
!! routine and the DSOLBT routine for the solution of the linear system
!! arising from chord iteration.
!! Communication with DSLSBT uses the following variables:
!!
!! WM
!!
!! : real work space containing the LU decomposition,
!! starting at WM(3).
!!
!! IWM
!!
!! : integer work space containing pivot information, starting at
!! IWM(21).  IWM also contains block structure parameters
!! MB = IWM(1) and NB = IWM(2).
!!
!! X
!!
!! : the right-hand side vector on input, and the solution vector
!! on output, of length N.
!!
!! TEM
!!
!! : vector of work space of length N, not used in this version.
!-----------------------------------------------------------------------
subroutine dslsbt(Wm,Iwm,X,Tem)
!
real(kind=dp) :: Wm(*)
integer       :: Iwm(*)
real(kind=dp) :: X(*)
real(kind=dp) :: Tem(*)
!
integer :: lblox, lpb, lpc, mb, nb
!
   mb = Iwm(1)
   nb = Iwm(2)
   lblox = mb*mb*nb
   lpb = 3 + lblox
   lpc = lpb + lblox
   call dsolbt(mb,nb,Wm(3),Wm(lpb),Wm(lpc),X,Iwm(21))

end subroutine dslsbt
