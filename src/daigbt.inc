!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This subroutine computes the initial value
!! of the vector YDOT satisfying
!!
!!      A * YDOT = g(t,y)
!!
!! when A is nonsingular.  It is called by DLSOIBT for
!! initialization only, when ISTATE = 0 .
!! DAIGBT returns an error flag IER:
!!
!!     IER  =  0  means DAIGBT was successful.
!!     IER .ge. 2 means RES returned an error flag IRES = IER.
!!     IER .lt. 0 means the A matrix was found to have a singular
!!                diagonal block (hence YDOT could not be solved for).
!!
!-----------------------------------------------------------------------
subroutine daigbt(res,adda,Neq,T,Y,Ydot,Mb,Nb,Pw,Ipvt,Ier)
!
external res, adda
integer               :: Neq(*)
real(kind=dp)         :: T
real(kind=dp)         :: Y(*)
real(kind=dp)         :: Ydot(*)
integer,intent(inout) :: Mb
integer,intent(inout) :: Nb
real(kind=dp)         :: Pw(*)
integer               :: Ipvt(*)
integer,intent(inout) :: Ier

integer :: i , lblox , lenpw , lpb , lpc
!
   lblox = Mb*Mb*Nb
   lpb = 1 + lblox
   lpc = lpb + lblox
   lenpw = 3*lblox

   do i = 1 , lenpw
      Pw(i) = 0.0D0
   enddo

   Ier = 1
   call res(Neq,T,Y,Pw,Ydot,Ier)
   if ( Ier>1 ) return
   call adda(Neq,T,Y,Mb,Nb,Pw(1),Pw(lpb),Pw(lpc))
   call ddecbt(Mb,Nb,Pw,Pw(lpb),Pw(lpc),Ipvt,Ier)
   if ( Ier==0 ) then
      call dsolbt(Mb,Nb,Pw,Pw(lpb),Pw(lpc),Ydot,Ipvt)
      return
   endif

   Ier = -Ier
end subroutine daigbt
