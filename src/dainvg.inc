!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This subroutine computes the initial value
!! of the vector YDOT satisfying
!!
!!      A * YDOT = g(t,y)
!!
!! when A is nonsingular.  It is called by DLSODI for
!! initialization only, when ISTATE = 0 .
!! DAINVG returns an error flag IER:
!!
!!     IER  =  0  means DAINVG was successful.
!!     IER .ge. 2 means RES returned an error flag IRES = IER.
!!     IER .lt. 0 means the a-matrix was found to be singular.
!!
!-----------------------------------------------------------------------
subroutine dainvg(res,adda,Neq,T,Y,Ydot,Miter,Ml,Mu,Pw,Ipvt,Ier)
!
external res, adda
integer , intent(inout) :: Neq
real(kind=dp) :: T
real(kind=dp) :: Y(*)
real(kind=dp) :: Ydot(*)
integer,intent(inout) :: Miter
integer,intent(inout) :: Ml
integer,intent(inout) :: Mu
real(kind=dp) :: Pw(*)
integer :: Ipvt(*)
integer , intent(inout) :: Ier

integer :: i , lenpw , mlp1 , nrowpw
!
   if ( Miter>=4 ) then
   !
   !  Band matrix case -----------------------------------------------------
   !
      nrowpw = 2*Ml + Mu + 1
      lenpw = Neq*nrowpw
      do i = 1 , lenpw
         Pw(i) = 0.0D0
      enddo
   !
      Ier = 1
      call res(Neq,T,Y,Pw,Ydot,Ier)
      if ( Ier>1 ) return
   !
      mlp1 = Ml + 1
      call adda(Neq,T,Y,Ml,Mu,Pw(mlp1),nrowpw)
      call dgbfa(Pw,nrowpw,Neq,Ml,Mu,Ipvt,Ier)
      if ( Ier==0 ) then
         call dgbsl(Pw,nrowpw,Neq,Ml,Mu,Ipvt,Ydot,0)
         return
      endif
   else
   !
   !  Full matrix case -----------------------------------------------------
   !
      lenpw = Neq*Neq
      do i = 1 , lenpw
         Pw(i) = 0.0D0
      enddo
   !
      Ier = 1
      call res(Neq,T,Y,Pw,Ydot,Ier)
      if ( Ier>1 ) return
   !
      call adda(Neq,T,Y,0,0,Pw,Neq)
      call dgefa(Pw,Neq,Neq,Ipvt,Ier)
      if ( Ier==0 ) then
         call dgesl(Pw,Neq,Neq,Ipvt,Ydot,0)
         return
      else
         Ier = -Ier
         return
      endif
   endif
   Ier = -Ier
end subroutine dainvg
