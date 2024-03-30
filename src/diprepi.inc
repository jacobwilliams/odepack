!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine serves as an interface between the driver and
!! Subroutine DPREPI.  Tasks performed here are:
!!
!!  * call DPREPI,
!!  * reset the required WM segment length LENWK,
!!  * move YH back to its final location (following WM in RWORK),
!!  * reset pointers for YH, SAVR, EWT, and ACOR, and
!!  * move EWT to its new position if ISTATE = 0 or 1.
!!
!! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
!! no trouble, and IPFLAG is the value of the DPREPI error flag IPPER
!! if there was trouble in Subroutine DPREPI.
!-----------------------------------------------------------------------
subroutine diprepi(Neq,Y,S,Rwork,Ia,Ja,Ic,Jc,Ipflag,res,jac,adda)
!
integer                        :: Neq(*)
real(kind=dp)                  :: Y(*)
real(kind=dp)                  :: S(*)
real(kind=dp), intent(inout)  :: Rwork(*)
integer                        :: Ia(*)
integer                        :: Ja(*)
integer                        :: Ic(*)
integer                        :: Jc(*)
integer,intent(inout)          :: Ipflag
external :: res
external :: jac
external :: adda

integer :: i, imax, lewtn, lyhd, lyhn

   Ipflag = 0

   !  Call DPREPI to do matrix preprocessing operations. -------------------
   call dprepi(Neq,Y,S,Rwork(dls1%lyh),Rwork(dls1%lsavf),Rwork(dls1%lewt), &
   & Rwork(dls1%lacor),Ia,Ja,Ic,Jc,Rwork(dls1%lwm),Rwork(dls1%lwm), &
   & Ipflag,res,jac,adda)

   dlss%lenwk = max(dlss%lreq,dlss%lwmin)
   if ( Ipflag < 0 ) return

   !  If DPREPI was successful, move YH to end of required space for WM. ---
   lyhn = dls1%lwm + dlss%lenwk
   if ( lyhn > dls1%lyh ) return
   lyhd = dls1%lyh - lyhn

   if ( lyhd /= 0 ) then
      imax = lyhn - 1 + dlss%lenyhm
      do i = lyhn, imax
         Rwork(i) = Rwork(i+lyhd)
      enddo
      dls1%lyh = lyhn
   endif

   !  Reset pointers for SAVR, EWT, and ACOR. ------------------------------
   dls1%lsavf = dls1%lyh + dlss%lenyh
   lewtn = dls1%lsavf + dls1%n
   dls1%lacor = lewtn + dls1%n

   if ( dlss%istatc /= 3 ) then
   !  If ISTATE = 1, move EWT (left) to its new position. ------------------
      if ( lewtn>dls1%lewt ) return
      do i = 1, dls1%n
         Rwork(i+lewtn-1) = Rwork(i+dls1%lewt-1)
      enddo
   endif

   dls1%lewt = lewtn

end subroutine diprepi
