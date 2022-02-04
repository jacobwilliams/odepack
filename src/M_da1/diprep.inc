!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine serves as an interface between the driver and
!! Subroutine DPREP.  It is called only if MITER is 1 or 2.
!! Tasks performed here are:
!!
!!  * call DPREP,
!!  * reset the required WM segment length LENWK,
!!  * move YH back to its final location (following WM in RWORK),
!!  * reset pointers for YH, SAVF, EWT, and ACOR, and
!!  * move EWT to its new position if ISTATE = 1.
!!
!! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
!! no trouble, and IPFLAG is the value of the DPREP error flag IPPER
!! if there was trouble in Subroutine DPREP.
!-----------------------------------------------------------------------
subroutine diprep(Neq,Y,Rwork,Ia,Ja,Ipflag,f,jac)

integer                     :: Neq(*)
real(kind=dp)               :: Y(*)
real(kind=dp),intent(inout) :: Rwork(*)
integer                     :: Ia(*)
integer                     :: Ja(*)
integer,intent(inout)       :: Ipflag
external                    :: f
external                    :: jac

integer :: i, imax, lewtn, lyhd, lyhn

!integer :: Iwork2(size(Rwork)*2) !SCRATCH SPACE

   Ipflag = 0
   ! ===============================================
   ! subroutine dprep (neq, y, yh, savf, ewt, ftem, ia, ja, wk, iwk, ipper, f, jac)
   ! integer          :: neq(*)
   ! double precision :: y(*)
   ! double precision :: yh(*)
   ! double precision :: savf(*)
   ! double precision :: ewt(*)
   ! double precision :: ftem(*)
   ! integer          :: ia(*)
   ! integer          :: ja(*)
   ! double precision :: wk(*)
   ! integer          :: iwk(*)
   ! integer          :: ipper
   ! external         :: f
   ! external         :: jac
   ! ===============================================
   !  Call DPREP to do matrix preprocessing operations. --------------------
   call dprep(Neq,Y, &
   & Rwork(dls1%lyh),Rwork(dls1%lsavf),Rwork(dls1%lewt), Rwork(dls1%lacor), &
   & Ia,Ja, &
   & Rwork(dls1%lwm), &
   & Rwork(dls1%lwm), &
   & Ipflag, &
   & f,jac)

   dlss%lenwk = max(dlss%lreq,dlss%lwmin)
   if ( Ipflag<0 ) return
   !  If DPREP was successful, move YH to end of required space for WM. ----
   lyhn = dls1%lwm + dlss%lenwk
   if ( lyhn>dls1%lyh ) return
   lyhd = dls1%lyh - lyhn
   if ( lyhd/=0 ) then
      imax = lyhn - 1 + dlss%lenyhm
      do i = lyhn , imax
         Rwork(i) = Rwork(i+lyhd)
      enddo
      dls1%lyh = lyhn
   endif
   !  Reset pointers for SAVF, EWT, and ACOR. ------------------------------
   dls1%lsavf = dls1%lyh + dlss%lenyh
   lewtn = dls1%lsavf + dls1%n
   dls1%lacor = lewtn + dls1%n
   if ( dlss%istatc/=3 ) then
   !  If ISTATE = 1, move EWT (left) to its new position. ------------------
      if ( lewtn>dls1%lewt ) return
      do i = 1 , dls1%n
         Rwork(i+lewtn-1) = Rwork(i+dls1%lewt-1)
      enddo
   endif
   dls1%lewt = lewtn
end subroutine diprep
