!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine checks for the presence of a root in the vicinity of
!! the current T, in a manner depending on the input flag JOB.  It calls
!! Subroutine DROOTS to locate the root as precisely as possible.
!!
!! In addition to variables described previously, DRCHEK
!! uses the following for communication:
!!
!! JOB
!!
!! : integer flag indicating type of call:
!!
!!  JOB | Description
!!  --- | ----------------------------------------------------------
!!   1  | means the problem is being initialized, and DRCHEK
!!      | is to look for a root at or very near the initial T.
!!      !
!!   2  | means a continuation call to the solver was just
!!      | made, and DRCHEK is to check for a root in the
!!      | relevant part of the step dlsr taken.
!!      !
!!   3  | means a successful step was just taken, and DRCHEK
!!      | is to look for a root in the interval of the step.
!!
!! G0
!!
!! : array of length NG, containing the value of g at T = T0.
!! G0 is input for JOB .ge. 2, and output in all cases.
!!
!! G1,GX
!!
!! : arrays of length NG for work space.
!!
!! IRT
!!
!! : completion flag:
!!
!!  IRT | Description
!!  --- | ----------------------------------------------------------
!!   0  | means no root was found.
!!  -1  | means JOB = 1 and a root was found too near to T.
!!   1  | means a legitimate root was found (JOB = 2 or 3).
!!
!! On return, T0 is the root location, and Y is the
!! corresponding solution vector.
!!
!! T0
!!
!! : value of T at one endpoint of interval of interest.  Only
!! roots beyond T0 in the direction of integration are sought.
!!
!! T0 is input if JOB .ge. 2, and output in all cases.
!!
!! T0 is updated by DRCHEK, whether a root is found or not.
!!
!! TLAST
!!
!! : dlsr value of T returned by the solver (input only).
!!
!! TOUTC
!!
!! : copy of TOUT (input only).
!!
!! IRFND
!!
!! : input flag showing whether the dlsr step taken had a root.
!! IRFND = 1 if it did, = 0 if not.
!!
!! ITASKC
!!
!! : copy of ITASK (input only).
!!
!! NGC
!!
!! : copy of NG (input only).
!-----------------------------------------------------------------------
subroutine drchek(Job,g_sub,Neq,Y,Yh,Nyh,G0,G1,Gx,Jroot,Irt)
!
integer, intent(in) :: Job
external g_sub
integer             :: Neq(*)
real(kind=dp)       :: Y(*)
integer,intent(in)  :: Nyh
real(kind=dp)       :: Yh(Nyh,*)
real(kind=dp)       :: G0(*)
real(kind=dp)       :: G1(*)
real(kind=dp)       :: Gx(*)
integer             :: Jroot(*)
integer,intent(out) :: Irt

real(kind=dp)       :: hming, t1, temp1, temp2, x
integer             :: i, iflag, jflag
logical             :: zroot

   Irt = 0
   do i = 1, dlsr%ngc
      Jroot(i) = 0
   enddo
   hming = (abs(dls1%tn)+abs(dls1%h))*dls1%uround*100.0D0

   select case (Job)
   case (2)

      if ( dlsr%irfnd/=0 ) then
         !  If a root was found on the previous step, evaluate G0 = g(T0). -------
         call dintdy(dlsr%t0,0,Yh,Nyh,Y,iflag)
         call g_sub(Neq,dlsr%t0,Y,dlsr%ngc,G0)
         dlsr%nge = dlsr%nge + 1
         zroot = .false.
         do i = 1, dlsr%ngc
            if ( abs(G0(i))<=0.0D0 ) zroot = .true.
         enddo
         if ( zroot ) then
            !  g has a zero at T0.  Look at g at T + (small increment). -------------
            temp1 = sign(hming,dls1%h)
            dlsr%t0 = dlsr%t0 + temp1
            if ( (dlsr%t0-dls1%tn)*dls1%h<0.0D0 ) then
               call dintdy(dlsr%t0,0,Yh,Nyh,Y,iflag)
            else
               temp2 = temp1/dls1%h
               do i = 1, dls1%n
                  Y(i) = Y(i) + temp2*Yh(i,2)
               enddo
            endif
            call g_sub(Neq,dlsr%t0,Y,dlsr%ngc,G0)
            dlsr%nge = dlsr%nge + 1
            zroot = .false.
            do i = 1, dlsr%ngc
               if ( abs(G0(i))<=0.0D0 ) then
                  Jroot(i) = 1
                  zroot = .true.
               endif
            enddo
            if ( zroot ) then
            !  g has a zero at T0 and also close to T0.  Return root. ---------------
               Irt = 1
               return
            endif
         endif
      endif
      !  G0 has no zero components.  Proceed to check relevant interval. ------
      if ( dls1%tn==dlsr%tlast ) then
         return
      endif
   case (3)
   case default
      !
      !  Evaluate g at initial T, and check for zero values. ------------------
      dlsr%t0 = dls1%tn
      call g_sub(Neq,dlsr%t0,Y,dlsr%ngc,G0)
      dlsr%nge = 1
      zroot = .false.
      do i = 1, dlsr%ngc
         if ( abs(G0(i))<=0.0D0 ) zroot = .true.
      enddo
      if ( zroot ) then
         !  g has a zero at T.  Look at g at T + (small increment). --------------
         temp2 = max(hming/abs(dls1%h),0.1D0)
         temp1 = temp2*dls1%h
         dlsr%t0 = dlsr%t0 + temp1
         do i = 1, dls1%n
            Y(i) = Y(i) + temp2*Yh(i,2)
         enddo
         call g_sub(Neq,dlsr%t0,Y,dlsr%ngc,G0)
         dlsr%nge = dlsr%nge + 1
         zroot = .false.
         do i = 1, dlsr%ngc
            if ( abs(G0(i))<=0.0D0 ) zroot = .true.
         enddo
         if ( zroot ) then
         !  g has a zero at T and also close to T.  Take error return. -----------
            Irt = -1
            return
         endif
      endif

      return
   endselect
   !
   !  Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
   if ( dlsr%itaskc/=2 .and. dlsr%itaskc/=3 .and. dlsr%itaskc/=5 ) then
      if ( (dlsr%toutc-dls1%tn)*dls1%h<0.0D0 ) then
         t1 = dlsr%toutc
         if ( (t1-dlsr%t0)*dls1%h<=0.0D0 ) then
            return
         endif
         call dintdy(t1,0,Yh,Nyh,Y,iflag)
      endif
      t1 = dls1%tn
      do i = 1, dls1%n
         Y(i) = Yh(i,1)
      enddo
   endif
   call g_sub(Neq,t1,Y,dlsr%ngc,G1)
   dlsr%nge = dlsr%nge + 1
   !  Call DROOTS to search for root in interval from T0 to T1. ------------
   jflag = 0
   do
      call droots(dlsr%ngc,hming,jflag,dlsr%t0,t1,G0,G1,Gx,x,Jroot)
      if ( jflag>1 ) then
         dlsr%t0 = x
         !X!call dcopy(dlsr%ngc,Gx,1,G0,1)
         G0(1:dlsr%ngc) = Gx(1:dlsr%ngc)!X!
         if ( jflag/=4 ) exit
         return
      else
         call dintdy(x,0,Yh,Nyh,Y,iflag)
         call g_sub(Neq,x,Y,dlsr%ngc,Gx)
         dlsr%nge = dlsr%nge + 1
      endif
   enddo
   !  Found a root.  Interpolate to X and return. --------------------------
   call dintdy(x,0,Yh,Nyh,Y,iflag)
   Irt = 1

end subroutine drchek
