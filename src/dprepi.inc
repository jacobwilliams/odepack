!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine performs preprocessing related to the sparse linear
!! systems that must be solved.
!!
!! The operations that are performed here are:
!!
!!  * compute sparseness structure of the iteration matrix
!!      P = A - con*J  according to MOSS,
!!  * compute grouping of column indices (MITER = 2),
!!  * compute a new ordering of rows and columns of the matrix,
!!  * reorder JA corresponding to the new ordering,
!!  * perform a symbolic LU factorization of the matrix, and
!!  * set pointers for segments of the IWK/WK array.
!!
!! In addition to variables described previously, DPREPI uses the
!! following for communication:
!!
!! YH
!!
!! : the history array.  Only the first column, containing the
!! current Y vector, is used.  Used only if MOSS .ne. 0.
!!
!! S
!!
!! : array of length NEQ, identical to YDOTI in the driver, used
!! only if MOSS .ne. 0.
!!
!! SAVR
!!
!! : a work array of length NEQ, used only if MOSS .ne. 0.
!!
!! EWT
!!
!! : array of length NEQ containing (inverted) error weights.
!! Used only if MOSS = 2 or 4 or if ISTATE = MOSS = 1.
!!
!! RTEM
!!
!! : a work array of length NEQ, identical to ACOR in the driver,
!! used only if MOSS = 2 or 4.
!!
!! WK
!!
!! : a real work array of length LENWK, identical to WM in
!! the driver.
!!
!! IWK
!!
!! : integer work array, assumed to occupy the same space as WK.
!!
!! LENWK
!!
!! : the length of the work arrays WK and IWK.
!!
!! ISTATC
!!
!! : a copy of the driver input argument ISTATE (= 1 on the
!! first call, = 3 on a continuation call).
!!
!! IYS
!!
!! : flag value from ODRV or CDRV.
!!
!! IPPER
!!
!! : output error flag, with the following values and meanings:
!!        =   0  no error.
!!        =  -1  insufficient storage for internal structure pointers.
!!        =  -2  insufficient storage for JGROUP.
!!        =  -3  insufficient storage for ODRV.
!!        =  -4  other error flag from ODRV (should never occur).
!!        =  -5  insufficient storage for CDRV.
!!        =  -6  other error flag from CDRV.
!!        =  -7  if the RES routine returned error flag IRES = IER = 2.
!!        =  -8  if the RES routine returned error flag IRES = IER = 3.
!-----------------------------------------------------------------------
subroutine dprepi(Neq,Y,S,Yh,Savr,Ewt,Rtem,Ia,Ja,Ic,Jc,Wk,Iwk,Ipper,res,jac,adda)
Use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer                     :: Neq(*)
real(kind=dp),intent(inout) :: Y(*)
real(kind=dp)               :: S(*)
real(kind=dp),intent(in)    :: Yh(*)
real(kind=dp),intent(inout) :: Savr(*)
real(kind=dp),intent(in)    :: Ewt(*)
real(kind=dp),intent(inout) :: Rtem(*)
integer,intent(in)          :: Ia(*)
integer,intent(in)          :: Ja(*)
integer,intent(in)          :: Ic(*)
integer,intent(in)          :: Jc(*)
real(kind=dp),intent(inout) :: Wk(*)
integer,intent(inout)       :: Iwk(*)
integer,intent(out)         :: Ipper
external res
external jac
external adda

real(kind=dp) :: erwt, fac, yj
integer :: i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, k, kamax, kamin, kcmax, kcmin, knew, ldif, lenigp, lenwk1, &
 & liwk, ljfo, maxg, np1, nzsut

   dlss%ibian = dlss%lrat*2
   dlss%ipian = dlss%ibian + 1
   np1 = dls1%n + 1
   dlss%ipjan = dlss%ipian + np1
   dlss%ibjan = dlss%ipjan - 1
   lenwk1 = dlss%lenwk - dls1%n
   liwk = dlss%lenwk*dlss%lrat
   if ( dlss%moss==0 ) liwk = liwk - dls1%n
   if ( dlss%moss==1 .or. dlss%moss==2 ) liwk = lenwk1*dlss%lrat
   if ( dlss%ipjan+dls1%n-1>liwk ) goto 600
   if ( dlss%moss/=0 ) then
      !
      if ( dlss%istatc/=3 ) then
         !  ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination.
         !  Initialize S with random nonzero elements for structure determination.
         do i = 1, dls1%n
            erwt = 1.0D0/Ewt(i)
            fac = 1.0D0 + 1.0D0/(i+1.0D0)
            Y(i) = Y(i) + fac*sign(erwt,Y(i))
            S(i) = 1.0D0 + fac*erwt
         enddo
         select case (dlss%moss)
         case (1)
            goto 100
         case (2)
            goto 200
         case (3)
            goto 300
         case (4)
            goto 400
         case default
         endselect
      endif
      !
      !  ISTATE = 3 and MOSS .ne. 0. Load Y from YH(*,1) and S from YH(*,2). --
      do i = 1, dls1%n
         Y(i) = Yh(i)
         S(i) = Yh(dls1%n+i)
      enddo
      select case (dlss%moss)
      case (1)
         goto 100
      case (2)
         goto 200
      case (3)
         goto 300
      case (4)
         goto 400
      case default
      endselect
   endif
   !
   !  MOSS = 0. Process user's IA,JA and IC,JC. ----------------------------
   knew = dlss%ipjan
   kamin = Ia(1)
   kcmin = Ic(1)
   Iwk(dlss%ipian) = 1
   do j = 1, dls1%n
      do i = 1, dls1%n
         Iwk(liwk+i) = 0
      enddo
      kamax = Ia(j+1) - 1
      if ( kamin<=kamax ) then
         do k = kamin, kamax
            i = Ja(k)
            Iwk(liwk+i) = 1
            if ( knew>liwk ) goto 600
            Iwk(knew) = i
            knew = knew + 1
         enddo
      endif
      kamin = kamax + 1
      kcmax = Ic(j+1) - 1
      if ( kcmin<=kcmax ) then
         do k = kcmin, kcmax
            i = Jc(k)
            if ( Iwk(liwk+i)==0 ) then
               if ( knew>liwk ) goto 600
               Iwk(knew) = i
               knew = knew + 1
            endif
         enddo
      endif
      Iwk(dlss%ipian+j) = knew + 1 - dlss%ipjan
      kcmin = kcmax + 1
   enddo
   goto 500
   !
   !  MOSS = 1. Compute structure from user-supplied Jacobian routine JAC. -
   !  A dummy call to RES allows user to create temporaries for use in JAC.
    100  continue
   ier = 1
   call res(Neq,dls1%tn,Y,S,Savr,ier)
   if ( ier>1 ) goto 1000
   do i = 1, dls1%n
      Savr(i) = 0.0D0
      Wk(lenwk1+i) = 0.0D0
   enddo
   k = dlss%ipjan
   Iwk(dlss%ipian) = 1
   do j = 1, dls1%n
      call adda(Neq,dls1%tn,Y,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk(lenwk1+1))
      call jac(Neq,dls1%tn,Y,S,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Savr)
      do i = 1, dls1%n
         ljfo = lenwk1 + i
         if ( Wk(ljfo)==0.0D0 ) then
            if ( Savr(i)==0.0D0 ) cycle
            Savr(i) = 0.0D0
         else
            Wk(ljfo) = 0.0D0
            Savr(i) = 0.0D0
         endif
         if ( k>liwk ) goto 600
         Iwk(k) = i
         k = k + 1
      enddo
      Iwk(dlss%ipian+j) = k + 1 - dlss%ipjan
   enddo
   goto 500
   !
   !  MOSS = 2. Compute structure from results of N + 1 calls to RES. ------
    200  continue
   do i = 1, dls1%n
      Wk(lenwk1+i) = 0.0D0
   enddo
   k = dlss%ipjan
   Iwk(dlss%ipian) = 1
   ier = -1
   if ( dls1%miter==1 ) ier = 1
   call res(Neq,dls1%tn,Y,S,Savr,ier)
   if ( ier>1 ) goto 1000
   do j = 1, dls1%n
      call adda(Neq,dls1%tn,Y,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk(lenwk1+1))
      yj = Y(j)
      erwt = 1.0D0/Ewt(j)
      Y(j) = yj + sign(erwt,yj)
      call res(Neq,dls1%tn,Y,S,Rtem,ier)
      if ( ier>1 ) return
      Y(j) = yj
      do i = 1, dls1%n
         ljfo = lenwk1 + i
         if ( Wk(ljfo)/=0.0D0 ) then
            Wk(ljfo) = 0.0D0
         elseif ( Rtem(i)==Savr(i) ) then
            cycle
         endif
         if ( k>liwk ) goto 600
         Iwk(k) = i
         k = k + 1
      enddo
      Iwk(dlss%ipian+j) = k + 1 - dlss%ipjan
   enddo
   goto 500
   !
   !  MOSS = 3. Compute structure from the user's IA/JA and JAC routine. ---
   !  A dummy call to RES allows user to create temporaries for use in JAC.
    300  continue
   ier = 1
   call res(Neq,dls1%tn,Y,S,Savr,ier)
   if ( ier>1 ) goto 1000
   do i = 1, dls1%n
      Savr(i) = 0.0D0
   enddo
   knew = dlss%ipjan
   kamin = Ia(1)
   Iwk(dlss%ipian) = 1
   do j = 1, dls1%n
      call jac(Neq,dls1%tn,Y,S,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Savr)
      kamax = Ia(j+1) - 1
      if ( kamin<=kamax ) then
         do k = kamin, kamax
            i = Ja(k)
            Savr(i) = 0.0D0
            if ( knew>liwk ) goto 600
            Iwk(knew) = i
            knew = knew + 1
         enddo
      endif
      kamin = kamax + 1
      do i = 1, dls1%n
         if ( Savr(i)/=0.0D0 ) then
            Savr(i) = 0.0D0
            if ( knew>liwk ) goto 600
            Iwk(knew) = i
            knew = knew + 1
         endif
      enddo
      Iwk(dlss%ipian+j) = knew + 1 - dlss%ipjan
   enddo
   goto 500
   !
   !  MOSS = 4. Compute structure from user's IA/JA and N + 1 RES calls. ---
    400  continue
   knew = dlss%ipjan
   kamin = Ia(1)
   Iwk(dlss%ipian) = 1
   ier = -1
   if ( dls1%miter==1 ) ier = 1
   call res(Neq,dls1%tn,Y,S,Savr,ier)
   if ( ier>1 ) goto 1000
   do j = 1, dls1%n
      yj = Y(j)
      erwt = 1.0D0/Ewt(j)
      Y(j) = yj + sign(erwt,yj)
      call res(Neq,dls1%tn,Y,S,Rtem,ier)
      if ( ier>1 ) return
      Y(j) = yj
      kamax = Ia(j+1) - 1
      if ( kamin<=kamax ) then
         do k = kamin, kamax
            i = Ja(k)
            Rtem(i) = Savr(i)
            if ( knew>liwk ) goto 600
            Iwk(knew) = i
            knew = knew + 1
         enddo
      endif
      kamin = kamax + 1
      do i = 1, dls1%n
         if ( Rtem(i)/=Savr(i) ) then
            if ( knew>liwk ) goto 600
            Iwk(knew) = i
            knew = knew + 1
         endif
      enddo
      Iwk(dlss%ipian+j) = knew + 1 - dlss%ipjan
   enddo
   !
    500  continue
   if ( dlss%moss/=0 .and. dlss%istatc/=3 ) then
   !  If ISTATE = 0 or 1 and MOSS .ne. 0, restore Y from YH. ---------------
      do i = 1, dls1%n
         Y(i) = Yh(i)
      enddo
   endif
   dlss%nnz = Iwk(dlss%ipian+dls1%n) - 1
   Ipper = 0
   dlss%ngp = 0
   lenigp = 0
   dlss%ipigp = dlss%ipjan + dlss%nnz
   if ( dls1%miter==2 ) then
   !
   !  Compute grouping of column indices (MITER = 2). ----------------------
   !
      maxg = np1
      dlss%ipjgp = dlss%ipjan + dlss%nnz
      dlss%ibjgp = dlss%ipjgp - 1
      dlss%ipigp = dlss%ipjgp + dls1%n
      iptt1 = dlss%ipigp + np1
      iptt2 = iptt1 + dls1%n
      dlss%lreq = iptt2 + dls1%n - 1
      if ( dlss%lreq>liwk ) goto 700
      call jgroup(dls1%n,Iwk(dlss%ipian),Iwk(dlss%ipjan),maxg,dlss%ngp,Iwk(dlss%ipigp),Iwk(dlss%ipjgp),Iwk(iptt1),Iwk(iptt2),ier)
      if ( ier/=0 ) goto 700
      lenigp = dlss%ngp + 1
   endif
   !
   !  Compute new ordering of rows/columns of Jacobian. --------------------
   dlss%ipr = dlss%ipigp + lenigp
   dlss%ipc = dlss%ipr
   dlss%ipic = dlss%ipc + dls1%n
   dlss%ipisp = dlss%ipic + dls1%n
   dlss%iprsp = (dlss%ipisp-2)/dlss%lrat + 2
   dlss%iesp = dlss%lenwk + 1 - dlss%iprsp
   if ( dlss%iesp<0 ) goto 800
   ibr = dlss%ipr - 1
   do i = 1, dls1%n
      Iwk(ibr+i) = i
   enddo
   dlss%nsp = liwk + 1 - dlss%ipisp
   call odrv(dls1%n,Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk,Iwk(dlss%ipr),Iwk(dlss%ipic),dlss%nsp,Iwk(dlss%ipisp),1,dlss%iys)
   if ( dlss%iys==11*dls1%n+1 ) then
   !
      Ipper = -4
      return
   else
      if ( dlss%iys/=0 ) goto 800
   !
   !  Reorder JAN and do symbolic LU factorization of matrix. --------------
      dlss%ipa = dlss%lenwk + 1 - dlss%nnz
      dlss%nsp = dlss%ipa - dlss%iprsp
      dlss%lreq = max(12*dls1%n/dlss%lrat,6*dls1%n/dlss%lrat+2*dls1%n+dlss%nnz) + 3
      dlss%lreq = dlss%lreq + dlss%iprsp - 1 + dlss%nnz
      if ( dlss%lreq>dlss%lenwk ) goto 900
      dlss%iba = dlss%ipa - 1
      do i = 1, dlss%nnz
         Wk(dlss%iba+i) = 0.0D0
      enddo
      dlss%ipisp = dlss%lrat*(dlss%iprsp-1) + 1
      call cdrv(dls1%n,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic),Iwk(dlss%ipian), &
       & Iwk(dlss%ipjan),Wk(dlss%ipa),Wk(dlss%ipa), &
       & Wk(dlss%ipa),dlss%nsp,Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,5,dlss%iys)
      dlss%lreq = dlss%lenwk - dlss%iesp
      if ( dlss%iys==10*dls1%n+1 ) goto 900
      if ( dlss%iys/=0 ) then
   !
         Ipper = -6
         dlss%lreq = dlss%lenwk
         return
      else
         ipil = dlss%ipisp
         ipiu = ipil + 2*dls1%n + 1
         dlss%nzu = Iwk(ipil+dls1%n) - Iwk(ipil)
         dlss%nzl = Iwk(ipiu+dls1%n) - Iwk(ipiu)
         if ( dlss%lrat<=1 ) then
            call adjlr(dls1%n,Iwk(dlss%ipisp),ldif)
            dlss%lreq = dlss%lreq + ldif
         endif
         if ( dlss%lrat==2 .and. dlss%nnz==dls1%n ) dlss%lreq = dlss%lreq + 1
         dlss%nsp = dlss%nsp + dlss%lreq - dlss%lenwk
         dlss%ipa = dlss%lreq + 1 - dlss%nnz
         dlss%iba = dlss%ipa - 1
         Ipper = 0
         return
      endif
   endif
   !
    600  continue
   Ipper = -1
   dlss%lreq = 2 + (2*dls1%n+1)/dlss%lrat
   dlss%lreq = max(dlss%lenwk+1,dlss%lreq)
   return
   !
    700  continue
   Ipper = -2
   dlss%lreq = (dlss%lreq-1)/dlss%lrat + 1
   return
   !
    800  continue
   Ipper = -3
   call cntnzu(dls1%n,Iwk(dlss%ipian),Iwk(dlss%ipjan),nzsut)
   dlss%lreq = dlss%lenwk - dlss%iesp + (3*dls1%n+4*nzsut-1)/dlss%lrat + 1
   return
   !
    900  continue
   Ipper = -5
   return
   !
    1000 continue
   Ipper = -ier - 5
   dlss%lreq = 2 + (2*dls1%n+1)/dlss%lrat
end subroutine dprepi
