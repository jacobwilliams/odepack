!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine performs preprocessing related to the sparse linear
!! systems that must be solved if MITER = 1 or 2.
!!
!! The operations that are performed here are:
!!
!!  * compute sparseness structure of Jacobian according to MOSS,
!!  * compute grouping of column indices (MITER = 2),
!!  * compute a new ordering of rows and columns of the matrix,
!!  * reorder JA corresponding to the new ordering,
!!  * perform a symbolic LU factorization of the matrix, and
!!  * set pointers for segments of the IWK/WK array.
!!
!! In addition to variables described previously, DPREP uses the
!! following for communication:
!!
!! YH
!!
!! : the history array.  Only the first column, containing the
!! current Y vector, is used.  Used only if MOSS .ne. 0.
!!
!! SAVF
!!
!! : a work array of length NEQ, used only if MOSS .ne. 0.
!!
!! EWT
!!
!! : array of length NEQ containing (inverted) error weights.
!! Used only if MOSS = 2 or if ISTATE = MOSS = 1.
!!
!! FTEM
!!
!! : a work array of length NEQ, identical to ACOR in the driver,
!! used only if MOSS = 2.
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
!! : output error flag with the following values and meanings:
!!          0  no error.
!!         -1  insufficient storage for internal structure pointers.
!!         -2  insufficient storage for JGROUP.
!!         -3  insufficient storage for ODRV.
!!         -4  other error flag from ODRV (should never occur).
!!         -5  insufficient storage for CDRV.
!!         -6  other error flag from CDRV.
!-----------------------------------------------------------------------
subroutine dprep(Neq,Y,Yh,Savf,Ewt,Ftem,Ia,Ja,Wk,Iwk,Ipper,f,jac)
Use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)

integer                      :: Neq(*)
real(kind=dp),intent(inout)  :: Y(*)
real(kind=dp),intent(in)     :: Yh(*)
real(kind=dp),intent(inout)  :: Savf(*)
real(kind=dp),intent(in)     :: Ewt(*)
real(kind=dp)                :: Ftem(*)
integer,intent(in)           :: Ia(*)
integer,intent(in)           :: Ja(*)
real(kind=dp)                :: Wk(*)
integer,intent(inout)        :: Iwk(*)
integer,intent(out)          :: Ipper
external                     :: f
external                     :: jac

real(kind=dp) :: dq, dyj, erwt, fac, yj
integer :: i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k, kmax, kmin, knew, ldif, lenigp, liwk, maxg, np1, nzsut

dlss%ibian = dlss%lrat*2
dlss%ipian = dlss%ibian + 1
np1 = dls1%n + 1
dlss%ipjan = dlss%ipian + np1
dlss%ibjan = dlss%ipjan - 1
liwk = dlss%lenwk*dlss%lrat
if ( dlss%ipjan+dls1%n-1>liwk ) then
   call wrapup400()
   return
endif
if ( dlss%moss/=0 ) then
!
   if ( dlss%istatc/=3 ) then
!  ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
      do i = 1, dls1%n
         erwt = 1.0D0/Ewt(i)
         fac = 1.0D0 + 1.0D0/(i+1.0D0)
         Y(i) = Y(i) + fac*sign(erwt,Y(i))
      enddo
      if ( dlss%moss==1 ) goto 100
      if ( dlss%moss==2 ) goto 200
   endif
!
!  ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
   do i = 1, dls1%n
      Y(i) = Yh(i)
   enddo
   if ( dlss%moss==1 ) goto 100
   if ( dlss%moss==2 ) goto 200
endif
!
!  MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
knew = dlss%ipjan
kmin = Ia(1)
Iwk(dlss%ipian) = 1
do j = 1, dls1%n
   jfound = 0
   kmax = Ia(j+1) - 1
   if ( kmin<=kmax ) then
      do k = kmin, kmax
         i = Ja(k)
         if ( i==j ) jfound = 1
         if ( knew>liwk ) then
            call wrapup400()
            return
         endif
         Iwk(knew) = i
         knew = knew + 1
      enddo
      if ( jfound==1 ) goto 50
   endif
   if ( knew>liwk ) then
      call wrapup400()
      return
   endif
   Iwk(knew) = j
   knew = knew + 1
 50   continue
   Iwk(dlss%ipian+j) = knew + 1 - dlss%ipjan
   kmin = kmax + 1
enddo
goto 300
!
!  MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
!  A dummy call to F allows user to create temporaries for use in JAC. --
 100  continue
call f(Neq,dls1%tn,Y,Savf)
k = dlss%ipjan
Iwk(dlss%ipian) = 1
do j = 1, dls1%n
   if ( k>liwk )then
      call wrapup400()
      return
   endif
   Iwk(k) = j
   k = k + 1
   do i = 1, dls1%n
      Savf(i) = 0.0D0
   enddo
   call jac(Neq,dls1%tn,Y,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Savf)
   do i = 1, dls1%n
      if ( abs(Savf(i))>dlss%seth ) then
         if ( i/=j ) then
            if ( k>liwk ) then
               call wrapup400()
               return
            endif
            Iwk(k) = i
            k = k + 1
         endif
      endif
   enddo
   Iwk(dlss%ipian+j) = k + 1 - dlss%ipjan
enddo
goto 300
!
!  MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
 200  continue
k = dlss%ipjan
Iwk(dlss%ipian) = 1
call f(Neq,dls1%tn,Y,Savf)
do j = 1, dls1%n
   if ( k>liwk ) then
      call wrapup400()
      return
   endif
   Iwk(k) = j
   k = k + 1
   yj = Y(j)
   erwt = 1.0D0/Ewt(j)
   dyj = sign(erwt,yj)
   Y(j) = yj + dyj
   call f(Neq,dls1%tn,Y,Ftem)
   Y(j) = yj
   do i = 1, dls1%n
      dq = (Ftem(i)-Savf(i))/dyj
      if ( abs(dq)>dlss%seth ) then
         if ( i/=j ) then
            if ( k>liwk ) then
               call wrapup400()
               return
            endif
            Iwk(k) = i
            k = k + 1
         endif
      endif
   enddo
   Iwk(dlss%ipian+j) = k + 1 - dlss%ipjan
enddo
!
 300  continue
if ( dlss%moss/=0 .and. dlss%istatc==1 ) then
!  If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
   do i = 1, dls1%n
      Y(i) = Yh(i)
   enddo
endif
dlss%nnz = Iwk(dlss%ipian+dls1%n) - 1
lenigp = 0
dlss%ipigp = dlss%ipjan + dlss%nnz
if ( dls1%miter==2 ) then
!
!  Compute grouping of column indices (MITER = 2). ----------------------
   maxg = np1
   dlss%ipjgp = dlss%ipjan + dlss%nnz
   dlss%ibjgp = dlss%ipjgp - 1
   dlss%ipigp = dlss%ipjgp + dls1%n
   iptt1 = dlss%ipigp + np1
   iptt2 = iptt1 + dls1%n
   dlss%lreq = iptt2 + dls1%n - 1
   if ( dlss%lreq>liwk ) then
      call wrapup500()
      return
   endif
   call jgroup(dls1%n,Iwk(dlss%ipian),Iwk(dlss%ipjan),maxg,dlss%ngp,Iwk(dlss%ipigp),Iwk(dlss%ipjgp),Iwk(iptt1),Iwk(iptt2),ier)
   if ( ier/=0 ) then
      call wrapup500()
      return
   endif
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
if ( dlss%iesp<0 ) then
   call wrapup600()
   return
endif
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
   if ( dlss%iys/=0 ) then
       call wrapup600()
       return
   endif
!
!  Reorder JAN and do symbolic LU factorization of matrix. --------------
   dlss%ipa = dlss%lenwk + 1 - dlss%nnz
   dlss%nsp = dlss%ipa - dlss%iprsp
   dlss%lreq = max(12*dls1%n/dlss%lrat,6*dls1%n/dlss%lrat+2*dls1%n+dlss%nnz) + 3
   dlss%lreq = dlss%lreq + dlss%iprsp - 1 + dlss%nnz
   if ( dlss%lreq>dlss%lenwk ) then
       call wrapup700()
       return
   endif
   dlss%iba = dlss%ipa - 1
   do i = 1, dlss%nnz
      Wk(dlss%iba+i) = 0.0D0
   enddo
   dlss%ipisp = dlss%lrat*(dlss%iprsp-1) + 1

   call cdrv(dls1%n,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic), &
    & Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk(dlss%ipa),Wk(dlss%ipa), &
    & Wk(dlss%ipa),dlss%nsp,Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,5,dlss%iys)

   dlss%lreq = dlss%lenwk - dlss%iesp
   if ( dlss%iys==10*dls1%n+1 ) then
       call wrapup700()
       return
   endif
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

contains

subroutine wrapup400()
   Ipper = -1
   dlss%lreq = 2 + (2*dls1%n+1)/dlss%lrat
   dlss%lreq = max(dlss%lenwk+1,dlss%lreq)
end subroutine wrapup400

subroutine wrapup500()
   Ipper = -2
   dlss%lreq = (dlss%lreq-1)/dlss%lrat + 1
end subroutine wrapup500

subroutine wrapup600()
   Ipper = -3
   call cntnzu(dls1%n,Iwk(dlss%ipian),Iwk(dlss%ipjan),nzsut)
   dlss%lreq = dlss%lenwk - dlss%iesp + (3*dls1%n+4*nzsut-1)/dlss%lrat + 1
end subroutine wrapup600

subroutine wrapup700()
   Ipper = -5
end subroutine wrapup700

end subroutine dprep
