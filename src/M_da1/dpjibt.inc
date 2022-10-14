!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DPJIBT is called by DSTODI to compute and process the matrix
!! P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
!! and r = g(t,y) - A(t,y)*s.
!!
!! Here J is computed by the user-supplied
!! routine JAC if MITER = 1, or by finite differencing if MITER = 2.
!!
!! J is stored in WM, rescaled, and ADDA is called to generate P.
!!
!! P is then subjected to LU decomposition by DDECBT in preparation
!! for later solution of linear systems with P as coefficient matrix.
!!
!! In addition to variables described previously, communication
!! with DPJIBT uses the following:
!!
!! Y
!!
!! : array containing predicted values on entry.
!!
!! RTEM
!!
!! : work array of length N (ACOR in DSTODI).
!!
!! SAVR
!!
!! : array used for output only.  On output it contains the
!! residual evaluated at current values of t and y.
!!
!! S
!!
!! : array containing predicted values of dy/dt (SAVF in DSTODI).
!!
!! WM
!!
!! : real work space for matrices.  On output it contains the
!! LU decomposition of P.
!! Storage of matrix elements starts at WM(3).
!! WM also contains the following matrix-related data:
!! WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
!!
!! IWM
!!
!! : integer work space containing pivot information, starting at
!! IWM(21).  IWM also contains block structure parameters
!! MB = IWM(1) and NB = IWM(2).
!! EL0
!!
!! : EL(1) (input).
!!
!! IERPJ
!!
!! : output error flag.
!!         = 0 if no trouble occurred,
!!         = 1 if the P matrix was found to be unfactorable,
!!         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
!!
!! JCUR
!!
!! : output flag = 1 to indicate that the Jacobian matrix
!! (or approximation) is now current.
!!
!! This routine also uses the Common variables EL0, H, TN, UROUND,
!! MITER, N, NFE, and NJE.
!-----------------------------------------------------------------------
subroutine dpjibt(Neq,Y,Yh,Nyh,Ewt,Rtem,Savr,S,Wm,Iwm,res,jac,adda)
!
integer, dimension(*) :: Neq
real(kind=dp),intent(inout), dimension(*) :: Y
integer,intent(in) :: Nyh
real(kind=dp),intent(in), dimension(Nyh,*) :: Yh
real(kind=dp),intent(in), dimension(*) :: Ewt
real(kind=dp),intent(inout), dimension(*) :: Rtem
real(kind=dp),dimension(*) :: Savr
real(kind=dp),dimension(*) :: S
real(kind=dp),intent(inout), dimension(*) :: Wm
integer, dimension(*) :: Iwm
external :: res
external :: jac
external :: adda
!
real(kind=dp) :: con, fac, hl0, r, srur
integer :: i, ier, iia, iib, iic, ipa, ipb, ipc, ires, j, j1, j2, k, k1, lblox, lenp, lpb, lpc, mb, mbsq, mwid, nb
!
dls1%nje = dls1%nje + 1
hl0 = dls1%h*dls1%el0
dls1%ierpj = 0
dls1%jcur = 1
mb = Iwm(1)
nb = Iwm(2)
mbsq = mb*mb
lblox = mbsq*nb
lpb = 3 + lblox
lpc = lpb + lblox
lenp = 3*lblox
if ( dls1%miter==2 ) then
!
!  If MITER = 2, make 3*MB + 1 calls to RES to approximate J. -----------
   ires = -1
   call res(Neq,dls1%tn,Y,S,Savr,ires)
   dls1%nfe = dls1%nfe + 1
   if ( ires>1 ) then
!  Error return for IRES = 2 or IRES = 3 return from RES. ---------------
      dls1%ierpj = ires
      return
   else
      mwid = 3*mb
      srur = Wm(1)
      do i = 1, lenp
         Wm(2+i) = 0.0D0
      enddo
      do k = 1, 3
         do j = 1, mb
!          Increment Y(I) for group of column indices, and call RES. ----
            j1 = j + (k-1)*mb
            do i = j1, dls1%n, mwid
               r = max(srur*abs(Y(i)),0.01D0/Ewt(i))
               Y(i) = Y(i) + r
            enddo
            call res(Neq,dls1%tn,Y,S,Rtem,ires)
            dls1%nfe = dls1%nfe + 1
            if ( ires>1 ) then
               dls1%ierpj = ires
               return
            else
               do i = 1, dls1%n
                  Rtem(i) = Rtem(i) - Savr(i)
               enddo
               k1 = k
               do i = j1, dls1%n, mwid
!            Get Jacobian elements in column I (block-column K1). -------
                  Y(i) = Yh(i,1)
                  r = max(srur*abs(Y(i)),0.01D0/Ewt(i))
                  fac = -hl0/r
!            Compute and load elements PA(*,J,K1). ----------------------
                  iia = i - j
                  ipa = 2 + (j-1)*mb + (k1-1)*mbsq
                  do j2 = 1, mb
                     Wm(ipa+j2) = Rtem(iia+j2)*fac
                  enddo
                  if ( k1>1 ) then
!            Compute and load elements PB(*,J,K1-1). --------------------
                     iib = iia - mb
                     ipb = ipa + lblox - mbsq
                     do j2 = 1, mb
                        Wm(ipb+j2) = Rtem(iib+j2)*fac
                     enddo
                  endif
                  if ( k1<nb ) then
!            Compute and load elements PC(*,J,K1+1). --------------------
                     iic = iia + mb
                     ipc = ipa + 2*lblox + mbsq
                     do j2 = 1, mb
                        Wm(ipc+j2) = Rtem(iic+j2)*fac
                     enddo
                  endif
                  if ( k1==3 ) then
!            Compute and load elements PC(*,J,1). -----------------------
                     ipc = ipa - 2*mbsq + 2*lblox
                     do j2 = 1, mb
                        Wm(ipc+j2) = Rtem(j2)*fac
                     enddo
                  endif
                  if ( k1==nb-2 ) then
!            Compute and load elements PB(*,J,NB). ----------------------
                     iib = dls1%n - mb
                     ipb = ipa + 2*mbsq + lblox
                     do j2 = 1, mb
                        Wm(ipb+j2) = Rtem(iib+j2)*fac
                     enddo
                  endif
                  k1 = k1 + 3
               enddo
            endif
         enddo
      enddo
!  RES call for first corrector iteration. ------------------------------
      ires = 1
      call res(Neq,dls1%tn,Y,S,Savr,ires)
      dls1%nfe = dls1%nfe + 1
      if ( ires>1 ) then
         dls1%ierpj = ires
         return
      endif
   endif
else
!  If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
   ires = 1
   call res(Neq,dls1%tn,Y,S,Savr,ires)
   dls1%nfe = dls1%nfe + 1
   if ( ires>1 ) then
      dls1%ierpj = ires
      return
   else
      do i = 1, lenp
         Wm(i+2) = 0.0D0
      enddo
      call jac(Neq,dls1%tn,Y,S,mb,nb,Wm(3),Wm(lpb),Wm(lpc))
      con = -hl0
      do i = 1, lenp
         Wm(i+2) = Wm(i+2)*con
      enddo
   endif
endif
!  Add matrix A. --------------------------------------------------------
call adda(Neq,dls1%tn,Y,mb,nb,Wm(3),Wm(lpb),Wm(lpc))
!  Do LU decomposition on P. --------------------------------------------
call ddecbt(mb,nb,Wm(3),Wm(lpb),Wm(lpc),Iwm(21),ier)
if ( ier/=0 ) dls1%ierpj = 1

end subroutine dpjibt
