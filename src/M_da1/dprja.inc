!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DPRJA is called by DSTODA to compute and process the matrix
!! P = I - H*EL(1)*J, where J is an approximation to the Jacobian.
!!
!! Here J is computed by the user-supplied routine JAC if
!! MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
!!
!! J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
!! matrix norm consistent with the weighted max-norm on vectors given
!! by DMNORM) is computed, and J is overwritten by P.
!!
!! P is then
!! subjected to LU decomposition in preparation for later solution
!! of linear systems with P as coefficient matrix.  This is done
!! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!!
!! In addition to variables described previously, communication
!! with DPRJA uses the following:
!!
!! Y
!!
!! : array containing predicted values on entry.
!!
!! FTEM
!!
!! : work array of length N (ACOR in DSTODA).
!!
!! SAVF
!!
!! : array containing f evaluated at predicted y.
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
!! IWM(21).   IWM also contains the band parameters
!! ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!!
!! EL0
!!
!! : EL(1) (input).
!!
!! PDNORM
!!
!! : norm of Jacobian matrix. (Output).
!!
!! IERPJ
!!
!! : output error flag,  = 0 if no trouble, .gt. 0 if
!! P matrix found to be singular.
!!
!! JCUR
!!
!! : output flag = 1 to indicate that the Jacobian matrix
!! (or approximation) is now current.
!!
!! This routine also uses the Common variables EL0, H, TN, UROUND,
!! MITER, N, NFE, and NJE.
!-----------------------------------------------------------------------
subroutine dprja(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wm,Iwm,f,jac)
!
integer, dimension(*) :: Neq
real(kind=dp), intent(inout), dimension(*) :: Y
integer, intent(in) :: Nyh
real(kind=dp), intent(in), dimension(Nyh,*) :: Yh
real(kind=dp), dimension(*) :: Ewt
real(kind=dp), dimension(*) :: Ftem
real(kind=dp), dimension(*) :: Savf
real(kind=dp), intent(inout), dimension(*) :: Wm
integer, dimension(*) :: Iwm
external :: f
external :: jac
!
real(kind=dp) :: con, fac, hl0, r, r0, srur, yi, yj, yjj
integer :: i, i1, i2, ier, ii, j, j1, jj, lenp, mba, mband, meb1, meband, ml, ml3, mu, np1
!
dls1%nje = dls1%nje + 1
dls1%ierpj = 0
dls1%jcur = 1
hl0 = dls1%h*dls1%el0
select case (dls1%miter)
case (2)
!  If MITER = 2, make N calls to F to approximate J. --------------------
   fac = dmnorm(dls1%n,Savf,Ewt)
   r0 = 1000.0D0*abs(dls1%h)*dls1%uround*dls1%n*fac
   if ( r0==0.0D0 ) r0 = 1.0D0
   srur = Wm(1)
   j1 = 2
   do j = 1, dls1%n
      yj = Y(j)
      r = max(srur*abs(yj),r0/Ewt(j))
      Y(j) = Y(j) + r
      fac = -hl0/r
      call f(Neq,dls1%tn,Y,Ftem)
      do i = 1, dls1%n
         Wm(i+j1) = (Ftem(i)-Savf(i))*fac
      enddo
      Y(j) = yj
      j1 = j1 + dls1%n
   enddo
   dls1%nfe = dls1%nfe + dls1%n
case (3)
!  Dummy block only, since MITER is never 3 in this routine. ------------
   return
case (4)
!  If MITER = 4, call JAC and multiply by scalar. -----------------------
   ml = Iwm(1)
   mu = Iwm(2)
   ml3 = ml + 3
   mband = ml + mu + 1
   meband = mband + ml
   lenp = meband*dls1%n
   do i = 1, lenp
      Wm(i+2) = 0.0D0
   enddo
   call jac(Neq,dls1%tn,Y,ml,mu,Wm(ml3),meband)
   con = -hl0
   do i = 1, lenp
      Wm(i+2) = Wm(i+2)*con
   enddo
   call wrapup()
   return
case (5)
!  If MITER = 5, make MBAND calls to F to approximate J. ----------------
   ml = Iwm(1)
   mu = Iwm(2)
   mband = ml + mu + 1
   mba = min(mband,dls1%n)
   meband = mband + ml
   meb1 = meband - 1
   srur = Wm(1)
   fac = dmnorm(dls1%n,Savf,Ewt)
   r0 = 1000.0D0*abs(dls1%h)*dls1%uround*dls1%n*fac
   if ( r0==0.0D0 ) r0 = 1.0D0
   do j = 1, mba
      do i = j, dls1%n, mband
         yi = Y(i)
         r = max(srur*abs(yi),r0/Ewt(i))
         Y(i) = Y(i) + r
      enddo
      call f(Neq,dls1%tn,Y,Ftem)
      do jj = j, dls1%n, mband
         Y(jj) = Yh(jj,1)
         yjj = Y(jj)
         r = max(srur*abs(yjj),r0/Ewt(jj))
         fac = -hl0/r
         i1 = max(jj-mu,1)
         i2 = min(jj+ml,dls1%n)
         ii = jj*meb1 - ml + 2
         do i = i1, i2
            Wm(ii+i) = (Ftem(i)-Savf(i))*fac
         enddo
      enddo
   enddo
   dls1%nfe = dls1%nfe + mba
   call wrapup()
   return
case default
!  If MITER = 1, call JAC and multiply by scalar. -----------------------
   lenp = dls1%n*dls1%n
   do i = 1, lenp
      Wm(i+2) = 0.0D0
   enddo
   call jac(Neq,dls1%tn,Y,0,0,Wm(3),dls1%n)
   con = -hl0
   do i = 1, lenp
      Wm(i+2) = Wm(i+2)*con
   enddo
endselect
!  Compute norm of Jacobian. --------------------------------------------
dlsa%pdnorm = dfnorm(dls1%n,Wm(3),Ewt)/abs(hl0)
!  Add identity matrix. -------------------------------------------------
j = 3
np1 = dls1%n + 1
do i = 1, dls1%n
   Wm(j) = Wm(j) + 1.0D0
   j = j + np1
enddo
!  Do LU decomposition on P. --------------------------------------------
call dgefa(Wm(3),dls1%n,dls1%n,Iwm(21),ier)
if ( ier/=0 ) dls1%ierpj = 1
return

contains

subroutine wrapup()

   !  Compute norm of Jacobian. --------------------------------------------
   dlsa%pdnorm = dbnorm(dls1%n,Wm(ml+3),meband,ml,mu,Ewt)/abs(hl0)
   !  Add identity matrix. -------------------------------------------------
   ii = mband + 2

   do i = 1, dls1%n
      Wm(ii) = Wm(ii) + 1.0D0
      ii = ii + meband
   enddo

   !  Do LU decomposition of P. --------------------------------------------
   call dgbfa(Wm(3),meband,dls1%n,ml,mu,Iwm(21),ier)

   if ( ier/=0 ) dls1%ierpj = 1

end subroutine wrapup

end subroutine dprja
