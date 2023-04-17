!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!    dprepj(3f) - [M_odepack] Compute and process Newton iteration matrix.
!!
!!### DESCRIPTION
!!  DPREPJ is called by DSTODE to compute and process the matrix
!!  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
!!
!!  Here J is computed by the user-supplied routine JAC if
!!  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
!!
!!  If MITER = 3, a diagonal approximation to J is used.
!!
!!  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
!!  subjected to LU decomposition in preparation for later solution
!!  of linear systems with P as coefficient matrix.  This is done
!!  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!!
!!  In addition to variables described in DSTODE and DLSODE prologues,
!!  communication with DPREPJ uses the following:
!!
!!  Y
!!
!!  : array containing predicted values on entry.
!!
!!  FTEM
!!
!!  : work array of length N (ACOR in DSTODE).
!!
!!  SAVF
!!
!!  : array containing f evaluated at predicted y.
!!
!!  WM
!!
!!  : real work space for matrices.  On output it contains the
!!    inverse diagonal matrix if MITER = 3 and the LU decomposition
!!    of P if MITER is 1, 2 , 4, or 5.
!!    Storage of matrix elements starts at WM(3).
!!    WM also contains the following matrix-related data:
!!    WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
!!    WM(2) = H*EL0, saved for later use if MITER = 3.
!!
!!  IWM
!!
!!  : integer work space containing pivot information, starting at
!!    IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
!!    parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!!
!!  EL0
!!
!!  : EL(1) (input).
!!
!!  IERPJ
!!
!!  : output error flag,  = 0 if no trouble, .gt. 0 if
!!    P matrix found to be singular.
!!
!!  JCUR
!!
!!  : output flag = 1 to indicate that the Jacobian matrix
!!    (or approximation) is now current.
!!
!!  This routine also uses the COMMON variables EL0, H, TN, UROUND,
!!  MITER, N, NFE, and NJE.
!!
! ### BEGIN PROLOGUE  DPREPJ
! ### SUBSIDIARY
! ### TYPE      DOUBLE PRECISION (SPREPJ-S, DPREPJ-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  DGBFA, DGEFA, DVNORM
! ### COMMON BLOCKS    DLS001
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890504  Minor cosmetic changes.  (FNF)
!     19930809  Renamed to allow single/double precision versions. (ACH)
!     19010418  Reduced size of Common block /DLS001/. (ACH)
!     19031105  Restored 'own' variables to Common block /DLS001/, to
!               enable interrupt/restart feature. (ACH)
! ### END PROLOGUE  DPREPJ
! **End
!-----------------------------------------------------------------------
subroutine dprepj(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wm,Iwm,f,jac)
!
integer                     :: Neq(*)
real(kind=dp),intent(inout) :: Y(*)
integer, intent(in)         :: Nyh
real(kind=dp),intent(in)    :: Yh(Nyh,*)
real(kind=dp)               :: Ewt(*)
real(kind=dp)               :: Ftem(*)
real(kind=dp)               :: Savf(*)
real(kind=dp),intent(inout) :: Wm(*)
integer                     :: Iwm(*)
external                    :: f
external                    :: jac
!
real(kind=dp) :: con, di, fac, hl0, r, r0, srur, yi, yj, yjj
integer :: i, i1, i2, ier, ii, j, j1, jj, lenp, mba, mband, meb1, meband, ml, ml3, mu, np1
!
   dls1%nje = dls1%nje + 1
   dls1%ierpj = 0
   dls1%jcur = 1
   hl0 = dls1%h*dls1%el0
   select case (dls1%miter)
   case (2)
   !  If MITER = 2, make N calls to F to approximate J. --------------------
      fac = dvnorm(dls1%n,Savf,Ewt)
      r0 = 1000.0D0*abs(dls1%h)*dls1%uround*dls1%n*fac
      if ( r0==0.0D0 ) r0 = 1.0D0
      srur = Wm(1)
      j1 = 2
      do j = 1 , dls1%n
         yj = Y(j)
         r = max(srur*abs(yj),r0/Ewt(j))
         Y(j) = Y(j) + r
         fac = -hl0/r
         call f(Neq,dls1%tn,Y,Ftem)
         do i = 1 , dls1%n
            Wm(i+j1) = (Ftem(i)-Savf(i))*fac
         enddo
         Y(j) = yj
         j1 = j1 + dls1%n
      enddo
      dls1%nfe = dls1%nfe + dls1%n
   case (3)
   !  If MITER = 3, construct a diagonal approximation to J and P. ---------
      Wm(2) = hl0
      r = dls1%el0*0.1D0
      do i = 1 , dls1%n
         Y(i) = Y(i) + r*(dls1%h*Savf(i)-Yh(i,2))
      enddo
      call f(Neq,dls1%tn,Y,Wm(3))
      dls1%nfe = dls1%nfe + 1
      do i = 1 , dls1%n
         r0 = dls1%h*Savf(i) - Yh(i,2)
         di = 0.1D0*r0 - dls1%h*(Wm(i+2)-Savf(i))
         Wm(i+2) = 1.0D0
         if ( abs(r0)>=dls1%uround/Ewt(i) ) then
            if ( abs(di)==0.0D0 ) then
               dls1%ierpj = 1
               return
            else
               Wm(i+2) = 0.1D0*r0/di
            endif
         endif
      enddo
      return
   case (4)
   !  If MITER = 4, call JAC and multiply by scalar. -----------------------
      ml = Iwm(1)
      mu = Iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*dls1%n
      do i = 1 , lenp
         Wm(i+2) = 0.0D0
      enddo
      call jac(Neq,dls1%tn,Y,ml,mu,Wm(ml3),meband)
      con = -hl0
      do i = 1 , lenp
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
      fac = dvnorm(dls1%n,Savf,Ewt)
      r0 = 1000.0D0*abs(dls1%h)*dls1%uround*dls1%n*fac
      if ( r0==0.0D0 ) r0 = 1.0D0
      do j = 1 , mba
         do i = j , dls1%n , mband
            yi = Y(i)
            r = max(srur*abs(yi),r0/Ewt(i))
            Y(i) = Y(i) + r
         enddo
         call f(Neq,dls1%tn,Y,Ftem)
         do jj = j , dls1%n , mband
            Y(jj) = Yh(jj,1)
            yjj = Y(jj)
            r = max(srur*abs(yjj),r0/Ewt(jj))
            fac = -hl0/r
            i1 = max(jj-mu,1)
            i2 = min(jj+ml,dls1%n)
            ii = jj*meb1 - ml + 2
            do i = i1 , i2
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
      do i = 1 , lenp
         Wm(i+2) = 0.0D0
      enddo
      call jac(Neq,dls1%tn,Y,0,0,Wm(3),dls1%n)
      con = -hl0
      do i = 1 , lenp
         Wm(i+2) = Wm(i+2)*con
      enddo
   endselect
   !  Add identity matrix. -------------------------------------------------
   j = 3
   np1 = dls1%n + 1
   do i = 1 , dls1%n
      Wm(j) = Wm(j) + 1.0D0
      j = j + np1
   enddo
   !  Do LU decomposition on P. --------------------------------------------
   call dgefa(Wm(3),dls1%n,dls1%n,Iwm(21),ier)
   if ( ier/=0 ) dls1%ierpj = 1

contains

subroutine wrapup()
   !  Add identity matrix. -------------------------------------------------
   ii = mband + 2
   do i = 1 , dls1%n
      Wm(ii) = Wm(ii) + 1.0D0
      ii = ii + meband
   enddo
   !  Do LU decomposition of P. --------------------------------------------
   call dgbfa(Wm(3),meband,dls1%n,ml,mu,Iwm(21),ier)
   if ( ier/=0 ) dls1%ierpj = 1
end subroutine wrapup

end subroutine dprepj
