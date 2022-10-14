!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DPREPJI is called by DSTODI to compute and process the matrix
!! P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
!! where r = g(t,y) - A(t,y)*s.
!!
!! Here J is computed by the user-supplied
!! routine JAC if MITER = 1 or 4, or by finite differencing if MITER =
!! 2 or 5.
!!
!! J is stored in WM, rescaled, and ADDA is called to generate
!! P.
!!
!! P is then subjected to LU decomposition in preparation
!! for later solution of linear systems with P as coefficient
!! matrix.  This is done by DGEFA if MITER = 1 or 2, and by
!! DGBFA if MITER = 4 or 5.
!!
!! In addition to variables described previously, communication
!! with DPREPJI uses the following:
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
!! IWM(21).  IWM also contains the band parameters
!! ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!!
!! EL0
!!
!! : el(1) (input).
!!
!! IERPJ
!!
!! : output error flag.

!!         = 0 if no trouble occurred,
!!         = 1 if the P matrix was found to be singular,
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
subroutine dprepji(Neq,Y,Yh,Nyh,Ewt,Rtem,Savr,S,Wm,Iwm,res,jac,adda)

!
integer  :: Neq(*)
real(kind=dp),intent(inout) :: Y(*)
integer,intent(in)          :: Nyh
real(kind=dp),intent(in)    :: Yh(Nyh,*)
real(kind=dp),intent(in)    :: Ewt(*)
real(kind=dp)               :: Rtem(*)
real(kind=dp)               :: Savr(*)
real(kind=dp)               :: S(*)
real(kind=dp),intent(inout) :: Wm(*)
integer                     :: Iwm(*)

external adda
external jac
external res
!
real(kind=dp) :: con, fac, hl0, r, srur, yi, yj, yjj
integer :: i, i1, i2, ier, ii, ires, j, j1, jj, lenp, mba, mband, meb1, meband, ml, ml3, mu
!
   dls1%nje = dls1%nje + 1
   hl0 = dls1%h*dls1%el0
   dls1%ierpj = 0
   dls1%jcur = 1
   select case (dls1%miter)
   case (2)
   !  If MITER = 2, make N + 1 calls to RES to approximate J. --------------
      ires = -1
      call res(Neq,dls1%tn,Y,S,Savr,ires)
      dls1%nfe = dls1%nfe + 1
      if ( ires>1 ) then
   !  Error return for IRES = 2 or IRES = 3 return from RES. ---------------
         dls1%ierpj = ires
         return
      else
         srur = Wm(1)
         j1 = 2
         do j = 1, dls1%n
            yj = Y(j)
            r = max(srur*abs(yj),0.01D0/Ewt(j))
            Y(j) = Y(j) + r
            fac = -hl0/r
            call res(Neq,dls1%tn,Y,S,Rtem,ires)
            dls1%nfe = dls1%nfe + 1
            if ( ires>1 ) then
               dls1%ierpj = ires
               return
            else
               do i = 1, dls1%n
                  Wm(i+j1) = (Rtem(i)-Savr(i))*fac
               enddo
               Y(j) = yj
               j1 = j1 + dls1%n
            endif
         enddo
         ires = 1
         call res(Neq,dls1%tn,Y,S,Savr,ires)
         dls1%nfe = dls1%nfe + 1
         if ( ires>1 ) then
            dls1%ierpj = ires
            return
         endif
      endif
   case (3)
   !  Dummy section for MITER = 3
      return
   case (4)
   !  If MITER = 4, call RES, then JAC, and multiply by scalar. ------------
      ires = 1
      call res(Neq,dls1%tn,Y,S,Savr,ires)
      dls1%nfe = dls1%nfe + 1
      if ( ires>1 ) then
         dls1%ierpj = ires
         return
      else
         ml = Iwm(1)
         mu = Iwm(2)
         ml3 = ml + 3
         mband = ml + mu + 1
         meband = mband + ml
         lenp = meband*dls1%n
         do i = 1, lenp
            Wm(i+2) = 0.0D0
         enddo
         call jac(Neq,dls1%tn,Y,S,ml,mu,Wm(ml3),meband)
         con = -hl0
         do i = 1, lenp
            Wm(i+2) = Wm(i+2)*con
         enddo
         call wrapup()
         return
      endif
   case (5)
   !  If MITER = 5, make ML + MU + 2 calls to RES to approximate J. --------
      ires = -1
      call res(Neq,dls1%tn,Y,S,Savr,ires)
      dls1%nfe = dls1%nfe + 1
      if ( ires>1 ) then
         dls1%ierpj = ires
         return
      else
         ml = Iwm(1)
         mu = Iwm(2)
         ml3 = ml + 3
         mband = ml + mu + 1
         mba = min(mband,dls1%n)
         meband = mband + ml
         meb1 = meband - 1
         srur = Wm(1)
         do j = 1, mba
            do i = j, dls1%n, mband
               yi = Y(i)
               r = max(srur*abs(yi),0.01D0/Ewt(i))
               Y(i) = Y(i) + r
            enddo
            call res(Neq,dls1%tn,Y,S,Rtem,ires)
            dls1%nfe = dls1%nfe + 1
            if ( ires>1 ) then
               dls1%ierpj = ires
               return
            else
               do jj = j, dls1%n, mband
                  Y(jj) = Yh(jj,1)
                  yjj = Y(jj)
                  r = max(srur*abs(yjj),0.01D0/Ewt(jj))
                  fac = -hl0/r
                  i1 = max(jj-mu,1)
                  i2 = min(jj+ml,dls1%n)
                  ii = jj*meb1 - ml + 2
                  do i = i1, i2
                     Wm(ii+i) = (Rtem(i)-Savr(i))*fac
                  enddo
               enddo
            endif
         enddo
         ires = 1
         call res(Neq,dls1%tn,Y,S,Savr,ires)
         dls1%nfe = dls1%nfe + 1
         if ( ires<=1 ) then
             call wrapup()
             return
         endif
         dls1%ierpj = ires
         return
      endif
   case default
   !  If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
      ires = 1
      call res(Neq,dls1%tn,Y,S,Savr,ires)
      dls1%nfe = dls1%nfe + 1

      if ( ires>1 ) then
         dls1%ierpj = ires
         return
      else
         lenp = dls1%n*dls1%n
         do i = 1, lenp
            Wm(i+2) = 0.0D0
         enddo
         call jac(Neq,dls1%tn,Y,S,0,0,Wm(3),dls1%n)
         con = -hl0
         do i = 1, lenp
            Wm(i+2) = Wm(i+2)*con
         enddo
      endif

   endselect
   !  Add matrix A. --------------------------------------------------------
   call adda(Neq,dls1%tn,Y,0,0,Wm(3),dls1%n)
   !  Do LU decomposition on P. --------------------------------------------
   call dgefa(Wm(3),dls1%n,dls1%n,Iwm(21),ier)
   if ( ier/=0 ) dls1%ierpj = 1

contains

subroutine wrapup()
   !  Add matrix A. --------------------------------------------------------
   call adda(Neq,dls1%tn,Y,ml,mu,Wm(ml3),meband)
   !  Do LU decomposition of P. --------------------------------------------
   call dgbfa(Wm(3),meband,dls1%n,ml,mu,Iwm(21),ier)
   if ( ier/=0 ) dls1%ierpj = 1
end subroutine wrapup

end subroutine dprepji
