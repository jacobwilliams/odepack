!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DPRJS is called to compute and process the matrix
!! P = I - H*EL(1)*J, where J is an approximation to the Jacobian.
!! J is computed by columns, either by the user-supplied routine JAC
!! if MITER = 1, or by finite differencing if MITER = 2.
!!
!! Alternatively, if MITER = 3, a diagonal approximation to J is used.
!!
!! if MITER = 1 or 2, and if the existing value of the Jacobian
!! (as contained in P) is considered acceptable, then a new value of
!! P is reconstructed from the old value.
!!
!! In any case, when MITER
!! is 1 or 2, the P matrix is subjected to LU decomposition in CDRV.
!!
!! P and its LU decomposition are stored (separately) in WK.
!!
!! In addition to variables described previously, communication
!! with DPRJS uses the following:
!!
!! Y
!!
!! : array containing predicted values on entry.
!!
!! FTEM
!!
!! : work array of length N (ACOR in DSTODE).
!!
!! SAVF
!!
!! : array containing f evaluated at predicted y.
!!
!! WK
!!
!! : real work space for matrices.  On output it contains the
!! inverse diagonal matrix if MITER = 3, and P and its sparse
!! LU decomposition if MITER is 1 or 2.
!! Storage of matrix elements starts at WK(3).
!! WK also contains the following matrix-related data:
!! WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
!! WK(2) = H*EL0, saved for later use if MITER = 3.
!!
!! IWK
!!
!! : integer work space for matrix-related data, assumed to
!! be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
!! are assumed to have identical locations.
!!
!! EL0
!!
!! : EL(1) (input).
!!
!! IERPJ
!!
!! : output error flag (in Common).
!!       = 0 if no error.
!!       = 1  if zero pivot found in CDRV.
!!       = 2  if a singular matrix arose with MITER = 3.
!!       = -1 if insufficient storage for CDRV (should not occur here).
!!       = -2 if other error found in CDRV (should not occur here).
!!
!! JCUR
!!
!! : output flag showing status of (approximate) Jacobian matrix:
!!       = 1 to indicate that the Jacobian is now current, or
!!       = 0 to indicate that a saved value was used.
!! This routine also uses other variables in Common.
!-----------------------------------------------------------------------
subroutine dprjs(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wk,Iwk,f,jac)

integer                       :: Neq(*)
real(kind=dp),intent(inout)   :: Y(*)
integer,intent(in)            :: Nyh
real(kind=dp),intent(in)      :: Yh(Nyh,*)
real(kind=dp)                 :: Ewt(*)
real(kind=dp),intent(inout)   :: Ftem(*)
real(kind=dp)                 :: Savf(*)
real(kind=dp),intent(inout)   :: Wk(*)
integer                       :: Iwk(*)
external                      :: f
external                      :: jac

real(kind=dp) :: con, di, fac, hl0, pij, r, r0, rcon, rcont, srur
integer :: i, imul, j, jj, jmax, jmin, jok, k, kmax, kmin, ng

   hl0 = dls1%h*dls1%el0
   con = -hl0
   if ( dls1%miter==3 ) then
   !
   !  If MITER = 3, construct a diagonal approximation to J and P. ---------
      dls1%jcur = 1
      dls1%nje = dls1%nje + 1
      Wk(2) = hl0
      dls1%ierpj = 0
      r = dls1%el0*0.1D0
      do i = 1, dls1%n
         Y(i) = Y(i) + r*(dls1%h*Savf(i)-Yh(i,2))
      enddo
      call f(Neq,dls1%tn,Y,Wk(3))
      dls1%nfe = dls1%nfe + 1
      do i = 1, dls1%n
         r0 = dls1%h*Savf(i) - Yh(i,2)
         di = 0.1D0*r0 - dls1%h*(Wk(i+2)-Savf(i))
         Wk(i+2) = 1.0D0
         if ( abs(r0)>=dls1%uround/Ewt(i) ) then
            if ( abs(di)==0.0D0 ) then
               dls1%ierpj = 2
               return
            else
               Wk(i+2) = 0.1D0*r0/di
            endif
         endif
      enddo
      return
   else
   !  See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
      jok = 1
      if ( dls1%nst==0 .or. dls1%nst>=dlss%nslj+dlss%msbj ) jok = 0
      if ( dls1%icf==1 .and. abs(dls1%rc-1.0D0)<dlss%ccmxj ) jok = 0
      if ( dls1%icf==2 ) jok = 0
      if ( jok==1 ) then
   !
   !  If JOK = 1, reconstruct new P from old P. ----------------------------
         dls1%jcur = 0
         rcon = con/dlss%con0
         rcont = abs(con)/dlss%conmin
         if ( rcont<=dlss%rbig .or. dlss%iplost/=1 ) then
            kmin = Iwk(dlss%ipian)
            do j = 1, dls1%n
               kmax = Iwk(dlss%ipian+j) - 1
               do k = kmin, kmax
                  i = Iwk(dlss%ibjan+k)
                  pij = Wk(dlss%iba+k)
                  if ( i==j ) then
                     pij = pij - 1.0D0
                     if ( abs(pij)<dlss%psmall ) then
                        dlss%iplost = 1
                        dlss%conmin = min(abs(dlss%con0),dlss%conmin)
                     endif
                  endif
                  pij = pij*rcon
                  if ( i==j ) pij = pij + 1.0D0
                  Wk(dlss%iba+k) = pij
               enddo
               kmin = kmax + 1
            enddo
            call wrapup()
            return
         endif
      endif
   !
   !  MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
      dls1%jcur = 1
      dls1%nje = dls1%nje + 1
      dlss%nslj = dls1%nst
      dlss%iplost = 0
      dlss%conmin = abs(con)
      if ( dls1%miter==2 ) then
   !
   !  If MITER = 2, make NGP calls to F to approximate J and P. ------------
         fac = dvnorm(dls1%n,Savf,Ewt)
         r0 = 1000.0D0*abs(dls1%h)*dls1%uround*dls1%n*fac
         if ( r0==0.0D0 ) r0 = 1.0D0
         srur = Wk(1)
         jmin = Iwk(dlss%ipigp)
         do ng = 1, dlss%ngp
            jmax = Iwk(dlss%ipigp+ng) - 1
            do j = jmin, jmax
               jj = Iwk(dlss%ibjgp+j)
               r = max(srur*abs(Y(jj)),r0/Ewt(jj))
               Y(jj) = Y(jj) + r
            enddo
            call f(Neq,dls1%tn,Y,Ftem)
            do j = jmin, jmax
               jj = Iwk(dlss%ibjgp+j)
               Y(jj) = Yh(jj,1)
               r = max(srur*abs(Y(jj)),r0/Ewt(jj))
               fac = -hl0/r
               kmin = Iwk(dlss%ibian+jj)
               kmax = Iwk(dlss%ibian+jj+1) - 1
               do k = kmin, kmax
                  i = Iwk(dlss%ibjan+k)
                  Wk(dlss%iba+k) = (Ftem(i)-Savf(i))*fac
                  if ( i==jj ) Wk(dlss%iba+k) = Wk(dlss%iba+k) + 1.0D0
               enddo
            enddo
            jmin = jmax + 1
         enddo
         dls1%nfe = dls1%nfe + dlss%ngp
      else
   !
   !  If MITER = 1, call JAC, multiply by scalar, and add identity. --------
         kmin = Iwk(dlss%ipian)
         do j = 1, dls1%n
            kmax = Iwk(dlss%ipian+j) - 1
            do i = 1, dls1%n
               Ftem(i) = 0.0D0
            enddo
            call jac(Neq,dls1%tn,Y,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Ftem)
            do k = kmin, kmax
               i = Iwk(dlss%ibjan+k)
               Wk(dlss%iba+k) = Ftem(i)*con
               if ( i==j ) Wk(dlss%iba+k) = Wk(dlss%iba+k) + 1.0D0
            enddo
            kmin = kmax + 1
         enddo
      endif
   endif

   call wrapup()

contains

subroutine wrapup()
   !
   !  Do numerical factorization of P matrix. ------------------------------
   dlss%nlu = dlss%nlu + 1
   dlss%con0 = con
   dls1%ierpj = 0

   do i = 1, dls1%n
      Ftem(i) = 0.0D0
   enddo

   call cdrv(dls1%n,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic), &
    & Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk(dlss%ipa),Ftem,Ftem,dlss%nsp, &
    & Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,2,dlss%iys)

   if ( dlss%iys==0 ) return
   imul = (dlss%iys-1)/dls1%n
   dls1%ierpj = -2
   if ( imul==8 ) dls1%ierpj = 1
   if ( imul==10 ) dls1%ierpj = -1

end subroutine wrapup

end subroutine dprjs
