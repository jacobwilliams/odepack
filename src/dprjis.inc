!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! DPRJIS is called to compute and process the matrix
!! P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
!! where r = g(t,y) - A(t,y)*s.
!!
!! J is computed by columns, either by
!! the user-supplied routine JAC if MITER = 1, or by finite differencing
!! if MITER = 2.
!!
!! J is stored in WK, rescaled, and ADDA is called to
!! generate P.
!!
!! The matrix P is subjected to LU decomposition in CDRV.
!! P and its LU decomposition are stored separately in WK.
!!
!! In addition to variables described previously, communication
!! with DPRJIS uses the following:
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
!! : array containing r evaluated at predicted y. On output it
!! contains the residual evaluated at current values of t and y.
!!
!! S
!!
!! : array containing predicted values of dy/dt (SAVF in DSTODI).
!!
!! WK
!!
!! : real work space for matrices.  On output it contains P and
!! its sparse LU decomposition.  Storage of matrix elements
!! starts at WK(3).
!! WK also contains the following matrix-related data.
!! WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
!!
!! IWK
!!
!! : integer work space for matrix-related data, assumed to be
!! equivalenced to WK.  In addition,  WK(IPRSP) and IWK(IPISP)
!! are assumed to have identical locations.
!!
!! EL0
!!
!! : EL(1) (input).
!!
!! IERPJ
!!
!! : output error flag (in COMMON).
!!         =  0 if no error.
!!         =  1 if zero pivot found in CDRV.
!!         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
!!         = -1 if insufficient storage for CDRV (should not occur).
!!         = -2 if other error found in CDRV (should not occur here).
!!
!! JCUR
!!
!! : output flag = 1 to indicate that the Jacobian matrix
!! (or approximation) is now current.
!!
!! This routine also uses other variables in global structures.
!-----------------------------------------------------------------------
subroutine dprjis(Neq,Y,Yh,Nyh,Ewt,Rtem,Savr,S,Wk,Iwk,res,jac,adda)

!
integer, dimension(*) :: Neq
real(kind=dp), intent(inout), dimension(*) :: Y
integer, intent(in) :: Nyh
real(kind=dp), intent(in), dimension(Nyh,*) :: Yh
real(kind=dp), intent(in), dimension(*) :: Ewt
real(kind=dp), intent(inout), dimension(*) :: Rtem
real(kind=dp), dimension(*) :: Savr
real(kind=dp), dimension(*) :: S
real(kind=dp), intent(inout), dimension(*) :: Wk
integer, dimension(*) :: Iwk
external :: res
external :: jac
external :: adda
!
real(kind=dp) :: con, fac, hl0, r, srur

integer :: i, imul, ires, j, jj, jmax, jmin, k, kmax, kmin, ng
!
hl0 = dls1%h*dls1%el0
con = -hl0
dls1%jcur = 1
dls1%nje = dls1%nje + 1
if ( dls1%miter==2 ) then
!
!  If MITER = 2, make NGP + 1 calls to RES to approximate J and P. ------
   ires = -1
   call res(Neq,dls1%tn,Y,S,Savr,ires)
   dls1%nfe = dls1%nfe + 1
   if ( ires>1 ) then
!  Error return for IRES = 2 or IRES = 3 return from RES. ---------------
      dls1%ierpj = ires
      return
   else
      srur = Wk(1)
      jmin = Iwk(dlss%ipigp)
      do ng = 1, dlss%ngp
         jmax = Iwk(dlss%ipigp+ng) - 1
         do j = jmin, jmax
            jj = Iwk(dlss%ibjgp+j)
            r = max(srur*abs(Y(jj)),0.01D0/Ewt(jj))
            Y(jj) = Y(jj) + r
         enddo
         call res(Neq,dls1%tn,Y,S,Rtem,ires)
         dls1%nfe = dls1%nfe + 1
         if ( ires>1 ) then
            dls1%ierpj = ires
            return
         else
            do j = jmin, jmax
               jj = Iwk(dlss%ibjgp+j)
               Y(jj) = Yh(jj,1)
               r = max(srur*abs(Y(jj)),0.01D0/Ewt(jj))
               fac = -hl0/r
               kmin = Iwk(dlss%ibian+jj)
               kmax = Iwk(dlss%ibian+jj+1) - 1
               do k = kmin, kmax
                  i = Iwk(dlss%ibjan+k)
                  Rtem(i) = (Rtem(i)-Savr(i))*fac
               enddo
               call adda(Neq,dls1%tn,Y,jj,Iwk(dlss%ipian),Iwk(dlss%ipjan),Rtem)
               do k = kmin, kmax
                  i = Iwk(dlss%ibjan+k)
                  Wk(dlss%iba+k) = Rtem(i)
               enddo
            enddo
            jmin = jmax + 1
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
else
!
!  If MITER = 1, call RES, then call JAC and ADDA for each column. ------
   ires = 1
   call res(Neq,dls1%tn,Y,S,Savr,ires)
   dls1%nfe = dls1%nfe + 1
   if ( ires>1 ) then
      dls1%ierpj = ires
      return
   else
      kmin = Iwk(dlss%ipian)
      do j = 1, dls1%n
         kmax = Iwk(dlss%ipian+j) - 1
         do i = 1, dls1%n
            Rtem(i) = 0.0D0
         enddo
         call jac(Neq,dls1%tn,Y,S,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Rtem)
         do i = 1, dls1%n
            Rtem(i) = Rtem(i)*con
         enddo
         call adda(Neq,dls1%tn,Y,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Rtem)
         do k = kmin, kmax
            i = Iwk(dlss%ibjan+k)
            Wk(dlss%iba+k) = Rtem(i)
         enddo
         kmin = kmax + 1
      enddo
   endif
endif
!
!  Do numerical factorization of P matrix. ------------------------------
dlss%nlu = dlss%nlu + 1
dls1%ierpj = 0
do i = 1, dls1%n
   Rtem(i) = 0.0D0
enddo
call cdrv(dls1%n,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic), &
 & Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk(dlss%ipa),Rtem,Rtem,dlss%nsp, &
 & Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,2,dlss%iys)
if ( dlss%iys==0 ) return
imul = (dlss%iys-1)/dls1%n
dls1%ierpj = -2
if ( imul==8 ) dls1%ierpj = 1
if ( imul==10 ) dls1%ierpj = -1

end subroutine dprjis
