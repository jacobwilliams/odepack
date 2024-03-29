!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!
!!### DESCRIPTION
!!  DINTDY computes interpolated values of the K-th derivative of the
!!  dependent variable vector y, and stores it in DKY.  This routine
!!  is called within the package with K = 0 and T = TOUT, but may
!!  also be called by the user for any K up to the current order.
!!  (See detailed instructions in the usage documentation.)
!!
!!  The computed values in DKY are gotten by interpolation using the
!!  Nordsieck history array YH.  This array corresponds uniquely to a
!!  vector-valued polynomial of degree NQCUR or less, and DKY is set
!!  to the K-th derivative of this polynomial at T.
!!  The formula for DKY is:
!!```text
!!               q
!!   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
!!              j=K
!!```
!!  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
!!
!!  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
!!  communicated by COMMON.  The above sum is done in reverse order.
!!  IFLAG is returned negative if either K or T is out of bounds.
!!
! ### SUBSIDIARY
! ### PURPOSE  Interpolate solution derivatives.
! ### TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  XERRWD
! ### COMMON BLOCKS    DLS001
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890503  Minor cosmetic changes.  (FNF)
!     19930809  Renamed to allow single/double precision versions. (ACH)
!     20010418  Reduced size of Common block /DLS001/. (ACH)
!     20031105  Restored 'own' variables to Common block /DLS001/, to
!               enable interrupt/restart feature. (ACH)
!     20050427  Corrected roundoff decrement in TP. (ACH)
!-----------------------------------------------------------------------
subroutine dintdy(T,K,Yh,Nyh,Dky,Iflag)
!
real(kind=dp)               :: T
integer                     :: K
integer,intent(in)          :: Nyh
real(kind=dp),intent(in)    :: Yh(Nyh,*)
real(kind=dp),intent(inout) :: Dky(*)
integer,intent(out)         :: Iflag
!
real(kind=dp)     :: c , r , s , tp
integer           :: i , ic , j , jb , jb2 , jj , jj1 , jp1
character(len=80) :: msg
   !
   Iflag = 0
   if ( K<0 .or. K>dls1%nq ) then
      msg = 'DINTDY-  K (=I1) illegal      '
      call xerrwd(msg,30,51,0,1,K,0,0,0.0D0,0.0D0)
      Iflag = -1
      return
   endif
   tp = dls1%tn - dls1%hu - 100.0D0*dls1%uround*sign(abs(dls1%tn)+abs(dls1%hu),dls1%hu)
   if ( (T-tp)*(T-dls1%tn)>0.0D0 ) then
      msg = 'DINTDY-  T (=R1) illegal      '
      call xerrwd(msg,30,52,0,0,0,0,1,T,0.0D0)
      msg = '      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      call xerrwd(msg,60,52,0,0,0,0,2,tp,dls1%tn)
      Iflag = -2
      return
   endif
   !
   s = (T-dls1%tn)/dls1%h
   ic = 1

   if ( K/=0 ) then
      jj1 = dls1%l - K
      do jj = jj1 , dls1%nq
         ic = ic*jj
      enddo
   endif

   c = ic

   do i = 1 , dls1%n
      Dky(i) = c*Yh(i,dls1%l)
   enddo

   if ( K/=dls1%nq ) then
      jb2 = dls1%nq - K
      do jb = 1 , jb2
         j = dls1%nq - jb
         jp1 = j + 1
         ic = 1
         if ( K/=0 ) then
            jj1 = jp1 - K
            do jj = jj1 , j
               ic = ic*jj
            enddo
         endif
         c = ic
         do i = 1 , dls1%n
            Dky(i) = c*Yh(i,jp1) + s*Dky(i)
         enddo
      enddo
      if ( K==0 ) return
   endif

   r = dls1%h**(-K)

   do i = 1 , dls1%n
      Dky(i) = r*Dky(i)
   enddo

end subroutine dintdy
