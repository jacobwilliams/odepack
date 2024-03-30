!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine manages the solution of the linear system arising from
!! a chord iteration.  It is called if MITER .ne. 0.
!!
!! If MITER is 1 or 2, it calls CDRV to accomplish this.
!!
!! If MITER = 3 it updates the coefficient H*EL0 in the diagonal
!! matrix, and then computes the solution.
!! communication with DSOLSS uses the following variables:
!!
!! WK
!!
!! : real work space containing the inverse diagonal matrix if
!! MITER = 3 and the LU decomposition of the matrix otherwise.
!! Storage of matrix elements starts at WK(3).
!! WK also contains the following matrix-related data:
!!
!!         WK(1) = SQRT(UROUND) (not used here),
!!         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
!!
!! IWK
!!
!! : integer work space for matrix-related data, assumed to
!! be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
!! are assumed to have identical locations.
!!
!! X
!!
!! : the right-hand side vector on input, and the solution vector
!! on output, of length N.
!!
!! TEM
!!
!! : vector of work space of length N, not used in this version.
!!
!! IERSL
!!
!! : output flag (in Common).
!!
!!            IERSL = 0  if no trouble occurred.
!!            IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
!!                       This should never occur and is considered fatal.
!!            IERSL = 1  if a singular matrix arose with MITER = 3.
!!
!-----------------------------------------------------------------------
! This routine also uses other variables in global structures
!-----------------------------------------------------------------------
subroutine dsolss(Wk,Iwk,X,Tem)
!
real(kind=dp),intent(inout) :: Wk(*)
integer                     :: Iwk(*)
real(kind=dp),intent(inout) :: X(*)
real(kind=dp)               :: Tem(*)
!
real(kind=dp) :: di, hl0, phl0, r
integer :: i
   !
   dls1%iersl = 0
   select case (dls1%miter)
   case (3)
   !
      phl0 = Wk(2)
      hl0 = dls1%h*dls1%el0
      Wk(2) = hl0
      if ( hl0/=phl0 ) then
         r = hl0/phl0
         do i = 1, dls1%n
            di = 1.0D0 - r*(1.0D0-1.0D0/Wk(i+2))
            if ( abs(di)==0.0D0 ) then
               dls1%iersl = 1
               return
            else
               Wk(i+2) = 1.0D0/di
            endif
         enddo
      endif
   case default
      call cdrv(dls1%n,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic),Iwk(dlss%ipian),Iwk(dlss%ipjan),Wk(dlss%ipa),X,X,dlss%nsp,&
       & Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,4,dls1%iersl)
      if ( dls1%iersl/=0 ) dls1%iersl = -1
      return
   endselect

   do i = 1, dls1%n
      X(i) = Wk(i+2)*X(i)
   enddo

end subroutine dsolss
