!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### DESCRIPTION
!!
!!  This routine manages the solution of the linear system arising from
!!  a chord iteration.  It is called if MITER .ne. 0.
!!
!!  If MITER is 1 or 2, it calls DGESL to accomplish this.
!!
!!  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
!!  matrix, and then computes the solution.
!!
!!  If MITER is 4 or 5, it calls DGBSL.
!!
!!### OPTIONS
!!
!!  Communication with DSOLSY uses the following variables:
!!
!!  WM
!!
!!  : real work space containing the inverse diagonal matrix if
!!  MITER = 3 and the LU decomposition of the matrix otherwise.
!!  Storage of matrix elements starts at WM(3).
!!  WM also contains the following matrix-related data:
!!
!!          WM(1) = SQRT(UROUND) (not used here),
!!          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
!!
!!  IWM
!!
!!  : integer work space containing pivot information, starting at
!!  IWM(21), if MITER is 1, 2, 4, or 5.
!!
!!  IWM also contains band parameters ML = IWM(1)
!!  and MU = IWM(2) if MITER is 4 or 5.
!!
!!  X
!!
!!  : the right-hand side vector on input, and the solution vector
!!  on output, of length N.
!!
!!  TEM
!!
!!  : vector of work space of length N, not used in this version.
!!
!!  IERSL
!!
!!  : output flag (in COMMON).  IERSL = 0 if no trouble occurred.
!!
!!  IERSL
!!
!!  : 1 if a singular matrix arose with MITER = 3.
!!
!-----------------------------------------------------------------------
! ### SUBSIDIARY
! ### PURPOSE  ODEPACK linear system solver.
! ### TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  DGBSL, DGESL
! ### COMMON BLOCKS    DLS001
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890503  Minor cosmetic changes.  (FNF)
!     19930809  Renamed to allow single/double precision versions. (ACH)
!     20010418  Reduced size of Common block /DLS001/. (ACH)
!     20031105  Restored 'own' variables to Common block /DLS001/, to
!             enable interrupt/restart feature. (ACH)
! ### END PROLOGUE  DSOLSY
!-----------------------------------------------------------------------
!  This routine also uses the COMMON variables EL0, H, MITER, and N.
!-----------------------------------------------------------------------
subroutine dsolsy(Wm,Iwm,X,Tem)

real(kind=dp),intent(inout) :: Wm(*)
integer                     :: Iwm(*)
real(kind=dp),intent(inout) :: X(*)
real(kind=dp)               :: Tem(*)

real(kind=dp) :: di, hl0, phl0, r
integer :: i, meband, ml, mu

   dls1%iersl = 0
   select case (dls1%miter)
   case (3)
      phl0 = Wm(2)
      hl0 = dls1%h*dls1%el0
      Wm(2) = hl0
      if ( hl0/=phl0 ) then
         r = hl0/phl0
         do i = 1, dls1%n
            di = 1.0D0 - r*(1.0D0-1.0D0/Wm(i+2))
            if ( abs(di)==0.0D0 ) then
               dls1%iersl = 1
               return
            else
               Wm(i+2) = 1.0D0/di
            endif
         enddo
      endif
      do i = 1, dls1%n
         X(i) = Wm(i+2)*X(i)
      enddo
   case (4,5)
      ml = Iwm(1)
      mu = Iwm(2)
      meband = 2*ml + mu + 1
      call dgbsl(Wm(3),meband,dls1%n,ml,mu,Iwm(21),X,0)
   case default
      call dgesl(Wm(3),dls1%n,dls1%n,Iwm(21),X,0)
   endselect
end subroutine dsolsy
