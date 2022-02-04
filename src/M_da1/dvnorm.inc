!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!! dvnorm(3f) - [M_odepack] Weighted root-mean-square vector norm.
!!
!!### DESCRIPTION
!!  This function routine computes the weighted root-mean-square norm
!!  of the vector of length N contained in the array V, with weights
!!  contained in the array W of length N:
!!
!!        DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
!!
!-----------------------------------------------------------------------
! ### SUBSIDIARY
! ### TYPE      real(kind=dp) (SVNORM-S, DVNORM-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     19791129  DATE WRITTEN
!     19890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     19890503  Minor cosmetic changes.  (FNF)
!     19930809  Renamed to allow single/real(kind=dp) versions. (ACH)
!-----------------------------------------------------------------------
pure function dvnorm (n, v, w)
integer,intent(in)       :: n
real(kind=dp),intent(in) :: v(n)
real(kind=dp),intent(in) :: w(n)
real(kind=dp)            :: dvnorm

integer                  :: i
real(kind=dp)            :: sum

   sum = 0.0d0
   do i = 1,n
       sum = sum + (v(i)*w(i))**2
   enddo
   dvnorm = sqrt(sum/n)

end function dvnorm
