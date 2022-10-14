!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This subroutine constructs groupings of the column indices of
!! the Jacobian matrix, used in the numerical evaluation of the
!! Jacobian by finite differences.
!!
!!### Input:
!!
!! N
!!
!! : the order of the matrix.
!!
!! IA,JA
!!
!! : sparse structure descriptors of the matrix by rows.
!!
!! MAXG
!!
!! : length of available storage in the IGP array.
!!
!!### Output:
!!
!! NGRP
!!
!! : number of groups.
!!
!! JGP
!!
!! : array of length N containing the column indices by groups.
!!
!! IGP
!!
!! : pointer array of length NGRP + 1 to the locations in JGP
!! of the beginning of each group.
!!
!! IER
!!
!! : error indicator.  IER = 0 if no error occurred, or 1 if
!! MAXG was insufficient.
!!
!! INCL and JDONE are working arrays of length N.
!-----------------------------------------------------------------------
subroutine jgroup(N,Ia,Ja,Maxg,Ngrp,Igp,Jgp,Incl,Jdone,Ier)
!
integer , intent(in) :: N
integer , intent(in) , dimension(*) :: Ia
integer , intent(in) , dimension(*) :: Ja
integer , intent(in) :: Maxg
integer , intent(out) :: Ngrp
integer , intent(inout) , dimension(*) :: Igp
integer , intent(out) , dimension(*) :: Jgp
integer , intent(inout) , dimension(*) :: Incl
integer , intent(inout) , dimension(*) :: Jdone
integer , intent(out) :: Ier
!
integer :: i , j , k , kmax , kmin , ncol , ng
!
   Ier = 0
   do j = 1 , N
      Jdone(j) = 0
   enddo
   ncol = 1
   do ng = 1 , Maxg
      Igp(ng) = ncol
      do i = 1 , N
         Incl(i) = 0
      enddo
      TEST: do j = 1 , N
   ! Reject column J if it is already in a group.--------------------------
         if ( Jdone(j)/=1 ) then
            kmin = Ia(j)
            kmax = Ia(j+1) - 1
            do k = kmin , kmax
   ! Reject column J if it overlaps any column already in this group.------
               i = Ja(k)
               if ( Incl(i)==1 ) cycle TEST
            enddo
   ! Accept column J into group NG.----------------------------------------
            Jgp(ncol) = j
            ncol = ncol + 1
            Jdone(j) = 1
            do k = kmin , kmax
               i = Ja(k)
               Incl(i) = 1
            enddo
         endif
      enddo TEST
   ! Stop if this group is empty (grouping is complete).-------------------
      if ( ncol==Igp(ng) )then
         Ngrp = ng - 1
         return
      endif
   enddo
   ! Error return if not all columns were chosen (MAXG too small).---------
   if ( ncol<=N ) then
      Ier = 1
      return
   else
      ng = Maxg
   endif

   Ngrp = ng - 1
end subroutine jgroup
