!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### Name
!!  sro -- symmetric reordering of sparse symmetric matrix
!!
!!### Description
!!
!!    the nonzero entries of the matrix m are assumed to be stored
!!    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
!!    are stored if i ne j).
!!
!!    sro does not rearrange the order of the rows, but does move
!!    nonzeroes from one row to another to ensure that if m(i,j) will be
!!    in the upper triangle of m with respect to the new ordering, then
!!    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
!!    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
!!    stored in row j (and thus m(i,j) is not stored).
!!
!!
!!### additional parameters
!!
!!   q
!!
!!   : integer one-dimensional work array.  dimension = n
!!
!!   r
!!
!!   : integer one-dimensional work array.  dimension = number of
!!   nonzero entries in the upper triangle of m
!!
!!   dflag
!!
!!   : logical variable.  if dflag = .true., then store nonzero
!!   diagonal elements at the beginning of the row
!!
!-----------------------------------------------------------------------
subroutine sro(N,Ip,Ia,Ja,A,Q,R,Dflag)
integer,intent(in)          :: N
integer,intent(in)          :: Ip(*)
integer,intent(inout)       :: Ia(*)
integer,intent(inout)       :: Ja(*)
real(kind=dp),intent(inout) :: A(*)
integer,intent(inout)       :: Q(*)
integer,intent(inout)       :: R(*)
logical,intent(in)          :: Dflag
!
real(kind=dp) :: ak
integer       :: i, ilast, j, jak, jdummy, jmax, jmin, k
!
!--phase 1 -- find row in which to store each nonzero
!----initialize count of nonzeroes to be stored in each row
do i = 1 , N
   Q(i) = 0
enddo
!----for each nonzero element a(j)
do i = 1 , N
   jmin = Ia(i)
   jmax = Ia(i+1) - 1
   if ( jmin<=jmax ) then
      do j = jmin , jmax
!--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
         k = Ja(j)
         if ( Ip(k)<Ip(i) ) Ja(j) = i
         if ( Ip(k)>=Ip(i) ) k = i
         R(j) = k
!--------... and increment count of nonzeroes (=q(r(j)) in that row
         Q(k) = Q(k) + 1
      enddo
   endif
enddo
!
!--phase 2 -- find new ia and permutation to apply to (ja,a)
!----determine pointers to delimit rows in permuted (ja,a)
do i = 1 , N
   Ia(i+1) = Ia(i) + Q(i)
   Q(i) = Ia(i+1)
enddo
!
!----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
!----for each nonzero element (in reverse order)
ilast = 0
jmin = Ia(1)
jmax = Ia(N+1) - 1
j = jmax
do jdummy = jmin , jmax
   i = R(j)
   if ( .not.Dflag .or. Ja(j)/=i .or. i==ilast ) then
!
!------put (off-diagonal) nonzero in last unused location in row
      Q(i) = Q(i) - 1
      R(j) = Q(i)
   else
!
!------if dflag, then put diagonal nonzero at beginning of row
      R(j) = Ia(i)
      ilast = i
   endif
!
   j = j - 1
enddo
!
!--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
do j = jmin , jmax
   do while ( R(j)/=j )
      k = R(j)
      R(j) = R(k)
      R(k) = k
      jak = Ja(k)
      Ja(k) = Ja(j)
      Ja(j) = jak
      ak = A(k)
      A(k) = A(j)
      A(j) = ak
   enddo
enddo
!
end subroutine sro
