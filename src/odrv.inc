!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### Name
!!  odrv(3f) [M_odepack] - driver for sparse matrix reordering routines
!!
!!### Description
!!
!!   odrv finds a minimum degree ordering of the rows and columns
!!   of a matrix m stored in (ia,ja,a) format (see below).  for the
!!   reordered matrix, the work and storage required to perform
!!   gaussian elimination is (usually) significantly less.
!!
!!   note.. odrv and its subordinate routines have been modified to
!!   compute orderings for general matrices, not necessarily having any
!!   symmetry.  the miminum degree ordering is computed for the
!!   structure of the symmetric matrix  m + m-transpose.
!!   modifications to the original odrv module have been made in
!!   the coding in subroutine mdi, and in the initial comments in
!!   subroutines odrv and md.
!!
!!   if only the nonzero entries in the upper triangle of m are being
!!   stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
!!   with the diagonal entries placed first in each row.  this is to
!!   ensure that if m(i,j) will be in the upper triangle of m with
!!   respect to the new ordering, then m(i,j) is stored in row i (and
!!   thus m(j,i) is not stored),  whereas if m(i,j) will be in the
!!   strict lower triangle of m, then m(j,i) is stored in row j (and
!!   thus m(i,j) is not stored).
!!
!!
!!### storage of sparse matrices
!!
!!   the nonzero entries of the matrix m are stored row-by-row in the
!!   array a.  to identify the individual nonzero entries in each row,
!!   we need to know in which column each entry lies.  these column
!!   indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!!   ja(k) = j.  to identify the individual rows, we need to know where
!!   each row starts.  these row pointers are stored in the array ia.
!!   i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
!!   and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
!!   the first location following the last element in the last row.
!!   thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
!!   the nonzero entries in the i-th row are stored consecutively in
!!```text
!!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!!```
!!   and the corresponding column indices are stored consecutively in
!!```text
!!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!!```
!!   when the coefficient matrix is symmetric, only the nonzero entries
!!   in the upper triangle need be stored.  for example, the matrix
!!```text
!!             ( 1  0  2  3  0 )
!!             ( 0  4  0  0  0 )
!!         m = ( 2  0  5  6  0 )
!!             ( 3  0  6  7  8 )
!!             ( 0  0  0  8  9 )
!!```
!!   could be stored as
!!```text
!!            - 1  2  3  4  5  6  7  8  9 10 11 12 13
!!         ---+--------------------------------------
!!         ia - 1  4  5  8 12 14
!!         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
!!          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
!!```
!!   or (symmetrically) as
!!```text
!!            - 1  2  3  4  5  6  7  8  9
!!         ---+--------------------------
!!         ia - 1  4  5  7  9 10
!!         ja - 1  3  4  2  3  4  4  5  5
!!          a - 1  2  3  4  5  6  7  8  9          .
!!```
!!
!!### parameters
!!```text
!!    n    - order of the matrix
!!
!!    ia   - integer one-dimensional array containing pointers to delimit
!!           rows in ja and a.  dimension = n+1
!!
!!    ja   - integer one-dimensional array containing the column indices
!!           corresponding to the elements of a.  dimension = number of
!!           nonzero entries in (the upper triangle of) m
!!
!!    a    - real one-dimensional array containing the nonzero entries in
!!           (the upper triangle of) m, stored by rows.  dimension =
!!           number of nonzero entries in (the upper triangle of) m
!!
!!    p    - integer one-dimensional array used to return the permutation
!!           of the rows and columns of m corresponding to the minimum
!!           degree ordering.  dimension = n
!!
!!    ip   - integer one-dimensional array used to return the inverse of
!!           the permutation returned in p.  dimension = n
!!
!!    nsp  - declared dimension of the one-dimensional array isp.  nsp
!!           must be at least  3n+4k,  where k is the number of nonzeroes
!!           in the strict upper triangle of m
!!
!!    isp  - integer one-dimensional array used for working storage.
!!           dimension = nsp
!!
!!    path - integer path specification.  values and their meanings are -
!!             1  find minimum degree ordering only
!!             2  find minimum degree ordering and reorder symmetrically
!!                  stored matrix (used when only the nonzero entries in
!!                  the upper triangle of m are being stored)
!!             3  reorder symmetrically stored matrix as specified by
!!                  input permutation (used when an ordering has already
!!                  been determined and only the nonzero entries in the
!!                  upper triangle of m are being stored)
!!             4  same as 2 but put diagonal entries at start of each row
!!             5  same as 3 but put diagonal entries at start of each row
!!
!!    flag - integer error flag.  values and their meanings are -
!!               0    no errors detected
!!              9n+k  insufficient storage in md
!!             10n+1  insufficient storage in odrv
!!             11n+1  illegal path specification
!!```
!!
!-----------------------------------------------------------------------
subroutine odrv(N,Ia,Ja,A,P,Ip,Nsp,Isp,Path,Flag)

integer                        :: N
integer,        dimension(*)   :: Ia
integer,        dimension(*)   :: Ja
real(kind=dp),  dimension(*)   :: A
integer,        dimension(*)   :: P
integer,        dimension(*)   :: Ip
integer,        intent(in)     :: Nsp
integer,        dimension(*)   :: Isp
integer,        intent(in)     :: Path
integer,        intent(inout)  :: Flag

logical :: dflag
integer :: head, l, max, next, q, tmp, v

   !----initialize error flag and validate path specification
   Flag = 0

   if ( Path<1 .or. 5<Path ) then

      !  ** error -- illegal path specified
      Flag = 11*N + 1
      return

   else

      !----allocate storage and find minimum degree ordering
      if ( (Path-1)*(Path-2)*(Path-4)==0 ) then
         max = (Nsp-N)/2
         v = 1
         l = v + max
         head = l + max
         next = head + N
         if ( max<N ) then
            Flag = 10*N + 1
            return
         endif

         call md(N,Ia,Ja,max,Isp(v),Isp(l),Isp(head),P,Ip,Isp(v),Flag)

         !  ** error -- error detected in md
         if ( Flag/=0 ) return
      endif

      !----allocate storage and symmetrically reorder matrix
      if ( (Path-2)*(Path-3)*(Path-4)*(Path-5)==0 ) then
         tmp = (Nsp+1) - N
         q = tmp - (Ia(N+1)-1)
         if ( q<1 )then
            Flag = 10*N + 1
            return
         endif

         dflag = Path==4 .or. Path==5
         call sro(N,Ip,Ia,Ja,A,Isp(tmp),Isp(q),dflag)
      endif

      return
   endif

   !  ** error -- insufficient storage
   Flag = 10*N + 1

end subroutine odrv
