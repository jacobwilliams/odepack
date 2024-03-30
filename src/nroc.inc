!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!
!!               yale sparse matrix package - nonsymmetric codes
!!                    solving the system of equations mx = b
!!
!!### i.   calling sequences
!!
!!         the coefficient matrix can be processed by an ordering routine
!!    (e.g., to reduce fillin or ensure numerical stability) before using
!!    the remaining subroutines.  if no reordering is done, then set
!!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
!!    is used, then nroc should be used to reorder the coefficient matrix
!!
!!### the calling sequence is --
!!```text
!!        (       (matrix ordering))
!!        (nroc   (matrix reordering))
!!         nsfc   (symbolic factorization to determine where fillin will
!!                  occur during numeric factorization)
!!         nnfc   (numeric factorization into product ldu of unit lower
!!                  triangular matrix l, diagonal matrix d, and unit
!!                  upper triangular matrix u, and solution of linear
!!                  system)
!!         nnsc   (solution of linear system for additional right-hand
!!                  side using ldu factorization from nnfc)
!!```
!!    (if only one system of equations is to be solved, then the
!!    subroutine trk should be used.)
!!
!!### ii.  storage of sparse matrices
!!         the nonzero entries of the coefficient matrix m are stored
!!
!!    row-by-row in the array a.  to identify the individual nonzero
!!    entries in each row, we need to know in which column each entry
!!    lies.  the column indices which correspond to the nonzero entries
!!    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!!    ja(k) = j.  in addition, we need to know where each row starts and
!!    how long it is.  the index positions in ja and a where the rows of
!!    m begin are stored in the array ia.  i.e., if m(i,j) is the first
!!    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
!!    ia(i) = k.  moreover, the index in ja and a of the first location
!!    following the last element in the last row is stored in ia(n+1).
!!    thus, the number of entries in the i-th row is given by
!!    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!!    consecutively in
!!
!!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!!
!!    and the corresponding column indices are stored consecutively in
!!
!!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!!
!!    for example, the 5 by 5 matrix
!!```text
!!                ( 1. 0. 2. 0. 0.)
!!                ( 0. 3. 0. 0. 0.)
!!            m = ( 0. 4. 5. 6. 0.)
!!                ( 0. 0. 0. 7. 0.)
!!                ( 0. 0. 0. 8. 9.)
!!```
!!    would be stored as
!!```text
!!               - 1  2  3  4  5  6  7  8  9
!!            ---+--------------------------
!!            ia - 1  3  4  7  8 10
!!            ja - 1  3  2  2  3  4  4  4  5
!!             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!!```
!!
!!         the strict upper (lower) triangular portion of the matrix
!!    u (l) is stored in a similar fashion using the arrays  iu, ju, u
!!    (il, jl, l)  except that an additional array iju (ijl) is used to
!!    compress storage of ju (jl) by allowing some sequences of column
!!    (row) indices to used for more than one row (column)  (n.b., l is
!!    stored by columns).  iju(k) (ijl(k)) points to the starting
!!    location in ju (jl) of entries for the kth row (column).
!!    compression in ju (jl) occurs in two ways.  first, if a row
!!    (column) i was merged into the current row (column) k, and the
!!    number of elements merged in from (the tail portion of) row
!!    (column) i is the same as the final length of row (column) k, then
!!    the kth row (column) and the tail of row (column) i are identical
!!    and iju(k) (ijl(k)) points to the start of the tail.  second, if
!!    some tail portion of the (k-1)st row (column) is identical to the
!!    head of the kth row (column), then iju(k) (ijl(k)) points to the
!!    start of that tail portion.  for example, the nonzero structure of
!!    the strict upper triangular part of the matrix
!!```text
!!            d 0 x x x
!!            0 d 0 x x
!!            0 0 d x 0
!!            0 0 0 d x
!!            0 0 0 0 d
!!```
!!    would be represented as
!!```text
!!                - 1 2 3 4 5 6
!!            ----+------------
!!             iu - 1 4 6 7 8 8
!!             ju - 3 4 5 4
!!            iju - 1 2 4 3           .
!!```
!!    the diagonal entries of l and u are assumed to be equal to one and
!!    are not stored.  the array d contains the reciprocals of the
!!    diagonal entries of the matrix d.
!!
!!### iii. additional storage savings
!!         in nsfc, r and ic can be the same array in the calling
!!    sequence if no reordering of the coefficient matrix has been done.
!!
!!         in nnfc, r, c, and ic can all be the same array if no
!!    reordering has been done.  if only the rows have been reordered,
!!    then c and ic can be the same array.  if the row and column
!!    orderings are the same, then r and c can be the same array.  z and
!!    row can be the same array.
!!
!!         in nnsc or nntc, r and c can be the same array if no
!!    reordering has been done or if the row and column orderings are the
!!    same.  z and b can be the same array.  however, then b will be
!!    destroyed.
!!
!!### iv.  parameters
!!         following is a list of parameters to the programs.  names are
!!    uniform among the various subroutines.  class abbreviations are --
!!```text
!!       n - integer variable
!!       f - real variable
!!       v - supplies a value to a subroutine
!!       r - returns a result from a subroutine
!!       i - used internally by a subroutine
!!       a - array
!!
!! class - parameter
!! ------+----------
!! fva   - a     - nonzero entries of the coefficient matrix m, stored
!!       -           by rows.
!!       -           size = number of nonzero entries in m.
!! fva   - b     - right-hand side b.
!!       -           size = n.
!! nva   - c     - ordering of the columns of m.
!!       -           size = n.
!! fvra  - d     - reciprocals of the diagonal entries of the matrix d.
!!       -           size = n.
!! nr    - flag  - error flag.  values and their meanings are --
!!       -            0     no errors detected
!!       -            n+k   null row in a  --  row = k
!!       -           2n+k   duplicate entry in a  --  row = k
!!       -           3n+k   insufficient storage for jl  --  row = k
!!       -           4n+1   insufficient storage for l
!!       -           5n+k   null pivot  --  row = k
!!       -           6n+k   insufficient storage for ju  --  row = k
!!       -           7n+1   insufficient storage for u
!!       -           8n+k   zero pivot  --  row = k
!! nva   - ia    - pointers to delimit the rows of a.
!!       -           size = n+1.
!! nvra  - ijl   - pointers to the first element in each column in jl,
!!       -           used to compress storage in jl.
!!       -           size = n.
!! nvra  - iju   - pointers to the first element in each row in ju, used
!!       -           to compress storage in ju.
!!       -           size = n.
!! nvra  - il    - pointers to delimit the columns of l.
!!       -           size = n+1.
!! nvra  - iu    - pointers to delimit the rows of u.
!!       -           size = n+1.
!! nva   - ja    - column numbers corresponding to the elements of a.
!!       -           size = size of a.
!! nvra  - jl    - row numbers corresponding to the elements of l.
!!       -           size = jlmax.
!! nv    - jlmax - declared dimension of jl.  jlmax must be larger than
!!       -           the number of nonzeros in the strict lower triangle
!!       -           of m plus fillin minus compression.
!! nvra  - ju    - column numbers corresponding to the elements of u.
!!       -           size = jumax.
!! nv    - jumax - declared dimension of ju.  jumax must be larger than
!!       -           the number of nonzeros in the strict upper triangle
!!       -           of m plus fillin minus compression.
!! fvra  - l     - nonzero entries in the strict lower triangular portion
!!       -           of the matrix l, stored by columns.
!!       -           size = lmax.
!! nv    - lmax  - declared dimension of l.  lmax must be larger than
!!       -           the number of nonzeros in the strict lower triangle
!!       -           of m plus fillin  (il(n+1)-1 after nsfc).
!! nv    - n     - number of variables/equations.
!! nva   - r     - ordering of the rows of m.
!!       -           size = n.
!! fvra  - u     - nonzero entries in the strict upper triangular portion
!!       -           of the matrix u, stored by rows.
!!       -           size = umax.
!! nv    - umax  - declared dimension of u.  umax must be larger than
!!       -           the number of nonzeros in the strict upper triangle
!!       -           of m plus fillin  (iu(n+1)-1 after nsfc).
!! fra   - z     - solution x.
!!       -           size = n.
!!
!!```
!!       ----------------------------------------------------------------
!!
!!###  subroutine nroc
!!###  reorders rows of a, leaving row order unchanged
!!
!!
!!       input parameters.. n, ic, ia, ja, a
!!       output parameters.. ja, a, flag
!!
!!       parameters used internally..
!! nia   - p     - at the kth step, p is a linked list of the reordered
!!       -           column indices of the kth row of a.  p(n+1) points
!!       -           to the first entry in the list.
!!       -           size = n+1.
!! nia   - jar   - at the kth step,jar contains the elements of the
!!       -           reordered column indices of a.
!!       -           size = n.
!! fia   - ar    - at the kth step, ar contains the elements of the
!!       -           reordered row of a.
!!       -           size = n.
!!
!-----------------------------------------------------------------------
subroutine nroc(N,Ic,Ia,Ja,A,Jar,Ar,P,Flag)
!
!      real  a(*), ar(*)
integer , intent(in) :: N
integer , intent(in) , dimension(*) :: Ic
integer , intent(in) , dimension(*) :: Ia
integer , intent(inout) , dimension(*) :: Ja
real(kind=dp) , intent(inout) , dimension(*) :: A
integer , intent(inout) , dimension(*) :: Jar
real(kind=dp) , intent(inout) , dimension(*) :: Ar
integer , intent(inout) , dimension(*) :: P
integer , intent(out) :: Flag
!
integer :: i , j , jmax , jmin , k , newj
   !
   !   ******  for each nonempty row  *******************************
   do k = 1 , N
      jmin = Ia(k)
      jmax = Ia(k+1) - 1
      if ( jmin<=jmax ) then
         P(N+1) = N + 1
         !   ******  insert each element in the list  *********************
         do j = jmin , jmax
            newj = Ic(Ja(j))
            i = N + 1
            do while ( P(i)<newj )
               i = P(i)
            enddo
            if ( P(i)==newj ) then
               !
               !  ** error.. duplicate entry in a
               Flag = N + k
            return
            else
               P(newj) = P(i)
               P(i) = newj
               Jar(newj) = Ja(j)
               Ar(newj) = A(j)
            endif
         enddo
         !   ******  replace old row in ja and a  *************************
         i = N + 1
         do j = jmin , jmax
            i = P(i)
            Ja(j) = Jar(i)
            A(j) = Ar(i)
         enddo
      endif
   enddo
   Flag = 0
end subroutine nroc
