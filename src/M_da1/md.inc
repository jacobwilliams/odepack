!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!  md -- minimum degree algorithm (based on element model)
!!
!!### description
!!
!!    md finds a minimum degree ordering of the rows and columns of a
!!    general sparse matrix m stored in (ia,ja,a) format.
!!    when the structure of m is nonsymmetric, the ordering is that
!!    obtained for the symmetric matrix  m + m-transpose.
!!
!!
!!### additional parameters
!!
!!    MAX
!!
!!    : declared dimension of the one-dimensional arrays v and l.
!!    max must be at least  n+2k,  where k is the number of
!!    nonzeroes in the strict upper triangle of m + m-transpose
!!
!!    V
!!
!!    : integer one-dimensional work array.  dimension = max
!!
!!    L
!!
!!    : integer one-dimensional work array.  dimension = max
!!
!!    HEAD
!!
!!    : integer one-dimensional work array.  dimension = n
!!
!!    LAST
!!
!!    : integer one-dimensional array used to return the permutation
!!    of the rows and columns of m corresponding to the minimum
!!    degree ordering.  dimension = n
!!
!!    NEXT
!!
!!    : integer one-dimensional array used to return the inverse of
!!    the permutation returned in last.  dimension = n
!!
!!    MARK
!!
!!    : integer one-dimensional work array (may be the same as v).
!!    dimension = n
!!
!!    FLAG
!!
!!    : integer error flag.  values and their meanings are -
!!
!!             0     no errors detected
!!             9n+k  insufficient storage in md
!!
!!
!!### definitions of internal parameters
!!```text
!!    ---------+---------------------------------------------------------
!!    v(s)     - value field of list entry
!!    ---------+---------------------------------------------------------
!!    l(s)     - link field of list entry  (0 =) end of list)
!!    ---------+---------------------------------------------------------
!!    l(vi)    - pointer to element list of uneliminated vertex vi
!!    ---------+---------------------------------------------------------
!!    l(ej)    - pointer to boundary list of active element ej
!!    ---------+---------------------------------------------------------
!!    head(d)  - vj =) vj head of d-list d
!!             -  0 =) no vertex in d-list d
!!
!!
!!             -                  vi uneliminated vertex
!!             -          vi in ek           -       vi not in ek
!!    ---------+-----------------------------+---------------------------
!!    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
!!             -                             -  0 =) vi tail of d-list
!!    ---------+-----------------------------+---------------------------
!!    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
!!             --vk =) compute degree        - vj =) vj last in d-list
!!             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
!!             -  0 =) do not compute degree -
!!    ---------+-----------------------------+---------------------------
!!    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
!!
!!
!!             -                   vi eliminated vertex
!!             -      ei active element      -           otherwise
!!    ---------+-----------------------------+---------------------------
!!    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
!!             -       to be eliminated      -       to be eliminated
!!    ---------+-----------------------------+---------------------------
!!    last(vi) -  m =) size of ei = m        - undefined
!!    ---------+-----------------------------+---------------------------
!!    mark(vi) - -m =) overlap count of ei   - undefined
!!             -       with ek = m           -
!!             - otherwise nonnegative tag   -
!!             -       .lt. mark(vk)         -
!!```
!!
!-----------------------------------------------------------------------
subroutine md(N,Ia,Ja,Max,V,L,Head,Last,Next,Mark,Flag)

integer               :: N
integer               :: Ia(*)
integer               :: Ja(*)
integer               :: Max
integer               :: V(*)
integer               :: L(*)
integer,intent(inout) :: Head(*)
integer,intent(inout) :: Last(*)
integer,intent(inout) :: Next(*)
integer               :: Mark(*)
integer               :: Flag

integer :: dmin, ek, k, tag, tail, vk

equivalence (vk,ek)

   !----initialization
   tag = 0
   call mdi(N,Ia,Ja,Max,V,L,Head,Last,Next,Mark,tag,Flag)
   if ( Flag/=0 ) return
   !
   k = 0
   dmin = 1
   !
   !----while  k .lt. n  do
   do while ( k<N )
   !
   !------search for vertex of minimum degree
      do while ( Head(dmin)<=0 )
         dmin = dmin + 1
      enddo
   !
   !------remove vertex vk of minimum degree from degree list
      vk = Head(dmin)
      Head(dmin) = Next(vk)
      if ( Head(dmin)>0 ) Last(Head(dmin)) = -dmin
   !
   !------number vertex vk, adjust tag, and tag vk
      k = k + 1
      Next(vk) = -k
      Last(ek) = dmin - 1
      tag = tag + Last(ek)
      Mark(vk) = tag
   !
   !------form element ek from uneliminated neighbors of vk
      call mdm(vk,tail,V,L,Last,Next,Mark)
   !
   !------purge inactive elements and do mass elimination
      call mdp(k,ek,tail,V,L,Head,Last,Next,Mark)
   !
   !------update degrees of uneliminated vertices in ek
   !
      call mdu(ek,dmin,V,L,Head,Last,Next,Mark)
   enddo
   !
   !----generate inverse permutation from permutation
   do k = 1, N
      Next(k) = -Next(k)
      Last(Next(k)) = k
   enddo

end subroutine md
