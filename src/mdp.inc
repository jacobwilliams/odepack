!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!  mdp -- purge inactive elements and do mass elimination
!!
!-----------------------------------------------------------------------
subroutine mdp(K,Ek,Tail,V,L,Head,Last,Next,Mark)
!
integer,intent(inout) :: K
integer,intent(in)    :: Ek
integer,intent(inout) :: Tail
integer,intent(inout) :: V(*)
integer,intent(inout) :: L(*)
integer,intent(out)   :: Head(*)
integer,intent(inout) :: Last(*)
integer,intent(inout) :: Next(*)
integer,intent(inout) :: Mark(*)
!
integer :: es, evi, i, ilp, ilpmax, li, ls, lvi, s, tag, vi
integer :: free
   free=-1
   !
   !----initialize tag
   tag = Mark(Ek)
   !
   !----for each vertex vi in ek
   li = Ek
   ilpmax = Last(Ek)
   if ( ilpmax>0 ) then
      do ilp = 1, ilpmax
         i = li
         li = L(i)
         vi = V(li)
         !
         !------remove vi from degree list
         if ( Last(vi)/=0 ) then
            if ( Last(vi)>0 ) then
               Next(Last(vi)) = Next(vi)
            else
               Head(-Last(vi)) = Next(vi)
            endif
            if ( Next(vi)>0 ) Last(Next(vi)) = Last(vi)
         endif
         !
         !------remove inactive items from element list of vi
         ls = vi
         do
            s = ls
            ls = L(s)
            if ( ls==0 ) then
               !
               !------if vi is interior vertex, then remove from list and eliminate
               lvi = L(vi)
               if ( lvi/=0 ) then
                  !
                  !------else ...
                  !--------classify vertex vi
                  if ( L(lvi)/=0 ) then
                     !
                     !----------else mark vi to compute degree
                     Last(vi) = -Ek
                  else
                     evi = V(lvi)
                     if ( Next(evi)>=0 ) then
                        Last(vi) = -Ek
                     elseif ( Mark(evi)<0 ) then
                        !
                        !----------else if vi is duplicate vertex, then mark as such and adjust
                        !----------overlap count for corresponding element
                        Last(vi) = 0
                        Mark(evi) = Mark(evi) - 1
                     else
                       !
                       !----------if vi is prototype vertex, then mark as such, initialize
                       !----------overlap count for corresponding element, and move vi to end
                       !----------of boundary list
                        Last(vi) = evi
                        Mark(evi) = -1
                        L(Tail) = li
                        Tail = li
                        L(i) = L(li)
                        li = i
                     endif
                  endif
                  !
                  !--------insert ek in element list of vi
                  if(free > 0)then
                     V(free) = Ek
                     L(free) = L(vi)
                     L(vi) = free
                  else
                     stop '*mdp* is uninitialized but being used'
                  endif
               else
                  L(i) = L(li)
                  li = i

                  K = K + 1
                  Next(vi) = -K
                  Last(Ek) = Last(Ek) - 1
               endif
               exit
            else
               es = V(ls)
               if ( Mark(es)>=tag ) then
                  free = ls
                  L(s) = L(ls)
                  ls = s
               endif
            endif
         enddo
      enddo
   endif
   !
   !----terminate boundary list
   L(Tail) = 0

end subroutine mdp
