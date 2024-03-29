!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!  mdm -- form element from uneliminated neighbors of vk
!!
!-----------------------------------------------------------------------
subroutine mdm(Vk,Tail,V,L,Last,Next,Mark)

integer,intent(in)    :: Vk
integer,intent(inout) :: Tail
integer,intent(in)    :: V(*)
integer,intent(inout) :: L(*)
integer,intent(in)    :: Last(*)
integer,intent(in)    :: Next(*)
integer,intent(inout) :: Mark(*)

integer :: b, blp, blpmax, es, lb, ls, s, tag, vb, vs

equivalence (vs,es)
   !
   !----initialize tag and list of uneliminated neighbors
   tag = Mark(Vk)
   Tail = Vk
   !
   !----for each vertex/element vs/es in element list of vk
   ls = L(Vk)
   do
      s = ls
      if ( s==0 ) then
   !
   !----terminate list of uneliminated neighbors
         L(Tail) = 0
         exit
      else
         ls = L(s)
         vs = V(s)
         if ( Next(vs)<0 ) then
   !
   !------if es is active element, then ...
   !--------for each vertex vb in boundary list of element es
            lb = L(es)
            blpmax = Last(es)
            do blp = 1, blpmax
               b = lb
               lb = L(b)
               vb = V(b)
   !
   !----------if vb is untagged vertex, then tag and append to list of
   !----------uneliminated neighbors
               if ( Mark(vb)<tag ) then
                  Mark(vb) = tag
                  L(Tail) = b
                  Tail = b
               endif
            enddo
   !
   !--------mark es inactive
   !
            Mark(es) = tag
         else
   !
   !------if vs is uneliminated vertex, then tag and append to list of
   !------uneliminated neighbors
            Mark(vs) = tag
            L(Tail) = s
            Tail = s
         endif
      endif
   enddo

end subroutine mdm
