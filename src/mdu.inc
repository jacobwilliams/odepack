!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!  mdu -- update degrees of uneliminated vertices in ek
!!
!-----------------------------------------------------------------------
subroutine mdu(Ek,Dmin,V,L,Head,Last,Next,Mark)
!
integer,intent(in)    :: Ek
integer,intent(inout) :: Dmin
integer,intent(in)    :: V(*)
integer,intent(in)    :: L(*)
integer,intent(inout) :: Head(*)
integer,intent(inout) :: Last(*)
integer,intent(inout) :: Next(*)
integer,intent(inout) :: Mark(*)

integer :: b, blp, blpmax, dvi, es, evi, i, ilp, ilpmax, s, tag, vb, vi, vs

equivalence (vs,es)
   !
   !----initialize tag
   tag = Mark(Ek) - Last(Ek)
   !
   !----for each vertex vi in ek
   i = Ek
   ilpmax = Last(Ek)
   if ( ilpmax>0 ) then
      MAIN: do ilp = 1, ilpmax
         i = L(i)
         vi = V(i)
         if ( Last(vi)<0 ) then
   !
   !------if vi neither prototype nor duplicate vertex, then merge elements
   !------to compute degree
            tag = tag + 1
            dvi = Last(Ek)
   !
   !--------for each vertex/element vs/es in element list of vi
            s = L(vi)
            do
               s = L(s)
               if ( s==0 ) exit
               vs = V(s)
               if ( Next(vs)>=0 ) then
   !
   !----------if vs is uneliminated vertex, then tag and adjust degree
                  Mark(vs) = tag
                  dvi = dvi + 1
   !
   !----------if es is active element, then expand
   !------------check for outmatched vertex
               elseif ( Mark(es)<0 ) then
   !
   !------else if vi is outmatched vertex, then adjust overlaps but do not
   !------compute degree
                  Last(vi) = 0
                  Mark(es) = Mark(es) - 1
                  do
                     s = L(s)
                     if ( s==0 ) cycle MAIN
                     es = V(s)
                     if ( Mark(es)<0 ) Mark(es) = Mark(es) - 1
                  enddo
               else
   !
   !------------for each vertex vb in es
                  b = es
                  blpmax = Last(es)
                  do blp = 1, blpmax
                     b = L(b)
                     vb = V(b)
   !
   !--------------if vb is untagged, then tag and adjust degree
                     if ( Mark(vb)<tag ) then
                        Mark(vb) = tag
                        dvi = dvi + 1
                     endif
   !
                  enddo
               endif
            enddo
         elseif ( Last(vi)==0 ) then
            cycle
         else
   !
   !------else if vi is prototype vertex, then calculate degree by
   !------inclusion/exclusion and reset overlap count
            evi = Last(vi)
            dvi = Last(Ek) + Last(evi) + Mark(evi)
            Mark(evi) = 0
         endif
   !
   !------insert vi in appropriate degree list
         Next(vi) = Head(dvi)
         Head(dvi) = vi
         Last(vi) = -dvi
         if ( Next(vi)>0 ) Last(Next(vi)) = vi
         if ( dvi<Dmin ) Dmin = dvi
   !
      enddo MAIN
   endif

end subroutine mdu
