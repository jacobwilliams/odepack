!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!  mdi -- initialization

!-----------------------------------------------------------------------
subroutine mdi(N,Ia,Ja,Max,V,L,Head,Last,Next,Mark,Tag,Flag)
!
integer,intent(in)    :: N
integer,intent(in)    :: Ia(*)
integer,intent(in)    :: Ja(*)
integer,intent(in)    :: Max
integer,intent(inout) :: V(*)
integer,intent(inout) :: L(*)
integer,intent(inout) :: Head(*)
integer,intent(out)   :: Last(*)
integer,intent(inout) :: Next(*)
integer,intent(inout) :: Mark(*)
integer,intent(in)    :: Tag
integer               :: Flag
!
integer :: dvi, j, jmax, jmin, k, kmax, lvk, nextvi, sfs, vi, vj
   !
   !----initialize degrees, element lists, and degree lists
   do vi = 1 , N
      Mark(vi) = 1
      L(vi) = 0
      Head(vi) = 0
   enddo
   sfs = N + 1
   !
   !----create nonzero structure
   !----for each nonzero entry a(vi,vj)
   do vi = 1 , N
      jmin = Ia(vi)
      jmax = Ia(vi+1) - 1
      if ( jmin<=jmax ) then
         MINMAX: do j = jmin , jmax
            vj = Ja(j)
            if ( vj<vi ) then
   !
   !------if a(vi,vj) is in strict lower triangle
   !------check for previous occurrence of a(vj,vi)
               lvk = vi
               kmax = Mark(vi) - 1
               if ( kmax/=0 ) then
                  do k = 1 , kmax
                     lvk = L(lvk)
                     if ( V(lvk)==vj ) cycle MINMAX
                  enddo
               endif
            elseif ( vj==vi ) then
               cycle
            endif
   !----for unentered entries a(vi,vj)
            if ( sfs>=Max ) then
   !
   !  ** error-  insufficient storage
               Flag = 9*N + vi
               return
            else
   !
   !------enter vj in element list for vi
               Mark(vi) = Mark(vi) + 1
               V(sfs) = vj
               L(sfs) = L(vi)
               L(vi) = sfs
               sfs = sfs + 1
   !
   !------enter vi in element list for vj
               Mark(vj) = Mark(vj) + 1
               V(sfs) = vi
               L(sfs) = L(vj)
               L(vj) = sfs
               sfs = sfs + 1
            endif
         enddo MINMAX
      endif
   enddo
   !
   !----create degree lists and initialize mark vector
   do vi = 1 , N
      dvi = Mark(vi)
      Next(vi) = Head(dvi)
      Head(dvi) = vi
      Last(vi) = -dvi
      nextvi = Next(vi)
      if ( nextvi>0 ) Last(nextvi) = vi
      Mark(vi) = Tag
   enddo

end subroutine mdi
