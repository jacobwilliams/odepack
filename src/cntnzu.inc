!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine counts the number of nonzero elements in the strict
!! upper triangle of the matrix M + M(transpose), where the sparsity
!! structure of M is given by pointer arrays IA and JA.
!!
!! This is needed to compute the storage requirements for the
!! sparse matrix reordering operation in ODRV.

!-----------------------------------------------------------------------
subroutine cntnzu(N,Ia,Ja,Nzsut)
!
integer , intent(in)  :: N
integer , intent(in)  :: Ia(*)
integer , intent(in)  :: Ja(*)
integer , intent(out) :: Nzsut
!
integer :: ii , j , jj , jmax , jmin , k , kmax , kmin , num
!
num = 0
OUTER: do ii = 1 , N
   jmin = Ia(ii)
   jmax = Ia(ii+1) - 1
   if ( jmin<=jmax ) then
      MAIN: do j = jmin , jmax
         if ( Ja(j)<ii ) then
            jj = Ja(j)
            kmin = Ia(jj)
            kmax = Ia(jj+1) - 1
            if ( kmin<=kmax ) then
               do k = kmin , kmax
                  if ( Ja(k)==ii ) cycle MAIN
               enddo
            endif
         elseif ( Ja(j)==ii ) then
            cycle MAIN
         endif
         num = num + 1
      enddo MAIN
   endif
enddo OUTER
Nzsut = num
end subroutine cntnzu
