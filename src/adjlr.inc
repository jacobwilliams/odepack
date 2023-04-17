!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine computes an adjustment, LDIF, to the required
!! integer storage space in IWK (sparse matrix work space).
!!
!! It is called only if the word length ratio is LRAT = 1.
!! This is to account for the possibility that the symbolic LU phase
!! may require more storage than the numerical LU and solution phases.
!!
!-----------------------------------------------------------------------
subroutine adjlr(N,Isp,Ldif)

integer,intent(in)  :: N
integer,intent(in)  :: Isp(*)
integer,intent(out) :: Ldif

integer :: ip , jlmax , jumax , lnfc , lsfc , nzlu

   ip = 2*N + 1

   ! Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). -------
   !
   jlmax = Isp(ip)
   jumax = Isp(ip+ip)

   ! NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
   !
   nzlu = Isp(N+1) - Isp(1) + Isp(ip+N+1) - Isp(ip+1)
   lsfc = 12*N + 3 + 2*max(jlmax,jumax)
   lnfc = 9*N + 2 + jlmax + jumax + nzlu
   Ldif = max(0,lsfc-lnfc)

end subroutine adjlr
