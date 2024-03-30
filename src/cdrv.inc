!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!! cdrv(3f) - [M_odepack ] driver for solving sparse non-symmetric
!!            systems of linear equations
!!
!!### DESCRIPTION
!!   driver for subroutines for solving sparse nonsymmetric systems of
!!   linear equations (compressed pointer storage)
!!
!!### PARAMETERS
!!
!!    class abbreviations are--
!!       n - integer variable
!!       f - real variable
!!       v - supplies a value to the driver
!!       r - returns a result from the driver
!!       i - used internally by the driver
!!       a - array
!!
!!### class - parameter
!!
!!   the nonzero entries of the coefficient matrix m are stored
!!   row-by-row in the array a.  to identify the individual nonzero
!!   entries in each row, we need to know in which column each entry
!!   lies.  the column indices which correspond to the nonzero entries
!!   of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
!!   ja(k) = j.  in addition, we need to know where each row starts and
!!   how long it is.  the index positions in ja and a where the rows of
!!   m begin are stored in the array ia.  i.e., if m(i,j) is the first
!!   nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
!!   ia(i) = k.  moreover, the index in ja and a of the first location
!!   following the last element in the last row is stored in ia(n+1).
!!   thus, the number of entries in the i-th row is given by
!!   ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!!   consecutively in
!!```text
!!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!!```
!!   and the corresponding column indices are stored consecutively in
!!```text
!!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!!```
!!   for example, the 5 by 5 matrix
!!```text
!!                ( 1. 0. 2. 0. 0.)
!!                ( 0. 3. 0. 0. 0.)
!!            m = ( 0. 4. 5. 6. 0.)
!!                ( 0. 0. 0. 7. 0.)
!!                ( 0. 0. 0. 8. 9.)
!!```
!!   would be stored as
!!```text
!!               - 1  2  3  4  5  6  7  8  9
!!            ---+--------------------------
!!            ia - 1  3  4  7  8 10
!!            ja - 1  3  2  2  3  4  4  4  5
!!             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!!```
!!
!!```text
!!     nv    - n     - number of variables/equations.
!!     fva   - a     - nonzero entries of the coefficient matrix m, stored
!!           -           by rows.
!!           -           size = number of nonzero entries in m.
!!     nva   - ia    - pointers to delimit the rows in a.
!!           -           size = n+1.
!!     nva   - ja    - column numbers corresponding to the elements of a.
!!           -           size = size of a.
!!     fva   - b     - right-hand side b.  b and z can the same array.
!!           -           size = n.
!!     fra   - z     - solution x.  b and z can be the same array.
!!           -           size = n.
!!```
!!
!!   the rows and columns of the original matrix m can be
!!   reordered (e.g., to reduce fillin or ensure numerical stability)
!!   before calling the driver.  if no reordering is done, then set
!!   r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
!!   in the original order.
!!
!!   if the columns have been reordered (i.e.,  c(i).ne.i  for some
!!   i), then the driver will call a subroutine (nroc) which rearranges
!!   each row of ja and a, leaving the rows in the original order, but
!!   placing the elements of each row in increasing order with respect
!!   to the new ordering.  if  path.ne.1,  then nroc is assumed to have
!!   been called already.
!!
!!```text
!! nva   - r     - ordering of the rows of m.
!!       -           size = n.
!! nva   - c     - ordering of the columns of m.
!!       -           size = n.
!! nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
!!       -           ic(c(i)) = i  for i=1,...,n.
!!       -           size = n.
!!```
!!
!!   the solution of the system of linear equations is divided into
!!   three stages --
!!
!!```text
!!      nsfc -- the matrix m is processed symbolically to determine where
!!               fillin will occur during the numeric factorization.
!!      nnfc -- the matrix m is factored numerically into the product ldu
!!               of a unit lower triangular matrix l, a diagonal matrix
!!               d, and a unit upper triangular matrix u, and the system
!!               mx = b  is solved.
!!      nnsc -- the linear system  mx = b  is solved using the ldu
!!  or           factorization from nnfc.
!!      nntc -- the transposed linear system  mt x = b  is solved using
!!               the ldu factorization from nnf.
!!```
!!
!!   for several systems whose coefficient matrices have the same
!!   nonzero structure, nsfc need be done only once (for the first
!!   system).  then nnfc is done once for each additional system.  for
!!   several systems with the same coefficient matrix, nsfc and nnfc
!!   need be done only once (for the first system).  then nnsc or nntc
!!   is done once for each additional right-hand side.
!!```text
!! nv    - path  - path specification.  values and their meanings are --
!!       -           1  perform nroc, nsfc, and nnfc.
!!       -           2  perform nnfc only  (nsfc is assumed to have been
!!       -               done in a manner compatible with the storage
!!       -               allocation used in the driver).
!!       -           3  perform nnsc only  (nsfc and nnfc are assumed to
!!       -               have been done in a manner compatible with the
!!       -               storage allocation used in the driver).
!!       -           4  perform nntc only  (nsfc and nnfc are assumed to
!!       -               have been done in a manner compatible with the
!!       -               storage allocation used in the driver).
!!       -           5  perform nroc and nsfc.
!!```
!!
!!   various errors are detected by the driver and the individual
!!   subroutines.
!!
!!```text
!! nr    - flag  - error flag.  values and their meanings are --
!!       -             0     no errors detected
!!       -             n+k   null row in a  --  row = k
!!       -            2n+k   duplicate entry in a  --  row = k
!!       -            3n+k   insufficient storage in nsfc  --  row = k
!!       -            4n+1   insufficient storage in nnfc
!!       -            5n+k   null pivot  --  row = k
!!       -            6n+k   insufficient storage in nsfc  --  row = k
!!       -            7n+1   insufficient storage in nnfc
!!       -            8n+k   zero pivot  --  row = k
!!       -           10n+1   insufficient storage in cdrv
!!       -           11n+1   illegal path specification
!!```
!!
!!   working storage is needed for the factored form of the matrix
!!   m plus various temporary vectors.  the arrays isp and rsp should be
!!   equivalenced.  integer storage is allocated from the beginning of
!!   isp and real storage from the end of rsp.
!!
!!```text
!!        nv    - nsp   - declared dimension of rsp.  nsp generally must
!!                        be larger than  8n+2 + 2k  (where  k = (number of
!!                        nonzero entries in m)).
!!        nvira - isp   - integer working storage divided up into various arrays
!!                        needed by the subroutines.  isp and rsp should be
!!                        equivalenced.
!!                             size = lratio*nsp.
!!        fvira - rsp   - real working storage divided up into various arrays
!!                        needed by the subroutines.  isp and rsp should be
!!                        equivalenced.
!!                             size = nsp.
!!        nr    - esp   - if sufficient storage was available to perform the
!!                        symbolic factorization (nsfc), then esp is set to
!!                        the amount of excess storage provided (negative if
!!                        insufficient storage was available to perform the
!!                        numeric factorization (nnfc)).
!!```
!!
!-----------------------------------------------------------------------
subroutine cdrv(N,R,C,Ic,Ia,Ja,A,B,Z,Nsp,Isp,Rsp,Esp,Path,Flag)
!
integer               :: N
integer               :: R(*)
integer               :: C(*)
integer               :: Ic(*)
integer               :: Ia(*)
integer               :: Ja(*)
real(kind=dp)         :: A(*)
real(kind=dp)         :: B(*)
real(kind=dp)         :: Z(*)
integer,intent(inout) :: Isp(*)
integer,intent(in)    :: Nsp
real(kind=dp)         :: Rsp(*)
integer,intent(out)   :: Esp
integer,intent(in)    :: Path
integer,intent(inout) :: Flag
!
integer :: ar, d, i, ijl, iju, il, ira, irac, irl, iru, iu, j, jl, jlmax, jra, jrl, jru, ju, jumax, jutmp, l,  &
         & lmax, max, q, row, tmp, u, umax
integer , save :: lratio
!
!   set lratio equal to the ratio between the length of floating point
!   and integer array data.  e. g., lratio = 1 for (real, integer),
!   lratio = 2 for (double precision, integer)
!
data lratio/2/
!
if ( Path<1 .or. 5<Path ) then
!  ** error.. illegal path specification
   Flag = 11*N + 1
   return
else
! ### ***initialize and divide up temporary storage  *******************
   il = 1
   ijl = il + (N+1)
   iu = ijl + N
   iju = iu + (N+1)
   irl = iju + N
   jrl = irl + N
   jl = jrl + N
!
!   ******  reorder a if necessary, call nsfc if flag is set  ***********
   if ( (Path-1)*(Path-5)==0 ) then
      max = (lratio*Nsp+1-jl) - (N+1) - 5*N
      jlmax = max/2
      q = jl + jlmax
      ira = q + (N+1)
      jra = ira + N
      irac = jra + N
      iru = irac + N
      jru = iru + N
      jutmp = jru + N
      jumax = lratio*Nsp + 1 - jutmp
      Esp = max/lratio
      if ( jlmax<=0 .or. jumax<=0 )then
         !  ** error.. insufficient storage
         Flag = 10*N + 1
         return
      endif
!
      do i = 1 , N
         if ( C(i)/=i ) then
            ar = Nsp + 1 - N
            call nroc(N,Ic,Ia,Ja,A,Isp(il),Rsp(ar),Isp(iu),Flag)
            if ( Flag==0 ) exit
               !  ** error.. error detected in nroc, nsfc, nnfc, or nnsc
               return
         endif
      enddo
!
      call nsfc(N,R,Ic,Ia,Ja,jlmax,Isp(il),Isp(jl),Isp(ijl),jumax,Isp(iu),Isp(jutmp),Isp(iju),Isp(q),Isp(ira),Isp(jra),Isp(irac), &
       & Isp(irl),Isp(jrl),Isp(iru),Isp(jru),Flag)
      if ( Flag/=0 )then
         !  ** error.. error detected in nroc, nsfc, nnfc, or nnsc
         return
      endif
!   ******  move ju next to jl  *****************************************
      jlmax = Isp(ijl+N-1)
      ju = jl + jlmax
      jumax = Isp(iju+N-1)
      if ( jumax>0 ) then
         do j = 1 , jumax
            Isp(ju+j-1) = Isp(jutmp+j-1)
         enddo
      endif
   endif
!
!   ******  call remaining subroutines  *********************************
   jlmax = Isp(ijl+N-1)
   ju = jl + jlmax
   jumax = Isp(iju+N-1)
   l = (ju+jumax-2+lratio)/lratio + 1
   lmax = Isp(il+N) - 1
   d = l + lmax
   u = d + N
   row = Nsp + 1 - N
   tmp = row - N
   umax = tmp - u
   Esp = umax - (Isp(iu+N)-1)
!
   if ( (Path-1)*(Path-2)==0 ) then
      if ( umax<0 ) then
         !  ** error.. insufficient storage
         Flag = 10*N + 1
         return
      endif
      call nnfc(N,R,C,Ic,Ia,Ja,A,Z,B,lmax,Isp(il),Isp(jl),Isp(ijl),Rsp(l),Rsp(d),umax,Isp(iu),Isp(ju),Isp(iju),Rsp(u),Rsp(row),&
              & Rsp(tmp),Isp(irl),Isp(jrl),Flag)
      if ( Flag/=0 ) then
         !  ** error.. error detected in nroc, nsfc, nnfc, or nnsc
         return
      endif
   endif
!
   if ( (Path-3)==0 ) call nnsc(N,R,C,Isp(il),Isp(jl),Isp(ijl),Rsp(l),Rsp(d),Isp(iu),Isp(ju),Isp(iju),Rsp(u),Z,B,Rsp(tmp))
!
   if ( (Path-4)==0 ) call nntc(N,R,C,Isp(il),Isp(jl),Isp(ijl),Rsp(l),Rsp(d),Isp(iu),Isp(ju),Isp(iju),Rsp(u),Z,B,Rsp(tmp))
   return
endif
!
end subroutine cdrv
