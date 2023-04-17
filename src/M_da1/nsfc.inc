!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!###  subroutine nsfc
!!
!!###  symbolic ldu-factorization of nonsymmetric sparse matrix
!!      (compressed pointer storage)
!!
!!
!!       input variables.. n, r, ic, ia, ja, jlmax, jumax.

!!       output variables.. il, jl, ijl, iu, ju, iju, flag.
!!
!!       parameters used internally..
!!```text
!! nia   - q     - suppose  m*  is the result of reordering  m.  if
!!       -           processing of the ith row of  m*  (hence the ith
!!       -           row of  u) is being done,  q(j)  is initially
!!       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
!!       -           values need not be stored, each entry points to the
!!       -           next nonzero and  q(n+1)  points to the first.  n+1
!!       -           indicates the end of the list.  for example, if n=9
!!       -           and the 5th row of  m*  is
!!       -              0 x x 0 x 0 0 x 0
!!       -           then  q  will initially be
!!       -              a a a a 8 a a 10 5           (a - arbitrary).
!!       -           as the algorithm proceeds, other elements of  q
!!       -           are inserted in the list because of fillin.
!!       -           q  is used in an analogous manner to compute the
!!       -           ith column of  l.
!!       -           size = n+1.
!! nia   - ira,  - vectors used to find the columns of  m.  at the kth
!! nia   - jra,      step of the factorization,  irac(k)  points to the
!! nia   - irac      head of a linked list in  jra  of row indices i
!!       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
!!       -           indicates the end of the list.  ira(i)  (i.ge.k)
!!       -           points to the smallest j such that j .ge. k and
!!       -           m(i,j)  is nonzero.
!!       -           size of each = n.
!! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
!! nia   - jrl       of the factorization,  jrl(k)  points to the head
!!       -           of a linked list in  jrl  of column indices j
!!       -           such j .lt. k and  l(k,j)  is nonzero.  zero
!!       -           indicates the end of the list.  irl(j)  (j.lt.k)
!!       -           points to the smallest i such that i .ge. k and
!!       -           l(i,j)  is nonzero.
!!       -           size of each = n.
!! nia   - iru,  - vectors used in a manner analogous to  irl and jrl
!! nia   - jru       to find the columns of  u.
!!       -           size of each = n.
!!```
!!
!!### internal variables..
!!
!!```text
!!    jlptr - points to the last position used in  jl.
!!    juptr - points to the last position used in  ju.
!!    jmin,jmax - are the indices in  a or u  of the first and last
!!                elements to be examined in a given row.
!!                for example,  jmin=ia(k), jmax=ia(k+1)-1.
!!```
!!
!-----------------------------------------------------------------------
subroutine nsfc(N,R,Ic,Ia,Ja,Jlmax,Il,Jl,Ijl,Jumax,Iu,Ju,Iju,Q,Ira,Jra,Irac,Irl,Jrl,Iru,Jru,Flag)
!
integer,intent(in)    :: N
integer,intent(inout) :: R(*)
integer,intent(in)    :: Ic(*)
integer,intent(in)    :: Ia(*)
integer,intent(in)    :: Ja(*)
integer,intent(in)    :: Jlmax
integer,intent(inout) :: Il(*)
integer,intent(inout) :: Jl(*)
integer,intent(inout) :: Ijl(*)
integer,intent(in)    :: Jumax
integer,intent(inout) :: Iu(*)
integer,intent(inout) :: Ju(*)
integer,intent(inout) :: Iju(*)
integer,intent(inout) :: Q(*)
integer,intent(inout) :: Ira(*)
integer,intent(inout) :: Jra(*)
integer,intent(inout) :: Irac(*)
integer,intent(inout) :: Irl(*)
integer,intent(inout) :: Jrl(*)
integer,intent(inout) :: Iru(*)
integer,intent(inout) :: Jru(*)
integer,intent(out)   :: Flag
!
integer :: cend, i, i1, iak, irai, irll, irul, j, jaiak, jairai, jlmin, jlptr, jmax, jmin, jtmp, jumin, juptr, k, &
         & lasti, lastid, long, luk, m, np1, qm, rend, rk, vj
   !
   !   ******  initialize pointers  ****************************************
   np1 = N + 1
   jlmin = 1
   jlptr = 0
   Il(1) = 1
   jumin = 1
   juptr = 0
   Iu(1) = 1
   do k = 1, N
      Irac(k) = 0
      Jra(k) = 0
      Jrl(k) = 0
      Jru(k) = 0
   enddo
   !   ******  initialize column pointers for a  ***************************
   do k = 1, N
      rk = R(k)
      iak = Ia(rk)
      if ( iak>=Ia(rk+1) ) then
   !
   !  ** error.. null row in a
         Flag = N + rk
         return
      else
         jaiak = Ic(Ja(iak))
         if ( jaiak>k ) then
            !  ** error.. null pivot
            Flag = 5*N + k
            return
         endif
         Jra(k) = Irac(jaiak)
         Irac(jaiak) = k
         Ira(k) = iak
      endif
   enddo
   !
   !   ******  for each column of l and row of u  **************************
   ALL: do k = 1, N
   !
   !   ******  initialize q for computing kth column of l  *****************
      Q(np1) = np1
      luk = -1
      !   ******  by filling in kth column of a  ******************************
      vj = Irac(k)
      if ( vj/=0 ) then
         qm = np1
         do
            m = qm
            qm = Q(m)
            if ( qm>=vj ) then
               if ( qm==vj )then
                  !  ** error.. duplicate entry in a
                  Flag = 2*N + rk
                  return
            endif
               luk = luk + 1
               Q(m) = vj
               Q(vj) = qm
               vj = Jra(vj)
               if ( vj==0 ) exit
               qm = np1
            endif
         enddo
      endif
      !   ******  link through jru  *******************************************
      lastid = 0
      lasti = 0
      Ijl(k) = jlptr
      i = k
      LOOP: do
         i = Jru(i)
         if ( i==0 ) then
         !   ******  lasti is the longest column merged into the kth  ************
         !   ******  see if it equals the entire kth column  *********************
            qm = Q(np1)
            if ( qm/=k ) then
               !  ** error.. null pivot
               Flag = 5*N + k
               return
         endif
            LUKZ: if ( luk/=0 ) then
               if ( lastid/=luk ) then
               !   ******  if not, see if kth column can overlap the previous one  *****
                  KTHCOL: if ( jlmin<=jlptr ) then
                     qm = Q(qm)
                     do j = jlmin, jlptr
                        if ( Jl(j)<qm ) cycle
                        if ( Jl(j)==qm ) then
                           Ijl(k) = j
                           do i = j, jlptr
                              if ( Jl(i)/=qm ) exit KTHCOL
                              qm = Q(qm)
                              if ( qm>N ) exit LUKZ
                           enddo
                           jlptr = j - 1
                        endif
                        exit
                     enddo
                  endif KTHCOL
                  !   ******  move column indices from q to jl, update vectors  ***********
                  jlmin = jlptr + 1
                  Ijl(k) = jlmin
                  if ( luk/=0 ) then
                     jlptr = jlptr + luk
                     if ( jlptr>Jlmax ) then
                        !  ** error.. insufficient storage for jl
                        Flag = 3*N + k
                        return
                     else
                        qm = Q(np1)
                        do j = jlmin, jlptr
                           qm = Q(qm)
                           Jl(j) = qm
                        enddo
                     endif
                  endif
               else
               !   ******  if so, jl can be compressed  ********************************
                  irll = Irl(lasti)
                  Ijl(k) = irll + 1
                  if ( Jl(irll)/=k ) Ijl(k) = Ijl(k) - 1
               endif
            endif LUKZ
            Irl(k) = Ijl(k)
            Il(k+1) = Il(k) + luk
            !
            !   ******  initialize q for computing kth row of u  ********************
            Q(np1) = np1
            luk = -1
            !   ******  by filling in kth row of reordered a  ***********************
            rk = R(k)
            jmin = Ira(k)
            jmax = Ia(rk+1) - 1
            if ( jmin<=jmax ) then
               do j = jmin, jmax
                  vj = Ic(Ja(j))
                  qm = np1
                  do
                     m = qm
                     qm = Q(m)
                     if ( qm>=vj ) then
                        if ( qm==vj )then
                           !  ** error.. duplicate entry in a
                           Flag = 2*N + rk
                           return
                     endif
                        luk = luk + 1
                        Q(m) = vj
                        Q(vj) = qm
                        exit
                     endif
                  enddo
               enddo
            endif
            !   ******  link through jrl,  ******************************************
            lastid = 0
            lasti = 0
            Iju(k) = juptr
            i = k
            i1 = Jrl(k)
            do
               i = i1
               if ( i==0 ) then
               !   ******  update jrl(k) and irl(k)  ***********************************
                  if ( Il(k+1)>Il(k) ) then
                     j = Jl(Irl(k))
                     Jrl(k) = Jrl(j)
                     Jrl(j) = k
                  endif
                  !   ******  lasti is the longest row merged into the kth  ***************
                  !   ******  see if it equals the entire kth row  ************************
                  qm = Q(np1)
                  if ( qm/=k )then
                     !  ** error.. null pivot
                     Flag = 5*N + k
                     return
               endif
                  LUKY: if ( luk/=0 ) then
                     if ( lastid/=luk ) then
                     !   ******  if not, see if kth row can overlap the previous one  ********
                        KTHROW: if ( jumin<=juptr ) then
                           qm = Q(qm)
                           do j = jumin, juptr
                              if ( Ju(j)<qm ) cycle
                              if ( Ju(j)==qm ) then
                                 Iju(k) = j
                                 do i = j, juptr
                                    if ( Ju(i)/=qm ) exit KTHROW
                                    qm = Q(qm)
                                    if ( qm>N ) exit LUKY
                                 enddo
                                 juptr = j - 1
                              endif
                              exit
                           enddo
                        endif KTHROW
                        !   ******  move row indices from q to ju, update vectors  **************
                        jumin = juptr + 1
                        Iju(k) = jumin
                        if ( luk/=0 ) then
                           juptr = juptr + luk
                           if ( juptr>Jumax ) then
                           !  ** error.. insufficient storage for ju
                              Flag = 6*N + k
                              return
                           else
                              qm = Q(np1)
                              do j = jumin, juptr
                                 qm = Q(qm)
                                 Ju(j) = qm
                              enddo
                           endif
                        endif
                     else
                        !   ******  if so, ju can be compressed  ********************************
                        irul = Iru(lasti)
                        Iju(k) = irul + 1
                        if ( Ju(irul)/=k ) Iju(k) = Iju(k) - 1
                     endif
                  endif LUKY
                  Iru(k) = Iju(k)
                  Iu(k+1) = Iu(k) + luk
               !
                  !   ******  update iru, jru  ********************************************
                  i = k
                  exit LOOP
               else
                  i1 = Jrl(i)
                  qm = np1
                  jmin = Iru(i)
                  jmax = Iju(i) + Iu(i+1) - Iu(i) - 1
                  long = jmax - jmin
                  if ( long>=0 ) then
                     jtmp = Ju(jmin)
                     if ( jtmp/=k ) then
                     !   ******  update irl and jrl, *****************************************
                        long = long + 1
                        cend = Ijl(i) + Il(i+1) - Il(i)
                        Irl(i) = Irl(i) + 1
                        if ( Irl(i)<cend ) then
                           j = Jl(Irl(i))
                           Jrl(i) = Jrl(j)
                           Jrl(j) = i
                        endif
                     endif
                     if ( lastid<long ) then
                        lasti = i
                        lastid = long
                     endif
                     !   ******  and merge the corresponding rows into the kth row  **********
                     do j = jmin, jmax
                        vj = Ju(j)
                        do
                           m = qm
                           qm = Q(m)
                           if ( qm>=vj ) then
                              if ( qm/=vj ) then
                                 luk = luk + 1
                                 Q(m) = vj
                                 Q(vj) = qm
                                 qm = vj
                              endif
                              exit
                           endif
                        enddo
                     enddo
                  endif
               endif
            enddo
         else
            qm = np1
            jmin = Irl(i)
            jmax = Ijl(i) + Il(i+1) - Il(i) - 1
            long = jmax - jmin
            if ( long>=0 ) then
               jtmp = Jl(jmin)
               if ( jtmp/=k ) long = long + 1
               if ( jtmp==k ) R(i) = -R(i)
               if ( lastid<long ) then
                  lasti = i
                  lastid = long
               endif
               !   ******  and merge the corresponding columns into the kth column  ****
               do j = jmin, jmax
                  vj = Jl(j)
                  do
                     m = qm
                     qm = Q(m)
                     if ( qm>=vj ) then
                        if ( qm/=vj ) then
                           luk = luk + 1
                           Q(m) = vj
                           Q(vj) = qm
                           qm = vj
                        endif
                        exit
                     endif
                  enddo
               enddo
            endif
         endif
      enddo LOOP
      BACK: do
         i1 = Jru(i)
         if ( R(i)<0 ) then
            R(i) = -R(i)
         else
            rend = Iju(i) + Iu(i+1) - Iu(i)
            if ( Iru(i)<rend ) then
               j = Ju(Iru(i))
               Jru(i) = Jru(j)
               Jru(j) = i
            endif
         endif
         i = i1
         if ( i==0 ) then
         !
         !   ******  update ira, jra, irac  **************************************
            i = Irac(k)
            if ( i==0 ) cycle ALL
         else
            Iru(i) = Iru(i) + 1
         cycle BACK
         endif
         exit BACK
      enddo BACK
      INFINITE: do
         i1 = Jra(i)
         Ira(i) = Ira(i) + 1
         if ( Ira(i)<Ia(R(i)+1) ) then
            irai = Ira(i)
            jairai = Ic(Ja(irai))
            if ( jairai<=i ) then
               Jra(i) = Irac(jairai)
               Irac(jairai) = i
            endif
         endif
         i = i1
         if ( i==0 ) exit INFINITE
      enddo INFINITE
   enddo ALL
   !
   Ijl(N) = jlptr
   Iju(N) = juptr
   Flag = 0
end subroutine nsfc
