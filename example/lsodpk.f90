!>
!! Demonstration program for DLSODPK.
!! ODE system from ns-species interaction pde in 2 dimensions.
!! This is the version of 14 June 2001.
!!
!! This version is in double precision.
!-----------------------------------------------------------------------
!! This program solves a stiff ODE system that arises from a system
!! of partial differential equations.  The PDE system is a food web
!! population model, with predator-prey interaction and diffusion on
!! the unit square in two dimensions.  The dependent variable vector is
!!```text
!!         1   2        ns
!!   c = (c , c , ..., c  )
!!```
!! and the PDEs are as follows:
!!```text
!!     i               i      i
!!   dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
!!                     xx     yy        i
!!
!! where
!!                  i          ns         j
!!   f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
!!    i                       j=1
!!```
!! The number of species is ns = 2*np, with the first np being prey and
!! the last np being predators.  The coefficients a(i,j), b(i), d(i) are:
!!```text
!!   a(i,i) = -a  (all i)
!!   a(i,j) = -g  (i .le. np, j .gt. np)
!!   a(i,j) =  e  (i .gt. np, j .le. np)
!!   b(i) =  b*(1 + alpha*x*y)  (i .le. np)
!!   b(i) = -b*(1 + alpha*x*y)  (i .gt. np)
!!   d(i) = dprey  (i .le. np)
!!   d(i) = dpred  (i .gt. np)
!!```
!! The various scalar parameters are set in subroutine setpar.
!!
!! The boundary conditions are: normal derivative = 0.
!! A polynomial in x and y is used to set the initial conditions.
!!
!! The PDEs are discretized by central differencing on a mx by my mesh.
!!
!! The ODE system is solved by DLSODPK using method flag values
!! mf = 10, 21, 22, 23, 24, 29.  The final time is tmax = 10, except
!! that for mf = 10 it is tmax = 1.0d-3 because the problem is stiff,
!! and for mf = 23 and 24 it is tmax = 2 because the lack of symmetry
!! in the problem makes these methods more costly.
!!
!! Two preconditioner matrices are used.  One uses a fixed number of
!! Gauss-Seidel iterations based on the diffusion terms only.
!! The other preconditioner is a block-diagonal matrix based on
!! the partial derivatives of the interaction terms f only, using
!! block-grouping (computing only a subset of the ns by ns blocks).
!! For mf = 21 and 22, these two preconditioners are applied on
!! the left and right, respectively, and for mf = 23 and 24 the product
!! of the two is used as the one preconditioner matrix.
!! For mf = 29, the inverse of the product is applied.
!!
!! Two output files are written: one with the problem description and
!! and performance statistics on unit 6, and one with solution profiles
!! at selected output times (for mf = 22 only) on unit 8.
!!-----------------------------------------------------------------------
!! Note: In addition to the main program and 10 subroutines
!! given below, this program requires the LINPACK subroutines
!! DGEFA and DGESL, and the BLAS routine DAXPY.
!!-----------------------------------------------------------------------
!!### Reference:
!!     Peter N. Brown and Alan C. Hindmarsh,
!!     Reduced Storage Matrix Methods in Stiff ODE Systems,
!!     J. Appl. Math. & Comp., 31 (1989), pp. 40-91;
!!     Also LLNL Report UCRL-95088, Rev. 1, June 1987.
!-----------------------------------------------------------------------
program lsodpk
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: aa , alph , ax , ay , bb , dpred , dprey , dx ,   &
                    & dy , ee , gg , srur , uround
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: itmax , meshx , meshy , mp , mpsq , mq , mx , mxmp ,   &
 & mxns , my , ngrp , ngx , ngy , ns
integer , dimension(21) :: jgx , jgy
integer , dimension(50) :: jigx , jigy
integer , dimension(20) :: jxr , jyr
common /pcom0 / aa , ee , gg , bb , dprey , dpred
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
 & cox , coy , ns , mx , my , mxns
common /pcom2 / uround , srur , mp , mq , mpsq , itmax
common /pcom3 / meshx , meshy , ngx , ngy , ngrp , mxmp , jgx ,   &
 & jgy , jigx , jigy , jxr , jyr
!
real(kind=dp) :: atol , avdim , hu , rcfl , rcfn , rtol , t , tout
real(kind=dp) , dimension(576) :: cc
integer :: i , imod3 , iopt , iout , istate , itask , itol ,      &
 & jacflg , jpre , leniw , lenrw , mf , ncfl , ncfn ,     &
 & neq , nfe , nfldif , nfndif , nli , nlidif , nni ,     &
 & nnidif , nout , npe , nps , nqu , nsdif , nst
integer , dimension(67) :: iwork
integer , save :: liw , lrw
real(kind=dp) , dimension(5213) :: rwork
external fweb , jacbg , solsbg
!
! The problem Common blocks below allow for up to 20 species,
! up to a 50x50 mesh, and up to a 20x20 group structure.
!
! The dimension of cc below must be .ge. 2*neq, where neq = ns*mx*my.
! The dimension lrw of rwork must be .ge. 17*neq + ns*ns*ngrp + 61,
! and the dimension liw of iwork must be .ge. 35 + ns*ngrp.
data lrw/5213/ , liw/67/
!
      !open (unit=6, file='lsodpk.out', status='new')
      open (unit=8,file='ccout',status='new')
!
      ax = 1.0d0
      ay = 1.0d0
!
! Call setpar to set problem parameters.
      call setpar
!
! Set remaining problem parameters.
      neq = ns*mx*my
      mxns = mx*ns
      dx = ax/(mx-1)
      dy = ay/(my-1)
      do i = 1 , ns
         cox(i) = diff(i)/dx**2
         coy(i) = diff(i)/dy**2
      enddo
!
! Write heading.
      write (6,99001) ns , mx , my , neq
99001 format (' Demonstration program for DLSODPK package'//            &
             &' Food web problem with ns species, ns =',i4/             &
             &' Predator-prey interaction and diffusion on a 2-d square'&
            & //' Mesh dimensions (mx,my) =',                           &
             &2i4/' Total system size is neq =',i7//)
      write (6,99002) aa , ee , gg , bb , dprey , dpred , alph
99002 format (' Matrix parameters:  a =',d12.4,'   e =',d12.4,'   g =', &
            & d12.4/20x,' b =',d12.4//' Diffusion coefficients: dprey ='&
            & ,d12.4,'   dpred =',d12.4/' Rate parameter alpha =',      &
            & d12.4//)
!
! Set remaining method parameters.
      jpre = 3
      jacflg = 1
      iwork(3) = jpre
      iwork(4) = jacflg
      iopt = 0
      mp = ns
      mq = mx*my
      mpsq = ns*ns
      uround = dumach()
      srur = sqrt(uround)
      meshx = mx
      meshy = my
      mxmp = meshx*mp
      ngx = 2
      ngy = 2
      ngrp = ngx*ngy
      call gset(meshx,ngx,jgx,jigx,jxr)
      call gset(meshy,ngy,jgy,jigy,jyr)
      iwork(1) = mpsq*ngrp
      iwork(2) = mp*ngrp
      itmax = 5
      itol = 1
      rtol = 1.0d-5
      atol = rtol
      itask = 1
      write (6,99003) ngrp , ngx , ngy , itmax , rtol , atol
99003 format (' Preconditioning uses interaction-only block-diagonal',  &
             &' matrix'/                                                &
             &' with block-grouping, and Gauss-Seidel iterations'//     &
             &' Number of diagonal block groups = ngrp =',i4,           &
             &'   (ngx by ngy, ngx =',i2,'  ngy =',i2,                  &
             &' )'//' G-S preconditioner uses itmax iterations, itmax ='&
            & ,i3//' Tolerance parameters: rtol =',d10.2,'   atol =',   &
            & d10.2)
!
!
! Loop over mf values 10, 21, ..., 24, 29.
!
      do mf = 10 , 29
         if ( mf<=10 .or. mf>=21 ) then
            if ( mf<=24 .or. mf>=29 ) then
               write (6,99004) mf
99004          format (//80('-')//' Solution with mf =',                &
                      &i3//'   t       nstep  nfe  nni  nli  npe  nq',  &
                     & 4x,'h          avdim    ncf rate    lcf rate')
!
               t = 0.0d0
               tout = 1.0d-8
               nout = 18
               if ( mf==10 ) nout = 6
               if ( mf==23 .or. mf==24 ) nout = 10
               call cinit(cc)
               if ( mf==22 ) call outweb(t,cc,ns,mx,my,8)
               istate = 1
               nli = 0
               nni = 0
               ncfn = 0
               ncfl = 0
               nst = 0
!
! Loop over output times and call DLSODPK.
!
               do iout = 1 , nout
                  call dlsodpk(fweb,[neq],cc,t,tout,itol,[rtol],[atol],itask, &
                             & istate,iopt,rwork,lrw,iwork,liw,jacbg,   &
                             & solsbg,mf)
                  nsdif = iwork(11) - nst
                  nst = iwork(11)
                  nfe = iwork(12)
                  npe = iwork(13)
                  nnidif = iwork(19) - nni
                  nni = iwork(19)
                  nlidif = iwork(20) - nli
                  nli = iwork(20)
                  nfndif = iwork(22) - ncfn
                  ncfn = iwork(22)
                  nfldif = iwork(23) - ncfl
                  ncfl = iwork(23)
                  nqu = iwork(14)
                  hu = rwork(11)
                  avdim = 0.0d0
                  rcfn = 0.0d0
                  rcfl = 0.0d0
                  if ( nnidif>0 ) avdim = real(nlidif)/real(nnidif)
                  if ( nsdif>0 ) rcfn = real(nfndif)/real(nsdif)
                  if ( nnidif>0 ) rcfl = real(nfldif)/real(nnidif)
                  write (6,99005) t , nst , nfe , nni , nli , npe ,     &
                                & nqu , hu , avdim , rcfn , rcfl
99005             format (d10.2,i5,i6,3i5,i4,2d11.2,d10.2,d12.2)
                  imod3 = iout - 3*(iout/3)
                  if ( mf==22 .and. imod3==0 )                          &
                     & call outweb(t,cc,ns,mx,my,8)
                  if ( istate==2 ) then
                     if ( tout>0.9d0 ) tout = tout + 1.0d0
                     if ( tout<0.9d0 ) tout = tout*10.0d0
                  else
                     write (6,99006) t
99006                format (//' final time reached =',d12.4//)
                     exit
                  endif
               enddo
!
               nst = iwork(11)
               nfe = iwork(12)
               npe = iwork(13)
               lenrw = iwork(17)
               leniw = iwork(18)
               nni = iwork(19)
               nli = iwork(20)
               nps = iwork(21)
               if ( nni>0 ) avdim = real(nli)/real(nni)
               ncfn = iwork(22)
               ncfl = iwork(23)
               write (6,99007) lenrw , leniw , nst , nfe , npe , nps ,  &
                             & nni , nli , avdim , ncfn , ncfl
99007          format (//' Final statistics for this run:'/             &
                      &' rwork size =',i8,'   iwork size =',            &
                      &i6/' number of time steps            =',         &
                      &i5/' number of f evaluations         =',         &
                      &i5/' number of preconditioner evals. =',         &
                      &i5/' number of preconditioner solves =',         &
                      &i5/' number of nonlinear iterations  =',         &
                      &i5/' number of linear iterations     =',         &
                      &i5/' average subspace dimension  =',f8.4/i5,     &
                      &' nonlinear conv. failures,',i5,                 &
                      &' linear conv. failures')
            endif
         endif
!
      enddo
!------  end of main program for DLSODPK demonstration program ----------
end program lsodpk

subroutine setpar
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: aa , alph , ax , ay , bb , dpred , dprey , dx ,   &
 & dy , ee , gg
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: mx , mxns , my , ns
common /pcom0 / aa , ee , gg , bb , dprey , dpred
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
                    & cox , coy , ns , mx , my , mxns
!
integer :: i , j , np
!
      np = 3
      mx = 6
      my = 6
      aa = 1.0d0
      ee = 1.0d4
      gg = 0.5d-6
      bb = 1.0d0
      dprey = 1.0d0
      dpred = 0.5d0
      alph = 1.0d0
      ns = 2*np
      do j = 1 , np
         do i = 1 , np
            acoef(np+i,j) = ee
            acoef(i,np+j) = -gg
         enddo
         acoef(j,j) = -aa
         acoef(np+j,np+j) = -aa
         bcoef(j) = bb
         bcoef(np+j) = -bb
         diff(j) = dprey
         diff(np+j) = dpred
      enddo
!
end subroutine setpar

subroutine gset(m,ng,jg,jig,jr)
implicit none
!
integer :: m
integer :: ng
integer , dimension(*) :: jg
integer , dimension(*) :: jig
integer , dimension(*) :: jr
!
integer :: ig , j , len1 , mper , ngm1
!
      mper = m/ng
      do ig = 1 , ng
         jg(ig) = 1 + (ig-1)*mper
      enddo
      jg(ng+1) = m + 1
!
      ngm1 = ng - 1
      len1 = ngm1*mper
      do j = 1 , len1
         jig(j) = 1 + (j-1)/mper
      enddo

      len1 = len1 + 1
      do j = len1 , m
         jig(j) = ng
      enddo
!
      do ig = 1 , ngm1
         jr(ig) = 0.5d0 + (ig-0.5d0)*mper
      enddo
      jr(ng) = 0.5d0*(1+ngm1*mper+m)
!
end subroutine gset

subroutine cinit(cc)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) :: alph , ax , ay , dx , dy
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: mx , mxns , my , ns
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
 & cox , coy , ns , mx , my , mxns
!
real(kind=dp) , dimension(*) :: cc
!
real(kind=dp) :: argx , argy , x , y
integer :: i , ici , ioff , iyoff , jx , jy
!
      do jy = 1 , my
         y = (jy-1)*dy
         argy = 16.0d0*y*y*(ay-y)*(ay-y)
         iyoff = mxns*(jy-1)
         do jx = 1 , mx
            x = (jx-1)*dx
            argx = 16.0d0*x*x*(ax-x)*(ax-x)
            ioff = iyoff + ns*(jx-1)
            do i = 1 , ns
               ici = ioff + i
               cc(ici) = 10.0d0 + i*argx*argy
            enddo
         enddo
      enddo
end subroutine cinit

subroutine outweb(t,c,ns,mx,my,lun)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: t
integer :: ns
integer :: mx
integer :: my
real(kind=dp) , dimension(ns,mx,my) :: c
integer :: lun
!
integer :: i , jx , jy
!
      write (lun,99001) t
99001 format (/80('-')/30x,'At time t = ',d16.8/80('-'))
!
      do i = 1 , ns
         write (lun,99002) i
99002    format (' the species c(',i2,') values are:')
         do jy = my , 1 , -1
            write (lun,99003) (c(i,jx,jy),jx=1,mx)
99003       format (6(1x,g12.6))
         enddo
         write (lun,99004)
99004    format (80('-'),/)
      enddo
!
end subroutine outweb

subroutine fweb(neq,t,cc,cdot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) :: alph , ax , ay , dx , dy
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: mx , mxns , my , ns
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
 & cox , coy , ns , mx , my , mxns
!
integer       :: neq
real(kind=dp) :: t
real(kind=dp) :: cc(2*neq)
real(kind=dp) :: cdot(neq)
!
real(kind=dp) :: dcxli , dcxui , dcyli , dcyui , x , y
integer :: i , ic , ici , idxl , idxu , idyl , idyu , iyoff , jx , jy
!
      do jy = 1 , my
         y = (jy-1)*dy
         iyoff = mxns*(jy-1)
         idyu = mxns
         if ( jy==my ) idyu = -mxns
         idyl = mxns
         if ( jy==1 ) idyl = -mxns
         do jx = 1 , mx
            x = (jx-1)*dx
            ic = iyoff + ns*(jx-1) + 1
! Get interaction rates at one point (x,y).
            !write(*,*)'GOT HERE A:',size(cc),ic,neq,neq+ic
            call webr(x,y,t,cc(ic),cc(neq+ic))
            idxu = ns
            if ( jx==mx ) idxu = -ns
            idxl = ns
            if ( jx==1 ) idxl = -ns
            do i = 1 , ns
               ici = ic + i - 1
! Do differencing in y.
               dcyli = cc(ici) - cc(ici-idyl)
               dcyui = cc(ici+idyu) - cc(ici)
! Do differencing in x.
               dcxli = cc(ici) - cc(ici-idxl)
               dcxui = cc(ici+idxu) - cc(ici)
! Collect terms and load cdot elements.
               cdot(ici) = coy(i)*(dcyui-dcyli) + cox(i)*(dcxui-dcxli)  &
                         & + cc(neq+ici)
            enddo
         enddo
      enddo
end subroutine fweb

subroutine webr(x,y,t,c,rate)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: x
real(kind=dp) :: y
real(kind=dp) :: t
real(kind=dp) :: c(*)
real(kind=dp) :: rate(*)
!
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) :: alph , ax , ay , dx , dy
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: mx , mxns , my , ns
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
 & cox , coy , ns , mx , my , mxns
!
!
real(kind=dp) :: fac
integer :: i
!
      do i = 1 , ns
         rate(i) = 0.0d0
      enddo
      do i = 1 , ns
         call daxpy(ns,c(i),acoef(1,i),1,rate,1)
      enddo
      fac = 1.0d0 + alph*x*y
      do i = 1 , ns
         rate(i) = c(i)*(bcoef(i)*fac+rate(i))
      enddo
end subroutine webr

subroutine jacbg(f,neq,t,cc,ccsv,rewt,f0,f1,hl0,bd,ipbd,ier)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: itmax , meshx , meshy , mp , mpsq , mq , mxmp , ngrp , &
 & ngx , ngy
integer , dimension(21) :: jgx , jgy
integer , dimension(50) :: jigx , jigy
integer , dimension(20) :: jxr , jyr
real(kind=dp) :: srur , uround
common /pcom2 / uround , srur , mp , mq , mpsq , itmax
common /pcom3 / meshx , meshy , ngx , ngy , ngrp , mxmp , jgx ,   &
 & jgy , jigx , jigy , jxr , jyr
!
real(kind=dp) , dimension(*) :: bd
integer :: neq
real(kind=dp) , dimension(neq) :: ccsv , f0 , f1 , rewt
real(kind=dp) , dimension(2*neq) :: cc 
real(kind=dp) :: fac , r , r0
real(kind=dp) :: hl0 , t
integer :: i , ibd , idiag , if0 , if00 , ig , igx , igy , iip ,  &
 & j , jj , jx , jy , n
integer :: ier 
integer , dimension(*) :: ipbd
external f
!
      n = neq
!-----------------------------------------------------------------------
! Make mp calls to fbg to approximate each diagonal block of Jacobian.
! Here cc(neq+1),...,cc(2*neq) contains the base fb value.
! r0 is a minimum increment factor for the difference quotient.
!-----------------------------------------------------------------------
      fac = dvnorm(n,f0,rewt)
      r0 = 1000.0d0*abs(hl0)*uround*n*fac
      if ( r0==0.0d0 ) r0 = 1.0d0
      ibd = 0
      do igy = 1 , ngy
         jy = jyr(igy)
         if00 = (jy-1)*mxmp
         do igx = 1 , ngx
            jx = jxr(igx)
            if0 = if00 + (jx-1)*mp
            do j = 1 , mp
               jj = if0 + j
               r = max(srur*abs(cc(jj)),r0/rewt(jj))
               cc(jj) = cc(jj) + r
               fac = -hl0/r
               call fbg(neq,t,cc,jx,jy,f1)
               do i = 1 , mp
                  bd(ibd+i) = (f1(i)-cc(neq+if0+i))*fac
               enddo
               cc(jj) = ccsv(jj)
               ibd = ibd + mp
            enddo
         enddo
      enddo
!
! Add identity matrix and do LU decompositions on blocks. --------------
      ibd = 1
      iip = 1
      do ig = 1 , ngrp
         idiag = ibd
         do i = 1 , mp
            bd(idiag) = bd(idiag) + 1.0d0
            idiag = idiag + (mp+1)
         enddo
         call dgefa(bd(ibd),mp,mp,ipbd(iip),ier)
         if ( ier/=0 ) exit
         ibd = ibd + mpsq
         iip = iip + mp
      enddo
end subroutine jacbg

subroutine fbg(neq,t,cc,jx,jy,cdot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) :: alph , ax , ay , dx , dy
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: mx , mxns , my , ns
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
 & cox , coy , ns , mx , my , mxns
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: cc
integer :: jx
integer :: jy
real(kind=dp) , dimension(neq) :: cdot
!
integer :: iblok , ic
real(kind=dp) :: x , y
!
      iblok = jx + (jy-1)*mx
      y = (jy-1)*dy
      x = (jx-1)*dx
      ic = ns*(iblok-1) + 1
      call webr(x,y,t,cc(ic),cdot)
end subroutine fbg

subroutine solsbg(n,t,cc,f0,wk,hl0,bd,ipbd,v,lr,ier)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: itmax , meshx , meshy , mp , mpsq , mq , mxmp , ngrp , &
 & ngx , ngy
integer , dimension(21) :: jgx , jgy
integer , dimension(50) :: jigx , jigy
integer , dimension(20) :: jxr , jyr
real(kind=dp) :: srur , uround
common /pcom2 / uround , srur , mp , mq , mpsq , itmax
common /pcom3 / meshx , meshy , ngx , ngy , ngrp , mxmp , jgx ,   &
 & jgy , jigx , jigy , jxr , jyr
!
integer :: n
real(kind=dp) :: t
real(kind=dp) , dimension(n) :: cc
real(kind=dp) , dimension(n) :: f0
real(kind=dp) , dimension(n) :: wk
real(kind=dp) :: hl0
real(kind=dp) , dimension(*) :: bd
integer , dimension(*) :: ipbd
real(kind=dp) , dimension(n) :: v
integer :: lr
integer :: ier
!
integer :: ibd , ig0 , igm1 , igx , igy , iip , iv , jx , jy
!
      ier = 0
!
      if ( lr==0 .or. lr==1 .or. lr==3 ) call gs(n,hl0,v,wk)
      if ( lr==0 .or. lr==2 .or. lr==3 ) then
         iv = 1
         do jy = 1 , meshy
            igy = jigy(jy)
            ig0 = (igy-1)*ngx
            do jx = 1 , meshx
               igx = jigx(jx)
               igm1 = igx - 1 + ig0
               ibd = 1 + igm1*mpsq
               iip = 1 + igm1*mp
               call dgesl(bd(ibd),mp,mp,ipbd(iip),v(iv),0)
               iv = iv + mp
            enddo
         enddo
      endif
!
end subroutine solsbg

subroutine gs(n,hl0,z,x)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) , dimension(20,20) :: acoef
real(kind=dp) :: alph , ax , ay , dx , dy , srur , uround
real(kind=dp) , dimension(20) :: bcoef , cox , coy , diff
integer :: itmax , mp , mpsq , mq , mx , mxns , my , ns
common /pcom1 / ax , ay , acoef , bcoef , dx , dy , alph , diff , &
 & cox , coy , ns , mx , my , mxns
common /pcom2 / uround , srur , mp , mq , mpsq , itmax
!
integer :: n
real(kind=dp) :: hl0
real(kind=dp) , dimension(n) :: z
real(kind=dp) , dimension(n) :: x
!
real(kind=dp) , dimension(20) :: beta , beta2 , cof1 , gamma ,     &
                                    & gamma2
real(kind=dp) :: elamda
integer :: i , ic , ici , iter , iyoff , jx , jy
!
!-----------------------------------------------------------------------
! Write matrix as P = D - L - U.
! Load local arrays beta, beta2, gamma, gamma2, and cof1.
!-----------------------------------------------------------------------
      do i = 1 , ns
         elamda = 1.d0/(1.d0+2.d0*hl0*(cox(i)+coy(i)))
         beta(i) = hl0*cox(i)*elamda
         beta2(i) = 2.d0*beta(i)
         gamma(i) = hl0*coy(i)*elamda
         gamma2(i) = 2.d0*gamma(i)
         cof1(i) = elamda
      enddo
!-----------------------------------------------------------------------
! Begin iteration loop.
! Load array x with (D-inverse)*z for first iteration.
!-----------------------------------------------------------------------
      iter = 1
!
      do jy = 1 , my
         iyoff = mxns*(jy-1)
         do jx = 1 , mx
            ic = iyoff + ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = cof1(i)*z(ici)
               z(ici) = 0.d0
            enddo
         enddo
      enddo
      do
!-----------------------------------------------------------------------
! Calculate (I - (D-inverse)*L)-inverse * x.
!-----------------------------------------------------------------------
         jy = 1
         do jx = 2 , mx - 1
            ic = ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = x(ici) + beta(i)*x(ici-ns)
            enddo
         enddo
         jx = mx
         ic = ns*(jx-1)
         do i = 1 , ns
            ici = ic + i
            x(ici) = x(ici) + beta2(i)*x(ici-ns)
         enddo
         do jy = 2 , my - 1
            iyoff = mxns*(jy-1)
            jx = 1
            ic = iyoff
            do i = 1 , ns
               ici = ic + i
               x(ici) = x(ici) + gamma(i)*x(ici-mxns)
            enddo
            do jx = 2 , mx - 1
               ic = iyoff + ns*(jx-1)
               do i = 1 , ns
                  ici = ic + i
                  x(ici) = (x(ici)+beta(i)*x(ici-ns)) + gamma(i)        &
                         & *x(ici-mxns)
               enddo
            enddo
            jx = mx
            ic = iyoff + ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = (x(ici)+beta2(i)*x(ici-ns)) + gamma(i)          &
                      & *x(ici-mxns)
            enddo
         enddo
         jy = my
         iyoff = mxns*(jy-1)
         jx = 1
         ic = iyoff
         do i = 1 , ns
            ici = ic + i
            x(ici) = x(ici) + gamma2(i)*x(ici-mxns)
         enddo
         do jx = 2 , mx - 1
            ic = iyoff + ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = (x(ici)+beta(i)*x(ici-ns)) + gamma2(i)          &
                      & *x(ici-mxns)
            enddo
         enddo
         jx = mx
         ic = iyoff + ns*(jx-1)
         do i = 1 , ns
            ici = ic + i
            x(ici) = (x(ici)+beta2(i)*x(ici-ns)) + gamma2(i)*x(ici-mxns)
         enddo
!-----------------------------------------------------------------------
! Add increment x to z.
!-----------------------------------------------------------------------
         do i = 1 , n
            z(i) = z(i) + x(i)
         enddo
!
         if ( iter>=itmax ) exit
!-----------------------------------------------------------------------
! Calculate (D-inverse)*U*x.
!-----------------------------------------------------------------------
         iter = iter + 1
         jy = 1
         jx = 1
         ic = ns*(jx-1)
         do i = 1 , ns
            ici = ic + i
            x(ici) = beta2(i)*x(ici+ns) + gamma2(i)*x(ici+mxns)
         enddo
         do jx = 2 , mx - 1
            ic = ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = beta(i)*x(ici+ns) + gamma2(i)*x(ici+mxns)
            enddo
         enddo
         jx = mx
         ic = ns*(jx-1)
         do i = 1 , ns
            ici = ic + i
            x(ici) = gamma2(i)*x(ici+mxns)
         enddo
         do jy = 2 , my - 1
            iyoff = mxns*(jy-1)
            jx = 1
            ic = iyoff
            do i = 1 , ns
               ici = ic + i
               x(ici) = beta2(i)*x(ici+ns) + gamma(i)*x(ici+mxns)
            enddo
            do jx = 2 , mx - 1
               ic = iyoff + ns*(jx-1)
               do i = 1 , ns
                  ici = ic + i
                  x(ici) = beta(i)*x(ici+ns) + gamma(i)*x(ici+mxns)
               enddo
            enddo
            jx = mx
            ic = iyoff + ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = gamma(i)*x(ici+mxns)
            enddo
         enddo
         jy = my
         iyoff = mxns*(jy-1)
         jx = 1
         ic = iyoff
         do i = 1 , ns
            ici = ic + i
            x(ici) = beta2(i)*x(ici+ns)
         enddo
         do jx = 2 , mx - 1
            ic = iyoff + ns*(jx-1)
            do i = 1 , ns
               ici = ic + i
               x(ici) = beta(i)*x(ici+ns)
            enddo
         enddo
         jx = mx
         ic = iyoff + ns*(jx-1)
         do i = 1 , ns
            ici = ic + i
            x(ici) = 0.0d0
         enddo
      enddo
end subroutine gs
