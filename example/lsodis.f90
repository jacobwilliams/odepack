!>
!! Demonstration program for the DLSODIS package.
!! This is the version of 14 June 2001.
!!
!! This version is in double precision.
!!
!! This program solves a semi-discretized form of the Burgers equation,
!!```text
!!     u  = -(u*u/2)  + eta * u
!!      t           x          xx
!!
!!```
!! for  -1 .le. x .le. 1, t .ge. 0.
!!
!! Here eta = 0.05.
!!
!! Boundary conditions: u(-1,t) = u(1,t) and du/dx(-1,t) = du/dx(1,t).
!!
!! Initial profile: square wave
!!     u(0,x) = 0    for 1/2 .lt. abs(x) .le. 1
!!     u(0,x) = 1/2  for abs(x) = 1/2
!!     u(0,x) = 1    for 0 .le. abs(x) .lt. 1/2
!!
!! An ODE system is generated by a simplified Galerkin treatment
!! of the spatial variable x.
!!
!!### Reference:
!! R. C. Y. Chin, G. W. Hedstrom, and K. E. Karlsson,
!! A Simplified Galerkin Method for Hyperbolic Equations,
!! Math. Comp., vol. 33, no. 146 (April 1979), pp. 647-658.
!!
!! The problem is run with the DLSODIS package with a 12-node mesh,
!! for various appropriate values of the method flag mf.
!! Output is on unit lout, set to 6 in a data statement below.
!!
!!### Problem specific data:
!!      npts  = number of unknowns (npts = 0 mod 4)
!!      nnz   = number of non-zeros in Jacobian before fill in
!!      nnza  = number of non-zeros in Jacobian after fill in
!!      lrwk  = length of real work array (taking into account fill in)
!!      liwk  = length of integer work array
!!      ipia  = pointer to ia in iw (ia(j) = iw(ipia+j-1)
!!      ipja  = pointer to ja in iw (ja(j) = iw(ipja+j-1)
!!      ipic  = pointer to ic in iw array (ic(j) = iw(ipic+j-1))
!!      ipjc  = pointer to jc in iw array (jc(j) = iw(ipjc+j-1))
!-----------------------------------------------------------------------
program lsodis
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer , parameter  ::  npts = 12 , nnz = 3*npts ,               &
                     & nnza = 5*npts , lrwk = 20 + 3*nnza +     &
                     & 28*npts , liwk = 32 + 2*nnza + 2*npts ,  &
                     & ipia = 31 , ipja = 31 + npts + 1 ,       &
                     & ll = 32 + npts + nnz , ipic = ll ,       &
                     & ipjc = ll + npts + 1
!
real(kind=dp) :: eodsq , r4d
common /test1 / r4d , eodsq
!
real(kind=dp) , save :: a , b , eta , fourth , half , hun , one ,  &
 & tinit , tlast , zero
real(kind=dp) , dimension(2) , save :: atol , rtol
real(kind=dp) :: delta , errfac , t
real(kind=dp) :: elkup
integer :: i , io , istate , j , liw , lrw , lyh , meth , mf ,    &
 & miter , moss , n , n14 , n34 , nerr
integer , save :: itol , lout , nout
integer , dimension(lrwk) :: iw
real(kind=dp) , dimension(lrwk) :: rw
real(kind=dp) , dimension(4) , save :: tout
real(kind=dp) , dimension(npts) :: y , ydoti
external addasp , jacsp , res
!
! Pass problem parameters in the Common block test1.
!
! Set problem parameters and run parameters
data eta/0.05d0/ , a/ - 1.0d0/ , b/1.0d0/
data zero/0.0d0/ , fourth/0.25d0/ , half/0.5d0/ , one/1.0d0/ ,    &
 & hun/100.0d0/
data tinit/0.0d0/ , tlast/0.4d0/
data tout/0.10d0 , 0.20d0 , 0.30d0 , 0.40d0/
data lout/6/ , nout/4/
data itol/1/ , rtol/1.0d-3 , 1.0d-6/ , atol/1.0d-3 , 1.0d-6/
!
      nerr = 0
      lrw = lrwk
      liw = liwk
!
! Compute the mesh width delta and other parameters.
      delta = (b-a)/npts
      r4d = fourth/delta
      eodsq = eta/delta**2
      n14 = npts/4 + 1
      n34 = 3*(npts/4) + 1
      n = npts
!
! Set the initial profile (for output purposes only).
      do i = 1 , n
         y(i) = zero
      enddo
      y(n14) = half
      do i = n14 + 1 , n34 - 1
         y(i) = one
      enddo
      y(n34) = half
!
      write (lout,99001)
!
99001 format (20x,' Demonstration Program for DLSODIS')
      write (lout,99002) eta , a , b , tinit , tlast , n
99002 format (//10x,'-- Simplified Galerkin solution of ',              &
             &'Burgers equation --'///13x,                              &
             &'Diffusion coefficient is eta =',d10.2/13x,               &
             &'Uniform mesh on interval',d12.3,' to ',d12.3/13x,        &
             &'Periodic boundary conditions'/13x,                       &
             &'Initial data are as follows:'//20x,'t0 = ',d12.5/20x,    &
             &'tlast = ',d12.5/20x,'n  = ',i3//)
      write (lout,99003) (y(i),i=1,n)
!
99003 format (/'Initial profile:',/20(6d12.4/))
!
! Set the initial sparse data structures for coefficient matrix A
! and the Jacobian matrix C
      call struct(iw(ipia),iw(ipja),iw(ipic),iw(ipjc),n)
!
! The j loop is over error tolerances.
      do j = 1 , 2
!
! This method flag loop is for demonstration only.
         do moss = 0 , 4
            do meth = 1 , 2
               do miter = 1 , 2
                  mf = 100*moss + 10*meth + miter
!
! Set the initial profile.
                  do i = 1 , n
                     y(i) = zero
                  enddo
                  y(n14) = half
                  do i = n14 + 1 , n34 - 1
                     y(i) = one
                  enddo
                  y(n34) = half
!
                  t = tinit
                  istate = 0
!
                  write (lout,99004) itol , rtol(j) , atol(j) , mf
!
99004             format (///85('*')///'Run with itol =',i2,'  rtol =', &
                        & d12.2,'  atol =',d12.2,'   mf = ',i3//)
!
! output loop for each case
                  do io = 1 , nout
!
! Call DLSODIS and print results.
                     call dlsodis(res,addasp,jacsp,[n],y,ydoti,t,tout(io),&
                                & itol,rtol(j),atol(j),1,istate,0,rw,   &
                                & lrw,iw,liw,mf)
                     write (lout,99005) t , rw(11) , iw(14) ,           &
                                      & (y(i),i=1,n)
!
99005                format (' Output for time t =',d12.5,              &
                            &'  current h =',d12.5,'  current order =', &
                           & i2/20(6d12.4/))
!
! If istate is not 2 on return, print message and go to next case.
                     if ( istate/=2 ) then
                        write (lout,99006) mf , t , istate
!
99006                   format (/'Final time reached for mf = ',i3,     &
                               &' was t = ',d12.5/25x,                  &
                               &'at which istate = ',i2//)
                        nerr = nerr + 1
                        goto 10
                     endif
                  enddo
                  write (lout,99007) mf , iw(11) , iw(12) , iw(13) ,    &
                                   & iw(17) , iw(18) , iw(20) , iw(21)
!
99007             format (/'Final statistics for mf = ',i3,': ',i5,     &
                         &' steps,',i6,' res,',i6,' Jacobians,'/20x,    &
                         &' rw size =',i6,',    iw size =',i6/20x,i4,   &
                         &' extra res for each jac,',i4,' decomps')
!
! Estimate final error and print result.
                  lyh = iw(22)
                  errfac = elkup(n,y,rw(lyh),itol,rtol(j),atol(j),lout)
                  if ( errfac<hun ) then
                     write (lout,99008) errfac
99008                format ('Final output is correct to within ',d9.2, &
                            &'  times local error tolerance.'/)
                  else
                     write (lout,99009) errfac
99009                format ('Final output is wrong by ',d9.2,          &
                            &'  times local error tolerance.'/)
                     nerr = nerr + 1
                  endif
 10            enddo
            enddo
         enddo
      enddo
!
      write (lout,99010) nerr
99010 format (///85('*')                                                &
             &//'Run completed: number of errors encountered =',i3)
!
! end of main program for the DLSODIS demonstration program
end program lsodis

subroutine struct(ia,ja,ic,jc,n)
implicit none
!
integer , dimension(*) :: ia
integer , dimension(*) :: ja
integer , dimension(*) :: ic
integer , dimension(*) :: jc
integer :: n
!
integer :: jj , k , l , m
!
      write (6,99001)
99001 format ('Initial sparse data structures'/)
      k = 0
      do l = 1 , n
         ia(l) = (l-1)*3 + 1
         ic(l) = (l-1)*3 + 1
         do m = l , l + 2
            k = k + 1
            ja(k) = m - 1
            jc(k) = m - 1
         enddo
      enddo
      ia(n+1) = 3*n + 1
      ic(n+1) = 3*n + 1
      ja(1) = n
      jc(1) = n
      ja(k) = 1
      jc(k) = 1
!
      write (6,99002) (ia(jj),jj=1,n+1)
99002 format (' ia  ',15i4/10(5x,15i4/))
      write (6,99003) (ja(jj),jj=1,k)
99003 format (' ja  ',15i4/10(5x,15i4/))
      write (6,99004) (ic(jj),jj=1,n+1)
99004 format (' ic  ',15i4/10(5x,15i4/))
      write (6,99005) (jc(jj),jj=1,k)
99005 format (' jc  ',15i4/10(5x,15i4/))
      end subroutine struct

subroutine res(n,t,y,v,r,ires)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: eodsq , r4d
common /test1 / r4d , eodsq
!
integer :: n
real(kind=dp) :: t
real(kind=dp) , dimension(n) :: y
real(kind=dp) , dimension(n) :: v
real(kind=dp) , dimension(n) :: r
integer :: ires
!
real(kind=dp) :: fact1 , fact4
real(kind=dp) , save :: four , one , six
integer :: i
!
data one/1.0d0/ , four/4.0d0/ , six/6.0d0/
!
      call gfun(n,t,y,r)
      if ( ires==-1 ) return
!
      fact1 = one/six
      fact4 = four/six
!
      r(1) = r(1) - (fact4*v(1)+fact1*(v(2)+v(n)))
      do i = 2 , n - 1
         r(i) = r(i) - (fact4*v(i)+fact1*(v(i-1)+v(i+1)))
      enddo
      r(n) = r(n) - (fact4*v(n)+fact1*(v(1)+v(n-1)))
end subroutine res

subroutine gfun(n,t,y,g)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: eodsq , r4d
common /test1 / r4d , eodsq
!
integer :: n
real(kind=dp) :: t
real(kind=dp) , dimension(n) :: y
real(kind=dp) , dimension(n) :: g
!
integer :: i
real(kind=dp) , save :: two
!
data two/2.0d0/
!
      g(1) = r4d*(y(n)**2-y(2)**2) + eodsq*(y(2)-two*y(1)+y(n))
!
      do i = 2 , n - 1
         g(i) = r4d*(y(i-1)**2-y(i+1)**2)                               &
              & + eodsq*(y(i+1)-two*y(i)+y(i-1))
      enddo
!
      g(n) = r4d*(y(n-1)**2-y(1)**2) + eodsq*(y(1)-two*y(n)+y(n-1))
!
end subroutine gfun

subroutine addasp(n,t,y,j,ip,jp,pa)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: n
real(kind=dp) :: t
real(kind=dp) , dimension(n) :: y
integer :: j
integer , dimension(*) :: ip
integer , dimension(*) :: jp
real(kind=dp) , dimension(n) :: pa
!
real(kind=dp) :: fact1 , fact4
real(kind=dp) , save :: four , one , six
integer :: jm1 , jp1
!
data one/1.0d0/ , four/4.0d0/ , six/6.0d0/
!
! Compute the elements of A.
      fact1 = one/six
      fact4 = four/six
      jm1 = j - 1
      jp1 = j + 1
      if ( j==n ) jp1 = 1
      if ( j==1 ) jm1 = n
!
! Add the matrix A to the matrix pa (sparse).
      pa(j) = pa(j) + fact4
      pa(jp1) = pa(jp1) + fact1
      pa(jm1) = pa(jm1) + fact1
end subroutine addasp

subroutine jacsp(n,t,y,s,j,ip,jp,pdj)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: eodsq , r4d
common /test1 / r4d , eodsq
!
integer :: n
real(kind=dp) :: t
real(kind=dp) , dimension(n) :: y
real(kind=dp) , dimension(n) :: s
integer :: j
integer , dimension(*) :: ip
integer , dimension(*) :: jp
real(kind=dp) , dimension(n) :: pdj
!
real(kind=dp) :: diag , r2d
integer :: jm1 , jp1
real(kind=dp) , save :: two
!
data two/2.0d0/
!
      diag = -two*eodsq
      r2d = two*r4d
      jm1 = j - 1
      jp1 = j + 1
      if ( j==1 ) jm1 = n
      if ( j==n ) jp1 = 1
!
      pdj(jm1) = -r2d*y(j) + eodsq
      pdj(j) = diag
      pdj(jp1) = r2d*y(j) + eodsq
end subroutine jacsp

function elkup(n,y,ewt,itol,rtol,atol,lout)
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
real(kind=dp) :: elkup
!
integer :: n
real(kind=dp) , dimension(n) :: y
real(kind=dp) , dimension(n) :: ewt
integer :: itol
real(kind=dp) :: rtol
real(kind=dp) :: atol
integer :: lout
!
integer :: i
real(kind=dp) , dimension(12) , save :: y12
real(kind=dp) , dimension(120) :: y120
real(kind=dp) , dimension(16) , save :: y120a , y120b , y120c ,    &
 & y120d , y120e , y120f , y120g
real(kind=dp) , dimension(8) , save :: y120h
!
equivalence (y120a(1),y120(1))
equivalence (y120b(1),y120(17))
equivalence (y120c(1),y120(33))
equivalence (y120d(1),y120(49))
equivalence (y120e(1),y120(65))
equivalence (y120f(1),y120(81))
equivalence (y120g(1),y120(97))
equivalence (y120h(1),y120(113))
      data y12/1.60581860d-02 , 3.23063251d-02 , 1.21903380d-01 ,       &
         & 2.70943828d-01 , 4.60951522d-01 , 6.57571216d-01 ,           &
         & 8.25154453d-01 , 9.35644796d-01 , 9.90167557d-01 ,           &
         & 9.22421221d-01 , 5.85764902d-01 , 1.81112615d-01/
      data y120a/1.89009068d-02 , 1.63261891d-02 , 1.47080563d-02 ,     &
         & 1.39263623d-02 , 1.38901341d-02 , 1.45336989d-02 ,           &
         & 1.58129308d-02 , 1.77017162d-02 , 2.01886844d-02 ,           &
         & 2.32742221d-02 , 2.69677715d-02 , 3.12854037d-02 ,           &
         & 3.62476563d-02 , 4.18776225d-02 , 4.81992825d-02 ,           &
         & 5.52360652d-02/
      data y120b/6.30096338d-02 , 7.15388849d-02 , 8.08391507d-02 ,     &
         & 9.09215944d-02 , 1.01792784d-01 , 1.13454431d-01 ,           &
         & 1.25903273d-01 , 1.39131085d-01 , 1.53124799d-01 ,           &
         & 1.67866712d-01 , 1.83334757d-01 , 1.99502830d-01 ,           &
         & 2.16341144d-01 , 2.33816600d-01 , 2.51893167d-01 ,           &
         & 2.70532241d-01/
      data y120c/2.89693007d-01 , 3.09332757d-01 , 3.29407198d-01 ,     &
         & 3.49870723d-01 , 3.70676646d-01 , 3.91777421d-01 ,           &
         & 4.13124817d-01 , 4.34670077d-01 , 4.56364053d-01 ,           &
         & 4.78157319d-01 , 5.00000270d-01 , 5.21843218d-01 ,           &
         & 5.43636473d-01 , 5.65330432d-01 , 5.86875670d-01 ,           &
         & 6.08223037d-01/
      data y120d/6.29323777d-01 , 6.50129662d-01 , 6.70593142d-01 ,     &
         & 6.90667536d-01 , 7.10307235d-01 , 7.29467947d-01 ,           &
         & 7.48106966d-01 , 7.66183477d-01 , 7.83658878d-01 ,           &
         & 8.00497138d-01 , 8.16665158d-01 , 8.32133153d-01 ,           &
         & 8.46875019d-01 , 8.60868691d-01 , 8.74096465d-01 ,           &
         & 8.86545273d-01/
      data y120e/8.98206892d-01 , 9.09078060d-01 , 9.19160487d-01 ,     &
         & 9.28460742d-01 , 9.36989986d-01 , 9.44763554d-01 ,           &
         & 9.51800339d-01 , 9.58122004d-01 , 9.63751979d-01 ,           &
         & 9.68714242d-01 , 9.73031887d-01 , 9.76725449d-01 ,           &
         & 9.79811001d-01 , 9.82297985d-01 , 9.84186787d-01 ,           &
         & 9.85466039d-01/
      data y120f/9.86109629d-01 , 9.86073433d-01 , 9.85291781d-01 ,     &
         & 9.83673704d-01 , 9.81099057d-01 , 9.77414704d-01 ,           &
         & 9.72431015d-01 , 9.65919133d-01 , 9.57609585d-01 ,           &
         & 9.47193093d-01 , 9.34324619d-01 , 9.18631922d-01 ,           &
         & 8.99729965d-01 , 8.77242371d-01 , 8.50830623d-01 ,           &
         & 8.20230644d-01/
      data y120g/7.85294781d-01 , 7.46035145d-01 , 7.02662039d-01 ,     &
         & 6.55609682d-01 , 6.05541326d-01 , 5.53327950d-01 ,           &
         & 4.99999118d-01 , 4.46670394d-01 , 3.94457322d-01 ,           &
         & 3.44389410d-01 , 2.97337561d-01 , 2.53964948d-01 ,           &
         & 2.14705729d-01 , 1.79770169d-01 , 1.49170367d-01 ,           &
         & 1.22758681d-01/
      data y120h/1.00271052d-01 , 8.13689920d-02 , 6.56761515d-02 ,     &
         & 5.28075160d-02 , 4.23908624d-02 , 3.40811650d-02 ,           &
         & 2.75691506d-02 , 2.25853507d-02/
!
      if ( (n-12)*(n-120)/=0 ) then
!
! error return
         write (lout,99001) n
99001    format (/5x,'Illegal use of elkup for n =',i4)
         elkup = 1.0d3
         goto 99999
      elseif ( n==120 ) then
!
! Compute local error tolerance using correct y (n = 120).
         call dewset(n,itol,[rtol],[atol],y120,ewt)
!
! Invert ewt and replace y by the error, y - ytrue.
         do i = 1 , 120
            ewt(i) = 1.0d0/ewt(i)
            y(i) = y(i) - y120(i)
         enddo
      else
!
! Compute local error tolerance using correct y (n = 12).
         call dewset(n,itol,[rtol],[atol],y12,ewt)
!
! Invert ewt and replace y by the error, y - ytrue.
         do i = 1 , 12
            ewt(i) = 1.0d0/ewt(i)
            y(i) = y(i) - y12(i)
         enddo
      endif
!
! Find weighted norm of the error.
      elkup = dvnorm(n,y,ewt)
      return
99999 continue
end function elkup