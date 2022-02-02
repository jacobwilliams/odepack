!>
!! Demonstration program for the DLSODA package.
!! This is the version of 14 June 2001.
!!
!! This version is in double precision.
!!
!! The package is used to solve two simple problems,
!! one with a full Jacobian, the other with a banded Jacobian,
!! with the 2 appropriate values of jt in each case.
!! If the errors are too large, or other difficulty occurs,
!! a warning message is printed.  All output is on unit lout = 6.
!!-----------------------------------------------------------------------
program lsoda
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: atol , dtout0 , dtout1 , er , erm , ero , hu ,    &
 & rtol , t , tout , tsw
real(kind=dp) , save :: dtout , tout1
integer :: i , iopar , iopt , iout , istate , itask , itol , jt , &
 & leniw , lenrw , liw , lrw , mband , ml , mu , mused ,  &
 & neq , nerr , nfe , nfea , nje , nout , nqu , nst
integer , dimension(45) :: iwork
integer , save :: lout
real(kind=dp) , dimension(522) :: rwork
real(kind=dp) , dimension(25) :: y
external f1 , f2 , jac1 , jac2
!
data lout/6/ , tout1/16.921743d0/ , dtout/17.341162d0/
!
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-8
      lrw = 522
      liw = 45
      iopt = 0
!
! First problem
!
      neq = 2
      nout = 4
      write (lout,99001) neq , itol , rtol , atol
99001 format (/'Demonstration program for DLSODA package'////           &
             &' Problem 1:   Van der Pol oscillator:'/                  &
             &'              xdotdot - 20*(1 - x**2)*xdot + x = 0, ',   &
             &'   x(0) = 2, xdot(0) = 0'/' neq =',i2/' itol =',i3,      &
             &'   rtol =',d10.1,'   atol =',d10.1//)
!
      do jt = 1 , 2
         write (lout,99002) jt
99002    format (//' Solution with jt =',                               &
                &i3//'  t               x               xdot       meth'&
               & ,'   nq     h           tsw'//)
         t = 0.0d0
         y(1) = 2.0d0
         y(2) = 0.0d0
         itask = 1
         istate = 1
         dtout0 = 0.5d0*tout1
         dtout1 = 0.5d0*dtout
         tout = dtout0
         ero = 0.0d0
         do iout = 1 , nout
            call dlsoda(f1,[neq],y,t,tout,itol,[rtol],[atol],itask,istate,    &
                      & iopt,rwork,lrw,iwork,liw,jac1,jt)
            hu = rwork(11)
            tsw = rwork(15)
            nqu = iwork(14)
            mused = iwork(19)
            write (lout,99003) t , y(1) , y(2) , mused , nqu , hu , tsw
99003       format (d12.5,d16.5,d14.3,2i6,2d13.3)
            if ( istate<0 ) exit
            iopar = iout - 2*(iout/2)
            if ( iopar==0 ) then
               er = abs(y(1))
               ero = max(ero,er)
               if ( er>1.0d-2 ) then
                  write (lout,99010)
                  nerr = nerr + 1
               endif
            endif
            if ( iout==1 ) tout = tout + dtout0
            if ( iout>1 ) tout = tout + dtout1
         enddo
         if ( istate<0 ) nerr = nerr + 1
         nst = iwork(11)
         nfe = iwork(12)
         nje = iwork(13)
         lenrw = iwork(17)
         leniw = iwork(18)
         nfea = nfe
         if ( jt==2 ) nfea = nfe - neq*nje
         write (lout,99004) lenrw , leniw , nst , nfe , nfea , nje , ero
99004    format (//' Final statistics for this run:'/' rwork size =',i4,&
                &'   iwork size =',i4/' number of steps =',             &
                &i5/' number of f-s   =',i5/' (excluding J-s) =',       &
                &i5/' number of J-s   =',i5/' max. error at root =',    &
               & d10.2)
      enddo
!
! Second problem
!
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      atol = 1.0d-6
      nout = 5
      write (lout,99005) neq , ml , mu , itol , rtol , atol
99005 format (///80('-')///' Problem 2: ydot = A * y , where',          &
             &'  A is a banded lower triangular matrix'/                &
             &'            derived from 2-D advection PDE'/' neq =',i3, &
             &'   ml =',i2,'   mu =',i2/' itol =',i3,'   rtol =',d10.1, &
             &'   atol =',d10.1//)
      do jt = 4 , 5
         write (lout,99006) jt
99006    format (//' Solution with jt =',                               &
                &i3//'     t             max.err.     meth   ',         &
                &'nq      h            tsw'//)
         t = 0.0d0
         do i = 2 , neq
            y(i) = 0.0d0
         enddo
         y(1) = 1.0d0
         itask = 1
         istate = 1
         tout = 0.01d0
         ero = 0.0d0
         do iout = 1 , nout
            call dlsoda(f2,[neq],y,t,tout,itol,[rtol],[atol],itask,istate,    &
                      & iopt,rwork,lrw,iwork,liw,jac2,jt)
            call edit2(y,t,erm)
            hu = rwork(11)
            tsw = rwork(15)
            nqu = iwork(14)
            mused = iwork(19)
            write (lout,99007) t , erm , mused , nqu , hu , tsw
99007       format (d15.5,d14.3,2i6,2d14.3)
            if ( istate<0 ) exit
            er = erm/atol
            ero = max(ero,er)
            if ( er>1000.0d0 ) then
               write (lout,99010)
               nerr = nerr + 1
            endif
            tout = tout*10.0d0
         enddo
         if ( istate<0 ) nerr = nerr + 1
         nst = iwork(11)
         nfe = iwork(12)
         nje = iwork(13)
         lenrw = iwork(17)
         leniw = iwork(18)
         nfea = nfe
         if ( jt==5 ) nfea = nfe - mband*nje
         write (lout,99008) lenrw , leniw , nst , nfe , nfea , nje , ero
99008    format (//' Final statistics for this run:'/' rwork size =',i4,&
                &'   iwork size =',i4/' number of steps =',             &
                &i5/' number of f-s   =',i5/' (excluding J-s) =',       &
                &i5/' number of J-s   =',i5/' error overrun =',d10.2)
      enddo
      write (lout,99009) nerr
99009 format (///' Number of errors encountered =',i3)
99010 format (//' Warning: value at root exceeds 1.0d-2'//)
end program lsoda

subroutine f1(neq,t,y,ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: y
real(kind=dp) , dimension(neq) :: ydot
!
      ydot(1) = y(2)
      ydot(2) = 20.0d0*(1.0d0-y(1)*y(1))*y(2) - y(1)
end subroutine f1

subroutine jac1(neq,t,y,ml,mu,pd,nrowpd)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: y
integer :: ml
integer :: mu
real(kind=dp) , dimension(nrowpd,neq) :: pd
integer :: nrowpd
!
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -40.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 20.0d0*(1.0d0-y(1)*y(1))
end subroutine jac1

subroutine f2(neq,t,y,ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer , parameter  ::  ng = 5
real(kind=dp) , parameter  ::  alph1 = 1.0d0 , alph2 = 1.0d0
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: y
real(kind=dp) , dimension(neq) :: ydot
!
real(kind=dp) :: d
integer :: i , j , k
!
      do j = 1 , ng
         do i = 1 , ng
            k = i + (j-1)*ng
            d = -2.0d0*y(k)
            if ( i/=1 ) d = d + y(k-1)*alph1
            if ( j/=1 ) d = d + y(k-ng)*alph2
            ydot(k) = d
         enddo
      enddo
end subroutine f2

subroutine jac2(neq,t,y,ml,mu,pd,nrowpd)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: y
integer :: ml
integer :: mu
real(kind=dp) , dimension(nrowpd,neq) :: pd
integer :: nrowpd
!
real(kind=dp) , parameter  ::  alph1 = 1.0d0 , alph2 = 1.0d0
integer, parameter :: ng=5
integer :: j , mband , mu1 , mu2
!
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do j = 1 , neq
         pd(mu1,j) = -2.0d0
         pd(mu2,j) = alph1
         pd(mband,j) = alph2
      enddo
      do j = ng , neq , ng
         pd(mu2,j) = 0.0d0
      enddo
end subroutine jac2

subroutine edit2(y,t,erm)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer , parameter  ::  ng = 5
real(kind=dp) , parameter  ::  alph1 = 1.0d0 , alph2 = 1.0d0
!
real(kind=dp) , intent(in) , dimension(*) :: y
real(kind=dp) , intent(in) :: t
real(kind=dp) , intent(out) :: erm
!
real(kind=dp) :: a1 , a2 , er , ex , yt
integer :: i , j , k
!
      erm = 0.0d0
      if ( t==0.0d0 ) return
      ex = 0.0d0
      if ( t<=30.0d0 ) ex = exp(-2.0d0*t)
      a2 = 1.0d0
      do j = 1 , ng
         a1 = 1.0d0
         do i = 1 , ng
            k = i + (j-1)*ng
            yt = t**(i+j-2)*ex*a1*a2
            er = abs(y(k)-yt)
            erm = max(erm,er)
            a1 = a1*alph1/i
         enddo
         a2 = a2*alph2/j
      enddo
end subroutine edit2
