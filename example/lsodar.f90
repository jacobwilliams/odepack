!>
!! Demonstration program for the DLSODAR package.
!! This is the version of 14 June 2001.
!!
!! This version is in double precision.
!!
!! The DLSODAR package is used to solve two simple problems,
!! one nonstiff and one intermittently stiff.
!! If the errors are too large, or other difficulty occurs,
!! a warning message is printed.  All output is on unit lout = 6.
!-----------------------------------------------------------------------
program lsodar
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) , dimension(2) :: atol , y
real(kind=dp) :: er , ero , errt , rtol , t , tout , tzero , yt
external :: jac_dum
integer :: iopt , iout , istate , itask , itol , jt ,      &
 & kroot , leniw , lenrw , liw , lrw , neq , nerr , nfe , &
 & nfea , ng , nge , nje , nst
integer , dimension(22) :: iwork
integer , dimension(2) :: jroot
integer , save :: lout
real(kind=dp) , dimension(57) :: rwork
external f1 , f2 , gr1 , gr2 , jac2
!
      data lout/6/
!
      nerr = 0
!-----------------------------------------------------------------------
! First problem.
! The initial value problem is:
!   dy/dt = ((2*log(y) + 8)/t - 5)*y,  y(1) = 1,  1 .le. t .le. 6
! The solution is  y(t) = exp(-t**2 + 5*t - 4)
! The two root functions are:
!   g1 = ((2*log(y)+8)/t - 5)*y (= dy/dt)  (with root at t = 2.5),
!   g2 = log(y) - 2.2491  (with roots at t = 2.47 and 2.53)
!-----------------------------------------------------------------------
! Set all input parameters and print heading.
      neq = 1
      y(1) = 1.0d0
      t = 1.0d0
      tout = 2.0d0
      itol = 1
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      itask = 1
      istate = 1
      iopt = 0
      lrw = 44
      liw = 21
      jt = 2
      ng = 2
      write (lout,99001) itol , rtol , atol(1) , jt
99001 format (/' Demonstration program for DLSODAR package'////         &
             &' First problem'//                                        &
             &/' Problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1' &
            & //' Solution is  y(t) = exp(-t**2 + 5*t - 4)'//           &
             &' Root functions are:'/10x,                               &
             &' g1 = dy/dt  (root at t = 2.5)'/10x,                     &
             &' g2 = log(y) - 2.2491  (roots at t = 2.47 and t = 2.53)' &
            & //' itol =',i3,'   rtol =',d10.1,'   atol =',             &
             &d10.1//' jt =',i3///)
!
! Call DLSODAR in loop over tout values 2,3,4,5,6.
      ero = 0.0d0
      do iout = 1 , 5
         do
            call dlsodar(f1,[neq],y,t,tout,itol,[rtol],atol,itask,istate,   &
                       & iopt,rwork,lrw,iwork,liw,jac_dum,jt,gr1,ng,jroot)
!
! Print y and error in y, and print warning if error too large.
            yt = exp(-t*t+5.0d0*t-4.0d0)
            er = y(1) - yt
            write (lout,99002) t , y(1) , er
99002       format (' At t =',d15.7,5x,'y =',d15.7,5x,'error =',d12.4)
            if ( istate<0 ) goto 100
            er = abs(er)/(rtol*abs(y(1))+atol(1))
            ero = max(ero,er)
            if ( er>1000.0d0 ) then
               write (lout,99003)
99003          format (//' Warning: error exceeds 1000 * tolerance'//)
               nerr = nerr + 1
            endif
            if ( istate/=3 ) then
!
! If no root found, increment tout and loop back.
               tout = tout + 1.0d0
               exit
            else
!
! If a root was found, write results and check root location.
! Then reset istate to 2 and return to DLSODAR call.
               write (lout,99004) t , jroot(1) , jroot(2)
99004          format (/' Root found at t =',d15.7,5x,'jroot =',2i5)
               if ( jroot(1)==1 ) errt = t - 2.5d0
               if ( jroot(2)==1 .and. t<=2.5d0 ) errt = t - 2.47d0
               if ( jroot(2)==1 .and. t>2.5d0 ) errt = t - 2.53d0
               write (lout,99005) errt
99005          format (' Error in t location of root is',d12.4/)
               if ( abs(errt)>1.0d-3 ) then
                  write (lout,99006)
99006             format (//' Warning: root error exceeds 1.0d-3'//)
                  nerr = nerr + 1
               endif
               istate = 2
            endif
         enddo
      enddo
!
! Problem complete.  Print final statistics.
 100  continue
      if ( istate<0 ) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if ( jt==2 ) nfea = nfe - neq*nje
      write (lout,99007) lenrw , leniw , nst , nfe , nfea , nje , nge , &
                       & ero
99007 format (//' Final statistics for this run:'/' rwork size =',i4,   &
             &'   iwork size =',i4/' number of steps =',                &
             &i5/' number of f-s   =',i5/' (excluding j-s) =',          &
             &i5/' number of j-s   =',i5/' number of g-s   =',          &
             &i5/' error overrun =',d10.2)
!
!-----------------------------------------------------------------------
! Second problem (Van der Pol oscillator).
! The initial value problem is (after reduction of 2nd order ODE):
!   dy1/dt = y2,  dy2/dt = 100*(1 - y1**2)*y2 - y1,
!   y1(0) = 2,  y2(0) = 0,  0 .le. t .le. 200
! The root function is  g = y1.
! An analytic solution is not known, but the zeros of y1 are known
! to 15 figures for purposes of checking the accuracy.
!-----------------------------------------------------------------------
! Set tolerance parameters and print heading.
      itol = 2
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      atol(2) = 1.0d-4
      write (lout,99008) itol , rtol , atol(1) , atol(2)
99008 format (////80('*')//' Second problem (Van der Pol oscillator)'// &
             &' Problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1'&
            & /'            y1(0) = 2,  y2(0) = 0'//                    &
             &' Root function is  g = y1'//' itol =',i3,'   rtol =',    &
            & d10.1,'   atol =',2d10.1)
!
! Loop over jt = 1, 2.  Set remaining parameters and print jt.
      do jt = 1 , 2
         neq = 2
         y(1) = 2.0d0
         y(2) = 0.0d0
         t = 0.0d0
         tout = 20.0d0
         itask = 1
         istate = 1
         iopt = 0
         lrw = 57
         liw = 22
         ng = 1
         write (lout,99009) jt
99009    format (///' Solution with jt =',i2//)
!
! Call DLSODAR in loop over tout values 20,40,...,200.
         do iout = 1 , 10
            do
               call dlsodar(f2,[neq],y,t,tout,itol,[rtol],atol,itask,istate,&
                          & iopt,rwork,lrw,iwork,liw,jac2,jt,gr2,ng,    &
                          & jroot)
!
! Print y1 and y2.
               write (lout,99010) t , y(1) , y(2)
99010          format (' At t =',d15.7,5x,'y1 =',d15.7,5x,'y2 =',d15.7)
               if ( istate<0 ) goto 150
               if ( istate/=3 ) then
!
! If no root found, increment tout and loop back.
                  tout = tout + 20.0d0
                  exit
               else
!
! If a root was found, write results and check root location.
! Then reset istate to 2 and return to DLSODAR call.
                  write (lout,99011) t
99011             format (/' Root found at t =',d15.7)
                  kroot = int(t/81.2d0+0.5d0)
                  tzero = 81.17237787055d0 + (kroot-1)*81.41853556212d0
                  errt = t - tzero
                  write (lout,99012) errt
99012             format (' Error in t location of root is',d12.4//)
                  if ( abs(errt)>1.0d-1 ) then
                     write (lout,99013)
99013                format (//' Warning: root error exceeds 1.0d-1'//)
                     nerr = nerr + 1
                  endif
                  istate = 2
               endif
            enddo
         enddo
!
! Problem complete.  Print final statistics.
 150     continue
         if ( istate<0 ) nerr = nerr + 1
         nst = iwork(11)
         nfe = iwork(12)
         nje = iwork(13)
         nge = iwork(10)
         lenrw = iwork(17)
         leniw = iwork(18)
         nfea = nfe
         if ( jt==2 ) nfea = nfe - neq*nje
         write (lout,99014) lenrw , leniw , nst , nfe , nfea , nje , nge
99014    format (//' Final statistics for this run:'/'  rwork size =',  &
               & i4,'   iwork size =',i4/'  number of steps =',         &
                &i5/'  number of f-s   =',i5/'  (excluding j-s) =',     &
                &i5/'  number of j-s   =',i5/'  number of g-s   =',i5)
      enddo
!
      write (lout,99015) nerr
99015 format (///' Total number of errors encountered =',i3)
end program lsodar

subroutine f1(neq,t,y,ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(1) :: y
real(kind=dp) , dimension(1) :: ydot
!
      ydot(1) = ((2.0d0*log(y(1))+8.0d0)/t-5.0d0)*y(1)
end subroutine f1

subroutine gr1(neq,t,y,ng,groot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(1) :: y
integer :: ng
real(kind=dp) , dimension(2) :: groot
!
      groot(1) = ((2.0d0*log(y(1))+8.0d0)/t-5.0d0)*y(1)
      groot(2) = log(y(1)) - 2.2491d0
end subroutine gr1

subroutine f2(neq,t,y,ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(2) :: y
real(kind=dp) , dimension(2) :: ydot
!
      ydot(1) = y(2)
      ydot(2) = 100.0d0*(1.0d0-y(1)*y(1))*y(2) - y(1)
end subroutine f2

subroutine jac_dum(neq,t,y,ml,mu,pd,nrowpd)
!> dummy JAC routine that is a noop
implicit none
integer,parameter :: dp=kind(0.0d0)
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(2) :: y
integer :: ml
integer :: mu
real(kind=dp) , dimension(nrowpd,2) :: pd
integer :: nrowpd
end subroutine jac_dum

subroutine jac2(neq,t,y,ml,mu,pd,nrowpd)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(2) :: y
integer :: ml
integer :: mu
real(kind=dp) , dimension(nrowpd,2) :: pd
integer :: nrowpd
!
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -200.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 100.0d0*(1.0d0-y(1)*y(1))
end subroutine jac2

subroutine gr2(neq,t,y,ng,groot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(2) :: y
integer :: ng
real(kind=dp) , dimension(1) :: groot
!
      groot(1) = y(1)
end subroutine gr2
