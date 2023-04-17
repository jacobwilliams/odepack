program dlsoda_ex
use m_odepack
implicit none
external fex
external jdum
integer,parameter            ::  dp=kind(0.0d0)
real(kind=dp),dimension(3)   ::  atol,y
integer                      ::  iopt,iout,istate,itask,itol,jt,liw,lrw,neq
integer,dimension(23)        ::  iwork
real(kind=dp)                ::  rtol,t,tout
real(kind=dp),dimension(70)  ::  rwork

   call reference()

   neq = 3
   y(1) = 1.
   y(2) = 0.
   y(3) = 0.
   t = 0.
   tout = .4
   itol = 2
   rtol = 1.D-4
   atol(1) = 1.D-6
   atol(2) = 1.D-10
   atol(3) = 1.D-6
   itask = 1
   istate = 1
   iopt = 0
   lrw = 70
   liw = 23
   jt = 2
   do iout = 1,12
      call dlsoda(fex,[neq],y,t,tout,itol,[rtol],atol,itask,istate,iopt,   &
                & rwork,lrw,iwork,liw,jdum,jt)
      write (6,99010) t,y(1),y(2),y(3)
      99010 format (' At t =',d12.4,'   Y =',3D14.6)
      if ( istate<0 ) then
         write (6,99020) istate
         99020 format (///' Error halt.. ISTATE =',i3)
         stop 1
      else
         tout = tout*10.
      endif
   enddo
   write (6,99030) iwork(11),iwork(12),iwork(13),iwork(19),        &
                 & rwork(15)
   99030 format (/' No. steps =',i4,'  No. f-s =',i4,'  No. J-s =',        &
                &i4/' Method last used =',i2,'   Last switch was at t =',  &
               & d12.4)

end program dlsoda_ex

subroutine jdum()
implicit none
end subroutine jdum
 
subroutine fex(Neq,T,Y,Ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
 
integer                                   ::  Neq
real(kind=dp)                             ::  T
real(kind=dp),intent(in),dimension(3)     ::  Y
real(kind=dp),intent(inout),dimension(3)  ::  Ydot
 
   Ydot(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
   Ydot(3) = 3.D7*Y(2)*Y(2)
   Ydot(2) = -Ydot(1) - Ydot(3)
end subroutine fex

subroutine reference()
implicit none
integer :: i
character(len=:),allocatable :: original(:)
write(*,'(a)') '>PROCEDURE dlsoda'
original=[ CHARACTER(LEN=128) :: &
'']

original=[ CHARACTER(LEN=128) :: &
' The output of this program (on a CDC-7600 in single precision)',&
' is as follows:',&
'',&
'   At t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02',&
'   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02',&
'   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01',&
'   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01',&
'   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01',&
'   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01',&
'   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01',&
'   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01',&
'   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01',&
'   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01',&
'   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01',&
'   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00',&
'',&
'   No. steps = 361  No. f-s = 693  No. J-s =  64',&
'   Method last used = 2   Last switch was at t =  6.0092e-03',&
'']

   write(*,'(a)') (trim(original(i)),i=1,size(original))
   write(*,'(a)') ' >NOW IT PRODUCES'
end subroutine reference
