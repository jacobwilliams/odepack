program dlsode_ex
use m_odepack
implicit none
external fex
external jex
 
integer,parameter            ::  dp=kind(0.0d0)
real(kind=dp),dimension(3)   ::  atol,y
integer                      ::  iopt,iout,istate,itask,itol,liw,lrw,mf,neq
integer,dimension(23)        ::  iwork
real(kind=dp)                ::  rtol,t,tout
real(kind=dp),dimension(58)  ::  rwork

   call reference()
 
   neq = 3
   y(1) = 1.D0
   y(2) = 0.D0
   y(3) = 0.D0
   t = 0.D0
   tout = .4D0
   itol = 2
   rtol = 1.D-4
   atol(1) = 1.D-6
   atol(2) = 1.D-10
   atol(3) = 1.D-6
   itask = 1
   istate = 1
   iopt = 0
   lrw = 58
   liw = 23
   mf = 21
   do iout = 1,12
      call dlsode(fex,[neq],y,t,tout,itol,[rtol],atol,itask,istate,iopt,   &
                & rwork,lrw,iwork,liw,jex,mf)
      write (6,99010) t,y(1),y(2),y(3)
   99010 format (' At t =',d12.4,'   y =',3D14.6)
      if ( istate<0 ) then
         write (6,99020) istate
   99020 format (///' Error halt.. ISTATE =',i3)
         stop 1
      else
         tout = tout*10.D0
      endif
   enddo
   write (6,99030) iwork(11),iwork(12),iwork(13)
   99030 format (/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)

end program dlsode_ex
 
subroutine fex(Neq,T,Y,Ydot)
implicit none
integer,parameter                         ::  dp=kind(0.0d0)
 
integer                                   ::  Neq
real(kind=dp)                             ::  T
real(kind=dp),intent(in),dimension(3)     ::  Y
real(kind=dp),intent(inout),dimension(3)  ::  Ydot
 
   Ydot(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
   Ydot(3) = 3.D7*Y(2)*Y(2)
   Ydot(2) = -Ydot(1) - Ydot(3)
end subroutine fex
 
subroutine jex(Neq,T,Y,Ml,Mu,Pd,Nrpd)
implicit none
 
integer,parameter                              ::  dp=kind(0.0d0)
integer                                        ::  Neq
real(kind=dp)                                  ::  T
real(kind=dp),intent(in),dimension(3)          ::  Y
integer                                        ::  Ml
integer                                        ::  Mu
integer,intent(in)                             ::  Nrpd
real(kind=dp),intent(inout),dimension(Nrpd,3)  ::  Pd
 
   Pd(1,1) = -.04D0
   Pd(1,2) = 1.D4*Y(3)
   Pd(1,3) = 1.D4*Y(2)
   Pd(2,1) = .04D0
   Pd(2,3) = -Pd(1,3)
   Pd(3,2) = 6.D7*Y(2)
   Pd(2,2) = -Pd(1,2) - Pd(3,2)
end subroutine jex

subroutine reference()
implicit none
integer :: i
character(len=:),allocatable :: original(:)
write(*,'(a)') '>PROCEDURE dlsode'
original=[ CHARACTER(LEN=128) :: &
'     The output from this program (on a Cray-1 in single precision)',&
'     is as follows.',&
'',&
'     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02',&
'     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02',&
'     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01',&
'     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01',&
'     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01',&
'     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01',&
'     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01',&
'     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01',&
'     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01',&
'     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01',&
'     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01',&
'     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00',&
'']

   write(*,'(a)') (trim(original(i)),i=1,size(original))
   write(*,'(a)') ' >NOW IT PRODUCES'
end subroutine reference
