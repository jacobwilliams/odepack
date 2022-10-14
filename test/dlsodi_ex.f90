program dlsodi_ex
use m_odepack
implicit none
external aplusp
external dgbydy
external resid
 
integer,parameter            ::  dp=kind(0.0d0)
real(kind=dp),dimension(3)   ::  atol,y,ydoti
integer                      ::  iopt,iout,istate,itask,itol,liw,lrw,mf,neq
integer,dimension(23)        ::  iwork
real(kind=dp)                ::  rtol,t,tout
real(kind=dp),dimension(58)  ::  rwork

   call reference()
 
   neq = 3
   y(1) = 1.
   y(2) = 0.
   y(3) = 0.
   ydoti(1) = -.04
   ydoti(2) = .04
   ydoti(3) = 0.
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
   lrw = 58
   liw = 23
   mf = 21
   do iout = 1,12
      call dlsodi(resid,aplusp,dgbydy,[neq],y,ydoti,t,tout,itol,[rtol],    &
                & atol,itask,istate,iopt,rwork,lrw,iwork,liw,mf)
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
   write (6,99030) iwork(11),iwork(12),iwork(13)
   99030 format (/' No. steps =',i4,'  No. r-s =',i4,'  No. J-s =',i4)
     
end program dlsodi_ex
 
subroutine resid(Neq,T,Y,S,R,Ires)
implicit none
integer,parameter :: dp=kind(0.0d0)
 
integer                                 ::  Neq
real(kind=dp)                           ::  T
real(kind=dp),intent(in),dimension(3)   ::  Y
real(kind=dp),intent(in),dimension(3)   ::  S
real(kind=dp),intent(out),dimension(3)  ::  R
integer                                 ::  Ires
 
   R(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3) - S(1)
   R(2) = .04*Y(1) - 1.D4*Y(2)*Y(3) - 3.D7*Y(2)*Y(2) - S(2)
   R(3) = Y(1) + Y(2) + Y(3) - 1.
end subroutine resid
 
subroutine aplusp(Neq,T,Y,Ml,Mu,P,Nrowp)
implicit none
integer,parameter :: dp=kind(0.0d0)
 
integer                                         ::  Neq
real(kind=dp)                                   ::  T
real(kind=dp),dimension(3)                      ::  Y
integer                                         ::  Ml
integer                                         ::  Mu
integer,intent(in)                              ::  Nrowp
real(kind=dp),intent(inout),dimension(Nrowp,3)  ::  P
 
   P(1,1) = P(1,1) + 1.
   P(2,2) = P(2,2) + 1.
end subroutine aplusp
 
subroutine dgbydy(Neq,T,Y,S,Ml,Mu,P,Nrowp)
implicit none
integer,parameter                             ::  dp=kind(0.0d0)
integer                                       ::  Neq
real(kind=dp)                                 ::  T
real(kind=dp),intent(in),dimension(3)         ::  Y
real(kind=dp),dimension(3)                    ::  S
integer                                       ::  Ml
integer                                       ::  Mu
integer,intent(in)                            ::  Nrowp
real(kind=dp),intent(out),dimension(Nrowp,3)  ::  P
 
   P(1,1) = -.04
   P(1,2) = 1.D4*Y(3)
   P(1,3) = 1.D4*Y(2)
   P(2,1) = .04
   P(2,2) = -1.D4*Y(3) - 6.D7*Y(2)
   P(2,3) = -1.D4*Y(2)
   P(3,1) = 1.
   P(3,2) = 1.
   P(3,3) = 1.
end subroutine dgbydy

subroutine reference()
implicit none
integer :: i
character(len=:),allocatable :: original(:)
write(*,'(a)') '>PROCEDURE dlsodi'
original=[ CHARACTER(LEN=128) :: &
' The output of this program (on a CDC-7600 in single precision)',&
' is as follows:',&
'',&
'   At t =  4.0000e-01   Y =  9.851726e-01  3.386406e-05  1.479357e-02',&
'   At t =  4.0000e+00   Y =  9.055142e-01  2.240418e-05  9.446344e-02',&
'   At t =  4.0000e+01   Y =  7.158050e-01  9.184616e-06  2.841858e-01',&
'   At t =  4.0000e+02   Y =  4.504846e-01  3.222434e-06  5.495122e-01',&
'   At t =  4.0000e+03   Y =  1.831701e-01  8.940379e-07  8.168290e-01',&
'   At t =  4.0000e+04   Y =  3.897016e-02  1.621193e-07  9.610297e-01',&
'   At t =  4.0000e+05   Y =  4.935213e-03  1.983756e-08  9.950648e-01',&
'   At t =  4.0000e+06   Y =  5.159269e-04  2.064759e-09  9.994841e-01',&
'   At t =  4.0000e+07   Y =  5.306413e-05  2.122677e-10  9.999469e-01',&
'   At t =  4.0000e+08   Y =  5.494532e-06  2.197826e-11  9.999945e-01',&
'   At t =  4.0000e+09   Y =  5.129457e-07  2.051784e-12  9.999995e-01',&
'   At t =  4.0000e+10   Y = -7.170472e-08 -2.868188e-13  1.000000e+00',&
'',&
'   No. steps = 330  No. r-s = 404  No. J-s =  69',&
'',&
'']

   write(*,'(a)') (trim(original(i)),i=1,size(original))
   write(*,'(a)') ' >NOW IT PRODUCES'
end subroutine reference
