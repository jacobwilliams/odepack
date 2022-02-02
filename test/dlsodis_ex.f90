module c_test1
implicit none
integer,parameter,private :: dp=kind(0.0d0)
real(kind=dp),public      :: EODsq, R4D
integer,public            :: NM1
end module c_test1

program dlsodis_ex
use m_odepack
use c_test1
implicit none
integer,parameter              ::  dp=kind(0.0d0)
external                       ::  addasp
external                       ::  jacsp
external                       ::  resid
real(kind=dp),save             ::  atol,rtol
real(kind=dp)                  ::  delx,t,tout
integer                        ::  i,io,istate,nnzlu
integer,save                   ::  iopt,itask,itol,liw,lrw,mf,neq
integer,dimension(30)          ::  iw
real(kind=dp),dimension(1409)  ::  rw
real(kind=dp),dimension(40)    ::  y,ydoti

data itol/1/,rtol/1.0D-3/,atol/1.0D-3/,itask/1/,iopt/0/
data neq/40/,lrw/1409/,liw/30/,mf/121/

   call reference()

   delx = 2.0/neq
   R4D = 0.25/delx
   EODsq = 0.05/delx**2
   NM1 = neq - 1
   do i = 1,neq
      y(i) = 0.0
   enddo
   y(11) = 0.5
   do i = 12,30
      y(i) = 1.0
   enddo
   y(31) = 0.5
   t = 0.0
   tout = 0.1
   istate = 0
   do io = 1,4
      call dlsodis(resid,addasp,jacsp,[neq],y,ydoti,t,tout,itol,[rtol],    &
                 & [atol],itask,istate,iopt,rw,lrw,iw,liw,mf)
      write (6,99010) t,iw(11),rw(11)
   99010 format (' At t =',f5.2,'   No. steps =',i4,'    Last step =',     &
               & d12.4)
      if ( istate/=2 ) then
         write (6,99020) istate
   99020 format (///' Error halt.. ISTATE =',i3)
         stop 1
      else
         tout = tout + 0.1
      endif
   enddo
   write (6,99030) (y(i),i=1,neq)
   99030 format (/' Final solution values..'/8(5D12.4/))
   write (6,99040) iw(17),iw(18),iw(11),iw(12),iw(13)
   99040 format (/' Required RW size =',i5,'   IW size =',                 &
               & i4/' No. steps =',i4,'   No. r-s =',i4,'   No. J-s =',i4)
   nnzlu = iw(25) + iw(26) + neq
   write (6,99050) iw(19),nnzlu
   99050 format (' No. of nonzeros in P matrix =',i4,                      &
                &'   No. of nonzeros in LU =',i4)

end program dlsodis_ex

subroutine gfun(N,T,Y,G)
use c_test1
implicit none
integer,parameter                       ::  dp=kind(0.0d0)
integer,intent(in)                      ::  N
real(kind=dp)                           ::  T
real(kind=dp),intent(in),dimension(N)   ::  Y
real(kind=dp),intent(out),dimension(N)  ::  G
integer                                 ::  i
 
   G(1) = R4D*(Y(N)**2-Y(2)**2) + EODsq*(Y(2)-2.0*Y(1)+Y(N))
   do i = 2,NM1
      G(i) = R4D*(Y(i-1)**2-Y(i+1)**2) + EODsq*(Y(i+1)-2.0*Y(i)+Y(i-1))
   enddo
   G(N) = R4D*(Y(NM1)**2-Y(1)**2) + EODsq*(Y(1)-2.0*Y(N)+Y(NM1))
end subroutine gfun

subroutine resid(N,T,Y,S,R,Ires)
use c_test1
implicit none
integer,parameter                         :: dp=kind(0.0d0)
external                                  :: gfun
 
integer                                   :: N
real(kind=dp)                             :: T
real(kind=dp),dimension(N)                :: Y
real(kind=dp),intent(in),dimension(N)     :: S
real(kind=dp),intent(inout),dimension(N)  :: R
integer                                   :: Ires
 
integer                                   :: i
 
   call gfun(N,T,Y,R)
   R(1) = R(1) - (S(N)+4.0*S(1)+S(2))/6.0
   do i = 2,NM1
      R(i) = R(i) - (S(i-1)+4.0*S(i)+S(i+1))/6.0
   enddo
   R(N) = R(N) - (S(NM1)+4.0*S(N)+S(1))/6.0
end subroutine resid

subroutine addasp(N,T,Y,J,Ip,Jp,P)
implicit none
integer,parameter                         ::  dp=kind(0.0d0)
integer,intent(in)                        ::  N
real(kind=dp)                             ::  T
real(kind=dp),dimension(N)                ::  Y
integer,intent(in)                        ::  J
integer,dimension(*)                      ::  Ip
integer,dimension(*)                      ::  Jp
real(kind=dp),intent(inout),dimension(N)  ::  P
integer                                   ::  jm1,jp1
 
   jm1 = J - 1
   jp1 = J + 1
   if ( J==N ) jp1 = 1
   if ( J==1 ) jm1 = N
   P(J) = P(J) + (2.0/3.0)
   P(jp1) = P(jp1) + (1.0/6.0)
   P(jm1) = P(jm1) + (1.0/6.0)
end subroutine addasp

subroutine jacsp(N,T,Y,S,J,Ip,Jp,Pdj)
use c_test1
implicit none
integer,parameter                       ::  dp=kind(0.0d0)
 
integer,intent(in)                      ::  N
real(kind=dp)                           ::  T
real(kind=dp),intent(in),dimension(N)   ::  Y
real(kind=dp),dimension(N)              ::  S
integer,intent(in)                      ::  J
integer,dimension(*)                    ::  Ip
integer,dimension(*)                    ::  Jp
real(kind=dp),intent(out),dimension(N)  ::  Pdj
 
integer                                 ::  jm1,jp1
 
   jm1 = J - 1
   jp1 = J + 1
   if ( J==1 ) jm1 = N
   if ( J==N ) jp1 = 1
   Pdj(jm1) = -2.0*R4D*Y(J) + EODsq
   Pdj(J) = -2.0*EODsq
   Pdj(jp1) = 2.0*R4D*Y(J) + EODsq
end subroutine jacsp

subroutine reference()
implicit none
integer :: i
character(len=:),allocatable :: original(:)
write(*,'(a)') '>PROCEDURE dlsodis'
original=[ CHARACTER(LEN=128) :: &
' The output of this program (on a CDC-7600 in single precision)',&
' is as follows:',&
'',&
' At t = 0.10   No. steps =  15    Last step =  1.6863e-02',&
' At t = 0.20   No. steps =  19    Last step =  2.4101e-02',&
' At t = 0.30   No. steps =  22    Last step =  4.3143e-02',&
' At t = 0.40   No. steps =  24    Last step =  5.7819e-02',&
'',&
' Final solution values..',&
'  1.8371e-02  1.3578e-02  1.5864e-02  2.3805e-02  3.7245e-02',&
'  5.6630e-02  8.2538e-02  1.1538e-01  1.5522e-01  2.0172e-01',&
'  2.5414e-01  3.1150e-01  3.7259e-01  4.3608e-01  5.0060e-01',&
'  5.6482e-01  6.2751e-01  6.8758e-01  7.4415e-01  7.9646e-01',&
'  8.4363e-01  8.8462e-01  9.1853e-01  9.4500e-01  9.6433e-01',&
'  9.7730e-01  9.8464e-01  9.8645e-01  9.8138e-01  9.6584e-01',&
'  9.3336e-01  8.7497e-01  7.8213e-01  6.5315e-01  4.9997e-01',&
'  3.4672e-01  2.1758e-01  1.2461e-01  6.6208e-02  3.3784e-02',&
'',&
' Required RW size = 1409   IW size =  30',&
' No. steps =  24   No. r-s =  33   No. J-s =   8',&
' No. of nonzeros in P matrix = 120   No. of nonzeros in LU = 194',&
'',&
'']

   write(*,'(a)') (trim(original(i)),i=1,size(original))
   write(*,'(a)') ' >NOW IT PRODUCES'
end subroutine reference
