program dlsoibt_ex
use m_odepack
implicit none
external addabt
external jacbt
external resid
 
integer,parameter :: dp=kind(0.0d0)
real(kind=dp) :: atol,rtol,t,tout
integer :: i,io,iopt,istate,itask,itol,liw,lrw,mf,neq
integer,dimension(61) :: iwork
real(kind=dp),dimension(514) :: rwork
real(kind=dp),dimension(41) :: y,ydoti

   call reference()

   neq = 41
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
   itol = 1
   rtol = 1.0D-4
   atol = 1.0D-5
   itask = 1
   istate = 0
   iopt = 0
   lrw = 514
   liw = 61
   iwork(1) = 1
   iwork(2) = neq
   mf = 21
   do io = 1,4
      call dlsoibt(resid,addabt,jacbt,[neq],y,ydoti,t,tout,itol,[rtol],    &
                 & [atol],itask,istate,iopt,rwork,lrw,iwork,liw,mf)
      write (6,99010) t,iwork(11),iwork(12),iwork(13)
   99010 format (' At t =',f5.2,'   No. steps =',i4,'  No. r-s =',i4,      &
                &'  No. J-s =',i3)
      if ( istate/=2 ) then
         write (6,99020) istate
   99020 format (///' Error halt.. ISTATE =',i3)
         stop 1
      else
         tout = tout + 0.1
      endif
   enddo
   write (6,99030) (y(i),i=1,neq)
   99030 format (/' Final solution values..'/9(5D12.4/))

end program dlsoibt_ex
 
subroutine resid(N,T,Y,S,R,Ires)
implicit none
integer,parameter                       ::  dp=kind(0.0d0)
 
integer,intent(in)                      ::  N
real(kind=dp)                           ::  T
real(kind=dp),intent(in),dimension(N)   ::  Y
real(kind=dp),intent(in),dimension(N)   ::  S
real(kind=dp),intent(out),dimension(N)  ::  R
integer                                 ::  Ires
 
real(kind=dp),save                      ::  delx,eta
real(kind=dp)                           ::  eodsq
integer                                 ::  i,nm1
 
data eta/0.05/,delx/0.05/
   eodsq = eta/delx**2
   R(1) = eodsq*(Y(3)-2.0*Y(2)+Y(1)) - S(1)
   nm1 = N - 1
   do i = 2,nm1
      R(i) = (Y(i-1)**2-Y(i+1)**2)/(4.0*delx)                              &
           & + eodsq*(Y(i+1)-2.0*Y(i)+Y(i-1)) - (S(i-1)+4.0*S(i)+S(i+1))   &
           & /6.0
   enddo
   R(N) = eodsq*(Y(N-2)-2.0*Y(nm1)+Y(N)) - S(N)
end subroutine resid
 
subroutine addabt(N,T,Y,Mb,Nb,Pa,Pb,Pc)
implicit                                         none
integer,parameter                                ::    dp=kind(0.0d0)
 
integer,intent(in)                               ::    N
real(kind=dp)                                    ::    T
real(kind=dp),dimension(N)                       ::    Y
integer,intent(in)                               ::    Mb
integer,intent(in)                               ::    Nb
real(kind=dp),intent(inout),dimension(Mb,Mb,Nb)  ::    Pa
real(kind=dp),intent(inout),dimension(Mb,Mb,Nb)  ::    Pb
real(kind=dp),intent(inout),dimension(Mb,Mb,Nb)  ::    Pc
 
integer                                          ::    k,nm1
 
   Pa(1,1,1) = Pa(1,1,1) + 1.0
   nm1 = N - 1
   do k = 2,nm1
      Pa(1,1,k) = Pa(1,1,k) + (4.0/6.0)
      Pb(1,1,k) = Pb(1,1,k) + (1.0/6.0)
      Pc(1,1,k) = Pc(1,1,k) + (1.0/6.0)
   enddo
   Pa(1,1,N) = Pa(1,1,N) + 1.0
end subroutine addabt
 
subroutine jacbt(N,T,Y,S,Mb,Nb,Pa,Pb,Pc)
implicit none
integer,parameter                              ::    dp=kind(0.0d0)
integer,intent(in)                             ::    N
real(kind=dp)                                  ::    T
real(kind=dp),intent(in),dimension(N)          ::    Y
real(kind=dp),dimension(N)                     ::    S
integer,intent(in)                             ::    Mb
integer,intent(in)                             ::    Nb
real(kind=dp),intent(out),dimension(Mb,Mb,Nb)  ::    Pa
real(kind=dp),intent(out),dimension(Mb,Mb,Nb)  ::    Pb
real(kind=dp),intent(out),dimension(Mb,Mb,Nb)  ::    Pc
real(kind=dp),save                             ::    delx,eta
real(kind=dp)                                  ::    eodsq
integer                                        ::    k
 
data eta/0.05/,delx/0.05/
   eodsq = eta/delx**2
   Pa(1,1,1) = eodsq
   Pb(1,1,1) = -2.0*eodsq
   Pc(1,1,1) = eodsq
   do k = 2,N
      Pa(1,1,k) = -2.0*eodsq
      Pb(1,1,k) = -Y(k+1)*(0.5/delx) + eodsq
      Pc(1,1,k) = Y(k-1)*(0.5/delx) + eodsq
   enddo
   Pb(1,1,N) = eodsq
   Pc(1,1,N) = -2.0*eodsq
   Pa(1,1,N) = eodsq
end subroutine jacbt

subroutine reference()
implicit none
integer :: i
character(len=:),allocatable :: original(:)
write(*,'(a)') '>PROCEDURE dlsoibt'
original=[ CHARACTER(LEN=128) :: &
' The output of this program (on a CDC-7600 in single precision)',&
' is as follows:',&
'',&
' At t = 0.10   No. steps =  35  No. r-s =  45  No. J-s =  9',&
' At t = 0.20   No. steps =  43  No. r-s =  54  No. J-s = 10',&
' At t = 0.30   No. steps =  48  No. r-s =  60  No. J-s = 11',&
' At t = 0.40   No. steps =  51  No. r-s =  64  No. J-s = 12',&
'',&
' Final solution values..',&
'  1.2747e-02  1.1997e-02  1.5560e-02  2.3767e-02  3.7224e-02',&
'  5.6646e-02  8.2645e-02  1.1557e-01  1.5541e-01  2.0177e-01',&
'  2.5397e-01  3.1104e-01  3.7189e-01  4.3530e-01  5.0000e-01',&
'  5.6472e-01  6.2816e-01  6.8903e-01  7.4612e-01  7.9829e-01',&
'  8.4460e-01  8.8438e-01  9.1727e-01  9.4330e-01  9.6281e-01',&
'  9.7632e-01  9.8426e-01  9.8648e-01  9.8162e-01  9.6617e-01',&
'  9.3374e-01  8.7535e-01  7.8236e-01  6.5321e-01  5.0003e-01',&
'  3.4709e-01  2.1876e-01  1.2771e-01  7.3671e-02  5.0642e-02',&
'  5.4496e-02',&
'',&
'']

   write(*,'(a)') (trim(original(i)),i=1,size(original))
   write(*,'(a)') ' >NOW IT PRODUCES'
end subroutine reference
