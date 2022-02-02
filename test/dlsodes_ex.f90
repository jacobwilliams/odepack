program dlsodes_ex
use m_odepack
implicit none
external fex
external jex
 
integer,parameter  ::  dp=kind(0.0d0)
real(kind=dp)      ::  atol,rtol,t,tout
integer            ::  i,iopt,iout,istate,itask,itol,leniw,lenrw,  &
                       & mf,neq,nfe,nje,nlu,nnz,nnzlu,nst
integer,dimension(30)         ::  iwork
integer,save                  ::  liw,lrw
real(kind=dp),dimension(500)  ::  rwork
real(kind=dp),dimension(12)   ::  y
 
data lrw/500/,liw/30/

   call reference()

   neq = 12
   do i = 1,neq
      y(i) = 0.0D0
   enddo
   y(1) = 1.0D0
   t = 0.0D0
   tout = 0.1D0
   itol = 1
   rtol = 1.0D-4
   atol = 1.0D-6
   itask = 1
   istate = 1
   iopt = 0
   mf = 121
   do iout = 1,5
      call dlsodes(fex,[neq],y,t,tout,itol,[rtol],[atol],itask,istate,iopt,&
                 & rwork,lrw,iwork,liw,jex,mf)
      write (6,99010) t,iwork(11),rwork(11),(y(i),i=1,neq)
   99010 format (//' At t =',d11.3,4x,' No. steps =',i5,4x,' Last step =', &
               & d11.3/'  Y array =  ',4D14.5/13x,4D14.5/13x,4D14.5)
      if ( istate<0 ) then
         write (6,99020) istate
   99020 format (///' Error halt.. ISTATE =',i3)
         stop 1
      else
         tout = tout*10.0D0
      endif
   enddo
   lenrw = iwork(17)
   leniw = iwork(18)
   nst = iwork(11)
   nfe = iwork(12)
   nje = iwork(13)
   nlu = iwork(21)
   nnz = iwork(19)
   nnzlu = iwork(25) + iwork(26) + neq
   write (6,99030) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
   99030 format (//' Required RWORK size =',i4,'   IWORK size =',          &
                &i4/' No. steps =',i4,'   No. f-s =',i4,'   No. J-s =',i4, &
                &'   No. LU-s =',i4/' No. of nonzeros in J =',i5,          &
                &'   No. of nonzeros in LU =',i5)
         
end program dlsodes_ex
 
subroutine fex(Neq,T,Y,Ydot)
implicit none
 
integer,parameter :: dp=kind(0.0d0)
integer                                  ::  Neq
real(kind=dp)                            ::  T
real(kind=dp),intent(in),dimension(12)   ::  Y
real(kind=dp),intent(out),dimension(12)  ::  Ydot
 
real(kind=dp),save :: rk1,rk10,rk11,rk12,rk13,rk14,rk15, &
                     & rk16,rk17,rk2,rk3,rk4,rk5,rk6,rk7 ,&
                     & rk8,rk9
real,save          :: rk18,rk19,rk20
 
data rk1/0.1D0/,rk2/10.0D0/,rk3/50.0D0/,rk4/2.5D0/,rk5/0.1D0/,&
   & rk6/10.0D0/,rk7/50.0D0/,rk8/2.5D0/,rk9/50.0D0/,            &
   & rk10/5.0D0/,rk11/50.0D0/,rk12/50.0D0/,rk13/50.0D0/,        &
   & rk14/30.0D0/,rk15/100.0D0/,rk16/2.5D0/,rk17/100.0D0/,      &
   & rk18/2.5D0/,rk19/50.0D0/,rk20/50.0D0/
   Ydot(1) = -rk1*Y(1)
   Ydot(2) = rk1*Y(1) + rk11*rk14*Y(4) + rk19*rk14*Y(5) - rk3*Y(2)*Y(3)    &
           & - rk15*Y(2)*Y(12) - rk2*Y(2)
   Ydot(3) = rk2*Y(2) - rk5*Y(3) - rk3*Y(2)*Y(3) - rk7*Y(10)*Y(3)          &
           & + rk11*rk14*Y(4) + rk12*rk14*Y(6)
   Ydot(4) = rk3*Y(2)*Y(3) - rk11*rk14*Y(4) - rk4*Y(4)
   Ydot(5) = rk15*Y(2)*Y(12) - rk19*rk14*Y(5) - rk16*Y(5)
   Ydot(6) = rk7*Y(10)*Y(3) - rk12*rk14*Y(6) - rk8*Y(6)
   Ydot(7) = rk17*Y(10)*Y(12) - rk20*rk14*Y(7) - rk18*Y(7)
   Ydot(8) = rk9*Y(10) - rk13*rk14*Y(8) - rk10*Y(8)
   Ydot(9) = rk4*Y(4) + rk16*Y(5) + rk8*Y(6) + rk18*Y(7)
   Ydot(10) = rk5*Y(3) + rk12*rk14*Y(6) + rk20*rk14*Y(7) + rk13*rk14*Y(8)  &
            & - rk7*Y(10)*Y(3) - rk17*Y(10)*Y(12) - rk6*Y(10) - rk9*Y(10)
   Ydot(11) = rk10*Y(8)
   Ydot(12) = rk6*Y(10) + rk19*rk14*Y(5) + rk20*rk14*Y(7) - rk15*Y(2)*Y(12)&
            & - rk17*Y(10)*Y(12)
end subroutine fex
 
subroutine jex(Neq,T,Y,J,Ia,Ja,Pdj)
implicit none
!
integer,parameter                        ::  dp=kind(0.0d0)
integer                                  ::  Neq
real(kind=dp)                            ::  T
real(kind=dp),intent(in),dimension(12)   ::  Y
integer,intent(in)                       ::  J
integer,dimension(*)                     ::  Ia
integer,dimension(*)                     ::  Ja
real(kind=dp),intent(out),dimension(12)  ::  Pdj
!
real(kind=dp),save :: rk1,rk10,rk11,rk12,rk13,rk14,rk15, &
                     & rk16,rk17,rk2,rk3,rk4,rk5,rk6,rk7 ,&
                     & rk8,rk9
real,save :: rk18,rk19,rk20
!
data rk1/0.1D0/,rk2/10.0D0/,rk3/50.0D0/,rk4/2.5D0/,rk5/0.1D0/,&
   & rk6/10.0D0/,rk7/50.0D0/,rk8/2.5D0/,rk9/50.0D0/,            &
   & rk10/5.0D0/,rk11/50.0D0/,rk12/50.0D0/,rk13/50.0D0/,        &
   & rk14/30.0D0/,rk15/100.0D0/,rk16/2.5D0/,rk17/100.0D0/,      &
   & rk18/2.5D0/,rk19/50.0D0/,rk20/50.0D0/
   select case (J)
   case (2)
      Pdj(2) = -rk3*Y(3) - rk15*Y(12) - rk2
      Pdj(3) = rk2 - rk3*Y(3)
      Pdj(4) = rk3*Y(3)
      Pdj(5) = rk15*Y(12)
      Pdj(12) = -rk15*Y(12)
   case (3)
      Pdj(2) = -rk3*Y(2)
      Pdj(3) = -rk5 - rk3*Y(2) - rk7*Y(10)
      Pdj(4) = rk3*Y(2)
      Pdj(6) = rk7*Y(10)
      Pdj(10) = rk5 - rk7*Y(10)
   case (4)
      Pdj(2) = rk11*rk14
      Pdj(3) = rk11*rk14
      Pdj(4) = -rk11*rk14 - rk4
      Pdj(9) = rk4
   case (5)
      Pdj(2) = rk19*rk14
      Pdj(5) = -rk19*rk14 - rk16
      Pdj(9) = rk16
      Pdj(12) = rk19*rk14
   case (6)
      Pdj(3) = rk12*rk14
      Pdj(6) = -rk12*rk14 - rk8
      Pdj(9) = rk8
      Pdj(10) = rk12*rk14
   case (7)
      Pdj(7) = -rk20*rk14 - rk18
      Pdj(9) = rk18
      Pdj(10) = rk20*rk14
      Pdj(12) = rk20*rk14
   case (8)
      Pdj(8) = -rk13*rk14 - rk10
      Pdj(10) = rk13*rk14
      Pdj(11) = rk10
   case (9)
   case (10)
      Pdj(3) = -rk7*Y(3)
      Pdj(6) = rk7*Y(3)
      Pdj(7) = rk17*Y(12)
      Pdj(8) = rk9
      Pdj(10) = -rk7*Y(3) - rk17*Y(12) - rk6 - rk9
      Pdj(12) = rk6 - rk17*Y(12)
   case (11)
   case (12)
      Pdj(2) = -rk15*Y(2)
      Pdj(5) = rk15*Y(2)
      Pdj(7) = rk17*Y(10)
      Pdj(10) = -rk17*Y(10)
      Pdj(12) = -rk15*Y(2) - rk17*Y(10)
   case default
      Pdj(1) = -rk1
      Pdj(2) = rk1
   endselect
   
end subroutine jex

subroutine reference()
implicit none
integer :: i
character(len=:),allocatable :: original(:)
write(*,'(a)') '>PROCEDURE dlsodes'
original=[ CHARACTER(LEN=128) :: &
' The output of this program (on a Cray-1 in single precision)',&
' is as follows:',&
'',&
'',&
' At t =  1.000e-01     No. steps =   12     Last step =  1.515e-02',&
'  Y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07',&
'                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07',&
'                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06',&
'',&
'',&
' At t =  1.000e+00     No. steps =   33     Last step =  7.880e-02',&
'  Y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05',&
'                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05',&
'                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03',&
'',&
'',&
' At t =  1.000e+01     No. steps =   48     Last step =  1.239e+00',&
'  Y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05',&
'                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04',&
'                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01',&
'',&
'',&
' At t =  1.000e+02     No. steps =   91     Last step =  3.764e+00',&
'  Y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11',&
'                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07',&
'                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01',&
'',&
'',&
' At t =  1.000e+03     No. steps =  111     Last step =  4.156e+02',&
'  Y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14',&
'               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15',&
'                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01',&
'',&
'',&
' Required RWORK size = 442   IWORK size =  30',&
' No. steps = 111   No. f-s = 142   No. J-s =   2   No. LU-s =  20',&
' No. of nonzeros in J =   44   No. of nonzeros in LU =   50',&
'']

   write(*,'(a)') (trim(original(i)),i=1,size(original))
   write(*,'(a)') ' >NOW IT PRODUCES'
end subroutine reference
