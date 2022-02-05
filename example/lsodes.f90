!>
!! Demonstration program for the DLSODES package.
!! This is the version of 14 June 2001.
!!
!! This version is in double precision.
!!
!! The package is used for each of the relevant values of mf to solve
!! the problem ydot = A * y, where A is the 9 by 9 sparse matrix
!!```text
!!
!!               -4  1     1
!!                1 -4  1     1
!!                   1 -4        1
!!                        -4  1     1
!!       A =               1 -4  1     1
!!                            1 -4        1
!!                                 -4  1
!!                                  1 -4  1
!!                                     1 -4
!!
!!```
!! The initial conditions are  y(0) = (1, 2, 3, ..., 9).
!! Output is printed at t = 1, 2, and 3.
!! Each case is solved first with nominal (large) values of lrw and liw,
!! and then with values given by lenrw and leniw (optional outputs)
!! on the first run, as a check on these computed work array lengths.
!! If the errors are too large, or other difficulty occurs,
!! a warning message is printed.
!! All output is on unit lout, which is data-loaded to 6 below.
!-----------------------------------------------------------------------
program lsodes
use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) :: atol , erm , ero , hu , rtol , t , tout
integer :: i , igrid , iopt , iout , irun , istate , itask ,      &
               & itol , j , k , l , leniw , lenrw , liw , lrw , m ,     &
               & meth , mf , miter , moss , neq , nerr , nfe , nfea ,   &
               & ngp , nje , nlu , nnz , nout , nqu , nst , nzl , nzu
integer , dimension(10) :: ia
integer , dimension(90) :: iwork
integer , dimension(50) :: ja
integer , save :: lout
real(kind=dp) , dimension(1000) :: rwork
real(kind=dp) , dimension(9) :: y
external fdem , jdem
external ssout
external edit
!
equivalence (ia(1),iwork(31))
equivalence (ja(1),iwork(41))
data lout/6/
!
! Write heading and set fixed parameters.
      write (lout,                                                      &
            &"(/,'Demonstration problem for the DLSODES package',/,/)")
      nerr = 0
      igrid = 3
      neq = igrid**2
      t = 0.0d0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-5
      itask = 1
      iopt = 0
      do i = 1 , neq
         y(i) = i
      enddo
      ia(1) = 1
      k = 1
      do m = 1 , igrid
         do l = 1 , igrid
            j = l + (m-1)*igrid
            if ( m>1 ) then
               ja(k) = j - igrid
               k = k + 1
            endif
            if ( l>1 ) then
               ja(k) = j - 1
               k = k + 1
            endif
            ja(k) = j
            k = k + 1
            if ( l<igrid ) then
               ja(k) = j + 1
               k = k + 1
            endif
            ia(j+1) = k
         enddo
      enddo
      write (lout,99001) neq , t , rtol , atol , (y(i),i=1,neq)
99001 format (' neq =',i4,5x,'t0 =',f4.1,5x,'rtol =',d12.3,5x,'atol =', &
            & d12.3//' Initial y vector =  ',9f5.1)
!
! Loop over all relevant values of mf.
      do moss = 0 , 2
         do meth = 1 , 2
            do miter = 0 , 3
               if ( .not.((miter==0 .or. miter==3) .and. moss/=0) ) then
                  mf = 100*moss + 10*meth + miter
                  write (lout,99008)
         ! First run: nominal work array lengths, 3 output points.
                  irun = 1
                  lrw = 1000
                  liw = 90
                  nout = 3
                  do
                     write (lout,99002) mf , lrw , liw
99002                format (//'Run with mf =',i4,'.',5x,               &
                            &'Input work lengths lrw, liw =',2i6/)
                     do i = 1 , neq
                        y(i) = i
                     enddo
                     t = 0.0d0
                     tout = 1.0d0
                     istate = 1
                     ero = 0.0d0
         ! Loop over output points.  Do output and accuracy check at each.
                     do iout = 1 , nout
                        call dlsodes(fdem,[neq],y,t,tout,itol,[rtol],[atol],  &
                                   & itask,istate,iopt,rwork,lrw,iwork, &
                                   & liw,jdem,mf)
                        nst = iwork(11)
                        hu = rwork(11)
                        nqu = iwork(14)
                        call edit(y,iout,erm)
                        write (lout,99003) t , nst , hu , nqu , erm ,   &
                             & (y(i),i=1,neq)
99003                   format ('At t =',f5.1,3x,'nst =',i4,3x,'hu =',  &
                              & d12.3,3x,'nqu =',i3,3x,' max. err. =',  &
                              & d11.3/'  y array =    ',4d15.6/5d15.6)
                        if ( istate<0 ) exit
                        erm = erm/atol
                        ero = max(ero,erm)
                        if ( erm>100.0d0 ) then
                           write (lout,99004)
99004                      format (//                                   &
                              &' Warning: error exceeds 100 * tolerance'&
                             & //)
                           nerr = nerr + 1
                        endif
                        tout = tout + 1.0d0
                     enddo
                     if ( istate<0 ) nerr = nerr + 1
                     if ( irun==2 ) exit
         ! Print final statistics (first run only)
                     nst = iwork(11)
                     nfe = iwork(12)
                     nje = iwork(13)
                     lenrw = iwork(17)
                     leniw = iwork(18)
                     nnz = iwork(19)
                     ngp = iwork(20)
                     nlu = iwork(21)
                     nzl = iwork(25)
                     nzu = iwork(26)
                     nfea = nfe
                     if ( miter==2 ) nfea = nfe - ngp*nje
                     if ( miter==3 ) nfea = nfe - nje
                     write (lout,99005) lenrw , leniw , nst , nfe ,     &
                                      & nfea , nje , ero
99005                format (/'Final statistics for this run:'/         &
                            &' rwork size =',i4,'   iwork size =',      &
                            &i4/' number of steps =',                   &
                            &i5/' number of f-s   =',                   &
                            &i5/' (excluding J-s) =',                   &
                            &i5/' number of J-s   =',                   &
                            &i5/' error overrun =',d10.2)
                     if ( miter==1 .or. miter==2 ) write (lout,99006)   &
                        & nnz , ngp , nlu , nzl , nzu
99006                format (' number of nonzeros in J = ',             &
                            &i5/' number of J index groups =',          &
                            &i5/' number of LU decomp-s    =',          &
                            &i5/' nonzeros in strict lower factor =',   &
                            &i5/' nonzeros in strict upper factor =',i5)
                     if ( istate<0 ) exit
                     if ( miter==1 .or. miter==2 )                      &
                        & call ssout(neq,rwork(21),iwork,lout)
         ! Return for second run: minimal work array lengths, 1 output point.
                     irun = irun + 1
                     lrw = lenrw
                     liw = leniw
                     nout = 1
                  enddo
               endif
            enddo
         enddo
      enddo
!
      write (lout,99008)
      write (lout,99007) nerr
99007 format (//'Number of errors encountered =',i3)
99008 format (//80('*'))
end program lsodes

subroutine fdem(neq,t,y,ydot)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: y
real(kind=dp) , dimension(neq) :: ydot
!
integer :: i , j , l , m
integer , save :: igrid
!
data igrid/3/
      do i = 1 , neq
         ydot(i) = 0.0d0
      enddo
      do m = 1 , igrid
         do l = 1 , igrid
            j = l + (m-1)*igrid
            if ( m/=1 ) ydot(j-igrid) = ydot(j-igrid) + y(j)
            if ( l/=1 ) ydot(j-1) = ydot(j-1) + y(j)
            ydot(j) = ydot(j) - 4.0d0*y(j)
            if ( l/=igrid ) ydot(j+1) = ydot(j+1) + y(j)
         enddo
      enddo
end subroutine fdem

subroutine jdem(neq,t,y,j,ia,ja,pdj)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
integer :: neq
real(kind=dp) :: t
real(kind=dp) , dimension(neq) :: y
integer :: j
integer , dimension(*) :: ia
integer , dimension(*) :: ja
real(kind=dp) , dimension(neq) :: pdj
!
integer , save :: igrid
integer :: l , m
!
data igrid/3/
      m = (j-1)/igrid + 1
      l = j - (m-1)*igrid
      pdj(j) = -4.0d0
      if ( m/=1 ) pdj(j-igrid) = 1.0d0
      if ( l/=1 ) pdj(j-1) = 1.0d0
      if ( l/=igrid ) pdj(j+1) = 1.0d0
end subroutine jdem

subroutine edit(y,iout,erm)
implicit none
integer,parameter :: dp=kind(0.0d0)
!
real(kind=dp) , dimension(*) :: y
integer :: iout
real(kind=dp) :: erm
!
real(kind=dp) :: er
integer :: i
integer , save :: neq
real(kind=dp) , dimension(9,3) , save :: yex
!
data neq/9/
data yex/6.687279d-01 , 9.901910d-01 , 7.603061d-01 ,           &
 & 8.077979d-01 , 1.170226e+00 , 8.810605d-01 , 5.013331d-01 ,  &
 & 7.201389d-01 , 5.379644d-01 , 1.340488d-01 , 1.917157d-01 ,  &
 & 1.374034d-01 , 1.007882d-01 , 1.437868d-01 , 1.028010d-01 ,  &
 & 3.844343d-02 , 5.477593d-02 , 3.911435d-02 , 1.929166d-02 ,  &
 & 2.735444d-02 , 1.939611d-02 , 1.055981d-02 , 1.496753d-02 ,  &
 & 1.060897d-02 , 2.913689d-03 , 4.128975d-03 , 2.925977d-03/
      erm = 0.0d0
      do i = 1 , neq
         er = abs(y(i)-yex(i,iout))
         erm = max(erm,er)
      enddo
      end subroutine edit

subroutine ssout(neq,iwk,iwork,lout)
implicit none
!
integer :: neq
integer , dimension(*) :: iwk
integer , dimension(*) :: iwork
integer :: lout
!
integer :: i , i1 , i2 , ipian , ipjan , nnz
!
      ipian = iwork(23)
      ipjan = iwork(24)
      nnz = iwork(19)
      i1 = ipian
      i2 = i1 + neq
      write (lout,99001) (iwk(i),i=i1,i2)
99001 format (/' structure descriptor array ian ='/(20i4))
      i1 = ipjan
      i2 = i1 + nnz - 1
      write (lout,99002) (iwk(i),i=i1,i2)
99002 format (/' structure descriptor array jan ='/(20i4))
end subroutine ssout
