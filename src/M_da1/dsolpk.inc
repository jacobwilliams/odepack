!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This routine interfaces to one of DSPIOM, DSPIGMR, DPCG, DPCGS, or
!! DUSOL, for the solution of the linear system arising from a Newton
!! iteration.  It is called if MITER .ne. 0.
!! In addition to variables described elsewhere,
!! communication with DSOLPK uses the following variables:
!!
!! WM
!!
!! : real work space containing data for the algorithm
!! (Krylov basis vectors, Hessenberg matrix, etc.)
!! IWM
!! : integer work space containing data for the algorithm
!!
!! X
!!
!! : the right-hand side vector on input, and the solution vector
!! on output, of length N.
!!
!! IERSL
!!
!! : output flag (in Common):
!!       ERSL =  0 means no trouble occurred.
!!       ERSL =  1 means the iterative method failed to converge.
!!                 If the preconditioner is out of date, the step
!!                 is repeated with a new preconditioner.
!!                 Otherwise, the stepsize is reduced (forcing a
!!                 new evaluation of the preconditioner) and the
!!                 step is repeated.
!!       ERSL = -1 means there was a nonrecoverable error in the
!!                 iterative solver, and an error exit occurs.
!!
!! This routine also uses the Common variables TN, EL0, H, N, MITER,
!! DELT, EPCON, SQRTN, RSQRTN, MAXL, KMP, MNEWT, NNI, NLI, NPS, NCFL,
!! LOCWP, LOCIWP.
!-----------------------------------------------------------------------
subroutine dsolpk(Neq,Y,Savf,X,Ewt,Wm,Iwm,f,psol)

integer        :: Neq(*)
real(kind=dp)  :: Y(*)
real(kind=dp)  :: Savf(*)
real(kind=dp)  :: X(*)
real(kind=dp)  :: Ewt(*)
real(kind=dp)  :: Wm(*)
integer        :: Iwm(*)
external       :: f
external       :: psol

real(kind=dp) :: delta, hl0
integer :: iflag, lb, ldl, lgmr, lhes, liom, lp, lpcg, lq, lr, lv, lw, lwk, lz, maxlp1, npsl

   dls1%iersl = 0
   hl0 = dls1%h*dls1%el0
   delta = dlpk%delt*dlpk%epcon
   select case (dls1%miter)
   case (2)
   !-----------------------------------------------------------------------
   !  Use the SPIGMR algorithm to solve the linear system P*x = -f.
   !-----------------------------------------------------------------------
      maxlp1 = dlpk%maxl + 1
      lv = 1
      lb = lv + dls1%n*dlpk%maxl
      lhes = lb + dls1%n + 1
      lq = lhes + dlpk%maxl*maxlp1
      lwk = lq + 2*dlpk%maxl
      ldl = lwk + min(1,dlpk%maxl-dlpk%kmp)*dls1%n
      !X!call dcopy(dls1%n,X,1,Wm(lb),1)
      Wm(lb:lb+dls1%n-1) = X(1:dls1%n)!X!
      call dscal(dls1%n,dlpk%rsqrtn,Ewt,1)
      call dspigmr(Neq,dls1%tn,Y,Savf,Wm(lb),Ewt,dls1%n,dlpk%maxl,maxlp1, &
       & dlpk%kmp,delta,hl0,dlpk%jpre,dlpk%mnewt,f,psol,npsl,X,   &
       & Wm(lv),Wm(lhes),Wm(lq),lgmr,Wm(dlpk%locwp),Iwm(dlpk%lociwp),Wm(lwk),Wm(ldl),iflag)
      dlpk%nni = dlpk%nni + 1
      dlpk%nli = dlpk%nli + lgmr
      dlpk%nps = dlpk%nps + npsl
      call dscal(dls1%n,dlpk%sqrtn,Ewt,1)
      if ( iflag/=0 ) dlpk%ncfl = dlpk%ncfl + 1
      if ( iflag>=2 ) dls1%iersl = 1
      if ( iflag<0 ) dls1%iersl = -1
   case (3)
   !-----------------------------------------------------------------------
   !  Use DPCG to solve the linear system P*x = -f
   !-----------------------------------------------------------------------
      lr = 1
      lp = lr + dls1%n
      lw = lp + dls1%n
      lz = lw + dls1%n
      lwk = lz + dls1%n
      !X!call dcopy(dls1%n,X,1,Wm(lr),1)
      Wm(lr:lr+dls1%n-1) = X(1:dls1%n)!X!
      call dpcg(Neq,dls1%tn,Y,Savf,Wm(lr),Ewt,dls1%n,dlpk%maxl,delta,hl0, &
       & dlpk%jpre,dlpk%mnewt,f,psol,npsl,X,Wm(lp),Wm(lw),Wm(lz), &
       & lpcg,Wm(dlpk%locwp),Iwm(dlpk%lociwp),Wm(lwk),iflag)
      dlpk%nni = dlpk%nni + 1
      dlpk%nli = dlpk%nli + lpcg
      dlpk%nps = dlpk%nps + npsl
      if ( iflag/=0 ) dlpk%ncfl = dlpk%ncfl + 1
      if ( iflag>=2 ) dls1%iersl = 1
      if ( iflag<0 ) dls1%iersl = -1
   case (4)
   !-----------------------------------------------------------------------
   !  Use DPCGS to solve the linear system P*x = -f
   !-----------------------------------------------------------------------
      lr = 1
      lp = lr + dls1%n
      lw = lp + dls1%n
      lz = lw + dls1%n
      lwk = lz + dls1%n
      !X!call dcopy(dls1%n,X,1,Wm(lr),1)
      Wm(lr:lr+dls1%n-1) = X(1:dls1%n)!X!
      call dpcgs(Neq,dls1%tn,Y,Savf,Wm(lr),Ewt,dls1%n,dlpk%maxl,delta,hl0, &
       & dlpk%jpre,dlpk%mnewt,f,psol,npsl,X,Wm(lp),Wm(lw),Wm(lz),&
       & lpcg,Wm(dlpk%locwp),Iwm(dlpk%lociwp),Wm(lwk),iflag)
      dlpk%nni = dlpk%nni + 1
      dlpk%nli = dlpk%nli + lpcg
      dlpk%nps = dlpk%nps + npsl
      if ( iflag/=0 ) dlpk%ncfl = dlpk%ncfl + 1
      if ( iflag>=2 ) dls1%iersl = 1
      if ( iflag<0 ) dls1%iersl = -1
   case (5,6,7,8,9)
   !-----------------------------------------------------------------------
   !  Use DUSOL, which interfaces to PSOL, to solve the linear system
   !  (no Krylov iteration).
   !-----------------------------------------------------------------------
      lb = 1
      lwk = lb + dls1%n
      !X!call dcopy(dls1%n,X,1,Wm(lb),1)
      Wm(lb:lb+dls1%n-1) = X(1:dls1%n)!X!
      call dusol(Neq,dls1%tn,Y,Savf,Wm(lb),Ewt,dls1%n,delta,hl0,dlpk%mnewt, &
       &psol,npsl,X,Wm(dlpk%locwp),Iwm(dlpk%lociwp),Wm(lwk),iflag)
      dlpk%nni = dlpk%nni + 1
      dlpk%nps = dlpk%nps + npsl
      if ( iflag/=0 ) dlpk%ncfl = dlpk%ncfl + 1
      if ( iflag==3 ) dls1%iersl = 1
      if ( iflag<0 ) dls1%iersl = -1
   case default
   !-----------------------------------------------------------------------
   !  Use the SPIOM algorithm to solve the linear system P*x = -f.
   !-----------------------------------------------------------------------
      lv = 1
      lb = lv + dls1%n*dlpk%maxl
      lhes = lb + dls1%n
      lwk = lhes + dlpk%maxl*dlpk%maxl
      !X!call dcopy(dls1%n,X,1,Wm(lb),1)
      Wm(lb:lb+dls1%n-1) = X(1:dls1%n)!X!
      call dscal(dls1%n,dlpk%rsqrtn,Ewt,1)
      call dspiom(Neq,dls1%tn,Y,Savf,Wm(lb),Ewt,dls1%n,dlpk%maxl,dlpk%kmp, &
       & delta,hl0,dlpk%jpre,dlpk%mnewt,f,psol,npsl,X,Wm(lv),       &
       & Wm(lhes),Iwm,liom,Wm(dlpk%locwp),Iwm(dlpk%lociwp),Wm(lwk),iflag)
      dlpk%nni = dlpk%nni + 1
      dlpk%nli = dlpk%nli + liom
      dlpk%nps = dlpk%nps + npsl
      call dscal(dls1%n,dlpk%sqrtn,Ewt,1)
      if ( iflag/=0 ) dlpk%ncfl = dlpk%ncfl + 1
      if ( iflag>=2 ) dls1%iersl = 1
      if ( iflag<0 ) dls1%iersl = -1

   endselect

end subroutine dsolpk
