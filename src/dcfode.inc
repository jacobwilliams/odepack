!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### NAME
!!   dcfode(3f) - [M_odepack] Set ODE integrator coefficients.
!!
!!### SYNOPSIS
!!         subroutine dcfode (meth, elco, tesco)
!!         integer          :: meth
!!         double precision :: elco(13,12)
!!         double precision :: tesco(3,12)
!!
!!### DESCRIPTION
!!   DCFODE is called by the integrator routine to set coefficients
!!   needed there.
!!
!!   DCFODE is called once at the beginning of the problem,
!!   and is not called again unless and until METH is changed.
!!
!!### OPTIONS
!!
!!   METH
!!
!!   : The coefficients for the current method, as
!!   given by the value of METH, are set for all orders and saved.
!!   The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
!!   (A smaller value of the maximum order is also allowed.)
!!
!!   ELCO
!!
!!   : The ELCO array contains the basic method coefficients.
!!   The coefficients el(i), 1 .le. i .le. nq+1, for the method of
!!   order nq are stored in ELCO(i,nq).  They are given by a genetrating
!!   polynomial, i.e.,
!!       l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
!!   For the implicit Adams methods, l(x) is given by
!!       dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
!!   For the BDF methods, l(x) is given by
!!       l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
!!   where
!!       K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!!
!!   TESCO
!!
!!   : The TESCO array contains test constants used for the
!!   local error test and the selection of step size and/or order.
!!   At order nq, TESCO(k,nq) is used for the selection of step
!!   size at order nq - 1 if k = 1, at order nq if k = 2, and at order
!!   nq + 1 if k = 3.
!!
! ### BEGIN PROLOGUE  DCFODE
! ### SUBSIDIARY
! ### PURPOSE  Set ODE integrator coefficients.
! ### TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### DESCRIPTION
! ### SEE ALSO  DLSODE
! ### ROUTINES CALLED  (NONE)
! ### REVISION HISTORY  (YYMMDD)
!     791129  DATE WRITTEN
!     890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!     890503  Minor cosmetic changes.  (FNF)
!     930809  Renamed to allow single/double precision versions. (ACH)
! ### END PROLOGUE  DCFODE
!-----------------------------------------------------------------------
subroutine dcfode (meth, elco, tesco)
integer,intent(in)          :: meth
real(kind=dp),intent(inout) :: elco(13,12)
real(kind=dp),intent(out)   :: tesco(3,12)

integer          :: i, ib, nq, nqm1, nqp1
double precision :: agamq, fnq, fnqm1, pc(12), pint, ragq, rqfac, rq1fac, tsign, xpin

   select case(meth)
   case(1)
      elco(1,1) = 1.0d0
      elco(2,1) = 1.0d0
      tesco(1,1) = 0.0d0
      tesco(2,1) = 2.0d0
      tesco(1,2) = 1.0d0
      tesco(3,12) = 0.0d0
      pc(1) = 1.0d0
      rqfac = 1.0d0
      do nq = 2,12
         !-----------------------------------------------------------------------
         ! The PC array will contain the coefficients of the polynomial
         !     p(x) = (x+1)*(x+2)*...*(x+nq-1).
         ! Initially, p(x) = 1.
         !-----------------------------------------------------------------------
         rq1fac = rqfac
         rqfac = rqfac/nq
         nqm1 = nq - 1
         fnqm1 = nqm1
         nqp1 = nq + 1
         ! Form coefficients of p(x)*(x+nq-1). ----------------------------------
         pc(nq) = 0.0d0
         do ib = 1,nqm1
           i = nqp1 - ib
           pc(i) = pc(i-1) + fnqm1*pc(i)
         enddo
         pc(1) = fnqm1*pc(1)
         ! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
         pint = pc(1)
         xpin = pc(1)/2.0d0
         tsign = 1.0d0
         do i = 2,nq
           tsign = -tsign
           pint = pint + tsign*pc(i)/i
           xpin = xpin + tsign*pc(i)/(i+1)
         enddo
         ! Store coefficients in ELCO and TESCO. --------------------------------
         elco(1,nq) = pint*rq1fac
         elco(2,nq) = 1.0d0
         do i = 2,nq
           elco(i+1,nq) = rq1fac*pc(i)/i
         enddo
         agamq = rqfac*xpin
         ragq = 1.0d0/agamq
         tesco(2,nq) = ragq
         if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/nqp1
         tesco(3,nqm1) = ragq
      enddo

   case(2)
      pc(1) = 1.0d0
      rq1fac = 1.0d0
      do nq = 1,5
         !-----------------------------------------------------------------------
         ! The PC array will contain the coefficients of the polynomial
         !     p(x) = (x+1)*(x+2)*...*(x+nq).
         ! Initially, p(x) = 1.
         !-----------------------------------------------------------------------
         fnq = nq
         nqp1 = nq + 1
         ! Form coefficients of p(x)*(x+nq). ------------------------------------
         pc(nqp1) = 0.0d0
         do ib = 1,nq
            i = nq + 2 - ib
            pc(i) = pc(i-1) + fnq*pc(i)
         enddo
         pc(1) = fnq*pc(1)
         ! Store coefficients in ELCO and TESCO. --------------------------------
         do i = 1,nqp1
            elco(i,nq) = pc(i)/pc(2)
         enddo
         elco(2,nq) = 1.0d0
         tesco(1,nq) = rq1fac
         tesco(2,nq) = nqp1/elco(1,nq)
         tesco(3,nq) = (nq+2)/elco(1,nq)
         rq1fac = rq1fac/fnq
      enddo
    case default
       stop '*dcfode* unknown value for meth'
    end select
end subroutine dcfode
