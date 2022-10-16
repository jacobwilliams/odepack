!==================================================================================================================================!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!==================================================================================================================================!
!>
!!### Synopsis
!!
!! DLSOIBT: Livermore Solver for Ordinary differential equations given
!!          in Implicit form, with Block-Tridiagonal Jacobian treatment.
!!
!! DLSOIBT solves the initial value problem for linearly implicit
!! systems of first order ODEs,
!!```text
!!     A(t,y) * dy/dt = g(t,y),  where A(t,y) is a square matrix,
!!```
!! or, in component form,
!!```text
!!     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
!!        i,1      1                     i,NEQ      NEQ
!!
!!      =   g ( t, y, y ,..., y    )   ( i = 1,...,NEQ )
!!           i      1   2       NEQ
!!```
!! If A is singular, this is a differential-algebraic system.
!!
!! DLSOIBT is a variant version of the DLSODI package, for the case where
!! the matrices A, dg/dy, and d(A*s)/dy are all block-tridiagonal.
!!
!! This version is in double precision.
!!
!!-----------------------------------------------------------------------
!!### Summary of Usage.
!!
!! Communication between the user and the DLSOIBT package, for normal
!! situations, is summarized here.  This summary describes only a subset
!! of the full set of options available.  See the full description for
!! details, including optional communication, nonstandard options,
!! and instructions for special situations.  See also the example
!! problem (with program and output) following this summary.
!!
!! **A.** First, provide a subroutine of the form:
!!```fortran
!!               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
!!               DOUBLE PRECISION T, Y(*), S(*), R(*)
!!```
!! which computes the residual function
!!```text
!!     r = g(t,y)  -  A(t,y) * s
!!```
!! as a function of t and the vectors y and s.  (s is an internally
!! generated approximation to dy/dt.)  The arrays Y and S are inputs
!! to the RES routine and should not be altered.  The residual
!! vector is to be stored in the array R.  The argument IRES should be
!! ignored for casual use of DLSOIBT.  (For uses of IRES, see the
!! paragraph on RES in the full description below.)
!!
!! **B.** Next, identify the block structure of the matrices A = A(t,y) and
!! dr/dy.  DLSOIBT must deal internally with a linear combination, P, of
!! these two matrices.  The matrix P (hence both A and dr/dy) must have
!! a block-tridiagonal form with fixed structure parameters
!!
!!        MB = block size, MB .ge. 1, and
!!        NB = number of blocks in each direction, NB .ge. 4,
!!
!! with MB\*NB = NEQ.  In each of the NB block-rows of the matrix P
!! (each consisting of MB consecutive rows), the nonzero elements are
!! to lie in three consecutive MB by MB blocks.  In block-rows
!! 2 through NB - 1, these are centered about the main diagonal.
!!
!! in block-rows 1 and NB, they are the diagonal blocks and the two
!! blocks adjacent to the diagonal block.  (Thus block positions (1,3)
!! and (NB,NB-2) can be nonzero.)
!!
!! Alternatively, P (hence A and dr/dy) may be only approximately
!! equal to matrices with this form, and DLSOIBT should still succeed.
!! The block-tridiagonal matrix P is described by three arrays,
!! each of size MB by MB by NB:
!!
!!         PA = array of diagonal blocks,
!!         PB = array of superdiagonal (and one subdiagonal) blocks, and
!!         PC = array of subdiagonal (and one superdiagonal) blocks.
!!
!! Specifically, the three MB by MB blocks in the k-th block-row of P
!! are stored in (reading across):
!!
!!         PC(*,*,k) = block to the left of the diagonal block,
!!         PA(*,*,k) = diagonal block, and
!!         PB(*,*,k) = block to the right of the diagonal block,
!!
!! except for k = 1, where the three blocks (reading across) are
!!
!!         PA(*,*,1) (= diagonal block), PB(*,*,1), and PC(*,*,1),
!!
!! and k = NB, where they are
!!
!!         PB(*,*,NB), PC(*,*,NB), and PA(*,*,NB) (= diagonal block).
!!
!! (Each asterisk \* stands for an index that ranges from 1 to MB.)
!!
!! **C.** You must also provide a subroutine of the form:
!!```fortran
!!     SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
!!     DOUBLE PRECISION T, Y(*), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
!!```
!! which adds the nonzero blocks of the matrix A = A(t,y) to the
!! contents of the arrays PA, PB, and PC, following the structure
!! description in Paragraph B above.
!! T and the Y array are input and should not be altered.
!! Thus the affect of ADDA should be the following:
!!```fortran
!!     DO K = 1,NB
!!       DO J = 1,MB
!!         DO I = 1,MB
!!           PA(I,J,K) = PA(I,J,K) +
!!             ( (I,J) element of K-th diagonal block of A)
!!           PB(I,J,K) = PB(I,J,K) +
!!             ( (I,J) element of block in block position (K,K+1) of A,
!!             or in block position (NB,NB-2) if K = NB)
!!           PC(I,J,K) = PC(I,J,K) +
!!             ( (I,J) element of block in block position (K,K-1) of A,
!!             or in block position (1,3) if K = 1)
!!         ENDDO
!!       ENDDO
!!     ENDDO
!!```
!!
!! **D.** For the sake of efficiency, you are encouraged to supply the
!! Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)\*s
!! (s = a fixed vector) as above.  If dr/dy is being supplied,
!! use MF = 21, and provide a subroutine of the form:
!!```fortran
!!     SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
!!     DOUBLE PRECISION T, Y(*), S(*), PA(MB,MB,NB), PB(MB,MB,NB), &
!!    &                 PC(MB,MB,NB)
!!```
!! which computes dr/dy as a function of t, y, and s.  Here T, Y, and
!! S are inputs, and the routine is to load dr/dy into PA, PB, PC,
!! according to the structure description in Paragraph B above.
!!
!! That is, load the diagonal blocks into PA, the superdiagonal blocks
!! (and block (NB,NB-2) ) into PB, and the subdiagonal blocks (and
!! block (1,3) ) into PC.  The blocks in block-row k of dr/dy are to
!! be loaded into PA(\*,\*,k), PB(\*,\*,k), and PC(\*,\*,k).
!!
!! Only nonzero elements need be loaded, and the indexing
!! of PA, PB, and PC is the same as in the ADDA routine.
!!
!! Note that if A is independent of Y (or this dependence
!! is weak enough to be ignored) then JAC is to compute dg/dy.
!!
!! If it is not feasible to provide a JAC routine, use
!! MF = 22, and DLSOIBT will compute an approximate Jacobian
!! internally by difference quotients.
!!
!! **E.** Next decide whether or not to provide the initial value of the
!! derivative vector dy/dt.  If the initial value of A(t,y) is
!! nonsingular (and not too ill-conditioned), you may let DLSOIBT compute
!! this vector (ISTATE = 0).  (DLSOIBT will solve the system A\*s = g for
!! s, with initial values of A and g.)  If A(t,y) is initially
!! singular, then the system is a differential-algebraic system, and
!! you must make use of the particular form of the system to compute the
!! initial values of y and dy/dt.  In that case, use ISTATE = 1 and
!! load the initial value of dy/dt into the array YDOTI.
!!
!! The input array YDOTI and the initial Y array must be consistent with
!! the equations A\*dy/dt = g.  This implies that the initial residual
!! r = g(t,y) - A(t,y)\*YDOTI  must be approximately zero.
!!
!! **F.** Write a main program which calls Subroutine DLSOIBT once for
!! each point at which answers are desired.  This should also provide
!! for possible use of logical unit 6 for output of error messages by
!! DLSOIBT.  on the first call to DLSOIBT, supply arguments as follows:
!!
!! RES
!!
!! : name of user subroutine for residual function r.
!!
!! ADDA
!!
!! : name of user subroutine for computing and adding A(t,y).
!!
!! JAC
!!
!! : name of user subroutine for Jacobian matrix dr/dy
!! (MF = 21).  If not used, pass a dummy name.
!!
!! Note: the names for the RES and ADDA routines and (if used) the
!! JAC routine must be declared External in the calling program.
!!
!! NEQ
!!
!! : number of scalar equations in the system.
!!
!! Y
!!
!! : array of initial values, of length NEQ.
!!
!! YDOTI
!!
!! : array of length NEQ (containing initial dy/dt if ISTATE = 1).
!!
!! T
!!
!! : the initial value of the independent variable.
!!
!! TOUT
!!
!! : first point where output is desired (.ne. T).
!!
!! ITOL
!!
!! : 1 or 2 according as ATOL (below) is a scalar or array.
!!
!! RTOL
!!
!! : relative tolerance parameter (scalar).
!!
!! ATOL
!!
!! : absolute tolerance parameter (scalar or array).
!! the estimated local error in y(i) will be controlled so as
!! to be roughly less (in magnitude) than
!!```text
!!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!!```
!! Thus the local error test passes if, in each component,
!! either the absolute error is less than ATOL (or ATOL(i)),
!! or the relative error is less than RTOL.
!! Use RTOL = 0.0 for pure absolute error control, and
!! use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!! control.  Caution: Actual (global) errors may exceed these
!! local tolerances, so choose them conservatively.
!!
!! ITASK
!!
!! : 1 for normal computation of output values of y at t = TOUT.
!!
!! ISTATE
!!
!! : integer flag (input and output).  Set ISTATE = 1 if the
!! initial dy/dt is supplied, and 0 otherwise.
!!
!! IOPT
!!
!! : 0 to indicate no optional inputs used.
!!
!! RWORK
!!
!! : real work array of length at least:
!!
!!             22 + 9*NEQ + 3*MB*MB*NB        for MF = 21 or 22.
!!
!! LRW
!!
!! : declared length of RWORK (in user's dimension).
!!
!! IWORK
!!
!! : integer work array of length at least 20 + NEQ.
!! Input in IWORK(1) the block size MB and in IWORK(2) the
!! number NB of blocks in each direction along the matrix A.
!! These must satisfy  MB .ge. 1, NB .ge. 4, and MB\*NB = NEQ.
!!
!! LIW
!!
!! : declared length of IWORK (in user's dimension).
!!
!! MF
!!
!! : method flag.  Standard values are:
!!
!!          21 for a user-supplied Jacobian.
!!          22 for an internally generated Jacobian.
!!
!!  For other choices of MF, see the paragraph on MF in
!!  the full description below.
!!
!! Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
!! and possibly ATOL.
!!
!! **G.** The output from the first call (or any call) is:
!!
!! Y
!!
!! : array of computed values of y(t) vector.
!!
!! T
!!
!! : corresponding value of independent variable (normally TOUT).
!!
!! ISTATE
!!
!! : ISTATE values:
!!
!!  value | description
!!  ----- | -----------------------------------------------------------
!!     2  | if DLSOIBT was successful, negative otherwise.
!!    -1  | means excess work done on this call (check all inputs).
!!    -2  | means excess accuracy requested (tolerances too small).
!!    -3  | means illegal input detected (see printed message).
!!    -4  | means repeated error test failures (check all inputs).
!!    -5  | means repeated convergence failures (perhaps bad Jacobian
!!        | supplied or wrong choice of tolerances).
!!    -6  | means error weight became zero during problem. (Solution
!!        | component i vanished, and ATOL or ATOL(i) = 0.)
!!    -7  | cannot occur in casual use.
!!    -8  | means DLSOIBT was unable to compute the initial dy/dt.
!!        | In casual use, this means A(t,y) is initially singular.
!!        | Supply YDOTI and use ISTATE = 1 on the first call.
!!
!! If DLSOIBT returns ISTATE = -1, -4, or -5, then the output of
!! DLSOIBT also includes YDOTI = array containing residual vector
!! r = g - A \* dy/dt  evaluated at the current t, y, and dy/dt.
!!
!! **H.** To continue the integration after a successful return, simply
!! reset TOUT and call DLSOIBT again.  No other parameters need be reset.
!!
!!-----------------------------------------------------------------------
!!### Example Problem.
!!
!! The following is an example problem, with the coding needed
!! for its solution by DLSOIBT.  The problem comes from the partial
!! differential equation (the Burgers equation)
!!```text
!!   du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
!! on -1 .le. x .le. 1.  The boundary conditions are
!!   du/dx = 0  at x = -1 and at x = 1.
!! The initial profile is a square wave,
!!   u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.
!! The PDE is discretized in x by a simplified Galerkin method,
!! using piecewise linear basis functions, on a grid of 40 intervals.
!! The equations at x = -1 and 1 use a 3-point difference approximation
!! for the right-hand side.  The result is a system A * dy/dt = g(y),
!! of size NEQ = 41, where y(i) is the approximation to u at x = x(i),
!! with x(i) = -1 + (i-1)*delx, delx = 2/(NEQ-1) = .05.  The individual
!! equations in the system are
!!   dy(1)/dt = ( y(3) - 2*y(2) + y(1) ) * eta / delx**2,
!!   dy(NEQ)/dt = ( y(NEQ-2) - 2*y(NEQ-1) + y(NEQ) ) * eta / delx**2,
!! and for i = 2, 3, ..., NEQ-1,
!!   (1/6) dy(i-1)/dt + (4/6) dy(i)/dt + (1/6) dy(i+1)/dt
!!       = ( y(i-1)**2 - y(i+1)**2 ) / (4*delx)
!!         + ( y(i+1) - 2*y(i) + y(i-1) ) * eta / delx**2.
!!```
!! The following coding solves the problem with MF = 21, with output
!! of solution statistics at t = .1, .2, .3, and .4, and of the
!! solution vector at t = .4.  Here the block size is just MB = 1.
!!
!!```fortran
!!program dlsoibt_ex
!!use m_odepack
!!implicit none
!!external addabt
!!external jacbt
!!external resid
!!
!!integer,parameter :: dp=kind(0.0d0)
!!real(kind=dp) :: atol,rtol,t,tout
!!integer :: i,io,iopt,istate,itask,itol,liw,lrw,mf,neq
!!integer,dimension(61) :: iwork
!!real(kind=dp),dimension(514) :: rwork
!!real(kind=dp),dimension(41) :: y,ydoti
!!
!!   neq = 41
!!   do i = 1,neq
!!      y(i) = 0.0
!!   enddo
!!   y(11) = 0.5
!!   do i = 12,30
!!      y(i) = 1.0
!!   enddo
!!   y(31) = 0.5
!!   t = 0.0
!!   tout = 0.1
!!   itol = 1
!!   rtol = 1.0D-4
!!   atol = 1.0D-5
!!   itask = 1
!!   istate = 0
!!   iopt = 0
!!   lrw = 514
!!   liw = 61
!!   iwork(1) = 1
!!   iwork(2) = neq
!!   mf = 21
!!   do io = 1,4
!!      call dlsoibt(resid,addabt,jacbt,[neq],y,ydoti,t,tout,itol,[rtol],    &
!!                 & [atol],itask,istate,iopt,rwork,lrw,iwork,liw,mf)
!!      write (6,99010) t,iwork(11),iwork(12),iwork(13)
!!   99010 format (' At t =',f5.2,'   No. steps =',i4,'  No. r-s =',i4,      &
!!                &'  No. J-s =',i3)
!!      if ( istate/=2 ) then
!!         write (6,99020) istate
!!   99020 format (///' Error halt.. ISTATE =',i3)
!!         stop 1
!!      else
!!         tout = tout + 0.1
!!      endif
!!   enddo
!!   write (6,99030) (y(i),i=1,neq)
!!   99030 format (/' Final solution values..'/9(5D12.4/))
!!
!!end program dlsoibt_ex
!!
!!subroutine resid(N,T,Y,S,R,Ires)
!!implicit none
!!integer,parameter                       ::  dp=kind(0.0d0)
!!
!!integer,intent(in)                      ::  N
!!real(kind=dp)                           ::  T
!!real(kind=dp),intent(in),dimension(N)   ::  Y
!!real(kind=dp),intent(in),dimension(N)   ::  S
!!real(kind=dp),intent(out),dimension(N)  ::  R
!!integer                                 ::  Ires
!!
!!real(kind=dp),save                      ::  delx,eta
!!real(kind=dp)                           ::  eodsq
!!integer                                 ::  i,nm1
!!
!!data eta/0.05/,delx/0.05/
!!   eodsq = eta/delx**2
!!   R(1) = eodsq*(Y(3)-2.0*Y(2)+Y(1)) - S(1)
!!   nm1 = N - 1
!!   do i = 2,nm1
!!      R(i) = (Y(i-1)**2-Y(i+1)**2)/(4.0*delx)                              &
!!           & + eodsq*(Y(i+1)-2.0*Y(i)+Y(i-1)) - (S(i-1)+4.0*S(i)+S(i+1))   &
!!           & /6.0
!!   enddo
!!   R(N) = eodsq*(Y(N-2)-2.0*Y(nm1)+Y(N)) - S(N)
!!end subroutine resid
!!
!!subroutine addabt(N,T,Y,Mb,Nb,Pa,Pb,Pc)
!!implicit                                         none
!!integer,parameter                                ::    dp=kind(0.0d0)
!!
!!integer,intent(in)                               ::    N
!!real(kind=dp)                                    ::    T
!!real(kind=dp),dimension(N)                       ::    Y
!!integer,intent(in)                               ::    Mb
!!integer,intent(in)                               ::    Nb
!!real(kind=dp),intent(inout),dimension(Mb,Mb,Nb)  ::    Pa
!!real(kind=dp),intent(inout),dimension(Mb,Mb,Nb)  ::    Pb
!!real(kind=dp),intent(inout),dimension(Mb,Mb,Nb)  ::    Pc
!!
!!integer                                          ::    k,nm1
!!
!!   Pa(1,1,1) = Pa(1,1,1) + 1.0
!!   nm1 = N - 1
!!   do k = 2,nm1
!!      Pa(1,1,k) = Pa(1,1,k) + (4.0/6.0)
!!      Pb(1,1,k) = Pb(1,1,k) + (1.0/6.0)
!!      Pc(1,1,k) = Pc(1,1,k) + (1.0/6.0)
!!   enddo
!!   Pa(1,1,N) = Pa(1,1,N) + 1.0
!!end subroutine addabt
!!
!!subroutine jacbt(N,T,Y,S,Mb,Nb,Pa,Pb,Pc)
!!implicit none
!!integer,parameter                              ::    dp=kind(0.0d0)
!!integer,intent(in)                             ::    N
!!real(kind=dp)                                  ::    T
!!real(kind=dp),intent(in),dimension(N)          ::    Y
!!real(kind=dp),dimension(N)                     ::    S
!!integer,intent(in)                             ::    Mb
!!integer,intent(in)                             ::    Nb
!!real(kind=dp),intent(out),dimension(Mb,Mb,Nb)  ::    Pa
!!real(kind=dp),intent(out),dimension(Mb,Mb,Nb)  ::    Pb
!!real(kind=dp),intent(out),dimension(Mb,Mb,Nb)  ::    Pc
!!real(kind=dp),save                             ::    delx,eta
!!real(kind=dp)                                  ::    eodsq
!!integer                                        ::    k
!!
!!data eta/0.05/,delx/0.05/
!!   eodsq = eta/delx**2
!!   Pa(1,1,1) = eodsq
!!   Pb(1,1,1) = -2.0*eodsq
!!   Pc(1,1,1) = eodsq
!!   do k = 2,N
!!      Pa(1,1,k) = -2.0*eodsq
!!      Pb(1,1,k) = -Y(k+1)*(0.5/delx) + eodsq
!!      Pc(1,1,k) = Y(k-1)*(0.5/delx) + eodsq
!!   enddo
!!   Pb(1,1,N) = eodsq
!!   Pc(1,1,N) = -2.0*eodsq
!!   Pa(1,1,N) = eodsq
!!end subroutine jacbt
!!```
!!
!! The output of this program (on a CDC-7600 in single precision)
!! is as follows:
!!```text
!! At t = 0.10   No. steps =  35  No. r-s =  45  No. J-s =  9
!! At t = 0.20   No. steps =  43  No. r-s =  54  No. J-s = 10
!! At t = 0.30   No. steps =  48  No. r-s =  60  No. J-s = 11
!! At t = 0.40   No. steps =  51  No. r-s =  64  No. J-s = 12
!!
!! Final solution values..
!!  1.2747e-02  1.1997e-02  1.5560e-02  2.3767e-02  3.7224e-02
!!  5.6646e-02  8.2645e-02  1.1557e-01  1.5541e-01  2.0177e-01
!!  2.5397e-01  3.1104e-01  3.7189e-01  4.3530e-01  5.0000e-01
!!  5.6472e-01  6.2816e-01  6.8903e-01  7.4612e-01  7.9829e-01
!!  8.4460e-01  8.8438e-01  9.1727e-01  9.4330e-01  9.6281e-01
!!  9.7632e-01  9.8426e-01  9.8648e-01  9.8162e-01  9.6617e-01
!!  9.3374e-01  8.7535e-01  7.8236e-01  6.5321e-01  5.0003e-01
!!  3.4709e-01  2.1876e-01  1.2771e-01  7.3671e-02  5.0642e-02
!!  5.4496e-02
!!```
!!-----------------------------------------------------------------------
!!#### Full Description of User Interface to DLSOIBT.
!!
!! The user interface to DLSOIBT consists of the following parts.
!!
!! 1.   The call sequence to Subroutine DLSOIBT, which is a driver
!!      routine for the solver.  This includes descriptions of both
!!      the call sequence arguments and of user-supplied routines.
!!      Following these descriptions is a description of
!!      optional inputs available through the call sequence, and then
!!      a description of optional outputs (in the work arrays).
!!
!! 2.   Descriptions of other routines in the DLSOIBT package that may be
!!      (optionally) called by the user.  These provide the ability to
!!      alter error message handling, save and restore the internal
!!      Common, and obtain specified derivatives of the solution y(t).
!!
!! 3.   Descriptions of Common blocks to be declared in overlay
!!      or similar environments, or to be saved when doing an interrupt
!!      of the problem and continued solution later.
!!
!! 4.   Description of two routines in the DLSOIBT package, either of
!!      which the user may replace with his/her own version, if desired.
!!      These relate to the measurement of errors.
!!
!!-----------------------------------------------------------------------
!!#### Part 1.  Call Sequence.
!!
!! The call sequence parameters used for input only are
!!     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
!!     IOPT, LRW, LIW, MF,
!!
!! and those used for both input and output are
!!     Y, T, ISTATE, YDOTI.
!!
!! The work arrays RWORK and IWORK are also used for additional and
!! optional inputs and optional outputs.  (The term output here refers
!! to the return from Subroutine DLSOIBT to the user's calling program.)
!!
!! The legality of input parameters will be thoroughly checked on the
!! initial call for the problem, but not checked thereafter unless a
!! change in input parameters is flagged by ISTATE = 3 on input.
!!
!! The descriptions of the call arguments are as follows.
!!
!! RES
!!
!! : the name of the user-supplied subroutine which supplies
!! the residual vector for the ODE system, defined by
!!```text
!!            r = g(t,y) - A(t,y) * s
!!```
!! as a function of the scalar t and the vectors
!! s and y (s approximates dy/dt).  This subroutine
!! is to have the form
!!```fortran
!!              SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
!!              DOUBLE PRECISION T, Y(*), S(*), R(*)
!!```
!! where NEQ, T, Y, S, and IRES are input, and R and
!! IRES are output. Y, S, and R are arrays of length NEQ.
!!
!! On input, IRES indicates how DLSOIBT will use the
!! returned array R, as follows:
!!
!!             IRES = 1  means that DLSOIBT needs the full residual,
!!                       r = g - A\*s, exactly.
!!             IRES = -1 means that DLSOIBT is using R only to compute
!!                       the Jacobian dr/dy by difference quotients.
!!
!! The RES routine can ignore IRES, or it can omit some terms
!! if IRES = -1.  If A does not depend on y, then RES can
!! just return R = g when IRES = -1.  If g - A\*s contains other
!! additive terms that are independent of y, these can also be
!! dropped, if done consistently, when IRES = -1.
!!
!! The subroutine should set the flag IRES if it
!! encounters a halt condition or illegal input.
!! Otherwise, it should not reset IRES.  On output,
!!
!! IRES = 1 or -1 represents a normal return, and
!! DLSOIBT continues integrating the ODE.  Leave IRES
!! unchanged from its input value.
!!
!! IRES = 2 tells DLSOIBT to immediately return control
!! to the calling program, with ISTATE = 3.  This lets
!! the calling program change parameters of the problem
!! if necessary.
!!
!! IRES = 3 represents an error condition (for example, an
!! illegal value of y).  DLSOIBT tries to integrate the system
!! without getting IRES = 3 from RES.  If it cannot, DLSOIBT
!! returns with ISTATE = -7 or -1.
!!
!! On an DLSOIBT return with ISTATE = 3, -1, or -7, the
!! values of T and Y returned correspond to the last point
!! reached successfully without getting the flag IRES = 2 or 3.
!!
!! The flag values IRES = 2 and 3 should not be used to
!! handle switches or root-stop conditions.  This is better
!! done by calling DLSOIBT in a one-step mode and checking the
!! stopping function for a sign change at each step.
!!
!! If quantities computed in the RES routine are needed
!! externally to DLSOIBT, an extra call to RES should be made
!! for this purpose, for consistent and accurate results.
!! To get the current dy/dt for the S argument, use DINTDY.
!!
!! RES must be declared External in the calling
!! program. See note below for more about RES.
!!
!! ADDA
!!
!! : the name of the user-supplied subroutine which adds the
!! matrix A = A(t,y) to another matrix, P, stored in
!! block-tridiagonal form.  This routine is to have the form
!!```fortran
!!         SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
!!         DOUBLE PRECISION T, Y(*), PA(MB,MB,NB), PB(MB,MB,NB), &
!!        &                 PC(MB,MB,NB)
!!```
!! where NEQ, T, Y, MB, NB, and the arrays PA, PB, and PC
!! are input, and the arrays PA, PB, and PC are output.
!! Y is an array of length NEQ, and the arrays PA, PB, PC
!! are all MB by MB by NB.
!!
!! Here a block-tridiagonal structure is assumed for A(t,y),
!! and also for the matrix P to which A is added here,
!! as described in Paragraph B of the Summary of Usage above.
!! Thus the affect of ADDA should be the following:
!!```fortran
!!         DO K = 1,NB
!!           DO J = 1,MB
!!             DO I = 1,MB
!!               PA(I,J,K) = PA(I,J,K) +
!!                 ( (I,J) element of K-th diagonal block of A)
!!               PB(I,J,K) = PB(I,J,K) +
!!                 ( (I,J) element of block (K,K+1) of A,
!!                 or block (NB,NB-2) if K = NB)
!!               PC(I,J,K) = PC(I,J,K) +
!!                 ( (I,J) element of block (K,K-1) of A,
!!                 or block (1,3) if K = 1)
!!             ENDDO
!!           ENDDO
!!         ENDDO
!!```
!! ADDA must be declared External in the calling program.
!! See note below for more information about ADDA.
!!
!! JAC
!!
!! : the name of the user-supplied subroutine which supplies
!! the Jacobian matrix, dr/dy, where r = g - A\*s.  JAC is
!! required if MITER = 1.  Otherwise a dummy name can be
!! passed.  This subroutine is to have the form
!!```fortran
!!         SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
!!         DOUBLE PRECISION T, Y(*), S(*), PA(MB,MB,NB),
!!        &                 PB(MB,MB,NB), PC(MB,MB,NB)
!!```
!! where NEQ, T, Y, S, MB, NB, and the arrays PA, PB, and PC
!! are input, and the arrays PA, PB, and PC are output.
!! Y and S are arrays of length NEQ, and the arrays PA, PB, PC
!! are all MB by MB by NB.
!!
!! PA, PB, and PC are to be loaded with partial derivatives
!! (elements of the Jacobian matrix) on output, in terms of the
!! block-tridiagonal structure assumed, as described
!! in Paragraph B of the Summary of Usage above.
!!
!! That is, load the diagonal blocks into PA, the
!! superdiagonal blocks (and block (NB,NB-2) ) into PB, and
!! the subdiagonal blocks (and block (1,3) ) into PC.
!!
!! The blocks in block-row k of dr/dy are to be loaded into
!! PA(\*,\*,k), PB(\*,\*,k), and PC(\*,\*,k).
!!
!! Thus the affect of JAC should be the following:
!!
!!```fortran
!!          DO K = 1,NB
!!            DO J = 1,MB
!!              DO I = 1,MB
!!                PA(I,J,K) = ( (I,J) element of
!!                  K-th diagonal block of dr/dy)
!!                PB(I,J,K) = ( (I,J) element of block (K,K+1)
!!                  of dr/dy, or block (NB,NB-2) if K = NB)
!!                PC(I,J,K) = ( (I,J) element of block (K,K-1)
!!                  of dr/dy, or block (1,3) if K = 1)
!!              ENDDO
!!            ENDDO
!!          ENDDO
!!```
!! PA, PB, and PC are preset to zero by the solver,
!! so that only the nonzero elements need be loaded by JAC.
!! Each call to JAC is preceded by a call to RES with the same
!! arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
!! intermediate quantities shared by both calculations may be
!! saved in a user Common block by RES and not recomputed by JAC
!! if desired.  Also, JAC may alter the Y array, if desired.
!!
!! JAC need not provide dr/dy exactly.  A crude
!! approximation will do, so that DLSOIBT may be used when
!! A and dr/dy are not really block-tridiagonal, but are close
!! to matrices that are.
!!
!! JAC must be declared External in the calling program.
!! See note below for more about JAC.
!!
!! Note on RES, ADDA, and JAC:
!!
!! These subroutines may access user-defined quantities in
!! NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!! (dimensioned in the subroutines) and/or Y has length
!! exceeding NEQ(1).  However, these routines should not alter
!! NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
!! See the descriptions of NEQ and Y below.
!!
!! NEQ
!!
!! : the size of the system (number of first order ordinary
!! differential equations or scalar algebraic equations).
!! Used only for input.
!!
!! NEQ may be decreased, but not increased, during the problem.
!! If NEQ is decreased (with ISTATE = 3 on input), the
!! remaining components of Y should be left undisturbed, if
!! these are to be accessed in RES, ADDA, or JAC.
!!
!! Normally, NEQ is a scalar, and it is generally referred to
!! as a scalar in this user interface description.  However,
!! NEQ may be an array, with NEQ(1) set to the system size.
!!
!! (The DLSOIBT package accesses only NEQ(1).)  In either case,
!! this parameter is passed as the NEQ argument in all calls
!! to RES, ADDA, and JAC.  Hence, if it is an array,
!!
!! locations NEQ(2),... may be used to store other integer data
!! and pass it to RES, ADDA, or JAC.  Each such subroutine
!! must include NEQ in a Dimension statement in that case.
!!
!! Y
!!
!! : a real array for the vector of dependent variables, of
!! length NEQ or more.  Used for both input and output on the
!! first call (ISTATE = 0 or 1), and only for output on other
!! calls.  On the first call, Y must contain the vector of
!! initial values.  On output, Y contains the computed solution
!! vector, evaluated at t.  If desired, the Y array may be used
!! for other purposes between calls to the solver.
!!
!! This array is passed as the Y argument in all calls to RES,
!! ADDA, and JAC.  Hence its length may exceed NEQ,
!! and locations Y(NEQ+1),... may be used to store other real
!! data and pass it to RES, ADDA, or JAC.  (The DLSOIBT
!! package accesses only Y(1),...,Y(NEQ). )
!!
!! YDOTI
!!
!! : a real array for the initial value of the vector
!! dy/dt and for work space, of dimension at least NEQ.
!!
!! On input:
!!
!! If ISTATE = 0 then DLSOIBT will compute the initial value
!! of dy/dt, if A is nonsingular.  Thus YDOTI will
!! serve only as work space and may have any value.
!!
!! If ISTATE = 1 then YDOTI must contain the initial value
!! of dy/dt.
!!
!! If ISTATE = 2 or 3 (continuation calls) then YDOTI
!! may have any value.
!!
!! Note: If the initial value of A is singular, then
!! DLSOIBT cannot compute the initial value of dy/dt, so
!! it must be provided in YDOTI, with ISTATE = 1.
!!
!! On output, when DLSOIBT terminates abnormally with ISTATE =
!! -1, -4, or -5, YDOTI will contain the residual
!! r = g(t,y) - A(t,y)\*(dy/dt).  If r is large, t is near
!! its initial value, and YDOTI is supplied with ISTATE = 1,
!! there may have been an incorrect input value of
!! YDOTI = dy/dt, or the problem (as given to DLSOIBT)
!! may not have a solution.
!!
!! If desired, the YDOTI array may be used for other
!! purposes between calls to the solver.
!!
!! T
!!
!! : the independent variable.  On input, T is used only on the
!! first call, as the initial point of the integration.
!!
!! On output, after each call, T is the value at which a
!! computed solution y is evaluated (usually the same as TOUT).
!! On an error return, T is the farthest point reached.
!!
!! TOUT
!!
!! : the next value of t at which a computed solution is desired.
!! Used only for input.
!!
!! When starting the problem (ISTATE = 0 or 1), TOUT may be
!! equal to T for one call, then should .ne. T for the next
!! call.  For the initial T, an input value of TOUT .ne. T is
!! used in order to determine the direction of the integration
!! (i.e. the algebraic sign of the step sizes) and the rough
!! scale of the problem.  Integration in either direction
!! (forward or backward in t) is permitted.
!!
!! If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!! the first call (i.e. the first call with TOUT .ne. T).
!! Otherwise, TOUT is required on every call.
!!
!! If ITASK = 1, 3, or 4, the values of TOUT need not be
!! monotone, but a value of TOUT which backs up is limited
!! to the current internal T interval, whose endpoints are
!! TCUR - HU and TCUR (see optional outputs, below, for
!! TCUR and HU).
!!
!! ITOL
!!
!! : an indicator for the type of error control.  See
!! description below under ATOL.  Used only for input.
!!
!! RTOL
!!
!! : a relative error tolerance parameter, either a scalar or
!! an array of length NEQ.  See description below under ATOL.
!! Input only.
!!
!! ATOL
!!
!! : an absolute error tolerance parameter, either a scalar or
!! an array of length NEQ.  Input only.
!!
!! The input parameters ITOL, RTOL, and ATOL determine
!! the error control performed by the solver.  The solver will
!! control the vector E = (E(i)) of estimated local errors
!! in y, according to an inequality of the form
!!```text
!!                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
!!          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!!```
!! and the RMS-norm (root-mean-square norm) here is
!! RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!! is a vector of weights which must always be positive, and
!! the values of RTOL and ATOL should all be non-negative.
!! The following table gives the types (scalar/array) of
!! RTOL and ATOL, and the corresponding form of EWT(i).
!!
!!  ITOL |  RTOL   |   ATOL   |      EWT(i)
!!  ---- | ------- | -------- | -----------------------------
!!   1   | scalar  |  scalar  |  RTOL\*ABS(Y(i)) + ATOL
!!   2   | scalar  |  array   |  RTOL\*ABS(Y(i)) + ATOL(i)
!!   3   | array   |  scalar  |  RTOL(i)\*ABS(Y(i)) + ATOL
!!   4   | array   |  scalar  |  RTOL(i)\*ABS(Y(i)) + ATOL(i)
!!
!! When either of these parameters is a scalar, it need not
!! be dimensioned in the user's calling program.
!!
!! If none of the above choices (with ITOL, RTOL, and ATOL
!! fixed throughout the problem) is suitable, more general
!! error controls can be obtained by substituting
!! user-supplied routines for the setting of EWT and/or for
!! the norm calculation.  See Part 4 below.
!!
!! If global errors are to be estimated by making a repeated
!! run on the same problem with smaller tolerances, then all
!! components of RTOL and ATOL (i.e. of EWT) should be scaled
!! down uniformly.
!!
!! ITASK
!!
!! : an index specifying the task to be performed.
!! Input only.  ITASK has the following values and meanings.
!!
!!  value | description
!!  ----- | ---------------------------------------------------------
!!      1 | means normal computation of output values of y(t) at
!!        | t = TOUT (by overshooting and interpolating).
!!      2 | means take one step only and return.
!!      3 | means stop at the first internal mesh point at or
!!        | beyond t = TOUT and return.
!!      4 | means normal computation of output values of y(t) at
!!        | t = TOUT but without overshooting t = TCRIT.
!!        | TCRIT must be input as RWORK(1).  TCRIT may be equal to
!!        | or beyond TOUT, but not behind it in the direction of
!!        | integration.  This option is useful if the problem
!!        | has a singularity at or beyond t = TCRIT.
!!      5 | means take one step, without passing TCRIT, and return.
!!        | TCRIT must be input as RWORK(1).
!!
!! **Note**:  If ITASK = 4 or 5 and the solver reaches TCRIT
!! (within roundoff), it will return T = TCRIT (exactly) to
!! indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!! in which case answers at t = TOUT are returned first).
!!
!! ISTATE
!!
!! : an index used for input and output to specify the
!! state of the calculation.
!!
!! On input, the values of ISTATE are as follows.
!!
!!  value | description
!!  ----- | ---------------------------------------------------------
!!      0 | means this is the first call for the problem, and
!!        | DLSOIBT is to compute the initial value of dy/dt
!!        | (while doing other initializations).  See note below.
!!      1 | means this is the first call for the problem, and
!!        | the initial value of dy/dt has been supplied in
!!        | YDOTI (DLSOIBT will do other initializations).
!!        | See note below.
!!      2 | means this is not the first call, and the calculation
!!        | is to continue normally, with no change in any input
!!        | parameters except possibly TOUT and ITASK.
!!        | (If ITOL, RTOL, and/or ATOL are changed between calls
!!        | with ISTATE = 2, the new values will be used but not
!!        | tested for legality.)
!!      3 | means this is not the first call, and the
!!        | calculation is to continue normally, but with
!!        | a change in input parameters other than
!!        | TOUT and ITASK.  Changes are allowed in
!!        | NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, MB, NB,
!!        | and any of the optional inputs except H0.
!!        | (See IWORK description for MB and NB.)
!!
!! Note:  A preliminary call with TOUT = T is not counted
!! as a first call here, as no initialization or checking of
!! input is done.  (Such a call is sometimes useful for the
!! purpose of outputting the initial conditions.)
!! Thus the first call for which TOUT .ne. T requires
!! ISTATE = 0 or 1 on input.
!!
!! On output, ISTATE has the following values and meanings.
!!
!!  value | description
!!  ----- | ---------------------------------------------------------
!!      0 | or 1  means nothing was done; TOUT = t and
!!        | ISTATE = 0 or 1 on input.
!!      2 | means that the integration was performed successfully.
!!      3 | means that the user-supplied Subroutine RES signalled
!!        | DLSOIBT to halt the integration and return (IRES = 2).
!!        | Integration as far as T was achieved with no occurrence
!!        | of IRES = 2, but this flag was set on attempting the
!!        | next step.
!!     -1 | means an excessive amount of work (more than MXSTEP
!!        | steps) was done on this call, before completing the
!!        | requested task, but the integration was otherwise
!!        | successful as far as T.  (MXSTEP is an optional input
!!        | and is normally 500.)  To continue, the user may
!!        | simply reset ISTATE to a value .gt. 1 and call again
!!        | (the excess work step counter will be reset to 0).
!!        | In addition, the user may increase MXSTEP to avoid
!!        | this error return (see below on optional inputs).
!!     -2 | means too much accuracy was requested for the precision
!!        | of the machine being used.  This was detected before
!!        | completing the requested task, but the integration
!!        | was successful as far as T.  To continue, the tolerance
!!        | parameters must be reset, and ISTATE must be set
!!        | to 3.  The optional output TOLSF may be used for this
!!        | purpose.  (Note: If this condition is detected before
!!        | taking any steps, then an illegal input return
!!        | (ISTATE = -3) occurs instead.)
!!     -3 | means illegal input was detected, before taking any
!!        | integration steps.  See written message for details.
!!        | Note:  If the solver detects an infinite loop of calls
!!        | to the solver with illegal input, it will cause
!!        | the run to stop.
!!     -4 | means there were repeated error test failures on
!!        | one attempted step, before completing the requested
!!        | task, but the integration was successful as far as T.
!!        | The problem may have a singularity, or the input
!!        | may be inappropriate.
!!     -5 | means there were repeated convergence test failures on
!!        | one attempted step, before completing the requested
!!        | task, but the integration was successful as far as T.
!!        | This may be caused by an inaccurate Jacobian matrix.
!!     -6 | means EWT(i) became zero for some i during the
!!        | integration.  Pure relative error control (ATOL(i) = 0.0)
!!        | was requested on a variable which has now vanished.
!!        | The integration was successful as far as T.
!!     -7 | means that the user-supplied Subroutine RES set
!!        | its error flag (IRES = 3) despite repeated tries by
!!        | DLSOIBT to avoid that condition.
!!     -8 | means that ISTATE was 0 on input but DLSOIBT was unable
!!        | to compute the initial value of dy/dt.  See the
!!        | printed message for details.
!!
!! Note:  Since the normal output value of ISTATE is 2,
!! it does not need to be reset for normal continuation.
!! Similarly, ISTATE (= 3) need not be reset if RES told
!! DLSOIBT to return because the calling program must change
!! the parameters of the problem.
!!
!! Also, since a negative input value of ISTATE will be
!! regarded as illegal, a negative output value requires the
!! user to change it, and possibly other inputs, before
!! calling the solver again.
!!
!! IOPT
!!
!! : an integer flag to specify whether or not any optional
!! inputs are being used on this call.  Input only.
!! The optional inputs are listed separately below.
!!
!!       IOPT = 0 means no optional inputs are being used.
!!                Default values will be used in all cases.
!!       IOPT = 1 means one or more optional inputs are being used.
!!
!! RWORK
!!
!! : a real working array (double precision).
!! The length of RWORK must be at least
!!```text
!!             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
!!          NYH    = the initial value of NEQ,
!!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!!                   smaller value is given as an optional input),
!!          LENWM  = 3*MB*MB*NB + 2.
!!          (See MF description for the definition of METH.)
!! Thus if MAXORD has its default value and NEQ is constant,
!! this length is
!!             22 + 16*NEQ + 3*MB*MB*NB     for MF = 11 or 12,
!!             22 + 9*NEQ + 3*MB*MB*NB      for MF = 21 or 22.
!!```text
!! The first 20 words of RWORK are reserved for conditional
!! and optional inputs and optional outputs.
!!
!! The following word in RWORK is a conditional input:
!!```text
!!            RWORK(1) = TCRIT = critical value of t which the solver
!!                       is not to overshoot.  Required if ITASK is
!!                       4 or 5, and ignored otherwise.  (See ITASK.)
!!```
!! LRW
!!
!! : the length of the array RWORK, as declared by the user.
!! (This will be checked by the solver.)
!!
!! IWORK
!!
!! : an integer work array.  The length of IWORK must be at least
!! 20 + NEQ .  The first few words of IWORK are used for
!! additional and optional inputs and optional outputs.
!!
!! The following 2 words in IWORK are additional required
!! inputs to DLSOIBT:
!!            IWORK(1) = MB = block size
!!            IWORK(2) = NB = number of blocks in the main diagonal
!! These must satisfy  MB .ge. 1, NB .ge. 4, and MB\*NB = NEQ.
!!
!! LIW
!!
!! : the length of the array IWORK, as declared by the user.
!! (This will be checked by the solver.)
!!
!! Note:  The work arrays must not be altered between calls to DLSOIBT
!! for the same problem, except possibly for the additional and
!! optional inputs, and except for the last 3\*NEQ words of RWORK.
!! The latter space is used for internal scratch space, and so is
!! available for use by the user outside DLSOIBT between calls, if
!! desired (but not for use by RES, ADDA, or JAC).
!!
!! MF
!!
!! : the method flag.  used only for input.  The legal values of
!! MF are 11, 12, 21, and 22.
!!
!! MF has decimal digits METH and MITER: MF = 10\*METH + MITER.
!! METH indicates the basic linear multistep method:
!!
!!  METH | description
!!  ---- | -------------------------------------
!!     1 | means the implicit Adams method.
!!     2 | means the method based on Backward
!!       | Differentiation Formulas (BDFS).
!!
!! The BDF method is strongly preferred for stiff
!! problems, while the Adams method is preferred when the
!! problem is not stiff.  If the matrix A(t,y) is
!! nonsingular, stiffness here can be taken to mean that of
!! the explicit ODE system dy/dt = A-inverse \* g.  If A is
!! singular, the concept of stiffness is not well defined.
!!
!! If you do not know whether the problem is stiff, we
!! recommend using METH = 2.  If it is stiff, the advantage
!! of METH = 2 over METH = 1 will be great, while if it is
!! not stiff, the advantage of METH = 1 will be slight.
!! If maximum efficiency is important, some experimentation
!! with METH may be necessary.
!!
!! MITER indicates the corrector iteration method:
!!
!!  MITER | description
!!  ----- | -------------------------------------
!!      1 | means chord iteration with a user-supplied
!!        | block-tridiagonal Jacobian.
!!        |
!!      2 | means chord iteration with an internally
!!        | generated (difference quotient) block-
!!        | tridiagonal Jacobian approximation, using
!!        | 3\*MB+1 extra calls to RES per dr/dy evaluation.
!!
!! If MITER = 1, the user must supply a Subroutine JAC
!! (the name is arbitrary) as described above under JAC.
!!
!! For MITER = 2, a dummy argument can be used.
!!
!!-----------------------------------------------------------------------
!!#### Optional Inputs.
!!
!! The following is a list of the optional inputs provided for in the
!! call sequence.  (See also Part 2.)  For each such input variable,
!! this table lists its name as used in this documentation, its
!! location in the call sequence, its meaning, and the default value.
!! The use of any of these inputs requires IOPT = 1, and in that
!! case all of these inputs are examined.  A value of zero for any
!! of these optional inputs will cause the default value to be used.
!! Thus to use a subset of the optional inputs, simply preload
!! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
!! then set those of interest to nonzero values.
!!
!!  Name   | Location |     Meaning and Default Value
!!  ------ | -------- | -------------------------------------------------
!!  H0     | RWORK(5) | the step size to be attempted on the first step.
!!         |          | The default value is determined by the solver.
!!         |          |
!!  HMAX   | RWORK(6) | the maximum absolute step size allowed.
!!         |          | The default value is infinite.
!!         |          |
!!  HMIN   | RWORK(7) | the minimum absolute step size allowed.
!!         |          | The default value is 0.  (This lower bound is not
!!         |          | enforced on the final step before reaching TCRIT
!!         |          | when ITASK = 4 or 5.)
!!         |          |
!!  MAXORD | IWORK(5) | the maximum order to be allowed.  The default
!!         |          | value is 12 if METH = 1, and 5 if METH = 2.
!!         |          | If MAXORD exceeds the default value, it will
!!         |          | be reduced to the default value.
!!         |          | If MAXORD is changed during the problem, it may
!!         |          | cause the current order to be reduced.
!!         |          |
!!  MXSTEP | IWORK(6) | maximum number of (internally defined) steps
!!         |          | allowed during one call to the solver.
!!         |          | The default value is 500.
!!         |          |
!!  MXHNIL | IWORK(7) | maximum number of messages printed (per problem)
!!         |          | warning that T + H = T on a step (H = step size).
!!         |          | This must be positive to result in a non-default
!!         |          | value.  The default value is 10.
!!
!!-----------------------------------------------------------------------
!! Optional Outputs.
!!
!! As optional additional output from DLSOIBT, the variables listed
!! below are quantities related to the performance of DLSOIBT
!! which are available to the user.  These are communicated by way of
!! the work arrays, but also have internal mnemonic names as shown.
!! Except where stated otherwise, all of these outputs are defined
!! on any successful return from DLSOIBT, and on any return with
!! ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
!! input) or -8, they will be unchanged from their existing values
!! (if any), except possibly for TOLSF, LENRW, and LENIW.
!! On any error return, outputs relevant to the error will be defined,
!! as noted below.
!!
!!  Name  | Location |     Meaning
!!  ----- | -------- | ---------------------------------------------------
!!  HU    | RWORK(11)| the step size in t last used (successfully).
!!        |          |
!!  HCUR  | RWORK(12)| the step size to be attempted on the next step.
!!        |          |
!!  TCUR  | RWORK(13)| the current value of the independent variable
!!        |          | which the solver has actually reached, i.e. the
!!        |          | current internal mesh point in t.  On output, TCUR
!!        |          | will always be at least as far as the argument
!!        |          | T, but may be farther (if interpolation was done).
!!        |          |
!!  TOLSF | RWORK(14)| a tolerance scale factor, greater than 1.0,
!!        |          | computed when a request for too much accuracy was
!!        |          | detected (ISTATE = -3 if detected at the start of
!!        |          | the problem, ISTATE = -2 otherwise).  If ITOL is
!!        |          | left unaltered but RTOL and ATOL are uniformly
!!        |          | scaled up by a factor of TOLSF for the next call,
!!        |          | then the solver is deemed likely to succeed.
!!        |          | (The user may also ignore TOLSF and alter the
!!        |          | tolerance parameters in any other way appropriate.)
!!        |          |
!!  NST   | IWORK(11)| the number of steps taken for the problem so far.
!!        |          |
!!  NRE   | IWORK(12)| the number of residual evaluations (RES calls)
!!        |          | for the problem so far.
!!        |          |
!!  NJE   | IWORK(13)| the number of Jacobian evaluations (each involving
!!        |          | an evaluation of a and dr/dy) for the problem so
!!        |          | far.  This equals the number of calls to ADDA and
!!        |          | (if MITER = 1) to JAC, and the number of matrix
!!        |          | LU decompositions.
!!        |          |
!!  NQU   | IWORK(14)| the method order last used (successfully).
!!        |          |
!!  NQCUR | IWORK(15)| the order to be attempted on the next step.
!!        |          |
!!  IMXER | IWORK(16)| the index of the component of largest magnitude in
!!        |          | the weighted local error vector ( E(i)/EWT(i) ),
!!        |          | on an error return with ISTATE = -4 or -5.
!!        |          |
!!  LENRW | IWORK(17)| the length of RWORK actually required.
!!        |          | This is defined on normal returns and on an illegal
!!        |          | input return for insufficient storage.
!!        |          |
!!  LENIW | IWORK(18)| the length of IWORK actually required.
!!        |          | This is defined on normal returns and on an illegal
!!        |          | input return for insufficient storage.
!!        |          |
!!
!! The following two arrays are segments of the RWORK array which
!! may also be of interest to the user as optional outputs.
!! For each array, the table below gives its internal name,
!! its base address in RWORK, and its description.
!!
!! Name  | Base Address |    Description
!! ----- | ------------ | ------------------------------------------------
!! YH    | 21           | the Nordsieck history array, of size NYH by
!!       |              | (NQCUR + 1), where NYH is the initial value
!!       |              | of NEQ.  For j = 0,1,...,NQCUR, column j+1
!!       |              | of YH contains HCUR\*\*j/factorial(j) times
!!       |              | the j-th derivative of the interpolating
!!       |              | polynomial currently representing the solution,
!!       |              | evaluated at t = TCUR.
!!       |              |
!! ACOR  |  LENRW-NEQ+1 | array of size NEQ used for the accumulated
!!       |              | corrections on each step, scaled on output to
!!       |              | represent the estimated local error in y on
!!       |              | the last step.  This is the vector E in the
!!       |              | description of the error control.  It is
!!       |              | defined only on a return from DLSOIBT with
!!       |              | ISTATE = 2.
!!
!!-----------------------------------------------------------------------
!!#### Part 2.  Other Routines Callable.
!!
!! The following are optional calls which the user may make to
!! gain additional capabilities in conjunction with DLSOIBT.
!! (The routines XSETUN and XSETF are designed to conform to the
!! SLATEC error handling package.)
!!
!!  Form of Call                |    Function
!!  --------------------------- | ------------------------------------------
!!   CALL XSETUN(LUN)           | Set the logical unit number, LUN, for
!!                              | output of messages from DLSOIBT, if
!!                              | the default is not desired.
!!                              | The default value of LUN is 6.
!!                              |
!!   CALL XSETF(MFLAG)          | Set a flag to control the printing of
!!                              | messages by DLSOIBT.
!!                              | MFLAG = 0 means do not print. (Danger:
!!                              | This risks losing valuable information.)
!!                              | MFLAG = 1 means print (the default).
!!                              |
!!                              | Either of the above calls may be made at
!!                              | any time and will take effect immediately.
!!                              |
!!   CALL DSRCOM(RSAV,ISAV,JOB) | saves and restores the contents of
!!                              | the internal Common blocks used by
!!                              | DLSOIBT (see Part 3 below).
!!                              | RSAV must be a real array of length 218
!!                              | or more, and ISAV must be an integer
!!                              | array of length 37 or more.
!!                              |      JOB=1 means save Common into RSAV/ISAV.
!!                              |      JOB=2 means restore Common from RSAV/ISAV.
!!                              | DSRCOM is useful if one is
!!                              | interrupting a run and restarting
!!                              | later, or alternating between two or
!!                              | more problems solved with DLSOIBT.
!!                              |
!!   CALL DINTDY(,,,,,)         | Provide derivatives of y, of various
!!        (see below)           | orders, at a specified point t, if
!!                              | desired.  It may be called only after
!!                              | a successful return from DLSOIBT.
!!
!! The detailed instructions for using DINTDY are as follows.
!! The form of the call is:
!!```fortran
!!   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!!```
!! The input parameters are:
!!
!! T
!!
!! : value of independent variable where answers are desired
!! (normally the same as the t last returned by DLSOIBT).
!! For valid results, T must lie between TCUR - HU and TCUR.
!! (See optional outputs for TCUR and HU.)
!!
!! K
!!
!! : integer order of the derivative desired.  K must satisfy
!! 0 .le. K .le. NQCUR, where NQCUR is the current order
!! (see optional outputs).  The capability corresponding
!! to K = 0, i.e. computing y(t), is already provided
!! by DLSOIBT directly.  Since NQCUR .ge. 1, the first
!! derivative dy/dt is always available with DINTDY.
!!
!! RWORK(21)
!!
!! : the base address of the history array YH.
!!
!! NYH
!!
!! : column length of YH, equal to the initial value of NEQ.
!!
!! The output parameters are:
!!
!! DKY
!!
!! : a real array of length NEQ containing the computed value
!! of the K-th derivative of y(t).
!! IFLAG
!!
!! : integer flag, returned as 0 if K and T were legal,
!! -1 if K was illegal, and -2 if T was illegal.
!! On an error return, a message is also written.
!!
!!-----------------------------------------------------------------------
!!#### Part 3.  Save and Restore Program State
!!
!! If the solution of a given problem by DLSOIBT is to be interrupted
!! and then later continued, such as when restarting an interrupted run
!! or alternating between two or more problems, the user should save,
!! following the return from the last DLSOIBT call prior to the
!! interruption, the contents of the call sequence variables and the
!! internal state variables, and later restore these values before the
!! next DLSOIBT call for that problem.  To save and restore the values
!! , use Subroutine DSRCOM (see Part 2 above).
!!
!!-----------------------------------------------------------------------
!!#### Part 4.  Optionally Replaceable Solver Routines.
!!
!! Below are descriptions of two routines in the DLSOIBT package which
!! relate to the measurement of errors.  Either routine can be
!! replaced by a user-supplied version, if desired.  However, since such
!! a replacement may have a major impact on performance, it should be
!! done only when absolutely necessary, and only with great caution.
!! (Note: The means by which the package version of a routine is
!! superseded by the user's version may be system-dependent.)
!!
!!#### (a) DEWSET()
!!
!! The following subroutine is called just before each internal
!! integration step, and sets the array of error weights, EWT, as
!! described under ITOL/RTOL/ATOL above:
!!
!!         SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
!!
!! where NEQ, ITOL, RTOL, and ATOL are as in the DLSOIBT call sequence,
!! YCUR contains the current dependent variable vector, and
!! EWT is the array of weights set by DEWSET.
!!
!! If the user supplies this subroutine, it must return in EWT(i)
!! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
!! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
!! routine (see below), and also used by DLSOIBT in the computation
!! of the optional output IMXER, the diagonal Jacobian approximation,
!! and the increments for difference quotient Jacobians.
!!
!! In the user-supplied version of DEWSET, it may be desirable to use
!! the current values of derivatives of y.  Derivatives up to order NQ
!! are available from the history array YH, described above under
!! optional outputs.  In DEWSET, YH is identical to the YCUR array,
!! extended to NQ + 1 columns with a column length of NYH and scale
!! factors of H\*\*j/factorial(j).  On the first call for the problem,
!! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
!! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
!! can be obtained by including in DEWSET the statements:
!!```fortran
!!      DOUBLE PRECISION RLS
!!      COMMON /DLS001/ RLS(218),ILS(37)
!!      NQ = ILS(33)
!!      NST = ILS(34)
!!      H = RLS(212)
!!```
!! Thus, for example, the current value of dy/dt can be obtained as
!! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
!! unnecessary when NST = 0).
!!
!!#### (b) DVNORM()
!!
!! The following is a real function routine which computes the weighted
!! root-mean-square norm of a vector v:
!!```text
!!     D = DVNORM (N, V, W)
!!```
!! where:
!!```text
!!   N = the length of the vector,
!!   V = real array of length N containing the vector,
!!   W = real array of length N containing weights,
!!   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
!!```
!! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
!! EWT is as set by Subroutine DEWSET.
!!
!! If the user supplies this function, it should return a non-negative
!! value of DVNORM suitable for use in the error control in DLSOIBT.
!! None of the arguments should be altered by DVNORM.
!! For example, a user-supplied DVNORM routine might:
!!  - substitute a max-norm of (V(i)\*W(i)) for the RMS-norm, or
!!  - ignore some components of V in the norm, with the effect of
!!    suppressing the error control on those components of y.
!!
!!-----------------------------------------------------------------------
!!### Pedigree:
!!
!! Reference:
!!     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!!     Solvers, in Scientific Computing,  R. S. Stepleman et al. (Eds.),
!!     North-Holland, Amsterdam, 1983, pp. 55-64.

!! This is based on the the 18 November 2003 version of ODEPACK
!!-----------------------------------------------------------------------
!!
!! Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
!!                Center for Applied Scientific Computing, L-561
!!                Lawrence Livermore National Laboratory
!!                Livermore, CA 94551
!! and
!!                Charles S. Kenney
!! formerly at:   Naval Weapons Center
!!                China Lake, CA 93555
!!
! -----------------------------------------------------------------------
! ### REVISION HISTORY  (YYYYMMDD)
!  19840625  DATE WRITTEN
!  19870330  Major update: corrected comments throughout;
!            removed TRET from Common; rewrote EWSET with 4 loops;
!            fixed t test in INTDY; added Cray directives in STODI;
!            in STODI, fixed DELP init. and logic around PJAC call;
!            combined routines to save/restore Common;
!            passed LEVEL = 0 in error message calls (except run abort).
!  20010425  Major update: convert source lines to upper case;
!            added *DECK lines; changed from 1 to * in dummy dimensions;
!            changed names R1MACH/D1MACH to RUMACH/DUMACH;
!            renamed routines for uniqueness across single/double prec.;
!            converted intrinsic names to generic form;
!            removed ILLIN and NTREP (data loaded) from Common;
!            removed all 'own' variables from Common;
!            changed error messages to quoted strings;
!            replaced XERRWV/XERRWD with 1993 revised version;
!            converted prologues, comments, error messages to mixed case;
!            converted arithmetic IF statements to logical IF statements;
!            numerous corrections to prologues and internal comments.
!  20010507  Converted single precision source to double precision.
!  20020502  Corrected declarations in descriptions of user routines.
!  20031105  Restored 'own' variables to Common block, to enable
!            interrupt/restart feature.
!  20031112  Added SAVE statements for data-loaded constants.
!  20031117  Changed internal names NRE, LSAVR to NFE, LSAVF resp.
!
! -----------------------------------------------------------------------
!  Other routines in the DLSOIBT package.
!
!  In addition to Subroutine DLSOIBT, the DLSOIBT package includes the
!  following subroutines and function routines:
!   DAIGBT   computes the initial value of the vector
!              dy/dt = A-inverse * g
!   DINTDY   computes an interpolated value of the y vector at t = TOUT.
!   DSTODI   is the core integrator, which does one step of the
!            integration and the associated error control.
!   DCFODE   sets all method coefficients and test constants.
!   DEWSET   sets the error weight vector EWT before each step.
!   DVNORM   computes the weighted RMS-norm of a vector.
!   DSRCOM   is a user-callable routine to save and restore
!            the contents of the internal Common blocks.
!   DPJIBT   computes and preprocesses the Jacobian matrix
!            and the Newton iteration matrix P.
!   DSLSBT   manages solution of linear system in chord iteration.
!   DDECBT and DSOLBT   are routines for solving block-tridiagonal
!            systems of linear algebraic equations.
!   DGEFA and DGESL   are routines from LINPACK for solving full
!            systems of linear algebraic equations.
!   DDOT     is one of the basic linear algebra modules (BLAS).
!   XERRWD, XSETUN, XSETF, and IXSAV handle the printing of all
!            error messages and warnings.  XERRWD is machine-dependent.
!  Note:  DVNORM, DDOT, and IXSAV are function routines.
!  All the others are subroutines.
!
!!-----------------------------------------------------------------------
subroutine dlsoibt(res,adda,jac,Neq,Y,Ydoti,T,Tout,Itol,Rtol,Atol,Itask,Istate,Iopt,Rwork,Lrw,Iwork,Liw,Mf)

external adda
external jac
external res
real(kind=dp), dimension(*) :: Atol, Rtol, Y, Ydoti
real(kind=dp) :: atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, size, sum, tcrit, tdist, tnext, tol, tolsf, tp, w0
integer :: i, i1, i2, ier, iflag, imxer, ires, kgo, leniw, lenrw, lenwm, lp, lyd0, mb, nb
logical :: ihit
integer :: Iopt, Istate, Itask, Itol, Liw, Lrw, Mf
integer, intent(inout), dimension(Liw) :: Iwork
integer, dimension(2), save :: mord
character(60) :: msg
integer, save :: mxhnl0, mxstp0
integer, dimension(*) :: Neq
real(kind=dp), intent(inout), dimension(Lrw) :: Rwork
real(kind=dp), intent(inout) :: T, Tout
!
! -----------------------------------------------------------------------
!  The following internal Common block contains
!  (a) variables which are local to any subroutine but whose values must
!      be preserved between calls to the routine ("own" variables), and
!  (b) variables which are communicated between subroutines.
!  The block DLS001 is declared in subroutines DLSOIBT, DINTDY, DSTODI,
!  DPJIBT, and DSLSBT.
!  Groups of variables are replaced by dummy arrays in the Common
!  declarations in routines where those variables are not used.
! -----------------------------------------------------------------------
!
data mord(1), mord(2)/12, 5/, mxstp0/500/, mxhnl0/10/
ihit=.false.
! -----------------------------------------------------------------------
!  Block A.
!  This code block is executed on every call.
!  It tests ISTATE and ITASK for legality and branches appropriately.
!  If ISTATE .gt. 1 but the flag INIT shows that initialization has
!  not yet been done, an error return occurs.
!  If ISTATE = 0 or 1 and TOUT = T, return immediately.
! -----------------------------------------------------------------------
if ( Istate<0 .or. Istate>3 ) then
! -----------------------------------------------------------------------
!  Block I.
!  The following block handles all error returns due to illegal input
!  (ISTATE = -3), as detected before calling the core integrator.
!  First the error message routine is called.  If the illegal input
!  is a negative ISTATE, the run is aborted (apparent infinite loop).
! -----------------------------------------------------------------------
   msg = 'DLSOIBT- ISTATE (=I1) illegal.'
   call xerrwd(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
   if ( Istate>=0 ) goto 1400
!
   msg = 'DLSOIBT- Run aborted.. apparent infinite loop.    '
   call xerrwd(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
   goto 99999
else
   if ( Itask<1 .or. Itask>5 ) then
      msg = 'DLSOIBT- ITASK (=I1) illegal. '
      call xerrwd(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
      goto 1400
   else
      if ( Istate<=1 ) then
         dls1%init = 0
         if ( Tout==T ) return
      elseif ( dls1%init==0 ) then
         msg = 'DLSOIBT- ISTATE.gt.1 but DLSOIBT not initialized. '
         call xerrwd(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
         goto 1400
      elseif ( Istate==2 ) then
         goto 50
      endif
! -----------------------------------------------------------------------
!  Block B.
!  The next code block is executed for the initial call (ISTATE = 0 or 1)
!  or for a continuation call with parameter changes (ISTATE = 3).
!  It contains checking of all inputs and various initializations.
!
!  First check legality of the non-optional inputs NEQ, ITOL, IOPT,
!  MF, MB, and NB.
! -----------------------------------------------------------------------
      if ( Neq(1)<=0 ) then
         msg = 'DLSOIBT- NEQ (=I1) .lt. 1     '
         call xerrwd(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
         goto 1400
      else
         if ( Istate>1 ) then
            if ( Neq(1)>dls1%n ) then
               msg = 'DLSOIBT- ISTATE = 3 and NEQ increased (I1 to I2). '
               call xerrwd(msg,50,5,0,2,dls1%n,Neq(1),0,0.0D0,0.0D0)
               goto 1400
            endif
         endif
         dls1%n = Neq(1)
         if ( Itol<1 .or. Itol>4 ) then
            msg = 'DLSOIBT- ITOL (=I1) illegal.  '
            call xerrwd(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
            goto 1400
         elseif ( Iopt<0 .or. Iopt>1 ) then
            msg = 'DLSOIBT- IOPT (=I1) illegal.  '
            call xerrwd(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
            goto 1400
         else
            dls1%meth = Mf/10
            dls1%miter = Mf - 10*dls1%meth
            if ( dls1%meth<1 .or. dls1%meth>2 ) goto 900
            if ( dls1%miter<1 .or. dls1%miter>2 ) goto 900
            mb = Iwork(1)
            nb = Iwork(2)
            if ( mb<1 .or. mb>dls1%n ) goto 1000
            if ( nb<4 ) then
               msg = 'DLSOIBT- NB (=I1) .lt. 4 illegal.       '
               call xerrwd(msg,40,10,0,1,nb,0,0,0.0D0,0.0D0)
               goto 1400
            else
               if ( mb*nb/=dls1%n ) goto 1000
!  Next process and check the optional inputs. --------------------------
               if ( Iopt==1 ) then
                  dls1%maxord = Iwork(5)
                  if ( dls1%maxord<0 ) then
                     msg = 'DLSOIBT- MAXORD (=I1) .lt. 0  '
                     call xerrwd(msg,30,11,0,1,dls1%maxord,0,0,0.0D0,0.0D0)
                     goto 1400
                  else
                     if ( dls1%maxord==0 ) dls1%maxord = 100
                     dls1%maxord = min(dls1%maxord,mord(dls1%meth))
                     dls1%mxstep = Iwork(6)
                     if ( dls1%mxstep<0 ) then
                        msg = 'DLSOIBT- MXSTEP (=I1) .lt. 0  '
                        call xerrwd(msg,30,12,0,1,dls1%mxstep,0,0,0.0D0,0.0D0)
                        goto 1400
                     else
                        if ( dls1%mxstep==0 ) dls1%mxstep = mxstp0
                        dls1%mxhnil = Iwork(7)
                        if ( dls1%mxhnil<0 ) then
                           msg = 'DLSOIBT- MXHNIL (=I1) .lt. 0  '
                           call xerrwd(msg,30,13,0,1,dls1%mxhnil,0,0,0.0D0,0.0D0)
                           goto 1400
                        else
                           if ( dls1%mxhnil==0 ) dls1%mxhnil = mxhnl0
                           if ( Istate<=1 ) then
                              h0 = Rwork(5)
                              if ( (Tout-T)*h0<0.0D0 ) then
                                 msg = 'DLSOIBT- TOUT (=R1) behind T (=R2)      '
                                 call xerrwd(msg,40,14,0,0,0,0,2,Tout,T)
                                 msg = '      Integration direction is given by H0 (=R1)  '
                                 call xerrwd(msg,50,14,0,0,0,0,1,h0,0.0D0)
                                 goto 1400
                              endif
                           endif
                           hmax = Rwork(6)
                           if ( hmax<0.0D0 ) then
                              msg = 'DLSOIBT- HMAX (=R1) .lt. 0.0  '
                              call xerrwd(msg,30,15,0,0,0,0,1,hmax,0.0D0)
                              goto 1400
                           else
                              dls1%hmxi = 0.0D0
                              if ( hmax>0.0D0 ) dls1%hmxi = 1.0D0/hmax
                              dls1%hmin = Rwork(7)
                              if ( dls1%hmin<0.0D0 ) then
                                 msg = 'DLSOIBT- HMIN (=R1) .lt. 0.0  '
                                 call xerrwd(msg,30,16,0,0,0,0,1,dls1%hmin,0.0D0)
                                 goto 1400
                              endif
                           endif
                        endif
                     endif
                  endif
               else
                  dls1%maxord = mord(dls1%meth)
                  dls1%mxstep = mxstp0
                  dls1%mxhnil = mxhnl0
                  if ( Istate<=1 ) h0 = 0.0D0
                  dls1%hmxi = 0.0D0
                  dls1%hmin = 0.0D0
               endif
! -----------------------------------------------------------------------
!  Set work array pointers and check lengths LRW and LIW.
!  Pointers to segments of RWORK and IWORK are named by prefixing L to
!  the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
!  Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
! -----------------------------------------------------------------------
               dls1%lyh = 21
               if ( Istate<=1 ) dls1%nyh = dls1%n
               dls1%lwm = dls1%lyh + (dls1%maxord+1)*dls1%nyh
               lenwm = 3*mb*mb*nb + 2
               dls1%lewt = dls1%lwm + lenwm
               dls1%lsavf = dls1%lewt + dls1%n
               dls1%lacor = dls1%lsavf + dls1%n
               lenrw = dls1%lacor + dls1%n - 1
               Iwork(17) = lenrw
               dls1%liwm = 1
               leniw = 20 + dls1%n
               Iwork(18) = leniw
               if ( lenrw>Lrw ) then
                  msg = 'DLSOIBT- RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
                  call xerrwd(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                  goto 1400
               elseif ( leniw>Liw ) then
                  msg = 'DLSOIBT- IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
                  call xerrwd(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                  goto 1400
               else
!  Check RTOL and ATOL for legality. ------------------------------------
                  rtoli = Rtol(1)
                  atoli = Atol(1)
                  do i = 1, dls1%n
                     if ( Itol>=3 ) rtoli = Rtol(i)
                     if ( Itol==2 .or. Itol==4 ) atoli = Atol(i)
                     if ( rtoli<0.0D0 ) then
                        msg = 'DLSOIBT- RTOL(=I1) is R1 .lt. 0.0       '
                        call xerrwd(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
                        goto 1400
                     elseif ( atoli<0.0D0 ) then
                        msg = 'DLSOIBT- ATOL(=I1) is R1 .lt. 0.0       '
                        call xerrwd(msg,40,20,0,1,i,0,1,atoli,0.0D0)
                        goto 1400
                     endif
                  enddo
                  if ( Istate<=1 ) then
! -----------------------------------------------------------------------
!  Block C.
!  The next block is for the initial call only (ISTATE = 0 or 1).
!  It contains all remaining initializations, the call to DAIGBT
!  (if ISTATE = 1), and the calculation of the initial step size.
!  The error weights in EWT are inverted after being loaded.
! -----------------------------------------------------------------------
                     dls1%uround = epsilon(0.0d0)
                     dls1%tn = T
                     if ( Itask==4 .or. Itask==5 ) then
                        tcrit = Rwork(1)
                        if ( (tcrit-Tout)*(Tout-T)<0.0D0 ) goto 1200
                        if ( h0/=0.0D0 .and. (T+h0-tcrit)*h0>0.0D0 ) h0 = tcrit - T
                     endif
                     dls1%jstart = 0
                     Rwork(dls1%lwm) = sqrt(dls1%uround)
                     dls1%nhnil = 0
                     dls1%nst = 0
                     dls1%nfe = 0
                     dls1%nje = 0
                     dls1%nslast = 0
                     dls1%hu = 0.0D0
                     dls1%nqu = 0
                     dls1%ccmax = 0.3D0
                     dls1%maxcor = 3
                     dls1%msbp = 20
                     dls1%mxncf = 10
!  Compute initial dy/dt, if necessary, and load it and initial Y into YH
                     lyd0 = dls1%lyh + dls1%nyh
                     lp = dls1%lwm + 1
                     if ( Istate==1 ) then
!  Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
                        do i = 1, dls1%n
                           Rwork(i+dls1%lyh-1) = Y(i)
                           Rwork(i+lyd0-1) = Ydoti(i)
                        enddo
                     else
!  DLSOIBT must compute initial dy/dt (LYD0 points to YH(*,2)). ---------
                        call daigbt(res,adda,Neq,T,Y,Rwork(lyd0),mb,nb,Rwork(lp),Iwork(21),ier)
                        dls1%nfe = dls1%nfe + 1
                        if ( ier<0 ) then
!  DAIGBT failed because a diagonal block of A matrix was singular. -----
                           ier = -ier
                           msg = 'DLSOIBT- Attempt to initialize dy/dt failed:  Matrix A has a'
                           call xerrwd(msg,60,207,0,0,0,0,0,0.0D0,0.0D0)
                           msg = '      singular diagonal block, block no. = (I1)   '
                           call xerrwd(msg,50,207,0,1,ier,0,0,0.0D0,0.0D0)
                           Istate = -8
                           return
                        elseif ( ier>0 ) then
!  DAIGBT failed because RES set IRES to 2 or 3. ------------------------
                           msg = 'DLSOIBT- Attempt to initialize dy/dt failed       '
                           call xerrwd(msg,50,208,0,0,0,0,0,0.0D0,0.0D0)
                           msg = '      because residual routine set its error flag '
                           call xerrwd(msg,50,208,0,0,0,0,0,0.0D0,0.0D0)
                           msg = '      to IRES = (I1)'
                           call xerrwd(msg,20,208,0,1,ier,0,0,0.0D0,0.0D0)
                           Istate = -8
                           return
                        else
                           do i = 1, dls1%n
                              Rwork(i+dls1%lyh-1) = Y(i)
                           enddo
                        endif
                     endif
!  Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                     dls1%nq = 1
                     dls1%h = 1.0D0
                     call dewset(dls1%n,Itol,Rtol,Atol,Rwork(dls1%lyh),Rwork(dls1%lewt))
                     do i = 1, dls1%n
                        if ( Rwork(i+dls1%lewt-1)<=0.0D0 ) then
                           ewti = Rwork(dls1%lewt+i-1)
                           msg = 'DLSOIBT- EWT(I1) is R1 .le. 0.0         '
                           call xerrwd(msg,40,21,0,1,i,0,1,ewti,0.0D0)
                           goto 1400
                        else
                           Rwork(i+dls1%lewt-1) = 1.0D0/Rwork(i+dls1%lewt-1)
                        endif
                     enddo
! -----------------------------------------------------------------------
!  The coding below computes the step size, H0, to be attempted on the
!  first step, unless the user has supplied a value for this.
!  First check that TOUT - T differs significantly from zero.
!  A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
!  if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
!  so as to be between 100*UROUND and 1.0E-3.
!  Then the computed value H0 is given by..
!                                       NEQ
!    H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
!                                        1
!  where   w0      = MAX ( ABS(T), ABS(TOUT) ),
!          YDOT(i) = i-th component of initial value of dy/dt,
!          ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
!  The sign of H0 is inferred from the initial values of TOUT and T.
! -----------------------------------------------------------------------
                     if ( h0==0.0D0 ) then
                        tdist = abs(Tout-T)
                        w0 = max(abs(T),abs(Tout))
                        if ( tdist<2.0D0*dls1%uround*w0 ) then
                           msg = 'DLSOIBT- TOUT(=R1) too close to T(=R2) to start integration.'
                           call xerrwd(msg,60,22,0,0,0,0,2,Tout,T)
                           goto 1400
                        else
                           tol = Rtol(1)
                           if ( Itol>2 ) then
                              do i = 1, dls1%n
                                 tol = max(tol,Rtol(i))
                              enddo
                           endif
                           if ( tol<=0.0D0 ) then
                              atoli = Atol(1)
                              do i = 1, dls1%n
                                 if ( Itol==2 .or. Itol==4 ) atoli = Atol(i)
                                 ayi = abs(Y(i))
                                 if ( ayi/=0.0D0 ) tol = max(tol,atoli/ayi)
                              enddo
                           endif
                           tol = max(tol,100.0D0*dls1%uround)
                           tol = min(tol,0.001D0)
                           sum = dvnorm(dls1%n,Rwork(lyd0),Rwork(dls1%lewt))
                           sum = 1.0D0/(tol*w0*w0) + tol*sum**2
                           h0 = 1.0D0/sqrt(sum)
                           h0 = min(h0,tdist)
                           h0 = sign(h0,Tout-T)
                        endif
                     endif
!  Adjust H0 if necessary to meet HMAX bound. ---------------------------
                     rh = abs(h0)*dls1%hmxi
                     if ( rh>1.0D0 ) h0 = h0/rh
!  Load H with H0 and scale YH(*,2) by H0. ------------------------------
                     dls1%h = h0
                     do i = 1, dls1%n
                        Rwork(i+lyd0-1) = h0*Rwork(i+lyd0-1)
                     enddo
                     goto 200
                  else
!  If ISTATE = 3, set flag to signal parameter changes to DSTODI. -------
                     dls1%jstart = -1
                     if ( dls1%nq>dls1%maxord ) then
!  MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
                        do i = 1, dls1%n
                           Ydoti(i) = Rwork(i+dls1%lwm-1)
                        enddo
                     endif
!  Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
                     Rwork(dls1%lwm) = sqrt(dls1%uround)
                     if ( dls1%n/=dls1%nyh ) then
!  NEQ was reduced.  Zero part of YH to avoid undefined references. -----
                        i1 = dls1%lyh + dls1%l*dls1%nyh
                        i2 = dls1%lyh + (dls1%maxord+1)*dls1%nyh - 1
                        if ( i1<=i2 ) then
                           do i = i1, i2
                              Rwork(i) = 0.0D0
                           enddo
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
   endif
! -----------------------------------------------------------------------
!  Block D.
!  The next code block is for continuation calls only (ISTATE = 2 or 3)
!  and is to check stop conditions before taking a step.
! -----------------------------------------------------------------------
 50   continue
   dls1%nslast = dls1%nst
   select case (Itask)
   case (2)
      goto 100
   case (3)
      tp = dls1%tn - dls1%hu*(1.0D0+100.0D0*dls1%uround)
      if ( (tp-Tout)*dls1%h>0.0D0 ) then
         msg = 'DLSOIBT- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
         call xerrwd(msg,60,23,0,1,Itask,0,2,Tout,tp)
         goto 1400
      else
         if ( (dls1%tn-Tout)*dls1%h>=0.0D0 ) goto 300
         goto 100
      endif
   case (4)
      tcrit = Rwork(1)
      if ( (dls1%tn-tcrit)*dls1%h>0.0D0 ) goto 1100
      if ( (tcrit-Tout)*dls1%h<0.0D0 ) goto 1200
      if ( (dls1%tn-Tout)*dls1%h>=0.0D0 ) then
         call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
         if ( iflag/=0 ) goto 1300
         T = Tout
         goto 400
      endif
   case (5)
      tcrit = Rwork(1)
      if ( (dls1%tn-tcrit)*dls1%h>0.0D0 ) goto 1100
   case default
      if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
      call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
      if ( iflag/=0 ) goto 1300
      T = Tout
      goto 400
   endselect
   hmx = abs(dls1%tn) + abs(dls1%h)
   ihit = abs(dls1%tn-tcrit)<=100.0D0*dls1%uround*hmx
   if ( ihit ) goto 300
   tnext = dls1%tn + dls1%h*(1.0D0+4.0D0*dls1%uround)
   if ( (tnext-tcrit)*dls1%h>0.0D0 ) then
      dls1%h = (tcrit-dls1%tn)*(1.0D0-4.0D0*dls1%uround)
      if ( Istate==2 ) dls1%jstart = -2
   endif
endif
! -----------------------------------------------------------------------
!  Block E.
!  The next block is normally executed for all calls and contains
!  the call to the one-step core integrator DSTODI.
!
!  This is a looping point for the integration steps.
!
!  First check for too many steps being taken, update EWT (if not at
!  start of problem), check for too much accuracy being requested, and
!  check for H below the roundoff level in T.
! -----------------------------------------------------------------------
 100  continue
if ( (dls1%nst-dls1%nslast)>=dls1%mxstep ) then
! -----------------------------------------------------------------------
!  Block H.
!  The following block handles all unsuccessful returns other than
!  those for illegal input.  First the error message routine is called.
!  If there was an error test or convergence test failure, IMXER is set.
!  Then Y is loaded from YH and T is set to TN.
!  The optional outputs are loaded into the work arrays before returning.
! -----------------------------------------------------------------------
!  The maximum number of steps was taken before reaching TOUT. ----------
   msg = 'DLSOIBT- At current T (=R1), MXSTEP (=I1) steps   '
   call xerrwd(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
   msg = '      taken on this call before reaching TOUT     '
   call xerrwd(msg,50,201,0,1,dls1%mxstep,0,1,dls1%tn,0.0D0)
   Istate = -1
   goto 600
else
   call dewset(dls1%n,Itol,Rtol,Atol,Rwork(dls1%lyh),Rwork(dls1%lewt))
   do i = 1, dls1%n
      if ( Rwork(i+dls1%lewt-1)<=0.0D0 ) then
!  EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
         ewti = Rwork(dls1%lewt+i-1)
         msg = 'DLSOIBT- At T (=R1), EWT(I1) has become R2 .le. 0.'
         call xerrwd(msg,50,202,0,1,i,0,2,dls1%tn,ewti)
         Istate = -6
         goto 700
      else
         Rwork(i+dls1%lewt-1) = 1.0D0/Rwork(i+dls1%lewt-1)
      endif
   enddo
endif
 200  continue
tolsf = dls1%uround*dvnorm(dls1%n,Rwork(dls1%lyh),Rwork(dls1%lewt))
if ( tolsf<=1.0D0 ) then
   if ( (dls1%tn+dls1%h)==dls1%tn ) then
      dls1%nhnil = dls1%nhnil + 1
      if ( dls1%nhnil<=dls1%mxhnil ) then
         msg = 'DLSOIBT- Warning..Internal T (=R1) and H (=R2) are'
         call xerrwd(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
         msg = '      such that in the machine, T + H = T on the next step  '
         call xerrwd(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
         msg = '     (H = step size). Solver will continue anyway.'
         call xerrwd(msg,50,101,0,0,0,0,2,dls1%tn,dls1%h)
         if ( dls1%nhnil>=dls1%mxhnil ) then
            msg = 'DLSOIBT- Above warning has been issued I1 times.  '
            call xerrwd(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
            msg = '     It will not be issued again for this problem.'
            call xerrwd(msg,50,102,0,1,dls1%mxhnil,0,0,0.0D0,0.0D0)
         endif
      endif
   endif
! -----------------------------------------------------------------------
!      CALL DSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
!                  ADDA,JAC,DPJIBT,DSLSBT)
!  Note: SAVF in DSTODI occupies the same space as YDOTI in DLSOIBT.
! -----------------------------------------------------------------------
   call dstodi(Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dls1%lyh),Rwork(dls1%lewt), &
   & Ydoti,Rwork(dls1%lsavf),Rwork(dls1%lacor),           &
             & Rwork(dls1%lwm),Iwork(dls1%liwm),res,adda,jac,dpjibt,dslsbt)
   kgo = 1 - dls1%kflag
   select case (kgo)
   case (2)
!  KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
      msg = 'DLSOIBT- At T (=R1) and step size H (=R2), the    '
      call xerrwd(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
      msg = 'error test failed repeatedly or with ABS(H) = HMIN'
      call xerrwd(msg,50,204,0,0,0,0,2,dls1%tn,dls1%h)
      Istate = -4
      goto 500
   case (3)
!  KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
      msg = 'DLSOIBT- At T (=R1) and step size H (=R2), the    '
      call xerrwd(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      corrector convergence failed repeatedly     '
      call xerrwd(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      or with ABS(H) = HMIN   '
      call xerrwd(msg,30,205,0,0,0,0,2,dls1%tn,dls1%h)
      Istate = -5
      goto 500
   case (4)
   case (5)
!  IRES = 3 returned by RES, despite retries by DSTODI.------------------
      msg = 'DLSOIBT- At T (=R1) residual routine returned     '
      call xerrwd(msg,50,206,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      error IRES = 3 repeatedly.        '
      call xerrwd(msg,40,206,0,0,0,0,1,dls1%tn,0.0D0)
      Istate = -7
      goto 700
   case default
!
!  KGO = 1:success; 2:error test failure; 3:convergence failure;
!        4:RES ordered return; 5:RES returned error.
! -----------------------------------------------------------------------
!  Block F.
!  The following block handles the case of a successful return from the
!  core integrator (KFLAG = 0).  Test for stop conditions.
! -----------------------------------------------------------------------
      dls1%init = 1
      select case (Itask)
      case (2)
      case (3)
!  ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
         if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
      case (4)
!  ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
         if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) then
            hmx = abs(dls1%tn) + abs(dls1%h)
            ihit = abs(dls1%tn-tcrit)<=100.0D0*dls1%uround*hmx
            if ( .not.(ihit) ) then
               tnext = dls1%tn + dls1%h*(1.0D0+4.0D0*dls1%uround)
               if ( (tnext-tcrit)*dls1%h>0.0D0 ) then
                  dls1%h = (tcrit-dls1%tn)*(1.0D0-4.0D0*dls1%uround)
                  dls1%jstart = -2
               endif
               goto 100
            endif
         else
            call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
            T = Tout
            goto 400
         endif
      case (5)
!  ITASK = 5.  see if TCRIT was reached and jump to exit. ---------------
         hmx = abs(dls1%tn) + abs(dls1%h)
         ihit = abs(dls1%tn-tcrit)<=100.0D0*dls1%uround*hmx
      case default
!  ITASK = 1.  If TOUT has been reached, interpolate. -------------------
         if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
         call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
         T = Tout
         goto 400
      endselect
   endselect
else
   tolsf = tolsf*2.0D0
   if ( dls1%nst==0 ) then
      msg = 'DLSOIBT- At start of problem, too much accuracy   '
      call xerrwd(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      requested for precision of machine..  See TOLSF (=R1) '
      call xerrwd(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
      Rwork(14) = tolsf
      goto 1400
   else
!  Too much accuracy requested for machine precision. -------------------
      msg = 'DLSOIBT- At T (=R1), too much accuracy requested  '
      call xerrwd(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      for precision of machine..  See TOLSF (=R2) '
      call xerrwd(msg,50,203,0,0,0,0,2,dls1%tn,tolsf)
      Rwork(14) = tolsf
      Istate = -2
      goto 700
   endif
endif
! -----------------------------------------------------------------------
!  Block G.
!  The following block handles all successful returns from DLSOIBT.
!  If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
!  ISTATE is set to 2, and the optional outputs are loaded into the
!  work arrays before returning.
! -----------------------------------------------------------------------
 300  continue
do i = 1, dls1%n
   Y(i) = Rwork(i+dls1%lyh-1)
enddo
T = dls1%tn
if ( Itask==4 .or. Itask==5 ) then
   if ( ihit ) T = tcrit
endif
 400  continue
Istate = 2
if ( dls1%kflag==-3 ) Istate = 3
Rwork(11) = dls1%hu
Rwork(12) = dls1%h
Rwork(13) = dls1%tn
Iwork(11) = dls1%nst
Iwork(12) = dls1%nfe
Iwork(13) = dls1%nje
Iwork(14) = dls1%nqu
Iwork(15) = dls1%nq
return
!  Compute IMXER if relevant. -------------------------------------------
 500  continue
big = 0.0D0
imxer = 1
do i = 1, dls1%n
   size = abs(Rwork(i+dls1%lacor-1)*Rwork(i+dls1%lewt-1))
   if ( big<size ) then
      big = size
      imxer = i
   endif
enddo
Iwork(16) = imxer
!  Compute residual if relevant. ----------------------------------------
 600  continue
lyd0 = dls1%lyh + dls1%nyh
do i = 1, dls1%n
   Rwork(i+dls1%lsavf-1) = Rwork(i+lyd0-1)/dls1%h
   Y(i) = Rwork(i+dls1%lyh-1)
enddo
ires = 1
call res(Neq,dls1%tn,Y,Rwork(dls1%lsavf),Ydoti,ires)
dls1%nfe = dls1%nfe + 1
if ( ires>1 ) then
   msg = 'DLSOIBT- Residual routine set its flag IRES       '
   call xerrwd(msg,50,210,0,0,0,0,0,0.0D0,0.0D0)
   msg = '      to (I1) when called for final output.       '
   call xerrwd(msg,50,210,0,1,ires,0,0,0.0D0,0.0D0)
endif
goto 800
!  Set Y vector, T, and optional outputs. -------------------------------
 700  continue
do i = 1, dls1%n
   Y(i) = Rwork(i+dls1%lyh-1)
enddo
 800  continue
T = dls1%tn
Rwork(11) = dls1%hu
Rwork(12) = dls1%h
Rwork(13) = dls1%tn
Iwork(11) = dls1%nst
Iwork(12) = dls1%nfe
Iwork(13) = dls1%nje
Iwork(14) = dls1%nqu
Iwork(15) = dls1%nq
return
 900  continue
msg = 'DLSOIBT- MF (=I1) illegal.    '
call xerrwd(msg,30,8,0,1,Mf,0,0,0.0D0,0.0D0)
goto 1400
 1000 continue
msg = 'DLSOIBT- MB (=I1) or NB (=I2) illegal.  '
call xerrwd(msg,40,9,0,2,mb,nb,0,0.0D0,0.0D0)
goto 1400
 1100 continue
msg = 'DLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
call xerrwd(msg,60,24,0,0,0,0,2,tcrit,dls1%tn)
goto 1400
 1200 continue
msg = 'DLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
call xerrwd(msg,60,25,0,0,0,0,2,tcrit,Tout)
goto 1400
 1300 continue
msg = 'DLSOIBT- Trouble in DINTDY.  ITASK = I1, TOUT = R1'
call xerrwd(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
!
 1400 continue
Istate = -3
return
99999 continue
end subroutine dlsoibt
