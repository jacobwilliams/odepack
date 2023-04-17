!==================================================================================================================================!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!==================================================================================================================================!
!>
!!### Synopsis
!!
!! _DLSODAR_ solves the initial value problem for stiff or nonstiff
!! systems of first order ODEs of the form
!!```text
!!     dy/dt = f(t,y),  or, in component form,
!!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
!!```
!! with Automatic method switching for stiff and nonstiff problems,
!! and with Root-finding.
!!
!! At the same time, it locates the roots of any of a set of functions
!!```text
!!     g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).
!!```
!! This a variant version of the _DLSODE_ package.  It differs from it
!! in two ways:
!!
!! (a) It switches automatically between stiff and nonstiff methods.
!! This means that the user does not have to determine whether the
!! problem is stiff or not, and the solver will automatically choose the
!! appropriate method.  It always starts with the nonstiff method.
!!
!! (b) It finds the root of at least one of a set of constraint
!! functions g(i) of the independent and dependent variables.
!! It finds only those roots for which some g(i), as a function
!! of t, changes sign in the interval of integration.
!! It then returns the solution at the root, if that occurs
!! sooner than the specified stop condition, and otherwise returns
!! the solution according the specified stop condition.
!!
!!-----------------------------------------------------------------------
!!#### Summary of Usage.
!!
!! Communication between the user and the DLSODAR package, for normal
!! situations, is summarized here.  This summary describes only a subset
!! of the full set of options available.  See the full description for
!! details, including alternative treatment of the Jacobian matrix,
!! optional inputs and outputs, nonstandard options, and
!! instructions for special situations.  See also the example
!! problem (with program and output) following this summary.
!!
!!#### Step A.
!!
!! First provide a subroutine of the form:
!!```fortran
!!      SUBROUTINE F (NEQ, T, Y, YDOT)
!!      INTEGER NEQ
!!      DOUBLE PRECISION T, Y(*), YDOT(*)
!!```
!! which supplies the vector function f by loading YDOT(i) with f(i).
!!
!!#### Step B.
!!
!! Provide a subroutine of the form:
!!```fortran
!!      SUBROUTINE G (NEQ, T, Y, NG, GOUT)
!!      INTEGER NEQ
!!      DOUBLE PRECISION T, Y(*), GOUT(NG)
!!```
!! which supplies the vector function g by loading GOUT(i) with
!! g(i), the i-th constraint function whose root is sought.
!!
!!#### Step C.
!!
!! Write a main program which calls Subroutine DLSODAR once for
!! each point at which answers are desired.  This should also provide
!! for possible use of logical unit 6 for output of error messages by
!! DLSODAR.
!!
!!#### On the first call to DLSODAR, supply arguments as follows:
!!
!! F
!!
!! : name of subroutine for right-hand side vector f.
!! This name must be declared External in calling program.
!! NEQ
!!
!! : number of first order ODEs.
!! Y
!!
!! : array of initial values, of length NEQ.
!! T
!!
!! : the initial value of the independent variable.
!! TOUT
!!
!! : first point where output is desired (.ne. T).
!! ITOL
!!
!! : 1 or 2 according as ATOL (below) is a scalar or array.
!! RTOL
!!
!! : relative tolerance parameter (scalar).
!! ATOL
!!
!! : absolute tolerance parameter (scalar or array).
!! the estimated local error in y(i) will be controlled so as
!! to be less than
!!```text
!!       EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!!       EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!!```
!! Thus the local error test passes if, in each component,
!! either the absolute error is less than ATOL (or ATOL(i)),
!! or the relative error is less than RTOL.
!! Use RTOL = 0.0 for pure absolute error control, and
!! use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!! control.  _Caution_: actual (global) errors may exceed these
!! local tolerances, so choose them conservatively.
!!
!! ITASK
!!
!! : 1 for normal computation of output values of y at t = TOUT.
!!
!! ISTATE
!!
!! : integer flag (input and output).  Set ISTATE = 1.
!!
!! IOPT
!!
!! : 0 to indicate no optional inputs used.
!!
!! RWORK
!!
!! : real work array of length at least:
!!```text
!!      22 + NEQ * MAX(16, NEQ + 9) + 3*NG.
!!```
!! See also Paragraph F below.
!!
!! LRW
!!
!! : declared length of RWORK (in user's dimension).
!!
!! IWORK
!!
!! : integer work array of length at least  20 + NEQ.
!!
!! LIW
!!
!! : declared length of IWORK (in user's dimension).
!!
!! JAC
!!
!! : name of subroutine for Jacobian matrix.
!! Use a dummy name to a noop function.  See also Paragraph F below.
!!
!! JT
!!
!! : Jacobian type indicator.  Set JT = 2.
!! See also Paragraph F below.
!!
!! G
!!
!! : name of subroutine for constraint functions, whose
!! roots are desired during the integration.
!! This name must be declared External in calling program.
!!
!! NG
!!
!! : number of constraint functions g(i).  If there are none,
!! set NG = 0, and pass a dummy name for G.
!!
!! JROOT
!!
!! : integer array of length NG for output of root information.
!! See next paragraph.
!! Note that the main program must declare arrays Y, RWORK, IWORK,
!! JROOT, and possibly ATOL.
!!
!!#### Step D.
!!
!! The output from the first call (or any call) is:
!!
!! Y
!!
!! : array of computed values of y(t) vector.
!!
!! T
!!
!! : corresponding value of independent variable.  This is
!! TOUT if ISTATE = 2, or the root location if ISTATE = 3,
!! or the farthest point reached if DLSODAR was unsuccessful.
!! ISTATE = 2 or 3  if DLSODAR was successful, negative otherwise.
!!
!! ISTATE | Description
!! ------ | --------------------------------------------------------------
!!      2 | means no root was found, and TOUT was reached as desired.
!!      3 | means a root was found prior to reaching TOUT.
!!     -1 | means excess work done on this call (perhaps wrong JT).
!!     -2 | means excess accuracy requested (tolerances too small).
!!     -3 | means illegal input detected (see printed message).
!!     -4 | means repeated error test failures (check all inputs).
!!     -5 | means repeated convergence failures (perhaps bad Jacobian
!!        | supplied or wrong choice of JT or tolerances).
!!     -6 | means error weight became zero during problem. (Solution
!!        | component i vanished, and ATOL or ATOL(i) = 0.)
!!     -7 | means work space insufficient to finish (see messages).
!!        |
!!
!! JROOT
!!
!! : array showing roots found if ISTATE = 3 on return.
!! JROOT(i) = 1 if g(i) has a root at t, or 0 otherwise.
!!
!!#### Step E.
!!
!! To continue the integration after a successful return, proceed as follows:
!!
!! (a)  If ISTATE = 2 on return, reset TOUT and call DLSODAR again.

!! (b)  If ISTATE = 3 on return, reset ISTATE to 2, call DLSODAR again.
!!
!! In either case, no other parameters need be reset.
!!
!!#### Step F.
!!
!! **Notes:**
!!
!! If and when DLSODAR regards the problem as stiff, and
!! switches methods accordingly, it must make use of the NEQ by NEQ
!! Jacobian matrix, J = df/dy.  For the sake of simplicity, the
!! inputs to DLSODAR recommended in Paragraph C above cause DLSODAR to
!! treat J as a full matrix, and to approximate it internally by
!! difference quotients.
!!
!! Alternatively, J can be treated as a band
!! matrix (with great potential reduction in the size of the RWORK
!! array).
!!
!! Also, in either the full or banded case, the user can supply
!! J in closed form, with a routine whose name is passed as the JAC
!! argument.  These alternatives are described in the paragraphs on
!! RWORK, JAC, and JT in the full description of the call sequence below.
!!
!!-----------------------------------------------------------------------
!!### Example Problem.
!!
!! The following is a simple example problem, with the coding
!! needed for its solution by DLSODAR.  The problem is from chemical
!! kinetics, and consists of the following three rate equations:
!!```text
!!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!!     dy3/dt = 3.e7*y2**2
!!```
!! on the interval from t = 0.0 to t = 4.e10, with initial conditions
!! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
!!
!! In addition, we want to find the values of t, y1, y2, and y3 at which
!!
!!   1.  y1 reaches the value 1.e-4, and
!!   2.  y3 reaches the value 1.e-2.
!!
!! The following coding solves this problem with DLSODAR,
!! printing results at t = .4, 4., ..., 4.e10, and at the computed
!! roots.  It uses ITOL = 2 and ATOL much smaller for y2 than y1 or y3
!! because y2 has much smaller values.
!!
!! At the end of the run, statistical quantities of interest are
!! printed (see optional outputs in the full description below).
!!
!!```fortran
!!program dlsodar_ex
!!use m_odepack
!!implicit none
!!external fex
!!external gex
!!external jdum
!!integer,parameter :: dp=kind(0.0d0)
!!
!!real(kind=dp),dimension(3)   ::  atol,y
!!integer                      ::  iopt,iout,istate,itask,itol,jt,liw,lrw,neq,ng
!!integer,dimension(23)        ::  iwork
!!integer,dimension(2)         ::  jroot
!!real(kind=dp)                ::  rtol,t,tout
!!real(kind=dp),dimension(76)  ::  rwork
!!
!!   neq = 3
!!   y(1) = 1.
!!   y(2) = 0.
!!   y(3) = 0.
!!   t = 0.
!!   tout = .4
!!   itol = 2
!!   rtol = 1.D-4
!!   atol(1) = 1.D-6
!!   atol(2) = 1.D-10
!!   atol(3) = 1.D-6
!!   itask = 1
!!   istate = 1
!!   iopt = 0
!!   lrw = 76
!!   liw = 23
!!   jt = 2
!!   ng = 2
!!   do iout = 1,12
!!      do
!!         call dlsodar(fex,[neq],y,t,tout,itol,[rtol],atol,itask,istate,    &
!!                    & iopt,rwork,lrw,iwork,liw,jdum,jt,gex,ng,jroot)
!!         write (6,99010) t,y(1),y(2),y(3)
!!   99010 format (' At t =',d12.4,'   Y =',3D14.6)
!!         if ( istate<0 ) then
!!            write (6,99020) istate
!!   99020    format (///' Error halt.. ISTATE =',i3)
!!            stop 1
!!         elseif ( istate==2 ) then
!!            tout = tout*10.
!!            exit
!!         else
!!            write (6,99030) jroot(1),jroot(2)
!!   99030    format (5x,' The above line is a root,  JROOT =',2I5)
!!            istate = 2
!!         endif
!!      enddo
!!   enddo
!!   write (6,99040) iwork(11),iwork(12),iwork(13),iwork(10),        &
!!                 & iwork(19),rwork(15)
!!   99040 format (/' No. steps =',i4,'  No. f-s =',i4,'  No. J-s =',i4,     &
!!                &'  No. g-s =',i4/' Method last used =',i2,                &
!!                &'   Last switch was at t =',d12.4)
!!
!!end program dlsodar_ex
!!
!!subroutine jdum()
!!implicit none
!!end subroutine jdum
!!
!!subroutine fex(Neq,T,Y,Ydot)
!!implicit none
!!integer,parameter                         ::  dp=kind(0.0d0)
!!integer                                   ::  Neq
!!real(kind=dp)                             ::  T
!!real(kind=dp),intent(in),dimension(3)     ::  Y
!!real(kind=dp),intent(inout),dimension(3)  ::  Ydot
!!
!!   Ydot(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
!!   Ydot(3) = 3.D7*Y(2)*Y(2)
!!   Ydot(2) = -Ydot(1) - Ydot(3)
!!end subroutine fex
!!
!!subroutine gex(Neq,T,Y,Ng,Gout)
!!implicit none
!!
!!integer                                 ::  Neq
!!integer,parameter                       ::  dp=kind(0.0d0)
!!real(kind=dp)                           ::  T
!!real(kind=dp),intent(in),dimension(3)   ::  Y
!!integer                                 ::  Ng
!!real(kind=dp),intent(out),dimension(2)  ::  Gout
!!   Gout(1) = Y(1) - 1.D-4
!!   Gout(2) = Y(3) - 1.D-2
!!end subroutine gex
!!```
!!
!! The output of this program (on a CDC-7600 in single precision)
!! is as follows:
!!```text
!!   At t =  2.6400e-01   y =  9.899653e-01  3.470563e-05  1.000000e-02
!!        The above line is a root,  JROOT =    0    1
!!   At t =  4.0000e-01   Y =  9.851712e-01  3.386380e-05  1.479493e-02
!!   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02
!!   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01
!!   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01
!!   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01
!!   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01
!!   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01
!!   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01
!!   At t =  2.0745e+07   Y =  1.000000e-04  4.000395e-10  9.999000e-01
!!        The above line is a root,  JROOT =    1    0
!!   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01
!!   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01
!!   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01
!!   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00
!!
!!   No. steps = 361  No. f-s = 693  No. J-s =  64  No. g-s = 390
!!   Method last used = 2   Last switch was at t =  6.0092e-03
!!```
!!-----------------------------------------------------------------------
!!### Full Description of User Interface to DLSODAR.
!!
!! The user interface to DLSODAR consists of the following parts.
!!
!! 1.   The call sequence to Subroutine DLSODAR, which is a driver
!!      routine for the solver.  This includes descriptions of both
!!      the call sequence arguments and of user-supplied routines.
!!      Following these descriptions is a description of
!!      optional inputs available through the call sequence, and then
!!      a description of optional outputs (in the work arrays).
!!
!! 2.   Descriptions of other routines in the DLSODAR package that may be
!!      (optionally) called by the user.  These provide the ability to
!!      alter error message handling, save and restore the internal
!!      Common, and obtain specified derivatives of the solution y(t).
!!
!! 3.   Descriptions of Common blocks to be declared in overlay
!!      or similar environments, or to be saved when doing an interrupt
!!      of the problem and continued solution later.
!!
!! 4.   Description of a subroutine in the DLSODAR package,
!!      which the user may replace with his/her own version, if desired.
!!      this relates to the measurement of errors.
!!
!!-----------------------------------------------------------------------
!!### Part 1.  Call Sequence.
!!
!! The call sequence parameters used for input only are
!! F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC,
!! JT, G, and NG.
!!
!! Used only for output is JROOT,
!!
!! Those used for both input and output are Y, T, ISTATE.
!!
!! The work arrays RWORK and IWORK are also used for conditional and
!! optional inputs and optional outputs.  (The term output here refers
!! to the return from Subroutine DLSODAR to the user's calling program.)
!!
!! The legality of input parameters will be thoroughly checked on the
!! initial call for the problem, but not checked thereafter unless a
!! change in input parameters is flagged by ISTATE = 3 on input.
!!
!!## Arguments
!! The descriptions of the call arguments are as follows.
!!
!! F
!!
!! : the name of the user-supplied subroutine defining the
!! ODE system.  The system must be put in the first-order
!! form dy/dt = f(t,y), where f is a vector-valued function
!! of the scalar t and the vector y.  Subroutine F is to
!! compute the function f.  It is to have the form
!!```fortran
!!       SUBROUTINE F (NEQ, T, Y, YDOT)
!!       DOUBLE PRECISION T, Y(*), YDOT(*)
!!```
!! where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!! is output.  Y and YDOT are arrays of length NEQ.
!! Subroutine F should not alter Y(1),...,Y(NEQ).
!! F must be declared External in the calling program.
!!
!! Subroutine F may access user-defined quantities in
!! NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!! (dimensioned in F) and/or Y has length exceeding NEQ(1).
!! See the descriptions of NEQ and Y below.
!!
!! If quantities computed in the F routine are needed
!! externally to DLSODAR, an extra call to F should be made
!! for this purpose, for consistent and accurate results.
!! If only the derivative dy/dt is needed, use DINTDY instead.
!!
!! NEQ
!!
!! : the size of the ODE system (number of first order
!! ordinary differential equations).  Used only for input.
!! NEQ may be decreased, but not increased, during the problem.
!! If NEQ is decreased (with ISTATE = 3 on input), the
!! remaining components of Y should be left undisturbed, if
!! these are to be accessed in F and/or JAC.
!!
!! Normally, NEQ is a scalar, and it is generally referred to
!! as a scalar in this user interface description.  However,
!! NEQ may be an array, with NEQ(1) set to the system size.
!! (The DLSODAR package accesses only NEQ(1).)  In either case,
!! this parameter is passed as the NEQ argument in all calls
!! to F, JAC, and G.  Hence, if it is an array, locations
!! NEQ(2),... may be used to store other integer data and pass
!! it to F, JAC, and G.  Each such subroutine must include
!! NEQ in a Dimension statement in that case.
!!
!! Y
!!
!! : = a real array for the vector of dependent variables, of
!! NEQ or more.  Used for both input and output on the
!! first call (ISTATE = 1), and only for output on other calls.
!! On the first call, Y must contain the vector of initial
!! values.  On output, Y contains the computed solution vector,
!! evaluated at T.  If desired, the Y array may be used
!! for other purposes between calls to the solver.
!!
!! This array is passed as the Y argument in all calls to F,
!! JAC, and G.  Hence its length may exceed NEQ, and locations
!! Y(NEQ+1),... may be used to store other real data and
!! pass it to F, JAC, and G.  (The DLSODAR package accesses only
!! Y(1),...,Y(NEQ).)
!!
!! T
!!
!! : = the independent variable.  On input, T is used only on the
!! first call, as the initial point of the integration.
!! On output, after each call, T is the value at which a
!! computed solution y is evaluated (usually the same as TOUT).
!! If a root was found, T is the computed location of the
!! root reached first, on output.
!!
!! On an error return, T is the farthest point reached.
!!
!! TOUT
!!
!! : = the next value of t at which a computed solution is desired.
!! Used only for input.
!!
!! When starting the problem (ISTATE = 1), TOUT may be equal
!! to T for one call, then should .ne. T for the next call.
!! For the initial T, an input value of TOUT .ne. T is used
!! in order to determine the direction of the integration
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
!!      max-norm of ( E(i)/EWT(i) )   .le.   1,
!! where EWT = (EWT(i)) is a vector of positive error weights.
!! The values of RTOL and ATOL should all be non-negative.
!! The following table gives the types (scalar/array) of
!! RTOL and ATOL, and the corresponding form of EWT(i).
!!
!!  ITOL |  RTOL   |   ATOL   |      EWT(i)
!!  ---- |-------- | -------- | -----------------------------
!!   1   | scalar  |  scalar  |  RTOL\*ABS(Y(i)) + ATOL
!!   2   | scalar  |  array   |  RTOL\*ABS(Y(i)) + ATOL(i)
!!   3   | array   |  scalar  |  RTOL(i)\*ABS(Y(i)) + ATOL
!!   4   | array   |  array   |  RTOL(i)\*ABS(Y(i)) + ATOL(i)
!!
!! When either of these parameters is a scalar, it need not
!! be dimensioned in the user's calling program.
!!
!! If none of the above choices (with ITOL, RTOL, and ATOL
!! fixed throughout the problem) is suitable, more general
!! error controls can be obtained by substituting a
!! user-supplied routine for the setting of EWT.
!! See Part 4 below.
!!
!! If global errors are to be estimated by making a repeated
!! run on the same problem with smaller tolerances, then all
!! components of RTOL and ATOL (i.e. of EWT) should be scaled
!! down uniformly.
!!
!! ITASK
!!
!! : an index specifying the task to be performed.
!! input only.  ITASK has the following values and meanings.
!!
!! value |description
!! ----- |---------------------------------------------------------------
!!     1 |means normal computation of output values of y(t) at
!!       |t = TOUT (by overshooting and interpolating).
!!       |
!!     2 |means take one step only and return.
!!       |
!!     3 |means stop at the first internal mesh point at or
!!       |beyond t = TOUT and return.
!!       |
!!     4 |means normal computation of output values of y(t) at
!!       |t = TOUT but without overshooting t = TCRIT.
!!       |TCRIT must be input as RWORK(1).  TCRIT may be equal to
!!       |or beyond TOUT, but not behind it in the direction of
!!       |integration.  This option is useful if the problem
!!       |has a singularity at or beyond t = TCRIT.
!!       |
!!     5 |means take one step, without passing TCRIT, and return.
!!       |TCRIT must be input as RWORK(1).
!!       |
!!       |  Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!!       |  (within roundoff), it will return T = TCRIT (exactly) to
!!       |  indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!!       |  in which case answers at t = TOUT are returned first).
!!
!! ISTATE
!!
!! : an index used for input and output to specify the
!! the state of the calculation.
!!
!!#### On input, the values of ISTATE are as follows.
!!
!! value | description
!! ----- |---------------------------------------------------------------
!!     1 | means this is the first call for the problem
!!       | (initializations will be done).  See note below.
!!       |
!!     2 | means this is not the first call, and the calculation
!!       | is to continue normally, with no change in any input
!!       | parameters except possibly TOUT and ITASK.
!!       | (If ITOL, RTOL, and/or ATOL are changed between calls
!!       | with ISTATE = 2, the new values will be used but not
!!       | tested for legality.)
!!       |
!!     3 | means this is not the first call, and the
!!       | calculation is to continue normally, but with
!!       | a change in input parameters other than
!!       | TOUT and ITASK.  Changes are allowed in
!!       | NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
!!       | and any optional inputs except H0, MXORDN, and MXORDS.
!!       | (See IWORK description for ML and MU.)
!!       | In addition, immediately following a return with
!!       | ISTATE = 3 (root found), NG and G may be changed.
!!       | (But changing NG from 0 to .gt. 0 is not allowed.)
!!
!! Note:  A preliminary call with TOUT = T is not counted
!! as a first call here, as no initialization or checking of
!! input is done.  (Such a call is sometimes useful for the
!! purpose of outputting the initial conditions.)
!! Thus the first call for which TOUT .ne. T requires
!! ISTATE = 1 on input.
!!
!! On output, ISTATE has the following values and meanings.
!!
!!  value    | description
!!  -------- | --------------------------------------------------------
!!         1 | means nothing was done; TOUT = t and ISTATE = 1 on input.
!!           |
!!         2 | means the integration was performed successfully, and
!!           | no roots were found.
!!           |
!!         3 | means the integration was successful, and one or more
!!           | roots were found before satisfying the stop condition
!!           | specified by ITASK.  See JROOT.
!!           |
!!        -1 | means an excessive amount of work (more than MXSTEP
!!           | steps) was done on this call, before completing the
!!           | requested task, but the integration was otherwise
!!           | successful as far as T.  (MXSTEP is an optional input
!!           | and is normally 500.)  To continue, the user may
!!           | simply reset ISTATE to a value .gt. 1 and call again
!!           | (the excess work step counter will be reset to 0).
!!           | In addition, the user may increase MXSTEP to avoid
!!           | this error return (see below on optional inputs).
!!           |
!!        -2 | means too much accuracy was requested for the precision
!!           | of the machine being used.  This was detected before
!!           | completing the requested task, but the integration
!!           | was successful as far as T.  To continue, the tolerance
!!           | parameters must be reset, and ISTATE must be set
!!           | to 3.  The optional output TOLSF may be used for this
!!           | purpose.  (Note: If this condition is detected before
!!           | taking any steps, then an illegal input return
!!           | (ISTATE = -3) occurs instead.)
!!           |
!!        -3 | means illegal input was detected, before taking any
!!           | integration steps.  See written message for details.
!!           | Note:  If the solver detects an infinite loop of calls
!!           | to the solver with illegal input, it will cause
!!           | the run to stop.
!!           |
!!        -4 | means there were repeated error test failures on
!!           | one attempted step, before completing the requested
!!           | task, but the integration was successful as far as T.
!!           | The problem may have a singularity, or the input
!!           | may be inappropriate.
!!           |
!!        -5 | means there were repeated convergence test failures on
!!           | one attempted step, before completing the requested
!!           | task, but the integration was successful as far as T.
!!           | This may be caused by an inaccurate Jacobian matrix,
!!           | if one is being used.
!!           |
!!        -6 | means EWT(i) became zero for some i during the
!!           | integration.  Pure relative error control (ATOL(i)=0.0)
!!           | was requested on a variable which has now vanished.
!!           | The integration was successful as far as T.
!!           |
!!        -7 | means the length of RWORK and/or IWORK was too small to
!!           | proceed, but the integration was successful as far as T.
!!           | This happens when DLSODAR chooses to switch methods
!!           | but LRW and/or LIW is too small for the new method.
!!           |
!!
!!   Note:  Since the normal output value of ISTATE is 2,
!!   it does not need to be reset for normal continuation.
!!   Also, since a negative input value of ISTATE will be
!!   regarded as illegal, a negative output value requires the
!!   user to change it, and possibly other inputs, before
!!   calling the solver again.
!!
!! IOPT
!!
!! : an integer flag to specify whether or not any optional
!! inputs are being used on this call.  Input only.
!! The optional inputs are listed separately below.
!!
!!          IOPT = 0 means no optional inputs are being used.
!!                   Default values will be used in all cases.
!!          IOPT = 1 means one or more optional inputs are being used.
!!
!! RWORK
!!
!! : a real array (double precision) for work space, and (in the
!! first 20 words) for conditional and optional inputs and
!! optional outputs.
!!
!! As DLSODAR switches automatically between stiff and nonstiff
!! methods, the required length of RWORK can change during the
!! problem.  Thus the RWORK array passed to DLSODAR can either
!! have a static (fixed) length large enough for both methods,
!! or have a dynamic (changing) length altered by the calling
!! program in response to output from DLSODAR.
!!
!!#### Fixed Length Case
!!
!!   If the RWORK length is to be fixed, it should be at least
!!   max (LRN, LRS),
!!   where LRN and LRS are the RWORK lengths required when the
!!   current method is nonstiff or stiff, respectively.
!!
!!   The separate RWORK length requirements LRN and LRS are
!!   as follows:
!!
!!   If NEQ is constant and the maximum method orders have
!!   their default values, then
!!```text
!!             LRN = 20 + 16*NEQ + 3*NG,
!!             LRS = 22 + 9*NEQ + NEQ**2 + 3*NG           (JT = 1 or 2),
!!             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ + 3*NG   (JT = 4 or 5).
!!          Under any other conditions, LRN and LRS are given by:
!!             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ + 3*NG,
!!             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT + 3*NG,
!!          where
!!             NYH    = the initial value of NEQ,
!!             MXORDN = 12, unless a smaller value is given as an
!!                      optional input,
!!             MXORDS = 5, unless a smaller value is given as an
!!                      optional input,
!!             LMAT   = length of matrix work space:
!!             LMAT   = NEQ**2 + 2              if JT = 1 or 2,
!!             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.
!!```
!!#### Dynamic Length Case
!! If the length of RWORK is to be dynamic, then it should
!! be at least LRN or LRS, as defined above, depending on the
!! current method.  Initially, it must be at least LRN (since
!! DLSODAR starts with the nonstiff method).  On any return
!! from DLSODAR, the optional output MCUR indicates the current
!! method.  If MCUR differs from the value it had on the
!! previous return, or if there has only been one call to
!! DLSODAR and MCUR is now 2, then DLSODAR has switched
!! methods during the last call, and the length of RWORK
!! should be reset (to LRN if MCUR = 1, or to LRS if
!! MCUR = 2).  (An increase in the RWORK length is required
!! if DLSODAR returned ISTATE = -7, but not otherwise.)
!! After resetting the length, call DLSODAR with ISTATE = 3
!! to signal that change.
!!
!! LRW
!!
!! : the length of the array RWORK, as declared by the user.
!!          (This will be checked by the solver.)
!!
!! IWORK
!!
!! : an integer array for work space.
!! As DLSODAR switches automatically between stiff and nonstiff
!! methods, the required length of IWORK can change during
!! problem, between
!!             LIS = 20 + NEQ   and   LIN = 20,
!! respectively.  Thus the IWORK array passed to DLSODAR can
!! either have a fixed length of at least 20 + NEQ, or have a
!! dynamic length of at least LIN or LIS, depending on the
!! current method.  The comments on dynamic length under
!! RWORK above apply here.  Initially, this length need
!! only be at least LIN = 20.
!!
!! The first few words of IWORK are used for conditional and
!! optional inputs and optional outputs.
!!
!! The following 2 words in IWORK are conditional inputs:
!!```text
!!            IWORK(1) = ML
!!            IWORK(2) = MU
!!```
!! These are the lower and upper half-bandwidths, respectively, of the
!! banded Jacobian, excluding the main diagonal.  The band is defined by
!! the matrix locations (i,j) with i-ML .le. j .le. i+MU.  ML and MU must
!! satisfy  0 .le.  ML,MU  .le. NEQ-1.  These are required if JT is 4 or 5,
!! and ignored otherwise.  ML and MU may in fact be the band parameters
!! for a matrix to which df/dy is only approximately equal.
!!
!! LIW
!!
!! : the length of the array IWORK, as declared by the user.
!! (This will be checked by the solver.)
!!
!! **Note**: The base addresses of the work arrays must not be
!! altered between calls to DLSODAR for the same problem.
!! The contents of the work arrays must not be altered
!! between calls, except possibly for the conditional and
!! optional inputs, and except for the last 3\*NEQ words of RWORK.
!! The latter space is used for internal scratch space, and so is
!! available for use by the user outside DLSODAR between calls, if
!! desired (but not for use by F, JAC, or G).
!!
!! JAC
!!
!! : the name of the user-supplied routine to compute the
!! Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine
!! is optional, but if the problem is expected to be stiff much
!! of the time, you are encouraged to supply JAC, for the sake
!! of efficiency.  (Alternatively, set JT = 2 or 5 to have
!! DLSODAR compute df/dy internally by difference quotients.)
!! If and when DLSODAR uses df/dy, it treats this NEQ by NEQ
!! matrix either as full (JT = 1 or 2), or as banded (JT =
!! 4 or 5) with half-bandwidths ML and MU (discussed under
!! IWORK above).  In either case, if JT = 1 or 4, the JAC
!! routine must compute df/dy as a function of the scalar t
!! and the vector y.  It is to have the form
!!```fortran
!!       SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
!!       DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
!!```
!! where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!! PD is to be loaded with partial derivatives (elements of
!! the Jacobian matrix) on output.  PD must be given a first
!! dimension of NROWPD.  T and Y have the same meaning as in
!! Subroutine F.
!!
!! In the full matrix case (JT = 1), ML and MU are
!! ignored, and the Jacobian is to be loaded into PD in
!! columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
!!
!! In the band matrix case (JT = 4), the elements
!! within the band are to be loaded into PD in columnwise
!! manner, with diagonal lines of df/dy loaded into the rows
!! of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!! ML and MU are the half-bandwidth parameters (see IWORK).
!! The locations in PD in the two triangular areas which
!! correspond to nonexistent matrix elements can be ignored
!! or loaded arbitrarily, as they are overwritten by DLSODAR.
!!
!! JAC need not provide df/dy exactly.  A crude
!! approximation (possibly with a smaller bandwidth) will do.
!!
!! In either case, PD is preset to zero by the solver,
!! so that only the nonzero elements need be loaded by JAC.
!! Each call to JAC is preceded by a call to F with the same
!! arguments NEQ, T, and Y.  Thus to gain some efficiency,
!! intermediate quantities shared by both calculations may be
!! saved in a user Common block by F and not recomputed by JAC,
!! if desired.  Also, JAC may alter the Y array, if desired.
!! JAC must be declared External in the calling program.
!!
!! Subroutine JAC may access user-defined quantities in
!! NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!! (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
!! See the descriptions of NEQ and Y above.
!!
!! JT
!!
!! : Jacobian type indicator.  Used only for input.
!!
!! JT specifies how the Jacobian matrix df/dy will be
!! treated, if and when DLSODAR requires this matrix.
!! JT has the following values and meanings:
!!
!! value | description
!! ----- | --------------------------------------------------------------
!!    1  | means a user-supplied full (NEQ by NEQ) Jacobian.
!!    2  | means an internally generated (difference quotient) full
!!       | Jacobian (using NEQ extra calls to F per df/dy value).
!!    4  | means a user-supplied banded Jacobian.
!!    5  | means an internally generated banded Jacobian (using
!!       | ML+MU+1 extra calls to F per df/dy evaluation).
!!
!! If JT = 1 or 4, the user must supply a Subroutine JAC
!! (the name is arbitrary) as described above under JAC.
!! If JT = 2 or 5, a dummy argument can be used.
!!
!! G
!!
!! : the name of subroutine for constraint functions, whose
!! roots are desired during the integration.  It is to have
!! the form
!!```fortran
!!       SUBROUTINE G (NEQ, T, Y, NG, GOUT)
!!       DOUBLE PRECISION T, Y(*), GOUT(NG)
!!```
!! where NEQ, T, Y, and NG are input, and the array GOUT
!! is output.  NEQ, T, and Y have the same meaning as in
!! the F routine, and GOUT is an array of length NG.
!! For i = 1,...,NG, this routine is to load into GOUT(i)
!! the value at (T,Y) of the i-th constraint function g(i).
!! DLSODAR will find roots of the g(i) of odd multiplicity
!! (i.e. sign changes) as they occur during the integration.
!! G must be declared External in the calling program.
!!
!! Caution:  Because of numerical errors in the functions
!! g(i) due to roundoff and integration error, DLSODAR may
!! return false roots, or return the same root at two or more
!! nearly equal values of t.  If such false roots are
!! suspected, the user should consider smaller error tolerances
!! and/or higher precision in the evaluation of the g(i).
!!
!! If a root of some g(i) defines the end of the problem,
!! the input to DLSODAR should nevertheless allow integration
!! to a point slightly past that root, so that DLSODAR can
!! locate the root by interpolation.
!!
!! Subroutine G may access user-defined quantities in
!! NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
!! (dimensioned in G) and/or Y has length exceeding NEQ(1).
!! See the descriptions of NEQ and Y above.
!!
!! NG
!!
!! : number of constraint functions g(i).  If there are none,
!! set NG = 0, and pass a dummy name for G.
!!
!! JROOT
!!
!! : integer array of length NG.  Used only for output.
!! On a return with ISTATE = 3 (one or more roots found),
!! JROOT(i) = 1 if g(i) has a root at T, or JROOT(i) = 0 if not.
!!
!!-----------------------------------------------------------------------
!!
!!### Optional Inputs.
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
!! Name  | Location |     Meaning and Default Value
!! ----- | -------- | ---------------------------------------------------
!! H0    | RWORK(5) | the step size to be attempted on the first step.
!!       |          | The default value is determined by the solver.
!!       |          |
!! HMAX  | RWORK(6) | the maximum absolute step size allowed.
!!       |          | The default value is infinite.
!!       |          |
!! HMIN  | RWORK(7) | the minimum absolute step size allowed.
!!       |          | The default value is 0.  (This lower bound is not
!!       |          | enforced on the final step before reaching TCRIT
!!       |          | when ITASK = 4 or 5.)
!!       |          |
!! IXPR  | IWORK(5) | flag to generate extra printing at method switches.
!!       |          | IXPR = 0 means no extra printing (the default).
!!       |          | IXPR = 1 means print data on each switch.
!!       |          | T, H, and NST will be printed on the same logical
!!       |          | unit as used for error messages.
!!       |          |
!! MXSTEP| IWORK(6) | maximum number of (internally defined) steps
!!       |          | allowed during one call to the solver.
!!       |          | The default value is 500.
!!       |          |
!! MXHNIL| IWORK(7) | maximum number of messages printed (per problem)
!!       |          | warning that T + H = T on a step (H = step size).
!!       |          | This must be positive to result in a non-default
!!       |          | value.  The default value is 10.
!!       |          |
!! MXORDN| IWORK(8) | the maximum order to be allowed for the nonstiff
!!       |          | (Adams) method.  The default value is 12.
!!       |          | If MXORDN exceeds the default value, it will
!!       |          | be reduced to the default value.
!!       |          | MXORDN is held constant during the problem.
!!       |          |
!! MXORDS| IWORK(9) | the maximum order to be allowed for the stiff
!!       |          | (BDF) method.  The default value is 5.
!!       |          | If MXORDS exceeds the default value, it will
!!       |          | be reduced to the default value.
!!       |          | MXORDS is held constant during the problem.
!!
!!-----------------------------------------------------------------------
!!
!!### Optional Outputs.
!!
!! As optional additional output from DLSODAR, the variables listed
!! below are quantities related to the performance of DLSODAR
!! which are available to the user.  These are communicated by way of
!! the work arrays, but also have internal mnemonic names as shown.
!!
!! Except where stated otherwise, all of these outputs are defined
!! on any successful return from DLSODAR, and on any return with
!! ISTATE = -1, -2, -4, -5, or -6.
!!
!! On an illegal input return
!! (ISTATE = -3), they will be unchanged from their existing values
!! (if any), except possibly for TOLSF, LENRW, and LENIW.
!!
!! On any error return, outputs relevant to the error will be defined,
!! as noted below.
!!
!! Name  | Location |     Meaning
!! ----- |--------- | ---------------------------------------------------
!! HU    | RWORK(11)| the step size in t last used (successfully).
!!       |          |
!! HCUR  | RWORK(12)| the step size to be attempted on the next step.
!!       |          |
!! TCUR  | RWORK(13)| the current value of the independent variable
!!       |          | which the solver has actually reached, i.e. the
!!       |          | current internal mesh point in t.  On output, TCUR
!!       |          | will always be at least as far as the argument
!!       |          | T, but may be farther (if interpolation was done).
!!       |          |
!! TOLSF | RWORK(14)| a tolerance scale factor, greater than 1.0,
!!       |          | computed when a request for too much accuracy was
!!       |          | detected (ISTATE = -3 if detected at the start of
!!       |          | the problem, ISTATE = -2 otherwise).  If ITOL is
!!       |          | left unaltered but RTOL and ATOL are uniformly
!!       |          | scaled up by a factor of TOLSF for the next call,
!!       |          | then the solver is deemed likely to succeed.
!!       |          | (The user may also ignore TOLSF and alter the
!!       |          | tolerance parameters in any other way appropriate.)
!!       |          |
!! TSW   | RWORK(15)| the value of t at the time of the last method
!!       |          | switch, if any.
!!       |          |
!! NGE   | IWORK(10)| the number of g evaluations for the problem so far.
!!       |          |
!! NST   | IWORK(11)| the number of steps taken for the problem so far.
!!       |          |
!! NFE   | IWORK(12)| the number of f evaluations for the problem so far.
!!       |          |
!! NJE   | IWORK(13)| the number of Jacobian evaluations (and of matrix
!!       |          | LU decompositions) for the problem so far.
!!       |          |
!! NQU   | IWORK(14)| the method order last used (successfully).
!!       |          |
!! NQCUR | IWORK(15)| the order to be attempted on the next step.
!!       |          |
!! IMXER | IWORK(16)| the index of the component of largest magnitude in
!!       |          | the weighted local error vector ( E(i)/EWT(i) ),
!!       |          | on an error return with ISTATE = -4 or -5.
!!       |          |
!! LENRW | IWORK(17)| the length of RWORK actually required, assuming
!!       |          | that the length of RWORK is to be fixed for the
!!       |          | rest of the problem, and that switching may occur.
!!       |          | This is defined on normal returns and on an illegal
!!       |          | input return for insufficient storage.
!!       |          |
!! LENIW | IWORK(18)| the length of IWORK actually required, assuming
!!       |          | that the length of IWORK is to be fixed for the
!!       |          | rest of the problem, and that switching may occur.
!!       |          | This is defined on normal returns and on an illegal
!!       |          | input return for insufficient storage.
!!       |          |
!! MUSED | IWORK(19)| the method indicator for the last successful step:
!!       |          | 1 means Adams (nonstiff), 2 means BDF (stiff).
!!       |          |
!! MCUR  | IWORK(20)| the current method indicator:
!!       |          | 1 means Adams (nonstiff), 2 means BDF (stiff).
!!       |          | This is the method to be attempted
!!       |          | on the next step.  Thus it differs from MUSED
!!       |          | only if a method switch has just been made.
!!
!! The following two arrays are segments of the RWORK array which
!! may also be of interest to the user as optional outputs.
!! For each array, the table below gives its internal name,
!! its base address in RWORK, and its description.
!!
!! Name  | Base Address |    Description
!! ----- | ------------ | ---------------------------------------------------
!!       |              |
!! YH    | 21 + 3\*NG   |  the Nordsieck history array, of size NYH by
!!       |              | (NQCUR + 1), where NYH is the initial value
!!       |              | of NEQ.  For j = 0,1,...,NQCUR, column j+1
!!       |              | of YH contains HCUR\*\*j/factorial(j) times
!!       |              | the j-th derivative of the interpolating
!!       |              | polynomial currently representing the solution,
!!       |              | evaluated at t = TCUR.
!!       |              |
!! ACOR  |  LACOR       | array of size NEQ used for the accumulated
!!       | (from Common | corrections on each step, scaled on output
!!       |   as noted)  | to represent the estimated local error in y
!!       |              | on the last step.  This is the vector E in
!!       |              | the description of the error control.  It is
!!       |              | defined only on a successful return from
!!       |              | DLSODAR.  The base address LACOR is obtained by
!!       |              | including in the user's program the
!!       |              | following 2 lines:
!!       |              |    COMMON /DLS001/ RLS(218), ILS(37)
!!       |              |    LACOR = ILS(22)
!!
!!### Part 2.  Other Routines Callable.
!!
!! The following are optional calls which the user may make to
!! gain additional capabilities in conjunction with DLSODAR.
!! (The routines XSETUN and XSETF are designed to conform to the
!! SLATEC error handling package.)
!!
!! Form of Call                | Function
!! --------------------------- | -----------------------------------------
!!  CALL XSETUN(LUN)           | Set the logical unit number, LUN, for
!!                             | output of messages from DLSODAR, if
!!                             | the default is not desired.
!!                             | The default value of LUN is 6.
!!                             |
!!  CALL XSETF(MFLAG)          | Set a flag to control the printing of
!!                             | messages by DLSODAR.
!!                             | MFLAG = 0 means do not print. (Danger:
!!                             | This risks losing valuable information.)
!!                             | MFLAG = 1 means print (the default).
!!                             |
!!                             | Either of the above calls may be made at
!!                             | any time and will take effect immediately.
!!                             |
!!  CALL DSRCAR(RSAV,ISAV,JOB) | saves and restores the contents of
!!                             | the internal Common blocks used by
!!                             | DLSODAR (see Part 3 below).
!!                             | RSAV must be a real array of length 245
!!                             | or more, and ISAV must be an integer
!!                             | array of length 55 or more.
!!                             | JOB=1 means save Common into RSAV/ISAV.
!!                             | JOB=2 means restore Common from RSAV/ISAV.
!!                             |
!!                             |    DSRCAR is useful if one is
!!                             | interrupting a run and restarting
!!                             | later, or alternating between two or
!!                             | more problems solved with DLSODAR.
!!                             |
!!  CALL DINTDY(,,,,,)         | Provide derivatives of y, of various
!!       (see below)           | orders, at a specified point t, if
!!                             | desired.  It may be called only after
!!                             | a successful return from DLSODAR.
!!
!! The detailed instructions for using DINTDY are as follows.
!! The form of the call is:
!!```fortran
!!   LYH = 21 + 3*NG
!!   CALL DINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
!!```
!! The input parameters are:
!!
!! T
!!
!! : value of independent variable where answers are desired
!! (normally the same as the T last returned by DLSODAR).
!! For valid results, T must lie between TCUR - HU and TCUR.
!! (See optional outputs for TCUR and HU.)
!!
!! K
!!
!! : integer order of the derivative desired.  K must satisfy
!! 0 .le. K .le. NQCUR, where NQCUR is the current order
!! (see optional outputs).  The capability corresponding
!! to K = 0, i.e. computing y(t), is already provided
!! by DLSODAR directly.  Since NQCUR .ge. 1, the first
!! derivative dy/dt is always available with DINTDY.
!!
!! LYH
!!
!! : 21 + 3\*NG = base address in RWORK of the history array YH.
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
!!
!! IFLAG
!!
!! : integer flag, returned as 0 if K and T were legal,
!! -1 if K was illegal, and -2 if T was illegal.
!! On an error return, a message is also written.
!!-----------------------------------------------------------------------
!!### Part 3.  Save and Restore Program State
!!
!! If the solution of a given problem by DLSODAR is to be interrupted
!! and then later continued, such as when restarting an interrupted run
!! or alternating between two or more problems, the user should save,
!! following the return from the last DLSODAR call prior to the
!! interruption, the contents of the call sequence variables and the
!! internal state variables, and later restore these values before the
!! next DLSODAR call for that problem.  To save and restore,
!! use Subroutine DSRCAR (see Part 2 above).
!!
!!-----------------------------------------------------------------------
!!### Part 4.  Optionally Replaceable Solver Routines.
!!
!! Below is a description of a routine in the DLSODAR package which
!! relates to the measurement of errors, and can be
!! replaced by a user-supplied version, if desired.  However, since such
!! a replacement may have a major impact on performance, it should be
!! done only when absolutely necessary, and only with great caution.
!! (Note: The means by which the package version of a routine is
!! superseded by the user's version may be system-dependent.)
!!
!! (a) DEWSET.
!! The following subroutine is called just before each internal
!! integration step, and sets the array of error weights, EWT, as
!! described under ITOL/RTOL/ATOL above:
!!     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
!! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODAR call sequence,
!! YCUR contains the current dependent variable vector, and
!! EWT is the array of weights set by DEWSET.
!!
!! If the user supplies this subroutine, it must return in EWT(i)
!! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
!! in y(i) to.  The EWT array returned by DEWSET is passed to the
!! DMNORM routine, and also used by DLSODAR in the computation
!! of the optional output IMXER, and the increments for difference
!! quotient Jacobians.
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
!!     DOUBLE PRECISION RLS
!!     COMMON /DLS001/ RLS(218),ILS(37)
!!     NQ = ILS(33)
!!     NST = ILS(34)
!!     H = RLS(212)
!! Thus, for example, the current value of dy/dt can be obtained as
!! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
!! unnecessary when NST = 0).
!!-----------------------------------------------------------------------
!!### Pedigree:
!!
!! _DLSODAR_ is derived from the Livermore Solver for Ordinary
!! Differential Equations package ODEPACK, and is in double precision.
!!
!!### Authors
!!                Alan C. Hindmarsh,
!!                Center for Applied Scientific Computing, L-561
!!                Lawrence Livermore National Laboratory
!!                Livermore, CA 94551
!! and
!!                Linda R. Petzold
!!                Univ. of California at Santa Barbara
!!                Dept. of Computer Science
!!                Santa Barbara, CA 93106
!!
!!### References
!! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!!     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
!!     North-Holland, Amsterdam, 1983, pp. 55-64.
!! 2.  Linda R. Petzold, Automatic Selection of Methods for Solving
!!     Stiff and Nonstiff Systems of Ordinary Differential Equations,
!!     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
!! 3.  Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
!!     Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
!!     February 1980.
!!-----------------------------------------------------------------------
!! Other routines in the DLSODAR package.
!!
!! In addition to Subroutine DLSODAR, the DLSODAR package includes the
!! following subroutines and function routines:
!!
!!  DRCHEK
!!
!!  : does preliminary checking for roots, and serves as an
!!  interface between Subroutine DLSODAR and Subroutine DROOTS.
!!
!!  DROOTS
!!
!!  : finds the leftmost root of a set of functions.
!!
!!  DINTDY
!!
!!  : computes an interpolated value of the y vector at t = TOUT.
!!
!!  DSTODA
!!
!!  : is the core integrator, which does one step of the
!!  integration and the associated error control.
!!
!!  DCFODE
!!
!!  : sets all method coefficients and test constants.
!!
!!  DPRJA
!!
!!  : computes and preprocesses the Jacobian matrix J = df/dy
!!  and the Newton iteration matrix P = I - h\*l0\*J.
!!
!!  DSOLSY
!!
!!  : manages solution of linear system in chord iteration.
!!
!!  DEWSET
!!
!!  : sets the error weight vector EWT before each step.
!!
!!  DMNORM
!!
!!  : computes the weighted max-norm of a vector.
!!
!!  DFNORM
!!
!!  : computes the norm of a full matrix consistent with the
!!  weighted max-norm on vectors.
!!
!!  DBNORM
!!
!!  : computes the norm of a band matrix consistent with the
!!  weighted max-norm on vectors.
!!
!!  DSRCAR
!!
!!  : is a user-callable routine to save and restore
!!  the contents of the internal Common blocks.
!!
!!  DGEFA and DGESL
!!
!!  : are routines from LINPACK for solving full
!!  systems of linear algebraic equations.
!!
!!  DGBFA and DGBSL
!!
!!  : are routines from LINPACK for solving banded
!!  linear systems.
!!  DCOPY
!!
!!  : is one of the basic linear algebra modules (BLAS).
!!
!!  XERRWD, XSETUN, XSETF, and IXSAV
!!
!!  : handle the printing of all
!!  error messages and warnings.  XERRWD is machine-dependent.
!!
!! Note:  DMNORM, DFNORM, DBNORM, and IXSAV are
!! function routines.  All the others are subroutines.
!!
!!-----------------------------------------------------------------------
!
! ### REVISION HISTORY  (YYYYMMDD)
!  19811102  DATE WRITTEN
!  19820126  Fixed bug in tests of work space lengths;
!            minor corrections in main prologue and comments.
!  19820507  Fixed bug in RCHEK in setting HMING.
!  19870330  Major update: corrected comments throughout;
!            removed TRET from Common; rewrote EWSET with 4 loops;
!            fixed t test in INTDY; added Cray directives in STODA;
!            in STODA, fixed DELP init. and logic around PJAC call;
!            combined routines to save/restore Common;
!            passed LEVEL = 0 in error message calls (except run abort).
!  19970225  Fixed lines setting JSTART = -2 in Subroutine LSODAR.
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
!            numerous corrections to prologues and internal comments.
!  20010507  Converted single precision source to double precision.
!  20010613  Revised excess accuracy test (to match rest of ODEPACK).
!  20010808  Fixed bug in DPRJA (matrix in DBNORM call).
!  20020502  Corrected declarations in descriptions of user routines.
!  20031105  Restored 'own' variables to Common blocks, to enable
!            interrupt/restart feature.
!  20031112  Added SAVE statements for data-loaded constants.
!
!!-----------------------------------------------------------------------
subroutine dlsodar(f,Neq,Y,T,Tout,Itol,Rtol,Atol,Itask,Istate,Iopt,Rwork,Lrw,Iwork,Liw,jac,Jt,g,Ng,Jroot)
!
external f
external g
external jac
real(kind=dp), dimension(*) :: Atol, Rtol, Y
real(kind=dp) :: atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, size, sum, tcrit, tdist, tnext, tol, tolsf, tp, w0
integer :: i, i1, i2, iflag, imxer, irfp, irt, kgo, len1, len1c, len1n, len1s, len2, leniw, leniwc, lenrw, lenrwc, &
         & lenwm, lenyh, lf0, lyhnew, ml, mu
logical :: ihit
integer :: Iopt, Istate, Itask, Itol, Jt, Liw, Lrw, Ng
integer, intent(inout), dimension(Liw) :: Iwork
integer, dimension(Ng) :: Jroot
integer, dimension(2), save :: mord
character(60) :: msg
integer, save :: mxhnl0, mxstp0
integer, dimension(*) :: Neq
real(kind=dp), intent(inout), dimension(Lrw) :: Rwork
real(kind=dp), intent(inout) :: T, Tout
!
! -----------------------------------------------------------------------
!  The following three internal Common blocks contain
!  (a) variables which are local to any subroutine but whose values must
!      be preserved between calls to the routine ("own" variables), and
!  (b) variables which are communicated between subroutines.
!  The block DLS001 is declared in subroutines DLSODAR, DINTDY, DSTODA,
!  DPRJA, and DSOLSY.
!  The type(DLSA01)::DLSA is declared in subroutines DLSODAR, DSTODA, DPRJA.
!  The block DLSR01 is declared in subroutines DLSODAR, DRCHEK, DROOTS.
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
!  If ISTATE = 1 and TOUT = T, return immediately.
! -----------------------------------------------------------------------
if(Ng.gt.0) Jroot(ng)=0

if ( Istate<1 .or. Istate>3 ) then
! -----------------------------------------------------------------------
!  Block I.
!  The following block handles all error returns due to illegal input
!  (ISTATE = -3), as detected before calling the core integrator.
!  First the error message routine is called.  If the illegal input
!  is a negative ISTATE, the run is aborted (apparent infinite loop).
! -----------------------------------------------------------------------
   msg = 'DLSODAR-  ISTATE(=I1) illegal.'
   call xerrwd(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
   if ( Istate>=0 ) goto 1200
!
   msg = 'DLSODAR-  Run aborted.. apparent infinite loop.   '
   call xerrwd(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
   goto 99999
else
   if ( Itask<1 .or. Itask>5 ) then
      msg = 'DLSODAR-  ITASK (=I1) illegal.'
      call xerrwd(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
      goto 1200
   else
      dlsr%itaskc = Itask
      if ( Istate==1 ) then
         dls1%init = 0
         if ( Tout==T ) return
      elseif ( dls1%init==0 ) then
         msg = 'DLSODAR-  ISTATE.gt.1 but DLSODAR not initialized.'
         call xerrwd(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
         goto 1200
      elseif ( Istate==2 ) then
         goto 50
      endif
! -----------------------------------------------------------------------
!  Block B.
!  The next code block is executed for the initial call (ISTATE = 1),
!  or for a continuation call with parameter changes (ISTATE = 3).
!  It contains checking of all inputs and various initializations.
!
!  First check legality of the non-optional inputs NEQ, ITOL, IOPT,
!  JT, ML, MU, and NG.
! -----------------------------------------------------------------------
      if ( Neq(1)<=0 ) then
         msg = 'DLSODAR-  NEQ (=I1) .lt. 1    '
         call xerrwd(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
         goto 1200
      else
         if ( Istate/=1 ) then
            if ( Neq(1)>dls1%n ) then
               msg = 'DLSODAR-  ISTATE = 3 and NEQ increased (I1 to I2).'
               call xerrwd(msg,50,5,0,2,dls1%n,Neq(1),0,0.0D0,0.0D0)
               goto 1200
            endif
         endif
         dls1%n = Neq(1)
         if ( Itol<1 .or. Itol>4 ) then
            msg = 'DLSODAR-  ITOL (=I1) illegal. '
            call xerrwd(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
            goto 1200
         elseif ( Iopt<0 .or. Iopt>1 ) then
            msg = 'DLSODAR-  IOPT (=I1) illegal. '
            call xerrwd(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
            goto 1200
         elseif ( Jt==3 .or. Jt<1 .or. Jt>5 ) then
            msg = 'DLSODAR-  JT (=I1) illegal.   '
            call xerrwd(msg,30,8,0,1,Jt,0,0,0.0D0,0.0D0)
            goto 1200
         else
            dlsa%jtyp = Jt
            if ( Jt>2 ) then
               ml = Iwork(1)
               mu = Iwork(2)
               if ( ml<0 .or. ml>=dls1%n ) then
                  msg = 'DLSODAR-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2)'
                  call xerrwd(msg,50,9,0,2,ml,Neq(1),0,0.0D0,0.0D0)
                  goto 1200
               elseif ( mu<0 .or. mu>=dls1%n ) then
                  msg = 'DLSODAR-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2)'
                  call xerrwd(msg,50,10,0,2,mu,Neq(1),0,0.0D0,0.0D0)
                  goto 1200
               endif
            endif
            if ( Ng<0 ) then
               msg = 'DLSODAR-  NG (=I1) .lt. 0     '
               call xerrwd(msg,30,30,0,1,Ng,0,0,0.0D0,0.0D0)
               goto 1200
            else
               if ( Istate/=1 ) then
                  if ( dlsr%irfnd==0 .and. Ng/=dlsr%ngc ) then
                     msg = 'DLSODAR-  NG changed (from I1 to I2) illegally,   '
                     call xerrwd(msg,50,31,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      i.e. not immediately after a root was found.'
                     call xerrwd(msg,50,31,0,2,dlsr%ngc,Ng,0,0.0D0,0.0D0)
                     goto 1200
                  endif
               endif
               dlsr%ngc = Ng
!  Next process and check the optional inputs. --------------------------
               if ( Iopt==1 ) then
                  dlsa%ixpr = Iwork(5)
                  if ( dlsa%ixpr<0 .or. dlsa%ixpr>1 ) then
                     msg = 'DLSODAR-  IXPR (=I1) illegal. '
                     call xerrwd(msg,30,11,0,1,dlsa%ixpr,0,0,0.0D0,0.0D0)
                     goto 1200
                  else
                     dls1%mxstep = Iwork(6)
                     if ( dls1%mxstep<0 ) then
                        msg = 'DLSODAR-  MXSTEP (=I1) .lt. 0 '
                        call xerrwd(msg,30,12,0,1,dls1%mxstep,0,0,0.0D0,0.0D0)
                        goto 1200
                     else
                        if ( dls1%mxstep==0 ) dls1%mxstep = mxstp0
                        dls1%mxhnil = Iwork(7)
                        if ( dls1%mxhnil<0 ) then
                           msg = 'DLSODAR-  MXHNIL (=I1) .lt. 0 '
                           call xerrwd(msg,30,13,0,1,dls1%mxhnil,0,0,0.0D0,0.0D0)
                           goto 1200
                        else
                           if ( dls1%mxhnil==0 ) dls1%mxhnil = mxhnl0
                           if ( Istate==1 ) then
                              h0 = Rwork(5)
                              dlsa%mxordn = Iwork(8)
                              if ( dlsa%mxordn<0 ) then
                                 msg = 'DLSODAR-  MXORDN (=I1) .lt. 0 '
                                 call xerrwd(msg,30,28,0,1,dlsa%mxordn,0,0,0.0D0,0.0D0)
                                 goto 1200
                              else
                                 if ( dlsa%mxordn==0 ) dlsa%mxordn = 100
                                 dlsa%mxordn = min(dlsa%mxordn,mord(1))
                                 dlsa%mxords = Iwork(9)
                                 if ( dlsa%mxords<0 ) then
                                    msg = 'DLSODAR-  MXORDS (=I1) .lt. 0 '
                                    call xerrwd(msg,30,29,0,1,dlsa%mxords,0,0,0.0D0,0.0D0)
                                    goto 1200
                                 else
                                    if ( dlsa%mxords==0 ) dlsa%mxords = 100
                                    dlsa%mxords = min(dlsa%mxords,mord(2))
                                    if ( (Tout-T)*h0<0.0D0 ) then
                                       msg = 'DLSODAR-  TOUT (=R1) behind T (=R2)     '
                                       call xerrwd(msg,40,14,0,0,0,0,2,Tout,T)
                                       msg = '      Integration direction is given by H0 (=R1)  '
                                       call xerrwd(msg,50,14,0,0,0,0,1,h0,0.0D0)
                                       goto 1200
                                    endif
                                 endif
                              endif
                           endif
                           hmax = Rwork(6)
                           if ( hmax<0.0D0 ) then
                              msg = 'DLSODAR-  HMAX (=R1) .lt. 0.0 '
                              call xerrwd(msg,30,15,0,0,0,0,1,hmax,0.0D0)
                              goto 1200
                           else
                              dls1%hmxi = 0.0D0
                              if ( hmax>0.0D0 ) dls1%hmxi = 1.0D0/hmax
                              dls1%hmin = Rwork(7)
                              if ( dls1%hmin<0.0D0 ) then
                                 msg = 'DLSODAR-  HMIN (=R1) .lt. 0.0 '
                                 call xerrwd(msg,30,16,0,0,0,0,1,dls1%hmin,0.0D0)
                                 goto 1200
                              endif
                           endif
                        endif
                     endif
                  endif
               else
                  dlsa%ixpr = 0
                  dls1%mxstep = mxstp0
                  dls1%mxhnil = mxhnl0
                  dls1%hmxi = 0.0D0
                  dls1%hmin = 0.0D0
                  if ( Istate==1 ) then
                     h0 = 0.0D0
                     dlsa%mxordn = mord(1)
                     dlsa%mxords = mord(2)
                  endif
               endif
! -----------------------------------------------------------------------
!  Set work array pointers and check lengths LRW and LIW.
!  If ISTATE = 1, METH is initialized to 1 here to facilitate the
!  checking of work space lengths.
!  Pointers to segments of RWORK and IWORK are named by prefixing L to
!  the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
!  Segments of RWORK (in order) are denoted  G0, G1, GX, YH, WM,
!  EWT, SAVF, ACOR.
!  If the lengths provided are insufficient for the current method,
!  an error return occurs.  This is treated as illegal input on the
!  first call, but as a problem interruption with ISTATE = -7 on a
!  continuation call.  If the lengths are sufficient for the current
!  method but not for both methods, a warning message is sent.
! -----------------------------------------------------------------------
               if ( Istate==1 ) dls1%meth = 1
               if ( Istate==1 ) dls1%nyh = dls1%n
               dlsr%lg0 = 21
               dlsr%lg1 = dlsr%lg0 + Ng
               dlsr%lgx = dlsr%lg1 + Ng
               lyhnew = dlsr%lgx + Ng
               if ( Istate==1 ) dls1%lyh = lyhnew
               if ( lyhnew/=dls1%lyh ) then
!  If ISTATE = 3 and NG was changed, shift YH to its new location. ------
                  lenyh = dls1%l*dls1%nyh
                  if ( Lrw>=lyhnew-1+lenyh ) then
                     i1 = 1
                     if ( lyhnew>dls1%lyh ) i1 = -1
                     !X!call dcopy(lenyh,Rwork(dls1%lyh),i1,Rwork(lyhnew),i1)
                     Rwork(lyhnew:lyhnew+lenyh-1:i1)=Rwork(dls1%lyh:dls1%lyh+lenyh-1:i1) !X!
                     dls1%lyh = lyhnew
                  endif
               endif
               len1n = lyhnew - 1 + (dlsa%mxordn+1)*dls1%nyh
               len1s = lyhnew - 1 + (dlsa%mxords+1)*dls1%nyh
               dls1%lwm = len1s + 1

               lenwm=0
               if ( Jt<=2 ) lenwm = dls1%n*dls1%n + 2
               if ( Jt>=4 ) lenwm = (2*ml+mu+1)*dls1%n + 2

               len1s = len1s + lenwm
               len1c = len1n
               if ( dls1%meth==2 ) len1c = len1s
               len1 = max(len1n,len1s)
               len2 = 3*dls1%n
               lenrw = len1 + len2
               lenrwc = len1c + len2
               Iwork(17) = lenrw
               dls1%liwm = 1
               leniw = 20 + dls1%n
               leniwc = 20
               if ( dls1%meth==2 ) leniwc = leniw
               Iwork(18) = leniw
               if ( Istate==1 .and. Lrw<lenrwc ) then
                  msg = 'DLSODAR-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) '
                  call xerrwd(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                  goto 1200
               elseif ( Istate==1 .and. Liw<leniwc ) then
                  msg = 'DLSODAR-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) '
                  call xerrwd(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                  goto 1200
               else
                  if ( Istate==3 .and. Lrw<lenrwc ) goto 500
                  if ( Istate==3 .and. Liw<leniwc ) goto 600
                  dls1%lewt = len1 + 1
                  dlsa%insufr = 0
                  if ( Lrw<lenrw ) then
                     dlsa%insufr = 2
                     dls1%lewt = len1c + 1
                     msg = 'DLSODAR-  Warning.. RWORK length is sufficient for now, but '
                     call xerrwd(msg,60,103,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      may not be later.  Integration will proceed anyway.   '
                     call xerrwd(msg,60,103,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      Length needed is LENRW = I1, while LRW = I2.'
                     call xerrwd(msg,50,103,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                  endif
                  dls1%lsavf = dls1%lewt + dls1%n
                  dls1%lacor = dls1%lsavf + dls1%n
                  dlsa%insufi = 0
                  if ( Liw<leniw ) then
                     dlsa%insufi = 2
                     msg = 'DLSODAR-  Warning.. IWORK length is sufficient for now, but '
                     call xerrwd(msg,60,104,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      may not be later.  Integration will proceed anyway.   '
                     call xerrwd(msg,60,104,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      Length needed is LENIW = I1, while LIW = I2.'
                     call xerrwd(msg,50,104,0,2,leniw,Liw,0,0.0D0,0.0D0)
                  endif
!  Check RTOL and ATOL for legality. ------------------------------------
                  rtoli = Rtol(1)
                  atoli = Atol(1)
                  do i = 1, dls1%n
                     if ( Itol>=3 ) rtoli = Rtol(i)
                     if ( Itol==2 .or. Itol==4 ) atoli = Atol(i)
                     if ( rtoli<0.0D0 ) then
                        msg = 'DLSODAR-  RTOL(I1) is R1 .lt. 0.0       '
                        call xerrwd(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
                        goto 1200
                     elseif ( atoli<0.0D0 ) then
                        msg = 'DLSODAR-  ATOL(I1) is R1 .lt. 0.0       '
                        call xerrwd(msg,40,20,0,1,i,0,1,atoli,0.0D0)
                        goto 1200
                     endif
                  enddo
                  if ( Istate==1 ) then
! -----------------------------------------------------------------------
!  Block C.
!  The next block is for the initial call only (ISTATE = 1).
!  It contains all remaining initializations, the initial call to F,
!  and the calculation of the initial step size.
!  The error weights in EWT are inverted after being loaded.
! -----------------------------------------------------------------------
                     dls1%uround = epsilon(0.0d0)
                     dls1%tn = T
                     dlsa%tsw = T
                     dls1%maxord = dlsa%mxordn
                     if ( Itask==4 .or. Itask==5 ) then
                        tcrit = Rwork(1)
                        if ( (tcrit-Tout)*(Tout-T)<0.0D0 ) goto 1000
                        if ( h0/=0.0D0 .and. (T+h0-tcrit)*h0>0.0D0 ) h0 = tcrit - T
                     endif
                     dls1%jstart = 0
                     dls1%nhnil = 0
                     dls1%nst = 0
                     dls1%nje = 0
                     dls1%nslast = 0
                     dls1%hu = 0.0D0
                     dls1%nqu = 0
                     dlsa%mused = 0
                     dls1%miter = 0
                     dls1%ccmax = 0.3D0
                     dls1%maxcor = 3
                     dls1%msbp = 20
                     dls1%mxncf = 10
!  Initial call to F.  (LF0 points to YH(*,2).) -------------------------
                     lf0 = dls1%lyh + dls1%nyh
                     call f(Neq,T,Y,Rwork(lf0))
                     dls1%nfe = 1
!  Load the initial value vector in YH. ---------------------------------
                     do i = 1, dls1%n
                        Rwork(i+dls1%lyh-1) = Y(i)
                     enddo
!  Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                     dls1%nq = 1
                     dls1%h = 1.0D0
                     call dewset(dls1%n,Itol,Rtol,Atol,Rwork(dls1%lyh),Rwork(dls1%lewt))
                     do i = 1, dls1%n
                        if ( Rwork(i+dls1%lewt-1)<=0.0D0 ) then
                           ewti = Rwork(dls1%lewt+i-1)
                           msg = 'DLSODAR-  EWT(I1) is R1 .le. 0.0        '
                           call xerrwd(msg,40,21,0,1,i,0,1,ewti,0.0D0)
                           goto 1200
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
!  Then the computed value H0 is given by:
!
!    H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
!
!  where   w0     = MAX ( ABS(T), ABS(TOUT) ),
!          F      = the initial value of the vector f(t,y), and
!          norm() = the weighted vector norm used throughout, given by
!                   the DMNORM function routine, and weighted by the
!                   tolerances initially loaded into the EWT array.
!  The sign of H0 is inferred from the initial values of TOUT and T.
!  ABS(H0) is made .le. ABS(TOUT-T) in any case.
! -----------------------------------------------------------------------
                     if ( h0==0.0D0 ) then
                        tdist = abs(Tout-T)
                        w0 = max(abs(T),abs(Tout))
                        if ( tdist<2.0D0*dls1%uround*w0 ) then
                           msg = 'DLSODAR- TOUT(=R1) too close to T(=R2) to start integration.'
                           call xerrwd(msg,60,22,0,0,0,0,2,Tout,T)
                           goto 1200
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
                           sum = dmnorm(dls1%n,Rwork(lf0),Rwork(dls1%lewt))
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
                        Rwork(i+lf0-1) = h0*Rwork(i+lf0-1)
                     enddo
!
!  Check for a zero of g at T. ------------------------------------------
                     dlsr%irfnd = 0
                     dlsr%toutc = Tout
                     if ( dlsr%ngc==0 ) goto 200
                     call drchek(1,g,Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dlsr%lg0),Rwork(dlsr%lg1),Rwork(dlsr%lgx),Jroot,irt)
                     if ( irt==0 ) goto 200
                     msg = 'DLSODAR-  One or more components of g has a root  '
                     call xerrwd(msg,50,32,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      too near to the initial point.    '
                     call xerrwd(msg,40,32,0,0,0,0,0,0.0D0,0.0D0)
                     goto 1200
                  else
!  if ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
                     dls1%jstart = -1
                     if ( dls1%n/=dls1%nyh ) then
!  NEQ was reduced.  zero part of yh to avoid undefined references. -----
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
!  First, DRCHEK is called to check for a root within the dlsr%last step
!  taken, other than the dlsr%last root found there, if any.
!  If ITASK = 2 or 5, and y(TN) has not yet been returned to the user
!  because of an intervening root, return through Block G.
! -----------------------------------------------------------------------
 50   continue
   dls1%nslast = dls1%nst
!
   irfp = dlsr%irfnd
   if ( dlsr%ngc/=0 ) then
      if ( Itask==1 .or. Itask==4 ) dlsr%toutc = Tout
      call drchek(2,g,Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dlsr%lg0),Rwork(dlsr%lg1),Rwork(dlsr%lgx),Jroot,irt)
      if ( irt==1 ) then
         dlsr%irfnd = 1
         Istate = 3
         T = dlsr%t0
         goto 400
      endif
   endif
   dlsr%irfnd = 0
   if ( irfp==1 .and. dlsr%tlast/=dls1%tn .and. Itask==2 ) goto 300
!
   select case (Itask)
   case (2)
      goto 100
   case (3)
      tp = dls1%tn - dls1%hu*(1.0D0+100.0D0*dls1%uround)
      if ( (tp-Tout)*dls1%h>0.0D0 ) then
         msg = 'DLSODAR-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2) '
         call xerrwd(msg,60,23,0,1,Itask,0,2,Tout,tp)
         goto 1200
      else
         if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
         T = dls1%tn
         goto 300
      endif
   case (4)
      tcrit = Rwork(1)
      if ( (dls1%tn-tcrit)*dls1%h>0.0D0 ) goto 900
      if ( (tcrit-Tout)*dls1%h<0.0D0 ) goto 1000
      if ( (dls1%tn-Tout)*dls1%h>=0.0D0 ) then
         call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
         if ( iflag/=0 ) goto 1100
         T = Tout
         Istate = 2
         goto 400
      endif
   case (5)
      tcrit = Rwork(1)
      if ( (dls1%tn-tcrit)*dls1%h>0.0D0 ) goto 900
   case default
      if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
      call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
      if ( iflag/=0 ) goto 1100
      T = Tout
      Istate = 2
      goto 400
   endselect
   hmx = abs(dls1%tn) + abs(dls1%h)
   ihit = abs(dls1%tn-tcrit)<=100.0D0*dls1%uround*hmx
   if ( ihit ) T = tcrit
   if ( irfp==1 .and. dlsr%tlast/=dls1%tn .and. Itask==5 ) goto 300
   if ( ihit ) goto 300
   tnext = dls1%tn + dls1%h*(1.0D0+4.0D0*dls1%uround)
   if ( (tnext-tcrit)*dls1%h>0.0D0 ) then
      dls1%h = (tcrit-dls1%tn)*(1.0D0-4.0D0*dls1%uround)
      if ( Istate==2 .and. dls1%jstart>=0 ) dls1%jstart = -2
   endif
endif
! -----------------------------------------------------------------------
!  Block E.
!  The next block is normally executed for all calls and contains
!  the call to the one-step core integrator DSTODA.
!
!  This is a looping point for the integration steps.
!
!  First check for too many steps being taken, update EWT (if not at
!  start of problem), check for too much accuracy being requested, and
!  check for H below the roundoff level in T.
! -----------------------------------------------------------------------
 100  continue
if ( dls1%meth/=dlsa%mused ) then
   if ( dlsa%insufr==1 ) goto 500
   if ( dlsa%insufi==1 ) goto 600
endif
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
   msg = 'DLSODAR-  At current T (=R1), MXSTEP (=I1) steps  '
   call xerrwd(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
   msg = '      taken on this call before reaching TOUT     '
   call xerrwd(msg,50,201,0,1,dls1%mxstep,0,1,dls1%tn,0.0D0)
   Istate = -1
   goto 800
else
   call dewset(dls1%n,Itol,Rtol,Atol,Rwork(dls1%lyh),Rwork(dls1%lewt))
   do i = 1, dls1%n
      if ( Rwork(i+dls1%lewt-1)<=0.0D0 ) then
!  EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
         ewti = Rwork(dls1%lewt+i-1)
         msg = 'DLSODAR-  At T(=R1), EWT(I1) has become R2 .le. 0.'
         call xerrwd(msg,50,202,0,1,i,0,2,dls1%tn,ewti)
         Istate = -6
         goto 800
      else
         Rwork(i+dls1%lewt-1) = 1.0D0/Rwork(i+dls1%lewt-1)
      endif
   enddo
endif
 200  continue
tolsf = dls1%uround*dmnorm(dls1%n,Rwork(dls1%lyh),Rwork(dls1%lewt))
if ( tolsf<=1.0D0 ) then
   if ( (dls1%tn+dls1%h)==dls1%tn ) then
      dls1%nhnil = dls1%nhnil + 1
      if ( dls1%nhnil<=dls1%mxhnil ) then
         msg = 'DLSODAR-  Warning..Internal T(=R1) and H(=R2) are '
         call xerrwd(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
         msg = '      such that in the machine, T + H = T on the next step  '
         call xerrwd(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
         msg = '     (H = step size). Solver will continue anyway.'
         call xerrwd(msg,50,101,0,0,0,0,2,dls1%tn,dls1%h)
         if ( dls1%nhnil>=dls1%mxhnil ) then
            msg = 'DLSODAR-  Above warning has been issued I1 times. '
            call xerrwd(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
            msg = '     It will not be issued again for this problem.'
            call xerrwd(msg,50,102,0,1,dls1%mxhnil,0,0,0.0D0,0.0D0)
         endif
      endif
   endif
! -----------------------------------------------------------------------
!    CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,f,JAC,DPRJA,DSOLSY)
! -----------------------------------------------------------------------
   call dstoda(Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dls1%lyh),Rwork(dls1%lewt), &
    & Rwork(dls1%lsavf),Rwork(dls1%lacor),Rwork(dls1%lwm), &
             & Iwork(dls1%liwm),f,jac,dprja,dsolsy)
   kgo = 1 - dls1%kflag
   select case (kgo)
   case (2)
!  KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
      msg = 'DLSODAR-  At T(=R1), step size H(=R2), the error  '
      call xerrwd(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      test failed repeatedly or with ABS(H) = HMIN'
      call xerrwd(msg,50,204,0,0,0,0,2,dls1%tn,dls1%h)
      Istate = -4
      goto 700
   case (3)
!  KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
      msg = 'DLSODAR-  At T (=R1) and step size H (=R2), the   '
      call xerrwd(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      corrector convergence failed repeatedly     '
      call xerrwd(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      or with ABS(H) = HMIN   '
      call xerrwd(msg,30,205,0,0,0,0,2,dls1%tn,dls1%h)
      Istate = -5
      goto 700
   case default
! -----------------------------------------------------------------------
!  Block F.
!  The following block handles the case of a successful return from the
!  core integrator (KFLAG = 0).
!  If a method switch was just made, record TSW, reset MAXORD,
!  set JSTART to -1 to signal DSTODA to complete the switch,
!  and do extra printing of data if IXPR = 1.
!  Then call DRCHEK to check for a root within the dlsr%last step.
!  Then, if no root was found, check for stop conditions.
! -----------------------------------------------------------------------
      dls1%init = 1
      if ( dls1%meth/=dlsa%mused ) then
         dlsa%tsw = dls1%tn
         dls1%maxord = dlsa%mxordn
         if ( dls1%meth==2 ) dls1%maxord = dlsa%mxords
         if ( dls1%meth==2 ) Rwork(dls1%lwm) = sqrt(dls1%uround)
         dlsa%insufr = min(dlsa%insufr,1)
         dlsa%insufi = min(dlsa%insufi,1)
         dls1%jstart = -1
         if ( dlsa%ixpr/=0 ) then
            if ( dls1%meth==2 ) then
               msg = 'DLSODAR- A switch to the BDF (stiff) method has occurred    '
               call xerrwd(msg,60,105,0,0,0,0,0,0.0D0,0.0D0)
            endif
            if ( dls1%meth==1 ) then
               msg = 'DLSODAR- A switch to the Adams (nonstiff) method occurred   '
               call xerrwd(msg,60,106,0,0,0,0,0,0.0D0,0.0D0)
            endif
            msg = '     at T = R1,  tentative step size H = R2,  step NST = I1 '
            call xerrwd(msg,60,107,0,1,dls1%nst,0,2,dls1%tn,dls1%h)
         endif
      endif
!
      if ( dlsr%ngc/=0 ) then
         call drchek(3,g,Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dlsr%lg0),Rwork(dlsr%lg1),Rwork(dlsr%lgx),Jroot,irt)
         if ( irt==1 ) then
            dlsr%irfnd = 1
            Istate = 3
            T = dlsr%t0
            goto 400
         endif
      endif
!
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
                  if ( dls1%jstart>=0 ) dls1%jstart = -2
               endif
               goto 100
            endif
         else
            call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
            T = Tout
            Istate = 2
            goto 400
         endif
      case (5)
!  ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
         hmx = abs(dls1%tn) + abs(dls1%h)
         ihit = abs(dls1%tn-tcrit)<=100.0D0*dls1%uround*hmx
      case default
!  ITASK = 1.  If TOUT has been reached, interpolate. -------------------
         if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
         call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
         T = Tout
         Istate = 2
         goto 400
      endselect
   endselect
else
   tolsf = tolsf*2.0D0
   if ( dls1%nst==0 ) then
      msg = 'DLSODAR-  At start of problem, too much accuracy  '
      call xerrwd(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      requested for precision of machine..  See TOLSF (=R1) '
      call xerrwd(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
      Rwork(14) = tolsf
      goto 1200
   else
!  Too much accuracy requested for machine precision. -------------------
      msg = 'DLSODAR-  At T (=R1), too much accuracy requested '
      call xerrwd(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      for precision of machine..  See TOLSF (=R2) '
      call xerrwd(msg,50,203,0,0,0,0,2,dls1%tn,tolsf)
      Rwork(14) = tolsf
      Istate = -2
      goto 800
   endif
endif
! -----------------------------------------------------------------------
!  Block G.
!  The following block handles all successful returns from DLSODAR.
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
Istate = 2
 400  continue
Rwork(11) = dls1%hu
Rwork(12) = dls1%h
Rwork(13) = dls1%tn
Rwork(15) = dlsa%tsw
Iwork(11) = dls1%nst
Iwork(12) = dls1%nfe
Iwork(13) = dls1%nje
Iwork(14) = dls1%nqu
Iwork(15) = dls1%nq
Iwork(19) = dlsa%mused
Iwork(20) = dls1%meth
Iwork(10) = dlsr%nge
dlsr%tlast = T
return
!  RWORK length too small to proceed. -----------------------------------
 500  continue
msg = 'DLSODAR- At current T(=R1), RWORK length too small'
call xerrwd(msg,50,206,0,0,0,0,0,0.0D0,0.0D0)
msg = '      to proceed.  The integration was otherwise successful.'
call xerrwd(msg,60,206,0,0,0,0,1,dls1%tn,0.0D0)
Istate = -7
goto 800
!  IWORK length too small to proceed. -----------------------------------
 600  continue
msg = 'DLSODAR- At current T(=R1), IWORK length too small'
call xerrwd(msg,50,207,0,0,0,0,0,0.0D0,0.0D0)
msg = '      to proceed.  The integration was otherwise successful.'
call xerrwd(msg,60,207,0,0,0,0,1,dls1%tn,0.0D0)
Istate = -7
goto 800
!  Compute IMXER if relevant. -------------------------------------------
 700  continue
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
!  Set Y vector, T, and optional outputs. -------------------------------
 800  continue
do i = 1, dls1%n
   Y(i) = Rwork(i+dls1%lyh-1)
enddo
T = dls1%tn
Rwork(11) = dls1%hu
Rwork(12) = dls1%h
Rwork(13) = dls1%tn
Rwork(15) = dlsa%tsw
Iwork(11) = dls1%nst
Iwork(12) = dls1%nfe
Iwork(13) = dls1%nje
Iwork(14) = dls1%nqu
Iwork(15) = dls1%nq
Iwork(19) = dlsa%mused
Iwork(20) = dls1%meth
Iwork(10) = dlsr%nge
dlsr%tlast = T
return
 900  continue
msg = 'DLSODAR-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  '
call xerrwd(msg,60,24,0,0,0,0,2,tcrit,dls1%tn)
goto 1200
 1000 continue
msg = 'DLSODAR-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  '
call xerrwd(msg,60,25,0,0,0,0,2,tcrit,Tout)
goto 1200
 1100 continue
msg = 'DLSODAR-  Trouble in DINTDY. ITASK = I1, TOUT = R1'
call xerrwd(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
!
 1200 continue
Istate = -3
return
99999 continue
end subroutine dlsodar
