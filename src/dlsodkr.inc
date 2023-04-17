!==================================================================================================================================!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!==================================================================================================================================!
!>
!!### Synopsis
!!
!! DLSODKR: Livermore Solver for Ordinary Differential equations,
!!          with preconditioned Krylov iteration methods for the
!!          Newton correction linear systems, and with Rootfinding.
!!
!! DLSODKR solves the initial value problem for stiff or nonstiff
!! systems of first order ODEs,
!!```text
!!     dy/dt = f(t,y),  or, in component form,
!!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
!!```
!! At the same time, it locates the roots of any of a set of functions
!!```text
!!     g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).
!!```
!!
!!   This version is in double precision.
!!
!!-----------------------------------------------------------------------
!!### Introduction.
!!
!! This is a modification of the DLSODE package, and differs from it
!! in five ways:

!! 1.  It uses various preconditioned Krylov subspace iteration methods
!!     for the linear algebraic systems that arise in the case of stiff
!!     systems.  See the introductory notes below.
!! 2.  It does automatic switching between functional (fixpoint)
!!     iteration and Newton iteration in the corrector iteration.
!! 3.  It finds the root of at least one of a set of constraint
!!     functions g(i) of the independent and dependent variables.
!!     It finds only those roots for which some g(i), as a function
!!     of t, changes sign in the interval of integration.
!!     It then returns the solution at the root, if that occurs
!!     sooner than the specified stop condition, and otherwise returns
!!     the solution according the specified stop condition.
!! 4.  It supplies to JAC an input flag, JOK, which indicates whether
!!     JAC may (optionally) bypass the evaluation of Jacobian matrix data
!!     and instead process saved data (with the current value of scalar hl0).
!! 5.  It contains a new subroutine that calculates the initial step
!!     size to be attempted.
!!
!!### Introduction to the Krylov methods in DLSODKR:
!!
!! The linear systems that must be solved have the form
!!```text
!!   A * x  = b,  where  A = identity - hl0 * (df/dy) .
!!```
!! Here hl0 is a scalar, and df/dy is the Jacobian matrix of partial
!! derivatives of f (NEQ by NEQ).
!!
!! The particular Krylov method is chosen by setting the second digit,
!! MITER, in the method flag MF.
!! Currently, the values of MITER have the following meanings:
!!
!!  MITER  | description
!!  ------ | ------------------------------------------------------------
!!   1     | means the Scaled Preconditioned Incomplete
!!         | Orthogonalization Method (SPIOM).
!!         |
!!   2     | means an incomplete version of the preconditioned scaled
!!         | Generalized Minimal Residual method (SPIGMR).
!!         | This is the best choice in general.
!!         |
!!   3     | means the Preconditioned Conjugate Gradient method (PCG).
!!         | Recommended only when df/dy is symmetric or nearly so.
!!         |
!!   4     | means the scaled Preconditioned Conjugate Gradient method
!!         | (PCGS).  Recommended only when D-inverse \* df/dy \* D is
!!         | symmetric or nearly so, where D is the diagonal scaling
!!         | matrix with elements 1/EWT(i) (see RTOL/ATOL description).
!!         |
!!   9     | means that only a user-supplied matrix P (approximating A)
!!         | will be used, with no Krylov iteration done.  This option
!!         | allows the user to provide the complete linear system
!!         | solution algorithm, if desired.
!!
!! The user can apply preconditioning to the linear system A\*x = b,
!! by means of arbitrary matrices (the preconditioners).
!!
!! In the case of SPIOM and SPIGMR, one can apply left and right
!! preconditioners P1 and P2, and the basic iterative method is then
!! applied to the matrix (P1-inverse)\*A\*(P2-inverse) instead of to the
!! matrix A.  The product P1\*P2 should be an approximation to matrix A
!! such that linear systems with P1 or P2 are easier to solve than with
!! A.  Preconditioning from the left only or right only means using
!! P2 = identity or P1 = identity, respectively.
!!
!! In the case of the PCG and PCGS methods, there is only one
!! preconditioner matrix P (but it can be the product of more than one).
!! It should approximate the matrix A but allow for relatively
!! easy solution of linear systems with coefficient matrix P.
!! For PCG, P should be positive definite symmetric, or nearly so,
!! and for PCGS, the scaled preconditioner D-inverse \* P \* D
!! should be symmetric or nearly so.
!!
!! If the Jacobian J = df/dy splits in a natural way into a sum
!! J = J1 + J2, then one possible choice of preconditioners is
!! P1 = identity - hl0 \* J1  and  P2 = identity - hl0 \* J2
!! provided each of these is easy to solve (or approximately solve).
!!
!!-----------------------------------------------------------------------
!! Summary of Usage.
!!
!! Communication between the user and the DLSODKR package, for normal
!! situations, is summarized here.  This summary describes only a subset
!! of the full set of options available.  See the full description for
!! details, including optional communication, nonstandard options,
!! and instructions for special situations.  See also the demonstration
!! program distributed with this solver.
!!
!! **A.** First provide a subroutine of the form:
!!```fortran
!!               SUBROUTINE F (NEQ, T, Y, YDOT)
!!               DOUBLE PRECISION T, Y(*), YDOT(*)
!!```
!! which supplies the vector function f by loading YDOT(i) with f(i).
!!
!! **B.** Provide a subroutine of the form:
!!```fortran
!!               SUBROUTINE G (NEQ, T, Y, NG, GOUT)
!!               DOUBLE PRECISION T, Y(*), GOUT(NG)
!!```
!! which supplies the vector function g by loading GOUT(i) with
!! g(i), the i-th constraint function whose root is sought.
!!
!! **C.** Next determine (or guess) whether or not the problem is stiff.
!! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
!! whose real part is negative and large in magnitude, compared to the
!! reciprocal of the t span of interest.  If the problem is nonstiff,
!! use a method flag MF = 10.  If it is stiff, MF should be between 21
!! and 24, or possibly 29.  MF = 22 is generally the best choice.
!! Use 23 or 24 only if symmetry is present.  Use MF = 29 if the
!! complete linear system solution is to be provided by the user.
!! The following four parameters must also be set.
!!```text
!!  IWORK(1) = LWP  = length of real array WP for preconditioning.
!!  IWORK(2) = LIWP = length of integer array IWP for preconditioning.
!!  IWORK(3) = JPRE = preconditioner type flag:
!!                  = 0 for no preconditioning (P1 = P2 = P = identity)
!!                  = 1 for left-only preconditioning (P2 = identity)
!!                  = 2 for right-only preconditioning (P1 = identity)
!!                  = 3 for two-sided preconditioning (and PCG or PCGS)
!!  IWORK(4) = JACFLG = flag for whether JAC is called.
!!                    = 0 if JAC is not to be called,
!!                    = 1 if JAC is to be called.
!!```
!! Use JACFLG = 1 if JAC computes any nonconstant data for use in
!! preconditioning, such as Jacobian elements.
!! The arrays WP and IWP are work arrays under the user's control,
!! for use in the routines that perform preconditioning operations.
!!
!! **D.** If the problem is stiff, you must supply two routines that deal
!! with the preconditioning of the linear systems to be solved.
!!
!! These are as follows:
!!```fortran
!!     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY,V,HL0,JOK,WP,IWP,IER)
!!     DOUBLE PRECISION T, Y(*), YSV(*), REWT(*), FTY(*), V(*), HL0,WP(*)
!!     INTEGER IWP(*)
!!```
!! This routine must evaluate and preprocess any parts of the
!! Jacobian matrix df/dy involved in the preconditioners P1, P2, P.
!! The Y and FTY arrays contain the current values of y and f(t,y),
!! respectively, and YSV also contains the current value of y.
!! The array V is work space of length NEQ.
!!
!! JAC must multiply all computed Jacobian elements by the scalar
!! -HL0, add the identity matrix, and do any factorization
!! operations called for, in preparation for solving linear systems
!! with a coefficient matrix of P1, P2, or P.  The matrix P1\*P2 or P
!! should be an approximation to  identity - hl0 \* (df/dy).
!! JAC should return IER = 0 if successful, and IER .ne. 0 if not.
!! (If IER .ne. 0, a smaller time step will be tried.)
!!
!! JAC may alter Y and V, but not YSV, REWT, FTY, or HL0.
!!
!! The JOK argument can be ignored (or see full description below).
!!```fortran
!!     SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HL0, WP, IWP, B, LR, IER)
!!     DOUBLE PRECISION T, Y(*), FTY(*), WK(*), HL0, WP(*), B(*)
!!     INTEGER IWP(*)
!!```
!! This routine must solve a linear system with B as right-hand
!! side and one of the preconditioning matrices, P1, P2, or P, as
!! coefficient matrix, and return the solution vector in B.
!!
!! LR is a flag concerning left vs right preconditioning, input
!! to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
!!
!! In the case of the PCG or PCGS method, LR will be 3, and PSOL
!! should solve the system P\*x = B with the preconditioner matrix P.
!!
!! In the case MF = 29 (no Krylov iteration), LR will be 0,
!! and PSOL is to return in B the desired approximate solution
!! to A \* x = B, where A = identity - hl0 \* (df/dy).
!!
!! PSOL can use data generated in the JAC routine and stored in
!! WP and IWP.  WK is a work array of length NEQ.
!!
!! The argument HL0 is the current value of the scalar appearing
!! in the linear system.  If the old value, at the time of the last
!! JAC call, is needed, it must have been saved by JAC in WP.
!!
!! on return, PSOL should set the error flag IER as follows:
!!```text
!!       IER = 0 if PSOL was successful,
!!       IER .gt. 0 if a recoverable error occurred, meaning that the
!!              time step will be retried,
!!       IER .lt. 0 if an unrecoverable error occurred, meaning that the
!!              solver is to stop immediately.
!!```
!! **E.** Write a main program which calls Subroutine DLSODKR once for
!! each point at which answers are desired.  This should also provide
!! for possible use of logical unit 6 for output of error messages
!! by DLSODKR.  On the first call to DLSODKR, supply arguments as
!! follows:
!!
!! F
!!
!! : name of subroutine for right-hand side vector f.
!!   This name must be declared External in calling program.
!!
!! NEQ
!!
!! : number of first order ODEs.
!!
!! Y
!!
!! : array of initial values, of length NEQ.
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
!!   The estimated local error in y(i) will be controlled so as
!!   to be roughly less (in magnitude) than
!!             EWT(i) = RTOL\*ABS(Y(i)) + ATOL     if ITOL = 1, or
!!             EWT(i) = RTOL\*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!!   Thus the local error test passes if, in each component,
!!   either the absolute error is less than ATOL (or ATOL(i)),
!!   or the relative error is less than RTOL.
!!   Use RTOL = 0.0 for pure absolute error control, and
!!   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!!   control.  Caution: Actual (global) errors may exceed these
!!   local tolerances, so choose them conservatively.
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
!!        20 + 16*NEQ + 3*NG           for MF = 10,
!!        45 + 17*NEQ + 3*NG + LWP     for MF = 21,
!!        61 + 17*NEQ + 3*NG + LWP     for MF = 22,
!!        20 + 15*NEQ + 3*NG + LWP     for MF = 23 or 24,
!!        20 + 12*NEQ + 3*NG + LWP     for MF = 29.
!!```
!!
!! LRW
!!
!! : declared length of RWORK (in user's dimension).
!!
!! IWORK
!!
!! : integer work array of length at least:
!!```text
!!        30            for MF = 10,
!!        35 + LIWP     for MF = 21,
!!        30 + LIWP     for MF = 22, 23, 24, or 29.
!!```
!!
!! LIW
!!
!! : declared length of IWORK (in user's dimension).
!!
!! JAC,PSOL
!!
!! : names of subroutines for preconditioning.
!! These names must be declared External in the calling program.
!!
!! MF
!!
!! : method flag.  Standard values are:
!!
!!  value | description
!!  ----- | --------------------------------------------------------
!!    10  | for nonstiff (Adams) method.
!!    21  | for stiff (BDF) method, with preconditioned SIOM.
!!    22  | for stiff method, with preconditioned GMRES method.
!!    23  | for stiff method, with preconditioned CG method.
!!    24  | for stiff method, with scaled preconditioned CG method.
!!    29  | for stiff method, with user's PSOL routine only.
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
!!
!! Note that the main program must declare arrays Y, RWORK, IWORK,
!! JROOT, and possibly ATOL.
!!
!! **F.** The output from the first call (or any call) is:
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
!! : values and meanings:
!!
!!  value  | description
!!  ------ | -------------------------------------------------------
!!  2 or 3 | if DLSODKR was successful, negative otherwise.
!!       2 | means no root was found, and TOUT was reached as desired.
!!       3 | means a root was found prior to reaching TOUT.
!!      -1 | means excess work done on this call (perhaps wrong MF).
!!      -2 | means excess accuracy requested (tolerances too small).
!!      -3 | means illegal input detected (see printed message).
!!      -4 | means repeated error test failures (check all inputs).
!!      -5 | means repeated convergence failures (perhaps bad JAC
!!         | or PSOL routine supplied or wrong choice of MF or
!!         | tolerances, or this solver is inappropriate).
!!      -6 | means error weight became zero during problem. (Solution
!!         | component i vanished, and ATOL or ATOL(i) = 0.)
!!      -7 | means an unrecoverable error occurred in PSOL.
!!
!! JROOT
!!
!! : array showing roots found if ISTATE = 3 on return.
!!   JROOT(i) = 1 if g(i) has a root at T, or 0 otherwise.
!!
!! **G.** To continue the integration after a successful return, proceed
!! as follows:
!!
!! 1.  If ISTATE = 2 on return, reset TOUT and call DLSODKR again.
!! 2.  If ISTATE = 3 on return, reset ISTATE to 2 and call DLSODKR again.
!!
!! In either case, no other parameters need be reset.
!!
!!-----------------------------------------------------------------------
!!#### Full Description of User Interface to DLSODKR.
!!
!! The user interface to DLSODKR consists of the following parts.
!!
!! 1.   The call sequence to Subroutine DLSODKR, which is a driver
!!      routine for the solver.  This includes descriptions of both
!!      the call sequence arguments and of user-supplied routines.
!!      Following these descriptions is a description of
!!      optional inputs available through the call sequence, and then
!!      a description of optional outputs (in the work arrays).
!!
!! 2.   Descriptions of other routines in the DLSODKR package that may be
!!      (optionally) called by the user.  These provide the ability to
!!      alter error message handling, save and restore the internal
!!      Common, and obtain specified derivatives of the solution y(t).
!!
!! 3.   Descriptions of Common blocks to be declared in overlay
!!      or similar environments, or to be saved when doing an interrupt
!!      of the problem and continued solution later.
!!
!! 4.   Description of two routines in the DLSODKR package, either of
!!      which the user may replace with his/her own version, if desired.
!!      These relate to the measurement of errors.
!!
!!-----------------------------------------------------------------------
!!#### Part 1.  Call Sequence.
!!
!! The call sequence parameters used for input only are
!!  F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, PSOL, MF,
!!  G, and NG,
!!
!! that used only for output is  JROOT,
!! and those used for both input and output are
!!  Y, T, ISTATE.
!!
!! The work arrays RWORK and IWORK are also used for conditional and
!! optional inputs and optional outputs.  (The term output here refers
!! to the return from Subroutine DLSODKR to the user's calling program.)
!!
!! The legality of input parameters will be thoroughly checked on the
!! initial call for the problem, but not checked thereafter unless a
!! change in input parameters is flagged by ISTATE = 3 on input.
!!
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
!!           SUBROUTINE F (NEQ, T, Y, YDOT)
!!           DOUBLE PRECISION T, Y(*), YDOT(*)
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
!! externally to DLSODKR, an extra call to F should be made
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
!! these are to be accessed in the user-supplied routines.
!!
!! Normally, NEQ is a scalar, and it is generally referred to
!! as a scalar in this user interface description.  However,
!! NEQ may be an array, with NEQ(1) set to the system size.
!! (The DLSODKR package accesses only NEQ(1).)  In either case,
!! this parameter is passed as the NEQ argument in all calls
!! to the user-supplied routines.  Hence, if it is an array,
!! locations NEQ(2),... may be used to store other integer data
!! and pass it to the user-supplied routines. Each such routine
!! must include NEQ in a Dimension statement in that case.
!!
!! Y
!!
!! : a real array for the vector of dependent variables, of
!! length NEQ or more.  Used for both input and output on the
!! first call (ISTATE = 1), and only for output on other calls.
!! On the first call, Y must contain the vector of initial
!! values.  On output, Y contains the computed solution vector,
!! evaluated at T.  If desired, the Y array may be used
!! for other purposes between calls to the solver.
!!
!! This array is passed as the Y argument in all calls to F, G,
!! JAC, and PSOL.  Hence its length may exceed NEQ, and
!! locations Y(NEQ+1),... may be used to store other real data
!! and pass it to the user-supplied routines.
!! (The DLSODKR package accesses only Y(1),...,Y(NEQ).)
!!
!! T
!!
!! : the independent variable.  On input, T is used only on the
!! first call, as the initial point of the integration.
!! On output, after each call, T is the value at which a
!! computed solution y is evaluated (usually the same as TOUT).
!! If a root was found, T is the computed location of the
!! root reached first, on output.
!! On an error return, T is the farthest point reached.
!!
!! TOUT
!!
!! : the next value of t at which a computed solution is desired.
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
!!```text
!!                    RMS-norm of ( E(i)/EWT(i) )   .le.   1,
!!        where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!!```
!! and the RMS-norm (root-mean-square norm) here is
!! RMS-norm(v) = SQRT(sum v(i)\*\*2 / NEQ).  Here EWT = (EWT(i))
!! is a vector of weights which must always be positive, and
!! the values of RTOL and ATOL should all be non-negative.
!! The following table gives the types (scalar/array) of
!! RTOL and ATOL, and the corresponding form of EWT(i).
!!
!!  ITOL  | RTOL    | ATOL     | EWT(i)
!!  ----- | ------- | -------- | -----------------------------
!!    1   | scalar  |  scalar  |  RTOL\*ABS(Y(i)) + ATOL
!!    2   | scalar  |  array   |  RTOL\*ABS(Y(i)) + ATOL(i)
!!    3   | array   |  scalar  |  RTOL(i)\*ABS(Y(i)) + ATOL
!!    4   | array   |  array   |  RTOL(i)\*ABS(Y(i)) + ATOL(i)
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
!!  ----- | -------------------------------------------------------
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
!! Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!! (within roundoff), it will return T = TCRIT (exactly) to
!! indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!! in which case answers at T = TOUT are returned first).
!!
!! ISTATE
!!
!! : an index used for input and output to specify the
!! the state of the calculation.
!!
!! On input, the values of ISTATE are as follows.
!!
!!  value | description
!!  ----- | -------------------------------------------------------
!!      1 | means this is the first call for the problem
!!        | (initializations will be done).  See note below.
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
!!        | NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
!!        | and any of the optional inputs except H0.
!!        | In addition, immediately following a return with
!!        | ISTATE = 3 (root found), NG and G may be changed.
!!        | (But changing NG from 0 to .gt. 0 is not allowed.)
!!        |
!!        |  Note:  A preliminary call with TOUT = T is not counted
!!        |  as a first call here, as no initialization or checking of
!!        |  input is done.  (Such a call is sometimes useful for the
!!        |  purpose of outputting the initial conditions.)
!!        |  Thus the first call for which TOUT .ne. T requires
!!        |  ISTATE = 1 on input.
!!
!! On output, ISTATE has the following values and meanings.
!!
!!  value | description
!!  ----- | -------------------------------------------------------
!!      1 | means nothing was done; TOUT = T and ISTATE = 1 on input.
!!      2 | means the integration was performed successfully.
!!      3 | means the integration was successful, and one or more
!!        | roots were found before satisfying the stop condition
!!        | specified by ITASK.  See JROOT.
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
!!     -6 | means EWT(i) became zero for some i during the
!!        | integration.  Pure relative error control (ATOL(i)=0.0)
!!        | was requested on a variable which has now vanished.
!!        | The integration was successful as far as T.
!!     -7 | means the PSOL routine returned an unrecoverable error
!!        | flag (IER .lt. 0).  The integration was successful as
!!        | far as T.
!!
!!  **Note**:  Since the normal output value of ISTATE is 2,
!!  it does not need to be reset for normal continuation.
!!  Also, since a negative input value of ISTATE will be
!!  regarded as illegal, a negative output value requires the
!!  user to change it, and possibly other inputs, before
!!  calling the solver again.
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
!! : a real working array (double precision).
!!```text
!!      The length of RWORK must be at least
!!         20 + NYH*(MAXORD+1) + 3*NEQ + 3*NG + LENLS + LWP    where
!!      NYH    = the initial value of NEQ,
!!      MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!!               smaller value is given as an optional input),
!!      LENLS = length of work space for linear system (Krylov)
!!              method, excluding preconditioning:
!!        LENLS = 0                               if MITER = 0,
!!        LENLS = NEQ*(MAXL+3) + MAXL**2          if MITER = 1,
!!        LENLS = NEQ*(MAXL+3+MIN(1,MAXL-KMP))
!!             + (MAXL+3)*MAXL + 1                if MITER = 2,
!!        LENLS = 6*NEQ                           if MITER = 3 or 4,
!!        LENLS = 3*NEQ                           if MITER = 9.
!!      (See the MF description for METH and MITER, and the
!!      list of optional inputs for MAXL and KMP.)
!!      LWP = length of real user work space for preconditioning
!!      (see JAC/PSOL).
!!      Thus if default values are used and NEQ is constant,
!!      this length is:
!!         20 + 16*NEQ + 3*NG           for MF = 10,
!!         45 + 24*NEQ + 3*NG + LWP     for MF = 11,
!!         61 + 24*NEQ + 3*NG + LWP     for MF = 12,
!!         20 + 22*NEQ + 3*NG + LWP     for MF = 13 or 14,
!!         20 + 19*NEQ + 3*NG + LWP     for MF = 19,
!!         20 + 9*NEQ + 3*NG            for MF = 20,
!!         45 + 17*NEQ + 3*NG + LWP     for MF = 21,
!!         61 + 17*NEQ + 3*NG + LWP     for MF = 22,
!!         20 + 15*NEQ + 3*NG + LWP     for MF = 23 or 24,
!!         20 + 12*NEQ + 3*NG + LWP     for MF = 29.
!!```
!! The first 20 words of RWORK are reserved for conditional
!! and optional inputs and optional outputs.
!!
!! The following word in RWORK is a conditional input:
!!
!! RWORK(1) = TCRIT = critical value of t which the solver
!! is not to overshoot.  Required if ITASK is
!! 4 or 5, and ignored otherwise.  (See ITASK.)
!!
!! LRW
!!
!! : the length of the array RWORK, as declared by the user.
!! (This will be checked by the solver.)
!!
!! IWORK
!!
!! : an integer work array.  The length of IWORK must be at least
!!             30                 if MITER = 0 (MF = 10 or 20),
!!             30 + MAXL + LIWP   if MITER = 1 (MF = 11, 21),
!!             30 + LIWP          if MITER = 2, 3, 4, or 9.
!!          MAXL = 5 unless a different optional input value is given.
!!          LIWP = length of integer user work space for preconditioning
!!          (see conditional input list following).
!!          The first few words of IWORK are used for conditional and
!!          optional inputs and optional outputs.
!!
!!          The following 4 words in IWORK are conditional inputs,
!!          required if MITER .ge. 1:
!!          IWORK(1) = LWP  = length of real array WP for use in
!!                     preconditioning (part of RWORK array).
!!          IWORK(2) = LIWP = length of integer array IWP for use in
!!                     preconditioning (part of IWORK array).
!!                     The arrays WP and IWP are work arrays under the
!!                     user's control, for use in the routines that
!!                     perform preconditioning operations (JAC and PSOL).
!!          IWORK(3) = JPRE = preconditioner type flag:
!!                   = 0 for no preconditioning (P1 = P2 = P = identity)
!!                   = 1 for left-only preconditioning (P2 = identity)
!!                   = 2 for right-only preconditioning (P1 = identity)
!!                   = 3 for two-sided preconditioning (and PCG or PCGS)
!!          IWORK(4) = JACFLG = flag for whether JAC is called.
!!                   = 0 if JAC is not to be called,
!!                   = 1 if JAC is to be called.
!!                     Use JACFLG = 1 if JAC computes any nonconstant
!!                     data needed in preconditioning operations,
!!                     such as some of the Jacobian elements.
!!
!! LIW
!!
!! : the length of the array IWORK, as declared by the user.
!! (This will be checked by the solver.)
!!
!! Note:  The work arrays must not be altered between calls to DLSODKR
!! for the same problem, except possibly for the conditional and
!! optional inputs, and except for the last 3\*NEQ words of RWORK.
!! The latter space is used for internal scratch space, and so is
!! available for use by the user outside DLSODKR between calls, if
!! desired (but not for use by any of the user-supplied routines).
!!
!! JAC
!!
!! : the name of the user-supplied routine to compute any
!! Jacobian elements (or approximations) involved in the
!! matrix preconditioning operations (MITER .ge. 1).
!! It is to have the form
!!```fortran
!!        SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, &
!!       &                HL0, JOK, WP, IWP, IER)
!!        DOUBLE PRECISION T, Y(*), YSV(*), REWT(*), FTY(*), V(*), &
!!       &                 HL0, WP(*)
!!        INTEGER IWP(*)
!!```
!! This routine must evaluate and preprocess any parts of the
!! Jacobian matrix df/dy used in the preconditioners P1, P2, P.
!!
!! The Y and FTY arrays contain the current values of y and
!! f(t,y), respectively, and the YSV array also contains
!! the current y vector.  The array V is work space of length
!! NEQ for use by JAC.  REWT is the array of reciprocal error
!! weights (1/EWT).  JAC must multiply all computed Jacobian
!! elements by the scalar -HL0, add the identity matrix, and do
!! any factorization operations called for, in preparation
!! for solving linear systems with a coefficient matrix of
!! P1, P2, or P.  The matrix P1\*P2 or P should be an &
!! approximation to  identity - hl0 \* (df/dy).  JAC should
!! return IER = 0 if successful, and IER .ne. 0 if not.
!!
!! (If IER .ne. 0, a smaller time step will be tried.)
!! The arrays WP (of length LWP) and IWP (of length LIWP)
!! are for use by JAC and PSOL for work space and for storage
!! of data needed for the solution of the preconditioner
!! linear systems.  Their lengths and contents are under the
!! user's control.
!!
!! The argument JOK is an input flag for optional use
!! by JAC in deciding whether to recompute Jacobian elements
!! or use saved values.  If JOK = -1, then JAC must compute
!! any relevant Jacobian elements (or approximations) used in
!! the preconditioners.  Optionally, JAC may also save these
!! elements for later reuse.  If JOK = 1, the integrator has
!! made a judgement (based on the convergence history and the
!! value of HL0) that JAC need not recompute Jacobian elements,
!! but instead use saved values, and the current value of HL0,
!! to reconstruct the preconditioner matrices, followed by
!! any required factorizations.  This may be cost-effective if
!! Jacobian elements are costly and storage is available.
!!
!! JAC may alter Y and V, but not YSV, REWT, FTY, or HL0.
!! JAC must be declared External in the calling program.
!!
!! Subroutine JAC may access user-defined quantities in
!! NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!! (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
!! See the descriptions of NEQ and Y above.
!!
!! PSOL
!!
!! : the name of the user-supplied routine for the
!! solution of preconditioner linear systems.
!! It is to have the form
!!```fortran
!!         SUBROUTINE PSOL (NEQ, T, Y, FTY, WK,HL0, WP,IWP, B, LR,IER)
!!         DOUBLE PRECISION T, Y(*), FTY(*), WK(*), HL0, WP(*), B(*)
!!         INTEGER IWP(*)
!!```
!! This routine must solve a linear system with B as right-hand
!! side and one of the preconditioning matrices, P1, P2, or P,
!! as coefficient matrix, and return the solution vector in B.
!! LR is a flag concerning left vs right preconditioning, input
!! to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
!! In the case of the PCG or PCGS method, LR will be 3, and PSOL
!! should solve the system P\*x = B with the preconditioner P.
!! In the case MITER = 9 (no Krylov iteration), LR will be 0,
!! and PSOL is to return in B the desired approximate solution
!! to A \* x = B, where A = identity - hl0 \* (df/dy).
!! PSOL can use data generated in the JAC routine and stored in
!! WP and IWP.
!!
!! The Y and FTY arrays contain the current values of y and
!! f(t,y), respectively.  The array WK is work space of length
!! NEQ for use by PSOL.
!!
!! The argument HL0 is the current value of the scalar appearing
!! in the linear system.  If the old value, as of the last
!! JAC call, is needed, it must have been saved by JAC in WP.
!!
!! On return, PSOL should set the error flag IER as follows:
!!```text
!!        IER = 0 if PSOL was successful,
!!        IER .gt. 0 on a recoverable error, meaning that the
!!               time step will be retried,
!!        IER .lt. 0 on an unrecoverable error, meaning that the
!!               solver is to stop immediately.
!!      PSOL may not alter Y, FTY, or HL0.
!!      PSOL must be declared External in the calling program.
!!           Subroutine PSOL may access user-defined quantities in
!!      NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
!!      (dimensioned in PSOL) and/or Y has length exceeding NEQ(1).
!!      See the descriptions of NEQ and Y above.
!!```
!! MF
!!
!! : the method flag.  Used only for input.  The legal values of
!! MF are 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, and 29.
!! MF has decimal digits METH and MITER: MF = 10\*METH + MITER.
!! METH indicates the basic linear multistep method:
!!
!!  METH | description
!!  ---- | ----------------------------------------------------
!!     1 | means the implicit Adams method.
!!     2 | means the method based on Backward
!!
!! Differentiation Formulas (BDFs).
!!
!! MITER indicates the corrector iteration method:

!!  MITER | description
!!  ----- | ----------------------------------------------------
!!  0     | means functional iteration (no linear system
!!        | is involved).
!!  1     | means Newton iteration with Scaled Preconditioned
!!        | Incomplete Orthogonalization Method (SPIOM)
!!        | for the linear systems.
!!  2     | means Newton iteration with Scaled Preconditioned
!!        | Incomplete Generalized Minimal Residual method
!!        | (SPIGMR) for the linear systems.
!!  3     | means Newton iteration with Preconditioned
!!        | Conjugate Gradient method (PCG)
!!        | for the linear systems.
!!  4     | means Newton iteration with scaled preconditioned
!!        | Conjugate Gradient method (PCGS)
!!        | for the linear systems.
!!  9     | means Newton iteration with only the
!!        | user-supplied PSOL routine called (no Krylov
!!        | iteration) for the linear systems.
!!        | JPRE is ignored, and PSOL is called with LR = 0.
!!
!! See comments in the introduction about the choice of MITER.
!! If MITER .ge. 1, the user must supply routines JAC and PSOL
!! (the names are arbitrary) as described above.
!! For MITER = 0, a dummy argument can be used.
!!
!! G
!!
!! : the name of subroutine for constraint functions, whose
!! roots are desired during the integration.  It is to have
!! the form
!!```fortran
!!          SUBROUTINE G (NEQ, T, Y, NG, GOUT)
!!          DOUBLE PRECISION T, Y(*), GOUT(NG)
!!```
!! where NEQ, T, Y, and NG are input, and the array GOUT
!! is output.  NEQ, T, and Y have the same meaning as in
!! the F routine, and GOUT is an array of length NG.
!! For i = 1,...,NG, this routine is to load into GOUT(i)
!! the value at (t,y) of the i-th constraint function g(i).
!! DLSODKR will find roots of the g(i) of odd multiplicity
!! (i.e. sign changes) as they occur during the integration.
!! G must be declared External in the calling program.
!!
!! Caution: Because of numerical errors in the functions
!! g(i) due to roundoff and integration error, DLSODKR may
!! return false roots, or return the same root at two or more
!! nearly equal values of t.  If such false roots are
!! suspected, the user should consider smaller error tolerances
!! and/or higher precision in the evaluation of the g(i).
!!
!! If a root of some g(i) defines the end of the problem,
!! the input to DLSODKR should nevertheless allow integration
!! to a point slightly past that root, so that DLSODKR can
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
!! JROOT(i) = 1 if g(i) has a root at t, or JROOT(i) = 0 if not.
!!-----------------------------------------------------------------------
!! Optional Inputs.
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
!!       |          |
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
!! DELT  | RWORK(8) | convergence test constant in Krylov iteration
!!       |          | algorithm.  The default is .05.
!!       |          |
!! MAXORD| IWORK(5) | the maximum order to be allowed.  The default
!!       |          | value is 12 if METH = 1, and 5 if METH = 2.
!!       |          | If MAXORD exceeds the default value, it will
!!       |          | be reduced to the default value.
!!       |          | If MAXORD is changed during the problem, it may
!!       |          | cause the current order to be reduced.
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
!! MAXL  | IWORK(8) | maximum number of iterations in the SPIOM, SPIGMR,
!!       |          | PCG, or PCGS algorithm (.le. NEQ).
!!       |          | The default is MAXL = MIN(5,NEQ).
!!       |          |
!! KMP   | IWORK(9) | number of vectors on which orthogonalization
!!       |          | is done in SPIOM or SPIGMR algorithm (.le. MAXL).
!!       |          | The default is KMP = MAXL.
!!       |          |
!!       |          | Note:  When KMP .lt. MAXL and MF = 22, the length
!!       |          | of RWORK must be defined accordingly.  See
!!       |          | the definition of RWORK above.
!!       |          |
!!
!!-----------------------------------------------------------------------
!! Optional Outputs.
!!
!! As optional additional output from DLSODKR, the variables listed
!! below are quantities related to the performance of DLSODKR
!! which are available to the user.  These are communicated by way of
!! the work arrays, but also have internal mnemonic names as shown.
!!
!! Except where stated otherwise, all of these outputs are defined
!! on any successful return from DLSODKR, and on any return with
!! ISTATE = -1, -2, -4, -5, -6, or -7.  On an illegal input return
!! (ISTATE = -3), they will be unchanged from their existing values
!! (if any), except possibly for TOLSF, LENRW, and LENIW.
!!
!! On any error return, outputs relevant to the error will be defined,
!! as noted below.
!!
!!  Name  |  Location |     Meaning
!!  ----- | --------- | ---------------------------------------------------
!!  HU    |  RWORK(11)| the step size in t last used (successfully).
!!        |           |
!!  HCUR  |  RWORK(12)| the step size to be attempted on the next step.
!!        |           |
!!  TCUR  |  RWORK(13)| the current value of the independent variable
!!        |           | which the solver has actually reached, i.e. the
!!        |           | current internal mesh point in t.  On output, TCUR
!!        |           | will always be at least as far as the argument
!!        |           | T, but may be farther (if interpolation was done).
!!        |           |
!!  TOLSF |  RWORK(14)| a tolerance scale factor, greater than 1.0,
!!        |           | computed when a request for too much accuracy was
!!        |           | detected (ISTATE = -3 if detected at the start of
!!        |           | the problem, ISTATE = -2 otherwise).  If ITOL is
!!        |           | left unaltered but RTOL and ATOL are uniformly
!!        |           | scaled up by a factor of TOLSF for the next call,
!!        |           | then the solver is deemed likely to succeed.
!!        |           | (The user may also ignore TOLSF and alter the
!!        |           | tolerance parameters in any other way appropriate.)
!!        |           |
!!  NGE   |  IWORK(10)| the number of g evaluations for the problem so far.
!!        |           |
!!  NST   |  IWORK(11)| the number of steps taken for the problem so far.
!!        |           |
!!  NFE   |  IWORK(12)| the number of f evaluations for the problem so far.
!!        |           |
!!  NPE   |  IWORK(13)| the number of calls to JAC so far (for evaluation
!!        |           | of preconditioners).
!!        |           |
!!  NQU   |  IWORK(14)| the method order last used (successfully).
!!        |           |
!!  NQCUR |  IWORK(15)| the order to be attempted on the next step.
!!        |           |
!!  IMXER |  IWORK(16)| the index of the component of largest magnitude in
!!        |           | the weighted local error vector ( E(i)/EWT(i) ),
!!        |           | on an error return with ISTATE = -4 or -5.
!!        |           |
!!  LENRW |  IWORK(17)| the length of RWORK actually required.
!!        |           | This is defined on normal returns and on an illegal
!!        |           | input return for insufficient storage.
!!        |           |
!!  LENIW |  IWORK(18)| the length of IWORK actually required.
!!        |           | This is defined on normal returns and on an illegal
!!        |           | input return for insufficient storage.
!!        |           |
!!  NNI   |  IWORK(19)| number of nonlinear iterations so far (each of
!!        |           | which calls an iterative linear solver).
!!        |           |
!!  NLI   |  IWORK(20)| number of linear iterations so far.
!!        |           | Note: A measure of the success of algorithm is
!!        |           | the average number of linear iterations per
!!        |           | nonlinear iteration, given by NLI/NNI.
!!        |           | If this is close to MAXL, MAXL may be too small.
!!        |           |
!!  NPS   |  IWORK(21)| number of preconditioning solve operations
!!        |           | (PSOL calls) so far.
!!        |           |
!!  NCFN  |  IWORK(22)| number of convergence failures of the nonlinear
!!        |           | (Newton) iteration so far.
!!        |           | Note: A measure of success is the overall
!!        |           | rate of nonlinear convergence failures, NCFN/NST.
!!        |           |
!!  NCFL  |  IWORK(23)| number of convergence failures of the linear
!!        |           | iteration so far.
!!        |           | Note: A measure of success is the overall
!!        |           | rate of linear convergence failures, NCFL/NNI.
!!        |           |
!!  NSFI  |  IWORK(24)| number of functional iteration steps so far.
!!        |           | Note: A measure of the extent to which the
!!        |           | problem is nonstiff is the ratio NSFI/NST.
!!        |           |
!!  NJEV  |  IWORK(25)| number of JAC calls with JOK = -1 so far
!!        |           | (number of evaluations of Jacobian data).
!!        |           |
!!
!! The following two arrays are segments of the RWORK array which
!! may also be of interest to the user as optional outputs.
!! For each array, the table below gives its internal name,
!! its base address in RWORK, and its description.
!!
!!  Name  | Base Address  |   Description
!!  ----- | ------------- | -----------------------------------------------
!!  YH    | 21 + 3\*NG    | the Nordsieck history array, of size NYH by
!!        |               | (NQCUR + 1), where NYH is the initial value
!!        |               | of NEQ.  For j = 0,1,...,NQCUR, column j+1
!!        |               | of YH contains HCUR\*\*j/factorial(j) times
!!        |               | the j-th derivative of the interpolating
!!        |               | polynomial currently representing the solution,
!!        |               | evaluated at t = TCUR.
!!        |               |
!!  ACOR  |  LENRW-NEQ+1  | array of size NEQ used for the accumulated
!!        |               | corrections on each step, scaled on output
!!        |               | to represent the estimated local error in y
!!        |               | on the last step.  This is the vector E in
!!        |               | the description of the error control.  It is
!!        |               | defined only on a successful return from
!!        |               | DLSODKR.
!!
!!-----------------------------------------------------------------------
!!#### Part 2.  Other Routines Callable.
!!
!! The following are optional calls which the user may make to
!! gain additional capabilities in conjunction with DLSODKR.
!! (The routines XSETUN and XSETF are designed to conform to the
!! SLATEC error handling package.)
!!
!!  Form of Call              |      Function
!!  ------------------------- | ------------------------------------------
!!  CALL XSETUN(LUN)          | Set the logical unit number, LUN, for
!!                            | output of messages from DLSODKR, if
!!                            | the default is not desired.
!!                            | The default value of LUN is 6.
!!                            |
!!  CALL XSETF(MFLAG)         | Set a flag to control the printing of
!!                            | messages by DLSODKR.
!!                            | MFLAG = 0 means do not print. (Danger:
!!                            | This risks losing valuable information.)
!!                            | MFLAG = 1 means print (the default).
!!                            |
!!                            | Either of the above calls may be made at
!!                            | any time and will take effect immediately.
!!                            |
!!  CALL DSRCKR(RSAV,ISAV,JOB)| saves and restores the contents of
!!                            | the internal Common blocks used by
!!                            | DLSODKR (see Part 3 below).
!!                            | RSAV must be a real array of length 228
!!                            | or more, and ISAV must be an integer
!!                            | array of length 63 or more.
!!                            |     JOB=1 means save Common into RSAV/ISAV.
!!                            |     JOB=2 means restore Common from RSAV/ISAV.
!!                            |
!!                            | DSRCKR is useful if one is
!!                            | interrupting a run and restarting
!!                            | later, or alternating between two or
!!                            | more problems solved with DLSODKR.
!!                            |
!!  CALL DINTDY(,,,,,)        | Provide derivatives of y, of various
!!       (see below)          | orders, at a specified point t, if
!!                            | desired.  It may be called only after
!!                            | a successful return from DLSODKR.
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
!! (normally the same as the T last returned by DLSODKR).
!! For valid results, T must lie between TCUR - HU and TCUR.
!! (See optional outputs for TCUR and HU.)
!! K
!!
!! : integer order of the derivative desired.  K must satisfy
!! 0 .le. K .le. NQCUR, where NQCUR is the current order
!! (see optional outputs).  The capability corresponding
!! to K = 0, i.e. computing y(T), is already provided
!! by DLSODKR directly.  Since NQCUR .ge. 1, the first
!! derivative dy/dt is always available with DINTDY.
!! LYH
!!
!! : 21 + 3\*NG = base address in RWORK of the history array YH.
!! NYH
!!
!! : column length of YH, equal to the initial value of NEQ.
!!
!! The output parameters are:
!!
!! DKY
!!
!! : a real array of length NEQ containing the computed value
!!             of the K-th derivative of y(t).
!! IFLAG
!!
!! : integer flag, returned as 0 if K and T were legal,
!! -1 if K was illegal, and -2 if T was illegal.
!! On an error return, a message is also written.
!!-----------------------------------------------------------------------
!!#### Part 3.  Common Blocks.
!!
!! If DLSODKR is to be used in an overlay situation, the user
!! must declare, in the primary overlay, the variables in:
!!```text
!!   (1) the call sequence to DLSODKR, and
!!   (2) the four internal Common blocks
!!         /DLS001/  of length  255  (218 double precision words
!!                      followed by 37 integer words),
!!         /DLS002/  of length   5  (1 double precision word
!!                      followed by  4 integer words),
!!         /DLPK01/  of length  17  (4 double precision words
!!                      followed by 13 integer words),
!!         /DLSR01/  of length  14     (5 double precision words
!!                      followed by  9 integer words).
!!```
!!
!! If DLSODKR is used on a system in which the contents of internal
!! Common blocks are not preserved between calls, the user should
!! declare the above Common blocks in the calling program to insure
!! that their contents are preserved.
!!
!! If the solution of a given problem by DLSODKR is to be interrupted
!! and then later continued, such as when restarting an interrupted run
!! or alternating between two or more problems, the user should save,
!! following the return from the last DLSODKR call prior to the
!! interruption, the contents of the call sequence variables and the
!! internal Common blocks, and later restore these values before the
!! next DLSODKR call for that problem.  To save and restore the Common
!! blocks, use Subroutine DSRCKR (see Part 2 above).
!!
!!-----------------------------------------------------------------------
!!#### Part 4.  Optionally Replaceable Solver Routines.
!!
!! Below are descriptions of two routines in the DLSODKR package which
!! relate to the measurement of errors.  Either routine can be
!! replaced by a user-supplied version, if desired.  However, since such
!! a replacement may have a major impact on performance, it should be
!! done only when absolutely necessary, and only with great caution.
!! (Note: The means by which the package version of a routine is
!! superseded by the user's version may be system-dependent.)
!!
!! (a) DEWSET()
!!
!! The following subroutine is called just before each internal
!! integration step, and sets the array of error weights, EWT, as
!! described under ITOL/RTOL/ATOL above:
!!```fortran
!!     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
!!```
!! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODKR call sequence,
!! YCUR contains the current dependent variable vector, and
!! EWT is the array of weights set by DEWSET.
!!
!! If the user supplies this subroutine, it must return in EWT(i)
!! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
!! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
!! routine (see below), and also used by DLSODKR in the computation
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
!!     DOUBLE PRECISION RLS
!!     COMMON /DLS001/ RLS(218),ILS(37)
!!     NQ = ILS(33)
!!     NST = ILS(34)
!!     H = RLS(212)
!!```
!! Thus, for example, the current value of dy/dt can be obtained as
!! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
!! unnecessary when NST = 0).
!!
!! (b) DVNORM()
!! The following is a real function routine which computes the weighted
!! root-mean-square norm of a vector v:
!!```text
!!     D = DVNORM (N, V, W)
!! where:
!!   N = the length of the vector,
!!   V = real array of length N containing the vector,
!!   W = real array of length N containing weights,
!!   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
!!```
!! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
!! EWT is as set by Subroutine DEWSET.
!!
!! If the user supplies this function, it should return a non-negative
!! value of DVNORM suitable for use in the error control in DLSODKR.
!! None of the arguments should be altered by DVNORM.
!! For example, a user-supplied DVNORM routine might:
!!
!!   - substitute a max-norm of (V(i)\*W(i)) for the RMS-norm, or
!!   - ignore some components of V in the norm, with the effect of
!!     suppressing the error control on those components of y.
!!
!!-----------------------------------------------------------------------
!!### Pedigree:
!!
!!   This is the 18 November 2003 version of
!!   DLSODKR is derived from the Livermore Solver for Ordinary
!!   Differential Equations package ODEPACK,
!!
!!-----------------------------------------------------------------------
!! References:
!! 1.  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage Matrix
!!     Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31 (1989),
!!     pp. 40-91; also  L.L.N.L. Report UCRL-95088, Rev. 1, June 1987.
!! 2.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!!     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
!!     North-Holland, Amsterdam, 1983, pp. 55-64.
!!-----------------------------------------------------------------------
!! Authors:       Alan C. Hindmarsh and Peter N. Brown
!!                Center for Applied Scientific Computing, L-561
!!                Lawrence Livermore National Laboratory
!!                Livermore, CA 94551
!!-----------------------------------------------------------------------
!
! ### REVISION HISTORY  (YYYYMMDD)
!  19900117  DATE WRITTEN
!  19900503  Added iteration switching (functional/Newton).
!  19900802  Added flag for Jacobian-saving in user preconditioner.
!  19900910  Added new initial stepsize routine LHIN.
!  19901019  Corrected LHIN - y array restored.
!  19910909  Changed names STOPK to STOKA, PKSET to SETPK;
!            removed unused variables in driver declarations;
!            minor corrections to main prologue.
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
!  20020502  Corrected declarations in descriptions of user routines.
!  20030603  Corrected duplicate type declaration for DUMACH.
!  20031105  Restored 'own' variables to Common blocks, to enable
!            interrupt/restart feature.
!  20031112  Added SAVE statements for data-loaded constants.
!  20031117  Changed internal name NPE to NJE.
!
! -----------------------------------------------------------------------
!  Other routines in the DLSODKR package.
!
!  In addition to Subroutine DLSODKR, the DLSODKR package includes the
!  following subroutines and function routines:
!   DLHIN    calculates a step size to be attempted initially.
!   DRCHEK   does preliminary checking for roots, and serves as an
!            interface between Subroutine DLSODKR and Subroutine DROOTS.
!   DROOTS   finds the leftmost root of a set of functions.
!   DINTDY   computes an interpolated value of the y vector at t = TOUT.
!   DEWSET   sets the error weight vector EWT before each step.
!   DVNORM   computes the weighted RMS-norm of a vector.
!   DSTOKA   is the core integrator, which does one step of the
!            integration and the associated error control.
!   DCFODE   sets all method coefficients and test constants.
!   DSETPK   interfaces between DSTOKA and the JAC routine.
!   DSOLPK   manages solution of linear system in Newton iteration.
!   DSPIOM   performs the SPIOM algorithm.
!   DATV     computes a scaled, preconditioned product (I-hl0*J)*v.
!   DORTHOG  orthogonalizes a vector against previous basis vectors.
!   DHEFA    generates an LU factorization of a Hessenberg matrix.
!   DHESL    solves a Hessenberg square linear system.
!   DSPIGMR  performs the SPIGMR algorithm.
!   DHEQR    generates a QR factorization of a Hessenberg matrix.
!   DHELS    finds the least squares solution of a Hessenberg system.
!   DPCG     performs preconditioned conjugate gradient algorithm (PCG).
!   DPCGS    performs the PCGS algorithm.
!   DATP     computes the product A*p, where A = I - hl0*df/dy.
!   DUSOL    interfaces to the user's PSOL routine (MITER = 9).
!   DSRCKR   is a user-callable routine to save and restore
!            the contents of the internal Common blocks.
!   DAXPY, DCOPY, DDOT, DNRM2, and DSCAL   are basic linear
!            algebra modules (from the BLAS collection).
!   XERRWD, XSETUN, XSETF, and IXSAV handle the printing of all
!            error messages and warnings.  XERRWD is machine-dependent.
!  Note:  DVNORM, DDOT, DNRM2, and IXSAV are function
!  routines.  All the others are subroutines.
!
!------------------------------------------------------------------------
subroutine dlsodkr(f,Neq,Y,T,Tout,Itol,Rtol,Atol,Itask,Istate,Iopt,Rwork,Lrw,Iwork,Liw,jac,psol,Mf,g,Ng,Jroot)

external f
external g
external jac
external psol
real(kind=dp), dimension(*) :: Atol, Rtol, Y
real(kind=dp) :: atoli, avdim, big, ewti, h0, hmax, hmx, rcfl, rcfn, rh, rtoli, size, tcrit, tnext, tolsf, tp
integer :: i, i1, i2, ier, iflag, imxer, irfp, irt, kgo, leniw, leniwk, lenrw, lenwk, lenwm, lenyh, lf0, liwp,     &
         & lwp, lyhnew, ncfl0, ncfn0, niter, nli0, nni0, nnid, nstd, nwarn
logical :: ihit, lavd, lcfl, lcfn, lwarn
integer :: Iopt, Istate, Itask, Itol, Liw, Lrw, Mf, Ng
integer, intent(inout), dimension(Liw) :: Iwork
integer, dimension(*) :: Neq
integer :: Jroot(*)
integer, dimension(2), save :: mord
character(60) :: msg
integer, save :: mxhnl0, mxstp0
real(kind=dp), intent(inout), dimension(Lrw) :: Rwork
real(kind=dp), intent(inout) :: T, Tout
!
! -----------------------------------------------------------------------
!  The following four internal Common blocks contain
!  (a) variables which are local to any subroutine but whose values must
!      be preserved between calls to the routine ("own" variables), and
!  (b) variables which are communicated between subroutines.
!  The block DLS001 is declared in subroutines DLSODKR, DINTDY,
!  DSTOKA, DSOLPK, and DATV.
!  The block DLS002 is declared in subroutines DLSODKR and DSTOKA.
!  The block DLSR01 is declared in subroutines DLSODKR, DRCHEK, DROOTS.
!  The block DLPK01 is declared in subroutines DLSODKR, DSTOKA, DSETPK,
!  and DSOLPK.
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
if(ng.ne.0) jroot(:ng) = 0

if ( Istate<1 .or. Istate>3 ) then
! -----------------------------------------------------------------------
!  Block I.
!  The following block handles all error returns due to illegal input
!  (ISTATE = -3), as detected before calling the core integrator.
!  First the error message routine is called.  If the illegal input
!  is a negative ISTATE, the run is aborted (apparent infinite loop).
! -----------------------------------------------------------------------
   msg = 'DLSODKR-  ISTATE(=I1) illegal.'
   call xerrwd(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
   if ( Istate>=0 ) goto 1000
!
   msg = 'DLSODKR-  Run aborted.. apparent infinite loop.   '
   call xerrwd(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
   goto 99999
else
   if ( Itask<1 .or. Itask>5 ) then
      msg = 'DLSODKR-  ITASK (=I1) illegal.'
      call xerrwd(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
      goto 1000
   else
      dlsr%itaskc = Itask
      if ( Istate==1 ) then
         dls1%init = 0
         if ( Tout==T ) return
      elseif ( dls1%init==0 ) then
         msg = 'DLSODKR- ISTATE.gt.1 but DLSODKR not initialized. '
         call xerrwd(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
         goto 1000
      elseif ( Istate==2 ) then
         goto 50
      endif
! -----------------------------------------------------------------------
!  Block B.
!  The next code block is executed for the initial call (ISTATE = 1),
!  or for a continuation call with parameter changes (ISTATE = 3).
!  It contains checking of all inputs and various initializations.
!
!  First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF,
!  and NG.
! -----------------------------------------------------------------------
      if ( Neq(1)<=0 ) then
         msg = 'DLSODKR-  NEQ (=I1) .lt. 1    '
         call xerrwd(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
         goto 1000
      else
         if ( Istate/=1 ) then
            if ( Neq(1)>dls1%n ) then
               msg = 'DLSODKR-  ISTATE = 3 and NEQ increased (I1 to I2).'
               call xerrwd(msg,50,5,0,2,dls1%n,Neq(1),0,0.0D0,0.0D0)
               goto 1000
            endif
         endif
         dls1%n = Neq(1)
         if ( Itol<1 .or. Itol>4 ) then
            msg = 'DLSODKR-  ITOL (=I1) illegal. '
            call xerrwd(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
            goto 1000
         elseif ( Iopt<0 .or. Iopt>1 ) then
            msg = 'DLSODKR-  IOPT (=I1) illegal. '
            call xerrwd(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
            goto 1000
         else
            dls1%meth = Mf/10
            dls1%miter = Mf - 10*dls1%meth
            if ( dls1%meth<1 .or. dls1%meth>2 ) goto 600
            if ( dls1%miter<0 ) goto 600
            if ( dls1%miter>4 .and. dls1%miter<9 ) goto 600
            if ( dls1%miter>=1 ) dlpk%jpre = Iwork(3)
            dlpk%jacflg = 0
            if ( dls1%miter>=1 ) dlpk%jacflg = Iwork(4)
            if ( Ng<0 ) then
               msg = 'DLSODKR-  NG (=I1) .lt. 0     '
               call xerrwd(msg,30,30,0,1,Ng,0,0,0.0D0,0.0D0)
               goto 1000
            else
               if ( Istate/=1 ) then
                  if ( dlsr%irfnd==0 .and. Ng/=dlsr%ngc ) then
                     msg = 'DLSODKR-  NG changed (from I1 to I2) illegally,   '
                     call xerrwd(msg,50,31,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      i.e. not immediately after a root was found.'
                     call xerrwd(msg,50,31,0,2,dlsr%ngc,Ng,0,0.0D0,0.0D0)
                     goto 1000
                  endif
               endif
               dlsr%ngc = Ng
!  Next process and check the optional inputs. --------------------------
               if ( Iopt==1 ) then
                  dls1%maxord = Iwork(5)
                  if ( dls1%maxord<0 ) then
                     msg = 'DLSODKR-  MAXORD (=I1) .lt. 0 '
                     call xerrwd(msg,30,11,0,1,dls1%maxord,0,0,0.0D0,0.0D0)
                     goto 1000
                  else
                     if ( dls1%maxord==0 ) dls1%maxord = 100
                     dls1%maxord = min(dls1%maxord,mord(dls1%meth))
                     dls1%mxstep = Iwork(6)
                     if ( dls1%mxstep<0 ) then
                        msg = 'DLSODKR-  MXSTEP (=I1) .lt. 0 '
                        call xerrwd(msg,30,12,0,1,dls1%mxstep,0,0,0.0D0,0.0D0)
                        goto 1000
                     else
                        if ( dls1%mxstep==0 ) dls1%mxstep = mxstp0
                        dls1%mxhnil = Iwork(7)
                        if ( dls1%mxhnil<0 ) then
                           msg = 'DLSODKR-  MXHNIL (=I1) .lt. 0 '
                           call xerrwd(msg,30,13,0,1,dls1%mxhnil,0,0,0.0D0,0.0D0)
                           goto 1000
                        else
                           if ( dls1%mxhnil==0 ) dls1%mxhnil = mxhnl0
                           if ( Istate==1 ) then
                              h0 = Rwork(5)
                              if ( (Tout-T)*h0<0.0D0 ) then
                                 msg = 'DLSODKR-  TOUT (=R1) behind T (=R2)     '
                                 call xerrwd(msg,40,14,0,0,0,0,2,Tout,T)
                                 msg = '      Integration direction is given by H0 (=R1)  '
                                 call xerrwd(msg,50,14,0,0,0,0,1,h0,0.0D0)
                                 goto 1000
                              endif
                           endif
                           hmax = Rwork(6)
                           if ( hmax<0.0D0 ) then
                              msg = 'DLSODKR-  HMAX (=R1) .lt. 0.0 '
                              call xerrwd(msg,30,15,0,0,0,0,1,hmax,0.0D0)
                              goto 1000
                           else
                              dls1%hmxi = 0.0D0
                              if ( hmax>0.0D0 ) dls1%hmxi = 1.0D0/hmax
                              dls1%hmin = Rwork(7)
                              if ( dls1%hmin<0.0D0 ) then
                                 msg = 'DLSODKR-  HMIN (=R1) .lt. 0.0 '
                                 call xerrwd(msg,30,16,0,0,0,0,1,dls1%hmin,0.0D0)
                                 goto 1000
                              else
                                 dlpk%maxl = Iwork(8)
                                 if ( dlpk%maxl==0 ) dlpk%maxl = 5
                                 dlpk%maxl = min(dlpk%maxl,dls1%n)
                                 dlpk%kmp = Iwork(9)
                                 if ( dlpk%kmp==0 .or. dlpk%kmp>dlpk%maxl ) dlpk%kmp = dlpk%maxl
                                 dlpk%delt = Rwork(8)
                                 if ( dlpk%delt==0.0D0 ) dlpk%delt = 0.05D0
                              endif
                           endif
                        endif
                     endif
                  endif
               else
                  dls1%maxord = mord(dls1%meth)
                  dls1%mxstep = mxstp0
                  dls1%mxhnil = mxhnl0
                  if ( Istate==1 ) h0 = 0.0D0
                  dls1%hmxi = 0.0D0
                  dls1%hmin = 0.0D0
                  dlpk%maxl = min(5,dls1%n)
                  dlpk%kmp = dlpk%maxl
                  dlpk%delt = 0.05D0
               endif
! -----------------------------------------------------------------------
!  Set work array pointers and check lengths LRW and LIW.
!  Pointers to segments of RWORK and IWORK are named by prefixing L to
!  the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
!  RWORK segments (in order) are denoted  G0, G1, GX, YH, WM,
!  EWT, SAVF, SAVX, ACOR.
! -----------------------------------------------------------------------
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
               dls1%lwm = dls1%lyh + (dls1%maxord+1)*dls1%nyh
               if ( dls1%miter==0 ) lenwk = 0
               if ( dls1%miter==1 ) lenwk = dls1%n*(dlpk%maxl+2) + dlpk%maxl*dlpk%maxl
               if ( dls1%miter==2 ) lenwk = dls1%n*(dlpk%maxl+2+min(1,dlpk%maxl-dlpk%kmp)) + (dlpk%maxl+3)*dlpk%maxl + 1
               if ( dls1%miter==3 .or. dls1%miter==4 ) lenwk = 5*dls1%n
               if ( dls1%miter==9 ) lenwk = 2*dls1%n
               lwp = 0
               if ( dls1%miter>=1 ) lwp = Iwork(1)
               lenwm = lenwk + lwp
               dlpk%locwp = lenwk + 1
               dls1%lewt = dls1%lwm + lenwm
               dls1%lsavf = dls1%lewt + dls1%n
               dlpk%lsavx = dls1%lsavf + dls1%n
               dls1%lacor = dlpk%lsavx + dls1%n
               if ( dls1%miter==0 ) dls1%lacor = dls1%lsavf + dls1%n
               lenrw = dls1%lacor + dls1%n - 1
               Iwork(17) = lenrw
               dls1%liwm = 31
               leniwk = 0
               if ( dls1%miter==1 ) leniwk = dlpk%maxl
               liwp = 0
               if ( dls1%miter>=1 ) liwp = Iwork(2)
               leniw = 30 + leniwk + liwp
               dlpk%lociwp = leniwk + 1
               Iwork(18) = leniw
               if ( lenrw>Lrw ) then
                  msg = 'DLSODKR-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) '
                  call xerrwd(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                  goto 1000
               elseif ( leniw>Liw ) then
                  msg = 'DLSODKR-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) '
                  call xerrwd(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                  goto 1000
               else
!  Check RTOL and ATOL for legality. ------------------------------------
                  rtoli = Rtol(1)
                  atoli = Atol(1)
                  do i = 1, dls1%n
                     if ( Itol>=3 ) rtoli = Rtol(i)
                     if ( Itol==2 .or. Itol==4 ) atoli = Atol(i)
                     if ( rtoli<0.0D0 ) then
                        msg = 'DLSODKR-  RTOL(I1) is R1 .lt. 0.0       '
                        call xerrwd(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
                        goto 1000
                     elseif ( atoli<0.0D0 ) then
                        msg = 'DLSODKR-  ATOL(I1) is R1 .lt. 0.0       '
                        call xerrwd(msg,40,20,0,1,i,0,1,atoli,0.0D0)
                        goto 1000
                     endif
                  enddo
!  Load SQRT(N) and its reciprocal in Common. ---------------------------
                  dlpk%sqrtn = sqrt(real(dls1%n))
                  dlpk%rsqrtn = 1.0D0/dlpk%sqrtn
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
                     if ( Itask==4 .or. Itask==5 ) then
                        tcrit = Rwork(1)
                        if ( (tcrit-Tout)*(Tout-T)<0.0D0 ) goto 800
                        if ( h0/=0.0D0 .and. (T+h0-tcrit)*h0>0.0D0 ) h0 = tcrit - T
                     endif
                     dls1%jstart = 0
                     dls1%nhnil = 0
                     dls1%nst = 0
                     dls1%nje = 0
                     dls1%nslast = 0
                     nli0 = 0
                     nni0 = 0
                     ncfn0 = 0
                     ncfl0 = 0
                     nwarn = 0
                     dls1%hu = 0.0D0
                     dls1%nqu = 0
                     dls1%ccmax = 0.3D0
                     dls1%maxcor = 3
                     dls1%msbp = 20
                     dls1%mxncf = 10
                     dlpk%nni = 0
                     dlpk%nli = 0
                     dlpk%nps = 0
                     dlpk%ncfn = 0
                     dlpk%ncfl = 0
                     dls%nsfi = 0
                     dls%njev = 0
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
                           msg = 'DLSODKR-  EWT(I1) is R1 .le. 0.0        '
                           call xerrwd(msg,40,21,0,1,i,0,1,ewti,0.0D0)
                           goto 1000
                        else
                           Rwork(i+dls1%lewt-1) = 1.0D0/Rwork(i+dls1%lewt-1)
                        endif
                     enddo
                     if ( h0==0.0D0 ) then
!  Call DLHIN to set initial step size H0 to be attempted. --------------
                        call dlhin(Neq,dls1%n,T,Rwork(dls1%lyh),Rwork(lf0), &
                        & f,Tout,dls1%uround,Rwork(dls1%lewt),Itol,Atol,Y,     &
                                 & Rwork(dls1%lacor),h0,niter,ier)
                        dls1%nfe = dls1%nfe + niter
                        if ( ier/=0 ) then
                           msg = 'DLSODKR- TOUT(=R1) too close to T(=R2) to start integration.'
                           call xerrwd(msg,60,22,0,0,0,0,2,Tout,T)
                           goto 1000
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
!  Check for a zero of g at T. ------------------------------------------
                     dlsr%irfnd = 0
                     dlsr%toutc = Tout
                     if ( dlsr%ngc==0 ) goto 200
                     call drchek(1,g,Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dlsr%lg0),Rwork(dlsr%lg1),Rwork(dlsr%lgx),Jroot,irt)
                     if ( irt==0 ) goto 200
                     msg = 'DLSODKR-  One or more components of g has a root  '
                     call xerrwd(msg,50,32,0,0,0,0,0,0.0D0,0.0D0)
                     msg = '      too near to the initial point.    '
                     call xerrwd(msg,40,32,0,0,0,0,0,0.0D0,0.0D0)
                     goto 1000
                  else
!  If ISTATE = 3, set flag to signal parameter changes to DSTOKA.--------
                     dls1%jstart = -1
                     if ( dls1%nq>dls1%maxord ) then
!  MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
                        do i = 1, dls1%n
                           Rwork(i+dls1%lsavf-1) = Rwork(i+dls1%lwm-1)
                        enddo
                     endif
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
   nli0 = dlpk%nli
   nni0 = dlpk%nni
   ncfn0 = dlpk%ncfn
   ncfl0 = dlpk%ncfl
   nwarn = 0
   select case (Itask)
   case (2)
      goto 100
   case (3)
      tp = dls1%tn - dls1%hu*(1.0D0+100.0D0*dls1%uround)
      if ( (tp-Tout)*dls1%h>0.0D0 ) then
         msg = 'DLSODKR-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2) '
         call xerrwd(msg,60,23,0,1,Itask,0,2,Tout,tp)
         goto 1000
      else
         if ( (dls1%tn-Tout)*dls1%h>=0.0D0 ) goto 300
         goto 100
      endif
   case (4)
      tcrit = Rwork(1)
      if ( (dls1%tn-tcrit)*dls1%h>0.0D0 ) goto 700
      if ( (tcrit-Tout)*dls1%h<0.0D0 ) goto 800
      if ( (dls1%tn-Tout)*dls1%h>=0.0D0 ) then
         call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
         if ( iflag/=0 ) goto 900
         T = Tout
         Istate = 2
         goto 400
      endif
   case (5)
      tcrit = Rwork(1)
      if ( (dls1%tn-tcrit)*dls1%h>0.0D0 ) goto 700
   case default
      if ( (dls1%tn-Tout)*dls1%h<0.0D0 ) goto 100
      call dintdy(Tout,0,Rwork(dls1%lyh),dls1%nyh,Y,iflag)
      if ( iflag/=0 ) goto 900
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
      if ( Istate==2 ) dls1%jstart = -2
   endif
endif
! -----------------------------------------------------------------------
!  Block E.
!  The next block is normally executed for all calls and contains
!  the call to the one-step core integrator DSTOKA.
!
!  This is a looping point for the integration steps.
!
!  First check for too many steps being taken,
!  check for poor Newton/Krylov method performance, update EWT (if not
!  at start of problem), check for too much accuracy being requested,
!  and check for H below the roundoff level in T.
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
   msg = 'DLSODKR-  At current T (=R1), MXSTEP (=I1) steps  '
   call xerrwd(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
   msg = '      taken on this call before reaching TOUT     '
   call xerrwd(msg,50,201,0,1,dls1%mxstep,0,1,dls1%tn,0.0D0)
   Istate = -1
   goto 500
else
   nstd = dls1%nst - dls1%nslast
   nnid = dlpk%nni - nni0
   if ( nstd>=10 .and. nnid/=0 ) then
      avdim = real(dlpk%nli-nli0)/real(nnid)
      rcfn = real(dlpk%ncfn-ncfn0)/real(nstd)
      rcfl = real(dlpk%ncfl-ncfl0)/real(nnid)
      lavd = avdim>(dlpk%maxl-0.05D0)
      lcfn = rcfn>0.9D0
      lcfl = rcfl>0.9D0
      lwarn = lavd .or. lcfn .or. lcfl
      if ( lwarn ) then
         nwarn = nwarn + 1
         if ( nwarn<=10 ) then
            if ( lavd ) then
               msg = 'DLSODKR- Warning. Poor iterative algorithm performance seen '
               call xerrwd(msg,60,111,0,0,0,0,0,0.0D0,0.0D0)
            endif
            if ( lavd ) then
               msg = '      at T = R1 by average no. of linear iterations = R2    '
               call xerrwd(msg,60,111,0,0,0,0,2,dls1%tn,avdim)
            endif
            if ( lcfn ) then
               msg = 'DLSODKR- Warning. Poor iterative algorithm performance seen '
               call xerrwd(msg,60,112,0,0,0,0,0,0.0D0,0.0D0)
            endif
            if ( lcfn ) then
               msg = '      at T = R1 by nonlinear convergence failure rate = R2  '
               call xerrwd(msg,60,112,0,0,0,0,2,dls1%tn,rcfn)
            endif
            if ( lcfl ) then
               msg = 'DLSODKR- Warning. Poor iterative algorithm performance seen '
               call xerrwd(msg,60,113,0,0,0,0,0,0.0D0,0.0D0)
            endif
            if ( lcfl ) then
               msg = '      at T = R1 by linear convergence failure rate = R2     '
               call xerrwd(msg,60,113,0,0,0,0,2,dls1%tn,rcfl)
            endif
         endif
      endif
   endif
   call dewset(dls1%n,Itol,Rtol,Atol,Rwork(dls1%lyh),Rwork(dls1%lewt))
   do i = 1, dls1%n
      if ( Rwork(i+dls1%lewt-1)<=0.0D0 ) then
!  EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
         ewti = Rwork(dls1%lewt+i-1)
         msg = 'DLSODKR-  At T(=R1), EWT(I1) has become R2 .le. 0.'
         call xerrwd(msg,50,202,0,1,i,0,2,dls1%tn,ewti)
         Istate = -6
         goto 500
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
         msg = 'DLSODKR-  Warning.. Internal T(=R1) and H(=R2) are'
         call xerrwd(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
         msg = '      such that in the machine, T + H = T on the next step  '
         call xerrwd(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
         msg = '     (H = step size). Solver will continue anyway.'
         call xerrwd(msg,50,101,0,0,0,0,2,dls1%tn,dls1%h)
         if ( dls1%nhnil>=dls1%mxhnil ) then
            msg = 'DLSODKR-  Above warning has been issued I1 times. '
            call xerrwd(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
            msg = '     It will not be issued again for this problem.'
            call xerrwd(msg,50,102,0,1,dls1%mxhnil,0,0,0.0D0,0.0D0)
         endif
      endif
   endif
! -----------------------------------------------------------------------
!      CALL DSTOKA(NEQ,Y,YH,NYH,YH,EWT,SAVF,SAVX,ACOR,WM,IWM,f,JAC,PSOL)
! -----------------------------------------------------------------------
   call dstoka(Neq,Y,Rwork(dls1%lyh),dls1%nyh,Rwork(dls1%lyh),Rwork(dls1%lewt), &
   & Rwork(dls1%lsavf),Rwork(dlpk%lsavx),Rwork(dls1%lacor)&
             & ,Rwork(dls1%lwm),Iwork(dls1%liwm),f,jac,psol)
   kgo = 1 - dls1%kflag
   select case (kgo)
   case (2)
!  KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
      msg = 'DLSODKR- At T(=R1) and step size H(=R2), the error'
      call xerrwd(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      test failed repeatedly or with ABS(H) = HMIN'
      call xerrwd(msg,50,204,0,0,0,0,2,dls1%tn,dls1%h)
      Istate = -4
!  Compute IMXER if relevant. -------------------------------------------
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
      goto 500
   case (3)
!  KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
      msg = 'DLSODKR-  At T (=R1) and step size H (=R2), the   '
      call xerrwd(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      corrector convergence failed repeatedly     '
      call xerrwd(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      or with ABS(H) = HMIN   '
      call xerrwd(msg,30,205,0,0,0,0,2,dls1%tn,dls1%h)
      Istate = -5
      goto 500
   case (4)
!  KFLAG = -3.  Unrecoverable error from PSOL. --------------------------
      msg = 'DLSODKR-  At T (=R1) an unrecoverable error return'
      call xerrwd(msg,50,206,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      was made from Subroutine PSOL     '
      call xerrwd(msg,40,206,0,0,0,0,1,dls1%tn,0.0D0)
      Istate = -7
      goto 500
   case default
! -----------------------------------------------------------------------
!  Block F.
!  The following block handles the case of a successful return from the
!  core integrator (KFLAG = 0).
!  Call DRCHEK to check for a root within the dlsr%last step.
!  Then, if no root was found, check for stop conditions.
! -----------------------------------------------------------------------
      dls1%init = 1
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
                  dls1%jstart = -2
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
      msg = 'DLSODKR-  At start of problem, too much accuracy  '
      call xerrwd(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      requested for precision of machine..  See TOLSF (=R1) '
      call xerrwd(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
      Rwork(14) = tolsf
      goto 1000
   else
!  Too much accuracy requested for machine precision. -------------------
      msg = 'DLSODKR-  At T (=R1), too much accuracy requested '
      call xerrwd(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
      msg = '      for precision of machine..  See TOLSF (=R2) '
      call xerrwd(msg,50,203,0,0,0,0,2,dls1%tn,tolsf)
      Rwork(14) = tolsf
      Istate = -2
      goto 500
   endif
endif
! -----------------------------------------------------------------------
!  Block G.
!  The following block handles all successful returns from DLSODKR.
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
Iwork(11) = dls1%nst
Iwork(12) = dls1%nfe
Iwork(13) = dls1%nje
Iwork(14) = dls1%nqu
Iwork(15) = dls1%nq
Iwork(19) = dlpk%nni
Iwork(20) = dlpk%nli
Iwork(21) = dlpk%nps
Iwork(22) = dlpk%ncfn
Iwork(23) = dlpk%ncfl
Iwork(24) = dls%nsfi
Iwork(25) = dls%njev
Iwork(10) = dlsr%nge
dlsr%tlast = T
return
!  Set Y vector, T, and optional outputs. -------------------------------
 500  continue
do i = 1, dls1%n
   Y(i) = Rwork(i+dls1%lyh-1)
enddo
T = dls1%tn
Rwork(11) = dls1%hu
Rwork(12) = dls1%h
Rwork(13) = dls1%tn
Iwork(11) = dls1%nst
Iwork(12) = dls1%nfe
Iwork(13) = dls1%nje
Iwork(14) = dls1%nqu
Iwork(15) = dls1%nq
Iwork(19) = dlpk%nni
Iwork(20) = dlpk%nli
Iwork(21) = dlpk%nps
Iwork(22) = dlpk%ncfn
Iwork(23) = dlpk%ncfl
Iwork(24) = dls%nsfi
Iwork(25) = dls%njev
Iwork(10) = dlsr%nge
dlsr%tlast = T
return
 600  continue
msg = 'DLSODKR-  MF (=I1) illegal.   '
call xerrwd(msg,30,8,0,1,Mf,0,0,0.0D0,0.0D0)
goto 1000
 700  continue
msg = 'DLSODKR-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  '
call xerrwd(msg,60,24,0,0,0,0,2,tcrit,dls1%tn)
goto 1000
 800  continue
msg = 'DLSODKR-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  '
call xerrwd(msg,60,25,0,0,0,0,2,tcrit,Tout)
goto 1000
 900  continue
msg = 'DLSODKR-  Trouble in DINTDY. ITASK = I1, TOUT = R1'
call xerrwd(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
!
 1000 continue
Istate = -3
return
99999 continue
end subroutine dlsodkr
