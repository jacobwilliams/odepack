<!--
see more recent:

   https://computing.llnl.gov/projects/odepack

might have test cases:

   https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html

markdown references

   https://michelf.ca/projects/php-markdown/extra
-->

![ODEPACK](docs/images/odepack.png)
============

# module ODEPACK

This is a WIP(Work In Progress) to evaluate the effort required and
feasibility of updating older Fortran code from the netlib repository
using a combination of the commercial plusFORT/spag software and GNU
utilities as well as conventional manual refactoring.

Many thanks to plusFORT for making an evaluation copy available for
several months to the Fortran Community. The plusFORT tools were crucial
to this study and unmatched in my experience for removing deprecated
syntax from pre-f2003 code.

This began with [ODEPACK](https://computing.llnl.gov/projects/odepack)

the ODEPACK package has been selected as a significant code that is well
documented and structured for its vintage, and available on the netlib
site and covered by a public domain license.

#### preliminary targets for phase I:

 + remove obsolescent syntax (trying plusFORT and spag)
 + able to build using fpm(1)  (the Fortran Package Manager)
   - in debug mode (ongoing)
   - with ifort, gfortran, nvfortran
 + text viewable in ford(1) and extractable as markup that can run through pandoc(1) (ongoing)
 + available on github (or equivalent)
 + no common blocks
 + no equivalences
 + build as a (single?) module M_odepack (ongoing)
 + complete set of unit tests (tests currently only use the original examples)

The originally bundled subset of BLAS/LAPACK routines are
being included in the module.  In production, this might not be done in
order to be able to easily call external optimized versions.

The biggest hinderance is some storage used for both INTEGER and DOUBLE
PRECISION values.

One take-away is how critical unit tests are to enable rapid development
(which so far this package does not have)

The initial pass was done just using the original sample programs as
unit tests. This may have allowed for introduction of errors as this
is a WIP but the original samples run with the same output as the original.

plusFORT was invaluable and reduced the effort by an estimated 85
percent. The results have been encouraging enough to inspire completing
the transformation.

### Phase II ?

Some of the goals of phase I are incomplete, but The results of this first
pass were significant enough that this project will hopefully continue.

A complete unit test suite is required to allow development to proceed
rapidily. Contributions, especially from current ODEPACK users are
particularly welcome.

Another major issue is the remaining non-standard code.  Non-standard (but
at the time de-facto fortran standard) such as equivalencing different
types, creating scratch space that is used as different numeric types,
and treating scalars as arrays and vice-versa as well as passing the same
arrays or values multiple times are the most time-consuming usages to
correct to standard-conforming, particularly since spag(1) had already
done an excellent job with updating the pre-f2003 code. spag(1) is not
(currently?) sufficient by itself to automate the additional refactoring
desired, which includes using post-f95 features and code restructuring,
so the remaining work requires manual recoding.

The type-mismatch issues have not been eliminated enough to include all
the routines in the module, so  those in the files "M_da1/dprep.inc"
"M_da1/dainvgs.inc" "M_da1/dprepi.inc" "M_da1/dstodi.inc" and
"M_da1_/dstode.inc" still require being built without an interface
definition.

### building with fpm(1)

This version of ODEPACK already builds with an included make(1) file
and as an fpm(1) package with the current options:
```bash
 fpm run                     --compiler nvfortran --example '*'
 fpm run --profile release   --compiler ifort     --example '*'
 # gfortran for production
 fpm run --profile release --flag -fallow-argument-mismatch --compiler gfortran  --example '*'
 # gfortran for debug
 fpm run --profile debug --flag -fallow-argument-mismatch --flag -std=f2018 --compiler gfortran  --example '*' --verbose
```
### building with make(1)
```bash

cd src
# gfortran
make clean
make gfortran
make run
make test

# ifort
make clean
make ifort
make run
make test

# nvfortran
make clean
make nvfortran
make run
make test
```
### To rebuild the documentation with ford(1), use

```bash
 ford ford.md
```

The code is far more readable after having been refactored by
a combination of using spag(1) from the plusFORT package and manual
editing, and is believed as useable as the original.

There are a few notes in src/M_odepack.f90 concerning continuing issues.

Current users of ODEPACK are encouraged to try this version and provide
feedback.

Hopefully as a community we can complete creating a new maintained
production-quality version of this venerable and still-valuable package.

### Documentation

The ongoing API documentation for the current `master` branch can
be found [here](https://jacobwilliams.github.io/odepack/).
This is generated by processing the source files with
[FORD](https://github.com/Fortran-FOSS-Programmers/ford).

#### Solvers for explicitly given systems.
Links to the solver-specific documenation for the main procedures (as described below):

   + [dlsoda](https://jacobwilliams.github.io/odepack/proc/dlsoda.html)
   + [dlsodar](https://jacobwilliams.github.io/odepack/proc/dlsodar.html)
   + [dlsode](https://jacobwilliams.github.io/odepack/proc/dlsode.html)
   + [dlsodes](https://jacobwilliams.github.io/odepack/proc/dlsodes.html)
   + [dlsodi](https://jacobwilliams.github.io/odepack/proc/dlsodi.html)
   + [dlsodis](https://jacobwilliams.github.io/odepack/proc/dlsodis.html)
   + [dlsodkr](https://jacobwilliams.github.io/odepack/proc/dlsodkr.html)
   + [dlsodpk](https://jacobwilliams.github.io/odepack/proc/dlsodpk.html)
   + [dlsoibt](https://jacobwilliams.github.io/odepack/proc/dlsoibt.html)

Here is an overview primarily from the original documentation  ...

--------------------------------------------------------------------------------

## Brief Description of ODEPACK - A Systematized Collection of ODE Solvers (Double Precision Version)

```text
Alan C. Hindmarsh
Center for Applied Scientific Computing, L-561
Lawrence Livermore National Laboratory
Livermore, CA 94551, U.S.A.

20 June 2001
```

Work performed under the auspices of the U.S. Department of Energy
by the Lawrence Livermore National Laboratory under contract
No. W-7405-Eng-48, and supported (formerly) by the DOE Office of
Energy Research, Applied Mathematical Sciences Research Program.

---------------------------------------------------------------------------

   ODEPACK is a collection of Fortran solvers for the initial value
problem for ordinary differential equation systems.  It consists of nine
solvers, namely a basic solver called LSODE and eight variants of it --
LSODES, LSODA, LSODAR, LSODPK, LSODKR, LSODI, LSOIBT, and LSODIS.
The collection is suitable for both stiff and nonstiff systems.  It
includes solvers for systems given in explicit form, dy/dt = f(t,y),
and also solvers for systems given in linearly implicit form,
A(t,y) dy/dt = g(t,y).  Two of the solvers use general sparse matrix
solvers for the linear systems that arise.  Two others use iterative
(preconditioned Krylov) methods instead of direct methods for these
linear systems.  The most recent addition is LSODIS, which solves
implicit problems with general sparse treatment of all matrices involved.

   The ODEPACK solvers are written in standard Fortran 77, with a few
exceptions, and with minimal machine dependencies.  There are separate
double and single precision versions of ODEPACK.  The actual solver
names are those given above with a prefix of D- or S- for the double
or single precision version, respectively, i.e. DLSODE/SLSODE, etc.
Each solver consists of a main driver subroutine having the same name
as the solver and some number of subordinate routines.  For each
solver, there is also a demonstration program, which solves one or two
simple problems in a somewhat self-checking manner.

   Recently, the ODEPACK solvers were upgraded to improve their
portability in numerous ways.  Among the improvements are (a) renaming
of routines and Common blocks to distinguish double and single
precision versions, (b) use of generic intrinsic function names, (c)
elimination of the Block Data subprogram, (d) use of a portable
routine to set the unit roundoff, and (e) passing of quoted strings to
the error message handler.  In addition, the prologue and internal
comments were reformatted, and use mixed upper/lower case.  Numerous
minor corrections and improvements were also made.

   The above upgrade operations were applied to LSODE earlier than they
were to the rest of ODEPACK, and the two upgrades were done somewhat
independently.  As a result, some differences will be apparent in the
source files of LSODE and the other solvers -- primarily in the
formatting of the comment line prologue of the main driver routine.
In Subroutines DLSODE/SLSODE and their subordinate routines, the
prologue was written in "SLATEC format", while for the other solvers a
more relaxed style was used.  The differences are entirely cosmetic,
however, and do not affect performance.

   Documentation on the usage of each solver is provided in the
initial block of comment lines in the source file, which (in most
cases) includes a simple example.  A demonstration program (in
separate double/single precision versions) is also available.

   What follows is a summary of the capabilities of ODEPACK, comments
about usage documentation, and notes about installing the collection.
For additional documentation on ODEPACK, see also the papers [1], [2]
(for LSODE), and [3] (for LSODPK and LSODKR), and in the references
cited there.  (However, the document [2] does not reflect the upgrade
operations described above.)

### References:

1.  A. C. Hindmarsh, "[ODEPACK, A Systematized Collection of ODE Solvers](https://computation.llnl.gov/casc/nsde/pubs/u88007.pdf),"
     in Scientific Computing, R. S. Stepleman et al. (eds.), North-Holland,
     Amsterdam, 1983 (vol. 1 of IMACS Transactions on Scientific Computation),
     pp. 55-64.

2.  K. Radhakrishnan and A. C. Hindmarsh, "[Description and Use of LSODE,
     the Livermore Solver for Ordinary Differential Equations](https://computation.llnl.gov/casc/nsde/pubs/u113855.pdf)," LLNL
     report UCRL-ID-113855, December 1993.

3.  P. N. Brown and A. C. Hindmarsh, "[Reduced Storage Matrix Methods
     in Stiff ODE Systems](https://computation.llnl.gov/casc/nsde/pubs/u95088.pdf)," J. Appl. Math. & Comp., 31 (1989), pp.40-91.

---------------------------------------------------------------------------

### I. Summary of the ODEPACK Solvers

#### A. Solvers for explicitly given systems.

For each of the following solvers, it is assumed that the ODEs are
given explicitly, so that the system can be written in the form
dy/dt = f(t,y), where y is the vector of dependent variables, and t is
the independent variable.

1. LSODE (Livermore Solver for Ordinary Differential Equations) is the
   basic solver of the collection.  It solves stiff and nonstiff systems
   of the form dy/dt = f.  In the stiff case, it treats the Jacobian
   matrix df/dy as either a dense (full) or a banded matrix, and as either
   user-supplied or internally approximated by difference quotients.
   It uses Adams methods (predictor-corrector) in the nonstiff case,
   and Backward Differentiation Formula (BDF) methods (the Gear methods)
   in the stiff case.  The linear systems that arise are solved by direct
   methods (LU factor/solve).  LSODE supersedes the older GEAR and GEARB
   packages, and reflects a complete redesign of the user interface
   and internal organization, with some algorithmic improvements.

2. LSODES, written jointly with A. H. Sherman, solves systems dy/dt = f
   and in the stiff case treats the Jacobian matrix in general sparse
   form.  It determines the sparsity structure on its own, or optionally
   accepts this information from the user.  It then uses parts of the
   Yale Sparse Matrix Package (YSMP) to solve the linear systems that
   arise, by a sparse (direct) LU factorization/backsolve method.
   LSODES supersedes, and improves upon, the older GEARS package.

3. LSODA, written jointly with L. R. Petzold, solves systems dy/dt = f
   with a dense or banded Jacobian when the problem is stiff, but it
   automatically selects between nonstiff (Adams) and stiff (BDF)
   methods.  It uses the nonstiff method initially, and dynamically
   monitors data in order to decide which method to use.

4. LSODAR, also written jointly with L. R. Petzold, is a variant of
   LSODA with a rootfinding capability added.  Thus it solves problems
   dy/dt = f with dense or banded Jacobian and automatic method
   selection, and at the same time, it finds the roots of any of a
   set of given functions of the form g(t,y).  This is often useful
   for finding stop conditions, or for finding points at which a switch
   is to be made in the function f.

5. LSODPK, written jointly with Peter N. Brown, is a variant of LSODE
   in which the direct solvers for the linear systems have been replaced
   by a selection of four preconditioned Krylov (iterative) solvers.
   The user must supply a pair of routine to evaluate, preprocess, and
   solve the (left and/or right) preconditioner matrices.  LSODPK also
   includes an option for a user-supplied linear system solver to be used
   without Krylov iteration.

6. LSODKR is a variant of LSODPK with the addition of the same
   rootfinding capability as in LSODAR, and also of automatic switching
   between functional and Newton iteration.  The nonlinear iteration
   method-switching differs from the method-switching in LSODA and LSODAR,
   but provides similar savings by using the cheaper method in the non-stiff
   regions of the problem.  LSODKR also improves on the Krylov methods in
   LSODPK by offering the option to save and reuse the approximate Jacobian
   data underlying the preconditioner.


#### B. Solvers for linearly implicit systems.

The following solvers treat systems in the linearly implicit form
A(t,y) dy/dt = g(t,y), A = a square matrix, i.e. with the derivative
dy/dt implicit, but linearly so.  These solvers allow A to be
singular, in which case the system is a differential-algebraic
equation (DAE) system.  In that case, the user must be very careful
to supply a well-posed problem with consistent initial conditions.

7. LSODI, written jointly with J. F. Painter, solves linearly implicit
   systems in which the matrices involved (A, dg/dy, and d(A dy/dt)/dy)
   are all assumed to be either dense or banded.  LSODI supersedes the
   older GEARIB solver and improves upon it in numerous ways.

8. LSOIBT, written jointly with C. S. Kenney, solves linearly implicit
   systems in which the matrices involved are all assumed to be
   block-tridiagonal.  Linear systems are solved by the LU method.

9. LSODIS, written jointly with S. Balsdon, solves linearly implicit
   systems in which the matrices involved are all assumed to be sparse.
   Like LSODES, LSODIS either determines the sparsity structure or
   accepts it from the user, and uses parts of the Yale Sparse Matrix
   Package to solve the linear systems that arise, by a direct method.

---------------------------------------------------------------------------

### II. Usage Documentation

   Each of the solvers in the ODEPACK collection is headed by a
user-callable driver subroutine, with the same name as the solver
(DLSODE, etc.).  The call sequence of the driver routine includes the
names of one or more user-supplied subroutines that define the ODE
system, and various other problem and solution parameters.  Complete
user documentation is given in the initial block of comment lines
(the prologue) in the driver routine.  In each case, this prologue is
organized as follows:

 * Summary of Usage (short, for standard modes of use)
 * Example Problem, with code and output (except for LSODPK and LSODKR)
 * Full Description of User Interface, further divided as follows:
      1. Call sequence description (including optional inputs/outputs)
      2. Optionally callable routines
      3. Descriptions of internal Common blocks
      4. Optionally user-replaceable routines
 * Revision History, showing date written and dates of revisions
 * Other Routines, a list of all subordinate routines for the solver

First-time users should read only the Summary of Usage and look at the
the Example Problem (or demonstration program), then later refer to the
Full Description if and when more details or nonstandard options are needed.

---------------------------------------------------------------------------

### III. Installation Notes

##### Use of supplied matrix procedures

THe src/M_matrix/ directory  includes modified versions of routines from the
LINPACK and BLAS collections that are needed by the solvers (and by two
of the demonstration programs), for the solution of dense and banded
linear systems and associated basic linear algebra operations.
These routine are:

       _From LINPACK_ :  DGEFA, DGESL, DGBFA, DGBSL
       _From the BLAS_: DAXPY, DCOPY, DDOT, DSCAL, DNRM2, IDAMAX

If your computer system already has these routines, and especially if it
has machine-optimized versions, the copies provided in the M_module
module should probably not be called by your program if high performance
is required.

##### The message package

The second auxiliary source file directory M_da1/ includes a set of five
routines -- XERRWD, XSETUN, XSETF, IXSAV, IUMACH -- which handle error
messages from the solvers.

These routines are slated for elimination and replacement with a 
more intuitive interface.

This set is in fact a reduced version (sufficient for the needs of
ODEPACK) of a much larger error handling package from the SLATEC Library.
If your computer system already has the full SLATEC error handler, the
version provided here can be ignored.  If the reduced version is used,
its machine-dependent features should be checked first; see comments in
Subroutine XERRWD.

##### Non-standard code

4. ODEPACK contains a few instances where the Fortran Standard is violated:

   (a) In various places in the LSODES and LSODIS solvers, a call to a
   subroutine has a subscripted real array as an argument where the
   subroutine called expects an integer array.  Calls of this form
   occur in Subroutine DLSODES (to DSTODE), in DIPREP (to DPREP),
   in Subroutine DLSODIS (to DSTODI), and in DIPREPI (to DPREPI).
   Another such call occurs in the DLSODES demonstration program,
   from the main program to Subroutine SSOUT.  

   This is done in order
   to use work space in an efficient manner, as the same space is
   sometimes used for real work space and sometimes for integer work
   space.  If your compiler does not accept this feature, one possible
   way to get the desired result is to compile the called routines
   and calling routines in separate jobs, and then combine the binary
   modules in an appropriate manner. 

   If this procedure is still not
   acceptable under your system, it will be necessary to radically
   alter the structure of the array RWORK within the LSODES or LSODIS
   solver package.  (See also Note 5 below.)

   (b) Each ODEPACK solver treats the arguments NEQ, Y, RTOL, and ATOL
   as arrays, even though the length may be only 1.  Moreover,
   except for Y, the usage instructions say that these arguments
   may be either arrays or scalars.  If your system does not allow
   such a mismatch, then the documentation of these arguments
   should be changed accordingly.

5. For maximum storage economy, the LSODES and LSODIS solvers make use
of the real to integer wordlength ratio.  This is assumed to be an
integer L such that if a real array R and an integer array M occupy
the same space in memory, R(1) having the same bit address as M(1),
then R(I) has the same address as M((I-1)*L+1).  This ratio L is
usually 2 for double precision, and this is the value used in the
double precision version supplied.  If this value is incorrect, it
needs to be changed in two places:

  (a) The integer LENRAT is DATA-loaded in Subroutines DLSODES and
      DLSODIS to this ratio, shortly below the prologue.

  (b) The integer LRATIO is DATA-loaded in Subroutine CDRV to this
      ratio, shortly below the prologue of that routine.

(See comments in both places.)  If the ratio is not an integer, use
the greatest integer not exceeding the ratio.

6. For installation of ODEPACK on a Cray computer, the source files
supplied include compiler directives for the CFT compiler.  These have
the form CDIR$ IVDEP and occur prior to certain loops that involve
subscript shifts (and would otherwise not be vectorized).  These
directives are (or should be) treated as comments by any other compiler.

7. On first obtaining ODEPACK, the demonstration programs should be
compiled and executed prior to any other use of the solvers.  In most
cases, these exercise all of the major method options in each solver,
and are self-checking.  (In the case of LSODPK and LSODKR, the
demonstration programs are not self-checking, and for LSODKR only one
major method option is used.)  In any case, the output can be compared
with the sample output supplied, which was generated from the double
precision version of ODEPACK on a 32-bit computer.  When comparing your
output with that supplied, differences of 10-20% in the final values of
the various statistical counters can be attributed to differences in
the roundoff properties of different computer systems.

8. If some subset of the whole ODEPACK collection is desired, without
unneeded routines, the appropriate routines must be extracted
accordingly.  The following lists give the routines needed for the
double precision version of each solver.

    The DLSODE solver consists of the routines
  DLSODE, DINTDY, DSTODE, DCFODE, DPREPJ, DSOLSY, DEWSET, DVNORM, DSRCOM,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

    The DLSODES solver consists of the routines
  DLSODES, DIPREP, DPREP, JGROUP, ADJLR, CNTNZU, DINTDY, DSTODE, DCFODE,
  DPRJS, DSOLSS, DEWSET, DVNORM, DSRCMS,
  ODRV, MD, MDI, MDM, MDP, MDU, SRO,
  CDRV, NROC, NSFC, NNFC, NNSC, NNTC,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

    The DLSODA solver consists of the routines
  DLSODA, DINTDY, DSTODA, DCFODE, DPRJA, DSOLSY, DEWSET,
  DMNORM, DFNORM, DBNORM, DSRCMA,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODAR solver consists of the routines
  DLSODAR, DRCHEK, DROOTS, DINTDY, DSTODA, DCFODE, DPRJA, DSOLSY, DEWSET,
  DMNORM, DFNORM, DBNORM, DSRCAR,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, DCOPY, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODPK solver consists of the routines
  DLSODPK, DINTDY, DEWSET, DVNORM, DSTODPK, DCFODE, DPKSET, DSOLPK,
  DSPIOM, DATV, DORTHOG, DHEFA, DHESL, DSPIGMR, DHEQR, DHELS,
  DPCG, DPCGS, DATP, DUSOL, DSRCPK,
  DAXPY, DSCAL, DCOPY, DDOT, DNRM2, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODKR solver consists of the routines
  DLSODKR, DRCHEK, DROOTS, DLHIN, DINTDY, DEWSET, DVNORM, DSTOKA,
  DCFODE, DSETPK, DSOLPK, DSPIOM, DATV, DORTHOG, DHEFA, DHESL, DSPIGMR,
  DHEQR, DHELS, DPCG, DPCGS, DATP, DUSOL, DSRCKR,
  DAXPY, DSCAL, DCOPY, DDOT, DNRM2, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODI solver consists of the routines
  DLSODI, DAINVG, DINTDY, DSTODI, DCFODE, DPREPJI, DSOLSY, DEWSET,
  DVNORM, DSRCOM,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSOIBT solver consists of the routines
  DLSOIBT, DAIGBT, DINTDY, DSTODI, DCFODE, DPJIBT, DSLSBT, DEWSET,
  DVNORM, DSRCOM, DDECBT, DSOLBT,
  DGEFA, DGESL, DAXPY, DSCAL, DDOT, IDAMAX,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODIS solver consists of the routines
  DLSODIS, DAINVGS, DIPREPI, DPREPI, JGROUP, ADJLR, CNTNZU, DINTDY,
  DSTODI, DCFODE, DPRJIS, DSOLSS, DEWSET, DVNORM, DSRCMS,
  ODRV, MD, MDI, MDM, MDP, MDU, SRO,
  CDRV, NROC, NSFC, NNFC, NNSC, NNTC,
  XERRWD, XSETUN, XSETF, IXSAV, IUMACH
