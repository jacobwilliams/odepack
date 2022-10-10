module M_odepack
!
! dcopy can be replaced with array syntax
! SPAG removes save attribute and type from f90 function declaration
! SPAG removes external attribute
! determine what can be private and public
! many equality tests for real values; could replace
! simplify checkpoint capability, maybe TRANSFER()
! look for SAVEs that should be parameters
! remaining GOTOs
! replace dnrm2 with intrinsic norm2(3f) ?
! in production probably would make matrix routines separate so could use alternate library of BLAS matrix procedures
! same arrays passed on single call all over the place
! NOTE: DLSODKR did not initialize JROOT and neither did lsodkr.f90 example
!
implicit none
integer,parameter :: dp=kind(0.0d0)

private
!----------------------------------------------------------------------------------------------------------------------------------!
! da1 routines
PUBLIC  :: DEWSET
PUBLIC  :: DINTDY
PUBLIC  :: DSRCAR
PUBLIC  :: DSRCKR
PUBLIC  :: DSRCMA
PUBLIC  :: DSRCMS
PUBLIC  :: DSRCOM
PUBLIC  :: DSRCPK
PUBLIC  :: DUMACH ! COMPUTE THE UNIT ROUNDOFF OF THE MACHINE.
PUBLIC  :: DVNORM
!-PUBLIC  :: DSTODE

!-- maybe private once modularized
public :: dmnorm
public :: dpjibt
public :: dprepj
public :: dprepji
public :: dsolss
public :: dprja
public :: dprjis
public :: dprjs
public :: dslsbt
public :: dsolsy
public :: adjlr
public :: cdrv
public :: cntnzu
public :: daigbt
public :: dainvg
public :: dcfode ! Set ODE integrator coefficients.
public :: diprep
public :: diprepi
public :: dlhin
public :: drchek
public :: dstoda
public :: dstodpk
public :: dstoka
public :: jgroup
public :: odrv
!-- see nothing in documentation about being public
private :: datp
private :: datv
private :: dbnorm
private :: ddecbt
private :: dfnorm
private :: dhefa
private :: dhels
private :: dheqr
private :: dhesl
private :: dorthog
private :: dpcg
private :: dpcgs
private :: dpkset
private :: droots
private :: dsetpk
private :: dsolbt
private :: dsolpk
private :: dspigmr
private :: dspiom
private :: dusol
private :: md
private :: mdi
private :: mdm
private :: mdp
private :: mdu
private :: nnfc
private :: nnsc
private :: nntc
private :: nroc
private :: nsfc
private :: sro
private :: dnrm2
!-private :: dstodi
!-private :: dprepi
!-private :: dainvgs
!-private :: dprep
!----------------------------------------------------------------------------------------------------------------------------------!
! matrix routines
PUBLIC  :: DAXPY
PUBLIC  :: XSETUN
PUBLIC  :: XSETF
PUBLIC  :: DSCAL
PUBLIC  :: DCOPY
PUBLIC  :: XERRWD
PUBLIC  :: DGEFA
PUBLIC  :: DGESL
! should any matrix procedures be private?
private :: idamax
private :: ixsav
private :: ddot
private :: dgbsl
private :: dgbfa
! main routines
PUBLIC :: DLSODIS
PUBLIC :: DLSODI
PUBLIC :: DLSOIBT
PUBLIC :: DLSODKR
PUBLIC :: DLSODPK
PUBLIC :: DLSODE
PUBLIC :: DLSODA
PUBLIC :: DLSODES
PUBLIC :: DLSODAR
!----------------------------------------------------------------------------------------------------------------------------------!
! state save and restore

type dls002
   real(kind=dp) :: stifr
   integer       :: newt, nsfi, nslj, njev
end type dls002
type(dls002),public,save :: dls

type dlss01
   real(kind=dp) :: con0, conmin, ccmxj, psmall, rbig, seth
   integer       :: iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp
   integer       :: ipa, lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj, nslj, ngp, nlu, nnz, nsp, nzl, nzu
end type dlss01
type(dlss01),public,save :: dlss

type dlsa01
   real(kind=dp) :: tsw, cm1(12), cm2(5), pdest, pdlast, ratio, pdnorm
   integer       :: insufr, insufi, ixpr, icount, irflag, jtyp, mused, mxordn, mxords
end type dlsa01
type(dlsa01),public,save :: dlsa

type dlsr01
   real(kind=dp) :: alpha, x2, t0, tlast, toutc
   integer       :: lg0, lg1, lgx, imax, last, irfnd, itaskc, ngc, nge
end type dlsr01
type(dlsr01),public,save :: dlsr

type dlpk01
   real(kind=dp) :: delt, epcon, sqrtn, rsqrtn
   integer       :: jpre, jacflg, locwp, lociwp, lsavx, kmp, maxl, mnewt, nni, nli, nps, ncfn, ncfl
   end type dlpk01
type(dlpk01),public,save :: dlpk

type dls001
   real(kind=dp) :: conit, crate, el(13), elco(13,12), hold, rmax, tesco(3,12), ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer       :: init,   mxstep,  mxhnil,  nhnil,  nslast,  nyh,   ialth,   ipup,    lmax,    meo
   integer       :: nqnyh,  nslp,    icf,     ierpj,  iersl,   jcur,  jstart,  kflag,   l,       lyh
   integer       :: lewt,   lacor,   lsavf,   lwm,    liwm,    meth,  miter,   maxord,  maxcor,  msbp
   integer       :: mxncf,  n,       nq,      nst,    nfe,     nje,   nqu
end type dls001
type(dls001),public,save :: dls1

private :: return_dls1_real
private :: return_dls1_int
private :: set_dls1_real
private :: set_dls1_int

contains
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!

function return_dls1_real()
real(kind=dp),allocatable :: return_dls1_real(:)

   return_dls1_real=[ dls1%conit, dls1%crate, dls1%el(:), dls1%elco(:,:), dls1%hold, &
                    & dls1%rmax, dls1%tesco(:,:), dls1%ccmax, dls1%el0, dls1%h, dls1%hmin, &
                    & dls1%hmxi, dls1%hu, dls1%rc, dls1%tn, dls1%uround ]

end function return_dls1_real

function return_dls1_int()
integer,allocatable :: return_dls1_int(:)

   return_dls1_int=[ dls1%init, dls1%mxstep, dls1%mxhnil, dls1%nhnil, dls1%nslast, dls1%nyh, dls1%ialth, &
                   & dls1%ipup, dls1%lmax, dls1%meo, dls1%nqnyh, dls1%nslp, dls1%icf, dls1%ierpj, &
                   & dls1%iersl, dls1%jcur, dls1%jstart, dls1%kflag, dls1%l, dls1%lyh, dls1%lewt, &
                   & dls1%lacor, dls1%lsavf, dls1%lwm, dls1%liwm, dls1%meth, dls1%miter, dls1%maxord, &
                   & dls1%maxcor, dls1%msbp, dls1%mxncf, dls1%n, dls1%nq, dls1%nst, dls1%nfe, &
                   & dls1%nje, dls1%nqu ]

end function return_dls1_int

subroutine set_dls1_real(r)
real(kind=dp),intent(in) :: r(218)
   dls1%conit        =  r(1)
   dls1%crate        =  r(2)
   dls1%el           =  r(3:15)
   dls1%elco         =  reshape(r(16:16+156-1),[13,12])
   dls1%hold         =  r(172)
   dls1%rmax         =  r(173)
   dls1%tesco        =  reshape(r(174:174+36-1),[3,12])
   dls1%ccmax        =  r(210)
   dls1%el0          =  r(211)
   dls1%h            =  r(212)
   dls1%hmin         =  r(213)
   dls1%hmxi         =  r(214)
   dls1%hu           =  r(215)
   dls1%rc           =  r(216)
   dls1%tn           =  r(217)
   dls1%uround       =  r(218)

end subroutine set_dls1_real

subroutine set_dls1_int(int)
integer,intent(in) :: int(37)
   dls1%init    =  int(1)
   dls1%mxstep  =  int(2)
   dls1%mxhnil  =  int(3)
   dls1%nhnil   =  int(4)
   dls1%nslast  =  int(5)
   dls1%nyh     =  int(6)
   dls1%ialth   =  int(7)
   dls1%ipup    =  int(8)
   dls1%lmax    =  int(9)
   dls1%meo     =  int(10)
   dls1%nqnyh   =  int(11)
   dls1%nslp    =  int(12)
   dls1%icf     =  int(13)
   dls1%ierpj   =  int(14)
   dls1%iersl   =  int(15)
   dls1%jcur    =  int(16)
   dls1%jstart  =  int(17)
   dls1%kflag   =  int(18)
   dls1%l       =  int(19)
   dls1%lyh     =  int(20)
   dls1%lewt    =  int(21)
   dls1%lacor   =  int(22)
   dls1%lsavf   =  int(23)
   dls1%lwm     =  int(24)
   dls1%liwm    =  int(25)
   dls1%meth    =  int(26)
   dls1%miter   =  int(27)
   dls1%maxord  =  int(28)
   dls1%maxcor  =  int(29)
   dls1%msbp    =  int(30)
   dls1%mxncf   =  int(31)
   dls1%n       =  int(32)
   dls1%nq      =  int(33)
   dls1%nst     =  int(34)
   dls1%nfe     =  int(35)
   dls1%nje     =  int(36)
   dls1%nqu     =  int(37)

end subroutine set_dls1_int
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
! da1 routines
include "M_da1/dumach.inc"
include "M_da1/dcfode.inc"
include "M_da1/dmnorm.inc"
include "M_da1/dfnorm.inc"
include "M_da1/dbnorm.inc"
include "M_da1/dvnorm.inc"
include "M_da1/sro.inc"
include "M_da1/adjlr.inc"
include "M_da1/dewset.inc"
include "M_da1/jgroup.inc"
include "M_da1/odrv.inc"
include "M_da1/dhels.inc"
include "M_da1/cntnzu.inc"
include "M_da1/daigbt.inc"
include "M_da1/ddecbt.inc"
include "M_da1/dslsbt.inc"
include "M_da1/dsolbt.inc"
include "M_da1/md.inc"
include "M_da1/mdi.inc"
include "M_da1/mdm.inc"
include "M_da1/mdp.inc"
include "M_da1/mdu.inc"
include "M_da1/nnsc.inc"
include "M_da1/nntc.inc"
include "M_da1/cdrv.inc"
include "M_da1/dainvg.inc"
include "M_da1/dheqr.inc"
include "M_da1/dhesl.inc"
include "M_da1/dorthog.inc"
include "M_da1/dusol.inc"
include "M_da1/dhefa.inc"
include "M_da1/dpcg.inc"
include "M_da1/dspiom.inc"
include "M_da1/dpcgs.inc"
include "M_da1/nsfc.inc"
include "M_da1/dspigmr.inc"
include "M_da1/nroc.inc"
include "M_da1/nnfc.inc"
include "M_da1/dintdy.inc"
include "M_da1/droots.inc"
include "M_da1/dsolss.inc"
include "M_da1/dsolsy.inc"
include "M_da1/dsrcar.inc"
include "M_da1/dsrckr.inc"
include "M_da1/dsrcma.inc"
include "M_da1/dsrcms.inc"
include "M_da1/dsrcom.inc"
include "M_da1/dsrcpk.inc"
include "M_da1/dlhin.inc"
include "M_da1/diprep.inc"
include "M_da1/diprepi.inc"
include "M_da1/dprepji.inc"
include "M_da1/drchek.inc"
include "M_da1/dsetpk.inc"
include "M_da1/datv.inc"
include "M_da1/dprjis.inc"
include "M_da1/dprja.inc"
include "M_da1/dpjibt.inc"
include "M_da1/dpkset.inc"
include "M_da1/dsolpk.inc"
include "M_da1/dprjs.inc"
include "M_da1/datp.inc"
include "M_da1/dprepj.inc"
include "M_da1/dstoda.inc"
include "M_da1/dstoka.inc"
include "M_da1/dstodpk.inc"
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
! main routines
include "M_main/dlsodis.inc"
include "M_main/dlsodi.inc"
include "M_main/dlsoibt.inc"
include "M_main/dlsodkr.inc"
include "M_main/dlsodpk.inc"
include "M_main/dlsode.inc"
include "M_main/dlsoda.inc"
include "M_main/dlsodes.inc"
include "M_main/dlsodar.inc"
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
! matrix routines
include "M_matrix/daxpy.inc"
include "M_matrix/xsetun.inc"
include "M_matrix/xsetf.inc"
include "M_matrix/dscal.inc"
include "M_matrix/idamax.inc"
include "M_matrix/dcopy.inc"
include "M_matrix/ixsav.inc"
include "M_matrix/ddot.inc"
include "M_matrix/xerrwd.inc"
include "M_matrix/dgefa.inc"
include "M_matrix/dgesl.inc"
include "M_matrix/dgbsl.inc"
include "M_matrix/dgbfa.inc"
include "M_matrix/dnrm2.inc"
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!-include "M_da1/dprep.inc"
!-include "M_da1/dainvgs.inc"
!-include "M_da1/dprepi.inc"
!-include "M_da1/dstodi.inc"
!-include "M_da1_/dstode.inc"
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
end module M_odepack
!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!- TYPE MISMATCH
include "M_da1/dprep.inc"
include "M_da1/dainvgs.inc"
include "M_da1/dprepi.inc"
include "M_da1/dstodi.inc"
include "M_da1/dstode.inc"
