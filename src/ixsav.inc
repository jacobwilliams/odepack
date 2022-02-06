!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!!### Name
!! ixsav(3f) - [M_odesave::matrix] Save and recall error message control parameters.
!!
!!### Synopsis
!!    function ixsav(Ipar,Ivalue,Iset)
!!    integer            :: ixsav
!!
!!    integer,intent(in) :: Ipar
!!    integer,intent(in) :: Ivalue
!!    logical,intent(in) :: Iset
!!
!!### Description
!!   IXSAV always recalls the _current_ value of one of two error message
!!   control flags.
!!
!!   It also can set new values to return for subsequent queries if
!!   Iset=.true.
!!
!!   This is a modification of the SLATEC library routine J4SAVE.
!!
!!### Saved local variables..
!!
!!   LUNIT
!!
!!   :   Logical unit number for messages.
!!
!!   MESFLG
!!
!!   :   Print control flag..
!!
!!           1 means print all messages (the default).
!!           0 means no message printing.
!!
!!### Options
!!
!!   IPAR
!!
!!   :   Parameter indicator (1 for LUNIT, 2 for MESFLG).
!!
!!   IVALUE
!!
!!   :   The value to be set for the parameter, if ISET = .TRUE.
!!
!!   ISET
!!
!!   :   Logical flag to indicate whether to read or write.
!!       If ISET = .TRUE., the parameter will be given
!!       the value IVALUE.  If ISET = .FALSE., the parameter
!!       will be unchanged, and IVALUE is a dummy argument.
!!
!!### Returns
!!
!!   IXSAV
!!   : The (old) value of the parameter.
!!
! -----------------------------------------------------------------------
! ### SUBSIDIARY
! ### PURPOSE  Save and recall error message control parameters.
! ### CATEGORY  R3C
! ### TYPE      ALL (IXSAV-A)
! ### AUTHOR  Hindmarsh, Alan C., (LLNL)
! ### DESCRIPTION
! ### SEE ALSO  XERRWD, XERRWV
! ### REVISION HISTORY  (YYMMDD)
!     19921118   DATE WRITTEN
!     19930329   Modified prologue to SLATEC format. (FNF)
!     19930922   Minor cosmetic changes. (FNF)
!     20220125   Changed to use ISO_FORTRAN_ENV (JSU)
!
!  Subroutines called by IXSAV.. None
!  Function routine called by IXSAV.. None
! -----------------------------------------------------------------------
function ixsav(Ipar,Ivalue,Iset)
use,intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                         stdout=>output_unit, &
                                         stderr=>error_unit

integer            :: ixsav

integer,intent(in) :: Ipar
integer,intent(in) :: Ivalue
logical,intent(in) :: Iset

integer,save :: lunit = stdout
integer,save :: mesflg = 1

   select case(Ipar)

   case(1) ! lunit
      ixsav = lunit
      if ( Iset ) lunit = Ivalue

   case(2) ! mesflg
      ixsav = mesflg
      if ( Iset ) mesflg = Ivalue

   case default
      stop '*ixsav* unknown value for Ipar'

   end select

end function ixsav
