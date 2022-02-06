!----------------------------------------------------------------------------------------------------------------------------------!
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()!
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!! This subroutine computes the initial value of the vector YDOT
!! satisfying
!!
!!      A * YDOT = g(t,y)
!!
!! when A is nonsingular.  It is called by DLSODIS for initialization
!! only, when ISTATE = 0.  The matrix A is subjected to LU
!! decomposition in CDRV.  Then the system A*YDOT = g(t,y) is solved
!! in CDRV.
!!
!! In addition to variables described previously, communication
!! with DAINVGS uses the following:
!!
!! Y
!!
!! : array of initial values.
!!
!! WK
!!
!! : real work space for matrices.  On output it contains A and
!!   its LU decomposition.  The LU decomposition is not entirely
!!   sparse unless the structure of the matrix A is identical to
!!   the structure of the Jacobian matrix dr/dy.
!!   Storage of matrix elements starts at WK(3).
!!   WK(1) = SQRT(UROUND), not used here.
!!
!! IWK
!!
!! : integer work space for matrix-related data, assumed to
!!   be equivalenced to WK.  In addition, WK(IPRSP) and WK(IPISP)
!!   are assumed to have identical locations.
!!
!! TEM
!!
!! : vector of work space of length N (ACOR in DSTODI).
!!
!! YDOT
!!
!! : output vector containing the initial dy/dt. YDOT(i) contains
!!   dy(i)/dt when the matrix A is non-singular.
!!
!! IER
!!
!! : output error flag with the following values and meanings:
!!       = 0  if DAINVGS was successful.
!!       = 1  if the A-matrix was found to be singular.
!!       = 2  if RES returned an error flag IRES = IER = 2.
!!       = 3  if RES returned an error flag IRES = IER = 3.
!!       = 4  if insufficient storage for CDRV (should not occur here).
!!       = 5  if other error found in CDRV (should not occur here).
!-----------------------------------------------------------------------
subroutine dainvgs(Neq,T,Y,Wk,Iwk,Tem,Ydot,Ier,res,adda)
Use M_odepack
implicit none
integer,parameter :: dp=kind(0.0d0)

integer                     :: Neq
real(kind=dp)               :: T
real(kind=dp)               :: Y(*)
real(kind=dp)               :: Wk(*)
integer                     :: Iwk(*)
real(kind=dp),intent(inout) :: Tem(*)
real(kind=dp)               :: Ydot(*)
integer, intent(inout)     :: Ier
external                    :: res
external                    :: adda

integer :: i, imul, j, k, kmax, kmin

   do i = 1, dlss%nnz
      Wk(dlss%iba+i) = 0.0D0
   enddo

   Ier = 1
   call res(Neq,T,Y,Wk(dlss%ipa),Ydot,Ier)
   if ( Ier>1 )then
      return
   endif

   kmin = Iwk(dlss%ipian)

   do j = 1, Neq
      kmax = Iwk(dlss%ipian+j) - 1

      do k = kmin, kmax
         i = Iwk(dlss%ibjan+k)
         Tem(i) = 0.0D0
      enddo

      call adda(Neq,T,Y,j,Iwk(dlss%ipian),Iwk(dlss%ipjan),Tem)

      do k = kmin, kmax
         i = Iwk(dlss%ibjan+k)
         Wk(dlss%iba+k) = Tem(i)
      enddo

      kmin = kmax + 1
   enddo

   dlss%nlu = dlss%nlu + 1
   Ier = 0

   do i = 1, Neq
      Tem(i) = 0.0D0
   enddo

   ! Numerical factorization of matrix A. ---------------------------------
   call cdrv(Neq,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic),Iwk(dlss%ipian), &
    & Iwk(dlss%ipjan),Wk(dlss%ipa),Tem,Tem,dlss%nsp, &
    & Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,2,dlss%iys)

   if ( dlss%iys==0 ) then
      !
      ! Solution of the linear system. ---------------------------------------
      call cdrv(Neq,Iwk(dlss%ipr),Iwk(dlss%ipc),Iwk(dlss%ipic),Iwk(dlss%ipian), &
       & Iwk(dlss%ipjan),Wk(dlss%ipa),Ydot,Ydot,dlss%nsp,        &
       & Iwk(dlss%ipisp),Wk(dlss%iprsp),dlss%iesp,4,dlss%iys)
      if ( dlss%iys/=0 ) Ier = 5

   else
      imul = (dlss%iys-1)/Neq
      Ier = 5
      if ( imul==8 ) Ier = 1
      if ( imul==10 ) Ier = 4
   endif

end subroutine dainvgs
