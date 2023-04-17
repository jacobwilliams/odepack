
c-----------------------------------------------------------------------
c Demonstration program for the DLSODIS package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c This program solves a semi-discretized form of the Burgers equation,
c
c     u  = -(u*u/2)  + eta * u
c      t           x          xx
c
c for  -1 .le. x .le. 1, t .ge. 0.
c Here eta = 0.05.
c Boundary conditions: u(-1,t) = u(1,t) and du/dx(-1,t) = du/dx(1,t).
c Initial profile: square wave
c     u(0,x) = 0    for 1/2 .lt. abs(x) .le. 1
c     u(0,x) = 1/2  for abs(x) = 1/2
c     u(0,x) = 1    for 0 .le. abs(x) .lt. 1/2
c
c An ODE system is generated by a simplified Galerkin treatment
c of the spatial variable x.
c
c Reference:
c R. C. Y. Chin, G. W. Hedstrom, and K. E. Karlsson,
c A Simplified Galerkin Method for Hyperbolic Equations,
c Math. Comp., vol. 33, no. 146 (April 1979), pp. 647-658.
c
c The problem is run with the DLSODIS package with a 12-node mesh,
c for various appropriate values of the method flag mf.
c Output is on unit lout, set to 6 in a data statement below.
c
c Problem specific data:
c npts  = number of unknowns (npts = 0 mod 4)
c nnz   = number of non-zeros in Jacobian before fill in
c nnza  = number of non-zeros in Jacobian after fill in
c lrwk  = length of real work array (taking into account fill in)
c liwk  = length of integer work array
c ipia  = pointer to ia in iw (ia(j) = iw(ipia+j-1)
c ipja  = pointer to ja in iw (ja(j) = iw(ipja+j-1)
c ipic  = pointer to ic in iw array (ic(j) = iw(ipic+j-1))
c ipjc  = pointer to jc in iw array (jc(j) = iw(ipjc+j-1))
c-----------------------------------------------------------------------
      parameter (npts = 12, nnz = 3*npts, nnza = 5*npts)
      parameter (lrwk = 20+3*nnza+28*npts)
      parameter (liwk = 32+2*nnza+2*npts)
      parameter (ipia = 31, ipja = 31+npts+1)
      parameter (ll = 32+npts+nnz, ipic = ll,ipjc = ll+npts+1)
c
      external res, addasp, jacsp
      integer i, io, istate, itol, iw, j, lout, lrw, liw, lyh,
     1   meth, miter, mf, moss, n, nout, nerr, n14, n34
      double precision eta, delta, r4d, eodsq, a, b, t, tout, tlast,
     1    tinit, errfac
      double precision atol, rtol, rw, y, ydoti, elkup
      double precision zero, fourth, half, one, hun
      dimension y(npts), ydoti(npts), tout(4), atol(2), rtol(2)
      dimension rw(lrwk), iw(liwk)
c Pass problem parameters in the Common block test1.
      common /test1/ r4d, eodsq
c
c Set problem parameters and run parameters
      data eta/0.05d0/, a/-1.0d0/, b/1.0d0/
      data zero/0.0d0/, fourth/0.25d0/, half/0.5d0/, one/1.0d0/,
     1   hun/100.0d0/
      data tinit/0.0d0/, tlast/0.4d0/
      data tout/0.10d0, 0.20d0, 0.30d0, 0.40d0/
      data lout/6/, nout/4/
      data itol/1/, rtol/1.0d-3, 1.0d-6/, atol/1.0d-3, 1.0d-6/
c
      nerr = 0
      lrw = lrwk
      liw = liwk
c
c Compute the mesh width delta and other parameters.
      delta = (b - a)/npts
      r4d = fourth/delta
      eodsq = eta/delta**2
      n14 = npts/4 + 1
      n34 = 3 * (npts/4) + 1
      n = npts
c
c Set the initial profile (for output purposes only).
      do 10 i = 1,n
   10   y(i) = zero
      y(n14) = half
      do 20 i = n14+1,n34-1
   20   y(i) = one
      y(n34) = half
c
      write(lout,1000)
      write(lout,1100) eta,a,b,tinit,tlast,n
      write(lout,1200) (y(i), i=1,n)
c
c Set the initial sparse data structures for coefficient matrix A
c and the Jacobian matrix C
      call struct(iw(ipia), iw(ipja), iw(ipic), iw(ipjc), n)
c
c The j loop is over error tolerances.
      do 200 j = 1,2
c
c This method flag loop is for demonstration only.
      do 200 moss = 0,4
        do 100 meth = 1,2
          do 100 miter = 1,2
   35       mf = 100*moss + 10*meth + miter
c
c Set the initial profile.
            do 40 i = 1,n
   40         y(i) = zero
            y(n14) = half
            do 50 i = n14+1,n34-1
   50         y(i) = one
            y(n34) = half
c
            t = tinit
            istate = 0
c
            write(lout,1500) itol, rtol(j), atol(j), mf
c
c output loop for each case
            do 80 io = 1,nout
c
c Call DLSODIS and print results.
              call dlsodis (res, addasp, jacsp, n, y, ydoti, t,
     1             tout(io), itol, rtol(j), atol(j), 1, istate, 0,
     2             rw, lrw, iw, liw, mf)
              write(lout,2000) t, rw(11), iw(14), (y(i), i=1,n)
c
c If istate is not 2 on return, print message and go to next case.
              if (istate .ne. 2) then
                write(lout,4000) mf, t, istate
                nerr = nerr + 1
                go to 100
                endif
   80       continue
            write(lout,3000) mf, iw(11), iw(12), iw(13),
     1                       iw(17), iw(18), iw(20), iw(21)
c
c Estimate final error and print result.
            lyh = iw(22)
            errfac = elkup(n, y, rw(lyh), itol, rtol(j), atol(j), lout)
            if (errfac .lt. hun) then
              write(lout,5000) errfac
            else
              write(lout,5001) errfac
              nerr = nerr + 1
              endif
  100     continue
  200   continue
c
      write(lout,6000) nerr
      stop
c
 1000 format(20x,' Demonstration Program for DLSODIS' )
 1100 format(//10x,'-- Simplified Galerkin solution of ',
     1       'Burgers equation --'///
     1       13x,'Diffusion coefficient is eta =',d10.2/
     1       13x,'Uniform mesh on interval',d12.3,' to ',d12.3/
     2       13x,'Periodic boundary conditions'/
     2       13x,'Initial data are as follows:'//20x,'t0 = ',d12.5/
     2       20x,'tlast = ',d12.5/20x,'n  = ',i3//)
c
 1200 format(/'Initial profile:',/20(6d12.4/))
c
 1500 format(///85('*')///'Run with itol =',i2,'  rtol =',d12.2,
     1       '  atol =',d12.2,'   mf = ',i3//)
c
 2000 format(' Output for time t =',d12.5,'  current h =',d12.5,
     1       '  current order =',i2/20(6d12.4/))
c
 3000 format(/'Final statistics for mf = ',i3,': ',
     1       i5,' steps,',i6,' res,',i6,' Jacobians,'/
     2       20x,' rw size =',i6,',    iw size =',i6/
     3       20x,i4,' extra res for each jac,',i4,' decomps')
c
 4000 format(/'Final time reached for mf = ',i3,
     1       ' was t = ',d12.5/25x,'at which istate = ',i2//)
 5000 format('Final output is correct to within ',d9.2,
     1       '  times local error tolerance.'/)
 5001 format('Final output is wrong by ',d9.2,
     1       '  times local error tolerance.'/)
 6000 format(///85('*')//
     1       'Run completed: number of errors encountered =',i3)
c
c end of main program for the DLSODIS demonstration program
      end

      subroutine struct(ia, ja, ic, jc, n)
c This subroutine computes the initial sparse data structure of
c the mass (ia,ja) and Jacobian (ic,jc) matrices.
c
      integer ia(*), ja(*), ic(*), jc(*), n,  jj, k, l, m
c
      write(6,1200)
      k = 0
      do 33 l = 1,n
         ia(l) = (l-1)*3+1
         ic(l) = (l-1)*3+1
         do 32 m = l,l+2
            k = k + 1
            ja(k) = m - 1
            jc(k) = m - 1
   32    continue
   33 continue
      ia(n+1) = 3*n + 1
      ic(n+1) = 3*n+1
      ja(1) = n
      jc(1) = n
      ja(k) = 1
      jc(k) = 1
c
      write(6,1300) (ia(jj),jj=1,n+1)
      write(6,1350) (ja(jj),jj=1,k)
      write(6,1400) (ic(jj),jj=1,n+1)
      write(6,1450) (jc(jj),jj=1,k)
      return
1200  format('Initial sparse data structures'/)
1300  format(' ia  ',15i4/10(5x,15i4/))
1350  format(' ja  ',15i4/10(5x,15i4/))
1400  format(' ic  ',15i4/10(5x,15i4/))
1450  format(' jc  ',15i4/10(5x,15i4/))
      end

      subroutine res (n, t, y, v, r, ires)
c This subroutine computes the residual vector
c   r = g(t,y) - A(t,y)*v .
c If ires = -1, only g(t,y) is returned in r, since A(t,y) does
c not depend on y.
c No changes need to be made to this routine if n is changed.
c
      integer n, ires,  i
      double precision t, y(n), v(n), r(n), r4d, eodsq, one, four, six,
     1   fact1, fact4
      common /test1/ r4d, eodsq
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      call gfun (n, t, y, r)
      if (ires .eq. -1) return
c
      fact1 = one/six
      fact4 = four/six
c
      r(1) = r(1) - (fact4*v(1) + fact1*(v(2) + v(n)))
      do 10 i = 2,n-1
         r(i) = r(i) - (fact4*v(i) + fact1*(v(i-1) + v(i+1)))
   10 continue
      r(n) = r(n) - (fact4*v(n) + fact1*(v(1) + v(n-1)))
      return
c end of subroutine res for the DLSODIS demonstration program
      end

      subroutine gfun (n, t, y, g)
c This subroutine computes the right-hand side function g(y,t).
c It uses r4d = 1/(4*delta) and eodsq = eta/delta**2
c from the Common block test1.
c
      integer n, i
      double precision t, y(n), g(n), r4d, eodsq, two
      common /test1/ r4d, eodsq
      data two/2.0d0/
c
      g(1) = r4d*(y(n)**2 - y(2)**2) + eodsq*(y(2) - two*y(1) + y(n))
c
      do 20 i = 2,n-1
        g(i) = r4d*(y(i-1)**2 - y(i+1)**2)
     1        + eodsq*(y(i+1) - two*y(i) + y(i-1))
   20   continue
c
      g(n) = r4d*(y(n-1)**2 - y(1)**2) + eodsq*(y(1)-two*y(n)+y(n-1))
c
      return
c end of subroutine gfun for the DLSODIS demonstration program
      end

      subroutine addasp (n, t, y, j, ip, jp, pa)
c This subroutine computes the sparse matrix A by columns, adds it to
c pa, and returns the sum in pa.
c The matrix A is periodic tridiagonal, of order n, with nonzero elements
c (reading across) of  1/6, 4/6, 1/6, with 1/6 in the lower left and
c upper right corners.
c
      integer n, j, ip(*), jp(*),  jm1, jp1
      double precision t, y(n), pa(n), fact1, fact4, one, four, six
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
c Compute the elements of A.
      fact1 = one/six
      fact4 = four/six
      jm1 = j - 1
      jp1 = j + 1
      if (j .eq. n) jp1 = 1
      if (j .eq. 1) jm1 = n
c
c Add the matrix A to the matrix pa (sparse).
      pa(j) = pa(j) + fact4
      pa(jp1) = pa(jp1) + fact1
      pa(jm1) = pa(jm1) + fact1
      return
c end of subroutine addasp for the DLSODIS demonstration program
      end

      subroutine jacsp (n, t, y, s, j, ip, jp, pdj)
c This subroutine computes the Jacobian dg/dy = d(g-A*s)/dy by
c columns in sparse matrix format.  Only nonzeros are loaded.
c It uses r4d = 1/(4*delta) and eodsq = eta/delta**2 from the Common
c block test1.
c
      integer n, j, ip(*), jp(*),  jm1, jp1
      double precision t, y(n), s(n), pdj(n), r4d, eodsq, two, diag, r2d
      common /test1/ r4d, eodsq
      data two/2.0d0/
c
      diag = -two*eodsq
      r2d = two*r4d
      jm1 = j - 1
      jp1 = j + 1
      if (j .eq. 1) jm1 = n
      if (j .eq. n) jp1 = 1
c
      pdj(jm1) = -r2d*y(j) + eodsq
      pdj(j) = diag
      pdj(jp1) = r2d*y(j) + eodsq
      return
c end of subroutine jacsp for the DLSODIS demonstration program
      end

      double precision function elkup (n, y, ewt, itol, rtol,atol, lout)
c This routine looks up approximately correct values of y at t = 0.4,
c ytrue = y12 or y120 depending on whether n = 12 or 120.
c These were obtained by running with very tight tolerances.
c The returned value is
c   elkup = norm of [ (y - ytrue) / (rtol*abs(ytrue) + atol) ].
c
      integer n, itol, lout, i
      double precision y(n), ewt(n), rtol, atol, y12(12), y120(120),
     1    y120a(16), y120b(16), y120c(16), y120d(16), y120e(16),
     2    y120f(16), y120g(16), y120h(8), dvnorm
      equivalence (y120a(1),y120(1)), (y120b(1),y120(17)),
     1      (y120c(1),y120(33)), (y120d(1),y120(49)),
     1      (y120e(1),y120(65)),
     1      (y120f(1),y120(81)), (y120g(1),y120(97)),
     1      (y120h(1),y120(113))
      data y12/
     1 1.60581860d-02, 3.23063251d-02, 1.21903380d-01, 2.70943828d-01,
     1 4.60951522d-01, 6.57571216d-01, 8.25154453d-01, 9.35644796d-01,
     1 9.90167557d-01, 9.22421221d-01, 5.85764902d-01, 1.81112615d-01/
      data y120a /
     1 1.89009068d-02, 1.63261891d-02, 1.47080563d-02, 1.39263623d-02,
     1 1.38901341d-02, 1.45336989d-02, 1.58129308d-02, 1.77017162d-02,
     1 2.01886844d-02, 2.32742221d-02, 2.69677715d-02, 3.12854037d-02,
     1 3.62476563d-02, 4.18776225d-02, 4.81992825d-02, 5.52360652d-02/
      data y120b /
     1 6.30096338d-02, 7.15388849d-02, 8.08391507d-02, 9.09215944d-02,
     1 1.01792784d-01, 1.13454431d-01, 1.25903273d-01, 1.39131085d-01,
     1 1.53124799d-01, 1.67866712d-01, 1.83334757d-01, 1.99502830d-01,
     1 2.16341144d-01, 2.33816600d-01, 2.51893167d-01, 2.70532241d-01/
      data y120c /
     1 2.89693007d-01, 3.09332757d-01, 3.29407198d-01, 3.49870723d-01,
     1 3.70676646d-01, 3.91777421d-01, 4.13124817d-01, 4.34670077d-01,
     1 4.56364053d-01, 4.78157319d-01, 5.00000270d-01, 5.21843218d-01,
     1 5.43636473d-01, 5.65330432d-01, 5.86875670d-01, 6.08223037d-01/
      data y120d /
     1 6.29323777d-01, 6.50129662d-01, 6.70593142d-01, 6.90667536d-01,
     1 7.10307235d-01, 7.29467947d-01, 7.48106966d-01, 7.66183477d-01,
     1 7.83658878d-01, 8.00497138d-01, 8.16665158d-01, 8.32133153d-01,
     1 8.46875019d-01, 8.60868691d-01, 8.74096465d-01, 8.86545273d-01/
      data y120e /
     1 8.98206892d-01, 9.09078060d-01, 9.19160487d-01, 9.28460742d-01,
     1 9.36989986d-01, 9.44763554d-01, 9.51800339d-01, 9.58122004d-01,
     1 9.63751979d-01, 9.68714242d-01, 9.73031887d-01, 9.76725449d-01,
     1 9.79811001d-01, 9.82297985d-01, 9.84186787d-01, 9.85466039d-01/
      data y120f /
     1 9.86109629d-01, 9.86073433d-01, 9.85291781d-01, 9.83673704d-01,
     1 9.81099057d-01, 9.77414704d-01, 9.72431015d-01, 9.65919133d-01,
     1 9.57609585d-01, 9.47193093d-01, 9.34324619d-01, 9.18631922d-01,
     1 8.99729965d-01, 8.77242371d-01, 8.50830623d-01, 8.20230644d-01/
      data y120g /
     1 7.85294781d-01, 7.46035145d-01, 7.02662039d-01, 6.55609682d-01,
     1 6.05541326d-01, 5.53327950d-01, 4.99999118d-01, 4.46670394d-01,
     1 3.94457322d-01, 3.44389410d-01, 2.97337561d-01, 2.53964948d-01,
     1 2.14705729d-01, 1.79770169d-01, 1.49170367d-01, 1.22758681d-01/
      data y120h /
     1 1.00271052d-01, 8.13689920d-02, 6.56761515d-02, 5.28075160d-02,
     1 4.23908624d-02, 3.40811650d-02, 2.75691506d-02, 2.25853507d-02/
c
      if ((n-12)*(n-120) .ne. 0) go to 300
      if (n .eq. 120) go to 100
c
c Compute local error tolerance using correct y (n = 12).
      call dewset(n, itol, rtol, atol, y12, ewt)
c
c Invert ewt and replace y by the error, y - ytrue.
      do 20  i = 1, 12
        ewt(i) = 1.0d0/ewt(i)
 20     y(i) = y(i) - y12(i)
      go to 200
c
c Compute local error tolerance using correct y (n = 120).
 100  call dewset( n, itol, rtol, atol, y120, ewt )
c
c Invert ewt and replace y by the error, y - ytrue.
      do 120  i = 1, 120
        ewt(i) = 1.0d0/ewt(i)
 120    y(i) = y(i) - y120(i)
c
c Find weighted norm of the error.
 200  elkup = dvnorm (n, y, ewt)
      return
c
c error return
 300  write(lout,400) n
      elkup = 1.0d3
 400  format(/5x,'Illegal use of elkup for n =',i4)
      return
c end of function elkup for the DLSODIS demonstration program
      end

