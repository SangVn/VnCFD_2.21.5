cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...File cdnozg.f                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program cdnozg

c        Computes the grid for a converging-diverging nozzle (CD Nozzle)
c        as described in AIAA-87-0355 by M.S. Liou.  Geometry and grid 
c        are in units of inches.
c
c-----------------------------------------------------------------------

      Implicit none

      Integer nim, njm, nkm

      Parameter ( nim  = 501 ) 
      Parameter ( njm  = 151 ) 
      Parameter ( nkm  =  31 ) 

      Integer idim
      Integer i, j, k
      Integer ni, nj, nk

      Real farea
      Real pi
      Real xx, ax, rx 
      Real spwall
      Real thetnz, thet
      Real sprmi, sprmax

      Real x (nim,njm,nkm)
      Real y (nim,njm,nkm)
      Real z (nim,njm,nkm)
      Real s (njm)

c-----------------------------------------------------------------------

c...Header. 

      write(*,*) ' '
      write(*,*) '--- cdnozg ---' 

c...Constants.

      pi = 4.0 * atan(1.0)

c...Specify dimension of domain.

      write(*,*) ' '
      write(*,*) 'Enter domain dimensions'
      write(*,*) '(1) 2D, (2) Axisymmetric, (3) Quasi-3D, (4) 3D'
      read (*,*) idim 

c...Specify the dimensions of the grid.  For 3D (idim=4), a 5 degree 
c...section of the nozzle is swept out.

      write(*,*) ' '

      if ( idim .eq. 4 ) then
        write(*,*) 'Enter ni, nj, nk'
        read (*,*) ni, nj, nk
        thetnz = 5.0 * pi / 180.0 
      else
        nk = 1
        write(*,*) 'Enter ni, nj'
        read (*,*) ni, nj
      endif

c...Specify the wall spacing.

      write(*,*) ' '
      write(*,*) 'Enter the wall spacing (inches)'
      read (*,*) spwall

c...Generate the grid. 

      open ( unit=7, file='rx.d', form='formatted', status='unknown' )

      sprmax = 0.0

      do  i = 1, ni

        xx = 10.0 * ( i - 1.0 ) / ( ni - 1.0 )
        ax = farea( xx )
        rx = sqrt( ax / pi )
        write(7,'(2(2x,f12.5))') xx, rx

c...    Compute the distribution of radial grid points.

        s(1)  = 0.0
        if     ( idim .eq. 1 ) then
          s(nj) = ax
          call sptanh ( nj, s, spwall, spwall, 0, 0 )
        elseif ( idim .eq. 2 ) then
          s(nj) = rx
          call sptanh ( nj, s, spwall, spwall, 0, 0 )
        elseif ( idim .eq. 3 ) then
          s(nj) = 1.0
          call speven ( nj, s )
        elseif ( idim .eq. 4 ) then
          s(nj) = rx
          call sptanh ( nj, s, spwall, spwall, 0, 0 )
        endif
    
        call gsprat ( s, nj, sprmi )
        if ( sprmi .gt. sprmax )  sprmax = sprmi

c...    Load grid.

        if     ( idim .eq. 1  .or.  idim .eq. 2 ) then
          do  j = 1, nj
            x(i,j,1) = xx
            y(i,j,1) = s(j)
            z(i,j,1) = 1.0
          enddo
        elseif ( idim .eq. 3 ) then
          do  j = 1, nj
            x(i,j,1) = xx
            y(i,j,1) = s(j)
            z(i,j,1) = rx
          enddo
        elseif ( idim .eq. 4 ) then
          do  j = 1, nj
            do  k = 1, nk
              thet     = thetnz * ( k - 1.0 ) / ( nk - 1.0 )
              x(i,j,k) = xx
              y(i,j,k) = s(j) * cos(thet)
              z(i,j,k) = s(j) * sin(thet)
            enddo
          enddo
        endif

      enddo 

c...Write out grid quality statements.

      write (*,*) ' '
      write (*,*) 'sprmax = ', sprmax

c...Write out the grid as a Plot3d file.

      write (*,*) ' '
      write (*,*) 'Writing out grid (cdnoz.txt)'

c...  open ( unit=8, file='cdnoz.z', form='unformatted', status='unknown' )
c...  write(8) c cdnoz.z --zip --binary

      open(unit=8,file='cdnoz.txt',status='unknown',form='formatted')


      write (8,*) ni, nj, nk
      write (8,*) ((( x(i,j,k), i=1,ni), j=1,nj), k=1,nk),
     &          ((( y(i,j,k), i=1,ni), j=1,nj), k=1,nk),
     &          ((( z(i,j,k), i=1,ni), j=1,nj), k=1,nk)

c...End Statement.

      write (*,*) '     '
      write (*,*) '--- End of cdnozg ---'
      write (*,*) '     '

c-----------------------------------------------------------------------

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function farea ( x )

c...Function for area of converging-divirging duct from aiaa-87-0355, 
c...M.S. liou (reference 75).
c
c-----------------------------------------------------------------------

      Implicit none
      Real farea, x, pi

c-----------------------------------------------------------------------

      pi = 4.0 * atan(1.0)

      if ( x .lt. 5.0 ) then
        farea = 1.75 - 0.75 * cos( ( 0.2 * x - 1.0 ) * pi )
      else
        farea = 1.25 - 0.25 * cos( ( 0.2 * x - 1.0 ) * pi )
      endif

c-----------------------------------------------------------------------

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...file  gsprat.f  ( grid spacing ratio )                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gsprat ( s, n, sprmax )

c     Computes the maximum ratio between two grid points for a curve.
c   
c-----------------------------------------------------------------------

      Integer i, n
      Real s(*), spa, spb, spr, sprmax

c-----------------------------------------------------------------------

      spa    = s(2) - s(1)
      spb    = s(3) - s(2)
      sprmax = max( spa, spb ) / min( spa, spb ) - 1.0

      do  i = 3, n - 1
        spa    = s(i) - s(i-1)
        spb    = s(i+1) - s(i)
        spr    = max( spa, spb ) / min( spa, spb ) - 1.0
        sprmax = max( sprmax, spr )
      enddo

c-----------------------------------------------------------------------

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...file  speven.f  ( evenly spaced grid )                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine speven ( n, s )

c     Distributes s evenly given number of points, n, and the boundary 
c     values, s(1) and s(n).
c
c     n       Number of points
c     s       Parametric variable
c
c-----------------------------------------------------------------------

      Implicit none

c.....Argument variables.

      Integer n

      Real s (n)

c.....Local variables.

      Integer i

c-----------------------------------------------------------------------

c.....Perform a uniform distribution.

      do  i = 2, n-1
        s(i) = s(1) + ( i - 1. ) * ( s(n) - s(1) ) / ( n - 1. )
      enddo

c-----------------------------------------------------------------------

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   The subroutine sptanh is a front end for the routine vinkiter and  c
c   vink                                                               c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sptanh ( lmax, s, dsae, dsbe, kase, nlast )

c      implicit double precision (a-h,o-z)

      dimension s(lmax)
      integer mxit

      smin = s(1)
      smax = s(lmax)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c.....VINKITER
c.....Vinokur's function creates a distribution of grid points which 
c.....satisfy a specified derivative condition, but we require a delta-s
c.....constraint instead. These two values are equivalent only to first 
c.....order, and hence, we resort to an iterative procedure to obtain  
c.....more accurate delta-s's. Up to ten iterative sweeps are
c.....made. The first guess sets ds/dxi = delta-s. The next guess 
c.....recalculates ds/dxi using the leading term in the truncation error
c.....(d2s/d(xi)**2). The next eight iterations use a 2-d secant 
c.....algorithm to home in on the ds/dxi's at both ends which will give 
c.....the correct delta-s.  In the cases where a single-sided 
c.....stretching function is required, (kase = 1 or 2) a secant 
c.....algorithm in 1-d is applied instead. 
c.....Inputs:
c.....   kase = 0 -> delta-s specified on both ends so use a 2D secant 
c.....               method to arrive at the values of dsa and dsb 
c.....               that will satisfy ds1e and ds2e within roundoff.
c.....   kase = 1 -> delta-s specified at the last endpoint only, so 
c.....               use a 1-D secant method to arrive at the value 
c.....               of dsb that will satisfy dsbe within roundoff.
c.....   kase = 2 -> delta-s specified at the first endpoint only, so 
c.....               use a 1-D secant method to arrive at the value 
c.....               of dsa which will satisfy ds1e within roundoff.
c
c      subroutine vinkiter( s, lmax, smin, smax, dsae, dsbe, 
c     &                     kase, nlast )
c
c#ifdef DP
c      implicit double precision (a-h,o-z)
c#endif

c      dimension s(*)
c      integer mxit

c#ifdef DP
c      tolmin2 = 1.0d-14
c#else
      tolmin2 = 1.0e-07
c#endif
      tolmin  = tolmin2 * tolmin2

      mxit = 20

      if ( kase .eq. 0 ) then
  
c........Initial guess, an = dsae, bn = dsbe.

         an2 = dsae 
         bn2 = dsbe  
         call vink(s,lmax,smin,smax,an2,bn2,esa,esb,kase )
         fn2 = esa / dsae - 1
         gn2 = esb / dsbe - 1

c........Second guess, calculate ds1 and ds2 from a truncated 
c........Taylor Series.

         dssa = 2.*s(   1)-5.*s(     2)+4.*s(     3)   -s(     4)
         dssb = 2.*s(lmax)-5.*s(lmax-1)+4.*s(lmax-2)   -s(lmax-3)
         an1 = dsae-0.5*dssa
         bn1 = dsbe+0.5*dssb
         call vink(s,lmax,smin,smax,an1,bn1,esa,esb,kase )
         fn1 = esa/dsae - 1
         gn1 = esb/dsbe - 1
         an  = an1
         bn  = bn1

c........3rd thru nth guesses, use 2-d secant method.

         do 10 n = 3, mxit 
c...........Calculate offset derivatives.
            call vink(s,lmax,smin,smax,an2,bn1,esa21,esb21,kase )
            call vink(s,lmax,smin,smax,an1,bn2,esa12,esb12,kase )
            fa = ( esa - esa21 )/dsae
            fb = ( esa - esa12 )/dsae
            ga = ( esb - esb21 )/dsbe
            gb = ( esb - esb12 )/dsbe
            den = fa*gb - fb*ga
            dela = -(an1 - an2)
            delb = -(bn1 - bn2)
c...........Stick with last guess if approaching roundoff.
            if ( abs(den) .lt. tolmin)then
               call vink(s,lmax,smin,smax,an,bn,esa,esb,kase )
               return 
            end if
c...........Calculate next distribution.
            an = an1 + dela*(  gb*fn1 - fb*gn1 )/den
            bn = bn1 + delb*(- ga*fn1 + fa*gn1 )/den
            call vink(s,lmax,smin,smax,an,bn,esa,esb,kase )
            fn = esa/dsae - 1
            gn = esb/dsbe - 1
c...........Update n, n-1, n-2 and continue.
            an2 = an1
            bn2 = bn1
            an1 = an
            bn1 = bn
            fn1 = fn
            gn1 = gn
   10    continue

      else if ( kase .eq. 1 ) then
  
c........Initial guess, bn = dsbe.

         bn2 = dsbe 
         call vink(s,lmax,smin,smax,dsae,bn2,esa,esb,kase )
         fn2 = esb/dsbe - 1

c........Second guess, calculate ds1 and ds2 from a truncated 
c........Taylor Series.

         dssb = 2.*s(lmax)-5.*s(lmax-1)+4.*s(lmax-2)   -s(lmax-3)
         bn1 = dsbe-0.5*dssb
         call vink(s,lmax,smin,smax,dsae,bn1,esa,esb,kase )
         fn1 = esb/dsbe - 1
         bn = bn1

c........3rd thru nth guesses, use 1-d secant method.

         do 20 n = 3, mxit 
c...........Stick with last guess if approaching roundoff.
            den = fn1-fn2
            if ( abs(den) .lt. tolmin2)then
               call vink(s,lmax,smin,smax,dsae,bn,esa,esb,kase )
               return
            end if
c...........Calculate next distribution.
            bn = bn1 - fn1*(bn1-bn2)/den 
            call vink(s,lmax,smin,smax,dsae,bn,esa,esb,kase )
            fn = esb/dsbe - 1
c...........Update n, n-1, n-2 and continue.
            bn2 = bn1
            bn1 = bn
            fn2 = fn1
            fn1 = fn
   20    continue

      else if ( kase .eq. 2 ) then
  
c........Initial guess, an = dsae.

         an2 = dsae 
         call vink(s,lmax,smin,smax,an2,dsbe,esa,esb,kase )
         fn2 = esa/dsae - 1

c........Second guess, calculate ds1 and ds2 from a truncated 
c........Taylor Series.

         dssa = 2.*s(   1)-5.*s(     2)+4.*s(     3)   -s(     4)
         an1 = dsae-0.5*dssa
         call vink(s,lmax,smin,smax,an1,dsbe,esa,esb,kase )
         fn1 = esa/dsae - 1
         an = an1

c........3rd thru nth guesses, use 1-d secant method.

         do 30 n = 3, mxit 
c...........Stick with last guess if approaching roundoff.
            den = fn1-fn2
            if(abs(den) .lt. tolmin2)then
               call vink(s,lmax,smin,smax,an,dsbe,esa,esb,kase )
               return
            end if
c...........Calculate next distribution.
            an = an1 - fn1*(an1-an2)/den 
            call vink(s,lmax,smin,smax,an,dsbe,esa,esb,kase )
            fn = esa/dsae - 1
c...........Update n, n-1, n-2 and continue.
            an2 = an1
            an1 = an
            fn2 = fn1
            fn1 = fn
   30    continue

      end if

      return 
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     
c.....VINK
c.....Distributes points on a curve so that specified derivatives
c.....at the edges are satisfied (from NASA CR 3313 by Vinokur, 1980).
c.....Vinokur's algorithm is designed to distribute a specified 
c.....number of points along a curve, given the number of points, 
c.....the length of the curve, and the deriviative conditions at 
c.....both ends of the curve. In CFD applications, the user would
c.....usually rather specify the delta-s's at the ends of the curve, 
c.....which are equivalent to the derivatives only to first order. 
c.....Therefore, the user may wish to apply this algorithm iteratively
c.....to obtain an exact delta-s specification. Subroutine 
c.....vinkiter will iterate on this scheme until the proper delta-s 
c.....constraints are satisfied. 
c.....Inputs:
c.....   lmax -> number of points on the curve
c.....   smin, smax -> beginning and end values of s
c.....   ds1,  ds2  -> the derivative end conditions input into 
c.....                 Vinokur's function
c.....   kase = 0 -> satisfy delta-s on both ends
c.....        = 1 -> satisfy delta-s only at xi=ximax     
c.....        = 2 -> satisfy delta-s only at xi=ximin
c.....Outputs:
c.....   s(xi) -> resulting s distribution from Vinokur's function
c.....   es1 -> ( s(ximin+1)-s(ximin) ) <-  calculated delta-s
c.....   es2 -> ( s(ximax)-s(ximax-1) ) <-
c.....Note:
c.....   Additionally, this version uses the approximate inverse 
c.....   solution for y=sin(x)/x and y=sinh(x)/x rather than a Newton 
c.....   iteration. The approximate solution was also taken from 
c.....   NASA CR 3313.

      subroutine vink( s, lmax, smin, smax, ds1, ds2, es1, es2, kase )

c#ifdef DP
c      implicit double precision (a-h,o-z)
c#endif

      dimension s(*), d1(4,2), d2(4,2)

      pi = 4.0 * atan(1.0)

c.....Initialization.

      sdel = smax - smin
      s0 = sdel / real( lmax-1 ) / ds1
      s1 = sdel / real( lmax-1 ) / ds2
      b = sqrt( s0 * s1 )
      a = sqrt( s0 / s1 )
      if ( kase .eq. 1 ) then
         b = s1
      else if ( kase .eq. 2 ) then
         b = s0
      endif

c.....Calculate x based on value of B.

      if ( b .lt. 1.0 ) then

c........x is real.

         if ( b .lt. 0.26938972 ) then
            x = ( 1.0 
     &          - b 
     &          + b**2 
     &          - ( 1.0 + PI ** 2.0 / 6.0 ) * b**3
     &          +  6.794732 * b**4 
     &          - 13.205501 * b**5  
     &          + 11.726095 * b**6 ) 
     &        * PI
         else
            c = 1.0 - b
            x = ( 1.0
     &          + 0.15 * c 
     &          + 0.057321429 * c**2 
     &          + 0.048774238 * c**3
     &          - 0.053337753 * c**4 
     &          + 0.075845134 * c**5 )
     &        * sqrt( 6.0 * c )
         endif

      else if ( b .eq. 1.0 ) then

c........x is zero.

         x = 0.0

      else

c........x is imaginary.

         if ( b .lt. 2.7829681 ) then
            c = b - 1.0
            x = ( 1.0
     &          - 0.15 * c 
     &          + 0.0573214290 * c**2 
     &          - 0.0249072950 * c**3
     &          + 0.0077424461 * c**4 
     &          - 0.0010794123 * c**5 )
     &        * sqrt( 6.0 * c )
         else
            v = log(b)
            w = 1.0 / b - 0.028527431
            x = v + ( 1.0 + 1.0 / v ) * log( 2.0*v ) - 0.02041793
     &        + 0.24902722 * w 
     &        + 1.9496443 * w**2
     &        - 2.6294547 * w**3 
     &        + 8.56795911 * w**4
         endif

      endif

c.....Distribute points along edge.

      if ( kase .eq. 1 .or. kase .eq. 2 ) then
         s(1   ) = 0.0
         s(lmax) = sdel
         do 9 i = 2, lmax-1
            j = lmax + 1 - i
            xi = real(i-1) / (lmax-1)
            if ( b .gt. 1.0001 ) then
               u1 = 1.0 + tanh( x/2.0 * ( xi - 1.0 ) )
     &                  / tanh( x/2.0 )
            else if ( b .lt. 0.9999 ) then
               u1 = 1.0 + tan ( x/2.0 * ( xi - 1.0 ) )
     &                  / tan ( x/2.0 )
            else
               u1 = xi * ( 1.0 - 0.5 * ( b - 1.0 ) * ( 1.0 - xi ) 
     &            * ( 2.0 - xi ) )
            endif
            u2 = sinh( xi*x ) / sinh( x )
            if ( kase .eq. 1 ) then
               fact = abs( ds1 )
               s(j) = ( ( 1.0 - fact ) * ( 1.0 - u1 ) 
     &              + fact * ( 1.0 - u2 ) ) * sdel
            else if ( kase .eq. 2 ) then
               fact = abs( ds2 )
               s(i) = ( ( 1.0 - fact ) * u1 + fact * u2 ) * sdel
            endif
    9    continue
      else
         do 5 i = 1, lmax
            xi = real(i-1) / real(lmax-1)
            cnum = x * ( xi-0.5 )
            cden = x / 2.0
            if ( b .lt. 0.9999 ) then
               cc = tan(cnum) / tan(cden)
               u = 0.5 * ( 1.0 + cc )
            else if ( b .ge. 0.9999 .and. b .le. 1.0001 ) then
               u = xi * ( 1.0 + 2.0 * ( b - 1.0 ) * ( xi - 0.5 ) 
     &           * ( 1.0 - xi ) )
            else if ( b .gt. 1.0001 ) then
               cc = tanh(cnum) / tanh(cden)
               u = 0.5 * ( 1.0 + cc )
            endif
            s(i) = u * sdel / ( a + ( 1.0 - a ) * u )
    5    continue
      endif

      do 8 l = 1, lmax
         s(l) = s(l) + smin
    8 continue
      es1 = s(   2) - s(     1)
      es2 = s(lmax) - s(lmax-1)

      return
      end


