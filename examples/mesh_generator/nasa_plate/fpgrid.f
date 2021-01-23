      Program fpgrid

c     Creates a grid for a laminar boundary layer on a flat plate.

      parameter ( nimax = 200 )
      parameter ( njmax = 100 )

      dimension x (nimax,njmax)
      dimension y (nimax,njmax)

      dimension eta (njmax)

c.....Get grid parameters.

      fplen = 1.0
      xin   = - 0.25

      write(*,*) ' '
      write(*,*) 'Enter axial grid spacing (ft)' !0.025
      read(*,*) spax

      write(*,*) ' '
      write(*,*) 'Enter wall spacing (in terms of eta)' !0.4
      read (*,*) spwall

c.....Create an evenly-spaced axial grid.

      nia = 0.25 / spax + 1 
      nib = 4 * nia - 3
      ni  = nia + nib - 1 

      write(*,*) ' '
      write(*,*) 'Number of axial grid points: ', ni
      write(*,*) 'Leading edge is at i = ', nia

      x(1,1) = xin
      do  i = 2, nia
        x(i,1) = x(i-1,1) + spax
      enddo
      x(nia,1) = 0.0

      do  i = 2, nib
        ix = nia + i - 1
        x(ix,1) = x(ix-1,1) + spax
      enddo
      x(ni,1) = fplen

      do  i = 1, ni
        write(*,*) i, x(i,1)
      enddo

c.....Create the eta distribution.

      eta(1) = 0.0
      ibl    = 0

      do  j = 2, njmax
        if ( ibl .eq. 0 ) then
          eta(j) = eta(j-1) + spwall
          if ( eta(j) .ge. 4.0 )  ibl = 1 
        else
          eta(j) = eta(j-1) + 1.1 * ( eta(j-1) - eta(j-2) )
        endif
        if ( eta(j) .gt. 50 )  go to 10
      enddo

      write(*,*) ' '
      write(*,*) 'Warning: njmax reached'

   10 continue

      nj = j

      write(*,*) ' '
      write(*,*) 'Number of wall normal grid points: ', nj

      do  j = 1, nj
        write(*,*) j, eta(j)
      enddo

c.....Create the x-coordinates for entire grid.

      do  i = 1, ni
        do  j = 2, nj
          x(i,j) = x(i,1)
        enddo
      enddo

c.....Create the y grid coordinates starting at about x = 0.25.

      istart = ni * 2 / 5 + 1 

      ufs    = 129.6974
      dynvis = 6.4849819E-04

      do  i = istart, ni
        do  j = 1, nj
          y(i,j) = eta(j) / sqrt( 0.5 * ufs / ( dynvis * x(i,j) ) )
        enddo
      enddo

c.....For i = 1, istart-1, use the y-coordinates at istart.

      do  i = 1, istart
        do  j = 1, nj
          y(i,j) = y(istart,j)
        enddo
      enddo

c.....Write out the grid to a Plot3d file.

c.....open ( unit=7, file='fplam.d', form='formatted', status='unknown' )
      open ( unit=7, file='fplam.txt', form='formatted' )
      write(7, *) ni, nj
      write(7, *) (( x(i,j), i=1,ni), j=1,nj),
     &         (( y(i,j), i=1,ni), j=1,nj)

      stop
      end
