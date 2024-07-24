      integer function ndex(x0, x, n)
      integer n, ndexd
      real*8 x0, x(n)
      ndex = ndexd(x0, x, n)
      return
      end

      function ndexd(x0, x, n)
      implicit real*8(a - h,o - z)
      dimension x(n)

c This function finds the value of ndexd such that x(ndexd) <=
c x0 < x(ndexd+1). If x0 < x(1) or x0 >= x(n), 0 and n are returned,
c respectively. If x(i) is decreasing instead of increasing, it finds
c i such that x(i) >=  x0 > x(i+1).
c The search time scales roughly as log2(n).

      ndexd = 1
      if (n .le. 1) return
      ndexd = 0

      if (x(2) .lt. x(1)) go to 30

      if (x0 .lt. x(1)) return

      if (x0 .gt. x(n)) then
        ndexd = n
        return
      end if






c x(i) increasing.

      ilow = 1
      ihigh = n

      do 10 ntries = 1, n
        ihalf = (ihigh + ilow) / 2

        if (x0 .ge. x(ihalf)) then
          ilow = ihalf
        else
          ihigh = ihalf
        end if

        if ((ihigh - ilow) .eq. 1) then
          ndexd = ilow
          return
        end if

  10  continue

  20  format(' routine NDEXD unable to to find correct range index.')
      write(7,20)
      write(7,*) ' n: ', n
      write(7,'(a4,e12.5,a4)') ' x0:', x0, '  x:'
      write(7,'(e12.5)') (x(i), i = 1, n)
      stop


c -------------
  30  continue
c -------------

      if (x0 .gt. x(1)) return

      if (x0 .le. x(n)) then
        ndexd = n
        return
      end if







c x(i) decreasing.
      ilow = 1
      ihigh = n

      do 40 ntries = 1, n
        ihalf = (ihigh + ilow) / 2

        if (x0 .gt. x(ihalf)) then
          ihigh = ihalf
        else
          ilow = ihalf
        end if

        if ((ihigh - ilow) .eq. 1) then
          ndexd = ilow
          return
        end if

  40  continue

      write(6,20)
      write(7,*) ' n: ', n
      write(7,'(a4,e12.5,a4)') ' x0:', x0, '  x:'
      write(7,'(e12.5)') (x(i), i = 1, n)

      end

      function ndexs(x0, x, n)
      implicit real*4 (a - h,o - z)
      dimension x(n)

c This function finds the value of ndexs such that x(ndexs) <=
c x0 < x(ndexs+1). If x0 < x(1) or x0 >= x(n), 0 and n are returned,
c respectively. If x(i) is decreasing instead of increasing, it finds
c i such that x(i) >=  x0 > x(i+1).
c The search time scales roughly as log2(n).

      ndexs = 1
      if (n .le. 1) return
      ndexs = 0

      if (x(2) .lt. x(1)) go to 30

      if (x0 .lt. x(1)) return

      if (x0 .gt. x(n)) then
        ndexs = n
        return
      end if






c x(i) increasing.

      ilow = 1
      ihigh = n

      do 10 ntries = 1, n
        ihalf = (ihigh + ilow) / 2

        if (x0 .ge. x(ihalf)) then
          ilow = ihalf
        else
          ihigh = ihalf
        end if

        if ((ihigh - ilow) .eq. 1) then
          ndexs = ilow
          return
        end if

  10  continue

  20  format(' routine NDEXS unable to to find correct range index.')
      write(7,20)
      write(7,*) ' n: ', n
      write(7,'(a4,e12.5,a4)') ' x0:', x0, '  x:'
      write(7,'(e12.5)') (x(i), i = 1, n)
      stop


c -------------
  30  continue
c -------------

      if (x0 .gt. x(1)) return

      if (x0 .le. x(n)) then
        ndexs = n
        return
      end if







c x(i) decreasing.
      ilow = 1
      ihigh = n

      do 40 ntries = 1, n
        ihalf = (ihigh + ilow) / 2

        if (x0 .gt. x(ihalf)) then
          ihigh = ihalf
        else
          ilow = ihalf
        end if

        if ((ihigh - ilow) .eq. 1) then
          ndexs = ilow
          return
        end if

  40  continue

      write(6,20)
      write(7,*) ' n: ', n
      write(7,'(a4,e12.5,a4)') ' x0:', x0, '  x:'
      write(7,'(e12.5)') (x(i), i = 1, n)

      end

      function ndexi(x0, x, n)
      integer n, x0, x(n)

c This function finds the value of ndexi such that x(ndexi) <=
c x0 < x(ndexi+1). If x0 < x(1) or x0 >= x(n), 0 and n are returned,
c respectively. If x(i) is decreasing instead of increasing, it finds
c i such that x(i) >=  x0 > x(i+1).
c The search time scales roughly as log2(n).

      ndexi = 1
      if (n .le. 1) return
      ndexi = 0

      if (x(2) .lt. x(1)) go to 30

      if (x0 .lt. x(1)) return

      if (x0 .gt. x(n)) then
        ndexi = n
        return
      end if






c x(i) increasing.

      ilow = 1
      ihigh = n

      do 10 ntries = 1, n
        ihalf = (ihigh + ilow) / 2

        if (x0 .ge. x(ihalf)) then
          ilow = ihalf
        else
          ihigh = ihalf
        end if

        if ((ihigh - ilow) .eq. 1) then
          ndexi = ilow
          return
        end if

  10  continue

  20  format(' routine NDEXI unable to to find correct range index.')
      write(7,20)
      write(7,*) ' n: ', n
      write(7,'(a4,i10,a4)') ' x0:', x0, '  x:'
      write(7,'(i10)') (x(i), i = 1, n)
      stop


c -------------
  30  continue
c -------------

      if (x0 .gt. x(1)) return

      if (x0 .le. x(n)) then
        ndexi = n
        return
      end if







c x(i) decreasing.
      ilow = 1
      ihigh = n

      do 40 ntries = 1, n
        ihalf = (ihigh + ilow) / 2

        if (x0 .gt. x(ihalf)) then
          ihigh = ihalf
        else
          ilow = ihalf
        end if

        if ((ihigh - ilow) .eq. 1) then
          ndexi = ilow
          return
        end if

  40  continue

      write(6,20)
      write(7,*) ' n: ', n
      write(7,'(a4,i10,a4)') ' x0:', x0, '  x:'
      write(7,'(i10)') (x(i), i = 1, n)

      end
