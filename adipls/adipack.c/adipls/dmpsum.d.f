      subroutine dmpsum(ids,iw)
c
c  dumps complete summary on d/s ids onto d/s iw
c
c  modified 13/8/87, to include iw as argument.
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      character a*132
c
      rewind ids
      write(iw,100)
   10 read(ids,110,end=50,err=50) a
      write(iw,110) a
      go to 10
   50 continue
      return
  100 format(1h1)
  110 format(a)
      end
c
c
