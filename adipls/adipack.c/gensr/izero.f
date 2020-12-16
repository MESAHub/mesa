      subroutine  izero(ia,nn)
c
c
c     sets ia(n)=0, n=1,nn
c
c  Dated 13/8/91
c
      dimension ia(nn)
      do 1 n=1,nn
    1 ia(n)=0
      return
c
      end
