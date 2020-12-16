      integer function intgpt(x)
      double precision x
c  finds integer part of x, i.e. largest integer  .le. x
c  (note: differs from standard fortran function int for negative x)
c
      ix=int(x)
      if(x.lt.0.and.ix.ne.x) ix=ix-1
c
      intgpt=ix
      return
      end
