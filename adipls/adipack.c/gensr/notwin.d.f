      logical function notwin(xw1,xw2,x) 
c
c  windowing function
c  returns true, if x is not in window defined by xw1 and xw2
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      if(xw1.gt.xw2) then
        notwin = .false.
      else
        notwin = x.lt.xw1.or.x.gt.xw2
      end if
      return
      end
