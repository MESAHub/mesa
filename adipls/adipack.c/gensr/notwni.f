      logical function notwni(kw1,kw2,k) 
c
c  windowing function
c  returns true, if k is not in window defined by kw1 and kw2
c
      if(kw1.gt.kw2) then
        notwni = .false.
      else
        notwni = k.lt.kw1.or.k.gt.kw2
      end if
      return
      end
