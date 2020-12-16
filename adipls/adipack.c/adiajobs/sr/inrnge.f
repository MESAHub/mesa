      integer function inrnge(n,nrange,krange,init)
c
c  sets flag for inclusion of data read in, based on range values
c  set into the array nrange(1-krange).
c
c  a single record is specified by its number. a range of records, n1
c  to n2, is indicated by setting nrange(i) = n1, nrange(i+1) = -n2.
c  thus nrange = 2,4,5,-8,10,.. selects records with n = 2,4,5,6,7,8,10,..
c
c  nrange(1) .lt. 0 is used to indicate that no ranges are set in
c  nrange.
c
c  the ranges specified must increase monotonically. failure of this
c  (i.e. if the last values of nrange are set to 0), indicates
c  the end of nrange.
c
c  returns inrnge = 1 if n is in range.
c  returns inrnge = 0 if n is not in range.
c  returns inrnge = -1 if end is reached of nrange.
c
c  init must be 1 in first call with new nrange, and 0 otherwise.
c
c  .................................................................
c
      dimension nrange(*)
c
      save
c
      data islrd /0/
c
c  test for initialization
c
      if(init.ne.1) go to 20
      islrd=0
      if(nrange(1).lt.0) go to 20
      nr1=nrange(1)
      nr2=nrange(2)
      if(nr2.lt.0) go to 15
c
      islrd=1
      nr2=nr1
      kr=1
      go to 20
c
   15 nr2=-nr2
      if(nr2.lt.nr1) go to 20
      islrd=1
      kr=2
c
c  test for n in range, if islrd = 1
c
   20 if(islrd) 70,60,22
c
   22 continue
      if(n.ge.nr1) go to 25
      inrnge=0
      go to 80
c
   25 if(n.gt.nr2) go to 30
      inrnge=1
      go to 80
c
c  set new range
c
   30 if(kr.ge.krange) go to 50
      kr=kr+1
      nr1=nrange(kr)
c
c  test for new range in nrange
c
      if(nr1.le.nr2) go to 50
c
c  test for single or several points
c
      if(kr.eq.krange) go to 37
      nr2=nrange(kr+1)
      if(nr2.lt.0) go to 40
c
   37 nr2=nr1
      go to 22
c
   40 nr2=-nr2
      if(nr2.lt.nr1) go to 50
      kr=kr+1
      go to 22
c
c  end of nrange has been reached
c
   50 inrnge=-1
      islrd=-1
      go to 80
c
c  range test not applied
c
   60 inrnge=1
      go to 80
c
c  end reached in nrange
c
   70 inrnge=-1
c
   80 return
      end
