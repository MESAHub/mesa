      subroutine decicd(icode,iel,iekin)
c
c  decodes icode (for adiabatic pulsation summary) and sets
c  iel and iekin
c
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      common /cicode/ ic(9)
c
      save
c
      data nmax/9/
c
      j=icode
      do 10 n=1,nmax
      if(j.eq.0) go to 5
      j1=j/10
      ic(n)=j-10*j1
      j=j1
      go to 10
c
    5 ic(n)=0
   10 continue
c
      iel=1
      if(ic(3)+ic(4).gt.0) iel=2
      if(ic(1).eq.0) then
        iekin=iel+5
      else
        iekin=iel+6
      end if
      return
      end
