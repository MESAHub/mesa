      subroutine setssm(cs,ics,ss,iss,ssmod,irmod)
c
c  sets short summary from grand summary
c
c  original version 12/7/1985
c
c  modified 26/9/86, to test for new model. for new model (or
c  first call) ssmod is set to model record for short summary,
c  and irmod is set to 1. otherwise irmod is returned as 0.
c
c  modified 16/2/89 to ensure that cyclic frequency is set from 
c  corrected eigenfrequency, when Cowling approximation is used.
c  Note that cs(25) contains period based on uncorrected 
c  eigenfrequency.
c  .................................................................
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      dimension cs(*),ics(*),ss(*),iss(*),ssmod(*)
      dimension csmodp(5)
c
      save csmodp
      data csmodp /-1.,-1.,-1.,-1.,-1./
c
      ss(1)=cs(18)
      ss(2)=cs(19)
c
c  set ss(3) to corrected sigma**2, if cowling approximation is used
c
      ss(3)=cs(21)
      if(ss(3).le.0) ss(3)=cs(20)
c
      ss(4)=cs(24)
      ss(5)=cs(27)
      if(ss(5).le.0) ss(5)=sqrt(cs(21)/cs(20))*16.666667/cs(25)
      iss(1)=ics(5)
      iss(2)=ics(6)
c
c  test for reset model
c
      irmod=0
c
      do 20 i=2,5
      if(abs(csmodp(i)/cs(i)-1).gt.1.e-5) irmod=1
   20 continue
c
      if(irmod.eq.1) then
        do 25 i=2,5
        csmodp(i)=cs(i)
   25   ssmod(i+1)=cs(i)
        ssmod(1)=-1
      end if
c
      return
      end
