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
c  Modified 2/7/95 to include xmod in model record, ssmod(2)
c  .................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension cs(*),ics(*),ss(*),iss(*),ssmod(*)
      dimension csmodp(5)
c
      save 
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
      if(ss(5).le.0) ss(5)=16.666666667d0/cs(25)
      iss(1)=ics(5)
      iss(2)=ics(6)
c
c  test for reset model
c
      irmod=0
c
      do 20 i=2,5
      if(cs(i).ne.0) then
        if(abs(csmodp(i)/cs(i)-1).gt.1.e-5) irmod=1
      end if
   20 continue
c
      if(irmod.eq.1) then
        do 25 i=2,5
        csmodp(i)=cs(i)
   25   ssmod(i+1)=cs(i)
        ssmod(1)=-1
	ssmod(2)=cs(1)
      end if
c
      return
      end
