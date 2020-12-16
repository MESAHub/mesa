c
      subroutine cmparr(ia,ib,n,aeq,is)
c  sets logical aeq to .true. if ia(i) = ib(i) for i = 1,n,
c  otherwise to false. if is = 1 ia is stored in ib after comparison
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical aeq,store
      dimension ia(n),ib(n)
c
      aeq=.true.
      store=is.eq.1
      do 10 i=1,n
      if(ia(i).ne.ib(i)) aeq=.false.
      if(store) ib(i)=ia(i)
   10 continue
      return
      end
c
c
