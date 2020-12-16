      subroutine sort(a,b,nn)
c  sorts the array a increasingly and makes the same rearrange-
c  ments on b, so that a dependence b=f(a) is kept
c
c   
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 15/7/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      dimension a(nn),b(nn),ip(50)
c  initial setting
      ip(1)=nn+1
      ip(2)=1
      k=2
c
    5 n1=ip(k)
      n2=ip(k-1)-1
      k=k-1
c  check if n2-n1>1
    7 if(n2-n1-1) 60,50,10
   10 continue
c
c  sorting in elements smaller and greater than the
c  aritmethical mean
c
      sum=0.d0
      do 12 i=n1,n2
   12 sum=sum+a(i)
      amean=sum/dfloat(n2-n1+1)
      ii=n1-1
      jj=n2+1
      go to 16
c
c  seek exchange element from below
   14 ii=ii+1
      if(a(ii)-amean) 14,14,15
   15 if(ii.gt.jj) go to 40
c  exchange ii with jj
      aa=a(ii)
      a(ii)=a(jj)
      a(jj)=aa
c
      bb=b(ii)
      b(ii)=b(jj)
      b(jj)=bb
c
c  seek wrong element from below
   16 ii=ii+1
      if(ii.eq.n2) go to 40
      if(a(ii)-amean) 16,16,17
   17 if(ii.ge.jj-1) go to 40
c
c  seek exchange element from above
c
   18 jj=jj-1
      if(jj.eq.n1) go to 40
      if(amean-a(jj)) 18,19,19
c  exchange jj with ii
   19 if(ii.gt.jj) go to 40
      aa=a(jj)
      a(jj)=a(ii)
      a(ii)=aa
c
      bb=b(jj)
      b(jj)=b(ii)
      b(ii)=bb
c
c  seek wrong element from above
   20 jj=jj-1
      if(amean-a(jj)) 20,21,21
   21 if(ii.eq.jj) go to 40
      go to 14
c
c
   40 if(a(ii).gt.amean) ii=ii-1
c
c  check if all elements in the interval are equal
      if(ii.eq.n2) go to 60
      if(ii.lt.n1) go to 60
c  now a(n1-n2) is sorted so that elements .le. amean is before
c  or at ii, the rest after
c
      k=k+1
      ip(k)=ii+1
      n2=ii
      go to 7
c
c  the case of two elements
   50 if(a(n1).le.a(n2)) go to 60
c  exchange n1 and n2
      aa=a(n1)
      a(n1)=a(n2)
      a(n2)=aa
c
      bb=b(n1)
      b(n1)=b(n2)
      b(n2)=bb
c
c
   60 if(n2.eq.nn) return
      go to 5
      end
