      subroutine anldet(ds1,c,x,y,iy,aa,ia,nn,nfit)
c
c  analyzes matching determinant for slowly and rapidly varying parts
c  of solution
c
c  modified 2/4/1985 to make ds1 double precision, to take
c  into account corresponding change in s/r adipls.
c
c  restored 11/7/1985 to single precision.
c
c  modified 13/8/87 to standardize output
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
c  ....................................................................
c
      implicit double precision (a-h,o-z)
      dimension ds1(4,4),c(4),x(nn),y(iy,nn),aa(ia,nn)
      common/rhsdat/ el,ell,alb,els
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  modify determinant
      fct=alb*aa(5,nfit)/(els*aa(1,nfit))
      cc=(ds1(3,1)-fct*ds1(1,1))/(fct*ds1(1,2)-ds1(3,2))
      cs=(ds1(3,3)-fct*ds1(1,3))/(fct*ds1(1,4)-ds1(3,4))
      do 38 i=1,4
      ds1(i,1)=ds1(i,1)+cc*ds1(i,2)
      ds1(i,3)=ds1(i,3)+cs*ds1(i,4)
      dd=ds1(i,2)
      ds1(i,2)=-ds1(i,3)
   38 ds1(i,3)=dd
c
      do 39 j=1,4
      do 39 k=3,4
   39 ds1(k,j)=ds1(k,j)-fct*ds1(k-2,j)
c  subdeterminants
      detr=ds1(1,1)*ds1(2,2)-ds1(1,2)*ds1(2,1)
      dets=ds1(3,3)*ds1(4,4)-ds1(3,4)*ds1(4,3)
      if(istdpr.gt.0) write(istdpr,315) detr,dets
  315 format(/'  detr= ',1pe13.5,'   dets =',e13.5)
      if(istdpr.gt.0) write(istdpr,140) fct,cc,cs
  140 format(/' determinant modified with fct =',1pe13.5,
     *  '  cc =',e13.5,'  cs =',e13.5)
c  reset coefficients
      c(2)=c(2)-cc*c(1)
      c(4)=c(4)-cs*c(3)
      if(istdpr.gt.0) write(istdpr,320) ((ds1(i,j),j=1,4),i=1,4)
  320 format(//' after modification determinant is'/(1p4e13.5))
c  modify solution
      do 41 n=1,nfit
      do 41 i=1,4
   41 y(i,n)=y(i,n)+cc*y(i+4,n)
c
      nnp1=nn+1
      nfit1=nfit+1
      do 42 n=nfit1,nnp1
      do 42 i=1,4
   42 y(i,n)=y(i,n)+cs*y(i+4,n)
      return
      end
