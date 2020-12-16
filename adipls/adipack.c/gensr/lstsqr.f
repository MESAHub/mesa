      subroutine lstsqr(n,x,y,a,b,rms,sa,sb,idiag)
c
c  fits least square line
c
c    y = a + b*x
c
c  to data (x(i),y(i)),i=1,n, and computes scatter sa and sb
c  in coefficients of line
c
c  ....................................................................
c
      dimension x(n),y(n)
   15 s1=0.
      s2=0.
      s3=0.
      s4=0.
      do 6 i=1,n
      xx=x(i)
      yy=y(i)
      s1=s1+xx*yy
      s2=s2+xx
      s3=s3+yy
    6 s4=s4+xx*xx
      b=(s1*n-s2*s3)/(s4*n-s2*s2)
      a=s3/n-b*s2/n
      if(idiag.ge.1) write(6,100) a,b
      if(idiag.ge.2) write(6,105)
c
      sum=0.0
      s2=s2/n
      s3=0.0
c
      do 8 i=1,n
      r=y(i)-a-b*x(i)
      s3=s3+(x(i)-s2)**2
      sum=sum+r*r
      if(idiag.ge.2) write(6,110) i,x(i),y(i),r
    8 continue
c
      s1=sum/(n-2)
      sb=sqrt(s1/s3)
      sa=sqrt(s1*(1/float(n)+s2*s2/s3))
c..      write(6,*) 's1, s2, s3, sa', s1, s2, s3, sa
      rms=sqrt(sum/n)
      if(idiag.ge.1) write(6,120) rms,sa,sb
      return
  100 format(///' the line is y = a + b*x, where a =',1pe13.5,
     *  ', and b =',e13.5)
  105 format(//' n,x,y(measured),y(measured)-y(line):'/)
  110 format(i5,1p2e13.5,3x,e11.3)
  120 format(//' root mean square residuum =',1pe13.5//
     *  ' error in a =',e13.5/' error in b =',e13.5)
      end
