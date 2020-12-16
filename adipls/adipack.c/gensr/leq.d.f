      subroutine leq(a,b,nn,mm,ia,ib,err)   
c   
c     this routine will find the inverse of a matrix by the method of   
c     partial pivoting and gaussian elimination if b is set equal to
c     the identity matrix   
c   
c     nn - dimension of segment of a to be used  
c     mm - number of right hand columns of b to be used  
c     ia - the total number of rows in large array a
c     ib - the total number of rows in large array b
c     the matrix equation is    ax=b
c     err = det(a)  
c     if mm = 0 leq calculates err = det(a) 
c   
c  note on modification on 23/1 1985:   
c   
c  previously, the routine contained the statements 
c   
c      do 14 k=i2,m1,ib 
c      b(k)=b(k)+b(i1)*r
c   14 i1=i1+ib 
c   
c  this caused problems, on some computers, for m = ib = 1. 
c  then m1 = 1 and the loop was skipped when i2 .gt. 1. 
c  this has been changed.   
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
      dimension a(100),b(100)   
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
10000 n=nn  
      m=mm  
      err=0.0d0 
      detsc=1.0d0   
c   
c     treat n .le. 1 separately 
      if(n-1) 280,281,285   
c     no equation   
  280 if(istdpr.gt.0) write(istdpr,295) n
      return
c     n = 1 
  281 err=a(1)  
      if(m.le.0) return 
      ai=1.d0/err   
      m1=ib*m   
      do 282 j=1,m1,ib  
  282 b(j)=b(j)*ai  
      return
c   
c   
c     find maximum element in each row and divide the row by it 
  285 n1=ia*n   
      ia1=ia+1  
      m1=ib*m   
      do 1 i=1,n
      r= abs(a(i))  
      do 2 j=  1,n1,ia  
      ij=j+i-1  
    2 r=max(r, abs(a(ij)))
      if(r)31,30,31 
   30 if(istdpr.gt.0) write(istdpr,298)i 
      return
   31 do 3 j=1,n1,ia
      ij=j+i-1  
    3 a(ij)=a(ij)/r 
      if(m.eq.0) go to 1
      do 4 j=1,m1,ib
      ij=j+i-1  
    4 b(ij)=b(ij)/r 
    1 detsc=detsc*r 
c   
c   
c     find maximum element in the i'th column   
      n2=n-1
      do 5 i=1,n2   
      ialow=(i-1)*ia+i  
      iaup=(i-1)*ia+n   
      r= abs(a(ialow))  
      ialow1=ialow+1
      imax=ialow
      do 6 j=ialow1,iaup
      if(r- abs(a(j)))7,6,6 
    7 imax=j
      r= abs(a(j))  
    6 continue  
      if(imax-ialow)8,8,9   
c     replace the i'th row with the row that has the maximum element in 
c         the respective column and put the i'th row in its place   
    9 im=imax   
   72 if(im-ia)70,70,71 
   71 im=im-ia  
      go to 72  
   70 do 10 j=1,n1,ia   
      jj=i+j-1 
      ji=im+j-1
      r=a(jj)  
      a(jj)=a(ji)  
   10 a(ji)=r  
c     change sign of determinant   
      detsc=-detsc 
c  
      if(m.eq.0) go to 8   
      do 11 j=1,m1,ib  
      jj=i+j-1 
      ji=im+j-1
      r=b(jj)  
      b(jj)=b(ji)  
   11 b(ji)=r  
c     multiply the i'th row by (the negative of each i'th column element   
c       below the diagonal divided by the diagonal element) and add the
c     resulting row to the respective row of the element used  
    8 iaup1=iaup-1 
c  
c  
      do 12 j=ialow,iaup1  
      if(a(ialow))32,33,32 
   33 joy=i
      if(a(ialow1))81,82,81
   82 if(istdpr.gt.0) write(istdpr,299)joy,joy  
      return   
   81 if(istdpr.gt.0) write(istdpr,297)joy,joy  
      do 34 k=1,n1,ia  
      jj=joy+k-1   
      ji=joy+k 
      if(joy+1-n)35,36,36  
   35 if(istdpr.gt.0) write(istdpr,296) 
      return   
   36 r=a(jj)  
      a(jj)=a(ji)  
   34 a(ji)=r  
c     change sign of determinant   
      detsc=-detsc 
c  
      if(m.eq.0) go to 8   
      do 37 k=1,m1,ib  
      jj=joy+k-1   
      ji=joy+k 
      r=b(jj)  
      b(jj)=b(ji)  
   37 b(ji)=r  
      go to 8  
   32 j1=j+1   
      r=-a(j1)/a(ialow)
      i1=ialow 
      do 13 k=j1,n1,ia 
      a(k)=a(k)+a(i1)*r
   13 i1=i1+ia 
c  
c  loop to reset b has been modified, 25/1/1985.   
c  
      if(m.eq.0) go to 12  
      i1=i 
      i2=j-ialow+i+1   
      do 14 k=1,m1,ib  
      b(i2)=b(i2)+b(i1)*r  
      i1=i1+ib 
   14 i2=i2+ib 
   12 continue 
c  
c  
    5 continue 
c  
c  
c     the matrix is now in triangular form 
c     first set err=1.0d0  
      err=1.0d0
c     if(any diagonal element of a is zero x cannot be solved for  
      do 15 i=1,n  
      idiag=(i-1)*ia+i 
      err=err*a(idiag) 
      if(err) 15,16,15 
   16 if(istdpr.gt.0) write(istdpr,299)i,i  
      return   
   15 continue 
c     scale determinant
      err=err*detsc
c  
      if(m.eq.0) return
c     find solution to ax=b
      do 18 k=1,m  
      ka=(n-1)*ia+n
      kb=(k-1)*ib+n
      b(kb)=b(kb)/a(ka)
      do 19 l=1,n2 
      i=n-l
      r=0.0d0  
      imax=i+1 
      do 20 j=imax,n   
      jj=i+n+1-j   
      ja=(jj-1)*ia+i   
      jb=(k-1)*ib+jj   
   20 r=r+a(ja)*b(jb)  
      la=(i-1)*ia+i
      lb=(k-1)*ib+i
   19 b(lb)=(b(lb)-r)/a(la)
   18 continue 
c  
c  
      return   
c  
  295 format(///20h leq called with n =,i4)
  296 format(///48h the row cannot be changed with the row below it  , 
     .40h because it it the last row in the array   /  
     .55h the solution, matrix x, cannot be found by this method  )
  297 format(5h1  a(,i4,1h,,i4,13h) equals zero  / 
     .47h   try switching this row with the row below it   ,   
     .48h and go back to statement number 8 and try again  )   
  298 format(26h1  all the elements in row,i5,20h  are zero therefore, 
     .55h the solution, matrix x, cannot be found by this method ) 
  299 format(50h1  the solution, matrix x, cannot be found by this  ,  
     .57h method because there is a zero array element in the main ,   
     .9h diagonal  / 30x,2ha(,i4,1h,,i4,8h) = zero  )  
      end  
      subroutine leqsd(a,b,n,m,ia,ib,det,isa,isb)
      implicit double precision (a-h,o-z)
      dimension a(ia,1),b(ib,1)
c
      call leq(a,b,n,m,ia,ib,det)
      return
      end
