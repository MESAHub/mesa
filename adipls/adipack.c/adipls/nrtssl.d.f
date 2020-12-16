      subroutine nrtssl(x,y,zk,ap,aq,rhs,bcs,ii,kk,ki,nn,
     *  id,ndpr,idiag)
c
c  tests solution from nrk. currently only test of right hand side
c  subroutine rhs is implemented. also
c  assumes that there are no integral constraints.
c
c  test is printed at ndpr points. if idiag = 1 additional details
c  are printed about solution and equations
c
c  common /work/ is used as work area.
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
c  Modified 26/4/99, fixing initialization of if0 and ifd0
c
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(id,1),zk(1),ap(1),aq(1)
      common/work/ w(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      external rhs,bcs
c
c  storage parameters, etc
c
      ik=ii+kk
      nfd=ik*ii
      nerr=4*ii
      if=1
      ifd=if+ii+ii
      ih=ifd+nfd+nfd
      ihd=ih+1
      ierr=ihd+2
      in0=0
      in1=1
      if0=if
      ifd0=ifd
c
c  step through solution
c
      idpr=0
      if(ndpr.le.1) go to 10
      idpr=max0(1,(nn-1)/(ndpr-1))
      if(istdpr.gt.0) write(istdpr,100)
c
   10 nw=ierr-1
c
      do 40 n=1,nn
c
c  set storage indices at n
c
      if0=if1
      ifd0=ifd1
      if1=if+in1*ii
      ifd1=ifd+in1*nfd
c
c  zero fd
c
      do 12 i=1,nfd
   12 w(ifd1+i-1)=0
c
      call rhs(x(n),y(1,n),zk,ap,aq,w(if1),w(ifd1),w(ih),w(ihd),
     *  ii,n,1)
c
      do 15 i=1,ii
   15 w(nw+i)=w(if1+i-1)
c
      if(n.eq.1) go to 35
c
c  set difference equations
c
      np=n-1
      dx=x(n)-x(np)
      nw1=nw+ii
      nw2=nw1+ii
      nw3=nw2+ii
c
      do 20 i=1,ii
      i1=i-1
      w(nw1+i)=(y(i,n)-y(i,np))/dx
      w(nw2+i)=0.5*(w(if1+i1)+w(if0+i1))
   20 w(nw3+i)=(w(nw1+i)-w(nw2+i))/(abs(w(nw2+i))+1.e-10)
c
c  test for output
c
      if(idpr.eq.0) go to 35
      if(mod(n-1,idpr).eq.0.and.istdpr.gt.0) 
     *  write(istdpr,110) n,x(n),(w(nw3+i),i=1,ii)
c
c  permute storage indices, increment nw
c
   35 i=in0
      in0=in1
      in1=i
   40 nw=nw+nerr
c
      if(idiag.le.0) return
c
c  additional output
c
      call printw(x,y,ii,nn,id,idpr,'y')
      call printw(x,w(ierr),ii,nn,nerr,idpr,'f')
      call printw(x,w(ierr+ii),ii,nn,nerr,idpr,'dy/dx')
      call printw(x,w(ierr+ii+ii),ii,nn,nerr,idpr,'f combin.')
c
      return
c
  100 format(///' test solution from nrk:'/
     *  ' *************************'//
     *  ' n, x, relative difference dy/dx - f:'/)
  110 format(i4,1pe13.5,3x,8e12.4/
     *  (20x,8e12.4))
      end
      subroutine printw(x,y,ii,nn,id,idpr,title)
c
c  prints array y(i,n),i=1,ii, n=1,nn against x(n), with step
c  idpr in n
c
      implicit double precision (a-h,o-z)
      character*(*) title
      dimension x(1),y(id,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(idpr.le.0) return
c
      if(istdpr.gt.0) then
        write(istdpr,100) title
        do 10 n=1,nn,idpr
   10   write(istdpr,110) n,x(n),(y(i,n),i=1,ii)
      end if
      return
  100 format(///' n, x, ',a8/)
  110 format(i4,1pe13.5,3x,8e12.4/
     *  (20x,8e12.4))
      end
