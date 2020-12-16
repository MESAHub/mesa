      subroutine rseta4(x,aa,nn,data,iaa,iresa4)   
c  
c  For iresa4 = 1:
c  Reset A4 near centre, to correct for problems near
c  end of hydrogen burning
c  Resetting is only applied if A4/x**2 is non-monotonic
c  near centre
c
c  For iresa4 ge 10:
c  Reset A4 from A2 and derivative of density
c
c  Original version: 25/7/92
c  
c  Modified 8/6/03, suppressing resetting in convective core. In this
c  case, also reset data(6) to ensure consistency.
c
c  Modified 11/9/12 to include option of iresa4 .ge. 10
c
      implicit double precision (a-h, o-z)
      include '../adipr.incl'
      parameter(iw = 5)
      dimension x(*),aa(iaa,*),data(*)
      dimension w(iw,nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) write(istdpr,'(/'' Entering rseta4''/)')
c
      if(iresa4.ge.10) go to 30
c
c  old resetting, for iresa4 = 1
c
c  reference values
c
      axc=data(6)-data(5)
      nref=5
      axref=aa(4,nref)/x(nref)**2
c
      ireset=0
      do 10 n=2,nref-1
      ax=aa(4,n)/x(n)**2
      if((axc-ax)*(ax-axref).lt.0) ireset=1
   10 continue
c
      if(ireset.eq.1) then
c
c  suppress resetting in case of convective core
c
	if(axref.lt.0) then
	  write(istdou,105) 
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,105)
	  datan6=data(5)+axref
	  write(istdou,107) data(6), datan6
	  if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *      write(istdpr,107) data(6), datan6
	  data(6)=datan6
	  return
        end if
c
	write(istdpr,110)
	axcoef=(axref-axc)/x(nref)**2
	do 20 n=2,nref-1
	x2=x(n)**2
	ax=aa(4,n)/x2
	axnew=axc+axcoef*x2
	write(istdpr,115) n, x(n), ax, axnew
   20   aa(4,n)=axnew*x2
c
      end if
c  
      return   
c
c  --------------------------------------------------------------
c
c  resetting from numerical derivative of rho
c
   30 continue
c
c  set logarithm of dimensionless rho
c
      do n=1,nn
       w(1,n)=log(aa(1,n)*aa(5,n))
      end do
c
      call derive(x,w(1,1),w(2,1),nn,iw,iw,1,1)
c
c  set new A4
c
      do n=1,nn
        a4_new=-aa(2,n)-x(n)*w(2,n)
        if(n.ge.nn-3.and.a4_new.le.0) then
          write(istdou,'('' **** Warning. At n, x ='',i5,f12.7,
     *      '' A4_new = '',1pe13.5,'' is negative''/
     *      ''      Replace by old A_4 ='',e13.5)') n, x(n),a4_new,
     *      aa(4,n)
          if(istdpr.gt.0.and.istdpr.ne.istdou)
     *      write(istdou,'('' **** Warning. At n, x ='',i5,f12.7,
     *      '' A4_new = '',1pe13.5,'' is negative''/
     *      ''      Replace by old A_4 ='',e13.5)') n, x(n),a4_new,
     *      aa(4,n)
            a4_new=aa(4,n)
        end if
        aa(4,n)=a4_new
c..        write(97,'(0pf12.8,1p4e13.5)') x(n),aa(4,n),a4_new,aa(2,n),
c..     *    w(2,n)
      end do
      return
  105 format(//
     *  ' Resetting of A4 near centre suppressed in convective core')
  107 format(/ 'Reset data(6). Old, new values =',1p2e13.5)
  110 format(//' Reset A4 near centre. ',
     *   ' n, x, old, new values of A4/x**2:'/)
  115 format(i5,1p3e13.5)
      end  
