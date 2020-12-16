      subroutine adirhs(x,fd,finh,if,n)
c
c  right hand side subroutine for adiabatic oscillations
c  *****************************************************
c
c  modified 13/8/87 to standardize output
c
c  modified 19/2/88 to include cowling approximation, and variable
c  lambda, for the radial case.
c
c  modified 16/12/93, to enable integration of inhomogeneous equation
c  for polytropic case, when 0 lt npol lt 1
c
c  Modified 20/1/95, to include CR formulation for equilibrium
c  model with turbulent pressure, for iturpr = 2
c  Note: certain aspects of this formulation still need checking
c  (e.g. radial equations). Also, the code can be streamlined and
c  made more transparent.
c
c  Modified 2/10/97, to fix various errors in connection with
c  the treatment of turbulent pressure for iturpr = 1 and 2.
c
c  Modified 17/2/98, to allow use also g/(g tilde) when
c  iturpr = 8 (in cases where aa(6,.) is set to g/(g tilde) in
c  input model)
c
c  Modified 31/10/02 to include option for first-order rotational
c  effects (flagged by irotsl = 1) in solution, according to
c  Soufi et al. formalism.
c
c  Modified 9/7/05, changing only name for consistency with evolution code
c
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      include 'adipls.c.d.incl'
      parameter (iaa = 10, iaa1 = 10, iy = 8)
      dimension fd(if,4),finh(*)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8),
     *  aa(iaa,1)
      common/rotdat/ em, irotsl
      common/cdiagn/ idgrhs,idgrh1
      common/csexpc/ y11, y13, y21, y23, y31, y33, y41, y43,
     *  yt11, yt13, yt21, yt23, yt31, yt33, yt41, yt43, 
     *  yt241, yt243, ytt41, ytt43 
      common/sysord/ ii
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data ioprhs /0/
c
      save
c
c  hardwired flag for diagnostics
c
      idiag=0
      idgrhs=0
c
      if(idiag.gt.0.and.ioprhs.eq.0) then
	open(87,file='ttt.rhs.out',status='unknown')
	ioprhs=1
	write(87,*) el, alb, sig
      end if
c
c  test for plane-parallel case
c
      if(idgrhs.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,*) 'Entering rhs at x, n =',x, n,
     *  '  el, el1, els, sig, ii, iplneq =', 
     *  el, el1, els, sig, ii, iplneq
      if(iplneq.eq.1) then
        fd(1,1)=aa(2,n)-el1
        fd(1,2)=1.d0-aa(2,n)/els
        fd(2,1)=ell-aa(4,n)*els
        fd(2,2)=aa(4,n)-el1
	finh(1)=0.d0
	finh(2)=0.d0
        return
c
      end if
c
   10 xi=1.d0/x
      vgt=aa(2,n)
      at=aa(4,n)
      u=aa(5,n)
c
c  set untilded variables for use when iturpr = 1 or iturpr = 8
c
      if(iturpr.eq.1.or.iturpr.eq.8) then
        ggt=aa(10,n)
        vg=ggt*vgt
        a=at-(vg-vgt)
        if(idgrhs.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    'xi, at, u, ggt, vg, a',xi, at, u, ggt, vg, a
      else
        ggt=1.d0
	vg=vgt
	a=at
      end if
c
c  set Vg0, Vg1, A0, A1 for iturpr = 2
c
      if(iturpr.eq.2) then
	vg1=aa(2,n)
	vg0=vg1+aa(6,n)
	a1=aa(4,n)
	a0=a1-aa(6,n)
      end if
c
      if(el.le.1.e-6) go to 20
c
c  nonradial case
c
      etat=els*aa(1,n)
      eta=ggt*etat
      vge=vg/eta
      aet=at*etat
c
      if(iturpr.ne.2) then
        fd(1,1)=(vgt-el1-2.d0)*xi
        fd(1,2)=(1.d0-vge)*xi
        fd(2,1)=(ell-aet)*xi
        fd(2,2)=(at-el1-1.d0)*xi
      else
        fd(1,1)=(vg0-el1-2.d0)*xi
        fd(1,2)=(1.d0-vg1/eta)*xi
        fd(2,1)=(ell-eta*a0)*xi
        fd(2,2)=(a1-el1-1.d0)*xi
      end if
c
      finh(1)=0.d0
      finh(2)=0.d0
c
      if(icow.ne.0) return
c
      if(iturpr.ne.2) then
        fd(1,3)=-vg*xi
        fd(1,4)=0.d0
        fd(2,3)=eta*at*xi
        fd(2,4)=0.d0
        fd(3,1)=0.d0
        fd(3,2)=0.d0
        fd(3,3)=(1.d0-el1)*xi
        fd(3,4)=xi
        fd(4,1)=-alb*u*at*xi
        fd(4,2)=-alb*u*vge*xi
        fd(4,3)=(ell+u*(a-2.d0)+(1.d0-alb)*u*vg)*xi
        fd(4,4)=(2.d0-el1-2.d0*u)*xi
      else
        fd(1,3)=-vg1*xi
        fd(1,4)=0.d0
        fd(2,3)=eta*a1*xi
        fd(2,4)=0.d0
        fd(3,1)=0.d0
        fd(3,2)=0.d0
        fd(3,3)=(1.d0-el1)*xi
        fd(3,4)=xi
        fd(4,1)=-alb*u*a0*xi
        fd(4,2)=-alb*u*(vg1/eta)*xi
        fd(4,3)=(ell+u*(a1-2.d0)+(1.d0-alb)*u*vg1)*xi
        fd(4,4)=(2.d0-el1-2.d0*u)*xi
      end if
c
      finh(5)=0.d0
      finh(6)=0.d0
c
c  test for setting inhomogeneous terms
c
      if(anres.eq.0) then
	finh(3)=0.d0
	finh(4)=0.d0
	finh(7)=0.d0
	finh(8)=0.d0
      else
	t=1.d0-x
	tn=t**anres
	finh(3)=xi*yt41*tn
	finh(4)=(fd(4,4)+anres/t)*yt41*tn
	finh(7)=xi*yt43*tn
	finh(8)=(fd(4,4)+anres/t)*yt43*tn
      end if
c
c  test for including contributions from rotation
c
      if(irotsl.eq.1.and.em.ne.0) then
        sigom2=1.5d0*(aa(10,n)-1.d0)*aa(1,n)
        sigom=sqrt(sigom2)
	sighat=sqrt(sig)+em*sigom
        sigfct=sighat*sighat/sig
	alpha=2.d0*em*sigom/sighat
	ellfct=ell/(ell-alpha)
	h1=ellfct*alpha
	zeta=ellfct*etat/sigfct
c
        fd(1,1)=fd(1,1)+h1*xi
	fd(1,2)=(zeta-vgt)*xi/etat
	fd(2,1)=(ell*sigfct-etat*(at+h1*h1/zeta))*xi
	fd(2,2)=fd(2,2)-h1*xi
        fd(4,3)=(ell+u*(a-2.d0)+(1.d0/ggt-alb)*u*vg)*xi
      end if
c
      if(idiag.gt.0.and.x.gt.0.95) then
        t = 1.d0-x
        write(87,'(f15.10,1p6e15.7)') t,(fd(4,i),i=1,4),finh(4),finh(8)
      end if
      if(idgrh1.eq.1) then
c..        write(68,90010) n, x,((fd(i,j),j=1,4),i=1,4)
c..90010   format(/' n, x =',i5,f13.8,'   f:'/(1p4e13.5))
      end if
      return
c
c  radial case, depending on icow
c
   20 qxt=aa(1,n)
      qx=ggt*qxt
      if(iturpr.ne.2) then
        fd(1,1)=(vgt-2.d0)*xi
        fd(1,2)=-vg*sig/(qx*x*x)
        if(icow.eq.0) then
          fd(2,1)=1.d0-(qxt*at-alb*qx*u)/sig
        else
          fd(2,1)=1.d0-qxt*at/sig
        end if
        fd(2,2)=at*xi
      else
        fd(1,1)=(vg0-2.d0)*xi
        fd(1,2)=-vg1*sig/(qx*x*x)
        if(icow.eq.0) then
          fd(2,1)=1.d0-(qxt*a0-alb*qx*u)/sig
        else
          fd(2,1)=1.d0-qxt*a0/sig
        end if
        fd(2,2)=a1*xi
      end if
c
      finh(1)=0.d0
      finh(2)=0.d0
c
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,20090) 
     *  n,x,((fd(i,j),j=1,2),i=1,2)
20090 format(i5,f12.6,1p4e13.5)
      return
      end
