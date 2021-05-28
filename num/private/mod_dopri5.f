! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
      module mod_dopri5
      use const_def, only: dp
      use math_lib      
      
      contains

      subroutine do_dopri5(
     &         n,fcn,x,y,xend,h,max_step_size,max_steps,
     &         rtol,atol,itol,solout,iout,work,lwork,iwork,liwork,
     &         lrpar,rpar,lipar,ipar,lout,idid)
c *** *** *** *** *** *** *** *** *** *** *** *** ***
c          declarations 
c *** *** *** *** *** *** *** *** *** *** *** *** ***
      implicit real(dp) (a-h,o-z)
      integer, intent(in) :: n ! the dimension of the system
      interface
#include "num_fcn.dek"
#include "num_solout.dek"
      end interface
      real(dp), intent(inout) :: x 
      real(dp), intent(inout), pointer :: y(:) ! (n)
      real(dp), intent(in) :: xend, h, max_step_size
      real(dp), intent(in) :: rtol(*)
      real(dp), intent(in) :: atol(*)
      integer, intent(in) :: itol, iout, liwork, lwork, max_steps
      integer, intent(inout) :: iwork(liwork)
      real(dp), intent(inout) :: work(lwork)
      integer, intent(in) :: lrpar, lipar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      integer, intent(in)  :: lout
      integer, intent(out)  :: idid

      logical arret
c *** *** *** *** *** *** ***
c        setting the parameters 
c *** *** *** *** *** *** ***
      nfcn=0
      nstep=0
      naccpt=0
      nrejct=0
      arret=.false.
c -------- nmax , the maximal number of steps ----- 
      if(max_steps.eq.0)then
         nmax=100000
      else
         nmax=max_steps
         if(nmax.le.0)then
            if (lout.gt.0) write(lout,*)
     &          ' wrong input max_steps=',max_steps
            arret=.true.
         end if
      end if
c -------- meth   coefficients of the method
      if(iwork(2).eq.0)then
         meth=1
      else
         meth=iwork(2)
         if(meth.le.0.or.meth.ge.4)then
            if (lout.gt.0) write(lout,*)
     &          ' curious input iwork(2)=',iwork(2)
            arret=.true.
         end if
      end if  
c -------- nstiff   parameter for stiffness detection  
      nstiff=iwork(4) 
      if (nstiff.eq.0) nstiff=1000
      if (nstiff.lt.0) nstiff=nmax+10
c -------- nrdens   number of dense output components
      nrdens=iwork(5)
      if(nrdens.lt.0.or.nrdens.gt.n)then
         if (lout.gt.0) write(lout,*)
     &           ' curious input iwork(5)=',iwork(5)
         arret=.true.
      else
            if(nrdens.gt.0.and.iout.lt.2)then
               if (lout.gt.0) write(lout,*)
     &      ' warning: put iout=2 for dense output '
            end if 
            if (nrdens.eq.n) then
                do 16 i=1,nrdens
  16            iwork(20+i)=i 
            end if
      end if
c -------- uround   smallest number satisfying 1.d0+uround>1.d0  
      if(work(1).eq.0.d0)then
         uround=2.3d-16
      else
         uround=work(1)
         if(uround.le.1.d-35.or.uround.ge.1.d0)then
            if (lout.gt.0) write(lout,*)
     &        ' which machine do you have? your uround was:',work(1)
            arret=.true.
         end if
      end if
c -------  safety factor -------------
      if(work(2).eq.0.d0)then
         safe=0.9d0
      else
         safe=work(2)
         if(safe.ge.1.d0.or.safe.le.1.d-4)then
            if (lout.gt.0) write(lout,*)
     &          ' curious input for safety factor work(2)=',work(2)
            arret=.true.
         end if
      end if
c -------  fac1,fac2     parameters for step size selection
      if(work(3).eq.0.d0)then
         fac1=0.2d0
      else
         fac1=work(3)
      end if
      if(work(4).eq.0.d0)then
         fac2=10.d0
      else
         fac2=work(4)
      end if
c --------- beta for step control stabilization -----------
      if(work(5).eq.0.d0)then
         beta=0.04d0
      else
         if(work(5).lt.0.d0)then
            beta=0.d0
         else
            beta=work(5)
            if(beta.gt.0.2d0)then
               if (lout.gt.0) write(lout,*)
     &          ' curious input for beta: work(5)=',work(5)
            arret=.true.
         end if
         end if
      end if
c -------- maximal step size
      if(max_step_size.eq.0.d0)then
         hmax=xend-x
      else
         hmax=max_step_size
      end if
c ------- prepare the entry-points for the arrays in work -----
      iey1=21
      iek1=iey1+n
      iek2=iek1+n
      iek3=iek2+n
      iek4=iek3+n
      iek5=iek4+n
      iek6=iek5+n
      ieys=iek6+n
      ieco=ieys+n
c ------ total storage requirement -----------
      istore=ieys+(3+5*nrdens)-1
      if(istore.gt.lwork)then
        if (lout.gt.0) write(lout,*)
     &   ' insufficient storage for work, min. lwork=',istore
        arret=.true.
      end if
      icomp=21
      istore=icomp+nrdens-1
      if(istore.gt.liwork)then
        if (lout.gt.0) write(lout,*)
     &   ' insufficient storage for iwork, min. liwork=',istore
        arret=.true.
      end if
c ------ when a fail has occured, we return with idid=-1
      if (arret) then
         idid=-1
         return
      end if
c -------- call to core integrator ------------
      call dopcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,lout,
     &   solout,iout,idid,nmax,uround,meth,nstiff,safe,beta,fac1,fac2,
     &   work(iey1),work(iek1),work(iek2),work(iek3),work(iek4),
     &   work(iek5),work(iek6),work(ieys),work(ieco),iwork(icomp),
     &   nrdens,lrpar,rpar,lipar,ipar,nfcn,nstep,naccpt,nrejct)
      work(7)=h
      iwork(17)=nfcn
      iwork(18)=nstep
      iwork(19)=naccpt
      iwork(20)=nrejct
c ----------- return -----------
      return
      end subroutine do_dopri5
c
c
c
c  ----- ... and here is the core integrator  ----------
c
      subroutine dopcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,lout,
     &   solout,iout,idid,nmax,uround,meth,nstiff,safe,beta,fac1,fac2,
     &   y1,k1,k2,k3,k4,k5,k6,ysti,rwork,icomp,nrd,lrpar,rpar,lipar,ipar,
     &   nfcn,nstep,naccpt,nrejct)
c ----------------------------------------------------------
c     core integrator for dopri5
c     parameters same as in dopri5 with workspace added 
c ---------------------------------------------------------- 
c         declarations 
c ---------------------------------------------------------- 
      implicit real(dp) (a-h,o-z)

      real(dp) k1(n),k2(n),k3(n),k4(n),k5(n),k6(n)
      dimension y(n),y1(n),ysti(n),atol(*),rtol(*)
      dimension icomp(nrd),iwork(nrd+1)
      logical reject,last 
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      !common /condo5/xold,hout
         
      interface
#include "num_fcn.dek"
#include "num_solout.dek"
      end interface
      
      real(dp), target :: rwork(3+5*nrd)
      real(dp), pointer :: cont(:)
      cont => rwork(3:3+5*nrd)


c *** *** *** *** *** *** ***
c  initialisations
c *** *** *** *** *** *** *** 
      if (meth.eq.1) call cdopri(c2,c3,c4,c5,e1,e3,e4,e5,e6,e7,
     &                    a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,
     &                    a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,
     &                    d1,d3,d4,d5,d6,d7)
      facold=1.d-4  
      expo1=0.2d0-beta*0.75d0
      facc1=1.d0/fac1
      facc2=1.d0/fac2
      posneg=sign(1.d0,xend-x)  
      nonsti=0
c --- initial preparations   
      atoli=atol(1)
      rtoli=rtol(1)    
      last=.false. 
      hlamb=0.d0
      iasti=0
      call fcn(n,x,h,y,k1,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; idid=-5; return; end if
      hmax=abs(hmax)     
      iord=5  
      if (h.eq.0.d0) h=hinit(n,fcn,x,y,xend,posneg,k1,k2,k3,iord,
     &                       hmax,atol,rtol,itol,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; idid=-5; return; end if
      nfcn=nfcn+2
      reject=.false.
      xold=x
      if (iout.ne.0) then 
          irtrn=1
          hout=h
          rwork(1) = xold
          rwork(2) = hout
          iwork(1) = nrd
          iwork(2:nrd+1) = icomp(1:nrd)

         if (iout.ge.2) then 
            do 443 j=1,nrd
            i=icomp(j)
            cont(j)=y(i)
            cont(nrd+j)=0
            cont(2*nrd+j)=0
            cont(3*nrd+j)=0
 443        continue
         end if

          call solout(naccpt+1,xold,x,n,y,rwork,iwork,contd5,lrpar,rpar,lipar,ipar,irtrn)
          if (irtrn.lt.0) goto 79
      else
          irtrn=0
      end if
c --- basic integration step  
   1  continue
      if (nstep.gt.nmax) goto 78
      if (0.1d0*abs(h).le.abs(x)*uround)goto 77
      if ((x+1.01d0*h-xend)*posneg.gt.0.d0) then
         h=xend-x 
         last=.true.
      end if
      nstep=nstep+1
c --- the first 6 stages
      if (irtrn.ge.2) then
         call fcn(n,x,h,y,k1,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      end if
      do 22 i=1,n 
  22  y1(i)=y(i)+h*a21*k1(i)
      call fcn(n,x+c2*h,h,y1,k2,lrpar,rpar,lipar,ipar,ierr)      
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 23 i=1,n 
  23  y1(i)=y(i)+h*(a31*k1(i)+a32*k2(i))
      call fcn(n,x+c3*h,h,y1,k3,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 24 i=1,n 
  24  y1(i)=y(i)+h*(a41*k1(i)+a42*k2(i)+a43*k3(i))
      call fcn(n,x+c4*h,h,y1,k4,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 25 i=1,n 
  25  y1(i)=y(i)+h*(a51*k1(i)+a52*k2(i)+a53*k3(i)+a54*k4(i))      
      call fcn(n,x+c5*h,h,y1,k5,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 26 i=1,n 
  26  ysti(i)=y(i)+h*(a61*k1(i)+a62*k2(i)+a63*k3(i)+a64*k4(i)+a65*k5(i))
      xph=x+h
      call fcn(n,xph,h,ysti,k6,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 27 i=1,n 
  27  y1(i)=y(i)+h*(a71*k1(i)+a73*k3(i)+a74*k4(i)+a75*k5(i)+a76*k6(i))  
      call fcn(n,xph,h,y1,k2,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      if (iout.ge.2) then 
            do 40 j=1,nrd
            i=icomp(j)
            cont(4*nrd+j)=h*(d1*k1(i)+d3*k3(i)+d4*k4(i)+d5*k5(i)
     &                   +d6*k6(i)+d7*k2(i)) 
  40        continue
      end if
      do 28 i=1,n 
  28  k4(i)=(e1*k1(i)+e3*k3(i)+e4*k4(i)+e5*k5(i)+e6*k6(i)+e7*k2(i))*h
      nfcn=nfcn+6 
c --- error estimation  
      err=0.d0
      if (itol.eq.0) then   
        do 41 i=1,n 
        sk=atoli+rtoli*max(abs(y(i)),abs(y1(i)))
  41    err=err+pow2(k4(i)/sk)
      else
        do 42 i=1,n 
        sk=atol(i)+rtol(i)*max(abs(y(i)),abs(y1(i)))
  42    err=err+pow2(k4(i)/sk)
      end if
      err=sqrt(err/n)  
c --- computation of hnew
      fac11=pow(err,expo1)
c --- lund-stabilization
      fac=fac11/pow(facold,beta)
c --- we require  fac1 <= hnew/h <= fac2
      fac=max(facc2,min(facc1,fac/safe))
      hnew=h/fac  
      if(err.le.1.d0)then
c --- step is accepted  
         facold=max(err,1.0d-4)
         naccpt=naccpt+1
c ------- stiffness detection
         if (mod(naccpt,nstiff).eq.0.or.iasti.gt.0) then
            stnum=0.d0
            stden=0.d0
            do 64 i=1,n 
               stnum=stnum+pow2(k2(i)-k6(i))
               stden=stden+pow2(y1(i)-ysti(i))
 64         continue  
            if (stden.gt.0.d0) hlamb=h*sqrt(stnum/stden) 
            if (hlamb.gt.3.25d0) then
               nonsti=0
               iasti=iasti+1  
               if (iasti.eq.15) then
                  if (lout.gt.0) write (lout,*) 
     &               ' the problem seems to become stiff at x = ',x   
                  if (lout.lt.0) goto 76
               end if
            else
               nonsti=nonsti+1  
               if (nonsti.eq.6) iasti=0
            end if
         end if 
         if (iout.ge.2) then 
            do 43 j=1,nrd
            i=icomp(j)
            yd0=y(i)
            ydiff=y1(i)-yd0
            bspl=h*k1(i)-ydiff 
            cont(j)=y(i)
            cont(nrd+j)=ydiff
            cont(2*nrd+j)=bspl
            cont(3*nrd+j)=-h*k2(i)+ydiff-bspl
  43        continue
         end if
         do 44 i=1,n
         k1(i)=k2(i)
  44     y(i)=y1(i)
         xold=x
         x=xph
         if (iout.ne.0) then
            irtrn=1
            hout=h
            rwork(1) = xold
            rwork(2) = hout
            iwork(1) = nrd
            iwork(2:nrd+1) = icomp(1:nrd)
            call solout(naccpt+1,xold,x,n,y,rwork,iwork,contd5,lrpar,rpar,lipar,ipar,irtrn)
            if (irtrn.lt.0) goto 79
         end if 
c ------- normal exit
         if (last) then
            h=hnew
            idid=1
            return
         end if
         if(abs(hnew).gt.hmax)hnew=posneg*hmax  
         if(reject)hnew=posneg*min(abs(hnew),abs(h))
         reject=.false. 
      else  
c --- step is rejected   
         hnew=h/min(facc1,fac11/safe)
         reject=.true.  
         if(naccpt.ge.1)nrejct=nrejct+1   
         last=.false.
      end if
      h=hnew
      goto 1
c --- fail exit
  76  continue
      idid=-4
      return
  77  continue
      if (lout.gt.0) write(lout,979)x   
      if (lout.gt.0) write(lout,*)' step size too small, h=',h
      idid=-3
      return
  78  continue
      if (lout.gt.0) write(lout,979)x   
      if (lout.gt.0) write(lout,*)
     &     ' more than nmax =',nmax,'steps are needed' 
      idid=-2
      return
  79  continue
      !if (lout.gt.0) write(lout,979)x
 979  format(' exit of dopri5 at x=',e18.4) 
      idid=2
      return
      end subroutine dopcor
c
      function hinit(n,fcn,x,y,xend,posneg,f0,f1,y1,iord,
     &                 hmax,atol,rtol,itol,lrpar,rpar,lipar,ipar,ierr)
c ----------------------------------------------------------
c ----  computation of an initial step size guess
c ----------------------------------------------------------
      implicit real(dp) (a-h,o-z)
      dimension y(n),y1(n),f0(n),f1(n),atol(*),rtol(*)
      integer, intent(in) :: lrpar, lipar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      
      interface
#include "num_fcn.dek"
      end interface

c ---- compute a first guess for explicit euler as
c ----   h = 0.01 * norm (y0) / norm (f0)
c ---- the increment for explicit euler is small
c ---- compared to the solution
      dnf=0.0d0
      dny=0.0d0 
      atoli=atol(1)
      rtoli=rtol(1)    
      if (itol.eq.0) then   
        do 10 i=1,n 
        sk=atoli+rtoli*abs(y(i))
        dnf=dnf+pow2(f0(i)/sk)
  10    dny=dny+pow2(y(i)/sk)
      else
        do 11 i=1,n 
        sk=atol(i)+rtol(i)*abs(y(i))
        dnf=dnf+pow2(f0(i)/sk)
  11    dny=dny+pow2(y(i)/sk)
      end if
      if (dnf.le.1.d-10.or.dny.le.1.d-10) then
         h=1.0d-6
      else
         h=sqrt(dny/dnf)*0.01d0 
      end if
      h=min(h,hmax)
      h=sign(h,posneg) 
c ---- perform an explicit euler step
      do 12 i=1,n
  12  y1(i)=y(i)+h*f0(i)
      call fcn(n,x+h,h,y1,f1,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; idid=-5; return; end if
c ---- estimate the second derivative of the solution
      der2=0.0d0
      if (itol.eq.0) then   
        do 15 i=1,n 
        sk=atoli+rtoli*abs(y(i))
  15    der2=der2+pow2((f1(i)-f0(i))/sk)
      else
        do 16 i=1,n 
        sk=atol(i)+rtol(i)*abs(y(i))
  16    der2=der2+pow2((f1(i)-f0(i))/sk)
      end if
      der2=sqrt(der2)/h
c ---- step size is computed such that
c ----  h**iord * max ( norm (f0), norm (der2)) = 0.01
      der12=max(abs(der2),sqrt(dnf))
      if (der12.le.1.d-15) then
         h1=max(1.0d-6,abs(h)*1.0d-3)
      else
         h1=pow(0.01d0/der12,1.d0/iord) 
      end if
      h=min(100*abs(h),h1,hmax)
      hinit=sign(h,posneg)  
      return
      end function hinit
c

      real(dp) function contd5(ii,x,rwork,iwork,ierr)
! ----------------------------------------------------------
!     this function can be used for continuous output in connection
!     with the output-subroutine for dopri5. it provides an
!     approximation to the ii-th component of the solution at x.
! ----------------------------------------------------------
      integer, intent(in) :: ii ! result is interpolated approximation of y(i) at x=s.
      real(dp), intent(in) :: x ! interpolation x value (between xold and x).
      real(dp), intent(inout), target :: rwork(*)
      integer, intent(inout), target :: iwork(*)
      integer, intent(out) :: ierr
      
      real(dp) :: xold, h, s, s1
      integer :: nd, i, j
      real(dp), pointer :: con(:)
      integer, pointer :: icomp(:)
      
      !dimension con(5*nd),icomp(nd)
      !common /condo5/xold,h
   
      nd = iwork(1)
      icomp => iwork(2:nd+1)
      xold = rwork(1)
      h = rwork(2)
      con => rwork(3:2+5*nd)
      
c ----- compute place of ii-th component 
      i=0 
      do 5 j=1,nd 
      if (icomp(j).eq.ii) i=j
   5  continue
      if (i.eq.0) then
         contd5 = 0
         ierr = -1
         return
      end if  
      ierr=0
      theta=(x-xold)/h
      theta1=1.d0-theta
      contd5=con(i)+theta*(con(nd+i)+theta1*(con(2*nd+i)+theta*
     &           (con(3*nd+i)+theta1*con(4*nd+i))))
      end function contd5
c

      subroutine cdopri(c2,c3,c4,c5,e1,e3,e4,e5,e6,e7,
     &                    a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,
     &                    a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,
     &                    d1,d3,d4,d5,d6,d7)
! ----------------------------------------------------------
!     runge-kutta coefficients of dormand and prince (1980)
! ----------------------------------------------------------
      implicit real(dp) (a-h,o-z)
      c2=0.2d0
      c3=0.3d0
      c4=0.8d0
      c5=8.d0/9.d0
      a21=0.2d0
      a31=3.d0/40.d0
      a32=9.d0/40.d0
      a41=44.d0/45.d0
      a42=-56.d0/15.d0
      a43=32.d0/9.d0
      a51=19372.d0/6561.d0
      a52=-25360.d0/2187.d0
      a53=64448.d0/6561.d0
      a54=-212.d0/729.d0
      a61=9017.d0/3168.d0
      a62=-355.d0/33.d0
      a63=46732.d0/5247.d0
      a64=49.d0/176.d0
      a65=-5103.d0/18656.d0
      a71=35.d0/384.d0
      a73=500.d0/1113.d0
      a74=125.d0/192.d0
      a75=-2187.d0/6784.d0
      a76=11.d0/84.d0
      e1=71.d0/57600.d0
      e3=-71.d0/16695.d0
      e4=71.d0/1920.d0
      e5=-17253.d0/339200.d0
      e6=22.d0/525.d0
      e7=-1.d0/40.d0  
c ---- dense output of shampine (1986)
      d1=-12715105075.d0/11282082432.d0
      d3=87487479700.d0/32700410799.d0
      d4=-10690763975.d0/1880347072.d0
      d5=701980252875.d0/199316789632.d0
      d6=-1453857185.d0/822651844.d0
      d7=69997945.d0/29380423.d0
      return
      end subroutine cdopri


      end module mod_dopri5


