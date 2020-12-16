! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      module mod_dop853
      use const_def, only: dp
      use math_lib
      
      contains

      subroutine do_dop853(
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
      integer, intent(inout) :: iwork(liwork) ! should be 0 on entry
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
     &       ' warning: put iout=2 for dense output '
         end if 
         if (nrdens.eq.n) then
            do i=1,nrdens
               iwork(i+20)=i
            end do
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
         fac1=0.333d0
      else
         fac1=work(3)
      end if
      if(work(4).eq.0.d0)then
         fac2=6.d0
      else
         fac2=work(4)
      end if
c --------- beta for step control stabilization -----------
      if(work(5).eq.0.d0)then
         beta=0.0d0
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
      iek1=21
      iek2=iek1+n
      iek3=iek2+n
      iek4=iek3+n
      iek5=iek4+n
      iek6=iek5+n
      iek7=iek6+n
      iek8=iek7+n
      iek9=iek8+n
      iek10=iek9+n
      iey1=iek10+n
      ieco=iey1+n
c ------ total storage requirement -----------
      istore=ieco+(3+8*nrdens)-1
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
c -------- when a fail has occured, we return with idid=-1
      if (arret) then
         idid=-1
         return
      end if
c -------- call to core integrator ------------
      call dp86co(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,lout,
     &   solout,iout,idid,nmax,uround,meth,nstiff,safe,beta,fac1,fac2,
     &   work(iek1),work(iek2),work(iek3),work(iek4),work(iek5),
     &   work(iek6),work(iek7),work(iek8),work(iek9),work(iek10),
     &   work(iey1),work(ieco),iwork(icomp),nrdens,lrpar,rpar,lipar,ipar,
     &   nfcn,nstep,naccpt,nrejct)
      work(7)=h
      iwork(17)=nfcn
      iwork(18)=nstep
      iwork(19)=naccpt
      iwork(20)=nrejct
c ----------- return -----------
      return
      end subroutine do_dop853
c
c
c
c  ----- ... and here is the core integrator  ----------
c
      subroutine dp86co(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,lout,
     &   solout,iout,idid,nmax,uround,meth,nstiff,safe,beta,fac1,fac2,
     &   k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,y1,rwork,icomp,nrd,lrpar,rpar,lipar,ipar,
     &   nfcn,nstep,naccpt,nrejct)
c ----------------------------------------------------------
c     core integrator for dop853
c     parameters same as in dop853 with workspace added 
c ---------------------------------------------------------- 
c         declarations 
c ---------------------------------------------------------- 
      implicit real(dp) (a-h,o-z)
      parameter (
     &  c2  = 0.526001519587677318785587544488d-01,
     &  c3  = 0.789002279381515978178381316732d-01,
     &  c4  = 0.118350341907227396726757197510d+00,
     &  c5  = 0.281649658092772603273242802490d+00,
     &  c6  = 0.333333333333333333333333333333d+00,
     &  c7  = 0.25d+00,
     &  c8  = 0.307692307692307692307692307692d+00,
     &  c9  = 0.651282051282051282051282051282d+00,
     &  c10 = 0.6d+00,
     &  c11 = 0.857142857142857142857142857142d+00,
     &  c14 = 0.1d+00,
     &  c15 = 0.2d+00,
     &  c16 = 0.777777777777777777777777777778d+00)
      parameter (
     &  b1 =   5.42937341165687622380535766363d-2,
     &  b6 =   4.45031289275240888144113950566d0,
     &  b7 =   1.89151789931450038304281599044d0,
     &  b8 =  -5.8012039600105847814672114227d0,
     &  b9 =   3.1116436695781989440891606237d-1,
     &  b10 = -1.52160949662516078556178806805d-1,
     &  b11 =  2.01365400804030348374776537501d-1,
     &  b12 =  4.47106157277725905176885569043d-2)
      parameter (
     &  bhh1 = 0.244094488188976377952755905512d+00,
     &  bhh2 = 0.733846688281611857341361741547d+00,
     &  bhh3 = 0.220588235294117647058823529412d-01)
      parameter (
     &  er 1 =  0.1312004499419488073250102996d-01,
     &  er 6 = -0.1225156446376204440720569753d+01,
     &  er 7 = -0.4957589496572501915214079952d+00,
     &  er 8 =  0.1664377182454986536961530415d+01,
     &  er 9 = -0.3503288487499736816886487290d+00,
     &  er10 =  0.3341791187130174790297318841d+00,
     &  er11 =  0.8192320648511571246570742613d-01,
     &  er12 = -0.2235530786388629525884427845d-01)
      parameter (
     &  a21 =    5.26001519587677318785587544488d-2,
     &  a31 =    1.97250569845378994544595329183d-2,
     &  a32 =    5.91751709536136983633785987549d-2,
     &  a41 =    2.95875854768068491816892993775d-2,
     &  a43 =    8.87627564304205475450678981324d-2,
     &  a51 =    2.41365134159266685502369798665d-1,
     &  a53 =   -8.84549479328286085344864962717d-1,
     &  a54 =    9.24834003261792003115737966543d-1,
     &  a61 =    3.7037037037037037037037037037d-2,
     &  a64 =    1.70828608729473871279604482173d-1,
     &  a65 =    1.25467687566822425016691814123d-1,
     &  a71 =    3.7109375d-2,
     &  a74 =    1.70252211019544039314978060272d-1,
     &  a75 =    6.02165389804559606850219397283d-2,
     &  a76 =   -1.7578125d-2)
      parameter (
     &  a81 =    3.70920001185047927108779319836d-2,
     &  a84 =    1.70383925712239993810214054705d-1,
     &  a85 =    1.07262030446373284651809199168d-1,
     &  a86 =   -1.53194377486244017527936158236d-2,
     &  a87 =    8.27378916381402288758473766002d-3,
     &  a91 =    6.24110958716075717114429577812d-1,
     &  a94 =   -3.36089262944694129406857109825d0,
     &  a95 =   -8.68219346841726006818189891453d-1,
     &  a96 =    2.75920996994467083049415600797d1,
     &  a97 =    2.01540675504778934086186788979d1,
     &  a98 =   -4.34898841810699588477366255144d1,
     &  a101 =   4.77662536438264365890433908527d-1,
     &  a104 =  -2.48811461997166764192642586468d0,
     &  a105 =  -5.90290826836842996371446475743d-1,
     &  a106 =   2.12300514481811942347288949897d1,
     &  a107 =   1.52792336328824235832596922938d1,
     &  a108 =  -3.32882109689848629194453265587d1,
     &  a109 =  -2.03312017085086261358222928593d-2)
      parameter (
     &  a111 =  -9.3714243008598732571704021658d-1,
     &  a114 =   5.18637242884406370830023853209d0,
     &  a115 =   1.09143734899672957818500254654d0,
     &  a116 =  -8.14978701074692612513997267357d0,
     &  a117 =  -1.85200656599969598641566180701d1,
     &  a118 =   2.27394870993505042818970056734d1,
     &  a119 =   2.49360555267965238987089396762d0,
     &  a1110 = -3.0467644718982195003823669022d0,
     &  a121 =   2.27331014751653820792359768449d0,
     &  a124 =  -1.05344954667372501984066689879d1,
     &  a125 =  -2.00087205822486249909675718444d0,
     &  a126 =  -1.79589318631187989172765950534d1,
     &  a127 =   2.79488845294199600508499808837d1,
     &  a128 =  -2.85899827713502369474065508674d0,
     &  a129 =  -8.87285693353062954433549289258d0,
     &  a1210 =  1.23605671757943030647266201528d1,
     &  a1211 =  6.43392746015763530355970484046d-1)
      parameter (
     &  a141 =  5.61675022830479523392909219681d-2,
     &  a147 =  2.53500210216624811088794765333d-1,
     &  a148 = -2.46239037470802489917441475441d-1,
     &  a149 = -1.24191423263816360469010140626d-1,
     &  a1410 =  1.5329179827876569731206322685d-1,
     &  a1411 =  8.20105229563468988491666602057d-3,
     &  a1412 =  7.56789766054569976138603589584d-3,
     &  a1413 = -8.298d-3)
      parameter (
     &  a151 =  3.18346481635021405060768473261d-2,
     &  a156 =  2.83009096723667755288322961402d-2,
     &  a157 =  5.35419883074385676223797384372d-2,
     &  a158 = -5.49237485713909884646569340306d-2,
     &  a1511 = -1.08347328697249322858509316994d-4,
     &  a1512 =  3.82571090835658412954920192323d-4,
     &  a1513 = -3.40465008687404560802977114492d-4,
     &  a1514 =  1.41312443674632500278074618366d-1,
     &  a161 = -4.28896301583791923408573538692d-1,
     &  a166 = -4.69762141536116384314449447206d0,
     &  a167 =  7.68342119606259904184240953878d0,
     &  a168 =  4.06898981839711007970213554331d0,
     &  a169 =  3.56727187455281109270669543021d-1,
     &  a1613 = -1.39902416515901462129418009734d-3,
     &  a1614 =  2.9475147891527723389556272149d0,
     &  a1615 = -9.15095847217987001081870187138d0)
      parameter (
     &  d41  = -0.84289382761090128651353491142d+01,
     &  d46  =  0.56671495351937776962531783590d+00,
     &  d47  = -0.30689499459498916912797304727d+01,
     &  d48  =  0.23846676565120698287728149680d+01,
     &  d49  =  0.21170345824450282767155149946d+01,
     &  d410 = -0.87139158377797299206789907490d+00,
     &  d411 =  0.22404374302607882758541771650d+01,
     &  d412 =  0.63157877876946881815570249290d+00,
     &  d413 = -0.88990336451333310820698117400d-01,
     &  d414 =  0.18148505520854727256656404962d+02,
     &  d415 = -0.91946323924783554000451984436d+01,
     &  d416 = -0.44360363875948939664310572000d+01)
      parameter (
     &  d51  =  0.10427508642579134603413151009d+02,
     &  d56  =  0.24228349177525818288430175319d+03,
     &  d57  =  0.16520045171727028198505394887d+03,
     &  d58  = -0.37454675472269020279518312152d+03,
     &  d59  = -0.22113666853125306036270938578d+02,
     &  d510 =  0.77334326684722638389603898808d+01,
     &  d511 = -0.30674084731089398182061213626d+02,
     &  d512 = -0.93321305264302278729567221706d+01,
     &  d513 =  0.15697238121770843886131091075d+02,
     &  d514 = -0.31139403219565177677282850411d+02,
     &  d515 = -0.93529243588444783865713862664d+01,
     &  d516 =  0.35816841486394083752465898540d+02)
      parameter (
     &  d61 =  0.19985053242002433820987653617d+02,
     &  d66 = -0.38703730874935176555105901742d+03,
     &  d67 = -0.18917813819516756882830838328d+03,
     &  d68 =  0.52780815920542364900561016686d+03,
     &  d69 = -0.11573902539959630126141871134d+02,
     &  d610 =  0.68812326946963000169666922661d+01,
     &  d611 = -0.10006050966910838403183860980d+01,
     &  d612 =  0.77771377980534432092869265740d+00,
     &  d613 = -0.27782057523535084065932004339d+01,
     &  d614 = -0.60196695231264120758267380846d+02,
     &  d615 =  0.84320405506677161018159903784d+02,
     &  d616 =  0.11992291136182789328035130030d+02)
      parameter (
     &  d71  = -0.25693933462703749003312586129d+02,
     &  d76  = -0.15418974869023643374053993627d+03,
     &  d77  = -0.23152937917604549567536039109d+03,
     &  d78  =  0.35763911791061412378285349910d+03,
     &  d79  =  0.93405324183624310003907691704d+02,
     &  d710 = -0.37458323136451633156875139351d+02,
     &  d711 =  0.10409964950896230045147246184d+03,
     &  d712 =  0.29840293426660503123344363579d+02,
     &  d713 = -0.43533456590011143754432175058d+02,
     &  d714 =  0.96324553959188282948394950600d+02,
     &  d715 = -0.39177261675615439165231486172d+02,
     &  d716 = -0.14972683625798562581422125276d+03)
      real(dp) y(n),y1(n),k1(n),k2(n),k3(n),k4(n),k5(n),k6(n)
      real(dp) k7(n),k8(n),k9(n),k10(n),atol(*),rtol(*)     
      dimension icomp(nrd),iwork(nrd+1)
      logical reject,last 
      integer ierr
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      
      interface
#include "num_fcn.dek"
#include "num_solout.dek"
      end interface
      
      real(dp), target :: rwork(3+8*nrd)
      real(dp), pointer :: cont(:)
      cont => rwork(3:3+8*nrd)
         
c *** *** *** *** *** *** ***
c  initialisations
c *** *** *** *** *** *** *** 
      facold=1.d-4  
      expo1=1.d0/8.d0-beta*0.2d0
      facc1=1.d0/fac1
      facc2=1.d0/fac2
      posneg=sign(1.d0,xend-x) 
c --- initial preparations   
      atoli=atol(1)
      rtoli=rtol(1)    
      last=.false. 
      hlamb=0.d0
      iasti=0
      nonsti=0
      ierr=0
      call fcn(n,x,h,y,k1,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; idid=-5; return; end if
      hmax=abs(hmax)     
      iord=8  
      if (h.eq.0.d0) h=hinit(n,fcn,x,y,xend,posneg,k1,k2,k3,iord,
     &                       hmax,atol,rtol,itol,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; idid=-5; return; end if
      nfcn=nfcn+2
      reject=.false.
      xold=x
      irtrn=1 
      if (iout.ge.1) then 
          hout=1.d0
          rwork(1) = xold
          rwork(2) = hout
          iwork(1) = nrd
          iwork(2:nrd+1) = icomp(1:nrd)

            do 662 j=1,nrd
               i=icomp(j)
               cont(j)=y(i)
               cont(j+nrd)=0
               cont(j+nrd*2)=0
               cont(j+nrd*3)=0
               cont(j+nrd*4)=0
               cont(j+nrd*5)=0
               cont(j+nrd*6)=0
               cont(j+nrd*7)=0
  662       continue 

          call solout(naccpt+1,xold,x,n,y,rwork,iwork,contd8,lrpar,rpar,lipar,ipar,irtrn)
          if (irtrn.lt.0) goto 79
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
c --- the twelve stages
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
  24  y1(i)=y(i)+h*(a41*k1(i)+a43*k3(i))  
      call fcn(n,x+c4*h,h,y1,k4,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 25 i=1,n 
  25  y1(i)=y(i)+h*(a51*k1(i)+a53*k3(i)+a54*k4(i))
      call fcn(n,x+c5*h,h,y1,k5,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 26 i=1,n 
  26  y1(i)=y(i)+h*(a61*k1(i)+a64*k4(i)+a65*k5(i))
      call fcn(n,x+c6*h,h,y1,k6,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 27 i=1,n 
  27  y1(i)=y(i)+h*(a71*k1(i)+a74*k4(i)+a75*k5(i)+a76*k6(i))
      call fcn(n,x+c7*h,h,y1,k7,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 28 i=1,n 
  28  y1(i)=y(i)+h*(a81*k1(i)+a84*k4(i)+a85*k5(i)+a86*k6(i)+a87*k7(i))  
      call fcn(n,x+c8*h,h,y1,k8,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 29 i=1,n 
  29  y1(i)=y(i)+h*(a91*k1(i)+a94*k4(i)+a95*k5(i)+a96*k6(i)+a97*k7(i)
     &   +a98*k8(i))
      call fcn(n,x+c9*h,h,y1,k9,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 30 i=1,n 
  30  y1(i)=y(i)+h*(a101*k1(i)+a104*k4(i)+a105*k5(i)+a106*k6(i)
     &   +a107*k7(i)+a108*k8(i)+a109*k9(i))
      call fcn(n,x+c10*h,h,y1,k10,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      do 31 i=1,n 
  31  y1(i)=y(i)+h*(a111*k1(i)+a114*k4(i)+a115*k5(i)+a116*k6(i)
     &   +a117*k7(i)+a118*k8(i)+a119*k9(i)+a1110*k10(i))
      call fcn(n,x+c11*h,h,y1,k2,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      xph=x+h
      do 32 i=1,n 
  32  y1(i)=y(i)+h*(a121*k1(i)+a124*k4(i)+a125*k5(i)+a126*k6(i)
     &   +a127*k7(i)+a128*k8(i)+a129*k9(i)+a1210*k10(i)+a1211*k2(i))
      call fcn(n,xph,h,y1,k3,lrpar,rpar,lipar,ipar,ierr)
      if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
      nfcn=nfcn+11
      do 35 i=1,n 
      k4(i)=b1*k1(i)+b6*k6(i)+b7*k7(i)+b8*k8(i)+b9*k9(i)
     &   +b10*k10(i)+b11*k2(i)+b12*k3(i)
  35  k5(i)=y(i)+h*k4(i)
c --- error estimation  
      err=0.d0
      err2=0.d0
      if (itol.eq.0) then   
        do 41 i=1,n 
        sk=atoli+rtoli*max(abs(y(i)),abs(k5(i)))
        erri=k4(i)-bhh1*k1(i)-bhh2*k9(i)-bhh3*k3(i)
        err2=err2+(erri/sk)**2
        erri=er1*k1(i)+er6*k6(i)+er7*k7(i)+er8*k8(i)+er9*k9(i)
     &      +er10*k10(i)+er11*k2(i)+er12*k3(i)
  41    err=err+(erri/sk)**2
      else
        do 42 i=1,n 
        sk=atol(i)+rtol(i)*max(abs(y(i)),abs(k5(i)))
        erri=k4(i)-bhh1*k1(i)-bhh2*k9(i)-bhh3*k3(i)
        err2=err2+(erri/sk)**2
        erri=er1*k1(i)+er6*k6(i)+er7*k7(i)+er8*k8(i)+er9*k9(i)
     &      +er10*k10(i)+er11*k2(i)+er12*k3(i)
  42    err=err+(erri/sk)**2
      end if  
      deno=err+0.01d0*err2
      if (deno.le.0.d0) deno=1.d0
      err=abs(h)*err*sqrt(1.d0/(n*deno))
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
         call fcn(n,xph,h,k5,k4,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
         nfcn=nfcn+1
c ------- stiffness detection
         if (mod(naccpt,nstiff).eq.0.or.iasti.gt.0) then
            stnum=0.d0
            stden=0.d0
            do 64 i=1,n 
               stnum=stnum+(k4(i)-k3(i))**2
               stden=stden+(k5(i)-y1(i))**2
 64         continue  
            if (stden.gt.0.d0) hlamb=abs(h)*sqrt(stnum/stden) 
            if (hlamb.gt.6.1d0) then
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
c ------- final preparation for dense output
         if (iout.ge.2) then
c ----    save the first function evaluations   
            do 62 j=1,nrd
               i=icomp(j)
               cont(j)=y(i)
               ydiff=k5(i)-y(i)
               cont(j+nrd)=ydiff
               bspl=h*k1(i)-ydiff
               cont(j+nrd*2)=bspl
               cont(j+nrd*3)=ydiff-h*k4(i)-bspl
               cont(j+nrd*4)=d41*k1(i)+d46*k6(i)+d47*k7(i)+d48*k8(i)
     &                  +d49*k9(i)+d410*k10(i)+d411*k2(i)+d412*k3(i)
               cont(j+nrd*5)=d51*k1(i)+d56*k6(i)+d57*k7(i)+d58*k8(i)
     &                  +d59*k9(i)+d510*k10(i)+d511*k2(i)+d512*k3(i)
               cont(j+nrd*6)=d61*k1(i)+d66*k6(i)+d67*k7(i)+d68*k8(i)
     &                  +d69*k9(i)+d610*k10(i)+d611*k2(i)+d612*k3(i)
               cont(j+nrd*7)=d71*k1(i)+d76*k6(i)+d77*k7(i)+d78*k8(i)
     &                  +d79*k9(i)+d710*k10(i)+d711*k2(i)+d712*k3(i)
   62       continue 
c ---     the next three function evaluations
            do 51 i=1,n 
  51           y1(i)=y(i)+h*(a141*k1(i)+a147*k7(i)+a148*k8(i)
     &            +a149*k9(i)+a1410*k10(i)+a1411*k2(i)+a1412*k3(i)
     &            +a1413*k4(i))
            call fcn(n,x+c14*h,h,y1,k10,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
            do 52 i=1,n 
  52           y1(i)=y(i)+h*(a151*k1(i)+a156*k6(i)+a157*k7(i)
     &            +a158*k8(i)+a1511*k2(i)+a1512*k3(i)+a1513*k4(i)
     &            +a1514*k10(i))
            call fcn(n,x+c15*h,h,y1,k2,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
            do 53 i=1,n 
  53           y1(i)=y(i)+h*(a161*k1(i)+a166*k6(i)+a167*k7(i)
     &            +a168*k8(i)+a169*k9(i)+a1613*k4(i)+a1614*k10(i)
     &            +a1615*k2(i))
            call fcn(n,x+c16*h,h,y1,k3,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then; hnew=h/facc1; h=hnew; goto 1; end if
            nfcn=nfcn+3 
c ---     final preparation
            do 63 j=1,nrd
               i=icomp(j)
               cont(j+nrd*4)=h*(cont(j+nrd*4)+d413*k4(i)+d414*k10(i)
     &            +d415*k2(i)+d416*k3(i))
               cont(j+nrd*5)=h*(cont(j+nrd*5)+d513*k4(i)+d514*k10(i)
     &            +d515*k2(i)+d516*k3(i))
               cont(j+nrd*6)=h*(cont(j+nrd*6)+d613*k4(i)+d614*k10(i)
     &            +d615*k2(i)+d616*k3(i))
               cont(j+nrd*7)=h*(cont(j+nrd*7)+d713*k4(i)+d714*k10(i)
     &            +d715*k2(i)+d716*k3(i))
  63        continue
            hout=h
         end if
         do 67 i=1,n
         k1(i)=k4(i)
  67     y(i)=k5(i)
         xold=x
         x=xph
         if (iout.ge.1) then
            rwork(1) = xold
            rwork(2) = hout
            iwork(1) = nrd
            iwork(2:nrd+1) = icomp(1:nrd)
            call solout(naccpt+1,xold,x,n,y,rwork,iwork,contd8,lrpar,rpar,lipar,ipar,irtrn)
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
      if (lout.gt.0) write(lout,979) x   
      if (lout.gt.0) write(lout,*)' step size too small, h=',h
      idid=-3
      return
  78  continue
      if (lout.gt.0) write(lout,979) x   
      if (lout.gt.0) write(lout,*)
     &     ' more than nmax =',nmax,'steps are needed' 
      idid=-2
      return
  79  continue
      !if (lout.gt.0) write(lout,979) x
 979  format(' exit of dop853 at x=',e18.4) 
      idid=2
      return
      end subroutine dp86co
c
      function hinit(n,fcn,x,y,xend,posneg,f0,f1,y1,iord,
     &                       hmax,atol,rtol,itol,lrpar,rpar,lipar,ipar,ierr)
c ----------------------------------------------------------
c ----  computation of an initial step size guess
c ----------------------------------------------------------
      implicit real(dp) (a-h,o-z)
      
      real(dp), intent(inout) :: x 
      real(dp), intent(inout), dimension(:) :: y, f0, f1 ! (n)
      real(dp), intent(in) :: xend
      real(dp), intent(in) :: rtol(*)
      real(dp), intent(in) :: atol(*)
      integer, intent(in) :: itol
      integer, intent(in) :: lrpar, lipar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      integer, intent(out)  :: ierr
      
      dimension y1(n)
      
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
        dnf=dnf+(f0(i)/sk)**2
  10    dny=dny+(y(i)/sk)**2 
      else
        do 11 i=1,n 
        sk=atol(i)+rtol(i)*abs(y(i))
        dnf=dnf+(f0(i)/sk)**2
  11    dny=dny+(y(i)/sk)**2 
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
  15    der2=der2+((f1(i)-f0(i))/sk)**2   
      else
        do 16 i=1,n 
        sk=atol(i)+rtol(i)*abs(y(i))
  16    der2=der2+((f1(i)-f0(i))/sk)**2   
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



      real(dp) function contd8(ii,x,rwork,iwork,ierr)
! ----------------------------------------------------------
!     this function can be used for continuous output in connection
!     with the output-subroutine for dop853. it provides an
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
      
         nd = iwork(1)
         icomp => iwork(2:nd+1)
         xold = rwork(1)
         h = rwork(2)
         con => rwork(3:3+8*nd)
         
         ! ----- compute place of ii-th component 
         i=0 
         do j=1,nd 
            if (icomp(j).eq.ii) then
               i=j; exit
            end if
         end do
         if (i.eq.0) then
            contd8 = 0
            ierr = -1
            return
         end if
         ierr=0
         s=(x-xold)/h
         s1=1.d0-s
         conpar=con(i+nd*4)+s*(con(i+nd*5)+s1*(con(i+nd*6)+s*con(i+nd*7)))
         contd8=con(i)+s*(con(i+nd)+s1*(con(i+nd*2)+s*(con(i+nd*3)+s1*conpar)))
         
      end function contd8


      end module mod_dop853
