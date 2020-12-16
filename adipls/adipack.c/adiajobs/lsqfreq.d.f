      program main
c
c  dog-like least square fit to oscillation frequencies
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
c  Modified 16/8/93, to output residual from combined fit to
c  file ttt.lsqfreq.resid
c
      implicit double precision (a-h, o-z)
      character*280 form, fin, fout
      dimension cs(50),x(500),y(500),ll(100),alpha(100),beta(100),
     *  gamma(100),cc(10),frobs(3000),lnobs(2,3000),lval(100),
     *  ccl(11,50),frq12(2,50),form(5),ccal(2),ccbt(2),ccgm(2),
     *  ls(500),ns(500),freqs(500),xs(500)
c..      namelist/exec/ icase,nmodel,
c..     *  lw1,lw2,nw1,nw2,xw1,xw2,frqw1,frqw2,
c..     *  x0,kpol,lmaxft,
c..     *  iobsfr,irdfrq,irdiag
      data form 
     *  /'unformatted','unformatted','formatted','unformatted',
     *  'unformatted'/
c
c  icase: type of modes used for input.
c  icase = 1: grand summary.
c  icase = 2: short summary.
c  icase = 3: observed frequencies.
c  icase = 4: grand summary, use straight eigenfrequency
c  icase = 5: grand summary, use Richardson extrapolated eigenfrequency
      icase=1
c  nmodel: if nmodel .ge. 1, consider only model no nmodel
      nmodel=0
c  lw1, lw2, nw1, nw2, xw1, xw2, frqw1, frqw2: window in l, order,
c  n + l/2 and frequency (in muhz)
      lw1=0
      lw2=1000
      nw1=1
      nw2=100
      xw1=0
      xw2=100
      frqw1=2400
      frqw2=4200
      x0=22
c  kpol: degree of polynomial fitted to frequencies
      kpol=1
c  lmaxft: maximum l-value used in least-squares fit to coefficients
c     if lmaxft .le. 0, use all values available
      lmaxft = 0
c  iobsfr: obsolete. Kept for consistency with old input files.
      iobsfr=0
c  irdfrq, irdiag: for iobsfr = 1 read frequencies only if irdfrq .ne. 0
c          and frequencies have not already been read. for irdiag = 1
c          print frequencies
      irdfrq=0
      irdiag=0
c
c  initialize flag for first run
c
      ifsrun=1
c
c  open data sets
c
      write(6,*) 'Enter input file'
      read(5,'(a)') fin
      write(6,*) 'Enter output file'
      read(5,'(a)') fout
c
    5 continue
c
c..    5 read(5,exec,err=90,end=90)
c
      write(6,*) 'icase, nmodel?'
      write(6,*) icase, nmodel
      read(5,*,err=90,end=90) icase, nmodel
      write(6,*) 'lw1,lw2,nw1,nw2,xw1,xw2,frqw1,frqw2?'
      write(6,*) lw1,lw2,nw1,nw2,xw1,xw2,frqw1,frqw2
      read(5,*) lw1,lw2,nw1,nw2,xw1,xw2,frqw1,frqw2
      write(6,*) 'x0,kpol,lmaxft?'
      write(6,*) x0,kpol,lmaxft
      read(5,*) x0,kpol,lmaxft
      write(6,*) 'iobsfr,irdfrq,irdiag?'
      write(6,*) iobsfr,irdfrq,irdiag
      read(5,*) iobsfr,irdfrq,irdiag
c
c..      write(6,100)
c..      write(6,exec)
c
      write(6,*) 'icase, nmodel'
      write(6,*) icase, nmodel
      write(6,*) 'lw1,lw2,nw1,nw2,xw1,xw2,frqw1,frqw2'
      write(6,*) lw1,lw2,nw1,nw2,xw1,xw2,frqw1,frqw2
      write(6,*) 'x0,kpol,lmaxft'
      write(6,*) x0,kpol,lmaxft
      write(6,*) 'iobsfr,irdfrq,irdiag'
      write(6,*) iobsfr,irdfrq,irdiag
c
c  set data sets
c
      open(2,file=fin,status='old',form=form(icase))
      open(10,file=fout,status='unknown')
      open(15,file='ttt.lsqfreq.sum',status='unknown')
      open(16,file='ttt.lsqfreq.resid',status='unknown')
c
      lfin=length(fin)
      write(6,101) fin(1:lfin), nmodel, icase
      write(10,101) fin(1:lfin), nmodel, icase
      write(16,102) fin(1:lfin), nmodel, icase
c
c  read modes (assumed ordered with l first)
      l1=0
      lp=-1
      n=0
      n1=0
      ktitle=1
c
c  test for using observed frequencies
c
      if(iobsfr.ne.1) go to 7
c
c  read observed frequencies
c
      if(irdfrq.eq.0.and.ifsrun.eq.0) go to 6
c
      call rdofrq
c
      ifsrun=0
c
    6 no=0
      go to 15
c
c  prepare for reading theoretical frequencies
c
    7 rewind 2
      go to 10
c
c  entry point for new l-value
c  ***************************
c
    8 if(iobsfr.eq.1) go to 17
      go to 12
c
c  step in theoretical frequencies
c
   10 iend=1
      call rdfrqm(icase,2,nmodel,cs,l,nord,sig,frq,ekin,ierr)
      if(ierr.gt.0) go to 20
c
c  test for model line. skip for now
c
      if(l.lt.0) go to 10
c
      xx=float(nord)+float(l)/2
c
      if(nord.lt.nw1.or.nord.gt.nw2.or.xx.lt.xw1.or.xx.gt.xw2.or.
     *  frq.lt.frqw1.or.frq.gt.frqw2.or.l.lt.lw1) go to 10
c
c  test for end when l exceeds lw2
c
      if(l.gt.lw2) go to 20
c
c  test for writing title
c  this has to be fixed up, dependent on data set case
c
      if(ktitle.ne.1) go to 11
      ktitle=0
      write(6,105) (cs(i),i=2,5),x0
      write(10,105) (cs(i),i=2,5),x0
      write(16,106) (cs(i),i=2,5),x0
c
c  test for setting lp
c
   11 if(lp.eq.-1) lp=l
c
      iend=0
      if(l.ne.lp) go to 20
c
c  continue after analysis
c
   12 lx=l
      n=n+1
      x(n)=float(nord)+float(l)/2-x0
      y(n)=frq
      n1=n1+1
      ls(n1)=l
      ns(n1)=nord
      freqs(n1)=frq
      xs(n1)=x(n)
      go to 10
c
c  step in observed frequencies
c
   15 iend=1
c
c  test for printing title
c
      if(ktitle.eq.1) go to 16
      write(6,110) x0
      write(10,110) x0
c
   16 no=no+1
      if(no.gt.nnobs) go to 20
      l=lnobs(1,no)
      nord=lnobs(2,no)
      frq=frobs(no)
      xx=nord+float(l)/2
      if(nord.lt.nw1.or.nord.gt.nw2.or.xx.lt.xw1.or.xx.gt.xw2.or.
     *  frq.lt.frqw1.or.frq.gt.frqw2.or.l.lt.lw1) go to 15
c
c  test for end when l exceeds lw2
c
      if(l.gt.lw2) go to 20
c
c  test for setting lp
c
      if(lp.eq.-1) lp=l
c
      iend=0
      if(l.ne.lp) go to 20
c
c  continue after analysis
c
   17 lx=l
      n=n+1
      x(n)=lnobs(2,no)+float(l)/2-x0
      y(n)=frobs(no)
      go to 15
c
c  new l
c
   20 if(iend.eq.1) go to 30
      lp=l
      if(n.eq.0) go to 8
c
c  least square analysis
c  *********************
c
   30 write(6,112) lx
      call lsqpol(x,y,n,1,1,kpol,cc,rmsres,1)
      l1=l1+1
      ll(l1)=lx
      alpha(l1)=cc(1)
      beta(l1)=cc(2)
      gamma(l1)=cc(3)
      kpol1=kpol+1
      kpol2=kpol1+1
      do 35 i=1,kpol1
   35 ccl(i,l1)=cc(i)
      ccl(kpol2,l1)=rmsres
c
c  actual range in frequency used in analysis
c
      frq12(1,l1)=y(1)
      frq12(2,l1)=y(n)
c
      n=0
      if(iend.eq.0) go to 8
c
c  final output and analysis
c  *************************
c
      nn1=n1
c
      write(6,115) kpol,nw1,nw2,xw1,xw2,frqw1,frqw2
      write(10,115) kpol,nw1,nw2,xw1,xw2,frqw1,frqw2
      do 40 i=1,l1
      yy=ll(i)*(ll(i)+1)
      x(i)=yy
      write(6,120)  ll(i),yy,(frq12(j,i),j=1,2),(ccl(j,i),j=1,kpol2)
   40 write(10,120) ll(i),yy,(frq12(j,i),j=1,2),(ccl(j,i),j=1,kpol2)
c
c  output alpha, beta, gamma and epsilon
c
      write(10,122)
      do 45 i=1,l1
      epsil=ccl(1,i)/ccl(2,i)-x0
   45 write(10,124) ll(i),(ccl(j,i),j=1,3),epsil
c
c  set number of values included in fit to coefficients
c
      if(lmaxft.le.0) then
        lf1=l1
      else
        lf1=lmaxft+1
      end if
c
      lf=lf1-1
c
      write(6,130) lf
      call lsqpol(x,alpha,lf1,1,1,1,ccal,rmsres,1)
      write(10,150) lf,ccal(1),ccal(2)
      write(6,150) lf,ccal(1),ccal(2)
      al0=ccal(1)
c
      write(6,140) lf
      call lsqpol(x,beta,lf1,1,1,1,ccbt,rmsres,1)
      write(10,160) ccbt(1),ccbt(2)
      write(6,160) ccbt(1),ccbt(2)
      bet0=ccbt(1)
c
      if(kpol.ge.2) then
        write(6,165) lf
        call lsqpol(x,gamma,lf1,1,1,1,ccgm,rmsres,1)
        write(10,167) lf,ccgm(1),ccgm(2)
        write(6,167) lf,ccgm(1),ccgm(2)
      end if
c
      epsi=al0/bet0-x0
      write(6,170) epsi
      write(10,170) epsi
c
c  test for setting extrapolated values
c
      if(lmaxft.gt.0) then
        ll1=(lmaxft+1)*(lmaxft+2)
        alphl=ccal(1)+ll1*ccal(2)
        betal=ccbt(1)+ll1*ccbt(2)
        gamml=ccgm(1)+ll1*ccgm(2)
        epsil=alphl/betal-x0
        write(10,180) lf1,alphl,betal,gamml,epsil
      end if
c
      write(15,185) icase, ccal(1), ccbt(1), -ccal(2), -6.*ccal(2),
     *  epsi, fin(1:lfin)
c
c  set and output residuals
c
      write(16,187)
      do 60 n=1,nn1
      lln=ls(n)*(ls(n)+1)
      alphal=ccal(1)+lln*ccal(2)
      betal=ccbt(1)+lln*ccbt(2)
      gammal=ccgm(1)+lln*ccgm(2)
      dx=xs(n)
      frqfit=alphal+dx*(betal+dx*gammal)
  60  write(16,190) ls(n),ns(n),freqs(n),freqs(n)-frqfit
c
      go to 5
c
   90 continue
      stop
  100 format(1h1)
  101 format(' Input file: ',a/' nmodel =',i4/' icase =',i3/)
  102 format('# Input file: ',a/'# nmodel =',i4/'# icase =',i3/'#')
  105 format(///' analysis of theoretical frequencies'//
     *  ' the model has mass =',1pe13.5,'  radius =',e13.5,
     *  ' pc =',e13.5,'  rhoc =',e13.5//' analysis uses x0 =',
     *  0pf8.2)
  106 format('#'/'# analysis of theoretical frequencies'/'#'/
     *  '# the model has mass =',1pe13.5,'  radius =',e13.5/
     *  '# pc =',e13.5,'  rhoc =',e13.5/'#'/'# analysis uses x0 =',
     *  0pf8.2)
  110 format(///
     *  ' analysis of observed frequencies'//' analysis uses x0 =',f8.2)
  112 format(////' l =',i4/1x,7(1h*))
  115 format(//' results of least squares analysis ',
     *  'with polynomial of degree',i3/
     *  ' including modes with'/
     *  ' order     between',i10,' and ',i10/
     *  ' n + l/2   between',f10.2,' and ',f10.2/
     *  ' frequency between',f10.2,' and ',f10.2,' muhz'//
     *  ' l,l(l+1), actual nu1, nu2 (muhz),',
     *  ' alpha,beta,gamma,...,rms residual:'/)
  120 format(i4,f10.1,2f10.2,7f12.5)
  122 format(//' l, alpha, beta, gamma, epsilon = alpha/beta - x0'/)
  124 format(i5,f8.2,f8.2,f8.4,f10.4)
  130 format(//' least square analysis on alpha, l = 0 - ',i3,':'//)
  140 format(//' least square analysis on beta, l = 0 - ',i3,':'//)
  150 format(///' least squares fits to alpha and beta,',
     *  ' l = 0 - ',i3//
     *  ' alpha =',f12.3,' + l(l+1)*',f10.5)
  160 format(//' beta  =',f12.3,' + l(l+1)*',f10.5)
  165 format(//' least square analysis on gamma, l = 0 - ',i3,':'//)
  167 format(///' least squares fit to gamma, l = 0 - ',i3//
     *  ' gamma =',f12.5,' + l(l+1)*',f10.5)
  170 format(//' epsilon = alpha0/beta0 - x0 =',f10.5)
  180 format(//' extrapolated values'/
     *  ' l, alpha, beta, gamma, epsilon = alpha/beta - x0'//
     *  i5,f8.2,f8.2,f8.4,f10.4)
  185 format(
     *  '#  case, nu_0, Delta nu_0, D_0, 6*D_0, epsilon_0, modeset:'/
     *  i3,f12.5,f12.6,3f12.8,2x,a)
  187 format('#'/'# l, order, frequency, residual (microHz)'/'#')
  190 format(2i5,f11.3,f11.5)
      end
      subroutine rdofrq
c  dummy subroutine
      return
      end
