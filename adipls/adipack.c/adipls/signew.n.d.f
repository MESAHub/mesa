      double precision function signew(itrsig,nsig,iselct,
     *  inomde,sigp,dfsig,ids,el,icry)
c
c  set new trial squared frequency.
c  ********************************
c
c  for itrsig = 1 sig is found from previous value sigp and dfsig.
c  meaning of dfsig depends on nsig:
c   nsig = 1: sig = sigp + dfsig
c   nsig = 2: dfsig is increment in frequency (i.e. in sqrt(sig))
c   nsig = 3: dfsig is increment in 1/(frequency) (i.e. in 1/sqrt(sig))
c
c  for itrsig .gt. 1 sig is found from grand summary (itrsig = 2
c  or 3), short summary (itrsig = 4 or 5) or cyclic frequencies from
c  file on the form of observed data (itrsig = 6 or 7), 
c  residing on d/s ids. 
c
c  for itrsig .lt. -1, sig is found as for abs(itrsig), but using
c  single precision files.
c
c  if iselct .ne. -1 inomde determines which mode is used.
c  if itrsig is even mode no inomde is taken, and el is reset to the value for
c  this mode.
c  if itrsig is odd the mode with the given value of el and of
c  order inomde is used.
c
c  if iselct = -1 the next mode on the file is considered.
c  if itrsig is even this mode is used, and el is reset to the value for
c  this mode.
c  if itrsig is odd the mode is only used if its value of el agrees
c  with the value of el in the argument list. otherwise icry is set
c  to -1.
c
c  The status of the call is returned in icry:
c  icry = -1: A mode with the specified order and degree was not found,
c             or the specified mode was not in the prescribed window.
c             (It makes sense to try the next mode).
c  icry = -2: Number of modes on input file exceeded. Terminate for
c             this case.
c  icry = -3: End of file reached on trial-frequency dataset
c             this case.
c
c  modified 12/7/1985 to let short summary have length 7 words
c
c  modified 13/8/87 to standardize output
c
c  modified 18/8/90 to include option of using observed data.
c
c  Modified 2/7/95, to include windowing in trial frequency and
c  degree, as determined by parameters in common/cwindw/
c
c  Modified 28/4/99, to allow real (non-integer) degree read from
c  file of observed frequencies (replacing rdfreq by rdfrqr)
c
c  ................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical notwin
      dimension gs(50)
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
      common/rhsdat/ elcom,ell,alb,els,el1,sigcom,anres,perfac,data(8)
      common/ccgrav/ cgrav
      common/cwindw/ eltrw1, eltrw2, sgtrw1, sgtrw2
      common/cmdtst/ iordtr,icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *  iorwn1, iorwn2, frlwn1, frlwn2
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data idsp, ninit/-1, 10 000 000/
c
c..      write(6,*) 'Entering signew with itrsig, inomde, el =',
c..     *  itrsig, inomde, el
c..      write(6,*) 'ninit =',ninit
c
      itrsga=iabs(itrsig)
c
      icnrew=0
      icry=1
c
      if(itrsga.eq.1) then
c
c  find sig from sigp and dfsig
c
        if(nsig.eq.1) then
c
c  step in sig
c
          sign=sigp+dfsig
c
        else if(nsig.eq.2.or.nsig.eq.3) then
          isig=1
          if(sigp.lt.0) isig=-1
          sig=sqrt(isig*sigp)
          if(nsig.eq.2) then
c
c  step in sqrt(sig)
c
            sig=sig+dfsig
          else
c
c  step in 1/sqrt(sig)
c
            sig=sig/(1+dfsig*sig)
          end if
c
          sign=isig*sig*sig
        else
c
c  nsig outside allowed range. 
c
          write(istdou,100) nsig
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,100) nsig
          icry=-1
          return
        end if
c
c  test for sig in window
c
        if(notwin(sgtrw1,sgtrw2,sign)) then
          if(istdpr.gt.0) write(istdpr,102) sign, sgtrw1, sgtrw2
          icry = -1
        else
          signew = sign
        end if
c
        return
c
      end if
        
c
c                    -----------------------------------
c
c  find sig from summary on d/s ids
c  ********************************
c
c  test for itrsig in range
c
      if(itrsga.gt.7) then
c
c  write diagnostics
c
        icry=-1
        write(istdou,130) itrsig
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,130) itrsig
        return
c
      end if
c
c  determine case
c
      itr=itrsga/2
      itrs=isign(itr,itrsig)
      imd=mod(itrsga,2)
c
c  test for new d/s or otherwise starting from beginning
c
      if(ids.ne.idsp.or.(imd.eq.0.and.inomde.lt.ninit)) then
	if(istdpr.gt.0) write(istdpr,*) 'Rewinding on entry'
        rewind ids
	icnrew=1
        ninit=0
        ordn=epsofl
        eln=-1
        idsp=ids
        if(itr.eq.3) then
c
c  set factor for going from cyclic frequencies (in microhz) to
c  dimensionless squared frequency
c
          sigfct=3.947842d-11*data(2)**3/(cgrav*data(1))
        end if
      end if
c
      nstop=ninit-1
      n=ninit
c
c  test for reading next mode on file
c
   40 if(n.le.0.or.iselct.eq.-1) then
c
        n=n+1
   42   call rdfrqr(itrs,ids,gs,elrn,nord,sign,frqn,ekin,ierr)
c
c  test for model record in short summary
c
        if(elrn.lt.0) go to 42
c
        ordn=nord
        if(itr.eq.1) then
          eln=gs(18)
        else if(itr.eq.2) then
          eln=gs(1)
        else 
          eln=elrn
          sign=sigfct*frqn*frqn
        end if
      end if
c
c  test for taking next mode on file
c
      if(iselct.eq.-1) then
c
c  use next mode on file. test for test on value of el
c
        ninit=n
        if(imd.eq.1.and.abs(eln-el).gt.1.e-10) then
c
c  wrong el
c
          icry=-1
c
c  Test on el and sig in window
c
        else if(notwin(sgtrw1,sgtrw2,sign)) then
          write(istdou,102) sign, sgtrw1, sgtrw2
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,102) 
     *      sign, sgtrw1, sgtrw2
          icry = -1
        else if(notwin(eltrw1,eltrw2,eln)) then
          write(istdou,104) eln, eltrw1, eltrw2
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,104) 
     *      eln, eltrw1, eltrw2
          icry = -1
        else
          signew = sign
	  iordtr = ordn
          if(imd.eq.0) el = eln
        end if
        return
      end if
c
c  otherwise test for correct mode number or order and degree
c
   50 if((imd.eq.0.and.n.eq.inomde).or.(imd.eq.1.and.
     *  abs(el-eln).lt.epsprc.and.abs(inomde-ordn).lt.0.1)) then
c
c  correct target mode found. Test on el and sig in window
c
        if(notwin(sgtrw1,sgtrw2,sign)) then
          write(istdou,102) sign, sgtrw1, sgtrw2
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,102) 
     *      sign, sgtrw1, sgtrw2
          icry = -1
        else if(notwin(eltrw1,eltrw2,eln)) then
          write(istdou,104) eln, eltrw1, eltrw2
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,104) 
     *      eln, eltrw1, eltrw2
          icry = -1
        else
          signew = sign
	  iordtr = ordn
          if(imd.eq.0) el = eln
        end if
c
        ninit=n
        return
c
c  test whether whole d/s has been scanned
c
      else if(n.eq.nstop) then
        if(imd.eq.0) then
          write(istdou,110) inomde,n,ids
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) 
     *      inomde,n,ids
	  icry=-2
        else
          write(istdou,120) el,inomde,ids
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) 
     *      el,inomde,ids
          icry=-1
        end if
        ninit=n
        return
c
      else
c
c  read next mode
c
   55   call rdfrqr(itrs,ids,gs,elrn,nord,sign,frqn,ekin,ierr)
c
c  test for model record in short summary
c
        if(elrn.lt.0) go to 55
c
        if(ierr.gt.0) go to 75
        n=n+1
        ordn=nord
        if(itr.eq.1) then
          eln=gs(18)
        else if(itr.eq.2) then
          eln=gs(1)
        else 
          eln=elrn
          sign=sigfct*frqn*frqn
        end if
c
        go to 50
c
      end if
c
c  end of d/s ids reached
c
   75 if(ierr.eq.2) then
        write(istdou,140) ids
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,140) ids
        icry=-3
c
c  test whether scan started at beginning of d/s
c
      else if(ninit.le.1) then
        if(imd.eq.0) then
          write(istdou,110) inomde,n,ids
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) 
     *      inomde,n,ids
          icry=-2
        else
          write(istdou,120) el,inomde,ids
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) 
     *      el,inomde,ids
          icry=-1
        end if
        ninit=n
        return
c
c  test for excessive number of rewinds (should not occur)
c
      else if(icnrew.eq.1) then
        write(istdou,145) ids
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,145) ids
        icry=-2
        return
c
c  rewind d/s and start from beginning
c
      else
        rewind ids
	if(istdpr.gt.0) 
     *  write(istdpr,*) 'Rewind after reaching end of file'
	icnrew=icnrew+1
        ninit=0
        n=0
        go to 40
      end if
c
  100 format(//1x,10(1h*),' nsig =',i10,' not allowed in signew')
  102 format(/' ***** Error in signew. sigma**2 = ',1pe11.3,
     *  ' is not in window ',2e11.3)
  104 format(/' ***** Error in signew. degree = ',1pe11.3,
     *  ' is not in window ',2e11.3)
  110 format(//' ***** Failure in signew. inomde =',i5,
     *  ' .gt. number of modes =',i5,' on d/s',i3)
  120 format(//' ***** Failure in signew. no mode with l =',
     *  f10.3,' and order =',i5,' on d/s',i3)
  130 format(//' ***** Error in signew. itrsig =',i10,
     *  ' has not been implemented')
  140 format(//' ***** Error in call of s/r rdfrqr on d/s',i3,
     *  ' from s/r signew.'/
     *         '       Execution terminated.')
  145 format(//' ***** Error. Excessive number of rewinds of d/s',i3,
     *  ' from s/r signew.'/
     *         '       Execution terminated.')
      end
