      subroutine readml(ifind,xmod,mdintg,ids,idsp,in,data,x,aa,ia,nn,
     *  nprmod,ivar,mlname,nwmod,icry)
c
c  reads model xmod on d/s ids and takes every in-th point
c
c  Returns as diagnostics:
c  nwmod: nwmod = 1 if new model has been set, 0 if no model was set.
c  icry : icry =  1 for succesful read.
c         icry = -1 if read or taking subset of points failed.
c         icry =  2 for attempted interpolation between models
c                   in sequence, with non-integral xmod
c                   (not yet implemented)
c
c  modified 13/8/87 to standardize output
c
c  modified 12/2/89 to include option for Richardson extrapolation
c  this requires that total number of points, excluding the centre
c  and a possibly singular surface, is odd.
c
c  modified 18/2/90: correct setting of turbulent pressure ratio in aa(6,.)
c     Previously aa(6,n) was set to 1/(x**3*aa(1,n)) for x .ge. 0.999.
c     Now aa(6,n) is set to 1, unless iturpr = 1
c  modified 19/1/95: read 6 columns if iturpr = 2, to include
c     CR turbulent pressure corrections. Also move gtilde/g
c     (for iturpr = 1) from aa(6,n) to aa(10,n)
c  modified 2/2/95, (partially revoking previous modification)
c     to check for number of variables in aa, depending on data(8).
c     (Also move open of model file to this routine).
c
c  modified 11/10/95, to keep innermost point fixed when resetting
c     model to allow for Richardson extrapolation, in case
c     of envelope model
c
c  modified 16/2/98, to allow reading in general amdl format,
c  depending on the value of data(8)
c
c  modified 17/2/98, to set iturpr to 8, if model contains 
c  g/(g tilde) in aa(6,.), as flagged by data(8).
c
c  modified 7/7/05, adding ivar to the argument list (as it should be
c  for consistency with calling routine!)
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
      logical sincen, sinsur
      dimension data(8),x(1),aa(ia,nn),mlname(4)
      common/worksp/ aa1(iaa1,1)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data iblank/8h        /
      data inp /-1/
      data amsun /1.989d33/
c
c  open file for model input
c
      if(ids.ne.idsp) then 
	if(idsp.gt.0) close(idsp)
	call openf(ids,'o','u')
c
c  read start of file, to check for number of variables
c
        read(ids) nmod,nn,data
        close(ids)
	call openf(ids,'o','u')
	idata8 = idint(data(8)+0.1)
	if(idata8.ge.100) then
	  ivar = 8
        else if(idata8.ge.10) then
	  ivar = 6
        else
	  ivar = 5
        end if
c
	ivarrd=ivar+1
c
	idsp=ids
c
      end if
c
      icry=1
      nwmod=0
c
c  test for non-integer xmod
c
      nummod=xmod
      frmod=xmod-nummod
      if(frmod.gt.0.and.ifind.eq.2) icry=2
c
      if(ifind.eq.2) go to 201
      if(ifind) 2,205,2
  201 if(ids-idsp) 205,203,205
c  read model
    2 read(ids,end=90,err=90) nmod,nn,data,
     *   ((aa1(i,n),i=1,ivarrd),n=1,nn)
c
      nwmod=1
c
      nrd=nrd+1
      idsp=ids
      if(nummod.eq.0.or.ifind.ne.2) go to 215
  203 if(nrd-nummod) 2,210,205
c
  205 rewind ids
      nrd=0
      go to 2
c
c  test whether new model has been read
c
  210 if(nwmod.eq.1) go to 215
c
c  no model read. test for same in
c
      if(in.eq.inp) return
c
c  otherwise have to re-read model
c
      go to 205
c
c   end of reading model
c   ********************
c
c  test for singular centre and/or surface
c
  215 sincen=aa1(1,1).eq.0
      sinsur=data(7).ge.0
      nsin=0
      if(sincen) nsin=nsin+1
      if(sinsur) nsin=nsin+1
c
c  test for inclusion of g/(g tilde)
c
	if(mod(idata8/10,10).eq.2) then
	  iggt=1
	  iturpr=8
	  write(istdou,105) aa1(7,nn)
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,105) 
     *      aa1(7,nn)
        else
	  iggt=0
        end if
c
c  take every in-th point in model
c  if iriche = 1, in addition need to check that number of 
c  non-singular points is odd
c
      inp=in
      if(in.eq.1) then
c
c  test for number of nonsingular points
c     
        if(iriche.ne.1.or.mod(nn-nsin,2).eq.1) then
          nshift=0
        else
          nshift=1
        end if
        nnr=nn
        if(nshift.ne.0) then
	  write(istdou,111) nn,nn-nshift
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,111) 
     *      nn,nn-nshift
          nn=nn-nshift
        end if
        if(sincen) then
          x(1)=aa1(1,1)
          do 21 i=1,ivar
   21     aa(i,1)=aa1(i+1,1)
          do 22 n=2,nnr
          n1=n+nshift
          x(n)=aa1(1,n1)
          do 22 i=1,ivar
   22     aa(i,n)=aa1(i+1,n1)
        else
          do 25 n=1,nnr
          if(n.eq.1) then
            n1=1
          else
            n1=n+nshift
	  end if
          x(n)=aa1(1,n1)
          do 25 i=1,ivar
   25     aa(i,n)=aa1(i+1,n1)
        end if
c
      else
c
c  take every in-th point
c
        nnr=nn
        if(sincen) then
          nn=(nn-2)/in+2
        else
          nn=(nn-1)/in+1
        end if
c
c  test for reducing for Richardson extrapolation
c
        if(iriche.eq.1.and.mod(nn-nsin,2).eq.0) then
	  write(istdou,111) nn,nn-1
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,111) nn,nn-1
	  nn=nn-1
        end if
c
        if(sincen) then
          x(1)=aa1(1,1)
          do 42 i=1,ivar
   42     aa(i,1)=aa1(i+1,1)
          n1=nnr-in*(nn-1)
          nstart=2
        else
          n1=nnr-in*(nn+1)
          nstart=1
        end if
c
        do 45 n=nstart,nn
        n1=n1+in
        x(n)=aa1(1,n1)
        do 45 i=1,ivar
   45   aa(i,n)=aa1(i+1,n1)
c
        if(n1.ne.nnr) then
          write(istdou,112) in, nnr, nn
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,112) 
     *      in, nnr, nn
          icry=-1
          go to 95
        end if
c
      end if
c
c  set g/gtilde (=1 in models without turbulent pressure)
c
      if(iturpr.eq.1) then
        do 55 n=1,nn
        if(x(n).lt.0.999) then
          ggt=1
        else
          ggt=1./(x(n)*x(n)*x(n)*aa(1,n))
        end if
   55   aa(10,n)=ggt
      else if(iggt.eq.1) then
	do 56 n=1,nn
   56   aa(10,n)=aa(6,n)
      else
	do 57 n=1,nn
   57   aa(10,n)=1
      end if
c
c  reset xmod, if zero, to number of actual model
c
      if(xmod.lt.1) xmod=nmod
c
      write(istdou,'(/'' Finish reading model with M ='',f10.5,
     *  '' Msun, R ='',1pe13.5,'' cm''/)') data(1)/amsun, data(2)
      if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *  write(istdpr,'(/'' Finish reading model with M ='',f10.5,
     *  '' Msun, R ='',1pe13.5,'' cm''/)') data(1)/amsun, data(2)
c
c  print model?
c
      if(nprmod.gt.0.and.istdpr.gt.0) then
        ndprmd=max0(1,nn/nprmod)
        write(istdpr,128) xmod,data,(n,x(n),(aa(i,n),i=1,6),
     *    n=1,nn,ndprmd)
      end if
c
      return
c
   90 write(istdou,150) nummod,ids,nrd
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,150) 
     *  nummod,ids,nrd
      icry=-1
   95 return
  105 format(//' Model contains g/(g tilde) in aa(6,.).'/
     *         ' Surface value =',1pe13.5/
     *         ' iturpr reset to 8')
  111 format(//' **** Warning in s/r readml. nn reduced from',i5,
     *  ' to',i5)
  112 format(//' **** error in s/r readml in taking every',i4,
     *  '-th point. nnr, nn =',2i5)
  115 format(4a8)
  120 format(//' name of model is: ',4a8)
  128 format(///' model no',f15.8//' data:',1p8e13.5//' n,x,aa(1-6):'
     *  //(i5,0pf10.7,1p6e13.5))
  150 format(//1x,10(1h*),' model no',i4,' not on d/s',i4,
     *  '. last model was no',i4)
      end
