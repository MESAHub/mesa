      subroutine res_adimod(x_in, aa_in, data_in, nn_in, ivar_in, in, 
     *  data, x, aa, iaa_in, iaa, nn, icry)
c
c  Resets model provided as arguments in x_in, aa_in, data_in, for
c  pulsation calculation, possibly setting every in-th point. Also
c  resets number of points if required by Richardson extrapolation.
c
c  Actions correspond to those taken in s/r readml.
c
c  Original version: 14/2/06.
c
      implicit double precision (a-h,o-z)
      logical sincen, sinsur
      dimension data_in(8),x_in(1),aa_in(iaa_in,1),
     *  data(8),x(1),aa(iaa,nn)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/comgrp/ isprtp, irotcp, omgrtp(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data amsun /1.989d33/
c
      save
c
      if(istdpr.gt.0) write(istdpr,'(/'' Enter res_adimod''/)')
c
      icry = 0
c
      do i=1,8
	data(i)=data_in(i)
      end do
c
c  set flag for rotation or turbulent pressure setting
c
      idata8 = idint(data_in(8)+0.1)
c
c  test for inclusion of g/(g tilde)
c
      if(mod(idata8/10,10).eq.2) then
        iggt=1
        iturpr=8
        write(istdou,105) aa_in(6,nn_in)
        if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *    write(istdpr,105) aa_in(6,nn_in)
      else
        iggt=0
      end if
c
c  test for singular centre and/or surface
c
      sincen=x_in(1).eq.0
      sinsur=data_in(7).ge.0
      nsin=0
      if(sincen) nsin=nsin+1
      if(sinsur) nsin=nsin+1
c
c  Test for taking every in-th point in model
c  if iriche = 1, in addition need to check that number of 
c  non-singular points is odd
c
      if(in.le.1) then
c
c  Take every point
c
c  test for number of nonsingular points
c     
        if(iriche.ne.1.or.mod(nn_in-nsin,2).eq.1) then
          nshift=0
        else
          nshift=1
        end if
        if(nshift.ne.0) then
          write(istdou,111) nn_in,nn_in-nshift
          if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *      write(istdpr,111) nn_in,nn_in-nshift
        end if
        nn=nn_in-nshift
        if(sincen) then
          x(1)=x_in(1)
          do i=1,ivar_in
            aa(i,1)=aa_in(i,1)
	  end do
          do n=2,nn
            n1=n+nshift
            x(n)=x_in(n1)
	    if(isprtp.ne.0) omgrtp(n)=omgrtp(n1)
            do i=1,ivar_in
              aa(i,n)=aa_in(i,n1)
	    end do
	  end do
        else
          do n=1,nn
            if(n.eq.1) then
              n1=1
            else
              n1=n+nshift
            end if
            x(n)=x_in(n1)
	    if(isprtp.ne.0) omgrtp(n)=omgrtp(n1)
            do i=1,ivar_in
              aa(i,n)=aa_in(i,n1)
	    end do
	  end do
        end if
c
      else
c
c  take every in-th point
c
        if(sincen) then
          nn=(nn_in-2)/in+2
        else
          nn=(nn_in-1)/in+1
        end if
c
c  test for reducing for Richardson extrapolation
c
        if(iriche.eq.1.and.mod(nn-nsin,2).eq.0) then
          write(istdou,111) nn,nn-1
          if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *      write(istdpr,111) nn,nn-1
          nn=nn-1
        end if
c
        if(sincen) then
          x(1)=x_in(1)
          do i=1,ivar_in
            aa(i,1)=aa_in(i,1)
	  end do
          n1=nn_in-in*(nn-1)
          nstart=2
        else
          n1=nn_in-in*(nn+1)
          nstart=1
        end if
c
        do n=nstart,nn
          n1=n1+in
          x(n)=x_in(n1)
	  if(isprtp.ne.0) omgrtp(n)=omgrtp(n1)
          do i=1,ivar_in
            aa(i,n)=aa_in(i,n1)
	  end do
	end do
c
        if(n1.ne.nn_in) then
          write(istdou,112) in, nn_in, nn
          if(istdou.ne.istdpr.and.istdpr.gt.0)
     *      write(istdpr,112) in, nn_in, nn
          icry=-1
          go to 95
        end if
c
      end if
c
c  set g/gtilde (=1 in models without turbulent pressure)
c
      if(iturpr.eq.1) then
        do n=1,nn
          if(x(n).lt.0.999) then
            ggt=1
          else
            ggt=1./(x(n)*x(n)*x(n)*aa(1,n))
          end if
          aa(10,n)=ggt
	end do
      else if(iggt.eq.1) then
        do n=1,nn
          aa(10,n)=aa(6,n)
	end do
      else
        do n=1,nn
          aa(10,n)=1
	end do
      end if
c
   95 write(istdou,'(/'' Finish setting model with M ='',f10.5,
     *  '' Msun, R ='',1pe13.5,'' cm''/)') data(1)/amsun, data(2)
      if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *  write(istdpr,'(/'' Finish setting model with M ='',f10.5,
     *  '' Msun, R ='',1pe13.5,'' cm''/)') data(1)/amsun, data(2)
      return
c
  105 format(//' Model contains g/(g tilde) in aa(6,.).'/
     *         ' Surface value =',1pe13.5/
     *         ' iturpr reset to 8')
  111 format(//' **** Warning in s/r res_adimod. nn reduced from',i5,
     *  ' to',i5)
  112 format(//' **** error in s/r res_adimod in taking every',i4,
     *  '-th point. nn_in, nn =',2i5)
      end
