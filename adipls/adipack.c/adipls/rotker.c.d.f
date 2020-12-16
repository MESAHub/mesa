      subroutine rotker(idsrkr,x,y,aa,el,sig,iy,ia,nn,nprtkr)
c
c  calculates and outputs rotational kernel and beta,
c  to unit idsrkr, if open
c 
c  Also possibly calculates rotational splitting for angular velocity
c  in common/comgrp/. This is flagged by irotkr (in common/coutpt/) .gt. 10
c  The rotational splitting and beta are stored in common/crot_split/
c  The splitting is given in terms of cyclic frequency, in microHz
c
c  Modified 13/8/87 to standardize output
c
c  Modified 21/7/94, adding unit number as argument
c
c  Modified 3/8/05, allowing for computation of rotational splitting
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
      parameter(iwork=10*nnmax)
      logical opnfil, nscfil
      dimension x(1),y(iy,1),aa(ia,1)
      dimension ea(2,3),v(2)
      dimension u(2,nnmax), h(nnmax)
      dimension del1(nnmax),del2(nnmax),del3(nnmax)
      dimension xis(2,nnmax)
      dimension del1as(nnmax), del2as(nnmax), del3as(nnmax)
      real mass, rad, gconst, oms2
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtk1,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
      common/csumma/ cs(50)
      common/comgrp/ isprtp, irotcp, omgrtp(1)
      common/crot_split/ beta, split(4)
      common/worksp/ w(4,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c   
      common/work/work(iwork) 

      save
c
      idiag=1
      pi=4.d0*atan(1.d0)
c
      iw=4

      mass = cs(2)/1000
      rad = cs(3)/100
      gconst = 6.67e-11

c
c Set conversion factor      
c
      oms2 = gconst*mass/(rad*rad*rad)
      oms=sqrt(oms2)

c
c  possibly return in radial case, if second-order effects are
c  not included
c
      if(el.lt.1.e-6) then
	if(irotkr.gt.10) then
	  beta=0.d0
	  split(1)=0.d0
	  if(irotkr.gt.20) then
	    del3(nn)=0.d0
            del2(nn)=0.d0
            go to 777
	  end if
        end if
        return
      end if
c
c  set integrands
c
      ell=el*(el+1)
      fct=1.-2./ell
      do 10 n=1,nn
      w(4,n)=x(n)*x(n)*aa(1,n)*aa(5,n)
      w(3,n)=y(2,n)*y(2,n)/ell
      w(1,n)=w(4,n)*(y(1,n)*y(1,n)+w(3,n))
      w(2,n)=y(1,n)-y(2,n)/ell
   10 w(2,n)=w(4,n)*(w(2,n)*w(2,n)+fct*w(3,n))
c  integrate
      do 20 k=1,2
      k2=k+2
   20 call vinta(x,w(k,1),w(k2,1),nn,iw,iw)
c  set beta and kernel
      beta=w(4,nn)/w(3,nn)
      fct=1./w(4,nn)
      do 30 n=1,nn
   30 w(1,n)=fct*w(2,n)
c  output
      if(istdpr.gt.0) write(istdpr,100) beta




      if(nprtkr.gt.1.and.istdpr.gt.0) then
        nd=max0(1,(nn-1)/(nprtkr-1))
        write(istdpr,110) (n,x(n),w(1,n),n=1,nn,nd)
      end if
c  file result
      cs(36)=beta
      if(opnfil(idsrkr).and.nscfil(idsrkr)) 
     *  write(idsrkr) cs,nn,(x(n),w(1,n),n=1,nn)

c
c  test for computing splitting
c
      if(irotkr.gt.10) then
	do n=1,nn
	  w(2,n)=w(1,n)*omgrtp(n)
        end do
        call vinta(x,w(2,1),w(3,1),nn,iw,iw)
	split(1)=1.d6*beta*w(3,nn)/(2.d0*pi)

        if(istdpr.gt.0) write(istdpr,*) 'split:'
        if(istdpr.gt.0) write(istdpr,*) split(1)

c
        if(idiag.gt.0)
     *    write(istdpr,'(//'' n, x, ker, omega, integral:''/
     *    (i5,0pf10.5,1p3e13.5))') 
     *    (n, x(n),w(1,n),omgrtp(n),w(3,n),n=1,nn)

      end if

c
c  when not including second-order terms, this is the end
c

      if(irotkr.le.20) return

c
c  Kara Burke second-order formulation
c
c call subroutine to calculate h(x)

      call uhx(x,u,aa,omgrtp,ia,nn)

c call subroutine to calculate xi_s and eta_s

      call sph(x,y,xis,aa,omgrtp,el,sig,iy,ia,nn)

c call subroutine to calculate 2nd order direct terms dependent on m^2
      if (el.gt.1.e-6) then
        call delta2(x,y,del2,del2as,xis,aa,omgrtp,sig,iy,ia,nn,el)
      
c call subroutine to calculate 2nd order distortion terms

        call delta3(x,y,del3,del3as,u,aa,omgrtp,sig,iy,ia,nn,el)
      end if

c
c  entry point for radial modes
c
c call subroutine to calculate 2nd order direct terms independent of m
c
  777 continue
      call delta1(x,y,del1,del1as,aa,omgrtp,sig,iy,ia,nn,el)
      
      freq = sqrt(sig)

c
c  no do loop here???
c..     do n=1,nn
       write(istdpr,*) '#D# omgrtp(nn) etc:', omgrtp(nn), freq,
     *  del1(nn), oms
       split(2) = omgrtp(nn)*omgrtp(nn)/freq*del1(nn)/oms
     *        + omgrtp(nn)*omgrtp(nn)/oms2*freq*del3(nn)*ell/(2.*el-1.)
     *        /(2.*el+3.)*oms

       split(2) = split(2)*1.d6/(2.d0*pi)


       split(3) = omgrtp(nn)*omgrtp(nn)/freq*del2(nn)/oms
     *        - omgrtp(nn)*omgrtp(nn)/oms2*freq*del3(nn)*3./(2.*el-1.)
     *        /(2.*el+3.)*oms

       split(3) = split(3)*1.d6/(2.d0*pi)
c..      end do

      if(istdpr.gt.0) write(istdpr,*) 'split2:'
      if(istdpr.gt.0) write(istdpr,*) split(2), split(3)


      if(idiag.gt.0) then

         freq = sqrt(sig)

c     OUTPUT 2ND ORDER ROTATIONAL SPLITTING RESULTS

         open(77,file="2ndorder")
         write(77,1042) el,cs(19),freq,1.-beta,del1(nn),
     *        del2(nn),del3(nn)

         open(78,file="2ndorderas")
         write(78,1042) el,cs(19),freq,1.-1.,del1as(nn),
     *        del2as(nn),del3as(nn) 

      end if

      return


12345 format(i5,2x,1pe22.15)

 3454 format(5(2x,1pe22.15))

  100 format(///' results on rotational splitting'//
     *  ' beta =',f12.7)
  110 format(/' n, x, kernel:'//(i5,0pf10.6,1pe13.5))
      
 1042 format(f4.1,1x,f5.1,1x,f10.2,1x,1p4e11.3,0p)

      
      end


