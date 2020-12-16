      subroutine instruct
c-----VERSION of October 5, 1995-----------------------------------------
c-----------------------------------------------------------------------

c          This subroutine contains instructions for using the subroutine
c     OPACGN93( z, xh, t6, R ) and OPAC(izi,mzin,xh,t6,r). 
c          The purpose of these subroutines is to perform 3 or 4 variable
c     interpolation on log10(kappa).  The opacity tables to be interpolated
c     are known to have somewhat random numerical errors of a few percent.
c     Consequently, adjusting the data prior to performing the interpolation is
c     justified at this level.  These codes are set-up to read the original
c     (unsmoothed) tabular data, this data is then passed through a smoothing
c     filter; using a set of routines developed by Mike Seaton (see M. J. Seaton
c     MNRAS 265,L25(1993)). It is the adjusted data that is actually used in
c     carrying out the interpolations in OPACGN93 and OPAC.  The interpolation
c     method, described below,  also produces smooth values for logK and its
c     first derivatives. The initial adjustment step helps improve the smoothness
c     of the OPACGn93 and OPAC  output,particularly at the smallest values of R.
c     The medium to large R output is only slightly effected by this step.  It
c     takes only a few seconds to carryout the initial data smoothing step. This
c     step can be skipped by setting the parameter ismdata in subroutine
c     readco =1.

c          The interpolation variables are :
c
c          z       The metallicity, Z
c          xh      The hydrogen mass fraction, X
c          t6      The temperature in millions of degrees Kelvin, T6
c          r       =rho(gm/cc)/T6**3, R
c
c          Additional input to OPAC is:

c          izi     Keeps or recalculates table indices. The value 0 causes
c                  the table indices to be recalculated.  A value other than 0
c                  causes the previous indices to be used.

c          mzin    The integer value of i of the Z value to use.  The 
c                  choices are:
c                  1=0.0  2=0.0001 3=0.0003 4=0.001 5=0.002 6=0.004 7=0.01
c                  8=0.02 9=0.03  10=0.04  11=0.06 12=0.08 13=0.1 
c
c
c          An interpolation between overlapping quadratics is used to obtain 
c     smoothed results.  A 4x4 grid in logT6 and logR is used to interpolate
c     in four different 3x3 sub-grids. Linear interpolation between quadratic
c     fits in these different sub-grids gives smoothed results in both log T6
c     and Log R directions. Compared to ealier versions of this code, the
c     interpolation in Z is in Kappa vs. Z; not log Kappa vs. Z.
c     The overlapping quadratic procedure produces results that are  smooth,
c     similar to bicubic spline interpolation, but require storage of only local
c     information.
c
c          The code OPACGN93 performs interpolation in Z, X, T6, and R. It calls
c     the subroutine OPAC at each Z.  If you are working with a fixed Z, as
c     listed above, it is more efficient to call OPAC directly. In this case use
c     izi=0 and the appropiate value of mzin. The  opacity data will be read
c     from unit 2 in the subroutine readco.  you will need to have the file
c     GN93hz available on disk.
c         Each of the individual tables in the GN93hz file cover 70 temperatures
c     in the range logT=3.75[T6=0.0056341325]( referred to as temperature 1) to
c     logT=8.7[T6=500]. and 19 values of log R in the range -8 (referred to as 1)
c     to +1; at half-integer steps.  (NOTE: earlier tables were explicitly in
c     terms of T6. For convenience the present tables tabulate log Kappa vs logT. The
c     interpolation however still uses T6 for the temperature variable)
c     For specialized problems, if storage space is a problem, a reduced set of
c     data can be input .  This requires a recompilation with altered parameter
c     values.  In order to limit the range of R, set the parameter nrb= index of
c     lowest value of log R to use(count from log R=-8).  Then set the parameter
c     nre to the index of the largest value of log R to use.  (NOTE: nre-nrb must
c     be at least 5). To ignore the first ntb-1 values of T6 (starting from
c     temperature 1) set the parameter ntb to desired value.  No other parameters
c     should be modified.

c     ***CAUTION***
c         As a result of the mixing procedure used to calculate the data a few
c     X=0.0, low T-small R, table values fell outside the range of T and R accessible
c     from the X=0.35 data directly calculated for this purpose.  These T-R 
c     locations are filled in with 9.99 (or for diagnostic purposes in some cases
c     larger values.  At the same locations the derivatives are set to 99.9.  When
c     T-R falls in this region a message is issued by the interpolation code to 
c     inform the user of this situation.  Presumable very few users will have
c     applications that take them into this region.
c     
c          Your routine that performs the call to OPAC should include the
c      statement:
c
c         common/e/ opact,dopact,dopacr,dopactd
c
c         These variables have the following meanings:
c        
c         OPACT        Is the Log of the Rosseland mean opacity: Log(kappa)  
c         DOPACT      Is Dlog(kappa)/Dlog(T6)   at constant R
c         DOPACR      Is Dlog(kappa)/Dlog(R)    at constant T
c         DOPACTD     Is Dlog(kappa)/Dlog(T6)   at constant Rho
c
      dum=0.0
      return
      end

C*********************************************************************
      subroutine opacgn93 (z,xh,t6,r,filename)
c.....The purpose of this subroutine is to interpolate the data along Z
      save
      character (len=*) :: filename
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/a/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)
c..... OPACT- opacity obtained from a quadraric interpolation at
c      fixed log T6 at three values of log R; followed by quadratic
c      interpolation along log T6. Results smoothed bt mixing
c      overlapping quadratics.
c..... DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
c              at fixed R  
c..... DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
c..... DOPACTD- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
c               at fixed rho
      common/e/ opact,dopact,dopacr,dopactd
      common/ee/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/eee/m1,zval
      real kapz,kapz1,kapz2
      dimension kapz(mz),dkapdtr(mz),dkapdrt(mz)
      zval=z
      zzl=z   ! use zzl=log10(.0001+z) for log interpolation
      if(itime .ne. 12345678) then
      do i=1,mz
      zza(i)=za(i)  ! use zza=log10(0.0001+za(i)) for log interpolation
      enddo
      endif
      if(itime .ne. 12345678) then
        itime=12345678
      endif

      do i=1,mz
        if(abs(z-za(i)) .lt. 1.e-7 ) then 
          izz=i
          call opac (0,izz,xh,t6,r,filename)
          if (opact .gt. 9.0) write (*,'(" logK > 9.0, X=",f7.5," Z=",
     x    f7.5," T6=",f10.5," R=",e12.4)') xh,z,t6,r
          return
        endif
      enddo

      ilo=2
      ihi=mz
    8 if(ihi-ilo .gt. 1) then
      imd=(ihi+ilo)/2
      if(z .le. za(imd)+1.e-7) then
        ihi=imd
      else
        ilo=imd
        endif
        go to 8
      endif
      i=ihi
      m1=i-2
      m2=i-1
      m3=i
      m4=i+1
      mfm=m4
c.....check whether Z is near a table limit 
      if((z .le. za(2)+1.e-7) .or. (z .ge. za(mz-1))) mfm=m3
c.....  Check if Z+X interpolation sums exceed unity at needed indices.
c       If so, backup to lower Z indices to perform interpolation.
c       This should work OK, due to density of Z-grid points and the 
c       slow Z variation(except at very small Z)
      if(xh+za(mfm) .gt. 1.) mfm=m3 
        if(xh+za(mfm) .gt. 1.) then
            if(m1 .le. 1) then
              write(*,'("special case: X,Z location not covered by"
     x                 ," logic")')
              stop
            endif 
          m1=m1-1
          m2=m2-1
          m3=m3-1
        endif
c
      izi=0
      do iz=m1,mfm
        izz=iz
        call opac(izi,izz,xh,t6,r,filename)
          if (opact .gt. 9.0) write (*,'(" logK > 9.0, X=",f7.5," Z=",
     x    f7.5," T6=",f10.5," R=",e12.4)') xh,z,t6,r
        izi=1
        kapz(iz)=10.**opact ! converts logK to K
        dkapdtr(iz)=dopact
        dkapdrt(iz)=dopacr
      enddo
      is=0
      iw=1
      kapz1=quad(is,iw,zzl,kapz(m1),kapz(m2),kapz(m3)
     x ,zza(m1),zza(m2),zza(m3))
      is=1
      dkapz1=quad(is,iw,zzl,dkapdtr(m1),dkapdtr(m2),dkapdtr(m3)
     x ,zza(m1),zza(m2),zza(m3))
      dkapz3=quad(is,iw,zzl,dkapdrt(m1),dkapdrt(m2),dkapdrt(m3)
     x ,zza(m1),zza(m2),zza(m3))
      if (mfm .eq. m3) then
        opact=log10(kapz1)   ! converts K to logK
        dopact=dkapz1
        dopacr=dkapz3
        dopactd=-3.*dopacr+dopact
        is=0

        return

      endif
      is=0
      iw=2
      kapz2=quad(is,iw,zzl,kapz(m2),kapz(m3),kapz(m4)
     x ,zza(m2),zza(m3),zza(m4))
      is=1
      dkapz2=quad(is,iw,zzl,dkapdtr(m2),dkapdtr(m3),dkapdtr(m4)
     x ,zza(m2),zza(m3),zza(m4))
      dkapz4=quad(is,iw,zzl,dkapdrt(m2),dkapdrt(m3),dkapdrt(m4)
     x ,zza(m2),zza(m3),zza(m4))
      dix=(zza(m3)-zzl)*dfsz(m3)
      opact=log10(kapz1*dix+kapz2*(1.-dix))   ! converts K to logK
      dopact=dkapz1*dix+dkapz2*(1.-dix)
      dopacr=dkapz3*dix+dkapz4*(1.-dix)
      dopactd=-3.*dopacr+dopact
      is=0
      return 
      end
c***********************************************************************
c
      subroutine opac (izi,mzin,xh,t6,r,filename)
c..... The purpose of this subroutine is to interpolate log kappa
c      in in X, T6, R
c        izi=0 recalulate table indices to use; =1 keep previous
c        mzin=index of za(i) in block data. za(i) are metallicities
c        t6=T6=temperature in millions of degrees kelvin
c        r=R=density(g/cm**3)/T6**3
c..... to use opac insert common/e/ in the calling routine.
c      This common contains interpolated values for kappa and its
c      first derivities.
c
      save
      integer w
      character (len=*) :: filename
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/aa/ q(4),h(4),xxh
      common/a/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)
      common/b/ itab(mx,mz),nta(nr),x(mx,mz),y(mx,mz),
     . zz(mx,mz)
      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/ee/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/eee/m1,zval
c..... OPACT- opacity obtained from a quadraric interpolation at
c      fixed log T6 at three values of log R; followed by quadratic
c      interpolation along log T6. Results smoothed bt mixing
c      overlapping quadratics.
c..... DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
c              at fixed R  
c..... DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
c..... DOPACTD- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics
c               at fixed rho
      common/e/ opact,dopact,dopacr,dopactd
c
      iop=1   ! provides smoothed interpolations; iop=0 gives no smoothing
      mzz=mzin
      z=za(mzz)

      if(nre .lt. 6) go to 65
c
      if((izi .eq. 0) .and. (z+xh-1.e-6 .gt. 1 )) go to 61
      if((izi .ne. 0) .and. (zval+xh-1.e-6 .gt. 1 )) go to 61
      xxh=xh
      xxi=xh
      t6i=t6
      ri=r
c
      xxx=log10(.005+xh)
      slt=log10(t6)
      slr=log10(r)
c
      if(itime .ne. 12345678) then
        itime=12345678
        do  i=1,mx
          xx(i)=log10(.005+xa(i))
        enddo
c..... this is the first time throught. Calculate the decadic
c      log of the perimeter points shifted by Z. m refers to
c      xa(m); the hydrogen table value.
 
c..... read the data files
        call readco(filename)
        xamx1=xa(mx-1)
        xxmx1=xx(mx-1)
        dfsxmx1=dfsx(mx-1)
      endif
      mxend=mx
      xa(mx)=1.-z
      xa(mx-1)=xamx1
      xx(mx-1)=xxmx1
      dfsx(mx-1)=dfsxmx1
        if (xa(mx) .lt. xa(mx-1)) then
          mxend=mx-1
          xa(mxend)=xa(mx)
        endif
        if (xh .ge. 0.8 ) then
          xx(mxend)=log10 (0.005+xa(mxend))
          dfsx(mxend)=1./(xx(mxend)-xx(mxend-1))
        endif
c
c
c..... Determine log R and log T6 grid points to use in the
c      interpolation.
      if((slt .lt. alt(1)).or.(slt .gt. alt(nt))) go to 62
      if((slr .lt. alr (1)).or.(slr .gt. alr(nre))) go to 62
c
c
c
      if (izi .eq. 0) then  ! freeze table indices
        ilo=2
        ihi=mx
    8   if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(xh .le. xa(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 8
        endif
        i=ihi
        mf=i-2
        mg=i-1
        mh=i
        mi=i+1
        mf2=mi
        if (xh .lt. 1.e-6) then
        mh=1
        mg=1
        mi=2
        mf2=1
        endif
        if((xh .le. xa(2)+1.e-7) .or. (xh .ge. xa(mx-2)-1.e-7)) mf2=mh
c
        ilo=2
        ihi=nre
   12     if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(slr .le. alr(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 12
          endif
        i=ihi
        l1=i-2
        l2=i-1
        l3=i
        l4=l3+1
c
        ilo=2
        ihi=nt
   11     if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(t6 .le. t6list(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 11
          endif
        i=ihi
        k1=i-2
        k2=i-1
        k3=i
        k4=k3+1
        l3s=l3+nrb-1
        k3s=k3+ntb-1
      endif

      kmin=0
      k1in=k1
      iadvance=0
      mfin=mf
      if ((mfin .eq. 1) .and. (xz(1,mzz,k1,l1) .gt. 9.)) then   ! no data
      do i=1,6
        if (xz(1,mzz,i,l1) .gt. 9.)  then
          if (xh .lt. .1) then
           kmin=i+1
          else

            if (iadvance .eq. 0) then
            iadvance=iadvance+1
            mf=mf+1 
            mg=mg+1
            mh=mh+1
            mi=mi+1
            mf2=mf2+1
            endif
          endif
        endif
      enddo
      if ((iadvance .eq. 0) .and. (k1 .le. kmin) .and.
     x    (slt .le. alt(kmin))) then
      k1=kmin
      if ((xz(1,mzz,kmin,l1+1) .lt. 9.) .and.
     x   ((slr+.01) .gt. alr(l1+1))) then
      l1=l1+1
      kmin=0
      k1=k1in
      do i=1,6
      if (xz(1,mzz,i,l1) .gt. 9.) kmin=i+1
      enddo
      if ((kmin .ne. 0) .and. (k1in .lt. kmin)) k1=kmin
      endif
      endif
      if ((slt+.001) .lt. alt(k1)) then
      opact=30.
      dopact=99.
      dopacr=99.
      return
      endif
      l2=l1+1
      l3=l2+1
      l4=l3+1
      l3s=l3+nrb-1
      k2=k1+1
      k3=k2+1
      k4=k3+1
      k3s=k3+ntb-1
      endif
        do i=14,18   ! allows jagged edge at high T,rho
          if((l3s .gt. i) .and. (k3s .gt. nta(i+1))) go to 62
        enddo
      do 123 m=mf,mf2
      ip=3
      iq=3
      ntlimit=nta(l3s)
      if((k3. eq. ntlimit) .or. (iop .eq. 0)) then 
        ip=2
        iq=2
      endif
      if(t6 .le. t6list(2)+1.e-7) ip=2

      if((l3 .eq. nre) .or. (iop .eq. 0)) then 
       iq=2
       ip=2
      endif
      if ((l4 .le.nr) .and. (xz(m,mzz,k3,l4) .eq. .0)) iq=2
      if(slr .le. alr(2)+1.e-7) iq=2
c
      is=0
c
c__________
      do ir=l1,l1+iq
        do it=k1,k1+ip
        opl(m,it,ir)=xz(m,mzz,it,ir)
        is=1
        enddo
      enddo
  123 continue
      if((zz(mg,mzin) .ne. zz(mf,mzin)) .or.
     x    (zz(mh,mzin) .ne. zz(mf,mzin))) then
      write(*,'("Z does not match Z in GN93hz files you are"
     x ," using")')
      stop
      endif
      if(z .ne. zz(mf,mzin)) go to 66
c                  with return
      is=0
      iw=1
      do 45 ir=l1,l1+iq
        do it=k1,k1+ip
        if (mf2 .eq. 1) then
        opk(it,ir)=opl(mf,it,ir)
        go to 46
        endif
        opk(it,ir)=quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir)
     x  ,opl(mh,it,ir),xx(mf),xx(mg),xx(mh))
        is=1
   46   continue
        enddo
   45 continue

      if (mi .eq. mf2) then  ! interpolate between quadratics
      is=0
      iw=1
       dixr=(xx(mh)-xxx)*dfsx(mh)
      do 47 ir=l1,l1+iq
        do it=k1,k1+ip
        opk2(it,ir)=quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir)
     x  ,opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
        opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.-dixr)
        is=1
        enddo
   47 continue
c     interpolate X between two overlapping quadratics
      endif

      is=0
c
c..... completed H,Z interpolation. Now interpolate T6 and log R on a
c      4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure
c      mixes overlapping quadratics to obtain smoothed derivatives.
c
c
      call t6rinterp(slr,slt)
      return

   61 write(*,'(" Mass fractions exceed unity")')
c                  with a return
      stop
   62 write(*,'(" T6/LogR outside of table range")')
c     write(*,'("slt,alt(1),alt(nt),slr,alr(1),alr(nre),l3s,i,k3s,
c    x nta(i+1)",6e12.5,4i5)') slt,alt(1),alt(nt),slr,alr(1),
c    x alr(nre),l3s,i,k3s,nta(i+1)
c                  with a return 
      stop
   64 write(*,'(" X not equal to zero: To run this case it
     .is necessary"/ "to recompile with parameter (mx=1)")')
      stop
   65 write(*,'("Too few R values; NRE+1-NRB < 6")')
      return
   66 write(*,'(" Z does not match Z in codata* files you are",
     . " using")')
      stop
      end

c***********************************************************************
      subroutine t6rinterp(slr,slt)
c     The purpose of this subroutine is to interpolate in logT6 and logR
      save
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/ee/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/aa/ q(4),h(4),xxh
      common/a/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)
      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/e/ opact,dopact,dopacr,dopactd
c
      iu=0
      is=0
      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3),
     x  alr(l1),alr(l2),alr(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quad(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4),
     x      alr(l2),alr(l3),alr(l4))
          endif
        is=1
      enddo
c
      is=0
      iw=1
c..... k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opact=quad(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
      dopact=dkap
      dkap1=dkap
        if (iq. eq. 3) then
c.....    k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          opactq=quad(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3))
          dkapq1=dkap
        endif
        if(ip .eq. 3) then
c.....    k and Dlog(k)/Dlog(T6) in lower-left 3x3.
          opact2=quad(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4))
          dkap2=dkap
c.....    k and Dlog(k)/Dlog(T6) smoothed in left 3x4
          dix=(alt(k3)-slt)*dfs(k3)
          dopact=dkap1*dix+dkap2*(1.-dix)
          opact=opact*dix+opact2*(1.-dix)
        endif
        if(iq .eq. 3) then
 
c.....    k and Dlog(k)/Dlog(T6) in upper-right 3x3.
          opactq2=quad(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
          dkapq2=dkap
          dopactq=dkapq1*dix+dkapq2*(1.-dix)
          opactq=opactq*dix+opactq2*(1.-dix)
        endif
c
      iu=0
      do lx=l1,l1+iq
        iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx),
     x  alt(k1),alt(k2),alt(k3))
          if(ip .eq. 3) then
            iw=2
            q(iu)=quad(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx),
     x      alt(k2),alt(k3),alt(k4))
          endif
        is=1
      enddo
c
      is=0
      iw=1
c..... k and Dlog(k)/Dlog(R) in lower-left 3x3
      opacr=quad(is,iw,slr,h(1),h(2),h(3),alr(l1),alr(l2),alr(l3))
      dopacr=dkap
        if(ip .eq. 3) then
          opacrq=quad(is,iw,slr,q(1),q(2),q(3),alr(l1),alr(l2),alr(l3))
c.....    k and Dlog(k)/Dlog(R) in upper-left 3x3.
          dkapq3=dkap
        endif
        if(iq .eq. 3) then
c.....    k and Dlog(k)/Dlog(R) in lower-right 3x3.
          opact2=quad(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4))
          dix2=(alr(l3)-slr)*dfsr(l3)
          dopacr=dopacr*dix2+dkap*(1.-dix2)
            if(ip .eq. 3) then
c.....        k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
              dopact=dopact*dix2+dopactq*(1.-dix2)
              opact=opact*dix2+opactq*(1-dix2)
            endif
         endif
        if(ip .eq. 3) then
c.....    k and Dlog(k)/Dlog(R) in upper-right 3x3.
          opacrq=quad(is,iw,slr,q(2),q(3),q(4),alr(l2),alr(l3),alr(l4))
            if(iq .eq. 3) then
c.....        Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
              dopacrq=dkapq3*dix2+dkap*(1.-dix2)
              dopacr=dopacr*dix+dopacrq*(1.-dix)
            endif
        endif
      dopactd=dopact-3.*dopacr
        if (opact .gt. 1.e+15) then
          write(*,'("Interpolation indices out of range",
     x              ";please report conditions.")') 
          stop
        endif
          if (opact .gt. 9.) then
            dopact=99.
            dopacr=99.
            dopactd=99.
          endif
      return
      end

c
c**********************************************************************
      subroutine readco(filename)
c..... The purpose of this subroutine is to read the data tables
      save
      parameter (ismdata=0)   ! modified
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      character*1 dumarra(250)
      character (len=*) :: filename
      common/aa/ q(4),h(4),xxh
      common/a/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)
      common/b/ itab(mx,mz),nta(nr),x(mx,mz),y(mx,mz),
     . zz(mx,mz)
      common/e/ opact,dopact,dopacr,dopactd
      common/ee/ opl(mx,nt,nr),xx(mx),zza(mz)
      common/alink/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),xzff(100,nr)  
      COMMON/CST/NRL,RLS,nset,tmax  ! modified



        if (itimeco .ne. 12345678) then
        do i=1,mx
          do j=1,mz 
            do k=1,nt
              do l=1,nr
                xz(i,j,k,l)=1.e+35
              enddo
            enddo
          enddo
        enddo
        itimeco=12345678
        endif
c
      close (2)
c..... read  tables
       open(2, FILE=trim(filename))
c
c      read header
      do i=1,240
      read (2,'(a)') dumarra(i)
      enddo
c
      do 3 m=1,mx
      do 2 i=1,n(m)
c
      read(2,'(f10.5)') dum
      read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4)')
     .itab(m,i),x(m,i),y(m,i),zz(m,i)
      read(2,'(f10.5)') dum,dum,dum
      read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
      read(2,'(f10.5)') dum
        do k=1,ntm
        read(2,'(f4.2,19f7.3)') alt(k),(xzf(k,l), l=1,nrm)
        alt(k)=alt(k)-6.
        if (isett6 .ne. 1234567) then
        t6listf(k)=10.**alt(k)
        t6arr(k)=t6listf(k)
        endif 
          do ll=1,nrm   ! modified
          xzff(k,ll)=xzf(k,ll)
          enddo
        enddo
        isett6=1234567

       if (ismdata .eq. 0) then
        tmax=10.   ! modified
        nset=65
        RLS=-8.
        nsm=1 
          RLE=1.
          nrlow=1
          nrhigh=2*(RLE-RLS)+1

        call opaltab    !modified
       endif

 1010  continue

      ll=1
      do 110 kk=1,nre
      alr(ll)=alrf(kk)
        do k=1,nt
        t6list(k)=t6listf(k+ntb-1)
        if(ismdata .eq. 0) then
c           Following skip required because, due to missing data,
c           the X=0  low T data cannot be smoothed
          if ((m  .eq. 1) .and. (k .le. 9)) then
            xz(m,i,k,ll)=xzf(k+ntb-1,kk)
          else
            xz(m,i,k,ll)=xzff(k+ntb-1,kk)
          endif
        else
         xz(m,i,k,ll)=xzf(k+ntb-1,kk)
        endif
        enddo
  110 ll=ll+1

    2 continue
    3 continue
 
      do 12 i=2,nt
   12 dfs(i)=1./(alt(i)-alt(i-1))
      do 13 i=2,nr
   13 dfsr(i)=1./(alr(i)-alr(i-1))
      do i=2,mx-1
      dfsx(i)=1./(xx(i)-xx(i-1))
      enddo
      do i=2,mz
      dfsz(i)=1./(zza(i)-zza(i-1))
      enddo
      return
      end
c
c***********************************************************************
      function quad(ic,i,x,y1,y2,y3,x1,x2,x3)
c..... this function performs a quadratic interpolation.
      save
      common/d/dkap
      dimension  xx(3),yy(3),xx12(30),xx13(30),xx23(30),xx1sq(30)
     . ,xx1pxx2(30)
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      yy(1)=y1
      yy(2)=y2
      yy(3)=y3
        if(ic .eq. 0) then
          xx12(i)=1./(xx(1)-xx(2))
          xx13(i)=1./(xx(1)-xx(3))
          xx23(i)=1./(xx(2)-xx(3))
          xx1sq(i)=xx(1)*xx(1)
          xx1pxx2(i)=xx(1)+xx(2)
        endif
      c3=(yy(1)-yy(2))*xx12(i)
      c3=c3-(yy(2)-yy(3))*xx23(i)
      c3=c3*xx13(i)
      c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3
      c1=yy(1)-xx(1)*c2-xx1sq(i)*c3
      dkap=c2+(x+x)*c3
      quad=c1+x*(c2+x*c3)
      return
      end
c
c***********************************************************************
      block data
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/aa/ q(4),h(4),xxh
      common/a/ mzz, xz(mx,mz,nt,nr),  
     . t6list(nt),alr(nr),n(mx),alt(nt),opk(nt,nr),opk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),dfsz(mz),a(3,mx),b(3),m,mf,xa(mx)
     . ,alrf(nrm),xzf(nt,nr),t6listf(ntm),za(mz)
      common/b/ itab(mx,mz),nta(nr),x(mx,mz),y(mx,mz),
     . zz(mx,mz)
      data (xa(i),i=1,mx-1)/0.0,0.1,0.2,0.35,0.5,.7,.8,.9,.95/
      data (za(i),i=1,mz)/.0,0.0001,.0003,.001,.002,.004,.01,.02,.03,
     x .04,.06,.08,.1/
      data (nta(i),i=1,nrm)/14*70,69,64,60,58,57/
      data (n(i),i=1,mx)/13,13,13,13,13,13,13,13,10,12/
      end
c
c***********************************************************************
      subroutine opaltab
C
C  CODE FOR FITTING AND SMOOTHING OPAL DATA. ADAPTED FROM A CODE
C     WRITTEN BY MIKE SEATON(obtained june 1993)
C
C     OPAL DATA.
C     ASSUMES FIRST T6=0.006, LAST T6=10.OR 0.04). Depending on position
C     in the table. 
C     USES RECTANGULAR ARRAY FOR VARIABLES T6 AND LOG10(R)
C
C     (1) NSM=NUMBER OF PASSES THROUGH SMOOTHING FILTER.
C     USE OF NSM=1 OR 2 IS RECOMMENDED.
C     NO SMOOTHING WITH NSM=0
C     (2) RANGE FOR LOG10(R),
C     RLS=FIRST VALUE, RLE=LAST VALE
C     (RLS MUST BE FIRST VALUYE IN TABLE)
C
C  SUBROUTINE INTERP
C     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF
C               SUBROUTINE INTERP
C     WHICH GIVES LOG(ROSS) AND TWO FIRST DERIVATIVES FOR ANY
C     VALUES OF LOG(T) AND LOG(RHO). SEE BELOW FOR FURTHER
C     EXPLANATION.
C
C  OUTPUT FOR THE CASE OF NSM.GT.0.
C     INTERP IS USED TO OBTAIN SMOOTHED DATA INTERPOLATED
C     BACK TO THE ORIGINAL OPAL MESH. TWO FILES ARE WRITTEN.
C
C
C  THE SUBROUTINES SPLINE AND SPLINT ARE ADAPTED FROM THOSE GIVE BY
C  W.H. Press, S.A. Teulolsky, W.T. Vettering and B.P. Flannery,
C  "Numerical Recipes in FORTRAN", 2nd edn., 1992, C.U.P.
C  OTHER REFERENCES ARE MADE TO METHODS DESCRIBED IN THAT BOOK.
C
      PARAMETER(IP=100,IPR=20)
      parameter (mx=10,mz=13,nrm=19,nrb=1,nre=19,nr=nrm+1-nrb
     . ,ntm=70,ntb=1,nt=ntm+1-ntb)
      DIMENSION U(IP),ROSSL(IP,IPR),V(IP),V2(IP)
      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      CHARACTER*1 HEAD(100)
      COMMON/CST/NRL,RLS,nset,tmax  ! modified
      common/alink/ N,NSM,nrlow,nrhigh,RLE,t6arr(100),xzff(100,nr)  
      LOGICAL IERR

C
      NRL=2*(RLE-RLS)+1
C
C     STORE LOG10(T) IN U AND LOG10(ROSS) IN ROSSL
C     CHECK FIRST VALUE OF T6
      T6=t6arr(1)
      do j=1,NRL
      ROSSL(1,j)=xzff(1,j)
      enddo

      if (abs(T6-.0056341325) .lt. 1.e-8) then
         U(1)=6.+LOG10(T6)
      ENDIF
C     SET ROSSL UP TO T6=t6arr(nset)
      I=1
    5 I=I+1
      T6=t6arr(I)
      do j=1,NRL
      ROSSL(I,j)=xzff(I,j)
      enddo
         U(I)=6+LOG10(T6)
         IF(T6.LT.tmax)GOTO 5
      N=I
      IF(N.GT.IP)THEN
         PRINT*,' REQUIRE PARAMETER IP OF AT LEAST ',N
         STOP
      ENDIF
C
C
C     DEFINE VARIABLES
C         X=20.0*(LOG10(T)-3.80)+1
C         Y=2.0*(LOG10(R)-RLS)+1
C     USE INDICES I=1 TO nset AND J=1 TO NRL
C     X AND Y ARE SUCH THAT, ON MESH-POINT (I,J), X=I AND Y=J
C     OBTAIN:-
C         F(I,J)=LOG10(ROSS)
C         FX(I,J)=dF/dX
C         FY(I,J)=dF/dY
C         FXY(I,J)=ddF/dXdY
C
C
C     FIRST GET F AND FX, INTERPOLATING FROM OPAL T6 TO
C     INTERVAL OF 0.05 IN LOG10(T).
      DO 40 J=1,NRL
C        FOR EACH LOG10(R), STORE LOG10(ROSS) IN V(I)
         DO 20 I=1,N
            V(I)=ROSSL(I,J)
   20    CONTINUE
C
C        GET FIRST DERIVATIVES AT END POINTS
C
C        GET SECOND DERIVATIVES FOR SPLINE FIT
         CALL SPLINE(U,V,N,V2)
C
C        INTERPOLATE TO LOG10(T)=FLT, FLT=3.8(0.05)8.0
         DO 30 I=1,nset ! modified
            FLT=3.75+0.05*I
            CALL SPLINT(U,V,N,V2,FLT,F(I,J),FX(I,J))
   30    CONTINUE
C
   40 CONTINUE
C
C
C  OPTION FOR SMOOTHING
      IF(NSM.GT.0)THEN
         DO 35 NS=1,NSM
            CALL SMOOTH
   35    CONTINUE
         CALL FITX
      ENDIF
C
C
C  GET FY AND FXY
      CALL FITY
C
C  THE ARRAYS F, FX, FY AND FXY ARE NOW STORED
C
C  CAN NOW DO INTERPOLATIONS USING
C       CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
C       INPUT IS FLT=LOG10(T), FLRHO=LOG10(RHO)
C       OUTPUT IS G=LOG10(ROSS)
C              DGDT=dG/d(LOG10(T))
C            DGDRHO=dG/d(LOG10(RHO))
C              IERR=.TRUE. IF INPUT FLT, FLRHO ARE OUT-OF-RANGE,
C                          ELSE IERR=.FALSE.
C
C INTERPOLATE BACK TO OPAL POINTS
      IF(NSM.GT.0)THEN
         do l=1,NRL
         xzff(1,l)=ROSSL(1,l)
         enddo

         DO 70 K=2,N
            FLT=U(K)
            DO 50 L=nrlow,nrhigh
               FLR=RLS+.5*(L-1)
               FLRHO=FLR-18.+3.*FLT
               CALL INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
               IF(IERR)THEN
               ENDIF
               V(L)=G
   50       CONTINUE
            T6=t6arr(K)
            do l=nrlow,nrhigh
            xzff(K,l)=V(l)

            enddo

   70    CONTINUE
      ENDIF
C
C
 1000 FORMAT('  SMOOTHED OPAL DATA')
 1100 FORMAT('  OPAL DATA, (SMOOTHED-ORIGINAL)')
 2000 FORMAT(100A1)
 2222 FORMAT(F8.3,20F7.3)
 6000 FORMAT(/' FIRST T6=',1P,E10.3,', SHOULD BE 0.006')
 6003 FORMAT(/' !!! OUT-OF-RANGE !!!'/' FLT=',1P,E10.3,', FLRHO=',E10.3,
     + ', FLR=',E10.3)
C
      END
C
C*********************************************************************
      SUBROUTINE SPLINE(X,Y,N,Y2)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
C
C     FIRST DERIVATIVES AT END POINTS USING CUBIC FIT
         YP1=((Y(3)-Y(1))*(X(2)-X(1))**2
     +   -(Y(2)-Y(1))*(X(3)-X(1))**2)/
     +   ((X(3)-X(1))*(X(2)-X(1))*(X(2)-X(3)))
         YPN=((Y(N-2)-Y(N))*(X(N-1)-X(N))**2
     +   -(Y(N-1)-Y(N))*(X(N-2)-X(N))**2)/
     +   ((X(N-2)-X(N))*(X(N-1)-X(N))*(X(N-1)-X(N-2)))
C
      Y2(1)=-0.5
      U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      QN=0.5
      UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
C*********************************************************************
      SUBROUTINE SPLINT(XA,YA,N,Y2A,X,Y,YP)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      YP=0.05*  (  (-YA(KLO)+YA(KHI))/H
     +   +      ( -(3*A**2-1)*Y2A(KLO)
     +            +(3*B**2-1)*Y2A(KHI) )*H/6. )
      RETURN
      END
C*********************************************************************
      SUBROUTINE FITY
C
C  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS
C  FY AND FXY
C
C
      COMMON/CST/NRL,RLS,nset,tmax  ! modified
C
      PARAMETER(IPR=20)
      DIMENSION A(IPR),B(IPR),AD(IPR),BD(IPR)
      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
C
      DO 30 I=1,nset   ! modified
         DO 10 J=1,NRL
            A(J)=F(I,J)
            B(J)=FX(I,J)
   10    CONTINUE
C
         CALL GETD(A,NRL,AD,AP1,APN)
         CALL GETD(B,NRL,BD,BP1,BPN)
C
         FY(I,1)=AP1
         FY(I,NRL)=APN
         FXY(I,1)=BP1
         FXY(I,NRL)=BPN
         DO 20 J=2,NRL-1
            FY(I,J)= -A(J)+A(J+1)-2.*AD(J)-AD(J+1)
            FXY(I,J)=-B(J)+B(J+1)-2.*BD(J)-BD(J+1)
   20    CONTINUE
   30 CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE FITX
C
C  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.
C  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.
C
C
      PARAMETER(IPR=20)
      DIMENSION A(85),D(85)
C
      COMMON/CST/NRL,RLS,nset,tmax  ! modified
      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
C
      DO 30 J=1,NRL
         DO 10 I=1,nset ! modified
            A(I)=F(I,J)
   10    CONTINUE
         CALL GETD(A,nset,D,AP1,APN)  ! modified
         FX(1,J)=AP1
         FX(nset,J)=APN   ! modified
         DO 20 I=2,nset-1  ! modified
            FX(I,J)=-A(I)+A(I+1)-2.*D(I)-D(I+1)
   20    CONTINUE
   30 CONTINUE
C
      RETURN
      END
C
C************************************************************************
      SUBROUTINE GETD(F,N,D,FP1,FPN)
C
C  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS
C  OF UNITY.
C
C
      DIMENSION F(N),D(N),T(85)
C
      FP1=(-11.*F(1)+18.*F(2)-9.*F(3)+2.*F(4))/6.
      FPN=(11.*F(N)-18.*F(N-1)+9.*F(N-2)-2.*F(N-3))/6.
C
      D(1)=-.5
      T(1)=.5*(-F(1)+F(2)-FP1)
C
      DO 10 J=2,N-1
         D(J)=-1./(4.+D(J-1))
         T(J)=-D(J)*(F(J-1)-2.*F(J)+F(J+1)-T(J-1))
   10 CONTINUE
C
      D(N)=(FPN+F(N-1)-F(N)-T(N-1))/(2.+D(N-1))
C
      DO 20 J=N-1,1,-1
         D(J)=D(J)*D(J+1)+T(J)
   20 CONTINUE
C
      RETURN
      END
C
C*********************************************************************
      SUBROUTINE INTERP(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
C
C  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE
C  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN
C  Numerical Recipes, PP. 118 TO 120
C
C
      PARAMETER(IPR=20)
      COMMON/CF/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      DIMENSION B(16)
      LOGICAL IERR
C
      COMMON/CST/NRL,RLS,nset,tmax  ! modified
C
C  EXTREME LIMITS ALLOWED ARE:-
C     (3.800-0.0125) TO (8.000+0.0125) FOR LOG10(T)
C     (RLS-0.125) TO (RLE+0.1254) FOR LOG10(R)
C     (ALLOWING FOR SMALL EXTRAPOLATIONS BEYOND TABULAR VALUES)
C
C  FUNCTION DEFINITIONS FOR CUBIC EXPANSION
C
      FF(S,T)=    B( 1)+T*(B( 2)+T*(B( 3)+T*B( 4)))
     +   +S*(     B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     +   +S*(     B( 9)+T*(B(10)+T*(B(11)+T*B(12)))
     +   +S*(     B(13)+T*(B(14)+T*(B(15)+T*B(16))) )))
C
      FFX(S,T)=   B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     +   +S*(  2*(B( 9)+T*(B(10)+T*(B(11)+T*B(12))))
     +   +S*(  3*(B(13)+T*(B(14)+T*(B(15)+T*B(16)))) ))
C
      FFY(S,T)=   B( 2)+S*(B( 6)+S*(B(10)+S*B(14)))
     +   +T*(  2*(B( 3)+S*(B( 7)+S*(B(11)+S*B(15))))
     +   +T*(  3*(B( 4)+S*(B( 8)+S*(B(12)+S*B(16)))) ))
C
      FFXY(S,T)=  B( 6)+T*(2*B( 7)+3*T*B( 8))
     +   +S*(   2*B(10)+T*(4*B(11)+6*T*B(12))
     +   +S*(   3*B(14)+T*(6*B(15)+9*T*B(16)) ))
C
C
      IERR=.FALSE.
C
      X=20.*(FLT-3.800)+1
      FLR=FLRHO+18.-3.*FLT
      Y=2*( FLR - RLS )+1
C
      IF(X.LT.2.)THEN
         IF(X.LT.0.75)THEN
            IERR=.TRUE.
         ELSE
            I=1
         ENDIF
      ELSEIF(X.GT.84)THEN
         IF(X.GT.85.25)THEN
            IERR=.TRUE.
         ELSE
            I=84
         ENDIF
      ELSE
         I=X
      ENDIF
      U=X-I
C
      IF(Y.LT.2.)THEN
         IF(Y.LT.0.75)THEN
            IERR=.TRUE.
         ELSE
            J=1
         ENDIF
      ELSEIF(Y.GT.NRL-1)THEN
         IF(Y.GT.NRL+.25)THEN
            IERR=.TRUE.
         ELSE
            J=NRL-1
         ENDIF
      ELSE
         J=Y
      ENDIF
      V=Y-J
C
      IF(IERR)THEN
         G=9.999
         DGDT=9.999
         DGDRHO=9.999
         RETURN
      ENDIF
C
C
C  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
      B(1)=F(I,J)
      B(2)=FY(I,J)
      B(3)=3*(-F(I,J)+F(I,J+1))-2*FY(I,J)-FY(I,J+1)
      B(4)=2*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1)
C
      B(5)=FX(I,J)
      B(6)=FXY(I,J)
      B(7)=3*(-FX(I,J)+FX(I,J+1))-2*FXY(I,J)-FXY(I,J+1)
      B(8)=2*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1)
C
      B(9)=3*(-F(I,J)+F(I+1,J))-2*FX(I,J)-FX(I+1,J)
      B(10)=3*(-FY(I,J)+FY(I+1,J))-2*FXY(I,J)-FXY(I+1,J)
      B(11)=9*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     + +6*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))
     + +4*FXY(I,J)
     + +3*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))
     + +2*(FXY(I,J+1)+FXY(I+1,J))
     + +FXY(I+1,J+1)
      B(12)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     + +4*(-FX(I,J)+FX(I,J+1))
     + +3*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))
     + +2*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))
     + -FXY(I+1,J)-FXY(I+1,J+1)
C
      B(13)=2*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J)
      B(14)=2*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J)
      B(15)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     + +4*(-FY(I,J)+FY(I+1,J))
     + +3*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))
     + +2*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))
     + -FXY(I+1,J+1)-FXY(I,J+1)
      B(16)=4*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     + +2*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)
     +    +FY(I,J)-FY(I+1,J)-FY(I+1,J+1)+FY(I,J+1))
     + +FXY(I,J)+FXY(I+1,J)+FXY(I+1,J+1)+FXY(I,J+1)
C
C  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),
C      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)
      G=FF(U,V)
      DGDT=20.*FFX(U,V)-6.*FFY(U,V)
      DGDRHO=2.*FFY(U,V)
C
C
      RETURN
      END

C*********************************************************************
      SUBROUTINE SMOOTH
C
C  THIS SUBROUTINE USES A 2-DIMENSIONAL GENERALISATION OF THE SMOOTHING
C  TECHNIQUES DESCRIBED ON PP. 644 TO 649 OF Numerical Recipes.
C
C  CONSIDER THE 25 POINTS DEFINED BY
C       I+n, n=-2,-1,0,1,2 AND J+m, m=-2,-1,0,1,2.
C  THE FUNCTION TO BE SMOOTHED IS FITTED TO A BI-CUBIC, INVOLVING
C  16 COEFFICIENTS, USING TECHNIQUES OF LEAST-SQUARES. THE SMOOTHED
C  FUNCTION (TEMPORARILY STORED IN FXY) IS GIVEN BY THE FITTED VALUE
C  AT THE POINT I AND J.
C
C  THE FITTING IS SHIFTED FOR POINTS CLOSE TO BOUNDARIES.
C
C
      PARAMETER(IPR=20)
C
      COMMON/CST/NRL,RLS,nset,tmax  ! modified
      COMMON/CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
C
      DIMENSION GAM(6)
      DATA GAM/+0.0073469388,-0.0293877551,-0.0416326531,
     +         +0.1175510204,+0.1665306122,+0.2359183673/
      DIMENSION BET(11)
      DATA BET/
     + -0.0048979592,-0.0661224490,-0.0293877551,+0.0195918367,
     +  0.2644897959,+0.1175510204,-0.0783673469,+0.0277551020,
     +  0.3746938776,+0.1665306122,-0.1110204082/
      DIMENSION ALP(11)
      DATA ALP/
     + -0.0844897959,-0.0048979592,+0.0073469388,+0.0012244898,
     +  0.3379591837,+0.0195918367,-0.0293877551,+0.4787755102,
     +  0.0277551020,-0.0416326531,-0.0069387755/
C
C
      DO 20 I=3,nset-2 
C
         J=1
         FXY(I,J)=
     +    ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     +   +ALP(2)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J+3)+F(I+2,J+3)
     +          +F(I-1,J+4)+F(I+1,J+4) )
     +   +ALP(3)*( F(I-2,J+2)+F(I+2,J+2) )
     +   +ALP(4)*( F(I-2,J+4)+F(I+2,J+4) )
     +   +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     +   +ALP(6)*( F(I-1,J+1)+F(I+1,J+1)+F(I-1,J+3)+F(I+1,J+3) )
     +   +ALP(7)*( F(I-1,J+2)+F(I+1,J+2) )
     +   +ALP(8)*  F(I  ,J  )
     +   +ALP(9)*( F(I  ,J+1)+F(I  ,J+3) )
     +   +ALP(10)* F(I  ,J+2) +ALP(11)*F(I  ,J+4)
C
         J=2
         FXY(I,J)=
     +    BET(1)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J+3)+F(I+2,J+3) )
     +   +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     +   +BET(3)*( F(I-2,J+1)+F(I+2,J+1) )
     +   +BET(4)*( F(I-2,J+2)+F(I+2,J+2)+F(I-1,J-1)+F(I+1,J-1)
     +            +F(I-1,J+3)+F(I+1,J+3) )
     +   +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     +   +BET(6)*( F(I-1,J+1)+F(I+1,J+1) )
     +   +BET(7)*( F(I-1,J+2)+F(I+1,J+2) )
     +   +BET(8)*( F(I  ,J-1)+F(I  ,J+3) )
     +   +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J+1) +BET(11)*F(I  ,J+2)
C
         DO 10 J=3,NRL-2
            FXY(I,J)=
     +         GAM(1)*( F(I-2,J-2)+F(I-2,J+2)+F(I+2,J-2)+F(I+2,J+2) )
     +        +GAM(2)*( F(I-2,J+1)+F(I-2,J-1)+F(I-1,J-2)+F(I-1,J+2)
     +                 +F(I+1,J-2)+F(I+1,J+2)+F(I+2,J-1)+F(I+2,J+1) )
     +        +GAM(3)*( F(I-2,J  )+F(I+2,J  )+F(I  ,J-2)+F(I  ,J+2) )
     +        +GAM(4)*( F(I-1,J-1)+F(I-1,J+1)+F(I+1,J-1)+F(I+1,J+1) )
     +        +GAM(5)*( F(I-1,J  )+F(I  ,J-1)+F(I  ,J+1)+F(I+1,J  ) )
     +        +GAM(6)*  F(I  ,J  )
   10    CONTINUE
C
         J=NRL-1
         FXY(I,J)=
     +     BET(1)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J-3)+F(I+2,J-3) )
     +    +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     +    +BET(3)*( F(I-2,J-1)+F(I+2,J-1) )
     +    +BET(4)*( F(I-2,J-2)+F(I+2,J-2)+F(I-1,J+1)+F(I+1,J+1)
     +             +F(I-1,J-3)+F(I+1,J-3) )
     +    +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     +    +BET(6)*( F(I-1,J-1)+F(I+1,J-1) )
     +    +BET(7)*( F(I-1,J-2)+F(I+1,J-2) )
     +    +BET(8)*( F(I  ,J+1)+F(I  ,J-3) )
     +    +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J-1) +BET(11)*F(I  ,J-2)
C
         J=NRL
         FXY(I,J)=
     +     ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     +    +ALP(2)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J-3)+F(I+2,J-3)
     +             +F(I-1,J-4)+F(I+1,J-4) )
     +    +ALP(3)*( F(I-2,J-2)+F(I+2,J-2) )
     +    +ALP(4)*( F(I-2,J-4)+F(I+2,J-4) )
     +    +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     +    +ALP(6)*( F(I-1,J-1)+F(I+1,J-1)+F(I-1,J-3)+F(I+1,J-3) )
     +    +ALP(7)*( F(I-1,J-2)+F(I+1,J-2) )
     +    +ALP(8)*  F(I  ,J  )
     +    +ALP(9)*( F(I  ,J-1)+F(I  ,J-3) )
     +    +ALP(10)* F(I  ,J-2) +ALP(11)*F(I  ,J-4)
C
   20 CONTINUE
C
      DO 40 I=3,nset-2   ! modified
         DO 30 J=1,NRL
            F(I,J)=FXY(I,J)
   30    CONTINUE
   40 CONTINUE
C
      RETURN
      END
C
