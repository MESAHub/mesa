! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************


! Modified by BP from the OPAL distributed version.
!
! most importantly, I'm filling in the high logT, high logR part of the tables.
! this is done by linearly extrapolating during the read.
! With this change, I get reasonable results out to the edges.


      subroutine init_opal_type2(z_in,xh_in)
         double precision, intent(in) :: z_in, xh_in
         real :: z, xh
         common/type2_recoin/ itimeco,mxzero,readco22_z
         z = z_in; xh = xh_in
         if (z /= readco22_z) itimeco = 0
      end subroutine init_opal_type2

! ***********************************************************************

      subroutine eval_opal_type2 (z_in,xh_in,xxc_in,xxo_in,t6_in,r_in,logKap,data_dir)
      common/type2_e/ opact,dopact,dopacr,dopactd
      double precision :: z_in,xh_in,xxc_in,xxo_in,t6_in,r_in,logKap
      character (len=*), intent(in) :: data_dir
      real :: z,xh,xxc,xxo,t6,r
      z = z_in; xh = xh_in; xxc = xxc_in; xxo = xxo_in; t6 = t6_in; r = r_in
 1    format(a40,1pe26.16)
!.....clip args to range of tables
      if (t6 < 0.0057e0) t6 = 0.0057e0
      if (t6 > 501.18e0) t6 = 501.18e0
      if (r < 1e-8) r = 1e-8
      if (r > 1e1) r = 1e1
      call opac2 (z,xh,xxc,xxo,t6,r,data_dir)
      logKap = opact
      end subroutine eval_opal_type2

! ***********************************************************************

      subroutine open_unit2(m,z,dir)
      use utils_lib, only: mesa_error
         integer, intent(IN) :: m
         real, intent(IN) :: z
         character (len=*) :: dir
         character (len=256) :: zstr, xstr, fname
         integer :: iz, ios
         iz = z * 1d3 + 0.5d0
         zstr = ''
         if (iz == 0) zstr = '000'
         if (iz == 1) zstr = '001'
         if (iz == 4) zstr = '004'
         if (iz == 10) zstr = '010'
         if (iz == 20) zstr = '020'
         if (iz == 30) zstr = '030'
         if (iz == 50) zstr = '050'
         if (iz == 100) zstr = '100'
         xstr = ''
         if (m == 1) xstr = '00'
         if (m == 2) xstr = '03'
         if (m == 3) xstr = '10'
         if (m == 4) xstr = '35'
         if (m == 5) xstr = '70'
         write(fname,'(5A)') trim(dir), '/Gz', trim(zstr), '.x', trim(xstr)
         open(2, FILE=trim(fname), action='read', status='old', iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open ', trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine open_unit2

! ***********************************************************************
      subroutine instruct2
      use utils_lib, only: mesa_error
!-----VERSION of November 20, 1995-----------------------------------------
!-----------------------------------------------------------------------

!          This subroutine contains instructions for using the subroutine
!     OPAC( z, xh, xxc, xxo, t6, r ).
!          The purpose of the subroutine OPAC is to perform 4 or 5 variable
!     interpolation on  log10(kappa).  The opacity tables to be interpolated
!     are known to have somewhat random numerical errors of a few percent.
!     Consequently, adjusting the data prior to performing the interpolation
!     is justified at this level.  The code first reads the original(unsmoothed)
!     tabular data, this data is then passed through a smoothing filter; using
!     a set of routines developed by Mike Seaton (see M.J. Seaton,MNRAS 265,
!     L25,1993). It is the adjusted data that is actually used in carrying out the
!     interpolations in OPAC.  Steps are taken, as described below to insure
!     that the interpolated data is also smooth. The initial adjustment step
!     helps improve the smoothness of the OPAC output,particularly at the
!     smallest values of R. The medium to large R output is only slightly
!     effected by this step. It takes only a few seconds to carryout the initial
!     data smoothing step.  This step can be skipped by setting the parameter
!     ismdata=1 in subroutine readco.

!          The interpolation variables are :

!          xh     The hydrogen mass fraction, X
!          xxc    The enhanced carbon mass fraction, delta Xc.
!                 The total C mass fraction, Xc, is the sum of the initial
!                 amount included in the metal mass fraction, Z, and
!                 delta Xc
!          xxo    The enhanced oxygen mass fraction, delta Xo
!          t6     The temperature in millions of degrees Kelvin, T6
!          r      =rho(gm/cc)/T6**3, R

!     Additional input to OPAC is:

!           z      The metallicity, Z

!          An interpolation between overlapping quadratics is used to obtain
!     smoothed results.  A 4x4 grid in logT6 and logR is used to interpolate
!     in four different 3x3 sub-grids. Linear interpolation between quadratic
!     fits in these different  sub-grids gives smoothed results in both log T6
!     and Log R directions. This procedure produces results that are similar
!     to bicubic spline interpolation, but require storage of only local
!     information.
!         Each of the individual tables in a file Gx**x**z covers 70 temperatures
!     in the range logT=3.75[T6=0.0056341325]( referred to as temperature 1) to
!     logT=8.7[T6=501.187] and 19 values of log R in the range -8 (referred to as 1)
!     to +1; at half-integer steps.  (NOTE: earlier tables were explicitly in
!     terms of T6. For convenience the present tables tabulate log Kappa vs logT. The
!     interpolation however still uses T6 for the temperature variable)
!     For specialized problems, if storage space is a problem, a reduced set of
!     data can be input .  This requires a recompilation with altered parameter
!     values.  In order to limit the range of R, set the parameter nrb= index of
!     lowest value of log R to use(count from log R=-8).  Then set the parameter
!     nre to the index of the largest value of log R to use.  (NOTE: nre-nrb must
!     be at least 5). To ignore the first ntb-1 values of T6 (starting from
!     temperature 1) set the parameter ntb to desired value.  No other parameters
!     should be modified.
!          A five variable interpolation is done when X is greater than 0.  In this
!     case it is assumed that variable values of X are also needed.  We have provided
!     sets of tables for X=0, 0.03, 0.1, 0.35 and 0.7, for each of the metallicities
!     0.0,0.001, 0.004,0.01,0.02,0.03,05 and .1.  The five sets of tables associated with
!     a particular value of Z (eg. Gx0z01, Gx03x01,Gx10z01, Gx35z01,
!     Gx70z01) should be placed in files named codataa, codatab, codatac, codatad,
!     and codatae, respectively.  To create storage for this data the routines must
!     be recompiled with the parameter mx=5.  Again if storage is a problem the T6
!     and log R ranges can be restricted.  This version of the interpolation routine
!     does not interpolate in Z.  Values of Z not in the table can be obtained by
!     interpolating the existing tables to produce similar tables for the Z of interest.
!     Interpolation in xh,xxo,xxc,t6, and r can be obtained as just described
!         A 4 variable interpolation in xxc, xxo, t6 and r is performed in the special
!     case  X=0.  The set of 60 data tables in (xxc,xxo) for a given Z that have
!     been provided for X=0 should be placed in a file called 'codataa'.  This file
!     will be read from unit 2 in the subroutine readco.  In this special case the
!     set of routines provided should be compiled with the parameter mx=1 (there are
!     4 occurrences).  (NOTE: The version of the code described above, intended
!     for X> 0, also  handles this case, but takes more storage space since mx=5).
!          If you want to work with a single value of X which is not zero, i.e.,
!     X=0.03, 0.1, 0.35, or 0.70, then compile the code with mx=1 but change the statement

!          data (xa(i), i=1,5 )/0.0, 0.03, 0.1, 0.35, 0.7/

!     If for example you want to use only the table for X=.70, then

!          data (xa(i), i=1,5 )/0.70, 0.03, 0.1, 0.35, 0.0/

!     You also need to place the tables for X=0.7 into a file named
!     codataa.
!     ***CAUTION***
!         As a result of the mixing procedure used to calculate the data a few
!     X=0.0, low T-small R, table values fell outside the range of T and R accessible
!     from the X=0.35 data directly calculated for this purpose.  These T-R
!     locations are filled in with 9.99 (or for diagnostic purposes in some cases
!     larger values.  At the same locations the derivatives are set to 99.9.  When
!     T-R falls in this region a message is issued by the interpolation code to
!     inform the user of this situation.  Presumable very few users will have
!     applications that take them into this region.

!          Your routine that performs the call to OPAC should include the
!      statement:

!         common/type2_e/ opact,dopact,dopacr,dopactd

!         These variables have the following meanings:

!         OPACT        Is the Log of the Rosseland mean opacity: Log(kappa)
!         DOPACT      Is Dlog(kappa)/Dlog(T6)   ! at constant R
!         DOPACTD     Is Dlog(kappa)/Dlog(T6)   ! at constant density
!         DOPACR      Is Dlog(kappa)/Dlog(R),   ! at constant T6

      dum=0.0
      end subroutine instruct2

! **********************************************************************
      subroutine opac2 (z,xh,xxc,xxo,t6,r,dir)
      use utils_lib, only: mesa_error
!..... The purpose of this subroutine is to interpolate log kappa
!          and obtain smooth derivatives.
!      in C/O abundance and in T6, R,i.e. (xc,xo,T6,R)
!      xxc=Xc=carbon mass fraction
!      xxo=Xo=oxygen mass abundance
!      t6=T6=temperature in millions of degrees kelvin
!      r=R=density(g/cm**3)/T6**3
!..... to use opac insert common/type2_e/ in the calling routine.
!      This common contains interpolated values for kappa and its
!      first derivities.

      save
      character (len=*) :: dir
      integer :: w
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/type2_aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/type2_aa/ q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx),nc,no
      common/type2_a/ co(mx,mc,mo,nt,nr),diag(mx,mc,nt,nr),
     &                index(101),t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr),
     &                dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm),t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/type2_b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
      common/type2_d/dkap
      common/type2_bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
!..... OPACT- opacity obtained from a quadraric interpolation at
!      fixed log T6 at three values of log R; followed by quadratic
!      interpolation along log T6. Results smoothed by mixing
!      overlapping quadratics.
!..... DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics.
!..... DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
      common/type2_e/ opact,dopact,dopacr,dopactd
      common/type2_recoin/ itimeco,mxzero,readco_z

!..... INDEX refers to the index,i, of xc(i) or xo(i); xc(i) and xo(i)
!      are the abundance grid points.
      iop=1   ! provides smoothed interpolations; iop=0 gives no smoothing
      if(nr < 6) GOTO 65
      if((xh > 1.e-6) .and. (mx < 4)) GOTO 64

!..... set-up C/O axis points
      xxco=xxc+xxo
      if(z+xh+xxco-1.e-6 > 1 ) GOTO 61
      zzz=z+0.001
      xxh=xh
       xxci=xxc
      xxoi=xxo
      xxi=xh
      t6i=t6
      ri=r

!..... convert xxc and xxo to logarithmic shifted by Z
      cxx=log10(zzz+xxc)
      oxx=log10(zzz+xxo)
      xxx=log10(0.005+xh)
      slt=log10(t6)
      slr=log10(r)

!..... set X indices
        ilo=2
        ihi=mx
    8   if(ihi-ilo > 1) then
          imd=(ihi+ilo)/2
            if(xh <= xa(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            end if
          GOTO 8
        end if
        i=ihi
        mf=i-2
        mg=i-1
        mh=i
        mi=i+1
        mf2=mi
        istep1=1
        if (mx > 1) then
        istep1=mx-1
        if((xh <= xa(2)+1.e-7) .or. (xh >= xa(istep1)-1.e-7)) mf2=mh
        end if

        if ((mx == 1) .or. (xh < 1.e-6)) then
          mf=1
          mg=1
          mh=1
          mi=2
          mf2=1
        end if

      if (itime(1) /= 12345678) then
      alr(1)=-8.+(nrb-1)*0.5
      do i=2,nr
        alr(i)=alr(i-1)+0.5
      end do
      alt(1)=-2.25+(ntb-1)*0.05
      do i=ntb+1,46
         alt(i)=alt(i-1)+.05
      end do
      ntd=47
      if ((ntb +1) > 47) ntd=ntb+1
      do i=ntd,68
        alt(i)=alt(i-1)+.1
      end do
      do i=68,70
        alt(i)=alt(i-1)+.2
      end do
      do i=1,nt
      t6list(i)=10.**alt(i)
      end do
      end if
        ilo=2
        ihi=nr
   12     if(ihi-ilo > 1) then
          imd=(ihi+ilo)/2
            if(slr <= alr(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            end if
          GOTO 12
          end if
        i=ihi
        l1=i-2
        l2=i-1
        l3=i
        l4=l3+1

        ilo=2
        ihi=nt
   11     if(ihi-ilo > 1) then
          imd=(ihi+ilo)/2
            if(t6 <= t6list(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            end if
          GOTO 11
          end if
        i=ihi
        k1=i-2
        k2=i-1
        k3=i
        k4=k3+1
        k3s=k3+ntb-1
        l3s=l3+nrb-1

!-----set-up indices to skip when low T&R data missing for X=0.
      kmin=0
      k1in=k1
      iadvance=0
      mfin=mf
      if ((mfin == 1) .and. (co(1,1,1,k1,l1) > 9.)) then! data missing
      do i=1,6
        if (co(1,1,1,i,l1) > 9.)  then
          if (xh <= .1) then
           kmin=i+1
          else

            if (iadvance == 0) then  ! sfift X index to avoid X=0.
            iadvance=iadvance+1
            mf=mf+1
            mg=mg+1
            mh=mh+1
            mi=mi+1
            mf2=mf2+1
            end if
          end if
        end if
      end do
      if ((iadvance == 0) .and. (k1 <= kmin) .and. (slt <= alt(kmin))) then
      k1=kmin
      if ((co(1,1,1,kmin,l1+1) < 9.) .and. ((slr+.01) > alr(l1+1))) then
      l1=l1+1
      kmin=0
      k1=k1in
      do i=1,6
      if (co(1,1,1,i,l1) > 9.) kmin=i+1
      end do
      if ((kmin /= 0) .and. (k1in < kmin)) k1=kmin
      end if
      end if
      if ((slt+.001) < alt(k1)) then
!      write (*,'("OPAL data not available for X=", f7.5," logT6=", f7.3, " logR=",f7.3)') xh,slt,slr
      opact=30.
      dopact=99.
      dopacr=99.
      dopactd=99.
      return
      end if
      l2=l1+1
      l3=l2+1
      l4=l3+1
      l3s=l3+nrb-1
      k2=k1+1
      k3=k2+1
      k4=k3+1
      k3s=k3+ntb-1
      end if
!-----end of check for missing data
      do m=mf,mf2
       if(mx >= 4) then
!.....  C and O  fractions determined by the ray through the origin that
!       also passes through the point (Xc,Xo). Specific interpolation
!       values determined by tabulated X values;i.e. xa(m).  Inter-
!       polation along the ray gives log (kappa(Xc,Xo)).  (Advantage
!       of method: keeps indices within table boundaries)
!      Subtract Z to prevent out-of-range C+O values for small X
       if(1.-xh-z > 1.e-6) then
          cmod=(1.-xa(m)-z)/(1.-xh-z)
       else
          cmod=0.
       end if
       xxc=cmod*xxci
       xxo=cmod*xxoi
       cxx=log10(zzz+xxc)
       oxx=log10(zzz+xxo)
       end if
!      ythism=z+xa(m)+xxc+xxo

         do i=1,mc
         xhe=1.-xa(m)-z
         nc=i
         no=i
         xc(i)=xcs(i)
         xo(i)=xos(i)
!          if(xcs(i) >= xhe-1.e-6) then
           if(xcs(i) > xhe) then
            xc(i)=xhe
            xo(i)=xhe
            exit
           end if
         end do

      if(itime(m) /= 12345678) then
      itime(m)=12345678
      mxzero=0
        do i=1,mx
          xx(i)=log10(0.005+xa(i))
          if (xa(i) == 0.0) mxzero=i
        end do
!  ... this is the first time through this m. Calculate the decadic
!      log of the perimeter points shifted by Z+0.001(to avoid divergence
!      at origin); m refers to xa(m); the hydrogen table value.

!      note that the nc-th elements are sometimes needed!
        do i=1,nc
          oxf(m,i)=log10(zzz+xo(i))
          cxf(m,i)=log10(zzz+xc(i))
          xcdf(m,i)=-xo(i)+xo(no)
          xodf(m,i)=-xc(i)+xc(nc)
          cxdf(m,i)=log10(zzz+xcdf(m,i))
          oxdf(m,i)=log10(zzz+xodf(m,i))
        end do

!      note that the nc-th elements are sometimes needed!
        do i=1,nc
          ox(i)=oxf(m,i)
          cx(i)=cxf(m,i)
          xcd(i)=xcdf(m,i)
          xod(i)=xodf(m,i)
          cxd(i)=cxdf(m,i)
          oxd(i)=oxdf(m,i)
        end do

!.....read the data files
      call readco2(z,dir)
      end if

        do i=1,nc
          ox(i)=oxf(m,i)
          cx(i)=cxf(m,i)
          xcd(i)=xcdf(m,i)
          xod(i)=xodf(m,i)
          cxd(i)=cxdf(m,i)
          oxd(i)=oxdf(m,i)
        end do

!..... Determine log R and log T6 grid points to use in the
!      interpolation.

      if((slt < alt(1)).or.(slt > alt(nt))) then
         !write(*,*) 'slt', slt, 'not in range', alt(1), 'to', alt(nt)
         GOTO 62
      end if

      if((slr < alr (1)).or.(slr > alr(nr))) then
         !write(*,*) 'slr', slr, 'not in range', alr(1), 'to', alr(nr)
         GOTO 62
      end if

      if (m == mf) then  !  calculate table indices

      if((mf2 /= mxzero) .and. (k3s > ntm)) then
         !write(*,*) '(mf2 /= mxzero) .and. (k3s > ntm)'
         !write(*,*) 'mf2', mf2, 'k3s', k3s
         GOTO 62
      end if

      do i=14,18 ! check for missing entries at right end of row
         if((l3s > i) .and. (k3s > nta(i+1))) then
            !write(*,*) '(l3s > i) .and. (k3s > nta(i+1))'
            GOTO 62
         end if
      end do
      ip=3
      iq=3
      ntlimit=nta(l3s)
      if((k3s == ntlimit) .or. (iop == 0)) then
        ip=2
        iq=2
      end if
      if(t6 <= t6list(2)+1.e-7 .or. iop == 0) ip=2

      if((l3 == nr) .or. (iop == 0)) then ! right edge of full table
        iq=2
        ip=2
      end if
      if(slr <= alr(2)+1.e-7 .or. iop == 0) iq=2
      end if

      xodp=max(-xxc+xc(nc),0.)
      xcdp=max(-xxo+xo(no),0.)
      is=0

      call cointerp(xxc,xxo)
      end do

        if((zz(mg,1) /= zz(mf,1)) .or. (zz(mh,1) /= zz(mf,1))) then
          write(*,*) 'xxc', xxc
          write(*,*) 'xxo', xxo
          write(*,*) 'zz(mg,1)', zz(mg,1)
          write(*,*) 'zz(mf,1)', zz(mf,1)
          write(*,*) 'zz(mh,1)', zz(mh,1)
          write(*,*) 'zz(mf,1)', zz(mf,1)
          write(*,'("Z does not match Z in codata* files you are using")')
          call mesa_error(__FILE__,__LINE__)
        end if
      if(z /= zz(mf,1)) GOTO 66
      xxc=xxci   ! restores input value; necessary if stop replaced
!                  with return
      xxo=xxoi   ! restores input value
      is=0
      iw=1
      do ir=l1,l1+iq
       do it=k1,k1+ip
         if((mx == 1) .or. (mf2 == 1)) then
           opk(it,ir)=opl(mf,it,ir)
           cycle
         end if
         opk(it,ir)=quad2(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir),opl(mh,it,ir),xx(mf),xx(mg),xx(mh))
         is=1
      end do
      end do

      if (mi == mf2) then  ! interpolate between quadratics
      is=0
      iw=1
       dixr=(xx(mh)-xxx)*dfsx(mh)
      do ir=l1,l1+iq
        do it=k1,k1+ip
        opk2(it,ir)=quad2(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir),opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
        opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.-dixr)
        is=1
        end do
      end do
!     interpolate X between two overlapping quadratics
      end if
      is=0

!..... completed H,C,O interpolation. Now interpolate T6 and log R on a
!      4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure
!      mixes overlapping quadratics to obtain smoothed derivatives.

      call t6rinterp2(slr,slt)
      return

   61 write(*,'(" Mass fractions exceed unity")')
      write(*,*) 'z', z, 'xh', xh, 'xxc', xxc, 'xxo', xxo, 'xxco', xxco, 'sum', z+xh+xxco
      xxc=xxci   ! restores input value; required if stop replaced
!                  with a return
      xxo=xxoi   ! restores input value
      call mesa_error(__FILE__,__LINE__)
   62 opact = 10.0
      !write(*,*) 'interpolation failed'
      !write(*,*) 't6', t6
      !write(*,*) 'r', r
      xxc=xxci   ! restores input value; required if stop replaced
!                  with a return
      xxo=xxoi   ! restores input value
      return
   64 write(*,'(" X not equal to zero: To run this case it is necessary"/ "to recompile with parameter (mx=5)")')
      call mesa_error(__FILE__,__LINE__)
   65 write(*,'("Too few R values; NRE+1-NRB < 6")')
      call mesa_error(__FILE__,__LINE__)
   66 write(*,'(" Z does not match Z in codata* files you are using")')
      call mesa_error(__FILE__,__LINE__)
      end subroutine opac2

! ***********************************************************************
      subroutine cointerp(xxc,xxo)
      use utils_lib, only: mesa_error
!     The purpose of this subroutine is to interpolate in C and O abund-
!     ances.
      save
      integer :: w
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/type2_aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/type2_aa/ q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx),nc,no
      common/type2_a/ co(mx,mc,mo,nt,nr),diag(mx,mc,nt,nr),
     &                index(101),t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr),
     &                dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm),t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/type2_bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
       is=0
      if(xxco < 1.e-6) then
        do ir=l1,l1+iq
          do it=k1,k1+ip
            opl(m,it,ir)=co(m,1,1,it,ir)
          end do
        end do
            is=1
            GOTO 123
      end if
!     include boundaries that could later cause division by 0!
      if(xxc > xcd(3)-1.e-6) then
!__________
      oxdp=log10(zzz+xodp)
!     handle possibility that xodp=0
      fac=max(min((oxx-ox(1))/max(oxdp-ox(1),1.e-6),1.),0.)
      do ir=l1,l1+iq
      do it=k1,k1+ip

!                    interpolation in region c1

!     include boundaries that could later cause division by 0!
      if(xxc > xcd(2)-1.e-6) then
      iw=1
      a(1,m)=quad2(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir),diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
      iw=iw+1
      a(2,m)=quad2(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir),diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
         do w=1,2
           b(w)=a(w,m)
         end do
!     handle possibility that xodp=0
           opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
      is=1
      cycle
      end if
!                    interpolation in region c2

      iw=1
      a(1,m)=quad2(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir),diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
      iw=iw+1
      a(2,m)=quad2(is,iw,cxx,co(m,n(m,2)-2,2,it,ir),co(m,n(m,2)-1,2,it,ir),diag(m,2,it,ir),cx(n(m,2)-2),cx(n(m,2)-1),cxd(2))
      iw=iw+1
      a(3,m)=quad2(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir),diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
        do w=1,3
          b(w)=a(w,m)
        end do
      iw=iw+1
      opl(m,it,ir)=quad2(is,iw,oxx,b(1),b(2),b(3),ox(1),ox(2),oxdp)
      is=1
      end do
      end do
      if(is == 1) GOTO 123
!__________
      end if

!                    interpolation in region c3 to c6
      is=0

      if(nc >= 5) then
!__________
      do i=4,nc-1
!     do not go beyond middle (where c3-c6 overlaps o3-o6), and
        if((xxc > xcd(i)-1.e-6) .and. (xxo > xo(i-1)-1.e-6) .and. (xcd(i-1) > xc(i-1))) then
      do ir=l1,l1+iq
      do it=k1,k1+ip
        oxdp=log10(zzz+xodp)
        iw=1
        m1=i-1
        m2=i-2
        a(1,m)=quad2(is,iw,cxx,co(m,n(m,m2)-2,m2,it,ir),co(m,n(m,m2)-1,m2
     &  ,it,ir),diag(m,m2,it,ir),cx(n(m,m2)-2),cx(n(m,m2)-1),cxd(m2))
        iw=iw+1
        a(2,m)=quad2(is,iw,cxx,co(m,n(m,m1)-2,m1,it,ir),co(m,n(m,m1)-1,m1
     &  ,it,ir),diag(m,m1,it,ir),cx(n(m,m1)-2),cx(n(m,m1)-1),cxd(m1))
        iw=iw+1
        a(3,m)=quad2(is,iw,cxx,diag(m,m2,it,ir),diag(m,m1,it,ir),
     &  diag(m,i,it,ir),cxd(m2),cxd(m1),cxd(i))
         do w=1,3
           b(w)=a(w,m)
         end do
        iw=iw+1
        opl(m,it,ir)=quad2(is,iw,oxx,b(1),b(2),b(3),ox(i-2),ox(i-1),oxdp)
        is=1
      end do
      end do
      if (is == 1) GOTO 123
        end if
      end do
!__________
      end if

      if (is == 1) GOTO 123

!     include boundaries that could later cause division by 0!
      if(xxo > xod(3)-1.e-6) then
!__________
      cxdp=log10(zzz+xcdp)
!     handle possibility that xcdp=0
      fac=max(min((cxx-cx(1))/max(cxdp-cx(1),1.e-6),1.),0.)
      do ir=l1,l1+iq
      do it=k1,k1+ip

!                    interpolation in region  o1

!     include boundaries that could later cause division by 0!
      if(xxo > xod(2)-1.e-6) then
      iw=1
      a(1,m)=quad2(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir),diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
      iw=iw+1
      a(2,m)=quad2(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir),diago(m,no-3,it,ir),oxd(1),oxd(2),oxd(3))
        do w=1,2
          b(w)=a(w,m)
        end do
!     handle possibility that xcdp=0
      opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
      is=1
      cycle
      end if
!                    interpolation in region  o2

      iw=1
      a(1,m)=quad2(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir),diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
      iw=iw+1
      a(2,m)=quad2(is,iw,oxx,co(m,2,n(m,2)-2,it,ir),co(m,2,n(m,2)-1,it,ir),diago(m,no-2,it,ir),ox(n(m,2)-2),ox(n(m,2)-1),oxd(2))
      iw=iw+1
      a(3,m)=quad2(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir),diago(m,nc-3,it,ir),oxd(1),oxd(2),oxd(3))
        do w=1,3
          b(w)=a(w,m)
        end do
      iw=iw+1
      opl(m,it,ir)=quad2(is,iw,cxx,b(1),b(2),b(3),cx(1),cx(2),cxdp)
      is=1
      end do
      end do
      if(is == 1) GOTO 123
!__________
      end if

!                    interpolation in region  o3 to o6
      is=0
      if(no >= 5) then
!__________
      do i=4,no-1
!     do not go beyond middle (where o3-o6 overlaps c3-c6), and
      if((xxo > xod(i)-1.e-6) .and. (xxc > xc(i-1)-1.e-6) .and. (xod(i-1) > xo(i-1)-1.e-6)) then
      do ir=l1,l1+iq
      do it=k1,k1+ip
      cxdp=log10(zzz+xcdp)
      iw=1
      m2=i-2
      m1=i-1
      a(1,m)=quad2(is,iw,oxx,co(m,m2,n(m,m2)-2,it,ir),co(m,m2,n(m,m2)-1,it,ir),
     &             diago(m,no-m2,it,ir),ox(n(m,m2)-2),ox(n(m,m2)-1),oxd(m2))
      iw=iw+1
      a(2,m)=quad2(is,iw,oxx,co(m,m1,n(m,m1)-2,it,ir),co(m,m1,n(m,m1)-1,it,ir),
     &             diago(m,no-m1,it,ir),ox(n(m,m1)-2),ox(n(m,m1)-1),oxd(m1))
      iw=iw+1
      a(3,m)=quad2(is,iw,oxx,diago(m,no-m2,it,ir),diago(m,no-m1,it,ir),diago(m,no-i,it,ir),oxd(m2),oxd(m1),oxd(i))
        do w=1,3
          b(w)=a(w,m)
        end do
      iw=iw+1
      opl(m,it,ir)=quad2(is,iw,cxx,b(1),b(2),b(3),cx(m2),cx(m1),cxdp)
      is=1
      end do
      end do
      if (is == 1) GOTO 123
      end if
      end do
!__________
      end if

      if (is == 1) GOTO 123

!.....find index of C grid.
   52 ie=100*xxc+1
      iei=index(ie)+1
!     must also allow index = nc, to avoid extrapolation
      if (iei > nc) iei=nc

        if(iei > 3) then
          i1=iei-2
          i2=iei-1
          i3=iei
        else
          i1=1
          i2=2
          i3=3
        end if

!.....find index of O grid
      ie=100*xxo+1
      iej=index(ie)+1
!     must also allow index = no, to avoid extrapolation
      if (iej > no) iej=no

        if(iej > 3) then
          j1=iej-2
          j2=iej-1
          j3=iej
        else
          j1=1
          j2=2
          j3=3
        end if

!     lower-O part of grid: interpolate C before O
      if(j3 < no .and. i3 <= n(m,j3) .and. (xxc < xcd(j3)+1.e-6 .or. xxc >= xxo)) then
      do ir=l1,l1+iq
      do it=k1,k1+ip
      iw=0
        do jx=j1,j1+2
          iw=iw+1
!     if i3=n(m,jx), then must replace cx(i3) with cxd(jx)
          a(iw,m)=quad2(is,iw,cxx,co(m,i1,jx,it,ir),co(m,i2,jx,it,ir),co(m,i3,jx,it,ir),cx(i1),cx(i2),min(cx(i3),cxd(jx)))
        end do
        do w=1,3
          b(w)=a(w,m)
        end do
      iw=iw+1
      opl(m,it,ir)=quad2(is,iw,oxx,b(1),b(2),b(3),ox(j1),ox(j2),ox(j3))
      is=1
      end do
      end do
!     else: high-O part of grid: must interpolate O before C
      else
       do ir=l1,l1+iq
       do it=k1,k1+ip
        iw=0
        do ix=i1,i1+2
          iw=iw+1
          if(j3 < n(m,ix)) then
            a(iw,m)=quad2(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir),co(m,ix,j3,it,ir),ox(j1),ox(j2),ox(j3))
          else
            a(iw,m)=quad2(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir),diago(m,no-ix,it,ir),ox(j1),ox(j2),oxd(ix))
          end if
        end do
        do w=1,3
          b(w)=a(w,m)
        end do
      iw=iw+1
      opl(m,it,ir)=quad2(is,iw,cxx,b(1),b(2),b(3),cx(i1),cx(i2),cx(i3))
      is=1
       end do
       end do
      end if
  123 continue
      end subroutine cointerp

! **********************************************************************
      subroutine t6rinterp2(slr,slt)
      use utils_lib, only: mesa_error
!     The purpose of this subroutine is to interpolate in logT6 and logR
      save
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/type2_aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/type2_aa/ q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx),nc,no
      common/type2_a/ co(mx,mc,mo,nt,nr),diag(mx,mc,nt,nr),
     &                index(101),t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr),
     &                dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm),t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/type2_d/dkap
      common/type2_bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
      common/type2_e/ opact,dopact,dopacr,dopactd

      is=0
      iu=0
      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quad2(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3),alr(l1),alr(l2),alr(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quad2(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4),alr(l2),alr(l3),alr(l4))
          end if
        is=1
      end do

      is=0
      iw=1
!..... k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opact=quad2(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
      dopact=dkap
      dkap1=dkap
        if (iq. eq. 3) then
!.....    k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          opactq=quad2(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3))
          dkapq1=dkap
          dopactq=0.d0
        end if
        if(ip == 3) then
!.....    k and Dlog(k)/Dlog(T6) in lower-left 3x3.
          opact2=quad2(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4))
          dkap2=dkap
!.....    k and Dlog(k)/Dlog(T6) smoothed in left 3x4
          dix=(alt(k3)-slt)*dfs(k3)
          dopact=dkap1*dix+dkap2*(1.-dix)
          opact=opact*dix+opact2*(1.-dix)

            if(iq == 3) then
!.....       k and Dlog(k)/Dlog(T6) in upper-right 3x3.
               opactq2=quad2(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
               dkapq2=dkap
               dopactq=dkapq1*dix+dkapq2*(1.-dix)
               opactq=opactq*dix+opactq2*(1.-dix)
            end if
         end if

      iu=0
      do lx=l1,l1+iq
        iw=1
        iu=iu+1
        h(iu)=quad2(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx),alt(k1),alt(k2),alt(k3))
          if(ip == 3) then
            iw=2
            q(iu)=quad2(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx),alt(k2),alt(k3),alt(k4))
          end if
        is=1
      end do

      is=0
      iw=1
!..... k and Dlog(k)/Dlog(R) in lower-left 3x3
      opacr=quad2(is,iw,slr,h(1),h(2),h(3),alr(l1),alr(l2),alr(l3))
      dopacr=dkap
        if(ip == 3) then
          opacrq=quad2(is,iw,slr,q(1),q(2),q(3),alr(l1),alr(l2),alr(l3))
!.....    k and Dlog(k)/Dlog(R) in upper-left 3x3.
          dopacrq=dkap
        end if
        if(iq == 3) then
!.....    k and Dlog(k)/Dlog(R) in lower-right 3x3.
          opact2=quad2(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4))
          dix2=(alr(l3)-slr)*dfsr(l3)
          dopacr=dopacr*dix2+dkap*(1.-dix2)
!.....        k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
              dopact=dopact*dix2+dopactq*(1.-dix2)
              opact=opact*dix2+opactq*(1.-dix2)
         end if
        if(ip == 3) then
         if(iq == 3) then
!.....    k and Dlog(k)/Dlog(R) in upper-right 3x3.
          opacrq=quad2(is,iw,slr,q(2),q(3),q(4),alr(l2),alr(l3),alr(l4))
!.....        Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
              dopacrq=dopacrq*dix2+dkap*(1.-dix2)
            end if
              dopacr=dopacr*dix+dopacrq*(1.-dix)
        end if
      dopactd=dopact-3.*dopacr
!        if (opact > 1.e+15) then
!          write(*,'("Interpolation indices out of range; please report conditions.")')
!          call mesa_error(__FILE__,__LINE__)
!        end if
      if (opact > 9) then
      opact=30.
      dopact=99.
      dopactr=99.
      dopactd=99.
      end if
      end subroutine t6rinterp2

! ***********************************************************************
      subroutine readco2(z,dir)
      use utils_lib, only: mesa_error
!..... The purpose of this subroutine is to read the data tables
      save
      character (len=*) :: dir
      parameter (ismdata=0)
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/type2_aa/ q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx),nc,no
      common/type2_a/ co(mx,mc,mo,nt,nr),diag(mx,mc,nt,nr),
     &                index(101),t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr),
     &                dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm),t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/type2_b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
      common/type2_alink/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(100),coff(100,nr)
      common/type2_CST/NRL,RLS,nset,tmax  ! modified
      common/type2_e/ opact,dopact,dopacr,dopacrd
      character(len=250) :: dumarra
      common/type2_recoin/ itimeco,mxzero,readco_z

      logical :: line_full

        if (itimeco /= 12345678) then
        readco_z = z
        do i=1,mx
          do j=1,mc
            do k=1,mo
              do l=1,nt
                do mq=1,nr
                  co(i,j,k,l,mq)=1.e+35
                end do
              end do
            end do
          end do
        end do
        do i=1,mx
          do j=1,mc
            do l=1,nt
              do mq=1,nr
                diag(i,j,l,mq)=1.e+35
                diago(i,j,l,mq)=1.e+35
              end do
            end do
          end do
        end do
        itimeco=12345678
        end if

      do j=1,nc-1
       do i=1,nc
         if(xcd(j) >= xc(i)) then
           n(m,j)=i+1
           ! BP
           !   the following line was in the original but causes an off-by-1 bug in a single case
           !if(xcd(j) < xc(i)+1.e-6) n(m,j)=i
           !   the bug only shows up for X=0.35 and Z=0.05
           !   the old version thought it should have 8 tables for dX2=0.0, but there are only 7.
           if(xcd(j) < xc(i)+1.e-6) then
             n(m,j)=i
             exit
           end if
         end if
       end do
      end do
      n(m,nc)=0

      close (2)

      call open_unit2(m,z,dir)

!      read header
      do i=1,240
         read(2,'(a)') dumarra(i:i)
      end do

      int=0
      do j=1,no-1
      do i=1,n(m,j)
        int=int+1
        read(2,'(f10.5)') dum
        read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)') itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
        xca(m,int)=min(xca(m,int),1.-x(m,int)-zz(m,int)-xoa(m,int))
        read(2,'(f10.5)') dum,dum,dum
        read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
        read(2,'(f10.5)') dum
          do k=1,ntm
            read(2,fmt='(f4.2)',advance='no') altin
            line_full = .true.
            do l=1,nrm
               read(2,fmt='(19f7.3)',advance='no',iostat=ioerror) cof(k,l)
               if (ioerror /= 0) then ! have reached end of data for this temperature
                  if (l < 3) then
                     write(*,*) 'lost in reading -- expected to find at least 2 entries on this line!'
                     call mesa_error(__FILE__,__LINE__)
                  end if
                  do ll=l,nrm
                     cof(k,ll) = 2. * cof(k,ll-1) - cof(k,ll-2)
                  end do
                  line_full = .false.
                  if (.false.) then
                     write(*,fmt='(f4.2)',advance='no') altin
                     do ll=1,nrm
                        write(*,fmt='(f7.3)',advance='no') cof(k,ll)
                     end do
                     write(*,*)
                  end if
                  exit
               end if
            end do

            if (line_full) read(2,'(f10.5)') dum ! advance to next line

            do ll=1,nrm   ! modified
              coff(k,ll)=cof(k,ll)
            end do

          end do

          if (isett6 /= 1234567) then
          do k=1,ntm
            t6arr(k)=t6list(k)
          end do
          end if
          isett6=1234567

       if (ismdata == 0) then
         if ((nrm /= nr) .or. (ntm /= nt)) then
           write (*,'("Not set up to smooth data with reduced T-Rho grid")')
           call mesa_error(__FILE__,__LINE__)
         end if
        tmax=10.   ! modified
        nset=67    ! 65 in earlier version
        RLS=-8.
        nsm=1
          RLE=1.0
          nrlow=1
          nrhigh=2*(RLE-RLS)+1

        call opaltab2    !modified
       end if

      ll=1
      do kk=nrb,nre
        alr(ll)=alrf(kk)
          do k=1,nt
            if (ismdata == 0) then
              if ((m == 1) .and. (k <= 9)) then
                co(m,i,j,k,ll)=cof(k+ntb-1,kk)
              else
                co(m,i,j,k,ll)=coff(k+ntb-1,kk)
              end if
            else
              co(m,i,j,k,ll)=coff(k+ntb-1,kk)
            end if
          end do
        ll=ll+1
      end do
      end do ! i
      end do ! j
      if(x(m,1) /= xa(m)) then
      write(*,'(" X in the codata? file does not match xa(m)")')
      call mesa_error(__FILE__,__LINE__)
      end if

      do i=1, nc-1
       do k=1,nt
        do l=1,nr
          diag(m,i,k,l)=co(m,n(m,i),i,k,l)
        end do
       end do
      end do

      do j=1,no-1
        int=int+1
        read(2,'(f10.5)') dum
        read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)') itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
        read(2,'(f10.5)') dum,dum,dum
        read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
        read(2,'(f10.5)') dum

          ! BP extrapolate to fill table
          do k=1,ntm
            read(2,fmt='(f4.2)',advance='no') altin
            line_full = .true.
            do l=1,nrm
               read(2,fmt='(19f7.3)',advance='no',iostat=ioerror) cof(k,l)
               if (ioerror /= 0) then ! have reached end of data for this temperature
                  if (l < 3) then
                     write(*,*) 'lost in reading -- expected to find at least 2 entries on this line!'
                     call mesa_error(__FILE__,__LINE__)
                  end if
                  do ll=l,nrm
                     cof(k,ll) = 2. * cof(k,ll-1) - cof(k,ll-2)
                  end do
                  line_full = .false.
                  if (.false.) then
                     write(*,fmt='(f4.2)',advance='no') altin
                     do ll=1,nrm
                        write(*,fmt='(f7.3)',advance='no') cof(k,ll)
                     end do
                     write(*,*)
                  end if
                  exit
               end if
            end do

            if (line_full) read(2,'(f10.5)') dum ! advance to next line

            do ll=1,nrm   ! modified
              coff(k,ll)=cof(k,ll)
            end do

          end do


!         do k=1,ntm
!           read (2,'(f4.2,19f7.3)') dum, (cof(k,l), l=1,nrm)
!     set up to smooth final "diago" opacity tables
!           do l=1,nrm
!              coff(k,l)=cof(k,l)
!           end do
!        end do


!     smooth final "diago" opacity tables too!
       if (ismdata == 0) then
        tmax=10.   ! modified
        nset=67 !65 in earlier version
        RLS=-8.
        nsm=1
          RLE=1.0
          nrlow=1
          nrhigh=2*(RLE-RLS)+1
        call opaltab2    !modified
        do k=3,NTEMP-2
          do ll=nrlow,nrhigh
!           Following skip required because, due to missing data,
!           the X=0  low T data cannot be smoothed
            if ((m == 1) .and. (k <= 9)) then
              cof(k,ll)=cof(k,ll)
            else
              cof(k,ll)=coff(k,ll)
            end if

          end do
        end do
       ll=1
       do kk=nrb,nre
         do k=1,nt
           diago(m,j,k,ll)=cof(k+ntb-1,kk)
         end do
       ll=ll+1
       end do
       end if
      end do

      do i=2,nt
        dfs(i)=1./(alt(i)-alt(i-1))
      end do
      do i=2,nr
       dfsr(i)=1./(alr(i)-alr(i-1))
      end do
      istep=-1
      if (mx > 1 ) then
        istep=1
        do i=2,mx,istep
          dfsx(i)=1./(xx(i)-xx(i-1))
        end do
      end if

      end subroutine readco2

! ***********************************************************************
      function quad2(ic,i,x,y1,y2,y3,x1,x2,x3)
      use utils_lib, only: mesa_error
!..... this function performs a quadratic interpolation.
      save
      common/type2_d/dkap
      common/type2_coquad/ xx(3),yy(3),xx12(30),xx13(30),xx23(30),xx1sq(30),xx1pxx2(30)
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      yy(1)=y1
      yy(2)=y2
      yy(3)=y3
        if(ic == 0) then
          if (xx(2) == xx(3)) then
             xx12(i)=1./(xx(1)-xx(2))
             xx13(i)=0
             xx23(i)=0
             xx1sq(i)=0
             xx1pxx2(i)=0
          else if (xx(1) == xx(2) .or. xx(1) == xx(3)) then
               write(*,*) 'divide by 0 in quad'
               write(*,*) 'x', x
               write(*,*) 'xx(1)', xx(1)
               write(*,*) 'xx(2)', xx(2)
               write(*,*) 'xx(3)', xx(3)
               write(*,*) 'yy(1)', yy(1)
               write(*,*) 'yy(2)', yy(2)
               write(*,*) 'yy(3)', yy(3)
               write(*,*)
               call mesa_error(__FILE__,__LINE__)
          else
             xx12(i)=1./(xx(1)-xx(2))
             xx13(i)=1./(xx(1)-xx(3))
             xx23(i)=1./(xx(2)-xx(3))
             xx1sq(i)=xx(1)*xx(1)
             xx1pxx2(i)=xx(1)+xx(2)
          end if
        end if
      c3=(yy(1)-yy(2))*xx12(i)
      c3=c3-(yy(2)-yy(3))*xx23(i)
      c3=c3*xx13(i)
      c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3
      c1=yy(1)-xx(1)*c2-xx1sq(i)*c3
      dkap=c2+(x+x)*c3
      quad2=c1+x*(c2+x*c3)
      end function quad2

! ***********************************************************************
      block data type2
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/type2_aa/ q(4),h(4),xcd(mc),xod(mc),xc(mc),xo(mo),xcs(mc),xos(mo),cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx),nc,no
      common/type2_a/ co(mx,mc,mo,nt,nr),diag(mx,mc,nt,nr),
     &                index(101),t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr),
     &                dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm),t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/type2_b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs),zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
      common/type2_aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc),opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/type2_recoin/ itimeco,mxzero,readco_z
      data itime/mx*0/
      data itimeco/0/
      data ( index(i),i=1,101)/1,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,
     &                         5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     &                         6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
     &                         7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
      data (xcs(i),i= 1,mc)/ 0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/
      data (xos(i),i= 1,mc)/0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/
      data (xa(i),i=1,5)/0.0,0.03,0.1,0.35,0.7/
      ! NOTE: nta tells how many logT values there are for a given column
      ! I've changed this to say that all columns have all logT.
      ! I've filled in the missing values by linear extrapolation during the read.   BP. July, 2007.
      data (nta(i),i=1,nrm)/19*70/ !/14*70,69,64,60,58,57/
      end

! ********************************************************************
      subroutine SPLINE2(X,Y,N,Y2)
      parameter (NMAX=100)
      dimension X(N),Y(N),Y2(N),U(NMAX)

!     FIRST DERIVATIVES AT END POINTS USING CUBIC FIT
         YP1=((Y(3)-Y(1))*(X(2)-X(1))**2-(Y(2)-Y(1))*(X(3)-X(1))**2)/((X(3)-X(1))*(X(2)-X(1))*(X(2)-X(3)))
         YPN=((Y(N-2)-Y(N))*(X(N-1)-X(N))**2-(Y(N-1)-Y(N))*(X(N-2)-X(N))**2)/((X(N-2)-X(N))*(X(N-1)-X(N))*(X(N-1)-X(N-2)))

      Y2(1)=-0.5
      U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      do I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      end do
      QN=0.5
      UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      do K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      end do
      end subroutine SPLINE2
! ********************************************************************
      subroutine SPLINT2(XA,YA,N,Y2A,X,Y,YP)
      use utils_lib, only: mesa_error
      dimension XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     if (KHI-KLO > 1) then
        K=(KHI+KLO)/2
        if(XA(K) > X) then
          KHI=K
        else
          KLO=K
        end if
      GOTO 1
      end if
      H=XA(KHI)-XA(KLO)
      if (H == 0.) then
         write(*,*) 'Bad XA input.'
         call mesa_error(__FILE__,__LINE__)
      end if
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      YP=0.05*((-YA(KLO)+YA(KHI))/H + (-(3*A**2-1)*Y2A(KLO)+(3*B**2-1)*Y2A(KHI))*H/6.)
      end subroutine SPLINT2
! ********************************************************************
      subroutine FITY2

!  THIS ROUTINE MAKES SPLINE FITS FOR F AND FX, AND OBTAINS
!  FY AND FXY


      common/type2_CST/NRL,RLS,nset,tmax  ! modified

      parameter(IPR=20)
      dimension A(IPR),B(IPR),AD(IPR),BD(IPR)
      common/type2_CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)

      do I=1,nset   ! modified
         do J=1,NRL
            A(J)=F(I,J)
            B(J)=FX(I,J)
         end do

         call getd2(A,NRL,AD,AP1,APN)
         call getd2(B,NRL,BD,BP1,BPN)

         FY(I,1)=AP1
         FY(I,NRL)=APN
         FXY(I,1)=BP1
         FXY(I,NRL)=BPN
         do J=2,NRL-1
            FY(I,J)= -A(J)+A(J+1)-2.*AD(J)-AD(J+1)
            FXY(I,J)=-B(J)+B(J+1)-2.*BD(J)-BD(J+1)
         end do
      end do

      end subroutine FITY2
! ********************************************************************
      subroutine FITX2

!  THIS ROUTINE IS USED ONLY AFTER SMOOTHING.
!  ITS FUNCTION IS TO RECOMPUTE FX USING SMOOTHED F.

      parameter(IPR=20)
      dimension A(85),D(85)

      common/type2_CST/NRL,RLS,nset,tmax  ! modified
      common/type2_CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)

      do J=1,NRL
         do I=1,nset ! modified
            A(I)=F(I,J)
         end do
         call getd2(A,nset,D,AP1,APN)  ! modified
         FX(1,J)=AP1
         FX(nset,J)=APN   ! modified
         do I=2,nset-1  ! modified
            FX(I,J)=-A(I)+A(I+1)-2.*D(I)-D(I+1)
         end do
      end do

      end subroutine FITX2

! ***********************************************************************
      subroutine GETD2(F,N,D,FP1,FPN)

!  SIMPLIFIED CODE FOR SPLINE COEFFICIENTS, FOR CASE OF INTERVALS
!  OF UNITY.

      dimension F(N),D(N),T(85)

      FP1=(-11.*F(1)+18.*F(2)-9.*F(3)+2.*F(4))/6.
      FPN=(11.*F(N)-18.*F(N-1)+9.*F(N-2)-2.*F(N-3))/6.

      D(1)=-.5
      T(1)=.5*(-F(1)+F(2)-FP1)

      do J=2,N-1
         D(J)=-1./(4.+D(J-1))
         T(J)=-D(J)*(F(J-1)-2.*F(J)+F(J+1)-T(J-1))
      end do

      D(N)=(FPN+F(N-1)-F(N)-T(N-1))/(2.+D(N-1))

      do J=N-1,1,-1
         D(J)=D(J)*D(J+1)+T(J)
      end do

      end subroutine GETD2

! ********************************************************************
      subroutine INTERP2(FLT,FLRHO,G,DGDT,DGDRHO,IERR)

!  GIVEN F,FX,FY AND FXY ON THE GRID POINTS, THIS ROUTINE
!  DOES BI-CUBIC INTERPOLATIONS USING METHODS DESCRIBED IN
!  Numerical Recipes, PP. 118 TO 120


      parameter(IPR=20)
      common/type2_CF/ F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      dimension B(16)
      logical :: IERR

      common/type2_CST/NRL,RLS,nset,tmax  ! modified

!  EXTREME LIMITS ALLOWED ARE:-
!     (3.800-0.0125) TO (8.000+0.0125) FOR LOG10(T)
!     (RLS-0.125) TO (RLE+0.1254) FOR LOG10(R)
!     (ALLOWING FOR SMALL EXTRAPOLATIONS BEYOND TABULAR VALUES)

!  FUNCTION DEFINITIONS FOR CUBIC EXPANSION

      FF(S,T)=    B( 1)+T*(B( 2)+T*(B( 3)+T*B( 4)))
     &   +S*(     B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     &   +S*(     B( 9)+T*(B(10)+T*(B(11)+T*B(12)))
     &   +S*(     B(13)+T*(B(14)+T*(B(15)+T*B(16))) )))

      FFX(S,T)=   B( 5)+T*(B( 6)+T*(B( 7)+T*B( 8)))
     &   +S*(  2*(B( 9)+T*(B(10)+T*(B(11)+T*B(12))))
     &   +S*(  3*(B(13)+T*(B(14)+T*(B(15)+T*B(16)))) ))

      FFY(S,T)=   B( 2)+S*(B( 6)+S*(B(10)+S*B(14)))
     &   +T*(  2*(B( 3)+S*(B( 7)+S*(B(11)+S*B(15))))
     &   +T*(  3*(B( 4)+S*(B( 8)+S*(B(12)+S*B(16)))) ))

      IERR=.false.

      X=20.*(FLT-3.800)+1
      FLR=FLRHO+18.-3.*FLT
      Y=2*( FLR - RLS )+1

      I=0
      if(X < 2.) then
         if(X < 0.75) then
            IERR=.true.
         else
            I=1
         end if
      else if(X > 84) then
         if(X > 85.25) then
            IERR=.true.
         else
            I=84
         end if
      else
         I=X
      end if
      U=X-I

      if(Y < 2.) then
         if(Y < 0.75) then
            IERR=.true.
         else
            J=1
         end if
      else if(Y > NRL-1) then
         if(Y > NRL+.25) then
            IERR=.true.
         else
            J=NRL-1
         end if
      else
         J=Y
      end if
      V=Y-J

      if (IERR) then
         G=9.999
         DGDT=9.999
         DGDRHO=9.999
          return
      end if


!  GIVEN FUNCTIONS AND DERIVATIVES AT GRID POINTS, COMPUTE COEFFICIENTS.
      B(1)=F(I,J)
      B(2)=FY(I,J)
      B(3)=3*(-F(I,J)+F(I,J+1))-2*FY(I,J)-FY(I,J+1)
      B(4)=2*(F(I,J)-F(I,J+1))+FY(I,J)+FY(I,J+1)

      B(5)=FX(I,J)
      B(6)=FXY(I,J)
      B(7)=3*(-FX(I,J)+FX(I,J+1))-2*FXY(I,J)-FXY(I,J+1)
      B(8)=2*(FX(I,J)-FX(I,J+1))+FXY(I,J)+FXY(I,J+1)

      B(9)=3*(-F(I,J)+F(I+1,J))-2*FX(I,J)-FX(I+1,J)
      B(10)=3*(-FY(I,J)+FY(I+1,J))-2*FXY(I,J)-FXY(I+1,J)
      B(11)=9*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     & +6*(FX(I,J)-FX(I,J+1)+FY(I,J)-FY(I+1,J))
     & +4*FXY(I,J)
     & +3*(FX(I+1,J)-FX(I+1,J+1)-FY(I+1,J+1)+FY(I,J+1))
     & +2*(FXY(I,J+1)+FXY(I+1,J))
     & +FXY(I+1,J+1)
      B(12)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     & +4*(-FX(I,J)+FX(I,J+1))
     & +3*(-FY(I,J)+FY(I+1,J)+FY(I+1,J+1)-FY(I,J+1))
     & +2*(-FX(I+1,J)+FX(I+1,J+1)-FXY(I,J)-FXY(I,J+1))
     & -FXY(I+1,J)-FXY(I+1,J+1)

      B(13)=2*(F(I,J)-F(I+1,J))+FX(I,J)+FX(I+1,J)
      B(14)=2*(FY(I,J)-FY(I+1,J))+FXY(I,J)+FXY(I+1,J)
      B(15)=6*(-F(I,J)+F(I+1,J)-F(I+1,J+1)+F(I,J+1))
     & +4*(-FY(I,J)+FY(I+1,J))
     & +3*(-FX(I,J)-FX(I+1,J)+FX(I+1,J+1)+FX(I,J+1))
     & +2*(FY(I+1,J+1)-FY(I,J+1)-FXY(I,J)-FXY(I+1,J))
     & -FXY(I+1,J+1)-FXY(I,J+1)
      B(16)=4*(F(I,J)-F(I+1,J)+F(I+1,J+1)-F(I,J+1))
     & +2*(FX(I,J)+FX(I+1,J)-FX(I+1,J+1)-FX(I,J+1)
     &    +FY(I,J)-FY(I+1,J)-FY(I+1,J+1)+FY(I,J+1))
     & +FXY(I,J)+FXY(I+1,J)+FXY(I+1,J+1)+FXY(I,J+1)

!  GET G=LOG10(ROSS), DGDT=d LOG10(ROSS)/d LOG10(T),
!      DGDRHO=d LOG10(ROSS)/d LOG10(RHO)
      G=FF(U,V)
      DGDT=20.*FFX(U,V)-6.*FFY(U,V)
      DGDRHO=2.*FFY(U,V)

      end subroutine INTERP2

! ********************************************************************
      subroutine SMOOTH2

!  THIS subroutine USES A 2-DIMENSIONAL GENERALISATION OF THE SMOOTHING
!  TECHNIQUES DESCRIBED ON PP. 644 TO 649 OF Numerical Recipes.

!  CONSIDER THE 25 POINTS DEFINED BY
!       I+n, n=-2,-1,0,1,2 AND J+m, m=-2,-1,0,1,2.
!  THE FUNCTION TO BE SMOOTHED IS FITTED TO A BI-CUBIC, INVOLVING
!  16 COEFFICIENTS, USING TECHNIQUES OF LEAST-SQUARES. THE SMOOTHED
!  FUNCTION (TEMPORARILY STORED IN FXY) IS GIVEN BY THE FITTED VALUE
!  AT THE POINT I AND J.

!  THE FITTING IS SHIFTED FOR POINTS CLOSE TO BOUNDARIES.


      parameter(IPR=20)

      common/type2_CST/NRL,RLS,nset,tmax  ! modified
      common/type2_CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)

      dimension GAM(6)
      DATA GAM/+0.0073469388,-0.0293877551,-0.0416326531,
     &         +0.1175510204,+0.1665306122,+0.2359183673/
      dimension BET(11)
      DATA BET/
     & -0.0048979592,-0.0661224490,-0.0293877551,+0.0195918367,
     &  0.2644897959,+0.1175510204,-0.0783673469,+0.0277551020,
     &  0.3746938776,+0.1665306122,-0.1110204082/
      dimension ALP(11)
      DATA ALP/
     & -0.0844897959,-0.0048979592,+0.0073469388,+0.0012244898,
     &  0.3379591837,+0.0195918367,-0.0293877551,+0.4787755102,
     &  0.0277551020,-0.0416326531,-0.0069387755/


      do I=3,nset-2

         J=1
         FXY(I,J)=
     &    ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     &   +ALP(2)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J+3)+F(I+2,J+3)
     &          +F(I-1,J+4)+F(I+1,J+4) )
     &   +ALP(3)*( F(I-2,J+2)+F(I+2,J+2) )
     &   +ALP(4)*( F(I-2,J+4)+F(I+2,J+4) )
     &   +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     &   +ALP(6)*( F(I-1,J+1)+F(I+1,J+1)+F(I-1,J+3)+F(I+1,J+3) )
     &   +ALP(7)*( F(I-1,J+2)+F(I+1,J+2) )
     &   +ALP(8)*  F(I  ,J  )
     &   +ALP(9)*( F(I  ,J+1)+F(I  ,J+3) )
     &   +ALP(10)* F(I  ,J+2) +ALP(11)*F(I  ,J+4)

         J=2
         FXY(I,J)=
     &    BET(1)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J+3)+F(I+2,J+3) )
     &   +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     &   +BET(3)*( F(I-2,J+1)+F(I+2,J+1) )
     &   +BET(4)*( F(I-2,J+2)+F(I+2,J+2)+F(I-1,J-1)+F(I+1,J-1)
     &            +F(I-1,J+3)+F(I+1,J+3) )
     &   +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     &   +BET(6)*( F(I-1,J+1)+F(I+1,J+1) )
     &   +BET(7)*( F(I-1,J+2)+F(I+1,J+2) )
     &   +BET(8)*( F(I  ,J-1)+F(I  ,J+3) )
     &   +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J+1) +BET(11)*F(I  ,J+2)

         do J=3,NRL-2
            FXY(I,J)=
     &         GAM(1)*( F(I-2,J-2)+F(I-2,J+2)+F(I+2,J-2)+F(I+2,J+2) )
     &        +GAM(2)*( F(I-2,J+1)+F(I-2,J-1)+F(I-1,J-2)+F(I-1,J+2)
     &                 +F(I+1,J-2)+F(I+1,J+2)+F(I+2,J-1)+F(I+2,J+1) )
     &        +GAM(3)*( F(I-2,J  )+F(I+2,J  )+F(I  ,J-2)+F(I  ,J+2) )
     &        +GAM(4)*( F(I-1,J-1)+F(I-1,J+1)+F(I+1,J-1)+F(I+1,J+1) )
     &        +GAM(5)*( F(I-1,J  )+F(I  ,J-1)+F(I  ,J+1)+F(I+1,J  ) )
     &        +GAM(6)*  F(I  ,J  )
         end do

         J=NRL-1
         FXY(I,J)=
     &     BET(1)*( F(I-2,J+1)+F(I+2,J+1)+F(I-2,J-3)+F(I+2,J-3) )
     &    +BET(2)*( F(I-2,J  )+F(I+2,J  ) )
     &    +BET(3)*( F(I-2,J-1)+F(I+2,J-1) )
     &    +BET(4)*( F(I-2,J-2)+F(I+2,J-2)+F(I-1,J+1)+F(I+1,J+1)
     &             +F(I-1,J-3)+F(I+1,J-3) )
     &    +BET(5)*( F(I-1,J  )+F(I+1,J  ) )
     &    +BET(6)*( F(I-1,J-1)+F(I+1,J-1) )
     &    +BET(7)*( F(I-1,J-2)+F(I+1,J-2) )
     &    +BET(8)*( F(I  ,J+1)+F(I  ,J-3) )
     &    +BET(9)*F(I  ,J  ) +BET(10)*F(I  ,J-1) +BET(11)*F(I  ,J-2)

         J=NRL
         FXY(I,J)=
     &     ALP(1)*( F(I-2,J  )+F(I+2,J  ) )
     &    +ALP(2)*( F(I-2,J-1)+F(I+2,J-1)+F(I-2,J-3)+F(I+2,J-3)
     &             +F(I-1,J-4)+F(I+1,J-4) )
     &    +ALP(3)*( F(I-2,J-2)+F(I+2,J-2) )
     &    +ALP(4)*( F(I-2,J-4)+F(I+2,J-4) )
     &    +ALP(5)*( F(I-1,J  )+F(I+1,J  ) )
     &    +ALP(6)*( F(I-1,J-1)+F(I+1,J-1)+F(I-1,J-3)+F(I+1,J-3) )
     &    +ALP(7)*( F(I-1,J-2)+F(I+1,J-2) )
     &    +ALP(8)*  F(I  ,J  )
     &    +ALP(9)*( F(I  ,J-1)+F(I  ,J-3) )
     &    +ALP(10)* F(I  ,J-2) +ALP(11)*F(I  ,J-4)

      end do

      do I=3,nset-2   ! modified
         do J=1,NRL
            F(I,J)=FXY(I,J)
         end do
      end do

      end

! ********************************************************************
      subroutine opaltab2

!  CODE FOR FITTING AND SMOOTHING OPAL DATA. ADAPTED FROM A CODE
!     WRITTEN BY MIKE SEATON(obtained june 1993)

!     OPAL DATA.
!     ASSUMES FIRST T6=0.006, LAST T6=10.OR 0.04). Depending on position
!     in the table.
!     USES RECTANGULAR ARRAY FOR VARIABLES T6 AND LOG10(R)

!     (1) NSM=NUMBER OF PASSES THROUGH SMOOTHING FILTER.
!     USE OF NSM=1 OR 2 IS RECOMMENDED.
!     NO SMOOTHING WITH NSM=0
!     (2) RANGE FOR LOG10(R),
!     RLS=FIRST VALUE, RLE=LAST VALE
!     (RLS MUST BE FIRST VALUYE IN TABLE)

!  subroutine INTERP
!     AFTER PROCESSING, DATA ARE IN A FORM FOR USE OF
!               subroutine INTERP
!     WHICH GIVES LOG(ROSS) AND TWO FIRST DERIVATIVES FOR ANY
!     VALUES OF LOG(T) AND LOG(RHO). SEE BELOW FOR FURTHER
!     EXPLANATION.

!  OUTPUT FOR THE CASE OF NSM > 0.
!     INTERP IS USED TO OBTAIN SMOOTHED DATA INTERPOLATED
!     BACK TO THE ORIGINAL OPAL MESH. TWO FILES ARE WRITTEN.


!  THE subroutineS SPLINE AND SPLINT ARE ADAPTED FROM THOSE GIVE BY
!  W.H. Press, S.A. Teulolsky, W.T. Vettering and B.P. Flannery,
!  "Numerical Recipes in FORTRAN", 2nd ed., 1992, C.U.P.
!  OTHER REFERENCES ARE MADE TO METHODS DESCRIBED IN THAT BOOK.

      parameter(IP=100,IPR=20)
      dimension U(IP),ROSSL(IP,IPR),V(IP),V2(IP)
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb,ntabs=60,ntm=70,ntb=1,nt=ntm+1-ntb)
      common/type2_CF/F(85,IPR),FX(85,IPR),FY(85,IPR),FXY(85,IPR)
      character(len=100) :: HEAD
      common/type2_CST/NRL,RLS,nset,tmax  ! modified
      common/type2_alink/ N,NSM,nrlow,nrhigh,RLE,t6arr(100),coff(100,nr)
      logical :: IERR


      NRL=2*(RLE-RLS)+1

!     STORE LOG10(T) IN U AND LOG10(ROSS) IN ROSSL
!     CHECK FIRST VALUE OF T6
      T6=t6arr(1)
      do j=1,NRL
      ROSSL(1,j)=coff(1,j)
      end do
      U(1)=0.d0
      if( abs(T6 -.0056341325) < 1.e-8) then
         U(1)=6.+LOG10(T6)
      end if
!     SET ROSSL UP TO T6=t6arr(nset)
      I=1
    5 I=I+1
      T6=t6arr(I)
      do j=1,NRL
      ROSSL(I,j)=coff(I,j)
      end do
         U(I)=6+LOG10(T6)
         if(T6 < tmax)GOTO 5
      N=I
      if(N > IP) then
         write(*,*) ' REQUIRE parameter IP OF AT LEAST ',N
         stop 1
      end if


!     DEFINE VARIABLES
!         X=20.0*(LOG10(T)-3.80)+1
!         Y=2.0*(LOG10(R)-RLS)+1
!     USE INDICES I=1 TO nset AND J=1 TO NRL
!     X AND Y ARE SUCH THAT, ON MESH-POINT (I,J), X=I AND Y=J
!     OBTAIN:-
!         F(I,J)=LOG10(ROSS)
!         FX(I,J)=dF/dX
!         FY(I,J)=dF/dY
!         FXY(I,J)=ddF/dXdY


!     FIRST GET F AND FX, INTERPOLATING FROM OPAL T6 TO
!     INTERVAL OF 0.05 IN LOG10(T).
      do J=1,NRL
!        FOR EACH LOG10(R), STORE LOG10(ROSS) IN V(I)
         do I=1,N
            V(I)=ROSSL(I,J)
         end do

!        GET FIRST DERIVATIVES AT END POINTS

!        GET SECOND DERIVATIVES FOR SPLINE FIT
         call spline2(U,V,N,V2)

!        INTERPOLATE TO LOG10(T)=FLT, FLT=3.8(0.05)8.0
         do I=1,nset ! modified
            FLT=3.75+0.05*I
            call splint2(U,V,N,V2,FLT,F(I,J),FX(I,J))
         end do

      end do


!  OPTION FOR SMOOTHING
      if(NSM > 0) then
         do NS=1,NSM
            call smooth2
         end do
         call fitx2
      end if

!  GET FY AND FXY
      call fity2

!  THE ARRAYS F, FX, FY AND FXY ARE NOW STORED

!  CAN NOW DO INTERPOLATIONS USING
!       call interp2(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
!       INPUT IS FLT=LOG10(T), FLRHO=LOG10(RHO)
!       OUTPUT IS G=LOG10(ROSS)
!              DGDT=dG/d(LOG10(T))
!            DGDRHO=dG/d(LOG10(RHO))
!              IERR=.true. if INPUT FLT, FLRHO ARE OUT-OF-RANGE,
!                          elseIERR=.false.

! INTERPOLATE BACK TO OPAL POINTS
      if(NSM > 0) then
         do l=1,NRL
         coff(1,l)=ROSSL(1,l)
         end do

         do K=2,N
            FLT=U(K)
            do L=nrlow,nrhigh
               FLR=RLS+.5*(L-1)
               FLRHO=FLR-18.+3.*FLT
               call interp2(FLT,FLRHO,G,DGDT,DGDRHO,IERR)
               if(IERR) then
               end if
               V(L)=G
            end do
            T6=t6arr(K)
            do l=nrlow,nrhigh
            coff(K,l)=V(l)

            end do

         end do
      end if
      return

 1000 FORMAT('  SMOOTHED OPAL DATA')
 1100 FORMAT('  OPAL DATA, (SMOOTHED-ORIGINAL)')
 2000 FORMAT(100A1)
 2222 FORMAT(F8.3,20F7.3)
 6000 FORMAT(/' FIRST T6=',1P,E10.3,', SHOULD BE 0.006')
 6003 FORMAT(/' !!! OUT-OF-RANGE !!!'/' FLT=',1P,E10.3,', FLRHO=',E10.3,', FLR=',E10.3)

      end subroutine opaltab2
