!     based on EOS5_xtrin_H-He.f from the OPAL release

! ************************************************************************
      subroutine esac_hhe (xh,t6,r,iorder,irad,filename,results_out,pgas_out,prad_out,info)
!..... The purpose of this subroutine is to interpolate
!      the equation of state and its derivatives in X, T6, density
!        izi=0 recalculate table indices to use; =1 keep previous

!        xh=hydrogen mass fraction
!        ztab is metal fraction of the EOS5_data_H-He tables you are using.
!           included only for purpose of preventing mismatch
!        t6=T6=temperature in millions of degrees kelvin
!        r=rho=Rho=density(g/cm**3)
!..... to use esac insert common/eeos_hhe/ esact,eos(mv) in the calling routine.
!      This common contains the interpolated EOS values for the EOS

!..... eos(i) are obtained from a quadratic interpolation at
!      fixed T6 at three values of Rho; followed by quadratic
!      interpolation along T6. Results smoothed by mixing
!      overlapping quadratics.
!         definitions:

!            T6 is the temperature in units of 10**6 K

!            rho is the density in grams/cc
!            R=Rho/T6**3

!            eos(1) is the pressure in megabars (10**12dyne/cm**2)
!            eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
!            eos(3) is the entropy in units of energy/T6
!            eos(4) is dE/dRHO at constant T6
!            eos(5) is the specific heat, dE/dT6 at constant V.
!            eos(6) is dlogP/dlogRho at constant T6.
!                   Cox and Guil1 eq 9.82
!            eos(7) is dlogP/dlogT6 at conxtant Rho.
!                   Cox and Guil1 eq 9.81
!            eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
!            eos(9) is gamma2/(gamma2-1). Eqs. 9.88 Cox and Guili

!            iorder sets maximum index for eos(i);i.e., iorder=1
!                   gives just the pressure

!            irad  if =0 no radiation correction; if =1 adds radiation

!            index(i),i=1,10  sets order in which the equation of state
!            variables are stored in eos(i).  Above order corresponds
!            to block data statement:
!                 data (index(i),i=1,10)/1,2,3,4,5,6,7,8,9,10/.
!            If you, for example, only want to return gamma1: set iorder=1
!            and set: data (index(i),i=1,10)/8,2,3,4,5,6,7,1,9,10/

      save
      double precision :: results_out(iorder), pgas_out, prad_out
      integer :: info ! returned = 0 if AOK
      character (len=*) :: filename
      integer :: w
      real :: moles
      parameter (mx=6,mv=12,nr=169,nt=197)
      character :: blank*1
      common/lreadco_hhe/itime
      common/eeeos_hhe/ epl(mx,nt,nr),xx(mx)
      common/aaeos_hhe/ q(4),h(4),xxh
      common/aeos_hhe/  xz(mx,mv,nt,nr), t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx),dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos_hhe/ iri(mv),index(mv),nta(nr),zz(mx)
      common/bbeos_hhe/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/eeos_hhe/esact,eos(mv)
      dimension frac(7)
      data aprop/83.14511/

      pgas_out = -1
      prad_out = -1

      blank=' '
      ztab=0.0
      if (iorder < 1 .or. iorder > mv) then
         write(*,*) 'invalid iorder argument (must be between 1 and ', mv, ')', iorder
         stop 1
      end if

      if ((irad /= 0) .and. (irad /= 1)) then
         write (*,'(" Irad must be 0 or 1")')
         stop
      end if

      info = 0
      xxi=xh
      ri=r

      slt=t6
      slr=r

      if(itime /= 12345678) then
        itime=12345678
        do i=1,mv
          do j=1,mv
            if  (index(i) == j) iri(i)=j
          end do
        end do
        do  i=1,mx
          xx(i)=xa(i)
        end do

!..... read the data files
        call readcoeos_hhe(filename)
        z=zz(1)
        if (abs(ztab-z) > 0.00001) then
        write (*,'("requested z=",f10.6," EOS5_data_H-He is for z=", f10.6)') ztab,z
        stop
        end if

        if(z+xh-1.e-6 > 1 ) GOTO 61
      end if


!..... Determine T6,rho grid points to use in the
!      interpolation.
      if((slt > t6a(1)).or.(slt < t6a(nt))) GOTO 62
      if((slr < rho(1)).or.(slr > rho(nr))) GOTO 62


!     if (izi == 0) then  ! freeze table indices if not 0
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
        if (xh < 1.e-6) then
        mf=1
        mg=1
        mh=1
        mi=2
        mf2=1
        end if
        if((xh <= xa(2)+1.e-7) .or. (xh >= xa(mx-2)-1.e-7)) mf2=mh

        ilo=2
        ihi=nr
   12     if(ihi-ilo > 1) then
          imd=(ihi+ilo)/2
           if (slr == rho(imd)) then
           ihi=imd
           GOTO 13
           end if
            if(slr <= rho(imd)) then
              ihi=imd
            else
              ilo=imd
            end if
          GOTO 12
          end if
   13     i=ihi
        l1=i-2
        l2=i-1
        l3=i
        l4=l3+1
        iqu=3
        if(l4 > nr) iqu=2

        ilo=nt
        ihi=2
   11     if(ilo-ihi > 1) then
          imd=(ihi+ilo)/2
           if (t6 == t6list(1,imd)) then
           ilo=imd
           GOTO 14
           end if
            if(t6 <= t6list(1,imd)) then
              ihi=imd
            else
              ilo=imd
            end if
          GOTO 11
          end if
   14     i=ilo
        k1=i-2
        k2=i-1
        k3=i
        k4=k3+1
        ipu=3
        if (k4 > nt) ipu=2
      if (k3 == 0) then
      write (*,'(" ihi,ilo,imd",3i5)')
      end if

!     end if

!     check to determine if interpolation indices fall within
!     table boundaries.  choose largest allowed size.
      sum1=0.0
      sum2=0.0
      sum23=0.0
      sum33=0.0
      do m=mf,mf2
        do ir=l1,l1+1
          do it=k1,k1+1
            sum1=sum1+xz(m,1,it,ir)
          end do
        end do
        do ir=l1,l1+2
          do it=k1,k1+2
            sum2=sum2+xz(m,1,it,ir)
          end do
        end do
        if (ipu == 3) then
          do ir=l1,l1+2
            do it=k1,k1+ipu
              sum23=sum23+xz(m,1,it,ir)
            end do
          end do
        else
          sum23=2.e+30
        end if
        if (iqu == 3) then
          do ir=l1,l1+3
            do it=k1,k1+ipu
              sum33=sum33+xz(m,1,it,ir)
            end do
          end do
        else
          sum33=2.e+30
        end if
      end do
      iq=2
      ip=2
      if (sum2 > 1.e+30) then
        if (sum1 < 1.e+25 ) then
          k1=k3-3
          k2=k1+1
          k3=k2+1
          l1=l3-3
          l2=l1+1
          l3=l2+1
          GOTO 15
            else
          GOTO 65
        end if
      end if
      if (sum23 < 1.e+30) ip=3
      if (sum33 < 1.e+30) iq=3

      if(t6 >= t6list(1,2)+1.e-7) ip=2
      if(slr <= rho(2)+1.e-7) iq=2

      if((l3 == nr) .or. (k3 == nt)) then
         iq=2
         ip=2
      end if

   15 continue
      do iv=1,iorder
      do m=mf,mf2

      is=0

!__________
      do ir=l1,l1+iq
        do it=k1,k1+ip
        epl(m,it,ir)=xz(m,iv,it,ir)
        is=1
        end do
      end do
      end do
      if((zz(mg) /= zz(mf)) .or. (zz(mh) /= zz(mf))) then
        write(*,'("Z does not match Z in EOS5_data_H-He files you are using")')
        stop
      end if
      if(z /= zz(mf)) GOTO 66
      is=0
      iw=1
      do ir=l1,l1+iq
        do it=k1,k1+ip
          if (mf2 == 1) then
          esk(it,ir)=epl(mf,it,ir)
          GOTO 46
          end if
          esk(it,ir)=quadeos(is,iw,xh,epl(mf,it,ir),epl(mg,it,ir),epl(mh,it,ir),xx(mf),xx(mg),xx(mh))
          if(esk(it,ir) > 1.e+20) then
          write(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir,l3,k3,iq,ip
          write(*,'(3e12.4)')  (epl(ms,it,ir),ms=mf,mf+2)
          end if
          is=1
   46     continue
        end do
      end do

      if (mi == mf2) then  ! interpolate between quadratics
      is=0
      iw=1
       dixr=(xx(mh)-xh)*dfsx(mh)
      do ir=l1,l1+iq
        do it=k1,k1+ip
          esk2(it,ir)=quadeos(is,iw,xh,epl(mg,it,ir),epl(mh,it,ir),epl(mi,it,ir),xx(mg),xx(mh),xx(mi))
          if(esk(it,ir) > 1.e+20) then
          write(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir,l3,k3,iq,ip
          write(*,'(3e12.4)')  (epl(ms,it,ir),ms=mg,mg+2)
          end if
          esk(it,ir)=esk(it,ir)*dixr+esk2(it,ir)*(1.-dixr)
          is=1
        end do
      end do


      end if

      is=0

!..... completed X interpolation. Now interpolate T6 and rho on a
!      4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
!      mixes overlapping quadratics to obtain smoothed derivatives.


      call t6rinteos_hhe(slr,slt)
      eos(iv)=esact
      end do

      p0=t6*r
      eos(iri(1))=eos(iri(1))*p0   ! interpolated in p/po
      eos(iri(2))=eos(iri(2))*t6   ! interpolated in E/T6
      tmass=gmass(xh,z,moles,eground,fracz,frac)
      pgas_out = eos(iri(1)) * 1d12
      prad_out = 1d12*4d0/3d0*(1.8914785d-3)*t6**4
      if (irad == 1) then
         call radsub_hhe (irad,t6,r,moles,tmass)
      else
         eos(iri(5))=eos(iri(5))*moles*aprop/tmass
      end if

      do j=1,iorder
         results_out(j) = eos(j)
      end do

      return

   61 continue
      info = -61; return
      write(*,'(" Mass fractions exceed unity (61)")')
      stop 1
   62 continue
      info = -62; return
      write(*,'(" T6/LogR outside of table range (62)")')
      write(*,*) 'T6', t6, 'LogR', log10(r), 'LogT', log10(t6*1d6)
      stop 1
   65 continue
      info = -65; return
      write(*,'("T6/log rho in empty region od table (65)")')
      write (*,'("xh,t6,r=", 3e12.4)') xh,t6,r
      stop 1
   66 continue
      info = -66; return
      write(*,'(" Z does not match Z in EOS5_data* files you are using (66)")')
      write (*,'("mf,zz(mf)=",i5,e12.4)') mf,zz(mf)
      write (*,'("  iq,ip,k3,l3,xh,t6,r,z= ",4i5,4e12.4)') ip,iq,k3,l3,xh,t6,r,z
      stop 1
      return
      end subroutine esac_hhe

! **********************************************************************
      subroutine t6rinteos_hhe(slr,slt)
!     The purpose of this subroutine is to interpolate in T6 and rho
      save
      parameter (mx=6,mv=12,nr=169,nt=197)
      common/eeeos_hhe/ epl(mx,nt,nr),xx(mx)
      common/aaeos_hhe/ q(4),h(4),xxh
      common/aeos_hhe/ xz(mx,mv,nt,nr),t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx),dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/bbeos_hhe/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/eeos_hhe/esact,eos(mv)

      iu=0
      is=0

      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quadeos(is,iw,slr,esk(kx,l1),esk(kx,l2),esk(kx,l3),rho(l1),rho(l2),rho(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quadeos(is,iw,slr,esk(kx,l2),esk(kx,l3),esk(kx,l4),rho(l2),rho(l3),rho(l4))
          end if
        is=1
      end do

      is=0
      iw=1
!..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      esact=quadeos(is,iw,slt,h(1),h(2),h(3),t6a(k1),t6a(k2),t6a(k3))
        if (iq. eq. 3) then
!.....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          esactq=quadeos(is,iw,slt,q(1),q(2),q(3),t6a(k1),t6a(k2),t6a(k3))
        end if
        if(ip == 3) then
!.....    eos(i) in lower-left 3x3.
          esact2=quadeos(is,iw,slt,h(2),h(3),h(4),t6a(k2),t6a(k3),t6a(k4))
!.....    eos(i) smoothed in left 3x4
          dix=(t6a(k3)-slt)*dfs(k3)
          esact=esact*dix+esact2*(1.-dix)
!       end if   ! moved to loc a
        if(iq == 3) then

!.....     eos(i) in upper-right 3x3.
          esactq2=quadeos(is,iw,slt,q(2),q(3),q(4),t6a(k2),t6a(k3),t6a(k4))
          esactq=esactq*dix+esactq2*(1.-dix)
        end if
        end if  ! loc a

        if(iq == 3) then
          dix2=(rho(l3)-slr)*dfsr(l3)
            if(ip == 3) then
!.....        eos(i) smoothed in both log(T6) and log(R)
              esact=esact*dix2+esactq*(1-dix2)
            end if
        end if
        if (esact > 1.e+15) then
          write(*,'("Interpolation indices out of range; please report conditions.")')
          stop
        end if

      return
      end subroutine t6rinteos_hhe


! *********************************************************************
      subroutine readcoeos_hhe(filename)
!..... The purpose of this subroutine is to read the data tables
      save
      character (len=256) :: filename
      parameter (mx=6,mv=12,nr=169,nt=197)
      real :: moles
      character(len=1) :: blank
      common/aaeos_hhe/ q(4),h(4),xxh
      common/aeos_hhe/  xz(mx,mv,nt,nr),t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx),dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos_hhe/ iri(mv),index(mv),nta(nr),zz(mx)
      common/eeos_hhe/esact,eos(mv)
      common/eeeos_hhe/ epl(mx,nt,nr),xx(mx)
      common/eeeeos_hhe/moles(mx),xin(mx),tmass(mx),icycuse(mx,nr),amu_M(mx,nr),alogNe(mx,nr),rhogr(mx,nr),frac(mx,6),alogr(nr,nt)

      blank=' '

        if (itimeco /= 12345678) then
        do i=1,mx
          do j=1,mv
            do k=1,nt
              do l=1,nr
                xz(i,j,k,l)=1.e+35
              end do
            end do
          end do
        end do
        itimeco=12345678
        end if

      close (2)
!..... read  tables
!       call system (' gunzip EOS5_data_H-He')
!       open(2, FILE='EOS5_data_H-He')

      write(*,*) 'read ' // trim(filename)

       open(2, FILE=trim(filename), IOSTAT=ios)
      if (ios /= 0) then
         write(*,*) 'failed to open ', trim(filename)
         stop 1
      end if


      do m=1,mx
      read (2,'(3x,f6.4,3x,f12.9,11x,f10.7,17x,f10.7)') xin(m),zz(m),moles(m),tmass(m)
      read (2,'(21x,e14.7,4x,e14.7,3x,e11.4,3x,e11.4,3x,e11.4,4x,e11.4)') (frac(m,i),i=1,6)
      read (2,'(a)') blank
      do jcs=1,nr
      read (2,'(2i5,2f12.7,17x,e15.7)') numtot,icycuse(m,jcs),dum,dum,rhogr(m,jcs)
      if(numtot /= jcs) then
         write (*,'(" Data file incorrect: numtot,jcs= ",2i5)') numtot,jcs
         stop
      end if
      read(2,'(a)') blank
      read(2,'(a)') blank
      if (icycuse(m,jcs) < nta(jcs)) then
!         write (*,'("problem with data files: X=",f6.4," density=",e14.4)') xin(m), rhogr(m,jcs)
         write (*,*) 'problem with ', filename
         write (*,*) 'm', m
         write (*,*) 'jcs', jcs
         write (*,*) 'icycuse(m,jcs)', icycuse(m,jcs)
         write (*,*) 'nta(jcs)', nta(jcs)
         write (*,*) 'readcoeos_hhe'
         stop 1
      end if
      do  i=1,icycuse(m,jcs)
      if (i > nta(jcs)) then
         read (2,'(a)') blank
         GOTO 4
      end if
      ! change to read mu_M and logNe as part of xz table
      read (2,'(f11.6,1x,f6.4,e11.4,2e13.6,2e11.3,5f10.6)')
     &         t6list(jcs,i),(xz(m,index(iv),i,jcs),iv=10,11),(xz(m,index(iv),i,jcs),iv=1,9)
!       t6list(jcs,i),amu_M(jcs,i),alogNe(jcs,j),
!       (xz(m,index(iv),i,jcs),iv=1,9)
    4 continue
      end do
      read(2,'(a)') blank
      read(2,'(a)') blank
      read(2,'(a)') blank
      end do
      read(2,'(a)') blank
       end do

      do i=1,nt
         if(t6list(1,i) == 0.0) then
            write(*,'("readcoeos_hhe: Error:",i4,"-th T6 value is zero")') i
            stop
         end if
         t6a(i)=t6list(1,i)
      end do
      do i=2,nt
         dfs(i)=1./(t6a(i)-t6a(i-1))
      end do
      rho(1)=rhogr(1,1)
      do i=2,nr
         rho(i)=rhogr(1,i)
         dfsr(i)=1./(rho(i)-rho(i-1))
      end do
      do i=2,mx
         dfsx(i)=1./(xx(i)-xx(i-1))
      end do
      !call system ('gzip EOS5_data_H-He')
      end subroutine readcoeos_hhe


! share quadeos and gmass with the standard xtrin

! **********************************************************************
      subroutine radsub_hhe (irad,t6,density,moles,tmass)
      parameter (mx=6,mv=12,nr=169,nt=197)
      real :: moles,k,molenak,Na
      common/eeos_hhe/esact,eos(mv)
      common/beos_hhe/ iri(mv),index(mv),nta(nr),zz(mx)

      data Na/6.0221367e+23/, k/1.380658e-16/, unitf/0.9648530/, unitfold/0.965296/,
     &     c/2.9979245e+10/, sigma/5.67051e-5/, sigmac/1.8914785e-15/, sigmacc/1.8914785e-3/, aprop/83.14510/

! Physical constants
!       Na=6.0221367e+23
!       k =1.380658e-16 !   erg/degree K
!       Na*k=6.0221367E+23*1.380658e-16 erg/degree K=8.314510E+7 erg/degree K
!           =8.314510e+7*11604.5 erg/eV=0.9648575E+12 erg/eV
!           Define unitf= Na*k/e+12=0.9648575
!           unitf=0.9648530  ! obtained with latest physical constants
!           unitfold=0.965296   ! old units- still used in the EOS code
!           In these units energy/density is in units of Mb-CC/gm
!           Pressure is in units of E+12 bars=Mb
!       sigma is the Stefan-Boltzmann constant
!       sigma=5.67051E-5 !   erg /(s*cm**2*K**4)
!       c=2.99792458E+10 !   cm/s

!     rat=sigma/!    ! dyne/(cm**2*K**4)

!     rat=rat*1.e+24  !   Convert degrees K to units 10**6 K (T replaced with T6)
      rat=sigmacc
      molenak=moles*aprop  ! Mb-cc/unit T6

!-----Calculate EOS without radiation correction
      pt=eos(iri(1))
      et=eos(iri(2))
      st=eos(iri(3))
      dedrho=eos(iri(4))
      chir=eos(iri(6))*eos(iri(1))/pt
      chitt=(eos(iri(1))*eos(iri(7)))/pt
      cvtt=(eos(iri(5))*molenak/tmass)
      gam3pt_norad=pt*chitt/(cvtt*density*t6)
      gam1t_norad=chir+chitt*gam3pt_norad
      gam2pt_norad=gam1t_norad/gam3pt_norad
!---- End  no radiation calculation

      if (irad /= 0) then
!---- Calculate EOS with radiation calculation
      pr=4./3.*rat*t6**4   ! Mb
      er=3.*pr/density   ! Mb-cc/gm
      sr=4./3.*er/t6   ! Mb-cc/(gm-unit T6)
      pt=eos(iri(1))+pr
      et=eos(iri(2))+er
      st=eos(iri(3))+sr
      dedrho=eos(iri(4))-er/density
      chir=eos(iri(6))*eos(iri(1))/pt
      chitt=(eos(iri(1))*eos(iri(7))+4.*pr)/pt
!     gam1t(jcs,i)=(p(jcs,i)*gam1(jcs,i)+4./3.*pr)/pt(jcs,i)
!     gam2pt(jcs,i)=(gam2p(jcs,i)*p(jcs,i)+4.*pr)/pt(jcs,i)
!     gam3pt(jcs,i)=gam1t(jcs,i)/gam2pt(jcs,i)
      cvtt=(eos(iri(5))*molenak/tmass+4.*er/t6)
      gam3pt=pt*chitt/(cvtt*density*t6)                              ! DIRECT
      gam1t=chir+chitt*gam3pt !eq 16.16 Landau_Lifshitz (Stat. Mech) ! DIRECT
      gam2pt=gam1t/gam3pt                                            ! DIRECT

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6) /molenak
!-----Add difference between EOS with and without radiation.  cvtt
!       calculation is not accurate enough to give accurate results using
!       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
      eos(iri(8))=eos(iri(8))+gam1t-gam1t_norad
      eos(iri(9))=eos(iri(9))+gam2pt-gam2pt_norad
!     eos(iri(mv))=eos(iri(mv))+gam3pt-gam3pt_norad
      end if
!-----End EOS calculations with radiation
!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)  /molenak
      eos(iri(1))=pt
      eos(iri(2))=et
      eos(iri(3))=st
      eos(iri(4))=dedrho
      eos(iri(5))=cvtt
      eos(iri(6))=chir
      eos(iri(7))=chitt
      return
      end subroutine radsub_hhe
! **********************************************************************
      block data

      ! for EOS5_xtrin_H_He
      parameter (mx_hhe=6,mv_hhe=12,nr_hhe=169,nt_hhe=197)
      common/aaeos_hhe/ q_hhe(4),h_hhe(4),xxh_hhe
      common/aeos_hhe/  xz_hhe(mx_hhe,mv_hhe,nt_hhe,nr_hhe),t6list_hhe(nr_hhe,nt_hhe),rho_hhe(nr_hhe),
     &                  t6a_hhe(nt_hhe),esk_hhe(nt_hhe,nr_hhe),esk2_hhe(nt_hhe,nr_hhe),
     &                  dfsx_hhe(mx_hhe),dfs_hhe(nt_hhe),dfsr_hhe(nr_hhe),m_hhe,mf_hhe,xa_hhe(mx_hhe)
      common/beos_hhe/ iri_hhe(mv_hhe),index_hhe(mv_hhe),nta_hhe(nr_hhe),zz_hhe(mx_hhe)

      data (xa_hhe(i),i=1,mx_hhe)/0.0,0.2,0.4,0.6,0.8,1.0/

      data (index_hhe(i),i=1,mv_hhe)/1,2,3,4,5,6,7,8,9,10,11,12/
      data (nta_hhe(i),i=1,nr_hhe)/92*197,196,194,190,2*189,3*187,2*159,133,125,123,122,120,115,113,
     &                             107,102,95,87,83,74,69,67,62,56,54,52,51,
     &                             3*50,2*49,26*48,47,45,43,42,41,39,37,36,35,34,32,31,30,29,27,26/


      ! for opal_core
      parameter (mx=5,mv=12,nr=169,nt=197)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx),dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/ iri(mv),index(mv),nta(nr),zz(mx)
      data (xa(i),i=1,mx)/0.0,0.2,0.4,0.6,0.8/

      data (index(i),i=1,mv)/1,2,3,4,5,6,7,8,9,10,11,12/
      data (nta(i),i=1,nr)/87*197,7*191,190,2*189,185,179,170,2*149,133,125,123,122,120,115,113,107,102,
     &                     2*80,72,68,66,64,62,56,54,52,51,2*50,49,47,2*45,43,42,28*40,39,37,36,35,34,32,31,30,29,27,26/

      end
