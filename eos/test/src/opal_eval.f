

      subroutine test_get_pure_opal
         double precision ::  
     >         Z, X, rho, temp,
     >         logPgas, logE, logS, chiRho, chiT, 
     >         Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >         mu, logNe, gamma1, gamma3, grad_ad
         character (len=256) :: data_dir
         integer :: ierr
         logical :: include_radiation
         
         include 'formats.dek'
         
         data_dir = '../eosDT_builder/eos_input_data'
         Z = 0.02d0
         X = 0.72d0
         rho = 1d-3
         temp = 1d5
         include_radiation = .true.
         ierr = 0
         call get_pure_opal(
     >         Z, X, rho, temp, data_dir, include_radiation,
     >         logPgas, logE, logS, chiRho, chiT, 
     >         Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >         mu, logNe, gamma1, gamma3, grad_ad,
     >         ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_pure_opal'
            stop 1
         end if
         
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'rho', rho
         write(*,1) 'temp', temp
         write(*,*)
         write(*,1) 'logPgas', logPgas
         write(*,1) 'logE', logE
         write(*,1) 'logS', logS
         write(*,1) 'chiRho', chiRho
         write(*,1) 'chiT', chiT
         write(*,1) 'Cp', Cp
         write(*,1) 'Cv', Cv
         write(*,1) 'dE_dRho', dE_dRho
         write(*,1) 'dS_dT', dS_dT
         write(*,1) 'dS_dRho', dS_dRho
         write(*,1) 'mu', mu
         write(*,1) 'logNe', logNe
         write(*,1) 'gamma1', gamma1
         write(*,1) 'gamma3', gamma3
         write(*,1) 'grad_ad', grad_ad
         write(*,*)
         
      
      end subroutine test_get_pure_opal


      subroutine get_pure_opal(
     >         Z, X, den, temp, data_dir, include_radiation,
     >         logPgas, logE, logS, chiRho, chiT, 
     >         Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >         mu, logNe, gamma1, gamma3, grad_ad,
     >         ierr)
         use const_def
         use math_lib
         double precision, intent(in) :: Z ! the desired Z
         double precision, intent(in) :: X ! the desired X         
         double precision, intent(in) :: den ! the density            
         double precision, intent(in) :: temp ! the temperature
         character (len=*), intent(in) :: data_dir
         logical, intent(in) :: include_radiation
         double precision, intent(out) ::  
     >         logPgas, logE, logS, chiRho, chiT, 
     >         Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >         mu, logNe, gamma1, gamma3, grad_ad
         integer, intent(out) :: ierr
         
         integer, parameter :: num_results = 11
         double precision :: opal_results(num_results)
         double precision t6, r 
         integer :: irad
         double precision :: pout, dpoutdd, dpoutdt, eout, deoutdd,
     >      sout, dsoutdd, dsoutdt, pgas_out, prad_out

         double precision log_free_e, log_free_e0, log_free_e1, Prad
         
         include 'formats.dek'
         
         t6 = temp * 1d-6
         r = den
         if (include_radiation) then
            irad = 1
         else
            irad = 0
         end if
         call opal_eos (X, Z, t6, r, num_results, 
     >      irad, opal_results, pgas_out, prad_out, data_dir, ierr)
         if (ierr /= 0) then
            !return
            
            write(*,*) 'opal_eos returned ierr', ierr
            write(*,*) ' opal_only', opal_only
            write(*,*) ' alfa', alfa
            write(*,*) ' beta', beta
            write(*,*) ' X', X
            write(*,*) ' Z', Z
            write(*,*) 't6', t6
            write(*,*) ' r', r
            stop 1
         end if
         
         pout    = opal_results(1) * 1d12
         dpoutdd = opal_results(6) * pout / den
         dpoutdt = opal_results(7) * pout / temp
   
         eout    = opal_results(2) * 1d12
         ! NOTE: rather than use the opal tabulated value for cv,
         ! it is better to calculate cv from other things.
         deoutdd = opal_results(4) * 1d12
   
         sout    = opal_results(3) * 1d6
         dsoutdd = (deoutdd - pout / (den*den)) / temp
         dsoutdt = opal_results(5) * 1d6 / temp

         mu = opal_results(10)
         logNe = opal_results(11)
         
         if (include_radiation .and. Prad >= pout) then
            write(*,*) 'Prad >= pout'
            stop 'opal_scvh_driver'
         end if
         
         logPgas = log10(pgas_out)
         logE = log10(eout)
         logS = log10(sout)
         logRho = log10(den)
         chiRho = dpoutdd*den/pout
         chiT = dpoutdt*temp/pout
         gamma1 = opal_results(8)
         ! gam1 / (gam3 - 1) = gam2 / (gam2 - 1)
         ! gam2 / (gam2 - 1) = opal_results(9)
         gamma3 = 1 + opal_results(8)/opal_results(9)
         Cv = chiT * pout / (den * temp * opal_results(8)/opal_results(9)) ! C&G 9.93
         Cp = Cv*gamma1/chiRho
         grad_ad = 1/opal_results(9)
         dE_dRho = deoutdd
         dS_dRho = dsoutdd
         dS_dT = dsoutdt
         log_free_e = logNe - log10(avo) - logRho


      end subroutine get_pure_opal



!     based on EOS5_xtrin.f from the OPAL release
!..***********************************************************************
      subroutine opal_eos (xh_in,ztab_in,t6_in,r_in,iorder,irad,results_out,pgas_out,prad_out,data_dir,info)
!        xh=hydrogen mass fraction
!        ztab is metal fraction of the EOS5_data tables you are using.
!           included only for purpose of preventing mismatch
!        t6=T6=temperature in millions of degrees kelvin
!        r=rho=Rho=density(g/cm**3)
!        iorder=11  ! gives all 1st and 2nd order data. See instructions
!                  in esac.
!        irad=0    ! if 1, then adds radiation corrections
!        results   ! array of length iorder to hold the results
!            results(1) is the pressure in megabars (10**12dyne/cm**2)
!            results(2) is energy in 10**12 ergs/gm. Zero is zero T6
!            results(3) is the entropy in units of energy/T6
!            results(4) is dE/dRHO at constant T6
!            results(5) is the specific heat, dE/dT6 at constant V.
!            results(6) is dlogP/dlogRho at constant T6. 
!                   Cox and Guil1 eq 9.82
!            results(7) is dlogP/dlogT6 at constant Rho.
!                   Cox and Guil1 eq 9.81
!            results(8) is gamma1. Eqs. 9.88 Cox and Guili.
!            results(9) is gamma2/(gamma2-1). Eqs. 9.88 Cox and Guili
!            results(10) is mu_M
!            results(11) is log_Ne
      integer info ! returned = 0 if AOK
      character (len=256), intent(in) :: data_dir
      parameter (mx=5,mv=12,nr=169,nt=197)
      common/eeos/esact,eos(mv)
      double precision xh_in,ztab_in,t6_in,r_in,results_out(iorder),pgas_out,prad_out
      real xh,ztab,t6,r
      character (len=256) :: filename
      integer iorder, irad
      if (iorder .lt. 1 .or. iorder .gt. mv) then
         write(*,*) 'invalid iorder argument (must be between 1 and ', mv, ')', iorder
         stop 1
      endif
      xh = xh_in
      ztab = ztab_in
      t6 = t6_in
      r = r_in
      if (ztab == 0) then
         filename = trim(data_dir) // '/opal/EOS5_data_H-He'
      else if (ztab < 0.005d0) then
         filename = trim(data_dir) // '/opal/EOS5_data_00z'
      else if (ztab < 0.1d0) then
         write(filename,'(a,a,i1,a)') trim(data_dir), '/opal/EOS5_data_0', floor(100d0 * ztab + 0.5d0), 'z'
      else
         write(filename,'(a,a,i2,a)') trim(data_dir), '/opal/EOS5_data_', floor(100d0 * ztab + 0.5d0), 'z'
      end if
      
      if (ztab == 0) then
         write(*,*) 'Z=0.00 not yet supported'
         stop 1
         !call esac_hhe(xh,t6,r,iorder,irad,filename,results_out,pgas_out,prad_out,info)
         return
      end if
      
      call esac(xh,ztab,t6,r,iorder,irad,filename,pgas_out,prad_out,info)
      if (info /= 0) return
      do j=1,iorder
         results_out(j) = eos(j)
      end do
      end subroutine opal_eos
      
      
!..***********************************************************************
      subroutine esac (xh,ztab,t6,r,iorder,irad,filename,pgas_out,prad_out,info)
!..... The purpose of this subroutine is to interpolate 
!      the equation of state and its derivatives in X, T6, density

!        xh=hydrogen mass fraction
!        ztab is metal fraction of the EOS5_data tables you are using.
!           included only for purpose of preventing mismatch
!        t6=T6=temperature in millions of degrees kelvin
!        r=rho=Rho=density(g/cm**3)
!..... to use esac insert common/eeos/ esact,eos(11) in the calling routine.
!      This common contains the interpolated EOS values for the EOS
!
!..... eos(i) are obtained from a quadradic interpolation at
!      fixed T6 at three values of Rho; followed by quadratic
!      interpolation along T6. Results smoothed by mixing
!      overlapping quadratics.

!            eos(1) is the pressure in megabars (10**12dyne/cm**2)
!            eos(2) is energy in 10**12 ergs/gm. Zero is zero T6
!            eos(3) is the entropy in units of energy/T6
!            eos(4) is dE/dRHO at constant T6
!            eos(5) is the specific heat, dE/dT6 at constant V.
!            eos(6) is dlogP/dlogRho at constant T6. 
!                   Cox and Guil1 eq 9.82
!            eos(7) is dlogP/dlogT6 at constant Rho.
!                   Cox and Guil1 eq 9.81
!            eos(8) is gamma1. Eqs. 9.88 Cox and Guili.
!            eos(9) is gamma2/(gamma2-1). Eqs. 9.88 Cox and Guili
!            eos(10) is mu_M
!            eos(11) is log_Ne

!            iorder sets maximum index for eos(i);i.e., iorder=1
!                   gives just the pressure
!
!            irad  if =0 no radiation correction; if =1 adds radiation

!            index(i),i=1,mv  sets order in which the equation of state
!            variables are stored in eos(i).  Above order corresponds
!            to block data statement:
!                 data (index(i),i=1,11)/1,2,3,4,5,6,7,8,9,10,11/. 
!            If you, for example, only want to return gamma1: set iorder=1 
!            and set: data (index(i),i=1,11)/8,2,3,4,5,6,7,1,9,10,11/
!
!
      save
      double precision pgas_out, prad_out
      integer info ! returned = 0 if AOK
      character (len=*) filename
      real moles
      parameter (mx=5,mv=12,nr=169,nt=197)
      character blank*1
      common/lreadco/itime
      common/eeeos/ epl(mx,nt,nr),xx(mx)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/ iri(mv),index(mv),nta(nr),zz(mx)
      common/bbeos/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/eeos/esact,eos(mv)
      dimension frac(7)
      data aprop/83.14511/

      info = 0
      blank=' '
      if (iorder .gt. mv ) then
         write (*,*) ' iorder cannot exceed ', mv
      endif
      if ((irad .ne. 0) .and. (irad .ne. 1)) then
         write (*,*) ' Irad must be 0 or 1'
         stop 1
      endif

      xxi=xh
      ri=r
!
      slt=t6
      slr=r
!
      if(itime .ne. 12345678) then
        itime=12345678
        do i=1,mv
          do j=1,mv
            if  (index(i) .eq. j) iri(i)=j
          enddo
        enddo
        do  i=1,mx
          xx(i)=xa(i)
        enddo
!
!..... read the data files
        call readcoeos(filename)
        z=zz(1)
        if (ztab .ne. z) then
            write (*,'("requested z=",f10.6," EOS5_data is for z=",
     x             f10.6)')
     x         ztab,z      
            stop 1
        endif

        if(z+xh-1.e-6 .gt. 1 ) go to 61
      endif
!
!
!..... Determine T6,rho grid points to use in the
!      interpolation.
      if((slt .gt. t6a(1)).or.(slt .lt. t6a(nt))) go to 62
      if((slr .lt. rho(1)).or.(slr .gt. rho(nr))) go to 62
!
!
!
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
           mf=1
           mg=1
           mh=1
           mi=2
           mf2=1
        endif
        if((xh .le. xa(2)+1.e-7) .or. (xh .ge. xa(mx-2)-1.e-7)) mf2=mh
!
        ilo=2
        ihi=nr
   12     if(ihi-ilo .gt. 1) then
              imd=(ihi+ilo)/2
              if (slr .eq. rho(imd)) then
                 ihi=imd
                 go to 13
              endif
               if(slr .le. rho(imd)) then
                 ihi=imd
               else
                 ilo=imd
               endif
             go to 12
          endif
   13     i=ihi
        l1=i-2
        l2=i-1
        l3=i
        l4=l3+1
        iqu=3
        if(l4 .gt. nr) iqu=2
!
        ilo=nt
        ihi=2
   11     if(ilo-ihi .gt. 1) then
          imd=(ihi+ilo)/2
           if (t6 .eq. t6list(1,imd)) then
           ilo=imd
           go to 14
           endif
            if(t6 .le. t6list(1,imd)) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 11
          endif
   14     i=ilo
        k1=i-2
        k2=i-1
        k3=i
        k4=k3+1
        ipu=3
        if (k4 .gt. nt) ipu=2
      if (k3 .eq. 0) then
            write (*,'(" ihi,ilo,imd",3i5)')
      endif

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
          enddo
        enddo
        do ir=l1,l1+2
          do it=k1,k1+2
            sum2=sum2+xz(m,1,it,ir)
          enddo
        enddo
        if (ipu .eq. 3) then
          do ir=l1,l1+2
            do it=k1,k1+ipu
              sum23=sum23+xz(m,1,it,ir)
            enddo
          enddo
        else
          sum23=2.e+30
        endif
        if (iqu .eq. 3) then
          do ir=l1,l1+3
            do it=k1,k1+ipu
              sum33=sum33+xz(m,1,it,ir)
            enddo
          enddo
        else
          sum33=2.e+30
        endif
      enddo
      iq=2
      ip=2
      if (sum2 .gt. 1.e+30) then
        if (sum1 .lt. 1.e+25 ) then
          k1=k3-3
          k2=k1+1
          k3=k2+1
          l1=l3-3
          l2=l1+1
          l3=l2+1
          go to 15
        else
          go to 65
        endif
      endif
      if (sum23 .lt. 1.e+30) ip=3
      if (sum33 .lt. 1.e+30) iq=3

      if(t6 .ge. t6list(1,2)+1.e-7) ip=2
      if(slr .le. rho(2)+1.e-7) iq=2

      if((l3 .eq. nr) .or. (k3 .eq. nt)) then
         iq=2
         ip=2
      endif

   15 continue
      do 124 iv=1,iorder
      do 123 m=mf,mf2
         is=0
         do ir=l1,l1+iq
           do it=k1,k1+ip
              epl(m,it,ir)=xz(m,iv,it,ir)
              is=1
           enddo
         enddo
  123 continue
      if((zz(mg) .ne. zz(mf)) .or. (zz(mh) .ne. zz(mf))) then
        write(*,'("Z does not match Z in EOS5_data files you are"
     x ," using")')
        stop 1
      endif
      if(z .ne. zz(mf)) go to 66
      is=0
      iw=1
      do 45 ir=l1,l1+iq
        do it=k1,k1+ip
          if (mf2 .eq. 1) then
            esk(it,ir)=epl(mf,it,ir)
            go to 46
          endif
          esk(it,ir)=quadeos(is,iw,xh,epl(mf,it,ir),epl(mg,it,ir)
     x      ,epl(mh,it,ir),xx(mf),xx(mg),xx(mh))
          if(esk(it,ir) .gt. 1.e+20) then
            write(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir
     x         ,l3,k3,iq,ip
            write(*,'(3e12.4)')  (epl(ms,it,ir),ms=mf,mf+2)
          endif
          is=1
   46     continue
        enddo
   45 continue

      if (mi .eq. mf2) then  ! interpolate between quadratics
         is=0
         iw=1
         dixr=(xx(mh)-xh)*dfsx(mh)
         do 47 ir=l1,l1+iq
           do it=k1,k1+ip
             esk2(it,ir)=quadeos(is,iw,xh,epl(mg,it,ir),epl(mh,it,ir)
     x         ,epl(mi,it,ir),xx(mg),xx(mh),xx(mi))
             if(esk(it,ir) .gt. 1.e+20) then
             write(*,'(" problem it ir,l3,k3,iq,ip=", 6i5)') it,ir
     x         ,l3,k3,iq,ip
             write(*,'(3e12.4)')  (epl(ms,it,ir),ms=mg,mg+2)
             endif
             esk(it,ir)=esk(it,ir)*dixr+esk2(it,ir)*(1.-dixr)
             is=1
           enddo
   47 continue
   

      endif

      is=0
!
!..... completed X interpolation. Now interpolate T6 and rho on a
!      4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
!      mixes overlapping quadratics to obtain smoothed derivatives.
!
!
      call t6rinteos(slr,slt)
      eos(iv)=esact
  124 continue

      p0=t6*r
      eos(iri(1))=eos(iri(1))*p0   ! interpolated in p/po
      eos(iri(2))=eos(iri(2))*t6   ! interpolated in E/T6
      tmass=gmass(xh,z,moles,eground,fracz,frac)
      pgas_out = eos(iri(1)) * 1d12
      prad_out = 1d12*4d0/3d0*(1.8914785d-3)*t6*t6*t6*t6
      if (irad .eq. 1) then
         call radsub (irad,t6,r,moles,tmass)
      else
         eos(iri(5))=eos(iri(5))*moles*aprop/tmass
      endif
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
      write(*,'(" Z does not match Z in EOS5_data* files you are",
     . " using (66)")')
      write (*,'("mf,zz(mf)=",i5,e12.4)') mf,zz(mf)
      write (*,'("  iq,ip,k3,l3,xh,t6,r,z= ",4i5,4e12.4)')
     x ip,iq,k3,l3,xh,t6,r,z
      stop 1
      return
      end subroutine esac

!..*********************************************************************
      subroutine t6rinteos(slr,slt)
!     The purpose of this subroutine is to interpolate in T6 and rho
      save
      parameter (mx=5,mv=12,nr=169,nt=197)
      common/eeeos/ epl(mx,nt,nr),xx(mx)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/bbeos/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
      common/eeos/esact,eos(mv)
!
      iu=0
      is=0

      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quadeos(is,iw,slr,esk(kx,l1),esk(kx,l2),esk(kx,l3),
     x  rho(l1),rho(l2),rho(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quadeos(is,iw,slr,esk(kx,l2),esk(kx,l3),esk(kx,l4),
     x      rho(l2),rho(l3),rho(l4))
          endif
        is=1
      enddo
!
      is=0
      iw=1
!..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      esact=quadeos(is,iw,slt,h(1),h(2),h(3),t6a(k1),t6a(k2),t6a(k3))
        if (iq. eq. 3) then
!.....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          esactq=quadeos(is,iw,slt,q(1),q(2),q(3),t6a(k1),t6a(k2),
     x           t6a(k3))
        endif
        if(ip .eq. 3) then
!.....    eos(i) in lower-left 3x3.
          esact2=quadeos(is,iw,slt,h(2),h(3),h(4),t6a(k2),t6a(k3),
     x           t6a(k4))
!.....    eos(i) smoothed in left 3x4
          dix=(t6a(k3)-slt)*dfs(k3)
          esact=esact*dix+esact2*(1.-dix)
!       endif   ! moved to loc a
        if(iq .eq. 3) then
 
!.....     eos(i) in upper-right 3x3.
          esactq2=quadeos(is,iw,slt,q(2),q(3),q(4),t6a(k2),t6a(k3),
     x            t6a(k4))
          esactq=esactq*dix+esactq2*(1.-dix)
        endif
        endif  ! loc a
!
        if(iq .eq. 3) then
          dix2=(rho(l3)-slr)*dfsr(l3)
            if(ip .eq. 3) then
!.....        eos(i) smoothed in both log(T6) and log(R)
              esact=esact*dix2+esactq*(1-dix2)
            endif
        endif
        if (esact .gt. 1.e+15) then
          write(*,'("Interpolation indices out of range",
     x              ";please report conditions.")') 
          stop 1
        endif

      return
      end subroutine t6rinteos

!
!..********************************************************************
      subroutine readcoeos(filename)
!..... The purpose of this subroutine is to read the data tables
      save
      character (len=256) filename
      parameter (mx=5,mv=12,nr=169,nt=197)
      real moles
      character*1 blank
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/ iri(mv),index(mv),nta(nr),zz(mx)
      common/eeos/esact,eos(mv)
      common/eeeos/ epl(mx,nt,nr),xx(mx)
      common/eeeeos/moles(mx),xin(mx),tmass(mx),icycuse(mx,nr),
     x    rhogr(mx,nr),frac(mx,6),
     x    alogr(nr,nt)



      blank=' '


        if (itimeco .ne. 12345678) then
        do i=1,mx
          do j=1,mv 
            do k=1,nt
              do l=1,nr
                xz(i,j,k,l)=1.e+35
              enddo
            enddo
          enddo
        enddo
        itimeco=12345678
        endif
 
      close (2)
!..... read  tables
!       call system (' gunzip EOS5_data')
!       open(2, FILE='EOS5_data')
      
      write(*,*) 'read ' // trim(filename)
      
       open(2, FILE=trim(filename), IOSTAT=ios)
      if (ios .ne. 0) then
         write(*,*) 'failed to open ', trim(filename)
         stop 1
      end if
 
 
      do 3 m=1,mx
      read (2,'(3x,f6.4,3x,f12.9,11x,f10.7,17x,f10.7)')
     x  xin(m),zz(m),moles(m),tmass(m)
      read (2,'(21x,e14.7,4x,e14.7,3x,e11.4,3x,e11.4,3x,e11.4,
     x 4x,e11.4)') (frac(m,i),i=1,6)
      read (2,'(a)') blank
      do 2 jcs=1,nr
      read (2,'(2i5,2f12.7,17x,e15.7)') numtot,icycuse(m,jcs)
     x  ,dum,dum,rhogr(m,jcs)
      if(numtot .ne. jcs) then
         write (*,'(" Data file incorrect: numtot,jcs= ",2i5)')
     x     numtot,jcs
         stop 1
      endif
      read(2,'(a)') blank
      read(2,'(a)') blank
      if (icycuse(m,jcs) .lt. nta(jcs)) then
!         write (*,'("problem with data files: X=",f6.4," density=",
!     x      e14.4)') xin(m), rhogr(m,jcs)
         write (*,*) 'problem with ', filename
         stop 1
      endif
      do  i=1,icycuse(m,jcs)
      if (i .gt. nta(jcs)) then
         read (2,'(a)') blank
         go to 4
      endif
      ! change to read mu_M and logNe as part of xz table
      read (2,'(f11.6,1x,f6.4,e11.4,2e13.6,2e11.3,5f10.6)') 
     x t6list(jcs,i),(xz(m,index(iv),i,jcs),iv=10,11),
     x (xz(m,index(iv),i,jcs),iv=1,9)
    4 continue
      enddo
      read(2,'(a)') blank
      read(2,'(a)') blank
      read(2,'(a)') blank
    2 continue
      read(2,'(a)') blank
    3 continue
 
      do i=1,nt
         if(t6list(1,i) .eq. 0.0) then
            write(*,'("READCOEOS: Error:",i4,
     $           "-th T6 value is zero")') i
            stop 1
         endif
         t6a(i)=t6list(1,i)
      enddo
      do 12 i=2,nt
      dfs(i)=1./(t6a(i)-t6a(i-1))
   12 continue
      rho(1)=rhogr(1,1)
      do 13 i=2,nr
      rho(i)=rhogr(1,i)
      dfsr(i)=1./(rho(i)-rho(i-1))
   13 continue
      do i=2,mx
         dfsx(i)=1./(xx(i)-xx(i-1))
      enddo
!       call system ('gzip EOS5_data')
      
      return
      end subroutine readcoeos
!
!..*********************************************************************
      function quadeos(ic,i,x,y1,y2,y3,x1,x2,x3)
!..... this function performs a quadratic interpolation.
      save
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
      quadeos=c1+x*(c2+x*c3)
      return
      end function quadeos
!..****************************************************************8
      function gmass(x,z,amoles,eground,fracz,frac)
      dimension anum(6),frac(7), amas(7),eion(7)
      data (eion(i),i=1,6)/-3394.873554,-1974.86545,-1433.92718,
     x  -993.326315,-76.1959403,-15.29409/
      data (anum(i),i=1,6)/10.,8.,7.,6.,2.,1./
      xc=0.247137766
      xn=.0620782
      xo=.52837118
      xne=.1624188
      amas(7)=1.0079
      amas(6)=4.0026
      amas(5)=12.011
      amas(4)=14.0067
      amas(3)=15.9994
      amas(2)=20.179
      amas(1)=0.00054858
      fracz=z/(xc*amas(5)+xn*amas(4)+xo*amas(3)+xne*amas(2))
      xc2=fracz*xc
      xn2=fracz*xn
      xo2=fracz*xo
      xne2=fracz*xne 
      xh=x/amas(7)
      xhe=(1.-x -z)/amas(6)
      xtot=xh+xhe+xc2+xn2+xo2+xne2
      frac(6)=xh/xtot
      frac(5)=xhe/xtot
      frac(4)=xc2/xtot
      frac(3)=xn2/xtot
      frac(2)=xo2/xtot
      frac(1)=xne2/xtot
      eground=0.0
      amoles=0.0
      do i=1,6
      eground=eground+eion(i)*frac(i)
      amoles=amoles+(1.+anum(i))*frac(i)
      enddo
      anume=amoles-1.
      gmass=0.
      do i=2,7
      gmass=gmass+amas(i)*frac(i-1)
      enddo

      return
      end function gmass
!..*********************************************************************
      subroutine radsub (irad,t6,density,moles,tmass)
      parameter (mx=5,mv=12,nr=169,nt=197)
      real moles,k,molenak,Na
      common/eeos/esact,eos(mv)
      common/beos/ iri(mv),index(mv),nta(nr),zz(mx)

      data Na/6.0221367e+23/, k/1.380658e-16/, unitf/0.9648530/, 
     x unitfold/0.965296/, c/2.9979245e+10/, sigma/5.67051e-5/
     x , sigmac/1.8914785e-15/, sigmacc/1.8914785e-3/, aprop/83.14510/

!Physical constants
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
!       c=2.99792458E+10 !   cm/sec

!     rat=sigma/c    ! dyne/(cm**2*K**4)

!     rat=rat*1.e+24  !   Convert degrees K to units 10**6 K (T replaced with T6)
      rat=sigmacc

      molenak=moles*aprop  ! Mb-cc/unit T6

!..---Calculate EOS without radiation correction
      
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
!..-- End  no radiation calculation

      if (irad .ne. 0) then
!..-- Calculate EOS with radiation calculation
      pr=4./3.*rat*t6*t6*t6*t6   ! Mb 
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
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
!..---Add difference between EOS with and without radiation.  cvtt
!       calculation is not accurate enough to give accurate results using
!       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
      eos(iri(8))=eos(iri(8))+gam1t-gam1t_norad
      eos(iri(9))=eos(iri(9))+gam2pt-gam2pt_norad
!     eos(iri(mv))=eos(iri(mv))+gam3pt-gam3pt_norad
      endif
!..---End EOS calculations with radiation
!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
      eos(iri(1))=pt
      eos(iri(2))=et
      eos(iri(3))=st
      eos(iri(4))=dedrho
      eos(iri(5))=cvtt
      eos(iri(6))=chir
      eos(iri(7))=chitt
      return
      end subroutine radsub
   
!..*********************************************************************
      block data
      parameter (mx=5,mv=12,nr=169,nt=197)
      common/aaeos/ q(4),h(4),xxh
      common/aeos/  xz(mx,mv,nt,nr),  
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos/ iri(mv),index(mv),nta(nr),zz(mx)
      data (xa(i),i=1,mx)/0.0,0.2,0.4,0.6,0.8/

      data (index(i),i=1,mv)/1,2,3,4,5,6,7,8,9,10,11,12/
      data (nta(i),i=1,nr)/87*197,7*191,190,2*189,185,179,170,2*149,
     x  133,125,123,122,120,115,113,107,102,2*80,72,68,66,64,62,56,54,
     x  52,51,2*50,49,47,2*45,43,42,28*40,39,37,36,35,34,32,31,30,
     x  29,27,26/
      end

