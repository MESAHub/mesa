program main

   use interp_2d_lib_db, only: interp_mkbicub_db, interp_evbicub_db
   use interp_1d_lib, only: interp_pm, interp_value
   use interp_1d_def, only: pm_work_size
   use const_def
   use const_lib
   use colors_def
   use colors_lib
   use math_lib
   use utils_lib, only: mesa_error, mkdir, is_bad

   implicit none

   integer :: i,j,k,zone,denjmax,idum,ierr
   real*8, allocatable,dimension(:,:) :: &
      r,v,temp,den,kap,tempr,xm,smooth,tau,lum,n_bar,n_e
   real*8, allocatable,dimension(:) :: &
      t,m,dm,h,he,c,n,o,ne,na,mg,al,si,s,ar,ca,fe,ni
   real*8 :: dum,time,X,sum_tau,tauph,tau_extra,denmax,gdepos
   character*132 runname,filestr,fname,test_str
   character*256 line, my_mesa_dir

   real(dp), parameter :: &
      A_Fe56 = 56d0, lambda0 = 5169.02d-8, f = 0.023, Z_div_X_solar = 0.02293d0, &
      tau_sob_hi = 2d0, tau_sob_med = 1d0, tau_sob_lo = 0.2d0
   integer, parameter :: num_logRhos = 41, num_logTs = 117, iounit = 33, &
      max_lbol = 10000, n_colors = 5
   integer :: ilinx,iliny,ibcxmin,ibcxmax,ibcymin,ibcymax,ict(6), num_lbol, num_lbol_max
   real(dp) :: bcxmin(num_logTs), bcxmax(num_logTs), Ts(num_logTs)
   real(dp) :: bcymin(num_logRhos), bcymax(num_logRhos)
   real(dp) :: time_lbol(max_lbol), logL_lbol(max_lbol), t0, logL_lbol_max
   real(dp), pointer, dimension(:) :: logRhos, logTs, tau_sob_f1, tau_sob_values
   real(dp), pointer :: tau_sob_f(:,:,:)
   real(dp) :: tau_prev, tau_sob_prev, alfa, beta, L_div_Lsun, &
      v_sob_hi_tau, m_sob_hi_tau, r_sob_hi_tau, &
      v_sob_med_tau, m_sob_med_tau, r_sob_med_tau, &
      v_sob_lo_tau, m_sob_lo_tau, r_sob_lo_tau, &
      m_center, r_center, m_edge, r_edge, v_edge, &
      star_mass, mass_IB, skip, Lbol, logRho, logT, eta, &
      density, temperature, n_Fe, tau_sob, time_sec, fval(6), &
      rphot, mphot, log_g, Xsurf, Ysurf, Zsurf, Fe_H, Z_div_X, &
      bb_magU, bb_magB, bb_magV, bb_magR, bb_magI
   integer :: nm, num_models, cnt, k_phot, iday

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir,ierr)     
	if (ierr /= 0) then
	   write(*,*) 'const_init failed'
	   call mesa_error(__FILE__,__LINE__)
	end if        

   call math_init()

   call colors_init(1, &
      (/ trim(my_mesa_dir) // '/data/colors_data/blackbody_johnson.dat' /), &
      (/ n_colors /), ierr)
   if (ierr /= 0) then
      write(*,*) 'colors_init failed during initialization'
      return
   end if

   ! setup interpolation table for tau sob eta
   open(unit=iounit, file='FeII_5169_eta.dat', action='read')
   allocate(logRhos(num_logRhos), logTs(num_logTs), &
      tau_sob_f1(4*num_logRhos*num_logTs))
   tau_sob_f(1:4,1:num_logRhos,1:num_logTs) => &
      tau_sob_f1(1:4*num_logRhos*num_logTs)
   do j=1,num_logRhos
      do i=1,num_logTs
         read(iounit,*) density, temperature, eta
         logRho = log10(density)
         logRhos(j) = logRho
         logT = log10(temperature)
         if (j == 1) then
            Ts(i) = temperature
            logTs(i) = logT
         else if (logT /= logTs(i)) then
            write(*,*) 'bad T?', i, j, Ts(1), temperature, density, eta
            call mesa_error(__FILE__,__LINE__)
         end if
         tau_sob_f(1,j,i) = eta
      end do
   end do
   close(iounit)
   if (.false.) then
      do j=1,num_logTs
         if (is_bad(logTs(j))) write(*,*) 'logT', j, logTs(j)
      end do
      write(*,*)
      do j=1,num_logRhos
         if (is_bad(logRhos(j))) write(*,*) 'logRho', j, logRhos(j)
      end do
      write(*,*)
   end if
   ! just use "not a knot" bc's at edges of tables
   ibcxmin = 0; bcxmin(1:num_logTs) = 0
   ibcxmax = 0; bcxmax(1:num_logTs) = 0
   ibcymin = 0; bcymin(1:num_logRhos) = 0
   ibcymax = 0; bcymax(1:num_logRhos) = 0
   call interp_mkbicub_db( &
      logRhos, num_logRhos, logTs, num_logTs, tau_sob_f1, num_logRhos, &
      ibcxmin,bcxmin,ibcxmax,bcxmax, &
      ibcymin,bcymin,ibcymax,bcymax, &
      ilinx,iliny,ierr)
   if (ierr /= 0) then
      write(*,*) 'interp_mkbicub_db error'
      ierr = -1
      call mesa_error(__FILE__,__LINE__)
   end if

   do j=1,num_logRhos
      do i=1,num_logTs
         do k=1,4
            if (is_bad(tau_sob_f(k,j,i))) then
            write(*,*) 'tau_sob_f', i, j, k, tau_sob_f(k,j,i)
            end if
         end do
      end do
   end do

   ict = 0
   ict(1) = 1

   tauph=2d0/3d0
   tau_extra = 5d0

   if (iargc() > 0) then
      call GetArg(1,filestr)
   else
      filestr = 'mesa'
   end if
      
   fname = '../modmake/'//trim(filestr)//'.hyd'
   open(25,file=fname, status='old')
   write(*,*) 'read ' // trim(filestr)//'.hyd'
   read(25,*,end=333) dum,zone,m_center,r_center

   write(*,*) 'read ' // trim(filestr)//'.lbol'
   fname = trim(filestr)//'.lbol'
   open(21,file=fname, status='old')
   read(21,*) ! skip line
   num_lbol = 0
   num_lbol_max = 0
   time_lbol(1) = -1
   t0 = -99d0; logL_lbol_max = -99d0
   do i=1,10000
      read(21,*,end=334) time_lbol(num_lbol+1), skip, logL_lbol(num_lbol+1)
      if (time_lbol(num_lbol+1) < 0d0) cycle
      if (num_lbol > 1) then
         if (time_lbol(num_lbol+1) <= time_lbol(num_lbol)) cycle
      end if
      num_lbol = num_lbol+1
      if (logL_lbol(num_lbol) > logL_lbol_max) then
         logL_lbol_max = logL_lbol(num_lbol)
         t0 = time_lbol(num_lbol)
         num_lbol_max = num_lbol
      end if
   end do  

334 continue

   !write(*,*) 't0=runtime_logL_lbol_max, logL_lbol_max', t0, logL_lbol_max
   close(21)

   fname = trim(filestr)//'.swd'
   open(21,file=fname, status='old')
   num_models = 0
   do i=1,10000
      do j=1,zone
         read(21,*,end=111) time
      end do
      num_models = num_models+1
   end do

111    continue
   close(21)

   write(*,*) 'read ' // trim(filestr)//'.tt'
   fname = trim(filestr)//'.tt'
   open(21,file=fname, status='old')
  
   fname = trim(filestr)//'.lbol_lnuc.txt'
   open(24,file=fname, status='unknown')
  
   do i=1,1000
      read(21,'(a)',end=335) line
      if (len_trim(line) < 7) cycle
      if (line(1:6) /= '  time') cycle
      write(24,'(99a19)') 'days post max Lbol', 'log_Lbol (erg/s)', 'log_Lnuc (erg/s)'
      do j=1,10000
         read(21,*,end=335) time,&
            dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,gdepos
         if (time < t0-0.2d0) cycle
         !call interp_value(time_lbol, num_lbol, f1, time-t0, Lbol, ierr)
         !if (ierr /= 0) then
         !   write(*,*) 'failed in interp_value', time
         !   call mesa_error(__FILE__,__LINE__)
         !end if
         dum = interp_logLbol(time)
         write(24,'(99(1pe18.6,x))') time-t0, dum, log10(gdepos*1d50)
      end do
   end do  

335 continue
   close(21)
   close(24)

   fname = trim(filestr)//'.inner_boundary.txt'
   open(27,file=fname, status='unknown')

   fname = trim(filestr)//'.swd'
   open(21,file=fname, status='old')

   fname = '../modmake/'//trim(filestr)//'.abn'
   open(22,file=fname, status='old')

   fname = trim(filestr)//'.swd.ph'
   open(23,file=fname, status='unknown')

   fname = trim(filestr)//'.vel_feII'
   open(24,file=fname, status='unknown')

   nm = num_models
   allocate(r(zone,nm), xm(zone,nm), v(zone,nm), temp(zone,nm), den(zone,nm), &
      smooth(zone,nm), kap(zone,nm), tempr(zone,nm), tau(zone,nm), lum(zone,nm), &
      n_bar(zone,nm), n_e(zone,nm))

   allocate(t(nm),m(zone),dm(zone), &
      h(zone),he(zone),c(zone),n(zone),o(zone),ne(zone),na(zone),mg(zone), &
      al(zone),si(zone),s(zone),ar(zone),ca(zone),fe(zone),ni(zone))

   write(*,*) 'read ' // trim(filestr)//'.abn'
   do j=1,zone
      read(22,*,end=333) idum,m(j),dm(j),dum,h(j),he(j),c(j),n(j),o(j),&
      ne(j),na(j),mg(j),al(j),si(j),s(j),ar(j),ca(j),fe(j),dum,ni(j)
      if (idum /= j) then
         write(*,*) 'error in abn: idum /= j', idum, j
         call mesa_error(__FILE__,__LINE__)
      end if
   end do
   star_mass = m(zone)
   mass_IB = m(1)

   write(23,'(a)') '# photosphere.  bb_mag* are synthetic color magnitudes from mesa/colors.' &
     // '  see mesa.tt for stella color magnitudes from multiband rad-hydro.'
   write(23,'(a4,99(a18,x))') 'k', 't post max Lbol', 'radius(cm)', 'v(km/s)', 'mass(Msun)', &
      'logT', 'rho', 'kap', 'h', 'he', 'o', 'fe', 'tau_IB', 'Lbol', &
      'bb_magU', 'bb_magB', 'bb_magV', 'bb_magR', 'bb_magI', 'T_rad'

   write(24,'(a,3f8.3)') '# v where tau Sob for FeII 5169 is (low,medium,high) ', &
      tau_sob_lo, tau_sob_med, tau_sob_hi
   write(24,'(5x,99(a18,x))') &
      't post max Lbol', 'v_tau_lo', 'v_tau_med', 'v_tau_hi', 'rho', 'T', 'eta', &
      'm_tau_lo', 'r_tau_lo', 'm_tau_med', 'r_tau_med', 'm_tau_hi', 'r_tau_hi'
   write(27,'(99(a13))') 't-t0', 'tau', 'rho', 'T', 'r', 'v'

   write(*,*) 'read ' // trim(filestr)//'.swd'
   
   ! read data
   open(20,file=trim(filestr)//'.res', status='old')
   
   do nm=1,num_models

      do j=1,zone
         read(21,*) t(nm), idum, xm(j,nm),r(j,nm),v(j,nm),temp(j,nm),tempr(j,nm),&
            den(j,nm),dum,dum,dum,dum,kap(j,nm)
         if (idum /= j) then
            write(*,*) 'error in swd: idum /= j', idum, j
            call mesa_error(__FILE__,__LINE__)
         end if
         xm(j,nm)=10**xm(j,nm)
         r(j,nm)=10**r(j,nm)
         temp(j,nm)=10**temp(j,nm)
         if (tempr(j,nm) < 1d-6) then
            tempr(j,nm)=temp(j,nm)
         else
            tempr(j,nm)=10**tempr(j,nm)
         end if
         den(j,nm)=10**den(j,nm)*1d-6
      end do
      
      time = t(nm)
      if (time < 0.1d0) then
         iday = int(time*1d2 + 1d-6)
         write(test_str,'(a,i1)') '  OBS.TIME=        0.0', iday
      else if (time < 1d0) then
         iday = int(time*1d1 + 1d-6)
         write(test_str,'(a,i1)') '  OBS.TIME=        0.', iday
      else
         iday = int(time + 1d-6)
         if (iday < 10) then
            write(test_str,'(a,i1)') '  OBS.TIME=        ', iday
         else if (iday < 100) then
            write(test_str,'(a,i2)') '  OBS.TIME=       ', iday
         else
            write(test_str,'(a,i3)') '  OBS.TIME=      ', iday
         end if
      end if
      !write(*,*) 'zone', zone
      !write(*,*) 'test_str ' // trim(test_str)
      !write(*,*) 'len_trim(test_str)', len_trim(test_str)
      do i = 1, 500000
         read(20,'(a)') line
         if (line(1:len_trim(test_str)) /= trim(test_str)) cycle
         !write(*,*) 'i', i, line(1:len_trim(test_str))
         read(20,*) ! skip column headers
         do j=1,zone
            read(20,'(a)') line
            !write(*,*) i, j, trim(line)
            line(1:4) = '    '
            read(line,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,&
               lum(j,nm),dum,idum,n_bar(j,nm),n_e(j,nm)
         end do
         exit
      end do
      
   end do
   
   close(20)
   close(21)
   
   ! write results
   do nm=1,num_models
      time = t(nm)

      sum_tau=0.0d0
      tau(zone,nm) = sum_tau
      do j=zone,2,-1
         sum_tau=sum_tau+kap(j,nm)*den(j,nm)*(r(j,nm)-r(j-1,nm))
         tau(j-1,nm) = sum_tau
      end do
      
      write(27,'(99(1pd13.4))') time-t0, &
         tau(1,nm), den(1,nm), temp(1,nm), r(1,nm), v(1,nm)

      k_phot = 0
      if (sum_tau >= 1d0+tauph) then ! record photosphere information
         do j=zone,2,-1
            if(tau(j-1,nm) >= tauph) then
               k_phot = j
               beta = (tauph - tau(j,nm))/(tau(j-1,nm) - tau(j,nm))
               alfa = 1d0 - beta
               !call interp_value(time_lbol, num_lbol, f1, time, Lbol, ierr)
               !if (ierr /= 0) then
               !   write(*,*) 'failed in interp_value', time
               !   call mesa_error(__FILE__,__LINE__)
               !end if
               dum = interp_logLbol(time)
               Lbol = exp10(dum) 
               L_div_Lsun = Lbol/Lsun
               rphot = alfa*r(j,nm)+beta*r(j-1,nm)
               mphot = star_mass-(alfa*xm(j,nm)+beta*xm(j-1,nm))
               log_g = log10(standard_cgrav*mphot*Msun/(rphot*rphot))
               logT = log10(alfa*temp(j,nm)+beta*temp(j-1,nm))
               Xsurf = alfa*h(j)+beta*h(j-1)
               Ysurf = alfa*he(j)+beta*he(j-1)
               Zsurf = max(1d-99, min(1d0, 1d0 - (Xsurf + Ysurf)))
               Z_div_X = Zsurf/max(1d-99, Xsurf)
               Fe_H = log10(Z_div_X/Z_div_X_solar)
               bb_magU = get1_synthetic_color_abs_mag('bb_U')
               bb_magB = get1_synthetic_color_abs_mag('bb_B')
               bb_magV = get1_synthetic_color_abs_mag('bb_V')
               bb_magR = get1_synthetic_color_abs_mag('bb_R')
               bb_magI = get1_synthetic_color_abs_mag('bb_I')
               write(23,'(i5,99(1pe18.6,x))') j, time-t0, &
                  rphot, (alfa*v(j,nm)+beta*v(j-1,nm))*1d3, mphot, logT, &
                  alfa*den(j,nm)+beta*den(j-1,nm), alfa*kap(j,nm)+beta*kap(j-1,nm), &
                  Xsurf, Ysurf, alfa*o(j)+beta*o(j-1), alfa*fe(j)+beta*fe(j-1), &
                  sum_tau, Lbol, bb_magU, bb_magB, bb_magV, bb_magR, bb_magI, &
                  alfa*tempr(j,nm)+beta*tempr(j-1,nm)
               exit
            end if
         end do
      end if

      do k=1,zone
         smooth(k,nm) = den(k,nm)
      end do
      do j=1,5
         do k=2,zone-1
            smooth(k,nm) = sum(smooth(k-1:k+1,nm))/3d0
         end do
      end do
      tau_sob_prev = 0d0
      v_sob_lo_tau = 0d0
      v_sob_med_tau = 0d0
      v_sob_hi_tau = 0d0
      do k=zone-1,2,-1
         density = smooth(k,nm)
         n_Fe = density*avo*fe(k)/A_Fe56
         if (is_bad(n_Fe)) then
            write(*,*) 'n_Fe', k, n_Fe, density, avo, fe(k), A_Fe56
            call mesa_error(__FILE__,__LINE__)
         end if
         logRho = log10(density)
         logRho = min(logRhos(num_logRhos), max(logRhos(1), logRho))  
         logT = log10(temp(k,nm))
         logT = min(logTs(num_logTs), max(logTs(1), logT))
         ierr = 0
         call interp_evbicub_db( &
            logRho, logT, logRhos, num_logRhos, logTs, num_logTs, &
            ilinx, iliny, tau_sob_f1, num_logRhos, ict, fval, ierr)
         if (ierr /= 0) then
            write(*,*) 'logRho', k, logRho
            write(*,*) 'logT', k, logT
            write(*,*) 'interp failed in data_for_extra_profile_columns'
            call mesa_error(__FILE__,__LINE__)
         end if
         eta = fval(1)
         if (is_bad(eta)) then
            write(*,*) 'eta', k, eta, logT, logRho, n_Fe
            call mesa_error(__FILE__,__LINE__)
         end if
         time_sec = time*60*60*24
         tau_sob = pi*qe*qe/(me*clight)*n_Fe*eta*f*time_sec*lambda0
         
         if ((tau_sob > tau_sob_lo .or. k==2) .and. v_sob_lo_tau == 0d0) then
            alfa = (tau_sob_lo - tau_sob_prev)/(tau_sob - tau_sob_prev)
            beta = 1d0 - alfa
            v_sob_lo_tau = alfa*v(k,nm) + beta*v(k+1,nm)
            m_sob_lo_tau = star_mass - (alfa*xm(k,nm) + beta*xm(k+1,nm))
            r_sob_lo_tau = alfa*r(k,nm) + beta*r(k+1,nm)
         end if
         
         if ((tau_sob > tau_sob_med .or. k==2) .and. v_sob_med_tau == 0d0) then
            alfa = (tau_sob_med - tau_sob_prev)/(tau_sob - tau_sob_prev)
            beta = 1d0 - alfa
            v_sob_med_tau = alfa*v(k,nm) + beta*v(k+1,nm)
            m_sob_med_tau = star_mass - (alfa*xm(k,nm) + beta*xm(k+1,nm))
            r_sob_med_tau = alfa*r(k,nm) + beta*r(k+1,nm)
         end if
         
         if (tau_sob > tau_sob_hi .or. k==2) then
            alfa = (tau_sob_hi - tau_sob_prev)/(tau_sob - tau_sob_prev)
            beta = 1d0 - alfa
            v_sob_hi_tau = alfa*v(k,nm) + beta*v(k+1,nm)
            m_sob_hi_tau = star_mass - (alfa*xm(k,nm) + beta*xm(k+1,nm))
            r_sob_hi_tau = alfa*r(k,nm) + beta*r(k+1,nm)
            if (k > k_phot - 10) &
               write(24,'(i5,99(1pe18.4,x))') k,time-t0, &
                  v_sob_lo_tau*1d3, v_sob_med_tau*1d3, v_sob_hi_tau*1d3, &
                  10**logRho, 10**logT, eta, &
                  m_sob_lo_tau, r_sob_lo_tau, &
                  m_sob_med_tau, r_sob_med_tau, &
                  m_sob_hi_tau, r_sob_hi_tau
            exit
         end if 
               
         tau_sob_prev = tau_sob                
      end do

   end do

333 continue

   call save_day_post_Lbol_max(0.1d0,t0,zone,star_mass,mass_IB,'000.1')
   call save_day_post_Lbol_max(0.2d0,t0,zone,star_mass,mass_IB,'000.2')
   call save_day_post_Lbol_max(0.3d0,t0,zone,star_mass,mass_IB,'000.3')
   call save_day_post_Lbol_max(0.4d0,t0,zone,star_mass,mass_IB,'000.4')
   call save_day_post_Lbol_max(0.5d0,t0,zone,star_mass,mass_IB,'000.5')
   call save_day_post_Lbol_max(0.6d0,t0,zone,star_mass,mass_IB,'000.6')
   call save_day_post_Lbol_max(0.7d0,t0,zone,star_mass,mass_IB,'000.7')
   call save_day_post_Lbol_max(0.8d0,t0,zone,star_mass,mass_IB,'000.8')
   call save_day_post_Lbol_max(0.9d0,t0,zone,star_mass,mass_IB,'000.9')
   call save_day_post_Lbol_max(1d0,t0,zone,star_mass,mass_IB,'001')
   call save_day_post_Lbol_max(2d0,t0,zone,star_mass,mass_IB,'002')
   call save_day_post_Lbol_max(3d0,t0,zone,star_mass,mass_IB,'003')
   call save_day_post_Lbol_max(4d0,t0,zone,star_mass,mass_IB,'004')
   call save_day_post_Lbol_max(5d0,t0,zone,star_mass,mass_IB,'005')
   
   call save_day_post_Lbol_max(10d0,t0,zone,star_mass,mass_IB,'010')
   call save_day_post_Lbol_max(20d0,t0,zone,star_mass,mass_IB,'020')
   call save_day_post_Lbol_max(30d0,t0,zone,star_mass,mass_IB,'030')
   call save_day_post_Lbol_max(40d0,t0,zone,star_mass,mass_IB,'040')
   call save_day_post_Lbol_max(50d0,t0,zone,star_mass,mass_IB,'050')
   call save_day_post_Lbol_max(60d0,t0,zone,star_mass,mass_IB,'060')
   call save_day_post_Lbol_max(70d0,t0,zone,star_mass,mass_IB,'070')
   call save_day_post_Lbol_max(80d0,t0,zone,star_mass,mass_IB,'080')
   call save_day_post_Lbol_max(90d0,t0,zone,star_mass,mass_IB,'090')
   call save_day_post_Lbol_max(100d0,t0,zone,star_mass,mass_IB,'100')
   call save_day_post_Lbol_max(110d0,t0,zone,star_mass,mass_IB,'110')
   call save_day_post_Lbol_max(120d0,t0,zone,star_mass,mass_IB,'120')
   call save_day_post_Lbol_max(130d0,t0,zone,star_mass,mass_IB,'130')
   call save_day_post_Lbol_max(140d0,t0,zone,star_mass,mass_IB,'140')
   call save_day_post_Lbol_max(150d0,t0,zone,star_mass,mass_IB,'150')
   call save_day_post_Lbol_max(160d0,t0,zone,star_mass,mass_IB,'160')

   write(*,*) 'write ' // trim(filestr)//'.swd.ph'
   write(*,*) 'write ' // trim(filestr)//'.vel_feII'
   write(*,*) 'write ' // trim(filestr)//'.lbol_lnuc.txt'
   write(*,*) 'write ' // trim(filestr)//'.inner_boundary.txt'
   write(*,*) 'done'

   contains

   real(dp) function interp_logLbol(time)
      real(dp), intent(in) :: time ! time since start of run
      integer :: k
      real(dp) :: alfa, beta
      do k=1,num_lbol-1
         if (time_lbol(k) <= time .and. time < time_lbol(k+1)) then
            alfa = (time_lbol(k+1) - time)/(time_lbol(k+1) - time_lbol(k))
            beta = 1d0 - alfa
            interp_logLbol = alfa*logL_lbol(k) + beta*logL_lbol(k+1)
            return
         end if
      end do
      interp_logLbol = logL_lbol(num_lbol)
   end function interp_logLbol
   
   real(dp) function get1_synthetic_color_abs_mag(name) result(mag)
      character (len=*) :: name
      mag = get_abs_mag_by_name(name, logT, log_g, Fe_H, L_div_Lsun, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in get_abs_mag_by_id ' // trim(name), &
            time, logT, log_g, Fe_H, Lbol
         call mesa_error(__FILE__,__LINE__)
      end if
   end function get1_synthetic_color_abs_mag

   subroutine save_day_post_Lbol_max(day, t0, zone, star_mass, mass_IB, daystr)
      real(dp), intent(in) :: day, t0, star_mass, mass_IB
      integer, intent(in) :: zone
      character (len=*), intent(in) :: daystr
      
      integer, parameter :: io = 26, ncol = 36
      real(dp) :: tnm, t1, t2, alfa, beta, data1(ncol), data2(ncol)
      integer :: k, i, nm1, nm2
      character (len=132) :: fname
      include 'formats'
      
      if (t0 < 0d0) return
      
      tnm = day + t0 ! this is the desired run time
      nm1 = 0
      nm2 = 0
      do nm=1,num_models-1
         if (t(nm) <= tnm .and. tnm <= t(nm+1)) then
            nm1 = nm; nm2 = nm+1
            t1 = t(nm1); t2 = t(nm2)
            exit
         end if
      end do
      
      if (nm1 == 0 .or. nm2 == 0) then
         write(*,*) 'save_day_post_Lbol_max failed to find models for', day
         return
      end if

      alfa = (t2 - tnm)/(t2 - t1) ! fraction from 1st time
      beta = 1d0 - alfa
      
      if (alfa > 1d0 .or. alfa < 0d0) then
         write(*,*) 'save_day_post_Lbol_max failed for', day
         return
      end if
      
      write(fname,'(a)') trim(filestr)// '.day' // trim(daystr) // '_post_Lbol_max.data'
      open(io,file=trim(fname), status='unknown')

      write(io,'(a20,f25.1)') 'days post max Lbol', day
      write(io,'(a20,i25)') 'zones', zone
      write(io,'(a20,2(1p,e25.15))') 'inner boundary mass', mass_IB*msun, mass_IB
      write(io,'(a20,2(1p,e25.15))') 'total mass', star_mass*msun, star_mass
      write(io,*)
      write(io,'(99a25)') '', &
         'mass of cell (g)', & 
         'cell center m (g)', & 
         'cell center R (cm)', & 
         'cell center v (cm/s)', &
         'avg density', & 
         'radiation pressure', & 
         'avg temperature', & 
         'radiation temperature', & 
         'avg opacity', & 
         'tau', & 
         'outer edge m (g)', & 
         'outer edge r (cm)', & 
         'h1', & 
         'he3', & 
         'he4', & 
         'c12', & 
         'n14', & 
         'o16', & 
         'ne20', & 
         'na23', & 
         'mg24', & 
         'si28', & 
         's32', & 
         'ar36', & 
         'ca40', & 
         'ti44', & 
         'cr48', & 
         'cr60', & 
         'fe52', & 
         'fe54', & 
         'fe56', & 
         'co56', & 
         'ni56', &
         'luminosity', &
         'n_bar', &
         'n_e' 
      write(io,*)
      do j=1,zone  
         !read(io1,*) i1, data1(1:ncol)
         !read(io2,*) i2, data2(1:ncol)
         call get1_data(j,nm1,data1)
         call get1_data(j,nm2,data2)
         !if (i1 /= j .or. i2 /= j) then
         !   write(*,*) 'bad zone', i1, i2, j
         !   call mesa_error(__FILE__,__LINE__,'save_day50_post_Lbol_max')
         !end if
         write(io,'(i25)', advance = 'no') j
         do k=1,ncol
            !write(*,*) k, alfa*data1(k) + beta*data2(k), data1(k), data2(k)
            write(io,'(1p,e25.15)', advance = 'no') alfa*data1(k) + beta*data2(k)
         end do
         write(io,*)
      end do

      close(io)
      !close(io1)
      !close(io2)

      write(*,*) 'write ' // trim(fname)
      return
333    continue
      write(*,*) 'failed in save_day_post_Lbol_max'
      call mesa_error(__FILE__,__LINE__)
      
   end subroutine save_day_post_Lbol_max
      
   subroutine get1_data(j,nm,d)
      integer, parameter :: ncol = 36
      integer, intent(in) :: j,nm
      real(dp), intent(out) :: d(ncol)
      real(dp) :: m_edge, r_edge, v_edge
      if (j < zone) then
         m_edge = 0.5d0*(m(j) + m(j+1))*msun
         r_edge = 0.5d0*(r(j,nm) + r(j+1,nm))
         v_edge = 0.5d0*(v(j,nm) + v(j+1,nm))*1d8
      else
         m_edge = m(j)*msun
         r_edge = r(j,nm)
         v_edge = v(j,nm)*1d8
      end if
      d(1:ncol) = (/ &
         dm(j)*msun, m(j)*msun, r(j,nm), v_edge, den(j,nm), &
         crad*tempr(j,nm)**4/3, temp(j,nm), tempr(j,nm), kap(j,nm), tau(j,nm), &
         m_edge, r_edge, &
         h(j),0d0,he(j),c(j),n(j),o(j),ne(j),na(j),mg(j),si(j),s(j),ar(j), ca(j), &
         0d0, 0d0, 0d0, 0d0, 0d0, fe(j), 0d0, ni(j), lum(j,nm), n_bar(j,nm), n_e(j,nm) &
         /)
   end subroutine get1_data


end program main
