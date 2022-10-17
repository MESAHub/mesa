subroutine read_inlist_stella(ierr)
   integer, intent(out) :: ierr
   character (len = 256) :: filename, message
   integer :: unit
   filename = 'inlist_stella'
   
   write(*, *) 'read_inlist_stella'
   
   ! set defaults
   save_stella_data_when_terminate = .true.
   save_stella_data_filename = 'mesa'
   stella_num_points = 400
   stella_nz_extra = 0
   stella_min_surf_logRho = -100
   stella_min_velocity = 0
   stella_skip_inner_dm = 0
   stella_skip_inner_v_limit = 5d+07
   stella_mdot_years_for_wind = 1
   stella_mdot_for_wind = 0.25
   stella_v_wind = 10
   stella_show_headers = .false.
   
   open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', iostat = ierr)
   if (ierr /= 0) then
      write(*, *) 'Failed to open control namelist file ', trim(filename)
   else
      read(unit, nml = stella, iostat = ierr)
      close(unit)
      if (ierr /= 0) then
         write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
         write(*, '(a)') &
               'The following runtime error message might help you find the problem'
         write(*, *)
         open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
         read(unit, nml = stella)
         close(unit)
      end if
   end if

end subroutine read_inlist_stella


subroutine write_stella_data(s, ierr)
   type (star_info), pointer :: s
   integer, intent(out) :: ierr
   call write_stella_files(s, &
         save_stella_data_filename, &
         stella_min_surf_logRho, &
         stella_min_velocity, stella_skip_inner_dm, &
         stella_skip_inner_v_limit, &
         stella_num_points, stella_nz_extra, &
         stella_mdot_years_for_wind, stella_mdot_for_wind, &
         stella_v_wind, stella_show_headers, ierr)
end subroutine write_stella_data


subroutine write_stella_files(s, &
      fname, &
      min_surf_logRho, &
      min_velocity, &
      skip_inner_dm, &
      skip_inner_v_limit, &
      num_points, &
      nz_extra, &
      mdot_years_for_CSM, &
      mdot_for_CSM, &
      v_CSM, &
      show_headers, ierr)
   use chem_def
   use interp_1d_def
   use interp_1d_lib
   type (star_info), pointer :: s
   character (len = *), intent(in) :: fname
   real(dp), intent(in) :: &
         min_surf_logRho, &
         min_velocity, &
         skip_inner_dm, &
         skip_inner_v_limit, &
         mdot_years_for_CSM, &
         mdot_for_CSM, &
         v_CSM
   integer, intent(in) :: &
         num_points, &
         nz_extra
   logical, intent(in) :: show_headers
   integer, intent(out) :: ierr
   
   integer, parameter :: io_hyd = 34, io_abn = 35
   integer :: &
         i_h1, i_he3, i_he4, i_c12, i_n14, i_o16, i_ne20, &
         i_mg24, i_si28, i_s32, i_ar36, i_ca40, i_ti44, &
         i_cr48, i_cr60, i_fe52, i_fe54, i_fe56, i_co56, i_ni56, &
         nz, n, k_center, k, k_surf, nz_mod, kout, j, i
   real(dp) :: &
         radius_wind, di, m_center, r_center, rho_center, K_wind, &
         total_co56, total_ni56, f_na23, f_al27, f, &
         dr, r, mass, density, temperature, pressure, minx, m_prev, dm, &
         velocity, energy, opacity, entropy, tau, v_div_csound, &
         h1, he3, he4, c12, n14, o16, ne20, &
         mg24, si28, s32, ar36, ca40, ti44, &
         cr48, cr60, fe52, fe54, fe56, co56, ni56, &
         na23, al27, xsum, r_inner, vol, tot_ni56, tot_co56
   real(dp), pointer, dimension(:) :: &
         model_ks, interp_ks, mmid, rmid, taumid, work, &
         mmid_mod, rmid_mod, tau_mod, u_mod, rho_mod, &
         pressure_mod, temperature_mod, energy_mod, &
         entropy_mod, opacity_mod, v_div_csound_mod, &
         h1_mod, he3_mod, he4_mod, c12_mod, n14_mod, &
         o16_mod, ne20_mod, mg24_mod, si28_mod, s32_mod, &
         ar36_mod, ca40_mod, ti44_mod, cr48_mod, cr60_mod, &
         fe52_mod, fe54_mod, fe56_mod, co56_mod, ni56_mod, &
         mmid_out, rmid_out, tau_out, u_out, rho_out, &
         pressure_out, temperature_out, energy_out, &
         entropy_out, opacity_out, v_div_csound_out, &
         h1_out, he3_out, he4_out, c12_out, n14_out, &
         o16_out, ne20_out, mg24_out, si28_out, s32_out, &
         ar36_out, ca40_out, ti44_out, cr48_out, cr60_out, &
         fe52_out, fe54_out, fe56_out, co56_out, ni56_out
   
   include 'formats'
   
   ierr = 0
   
   i_h1 = s% net_iso(ih1)
   i_he3 = s% net_iso(ihe3)
   i_he4 = s% net_iso(ihe4)
   i_c12 = s% net_iso(ic12)
   i_n14 = s% net_iso(in14)
   i_o16 = s% net_iso(io16)
   i_ne20 = s% net_iso(ine20)
   i_mg24 = s% net_iso(img24)
   i_si28 = s% net_iso(isi28)
   i_s32 = s% net_iso(is32)
   i_ar36 = s% net_iso(iar36)
   i_ca40 = s% net_iso(ica40)
   i_ti44 = s% net_iso(iti44)
   i_cr48 = s% net_iso(icr48)
   i_cr60 = s% net_iso(icr60)
   i_fe52 = s% net_iso(ife52)
   i_fe54 = s% net_iso(ife54)
   i_fe56 = s% net_iso(ife56)
   i_co56 = s% net_iso(ico56)
   i_ni56 = s% net_iso(ini56)
   
   open(unit = io_hyd, file = trim(fname) // '.hyd', action = 'write', status = 'replace')
   open(unit = io_abn, file = trim(fname) // '.abn', action = 'write', status = 'replace')
   
   nz = s% nz
   n = num_points - nz_extra
   if (nz_extra > 0) then
      radius_wind = s% r(1) + v_CSM * 1d5 * mdot_years_for_CSM * secyer
   else
      radius_wind = 0
   end if
   
   k_center = nz
   do k = nz, 1, -1
      if (s% m(k) >= s% M_center + skip_inner_dm * Msun) then
         k_center = k
         exit
      end if
   end do
   do k = k_center, 1, -1
      if (s% u(k) > skip_inner_v_limit) then
         k_center = k
         exit
      end if
   end do
   k_surf = 1
   do k = 1, k_center
      if (s% lnd(k) >= min_surf_logRho * ln10) then
         k_surf = k
         exit
      end if
   end do
   
   tot_ni56 = dot_product(s% dm(1:k_center), s% xa(i_ni56, 1:k_center))
   tot_co56 = dot_product(s% dm(1:k_center), s% xa(i_co56, 1:k_center))
   !write(*,2) 'k_center', k_center, s% m(k_center)/Msun, (tot_co56 + tot_ni56)/Msun
   
   nz_mod = k_center - k_surf + 1
   
   allocate(&
         mmid(nz), rmid(nz), taumid(nz), work(nz * pm_work_size))
   k = nz_mod
   allocate(&
         model_ks(k), &
         mmid_mod(k), rmid_mod(k), tau_mod(k), u_mod(k), rho_mod(k), &
         pressure_mod(k), temperature_mod(k), energy_mod(k), &
         entropy_mod(k), opacity_mod(k), v_div_csound_mod(k), &
         h1_mod(k), he3_mod(k), he4_mod(k), c12_mod(k), n14_mod(k), &
         o16_mod(k), ne20_mod(k), mg24_mod(k), si28_mod(k), s32_mod(k), &
         ar36_mod(k), ca40_mod(k), ti44_mod(k), cr48_mod(k), cr60_mod(k), &
         fe52_mod(k), fe54_mod(k), fe56_mod(k), co56_mod(k), ni56_mod(k))
   allocate(&
         interp_ks(n), &
         mmid_out(n), rmid_out(n), tau_out(n), u_out(n), rho_out(n), &
         pressure_out(n), temperature_out(n), energy_out(n), &
         entropy_out(n), opacity_out(n), v_div_csound_out(n), &
         h1_out(n), he3_out(n), he4_out(n), c12_out(n), n14_out(n), &
         o16_out(n), ne20_out(n), mg24_out(n), si28_out(n), s32_out(n), &
         ar36_out(n), ca40_out(n), ti44_out(n), cr48_out(n), cr60_out(n), &
         fe52_out(n), fe54_out(n), fe56_out(n), co56_out(n), ni56_out(n))
   
   do k = 1, nz_mod
      model_ks(k) = k
   end do
   kout = 10
   di = (nz_mod - kout) * 1d0 / (n - kout)
   do k = 1, n - kout
      interp_ks(k) = 1 + (k - 1) * di
   end do
   do k = n - kout + 1, n
      interp_ks(k) = interp_ks(k - 1) + 1d0
   end do
   
   do k = 1, nz - 1
      mmid(k) = 0.5d0 * (s% m(k) + s% m(k + 1))
      rmid(k) = 0.5d0 * (s% r(k) + s% r(k + 1))
      taumid(k) = 0.5d0 * (s% tau(k) + s% tau(k + 1))
   end do
   mmid(nz) = 0.5d0 * (s% m(nz) + s% M_center)
   rmid(nz) = 0.5d0 * (s% r(nz) + s% R_center)
   taumid(nz) = s% tau(nz)
   
   ! vectors of model data at cell centers ordered from inner to outer
   do k = 1, nz_mod
      j = k_center - k + 1
      mmid_mod(k) = mmid(j)
      rmid_mod(k) = rmid(j)
      tau_mod(k) = taumid(j)
      u_mod(k) = s% u(j)
      rho_mod(k) = s% rho(j)
      pressure_mod(k) = s% Peos(j)
      temperature_mod(k) = s% T(j)
      energy_mod(k) = s% energy(j)
      entropy_mod(k) = s% entropy(j)
      opacity_mod(k) = s% opacity(j)
      v_div_csound_mod(k) = s% u(j) / s% csound(j)
      if (i_h1 > 0) h1_mod(k) = s% xa(i_h1, j)
      if (i_he3 > 0) he3_mod(k) = s% xa(i_he3, j)
      if (i_he4 > 0) he4_mod(k) = s% xa(i_he4, j)
      if (i_c12 > 0) c12_mod(k) = s% xa(i_c12, j)
      if (i_n14 > 0) n14_mod(k) = s% xa(i_n14, j)
      if (i_o16 > 0) o16_mod(k) = s% xa(i_o16, j)
      if (i_ne20 > 0) ne20_mod(k) = s% xa(i_ne20, j)
      if (i_mg24 > 0) mg24_mod(k) = s% xa(i_mg24, j)
      if (i_si28 > 0) si28_mod(k) = s% xa(i_si28, j)
      if (i_s32 > 0) s32_mod(k) = s% xa(i_s32, j)
      if (i_ar36 > 0) ar36_mod(k) = s% xa(i_ar36, j)
      if (i_ca40 > 0) ca40_mod(k) = s% xa(i_ca40, j)
      if (i_ti44 > 0) ti44_mod(k) = s% xa(i_ti44, j)
      if (i_cr48 > 0) cr48_mod(k) = s% xa(i_cr48, j)
      if (i_cr60 > 0) cr60_mod(k) = s% xa(i_cr60, j)
      if (i_fe52 > 0) fe52_mod(k) = s% xa(i_fe52, j)
      if (i_fe54 > 0) fe54_mod(k) = s% xa(i_fe54, j)
      if (i_fe56 > 0) fe56_mod(k) = s% xa(i_fe56, j)
      if (i_co56 > 0) co56_mod(k) = s% xa(i_co56, j)
      if (i_ni56 > 0) ni56_mod(k) = s% xa(i_ni56, j)
   end do
   
   call interp1(mmid_mod, mmid_out)
   call interp1(rmid_mod, rmid_out)
   call interp1(u_mod, u_out)
   call interp1(rho_mod, rho_out)
   call interp1(pressure_mod, pressure_out)
   call interp1(temperature_mod, temperature_out)
   call interp1(energy_mod, energy_out)
   call interp1(entropy_mod, entropy_out)
   call interp1(opacity_mod, opacity_out)
   call interp1(tau_mod, tau_out)
   call interp1(v_div_csound_mod, v_div_csound_out)
   if (i_h1 > 0) call interp1(h1_mod, h1_out)
   if (i_he3 > 0) call interp1(he3_mod, he3_out)
   if (i_he4 > 0) call interp1(he4_mod, he4_out)
   if (i_c12 > 0) call interp1(c12_mod, c12_out)
   if (i_n14 > 0) call interp1(n14_mod, n14_out)
   if (i_o16 > 0) call interp1(o16_mod, o16_out)
   if (i_ne20 > 0) call interp1(ne20_mod, ne20_out)
   if (i_mg24 > 0) call interp1(mg24_mod, mg24_out)
   if (i_si28 > 0) call interp1(si28_mod, si28_out)
   if (i_s32 > 0) call interp1(s32_mod, s32_out)
   if (i_ar36 > 0) call interp1(ar36_mod, ar36_out)
   if (i_ca40 > 0) call interp1(ca40_mod, ca40_out)
   if (i_ti44 > 0) call interp1(ti44_mod, ti44_out)
   if (i_cr48 > 0) call interp1(cr48_mod, cr48_out)
   if (i_cr60 > 0) call interp1(cr60_mod, cr60_out)
   if (i_fe52 > 0) call interp1(fe52_mod, fe52_out)
   if (i_fe54 > 0) call interp1(fe54_mod, fe54_out)
   if (i_fe56 > 0) call interp1(fe56_mod, fe56_out)
   if (i_co56 > 0) call interp1(co56_mod, co56_out)
   if (i_ni56 > 0) call interp1(ni56_mod, ni56_out)
   
   if (k_center == nz) then
      m_center = s% M_center
      r_center = s% R_center
   else
      m_center = s% m(k_center + 1)
      r_center = s% r(k_center + 1)
   end if
   rho_center = s% rho(k_center)
   
   if (show_headers) &
         write(io_hyd, '(12x,a)') '   n      mass IB        R IB           rho IB'
   write(io_hyd, '(1pe12.6,i5,3(3x,1pe12.6))') &
         0.0d0, n + nz_extra, m_center / Msun, r_center, rho_center
   
   if (show_headers) then
      write(io_hyd, '(99a24)') &
            "mass (Msun)", & !      # 1
            "radius (cm)", & !      # 2
            "density", & !          # 3
            "temperature", & !      # 4
            "velocity (cm/s)", & !  # 5
            "pressure", & !         # 6
            " ", & !                # 7
            "energy", & !           # 8
            "entropy", & !          # 9
            "opacity", & !          # 10
            "tau", & !              # 11
            "v_div_csound"   !      # 12
      write(io_abn, '(3a13,99a14)') &
            "mass", &
            "dm", &
            "r", &
            "h", &
            "he", &
            "c", &
            "n", &
            "o", &
            "ne", &
            "na", &
            "mg", &
            "al", &
            "si", &
            "s", &
            "ar", &
            "ca", &
            "fe", &
            "", &
            "co+ni"
   end if
   
   f_na23 = 0.044538489161177d0 ! Lodders 2009 ratio of na/n.
   f_al27 = 0.0766594958504439d0 ! Lodders 2009 ratio of al/n.
   
   total_co56 = 0
   total_ni56 = 0
   h1 = 0
   he3 = 0
   he4 = 0
   c12 = 0
   n14 = 0
   o16 = 0
   ne20 = 0
   mg24 = 0
   si28 = 0
   s32 = 0
   ar36 = 0
   ca40 = 0
   ti44 = 0
   cr48 = 0
   cr60 = 0
   fe52 = 0
   fe54 = 0
   fe56 = 0
   co56 = 0
   ni56 = 0
   na23 = 0
   al27 = 0
   minx = 1d-30
   
   ! measure ni56 + co56
   mass = m_center
   do i = 1, n
      m_prev = mass
      mass = mmid_out(i)
      dm = mass - m_prev
      
      if (i_h1 > 0) h1 = h1_out(i)
      if (i_he3 > 0) he3 = he3_out(i)
      if (i_he4 > 0) he4 = he4_out(i)
      if (i_c12 > 0) c12 = c12_out(i)
      if (i_n14 > 0) n14 = n14_out(i)
      if (i_o16 > 0) o16 = o16_out(i)
      if (i_ne20 > 0) ne20 = ne20_out(i)
      if (i_mg24 > 0) mg24 = mg24_out(i)
      if (i_si28 > 0) si28 = si28_out(i)
      if (i_s32 > 0) s32 = s32_out(i)
      if (i_ar36 > 0) ar36 = ar36_out(i)
      if (i_ca40 > 0) ca40 = ca40_out(i)
      if (i_ti44 > 0) ti44 = ti44_out(i)
      if (i_cr48 > 0) cr48 = cr48_out(i)
      if (i_cr60 > 0) cr60 = cr60_out(i)
      if (i_fe52 > 0) fe52 = fe52_out(i)
      if (i_fe54 > 0) fe54 = fe54_out(i)
      if (i_co56 > 0) co56 = co56_out(i)
      if (i_ni56 > 0) ni56 = ni56_out(i)
      na23 = f_na23 * n14
      al27 = f_al27 * n14
      n14 = n14 - na23 - al27
      
      xsum = h1 + he3 + he4 + c12 + n14 + o16 + ne20 + na23 + mg24 + al27 + si28 + &
            s32 + ar36 + ca40 + ti44 + cr48 + cr60 + fe52 + fe54 + fe56 + co56 + ni56
      if (xsum > 1.1d0 .or. xsum < 0.9d0 .or. is_bad(xsum)) then
         write(*, 2) 'bad raw xsum', i, xsum
         stop 'write_stella_files'
      end if
      h1 = h1 / xsum
      he3 = he3 / xsum
      he4 = he4 / xsum
      c12 = c12 / xsum
      n14 = n14 / xsum
      o16 = o16 / xsum
      ne20 = ne20 / xsum
      na23 = na23 / xsum
      mg24 = mg24 / xsum
      al27 = al27 / xsum
      si28 = si28 / xsum
      s32 = s32 / xsum
      ar36 = ar36 / xsum
      ca40 = ca40 / xsum
      ti44 = ti44 / xsum
      cr48 = cr48 / xsum
      cr60 = cr60 / xsum
      fe52 = fe52 / xsum
      fe54 = fe54 / xsum
      fe56 = fe56 / xsum
      co56 = co56 / xsum
      ni56 = ni56 / xsum
      xsum = h1 + he3 + he4 + c12 + n14 + o16 + ne20 + na23 + mg24 + al27 + si28 + &
            s32 + ar36 + ca40 + ti44 + cr48 + cr60 + fe52 + fe54 + fe56 + co56 + ni56
      if (xsum > 1.0001d0 .or. xsum < 0.9999d0) then
         write(*, 2) 'bad renormalized xsum', i, xsum
         stop 'write_stella_files'
      end if
      total_co56 = total_co56 + co56 * dm
      total_ni56 = total_ni56 + ni56 * dm
      if (is_bad(total_co56 + total_ni56)) then
         write(*, 2) 'total_co56 total_ni56', i, &
               total_co56, total_ni56, co56, ni56, dm
         stop 'write_stella_files'
      end if
   end do
   
   ! adjust mass fractions to get desired ni56 + co56
   f = (tot_co56 + tot_ni56) / (total_co56 + total_ni56)
   !write(*,1) 'f', f, (tot_co56 + tot_ni56)/Msun, (total_co56 + total_ni56)/Msun
   
   do i = 1, n
      o16_out(i) = o16_out(i) - (f - 1d0) * (co56_out(i) + ni56_out(i))
      if (o16_out(i) > 1d01 .or. o16_out(i) < 0d0) then
         write(*, 2) 'bad renormalized o16', i, o16_out(i)
         stop 'write_stella_files'
      end if
      co56_out(i) = f * co56_out(i)
      ni56_out(i) = f * ni56_out(i)
   end do
   
   r = 0
   mass = 0
   density = 0
   temperature = 0
   pressure = 0
   total_co56 = 0
   total_ni56 = 0
   
   mass = m_center
   do i = 1, n
      m_prev = mass
      mass = mmid_out(i)
      dm = mass - m_prev
      r = rmid_out(i)
      density = rho_out(i)
      temperature = temperature_out(i)
      velocity = max(min_velocity, u_out(i))
      pressure = pressure_out(i)
      energy = energy_out(i)
      entropy = entropy_out(i)
      opacity = opacity_out(i)
      tau = tau_out(i)
      v_div_csound = v_div_csound_out(i)
      write(io_hyd, '(i4,99(1pe24.16))') &
            i, mass / Msun, r, density, temperature, velocity, pressure, &
            0d0, energy, entropy, opacity, tau, v_div_csound
      
      if (i_h1 > 0) h1 = h1_out(i)
      if (i_he3 > 0) he3 = he3_out(i)
      if (i_he4 > 0) he4 = he4_out(i)
      if (i_c12 > 0) c12 = c12_out(i)
      if (i_n14 > 0) n14 = n14_out(i)
      if (i_o16 > 0) o16 = o16_out(i)
      if (i_ne20 > 0) ne20 = ne20_out(i)
      if (i_mg24 > 0) mg24 = mg24_out(i)
      if (i_si28 > 0) si28 = si28_out(i)
      if (i_s32 > 0) s32 = s32_out(i)
      if (i_ar36 > 0) ar36 = ar36_out(i)
      if (i_ca40 > 0) ca40 = ca40_out(i)
      if (i_ti44 > 0) ti44 = ti44_out(i)
      if (i_cr48 > 0) cr48 = cr48_out(i)
      if (i_cr60 > 0) cr60 = cr60_out(i)
      if (i_fe52 > 0) fe52 = fe52_out(i)
      if (i_fe54 > 0) fe54 = fe54_out(i)
      if (i_co56 > 0) co56 = co56_out(i)
      if (i_ni56 > 0) ni56 = ni56_out(i)
      na23 = f_na23 * n14
      al27 = f_al27 * n14
      n14 = n14 - na23 - al27
      
      xsum = h1 + he3 + he4 + c12 + n14 + o16 + ne20 + na23 + mg24 + al27 + si28 + &
            s32 + ar36 + ca40 + ti44 + cr48 + cr60 + fe52 + fe54 + fe56 + co56 + ni56
      if (xsum > 1.1d0 .or. xsum < 0.9d0) then
         write(*, 2) 'bad raw xsum', i, xsum
         stop 'write_stella_files'
      end if
      h1 = h1 / xsum
      he3 = he3 / xsum
      he4 = he4 / xsum
      c12 = c12 / xsum
      n14 = n14 / xsum
      o16 = o16 / xsum
      ne20 = ne20 / xsum
      na23 = na23 / xsum
      mg24 = mg24 / xsum
      al27 = al27 / xsum
      si28 = si28 / xsum
      s32 = s32 / xsum
      ar36 = ar36 / xsum
      ca40 = ca40 / xsum
      ti44 = ti44 / xsum
      cr48 = cr48 / xsum
      cr60 = cr60 / xsum
      fe52 = fe52 / xsum
      fe54 = fe54 / xsum
      fe56 = fe56 / xsum
      co56 = co56 / xsum
      ni56 = ni56 / xsum
      xsum = h1 + he3 + he4 + c12 + n14 + o16 + ne20 + na23 + mg24 + al27 + si28 + &
            s32 + ar36 + ca40 + ti44 + cr48 + cr60 + fe52 + fe54 + fe56 + co56 + ni56
      if (xsum > 1.0001d0 .or. xsum < 0.9999d0) then
         write(*, 2) 'bad renormalized xsum', i, xsum
         stop 'write_stella_files'
      end if
      total_co56 = total_co56 + co56 * dm
      total_ni56 = total_ni56 + ni56 * dm
      write(io_abn, '(i4,3(1pe13.6),99(1pe14.6))') &
            i, &
            mass / Msun, &
            dm / Msun, &
            r, &
            max(h1, minx), &
            max(he3 + he4, minx), &
            max(c12, minx), &
            max(n14, minx), &
            max(o16, minx), &
            max(ne20, minx), &
            max(na23, minx), &
            max(mg24, minx), &
            max(al27, minx), &
            max(si28, minx), &
            max(s32, minx), &
            max(ar36, minx), &
            max(ca40, minx), &
            max(fe52 + fe54 + fe56, minx), &
            minx, &
            max(co56 + ni56, minx)
   end do
   
   if (nz_extra > 0 .and. radius_wind > r) then
      velocity = v_CSM * 1d5 ! convert km/s to cm/s
      K_wind = mdot_for_CSM * Msun / secyer / (4 * pi * velocity)
      dr = (radius_wind - r) / (nz_extra - 1)
      do i = 1, nz_extra
         r_inner = r
         r = r_inner + dr
         vol = 4 * pi / 3 * (r * r * r - r_inner * r_inner * r_inner)
         density = K_wind / (r * r)
         dm = vol * density
         mass = mass + dm
         write(io_hyd, '(i4,99(1pe24.16))') &
               n + i, mass / Msun, r, density, temperature, velocity, pressure, &
               0d0, 0d0, 0d0, 0d0, 0d0, 0d0
         write(io_abn, '(i4,3(1pe13.6),99(1pe14.6))') &
               n + i, &
               mass / Msun, &
               dm / Msun, &
               r, &
               max(h1, minx), &
               max(he3 + he4, minx), &
               max(c12, minx), &
               max(n14, minx), &
               max(o16, minx), &
               max(ne20, minx), &
               max(na23, minx), &
               max(mg24, minx), &
               max(al27, minx), &
               max(si28, minx), &
               max(s32, minx), &
               max(ar36, minx), &
               max(ca40, minx), &
               max(fe52 + fe54 + fe56, minx), &
               minx, &
               max(co56 + ni56, minx)
      end do
   end if
   
   write(*, *)
   write(*, 1) 'stella total mass co56+ni56', (total_co56 + total_ni56) / Msun
   
   close(io_hyd)
   close(io_abn)
   
   deallocate(&
         model_ks, &
         mmid_mod, rmid_mod, tau_mod, u_mod, rho_mod, &
         pressure_mod, temperature_mod, energy_mod, &
         entropy_mod, opacity_mod, v_div_csound_mod, &
         h1_mod, he3_mod, he4_mod, c12_mod, n14_mod, &
         o16_mod, ne20_mod, mg24_mod, si28_mod, s32_mod, &
         ar36_mod, ca40_mod, ti44_mod, cr48_mod, cr60_mod, &
         fe52_mod, fe54_mod, fe56_mod, co56_mod, ni56_mod)
   deallocate(&
         interp_ks, &
         mmid_out, rmid_out, tau_out, u_out, rho_out, &
         pressure_out, temperature_out, energy_out, &
         entropy_out, opacity_out, v_div_csound_out, &
         h1_out, he3_out, he4_out, c12_out, n14_out, &
         o16_out, ne20_out, mg24_out, si28_out, s32_out, &
         ar36_out, ca40_out, ti44_out, cr48_out, cr60_out, &
         fe52_out, fe54_out, fe56_out, co56_out, ni56_out)
   deallocate(&
         mmid, rmid, taumid, work)

contains
   
   subroutine interp1(v_old, v_new)
      real(dp), pointer, dimension(:) :: v_old, v_new
      call interpolate_vector(&
            nz_mod, model_ks, n, interp_ks, v_old, v_new, &
            interp_pm, pm_work_size, work, 'write_stella_files', ierr)
      if (ierr /= 0) stop 'write_stella_files failed in interpolate_vector'
   end subroutine interp1


end subroutine write_stella_files


