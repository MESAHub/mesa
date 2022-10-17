! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the teerr of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras
   
   use star_lib
   use star_def
   use const_def
   use math_lib
   use auto_diff
   use utils_lib, only : mesa_error, is_bad
   
   implicit none
   
   !         x_ctrl(2) stop_m
   !         x_ctrl(3-6) magnetar - part 6 only
   !         x_ctrl(7) for setting stop_m in part2, fraction of He layer.
   !         x_ctrl(11-12) Ni mass - part 6 only
   !         x_ctrl(14) turn off RTI - part 6 only
   !         x_ctrl(15) length of time for Ni tail
   !         x_ctrl(16) dist below surface for stop_m
   !         x_ctrl(17) add csm - part 5 only
   !         x_ctrl(18-19) max timestep - part 6 only
   !         x_ctrl(21-23) alpha_MLT params
   !         x_ctrl(24-27) adjust kap factor - part 6 only
   !         x_ctrl(30) min_n_Fe for FeII velocity - part 6 only
   !         x_ctrl(31) tau_sob for FeII velocity - part 6 only
   !         x_ctrl(32-34) delta lgL timestep control
   !         x_ctrl(35-36) where to put Ni          num iters in x_integer_ctrl(2)
   !         x_ctrl(37) min center velocity for stella.
   !         x_ctrl(41-43) mass change - part 6 only
   !         x_ctrl(44) mass change - part 6 only
   !         x_ctrl(45-47) 1st smoothing of xa      num iters in x_integer_ctrl(3)
   !         x_ctrl(48-50) 2nd smoothing of xa      num iters in x_integer_ctrl(4)
   
   !         x_ctrl(71-75) extra energy deposition
   !         x_ctrl(98) forced stop_m
   !         x_ctrl(99) default stop_m
   
   include "test_suite_extras_def.inc"
   include 'stella/stella_def.inc'
   
   real(dp) :: &
      tp_photosphere, tp_L_eq_Lnuc, min_m_photosphere, initial_time, &
      initial_nico, initial_M_center, initial_he_core_mass, initial_mass, &
      start_m, stop_m
   
   integer, parameter :: num_logRhos = 41, num_logTs = 117, iounit = 33
   integer :: ilinx, iliny, ibcxmin, ibcxmax, ibcymin, ibcymax
   real(dp) :: bcxmin(num_logTs), bcxmax(num_logTs), Ts(num_logTs)
   real(dp) :: bcymin(num_logRhos), bcymax(num_logRhos)
   real(dp), pointer, dimension(:) :: logRhos, logTs, tau_sob_f1, &
      tau_sob_values, eta_i_values, n_Fe_values
   real(dp), pointer :: tau_sob_f(:, :, :)
   logical :: have_tau_sob_info_for_this_step

contains

include "test_suite_extras.inc"
   
   
   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      include 'stella/stella_controls.inc'
      if (ierr /= 0) return
      
      s% extras_startup => extras_startup
      s% extras_check_model => extras_check_model
      s% extras_start_step => extras_start_step
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns
      s% other_wind => low_density_wind_routine
      s% other_alpha_mlt => alpha_mlt_routine
   end subroutine extras_controls


include 'stella/stella.inc'
   
   
   subroutine alpha_mlt_routine(id, ierr)
      use chem_def, only : ih1
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k, h1
      real(dp) :: alpha_H, alpha_other, H_limit
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      alpha_H = s% x_ctrl(21)
      alpha_other = s% x_ctrl(22)
      H_limit = s% x_ctrl(23)
      h1 = s% net_iso(ih1)
      if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
      do k = 1, s% nz
         if (s% xa(h1, k) >= H_limit) then
            s% alpha_mlt(k) = alpha_H
         else
            s% alpha_mlt(k) = alpha_other
         end if
      end do
   end subroutine alpha_mlt_routine
   
   
   subroutine low_density_wind_routine(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
      real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: msum, lgrho_limit, lnd_limit
      integer :: k, i_lnd
      include 'formats'
      w = 0
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      i_lnd = s% i_lnd
      if (i_lnd <= 0) return
      lgrho_limit = -14.5d0 ! remove material from surface with density below this
      lnd_limit = ln10 * lgrho_limit
      msum = 0d0
      do k = 1, s% nz
         if (s% xh(i_lnd, k) >= lnd_limit) exit
         msum = msum + s% dm(k)
      end do
      if (msum == 0d0) return
      w = (msum / Msun) / (s% dt / secyer)
      write(*, 1) 'low_density_wind_routine lg(Mdot) msum/Msun', safe_log10(w), msum / Msun
   end subroutine low_density_wind_routine
   
   
   subroutine set_nico_mass(s, i_ni56, i_co56, new_ni, final_call, mass_ni56, ierr)
      use chem_def, only : io16
      type (star_info), pointer :: s
      integer, intent(in) :: i_ni56, i_co56
      logical, intent(in) :: final_call
      real(dp), intent(in) :: new_ni
      real(dp), intent(out) :: mass_ni56
      integer, intent(out) :: ierr
      
      integer :: i_o16, k, j, n, nz, species, jmax, kcut
      real(dp) :: old_nico, nico_change, &
         old_o16, new_o16, alfa_o16, alfa_nico, sum_dm, &
         xo16, dm, m00, mp1, mass_co56, mcut, &
         old_xnico, new_xnico, dxnico, x, f, old, new, xsum, &
         check_o16, check_ni56, check_co56, min_m, max_m, new_ni_frac
      include 'formats'
      ierr = 0
      i_o16 = s% net_iso(io16)
      if (i_o16 <= 0) then
         call mesa_error(__FILE__, __LINE__, 'need to have o16 in net for set_nico_mass')
      end if
      nz = s% nz
      species = s% species
      do k = 1, nz ! fixup abundances before revise
         do j = 1, species
            s% xa(j, k) = max(0d0, min(1d0, s% xa(j, k)))
         end do
         xsum = sum(s% xa(1:species, k))
         if (abs(xsum - 1d0) > 1d-12) then
            do j = 1, species
               s% xa(j, k) = s% xa(j, k) / xsum
            end do
         end if
      end do
      min_m = s% x_ctrl(35) * Msun + s% M_center
      max_m = s% x_ctrl(36) * Msun + s% he_core_mass * Msun
      
      if (s% u_flag) then
         do k = nz, 1, -1
            if (s% u(k) > s% x_ctrl(37)) then
               ! prepare for removal before give to Stella
               if (s% m(k) > min_m) then
                  min_m = s% m(k)
                  write(*, 1) 'increase min_m to ', min_m / Msun
                  exit
               end if
            end if
         end do
      end if
      
      write(*, 1) 'max_m min_m he_core_mass new_ni', &
         max_m / Msun, min_m / Msun, s% he_core_mass, new_ni
      if (max_m > 0d0) then
         do k = 1, nz ! replace ni56 + co56 by o16
            s% xa(i_o16, k) = s% xa(i_o16, k) + s% xa(i_ni56, k) + s% xa(i_co56, k)
            s% xa(i_co56, k) = 0d0
            s% xa(i_ni56, k) = 0d0
         end do
         max_m = min(max_m, s% m(1))
         if (min_m < 0) then
            min_m = max_m - new_ni * Msun
            mp1 = s% m(1)
            do k = 1, nz ! add ni56 from max_m to min_m
               m00 = mp1
               mp1 = m00 - s% dm(k)
               if (mp1 >= max_m) cycle
               if (m00 <= max_m .and. mp1 >= min_m) then
                  s% xa(1:species, k) = 0d0
                  s% xa(i_ni56, k) = 1d0
               else
                  dm = min(m00, max_m) - max(mp1, min_m)
                  f = max(0d0, min(1d0, 1d0 - dm / s% dm(k)))
                  do j = 1, species
                     s% xa(j, k) = f * s% xa(j, k)
                  end do
                  s% xa(i_ni56, k) = 0d0
                  s% xa(i_ni56, k) = 1d0 - sum(s% xa(1:species, k))
               end if
            end do
         else ! min_m >= 0
            min_m = max(min_m, s% M_center)
            if (min_m >= max_m) then
               write(*, 1) 'min_m >= max_m', min_m / Msun, max_m / Msun, s% M_center / Msun
               call mesa_error(__FILE__, __LINE__, 'set_nico_mass')
            end if
            new_ni_frac = new_ni * Msun / (max_m - min_m)
            write(*, 1) 'new_ni_frac min_m/Msun max_m/Msun', new_ni_frac, min_m / Msun, max_m / Msun
            mp1 = s% m(1)
            do k = 1, nz ! add ni56 from max_m to min_m
               m00 = mp1
               mp1 = m00 - s% dm(k)
               if (m00 <= min_m) exit
               if (mp1 >= max_m) cycle
               if (m00 <= max_m .and. mp1 >= min_m) then
                  dm = s% dm(k)
               else
                  dm = min(m00, max_m) - max(mp1, min_m)
               end if
               s% xa(i_ni56, k) = 0d0
               jmax = maxloc(s% xa(1:species, k), dim = 1)
               if (s% xa(jmax, k) < new_ni_frac) then
                  write(*, 2) 'too much ni for cell', k
                  call mesa_error(__FILE__, __LINE__, 'set_nico_mass')
               end if
               s% xa(i_ni56, k) = new_ni_frac * dm / s% dm(k)
               s% xa(jmax, k) = s% xa(jmax, k) - s% xa(i_ni56, k)
               if (is_bad(s% xa(i_ni56, k))) then
                  write(*, 2) 's% xa(i_ni56,k)', k, s% xa(i_ni56, k)
               end if
            end do
         end if
         mass_ni56 = dot_product(s% dm(1:nz), s% xa(i_ni56, 1:nz)) / Msun
         if (abs(mass_ni56 - new_ni) > 1d-6 * new_ni) then
            write(*, 1) 'bad new ni', new_ni, mass_ni56
            call mesa_error(__FILE__, __LINE__, 'set_nico_mass')
         end if
         write(*, 1) 'new ni56 mass', mass_ni56
         !call mesa_error(__FILE__,__LINE__,'set_nico_mass')
         return
      end if
      
      write(*, *) 'rescale Ni profile'
      do k = 1, nz
         do j = 1, species
            s% xa(j, k) = max(0d0, min(1d0, s% xa(j, k)))
         end do
         xsum = sum(s% xa(1:species, k))
         if (abs(xsum - 1d0) > 1d-12) then
            do j = 1, species
               s% xa(j, k) = s% xa(j, k) / xsum
            end do
         end if
      end do
      
      kcut = nz
      if (final_call .and. stella_skip_inner_dm > 0d0) then
         mcut = s% M_center + stella_skip_inner_dm * Msun
         do k = nz, 1, -1
            if (s% m(k) >= mcut .and. &
               (stella_skip_inner_v_limit < 0d0 .or. &
                  s% u(k) >= stella_skip_inner_v_limit)) then
               kcut = k
               exit
            end if
         end do
      end if
      
      do k = 1, nz ! replace co56 by ni56
         s% xa(i_ni56, k) = s% xa(i_ni56, k) + s% xa(i_co56, k)
         s% xa(i_co56, k) = 0d0
      end do
      ! only count cells 1:kcut
      old_o16 = dot_product(s% dm(1:kcut), s% xa(i_o16, 1:kcut)) / Msun
      mass_ni56 = dot_product(s% dm(1:kcut), s% xa(i_ni56, 1:kcut)) / Msun
      mass_co56 = dot_product(s% dm(1:kcut), s% xa(i_co56, 1:kcut)) / Msun
      old_nico = mass_ni56 + mass_co56
      nico_change = new_ni - old_nico
      if (new_ni < 0d0) then
         write(*, 1) 'new_ni', new_ni
         call mesa_error(__FILE__, __LINE__, 'bad mass set_nico_mass set_nico_mass')
      end if
      alfa_nico = new_ni / old_nico
      do k = 1, kcut
         xo16 = s% xa(i_o16, k)
         old_xnico = s% xa(i_ni56, k) + s% xa(i_co56, k)
         new_xnico = alfa_nico * old_xnico
         dxnico = new_xnico - old_xnico
         s% xa(i_co56, k) = 0d0
         s% xa(i_ni56, k) = new_xnico
         if (dxnico <= xo16) then
            s% xa(i_o16, k) = xo16 - dxnico
         else
            j = maxloc(s% xa(1:species, k), dim = 1)
            s% xa(j, k) = s% xa(j, k) - dxnico
            if (s% xa(j, k) <= 0d0) then
               write(*, 3) 'failed in set_nico_mass', j, k, s% xa(j, k), dxnico
               call mesa_error(__FILE__, __LINE__, 'set_nico_mass')
            end if
         end if
      end do
      check_ni56 = dot_product(s% dm(1:kcut), s% xa(i_ni56, 1:kcut)) / Msun
      if (abs(check_ni56 - new_ni) > 1d-10 * max(1d-10, new_ni)) then
         write(*, 1) 'check_ni56 - new_ni', check_ni56 - new_ni, check_ni56, new_ni
         write(*, 1) 'alfa_nico', alfa_nico
         call mesa_error(__FILE__, __LINE__, 'bad check mass in set_nico_mass')
      end if
      mass_ni56 = new_ni
      write(*, 1) 'revised mass Ni56', check_ni56
   end subroutine set_nico_mass
   
   
   subroutine create_sedona_Lbol_file()
      
      integer, parameter :: num_times = 205, num_freq = 2999, &   ! from 1st line of sedona data
         io1 = 34, io2 = 35
      integer :: j, k, ierr, n
      real(dp) :: Lbol, time0, freq_i00, Lnu_i00, time, freq_im1, Lnu_im1
      character (len = 256) :: file_in, file_out
      
      include 'formats'
      
      ierr = 0
      file_in = 'LOGS_comparison_apr10/sn_test_apr10_day20_no_mix_for_sedona.mod'
      file_out = 'LOGS_comparison_apr10/sn_test_apr10_day20_no_mix_for_sedona_Lbol.txt'
      
      open(unit = io1, file = trim(file_in), status = 'old', action = 'read', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) 'failed to open ' // trim(file_in)
         call mesa_error(__FILE__, __LINE__, 'create_sedona_Lbol_file')
         return
      end if
      open(io2, file = trim(file_out), action = 'write', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) 'failed to open ' // trim(file_out)
         return
      end if
      write(io2, '(a)') 'xxx day   log(Lbol)'
      
      read(io1, *) ! discard 1st line
      n = 0
      do j = 1, num_times
         read(io1, *) time0, freq_i00, Lnu_i00
         Lbol = 0d0
         do k = 1, num_freq - 1
            freq_im1 = freq_i00
            Lnu_im1 = Lnu_i00
            read(io1, *) time, freq_i00, Lnu_i00
            if (time /= time0) then
               write(*, *) 'time /= time0', time, time0
               call mesa_error(__FILE__, __LINE__, 'create_sedona_Lbol_file')
            end if
            Lbol = Lbol + 0.5d0 * (Lnu_i00 + Lnu_im1) * (freq_i00 - freq_im1)
         end do
         if (Lbol > 0d0) then
            write(io2, '(f8.2,f10.4)') time0 / secday, log10(Lbol)
            n = n + 1
         end if
      end do
      close(io1)
      close(io2)
      write(*, *) 'done'
      write(*, *) 'edit file to put number of times in 1st line', n
      write(*, *) trim(file_out)
      call mesa_error(__FILE__, __LINE__, 'create_sedona_Lbol_file')
   
   end subroutine create_sedona_Lbol_file
   
   
   subroutine test_numerical_diffusion(s)
      use chem_def, only : ih1
      type (star_info), pointer :: s
      integer :: i_h1, k
      include 'formats'
      i_h1 = s% net_iso(ih1)
      do k = 1, s% nz
         if (s% m(k) < 1.98d0 * Msun .or. s% m(k) > 2.02d0 * Msun) cycle
         s% xa(:, k) = 0d0
         s% xa(i_h1, k) = 1d0
      end do
   end subroutine test_numerical_diffusion
   
   
   subroutine extras_startup(id, restart, ierr)
      use chem_def, only : ini56, ico56, ih1, ihe4, io16
      use interp_2d_lib_db, only : interp_mkbicub_db
      use eos_lib, only : eosDT_get_T
      use eos_def
      use atm_lib, only : atm_L
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: xni56, xmax, &
         density, temperature, eta, logRho, logT, P_hse, Z, &
         logT_bnd1, logT_bnd2, logP_at_bnd1, logP_at_bnd2, &
         logT_guess, logT_result, logT_tol, logP_tol, &
         csm_mass, csm_mdot, windv, r0, rho0, T0, L0, r, f, rho, dm, dV, dq, &
         min_mass, max_mass, boxcar_mass, tot_h1, tot_he4
      integer :: i, j, k, kk, k1, k_max_v, ni56, co56, o16, he4, h1, &
         max_iter, eos_calls
      real(dp), dimension(num_eos_basic_results) :: &
         res, d_dlnd, d_dlnT
      real(dp), allocatable :: d_dxa(:, :)
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call test_suite_startup(s, restart, ierr)
      
      allocate(d_dxa(num_eos_d_dxa_results, s% species))
      
      xni56 = 0d0
      ni56 = s% net_iso(ini56)
      co56 = s% net_iso(ico56)
      o16 = s% net_iso(io16)
      he4 = s% net_iso(ihe4)
      h1 = s% net_iso(ih1)
      if (o16 <= 0 .or. he4 <= 0 .or. h1 <= 0) call mesa_error(__FILE__, __LINE__, 'missing o16, he4, or h1')
      
      if (s% eos_rq% logRho_min_OPAL_SCVH_limit > -12d0) then
         write(*, '(A)')
         write(*, *)'FIX: have set_logRho_OPAL_SCVH_limits too large'
         write(*, *) 'causes errors in HELM/OPAL blend partials'
         write(*, '(A)')
         call mesa_error(__FILE__, __LINE__, 'extras_startup')
      end if
      
      if (restart) then
         call unpack_extra_info(s)
      else
         
         initial_nico = 0
         min_m_photosphere = 1d99
         tp_photosphere = 0
         tp_L_eq_Lnuc = 0
         stop_m = 0
         initial_M_center = s% M_center
         initial_mass = s% star_mass
         initial_time = s% time
         initial_he_core_mass = s% he_core_mass
         call alloc_extra_info(s)
         
         if (s% x_integer_ctrl(1) == 0) then
            if (s% total_mass_for_inject_extra_ergs_sec > 0) then ! doing edep
               if (s% v_flag) then
                  do k = 1, s% nz
                     s% xh(s% i_v, k) = 0d0
                     s% v(k) = 0d0
                  end do
               end if
            end if
            return
         end if
         
         if (s% x_integer_ctrl(1) == 5 .and. &
            s% x_ctrl(17) > 0) call add_csm ! part 5, add csm
         
         if (s% x_ctrl(47) > 0d0 .and. s% x_integer_ctrl(3) > 0) then
            min_mass = s% x_ctrl(45) + s% M_center / Msun
            max_mass = s% x_ctrl(46) + s% he_core_mass
            boxcar_mass = s% x_ctrl(47)
            call smooth_xa_by_boxcar_mass(&
               s% id, min_mass, max_mass, boxcar_mass, s% x_integer_ctrl(3), ierr)
            if (ierr /= 0) return
         end if
         
         if (ni56 > 0 .and. co56 > 0) then
            if (s% x_ctrl(12) > 0) then
               call set_nico_mass(s, ni56, co56, s% x_ctrl(12), .false., xni56, ierr)
               if (ierr /= 0) return
            end if
            initial_nico = xni56
         end if
         
         if (s% x_ctrl(50) > 0d0 .and. s% x_integer_ctrl(4) > 0) then
            min_mass = s% x_ctrl(48) + s% M_center / Msun
            max_mass = s% x_ctrl(49) + s% he_core_mass
            boxcar_mass = s% x_ctrl(50)
            call smooth_xa_by_boxcar_mass(&
               s% id, min_mass, max_mass, boxcar_mass, s% x_integer_ctrl(4), ierr)
            if (ierr /= 0) return
         end if
         
         if (s% x_integer_ctrl(1) == 1) &
            s% cumulative_energy_error = 0d0 ! set to 0 at start of part1
         
         start_m = s% shock_mass
         if (start_m == 0d0) then ! use max v
            k_max_v = maxloc(s% u(1:s% nz), dim = 1)
            start_m = s% m(k_max_v) / Msun
            if (start_m < 2d0 * s% M_center / Msun) then
               write(*, 2) 'no shock; use max_v for start_m', k_max_v, start_m
            else
               start_m = s% M_center / Msun
               write(*, 2) 'no shock; use M_center for start_m', s% nz, start_m
            end if
         else
            write(*, 2) 'use shock for start_m', s% shock_k, s% shock_mass
         end if
         write(*, 2) 'M_center', s% nz, s% M_center / Msun
         
         write(*, 1) 's% x_ctrl(98)', s% x_ctrl(98)
         write(*, 2) 's% x_integer_ctrl(1)', s% x_integer_ctrl(1)
         
         if (s% x_ctrl(98) > 0d0) then
            stop_m = s% x_ctrl(98)
         else
            stop_m = 0d0
            select case(s% x_integer_ctrl(1))
            case (1)
               do k = 1, s% nz
                  if (s% xa(o16, k) > s% xa(he4, k)) then
                     stop_m = (0.15d0 * s% M_center + 0.85d0 * s% m(k)) / Msun
                     if (stop_m <= start_m) then
                        write(*, 2) '1st choice for stop_m too small', k, stop_m
                        stop_m = (0.05d0 * s% M_center + 0.95d0 * s% m(k)) / Msun
                        if (stop_m <= start_m) then
                           write(*, 2) '2nd choice for stop_m too small', k, stop_m
                           stop_m = (0.01d0 * s% M_center + 0.99d0 * s% m(k)) / Msun
                           if (stop_m <= start_m) then
                              write(*, 2) '3rd choice for stop_m too small', k, stop_m
                              stop_m = 1.001d0 * start_m
                           end if
                        end if
                     end if
                     exit
                  end if
               end do
               if (stop_m == 0d0) then
                  if (s% x_ctrl(99) <= 0d0) call mesa_error(__FILE__, __LINE__, 'failed to find stop_m')
                  stop_m = s% x_ctrl(99)
               end if
            case (2)
               k1 = 0
               do k = 1, s% nz
                  if (s% xa(he4, k) > s% xa(h1, k)) then
                     k1 = k
                     exit
                  end if
               end do
               if (k1 > 0) then
                  do k = k1 + 1, s% nz
                     if (s% xa(o16, k) > s% xa(he4, k)) then
                        f = s% x_ctrl(7)
                        stop_m = ((1d0 - f) * s% m(k) + f * s% m(k1)) / Msun
                        exit
                     end if
                  end do
               end if
               if (stop_m == 0d0) then
                  if (s% x_ctrl(99) <= 0d0) call mesa_error(__FILE__, __LINE__, 'failed to find stop_m')
                  stop_m = s% x_ctrl(99)
               end if
            case (3)
               do k = 1, s% nz
                  if (s% xa(he4, k) > s% xa(h1, k)) then
                     stop_m = 0.95d0 * s% m(k) / Msun + 0.05d0 * (s% star_mass - 0.1d0)
                     exit
                  end if
               end do
               if (stop_m == 0d0) call mesa_error(__FILE__, __LINE__, 'failed to find stop_m')
            case (4)
               stop_m = s% star_mass - s% x_ctrl(16)
            case (5)
               if (s% x_ctrl(16) > 0d0) &
                  stop_m = s% star_mass - s% x_ctrl(16)
            case (6)
            end select
         end if
      
      end if ! not restart
      
      write(*, 1) 's% x_ctrl(16)', s% x_ctrl(16)
      write(*, 1) 's% star_mass', s% star_mass
      write(*, 1) 'start_m', start_m
      write(*, 1) 'stop_m', stop_m
      
      !if (stop_m < start_m) call mesa_error(__FILE__,__LINE__,'bad stop_m')
      
      if (s% x_ctrl(16) > 0d0) &
         stop_m = min(stop_m, s% star_mass - s% x_ctrl(16))
      
      if (start_m > stop_m .and. stop_m > 0d0) then
         write(*, 1) 'start_m > stop_m', start_m, stop_m
         call mesa_error(__FILE__, __LINE__, 'extras_startup')
      end if
      
      if (stop_m > 0d0) then
         s% x_ctrl(2) = stop_m
         write(*, 1) 'stop when shock reaches', stop_m
         write(*, 1) 's% he_core_mass', s% he_core_mass
         write(*, 1) 's% co_core_mass', s% co_core_mass
         write(*, 1) 's% M_center/Msun', s% M_center / Msun
         write(*, 1) 's% star_mass', s% star_mass
         write(*, '(A)')
         !stop
      end if
      
      if (s% x_integer_ctrl(1) == 6 .and. s% x_logical_ctrl(6)) then
         ! setup interpolation table for tau sob eta_i
         open(unit = iounit, file = 'FeII_5169_eta.dat', action = 'read')
         allocate(logRhos(num_logRhos), logTs(num_logTs), &
            tau_sob_f1(4 * num_logRhos * num_logTs))
         tau_sob_f(1:4, 1:num_logRhos, 1:num_logTs) => &
            tau_sob_f1(1:4 * num_logRhos * num_logTs)
         do j = 1, num_logRhos
            do i = 1, num_logTs
               read(iounit, *) density, temperature, eta
               logRho = log10(density)
               logRhos(j) = logRho
               logT = log10(temperature)
               if (j == 1) then
                  Ts(i) = temperature
                  logTs(i) = logT
               else if (logT /= logTs(i)) then
                  write(*, 3) 'bad T?', i, j, Ts(1), temperature, density, eta
                  call mesa_error(__FILE__, __LINE__, 'table error?')
               end if
               tau_sob_f(1, j, i) = eta
            end do
         end do
         close(iounit)
         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(1:num_logTs) = 0
         ibcxmax = 0; bcxmax(1:num_logTs) = 0
         ibcymin = 0; bcymin(1:num_logRhos) = 0
         ibcymax = 0; bcymax(1:num_logRhos) = 0
         call interp_mkbicub_db(&
            logRhos, num_logRhos, logTs, num_logTs, tau_sob_f1, num_logRhos, &
            ibcxmin, bcxmin, ibcxmax, bcxmax, &
            ibcymin, bcymin, ibcymax, bcymax, &
            ilinx, iliny, ierr)
         if (ierr /= 0) then
            write(*, *) 'interp_mkbicub_db error'
            ierr = -1
            return
         end if
         !write(*,*) 'done with setup for tau_sob eta interpolation'
         allocate(&
            tau_sob_values(s% nz + 100), eta_i_values(s% nz + 100), n_Fe_values(s% nz + 100))
      end if
      
      have_tau_sob_info_for_this_step = .false.
   
   contains
      
      subroutine add_csm
         include 'formats'
         csm_mass = s% x_ctrl(17) * Msun
         csm_mdot = s% x_ctrl(16) * Msun / secyer
         do kk = 2, s% nz - 2
            if (s% m(1) - s% m(kk) < csm_mass) cycle
            r0 = s% r(kk)
            if (s% x_ctrl(28) > 0) then
               rho0 = s% x_ctrl(28)
            else
               rho0 = s% rho(kk)
            end if
            if (s% x_ctrl(29) > 0) then
               T0 = s% x_ctrl(29)
            else
               T0 = s% T(kk)
            end if
            L0 = s% L(1)
            windv = csm_mdot / (4 * pi * r0 * r0 * rho0) ! mdot velocity
            !windv = sqrt(2*s% cgrav(kk)*s% m(kk)/r0) ! escape velocity
            write(*, 1) 'old log(r(1)/Rsun), R/Rsun', log10(s% r(1) / Rsun), s% r(1) / Rsun
            write(*, 2) 'rho0, T0, r0, csm mass, csm v', &
               kk, rho0, T0, r0 / Rsun, (s% m(1) - s% m(kk)) / Msun, windv
            dm = sum(s% dm(1:kk - 1)) / (kk - 1)
            dq = dm / s% xmstar
            do k = 1, kk - 1
               s% dq(k) = dq
               s% q(k + 1) = s% q(k) - dq
               s% dm(k) = dm
               s% m(k + 1) = s% m(k) - dm
            end do
            logT_bnd1 = arg_not_provided
            logT_bnd2 = arg_not_provided
            logP_at_bnd1 = arg_not_provided
            logP_at_bnd2 = arg_not_provided
            max_iter = 100
            logT_tol = 1d-5
            logP_tol = 1d-8
            k = kk
            write(*, '(A)')
            write(*, 2) 'T rho P logR u', k, s% T(k), s% rho(k), s% Peos(k), &
               log10(s% r(k) / Rsun), s% u(k)
            do k = kk - 1, 1, -1
               r = s% r(k + 1) + 0.5d0 * (s% r(k + 1) - s% r(k + 2))
               f = r0 / r; f = f * f
               rho = rho0 * f ! rho proportional to 1/r^2
               dV = dm / rho
               r = s% r(k + 1)
               r = pow(dV / (4d0 * pi / 3d0) + r * r * r, 1d0 / 3d0)
               s% r(k) = r
               s% lnR(k) = log(s% r(k))
               s% xh(s% i_lnR, k) = s% lnR(k)
               s% rho(k) = rho
               s% lnd(k) = log(s% rho(k))
               s% xh(s% i_lnd, k) = s% lnd(k)
               
               s% u(k) = windv !* r/r0
               
               s% xh(s% i_u, k) = s% u(k)
               
               if (.true.) then ! set T to give P for HSE
                  r = s% r(k + 1)
                  P_hse = s% Peos(k + 1) - &
                     s% cgrav(k + 1) * s% m(k + 1) * (s% dm(k + 1) + s% dm(k)) / (8 * pi * r * r * r * r)
                  Z = max(0d0, min(1d0, 1d0 - (s% X(k) + s% Y(k))))
                  logT_guess = s% lnT(k + 1) / ln10
                  call eosDT_get_T(&
                     s% eos_handle, &
                     s% species, s% chem_id, s% net_iso, s% xa(:, k), &
                     s% lnd(k) / ln10, i_logPtot, log10(P_hse), &
                     logT_tol, logP_tol, max_iter, logT_guess, &
                     logT_bnd1, logT_bnd2, logP_at_bnd1, logP_at_bnd2, &
                     logT_result, res, d_dlnd, d_dlnT, &
                     d_dxa, eos_calls, ierr)
                  if (ierr /= 0) return
                  s% lnT(k) = logT_result * ln10
                  s% T(k) = exp(s% lnT(k))
                  if (.false. .and. s% T(k) > s% T(k + 1)) then ! can happen at base
                     s% T(k) = s% T(k + 1)
                     s% lnT(k) = s% lnT(k + 1)
                  end if
               else
                  s% T(k) = T0 * f
                  s% lnT(k) = log(s% T(k))
               end if
               s% xh(s% i_lnT, k) = s% lnT(k)
               
               if (.true.) then ! set to black body L
                  s% L(k) = atm_L(s% T(k), s% r(k))
               else
                  s% L(k) = L0
               end if
               s% xh(s% i_lum, k) = s% L(k)
               
               write(*, 2) 'T rho P logR u', k, s% T(k), s% rho(k), s% Peos(k), &
                  log10(s% r(k) / Rsun), s% u(k)
            end do
            write(*, '(A)')
            write(*, 1) 'new log(r(1)/Rsun), R/Rsun', log10(s% r(1) / Rsun), s% r(1) / Rsun
            write(*, '(A)')
            exit
         end do
      end subroutine add_csm
   
   end subroutine extras_startup
   
   
   subroutine extras_after_evolve(id, ierr)
      use chem_def, only : ini56, ico56
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: dt
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (s% x_integer_ctrl(1) == 5 .and. &
         save_stella_data_when_terminate) then
         call write_stella_data(s, ierr)
         if (ierr /= 0) return
      end if
      call test_suite_after_evolve(s, ierr)
   end subroutine extras_after_evolve
   
   
   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      use chem_def, only : ini56
      integer, intent(in) :: id
      integer :: ierr, k
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going
      have_tau_sob_info_for_this_step = .false.
   end function extras_check_model
   
   
   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 0
      if (s% x_integer_ctrl(1) == 6 .and. s% x_logical_ctrl(6)) &
         how_many_extra_history_columns = 6
   end function how_many_extra_history_columns
   
   
   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k, k_tau, k_lum
      integer :: iounit
      real(dp) :: tau_vel, alfa, tau_vel_FeII, &
         r_lum, m_lum, t_sum, dt_sum, dr, E_sum, Lnuc_sum
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (n == 0) return
      if (n /= 6) call mesa_error(__FILE__, __LINE__, 'bad num cols for data_for_extra_history_columns')
      names(1) = 'vel_FeII'
      names(2) = 'vel_FeII_cell'
      names(3) = 'vel_FeII_eta_i'
      names(4) = 'vel_FeII_n_Fe'
      names(5) = 'tau100_L'
      names(6) = 'tau100_logL'
      if (s% doing_relax) then
         vals = 0d0
         return
      end if
      call get_tau_sob_info_for_this_step(s)
      tau_vel_FeII = s% x_ctrl(31)
      k_tau = s% photosphere_cell_k
      tau_vel = s% photosphere_v
      if (.false. .and. tau_sob_values(1) >= tau_vel_FeII) then
         tau_vel = s% u(1)
         k_tau = 1
      else
         do k = 2, s% nz
            if (tau_sob_values(k) > tau_vel_FeII) then
               alfa = (tau_vel_FeII - tau_sob_values(k - 1)) / (tau_sob_values(k) - tau_sob_values(k - 1))
               tau_vel = s% u(k - 1) + (s% u(k) - s% u(k - 1)) * alfa
               k_tau = k
               exit
            end if
         end do
      end if
      vals(1) = tau_vel * 1d-5
      vals(2) = k_tau
      vals(3) = eta_i_values(k_tau)
      vals(4) = n_Fe_values(k_tau)
      vals(5:6) = 0
      do k = 2, s% nz
         if (s% tau(k) > 100) then
            alfa = (100 - s% tau(k - 1)) / (s% tau(k) - s% tau(k - 1))
            vals(5) = s% L(k - 1) + (s% L(k) - s% L(k - 1)) * alfa
            vals(6) = safe_log10(vals(5))
            exit
         end if
      end do
      open(newunit = iounit, file = 'tau100_L.data', &
         status = "unknown", form = 'formatted', position = "append")
      write(iounit, *) s% time / secday, vals(5), vals(6)
      close(iounit)
   
   end subroutine data_for_extra_history_columns
   
   
   subroutine get_tau_sob_info_for_this_step(s)
      use const_def, only : dp, avo, pi, qe, me, clight
      use chem_def, only : ife56
      use interp_2d_lib_db, only : interp_evbicub_db
      type (star_info), pointer :: s
      
      integer :: k, ict(6), fe56, nz, ierr, j
      real(dp) :: logRho, logT, tau_sob, eta_i, n_Fe, fval(6), time, smooth(s% nz)
      real(dp), parameter :: A_Fe56 = 56d0, lambda0 = 5169.02d-8, f = 0.023d0, t0 = 0d0 ! 8.8d4
      include 'formats'
      
      nz = s% nz
      if (.not. associated(tau_sob_values)) then
         allocate(tau_sob_values(nz + 100), eta_i_values(nz + 100), n_Fe_values(nz + 100))
         have_tau_sob_info_for_this_step = .false.
      else if (size(tau_sob_values, dim = 1) < nz) then
         deallocate(tau_sob_values, eta_i_values, n_Fe_values)
         allocate(tau_sob_values(nz + 100), eta_i_values(nz + 100), n_Fe_values(nz + 100))
         have_tau_sob_info_for_this_step = .false.
      end if
      
      if (have_tau_sob_info_for_this_step) return
      
      fe56 = s% net_iso(ife56)
      ict = 0
      ict(1) = 1
      !write(*,1) 'pi*qe*qe/(me*clight)', pi*qe*qe/(me*clight)
      do k = 1, nz
         smooth(k) = s% rho(k)
      end do
      do j = 1, 5
         do k = 2, nz - 1
            smooth(k) = sum(smooth(k - 1:k + 1)) / 3d0
         end do
      end do
      do k = 1, nz
         n_Fe = smooth(k) * avo * s% xa(fe56, k) / A_Fe56
         logRho = log10(smooth(k))
         logRho = min(logRhos(num_logRhos), max(logRhos(1), logRho))
         logT = min(logTs(num_logTs), max(logTs(1), s% lnT(k) / ln10))
         ierr = 0
         call interp_evbicub_db(&
            logRho, logT, logRhos, num_logRhos, logTs, num_logTs, &
            ilinx, iliny, tau_sob_f1, num_logRhos, ict, fval, ierr)
         if (ierr /= 0) then
            write(*, 2) 'logRho', k, s% lnd(k) / ln10
            write(*, 2) 'logT', k, s% lnT(k) / ln10
            call mesa_error(__FILE__, __LINE__, 'interp failed in data_for_extra_profile_columns')
         end if
         eta_i = fval(1)
         if (.false.) then
            write(*, 1) 'logT', logT
            write(*, 1) 'logTs(1)', logTs(1)
            write(*, 1) 'logTs(num_logTs)', logTs(num_logTs)
            write(*, 1) 'logRho', logRho
            write(*, 1) 'logRhos(1)', logRhos(1)
            write(*, 1) 'logRhos(num_logRhos)', logRhos(num_logRhos)
            write(*, 1) 'eta_i', eta_i
            write(*, 1) 's% time', s% time
            call mesa_error(__FILE__, __LINE__, 'data_for_extra_profile_columns')
         end if
         time = s% time + t0
         tau_sob = pi * qe * qe / (me * clight) * n_Fe * eta_i * f * time * lambda0
         tau_sob_values(k) = tau_sob
         eta_i_values(k) = eta_i
         n_Fe_values(k) = n_Fe
      end do
      
      do k = 1, nz
         smooth(k) = tau_sob_values(k)
      end do
      do j = 1, 5
         do k = 2, nz - 1
            smooth(k) = sum(smooth(k - 1:k + 1)) / 3d0
         end do
      end do
      do k = 1, nz
         tau_sob_values(k) = smooth(k)
      end do
      
      have_tau_sob_info_for_this_step = .true.
   
   end subroutine get_tau_sob_info_for_this_step
   
   
   integer function how_many_extra_profile_columns(id)
      use star_def, only : star_info
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 0
      if (s% x_integer_ctrl(1) == 3) &
         how_many_extra_profile_columns = 1
      if (s% x_integer_ctrl(1) == 6 .and. s% x_logical_ctrl(6)) &
         how_many_extra_profile_columns = 3
   end function how_many_extra_profile_columns
   
   
   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      use star_def, only : star_info, maxlen_profile_column_name
      use const_def, only : dp
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      real(dp) :: dtau, kap
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (n == 0) return
      !call mesa_error(__FILE__,__LINE__,'data_for_extra_profile_columns')
      if (s% x_integer_ctrl(1) == 3) then
         if (n /= 1) call mesa_error(__FILE__, __LINE__, 'bad num cols for data_for_extra_profile_columns')
         names(1) = 'du'
         vals(1, 1) = 0
         do k = 2, s% nz
            vals(k, 1) = 0 ! s% xtra1_array(k)
         end do
         return
      end if
      if (n /= 3) call mesa_error(__FILE__, __LINE__, 'bad num cols for data_for_extra_profile_columns')
      names(1) = 'dlnT4_dtau'
      names(2) = 'v_times_3_div_c'
      names(3) = 'log_ratio'
      if (s% doing_relax) then
         vals = 0d0
         return
      end if
      do k = 2, s% nz - 1
         kap = 0.5d0 * (s% opacity(k) + s% opacity(k - 1))
         dtau = -s% dm_bar(k) * kap / (4 * pi * s% r(k) * s% r(k))
         vals(k, 1) = (s% lnT(k) - s% lnT(k + 1)) / dtau
         vals(k, 2) = s% u_face_ad(k)%val * 3d0 / clight
         vals(k, 3) = safe_log10(vals(k, 1) / max(1d-99, vals(k, 2)))
      end do
      vals(s% nz, 1:3) = vals(s% nz - 1, 1:3)
      vals(1, 1:3) = vals(2, 1:3)
   end subroutine data_for_extra_profile_columns
   
   
   ! returns either keep_going or terminate.
   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      real(dp) :: alfa, beta, age_days, L_center, next_R_center, &
         start_ramp_up, end_ramp_up, &
         start_ramp_down, end_ramp_down, factor_end, factor_start
      
      include 'formats'
      extras_start_step = keep_going
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      if (s% x_integer_ctrl(1) == 1 .and. s% model_number >= 1000) &
         s% max_timestep = 0 ! turn off limit
      
      age_days = s% star_age * 365.25d0
      
      if (s% x_ctrl(14) > 0 .and. age_days >= s% x_ctrl(14) .and. &
         s% RTI_C > 0d0) then
         s% RTI_C = 0d0
         s% RTI_log_max_boost = 0d0
         s% RTI_m_full_boost = -1d0
         s% RTI_m_no_boost = 0d0
         s% dedt_RTI_diffusion_factor = 1d0
         s% composition_RTI_diffusion_factor = 1d0
         write(*, '(A)')
         write(*, 2) 'turn off RTI', s% model_number, age_days
         write(*, '(A)')
      end if
      
      L_center = s% x_ctrl(3)
      if (L_center > 0) then  ! magnetar
         start_ramp_up = s% x_ctrl(4)
         end_ramp_up = s% x_ctrl(5)
         start_ramp_down = s% x_ctrl(6)
         end_ramp_down = s% x_ctrl(7)
         if (age_days <= start_ramp_up) then
            s% L_center = 0d0
         else if (age_days < end_ramp_up) then
            s% L_center = &
               L_center * (age_days - start_ramp_up) / (end_ramp_up - start_ramp_up)
         else if (age_days < start_ramp_down) then
            s% L_center = L_center
         else if (age_days <= end_ramp_down) then
            s% L_center = &
               L_center * (end_ramp_down - age_days) / (end_ramp_down - start_ramp_down)
         else
            s% L_center = 0d0
         end if
      end if
      
      if (s% x_integer_ctrl(1) == 6) then
         
         if (s% x_ctrl(41) > 0d0) then
            if (s% rho(1) < s% x_ctrl(41) .and. age_days < s% x_ctrl(44)) then
               if (s% mass_change == 0) then
                  s% mass_change = s% x_ctrl(42)
                  s% mass_depth_for_L_surf = s% x_ctrl(43)
                  write(*, 2) 'turn on mass_change', s% model_number, s% mass_change
               end if
            else if (s% mass_change /= 0) then
               s% mass_change = 0d0
               s% mass_depth_for_L_surf = 0d0
               write(*, 2) 'turn off mass_change', s% model_number
            end if
         end if
         
         if (s% x_ctrl(24) > 0) then ! adjust kap factor
            start_ramp_up = s% x_ctrl(24)
            end_ramp_up = s% x_ctrl(25)
            factor_start = s% x_ctrl(26)
            factor_end = s% x_ctrl(27)
            if (age_days <= start_ramp_up) then
               s% opacity_factor = factor_start
            else if (age_days >= end_ramp_up) then
               s% opacity_factor = factor_end
            else
               s% opacity_factor = factor_start + &
                  (factor_end - factor_start) * &
                     (age_days - start_ramp_up) / (end_ramp_up - start_ramp_up)
               write(*, 2) 'opacity_factor age', s% model_number, s% opacity_factor, age_days
            end if
         end if
      
      end if
      
      if (s% x_ctrl(18) > 0d0 .and. s% x_ctrl(19) > 0d0 .and. &
         s% time > s% x_ctrl(18) .and. &
         abs(s% max_timestep - s% x_ctrl(19)) > 1d-6 * s% x_ctrl(19)) then
         write(*, 2) 'change max_timestep', s% model_number, &
            s% time, s% x_ctrl(18), s% x_ctrl(19), &
            s% max_timestep - s% x_ctrl(19)
         s% max_timestep = s% x_ctrl(19)
         !s% okay_to_remesh = .true.
      end if
      
      if (age_days >= s% x_ctrl(32)) then
         s% delta_lgL_limit = s% x_ctrl(33)
         s% delta_lgL_hard_limit = s% x_ctrl(34)
      end if
      
      if (s% x_logical_ctrl(1) .and. s% dt > 0d0) then
         s% v_center = &
            s% v_center - s% cgrav(s% nz) * s% M_center / (s% R_center * s% R_center)
         s% v_center = max(s% v_center, -s% x_ctrl(1) * clight)
         next_R_center = s% R_center + s% v_center * s% dt
         if (next_R_center < s% center_R_lower_limit) then
            s% v_center = (s% center_R_lower_limit - s% R_center) / s% dt
         end if
      end if
   
   end function extras_start_step
   
   
   ! returns either keep_going or terminate.
   integer function extras_finish_step(id)
      use chem_def, only : ih1, ico56, ini56
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      real(dp) :: age_days, log_L, shock_mass, log_Lnuc_burn, xmax, xni56
      integer :: k
      include 'formats'
      extras_finish_step = keep_going
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call store_extra_info(s)
      
      if (s% x_ctrl(2) <= 0) return
      shock_mass = s% shock_mass
      if (shock_mass >= s% x_ctrl(2)) then
         write(*, '(a,2f12.5)') 'shock has reached target location', &
            shock_mass, s% x_ctrl(2)
         extras_finish_step = terminate
         s% termination_code = t_extras_finish_step
         ! restore Ni+Co mass
         if (s% net_iso(ini56) > 0 .and. s% net_iso(ico56) > 0) then
            if (s% x_ctrl(12) > 0) then
               call set_nico_mass(&
                  s, s% net_iso(ini56), s% net_iso(ico56), s% x_ctrl(12), .true., xni56, ierr)
               if (ierr /= 0) return
            end if
         end if
      else if (shock_mass >= 0.9995d0 * s% x_ctrl(2)) then
         write(*, 1) 'shock has reached this fraction of target', &
            shock_mass / s% x_ctrl(2)
      end if
      
      if (s% x_integer_ctrl(1) == 5 .and. &
         s% model_number == save_stella_data_for_model_number) then
         call write_stella_data(s, ierr)
         if (ierr /= 0) return
      end if
   
   end function extras_finish_step
   
   
   ! routines for saving and restoring extra data so can do restarts
   
   ! put these defs at the top and delete from the following routines
   !integer, parameter :: extra_info_alloc = 1
   !integer, parameter :: extra_info_get = 2
   !integer, parameter :: extra_info_put = 3
   
   
   subroutine alloc_extra_info(s)
      integer, parameter :: extra_info_alloc = 1
      type (star_info), pointer :: s
      call move_extra_info(s, extra_info_alloc)
   end subroutine alloc_extra_info
   
   
   subroutine unpack_extra_info(s)
      integer, parameter :: extra_info_get = 2
      type (star_info), pointer :: s
      call move_extra_info(s, extra_info_get)
   end subroutine unpack_extra_info
   
   
   subroutine store_extra_info(s)
      integer, parameter :: extra_info_put = 3
      type (star_info), pointer :: s
      call move_extra_info(s, extra_info_put)
   end subroutine store_extra_info
   
   
   subroutine move_extra_info(s, op)
      integer, parameter :: extra_info_alloc = 1
      integer, parameter :: extra_info_get = 2
      integer, parameter :: extra_info_put = 3
      type (star_info), pointer :: s
      integer, intent(in) :: op
      
      integer :: i, j, num_ints, num_dbls, ierr
      
      include 'formats'
      
      i = 0
      ! call move_int or move_flg
      num_ints = i
      
      i = 0
      ! call move_dbl
      call move_dbl(initial_nico)
      call move_dbl(initial_M_center)
      call move_dbl(initial_mass)
      call move_dbl(initial_time)
      call move_dbl(initial_he_core_mass)
      call move_dbl(tp_photosphere)
      call move_dbl(tp_L_eq_Lnuc)
      call move_dbl(min_m_photosphere)
      call move_dbl(start_m)
      call move_dbl(stop_m)
      
      num_dbls = i
      
      if (op /= extra_info_alloc) return
      if (num_ints == 0 .and. num_dbls == 0) return
      
      ierr = 0
      call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in star_alloc_extras'
         write(*, *) 'alloc_extras num_ints', num_ints
         write(*, *) 'alloc_extras num_dbls', num_dbls
         call mesa_error(__FILE__, __LINE__, 'example_cccsn_IIp alloc failed')
         !call mesa_error(__FILE__,__LINE__)
      end if
   
   contains
      
      subroutine move_dbl(dbl)
         real(dp) :: dbl
         i = i + 1
         select case (op)
         case (extra_info_get)
            dbl = s% extra_work(i)
         case (extra_info_put)
            s% extra_work(i) = dbl
         end select
      end subroutine move_dbl
      
      subroutine move_int(int)
         integer :: int
         i = i + 1
         select case (op)
         case (extra_info_get)
            int = s% extra_iwork(i)
         case (extra_info_put)
            s% extra_iwork(i) = int
         end select
      end subroutine move_int
      
      subroutine move_flg(flg)
         logical :: flg
         i = i + 1
         select case (op)
         case (extra_info_get)
            flg = (s% extra_iwork(i) /= 0)
         case (extra_info_put)
            if (flg) then
               s% extra_iwork(i) = 1
            else
               s% extra_iwork(i) = 0
            end if
         end select
      end subroutine move_flg
   
   end subroutine move_extra_info


end module run_star_extras
      
