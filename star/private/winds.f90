! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************


      module winds

      use star_private_def
      use const_def
      use chem_def, only: ih1, ihe4, ic12, ic13, in14, io16
      use utils_lib, only: is_bad

      implicit none

      private
      public :: set_mdot

      contains


      subroutine set_mdot(s, L_phot, M_phot, T_phot, ierr)
         use chem_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: L_phot, M_phot, T_phot ! photosphere values (cgs)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         call do_set_mdot(s, L_phot, M_phot, T_phot, ierr)
         if (ierr /= 0) return
         if (s% use_other_adjust_mdot) then
            call s% other_adjust_mdot(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*, *) 'other_adjust_mdot'
               return
            end if
         end if         
      end subroutine set_mdot
      

      subroutine do_set_mdot(s, L_phot, M_phot, T_phot, ierr)
         use chem_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: L_phot, M_phot, T_phot ! photosphere values (cgs)
         integer, intent(out) :: ierr

         integer :: k, j, h1, he4, nz, base
         real(dp) :: max_ejection_mass, wind_mdot, wind, alfa, total_H, &
            X, Y, Z, w1, w2, T_high, T_low, L1, M1, R1, T1, &
            log_dtyr, log_dtyr_full_off, log_dtyr_full_on, beta, divisor, &
            center_h1, center_he4, surface_h1, surface_he4, mdot, xfer_ratio, &
            L_div_Ledd, full_off, full_on, max_boost, super_eddington_boost, &
            hot_wind, cool_wind, H_env_mass, H_He_env_mass, He_layer_mass
         character (len=strlen) :: message, cool_wind_scheme, hot_wind_scheme, scheme
         logical :: is_infalling, using_wind_scheme_mdot, use_other
         real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula

         logical, parameter :: dbg = .false.

         include 'formats'

         if (dbg) write(*,1) 'enter set_mdot mass_change', s% mass_change

         ierr = 0

         xfer_ratio = 1d0

         L1 = L_phot
         M1 = M_phot
         T1 = T_phot
         R1 = sqrt(L1/(pi*crad*clight*T1*T1*T1*T1)) ! assume L1 and T1 for photosphere

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         nz = s% nz
         wind_mdot = 0
         using_wind_scheme_mdot = .false.

         call eval_super_eddington_wind(s, L1, M1, R1, ierr)
         if (ierr /= 0) then
            if (dbg .or. s% report_ierr) write(*, *) 'set_mdot: eval_super_eddington_wind ierr'
            return
         end if

         if (s% super_eddington_wind_mdot > wind_mdot) then
            wind_mdot = s% super_eddington_wind_mdot
            if (dbg) write(*,1) 'super eddington wind lg_Mdot', &
               log10(s% super_eddington_wind_mdot/(Msun/secyer))
         end if

         mdot = eval_rlo_wind(s, L1/Lsun, R1/Rsun, T1, xfer_ratio, ierr) ! Msun/year
         mdot = mdot*Msun/secyer
         if (ierr /= 0) then
            if (dbg .or. s% report_ierr) write(*, *) 'set_mdot: eval_rlo_wind ierr'
            return
         end if
         s% doing_rlo_wind = (mdot /= 0)
         if (dbg) write(*,*) 's% doing_rlo_wind', s% doing_rlo_wind, mdot, wind_mdot

         if (s% doing_rlo_wind .and. mdot > wind_mdot) then
            wind_mdot = mdot
            if (dbg) write(*,1) 's% doing_rlo_wind mdot', wind_mdot
         end if

         if (h1 > 0) then
            center_h1 = s% xa(h1,nz)
            surface_h1 = s% xa(h1,1)
         else
            center_h1 = 0
            surface_h1 = 0
         end if
         if (he4 > 0) then
            center_he4 = s% xa(he4,nz)
            surface_he4 = s% xa(he4,1)
         else
            center_he4 = 0
            surface_he4 = 0
         end if

         if ((s% mass_change > 0 .and. wind_mdot == 0) .or. &
             (s% mass_change < 0 .and. -s% mass_change*Msun/secyer > wind_mdot)) then
            if (dbg) write(*,*) 'mass_change mdot', s% mass_change
            wind_mdot = -s% mass_change*Msun/secyer
            if (s% mass_change > 0 .and. xfer_ratio < 1d0) then
               write(*,1) 'almost full roche lobe: reduce mdot by xfer_ratio', xfer_ratio
               wind_mdot = wind_mdot*xfer_ratio
            end if
         end if

         !separate calls for hot and cool winds
         if(s% hot_wind_full_on_T < s% cool_wind_full_on_T)then
            ierr = -1
            write(*,*) ' *** set_mdot error: hot_wind_full_on_T < cool_wind_full_on_T '
            return
         endif

         if(T1 >= s% cool_wind_full_on_T)then !do hot_wind calculation
            call eval_wind_for_scheme(s% hot_wind_scheme,hot_wind)
         else
            hot_wind = 0d0
         endif

         if(T1 <= s% hot_wind_full_on_T)then
            if (center_h1 < 0.01d0 .and. center_he4 < s% RGB_to_AGB_wind_switch) then
               scheme = s% cool_wind_AGB_scheme
               if (dbg) &
                  write(*,1) 'using cool_wind_AGB_scheme: "' // trim(scheme) // '"', &
                    center_h1, center_he4, s% RGB_to_AGB_wind_switch
            else
               scheme = s% cool_wind_RGB_scheme
               if (dbg) &
                  write(*,*) 'using cool_wind_RGB_scheme: "' // trim(scheme) // '"'
            end if
            call eval_wind_for_scheme(scheme,cool_wind)
         else
            cool_wind = 0d0
         endif

         !now combine the contributions of hot and cool winds
         if(T1 >= s% hot_wind_full_on_T)then
            wind = hot_wind
         else if(T1 <= s% cool_wind_full_on_T)then
            wind = cool_wind
         else if(s% hot_wind_full_on_T == s% cool_wind_full_on_T)then
            wind = 0.5d0*(hot_wind + cool_wind)
         else ! blend
            divisor = s% hot_wind_full_on_T - s% cool_wind_full_on_T
            beta = min( (s% hot_wind_full_on_T - T1) / divisor, 1d0)
            alfa = 1d0 - beta
            wind = alfa*hot_wind + beta*cool_wind
         endif

         if (wind*Msun/secyer > abs(wind_mdot)) then
            using_wind_scheme_mdot = .true.
            if (dbg) write(*,1) 'use wind scheme mdot', wind*Msun/secyer, wind_mdot
            wind_mdot = wind*Msun/secyer
         end if

         if (dbg) write(*,1) 'wind_mdot 1', wind_mdot
         
         if (using_wind_scheme_mdot) then
            if (s% no_wind_if_no_rotation .and. .not. s% rotation_flag) then
               s% mstar_dot = 0
               if (s% trace_dt_control_mass_change) &
                  write(*,1) 'no_wind_if_no_rotation'
               return
            end if
            if (s% dt > 0 .and. s% dt < s% mass_change_full_on_dt) then
               if (s% dt <= s% mass_change_full_off_dt) then
                  s% mstar_dot = 0
                  if (s% trace_dt_control_mass_change .or. dbg) &
                     write(*,1) 'no wind: dt <= mass_change_full_off_dt'
                  return
               end if
               alfa = (s% dt - s% mass_change_full_off_dt)/ &
                        (s% mass_change_full_on_dt - s% mass_change_full_off_dt)
               if (s% trace_dt_control_mass_change .or. dbg) &
                  write(*,1) 'reduce wind: dt <= mass_change_full_on_dt', alfa
               wind_mdot = wind_mdot*alfa
            end if
         end if

         if (wind_mdot >= 0 .and. s% super_eddington_scaling_factor <= 0) then
            ! check for super eddington boost to wind
            L_div_Ledd = L1 / s% prev_Ledd
            full_off = s% wind_boost_full_off_L_div_Ledd
            if (L_div_Ledd > full_off) then
               full_on = s% wind_boost_full_on_L_div_Ledd
               max_boost = s% super_eddington_wind_max_boost
               if (L_div_Ledd >= full_on) then
                  super_eddington_boost = max_boost
               else
                  super_eddington_boost = &
                     1 + (max_boost-1)*(L_div_Ledd - full_off)/(full_on - full_off)
               end if
               wind_mdot = wind_mdot*super_eddington_boost
               if (s% trace_super_eddington_wind_boost .or. dbg) then
                  write(*,1) 'super eddington wind boost factor, L_div_Ledd', &
                     super_eddington_boost, L_div_Ledd
                  write(*,*)
               end if
            end if
         end if

         if (dbg) write(*,1) 'wind_mdot 2', wind_mdot

         if (wind_mdot >= 0 .and. s% min_wind > 0 .and. &
               wind_mdot < s% min_wind*Msun/secyer) then
            if (dbg) write(*,1) 'use s% min_wind', s% min_wind
            wind_mdot = s% min_wind*Msun/secyer
         end if

         if (dbg) write(*,1) 'wind_mdot 3', wind_mdot

         if (wind_mdot >= 0 .and. s% max_wind > 0 .and. &
               wind_mdot > s% max_wind*Msun/secyer) then
            if (dbg) write(*,1) 'use s% max_wind', s% max_wind
            wind_mdot = s% max_wind*Msun/secyer
         end if
         if (dbg) write(*,1) 'wind_mdot 4', wind_mdot

         if (wind_mdot >= 0) then
            if (s% starting_T_center > s% max_T_center_for_any_mass_loss) then
               if (dbg) write(*,1) 'starting_T_center > max_T_center_for_any_mass_loss', &
                        s% starting_T_center, s% max_T_center_for_any_mass_loss
               wind_mdot = 0
            else if (s% starting_T_center > s% max_T_center_for_full_mass_loss) then
               if (dbg) write(*,1) 'starting_T_center > max_T_center_for_full_mass_loss', &
                        s% starting_T_center, s% max_T_center_for_full_mass_loss
               wind_mdot = wind_mdot* &
                  (s% max_T_center_for_any_mass_loss - s% starting_T_center)/ &
                  (s% max_T_center_for_any_mass_loss - &
                     s% max_T_center_for_full_mass_loss)
            end if
         end if
         if (dbg) write(*,1) 'wind_mdot 5', wind_mdot

         if (wind_mdot >= 0) then
             H_env_mass = s% star_mass - s% he_core_mass
             H_He_env_mass = s% star_mass - s% co_core_mass
             He_layer_mass = s% he_core_mass - s% co_core_mass
             if (s% wind_H_envelope_limit > 0 .and. &
                   H_env_mass < s% wind_H_envelope_limit) then
                wind_mdot = 0
             else if (s% wind_H_He_envelope_limit > 0 .and. &
                   H_He_env_mass < s% wind_H_He_envelope_limit) then
                wind_mdot = 0
             else if (s% wind_He_layer_limit > 0 .and. &
                   He_layer_mass < s% wind_He_layer_limit) then
                wind_mdot = 0
             end if
         end if

         s% mstar_dot = -wind_mdot
         if (dbg) write(*,1) 'mstar_dot', s% mstar_dot
         if (s% mstar_dot < 0 .and. &
               (s% min_wind > 0 .or. &
                  using_wind_scheme_mdot .or. &
                     s% v_div_v_crit_avg_surf > 0.8d0)) then
            call rotation_enhancement(ierr)
            if (is_bad(s% rotational_mdot_boost)) then
               write(*,2) 'is_bad(s% rotational_mdot_boost)', s% model_number
               if (s% stop_for_bad_nums) stop 'winds: rotation_enhancement'
            end if
            if (ierr /= 0) then
               if (dbg .or. s% report_ierr) write(*, *) 'set_mdot: rotation_enhancement ierr'
               return
            end if
         end if

         s% explicit_mstar_dot = s% mstar_dot

         if (dbg) then
            write(*,1) 'final star_mdot', s% mstar_dot/(Msun/secyer)
            write(*,1) 'final lg abs s% mstar_dot/(Msun/secyer)', safe_log10(abs(s% mstar_dot/(Msun/secyer)))
            write(*,*)
         end if

         contains

           subroutine eval_wind_for_scheme(scheme,wind)
             character(len=strlen) :: scheme
             real(dp), intent(out) :: wind
             include 'formats'

             use_other = (s% use_other_wind .or. scheme == 'other')
             if ((.not. use_other) .and. len_trim(scheme) == 0) then
                wind = 0
             else
                wind = 4d-13*(L1*R1/M1)/(Lsun*Rsun/Msun) ! in Msun/year
                if (dbg) write(*,1) 'wind', wind
                if (wind <= 0 .or. is_bad_num(wind)) then
                   ierr = -1
                   write(*,*) 'bad value for wind :', wind,L1,R1,M1
                   if (dbg) stop 'debug: bad value for wind'
                   if (s% stop_for_bad_nums) stop 'winds'
                   return
                end if
                X = surface_h1
                Y = surface_he4
                Z = 1 - (X + Y)

                if (use_other) then
                   if (dbg) write(*,*) 'call other_wind'
                   call s% other_wind(s% id, L1, M1, R1, T1, X, Y, Z, wind, ierr)
                   if (ierr /= 0) return
                else if (scheme == 'Dutch') then
                   T_high = 11000
                   T_low = 10000
                   if (s% Dutch_scaling_factor == 0) then
                      wind = 0
                   else if (T1 <= T_low) then
                      call eval_lowT_Dutch(wind)
                   else if (T1 >= T_high) then
                      call eval_highT_Dutch(wind)
                   else ! transition
                      call eval_lowT_Dutch(w1)
                      call eval_highT_Dutch(w2)
                      alfa = (T1 - T_low)/(T_high - T_low)
                      wind = (1-alfa)*w1 + alfa*w2
                   end if
                   wind = s% Dutch_scaling_factor * wind
                else if (scheme == 'Reimers') then
                   wind = wind * s% Reimers_scaling_factor
                   if (dbg) then
                      write(*,1) 's% Reimers_scaling_factor', s% Reimers_scaling_factor
                      write(*,1) 'Reimers_wind', wind
                      write(*,1) 'L1/Lsun', L1/Lsun
                      write(*,1) 'R1/Rsun', R1/Rsun
                      write(*,1) 'M1/Msun', M1/Msun
                      write(*,1) 'Reimers_scaling_factorReimers_scaling_factor', s% Reimers_scaling_factor
                      write(*,1) 'wind', wind
                      write(*,1) 'log10 wind', log10(wind)
                      write(*,*)
                      stop 'debug: winds'
                   end if
                else if (scheme == 'Vink') then
                   call eval_Vink_wind(wind)
                   wind = wind * s% Vink_scaling_factor
                   if (dbg) write(*,1) 'Vink_wind', wind
                   if (dbg) write(*,1) 'Grafener_wind', wind
                else if (scheme == 'Blocker') then
                   call eval_blocker_wind(wind)
                   if (dbg) write(*,1) 'Blocker_wind', wind
                else if (scheme == 'de Jager') then
                   call eval_de_Jager_wind(wind)
                   wind = s% de_Jager_scaling_factor * wind
                   if (dbg) write(*,1) 'de_Jager_wind', wind
                else if (scheme == 'van Loon') then
                   call eval_van_Loon_wind(wind)
                   wind = s% van_Loon_scaling_factor * wind
                   if (dbg) write(*,1) 'van_Loon_wind', wind
                else if (scheme == 'Nieuwenhuijzen') then
                   call eval_Nieuwenhuijzen_wind(wind)
                   wind = s% Nieuwenhuijzen_scaling_factor * wind
                   if (dbg) write(*,1) 'Nieuwenhuijzen_wind', wind
                else
                   ierr = -1
                   write(*,*) 'unknown name for wind scheme : ' // trim(scheme)
                   if (dbg) stop 'debug: bad value for wind scheme'
                   return
                end if
             end if

           end subroutine eval_wind_for_scheme


         subroutine rotation_enhancement(ierr)
            use star_utils, only: eval_kh_timescale
            integer, intent(out) :: ierr
            ! as in Heger, Langer, and Woosley, 2000, ApJ, 528:368-396.  section 2.6
            ! Mdot = Mdot_no_rotation/(1 - Osurf/Osurf_crit)^mdot_omega_power
            ! where Osurf = angular velocity at surface
            !       Osurf_crit^2 = (1 - Gamma_edd)*G*M/R_equatorial^3
            !       Gamma_edd = kappa*L/(4 pi c G M), Eddington factor
            real(dp) :: enhancement, wind_mdot, &
               kh_timescale, mdot_lim, wind_mdot_prev, dmsfac, dmskhf, &
               wind_mdot_lim, v_div_v_crit_full_on, v_div_v_crit_full_off

            include 'formats'

            ierr = 0

            if (.not. s% rotation_flag) return
            if (s% mdot_omega_power <= 0) return
            if (s% mstar_dot >= 0) return

            wind_mdot = -s% mstar_dot

            kh_timescale = eval_kh_timescale(s% cgrav(1), M1, R1, L1)
            dmskhf = s% rotational_mdot_kh_fac
            dmsfac = s% rotational_mdot_boost_fac
            wind_mdot_lim = min(dmskhf*M1/kh_timescale, wind_mdot*dmsfac)

            enhancement = pow(max(1d-22, 1d0 - s% v_div_v_crit_avg_surf), -s% mdot_omega_power)
            if (s% max_rotational_mdot_boost > 0 .and. &
                  enhancement > s% max_rotational_mdot_boost) then
               enhancement = s% max_rotational_mdot_boost
            end if

            if (enhancement > s% lim_trace_rotational_mdot_boost) then
               if (dbg) write(*,1) &
                  'mdot rotation enhancement factor for mdot, v_div_v_crit_avg_surf', &
                  enhancement, s% v_div_v_crit_avg_surf
            end if

            if (wind_mdot*enhancement < wind_mdot_lim) then
               wind_mdot = wind_mdot*enhancement
               if (dbg) write(*,2) 'wind_mdot = wind_mdot*enhancement', &
                  s% model_number, enhancement, &
                  log10(wind_mdot/(Msun/secyer)), log10(wind_mdot_lim/(Msun/secyer))
            else
               enhancement = wind_mdot_lim/wind_mdot
               wind_mdot = wind_mdot_lim
               if (dbg) write(*,2) 'wind_mdot = wind_mdot_lim', &
                  s% model_number, log10(wind_mdot/(Msun/secyer))
            end if

            wind_mdot_prev = -s% mstar_dot_old
            if (wind_mdot > 0 .and. wind_mdot_prev > 0) &
               wind_mdot = min(wind_mdot, s% max_mdot_jump_for_rotation*wind_mdot_prev)

            ! recheck max_wind
            if (wind_mdot >= 0 .and. s% max_wind > 0 .and. &
                  wind_mdot > s% max_wind*Msun/secyer) then
               if (dbg) write(*,1) 'use s% max_wind', s% max_wind
               wind_mdot = s% max_wind*Msun/secyer
            end if

            s% mstar_dot = -wind_mdot
            s% rotational_mdot_boost = enhancement

         end subroutine rotation_enhancement


         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
               Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
               end if
            end if

            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z/Zsolar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if

            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30d0) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z/Zsolar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if

            w = alfa*w1 + (1 - alfa)*w2

            if (dbg) write(*,*) 'vink wind', w

         end subroutine eval_Vink_wind


         subroutine eval_blocker_wind(w)
            real(dp), intent(inout) :: w
            w = w * s% Blocker_scaling_factor * &
               4.83d-9 * pow(M1/Msun,-2.1d0) * pow(L1/Lsun,2.7d0)
            if (dbg) write(*,*) 'blocker wind', w
         end subroutine eval_blocker_wind


         subroutine eval_highT_Dutch(w)
            real(dp), intent(out) :: w
            include 'formats'
            if (surface_h1 < 0.4d0) then ! helium rich Wolf-Rayet star: Nugis & Lamers
               w = 1d-11 * pow(L1/Lsun,1.29d0) * pow(Y,1.7d0) * sqrt(Z)
               if (dbg) write(*,1) 'Dutch_wind = Nugis & Lamers', log10(wind)
            else
               call eval_Vink_wind(w)
            end if
         end subroutine eval_highT_Dutch


         subroutine eval_lowT_Dutch(w)
            real(dp), intent(out) :: w
            include 'formats'
            if (s% Dutch_wind_lowT_scheme == 'de Jager') then
               call eval_de_Jager_wind(w)
               if (dbg) write(*,1) 'Dutch_wind = de Jager', safe_log10(wind), T1, T_low, T_high
            else if (s% Dutch_wind_lowT_scheme == 'van Loon') then
               call eval_van_Loon_wind(w)
               if (dbg) write(*,1) 'Dutch_wind = van Loon', safe_log10(wind), T1, T_low, T_high
            else if (s% Dutch_wind_lowT_scheme == 'Nieuwenhuijzen') then
               call eval_Nieuwenhuijzen_wind(w)
               if (dbg) write(*,1) 'Dutch_wind = Nieuwenhuijzen', safe_log10(wind), T1, T_low, T_high
            else
               write(*,*) 'unknown value for Dutch_wind_lowT_scheme ' // &
                  trim(s% Dutch_wind_lowT_scheme)
               w = 0
            end if
         end subroutine eval_lowT_Dutch


         subroutine eval_de_Jager_wind(w)
            ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
            w = exp10(log10w)
            if (dbg) then
               write(*,1) 'de_Jager log10 wind', log10w
            end if
         end subroutine eval_de_Jager_wind


         subroutine eval_van_Loon_wind(w)
            ! van Loon et al. 2005, A&A, 438, 273
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -5.65d0 + 1.05d0*log10(L1/(1d4*Lsun)) - 6.3d0*log10(T1/35d2)
            w = exp10(log10w)
         end subroutine eval_van_Loon_wind


         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 + &
                     1.24d0*log10(L1/Lsun) + &
                     0.16d0*log10(M1/Msun) + &
                     0.81d0*log10(R1/Rsun)
            w = exp10(log10w)
            if (dbg) then
               write(*,1) 'Nieuwenhuijzen log10 wind', log10w
            end if
         end subroutine eval_Nieuwenhuijzen_wind

      end subroutine do_set_mdot


      subroutine eval_super_eddington_wind(s, L, M, R, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: L, M, R
         integer, intent(out) :: ierr

         real(dp) :: Ledd, Leff, vesc2
         include 'formats'

         ierr = 0
         s% super_eddington_wind_mdot = 0
         if (s% super_eddington_scaling_factor <= 0) return

         Ledd = s% prev_Ledd
         Leff = L/s% super_eddington_wind_Ledd_factor
         if (Leff <= Ledd) return
         vesc2 = s% cgrav(1)*M/R  ! GM/R vs. 2GM/R ?
         s% super_eddington_wind_mdot = s% super_eddington_scaling_factor*(Leff - Ledd)/vesc2
         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,'(a60,i12,1p2e12.4)') 'super eddington wind: lg_Mdot, L/Ledd', &
               s% model_number, log10(s% super_eddington_wind_mdot/(Msun/secyer)), L/Ledd

      end subroutine eval_super_eddington_wind


      real(dp) function eval_rlo_wind(s, L_surf, R, Teff, xfer_ratio, ierr) ! value in Msun/year
         type (star_info), pointer :: s
         real(dp), intent(in) :: L_surf, R, Teff ! Lsun, Rsun, K
         real(dp), intent(inout) :: xfer_ratio
         integer, intent(out) :: ierr
         real(dp) :: roche_lobe_radius ! Rsun
         real(dp) :: ratio, rho, p, grav, hp, scale_height, h, rho_exponent, rho_rl, rho_rl0, mdot
         include 'formats'
         ierr = 0
         eval_rlo_wind = 0
         if (s% rlo_scaling_factor <= 0) return
         if (L_surf < s% rlo_wind_min_L) return
         if (Teff > s% rlo_wind_max_Teff) return
         scale_height = s% rlo_wind_scale_height
         if (scale_height <= 0) then
            scale_height = s% Peos(1) / (s% cgrav(1)*s% m(1)*s% rho(1) / (s% r(1)**2)) / Rsun
         end if 
         roche_lobe_radius = s% rlo_wind_roche_lobe_radius
         ratio = R/roche_lobe_radius
         !write(*,2) 'R/roche_lobe_radius', s% model_number, ratio
         if (ratio < 1) then
            ! check for reduction in transfer ratio for almost full Roche lobe
            if (ratio < s% roche_lobe_xfer_full_on) return
            if (ratio > s% roche_lobe_xfer_full_off) then
               xfer_ratio = 0
               return
            end if
            xfer_ratio = (s% roche_lobe_xfer_full_off - ratio) / &
               (s% roche_lobe_xfer_full_off - s% roche_lobe_xfer_full_on)
            xfer_ratio = 0.5d0*(1 - cospi(xfer_ratio))
            return
         end if
         mdot = s% rlo_wind_base_mdot* &
            exp(min(6*ln10,(R - roche_lobe_radius)/scale_height))
         eval_rlo_wind = s% rlo_scaling_factor*mdot ! Msun/year
         
         !write(*,1) 's% rlo_wind_base_mdot', s% rlo_wind_base_mdot
         !write(*,1) 'R - roche_lobe_radius', R - roche_lobe_radius, R, roche_lobe_radius
         !write(*,1) 'scale_height', scale_height
         !write(*,1) 's% rlo_scaling_factor', s% rlo_scaling_factor
         !write(*,1) 'eval_rlo_wind, log eval_rlo_wind, R/R_L', log10(eval_rlo_wind), R/roche_lobe_radius

      end function eval_rlo_wind


      end module winds
