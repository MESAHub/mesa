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

      module report

      use star_private_def
      use chem_def
      use utils_lib
      use star_utils
      use num_lib, only: find0

      use const_def, only: avo, kerg, pi, clight, crad, Rsun, Lsun, Msun, &
         secyer, ln10, mev_amu, ev2erg, two_thirds

      implicit none

      contains
      
      
      subroutine set_phot_info(s)
         use atm_lib, only: atm_black_body_T
         type (star_info), pointer :: s
         real(dp) :: luminosity
         include 'formats'
         call get_phot_info(s, &
            s% photosphere_r, s% photosphere_m, s% photosphere_v, &
            s% photosphere_L, s% photosphere_T, s% photosphere_csound, &
            s% photosphere_opacity, s% photosphere_logg, &
            s% photosphere_column_density, s% photosphere_cell_k)
         s% photosphere_black_body_T = &
            atm_black_body_T(s% photosphere_L, s% photosphere_r)
         s% photosphere_r = s% photosphere_r/Rsun
         s% photosphere_m = s% photosphere_m/Msun
         s% photosphere_L = s% photosphere_L/Lsun
         s% L_phot = s% photosphere_L
         luminosity = s% L(1)
         if (is_bad(luminosity)) then
            write(*,2) 's% L(1)', s% model_number, s% L(1)
            write(*,2) 's% xh(s% i_lum,1)', s% model_number, s% xh(s% i_lum,1)
            stop 'set_phot_info'
            luminosity = 0d0
         end if
         if (s% Teff < 0 .or. is_bad(s% Teff)) s% Teff = s% photosphere_black_body_T
         s% L_surf = luminosity/Lsun
         s% log_surface_luminosity = log10(max(1d-99,luminosity/Lsun))
            ! log10(stellar luminosity in solar units)
         s% log_L_surf = s% log_surface_luminosity
         if (is_bad(s% L_surf)) then
            write(*,2) 's% L_surf', s% model_number, s% L_surf
            stop 'set_phot_info'
         end if
      end subroutine set_phot_info


      subroutine set_min_gamma1(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: vesc
         include 'formats'
         s% min_gamma1 = 1d99
         do k = s% nz, 1, -1
            if (s% q(k) > s% gamma1_limit_max_q) exit
            vesc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (s% u_flag) then
               if (s% u(k) > vesc*s% gamma1_limit_max_v_div_vesc) exit
            else if (s% v_flag) then
               if (s% v(k) > vesc*s% gamma1_limit_max_v_div_vesc) exit
            end if
            if (s% gamma1(k) < s% min_gamma1) s% min_gamma1 = s% gamma1(k)
         end do
      end subroutine set_min_gamma1


      subroutine set_power_info(s)
         use chem_def, only: category_name
         type (star_info), pointer :: s
         integer :: j, k, nz
         real(dp) :: eps_nuc
         include 'formats'
         nz = s% nz

         do j=1,num_categories
            s% L_by_category(j) = &
               dot_product(s% dm(1:nz), s% eps_nuc_categories(j,1:nz))/Lsun
            s% center_eps_burn(j) = center_value_eps_burn(j)
            if (is_bad(s% L_by_category(j))) then
               do k=1,nz
                  if (is_bad(s% eps_nuc_categories(j,k))) then
                     write(*,2) trim(category_name(j)) // ' eps_nuc logT', k, s% eps_nuc_categories(j,k), s% lnT(k)/ln10
                     if (s% stop_for_bad_nums) stop 'set_power_info'
                  end if
               end do
            end if
         end do  
                
         if (s% eps_nuc_factor == 0d0) then
            s% power_nuc_burn = 0d0
            s% power_nuc_neutrinos = 0d0
            s% power_nonnuc_neutrinos = 0d0
            s% power_neutrinos = 0d0
            s% power_h_burn = 0d0
            s% power_he_burn = 0d0
            s% power_z_burn = 0d0
            s% power_photo = 0d0
         else            
            ! better if set power_nuc_burn using eps_nuc instead of categories
            ! categories can be subject to numerical jitters at very high temperatures
            s% power_nuc_burn = 0d0
            do k=1,nz
               if (s% op_split_burn .and. s% T_start(k) >= s% op_split_burn_min_T) then
                  eps_nuc = s% burn_avg_epsnuc(k)
               else
                  eps_nuc = s% eps_nuc(k)
               end if
               s% power_nuc_burn = s% power_nuc_burn + eps_nuc*s% dm(k)
            end do
            s% power_nuc_burn = s% power_nuc_burn/Lsun            
            s% power_nuc_neutrinos = dot_product(s% dm(1:nz),s% eps_nuc_neu_total(1:nz))/Lsun
            s% power_h_burn = s% L_by_category(ipp) + s% L_by_category(icno)
            s% power_he_burn = s% L_by_category(i3alf)
            s% power_z_burn = s% power_nuc_burn - (s% power_h_burn + s% power_he_burn)
            s% power_photo = s% L_by_category(iphoto)
         end if
         
         if (s% non_nuc_neu_factor == 0d0) then
            s% power_nonnuc_neutrinos = 0d0
         else
            s% power_nonnuc_neutrinos = &
                 dot_product(s% dm(1:nz),s% non_nuc_neu(1:nz))/Lsun
         end if
         s% power_neutrinos = s% power_nuc_neutrinos + s% power_nonnuc_neutrinos
         s% L_nuc_burn_total = s% power_nuc_burn
         
         contains

         real(dp) function center_value_eps_burn(j)
            integer, intent(in) :: j
            real(dp) :: sum_x, sum_dq, dx, dq
            integer :: k
            sum_x = 0
            sum_dq = 0
            do k = s% nz, 1, -1
               dq = s% dq(k)
               dx = s% eps_nuc_categories(j,k)*dq
               if (sum_dq+dq >= s% center_avg_value_dq) then
                  sum_x = sum_x + dx*(s% center_avg_value_dq - sum_dq)/dq
                  sum_dq = s% center_avg_value_dq
                  exit
               end if
               sum_x = sum_x + dx
               sum_dq = sum_dq + dq
            end do
            center_value_eps_burn = sum_x/sum_dq
         end function center_value_eps_burn

      end subroutine set_power_info


      subroutine do_report(s, ierr)
         use rates_def, only: &
            i_rate, i_rate_dRho, i_rate_dT, std_reaction_Qs, std_reaction_neuQs
         use star_utils, only: get_phot_info
         use hydro_rotation, only: set_surf_avg_rotation_info
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: k, i, j, ic, nz, kcore, &
            h1, h2, he3, he4, c12, n14, o16, ne20, si28, co56, ni56, k_min
         real(dp) :: w1, radius, dr, dm, hpc, cur_m, cur_r, prev_r, tau_conv, &
            twoGmrc2, cur_h, prev_h, cur_he, non_fe_core_mass, nu_for_delta_Pg, &
            prev_he, cur_c, prev_c, v, mstar, pdg, pdg_prev, luminosity, &
            prev_m, cell_mass, wf, conv_time, mv, bminv, uminb, eps_nuc_sum, eps_cat_sum
         logical, parameter :: new_only = .false.
         integer, pointer :: net_iso(:)

         include 'formats'

         ierr = 0         
         nz = s% nz
         net_iso => s% net_iso

         h1 = net_iso(ih1)
         h2 = net_iso(ih2)
         he3 = net_iso(ihe3)
         he4 = net_iso(ihe4)
         c12 = net_iso(ic12)
         n14 = net_iso(in14)
         o16 = net_iso(io16)
         ne20 = net_iso(ine20)
         si28 = net_iso(isi28)
         co56 = net_iso(ico56)
         ni56 = net_iso(ini56)
         
         radius = s% r(1)  !  radius in cm
         s% log_surface_radius = log10(radius/Rsun)
            ! log10(stellar radius in solar units)
         s% log_center_density = center_value(s, s% lnd)/ln10
         s% log_max_temperature = maxval(s% lnT(1:nz))/ln10
         s% log_center_temperature = center_value(s, s% lnT)/ln10
         s% log_center_pressure = center_value(s, s% lnPeos)/ln10
         s% center_degeneracy = center_value(s, s% eta)

         s% center_eps_nuc = center_value(s, s% eps_nuc)

         s% d_center_eps_nuc_dlnT = center_value(s, s% d_epsnuc_dlnT)
         s% d_center_eps_nuc_dlnd = center_value(s, s% d_epsnuc_dlnd)
         s% center_non_nuc_neu = center_value(s, s% non_nuc_neu)

         s% center_gamma = center_value(s, s% gam)
         s% center_abar = center_value(s, s% abar)
         s% center_zbar = center_value(s, s% zbar)
         s% center_mu = center_value(s, s% mu)
         s% center_ye = center_value(s, s% ye)
         s% center_entropy = exp(center_value(s, s% lnS))*amu/kerg
         s% max_entropy = exp(maxval(s% lnS(1:nz)))*amu/kerg

         if (.not. s% rotation_flag) then
            s% omega(1:nz) = 0
            s% center_omega = 0
            s% center_omega_div_omega_crit = 0
         else
            s% center_omega = center_value(s, s% omega)
            s% center_omega_div_omega_crit = center_omega_div_omega_crit()
         end if

         s% log_surface_temperature = s% lnT(1)/ln10 ! log10(temperature at surface)
         s% log_surface_pressure = s% lnPeos(1)/ln10 ! log10(eos pressure at surface)
         s% log_surface_density = s% lnd(1)/ln10 ! log10(density at surface)
         
         luminosity = s% L(1)

         if (s% u_flag) then
            s% v_surf = s% u(1)
         else if (s% v_flag) then
            s% v_surf = s% v(1)
         else
            s% v_surf = 0d0
         end if

         call set_surf_avg_rotation_info(s)
         
         call set_min_gamma1(s)

         ! s% time is in seconds
         s% star_age = s% time/secyer
         s% time_years = s% time/secyer
         s% time_days = s% time/dble(60*60*24)
         if ( s% model_number <= 0 ) then
            s% star_age = 0d0
            s% time_days = 0d0
            s% time_years = 0d0
         end if
         
         ! s% dt is in seconds
         s% time_step = s% dt/secyer         ! timestep in years
         s% dt_years = s% dt/secyer
         s% dt_days = s% dt/dble(60*60*24)
         
         mstar = s% mstar
         s% star_mass = mstar/Msun             ! stellar mass in solar units

         s% kh_timescale = eval_kh_timescale(s% cgrav(1), mstar, radius, luminosity)/secyer
         ! kelvin-helmholtz timescale in years (about 1.6x10^7 for the sun)
         
         if (is_bad(s% kh_timescale)) then
            write(*,1) 's% kh_timescale', s% kh_timescale
            write(*,1) 's% cgrav(1)', s% cgrav(1)
            write(*,1) 'mstar', mstar
            write(*,1) 'radius', radius
            write(*,1) 'luminosity', luminosity
            stop 
         end if
         
         if (luminosity > 0d0) then
            s% nuc_timescale = 1d10*s% star_mass/(luminosity/Lsun)
         else
            s% nuc_timescale = 1d99
         end if
         ! nuclear timescale in years (e.g., about 10^10 years for sun)
         if (s% cgrav(1) <= 0d0) then
            s% dynamic_timescale = 1d99
         else
            s% dynamic_timescale = 2*pi*sqrt(radius*radius*radius/(s% cgrav(1)*mstar))
         end if

         if (h1 /= 0) then
            s% center_h1 = center_avg_x(s,h1)
            s% surface_h1 = surface_avg_x(s,h1)
         end if
         if (he3 /= 0) then
            s% center_he3 = center_avg_x(s,he3)
            s% surface_he3 = surface_avg_x(s,he3)
         end if
         if (he4 /= 0) then
            s% center_he4 = center_avg_x(s,he4)
            s% surface_he4 = surface_avg_x(s,he4)
         end if
         if (c12 /= 0) then
            s% center_c12 = center_avg_x(s,c12)
            s% surface_c12 = surface_avg_x(s,c12)
         end if
         if (n14 /= 0) then
            s% center_n14 = center_avg_x(s,n14)
            s% surface_n14 = surface_avg_x(s,n14)
         end if
         if (o16 /= 0) then
            s% center_o16 = center_avg_x(s,o16)
            s% surface_o16 = surface_avg_x(s,o16)
         end if
         if (ne20 /= 0) then
            s% center_ne20 = center_avg_x(s,ne20)
            s% surface_ne20 = surface_avg_x(s,ne20)
         end if
         if (si28 /= 0) then
            s% center_si28 = center_avg_x(s,si28)
         end if

         ! FYI profile stuff
         do k=1,nz

            s% entropy(k) = exp(s% lnS(k))/(avo*kerg)
            if (is_bad(s% entropy(k))) then
               ierr = -1
               write(*,2) 'report: s% entropy(k)', k, s% entropy(k)
               if (s% stop_for_bad_nums) stop 'report'
               return
            end if

            if (k == nz) then
               dr = s% r(k) - s% R_center
            else
               dr = s% r(k) - s% r(k+1)
            end if

            s% dr_div_csound(k) = dr/s% csound(k)
            if (is_bad(s% dr_div_csound(k))) then
               ierr = -1
               write(*,2) 'report: s% dr_div_csound(k)', &
                  k, s% dr_div_csound(k), dr, s% csound(k)
               if (s% stop_for_bad_nums) stop 'report'
               return
            end if
            
            if (s% u_flag) then
               v = s% u(k)
            else if (s% v_flag) then
               v = s% v(k)
            else
               v = 0d0
            end if

            if (is_bad(v)) then
               ierr = -1
               write(*,2) 'report: v', k, v
               if (s% stop_for_bad_nums) stop 'report'
               return
            end if
            
            s% v_div_csound(k) = v/s% csound_face(k)
            if (is_bad(s% v_div_csound(k))) then
               ierr = -1
               write(*,2) 'report: s% v_div_csound(k)', k, s% v_div_csound(k), &
                  v, s% csound_face(k)
               if (s% stop_for_bad_nums) stop 'report'
               return
            end if

         end do
         
         call set_phot_info(s)

         if (s% photosphere_r*Rsun >= s% r(1)) then
            s% photosphere_acoustic_r = sum(s% dr_div_csound(1:nz)) + &
               (s% photosphere_r*Rsun - s% r(1))/s% csound(1)
         else
            do k=2,nz
               if (s% photosphere_r*Rsun > s% r(k)) then
                  s% photosphere_acoustic_r = sum(s% dr_div_csound(k:nz)) + &
                     (s% photosphere_r*Rsun - s% r(k))/s% csound(k-1)
                  exit
               end if
            end do
         end if

         if (.not. s% get_delta_nu_from_scaled_solar) then
            s% delta_nu = 1d6/(2*s% photosphere_acoustic_r) ! microHz
         else
            s% delta_nu = &
               s% delta_nu_sun*sqrt(s% star_mass)*pow3(s% Teff/s% Teff_sun) / &
                  pow(s% L_phot,0.75d0)
         end if
         
         call get_mass_info(s, s% dm, ierr)
         if (failed('get_mass_info')) return
         
         s% nu_max = s% nu_max_sun*s% star_mass/ &
            (pow2(s% photosphere_r)*sqrt(max(0d0,s% Teff)/s% Teff_sun))
         s% acoustic_cutoff = &
            0.25d6/pi*s% grav(1)*sqrt(s% gamma1(1)*s% rho(1)/s% Peos(1))
         nu_for_delta_Pg = s% nu_max
         if (s% delta_Pg_mode_freq > 0) nu_for_delta_Pg = s% delta_Pg_mode_freq
         call get_delta_Pg(s, nu_for_delta_Pg, s% delta_Pg)

         if (s% rsp_flag) return

         call get_mixing_regions(s, ierr)
         if (failed('get_mixing_regions')) return

         call set_mass_conv_core
         call set_mass_semiconv_core
         call find_conv_mx_regions
         call find_mx_regions
         call get_burn_zone_info(s, ierr)
         if (failed('get_burn_zone_info')) return

         s% fe_core_infall = 0
         s% non_fe_core_infall = 0
         s% non_fe_core_rebound = 0
         if (s% u_flag) then
            k_min = minloc(s% u(1:nz),dim=1)
            if (k_min > 0) then
               s% max_infall_speed_mass = s% m(k_min)/Msun
            else
               s% max_infall_speed_mass = s% m(k_min)/Msun
            end if
            if (s% fe_core_mass > 0) then
               do k = 1, nz
                  if (s% m(k) > Msun*s% fe_core_mass) cycle
                  if (-s% u(k) > s% fe_core_infall) &
                     s% fe_core_infall = -s% u(k)
               end do
            end if
            non_fe_core_mass = s% he_core_mass
            if (non_fe_core_mass > 0) then
               do k = 1, nz
                  if (s% m(k) > Msun*non_fe_core_mass) cycle
                  if (s% m(k) < Msun*s% fe_core_mass) exit
                  if (-s% u(k) > s% non_fe_core_infall) &
                     s% non_fe_core_infall = -s% u(k)
                  if (s% u(k) > s% non_fe_core_rebound) then
                     s% non_fe_core_rebound = s% u(k)
                     !write(*,2) 's% non_fe_core_rebound', k, s% non_fe_core_rebound, s% m(k)/Msun
                  end if
               end do
            end if
         else if (s% v_flag) then
            k_min = minloc(s% v(1:nz),dim=1)
            if (k_min > 0) then
               s% max_infall_speed_mass = s% m(k_min)/Msun
            else
               s% max_infall_speed_mass = s% m(k_min)/Msun
            end if
            if (s% fe_core_mass > 0) then
               do k = 1, nz
                  if (s% m(k) > Msun*s% fe_core_mass) cycle
                  if (-s% v(k) > s% fe_core_infall) &
                     s% fe_core_infall = -s% v(k)
               end do
            end if
            non_fe_core_mass = max(s% he_core_mass, s% co_core_mass)
            if (non_fe_core_mass > 0) then
               do k = 1, nz
                  if (s% m(k) > Msun*non_fe_core_mass) cycle
                  if (s% m(k) < Msun*s% fe_core_mass) exit
                  if (-s% v(k) > s% non_fe_core_infall) &
                     s% non_fe_core_infall = -s% v(k)
                  if (s% v(k) > s% non_fe_core_rebound) &
                     s% non_fe_core_rebound = s% v(k)
               end do
            end if
         end if

         contains


         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed


         subroutine set_mass_conv_core
            integer :: j, nz, k
            real(dp) :: dm_limit
            include 'formats'
            s% mass_conv_core = 0
            dm_limit = s% conv_core_gap_dq_limit*s% xmstar
            nz = s% nz
            do j = 1, s% n_conv_regions
               ! ignore possible small gap at center
               if (s% cz_bot_mass(j) <= s% m(nz) + dm_limit) then
                  s% mass_conv_core = s% cz_top_mass(j)/Msun
                  ! jump over small gaps
                  do k = j+1, s% n_conv_regions
                     if (s% cz_bot_mass(k) - s% cz_top_mass(k-1) >= dm_limit) exit
                     s% mass_conv_core = s% cz_top_mass(k)/Msun
                  end do 
                  exit
               end if
            end do
         end subroutine set_mass_conv_core


         subroutine set_mass_semiconv_core
            integer :: k, ktop, nz
            real(dp) :: qb
            include 'formats'
            s% mass_semiconv_core = s% mass_conv_core
            nz = s% nz
            if (s% mixing_type(nz) /= convective_mixing) return
            ktop = 1
            do k=nz-1,1,-1
               if (s% mixing_type(k) == convective_mixing) then
                  ktop = k
                  exit
               end if
            end do
            if (ktop == 1 .or. s% mixing_type(ktop) /= semiconvective_mixing) return
            do k=ktop-1,1,-1
               if (s% mixing_type(k) /= semiconvective_mixing) then
                  qb = s% q(k) - s% cz_bdy_dq(k)
                  s% mass_semiconv_core = (qb*s% xmstar + s% M_center)/Msun
                  return
               end if
            end do
            s% mass_semiconv_core = s% star_mass
         end subroutine set_mass_semiconv_core


         real(dp) function volume_at_q(q)
            use interp_1d_def
            use interp_1d_lib
            real(dp), intent(in) :: q
            real(dp) :: vp2, vp1, v00, vm1
            integer, parameter :: n_old = 4, n_new = 1, nwork = pm_work_size
            real(dp) :: qlo, x_old(n_old), v_old(n_old), x_new(n_new), v_new(n_new)
            integer :: k, nz, k00, ierr
            real(dp), target :: work_ary(n_old*nwork)
            real(dp), pointer :: work(:)
            work => work_ary

            include 'formats'
            nz = s% nz

            if (q == 0d0) then
               volume_at_q = 0
               return
            end if

            if (q <= s% q(nz)) then
               volume_at_q = (q/s% q(nz))*four_thirds_pi*s% r(nz)*s% r(nz)*s% r(nz)
               return
            end if
            k00 = 1
            do k=nz-1,2,-1
               if (s% q(k) >= q) then
                  k00 = k; exit
               end if
            end do
            if (k00 == 1) then
               volume_at_q = four_thirds_pi*q*s% r(1)*s% r(1)*s% r(1)
               write(*,1) 'volume_at_q', volume_at_q
               return
            end if

            x_old(1) = 0
            if (k00+1 == nz) then
               v_old(1) = s% R_center*s% R_center*s% R_center
               qlo = 0
            else
               v_old(1) = s% r(k00+2)*s% r(k00+2)*s% r(k00+2)
               qlo = s% q(k00+2)
            end if

            x_old(2) = s% dq(k00+1)
            v_old(2) = s% r(k00+1)*s% r(k00+1)*s% r(k00+1)

            x_old(3) = x_old(2) + s% dq(k00)
            v_old(3) = s% r(k00)*s% r(k00)*s% r(k00)

            x_old(4) = x_old(3) + s% dq(k00-1)
            v_old(4) = s% r(k00-1)*s% r(k00-1)*s% r(k00-1)

            ierr = 0
            x_new(1) = q-qlo
            call interpolate_vector( &
               n_old, x_old, n_new, x_new, v_old, v_new, interp_pm, nwork, work, &
               'report volume_at_q', ierr)
            if (ierr /= 0) then
               write(*,*) 'volume_at_q: failed in interpolate_vector'
               volume_at_q = (v_old(2) + v_old(3))/2
               return
            end if

            volume_at_q = four_thirds_pi*v_new(1)

         end function volume_at_q


         real(dp) function center_omega_div_omega_crit()
            real(dp) :: sum_x, sum_dq, dx, dq
            integer :: k
            center_omega_div_omega_crit = 0
            if (.not. s% rotation_flag) return
            sum_x = 0
            sum_dq = 0
            do k = s% nz, 1, -1
               dq = s% dq(k)
               dx = dq*s% omega(k)/omega_crit(s,k)
               if (sum_dq+dq >= s% center_avg_value_dq) then
                  sum_x = sum_x+ dx*(s% center_avg_value_dq - sum_dq)/dq
                  sum_dq = s% center_avg_value_dq
                  exit
               end if
               sum_x = sum_x + dx
               sum_dq = sum_dq + dq
            end do
            center_omega_div_omega_crit = min(1d0,sum_x/sum_dq)
         end function center_omega_div_omega_crit


         subroutine find_conv_mx_regions
            real(dp) :: conv_mx1_dq, conv_mx2_dq, mx_dq
            integer :: i, ktop, kbot, conv_mx1, conv_mx2

            include 'formats'

            s% largest_conv_mixing_region = 0

            s% conv_mx1_top = 0
            s% conv_mx1_bot = 0
            s% conv_mx2_top = 0
            s% conv_mx2_bot = 0

            s% conv_mx1_top_r = 0
            s% conv_mx1_bot_r = 0
            s% conv_mx2_top_r = 0
            s% conv_mx2_bot_r = 0

            if (s% num_mixing_regions == 0) return

            conv_mx1 = 0
            conv_mx2 = 0
            conv_mx1_dq = 0
            conv_mx2_dq = 0

            do i = 1, s% num_mixing_regions
               if (s% mixing_region_type(i) /= convective_mixing) cycle
               ! mx_dq = s% q(s% mixing_region_top(i)) - s% q(s% mixing_region_bottom(i))
               mx_dq = sum(s% dq(s% mixing_region_top(i):s% mixing_region_bottom(i)))
               if (mx_dq > conv_mx1_dq) then
                  conv_mx2_dq = conv_mx1_dq; conv_mx1_dq = mx_dq
                  conv_mx2 = conv_mx1; conv_mx1 = i
               else if (mx_dq > conv_mx2_dq) then
                  conv_mx2_dq = mx_dq
                  conv_mx2 = i
               end if
            end do

            if (conv_mx1 > 0) then
               s% largest_conv_mixing_region = conv_mx1
               ktop = s% mixing_region_top(conv_mx1)
               kbot = s% mixing_region_bottom(conv_mx1)
               s% conv_mx1_top = s% q(ktop) - s% cz_bdy_dq(ktop)
               s% conv_mx1_bot = s% q(kbot) - s% cz_bdy_dq(kbot)
               s% conv_mx1_top_r = s% r(ktop)/Rsun
               s% conv_mx1_bot_r = s% r(kbot)/Rsun
            end if

            if (conv_mx2 > 0) then
               ktop = s% mixing_region_top(conv_mx2)
               kbot = s% mixing_region_bottom(conv_mx2)
               s% conv_mx2_top = s% q(ktop) - s% cz_bdy_dq(ktop)
               s% conv_mx2_bot = s% q(kbot) - s% cz_bdy_dq(kbot)
               s% conv_mx2_top_r = s% r(ktop)/Rsun
               s% conv_mx2_bot_r = s% r(kbot)/Rsun
            end if

         end subroutine find_conv_mx_regions


         subroutine find_mx_regions
            real(dp) :: mx1_dq, mx2_dq, mx_dq
            integer :: i, ktop, kbot, mx1_top_region, mx1_bottom_region, &
               mx2_top_region, mx2_bottom_region, &
               current_top_region, current_bottom_region, &
               current_top_point, current_bottom_point

            include 'formats'

            s% mx1_top = 0
            s% mx1_bot = 0
            s% mx2_top = 0
            s% mx2_bot = 0

            s% mx1_top_r = 0
            s% mx1_bot_r = 0
            s% mx2_top_r = 0
            s% mx2_bot_r = 0

            if (s% num_mixing_regions == 0) return

            mx1_top_region = 0
            mx1_bottom_region = 0
            mx1_dq = 0

            mx2_top_region = 0
            mx2_bottom_region = 0
            mx2_dq = 0

            i = 1
            do
               if (i > s% num_mixing_regions) exit
               current_top_region = i
               current_top_point = s% mixing_region_top(current_top_region)
               do
                  current_bottom_region = i
                  current_bottom_point = s% mixing_region_bottom(current_bottom_region)
                  i = i+1
                  if (i > s% num_mixing_regions) exit
                  if (s% mixing_region_top(i) /= current_bottom_point+1) exit
               end do
               ! mx_dq = s% q(current_top_point) - s% q(current_bottom_point)
               mx_dq = sum(s% dq(current_top_point:current_bottom_point))
               if (mx_dq > mx1_dq) then
                  mx2_dq = mx1_dq; mx1_dq = mx_dq
                  mx2_top_region = mx1_top_region
                  mx1_top_region = current_top_region
                  mx2_bottom_region = mx1_bottom_region
                  mx1_bottom_region = current_bottom_region
               else if (mx_dq > mx2_dq) then
                  mx2_dq = mx_dq
                  mx2_top_region = current_top_region
                  mx2_bottom_region = current_bottom_region
               end if
            end do

            if (mx1_top_region > 0) then
               ktop = s% mixing_region_top(mx1_top_region)
               kbot = s% mixing_region_bottom(mx1_bottom_region)
               s% mx1_top = s% q(ktop) - s% cz_bdy_dq(ktop)
               s% mx1_bot = s% q(kbot) - s% cz_bdy_dq(kbot)
               s% mx1_top_r = s% r(ktop)/Rsun
               s% mx1_bot_r = s% r(kbot)/Rsun
            end if

            if (mx2_top_region > 0) then
               ktop = s% mixing_region_top(mx2_top_region)
               kbot = s% mixing_region_bottom(mx2_bottom_region)
               s% mx2_top = s% q(ktop) - s% cz_bdy_dq(ktop)
               s% mx2_bot = s% q(kbot) - s% cz_bdy_dq(kbot)
               s% mx2_top_r = s% r(ktop)/Rsun
               s% mx2_bot_r = s% r(kbot)/Rsun
            end if

         end subroutine find_mx_regions


      end subroutine do_report


      subroutine get_burn_zone_info(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: burn_min1, burn_min2
         integer :: i, i_start
         include 'formats'
         burn_min1 = s% burn_min1; burn_min2 = s% burn_min2
         ! up to 3 zones where eps_nuc > burn_min1 erg/g/s
         ! for each zone have 4 numbers: start1, start2, end2, end1
         ! start1 is mass of inner edge where first goes > burn_min1 (or -20 if none such)
         ! start2 is mass of inner edge where first zone reaches burn_min2 erg/g/sec (or -20 if none such)
         ! end2 is mass of outer edge where first zone drops back below burn_min2 erg/g/s
         ! end1 is mass of outer edge where first zone ends (i.e. eps_nuc < burn_min1)
         ! similar for second and third zones
         i_start = s% nz
         do i=1,3
            call find_epsnuc_zone(s, i_start, &
               s% burn_zone_mass(1,i), s% burn_zone_mass(2,i), &
               s% burn_zone_mass(3,i), s% burn_zone_mass(4,i), &
               s% burn_min1, s% burn_min2, ierr)
         end do
      end subroutine get_burn_zone_info


      subroutine find_epsnuc_zone( &
            s, i_start, bzm_1, bzm_2, bzm_3, bzm_4, burn_min1, burn_min2, ierr)
         use const_def, only:Msun
         type (star_info), pointer :: s
         integer, intent(inout) :: i_start
         real(dp), intent(out) :: bzm_1, bzm_2, bzm_3, bzm_4
         real(dp), intent(in) :: burn_min1, burn_min2
         integer, intent(out) :: ierr

         real(dp), parameter :: null_zone = -20
         integer :: i, burn_zone
         real(dp) :: prev_m, prev_x, cur_m, cur_x
         ierr = 0
         bzm_1 = null_zone; bzm_2 = null_zone; bzm_3 = null_zone; bzm_4 = null_zone
         burn_zone = 0 ! haven't entered the zone yet
         if (i_start .ne. s% nz) then
            i = i_start+1
            prev_m = s% m(i)
            prev_x = s% eps_nuc(i)
         else ! keep the compiler happy
            prev_m = 0
            prev_x = 0
         end if
         do i = i_start, 1, -1
            cur_m = s% m(i)
            cur_x = s% eps_nuc(i)
            select case (burn_zone)
               case (0)
                  if ( cur_x .gt. burn_min2 ) then
                     if ( i .eq. s% nz ) then ! use star center as start of zone
                        bzm_2 = 0d0
                     else ! interpolate to estimate where rate reached burn_min1
                        bzm_2 = find0(prev_m, prev_x-burn_min2, cur_m, cur_x-burn_min2)
                     end if
                     bzm_1 = bzm_2
                     burn_zone = 2
                  elseif ( cur_x .gt. burn_min1 ) then
                     if ( i .eq. s% nz ) then ! use star center as start of zone
                        bzm_1 = 0d0
                     else ! interpolate to estimate where rate reached burn_min1
                        bzm_1 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     end if
                     burn_zone = 1
                  end if
               case (1) ! in the initial eps > burn_min1 region
                  if ( cur_x .gt. burn_min2 ) then
                     bzm_2 = find0(prev_m, prev_x-burn_min2, cur_m, cur_x-burn_min2)
                     burn_zone = 2
                  else if ( cur_x .lt. burn_min1 ) then
                     bzm_4 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     i_start = i
                     return
                  end if
               case (2) ! in the initial eps > burn_min2 region
                  if ( cur_x .lt. burn_min1 ) then
                     bzm_4 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     bzm_3 = bzm_4
                     i_start = i
                     return
                  end if
                  if ( cur_x .lt. burn_min2 ) then
                     bzm_3 = find0(prev_m, prev_x-burn_min2, cur_m, cur_x-burn_min2)
                     burn_zone = 3
                  end if
               case (3) ! in the final eps > burn_min1 region
                  if ( cur_x .lt. burn_min1 ) then
                     bzm_4 = find0(prev_m, prev_x-burn_min1, cur_m, cur_x-burn_min1)
                     i_start = i
                     return
                  end if
               case default
                  ierr = -1
                  write(*,*) 'error in find_eps_nuc_zone'
                  return
            end select
            prev_m = cur_m; prev_x = cur_x
         end do
         i_start = 0
         select case (burn_zone)
            case (0)
               return
            case (1)
               bzm_4 = cur_m
            case (2)
               bzm_3 = cur_m
               bzm_4 = cur_m
            case (3)
               bzm_4 = cur_m
            case default
               ierr = -1
               write(*,*) 'error in find_eps_nuc_zone'
               return
         end select
      end subroutine find_epsnuc_zone


      subroutine get_mass_info(s, cell_masses, ierr)
         type (star_info), pointer :: s
         real(dp), pointer, intent(in) :: cell_masses(:)
         integer, intent(out) :: ierr

         integer :: k, nz, j, nzlo, nzhi, kbdy, nzlo_prev
         real(dp) :: cell_mass
         integer, pointer :: net_iso(:)

         include 'formats'

         ierr = 0
         nz = s% nz
         net_iso => s% net_iso

         s% star_mass_h1=0d0
         s% star_mass_he3=0d0
         s% star_mass_he4=0d0
         s% star_mass_c12 = 0d0
         s% star_mass_n14 = 0d0
         s% star_mass_o16 = 0d0
         s% star_mass_ne20 = 0d0

         do k = 1, nz
            cell_mass = cell_masses(k)
            do j=1, s% species
               if (s% chem_id(j) == ih1) then
                  s% star_mass_h1 = s% star_mass_h1 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ihe3) then
                  s% star_mass_he3 = s% star_mass_he3 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ihe4) then
                  s% star_mass_he4 = s% star_mass_he4 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ic12) then
                  s% star_mass_c12 = s% star_mass_c12 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == in14) then
                  s% star_mass_n14 = s% star_mass_n14 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == io16) then
                  s% star_mass_o16 = s% star_mass_o16 + cell_mass*s% xa(j, k)
               else if (s% chem_id(j) == ine20) then
                  s% star_mass_ne20 = s% star_mass_ne20 + cell_mass*s% xa(j, k)
               end if
            end do
         end do

         s% star_mass_h1 = s% star_mass_h1 / Msun
         s% star_mass_he3 = s% star_mass_he3 / Msun
         s% star_mass_he4 = s% star_mass_he4 / Msun
         s% star_mass_c12 = s% star_mass_c12 / Msun
         s% star_mass_n14 = s% star_mass_n14 / Msun
         s% star_mass_o16 = s% star_mass_o16 / Msun
         s% star_mass_ne20 = s% star_mass_ne20 / Msun
         
         call get_core_info(s)
         call get_shock_info(s)

      end subroutine get_mass_info


      subroutine get_info_at_surface( &
            s, bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, bdy_time, &
            bdy_omega, bdy_omega_div_omega_crit)
         type (star_info), pointer :: s
         real(dp), intent(out) :: &
            bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, &
            bdy_omega, bdy_omega_div_omega_crit, bdy_time

         real(dp) :: bdy_omega_crit

         bdy_time = 0
         bdy_m = s% star_mass
         bdy_r = s% r(1)/Rsun
         bdy_lgT = s% lnT(1)/ln10
         bdy_lgRho = s% lnd(1)/ln10
         bdy_L = s% L(1)/Lsun
         if (s% v_flag) then
            bdy_v = s% v(1)
         else if (s% u_flag) then
            bdy_v = s% u_face_ad(1)%val
         else
            bdy_v = 0d0
         end if
         if (s% rotation_flag) then
            bdy_omega = s% omega_avg_surf
            bdy_omega_crit = s% omega_crit_avg_surf
            bdy_omega_div_omega_crit = s% w_div_w_crit_avg_surf
         else
            bdy_omega = 0
            bdy_omega_div_omega_crit = 0
         end if

      end subroutine get_info_at_surface


      subroutine get_mach1_location_info( &
            s, dbg, k_shock, v, r, &
            mach1_mass, &
            mach1_q, &
            mach1_radius, &
            mach1_velocity, &
            mach1_csound, &
            mach1_lgT, &
            mach1_lgRho, &
            mach1_lgP, &
            mach1_gamma1, &
            mach1_entropy, &
            mach1_tau, &
            mach1_k)
         type (star_info), pointer :: s
         integer, intent(in) :: k_shock
         logical, intent(in) :: dbg
         real(dp), intent(in), pointer :: v(:)
         real(dp), intent(in) :: r
         real(dp), intent(out) :: &
            mach1_mass, &
            mach1_q, &
            mach1_radius, &
            mach1_velocity, &
            mach1_csound, &
            mach1_lgT, &
            mach1_lgRho, &
            mach1_lgP, &
            mach1_gamma1, &
            mach1_entropy, &
            mach1_tau
         integer, intent(out) :: mach1_k

         integer :: k
         real(dp) :: alfa, beta

         include 'formats'

         k = k_shock
         if (r < s% R_center .or. r > s% r(1) .or. &
               k < 1 .or. k > s% nz .or. &
               .not. associated(v) .or. &
               .not. (s% v_flag .or. s% u_flag)) then
            mach1_mass = 0
            mach1_q = 0
            mach1_radius = 0
            mach1_velocity = 0
            mach1_csound = 0
            mach1_lgT = 0
            mach1_lgRho = 0
            mach1_lgP = 0
            mach1_gamma1 = 0
            mach1_entropy = 0
            mach1_tau = 0
            mach1_k = 0
            return
         end if
         
         mach1_radius = r/Rsun
         mach1_k = k
         if (k < s% nz) then
            alfa = (r - s% r(k))/(s% r(k+1) - s% r(k))
            beta = 1d0 - alfa
            mach1_mass = (alfa*s% m(k+1) + beta*s% m(k))/Msun
            mach1_q = alfa*s% q(k+1) + beta*s% q(k)
            mach1_velocity = alfa*v(k+1) + beta*v(k)
         else
            mach1_mass = s% m(k)/Msun
            mach1_q = s% q(k)
            mach1_velocity = v(k)
         end if
         mach1_csound = s% csound(k)
         mach1_lgT = s% lnT(k)/ln10
         mach1_lgRho = s% lnd(k)/ln10
         mach1_lgP = s% lnPeos(k)/ln10
         mach1_gamma1 = s% gamma1(k)
         mach1_entropy = s% entropy(k)
         mach1_tau = s% tau(k)

      end subroutine get_mach1_location_info


      subroutine get_core_info(s)
         type (star_info), pointer :: s

         integer :: k, j, A_max, h1, he4, c12, o16, ne20, si28, species, nz
         integer, pointer :: net_iso(:)
         logical :: have_he, have_co, have_one, have_fe
         real(dp) :: sumx, min_x
         integer, parameter :: &
            A_max_fe = 47, &
            A_max_si = 28, &
            A_max_one = 20, &
            A_max_co = 16, &
            A_max_he = 4

         include 'formats'

         net_iso => s% net_iso
         species = s% species
         nz = s% nz
         h1 = net_iso(ih1)
         he4 = net_iso(ihe4)
         c12 = net_iso(ic12)
         o16 = net_iso(io16)
         ne20 = net_iso(ine20)
         si28 = net_iso(isi28)
         min_x = s% min_boundary_fraction

         call clear_core_info(s% neutron_rich_core_k, &
            s% neutron_rich_core_mass, s% neutron_rich_core_radius, s% neutron_rich_core_lgT, &
            s% neutron_rich_core_lgRho, s% neutron_rich_core_L, s% neutron_rich_core_v, &
            s% neutron_rich_core_omega, s% neutron_rich_core_omega_div_omega_crit)

         do k=1,nz
            if (s% Ye(k) <= s% neutron_rich_core_boundary_Ye_max) then
               call set_core_info(s, k, s% neutron_rich_core_k, &
                  s% neutron_rich_core_mass, s% neutron_rich_core_radius, s% neutron_rich_core_lgT, &
                  s% neutron_rich_core_lgRho, s% neutron_rich_core_L, s% neutron_rich_core_v, &
                  s% neutron_rich_core_omega, s% neutron_rich_core_omega_div_omega_crit)
               exit
            end if
         end do

         call clear_core_info(s% fe_core_k, &
            s% fe_core_mass, s% fe_core_radius, s% fe_core_lgT, &
            s% fe_core_lgRho, s% fe_core_L, s% fe_core_v, &
            s% fe_core_omega, s% fe_core_omega_div_omega_crit)
         call clear_core_info(s% co_core_k, &
            s% co_core_mass, s% co_core_radius, s% co_core_lgT, &
            s% co_core_lgRho, s% co_core_L, s% co_core_v, &
            s% co_core_omega, s% co_core_omega_div_omega_crit)
         call clear_core_info(s% one_core_k, &
            s% one_core_mass, s% one_core_radius, s% one_core_lgT, &
            s% one_core_lgRho, s% one_core_L, s% one_core_v, &
            s% one_core_omega, s% one_core_omega_div_omega_crit)
         call clear_core_info(s% he_core_k, &
            s% he_core_mass, s% he_core_radius, s% he_core_lgT, &
            s% he_core_lgRho, s% he_core_L, s% he_core_v, &
            s% he_core_omega, s% he_core_omega_div_omega_crit)

         have_he = .false.
         have_co = .false.
         have_one = .false.
         have_fe = .false.

         do k=1,nz

            j = maxloc(s% xa(1:species,k),dim=1)
            A_max = chem_isos% Z_plus_N(s% chem_id(j))

            if (.not. have_fe) then
               if (s% fe_core_boundary_si28_fraction < 0) then
                  if (A_max >= A_max_fe) then
                     call set_core_info(s, k, s% fe_core_k, &
                        s% fe_core_mass, s% fe_core_radius, s% fe_core_lgT, &
                        s% fe_core_lgRho, s% fe_core_L, s% fe_core_v, &
                        s% fe_core_omega, s% fe_core_omega_div_omega_crit)
                     exit
                  end if
               else if (si28 /= 0) then
                  if (s% xa(si28,k) <= s% fe_core_boundary_si28_fraction) then
                     sumx = 0
                     do j=1,species
                        if (chem_isos% Z_plus_N(s% chem_id(j)) >= A_max_fe) &
                           sumx = sumx + s% xa(j,k)
                     end do
                     if (sumx >= min_x) then
                        call set_core_info(s, k, s% fe_core_k, &
                           s% fe_core_mass, s% fe_core_radius, s% fe_core_lgT, &
                           s% fe_core_lgRho, s% fe_core_L, s% fe_core_v, &
                           s% fe_core_omega, s% fe_core_omega_div_omega_crit)
                        exit
                     end if
                  end if
               end if
            end if

            if (.not. have_one) then
               if (s% one_core_boundary_he4_c12_fraction < 0) then
                  if (A_max >= A_max_one) then
                     call set_core_info(s, k, s% one_core_k, &
                        s% one_core_mass, s% one_core_radius, s% one_core_lgT, &
                        s% one_core_lgRho, s% one_core_L, s% one_core_v, &
                        s% one_core_omega, s% one_core_omega_div_omega_crit)
                     have_one = .true.
                  end if
               else if (he4 /= 0 .and. c12 /= 0 .and. o16 /= 0 .and. ne20 /=0) then
                  if (s% xa(he4,k)+s% xa(c12,k) <= s% one_core_boundary_he4_c12_fraction .and. &
                      s% xa(o16,k)+s% xa(ne20,k) >= min_x) then
                     call set_core_info(s, k, s% one_core_k, &
                        s% one_core_mass, s% one_core_radius, s% one_core_lgT, &
                        s% one_core_lgRho, s% one_core_L, s% one_core_v, &
                        s% one_core_omega, s% one_core_omega_div_omega_crit)
                     have_one = .true.
                  end if
               end if
            end if

            if (.not. have_co) then
               if (s% co_core_boundary_he4_fraction < 0) then
                  if (A_max >= A_max_co) then
                     call set_core_info(s, k, s% co_core_k, &
                        s% co_core_mass, s% co_core_radius, s% co_core_lgT, &
                        s% co_core_lgRho, s% co_core_L, s% co_core_v, &
                        s% co_core_omega, s% co_core_omega_div_omega_crit)
                     have_co = .true.
                  end if
               else if (he4 /= 0 .and. c12 /= 0 .and. o16 /= 0) then
                  if (s% xa(he4,k) <= s% co_core_boundary_he4_fraction .and. &
                      s% xa(c12,k)+s% xa(o16,k) >= min_x) then
                     call set_core_info(s, k, s% co_core_k, &
                        s% co_core_mass, s% co_core_radius, s% co_core_lgT, &
                        s% co_core_lgRho, s% co_core_L, s% co_core_v, &
                        s% co_core_omega, s% co_core_omega_div_omega_crit)
                     have_co = .true.
                  end if
               end if
            end if

            if (.not. have_he) then
               if (s% he_core_boundary_h1_fraction < 0) then
                  if (A_max >= A_max_he) then
                     call set_core_info(s, k, s% he_core_k, &
                        s% he_core_mass, s% he_core_radius, s% he_core_lgT, &
                        s% he_core_lgRho, s% he_core_L, s% he_core_v, &
                        s% he_core_omega, s% he_core_omega_div_omega_crit)
                     have_he = .true.
                  end if
               else if (h1 /= 0 .and. he4 /= 0) then
                  if (s% xa(h1,k) <= s% he_core_boundary_h1_fraction .and. &
                      s% xa(he4,k) >= min_x) then
                     call set_core_info(s, k, s% he_core_k, &
                        s% he_core_mass, s% he_core_radius, s% he_core_lgT, &
                        s% he_core_lgRho, s% he_core_L, s% he_core_v, &
                        s% he_core_omega, s% he_core_omega_div_omega_crit)
                     have_he = .true.
                  end if
               end if
            end if

         end do

      end subroutine get_core_info


      subroutine clear_core_info( &
            core_k, core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit)
         integer, intent(out) :: core_k
         real(dp), intent(out) :: &
            core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit
         core_k = 0
         core_m = 0
         core_r = 0
         core_lgT = 0
         core_lgRho = 0
         core_L = 0
         core_v = 0
         core_omega = 0
         core_omega_div_omega_crit = 0
      end subroutine clear_core_info


      subroutine set_core_info(s, k, &
            core_k, core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: core_k
         real(dp), intent(out) :: &
            core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_omega, core_omega_div_omega_crit

         integer :: j, jm1, j00
         real(dp) :: dm1, d00, qm1, q00, core_q, &
            core_lgP, core_g, core_X, core_Y, &
            core_scale_height, core_dlnX_dr, core_dlnY_dr, core_dlnRho_dr

         include 'formats'

         if (k == 1) then
            core_q = 1d0
         else
            jm1 = maxloc(s% xa(:,k-1), dim=1)
            j00 = maxloc(s% xa(:,k), dim=1)
            qm1 = s% q(k-1) - 0.5d0*s% dq(k-1) ! center of k-1
            q00 = s% q(k) - 0.5d0*s% dq(k) ! center of k
            dm1 = s% xa(j00,k-1) - s% xa(jm1,k-1)
            d00 = s% xa(j00,k) - s% xa(jm1,k)
            if (dm1*d00 > 0d0) then
               write(*,2) 'bad args for set_core_info', k, dm1, d00
               call mesa_error(__FILE__,__LINE__)
               core_q = 0.5d0*(qm1 + q00)
            else if (dm1 == 0d0 .and. d00 == 0d0) then
               core_q = 0.5d0*(qm1 + q00)
            else if (dm1 == 0d0) then
               core_q = qm1
            else if (d00 == 0d0) then
               core_q = q00
            else
               core_q = find0(qm1, dm1, q00, d00)
            end if
         end if

         call get_info_at_q(s, core_q, &
            core_k, core_m, core_r, core_lgT, core_lgRho, core_L, core_v, &
            core_lgP, core_g, core_X, core_Y, &
            core_scale_height, core_dlnX_dr, core_dlnY_dr, core_dlnRho_dr, &
            core_omega, core_omega_div_omega_crit)

      end subroutine set_core_info


      subroutine get_info_at_q(s, bdy_q, &
            kbdy, bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, &
            bdy_lgP, bdy_g, bdy_X, bdy_Y, &
            bdy_scale_height, bdy_dlnX_dr, bdy_dlnY_dr, bdy_dlnRho_dr, &
            bdy_omega, bdy_omega_div_omega_crit)

         type (star_info), pointer :: s
         real(dp), intent(in) :: bdy_q
         integer, intent(out) :: kbdy
         real(dp), intent(out) :: &
            bdy_m, bdy_r, bdy_lgT, bdy_lgRho, bdy_L, bdy_v, &
            bdy_lgP, bdy_g, bdy_X, bdy_Y, &
            bdy_scale_height, bdy_dlnX_dr, bdy_dlnY_dr, bdy_dlnRho_dr, &
            bdy_omega, bdy_omega_div_omega_crit

         real(dp) :: x, x0, x1, x2, alfa, beta, bdy_omega_crit
         integer :: k, ii, klo, khi

         include 'formats'

         bdy_m=0; bdy_r=0; bdy_lgT=0; bdy_lgRho=0; bdy_L=0; bdy_v=0
         bdy_lgP=0; bdy_g=0; bdy_X=0; bdy_Y=0
         bdy_scale_height=0; bdy_dlnX_dr=0; bdy_dlnY_dr=0; bdy_dlnRho_dr=0
         bdy_omega=0; bdy_omega_div_omega_crit=0
         kbdy = 0

         if (bdy_q <= 0) return
         k = k_for_q(s,bdy_q)
         if (k >= s% nz) then
            kbdy = s% nz
            return
         end if
         if (k <= 1) then
            bdy_m = s% star_mass
            bdy_r = s% r(1)/Rsun
            bdy_lgT = s% lnT(1)/ln10
            bdy_lgRho = s% lnd(1)/ln10
            bdy_L = s% L(1)/Lsun
            if (s% u_flag) then
               bdy_v = s% u(1)
            else if (s% v_flag) then
               bdy_v = s% v(1)
            end if
            bdy_lgP = s% lnPeos(1)/ln10
            bdy_g = s% grav(1)
            bdy_X = s% X(1)
            bdy_Y = s% Y(1)
            bdy_scale_height = s% scale_height(1)
            bdy_omega = s% omega(k)
            if (s% rotation_flag) then
               bdy_omega_crit = omega_crit(s,1)
               bdy_omega_div_omega_crit = min(1d0,bdy_omega/bdy_omega_crit)
            else
               bdy_omega_crit = 0
               bdy_omega_div_omega_crit = 0
            end if
            kbdy = 1
            return
         end if

         kbdy = k+1

         bdy_m = (s% M_center + s% xmstar*bdy_q)/Msun

         x = s% q(k-1) - bdy_q
         x0 = s% dq(k-1)/2
         x1 = s% dq(k)/2 + s% dq(k-1)
         x2 = s% dq(k+1)/2 + s% dq(k) + s% dq(k-1)

         alfa = max(0d0, min(1d0, (bdy_q - s% q(k+1))/s% dq(k)))

         bdy_lgT = interp3(s% lnT(k-1), s% lnT(k), s% lnT(k+1))/ln10
         bdy_lgRho = interp3(s% lnd(k-1), s% lnd(k), s% lnd(k+1))/ln10
         bdy_lgP = interp3(s% lnPeos(k-1), s% lnPeos(k), s% lnPeos(k+1))/ln10
         bdy_X = interp3(s% X(k-1), s% X(k), s% X(k+1))
         bdy_Y = interp3(s% Y(k-1), s% Y(k), s% Y(k+1))

         bdy_r = pow( &
            interp2(s% r(k)*s% r(k)*s% r(k), s% r(k+1)*s% r(k+1)*s% r(k+1)),one_third)/Rsun
         bdy_L = interp2(s% L(k), s% L(k+1))/Lsun
         bdy_g = interp2(s% grav(k), s% grav(k+1))
         bdy_scale_height = interp2(s% scale_height(k), s% scale_height(k+1))

         klo = k-1
         khi = k+1
         bdy_dlnX_dr = log(max(1d-99,max(1d-99,s% X(klo))/max(1d-99,s% X(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
         bdy_dlnY_dr = log(max(1d-99,max(1d-99,s% Y(klo))/max(1d-99,s% Y(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
         bdy_dlnRho_dr = (s% lnd(klo) - s% lnd(khi))/(s% rmid(klo) - s% rmid(khi))

         if (s% u_flag) then
            bdy_v = interp2(s% u(k), s% u(k+1))
         else if (s% v_flag) then
            bdy_v = interp2(s% v(k), s% v(k+1))
         end if
         if (s% rotation_flag) then
            bdy_omega = interp2(s% omega(k), s% omega(k+1))
            bdy_omega_crit = interp2(omega_crit(s,k), omega_crit(s,k+1))
            bdy_omega_div_omega_crit = min(1d0,bdy_omega/bdy_omega_crit)
         else
            bdy_omega = 0
            bdy_omega_crit = 0
            bdy_omega_div_omega_crit = 0
         end if

         contains

         real(dp) function interp3(f0, f1, f2)
            real(dp), intent(in) :: f0, f1, f2
            real(dp) :: fmin, fmax
            fmin = min(f0,f1,f2)
            fmax = max(f0,f1,f2)
            interp3 = (f1*(x-x0)*(x-x2)*(x0-x2)-&
                  (x-x1)*(f0*(x-x2)*(x1-x2) + (x-x0)*(x0-x1)*x2))/ &
                     ((x0-x1)*(x0-x2)*(-x1+x2))
            interp3 = min(fmax, max(fmin, interp3))
         end function interp3

         real(dp) function interp2(f0, f1)
            real(dp), intent(in) :: f0, f1
            interp2 = alfa*f0 + (1-alfa)*f1
         end function interp2

      end subroutine get_info_at_q


      subroutine get_mixing_regions(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: prev_type, cur_type, cur_top, n, k
         include 'formats'
         ierr = 0
         cur_type = s% mixing_type(1)
         cur_top = 1
         n = 0
         do k = 2, s% nz
            prev_type = cur_type
            cur_type = s% mixing_type(k)
            if (cur_type == prev_type .and. k < s% nz) cycle
            ! change of type from k-1 to k
            if (prev_type /= no_mixing) then
               n = n + 1
               s% mixing_region_type(n) = prev_type
               s% mixing_region_top(n) = cur_top
               if (k == s% nz) then
                  s% mixing_region_bottom(n) = k
               else
                  s% mixing_region_bottom(n) = k-1
               end if
               if (n == max_num_mixing_regions) exit
            end if
            cur_top = k
         end do
         s% num_mixing_regions = n
      end subroutine get_mixing_regions


      end module report
