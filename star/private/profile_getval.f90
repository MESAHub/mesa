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

      module profile_getval

      use star_private_def
      use star_profile_def
      use const_def
      use star_utils
      use utils_lib

      implicit none


      integer, parameter :: idel = 10000

      integer, parameter :: add_abundances = idel
      integer, parameter :: add_log_abundances = add_abundances + 1
      integer, parameter :: category_offset = add_log_abundances + 1
      integer, parameter :: abundance_offset = category_offset + idel
      integer, parameter :: log_abundance_offset = abundance_offset + idel
      integer, parameter :: xadot_offset = log_abundance_offset + idel
      integer, parameter :: xaprev_offset = xadot_offset + idel
      integer, parameter :: ionization_offset = xaprev_offset + idel
      integer, parameter :: typical_charge_offset = ionization_offset + idel
      integer, parameter :: edv_offset = typical_charge_offset + idel
      integer, parameter :: extra_diffusion_factor_offset = edv_offset + idel
      integer, parameter :: v_rad_offset = extra_diffusion_factor_offset + idel
      integer, parameter :: log_g_rad_offset = v_rad_offset + idel
      integer, parameter :: log_concentration_offset = log_g_rad_offset + idel
      integer, parameter :: diffusion_dX_offset = log_concentration_offset + idel
      integer, parameter :: diffusion_D_offset = diffusion_dX_offset + idel
      integer, parameter :: extra_offset = diffusion_D_offset + idel

      integer, parameter :: max_profile_offset = extra_offset + idel


      contains


      integer function do1_profile_spec( &
            iounit, n, i, string, buffer, report, ierr) result(spec)

         use utils_lib
         use utils_def
         use chem_def
         use chem_lib
         integer :: iounit, n, i, num, t

         character (len=*) :: string, buffer
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         integer :: id

         ierr = 0
         spec = -1

         id = do_get_profile_id(string)
         if (id > 0) then
            spec = id
            return
         end if

         select case(string)

            case ('xadot')
               call do1_nuclide(xadot_offset)

            case ('xaprev')
               call do1_nuclide(xaprev_offset)

            case ('ionization')
               call do1_nuclide(ionization_offset)

            case ('typical_charge')
               call do1_nuclide(typical_charge_offset)

            case ('edv')
               call do1_nuclide(edv_offset)

            case ('extra_diffusion_factor')
               call do1_nuclide(extra_diffusion_factor_offset)

            case ('v_rad')
               call do1_nuclide(v_rad_offset)

            case ('log_g_rad')
               call do1_nuclide(log_g_rad_offset)

            case ('log_concentration')
               call do1_nuclide(log_concentration_offset)

            case ('diffusion_dX')
               call do1_nuclide(diffusion_dX_offset)

            case ('diffusion_D')
               call do1_nuclide(diffusion_D_offset)

            case ('log') ! add log of abundance
               call do1_nuclide(log_abundance_offset)

            case ('extra')

               t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  ierr = -1; return
               end if
               read(string,fmt=*,iostat=ierr) num
               if (ierr /= 0 .or. num <= 0 .or. num > max_num_profile_extras) then
                  write(*,*) 'failed to find valid integer for extra: ' // trim(string)
                  ierr = -1
               end if
               spec = extra_offset + num

            case default

               id = chem_get_iso_id(string)
               if (id > 0) then
                  spec = abundance_offset + id
                  return
               end if
               id = rates_category_id(string)
               if (id > 0) then
                  spec = category_offset + id
                  return
               end if
               if (report) &
                  write(*,*) 'failed to recognize item for profile columns: ' // trim(string)
               ierr = -1

         end select


         contains


         subroutine do1_nuclide(offset)
            integer, intent(in) :: offset
            integer :: t, id
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token) then
               ierr = -1; return
            end if
            id = chem_get_iso_id(string)
            if (id > 0) then
               spec = offset + id
               return
            end if
            write(*,*) 'bad iso name: ' // trim(string)
            ierr = -1
         end subroutine do1_nuclide


      end function do1_profile_spec


      integer function get_profile_id(s, name) result(spec)
         use utils_lib, only: token
         use utils_def, only: name_token
         type (star_info), pointer :: s
         character (len=*), intent(in) :: name
         character (len=strlen) :: buffer, string
         integer :: i, n, iounit, ierr, t
         iounit = -1
         ierr = 0
         buffer = name
         n = len_trim(buffer) + 1
         buffer(n:n) = ' '
         i = 0
         t = token(iounit, n, i, buffer, string)
         if (t /= name_token) then
            spec = -1; return
         end if
         spec = do1_profile_spec(iounit, n, i, string, buffer, .false., ierr)
         if (ierr == 0) return
         ! check to see if it is one of the extra profile columns
         do i=1,s% num_extra_profile_cols
            if (name == s% extra_profile_col_names(i)) then
               spec = i + max_profile_offset
               return
            end if
         end do
         spec = -1
      end function get_profile_id


      real(dp) function get_profile_val(s, id, k)
         type (star_info), pointer :: s
         integer, intent(in) :: id, k
         integer :: ii, int_val
         logical :: int_flag
         if (id > max_profile_offset) then ! get from extras
            get_profile_val = s% extra_profile_col_vals(k, id - max_profile_offset)
            return
         end if
         call getval_for_profile(s, id, k, get_profile_val, int_flag, int_val)
         if (int_flag) get_profile_val = dble(int_val)
      end function get_profile_val






      subroutine getval_for_profile(s, c, k, val, int_flag, int_val)
         use chem_def
         use rates_def
         use ionization_def
         use ionization_lib, only: eval_typical_charge
         use rsp_def, only: rsp_WORK, rsp_WORKQ, rsp_WORKT, rsp_WORKC
         type (star_info), pointer :: s
         integer, intent(in) :: c, k
         real(dp), intent(out) :: val
         integer, intent(out) :: int_val
         logical, intent(inout) :: int_flag

         real(dp) :: cno, z, x, frac, eps, eps_alt, L_rad, L_edd, Pbar_00, Pbar_p1, &
            P_face, rho_face, dr, v, r00, rp1, v00, vp1, A00, Ap1, Amid, &
            r00_start, rp1_start, dr3, dr3_start, d_drL, d_drR, flxR, mmid, &
            d_dlnR00, d_dlnRp1, d_dv00, d_dvp1
         integer :: j, nz, ionization_k, klo, khi, i, ii, kk, ierr
         real(dp) :: ionization_res(num_ion_vals)
         real(dp) :: f, lgT, full_on, full_off, am_nu_factor, Lconv, conv_vel
         logical :: rsp_or_eturb
         include 'formats'

         if (s% rotation_flag) then
            full_on = s% D_mix_rotation_max_logT_full_on
            full_off = s% D_mix_rotation_min_logT_full_off
            lgT = s% lnT(k)/ln10
            if (lgT <= full_on) then
               f = 1d0
            else if (lgT >= full_off) then
               f = 0d0
            else                   ! lgT > full_on and < full_off
               f = (lgT - full_on) / (full_off - full_on)
            end if
            am_nu_factor = f*s% am_nu_factor
         else
            am_nu_factor = 1d0
         end if

         val = 0; int_val = 0; int_flag = .false.
         nz = s% nz
         ionization_k = 0

         int_flag = .false.
         rsp_or_eturb = s% RSP_flag .or. s% Eturb_flag

         if (c > extra_offset) then
            i = c - extra_offset
            val = s% profile_extra(k,i)
         else if (c > diffusion_D_offset) then
            i = c - diffusion_D_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% diffusion_D_self(ii,k)
         else if (c > diffusion_dX_offset) then
            i = c - diffusion_dX_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% diffusion_dX(ii,k)
         else if (c > log_concentration_offset) then
            i = c - log_concentration_offset
            ii = s% net_iso(i)
            if (ii > 0) val = get_log_concentration(s,ii,k)
         else if (c > log_g_rad_offset) then
            i = c - log_g_rad_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = safe_log10(s% g_rad(ii,k))
         else if (c > v_rad_offset) then
            i = c - v_rad_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% v_rad(ii,k)
         else if (c > extra_diffusion_factor_offset) then
            i = c - extra_diffusion_factor_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% extra_diffusion_factor(ii,k)
         else if (c > edv_offset) then
            i = c - edv_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% edv(ii,k)
         else if (c > typical_charge_offset) then
            i = c - typical_charge_offset
            ii = s% net_iso(i)
            if (ii > 0 .and. s% do_element_diffusion) val = s% typical_charge(ii,k)
         else if (c > ionization_offset) then
            i = c - ionization_offset
            ii = s% net_iso(i)
            val = eval_typical_charge( &
               i, s% abar(k), exp(s% lnfree_e(k)), &
               s% T(k), s% lnT(k)/ln10, s% rho(k), s% lnd(k)/ln10)
         else if (c > xaprev_offset) then
            i = c - xaprev_offset
            ii = s% net_iso(i)
            if (ii > 0) val = s% xa_start(ii,k)
         else if (c > xadot_offset) then
            i = c - xadot_offset
            ii = s% net_iso(i)
            if (ii > 0) val = s% xa(ii,k) - s% xa_start(ii,k)
         else if (c > log_abundance_offset) then
            i = c - log_abundance_offset
            ii = s% net_iso(i)
            if (ii > 0) then
               val = safe_log10(s% xa(ii,k))
            else
               val = -99d0
            end if
         else if (c > abundance_offset) then
            i = c - abundance_offset
            ii = s% net_iso(i)
            if (ii > 0) val = s% xa(ii,k)
         else if (c > category_offset) then
            i = c - category_offset
            val = s% eps_nuc_categories(i,k)
         else

         select case(c)
            case (p_zone)
               val = dble(k)
               int_val = k
               int_flag = .true.
            case (p_k)
               val = dble(k)
               int_val = k
               int_flag = .true.
            case (p_conv_L_div_L)
               if (s% L(k) > 0d0) val = s% L_conv(k)/s% L(k)
            case (p_log_conv_L_div_L)
               if (s% L(k) > 0d0) val = safe_log10(s% L_conv(k)/s% L(k))
            case (p_lum_erg_s)
               val = s% L(k)
            case (p_luminosity)
               val = s% L(k)/Lsun
            case (p_log_abs_lum_erg_s)
               val = safe_log10(abs(s% L(k)))
            case (p_lum_adv)
               val = get_Ladv(s,k)
            case (p_lum_plus_lum_adv)
               val = s% L(k) + get_Ladv(s,k)
            case (p_lum_rad)
               L_rad = get_Lrad(s,k)
               val = L_rad/Lsun
            case (p_lum_conv)
               L_rad = get_Lrad(s,k)
               val = (s% L(k) - L_rad)/Lsun
            case (p_lum_conv_MLT)
               L_rad = get_Lrad(s,k)
               val = s% L_conv(k)/Lsun

            !case (p_lum_rad_div_L_Edd_sub_fourPrad_div_PchiT)
            !   val = get_Lrad_div_Ledd(s,k) - 4*s% Prad(k)/(s% P(k)*s% chiT(k))
            case (p_lum_rad_div_L_Edd)
               val = get_Lrad_div_Ledd(s,k)
            case (p_lum_conv_div_lum_Edd)
               L_rad = get_Lrad(s,k)
               L_edd = get_Ledd(s,k)
               val = (s% L(k) - L_rad)/L_edd

            case (p_lum_conv_div_lum_rad)
               L_rad = get_Lrad(s,k)
               val = (s% L(k) - L_rad)/L_rad

            case (p_lum_rad_div_L)
               L_rad = get_Lrad(s,k)
               val = L_rad/max(1d0,s% L(k))
            case (p_lum_conv_div_L)
               L_rad = get_Lrad(s,k)
               val = (s% L(k) - L_rad)/max(1d0,s% L(k))

            case (p_log_Lrad)
               L_rad = get_Lrad(s,k)
               val = safe_log10(L_rad/Lsun)
            case (p_log_Lconv)
               L_rad = get_Lrad(s,k)
               val = safe_log10((s% L(k) - L_rad)/Lsun)
            case (p_log_Lconv_div_L)
               L_rad = get_Lrad(s,k)
               val = safe_log10((s% L(k) - L_rad)/s% L(k))

            case (p_log_Lrad_div_L)
               L_rad = get_Lrad(s,k)
               val = safe_log10(L_rad/s% L(k))
            case (p_log_Lrad_div_Ledd)
               val = safe_log10(get_Lrad_div_Ledd(s,k))

            case (p_log_g)
               val = safe_log10(s% grav(k))
            case (p_grav)
               val = s% grav(k)
            case (p_g_div_r)
               val = s% grav(k)/s% r(k)
            case (p_r_div_g)
               val = s% r(k)/s% grav(k)
            case (p_signed_log_eps_grav)
               val = s% eps_grav(k)
               val = sign(1d0,val)*log10(max(1d0,abs(val)))
            case (p_net_nuclear_energy)
               val = s% eps_nuc(k) - s% eps_nuc_neu_total(k) - s% non_nuc_neu(k)
               val = sign(1d0,val)*log10(max(1d0,abs(val)))
            case (p_eps_nuc_plus_nuc_neu)
               val = s% eps_nuc(k) + s% eps_nuc_neu_total(k)
            case (p_eps_nuc_minus_non_nuc_neu)
               val = s% eps_nuc(k) - s% non_nuc_neu(k)
            case (p_net_energy)
               val = s% eps_nuc(k) - s% non_nuc_neu(k) + s% eps_grav(k)
               val = sign(1d0,val)*log10(max(1d0,abs(val)))
            case (p_signed_log_power)
               val = s% L(k)
               val = sign(1d0,val)*log10(max(1d0,abs(val)))
            case (p_logL)
               val = safe_log10(max(1d-12,s% L(k)/Lsun))
            case (p_log_Ledd)
               val = safe_log10(get_Ledd(s,k)/Lsun)
            case (p_lum_div_Ledd)
               val = s% L(k)/get_Ledd(s,k)
            case (p_log_L_div_Ledd)
               val = safe_log10(max(1d-12,s% L(k)/get_Ledd(s,k)))
            case (p_log_abs_dvdt_div_v)
               if (s% u_flag) then
                  val = safe_log10(abs(s% du_dt(k))/max(1d-50,abs(s% u(k))))
               else if (s% v_flag) then
                  val = safe_log10(abs(s% dv_dt(k))/max(1d-50,abs(s% v(k))))
               end if
            case (p_log_abs_v)
               if (s% u_flag) then
                  val = safe_log10(abs(s% u(k)))
               else if (s% v_flag) then
                  val = safe_log10(abs(s% v(k)))
               end if
            
            case (p_superad_reduction_factor)
               val = s% superad_reduction_factor(k)
            case (p_gradT_excess_effect)
               val = s% gradT_excess_effect(k)
            case (p_diff_grads)
               val = s% gradr(k) - s% gradL(k) ! convective if this is > 0
            case (p_log_diff_grads)
               val = safe_log10(abs(s% gradr(k) - s% gradL(k)))
            case (p_v)
               if (s% u_flag) then
                  val = s% u(k)
               else if (s% v_flag) then
                  val = s% v(k)
               end if
            case (p_velocity)
               if (s% u_flag) then
                  val = s% u(k)
               else if (s% v_flag) then
                  val = s% v(k)
               end if
            case (p_vel_km_per_s)
               if (s% u_flag) then
                  val = s% u(k)*1d-5
               else if (s% v_flag) then
                  val = s% v(k)*1d-5
               end if
            case (p_v_div_r)
               if (s% u_flag) then
                  val = s% u_face(k)/s% r(k)
               else if (s% v_flag) then
                  val = s% v(k)/s% r(k)
               end if

            case (p_v_times_t_div_r)
               if (s% u_flag) then
                  val = s% u_face(k)*s% time/s% r(k)
               else if (s% v_flag) then
                  val = s% v(k)*s% time/s% r(k)
               end if
            case (p_radius)
               val = s% r(k)/Rsun
            case (p_radius_cm)
               val = s% r(k)
            case (p_radius_km)
               val = s% r(k)*1d-5
            case (p_rmid)
               val = s% rmid(k)/Rsun
            case (p_logR_cm)
               val = safe_log10(s% r(k))
            case (p_logR)
               val = safe_log10(s% r(k)/Rsun)

            case (p_q)
               val = s% q(k)
            case (p_log_q)
               val = safe_log10(s% q(k))
            case (p_dq)
               val = s% dq(k)
            case (p_log_dq)
               val = safe_log10(s% dq(k))
            case (p_mass)
               val = s% m(k)/Msun
            case (p_log_mass)
               val = safe_log10(s% m(k)/Msun)
            case (p_mass_grams)
               val = s% m(k)
            case (p_mmid)
               val = (s% M_center + s% xmstar*(s% q(k) - s% dq(k)/2))/Msun

            case (p_dm)
               val = s% dm(k)
            case (p_dm_bar)
               val = s% dm_bar(k)

            case (p_m_div_r)
               val = s% m(k)/s% r(k)
            case (p_dmbar_m_div_r)
               val = s% dm_bar(k)*s% m(k)/s% r(k)
            case (p_log_dmbar_m_div_r)
               val = safe_log10(s% dm_bar(k)*s% m(k)/s% r(k))

            case (p_m_grav)
               val = s% m_grav(k)/Msun
            case (p_m_grav_div_m_baryonic)
               val = s% m_grav(k)/s% m(k)
            case (p_mass_correction_factor)
               val = s% mass_correction(k)

            case (p_xr)
               val = (s% r(1) - s% r(k))/Rsun
            case (p_xr_cm)
               val = s% r(1) - s% r(k)
            case (p_xr_div_R)
               val = (s% r(1) - s% r(k))/s% r(1)
            case (p_log_xr)
               val = safe_log10((s% r(1) - s% r(k))/Rsun)
            case (p_log_xr_cm)
               val = safe_log10(s% r(1) - s% r(k))
            case (p_log_xr_div_R)
               val = safe_log10((s% r(1) - s% r(k))/s% r(1))

            case (p_x)
               val = s% X(k)
            case (p_log_x)
               val = safe_log10(s% X(k))
            case (p_y)
               val = s% Y(k)
            case (p_log_y)
               val = safe_log10(s% Y(k))
            case (p_z)
               val = s% Z(k)
            case (p_log_z)
               val = safe_log10(s% Z(k))
            case (p_xm)
               val = sum(s% dm(1:k-1))/Msun
            case (p_logxm)
               val = safe_log10(sum(s% dm(1:k-1))/Msun)
            case (p_xq)
               val = sum(s% dq(1:k-1))
            case (p_logxq)
               val = safe_log10(sum(s% dq(1:k-1)))
            case (p_logdq)
               val = safe_log10(s% dq(k))
            case (p_log_column_depth)
               val = safe_log10(s% xmstar*sum(s% dq(1:k-1))/(4*pi*s% r(k)*s% r(k)))
            case (p_log_radial_depth)
               val = safe_log10(s% r(1) - s% r(k))

            case (p_r_div_R)
               val = s% r(k)/s% r(1)
            case (p_log_dr)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10(val)
            case (p_dlogR)
               if (k == s% nz) then
                  val = s% lnR(k) - log(max(1d0,s% R_center))
               else
                  val = s% lnR(k) - s% lnR(k+1)
               end if
               val = val/ln10
            case (p_dr_div_rmid)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = val/s% rmid(k)
            case (p_log_dr_div_rmid)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10(val/s% rmid(k))
            case (p_log_acoustic_radius)
               val = safe_log10(sum(s% dr_div_csound(k:nz)))
            case (p_acoustic_radius)
               val = sum(s% dr_div_csound(k:nz))
            case (p_log_acoustic_depth)
               if (k > 1) &
                  val = sum(s% dr_div_csound(1:k-1))
               val = safe_log10(val)
            case (p_acoustic_depth)
               if (k > 1) &
                  val = sum(s% dr_div_csound(1:k-1))
            case (p_acoustic_r_div_R_phot)
               val = sum(s% dr_div_csound(k:nz))/s% photosphere_acoustic_r

            case (p_lnR_residual)
               if (s% i_dlnR_dt > 0) val = s% lnR_residual(k)
            case (p_log_lnR_residual)
               if (s% i_dlnR_dt > 0) &
                  val = safe_log10(abs(s% lnR_residual(k)))

            case (p_lnd_residual)
               if (s% i_dlnd_dt > 0) val = s% lnd_residual(k)
            case (p_log_lnd_residual)
               if (s% i_dlnd_dt > 0) &
                  val = safe_log10(abs(s% lnd_residual(k)))

            case (p_equL_residual)
               if (k > 1 .and. s% i_equL > 0) val = s% equL_residual(k)
            case (p_log_equL_residual)
               if (k > 1 .and. s% i_equL > 0) &
                  val = safe_log10(abs(s% equL_residual(k)))
               
            case (p_E_residual)
               if (s% i_dlnE_dt > 0) val = s% E_residual(k)
            case (p_log_E_residual)
               if (s% i_dlnE_dt > 0) val = safe_log10(s% E_residual(k))

            case (p_v_residual)
               if (s% i_dv_dt > 0 .and. k > 1) val = s% v_residual(k)
            case (p_log_v_residual)
               if (s% i_dv_dt > 0 .and. k > 1) &
                  val = safe_log10(abs(s% v_residual(k)))

            case (p_dvdt_residual)
               if (s% i_dv_dt > 0 .and. k > 1) val = s% v_residual(k)
            case (p_log_dvdt_residual)
               if (s% i_dv_dt > 0 .and. k > 1) &
                  val = safe_log10(abs(s% v_residual(k)))

            case (p_ergs_error)
               val = s% ergs_error(k)
            case (p_log_rel_E_err)
               val = safe_log10(abs(s% ergs_error(k)/s% total_energy_start))
            case (p_ergs_error_integral)
               val = sum(s% ergs_error(1:k))
            case (p_ergs_rel_error_integral)
               if (s% total_energy_end /= 0d0) &
                  val = sum(s% ergs_error(1:k))/s% total_energy_end
               
            case (p_cell_internal_energy_fraction)
               val = s% energy(k)*s% dm(k)/s% total_internal_energy_end
            case (p_cell_internal_energy_fraction_start)
               val = s% energy_start(k)*s% dm(k)/s% total_internal_energy_start
                  
            case (p_dr_div_R)
               if (k < s% nz) then
                  val = (s% r(k) - s% r(k+1))/s% r(1)
               else
                  val = (s% r(k) - s% r_center)/s% r(1)
               end if
            case (p_dRstar_div_dr)
               if (k < s% nz) then
                  val = (s% r(1) - s% R_center)/(s% r(k) - s% r(k+1))
               else
                  val = (s% r(1) - s% R_center)/(s% r(k) - s% r_center)
               end if
            case (p_log_dr_div_R)
               if (k < s% nz) then
                  val = (s% r(k) - s% r(k+1))/s% r(1)
               else
                  val = (s% r(k) - s% r_center)/s% r(1)
               end if
               val = safe_log10(val)
            
            case(p_t_rad)
               val = 1d0/(clight*s% opacity(k)*s% rho(k))
            case(p_log_t_rad)
               val = log10(1d0/(clight*s% opacity(k)*s% rho(k)))
            case (p_dt_cs_div_dr)
               if (k < s% nz) then
                  val = s% r(k) - s% r(k+1)
               else
                  val = s% r(k) - s% r_center
               end if
               val = s% dt*s% csound(k)/val
            case (p_log_dt_cs_div_dr)
               if (k < s% nz) then
                  val = s% r(k) - s% r(k+1)
               else
                  val = s% r(k) - s% r_center
               end if
               val = s% dt*s% csound(k)/val
               val = safe_log10(val)
            case (p_dr_div_cs)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = val/s% csound(k)
            case (p_log_dr_div_cs)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10(val/s% csound(k))

            case (p_dr_div_cs_yr)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = val/s% csound(k)/secyer
            case (p_log_dr_div_cs_yr)
               if (k == s% nz) then
                  val = s% r(k) - s% R_center
               else
                  val = s% r(k) - s% r(k+1)
               end if
               val = safe_log10(val/s% csound(k)/secyer)

            case (p_pgas_div_ptotal)
               val = s% Pgas(k)/s% P(k)
            case (p_prad_div_pgas)
               val = s% Prad(k)/s% Pgas(k)
            case(p_prad_div_pgas_div_L_div_Ledd)
               val = (s% Prad(k)/s% Pgas(k))/max(1d-12,s% L(k)/get_Ledd(s,k))
            case (p_pgas_div_p)
               val = s% Pgas(k)/s% P(k)

            case (p_cell_collapse_time)
               if (s% v_flag) then
                  if (k == s% nz) then
                     rp1 = s% R_center
                     vp1 = s% v_center
                  else
                     rp1 = s% r(k+1)
                     vp1 = s% v(k+1)
                  end if
                  r00 = s% r(k)
                  v00 = s% v(k)
                  if (vp1 > v00) val = (r00 - rp1)/(vp1 - v00)
               end if

            case (p_log_cell_collapse_time)
               if (s% v_flag) then
                  if (k == s% nz) then
                     rp1 = s% R_center
                     vp1 = s% v_center
                  else
                     rp1 = s% r(k+1)
                     vp1 = s% v(k+1)
                  end if
                  r00 = s% r(k)
                  v00 = s% v(k)
                  if (vp1 > v00) val = (r00 - rp1)/(vp1 - v00)
               end if
               val = safe_log10(val)

            case (p_dq_ratio)
               if (k == 1 .or. k == s% nz) then
                  val = 1
               else
                  val = s% dq(k-1)/s% dq(k)
               end if

            case (p_compression_gradient)
               if (k == 1) then
                  val = s% rho_start(1)*&
                     ((1/s% rho(1) - 1/s% rho_start(1)) - (1/s% rho(2) - 1/s% rho_start(2)))
               else
                  val = s% rho_start(k-1)*&
                     ((1/s% rho(k-1) - 1/s% rho_start(k-1)) - (1/s% rho(k) - 1/s% rho_start(k)))
               end if

            case (p_tau)
               val = s% tau(k)
            case (p_logtau)
               val = s% lntau(k)/ln10
            case (p_xtau)
               val = s% tau(nz) - s% tau(k)
            case (p_xlogtau)
               val = safe_log10(s% tau(nz) - s% tau(k))
            case (p_logtau_sub_xlogtau)
               val = safe_log10(s% tau(k)) - safe_log10(s% tau(nz) - s% tau(k))

            case (p_tau_eff)
               val = tau_eff(s,k)
            case (p_tau_eff_div_tau)
               val = tau_eff(s,k)/s% tau(k)

            case (p_kap_frac_lowT)
               val = s% kap_frac_lowT(k)
            case (p_kap_frac_highT)
               val = s% kap_frac_highT(k)
            case (p_kap_frac_Type2)
               val = s% kap_frac_Type2(k)
            case (p_kap_frac_Compton)
               val = s% kap_frac_Compton(k)
            case (p_kap_frac_op_mono)
               val = s% kap_frac_op_mono(k)
            case (p_log_opacity)
               val = safe_log10(s% opacity(k))
            case (p_extra_opacity_factor)
               val = s% extra_opacity_factor(k)
            case (p_log_kap_times_factor)
               val = safe_log10(s% opacity(k)*s% extra_opacity_factor(k))
            case (p_energy)
               val = s% energy(k)
            case (p_logM)
               val = safe_log10(s% m(k)/Msun)
            case (p_temperature)
               val = s% T(k)
            case (p_logT)
               val = s% lnT(k)/ln10

            case (p_logT_face)
               if (k == 1) then
                  val = safe_log10(s% T_surf)
               else
                  val = (s% dq(k-1)*s% lnT(k) + &
                         s% dq(k)*s% lnT(k-1))/(s% dq(k-1) + s% dq(k))/ln10
               end if
            case (p_logT_bb)
               val = safe_log10( &
                        pow(s% L(k)/(4*pi*s% r(k)*s% r(k)*boltz_sigma), 0.25d0))
            case (p_logT_face_div_logT_bb)
               if (k == 1) then
                  val = safe_log10(s% Teff)
               else
                  val = (s% dq(k-1)*s% lnT(k) + &
                         s% dq(k)*s% lnT(k-1))/(s% dq(k-1) + s% dq(k))/ln10
               end if
               val = val / safe_log10( &
                        pow(s% L(k)/(4*pi*s% r(k)*s% r(k)*boltz_sigma), 0.25d0))

            case (p_density)
               val = s% rho(k)
            case (p_rho)
               val = s% rho(k)
            case (p_logRho)
               val = s% lnd(k)/ln10
            case (p_pgas)
               val = s% Pgas(k)
            case (p_logPgas)
               val = s% lnPgas(k)/ln10
            case (p_prad)
               val = s% Prad(k)
            case (p_pressure)
               val = s% P(k)
            case (p_logP)
               val = s% lnP(k)/ln10
            case (p_logE)
               val = s% lnE(k)/ln10
            case (p_grada)
               val = s% grada(k)
            case (p_dE_dRho)
               val = s% dE_dRho(k)
            case (p_cv)
               val = s% Cv(k)
            case (p_cp)
               val = s% cp(k)

            case (p_thermal_time_to_surface)
               if (s% L(1) > 0) &
                  val = sum(s% dm(1:k)*s% cp(1:k)*s% T(1:k))/s% L(1)
            case (p_log_thermal_time_to_surface)
               if (s% L(1) > 0) then
                  val = sum(s% dm(1:k)*s% cp(1:k)*s% T(1:k))/s% L(1)
                  val = safe_log10(val)
               end if

            case (p_log_CpT)
               val = safe_log10(s% cp(k)*s% T(k))
            case (p_log_CpT_absMdot_div_L)
               val = safe_log10(s% cp(k)*s% T(k)*abs(s% mstar_dot)/max(1d-99,s% L(k)))
            case (p_logS)
               val = s% lnS(k)/ln10
            case (p_logS_per_baryon)
               val = s% lnS(k)/ln10 + log10(amu)
            case (p_gamma1)
               val = s% gamma1(k)
            case (p_gamma3)
               val = s% gamma3(k)
            case (p_eta)
               val = s% eta(k)
            case (p_theta_e)
               val = s% theta_e(k)
            case (p_gam)
               val = s% gam(k)
            case (p_mu)
               val = s% mu(k)

            case (p_eos_frac_OPAL_SCVH)
               val = s% eos_frac_OPAL_SCVH(k)
            case (p_eos_frac_HELM)
               val = s% eos_frac_HELM(k)
            case (p_eos_frac_Skye)
               val = s% eos_frac_Skye(k)
            case (p_eos_frac_PC)
               val = s% eos_frac_PC(k)
            case (p_eos_frac_FreeEOS)
               val = s% eos_frac_FreeEOS(k)
            case (p_eos_frac_CMS)
               val = s% eos_frac_CMS(k)
               
            case (p_log_c_div_tau)
               val = safe_log10(clight/s% tau(k))
            case (p_log_v_escape)
               val = safe_log10(sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k))))
            case (p_v_div_vesc)
               if (s% u_flag) then
                  val = s% u(k)
               else if (s% v_flag) then
                  val = s% v(k)
               end if
               val = val/sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            case (p_v_div_v_escape)
               if (s% u_flag) then
                  val = s% u(k)
               else if (s% v_flag) then
                  val = s% v(k)
               end if
               val = val/sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            case (p_v_div_cs)
               val = s% v_div_csound(k)
            case (p_v_div_csound)
               val = s% v_div_csound(k)
            case (p_log_csound)
               val = safe_log10(s% csound(k))
            case (p_csound)
               val = s% csound(k)
            case (p_csound_face)
               val = s% csound_face(k)
            case (p_scale_height)
               val = s% scale_height(k)/Rsun
            case (p_entropy)
               val = s% entropy(k)
            case (p_free_e)
               val = exp(s% lnfree_e(k))
            case (p_logfree_e)
               val = s% lnfree_e(k)/ln10
            case (p_chiRho)
               val = s% chiRho(k)
            case (p_chiT)
               val = s% chiT(k)
            case (p_QQ)
               val = s% QQ(k)

            case (p_phase)
               val = s% phase(k)
            case (p_latent_ddlnT)
               val = s% latent_ddlnT(k)
            case (p_latent_ddlnRho)
               val = s% latent_ddlnRho(k)
               
            case (p_chiRho_for_partials)
               val = s% chiRho_for_partials(k)
            case (p_chiT_for_partials)
               val = s% chiT_for_partials(k)
            case (p_rel_diff_chiRho_for_partials)
               val = (s% chiRho_for_partials(k) - s% chiRho(k))/s% chiRho(k)
            case (p_rel_diff_chiT_for_partials)
               val = (s% chiT_for_partials(k) - s% chiT(k))/s% chiT(k)

            case (p_x_mass_fraction_H)
               val = s% X(k)
            case (p_y_mass_fraction_He)
               val = s% Y(k)
            case (p_z_mass_fraction_metals)
               val = s% Z(k)

            case (p_abar)
               val = s% abar(k)
            case (p_zbar)
               val = s% zbar(k)
            case (p_z2bar)
               val = s% z2bar(k)
            case (p_ye)
               val = s% ye(k)
            case (p_opacity)
               val = s% opacity(k)
            case (p_dkap_dlnrho_face)
               val = interp_val_to_pt(s% d_opacity_dlnd,k,nz,s% dq,'p_dkap_dlnrho_face')
            case (p_dkap_dlnT_face)
               val = interp_val_to_pt(s% d_opacity_dlnT,k,nz,s% dq,'p_dkap_dlnT_face')

            case (p_eps_nuc)
               val = s% eps_nuc(k)
            case (p_signed_log_eps_nuc)
               val = s% eps_nuc(k)
               val = sign(1d0,val)*log10(max(1d0,abs(val)))               
            case (p_log_abs_eps_nuc)
               val = safe_log10(abs(s% eps_nuc(k)))
            case (p_d_epsnuc_dlnd)
               val = s% d_epsnuc_dlnd(k)
            case (p_d_lnepsnuc_dlnd)
               val = s% d_epsnuc_dlnd(k)/max(1d0,abs(s% eps_nuc(k)))
            case (p_d_epsnuc_dlnT)
               val = s% d_epsnuc_dlnT(k)
            case (p_d_lnepsnuc_dlnT)
               val = s% d_epsnuc_dlnT(k)/max(1d0,abs(s% eps_nuc(k)))

            case (p_deps_dlnd_face)
               val = interp_val_to_pt(s% d_epsnuc_dlnd,k,nz,s% dq,'p_deps_dlnd_face')
            case (p_deps_dlnT_face)
               val = interp_val_to_pt(s% d_epsnuc_dlnT,k,nz,s% dq,'p_deps_dlnT_face')
            case (p_eps_nuc_neu_total)
               val = s% eps_nuc_neu_total(k)
            case (p_non_nuc_neu)
               val = s% non_nuc_neu(k)
            case (p_nonnucneu_plas)
               val = s% nonnucneu_plas(k)
            case (p_nonnucneu_brem)
               val = s% nonnucneu_brem(k)
            case (p_nonnucneu_phot)
               val = s% nonnucneu_phot(k)
            case (p_nonnucneu_pair)
               val = s% nonnucneu_pair(k)
            case (p_nonnucneu_reco)
               val = s% nonnucneu_reco(k)

            case (p_log_irradiation_heat)
               val = safe_log10(s% irradiation_heat(k))
            case (p_cgrav_factor)
               val = s% cgrav(k)/standard_cgrav
            case (p_alpha_mlt)
               val = s% alpha_mlt(k)

            case (p_extra_jdot)
               val = s% extra_jdot(k)
            case (p_extra_omegadot)
               val = s% extra_omegadot(k)
            case (p_extra_heat)
               val = s% extra_heat(k)
            case (p_extra_grav)
               val = s% extra_grav(k)
            case (p_extra_L)
               val = dot_product(s% dm(k:s% nz),s% extra_heat(k:s% nz))/Lsun
            case (p_log_extra_L)
               val = safe_log10( &
                  dot_product(s% dm(k:s% nz),s% extra_heat(k:s% nz))/Lsun)

            case (p_log_abs_eps_grav_dm_div_L)
               val = safe_log10( &
                  abs(s% eps_grav(k))*s% dm(k)/max(1d0,abs(s% L(k))))

            case (p_eps_grav_composition_term)
               if (s% include_composition_in_eps_grav) &
                  val = s% eps_grav_composition_term(k)
                  
            case (p_eps_grav_plus_eps_mdot)
               val = s% eps_grav(k) + s% eps_mdot(k)
            case (p_ergs_eps_grav_plus_eps_mdot)
               val = (s% eps_grav(k) + s% eps_mdot(k))*s% dm(k)*s% dt
                  
            case (p_eps_mdot)
               val = s% eps_mdot(k)
            case (p_ergs_mdot)
               val = s% eps_mdot(k)*s% dm(k)*s% dt

            case (p_div_v)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     Ap1 = 4*pi*s% R_center*s% R_center
                  else
                     vp1 = s% v(k+1)
                     Ap1 = 4*pi*s% r(k+1)*s% r(k+1)
                  end if
                  val = (4*pi*s% r(k)*s% r(k)*s% v(k) - Ap1*vp1)*s% rho(k)/s% dm(k)
               end if

            case (p_d_v_div_r_dm)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     rp1 = s% R_center
                  else
                     vp1 = s% v(k+1)
                     rp1 = s% r(k+1)
                  end if
                  v00 = s% v(k)
                  r00 = s% r(k)
                  if (rp1 > 0) then
                     val = (v00/r00 - vp1/rp1)/s% dm(k)
                  end if
               end if

            case (p_d_v_div_r_dr)
               if (s% v_flag) then
                  if (k == s% nz) then
                     vp1 = s% V_center
                     rp1 = s% R_center
                  else
                     vp1 = s% v(k+1)
                     rp1 = s% r(k+1)
                  end if
                  v00 = s% v(k)
                  r00 = s% r(k)
                  if (rp1 > 0) then
                     val = 4*pi*s% rmid(k)*s% rmid(k)*s% rho(k)* &
                           (v00/r00 - vp1/rp1)/s% dm(k)
                  end if
               end if

            case (p_rho_times_r3)
               val = s% rho_face(k)*s% r(k)*s% r(k)*s% r(k)
            case (p_log_rho_times_r3)
               val = safe_log10(s% rho_face(k)*s% r(k)*s% r(k)*s% r(k))

            case(p_du)
               if (s% u_flag) then
                  if (k == s% nz) then
                     val = s% u(k)
                  else
                     val = s% u(k) - s% u(k+1)
                  end if
               end if

            case(p_P_face)
               if (s% u_flag) val = s% P_face(k)
            case(p_log_P_face)
               if (s% u_flag) val = safe_log10(s% P_face(k))

            case (p_hse_ratio)
               if (k > 1 .and. k < nz .and. s% cgrav(k) > 0d0) then
                  val = (s% P(k-1) - s% P(k))/(-s% cgrav(k)*s% m(k)*s% dm_bar(k)/(4d0*pi*pow4(s% r(k)))) - 1d0
               end if
            case (p_hse_ratio_gyre)
               if (k > 1 .and. k < nz .and. s% cgrav(k) > 0d0) then
                  Pbar_00 = (s% P(k-1)*s% dm(k) + s% P(k)*s% dm(k-1))/(s% dm(k) + s% dm(k-1))
                  Pbar_p1 = (s% P(k)*s% dm(k+1) + s% P(k+1)*s% dm(k))/(s% dm(k+1) + s% dm(k))
                  val = (Pbar_00 - Pbar_p1)/(-0.5d0*s% dm(k)*( &
                     s% cgrav(k)*s% m(k)/(4d0*pi*pow4(s% r(k))) + &
                     s% cgrav(k+1)*s% m(k+1)/(4d0*pi*pow4(s% r(k+1))))) - 1d0
               end if

            case (p_dPdr_div_grav)
               if (k > 1 .and. k < nz .and. s% cgrav(k) > 0d0 .and. s% RTI_flag) then
                  val = s% dPdr_info(k)/s% rho_face(k)
               end if

            case (p_gradP_div_rho)
               if (k > 1) val = 4*pi*s% r(k)*s% r(k)*(s% P(k-1) - s% P(k))/s% dm_bar(k)
            case (p_dlnP_dlnR)
               if (k > 1) val = log(s% P_face(k-1)/s% P_face(k)) / (s% lnR(k-1) - s% lnR(k))
            case (p_dlnRho_dlnR)
               if (k > 1) val = log(s% rho_face(k-1)/s% rho_face(k)) / (s% lnR(k-1) - s% lnR(k))

            case (p_dvdt_grav)
               val = -s% cgrav(k)*s% m(k)/(s% r(k)*s% r(k))
            case (p_dvdt_dPdm)
               if (k > 1) val = -4*pi*s% r(k)*s% r(k)*(s% P(k-1) - s% P(k))/s% dm_bar(k)

            case (p_dm_eps_grav)
               val = s% eps_grav(k)*s% dm(k)
            case (p_eps_grav)
               val = s% eps_grav(k)
               
            case (p_log_xm_div_delta_m)
               val = safe_log10((s% m(1) - s% m(k))/abs(s% dt*s% mstar_dot))
            case (p_xm_div_delta_m)
               val = (s% m(1) - s% m(k))/abs(s% dt*s% mstar_dot)

            case (p_env_eps_grav)
               val = -s% gradT_sub_grada(k)*s% grav(k)*s% mstar_dot*s% Cp(k)*s% T(k) / &
                        (4*pi*s% r(k)*s% r(k)*s% P(k))

            case (p_mlt_mixing_type)
               int_val = s% mlt_mixing_type(k)
               val = dble(int_val)
               int_flag = .true.
            case (p_mlt_mixing_length)
               val = s% mlt_mixing_length(k)
            case (p_mlt_Gamma)
               val = s% mlt_Gamma(k)
            case (p_mlt_Zeta)
               if (abs(s% gradr(k) - s% grada_face(k)) > 1d-20) &
                  val = (s% gradr(k) - s% gradT(k))/(s% gradr(k) - s% grada_face(k))
            case (p_mlt_Pturb)
               if (k < s% nz) &
                  val = s% mlt_Pturb_factor*s% rho(k)/3d0*(s% mlt_vc_start(k)**2 + s% mlt_vc_start(k+1)**2)/2d0

            case (p_grad_density)
               val = s% grad_density(k)
            case (p_grad_temperature)
               val = s% grad_temperature(k)

            case (p_gradL_sub_gradr)
               val = s% gradL(k) - s% gradr(k)
            case (p_grada_sub_gradr)
               val = s% grada_face(k) - s% gradr(k)

            case (p_gradL)
               val = s% gradL(k)
            case (p_sch_stable)
               if (s% grada(k) > s% gradr(k)) val = 1
            case (p_ledoux_stable)
               if (s% gradL(k) > s% gradr(k)) val = 1

            case (p_eps_nuc_start)
               val = s% eps_nuc_start(k)

            case (p_dominant_isoA_for_thermohaline)
               int_val = chem_isos% Z_plus_N(s% dominant_iso_for_thermohaline(k))
               int_flag = .true.
            case (p_dominant_isoZ_for_thermohaline)
               int_val = chem_isos% Z(s% dominant_iso_for_thermohaline(k))
               int_flag = .true.
            case (p_gradL_composition_term)
               val = s% gradL_composition_term(k)

            case (p_log_D_conv)
               if (s% mixing_type(k) == convective_mixing) then
                  val = safe_log10(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_leftover)
               if (s% mixing_type(k) == leftover_convective_mixing) then
                  val = safe_log10(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_semi)
               if (s% mixing_type(k) == semiconvective_mixing) then
                  val = safe_log10(s% D_mix_non_rotation(k))
               else
                  if (s% conv_vel_flag .and. s% conv_vel_ignore_semiconvection) then
                     val = safe_log10(s% mlt_D_semi(k))
                  else
                     val = -99
                  end if
               end if
            case (p_log_D_ovr)
               if (s% mixing_type(k) == overshoot_mixing) then
                  val = safe_log10(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_rayleigh_taylor)
               if(s% RTI_flag) then
                  val = safe_log10(s% eta_RTI(k))
               else
                  val =-99
               end if
            case (p_log_D_anon)
               if (s% mixing_type(k) == anonymous_mixing) then
                  val = safe_log10(s% D_mix_non_rotation(k))
               else
                  val = -99
               end if
            case (p_log_D_thrm)
               if (s% mixing_type(k) == thermohaline_mixing) then
                  val = safe_log10(s% D_mix_non_rotation(k))
               else
                  if (s% conv_vel_flag .and. s% conv_vel_ignore_thermohaline) then
                     val = safe_log10(s% mlt_D_thrm(k))
                  else
                     val = -99
                  end if
               end if

            case (p_log_D_minimum)
               if (s% mixing_type(k) == minimum_mixing) then
                  val = safe_log10(s% D_mix(k))
               else
                  val = -99
               end if
               
            case (p_log_lambda_RTI_div_Hrho)
               if (s% RTI_flag) val = safe_log10( &
                  sqrt(s% alpha_RTI(k))*s% r(k)/s% rho(k)*abs(s% dRhodr_info(k)))
            case (p_lambda_RTI)
               if (s% RTI_flag) val = sqrt(s% alpha_RTI(k))*s% r(k)
            case (p_dPdr_info)
               if (s% RTI_flag) val = s% dPdr_info(k)
            case (p_dRhodr_info)
               if (s% RTI_flag) val = s% dRhodr_info(k)
               
            case (p_source_plus_alpha_RTI)
               if (s% RTI_flag) val = s% source_plus_alpha_RTI(k)
            case (p_log_source_plus_alpha_RTI)
               if (s% RTI_flag) val = safe_log10(s% source_plus_alpha_RTI(k))
            case (p_log_source_RTI)
               if (s% RTI_flag) val = safe_log10(s% source_plus_alpha_RTI(k))
            case (p_source_minus_alpha_RTI)
               if (s% RTI_flag) val = s% source_minus_alpha_RTI(k)
            case (p_log_source_minus_alpha_RTI)
               if (s% RTI_flag) val = safe_log10(abs(s% source_minus_alpha_RTI(k)))

            case (p_dudt_RTI)
               if (s% RTI_flag) val = s% dudt_RTI(k)
            case (p_dedt_RTI)
               if (s% RTI_flag) val = s% dedt_RTI(k)

            case (p_eta_RTI)
               if (s% RTI_flag) val = s% eta_RTI(k)
            case (p_log_eta_RTI)
               if (s% RTI_flag) val = safe_log10(abs(s% eta_RTI(k)))
            case (p_boost_for_eta_RTI)
               if (s% RTI_flag) val = s% boost_for_eta_RTI(k)
            case (p_log_boost_for_eta_RTI)
               if (s% RTI_flag) val = safe_log10(abs(s% boost_for_eta_RTI(k)))

            case (p_alpha_RTI)
               if (s% RTI_flag) val = s% alpha_RTI(k)
            case (p_log_alpha_RTI)
               if (s% RTI_flag) val = safe_log10(s% alpha_RTI(k))
            case (p_log_etamid_RTI)
               if (s% RTI_flag) val = safe_log10(s% etamid_RTI(k))

            case (p_log_sig_RTI)
               if (s% RTI_flag) val = safe_log10(s% sig_RTI(k))
            case (p_log_sigmid_RTI)
               if (s% RTI_flag) val = safe_log10(s% sigmid_RTI(k))

            case (p_log_D_omega)
               if (s% D_omega_flag) val = safe_log10(s% D_omega(k))
               
            case (p_log_D_mix_non_rotation)
               val = safe_log10(s% D_mix_non_rotation(k))
            case (p_log_D_mix_rotation)
               val = safe_log10(s% D_mix(k) - s% D_mix_non_rotation(k))
            case (p_log_D_mix)
               val = safe_log10(s% D_mix(k))
            case (p_log_sig_mix)
               val = safe_log10(s% sig(k))
            case (p_log_sig_raw_mix)
               val = safe_log10(s% sig_raw(k))

            case (p_d_gradT_dlnd00)
               val = s% d_gradT_dlnd00(k)
            case (p_d_gradT_dlnT00)
               val = s% d_gradT_dlnT00(k)
            case (p_d_gradT_dlndm1)
               val = s% d_gradT_dlndm1(k)
            case (p_d_gradT_dlnTm1)
               val = s% d_gradT_dlnTm1(k)
            case (p_d_gradT_dlnR)
               val = s% d_gradT_dlnR(k)
            case (p_d_gradT_dL)
               val = s% d_gradT_dL(k)
            case (p_d_gradT_dln_cvpv0)
               val = s% d_gradT_dln_cvpv0(k)

            case (p_burn_avg_epsnuc)
               if (s% op_split_burn) val = s% burn_avg_epsnuc(k)
            case (p_log_burn_avg_epsnuc)
               if (s% op_split_burn) &
                  val = safe_log10(abs(s% burn_avg_epsnuc(k)))
            case (p_burn_num_iters)
               if (s% op_split_burn) then
                  int_val = s% burn_num_iters(k); val = dble(int_val)
               else
                  int_val = 0; val = 0
               end if
               int_flag = .true.

            case (p_conv_vel_div_mlt_vc)
               if (s% mlt_vc(k) > 0d0) val = s% conv_vel(k)/s% mlt_vc(k)
               
            case (p_conv_vel)
               val = s% conv_vel(k)
            case (p_dt_times_conv_vel_div_mixing_length)
               val = s% dt*s% conv_vel(k)/s% mlt_mixing_length(k)
            case (p_log_dt_times_conv_vel_div_mixing_length)
               val = safe_log10(s% dt*s% conv_vel(k)/s% mlt_mixing_length(k))
            case (p_log_conv_vel)
               val = safe_log10(s% conv_vel(k))
            case (p_conv_vel_div_L_vel)
               val = s% conv_vel(k)/max(1d0,get_L_vel(k))
            case (p_conv_vel_div_csound)
               val = s% conv_vel(k)/s% csound(k)
            case (p_log_tau_conv_yrs)
               if (s% conv_vel(k) > 1d-99) then
                  val = safe_log10(s% mlt_mixing_length(k)/(4*s% conv_vel(k)*secyer))
               else
                  val = -99
               end if
            case (p_mixing_type)
               val = dble(s% mixing_type(k))
               int_val = s% mixing_type(k)
               int_flag = .true.
            case (p_conv_mixing_type) ! OBSOLETE
               val = dble(s% mixing_type(k))
               int_val = s% mixing_type(k)
               int_flag = .true.
            case (p_log_mlt_D_mix)
               val = safe_log10(s% mlt_D(k))
            case (p_log_t_thermal)
               val = safe_log10(s% Cp(k)*s% T(k)*(s% m(1) - s% m(k))/s% L(k))
            case (p_log_cp_T_div_t_sound)
               val = safe_log10( &
                  s% Cp(k)*s% T(k)/(s% P(k)/(s% rho(k)*s% grav(k))/s% csound(k)))
            case (p_log_t_sound)
               val = safe_log10(s% P(k)/(s% rho(k)*s% grav(k))/s% csound(k))
            case (p_pressure_scale_height)
               val = s% P(k)/(s% rho(k)*s% grav(k))/Rsun
            case (p_pressure_scale_height_cm)
               val = s% P(k)/(s% rho(k)*s% grav(k))
            case (p_actual_gradT)
               val = s% actual_gradT(k)
            case (p_gradT_sub_actual_gradT)
               val = s% gradT(k) - s% actual_gradT(k)
            case (p_gradT)
               val = s% gradT(k)
            case (p_gradr)
               val = s% gradr(k)
            case (p_grada_sub_actual_gradT)
               val = s% grada_face(k) - s% actual_gradT(k)
            case (p_grada_sub_gradT)
               val = s% grada_face(k) - s% gradT(k)

            case (p_omega)
               val = if_rot(s% omega,k)

            case (p_log_omega)
               val = safe_log10(if_rot(s% omega,k))
            case (p_log_j_rot)
               val = safe_log10(if_rot(s% j_rot,k))
            case (p_log_J_inside)
               if (s% rotation_flag) then
                  val = safe_log10(dot_product(s% j_rot(k:s% nz), s% dm(k:s% nz)))
               else
                  val = -99.0d0
               end if
            case (p_log_J_div_M53)
               if (s% rotation_flag) then
                  val = safe_log10(&
                       dot_product(s% j_rot(k:s% nz), s% dm(k:s% nz)) * &
                       1d-50/pow(s% m(k)/Msun,5d0/3d0))
               else
                  val = -99.0d0
               end if

            case (p_shear)
               val = if_rot(s% omega_shear,k)
            case (p_log_abs_shear)
               if (s% rotation_flag) then
                  val = safe_log10(s% omega_shear(k))
                  if (is_bad(val)) then
                     write(*,2) 'val', k, val
                     write(*,2) 's% omega_shear(k)', k, s% omega_shear(k)
                     stop 'profile'
                  end if
               else
                  val = -99
               end if
            case (p_log_abs_dlnR_domega)
               if (s% rotation_flag) then
                  val = -safe_log10(s% omega_shear(k))
                  if (is_bad(val)) then
                     write(*,2) 'val', k, val
                     write(*,2) 's% omega_shear(k)', k, s% omega_shear(k)
                     stop 'profile'
                  end if
               else
                  val = -99
               end if
            case (p_i_rot)
               val = if_rot(s% i_rot,k)
            case (p_j_rot)
               val = if_rot(s% j_rot,k)
            case (p_v_rot)
               val = if_rot(s% omega,k)*if_rot(s% r_equatorial,k)*1d-5 ! km/sec
            case (p_fp_rot)
               val = if_rot(s% fp_rot,k, alt=1.0d0)
            case (p_ft_rot)
               val = if_rot(s% ft_rot,k, alt=1.0d0)
            case (p_ft_rot_div_fp_rot)
               if(s% rotation_flag) then
                  val = s% ft_rot(k)/s% fp_rot(k) 
               else
                  val = 1.0d0
               end if
            case (p_w_div_w_crit_roche)
               val = if_rot(s% w_div_w_crit_roche,k)
            case (p_w_div_w_crit_roche2)
               val = if_rot(s% xh(s% i_w_div_wc,:),k)
            case (p_log_am_nu_non_rot)
               val = safe_log10(if_rot(s% am_nu_non_rot,k))
            case (p_log_am_nu_rot)
               val = safe_log10(if_rot(s% am_nu_rot,k))
            case (p_log_am_nu)
               val = safe_log10(if_rot(s% am_nu_rot,k) + if_rot(s% am_nu_non_rot,k))

            case (p_r_polar)
               val = if_rot(s% r_polar,k, alt=s% r(k))/Rsun
            case (p_log_r_polar)
               val = safe_log10(if_rot(s% r_polar,k, alt=s% r(k))/Rsun)
            case (p_r_equatorial)
               val = if_rot(s% r_equatorial,k, alt=s% r(k))/Rsun
            case (p_log_r_equatorial)
               val = safe_log10(if_rot(s% r_equatorial,k, alt=s% r(k))/Rsun)
            case (p_r_e_div_r_p)
               if (s% rotation_flag) then
                   if(s% r_polar(k) > 1) val = s% r_equatorial(k)/s% r_polar(k)
               end if
            case (p_omega_crit)
               val = omega_crit(s,k)
            case (p_omega_div_omega_crit)
               if (s% rotation_flag) then
                  val = omega_crit(s,k)
                  if (val < 1d-50) then
                     val = 0
                  else
                     val = s% omega(k)/val
                  end if
               end if

            case (p_eps_WD_sedimentation)
               if (s% do_element_diffusion) val = s% eps_WD_sedimentation(k)
            case (p_log_eps_WD_sedimentation)
               if (s% do_element_diffusion) val = safe_log10(s% eps_WD_sedimentation(k))

            case (p_eps_diffusion)
               if (s% do_element_diffusion) val = s% eps_diffusion(k)
            case (p_log_eps_diffusion)
               if (s% do_element_diffusion) val = safe_log10(s% eps_diffusion(k))

            case (p_e_field)
               if (s% do_element_diffusion) val = s% E_field(k)
            case (p_log_e_field)
               if (s% do_element_diffusion) val = safe_log10(s% E_field(k))

            case (p_g_field_element_diffusion)
               if (s% do_element_diffusion) val = s% g_field_element_diffusion(k)
            case (p_log_g_field_element_diffusion)
               if (s% do_element_diffusion) &
                  val = safe_log10(s% g_field_element_diffusion(k))

            case (p_eE_div_mg_element_diffusion)
               if (s% do_element_diffusion) then
                  if ( s% g_field_element_diffusion(k) /= 0d0) then
                     val = qe * s% E_field(k)/(amu * s% g_field_element_diffusion(k))
                  else
                     val = 0d0
                  end if
               end if
            case (p_log_eE_div_mg_element_diffusion)
               if (s% do_element_diffusion) &
                  val = safe_log10(qe * s% E_field(k)/(amu * s% g_field_element_diffusion(k)))

            case (p_richardson_number)
               val = if_rot(s% richardson_number,k)
            case (p_am_domega_dlnR)
               val = if_rot(s% domega_dlnR,k)

            case (p_am_log_sig) ! == am_log_sig_omega
               val = safe_log10(if_rot(s% am_sig_omega,k))
            case (p_am_log_sig_omega)
               val = safe_log10(if_rot(s% am_sig_omega,k))
            case (p_am_log_sig_j)
               val = safe_log10(if_rot(s% am_sig_j,k))

            case (p_am_log_nu_omega)
               val = safe_log10(if_rot(s% am_nu_omega,k))
            case (p_am_log_nu_j)
               val = safe_log10(if_rot(s% am_nu_j,k))

            case (p_am_log_nu_rot)
               val = safe_log10(if_rot(s% am_nu_rot,k))
            case (p_am_log_nu_non_rot)
               val = safe_log10(if_rot(s% am_nu_non_rot,k))

            case (p_am_log_D_visc)
               if (s% am_nu_visc_factor >= 0) then
                  f = s% am_nu_visc_factor
               else
                  f = s% D_visc_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_visc,k))
            case (p_am_log_D_DSI)
               if (s% am_nu_DSI_factor >= 0) then
                  f = s% am_nu_DSI_factor
               else
                  f = s% D_DSI_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_DSI,k))
            case (p_am_log_D_SH)
               if (s% am_nu_SH_factor >= 0) then
                  f = s% am_nu_SH_factor
               else
                  f = s% D_SH_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_SH,k))
            case (p_am_log_D_SSI)
               if (s% am_nu_SSI_factor >= 0) then
                  f = s% am_nu_SSI_factor
               else
                  f = s% D_SSI_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_SSI,k))

            case (p_am_log_D_ES)
               if (s% am_nu_ES_factor >= 0) then
                  f = s% am_nu_ES_factor
               else
                  f = s% D_ES_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_ES,k))
            case (p_am_log_D_GSF)
               if (s% am_nu_GSF_factor >= 0) then
                  f = s% am_nu_GSF_factor
               else
                  f = s% D_GSF_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_GSF,k))
            case (p_am_log_D_ST)
               if (s% am_nu_ST_factor >= 0) then
                  f = s% am_nu_ST_factor
               else
                  f = s% D_ST_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% D_ST,k))
            case (p_am_log_nu_ST)
               if (s% am_nu_ST_factor >= 0) then
                  f = s% am_nu_ST_factor
               else
                  f = s% D_ST_factor
               end if
               val = safe_log10(am_nu_factor*f*if_rot(s% nu_ST,k))

            case (p_dynamo_log_B_r)
               val = safe_log10(if_rot(s% dynamo_B_r,k))
            case (p_dynamo_log_B_phi)
               val = safe_log10(if_rot(s% dynamo_B_phi,k))

            case (p_grada_face)
               val = s% grada_face(k)
            case (p_gradr_div_grada)
               val = s% gradr(k)/s% grada_face(k)
            case (p_gradr_sub_grada)
               val = s% gradr(k) - s% grada_face(k)
            case (p_gradT_sub_a)
               val = s% gradT(k) - s% grada_face(k)
            case (p_gradT_sub_grada)
               val = s% gradT(k) - s% grada_face(k)
            case (p_grad_superad)
               if (k > 1) val = s% grad_superad(k)
            case (p_grad_superad_actual)
               if (k > 1) val = s% grad_superad_actual(k)
            case (p_gradT_div_grada)
               val = s% gradT(k) / s% grada_face(k)
            case (p_gradr_sub_gradT)
               val = s% gradr(k) - s% gradT(k)
            case (p_gradT_sub_gradr)
               val = s% gradT(k) - s% gradr(k)
               
            case (p_gradT_rel_err)
               if (k > 1) then
                  val = (s% lnT(k-1) - s% lnT(k))/(s% lnP(k-1) - s% lnP(k))
                  val = (s% gradT(k) - val)/s% gradT(k)
               end if

            case (p_gradT_div_gradr)
               if (abs(s% gradr(k)) < 1d-99) then
                  val = 1d0
               else
                  val = s% gradT(k) / s% gradr(k)
               end if
            case (p_log_gradT_div_gradr)
               if (abs(s% gradr(k)) < 1d-99) then
                  val = 0d0
               else
                  val = safe_log10(s% gradT(k) / s% gradr(k))
               end if

            case (p_log_mlt_Gamma)
               val = safe_log10(s% mlt_Gamma(k))

            case (p_log_mlt_vc)
               val = safe_log10(s% mlt_vc(k))
            case (p_mlt_vc)
               val = s% mlt_vc(k)
            case (p_delta_r)
               val = s% r(k) - s% r_start(k)
            case (p_delta_L)
               val = s% L(k) - s% L_start(k)
            case (p_delta_cell_vol)
               if (k == s% nz) then
                  rp1 = s% R_center
                  rp1_start = s% R_center_old
               else
                  rp1 = s% r(k+1)
                  rp1_start = s% r_start(k+1)
               end if
               r00 = s% r(k)
               r00_start = s% r_start(k)
               dr3 = r00*r00*r00 - rp1*rp1*rp1
               dr3_start = r00_start*r00_start*r00_start - rp1_start*rp1_start*rp1_start
               val = 4d0/3d0*pi*(dr3 - dr3_start)
            case (p_delta_entropy)
               val = s% entropy(k) - exp(s% lnS_start(k))/(avo*kerg)
            case (p_delta_T)
               val = s% T(k) - s% T_start(k)
            case (p_delta_rho)
               val = s% rho(k) - exp(s% lnd_start(k))
            case (p_delta_eps_nuc)
               val = s% eps_nuc(k) - s% eps_nuc_start(k)
            case (p_delta_mu)
               val = s% mu(k) - s% mu_start(k)

            case (p_super_ad)
               val = max(0d0, s% gradT(k) - s% grada_face(k))

            case (p_accel_div_grav)
               if (s% v_flag) val = s% dv_dt(k)/s% grav(k)

            case (p_dlnd_dt_const_q)
               val = s% dlnd_dt_const_q(k)
            case (p_dlnT_dt_const_q)
               val = s% dlnT_dt_const_q(k)

            case (p_dlnd)
               val = s% dt*s% dlnd_dt(k)
            case (p_dlnT)
               val = s% dt*s% dlnT_dt(k)
            case (p_dlnR)
               val = s% dt*s% dlnR_dt(k)

            case (p_dlnd_dt)
               val = s% dlnd_dt(k)
            case (p_dlnT_dt)
               val = s% dlnT_dt(k)
            case (p_dlnR_dt)
               val = s% dlnR_dt(k)
            case (p_dr_dt)
               val = s% dlnR_dt(k)*s% r(k)
            case (p_dv_dt)
               if (s% v_flag) val = s% dv_dt(k)
            case (p_du_dt)
               if (s% u_flag) val = s% du_dt(k)

            case (p_del_entropy)
               val = s% entropy(k) - exp(s% lnS_start(k))/(avo*kerg)
            case (p_ds_from_eps_grav)
               val = -s% dt*s% eps_grav(k)/s% T(k)/(avo*kerg)

            case(p_dt_dm_eps_grav)
               val = s% eps_grav(k)*s% dm(k)*s% dt

            case(p_dm_de)
               val = s% dm(k)*(s% energy(k) - s% energy_start(k))
            case(p_dt_dL)
               if (k < s% nz) val = s% dt*(s% L(k) - s% L(k+1))

            case (p_signed_dlnd)
               val = 1d6*s% dlnd_dt(k)*s% dt
               val = sign(1d0,val)*log10(max(1d0,abs(val)))
            case (p_signed_dlnT)
               val = 1d6*s% dlnT_dt(k)*s% dt
               val = sign(1d0,val)*log10(max(1d0,abs(val)))
            case (p_cno_div_z)
               cno = s% xa(s% net_iso(ic12),k) + &
                     s% xa(s% net_iso(in14),k) + s% xa(s% net_iso(io16),k)
               z = 1 - (s% xa(s% net_iso(ih1),k) + s% xa(s% net_iso(ihe4),k))
               if (z > 1d-50) then
                  val = cno/z
               else
                  val = 0
               end if
            case (p_dE)
               val = s% energy(k) - s% energy_start(k)
            case (p_dr)
               if (k < s% nz) then
                  val = s% r(k) - s% r(k+1)
               else
                  val = s% r(k) - s% R_center
               end if
            case (p_dr_ratio)
               if (k == 1 .or. k == s% nz) then
                  val = 1
               else
                  val = (s% r(k-1) - s% r(k))/(s% r(k) - s% r(k+1))
               end if
            case (p_dv)
               if (.not. s% v_flag) then
                  val = 0
               else if (k < s% nz) then
                  val = s% v(k+1) - s% v(k)
               else
                  val = -s% v(k)
               end if
            case (p_dt_dv_div_dr)
               if (.not. s% v_flag) then
                  val = 0
               else if (k < s% nz) then
                  val = s% dt*(s% v(k+1) - s% v(k))/(s% r(k) - s% r(k+1))
               else
                  val = -s% dt*s% v(k)/s% r(k)
               end if

            case (p_dlog_h1_dlogP)
               val = get_dlogX_dlogP(ih1, k)
            case (p_dlog_he3_dlogP)
               val = get_dlogX_dlogP(ihe3, k)
            case (p_dlog_he4_dlogP)
               val = get_dlogX_dlogP(ihe4, k)
            case (p_dlog_c12_dlogP)
               val = get_dlogX_dlogP(ic12, k)
            case (p_dlog_c13_dlogP)
               val = get_dlogX_dlogP(ic13, k)
            case (p_dlog_n14_dlogP)
               val = get_dlogX_dlogP(in14, k)
            case (p_dlog_o16_dlogP)
               val = get_dlogX_dlogP(io16, k)
            case (p_dlog_ne20_dlogP)
               val = get_dlogX_dlogP(ine20, k)
            case (p_dlog_mg24_dlogP)
               val = get_dlogX_dlogP(img24, k)
            case (p_dlog_si28_dlogP)
               val = get_dlogX_dlogP(isi28, k)

            case (p_dlog_pp_dlogP)
               val = get_dlog_eps_dlogP(ipp, k)
            case (p_dlog_cno_dlogP)
               val = get_dlog_eps_dlogP(icno, k)
            case (p_dlog_3alf_dlogP)
               val = get_dlog_eps_dlogP(i3alf, k)

            case (p_dlog_burn_c_dlogP)
               val = get_dlog_eps_dlogP(i_burn_c, k)
            case (p_dlog_burn_n_dlogP)
               val = get_dlog_eps_dlogP(i_burn_n, k)
            case (p_dlog_burn_o_dlogP)
               val = get_dlog_eps_dlogP(i_burn_o, k)

            case (p_dlog_burn_ne_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ne, k)
            case (p_dlog_burn_na_dlogP)
               val = get_dlog_eps_dlogP(i_burn_na, k)
            case (p_dlog_burn_mg_dlogP)
               val = get_dlog_eps_dlogP(i_burn_mg, k)

            case (p_dlog_cc_dlogP)
               val = get_dlog_eps_dlogP(icc, k)
            case (p_dlog_co_dlogP)
               val = get_dlog_eps_dlogP(ico, k)
            case (p_dlog_oo_dlogP)
               val = get_dlog_eps_dlogP(ioo, k)

            case (p_dlog_burn_si_dlogP)
               val = get_dlog_eps_dlogP(i_burn_si, k)
            case (p_dlog_burn_s_dlogP)
               val = get_dlog_eps_dlogP(i_burn_s, k)
            case (p_dlog_burn_ar_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ar, k)
            case (p_dlog_burn_ca_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ca, k)
            case (p_dlog_burn_ti_dlogP)
               val = get_dlog_eps_dlogP(i_burn_ti, k)
            case (p_dlog_burn_cr_dlogP)
               val = get_dlog_eps_dlogP(i_burn_cr, k)
            case (p_dlog_burn_fe_dlogP)
               val = get_dlog_eps_dlogP(i_burn_fe, k)
            case (p_dlog_pnhe4_dlogP)
               val = get_dlog_eps_dlogP(ipnhe4, k)
            case (p_dlog_photo_dlogP)
               val = get_dlog_eps_dlogP(iphoto, k)
            case (p_dlog_other_dlogP)
               val = get_dlog_eps_dlogP(iother, k)

            case(p_d_u_div_rmid)
               if (s% u_flag .and. k > 1) &
                  val = s% u(k-1)/s% rmid(k-1) - s% u(k)/s% rmid(k)
            case(p_d_u_div_rmid_start)
               if (s% u_flag .and. k > 1) &
                  val = s% u(k-1)/s% rmid_start(k-1) - s% u(k)/s% rmid_start(k)

            case(p_Pturb)
               if (s% Eturb_flag) then
                  val = s% Eturb(k)*s% rho(k)
               else if (s% RSP_flag) then
                  val = s% Et(k)*s% rho(k)
               end if
            case(p_log_Pturb)
               if (s% Eturb_flag) then
                  val = safe_log10(s% Eturb(k)*s% rho(k))
               else if (s% RSP_flag) then
                  val = safe_log10(s% Et(k)*s% rho(k))
               end if
            case(p_Eturb)
               if (s% Eturb_flag) then
                  val = s% Eturb(k)
               else if (s% RSP_flag) then
                  val = s% Et(k)
               end if               
            case(p_log_Eturb)
               if (s% Eturb_flag) then
                  val = safe_log10(s% Eturb(k))
               else if (s% RSP_flag) then
                  val = safe_log10(s% Et(k))
               end if
            case(p_avQ)
               if (s% use_avQ_art_visc .or. s% RSP_flag) val = s% avQ(k)
            case(p_Hp_face)
               if (s% RSP_flag) val = s% Hp_face(k)
            case(p_Y_face)
               if (s% RSP_flag) val = s% Y_face(k)
            case(p_PII_face)
               if (s% RSP_flag) val = s% PII(k)
            case(p_Chi)
               if (s% RSP_flag) val = s% Chi(k)
            case(p_COUPL)
               if (rsp_or_eturb) val = s% COUPL(k)
            case(p_SOURCE)
               if (rsp_or_eturb) val = s% SOURCE(k)
            case(p_DAMP)
               if (rsp_or_eturb) val = s% DAMP(k)
            case(p_DAMPR)
               if (rsp_or_eturb) val = s% DAMPR(k)
            case(p_Eq)
               if (rsp_or_eturb) val = s% Eq(k)
            case(p_Uq)
               if (rsp_or_eturb) val = s% Uq(k)
            case(p_Lr)
               if (rsp_or_eturb) val = s% Lr(k)
            case(p_Lr_div_L)
               if (rsp_or_eturb) val = s% Lr(k)/s% L(k)
            case(p_Lc)
               if (rsp_or_eturb) val = s% Lc(k)
            case(p_Lc_div_L)
               if (rsp_or_eturb) val = s% Lc(k)/s% L(k)
            case(p_Lt)
               if (rsp_or_eturb) val = s% Lt(k)
            case(p_Lt_div_L)
               if (rsp_or_eturb) val = s% Lt(k)/s% L(k)
               
            case(p_rsp_Et)
               if (s% rsp_flag) val = s% Et(k)
            case(p_rsp_logEt)
               if (s% rsp_flag) &
                  val = safe_log10(s% Et(k))
            case(p_rsp_vt)
               if (s% rsp_flag) val = sqrt2*s% w(k)
            case(p_rsp_vt_div_cs)
               if (s% rsp_flag) val = sqrt2*s% w(k)/s% csound(k)
            case(p_rsp_Pt)
               if (s% rsp_flag) val = s% Pt(k)
            case(p_rsp_Eq)
               if (s% rsp_flag) val = s% Eq(k)
            case(p_rsp_PII_face)
               if (s% rsp_flag) val = s% PII(k)
            case(p_rsp_src_snk)
               if (s% rsp_flag) val = s% COUPL(k)
            case(p_rsp_src)
               if (s% rsp_flag) val = s% SOURCE(k)
            case(p_rsp_sink)
               if (s% rsp_flag) val = s% DAMP(k) + s% DAMPR(k)
            case(p_rsp_damp)
               if (s% rsp_flag) val = s% DAMP(k)
            case(p_rsp_dampR)
               if (s% rsp_flag) val = s% DAMPR(k)
            case(p_rsp_Hp_face)
               if (s% rsp_flag) val = s% Hp_face(k)
            case(p_rsp_Chi)
               if (s% rsp_flag) val = s% Chi(k)
            case(p_rsp_avQ)
               if (s% rsp_flag) val = s% avQ(k)
            case(p_rsp_erad)
               if (s% rsp_flag) val = s% erad(k)
            case(p_rsp_log_erad)
               if (s% rsp_flag) val = safe_log10(s% erad(k))
            case(p_rsp_log_dt_div_heat_exchange_timescale)
               if (s% rsp_flag) val = safe_log10(s% dt*clight*s% opacity(k)*s% rho(k))
            case(p_rsp_heat_exchange_timescale)
               if (s% rsp_flag) val = 1d0/(clight*s% opacity(k)*s% rho(k))
            case(p_rsp_log_heat_exchange_timescale)
               if (s% rsp_flag) &
                  val = safe_log10(1d0/(clight*s% opacity(k)*s% rho(k)))
            case(p_rsp_Y_face)
               if (s% rsp_flag) then
                  if (k > 1) then
                     val = s% Y_face(k)
                  else ! for plotting, use value at k=2
                     val = s% Y_face(2)
                  end if
               end if
            case(p_rsp_gradT)
               if (s% rsp_flag) then
                  if (k > 1) then ! Y is superadiabatic gradient
                     val = s% Y_face(k) + 0.5d0*(s% grada(k-1) + s% grada(k))
                  else ! for plotting, use value at k=2
                     val = s% Y_face(2) + 0.5d0*(s% grada(1) + s% grada(2))
                  end if
               end if
            case(p_rsp_Uq)
               if (s% rsp_flag) then
                  if (k > 1) then
                     val = s% Uq(k)
                  else ! for plotting, use value at k=2
                     val = s% Uq(2)
                  end if
               end if               
            case(p_rsp_Lr)
               if (s% rsp_flag) val = s% Fr(k)*4d0*pi*s% r(k)*s% r(k)
            case(p_rsp_Lr_div_L)
               if (s% rsp_flag) val = s% Fr(k)*4d0*pi*s% r(k)*s% r(k)/s% L(k)
            case(p_rsp_Lc)
               if (s% rsp_flag) then
                  val = s% Lc(k)
                  if (k > 1) then
                     val = s% Lc(k)
                  else ! for plotting, use value at k=2
                     val = s% Lc(2)
                  end if
               end if
            case(p_rsp_Lc_div_L)
               if (s% rsp_flag) then
                  if (k > 1) then
                     val = s% Lc(k)/s% L(k)
                  else ! for plotting, use value at k=2
                     val = s% Lc(2)/s% L(2)
                  end if
               end if
            case(p_rsp_Lt)
               if (s% rsp_flag) then
                  if (k > 1) then
                     val = s% Lt(k)
                  else ! for plotting, use value at k=2
                     val = s% Lt(2)
                  end if
               end if
            case(p_rsp_Lt_div_L)
               if (s% rsp_flag) then
                  if (k > 1) then
                     val = s% Lt(k)/s% L(k)
                  else ! for plotting, use value at k=2
                     val = s% Lt(2)/s% L(2)
                  end if
               end if

            case(p_rsp_WORK)
               if (s% rsp_flag) val = rsp_WORK(s,k)
            case(p_rsp_WORKQ)
               if (s% rsp_flag) val = rsp_WORKQ(s,k)
            case(p_rsp_WORKT)
               if (s% rsp_flag) val = rsp_WORKT(s,k)
            case(p_rsp_WORKC)
               if (s% rsp_flag) val = rsp_WORKC(s,k)
               
            case(p_conv_vel_residual)
               val = s% conv_vel_residual(k)
            case(p_log_conv_vel_residual)
               val = safe_log10(abs(s% conv_vel_residual(k)))
            case(p_dconv_vel_dt)
               val = s% dln_cvpv0_dt(k)*(s% conv_vel(k) + s% conv_vel_v0)

            case (p_total_energy) ! specific total energy at k
               val = eval_cell_section_total_energy(s,k,k)/s% dm(k)               
            case (p_total_energy_sign) ! specific total energy at k
               val = eval_cell_section_total_energy(s,k,k)
               if (val > 0d0) then
                  int_val = 1
               else if (val < 0d0) then
                  int_val = -1
               else
                  int_val = 0
               end if  
               val = dble(int_val)
               int_flag = .true.
            case (p_total_energy_integral) ! from surface down to k
               val = s%total_energy_integral_surface(k)
            case (p_total_energy_integral_outward) ! from center up to k
               val = s%total_energy_integral_center(k)
               
            case (p_cell_specific_IE)
               val = s% energy(k)
            case (p_cell_ie_div_star_ie)
               val = s% energy(k)*s% dm(k)/s% total_internal_energy_end
            case (p_log_cell_specific_IE)
               val = safe_log10(s% energy(k))
            case (p_log_cell_ie_div_star_ie)
               val = safe_log10(s% energy(k)*s% dm(k)/s% total_internal_energy_end)

            case (p_cell_specific_PE)
               val = cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)

            case (p_cell_specific_KE)
               val = cell_specific_KE(s,k,d_dv00,d_dvp1)

            case (p_cell_IE_div_IE_plus_KE)
               val = s% energy(k)/(s% energy(k) + cell_specific_KE(s,k,d_dv00,d_dvp1))
               
            case (p_cell_KE_div_IE_plus_KE)
               f = cell_specific_KE(s,k,d_dv00,d_dvp1)
               val = f/(s% energy(k) + f)

            case (p_dlnX_dr)
               klo = max(1,k-1)
               khi = min(nz,k+1)
               val = log(max(1d-99,max(1d-99,s% X(klo))/max(1d-99,s% X(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
            case (p_dlnY_dr)
               klo = max(1,k-1)
               khi = min(nz,k+1)
               val = log(max(1d-99,max(1d-99,s% Y(klo))/max(1d-99,s% Y(khi))))  &
                              /  (s% rmid(klo) - s% rmid(khi))
            case (p_dlnRho_dr)
               klo = max(1,k-1)
               khi = min(nz,k+1)
               val = (s% lnd(klo) - s% lnd(khi))/(s% rmid(klo) - s% rmid(khi))

            case (p_brunt_B)
               if (s% calculate_Brunt_N2) val = s% brunt_B(k)
            case (p_brunt_nonB)
               if (s% calculate_Brunt_N2) val = -s% gradT_sub_grada(k)
            case (p_log_brunt_B)
               val = log10(max(1d-99,s% brunt_B(k)))
            case (p_log_brunt_nonB)
               if (s% calculate_Brunt_N2) val = log10(max(1d-99,-s% gradT_sub_grada(k)))

            case (p_brunt_N2)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)
            case (p_brunt_N2_composition_term)
               if (s% calculate_Brunt_N2) val = s% brunt_N2_composition_term(k)
            case (p_brunt_N2_structure_term)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k) - s% brunt_N2_composition_term(k)
            case (p_log_brunt_N2_composition_term)
               if (s% calculate_Brunt_N2) val = &
                  safe_log10(s% brunt_N2_composition_term(k))
            case (p_log_brunt_N2_structure_term)
               if (s% calculate_Brunt_N2) val = &
                  safe_log10(s% brunt_N2(k) - s% brunt_N2_composition_term(k))

            case (p_brunt_A)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)*s% r(k)/s% grav(k)
            case (p_brunt_A_div_x2)
               x = s% r(k)/s% r(1)
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)*s% r(k)/s% grav(k)/x/x
            case (p_log_brunt_N2_dimensionless)
               if (s% calculate_Brunt_N2) val = &
                  safe_log10(s% brunt_N2(k)/(3*s% cgrav(1)*s% m_grav(1)/pow3(s% r(1))))
            case (p_brunt_N2_dimensionless)
               if (s% calculate_Brunt_N2) val = &
                  s% brunt_N2(k)/(3*s% cgrav(1)*s% m_grav(1)/pow3(s% r(1)))
            case (p_brunt_N_dimensionless)
               if (s% calculate_Brunt_N2) val = &
                  sqrt(max(0d0,s% brunt_N2(k))/(3*s% cgrav(1)*s% m_grav(1)/pow3(s% r(1))))
            case (p_brunt_N)
               if (s% calculate_Brunt_N2) val = sqrt(max(0d0,s% brunt_N2(k)))
            case (p_brunt_frequency) ! cycles per day
               if (s% calculate_Brunt_N2) val = &
                  (24d0*60d0*60d0/(2*pi))*sqrt(max(0d0,s% brunt_N2(k)))
            case (p_log_brunt_N)
               if (s% calculate_Brunt_N2) val = safe_log10(sqrt(max(0d0,s% brunt_N2(k))))
            case (p_log_brunt_N2)
               if (s% calculate_Brunt_N2) val = safe_log10(s% brunt_N2(k))

            case (p_brunt_nu) ! micro Hz
               if (s% calculate_Brunt_N2) val = s% brunt_N2(k)
               val = (1d6/(2*pi))*sqrt(max(0d0,val))
            case (p_log_brunt_nu) ! micro Hz
               if (s% calculate_Brunt_N2) &
                  val = safe_log10((1d6/(2*pi))*sqrt(max(0d0,s% brunt_N2(k))))

            case (p_lamb_S)
               val = sqrt(2d0)*s% csound_face(k)/s% r(k) ! for l=1
            case (p_lamb_S2)
               val = 2d0*pow2(s% csound_face(k)/s% r(k)) ! for l=1

            case (p_lamb_Sl1)
               val = (1d6/(2*pi))*sqrt(2d0)*s% csound_face(k)/s% r(k) ! microHz
            case (p_lamb_Sl2)
               val = (1d6/(2*pi))*sqrt(6d0)*s% csound_face(k)/s% r(k) ! microHz
            case (p_lamb_Sl3)
               val = (1d6/(2*pi))*sqrt(12d0)*s% csound_face(k)/s% r(k) ! microHz
            case (p_lamb_Sl10)
               val = (1d6/(2*pi))*sqrt(110d0)*s% csound_face(k)/s% r(k) ! microHz

            case (p_log_lamb_Sl1)
               val = safe_log10((1d6/(2*pi))*sqrt(2d0)*s% csound_face(k)/s% r(k)) ! microHz
            case (p_log_lamb_Sl2)
               val = safe_log10((1d6/(2*pi))*sqrt(6d0)*s% csound_face(k)/s% r(k)) ! microHz
            case (p_log_lamb_Sl3)
               val = safe_log10((1d6/(2*pi))*sqrt(12d0)*s% csound_face(k)/s% r(k)) ! microHz
            case (p_log_lamb_Sl10)
               val = safe_log10((1d6/(2*pi))*sqrt(110d0)*s% csound_face(k)/s% r(k)) ! microHz

            case (p_brunt_N_div_r_integral)
               if (s% calculate_Brunt_N2) val = get_brunt_N_div_r_integral(k)
            case (p_sign_brunt_N2)
               if (s% calculate_Brunt_N2) val = sign(1d0,s% brunt_N2(k))

            case (p_k_r_integral)
               if (s% calculate_Brunt_N2) val = get_k_r_integral(k,1,1d0)

            case (p_brunt_N2_sub_omega2)
               if (s% calculate_Brunt_N2) then
                  val = s% brunt_N2(k) - pow2(2*pi*s% nu_max/1d6)
                  if (val > 0d0) then
                     val = 1
                  else
                     val = 0
                  end if
               end if
            case (p_sl2_sub_omega2)
               if (s% calculate_Brunt_N2) then
                  val = 2*pow2(s% csound_face(k)/s% r(k)) - pow2(2*pi*s% nu_max/1d6)
                  if (val >= 0d0) then
                     val = 1
                  else
                     val = 0
                  end if
               end if

            case (p_cs_at_cell_bdy)
               val = s% csound_face(k)
            case (p_log_mdot_cs) ! log10(4 Pi r^2 csound rho / (Msun/year))
               val = safe_log10(4*pi*s% r(k)*s% r(k)*s% csound(k)*s% rho(k)/(Msun/secyer))
            case (p_log_mdot_v) ! log10(4 Pi r^2 v rho / (Msun/year))
               if (s% u_flag) then
                  val = safe_log10(4*pi*s% r(k)*s% r(k)*s% u_face(k)*s% rho(k)/(Msun/secyer))
               else if (s% v_flag) then
                  val = safe_log10(4*pi*s% r(k)*s% r(k)*s% v(k)*s% rho(k)/(Msun/secyer))
               end if
            case (p_log_L_div_CpTMdot)
               if (s% star_mdot == 0) then
                  val = 0
               else
                  val = safe_log10(s% L(k)/(s% cp(k)*s% T(k)*abs(s% star_mdot)*(Msun/secyer)))
               end if
            case (p_logR_kap)
               val = s% lnd(k)/ln10 - 3d0*s% lnT(k)/ln10 + 18d0
            case (p_logW)
               val = s% lnPgas(k)/ln10 - 4d0*s% lnT(k)/ln10
            case (p_logQ)
               val = s% lnd(k)/ln10 - 2d0*s% lnT(k)/ln10 + 12d0
            case (p_logV)
               val = s% lnd(k)/ln10 - 0.7d0*s% lnE(k)/ln10 + 20d0

            case (p_log_zFe)
               val = 0d0
               do j=1,s% species
                  if (chem_isos% Z(s% chem_id(j)) >= 24) val = val + s% xa(j,k)
               end do
               val = safe_log10(val)
            case (p_zFe)
               val = 0d0
               do j=1,s% species
                  if (chem_isos% Z(s% chem_id(j)) >= 24) val = val + s% xa(j,k)
               end do
            case(p_log_u_residual)
               if (s% u_flag) val = safe_log10(abs(s% u_residual(k)))
            case(p_u_residual)
               if (s% u_flag) val = s% u_residual(k)
            case(p_u)
               if (s% u_flag) val = s% u(k)
            case(p_u_face)
               if (s% u_flag) val = s% u_face(k)
            case (p_dPdr_dRhodr_info)
               if (s% RTI_flag) val = s% dPdr_dRhodr_info(k)
            case(p_signed_log_ergs_err)
               val = sign(safe_log10( &
                        abs(s% E_residual(k)*s% energy_start(k)*s% dm(k))), &
                        s% E_residual(k))
            case(p_RTI_du_diffusion_kick)
               if (s% u_flag) val = s% RTI_du_diffusion_kick(k)
            case(p_log_du_kick_div_du)
               if (s% u_flag .and. k > 1) then
                  if (abs(s% u_face(k)) > 1d0) &
                     val = safe_log10(abs(s% RTI_du_diffusion_kick(k)/s% u_face(k)))
               end if
               
            case(p_max_abs_xa_corr)
               val = s% max_abs_xa_corr(k)

            case default
               write(*,*) 'FATAL ERROR in profile_getval', c, k
               write(*,*) 'between ' // trim(profile_column_name(c-1)) // ' and ' // &
                  trim(profile_column_name(c+1)), c-1, c+1
               val = 0
               stop 'profile_getval'

         end select

         end if


         contains
         

         real(dp) function get_L_vel(k) result(v) ! velocity if L carried by convection
            integer, intent(in) :: k
            real(dp) :: rho_face
            integer :: j
            if (k == 1) then
               j = 2
            else
               j = k
            end if
            rho_face = interp_val_to_pt(s% rho,j,nz,s% dq,'profile get_L_vel')
            v = pow(max(1d0,s% L(k))/(4*pi*s% r(k)*s% r(k)*rho_face),1d0/3d0)
         end function get_L_vel


         real(dp) function get_k_r_integral(k_in, el, nu_factor)
            integer, intent(in) :: k_in
            integer, intent(in) :: el
            real(dp), intent(in) :: nu_factor
            real(dp) :: integral, integral_for_k, &
               cs2, r2, n2, sl2, omega2, L2, kr2, dr
            integer :: k, k1, k_inner, k_outer
            include 'formats'

            if (k_in == 1) then
               get_k_r_integral = 1
               return
            end if

            get_k_r_integral = 0
            L2 = el*(el+1)
            omega2 = pow2(1d-6*2*pi*s% nu_max*nu_factor)

            ! k_inner and k_outer are bounds of evanescent region

            ! k_outer is outermost k where Sl2 <= omega2 at k-1 and Sl2 > omega2 at k
            ! 1st find outermost where Sl2 <= omega2
            k1 = 0
            do k = 2, s% nz
               r2 = s% r(k)*s% r(k)
               cs2 = s% csound_face(k)*s% csound_face(k)
               sl2 = L2*cs2/r2
               if (sl2 <= omega2) then
                  k1 = k; exit
               end if
            end do
            if (k1 == 0) return
            ! then find next k where Sl2 >= omega2
            k_outer = 0
            do k = k1+1, s% nz
               r2 = s% r(k)*s% r(k)
               cs2 = s% csound_face(k)*s% csound_face(k)
               sl2 = L2*cs2/r2
               if (sl2 > omega2) then
                  k_outer = k; exit
               end if
            end do
            if (k_outer == 0) return
            if (k_in <= k_outer) then
               get_k_r_integral = 1
               return
            end if

            ! k_inner is next k where N2 >= omega2 at k+1 and N2 < omega2 at k
            k_inner = 0
            do k = k_outer+1, s% nz
               if (s% brunt_N2(k) >= omega2) then
                  k_inner= k; exit
               end if
            end do
            if (k_inner == 0) return
            if (k_in > k_inner) then
               get_k_r_integral = 1
               return
            end if

            integral = 0; integral_for_k = 0
            get_k_r_integral = 0
            do k = k_inner, k_outer, -1
               r2 = s% r(k)*s% r(k)
               cs2 = s% csound_face(k)*s% csound_face(k)
               n2 = s% brunt_N2(k)
               sl2 = L2*cs2/r2
               kr2 = (1 - n2/omega2)*(1 - Sl2/omega2)/cs2
               dr = s% rmid(k-1) - s% rmid(k)
               if (kr2 < 0 .and. omega2 < Sl2) integral = integral + sqrt(-kr2)*dr
               if (k == k_in) integral_for_k = integral
            end do
            if (integral < 1d-99) return
            get_k_r_integral = integral_for_k/integral

            if (is_bad(get_k_r_integral)) then
               write(*,2) 'get_k_r_integral', k_in, integral_for_k, integral
               stop 'get_k_r_integral'
            end if

         end function get_k_r_integral


         real(dp) function get_brunt_N_div_r_integral(k_in)
            integer, intent(in) :: k_in
            real(dp) :: integral, integral_for_k, dr
            integer :: k
            integral = 0
            integral_for_k = 0
            get_brunt_N_div_r_integral = 1
            if (k_in == 1) return
            get_brunt_N_div_r_integral = 0
            do k = s% nz, 2, -1
               dr = s% rmid(k-1) - s% rmid(k)
               if (s% brunt_N2(k) > 0) &
                  integral = integral + sqrt(s% brunt_N2(k))*dr/s% r(k)
               if (k == k_in) integral_for_k = integral
            end do
            if (integral < 1d-99) return
            get_brunt_N_div_r_integral = integral_for_k/integral
         end function get_brunt_N_div_r_integral


         real(dp) function get_dlogX_dlogP(j, k)
            integer, intent(in) :: j, k
            integer :: ii, i
            real(dp) :: val, x00, xm1, dlogP, dlogX
            include 'formats'
            get_dlogx_dlogp = 0
            if (k > 1) then
               ii = k
            else
               ii = 2
            end if
            i = s% net_iso(j)
            if (i == 0) return
            x00 = s% xa(i,ii)
            xm1 = s% xa(i,ii-1)
            if (x00 < 1d-20 .or. xm1 < 1d-20) return
            dlogP = (s% lnP(ii) - s% lnP(ii-1))/ln10
            if (dlogP <= 0d0) return
            dlogX = log10(x00/xm1)
            get_dlogX_dlogP = dlogX/dlogP
         end function get_dlogX_dlogP


         real(dp) function get_dlog_eps_dlogP(cat, k)
            integer, intent(in) :: cat, k
            integer :: ii
            real(dp) :: val, eps, epsm1, dlogP, dlog_eps
            get_dlog_eps_dlogP = 0
            if (k > 1) then
               ii = k
            else
               ii = 2
            end if
            eps = s% eps_nuc_categories(cat,ii)
            epsm1 = s% eps_nuc_categories(cat,ii-1)
            if (eps < 1d-3 .or. epsm1 < 1d-3) return
            dlogP = (s% lnP(ii) - s% lnP(ii-1))/ln10
            if (dlogP <= 0d0) return
            dlog_eps = log10(eps/epsm1)
            get_dlog_eps_dlogP = dlog_eps/dlogP
         end function get_dlog_eps_dlogP


         real(dp) function pt(v,k)
            integer, intent(in) :: k
            real(dp), pointer :: v(:)
            if (k == 1) then
               pt = v(k)
            else
               pt = (v(k)*s% dq(k-1) + v(k-1)*s% dq(k))/(s% dq(k-1) + s% dq(k))
            endif
         end function pt


         real(dp) function if_rot(v,k, alt)
            real(dp),dimension(:), intent(in) :: v
            integer, intent(in) :: k
            real(dp), optional, intent(in) :: alt
            if (s% rotation_flag) then
               if_rot = v(k)
            else
               if (present(alt)) then
                  if_rot = alt
               else
                  if_rot = 0
               end if
            endif
         end function if_rot

      end subroutine getval_for_profile

      end module profile_getval

