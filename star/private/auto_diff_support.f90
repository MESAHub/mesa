! ***********************************************************************
!
!   Copyright (C) 2020  Adam Jermyn, Bill Paxton & The MESA Team
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

      module auto_diff_support

      use star_private_def
      use const_def
      use auto_diff

      implicit none

      public

      contains

      type(auto_diff_real_18var_order1) function shift_p1(val_00) result(val_p1)
         type(auto_diff_real_18var_order1), intent(in) :: val_00
         val_p1%val = val_00%val
         ! make sure copy each partial before overwrite it in case output is same as input
         val_p1%d1Array(i_L_p1) = val_00%d1Array(i_L_00)
         val_p1%d1Array(i_L_00) = val_00%d1Array(i_L_m1)
         val_p1%d1Array(i_L_m1) = 0d0
         val_p1%d1Array(i_v_p1) = val_00%d1Array(i_v_00)
         val_p1%d1Array(i_v_00) = val_00%d1Array(i_v_m1)
         val_p1%d1Array(i_v_m1) = 0d0
         val_p1%d1Array(i_lnR_p1) = val_00%d1Array(i_lnR_00)
         val_p1%d1Array(i_lnR_00) = val_00%d1Array(i_lnR_m1)
         val_p1%d1Array(i_lnR_m1) = 0d0
         val_p1%d1Array(i_w_p1) = val_00%d1Array(i_w_00)
         val_p1%d1Array(i_w_00) = val_00%d1Array(i_w_m1)
         val_p1%d1Array(i_w_m1) = 0d0
         val_p1%d1Array(i_lnT_p1) = val_00%d1Array(i_lnT_00)
         val_p1%d1Array(i_lnT_00) = val_00%d1Array(i_lnT_m1)
         val_p1%d1Array(i_lnT_m1) = 0d0
         val_p1%d1Array(i_lnd_p1) = val_00%d1Array(i_lnd_00)
         val_p1%d1Array(i_lnd_00) = val_00%d1Array(i_lnd_m1)
         val_p1%d1Array(i_lnd_m1) = 0d0
      end function shift_p1

      type(auto_diff_real_18var_order1) function shift_m1(val_00) result(val_m1)
         type(auto_diff_real_18var_order1), intent(in) :: val_00
         val_m1%val = val_00%val
         ! make sure copy each partial before overwrite it in case output is same as input
         val_m1%d1Array(i_lnd_m1) = val_00%d1Array(i_lnd_00)
         val_m1%d1Array(i_lnd_00) = val_00%d1Array(i_lnd_p1)
         val_m1%d1Array(i_lnd_p1) = 0d0
         val_m1%d1Array(i_lnT_m1) = val_00%d1Array(i_lnT_00)
         val_m1%d1Array(i_lnT_00) = val_00%d1Array(i_lnT_p1)
         val_m1%d1Array(i_lnT_p1) = 0d0
         val_m1%d1Array(i_w_m1) = val_00%d1Array(i_w_00)
         val_m1%d1Array(i_w_00) = val_00%d1Array(i_w_p1)
         val_m1%d1Array(i_w_p1) = 0d0
         val_m1%d1Array(i_lnR_m1) = val_00%d1Array(i_lnR_00)
         val_m1%d1Array(i_lnR_00) = val_00%d1Array(i_lnR_p1)
         val_m1%d1Array(i_lnR_p1) = 0d0
         val_m1%d1Array(i_v_m1) = val_00%d1Array(i_v_00)
         val_m1%d1Array(i_v_00) = val_00%d1Array(i_v_p1)
         val_m1%d1Array(i_v_p1) = 0d0
         val_m1%d1Array(i_L_m1) = val_00%d1Array(i_L_00)
         val_m1%d1Array(i_L_00) = val_00%d1Array(i_L_p1)
         val_m1%d1Array(i_L_p1) = 0d0
      end function shift_m1

      subroutine unwrap(var, val, &
            dlnd_m1, dlnd_00, dlnd_p1, &
            dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, &
            dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, &
            dL_m1, dL_00, dL_p1)
         type(auto_diff_real_18var_order1), intent(in) :: var
         real(dp), intent(out) :: &
            val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1
         val = var%val
         dlnd_m1 = var%d1Array(i_lnd_m1)
         dlnd_00 = var%d1Array(i_lnd_00)
         dlnd_p1 = var%d1Array(i_lnd_p1)
         dlnT_m1 = var%d1Array(i_lnT_m1)
         dlnT_00 = var%d1Array(i_lnT_00)
         dlnT_p1 = var%d1Array(i_lnT_p1)
         dw_m1 = var%d1Array(i_w_m1)
         dw_00 = var%d1Array(i_w_00)
         dw_p1 = var%d1Array(i_w_p1)
         dlnR_m1 = var%d1Array(i_lnR_m1)
         dlnR_00 = var%d1Array(i_lnR_00)
         dlnR_p1 = var%d1Array(i_lnR_p1)
         dv_m1 = var%d1Array(i_v_m1)
         dv_00 = var%d1Array(i_v_00)
         dv_p1 = var%d1Array(i_v_p1)
         dL_m1 = var%d1Array(i_L_m1)
         dL_00 = var%d1Array(i_L_00)
         dL_p1 = var%d1Array(i_L_p1)
      end subroutine unwrap

      subroutine wrap(var, val, &
            dlnd_m1, dlnd_00, dlnd_p1, &
            dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, &
            dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, &
            dL_m1, dL_00, dL_p1)
         type(auto_diff_real_18var_order1), intent(out) :: var
         real(dp), intent(in) :: &
            val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1
          var%val = val
          var%d1Array(i_lnd_m1) = dlnd_m1
          var%d1Array(i_lnd_00) = dlnd_00
          var%d1Array(i_lnd_p1) = dlnd_p1
          var%d1Array(i_lnT_m1) = dlnT_m1
          var%d1Array(i_lnT_00) = dlnT_00
          var%d1Array(i_lnT_p1) = dlnT_p1
          var%d1Array(i_w_m1) = dw_m1
          var%d1Array(i_w_00) = dw_00
          var%d1Array(i_w_p1) = dw_p1
          var%d1Array(i_lnR_m1) = dlnR_m1
          var%d1Array(i_lnR_00) = dlnR_00
          var%d1Array(i_lnR_p1) = dlnR_p1
          var%d1Array(i_v_m1) = dv_m1
          var%d1Array(i_v_00) = dv_00
          var%d1Array(i_v_p1) = dv_p1
          var%d1Array(i_L_m1) = dL_m1
          var%d1Array(i_L_00) = dL_00
          var%d1Array(i_L_p1) = dL_p1
      end subroutine wrap


      ! The following routines turn regular star variables into auto_diff_real_18var_order1 variables.
      ! For independent variables this is a straightforward wrapping. For dependent variables like eos and kap
      ! outputs we also pull in information about their partials from the relevant module.
      
      !----------------------------------------------------------------------------------------------------
      !
      ! We begin with quantities defined on cells (T, rho, e_turb, kap, pressure, energy, entropy).
      ! We need these on cells k-1, k and k+1.
      !
      !----------------------------------------------------------------------------------------------------


      !! Wrap the temperature at cell k-1 with appropriate dependences on the (ordered) independent variables.
      function wrap_T_m1(s, k) result(T_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: T_m1
         integer, intent(in) :: k
         T_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            T_m1 % val = s%T(k-1)
            T_m1 % d1Array(i_lnT_m1) = s%T(k-1)
         end if
      end function wrap_T_m1

      !! Wrap the temperature at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_T_00(s, k) result(T_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: T_00
         integer, intent(in) :: k
         T_00 = 0d0 ! sets val and d1Array to 0
         T_00 % val = s%T(k)
         T_00 % d1Array(i_lnT_00) = s%T(k)
      end function wrap_T_00

      !! Wrap the temperature at cell k+1 with appropriate dependences on the (ordered) independent variables.
      function wrap_T_p1(s, k) result(T_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: T_p1
         integer, intent(in) :: k
         T_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            T_p1 % val = s%T(k+1)
            T_p1 % d1Array(i_lnT_p1) = s%T(k+1)
         end if
      end function wrap_T_p1


      !! Wrap the density at cell k-1 with appropriate dependences on the (ordered) independent variables.
      function wrap_d_m1(s, k) result(d_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: d_m1
         integer, intent(in) :: k
         d_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            d_m1 % val = s%rho(k-1)
            d_m1 % d1Array(i_lnd_m1) = s%rho(k-1)
         end if
      end function wrap_d_m1

      !! Wrap the density at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_d_00(s, k) result(d_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: d_00
         integer, intent(in) :: k
         d_00 = 0d0 ! sets val and d1Array to 0
         d_00 % val = s%rho(k)
         d_00 % d1Array(i_lnd_00) = s%rho(k)
      end function wrap_d_00

      !! Wrap the density at cell k+1 with appropriate dependences on the (ordered) independent variables.
      function wrap_d_p1(s, k) result(d_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: d_p1
         integer, intent(in) :: k
         d_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            d_p1 % val = s%rho(k+1)
            d_p1 % d1Array(i_lnd_p1) = s%rho(k+1)
         end if
      end function wrap_d_p1


      !! Wrap et at cell k-1 with appropriate dependences on the (ordered) independent variables.
      function wrap_w_m1(s, k) result(w_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: w_m1
         integer, intent(in) :: k
         w_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            w_m1 % val = s%w(k-1)
            w_m1 % d1Array(i_w_m1) = 1d0
         end if            
      end function wrap_w_m1

      !! Wrap et at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_w_00(s, k) result(w_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: w_00
         integer, intent(in) :: k
         w_00 = 0d0 ! sets val and d1Array to 0
         w_00 % val = s%w(k)
         w_00 % d1Array(i_w_00) = 1d0
      end function wrap_w_00

      !! Wrap et at cell k+1 with appropriate dependences on the (ordered) independent variables.
      function wrap_w_p1(s, k) result(w_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: w_p1
         integer, intent(in) :: k
         w_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            w_p1 % val = s%w(k+1)
            w_p1 % d1Array(i_w_p1) = 1d0
         end if
      end function wrap_w_p1


      !! Wrap opacity at cell k-1 with appropriate dependences on the (ordered) independent variables.
      function wrap_kap_m1(s, k) result(kap_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: kap_m1
         integer, intent(in) :: k
         kap_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            kap_m1 % val = s%opacity(k-1)
            kap_m1 % d1Array(i_lnd_m1) = s%d_opacity_dlnd(k-1)
            kap_m1 % d1Array(i_lnT_m1) = s%d_opacity_dlnT(k-1)
         end if
      end function wrap_kap_m1

      !! Wrap opacity at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_kap_00(s, k) result(kap_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: kap_00
         integer, intent(in) :: k
         kap_00 = 0d0 ! sets val and d1Array to 0
         kap_00 % val = s%opacity(k)
         kap_00 % d1Array(i_lnd_00) = s%d_opacity_dlnd(k)
         kap_00 % d1Array(i_lnT_00) = s%d_opacity_dlnT(k)
      end function wrap_kap_00

      !! Wrap opacity at cell k+1 with appropriate dependences on the (ordered) independent variables.
      function wrap_kap_p1(s, k) result(kap_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: kap_p1
         integer, intent(in) :: k
         kap_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            kap_p1 % val = s%opacity(k+1)
            kap_p1 % d1Array(i_lnd_p1) = s%d_opacity_dlnd(k+1)
            kap_p1 % d1Array(i_lnT_p1) = s%d_opacity_dlnT(k+1)
         end if
      end function wrap_kap_p1


      !! Wrap s for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_s_m1(s, k) result(s_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: s_m1
         integer, intent(in) :: k
         s_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            s_m1%val = s% entropy(k-1)
            s_m1%d1Array(i_lnd_m1) = s% dS_dRho_for_partials(k-1)*s% rho(k-1)
            s_m1%d1Array(i_lnT_m1) = s% dS_dT_for_partials(k-1)*s% T(k-1)
         end if   
      end function wrap_s_m1

      !! Wrap s for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_s_00(s, k) result(s_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: s_00
         integer, intent(in) :: k
         s_00 = 0d0 ! sets val and d1Array to 0
         s_00%val = s% entropy(k)
         s_00%d1Array(i_lnd_00) = s% dS_dRho_for_partials(k)*s% rho(k)
         s_00%d1Array(i_lnT_00) = s% dS_dT_for_partials(k)*s% T(k)
      end function wrap_s_00

      !! Wrap s for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_s_p1(s, k) result(s_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: s_p1
         integer, intent(in) :: k
         s_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            s_p1%val = s% entropy(k+1)
            s_p1%d1Array(i_lnd_p1) = s% dS_dRho_for_partials(k+1)*s% rho(k+1)
            s_p1%d1Array(i_lnT_p1) = s% dS_dT_for_partials(k+1)*s% T(k+1)
         end if
      end function wrap_s_p1


      !! Wrap internal energy for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_e_m1(s, k) result(e_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: e_m1
         integer, intent(in) :: k
         e_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            e_m1%val = s% energy(k-1)
            e_m1%d1Array(i_lnd_m1) = s% dE_dRho_for_partials(k-1)*s% rho(k-1)
            e_m1%d1Array(i_lnT_m1) = s% Cv_for_partials(k-1)*s% T(k-1)
         end if   
      end function wrap_e_m1

      !! Wrap internal energy for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_e_00(s, k) result(e_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: e_00
         integer, intent(in) :: k
         e_00 = 0d0 ! sets val and d1Array to 0
         e_00%val = s% energy(k)
         e_00%d1Array(i_lnd_00) = s% dE_dRho_for_partials(k)*s% rho(k)
         e_00%d1Array(i_lnT_00) = s% Cv_for_partials(k)*s% T(k)
      end function wrap_e_00

      !! Wrap internal energy for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_e_p1(s, k) result(e_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: e_p1
         integer, intent(in) :: k
         e_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            e_p1%val = s% energy(k+1)
            e_p1%d1Array(i_lnd_p1) = s% dE_dRho_for_partials(k+1)*s% rho(k+1)
            e_p1%d1Array(i_lnT_p1) = s% Cv_for_partials(k+1)*s% T(k+1)
         end if
      end function wrap_e_p1


      !! Wrap P for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_p_m1(s, k) result(p_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: p_m1
         integer, intent(in) :: k
         p_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            p_m1%val = s% P(k-1)
            p_m1%d1Array(i_lnd_m1) = s%P(k-1) * s% chiRho_for_partials(k-1)
            p_m1%d1Array(i_lnT_m1) = s%P(k-1) * s% chiT_for_partials(k-1)
         end if   
      end function wrap_p_m1

      !! Wrap P for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_p_00(s, k) result(p_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: p_00
         integer, intent(in) :: k
         p_00 = 0d0 ! sets val and d1Array to 0
         p_00%val = s% P(k)
         p_00%d1Array(i_lnd_00) = s%P(k) * s% chiRho_for_partials(k)
         p_00%d1Array(i_lnT_00) = s%P(k) * s% chiT_for_partials(k)
      end function wrap_p_00

      !! Wrap P for cell k+1 with appropriate dependences on the (ordered) independent variables
      function wrap_p_p1(s, k) result(p_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: p_p1
         integer, intent(in) :: k
         p_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            p_p1%val = s% P(k+1)
            p_p1%d1Array(i_lnd_p1) = s%P(k+1) * s% chiRho_for_partials(k+1)
            p_p1%d1Array(i_lnT_p1) = s%P(k+1) * s% chiT_for_partials(k+1)
         end if   
      end function wrap_p_p1


      !! Wrap ChiRho for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_ChiRho_m1(s, k) result(ChiRho_m1)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: ChiRho_m1
         integer, intent(in) :: k
         ChiRho_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            ChiRho_m1%val = s% ChiRho(k-1)
            ChiRho_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiRho,k-1)
            ChiRho_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiRho,k-1)
         end if
      end function wrap_ChiRho_m1

      !! Wrap ChiRho for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_ChiRho_00(s, k) result(ChiRho_00)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: ChiRho_00
         integer, intent(in) :: k
         ChiRho_00 = 0d0 ! sets val and d1Array to 0
         ChiRho_00%val = s% ChiRho(k)
         ChiRho_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiRho,k)
         ChiRho_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiRho,k)
      end function wrap_ChiRho_00

      !! Wrap ChiRho for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_ChiRho_p1(s, k) result(ChiRho_p1)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: ChiRho_p1
         integer, intent(in) :: k
         ChiRho_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s% nz) then
            ChiRho_p1%val = s% ChiRho(k+1)
            ChiRho_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiRho,k+1)
            ChiRho_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiRho,k+1)
         end if
      end function wrap_ChiRho_p1

      !! Wrap ChiT for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_ChiT_m1(s, k) result(ChiT_m1)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: ChiT_m1
         integer, intent(in) :: k
         ChiT_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            ChiT_m1%val = s% ChiT(k-1)
            ChiT_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiT,k-1)
            ChiT_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiT,k-1)
         end if
      end function wrap_ChiT_m1

      !! Wrap ChiT for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_ChiT_00(s, k) result(ChiT_00)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: ChiT_00
         integer, intent(in) :: k
         ChiT_00 = 0d0 ! sets val and d1Array to 0
         ChiT_00%val = s% ChiT(k)
         ChiT_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiT,k)
         ChiT_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiT,k)
      end function wrap_ChiT_00

      !! Wrap ChiT for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_ChiT_p1(s, k) result(ChiT_p1)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: ChiT_p1
         integer, intent(in) :: k
         ChiT_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s% nz) then
            ChiT_p1%val = s% ChiT(k+1)
            ChiT_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiT,k+1)
            ChiT_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiT,k+1)
         end if
      end function wrap_ChiT_p1


      !! Wrap Cp for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_Cp_m1(s, k) result(Cp_m1)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: Cp_m1
         integer, intent(in) :: k
         Cp_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            Cp_m1%val = s% Cp(k-1)
            Cp_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_Cp,k-1)
            Cp_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_Cp,k-1)
         end if   
      end function wrap_Cp_m1

      !! Wrap Cp for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_Cp_00(s, k) result(Cp_00)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: Cp_00
         integer, intent(in) :: k
         Cp_00 = 0d0 ! sets val and d1Array to 0
         Cp_00%val = s% Cp(k)
         Cp_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_Cp,k)
         Cp_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_Cp,k)
      end function wrap_Cp_00

      !! Wrap Cp for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_Cp_p1(s, k) result(Cp_p1)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: Cp_p1
         integer, intent(in) :: k
         Cp_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s% nz) then
            Cp_p1%val = s% Cp(k+1)
            Cp_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_Cp,k+1)
            Cp_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_Cp,k+1)
         end if   
      end function wrap_Cp_p1


      !! Wrap gamma1 for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_gamma1_m1(s, k) result(gamma1_m1)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: gamma1_m1
         integer, intent(in) :: k
         gamma1_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            gamma1_m1%val = s% gamma1(k-1)
            gamma1_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_gamma1,k-1)
            gamma1_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_gamma1,k-1)
         end if   
      end function wrap_gamma1_m1

      !! Wrap gamma1 for cell k with appropriate dependences on the (ordered) independent variables
      function wrap_gamma1_00(s, k) result(gamma1_00)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: gamma1_00
         integer, intent(in) :: k
         gamma1_00 = 0d0 ! sets val and d1Array to 0
         gamma1_00%val = s% gamma1(k)
         gamma1_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_gamma1,k)
         gamma1_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_gamma1,k)
      end function wrap_gamma1_00

      !! Wrap gamma1 for cell k-1 with appropriate dependences on the (ordered) independent variables
      function wrap_gamma1_p1(s, k) result(gamma1_p1)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: gamma1_p1
         integer, intent(in) :: k
         gamma1_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s% nz) then
            gamma1_p1%val = s% gamma1(k+1)
            gamma1_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_gamma1,k+1)
            gamma1_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_gamma1,k+1)
         end if   
      end function wrap_gamma1_p1


      !----------------------------------------------------------------------------------------------------
      !
      ! Next we handle quantities defined on faces (r, v, L). We need these on faces k+1, k and k-1.
      !
      !----------------------------------------------------------------------------------------------------


      !! Wrap L at cell k-1 with appropriate dependences on the (ordered) independent variables.
      function wrap_L_m1(s, k) result(L_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: L_m1
         integer, intent(in) :: k
         L_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            L_m1 % val = s%L(k-1)
            L_m1 % d1Array(i_L_m1) = 1d0
         end if
      end function wrap_L_m1

      !! Wrap L at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_L_00(s, k) result(L_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: L_00
         integer, intent(in) :: k
         L_00 = 0d0 ! sets val and d1Array to 0
         L_00 % val = s%L(k)
         L_00 % d1Array(i_L_00) = 1d0
      end function wrap_L_00

      !! Wrap L at cell k+1 with appropriate dependences on the (ordered) independent variables.
      function wrap_L_p1(s, k) result(L_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: L_p1
         integer, intent(in) :: k
         L_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            L_p1 % val = s%L(k+1)
            L_p1 % d1Array(i_L_p1) = 1d0
         else
            L_p1 %val = s%L_center
            ! L_center is a constant
         end if
      end function wrap_L_p1


      !! Wrap r at cell k-1 with appropriate dependences on the (ordered) independent variables.
      function wrap_r_m1(s, k) result(r_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: r_m1
         integer, intent(in) :: k
         r_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            r_m1 % val = s%r(k-1)
            r_m1 % d1Array(i_lnR_m1) = s%r(k-1)
         end if
      end function wrap_r_m1

      !! Wrap r at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_r_00(s, k) result(r_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: r_00
         integer, intent(in) :: k
         r_00 = 0d0 ! sets val and d1Array to 0
         r_00 % val = s%r(k)
         r_00 % d1Array(i_lnR_00) = s%r(k)
      end function wrap_r_00

      !! Wrap r at cell k+1 with appropriate dependences on the (ordered) independent variables.
      function wrap_r_p1(s, k) result(r_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: r_p1
         integer, intent(in) :: k
         r_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            r_p1 % val = s%r(k+1)
            r_p1 % d1Array(i_lnR_p1) = s%r(k+1)
         else
            r_p1 % val = s%r_center
         end if
      end function wrap_r_p1


      !! Wrap v at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_v_m1(s, k) result(v_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: v_m1
         integer, intent(in) :: k
         v_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            v_m1 % val = s%v(k-1)
            v_m1 % d1Array(i_v_m1) = 1d0
         end if
      end function wrap_v_m1

      !! Wrap v at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_v_00(s, k) result(v_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: v_00
         integer, intent(in) :: k
         v_00 = 0d0 ! sets val and d1Array to 0
         v_00 % val = s%v(k)
         v_00 % d1Array(i_v_00) = 1d0
      end function wrap_v_00

      !! Wrap v at cell k with appropriate dependences on the (ordered) independent variables.
      function wrap_v_p1(s, k) result(v_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: v_p1
         integer, intent(in) :: k
         v_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            v_p1 % val = s%v(k+1)
            v_p1 % d1Array(i_v_p1) = 1d0
         else
            v_p1 % val = s%v_center
            ! v_center is a constant
         end if
      end function wrap_v_p1

      ! u replaces v
      function wrap_u_m1(s, k) result(v_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: v_m1
         integer, intent(in) :: k
         v_m1 = 0d0 ! sets val and d1Array to 0
         if (k > 1) then
            v_m1 % val = s%u(k-1)
            v_m1 % d1Array(i_v_m1) = 1d0
         end if
      end function wrap_u_m1

      function wrap_u_00(s, k) result(v_00)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: v_00
         integer, intent(in) :: k
         v_00 = 0d0 ! sets val and d1Array to 0
         v_00 % val = s%u(k)
         v_00 % d1Array(i_v_00) = 1d0
      end function wrap_u_00

      function wrap_u_p1(s, k) result(v_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_18var_order1) :: v_p1
         integer, intent(in) :: k
         v_p1 = 0d0 ! sets val and d1Array to 0
         if (k < s%nz) then
            v_p1 % val = s%u(k+1)
            v_p1 % d1Array(i_v_p1) = 1d0
         else
            v_p1 % val = 0d0
            ! v_center is a constant
         end if
      end function wrap_u_p1


end module auto_diff_support
