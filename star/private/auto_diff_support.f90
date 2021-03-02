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
      
      ! current use of xtra's
      ! xtra1 is ln_cvpv0
      ! xtra2 is w_div_wc

      public

      contains

      type(auto_diff_real_star_order1) function shift_p1(val_00) result(val_p1)
         type(auto_diff_real_star_order1), intent(in) :: val_00
         integer :: j
         val_p1%val = val_00%val
         do j=auto_diff_star_num_vars-2,1,-3 ! p1 gets 00, 00 gets m1, m1 gets 0d0
            val_p1%d1Array(j+2) = val_00%d1Array(j+1)
            val_p1%d1Array(j+1) = val_00%d1Array(j)
            val_p1%d1Array(j) = 0d0
         end do
      end function shift_p1

      type(auto_diff_real_star_order1) function shift_m1(val_00) result(val_m1)
         type(auto_diff_real_star_order1), intent(in) :: val_00
         integer :: j
         val_m1%val = val_00%val
         do j=1,auto_diff_star_num_vars,3 ! m1 gets 00, 00 gets p1, p1 gets 0d0
            val_m1%d1Array(j) = val_00%d1Array(j+1)
            val_m1%d1Array(j+1) = val_00%d1Array(j+2)
            val_m1%d1Array(j+2) = 0d0
         end do
      end function shift_m1

      subroutine unwrap(var, val, &
            dlnd_m1, dlnd_00, dlnd_p1, &
            dlnT_m1, dlnT_00, dlnT_p1, &
            detrb_m1, detrb_00, detrb_p1, &
            dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, &
            dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dHp_m1, dHp_00, dHp_p1)
         type(auto_diff_real_star_order1), intent(in) :: var
         real(dp), intent(out) :: &
            val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            detrb_m1, detrb_00, detrb_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dHp_m1, dHp_00, dHp_p1
         val = var%val
         dlnd_m1 = var%d1Array(i_lnd_m1)
         dlnd_00 = var%d1Array(i_lnd_00)
         dlnd_p1 = var%d1Array(i_lnd_p1)
         dlnT_m1 = var%d1Array(i_lnT_m1)
         dlnT_00 = var%d1Array(i_lnT_00)
         dlnT_p1 = var%d1Array(i_lnT_p1)
         detrb_m1 = var%d1Array(i_etrb_m1)
         detrb_00 = var%d1Array(i_etrb_00)
         detrb_p1 = var%d1Array(i_etrb_p1)
         dlnR_m1 = var%d1Array(i_lnR_m1)
         dlnR_00 = var%d1Array(i_lnR_00)
         dlnR_p1 = var%d1Array(i_lnR_p1)
         dv_m1 = var%d1Array(i_v_m1)
         dv_00 = var%d1Array(i_v_00)
         dv_p1 = var%d1Array(i_v_p1)
         dL_m1 = var%d1Array(i_L_m1)
         dL_00 = var%d1Array(i_L_00)
         dL_p1 = var%d1Array(i_L_p1)
         dxtra1_m1 = var%d1Array(i_xtra1_m1)
         dxtra1_00 = var%d1Array(i_xtra1_00)
         dxtra1_p1 = var%d1Array(i_xtra1_p1)
         dxtra2_m1 = var%d1Array(i_xtra2_m1)
         dxtra2_00 = var%d1Array(i_xtra2_00)
         dxtra2_p1 = var%d1Array(i_xtra2_p1)
         dHp_m1 = var%d1Array(i_Hp_m1)
         dHp_00 = var%d1Array(i_Hp_00)
         dHp_p1 = var%d1Array(i_Hp_p1)
      end subroutine unwrap

      subroutine wrap(var, val, &
            dlnd_m1, dlnd_00, dlnd_p1, &
            dlnT_m1, dlnT_00, dlnT_p1, &
            detrb_m1, detrb_00, detrb_p1, &
            dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, &
            dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dHp_m1, dHp_00, dHp_p1)
         type(auto_diff_real_star_order1), intent(out) :: var
         real(dp), intent(in) :: &
            val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            detrb_m1, detrb_00, detrb_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dHp_m1, dHp_00, dHp_p1
         var%val = val
         var%d1Array(i_lnd_m1) = dlnd_m1
         var%d1Array(i_lnd_00) = dlnd_00
         var%d1Array(i_lnd_p1) = dlnd_p1
         var%d1Array(i_lnT_m1) = dlnT_m1
         var%d1Array(i_lnT_00) = dlnT_00
         var%d1Array(i_lnT_p1) = dlnT_p1
         var%d1Array(i_etrb_m1) = detrb_m1
         var%d1Array(i_etrb_00) = detrb_00
         var%d1Array(i_etrb_p1) = detrb_p1
         var%d1Array(i_lnR_m1) = dlnR_m1
         var%d1Array(i_lnR_00) = dlnR_00
         var%d1Array(i_lnR_p1) = dlnR_p1
         var%d1Array(i_v_m1) = dv_m1
         var%d1Array(i_v_00) = dv_00
         var%d1Array(i_v_p1) = dv_p1
         var%d1Array(i_L_m1) = dL_m1
         var%d1Array(i_L_00) = dL_00
         var%d1Array(i_L_p1) = dL_p1
         var%d1Array(i_xtra1_m1) = dxtra1_m1
         var%d1Array(i_xtra1_00) = dxtra1_00
         var%d1Array(i_xtra1_p1) = dxtra1_p1
         var%d1Array(i_xtra2_m1) = dxtra2_m1
         var%d1Array(i_xtra2_00) = dxtra2_00
         var%d1Array(i_xtra2_p1) = dxtra2_p1
         var%d1Array(i_Hp_m1) = dHp_m1
         var%d1Array(i_Hp_00) = dHp_00
         var%d1Array(i_Hp_p1) = dHp_p1
      end subroutine wrap


      ! The following routines turn regular star variables into auto_diff_real_star_order1 variables.
      ! For independent variables this is a straightforward wrapping. For dependent variables like eos and kap
      ! outputs we also pull in information about their partials from the relevant module.
      
      !----------------------------------------------------------------------------------------------------
      !
      ! We begin with quantities defined on cells (T, rho, e_turb, kap, pressure, energy, entropy).
      ! We need these on cells k-1, k and k+1.
      !
      !----------------------------------------------------------------------------------------------------


      function wrap_T_m1(s, k) result(T_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: T_m1
         integer, intent(in) :: k
         T_m1 = 0d0 
         if (k > 1) then
            T_m1 % val = s%T(k-1)
            if (s% solver_use_lnT) then
               T_m1 % d1Array(i_lnT_m1) = s%T(k-1)
            else
               T_m1 % d1Array(i_lnT_m1) = 1d0
            end if
         end if
      end function wrap_T_m1

      function wrap_T_00(s, k) result(T_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: T_00
         integer, intent(in) :: k
         T_00 = 0d0 
         T_00 % val = s%T(k)
         if (s% solver_use_lnT) then
            T_00 % d1Array(i_lnT_00) = s%T(k)
         else
            T_00 % d1Array(i_lnT_00) = 1d0
         end if
      end function wrap_T_00

      function wrap_T_p1(s, k) result(T_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: T_p1
         integer, intent(in) :: k
         T_p1 = 0d0 
         if (k < s%nz) then
            T_p1 % val = s%T(k+1)
            if (s% solver_use_lnT) then
               T_p1 % d1Array(i_lnT_p1) = s%T(k+1)
            else
               T_p1 % d1Array(i_lnT_p1) = 1d0
            end if
         end if
      end function wrap_T_p1

      function wrap_lnT_m1(s, k) result(lnT_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnT_m1
         integer, intent(in) :: k
         lnT_m1 = 0d0 
         if (k > 1) then
            lnT_m1 % val = s%lnT(k-1)
            if (s% solver_use_lnT) then
               lnT_m1 % d1Array(i_lnT_m1) = 1d0
            else
               lnT_m1 % d1Array(i_lnT_m1) = 1d0/s% T(k-1)
            end if
         end if
      end function wrap_lnT_m1

      function wrap_lnT_00(s, k) result(lnT_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnT_00
         integer, intent(in) :: k
         lnT_00 = 0d0 
         lnT_00 % val = s%lnT(k)
         if (s% solver_use_lnT) then
            lnT_00 % d1Array(i_lnT_00) = 1d0
         else
            lnT_00 % d1Array(i_lnT_00) = 1d0/s% T(k)
         end if
      end function wrap_lnT_00

      function wrap_lnT_p1(s, k) result(lnT_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnT_p1
         integer, intent(in) :: k
         lnT_p1 = 0d0 
         if (k < s%nz) then
            lnT_p1 % val = s%lnT(k+1)
            if (s% solver_use_lnT) then
               lnT_p1 % d1Array(i_lnT_p1) = 1d0
            else
               lnT_p1 % d1Array(i_lnT_p1) = 1d0/s% T(k+1)
            end if
         end if
      end function wrap_lnT_p1

      function wrap_dxh_lnT(s, k) result(dxh_lnT)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: dxh_lnT
         integer, intent(in) :: k
         dxh_lnT = 0d0 
         dxh_lnT % val = s%dxh_lnT(k)
         dxh_lnT % d1Array(i_lnT_00) = 1d0
      end function wrap_dxh_lnT

      function wrap_d_m1(s, k) result(d_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: d_m1
         integer, intent(in) :: k
         d_m1 = 0d0 
         if (k > 1) then
            d_m1 % val = s%rho(k-1)
            if (s% solver_use_lnd) then
               d_m1 % d1Array(i_lnd_m1) = s%rho(k-1)
            else
               d_m1 % d1Array(i_lnd_m1) = 1d0
            end if
         end if
      end function wrap_d_m1

      function wrap_d_00(s, k) result(d_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: d_00
         integer, intent(in) :: k
         d_00 = 0d0 
         d_00 % val = s%rho(k)
         if (s% solver_use_lnd) then
            d_00 % d1Array(i_lnd_00) = s%rho(k)
         else
            d_00 % d1Array(i_lnd_00) = 1d0
         end if
      end function wrap_d_00

      function wrap_d_p1(s, k) result(d_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: d_p1
         integer, intent(in) :: k
         d_p1 = 0d0 
         if (k < s%nz) then
            d_p1 % val = s%rho(k+1)
            if (s% solver_use_lnd) then
               d_p1 % d1Array(i_lnd_p1) = s%rho(k+1)
            else
               d_p1 % d1Array(i_lnd_p1) = 1d0
            end if
         end if
      end function wrap_d_p1

      function wrap_lnd_m1(s, k) result(lnd_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnd_m1
         integer, intent(in) :: k
         lnd_m1 = 0d0 
         if (k > 1) then
            lnd_m1 % val = s%lnd(k-1)
            if (s% solver_use_lnd) then
               lnd_m1 % d1Array(i_lnd_m1) = 1d0
            else
               lnd_m1 % d1Array(i_lnd_m1) = 1d0/s% rho(k-1)
            end if
         end if
      end function wrap_lnd_m1

      function wrap_lnd_00(s, k) result(lnd_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnd_00
         integer, intent(in) :: k
         lnd_00 = 0d0 
         lnd_00 % val = s%lnd(k)
         if (s% solver_use_lnd) then
            lnd_00 % d1Array(i_lnd_00) = 1d0
         else
            lnd_00 % d1Array(i_lnd_00) = 1d0/s% rho(k)
         end if
      end function wrap_lnd_00

      function wrap_lnd_p1(s, k) result(lnd_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnd_p1
         integer, intent(in) :: k
         lnd_p1 = 0d0 
         if (k < s%nz) then
            lnd_p1 % val = s%lnd(k+1)
            if (s% solver_use_lnd) then
               lnd_p1 % d1Array(i_lnd_p1) = 1d0
            else
               lnd_p1 % d1Array(i_lnd_p1) = 1d0/s% rho(k+1)
            end if
         end if
      end function wrap_lnd_p1

      function wrap_dxh_lnd(s, k) result(dxh_lnd)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: dxh_lnd
         integer, intent(in) :: k
         dxh_lnd = 0d0 
         dxh_lnd % val = s%dxh_lnd(k)
         dxh_lnd % d1Array(i_lnd_00) = 1d0
      end function wrap_dxh_lnd

      function wrap_etrb_m1(s, k) result(etrb_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: etrb_m1
         integer, intent(in) :: k
         etrb_m1 = 0d0 
         if (k > 1) then
            etrb_m1 % val = s%etrb(k-1)
            etrb_m1 % d1Array(i_etrb_m1) = 1d0
         end if            
      end function wrap_etrb_m1

      function wrap_etrb_00(s, k) result(etrb_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: etrb_00
         integer, intent(in) :: k
         etrb_00 = 0d0 
         etrb_00 % val = s%etrb(k)
         etrb_00 % d1Array(i_etrb_00) = 1d0
      end function wrap_etrb_00

      function wrap_etrb_p1(s, k) result(etrb_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: etrb_p1
         integer, intent(in) :: k
         etrb_p1 = 0d0 
         if (k < s%nz) then
            etrb_p1 % val = s%etrb(k+1)
            etrb_p1 % d1Array(i_etrb_p1) = 1d0
         end if
      end function wrap_etrb_p1

      function wrap_dxh_etrb(s, k) result(dxh_etrb)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: dxh_etrb
         integer, intent(in) :: k
         dxh_etrb = 0d0 
         dxh_etrb % val = s%dxh_etrb(k)
         dxh_etrb % d1Array(i_etrb_00) = 1d0
      end function wrap_dxh_etrb

      function wrap_w_m1(s, k) result(w_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: w_m1
         integer, intent(in) :: k
         w_m1 = 0d0 
         if (k > 1) then
            w_m1 % val = s%w(k-1)
            w_m1 % d1Array(i_etrb_m1) = 0.5d0/s%w(k-1)
         end if            
      end function wrap_w_m1

      function wrap_w_00(s, k) result(w_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: w_00
         integer, intent(in) :: k
         w_00 = 0d0 
         w_00 % val = s%w(k)
         w_00 % d1Array(i_etrb_00) = 0.5d0/s%w(k)
      end function wrap_w_00

      function wrap_w_p1(s, k) result(w_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: w_p1
         integer, intent(in) :: k
         w_p1 = 0d0 
         if (k < s%nz) then
            w_p1 % val = s%w(k+1)
            w_p1 % d1Array(i_etrb_p1) = 0.5d0/s%w(k+1)
         end if
      end function wrap_w_p1

      function wrap_kap_m1(s, k) result(kap_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: kap_m1
         integer, intent(in) :: k
         kap_m1 = 0d0 
         if (k > 1) then
            kap_m1 % val = s%opacity(k-1)
            if (s% solver_use_lnd) then
               kap_m1 % d1Array(i_lnd_m1) = s%d_opacity_dlnd(k-1)
            else
               kap_m1 % d1Array(i_lnd_m1) = s%d_opacity_dlnd(k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               kap_m1 % d1Array(i_lnT_m1) = s%d_opacity_dlnT(k-1)
            else
               kap_m1 % d1Array(i_lnT_m1) = s%d_opacity_dlnT(k-1)/s% T(k-1)
            end if
         end if
      end function wrap_kap_m1

      function wrap_kap_00(s, k) result(kap_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: kap_00
         integer, intent(in) :: k
         kap_00 = 0d0 
         kap_00 % val = s%opacity(k)
         if (s% solver_use_lnd) then
            kap_00 % d1Array(i_lnd_00) = s%d_opacity_dlnd(k)
         else
            kap_00 % d1Array(i_lnd_00) = s%d_opacity_dlnd(k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            kap_00 % d1Array(i_lnT_00) = s%d_opacity_dlnT(k)
         else
            kap_00 % d1Array(i_lnT_00) = s%d_opacity_dlnT(k)/s% T(k)
         end if
      end function wrap_kap_00

      function wrap_kap_p1(s, k) result(kap_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: kap_p1
         integer, intent(in) :: k
         kap_p1 = 0d0 
         if (k < s%nz) then
            kap_p1 % val = s%opacity(k+1)
            if (s% solver_use_lnd) then
               kap_p1 % d1Array(i_lnd_p1) = s%d_opacity_dlnd(k+1)/s% rho(k+1)
            else
               kap_p1 % d1Array(i_lnd_p1) = s%d_opacity_dlnd(k+1)
            end if
            if (s% solver_use_lnT) then
               kap_p1 % d1Array(i_lnT_p1) = s%d_opacity_dlnT(k+1)
            else
               kap_p1 % d1Array(i_lnT_p1) = s%d_opacity_dlnT(k+1)/s% T(k+1)
            end if
         end if
      end function wrap_kap_p1

      function wrap_s_m1(s, k) result(s_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: s_m1
         integer, intent(in) :: k
         s_m1 = 0d0 
         if (k > 1) then
            s_m1%val = s% entropy(k-1)
            if (s% solver_use_lnd) then
               s_m1%d1Array(i_lnd_m1) = s% dS_dRho_for_partials(k-1)*s% rho(k-1)
            else
               s_m1%d1Array(i_lnd_m1) = s% dS_dRho_for_partials(k-1)
            end if
            if (s% solver_use_lnT) then
               s_m1%d1Array(i_lnT_m1) = s% dS_dT_for_partials(k-1)*s% T(k-1)
            else
               s_m1%d1Array(i_lnT_m1) = s% dS_dT_for_partials(k-1)
            end if
         end if   
      end function wrap_s_m1

      function wrap_s_00(s, k) result(s_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: s_00
         integer, intent(in) :: k
         s_00 = 0d0 
         s_00%val = s% entropy(k)
         if (s% solver_use_lnd) then
            s_00%d1Array(i_lnd_00) = s% dS_dRho_for_partials(k)*s% rho(k)
         else
            s_00%d1Array(i_lnd_00) = s% dS_dRho_for_partials(k)
         end if
         if (s% solver_use_lnT) then
            s_00%d1Array(i_lnT_00) = s% dS_dT_for_partials(k)*s% T(k)
         else
            s_00%d1Array(i_lnT_00) = s%dS_dT_for_partials(k)
         end if
      end function wrap_s_00

      function wrap_s_p1(s, k) result(s_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: s_p1
         integer, intent(in) :: k
         s_p1 = 0d0 
         if (k < s%nz) then
            s_p1%val = s% entropy(k+1)
            if (s% solver_use_lnd) then
               s_p1%d1Array(i_lnd_p1) = s% dS_dRho_for_partials(k+1)*s% rho(k+1)
            else
               s_p1%d1Array(i_lnd_p1) = s% dS_dRho_for_partials(k+1)
            end if
            if (s% solver_use_lnT) then
               s_p1%d1Array(i_lnT_p1) = s% dS_dT_for_partials(k+1)*s% T(k+1)
            else
               s_p1 % d1Array(i_lnT_p1) = s%dS_dT_for_partials(k+1)
            end if
         end if
      end function wrap_s_p1

      function wrap_e_m1(s, k) result(e_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: e_m1
         integer, intent(in) :: k
         e_m1 = 0d0 
         if (k > 1) then
            e_m1%val = s% energy(k-1)
            if (s% solver_use_lnd) then
               e_m1%d1Array(i_lnd_m1) = s% dE_dRho_for_partials(k-1)*s% rho(k-1)
            else
               e_m1%d1Array(i_lnd_m1) = s% dE_dRho_for_partials(k-1)
            end if
            if (s% solver_use_lnT) then
               e_m1%d1Array(i_lnT_m1) = s% Cv_for_partials(k-1)*s% T(k-1)
            else
               e_m1%d1Array(i_lnT_m1) = s% Cv_for_partials(k-1)
            end if
         end if   
      end function wrap_e_m1

      function wrap_e_00(s, k) result(e_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: e_00
         integer, intent(in) :: k
         e_00 = 0d0 
         e_00%val = s% energy(k)
         if (s% solver_use_lnd) then
            e_00%d1Array(i_lnd_00) = s% dE_dRho_for_partials(k)*s% rho(k)
         else
            e_00%d1Array(i_lnd_00) = s% dE_dRho_for_partials(k)
         end if
         if (s% solver_use_lnT) then
            e_00%d1Array(i_lnT_00) = s% Cv_for_partials(k)*s% T(k)
         else
            e_00%d1Array(i_lnT_00) = s% Cv_for_partials(k)
         end if
      end function wrap_e_00

      function wrap_e_p1(s, k) result(e_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: e_p1
         integer, intent(in) :: k
         e_p1 = 0d0 
         if (k < s%nz) then
            e_p1%val = s% energy(k+1)
            e_p1%d1Array(i_lnd_p1) = s% dE_dRho_for_partials(k+1)*s% rho(k+1)
            if (s% solver_use_lnT) then
               e_p1%d1Array(i_lnT_p1) = s% Cv_for_partials(k+1)*s% T(k+1)
            else
               e_p1%d1Array(i_lnT_p1) = s% Cv_for_partials(k+1)
            end if
         end if
      end function wrap_e_p1

      function wrap_p_m1(s, k) result(p_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: p_m1
         integer, intent(in) :: k
         p_m1 = 0d0 
         if (k > 1) then
            p_m1%val = s% P(k-1)
            if (s% solver_use_lnd) then
               p_m1%d1Array(i_lnd_m1) = s%P(k-1) * s% chiRho_for_partials(k-1)
            else
               p_m1%d1Array(i_lnd_m1) = s%P(k-1) * s% chiRho_for_partials(k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               p_m1%d1Array(i_lnT_m1) = s%P(k-1) * s% chiT_for_partials(k-1)
            else
               p_m1%d1Array(i_lnT_m1) = s%P(k-1) * s% chiT_for_partials(k-1)/s% T(k-1)
            end if
         end if   
      end function wrap_p_m1

      function wrap_p_00(s, k) result(p_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: p_00
         integer, intent(in) :: k
         p_00 = 0d0 
         p_00%val = s% P(k)
         if (s% solver_use_lnd) then
            p_00%d1Array(i_lnd_00) = s%P(k) * s% chiRho_for_partials(k)
         else
            p_00%d1Array(i_lnd_00) = s%P(k) * s% chiRho_for_partials(k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            p_00%d1Array(i_lnT_00) = s%P(k) * s% chiT_for_partials(k)
         else
            p_00%d1Array(i_lnT_00) = s%P(k) * s% chiT_for_partials(k)/s% T(k)
         end if
      end function wrap_p_00

      function wrap_p_p1(s, k) result(p_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: p_p1
         integer, intent(in) :: k
         p_p1 = 0d0 
         if (k < s%nz) then
            p_p1%val = s% P(k+1)
            if (s% solver_use_lnd) then
               p_p1%d1Array(i_lnd_p1) = s%P(k+1) * s% chiRho_for_partials(k+1)
            else
               p_p1%d1Array(i_lnd_p1) = s%P(k+1) * s% chiRho_for_partials(k+1)/s% rho(k+1)
            end if
            if (s% solver_use_lnT) then
               p_p1%d1Array(i_lnT_p1) = s%P(k+1) * s% chiT_for_partials(k+1)
            else
               p_p1%d1Array(i_lnT_p1) = s%P(k+1) * s% chiT_for_partials(k+1)/s% T(k+1)
            end if
         end if   
      end function wrap_p_p1

      function wrap_lnP_m1(s, k) result(lnP_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnP_m1
         integer, intent(in) :: k
         lnP_m1 = 0d0 
         if (k > 1) then
            lnP_m1%val = s% lnP(k-1)
            if (s% solver_use_lnd) then
               lnP_m1%d1Array(i_lnd_m1) = s% chiRho_for_partials(k-1)
            else
               lnP_m1%d1Array(i_lnd_m1) = s% chiRho_for_partials(k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               lnP_m1%d1Array(i_lnT_m1) = s% chiT_for_partials(k-1)
            else
               lnP_m1%d1Array(i_lnT_m1) = s% chiT_for_partials(k-1)/s% T(k-1)
            end if
         end if   
      end function wrap_lnP_m1

      function wrap_lnP_00(s, k) result(lnP_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnP_00
         integer, intent(in) :: k
         lnP_00 = 0d0 
         lnP_00%val = s% lnP(k)
         if (s% solver_use_lnd) then
            lnP_00%d1Array(i_lnd_00) = s% chiRho_for_partials(k)
         else
            lnP_00%d1Array(i_lnd_00) = s% chiRho_for_partials(k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            lnP_00%d1Array(i_lnT_00) = s% chiT_for_partials(k)
         else
            lnP_00%d1Array(i_lnT_00) = s% chiT_for_partials(k)/s% T(k)
         end if
      end function wrap_lnP_00

      function wrap_lnP_p1(s, k) result(lnP_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnP_p1
         integer, intent(in) :: k
         lnP_p1 = 0d0 
         if (k < s%nz) then
            lnP_p1%val = s% lnP(k+1)
            if (s% solver_use_lnd) then
               lnP_p1%d1Array(i_lnd_p1) = s% chiRho_for_partials(k+1)
            else
               lnP_p1%d1Array(i_lnd_p1) = s% chiRho_for_partials(k+1)/s% rho(k+1)
            end if
            if (s% solver_use_lnT) then
               lnP_p1%d1Array(i_lnT_p1) = s% chiT_for_partials(k+1)
            else
               lnP_p1%d1Array(i_lnT_p1) = s% chiT_for_partials(k+1)/s% T(k+1)
            end if
         end if   
      end function wrap_lnP_p1

      function wrap_ChiRho_m1(s, k) result(ChiRho_m1)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiRho_m1
         integer, intent(in) :: k
         ChiRho_m1 = 0d0 
         if (k > 1) then
            ChiRho_m1%val = s% ChiRho(k-1)
            if (s% solver_use_lnd) then
               ChiRho_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiRho,k-1)
            else
               ChiRho_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiRho,k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               ChiRho_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiRho,k-1)
            else
               ChiRho_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiRho,k-1)/s% T(k-1)
            end if
         end if
      end function wrap_ChiRho_m1

      function wrap_ChiRho_00(s, k) result(ChiRho_00)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiRho_00
         integer, intent(in) :: k
         ChiRho_00 = 0d0 
         ChiRho_00%val = s% ChiRho(k)
         if (s% solver_use_lnd) then
            ChiRho_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiRho,k)
         else
            ChiRho_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiRho,k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            ChiRho_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiRho,k)
         else
            ChiRho_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiRho,k)/s% T(k)
         end if
      end function wrap_ChiRho_00

      function wrap_ChiRho_p1(s, k) result(ChiRho_p1)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiRho_p1
         integer, intent(in) :: k
         ChiRho_p1 = 0d0 
         if (k < s% nz) then
            ChiRho_p1%val = s% ChiRho(k+1)
            if (s% solver_use_lnd) then
               ChiRho_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiRho,k+1)
            else
               ChiRho_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiRho,k+1)/s% rho(k+1)
            end if
            if (s% solver_use_lnT) then
               ChiRho_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiRho,k+1)
            else
               ChiRho_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiRho,k+1)/s% T(k+1)
            end if
         end if
      end function wrap_ChiRho_p1

      function wrap_ChiT_m1(s, k) result(ChiT_m1)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiT_m1
         integer, intent(in) :: k
         ChiT_m1 = 0d0 
         if (k > 1) then
            ChiT_m1%val = s% ChiT(k-1)
            if (s% solver_use_lnd) then
               ChiT_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiT,k-1)
            else
               ChiT_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiT,k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               ChiT_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiT,k-1)
            else
               ChiT_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiT,k-1)/s% T(k-1)
            end if
         end if
      end function wrap_ChiT_m1

      function wrap_ChiT_00(s, k) result(ChiT_00)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiT_00
         integer, intent(in) :: k
         ChiT_00 = 0d0 
         ChiT_00%val = s% ChiT(k)
         if (s% solver_use_lnd) then
            ChiT_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiT,k)
         else
            ChiT_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiT,k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            ChiT_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiT,k)
         else
            ChiT_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiT,k)/s% T(k)
         end if
      end function wrap_ChiT_00

      function wrap_ChiT_p1(s, k) result(ChiT_p1)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiT_p1
         integer, intent(in) :: k
         ChiT_p1 = 0d0 
         if (k < s% nz) then
            ChiT_p1%val = s% ChiT(k+1)
            if (s% solver_use_lnd) then
               ChiT_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiT,k+1)
            else
               ChiT_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiT,k+1)/s% rho(k+1)
            end if
            if (s% solver_use_lnT) then
               ChiT_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiT,k+1)
            else
               ChiT_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiT,k+1)/s% T(k+1)
            end if
         end if
      end function wrap_ChiT_p1

      function wrap_Cp_m1(s, k) result(Cp_m1)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Cp_m1
         integer, intent(in) :: k
         Cp_m1 = 0d0 
         if (k > 1) then
            Cp_m1%val = s% Cp(k-1)
            if (s% solver_use_lnd) then
               Cp_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_Cp,k-1)
            else
               Cp_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_Cp,k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               Cp_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_Cp,k-1)
            else
               Cp_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_Cp,k-1)/s% T(k-1)
            end if
         end if   
      end function wrap_Cp_m1

      function wrap_Cp_00(s, k) result(Cp_00)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Cp_00
         integer, intent(in) :: k
         Cp_00 = 0d0 
         Cp_00%val = s% Cp(k)
         if (s% solver_use_lnd) then
            Cp_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_Cp,k)
         else
            Cp_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_Cp,k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            Cp_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_Cp,k)
         else
            Cp_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_Cp,k)/s% T(k)
         end if
      end function wrap_Cp_00

      function wrap_Cp_p1(s, k) result(Cp_p1)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Cp_p1
         integer, intent(in) :: k
         Cp_p1 = 0d0 
         if (k < s% nz) then
            Cp_p1%val = s% Cp(k+1)
            if (s% solver_use_lnd) then
               Cp_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_Cp,k+1)
            else
               Cp_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_Cp,k+1)/s% rho(k+1)
            end if
            if (s% solver_use_lnT) then
               Cp_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_Cp,k+1)
            else
               Cp_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_Cp,k+1)/s% T(k+1)
            end if
         end if   
      end function wrap_Cp_p1

      function wrap_gamma1_m1(s, k) result(gamma1_m1)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: gamma1_m1
         integer, intent(in) :: k
         gamma1_m1 = 0d0 
         if (k > 1) then
            gamma1_m1%val = s% gamma1(k-1)
            if (s% solver_use_lnd) then
               gamma1_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_gamma1,k-1)
            else
               gamma1_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_gamma1,k-1)/s% rho(k-1)
            end if
            if (s% solver_use_lnT) then
               gamma1_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_gamma1,k-1)
            else
               gamma1_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_gamma1,k-1)/s% T(k-1)
            end if
         end if   
      end function wrap_gamma1_m1

      function wrap_gamma1_00(s, k) result(gamma1_00)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: gamma1_00
         integer, intent(in) :: k
         gamma1_00 = 0d0 
         gamma1_00%val = s% gamma1(k)
         if (s% solver_use_lnd) then
            gamma1_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_gamma1,k)
         else
            gamma1_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_gamma1,k)/s% rho(k)
         end if
         if (s% solver_use_lnT) then
            gamma1_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_gamma1,k)
         else
            gamma1_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_gamma1,k)/s% T(k)
         end if
      end function wrap_gamma1_00

      function wrap_gamma1_p1(s, k) result(gamma1_p1)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: gamma1_p1
         integer, intent(in) :: k
         gamma1_p1 = 0d0 
         if (k < s% nz) then
            gamma1_p1%val = s% gamma1(k+1)
            if (s% solver_use_lnd) then
               gamma1_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_gamma1,k+1)
            else
               gamma1_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_gamma1,k+1)/s% rho(k+1)
            end if
            if (s% solver_use_lnT) then
               gamma1_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_gamma1,k+1)
            else
               gamma1_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_gamma1,k+1)/s% T(k+1)
            end if
         end if   
      end function wrap_gamma1_p1

      function wrap_L_m1(s, k) result(L_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: L_m1
         integer, intent(in) :: k
         L_m1 = 0d0 
         if (k > 1) then
            L_m1 % val = s%L(k-1)
            L_m1 % d1Array(i_L_m1) = 1d0
         end if
      end function wrap_L_m1

      function wrap_L_00(s, k) result(L_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: L_00
         integer, intent(in) :: k
         L_00 = 0d0 
         L_00 % val = s%L(k)
         L_00 % d1Array(i_L_00) = 1d0
      end function wrap_L_00

      function wrap_L_p1(s, k) result(L_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: L_p1
         integer, intent(in) :: k
         L_p1 = 0d0 
         if (k < s%nz) then
            L_p1 % val = s%L(k+1)
            L_p1 % d1Array(i_L_p1) = 1d0
         else
            L_p1 %val = s%L_center
            ! L_center is a constant
         end if
      end function wrap_L_p1

      function wrap_r_m1(s, k) result(r_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: r_m1
         integer, intent(in) :: k
         r_m1 = 0d0 
         if (k > 1) then
            r_m1 % val = s%r(k-1)
            if (s% solver_use_lnR) then
               r_m1 % d1Array(i_lnR_m1) = s%r(k-1)
            else
               r_m1 % d1Array(i_lnR_m1) = 1d0
            end if
         end if
      end function wrap_r_m1

      function wrap_r_00(s, k) result(r_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: r_00
         integer, intent(in) :: k
         r_00 = 0d0 
         r_00 % val = s%r(k)
         if (s% solver_use_lnR) then
            r_00 % d1Array(i_lnR_00) = s%r(k)
         else
            r_00 % d1Array(i_lnR_00) = 1d0
         end if
      end function wrap_r_00

      function wrap_r_p1(s, k) result(r_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: r_p1
         integer, intent(in) :: k
         r_p1 = 0d0 
         if (k < s%nz) then
            r_p1 % val = s%r(k+1)
            if (s% solver_use_lnR) then
               r_p1 % d1Array(i_lnR_p1) = s%r(k+1)
            else
               r_p1 % d1Array(i_lnR_p1) = 1d0
            end if
         else
            r_p1 % val = s%r_center
         end if
      end function wrap_r_p1

      function wrap_lnR_m1(s, k) result(lnR_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnR_m1
         integer, intent(in) :: k
         lnR_m1 = 0d0 
         if (k > 1) then
            lnR_m1 % val = s%lnR(k-1)
            if (s% solver_use_lnR) then
               lnR_m1 % d1Array(i_lnR_m1) = 1d0
            else
               lnR_m1 % d1Array(i_lnR_m1) = 1d0/s% r(k-1)
            end if
         end if
      end function wrap_lnR_m1

      function wrap_lnR_00(s, k) result(lnR_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnR_00
         integer, intent(in) :: k
         lnR_00 = 0d0 
         lnR_00 % val = s%lnR(k)
         if (s% solver_use_lnR) then
            lnR_00 % d1Array(i_lnR_00) = 1d0
         else
            lnR_00 % d1Array(i_lnR_00) = 1d0/s% r(k)
         end if
      end function wrap_lnR_00

      function wrap_lnR_p1(s, k) result(lnR_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnR_p1
         integer, intent(in) :: k
         lnR_p1 = 0d0 
         if (k < s%nz) then
            lnR_p1 % val = s%lnR(k+1)
            if (s% solver_use_lnR) then
               lnR_p1 % d1Array(i_lnR_p1) = 1d0
            else
               lnR_p1 % d1Array(i_lnR_p1) = 1d0/s% r(k+1)
            end if
         else
            lnR_p1 % val = log(max(1d0,s%r_center))
         end if
      end function wrap_lnR_p1

      function wrap_dxh_lnR(s, k) result(dxh_lnR)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: dxh_lnR
         integer, intent(in) :: k
         dxh_lnR = 0d0 
         dxh_lnR % val = s%dxh_lnR(k) 
         dxh_lnR % d1Array(i_lnR_00) = 1d0
      end function wrap_dxh_lnR

      function wrap_v_m1(s, k) result(v_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_m1
         integer, intent(in) :: k
         v_m1 = 0d0 
         if (k > 1) then
            v_m1 % val = s%v(k-1)
            v_m1 % d1Array(i_v_m1) = 1d0
         end if
      end function wrap_v_m1

      function wrap_v_00(s, k) result(v_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_00
         integer, intent(in) :: k
         v_00 = 0d0 
         v_00 % val = s%v(k)
         v_00 % d1Array(i_v_00) = 1d0
      end function wrap_v_00

      function wrap_v_p1(s, k) result(v_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_p1
         integer, intent(in) :: k
         v_p1 = 0d0 
         if (k < s%nz) then
            v_p1 % val = s%v(k+1)
            v_p1 % d1Array(i_v_p1) = 1d0
         else
            v_p1 % val = s%v_center
            ! v_center is a constant
         end if
      end function wrap_v_p1
      
      function wrap_opt_time_center_r_m1(s, k) result(r_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: r_tc
         r_tc = wrap_r_m1(s,k)
         if (s% using_velocity_time_centering) then
            if (k > 1) r_tc = 0.5d0*(r_tc + s% r_start(k-1))
         end if
      end function wrap_opt_time_center_r_m1
            
      function wrap_opt_time_center_r_00(s, k) result(r_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: r_tc
         r_tc = wrap_r_00(s,k)
         if (s% using_velocity_time_centering) &
            r_tc = 0.5d0*(r_tc + s% r_start(k))
      end function wrap_opt_time_center_r_00
            
      function wrap_opt_time_center_r_p1(s, k) result(r_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: r_tc
         r_tc = wrap_r_p1(s,k)
         if (s% using_velocity_time_centering) then
            if (k < s% nz) then
               r_tc = 0.5d0*(r_tc + s% r_start(k+1))
            else
               r_tc = 0.5d0*(r_tc + s% r_center)
            end if
         end if
      end function wrap_opt_time_center_r_p1      
      
      function wrap_opt_time_center_v_m1(s, k) result(v_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: v_tc
         v_tc = wrap_v_m1(s,k)
         if (s% using_velocity_time_centering) then
            if (s% v_flag) then
               if (k > 1) then
                  v_tc = 0.5d0*(v_tc + s% v_start(k-1))
               end if
            else
               stop 'fix wrap_opt_time_center_v'
            end if
         end if
      end function wrap_opt_time_center_v_m1      
      
      function wrap_opt_time_center_v_00(s, k) result(v_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: v_tc
         v_tc = wrap_v_00(s,k)
         if (s% using_velocity_time_centering) then
            if (s% v_flag) then
               v_tc = 0.5d0*(v_tc + s% v_start(k))
            else
               stop 'fix wrap_opt_time_center_v'
            end if
         end if
      end function wrap_opt_time_center_v_00     
      
      function wrap_opt_time_center_v_p1(s, k) result(v_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: v_tc
         v_tc = wrap_v_p1(s,k)
         if (s% using_velocity_time_centering) then
            if (s% v_flag) then
               if (k < s% nz) then
                  v_tc = 0.5d0*(v_tc + s% v_start(k+1))
               else
                  v_tc = 0.5d0*(v_tc + s% v_center)
               end if
            else
               stop 'fix wrap_opt_time_center_v'
            end if
         end if
      end function wrap_opt_time_center_v_p1

      function wrap_u_m1(s, k) result(v_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_m1
         integer, intent(in) :: k
         v_m1 = 0d0 
         if (k > 1) then
            v_m1 % val = s%u(k-1)
            v_m1 % d1Array(i_v_m1) = 1d0
         end if
      end function wrap_u_m1

      function wrap_u_00(s, k) result(v_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_00
         integer, intent(in) :: k
         v_00 = 0d0 
         v_00 % val = s%u(k)
         v_00 % d1Array(i_v_00) = 1d0
      end function wrap_u_00

      function wrap_u_p1(s, k) result(v_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_p1
         integer, intent(in) :: k
         v_p1 = 0d0 
         if (k < s%nz) then
            v_p1 % val = s%u(k+1)
            v_p1 % d1Array(i_v_p1) = 1d0
         else
            v_p1 % val = 0d0
            ! v_center is a constant
         end if
      end function wrap_u_p1

      function wrap_Hp_m1(s, k) result(Hp_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Hp_m1
         integer, intent(in) :: k
         Hp_m1 = 0d0 
         if (k > 1) then
            Hp_m1 % val = s% Hp_face(k-1)
            Hp_m1 % d1Array(i_Hp_m1) = 1d0
         end if            
      end function wrap_Hp_m1

      function wrap_Hp_00(s, k) result(Hp_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Hp_00
         integer, intent(in) :: k
         Hp_00 = 0d0 
         Hp_00 % val = s% Hp_face(k)
         Hp_00 % d1Array(i_Hp_00) = 1d0
      end function wrap_Hp_00

      function wrap_Hp_p1(s, k) result(Hp_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Hp_p1
         integer, intent(in) :: k
         Hp_p1 = 0d0 
         if (k < s%nz) then
            Hp_p1 % val = s% Hp_face(k+1)
            Hp_p1 % d1Array(i_Hp_p1) = 1d0
         end if
      end function wrap_Hp_p1


      function wrap_xtra1_m1(s, k) result(xtra1_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra1_m1
         integer, intent(in) :: k
         xtra1_m1 = 0d0 
         if (k > 1) then ! s%w(k-1)
            xtra1_m1 % d1Array(i_xtra1_m1) = 1d0
         end if            
      end function wrap_xtra1_m1

      function wrap_xtra1_00(s, k) result(xtra1_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra1_00
         integer, intent(in) :: k
         xtra1_00 = 0d0 
         xtra1_00 % val = 0d0
         xtra1_00 % d1Array(i_xtra1_00) = 1d0
      end function wrap_xtra1_00

      function wrap_xtra1_p1(s, k) result(xtra1_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra1_p1
         integer, intent(in) :: k
         xtra1_p1 = 0d0 
         if (k < s%nz) then
            xtra1_p1 % val = 0d0 
            xtra1_p1 % d1Array(i_xtra1_p1) = 1d0
         end if
      end function wrap_xtra1_p1


      function wrap_xtra2_m1(s, k) result(xtra2_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra2_m1
         integer, intent(in) :: k
         xtra2_m1 = 0d0 
         if (k > 1) then
            xtra2_m1 % val = 0d0
            xtra2_m1 % d1Array(i_xtra2_m1) = 1d0
         end if            
      end function wrap_xtra2_m1

      function wrap_xtra2_00(s, k) result(xtra2_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra2_00
         integer, intent(in) :: k
         xtra2_00 = 0d0 
         xtra2_00 % val = 0d0 
         xtra2_00 % d1Array(i_xtra2_00) = 1d0
      end function wrap_xtra2_00

      function wrap_xtra2_p1(s, k) result(xtra2_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra2_p1
         integer, intent(in) :: k
         xtra2_p1 = 0d0 
         if (k < s%nz) then
            xtra2_p1 % val = 0d0 
            xtra2_p1 % d1Array(i_xtra2_p1) = 1d0
         end if
      end function wrap_xtra2_p1


end module auto_diff_support
