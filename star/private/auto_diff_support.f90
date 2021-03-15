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
            dw_m1, dw_00, dw_p1, &
            dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, &
            dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dxtra3_m1, dxtra3_00, dxtra3_p1)
         type(auto_diff_real_star_order1), intent(in) :: var
         real(dp), intent(out) :: &
            val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dxtra3_m1, dxtra3_00, dxtra3_p1
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
         dxtra1_m1 = var%d1Array(i_xtra1_m1)
         dxtra1_00 = var%d1Array(i_xtra1_00)
         dxtra1_p1 = var%d1Array(i_xtra1_p1)
         dxtra2_m1 = var%d1Array(i_xtra2_m1)
         dxtra2_00 = var%d1Array(i_xtra2_00)
         dxtra2_p1 = var%d1Array(i_xtra2_p1)
         dxtra3_m1 = var%d1Array(i_xtra3_m1)
         dxtra3_00 = var%d1Array(i_xtra3_00)
         dxtra3_p1 = var%d1Array(i_xtra3_p1)
      end subroutine unwrap

      subroutine wrap(var, val, &
            dlnd_m1, dlnd_00, dlnd_p1, &
            dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, &
            dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, &
            dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dxtra3_m1, dxtra3_00, dxtra3_p1)
         type(auto_diff_real_star_order1), intent(out) :: var
         real(dp), intent(in) :: &
            val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            dw_m1, dw_00, dw_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1, &
            dxtra1_m1, dxtra1_00, dxtra1_p1, &
            dxtra2_m1, dxtra2_00, dxtra2_p1, &
            dxtra3_m1, dxtra3_00, dxtra3_p1
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
         var%d1Array(i_xtra1_m1) = dxtra1_m1
         var%d1Array(i_xtra1_00) = dxtra1_00
         var%d1Array(i_xtra1_p1) = dxtra1_p1
         var%d1Array(i_xtra2_m1) = dxtra2_m1
         var%d1Array(i_xtra2_00) = dxtra2_00
         var%d1Array(i_xtra2_p1) = dxtra2_p1
         var%d1Array(i_xtra3_m1) = dxtra3_m1
         var%d1Array(i_xtra3_00) = dxtra3_00
         var%d1Array(i_xtra3_p1) = dxtra3_p1
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
            T_m1 % d1Array(i_lnT_m1) = s%T(k-1)
         end if
      end function wrap_T_m1

      function wrap_T_00(s, k) result(T_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: T_00
         integer, intent(in) :: k
         T_00 = 0d0 
         T_00 % val = s%T(k)
         T_00 % d1Array(i_lnT_00) = s%T(k)
      end function wrap_T_00

      function wrap_T_p1(s, k) result(T_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: T_p1
         integer, intent(in) :: k
         T_p1 = 0d0 
         if (k < s%nz) then
            T_p1 % val = s%T(k+1)
            T_p1 % d1Array(i_lnT_p1) = s%T(k+1)
         end if
      end function wrap_T_p1

      function wrap_lnT_m1(s, k) result(lnT_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnT_m1
         integer, intent(in) :: k
         lnT_m1 = 0d0 
         if (k > 1) then
            lnT_m1 % val = s%lnT(k-1)
            lnT_m1 % d1Array(i_lnT_m1) = 1d0
         end if
      end function wrap_lnT_m1

      function wrap_lnT_00(s, k) result(lnT_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnT_00
         integer, intent(in) :: k
         lnT_00 = 0d0 
         lnT_00 % val = s%lnT(k)
         lnT_00 % d1Array(i_lnT_00) = 1d0
      end function wrap_lnT_00

      function wrap_lnT_p1(s, k) result(lnT_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnT_p1
         integer, intent(in) :: k
         lnT_p1 = 0d0 
         if (k < s%nz) then
            lnT_p1 % val = s%lnT(k+1)
            lnT_p1 % d1Array(i_lnT_p1) = 1d0
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
            d_m1 % d1Array(i_lnd_m1) = s%rho(k-1)
         end if
      end function wrap_d_m1

      function wrap_d_00(s, k) result(d_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: d_00
         integer, intent(in) :: k
         d_00 = 0d0 
         d_00 % val = s%rho(k)
         d_00 % d1Array(i_lnd_00) = s%rho(k)
      end function wrap_d_00

      function wrap_d_p1(s, k) result(d_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: d_p1
         integer, intent(in) :: k
         d_p1 = 0d0 
         if (k < s%nz) then
            d_p1 % val = s%rho(k+1)
            d_p1 % d1Array(i_lnd_p1) = s%rho(k+1)
         end if
      end function wrap_d_p1

      function wrap_lnd_m1(s, k) result(lnd_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnd_m1
         integer, intent(in) :: k
         lnd_m1 = 0d0 
         if (k > 1) then
            lnd_m1 % val = s%lnd(k-1)
            lnd_m1 % d1Array(i_lnd_m1) = 1d0
         end if
      end function wrap_lnd_m1

      function wrap_lnd_00(s, k) result(lnd_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnd_00
         integer, intent(in) :: k
         lnd_00 = 0d0 
         lnd_00 % val = s%lnd(k)
         lnd_00 % d1Array(i_lnd_00) = 1d0
      end function wrap_lnd_00

      function wrap_lnd_p1(s, k) result(lnd_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnd_p1
         integer, intent(in) :: k
         lnd_p1 = 0d0 
         if (k < s%nz) then
            lnd_p1 % val = s%lnd(k+1)
            lnd_p1 % d1Array(i_lnd_p1) = 1d0
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

      function wrap_w_m1(s, k) result(w_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: w_m1
         integer, intent(in) :: k
         w_m1 = 0d0 
         if (k > 1) then
            w_m1 % val = s%w(k-1)
            w_m1 % d1Array(i_w_m1) = 1d0
         end if            
      end function wrap_w_m1

      function wrap_w_00(s, k) result(w_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: w_00
         integer, intent(in) :: k
         w_00 = 0d0 
         w_00 % val = s%w(k)
         w_00 % d1Array(i_w_00) = 1d0
      end function wrap_w_00

      function wrap_w_p1(s, k) result(w_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: w_p1
         integer, intent(in) :: k
         w_p1 = 0d0 
         if (k < s%nz) then
            w_p1 % val = s%w(k+1)
            w_p1 % d1Array(i_w_p1) = 1d0
         end if
      end function wrap_w_p1

      function wrap_dxh_w(s, k) result(dxh_w)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: dxh_w
         integer, intent(in) :: k
         dxh_w = 0d0 
         dxh_w % val = s%dxh_w(k)
         dxh_w % d1Array(i_w_00) = 1d0
      end function wrap_dxh_w
      
      real(dp) function get_etrb(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         get_etrb = pow2(s% w(k))
      end function get_etrb
      
      real(dp) function get_etrb_start(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         get_etrb_start = pow2(s% w_start(k))
      end function get_etrb_start
      
      real(dp) function get_w(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         get_w = s% w(k)
      end function get_w
      
      real(dp) function get_w_start(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         get_w_start = s% w_start(k)
      end function get_w_start

      function wrap_etrb_m1(s, k) result(etrb_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: etrb_m1
         integer, intent(in) :: k
         etrb_m1 = pow2(wrap_w_m1(s,k)) 
      end function wrap_etrb_m1

      function wrap_etrb_00(s, k) result(etrb_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: etrb_00
         integer, intent(in) :: k
         etrb_00 = pow2(wrap_w_00(s,k)) 
      end function wrap_etrb_00

      function wrap_etrb_p1(s, k) result(etrb_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: etrb_p1
         integer, intent(in) :: k
         etrb_p1 = pow2(wrap_w_p1(s,k)) 
      end function wrap_etrb_p1

      function wrap_kap_m1(s, k) result(kap_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: kap_m1
         integer, intent(in) :: k
         kap_m1 = 0d0 
         if (k > 1) then
            kap_m1 % val = s%opacity(k-1)
            kap_m1 % d1Array(i_lnd_m1) = s%d_opacity_dlnd(k-1)
            kap_m1 % d1Array(i_lnT_m1) = s%d_opacity_dlnT(k-1)
         end if
      end function wrap_kap_m1

      function wrap_kap_00(s, k) result(kap_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: kap_00
         integer, intent(in) :: k
         kap_00 = 0d0 
         kap_00 % val = s%opacity(k)
         kap_00 % d1Array(i_lnd_00) = s%d_opacity_dlnd(k)
         kap_00 % d1Array(i_lnT_00) = s%d_opacity_dlnT(k)
      end function wrap_kap_00

      function wrap_kap_p1(s, k) result(kap_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: kap_p1
         integer, intent(in) :: k
         kap_p1 = 0d0 
         if (k < s%nz) then
            kap_p1 % val = s%opacity(k+1)
            kap_p1 % d1Array(i_lnd_p1) = s%d_opacity_dlnd(k+1)/s% rho(k+1)
            kap_p1 % d1Array(i_lnT_p1) = s%d_opacity_dlnT(k+1)
         end if
      end function wrap_kap_p1

      function wrap_s_m1(s, k) result(s_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: s_m1
         integer, intent(in) :: k
         s_m1 = 0d0 
         if (k > 1) then
            s_m1%val = s% entropy(k-1)
            s_m1%d1Array(i_lnd_m1) = s% dS_dRho_for_partials(k-1)*s% rho(k-1)
            s_m1%d1Array(i_lnT_m1) = s% dS_dT_for_partials(k-1)*s% T(k-1)
         end if   
      end function wrap_s_m1

      function wrap_s_00(s, k) result(s_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: s_00
         integer, intent(in) :: k
         s_00 = 0d0 
         s_00%val = s% entropy(k)
         s_00%d1Array(i_lnd_00) = s% dS_dRho_for_partials(k)*s% rho(k)
         s_00%d1Array(i_lnT_00) = s% dS_dT_for_partials(k)*s% T(k)
      end function wrap_s_00

      function wrap_s_p1(s, k) result(s_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: s_p1
         integer, intent(in) :: k
         s_p1 = 0d0 
         if (k < s%nz) then
            s_p1%val = s% entropy(k+1)
            s_p1%d1Array(i_lnd_p1) = s% dS_dRho_for_partials(k+1)*s% rho(k+1)
            s_p1%d1Array(i_lnT_p1) = s% dS_dT_for_partials(k+1)*s% T(k+1)
         end if
      end function wrap_s_p1

      function wrap_e_m1(s, k) result(e_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: e_m1
         integer, intent(in) :: k
         e_m1 = 0d0 
         if (k > 1) then
            e_m1%val = s% energy(k-1)
            e_m1%d1Array(i_lnd_m1) = s% dE_dRho_for_partials(k-1)*s% rho(k-1)
            e_m1%d1Array(i_lnT_m1) = s% Cv_for_partials(k-1)*s% T(k-1)
         end if   
      end function wrap_e_m1

      function wrap_e_00(s, k) result(e_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: e_00
         integer, intent(in) :: k
         e_00 = 0d0 
         e_00%val = s% energy(k)
         e_00%d1Array(i_lnd_00) = s% dE_dRho_for_partials(k)*s% rho(k)
         e_00%d1Array(i_lnT_00) = s% Cv_for_partials(k)*s% T(k)
      end function wrap_e_00

      function wrap_e_p1(s, k) result(e_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: e_p1
         integer, intent(in) :: k
         e_p1 = 0d0 
         if (k < s%nz) then
            e_p1%val = s% energy(k+1)
            e_p1%d1Array(i_lnd_p1) = s% dE_dRho_for_partials(k+1)*s% rho(k+1)
            e_p1%d1Array(i_lnT_p1) = s% Cv_for_partials(k+1)*s% T(k+1)
         end if
      end function wrap_e_p1

      function wrap_Peos_m1(s, k) result(Peos_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Peos_m1
         integer, intent(in) :: k
         Peos_m1 = 0d0 
         if (k > 1) then
            Peos_m1%val = s% Peos(k-1)
            Peos_m1%d1Array(i_lnd_m1) = s%Peos(k-1) * s% chiRho_for_partials(k-1)
            Peos_m1%d1Array(i_lnT_m1) = s%Peos(k-1) * s% chiT_for_partials(k-1)
         end if   
      end function wrap_Peos_m1

      function wrap_Peos_00(s, k) result(Peos_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Peos_00
         integer, intent(in) :: k
         Peos_00 = 0d0 
         Peos_00%val = s% Peos(k)
         Peos_00%d1Array(i_lnd_00) = s%Peos(k) * s% chiRho_for_partials(k)
         Peos_00%d1Array(i_lnT_00) = s%Peos(k) * s% chiT_for_partials(k)
      end function wrap_Peos_00

      function wrap_Peos_p1(s, k) result(Peos_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Peos_p1
         integer, intent(in) :: k
         Peos_p1 = 0d0 
         if (k < s%nz) then
            Peos_p1%val = s% Peos(k+1)
            Peos_p1%d1Array(i_lnd_p1) = s%Peos(k+1) * s% chiRho_for_partials(k+1)
            Peos_p1%d1Array(i_lnT_p1) = s%Peos(k+1) * s% chiT_for_partials(k+1)
         end if   
      end function wrap_Peos_p1

      function wrap_lnPeos_m1(s, k) result(lnPeos_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnPeos_m1
         integer, intent(in) :: k
         lnPeos_m1 = 0d0 
         if (k > 1) then
            lnPeos_m1%val = s% lnPeos(k-1)
            lnPeos_m1%d1Array(i_lnd_m1) = s% chiRho_for_partials(k-1)
            lnPeos_m1%d1Array(i_lnT_m1) = s% chiT_for_partials(k-1)
         end if   
      end function wrap_lnPeos_m1

      function wrap_lnPeos_00(s, k) result(lnPeos_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnPeos_00
         integer, intent(in) :: k
         lnPeos_00 = 0d0 
         lnPeos_00%val = s% lnPeos(k)
         lnPeos_00%d1Array(i_lnd_00) = s% chiRho_for_partials(k)
         lnPeos_00%d1Array(i_lnT_00) = s% chiT_for_partials(k)
      end function wrap_lnPeos_00

      function wrap_lnPeos_p1(s, k) result(lnPeos_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnPeos_p1
         integer, intent(in) :: k
         lnPeos_p1 = 0d0 
         if (k < s%nz) then
            lnPeos_p1%val = s% lnPeos(k+1)
            lnPeos_p1%d1Array(i_lnd_p1) = s% chiRho_for_partials(k+1)
            lnPeos_p1%d1Array(i_lnT_p1) = s% chiT_for_partials(k+1)
         end if   
      end function wrap_lnPeos_p1

      function wrap_ChiRho_m1(s, k) result(ChiRho_m1)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiRho_m1
         integer, intent(in) :: k
         ChiRho_m1 = 0d0 
         if (k > 1) then
            ChiRho_m1%val = s% ChiRho(k-1)
            ChiRho_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiRho,k-1)
            ChiRho_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiRho,k-1)
         end if
      end function wrap_ChiRho_m1

      function wrap_ChiRho_00(s, k) result(ChiRho_00)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiRho_00
         integer, intent(in) :: k
         ChiRho_00 = 0d0 
         ChiRho_00%val = s% ChiRho(k)
         ChiRho_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiRho,k)
         ChiRho_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiRho,k)
      end function wrap_ChiRho_00

      function wrap_ChiRho_p1(s, k) result(ChiRho_p1)
         use eos_def, only: i_ChiRho
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiRho_p1
         integer, intent(in) :: k
         ChiRho_p1 = 0d0 
         if (k < s% nz) then
            ChiRho_p1%val = s% ChiRho(k+1)
            ChiRho_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiRho,k+1)
            ChiRho_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiRho,k+1)
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
            ChiT_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_ChiT,k-1)
            ChiT_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_ChiT,k-1)
         end if
      end function wrap_ChiT_m1

      function wrap_ChiT_00(s, k) result(ChiT_00)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiT_00
         integer, intent(in) :: k
         ChiT_00 = 0d0 
         ChiT_00%val = s% ChiT(k)
         ChiT_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_ChiT,k)
         ChiT_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_ChiT,k)
      end function wrap_ChiT_00

      function wrap_ChiT_p1(s, k) result(ChiT_p1)
         use eos_def, only: i_ChiT
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: ChiT_p1
         integer, intent(in) :: k
         ChiT_p1 = 0d0 
         if (k < s% nz) then
            ChiT_p1%val = s% ChiT(k+1)
            ChiT_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_ChiT,k+1)
            ChiT_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_ChiT,k+1)
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
            Cp_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_Cp,k-1)
            Cp_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_Cp,k-1)
         end if   
      end function wrap_Cp_m1

      function wrap_Cp_00(s, k) result(Cp_00)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Cp_00
         integer, intent(in) :: k
         Cp_00 = 0d0 
         Cp_00%val = s% Cp(k)
         Cp_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_Cp,k)
         Cp_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_Cp,k)
      end function wrap_Cp_00

      function wrap_Cp_p1(s, k) result(Cp_p1)
         use eos_def, only: i_Cp
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: Cp_p1
         integer, intent(in) :: k
         Cp_p1 = 0d0 
         if (k < s% nz) then
            Cp_p1%val = s% Cp(k+1)
            Cp_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_Cp,k+1)
            Cp_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_Cp,k+1)
         end if   
      end function wrap_Cp_p1

      function wrap_grad_ad_m1(s, k) result(grad_ad_m1)
         use eos_def, only: i_grad_ad
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: grad_ad_m1
         integer, intent(in) :: k
         grad_ad_m1 = 0d0 
         if (k > 1) then
            grad_ad_m1%val = s% grada(k-1)
            grad_ad_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_grad_ad,k-1)
            grad_ad_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_grad_ad,k-1)
         end if   
      end function wrap_grad_ad_m1

      function wrap_grad_ad_00(s, k) result(grad_ad_00)
         use eos_def, only: i_grad_ad
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: grad_ad_00
         integer, intent(in) :: k
         grad_ad_00 = 0d0 
         grad_ad_00%val = s% grada(k)
         grad_ad_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_grad_ad,k)
         grad_ad_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_grad_ad,k)
      end function wrap_grad_ad_00

      function wrap_grad_ad_p1(s, k) result(grad_ad_p1)
         use eos_def, only: i_grad_ad
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: grad_ad_p1
         integer, intent(in) :: k
         grad_ad_p1 = 0d0 
         if (k < s% nz) then
            grad_ad_p1%val = s% grada(k+1)
            grad_ad_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_grad_ad,k+1)
            grad_ad_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_grad_ad,k+1)
         end if   
      end function wrap_grad_ad_p1

      function wrap_gamma1_m1(s, k) result(gamma1_m1)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: gamma1_m1
         integer, intent(in) :: k
         gamma1_m1 = 0d0 
         if (k > 1) then
            gamma1_m1%val = s% gamma1(k-1)
            gamma1_m1%d1Array(i_lnd_m1) = s% d_eos_dlnd(i_gamma1,k-1)
            gamma1_m1%d1Array(i_lnT_m1) = s% d_eos_dlnT(i_gamma1,k-1)
         end if   
      end function wrap_gamma1_m1

      function wrap_gamma1_00(s, k) result(gamma1_00)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: gamma1_00
         integer, intent(in) :: k
         gamma1_00 = 0d0 
         gamma1_00%val = s% gamma1(k)
         gamma1_00%d1Array(i_lnd_00) = s% d_eos_dlnd(i_gamma1,k)
         gamma1_00%d1Array(i_lnT_00) = s% d_eos_dlnT(i_gamma1,k)
      end function wrap_gamma1_00

      function wrap_gamma1_p1(s, k) result(gamma1_p1)
         use eos_def, only: i_gamma1
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: gamma1_p1
         integer, intent(in) :: k
         gamma1_p1 = 0d0 
         if (k < s% nz) then
            gamma1_p1%val = s% gamma1(k+1)
            gamma1_p1%d1Array(i_lnd_p1) = s% d_eos_dlnd(i_gamma1,k+1)
            gamma1_p1%d1Array(i_lnT_p1) = s% d_eos_dlnT(i_gamma1,k+1)
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
            r_m1 % d1Array(i_lnR_m1) = s%r(k-1)
         end if
      end function wrap_r_m1

      function wrap_r_00(s, k) result(r_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: r_00
         integer, intent(in) :: k
         r_00 = 0d0 
         r_00 % val = s%r(k)
         r_00 % d1Array(i_lnR_00) = s%r(k)
      end function wrap_r_00

      function wrap_r_p1(s, k) result(r_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: r_p1
         integer, intent(in) :: k
         r_p1 = 0d0 
         if (k < s%nz) then
            r_p1 % val = s%r(k+1)
            r_p1 % d1Array(i_lnR_p1) = s%r(k+1)
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
            lnR_m1 % d1Array(i_lnR_m1) = 1d0
         end if
      end function wrap_lnR_m1

      function wrap_lnR_00(s, k) result(lnR_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnR_00
         integer, intent(in) :: k
         lnR_00 = 0d0 
         lnR_00 % val = s%lnR(k)
         lnR_00 % d1Array(i_lnR_00) = 1d0
      end function wrap_lnR_00

      function wrap_lnR_p1(s, k) result(lnR_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: lnR_p1
         integer, intent(in) :: k
         lnR_p1 = 0d0 
         if (k < s%nz) then
            lnR_p1 % val = s%lnR(k+1)
            lnR_p1 % d1Array(i_lnR_p1) = 1d0
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

      function wrap_u_face_m1(s, k) result(v_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_m1
         integer, intent(in) :: k
         v_m1 = 0
         if (k > 1) v_m1 = shift_m1(s% u_face_ad(k-1))
      end function wrap_u_face_m1

      function wrap_u_face_00(s, k) result(v_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_00
         integer, intent(in) :: k
         v_00 = 0
         if (k > 1) v_00 = s% u_face_ad(k)
      end function wrap_u_face_00

      function wrap_u_face_p1(s, k) result(v_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_p1
         integer, intent(in) :: k
         v_p1 = 0
         if (k < s% nz) v_p1 = shift_p1(s% u_face_ad(k+1))
      end function wrap_u_face_p1

      function wrap_v_m1(s, k) result(v_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_m1
         integer, intent(in) :: k
         if (s% u_flag) then
            v_m1 = wrap_u_face_m1(s,k)
            return
         end if
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
         if (s% u_flag) then
            v_00 = wrap_u_face_00(s,k)
            return
         end if
         v_00 = 0d0 
         v_00 % val = s%v(k)
         v_00 % d1Array(i_v_00) = 1d0
      end function wrap_v_00

      function wrap_v_p1(s, k) result(v_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: v_p1
         integer, intent(in) :: k
         if (s% u_flag) then
            v_p1 = wrap_u_face_p1(s,k)
            return
         end if
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
         v_tc = 0
         if (k == 1) return
         if (s% v_flag) then
            v_tc = wrap_v_m1(s,k) 
            if (s% using_velocity_time_centering) &
               v_tc = 0.5d0*(v_tc + s% v_start(k-1))
         else if (s% u_flag) then
            v_tc = wrap_u_face_m1(s,k)
            if (s% using_velocity_time_centering) &
               v_tc = 0.5d0*(v_tc + s% u_face_start(k-1))
         end if
      end function wrap_opt_time_center_v_m1      
      
      function wrap_opt_time_center_v_00(s, k) result(v_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: v_tc
         v_tc = 0
         if (s% v_flag) then
            v_tc = wrap_v_00(s,k)
            if (s% using_velocity_time_centering) &
               v_tc = 0.5d0*(v_tc + s% v_start(k))
         else if (s% u_flag) then
            v_tc = wrap_u_face_00(s,k)
            if (s% using_velocity_time_centering) &
               v_tc = 0.5d0*(v_tc + s% u_face_start(k))
         end if
      end function wrap_opt_time_center_v_00     
      
      function wrap_opt_time_center_v_p1(s, k) result(v_tc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: v_tc
         v_tc = 0
         if (k == s% nz) return
         if (s% v_flag) then
            v_tc = wrap_v_p1(s,k) 
            if (s% using_velocity_time_centering) &
               v_tc = 0.5d0*(v_tc + s% v_start(k+1))
         else if (s% u_flag) then
            v_tc = wrap_u_face_p1(s,k)
            if (s% using_velocity_time_centering) &
               v_tc = 0.5d0*(v_tc + s% u_face_start(k+1))
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


      function wrap_xtra1_m1(s, k) result(xtra1_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra1_m1
         integer, intent(in) :: k
         xtra1_m1 = 0d0 
         if (k > 1) then 
            xtra1_m1% val = s% xtra1_array(k)
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


      function wrap_xtra3_m1(s, k) result(xtra3_m1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra3_m1
         integer, intent(in) :: k
         xtra3_m1 = 0d0 
         if (k > 1) then
            xtra3_m1 % val = 0d0
            xtra3_m1 % d1Array(i_xtra3_m1) = 1d0
         end if            
      end function wrap_xtra3_m1

      function wrap_xtra3_00(s, k) result(xtra3_00)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra3_00
         integer, intent(in) :: k
         xtra3_00 = 0d0 
         xtra3_00 % val = 0d0
         xtra3_00 % d1Array(i_xtra3_00) = 1d0
      end function wrap_xtra3_00

      function wrap_xtra3_p1(s, k) result(xtra3_p1)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: xtra3_p1
         integer, intent(in) :: k
         xtra3_p1 = 0d0 
         if (k < s%nz) then
            xtra3_p1 % val = 0d0
            xtra3_p1 % d1Array(i_xtra3_p1) = 1d0
         end if
      end function wrap_xtra3_p1


end module auto_diff_support
