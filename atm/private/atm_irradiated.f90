! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
!
! ***********************************************************************

module atm_irradiated
   
   ! Uses
   
   use const_def
   use math_lib
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: eval_irradiated
   
   ! Procedures

contains
   
   ! Evaluate irradiated atmosphere data with fixed, uniform opacity
   
   subroutine eval_irradiated(&
      L, R, M, cgrav, T_eq, P_surf, kap_guess, kap_v, gamma, &
      eos_proc, kap_proc, errtol, max_iters, skip_partials, &
      Teff, kap, tau, &
      lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
      ierr)
      
      use atm_def, only : atm_eos_iface, atm_kap_iface
      use atm_utils, only : eval_Teff_g
      use eos_def, only : num_eos_basic_results, i_chiRho, i_chiT
      
      real(dp), intent(in) :: L
      real(dp), intent(in) :: R
      real(dp), intent(in) :: M
      real(dp), intent(in) :: cgrav
      real(dp), intent(in) :: T_eq
      real(dp), intent(in) :: P_surf
      real(dp), intent(in) :: kap_guess
      real(dp), intent(in) :: kap_v
      real(dp), intent(in) :: gamma
      procedure(atm_eos_iface) :: eos_proc
      procedure(atm_kap_iface) :: kap_proc
      real(dp), intent(in) :: errtol
      integer, intent(in) :: max_iters
      logical, intent(in) :: skip_partials
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: kap
      real(dp), intent(out) :: tau
      real(dp), intent(out) :: lnT
      real(dp), intent(out) :: dlnT_dL
      real(dp), intent(out) :: dlnT_dlnR
      real(dp), intent(out) :: dlnT_dlnM
      real(dp), intent(out) :: dlnT_dlnkap
      integer, intent(out) :: ierr
      
      real(dp) :: T_int
      real(dp) :: g
      real(dp) :: lnP
      integer :: iters
      real(dp) :: lnRho
      real(dp) :: res(num_eos_basic_results)
      real(dp) :: dres_dlnRho(num_eos_basic_results)
      real(dp) :: dres_dlnT(num_eos_basic_results)
      real(dp) :: dlnkap_dlnRho
      real(dp) :: dlnkap_dlnT
      real(dp) :: kap_prev
      real(dp) :: err
      real(dp) :: chiRho
      real(dp) :: chiT
      real(dp) :: dlnkap_dlnT_P
      
      include 'formats'
      
      ierr = 0
      
      ! Sanity checks
      
      if (L <= 0._dp .OR. R <= 0._dp .OR. M <= 0._dp) then
         ierr = -1
         return
      endif
      
      ! Evaluate the 'interior' temperature & gravity
      
      call eval_Teff_g(L, R, M, cgrav, T_int, g)
      
      ! Evaluate atmosphere data using kap_guess as the opacity
      
      kap = kap_guess
      
      call eval_data(&
         T_int, g, L, T_eq, P_surf, kap, kap_v, gamma, skip_partials, &
         tau, lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         ierr)
      
      ! Iterate to find a consistent opacity
      
      lnP = log(P_surf)
      
      iterate_loop : do iters = 1, max_iters
         
         ! Calculate the density & eos results
         
         call eos_proc(&
            lnP, lnT, &
            lnRho, res, dres_dlnRho, dres_dlnT, &
            ierr)
         if (ierr /= 0) then
            write(*, *) 'Call to eos_proc failed in eval_T_tau_uniform'
            return
         end if
         
         ! Update the opacity
         
         kap_prev = kap
         
         call kap_proc(&
            lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            ierr)
         if (ierr /= 0) then
            write(*, *) 'Call to kap_proc failed in eval_T_tau_uniform'
            return
         end if
         
         ! Check for convergence
         
         err = abs(kap_prev - kap) / (errtol + errtol * kap)
         
         if (err < 1._dp) exit iterate_loop
         
         kap = kap_prev + 0.5_dp * (kap - kap_prev) ! under correct
         
         ! Re-evaluate atmosphere data
         
         call eval_data(&
            T_int, g, L, T_eq, P_surf, kap, kap_v, gamma, skip_partials, &
            tau, lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            ierr)
      
      end do iterate_loop
      
      if (max_iters > 0 .AND. iters > max_iters) then
         write(*, *) 'Exceeded max_iters iterations in eval_irradiated'
         ierr = -1
         return
      end if
      
      ! If necessary, fix up the partials to account for the implicit
      ! dependence of the opacity on the final T
      
      if (max_iters > 0 .AND. .NOT. skip_partials) then
         
         chiRho = res(i_chiRho)
         chiT = res(i_chiT)
         
         dlnkap_dlnT_P = dlnkap_dlnT - dlnkap_dlnRho * chiT / chiRho
         
         dlnT_dL = dlnT_dL / (1._dp - dlnkap_dlnT_P * dlnT_dlnkap)
         dlnT_dlnR = dlnT_dlnR / (1._dp - dlnkap_dlnT_P * dlnT_dlnkap)
         dlnT_dlnM = dlnT_dlnM / (1._dp - dlnkap_dlnT_P * dlnT_dlnkap)
         
         dlnT_dlnkap = 0._dp
      
      endif
      
      ! Set the effective temperature. This is equal to T_int, because
      ! irradiation has no effect on the *net* flux emerging from the
      ! atmosphere
      
      ! Teff = T_int
      
      ! Finish
      
      return
   
   end subroutine eval_irradiated
   
   !****
   
   ! Evaluate atmosphere data
   
   subroutine eval_data(&
      T_int, g, L, T_eq, P_surf, kap, kap_v, gamma, skip_partials, &
      tau, lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
      ierr)
      
      use atm_utils, only : eval_E2
      
      real(dp), intent(in) :: T_int
      real(dp), intent(in) :: g
      real(dp), intent(in) :: L
      real(dp), intent(in) :: T_eq
      real(dp), intent(in) :: P_surf
      real(dp), intent(in) :: kap
      real(dp), intent(in) :: kap_v
      real(dp), intent(in) :: gamma
      logical, intent(in) :: skip_partials
      real(dp), intent(out) :: tau
      real(dp), intent(out) :: lnT
      real(dp), intent(out) :: dlnT_dL
      real(dp), intent(out) :: dlnT_dlnR
      real(dp), intent(out) :: dlnT_dlnM
      real(dp), intent(out) :: dlnT_dlnkap
      integer, intent(out) :: ierr
      
      real(dp) :: gamma_eff
      real(dp) :: x
      real(dp) :: E2
      real(dp) :: dE2_dx
      real(dp) :: f_1
      real(dp) :: f_2
      real(dp) :: T4_int
      real(dp) :: T4_eq
      real(dp) :: T4
      real(dp) :: df_1_dtau
      real(dp) :: df_2_dtau
      real(dp) :: df_1_dx
      real(dp) :: df_2_dx
      real(dp) :: dlnT_dlnT_int
      real(dp) :: dlnT_dlntau
      real(dp) :: dlnT_dlnx
      real(dp) :: dlnT_int_dlnL
      real(dp) :: dlnT_int_dlnR
      real(dp) :: dlnT_int_dlnM
      real(dp) :: dlnT_int_dlnkap
      real(dp) :: dlntau_dlnL
      real(dp) :: dlntau_dlnR
      real(dp) :: dlntau_dlnM
      real(dp) :: dlntau_dlnkap
      real(dp) :: dlnx_dlnL
      real(dp) :: dlnx_dlnR
      real(dp) :: dlnx_dlnM
      real(dp) :: dlnx_dlnkap
      
      ! Calculate the optical depth corresponding to P_surf
      
      tau = P_surf * kap / g
      
      ! Evaluate irradiation terms [cf. eq. 6 of Guillot & Havel (2011,
      ! A&A 527, A20)]
      
      if (gamma > 0._dp) then
         gamma_eff = gamma
      else
         gamma_eff = kap_v / kap
      endif
      
      x = gamma_eff * tau
      
      call eval_E2(x, E2, dE2_dx, ierr)
      if (ierr /= 0) return
      
      f_1 = 2._dp * (1._dp + (x / 2._dp - 1._dp) * exp(-x)) * tau / (3._dp * x)
      f_2 = 2._dp * x * (1 - tau * tau / 2._dp) * E2 / (3._dp * tau)
      
      ! Evaluate the temperature
      
      T4_int = T_int * T_int * T_int * T_int
      T4_eq = T_eq * T_eq * T_eq * T_eq
      
      T4 = 0.75_dp * (T4_int * (tau + 2._dp / 3._dp) + T4_eq * (f_1 + f_2 + 2._dp / 3._dp))
      
      lnT = 0.25_dp * log(T4)
      
      ! Set up partials
      
      if (.NOT. skip_partials) then
         
         df_1_dtau = (2._dp + (x - 2._dp) * exp(-x)) / (3._dp * x)
         df_2_dtau = -x * (2._dp + tau * tau) * E2 / (3._dp * tau * tau)
         
         df_1_dx = -(2._dp - (2._dp + (2._dp - x) * x) * exp(-x)) * tau / (3._dp * x * x)
         df_2_dx = -(tau * tau - 2._dp) * (E2 + x * dE2_dx) / (3._dp * tau)
         
         dlnT_dlnT_int = 0.75_dp * (tau + 2._dp / 3._dp) * T4_int / T4
         dlnT_dlntau = 0.75_dp * (T4_int + T4_eq * (df_1_dtau + df_2_dtau)) * tau / (4._dp * T4)
         dlnT_dlnx = 0.75_dp * T4_eq * (df_1_dx + df_2_dx) * x / (4._dp * T4)
         
         dlnT_int_dlnL = 0.25_dp
         dlnT_int_dlnR = -0.5_dp
         dlnT_int_dlnM = 0._dp
         dlnT_int_dlnkap = 0._dp
         
         dlntau_dlnL = 0._dp
         dlntau_dlnR = 2._dp
         dlntau_dlnM = -1._dp
         dlntau_dlnkap = 1._dp
         
         dlnx_dlnL = 0._dp
         dlnx_dlnR = 2._dp
         dlnx_dlnM = -1._dp
         if (gamma > 0._dp) then
            dlnx_dlnkap = 1._dp
         else
            dlnx_dlnkap = 0._dp
         endif
         
         dlnT_dL = (dlnT_dlnT_int * dlnT_int_dlnL + &
            dlnT_dlntau * dlntau_dlnL + &
            dlnT_dlnx * dlnx_dlnL) / L
         dlnT_dlnR = dlnT_dlnT_int * dlnT_int_dlnR + &
            dlnT_dlntau * dlntau_dlnR + &
            dlnT_dlnx * dlnx_dlnR
         dlnT_dlnM = dlnT_dlnT_int * dlnT_int_dlnM + &
            dlnT_dlntau * dlntau_dlnM + &
            dlnT_dlnx * dlnx_dlnM
         dlnT_dlnkap = dlnT_dlnT_int * dlnT_int_dlnkap + &
            dlnT_dlntau * dlntau_dlnkap + &
            dlnT_dlnx * dlnx_dlnkap
      
      else
         
         dlnT_dL = 0._dp
         dlnT_dlnR = 0._dp
         dlnT_dlnM = 0._dp
         dlnT_dlnkap = 0._dp
      
      endif
      
      ! Finish
      
      return
   
   end subroutine eval_data

end module atm_irradiated
