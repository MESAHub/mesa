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

module atm_T_tau_relations
   
   ! Uses
   
   use const_def
   use math_lib
   use utils_lib, only : mesa_error
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: get_T_tau_base
   public :: eval_T_tau
   public :: eval_T_tau_dq_dtau
   
   ! Procedures

contains
   
   subroutine get_T_tau_base (id, tau_base, ierr)
      
      use atm_def, only : &
         ATM_T_TAU_EDDINGTON, &
         ATM_T_TAU_SOLAR_HOPF, &
         ATM_T_TAU_KRISHNA_SWAMY, &
         ATM_T_TAU_TRAMPEDACH_SOLAR
      
      integer, intent(in) :: id
      real(dp), intent(out) :: tau_base
      integer, intent(out) :: ierr
      
      ierr = 0
      
      ! Get the base optical depth
      
      select case (id)
      case (ATM_T_TAU_EDDINGTON)
         tau_base = 2._dp / 3._dp
      case (ATM_T_TAU_SOLAR_HOPF)
         tau_base = 0.4116433502_dp
      case (ATM_T_TAU_KRISHNA_SWAMY)
         tau_base = 0.3121563_dp
      case (ATM_T_TAU_TRAMPEDACH_SOLAR)
         tau_base = 0.5147929058057147_dp
      case default
         write(*, *) 'Invalid id in get_T_tau_base: ', id
         call mesa_error(__FILE__, __LINE__)
      end select
      
      ! Finish
      
      return
   
   end subroutine get_T_tau_base
   
   !****
   
   subroutine eval_T_tau (id, tau, Teff, lnT, ierr)
      
      use atm_def, only : &
         ATM_T_TAU_EDDINGTON, &
         ATM_T_TAU_SOLAR_HOPF, &
         ATM_T_TAU_KRISHNA_SWAMY, &
         ATM_T_TAU_TRAMPEDACH_SOLAR
      
      integer, intent(in) :: id
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      integer, intent(out) :: ierr
      
      ierr = 0
      
      ! Evaluate the T-tau relation
      
      select case (id)
      case (ATM_T_TAU_EDDINGTON)
         call eval_Eddington(tau, Teff, lnT)
      case (ATM_T_TAU_SOLAR_HOPF)
         call eval_solar_Hopf(tau, Teff, lnT)
      case (ATM_T_TAU_KRISHNA_SWAMY)
         call eval_Krishna_Swamy(tau, Teff, lnT)
      case (ATM_T_TAU_TRAMPEDACH_SOLAR)
         call eval_Trampedach_solar(tau, Teff, lnT)
      case default
         write(*, *) 'Invalid id in eval_T_tau: ', id
         call mesa_error(__FILE__, __LINE__)
      end select
      
      ! Finish
      
      return
   
   end subroutine eval_T_tau
   
   !****
   
   subroutine eval_Eddington (tau, Teff, lnT)
      
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      
      real(dp) :: Teff4
      real(dp) :: T4
      
      ! Evaluate the Eddington T-tau relation
      
      Teff4 = Teff * Teff * Teff * Teff
      T4 = 0.75d0 * Teff4 * (tau + two_thirds)
      
      lnT = log(T4) * 0.25d0
      
      ! Finish
      
      return
   
   end subroutine eval_Eddington
   
   !****
   
   subroutine eval_solar_Hopf (tau, Teff, lnT)
      
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      
      real(dp), parameter :: Q1 = 1.0361_dp
      real(dp), parameter :: Q2 = -0.3134_dp
      real(dp), parameter :: Q3 = 2.44799995_dp
      real(dp), parameter :: Q4 = -0.29589999_dp
      real(dp), parameter :: Q5 = 30._dp
      
      real(dp) :: e1
      real(dp) :: e2
      real(dp) :: Teff4
      real(dp) :: T4
      
      ! Evaluate the T-tau relation for an approximate Hopf function
      ! tuned to solar data; see MESA II, Sec. A.5. This is essentially
      ! equivalent to the fit given by Sonoi et al. (2019, A&A, 621,
      ! 84)
      
      e1 = exp(-Q3 * tau)
      e2 = exp(-Q5 * tau)
      
      Teff4 = Teff * Teff * Teff * Teff
      T4 = 0.75d0 * Teff4 * (tau + Q1 + Q2 * e1 + Q4 * e2)
      
      lnT = log(T4) * 0.25d0
      
      ! Finish
      
      return
   
   end subroutine eval_solar_Hopf
   
   !****
   
   subroutine eval_Krishna_Swamy (tau, Teff, lnT)
      
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      
      real(dp), parameter :: Q1 = 1.39_dp
      real(dp), parameter :: Q2 = -0.815_dp
      real(dp), parameter :: Q3 = 2.54_dp
      real(dp), parameter :: Q4 = -0.025_dp
      real(dp), parameter :: Q5 = 30._dp
      
      real(dp) :: e1
      real(dp) :: e2
      real(dp) :: Teff4
      real(dp) :: T4
      
      ! Evaluate the T-tau relation from Krishna-Swamy (1966, ApJ 145,
      ! 174–194)
      
      e1 = exp(-Q3 * tau)
      e2 = exp(-Q5 * tau)
      
      Teff4 = Teff * Teff * Teff * Teff
      T4 = 0.75d0 * Teff4 * (tau + Q1 + Q2 * e1 + Q4 * e2)
      
      lnT = log(T4) * 0.25d0
      
      ! Finish
      
      return
   
   end subroutine eval_Krishna_Swamy
   
   !****
   
   subroutine eval_Trampedach_solar (tau, Teff, lnT)
      
      real(dp), intent(in) :: tau
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      
      real(dp), parameter :: c0 = 0.6887302005929656_dp
      real(dp), parameter :: c1 = 0.0668697860833449_dp
      real(dp), parameter :: a = 0.9262126497691250_dp
      real(dp), parameter :: v = 0.7657856893402466_dp
      real(dp), parameter :: b = 0.1148742902769433_dp
      
      real(dp) :: x
      real(dp) :: Teff4
      real(dp) :: T4
      
      ! Evaluate the T-tau relation by Ball (2021, RNAAS 5, 7),
      ! which is a fit to the solar simulation by Trampedach
      ! et al. (2014, MNRAS 442, 805–820)
      
      x = log10(tau)
      
      if (x >= 0.07407427) then
         write(*, *) 'WARNING: evaluating Trampedach_solar T-tau relation beyond valid region (log10(tau) < 0.0741):', x
      end if
      
      Teff4 = Teff * Teff * Teff * Teff
      T4 = 0.75d0 * Teff4 * (tau + c0 + c1 * (x - b) + v * exp((x - a) / v))
      
      lnT = log(T4) * 0.25d0
      
      ! Finish
      
      return
   
   end subroutine eval_Trampedach_solar
   
   !****
   
   subroutine eval_T_tau_dq_dtau (id, tau, dq_dtau, ierr)
      
      use atm_def, only : &
         ATM_T_TAU_EDDINGTON, &
         ATM_T_TAU_SOLAR_HOPF, &
         ATM_T_TAU_KRISHNA_SWAMY, &
         ATM_T_TAU_TRAMPEDACH_SOLAR
      
      integer, intent(in) :: id
      real(dp), intent(in) :: tau
      real(dp), intent(out) :: dq_dtau
      integer, intent(out) :: ierr
      
      ierr = 0
      
      ! For a T(τ) relation of the standard form,
      !
      !    (T/Teff)⁴ = (3/4)[τ + q(τ)]
      !
      ! evaluate the derivative q'(τ).
      
      select case (id)
      case (ATM_T_TAU_EDDINGTON)
         call eval_Eddington_dq_dtau(tau, dq_dtau)
      case (ATM_T_TAU_SOLAR_HOPF)
         call eval_solar_Hopf_dq_dtau(tau, dq_dtau)
      case (ATM_T_TAU_KRISHNA_SWAMY)
         call eval_Krishna_Swamy_dq_dtau(tau, dq_dtau)
      case (ATM_T_TAU_TRAMPEDACH_SOLAR)
         call eval_Trampedach_solar_dq_dtau(tau, dq_dtau)
      case default
         write(*, *) 'Invalid id in eval_T_tau_dq_dtau: ', id
         call mesa_error(__FILE__, __LINE__)
      end select
      
      ! Finish
      
      return
   
   end subroutine eval_T_tau_dq_dtau
   
   !****
   
   subroutine eval_Eddington_dq_dtau (tau, dq_dtau)
      
      real(dp), intent(in) :: tau
      real(dp), intent(out) :: dq_dtau
      
      ! Evaluate the Eddington q'(τ)
      
      dq_dtau = 0.0_dp
      
      ! Finish
      
      return
   
   end subroutine eval_Eddington_dq_dtau
   
   !****
   
   subroutine eval_solar_Hopf_dq_dtau (tau, dq_dtau)
      
      real(dp), intent(in) :: tau
      real(dp), intent(out) :: dq_dtau
      
      real(dp), parameter :: Q1 = 1.0361_dp
      real(dp), parameter :: Q2 = -0.3134_dp
      real(dp), parameter :: Q3 = 2.44799995_dp
      real(dp), parameter :: Q4 = -0.29589999_dp
      real(dp), parameter :: Q5 = 30._dp
      
      real(dp) :: e1
      real(dp) :: e2
      
      ! Evaluate q'(τ) for an approximate Hopf function
      ! tuned to solar data; see MESA II, Sec. A.5. This is essentially
      ! equivalent to the fit given by Sonoi et al. (2019, A&A, 621,
      ! 84)
      
      e1 = exp(-Q3 * tau)
      e2 = exp(-Q5 * tau)
      
      dq_dtau = - Q2 * Q3 * e1 - Q4 * Q5 * e2
      
      ! Finish
      
      return
   
   end subroutine eval_solar_Hopf_dq_dtau
   
   !****
   
   subroutine eval_Krishna_Swamy_dq_dtau (tau, dq_dtau)
      
      real(dp), intent(in) :: tau
      real(dp), intent(out) :: dq_dtau
      
      real(dp), parameter :: Q1 = 1.39_dp
      real(dp), parameter :: Q2 = -0.815_dp
      real(dp), parameter :: Q3 = 2.54_dp
      real(dp), parameter :: Q4 = -0.025_dp
      real(dp), parameter :: Q5 = 30._dp
      
      real(dp) :: e1
      real(dp) :: e2
      
      ! Evaluate q'(τ) for Krishna-Swamy (1966, ApJ 145, 174–194)
      
      e1 = exp(-Q3 * tau)
      e2 = exp(-Q5 * tau)
      
      dq_dtau = - Q2 * Q3 * e1 - Q4 * Q5 * e2
      
      ! Finish
      
      return
   
   end subroutine eval_Krishna_Swamy_dq_dtau
   
   !****
   
   subroutine eval_Trampedach_solar_dq_dtau (tau, dq_dtau)
      
      real(dp), intent(in) :: tau
      real(dp), intent(out) :: dq_dtau
      
      real(dp), parameter :: c1 = 0.0668697860833449_dp
      real(dp), parameter :: a = 0.9262126497691250_dp
      real(dp), parameter :: v = 0.7657856893402466_dp
      real(dp), parameter :: b = 0.1148742902769433_dp
      real(dp), parameter :: w = 0.0514999047169869_dp
      
      real(dp) :: x
      
      ! Evaluate q'(τ) for Ball (2021, RNAAS ...)
      ! using fit to solar simulation by
      ! Trampedach et al. (2014, MNRAS 442, 805–820)
      
      x = log10(tau)
      
      dq_dtau = (c1 + exp((x - a) / v)) / (1._dp + exp((x - b) / w)) / tau * iln10
      
      ! Finish
      
      return
   
   end subroutine eval_Trampedach_solar_dq_dtau

end module atm_T_tau_relations
