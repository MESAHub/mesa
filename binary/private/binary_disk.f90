! ***********************************************************************
!
!   Copyright (C) 2025  Philip Mocz, Mathieu Renzo & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module binary_disk

   use const_def, only: dp, pi, pi2, two_thirds, standard_cgrav, Msun, Rsun, secyer, crad, boltzm, clight, mp
   use star_lib
   use star_def
   use math_lib
   use binary_def

   implicit none

   ! This module includes functions for mass transfer and L2 mass loss for binary systems with a thin disk
   ! TODO: switch Kramers opacity module to real one in eval_L2_mass_loss_fraction()

contains

   subroutine eval_L2_mass_loss_fraction(donor_mass, accretor_mass, mass_transfer_rate, orbital_separation, &
                                         disk_alpha, disk_mu, &
                                         fL2, ierr)
      ! Calculate the (outer) L2 mass-loss fraction
      ! according to Lu et al. (2023, MNRAS 519, 1409) "On rapid binary mass transfer -I. Physical model"
      real(dp), intent(in) :: donor_mass         ! [M_sun]
      real(dp), intent(in) :: accretor_mass      ! [M_sun]
      real(dp), intent(in) :: mass_transfer_rate ! [M_sun/yr]
      real(dp), intent(in) :: orbital_separation ! [R_sun]
      real(dp), intent(in) :: disk_alpha         ! disk alpha viscosity parameter (dimensionless)
      real(dp), intent(in) :: disk_mu            ! disk mean molecular weight (dimensionless)
      real(dp), intent(out) :: fL2               ! L2 mass-loss fraction (dimensionless)
      integer, intent(out) :: ierr

      real(dp), parameter :: eps_small = 1d-12  ! a very small number
      real(dp), parameter :: tol = 1d-8  ! fractional tolerance for bisection method
      real(dp) :: M2, M1dot, a, q, log_q
      real(dp) :: xL1, xL2, mu
      real(dp) :: Rd_over_a, PhiL1_dimless, PhiL2_dimless, PhiRd_dimless
      real(dp) :: GM2, Rd, Phi_units, PhiL1, PhiL2, PhiRd, omega_K
      real(dp) :: c1, c2, c3, c4

      integer :: i
      integer, parameter :: n_the = 50  ! number of grid points for disk thickness search
      real(dp), parameter :: the_grid_min = 0.1_dp
      real(dp), parameter :: the_grid_max = 1.0_dp
      real(dp) :: the_grid(n_the)  ! grid for disk thickness
      real(dp) :: T_arr(n_the)
      real(dp), parameter :: T_floor = 3.0d3  ! [K] the minimum value for disk temperature solution
      real(dp) :: T, T_max, the, the_left, the_right, the_min, the_max, separation_factor, dlogthe
      real(dp) :: T_left, T_right, f1, f1_left, f2, f2_left, f2_right, f, f_left
      real(dp) :: logT_arr(n_the), logthe_grid(n_the)

      ierr = 0

      ! key parameters
      M2 = accretor_mass * Msun
      M1dot = mass_transfer_rate * Msun/secyer
      a = orbital_separation * Rsun
      q = accretor_mass / donor_mass  ! mass ratio M2/M1

      log_q = log10(q)

      ! positions of Lagrangian points (based on analytic fits)
      xL1 = -0.0355_dp * log_q**2 + 0.251_dp * abs(log_q) + 0.500_dp  ! [a = SMA]
      xL2 = 0.0756_dp * log_q**2 - 0.424_dp * abs(log_q) + 1.699_dp  ! [a]

      if (log_q > 0.0_dp) then  ! m2 is more massive
         xL1 = 1.0_dp - xL1
         xL2 = 1.0_dp - xL2
      end if
      mu = q / (1.0_dp + q)

      ! outer disk radius
      Rd_over_a = pow4(1.0_dp - xL1) / mu
      ! relavent potential energies
      PhiL1_dimless = -((1.0_dp - mu)/abs(xL1) + mu/abs(1.0_dp - xL1) + 0.5*(xL1 - mu)**2)  ! [G(M1+M2)/a]
      PhiL2_dimless = -((1.0_dp - mu)/abs(xL2) + mu/abs(1.0_dp - xL2) + 0.5*(xL2 - mu)**2)  ! [G(M1+M2)/a]
      PhiRd_dimless = -(1.0_dp - mu + mu/Rd_over_a + 0.5*(1.0_dp - mu)**2)

      GM2 = standard_cgrav * M2

      Rd = Rd_over_a*a
      Phi_units = standard_cgrav * (M2 / mu) / a
      PhiL1 = PhiL1_dimless * Phi_units
      PhiL2 = PhiL2_dimless * Phi_units
      PhiRd = PhiRd_dimless * Phi_units
      ! Keplerian frequency at Rd
      omega_K = sqrt(GM2 / Rd**3)


      ! constants involved in numerical solutions
      c1 = 2.0_dp * pi * crad * disk_alpha * Rd / (3.0_dp * omega_K * M1dot)
      c2 = boltzm * Rd / (GM2 * disk_mu * mp)
      c3 = 8.0_dp * pi**2 * crad * disk_alpha * clight * Rd**2 / (M1dot**2 * omega_K)
      c4 = 2.0_dp * pi * disk_mu * crad * disk_alpha * omega_K * mp * Rd**3 / (boltzm * M1dot)

      ! Create logarithmically spaced grid for grid search for disk thickness [only used at the beginning]
      do i = 1, n_the
         the_grid(i) = 10.0_dp**(log10(the_grid_min) + (i - 1) * (log10(the_grid_max) - log10(the_grid_min)) / (n_the - 1))
      end do

      ! only T < T_max is possible to calculate
      T_max = pow(4.0_dp / (27.0_dp * c1**2 * c2), 1.0_dp / 9.0_dp)

      do i = 1, n_the
         T_arr(i) = 0.0_dp
      end do

      do i = 1, n_the
         the = the_grid(i)
         ! use bisection method
         T_left = 0.1_dp * min(the**2 / c2, T_max)
         f1_left = f1_the_T_fL2(the, T_left, 0.0_dp, c1, c2)
         T_right = T_max
         do while (abs((T_left - T_right) / T_right) > tol)
            T = (T_left + T_right) / 2.0_dp
            f1 = f1_the_T_fL2(the, T, 0.0_dp, c1, c2)
            if (f1 * f1_left > 0) then
               T_left = T
               f1_left = f1
            else
               T_right = T
            end if
         end do
         T_arr(i) = (T_left + T_right) / 2.0_dp
      end do
      ! now we have obtained numerical relation between the and T
      do i = 1, n_the
         logT_arr(i) = log10(T_arr(i))
         logthe_grid(i) = log10(the_grid(i))
      end do
      dlogthe = logthe_grid(2) - logthe_grid(1)

      ! bisection to find the numerical solution to f2(the, T, fL2=0)=0
      the_right = 1.0_dp
      f2_right = f2_the_T_fL2(the_right, T_the_nofL2(the_right, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
                              c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
      separation_factor = 0.95_dp
      the_left = separation_factor * the_right
      f2_left = f2_the_T_fL2(the_left, T_the_nofL2(the_left, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
                             c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
      do while (f2_left * f2_right > 0)
         ! need to decrease the_left
         the_right = the_left
         f2_right = f2_left
         the_left = the_left * separation_factor
         f2_left = f2_the_T_fL2(the_left, T_the_nofL2(the_left, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
                                c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
      end do
      ! now the solution is between the_left and the_right
      do while (abs((the_left - the_right) / the_right) > tol)
         the = (the_left + the_right) / 2.0_dp
         f2 = f2_the_T_fL2(the, T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
                           c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
         if (f2 * f2_left > 0.0_dp) then
            the_left = the
            f2_left = f2
         else
            the_right = the
         end if
       end do
      ! solution
      the = (the_left + the_right) / 2.0_dp
      T = T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe)
      the_max = sqrt(3.0_dp / 8.0_dp * c2 * T + 1.0_dp / 4.0_dp * (PhiL2 - PhiRd) / (GM2 / Rd) - 1.0_dp / 8.0_dp)

      if (the < the_max) then
         ! return a tiny numner
         fL2 = eps_small
      else
         the_min = ( 1.0_dp / 2.0_dp * sqrt((PhiL2 - PhiRd) / (GM2 / Rd) - 1.0_dp / 2.0_dp) )  ! corresponding to fL2=1, T=0
         ! need to find the maximum corresponding to fL2=0
         ! this is given by the intersection between T_the(the), T_the_nofL2(the)
         the_left = the_min
         the_right = 1.0_dp
         f_left = T_the(the_left, c2, PhiL2, PhiRd, GM2, Rd) - T_the_nofL2(the_left, logthe_grid, logT_arr, n_the, dlogthe)
         do while (abs((the_left - the_right) / the_right) > tol)
            the = (the_left + the_right) / 2.0_dp
            f = T_the(the, c2, PhiL2, PhiRd, GM2, Rd) - T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe)
            if (f * f_left > 0.0_dp) then
               the_left = the
               f_left = f
            else
               the_right = the
            end if
         end do
         the_max = (the_left + the_right) / 2.0_dp  ! this corresponds to fL2=0

         ! --- numerical solution for f2(the, T, fL2)=0 under non-zero fL2

         ! -- do not use exactly the_min (corresponding to T = 0, bc. kap table breaks down)
         ! -- define another the_min based on T_floor (kap table won't be a problem)
         the_min = sqrt( 3.0_dp / 8.0_dp * c2 * T_floor &
                       + 1.0_dp / 4.0_dp * (PhiL2 - PhiRd) / (GM2 / Rd) - 1.0_dp / 8.0_dp )
         the_left = the_min
         f2_left = f2_the_T_fL2(the_left, &
                                T_the(the_left, c2, PhiL2, PhiRd, GM2, Rd), &
                                fL2_the(the_left, c1, c2, PhiL2, PhiRd, GM2, Rd), &
                                c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
         the_right = the_max / (1.0_dp + eps_small)
         ! bisection again
         do while (abs((the_left - the_right) / the_right) > tol)
            the = (the_left + the_right) / 2
            f2 = f2_the_T_fL2(the, T_the(the, c2, PhiL2, PhiRd, GM2, Rd), fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd), &
                              c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
            if (f2 * f2_left > 0.0_dp) then
               the_left = the
               f2_left = f2
            else
               the_right = the
            end if
         end do
         ! solution
         the = (the_left + the_right) / 2.0_dp
         fL2 = fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd)
      end if

   end subroutine eval_L2_mass_loss_fraction


   ! Helper Functions
   real(dp) function f1_the_T_fL2(the, T, fL2, c1, c2)
      real(dp), intent(in) :: the, T, fL2, c1, c2
      f1_the_T_fL2 = c1 * T**4 * the**3 / (1.0_dp - fL2) - the**2 + c2 * T
   end function f1_the_T_fL2

   real(dp) function  kap(rho, T)
      ! simplified Kramers rule (cgs; approximate)
      real(dp), intent(in) :: rho, T
      kap = 0.34 + 3.0e24 * rho * pow(T, -3.5_dp)
   end function kap

   real(dp) function f2_the_T_fL2(the, T, fL2, &
                                  c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
      real(dp), intent(in) :: the, T, fL2
      real(dp), intent(in) :: c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2
      real(dp) :: x, U_over_P, rho
      x = c4 * (T * the) ** 3 / (1.0_dp - fL2)
      U_over_P = (1.5_dp + x) / (1.0_dp + 1.0_dp / 3.0_dp * x)
      rho = (1.0_dp - fL2) * M1dot / (2.0_dp * pi * disk_alpha * omega_K * Rd**3) / the**3
      f2_the_T_fL2 = &
         7.0_dp / 4.0_dp &
         - (1.5_dp * U_over_P + c3 * T**4 / kap(rho, T) / (1.0_dp - fL2) ** 2) * the**2 &
         - PhiRd / (GM2 / Rd) &
         + (PhiL1 - fL2 * PhiL2) / (GM2 / Rd) / (1.0_dp - fL2)
   end function f2_the_T_fL2

   real(dp) function T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe)
      real(dp), intent(in) :: the
      integer, intent(in) :: n_the
      real(dp), intent(in) :: logthe_grid(n_the), logT_arr(n_the)
      real(dp), intent(in) :: dlogthe
      real(dp) :: logthe, slope, logT
      integer :: i_the
      ! under the assumption fL2=0
      logthe = log10(the)
      if (logthe > logthe_grid(n_the-1)) then
         ! use analytic extrapolation
         T_the_nofL2 = 10 ** (logT_arr(n_the-1) - 0.25_dp * (logthe - logthe_grid(n_the-1)))
      else if (logthe < logthe_grid(1)) then
         ! analytic extrapolation
         T_the_nofL2 = 10 ** (logT_arr(1) + 2.0_dp * (logthe - logthe_grid(1)))
      else
         i_the = floor((logthe - logthe_grid(1)) / dlogthe) + 1
         slope = (logT_arr(i_the + 1) - logT_arr(i_the)) / dlogthe
         logT = logT_arr(i_the) + (logthe - logthe_grid(i_the)) * slope
         T_the_nofL2 = 10 ** logT
      end if
   end function T_the_nofL2

   real(dp) function T_the(the, c2, PhiL2, PhiRd, GM2, Rd)
      real(dp), intent(in) :: the
      real(dp), intent(in) :: c2, PhiL2, PhiRd, GM2, Rd
      ! only for non-zero fL2
      T_the = (8.0_dp * the**2 + 1.0_dp - 2.0_dp * (PhiL2 - PhiRd) / (GM2 / Rd)) / (3.0_dp * c2)
   end function T_the

   real(dp) function fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd)
      real(dp), intent(in) :: the
      real(dp), intent(in) :: c1, c2, PhiL2, PhiRd, GM2, Rd
      real(dp) :: T
      ! only for non-zero fL2
      T = T_the(the, c2, PhiL2, PhiRd, GM2, Rd)
      fL2_the =  1.0_dp - c1 * T**4 * the**3 / (the**2 - c2 * T)
   end function fL2_the


end module binary_disk
