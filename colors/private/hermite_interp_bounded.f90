! ***********************************************************************
!
!   Copyright (C) 2026  Niall Miller & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
! ***********************************************************************

module hermite_interp_bounded
   use const_def, only: dp
   use colors_def, only: Colors_General_Info
   use colors_utils, only: simpson_integration
   use hermite_interp, only: construct_sed_hermite
   use linear_interp, only: construct_sed_linear
   use utils_lib, only: is_inf, is_nan, mesa_error

   implicit none

   private
   public :: construct_sed_hermite_bounded

   ! Accept negligible Hermite undershoot when both its integrated magnitude
   ! and its peak amplitude are below 0.1% of the positive SED.
   real(dp), parameter :: negative_flux_fraction_tol = 1.0d-3

contains

   subroutine construct_sed_hermite_bounded(rq, teff, log_g, metallicity, R, d, &
                                            stellar_model_dir, wavelengths, fluxes)
      type(Colors_General_Info), intent(inout) :: rq
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      character(len=*), intent(in) :: stellar_model_dir
      real(dp), allocatable, intent(out) :: wavelengths(:), fluxes(:)

      real(dp), allocatable :: interp_wave(:), interp_flux(:)
      real(dp) :: integrated_flux
      logical :: interp_ok

      call construct_sed_hermite(rq, teff, log_g, metallicity, R, d, &
                                 stellar_model_dir, interp_wave, interp_flux)

      ! Accept Hermite whenever it is finite and physically usable.  Negligible
      ! negative undershoots are projected onto the physical flux bound.
      interp_ok = valid_sed(interp_wave, interp_flux, .false.)
      if (interp_ok) call repair_negative_flux(interp_wave, interp_flux, interp_ok)

      if (interp_ok) then
         call simpson_integration(interp_wave, interp_flux, integrated_flux)
         interp_ok = .not. is_nan(integrated_flux) .and. &
                     .not. is_inf(integrated_flux) .and. integrated_flux > 0.0_dp
      end if

      if (interp_ok) then
         call move_alloc(interp_wave, wavelengths)
         call move_alloc(interp_flux, fluxes)
         return
      end if

      ! Hermite was unusable.  Fall back to the positivity-preserving linear
      ! interpolant so a Colors failure does not terminate the MESA run.
      call construct_sed_linear(rq, teff, log_g, metallicity, R, d, &
                                stellar_model_dir, interp_wave, interp_flux)

      interp_ok = valid_sed(interp_wave, interp_flux, .true.)
      if (interp_ok) then
         call simpson_integration(interp_wave, interp_flux, integrated_flux)
         interp_ok = .not. is_nan(integrated_flux) .and. &
                     .not. is_inf(integrated_flux) .and. integrated_flux > 0.0_dp
      end if

      if (.not. interp_ok) then
         write (*, *) 'colors: Hermite and linear interpolation returned invalid SEDs'
         call mesa_error(__FILE__, __LINE__)
      end if

      call move_alloc(interp_wave, wavelengths)
      call move_alloc(interp_flux, fluxes)
   end subroutine construct_sed_hermite_bounded

   subroutine repair_negative_flux(wavelengths, fluxes, repair_ok)
      real(dp), intent(in) :: wavelengths(:)
      real(dp), intent(inout) :: fluxes(:)
      logical, intent(out) :: repair_ok

      real(dp), allocatable :: flux_part(:)
      real(dp) :: negative_integral, positive_integral
      real(dp) :: negative_peak, positive_peak
      integer :: n

      repair_ok = .false.
      n = size(fluxes)

      if (n < 2 .or. size(wavelengths) /= n) return

      ! Nothing needs repairing.  Zero flux is already physically safe.
      if (.not. any(fluxes < 0.0_dp)) then
         repair_ok = any(fluxes > 0.0_dp)
         return
      end if

      ! A useful SED must retain a positive interior.
      if (.not. any(fluxes > 0.0_dp)) return

      allocate (flux_part(n))

      flux_part = max(0.0_dp, -fluxes)
      call simpson_integration(wavelengths, flux_part, negative_integral)
      negative_peak = maxval(flux_part)

      flux_part = max(0.0_dp, fluxes)
      call simpson_integration(wavelengths, flux_part, positive_integral)
      positive_peak = maxval(flux_part)

      deallocate (flux_part)

      if (is_nan(negative_integral) .or. is_inf(negative_integral)) return
      if (is_nan(positive_integral) .or. is_inf(positive_integral)) return
      if (is_nan(negative_peak) .or. is_inf(negative_peak)) return
      if (is_nan(positive_peak) .or. is_inf(positive_peak)) return
      if (positive_integral <= 0.0_dp .or. positive_peak <= 0.0_dp) return

      if (negative_integral/positive_integral > negative_flux_fraction_tol) return
      if (negative_peak/positive_peak > negative_flux_fraction_tol) return

      ! The undershoot is negligible; project only the unphysical samples onto
      ! the physical constraint f_lambda >= 0.
      where (fluxes < 0.0_dp) fluxes = 0.0_dp

      repair_ok = all(fluxes >= 0.0_dp) .and. any(fluxes > 0.0_dp)
   end subroutine repair_negative_flux

   logical function valid_sed(wavelengths, fluxes, require_positive_flux)
      real(dp), intent(in) :: wavelengths(:), fluxes(:)
      logical, intent(in) :: require_positive_flux

      valid_sed = size(wavelengths) >= 2 .and. size(wavelengths) == size(fluxes)
      if (.not. valid_sed) return

      valid_sed = .not. any(is_nan(wavelengths)) .and. .not. any(is_inf(wavelengths)) .and. &
                  .not. any(is_nan(fluxes)) .and. .not. any(is_inf(fluxes))
      if (.not. valid_sed) return

      valid_sed = all(wavelengths > 0.0_dp) .and. &
                  all(wavelengths(2:) > wavelengths(:size(wavelengths) - 1))
      if (.not. valid_sed) return

      if (require_positive_flux) valid_sed = all(fluxes > 0.0_dp)
   end function valid_sed

end module hermite_interp_bounded
