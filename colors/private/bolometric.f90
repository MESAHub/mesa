! ***********************************************************************
!
!   Copyright (C) 2025  Niall Miller & The MESA Team
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

module bolometric

   use const_def, only: dp
   use colors_utils, only: romberg_integration, load_lookup_table
   use hermite_interp, only: construct_sed_hermite

   implicit none

   private
   public :: calculate_bolometric

contains

   !****************************
   ! Calculate Bolometric Photometry Using Multiple SEDs
   !****************************
   subroutine calculate_bolometric(teff, log_g, metallicity, R, d, bolometric_magnitude, &
                                   bolometric_flux, wavelengths, fluxes, sed_filepath)
      real(dp), intent(in) :: teff, log_g, metallicity, R, d
      character(len=*), intent(in) :: sed_filepath
      real(dp), intent(out) :: bolometric_magnitude, bolometric_flux
      real(dp), dimension(:), allocatable, intent(out) :: wavelengths, fluxes

      real(dp), allocatable :: lu_logg(:), lu_meta(:), lu_teff(:)
      character(len=100), allocatable :: file_names(:)
      REAL, dimension(:, :), allocatable :: lookup_table
      character(len=256) :: lookup_file

      lookup_file = trim(sed_filepath)//'/lookup_table.csv'

      call load_lookup_table(lookup_file, lookup_table, file_names, lu_logg, lu_meta, lu_teff)

      call construct_sed_hermite(teff, log_g, metallicity, R, d, file_names, &
                                 lu_teff, lu_logg, lu_meta, sed_filepath, &
                                 wavelengths, fluxes)

      ! Calculate bolometric flux and magnitude
      call calculate_bolometric_phot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
   end subroutine calculate_bolometric

   !****************************
   ! Calculate Bolometric Magnitude and Flux
   !****************************
   subroutine calculate_bolometric_phot(wavelengths, fluxes, bolometric_magnitude, bolometric_flux)
      real(dp), dimension(:), intent(inout) :: wavelengths, fluxes
      real(dp), intent(out) :: bolometric_magnitude, bolometric_flux
      integer :: i

      ! Validate inputs and replace invalid wavelengths with 0
      do i = 1, size(wavelengths) - 1
         if (wavelengths(i) <= 0.0d0 .or. fluxes(i) < 0.0d0) then
            fluxes(i) = 0.0d0  ! Replace invalid wavelength with 0
         end if
      end do

      ! Call Romberg integration
      call romberg_integration(wavelengths, fluxes, bolometric_flux)

      ! Validate integration result
      if (bolometric_flux <= 0.0d0) then
         print *, "Error: Flux integration resulted in non-positive value."
         bolometric_magnitude = 99.0d0
         return
      end if

      ! Calculate bolometric magnitude
      if (bolometric_flux <= 0.0d0) then
         print *, "Error: Flux integration resulted in non-positive value."
         bolometric_magnitude = 99.0d0
         return
      else if (bolometric_flux < 1.0d-10) then
         print *, "Warning: Flux value is very small, precision might be affected."
      end if

      bolometric_magnitude = flux_to_magnitude(bolometric_flux)
   end subroutine calculate_bolometric_phot

   !****************************
   ! Convert Flux to Magnitude
   !****************************
   real(dp) function flux_to_magnitude(flux)
      real(dp), intent(in) :: flux
      if (flux <= 0.0d0) then
         print *, "Error: Flux must be positive to calculate magnitude."
         flux_to_magnitude = 99.0d0  ! Return an error value
      else
         flux_to_magnitude = -2.5d0*log10(flux)
      end if
   end function flux_to_magnitude

end module bolometric
