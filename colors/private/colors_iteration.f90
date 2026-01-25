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

module colors_iteration

   use const_def, only: dp, mesa_dir
   use colors_def
   use colors_utils, only: remove_dat
   use bolometric, only: calculate_bolometric
   use synthetic, only: calculate_synthetic

   implicit none
   private
   public :: write_iteration_colors, open_iteration_file, close_iteration_file

contains

   subroutine open_iteration_file(colors_handle, ierr)
      integer, intent(in) :: colors_handle
      integer, intent(out) :: ierr
      type(Colors_General_Info), pointer :: cs
      character(len=512) :: filename
      character(len=100) :: filter_name
      integer :: i

      ierr = 0
      call get_colors_ptr(colors_handle, cs, ierr)
      if (ierr /= 0) return

      if (cs%iteration_file_open) return  ! already open

      filename = trim(cs%colors_results_directory)//'/iteration_colors.data'
      open(newunit=cs%iteration_output_unit, file=trim(filename), &
           status='replace', action='write', iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'Error opening iteration colors file: ', trim(filename)
         return
      end if

      ! Write header
      write(cs%iteration_output_unit, '(a)', advance='no') &
         '#  model  iter       star_age             dt           Teff'
      write(cs%iteration_output_unit, '(a)', advance='no') &
         '          log_g              R       Mag_bol      Flux_bol'
      do i = 1, num_color_filters
         filter_name = trim(remove_dat(color_filter_names(i)))
         write(cs%iteration_output_unit, '(2x,a14)', advance='no') trim(filter_name)
      end do
      write(cs%iteration_output_unit, *)

      cs%iteration_file_open = .true.

   end subroutine open_iteration_file


   subroutine write_iteration_colors( &
         colors_handle, model_number, iter, star_age, dt, &
         t_eff, log_g, R, metallicity, ierr)

      integer, intent(in) :: colors_handle, model_number, iter
      real(dp), intent(in) :: star_age, dt, t_eff, log_g, R, metallicity
      integer, intent(out) :: ierr

      type(Colors_General_Info), pointer :: cs
      real(dp) :: bolometric_magnitude, bolometric_flux, interpolation_radius
      real(dp) :: magnitude, d, zero_point
      character(len=256) :: sed_filepath
      real(dp), dimension(:), allocatable :: wavelengths, fluxes
      integer :: i, iounit
      logical :: make_sed

      ierr = 0
      call get_colors_ptr(colors_handle, cs, ierr)
      if (ierr /= 0) return

      ! Check if per-iteration colors is enabled
      if (.not. cs%colors_per_iteration) return
      if (.not. cs%use_colors) return

      ! Verify data was loaded at initialization
      if (.not. cs%lookup_loaded) then
         write(*,*) 'colors_iteration error: lookup table not loaded'
         ierr = -1
         return
      end if
      if (.not. cs%filters_loaded) then
         write(*,*) 'colors_iteration error: filter data not loaded'
         ierr = -1
         return
      end if

      ! Open file if needed
      if (.not. cs%iteration_file_open) then
         call open_iteration_file(colors_handle, ierr)
         if (ierr /= 0) return
      end if

      iounit = cs%iteration_output_unit

      ! Get paths and settings from colors_settings
      d = cs%distance
      sed_filepath = trim(mesa_dir)//cs%stellar_atm
      make_sed = .false.  ! Don't write individual SEDs for iteration output

      write(*,*) 'hello, it is I, newton solver'

      ! Calculate bolometric magnitude using cached lookup table
      ! Must pass the cached lookup arrays for atmosphere interpolation
      call calculate_bolometric(t_eff, log_g, metallicity, R, d, &
                                bolometric_magnitude, bolometric_flux, wavelengths, fluxes, &
                                sed_filepath, interpolation_radius, &
                                cs%lu_file_names, cs%lu_teff, cs%lu_logg, cs%lu_meta)

      ! Write basic data
      write(iounit, '(i8, i6)', advance='no') model_number, iter
      write(iounit, '(6(1pe15.7))', advance='no') &
         star_age, dt, t_eff, log_g, R, bolometric_magnitude
      write(iounit, '(1pe15.7)', advance='no') bolometric_flux

      ! Calculate and write each filter magnitude using cached filter data
      do i = 1, num_color_filters
         if (t_eff >= 0 .and. metallicity >= 0 .and. &
             allocated(wavelengths) .and. allocated(fluxes)) then

            ! Select precomputed zero-point based on magnitude system
            select case (trim(cs%mag_system))
            case ('VEGA', 'Vega', 'vega')
               zero_point = cs%filters(i)%vega_zero_point
            case ('AB', 'ab')
               zero_point = cs%filters(i)%ab_zero_point
            case ('ST', 'st')
               zero_point = cs%filters(i)%st_zero_point
            case default
               zero_point = cs%filters(i)%vega_zero_point
            end select

            ! Calculate synthetic magnitude using cached filter data and precomputed zero-point
            ! Arguments: Teff, log_g, metallicity, ierr, SED wavelengths, SED fluxes,
            !            filter wavelengths, filter transmission, zero_point,
            !            filter_name, make_sed_file, output_directory
            magnitude = calculate_synthetic(t_eff, log_g, metallicity, ierr, &
                                                             wavelengths, fluxes, &
                                                             cs%filters(i)%wavelengths, &
                                                             cs%filters(i)%transmission, &
                                                             zero_point, &
                                                             color_filter_names(i), &
                                                             make_sed, cs%sed_per_model, &
                                                             cs%colors_results_directory, model_number)
               



            if (ierr /= 0) magnitude = -99.0_dp
         else
            magnitude = -99.0_dp
         end if

         write(iounit, '(1pe15.7)', advance='no') magnitude
      end do

      write(iounit, *)  ! newline

      ! Clean up allocated arrays from calculate_bolometric
      if (allocated(wavelengths)) deallocate(wavelengths)
      if (allocated(fluxes)) deallocate(fluxes)

   end subroutine write_iteration_colors


   subroutine close_iteration_file(colors_handle, ierr)
      integer, intent(in) :: colors_handle
      integer, intent(out) :: ierr
      type(Colors_General_Info), pointer :: cs

      ierr = 0
      call get_colors_ptr(colors_handle, cs, ierr)
      if (ierr /= 0) return

      if (cs%iteration_file_open) then
         close(cs%iteration_output_unit)
         cs%iteration_file_open = .false.
      end if
   end subroutine close_iteration_file

end module colors_iteration