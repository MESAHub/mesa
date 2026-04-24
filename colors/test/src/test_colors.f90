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

! unit test for the MESA colors module.
!
! compares synthetic magnitudes (Vega system, Kurucz2003, Johnson filters) and
! sampled SED flux values against a reference test_output file generated from
! a known-good state of the code.  run via ./ck, which invokes ./rn and
! compares stdout against test_output using diff -b.
!
! group 1 - representative stellar types:
!   solar      Teff = 5778 K,  log g = 4.44,  [M/H] = 0.0
!   hot_ms     Teff = 15000 K, log g = 4.00,  [M/H] = 0.0
!   cool_giant Teff = 4000 K,  log g = 2.00,  [M/H] = 0.0
!
! group 2 - grid sweeps (exercises each table dimension independently):
!   vary [M/H]:  Teff = 5778, log g = 4.44, [M/H] = -2.0, -1.0, 0.0, 0.5
!   vary log g:  Teff = 5778, [M/H] = 0.0,  log g = 1.0, 2.5, 4.0, 5.0
!   vary Teff:   log g = 4.0, [M/H] = 0.0,  Teff = 3500, 6000, 10000, 20000

program test_colors

   use const_lib,  only: const_init
   use math_lib,   only: math_init
   use colors_lib, only: &
      colors_init, colors_shutdown, &
      alloc_colors_handle_using_inlist, free_colors_handle, colors_ptr, &
      colors_setup_tables, colors_setup_hooks, how_many_colors_history_columns, &
      data_for_colors_history_columns, calculate_bolometric
   use colors_def, only: Colors_General_Info
   use const_def,  only: dp, rsun, boltz_sigma
   use utils_lib,  only: mesa_error
   use colors_utils, only: resolve_path

   implicit none

   ! -----------------------------------------------------------------------
   ! group 1: representative stellar types
   ! -----------------------------------------------------------------------

   integer, parameter :: n_cases = 3

   real(dp), parameter :: test_teff(n_cases) = [5778d0,   15000d0,  4000d0  ]
   real(dp), parameter :: test_logg(n_cases) = [4.44d0,   4.0d0,    2.0d0   ]
   real(dp), parameter :: test_meta(n_cases) = [0.0d0,    0.0d0,    0.0d0   ]
   real(dp), parameter :: test_R(n_cases)    = [rsun,     5d0*rsun, 20d0*rsun]

   character(len=12), parameter :: labels(n_cases) = &
      ['solar       ', 'hot_ms      ', 'cool_giant  ']

   ! -----------------------------------------------------------------------
   ! group 2a: fixed Teff=5778, fixed log g=4.44, varying [M/H]
   ! -----------------------------------------------------------------------

   integer, parameter :: n_meta = 4

   real(dp), parameter :: sweep_meta(n_meta) = [-2.0d0, -1.0d0, 0.0d0, 0.5d0]

   ! -----------------------------------------------------------------------
   ! group 2b: fixed Teff=5778, fixed [M/H]=0.0, varying log g
   ! -----------------------------------------------------------------------

   integer, parameter :: n_logg = 4

   real(dp), parameter :: sweep_logg(n_logg) = [1.0d0, 2.5d0, 4.0d0, 5.0d0]

   ! -----------------------------------------------------------------------
   ! group 2c: fixed log g=4.0, fixed [M/H]=0.0, varying Teff
   ! -----------------------------------------------------------------------

   integer, parameter :: n_teff = 4

   real(dp), parameter :: sweep_teff(n_teff) = [3500d0, 6000d0, 10000d0, 20000d0]

   ! -----------------------------------------------------------------------
   ! shared
   ! -----------------------------------------------------------------------

   ! 10 parsecs in cm -> absolute magnitudes
   real(dp), parameter :: d_10pc = 3.0857d19

   ! number of sampled SED points printed for the SED comparison
   integer, parameter :: n_sed_samples = 20

   ! require the integrated SED flux to stay close to sigma*T^4 scaled to
   ! the test distance. large deficits usually mean the wavelength coverage
   ! is missing too much UV or IR flux.
   real(dp), parameter :: bol_flux_rel_tol = 5d-3

   character(len=32) :: my_mesa_dir
   integer :: handle, ierr, n_cols, i, j, k
   integer :: model_num
   type(Colors_General_Info), pointer :: cs
   character(len=80), allocatable :: col_names(:)
   real(dp), allocatable :: col_vals(:)
   real(dp), allocatable :: wavelengths(:), fluxes(:)
   real(dp) :: bol_mag, bol_flux, interp_rad
   character(len=256) :: sed_filepath
   integer :: n_wav, stride

   ! -----------------------------------------------------------------------
   ! module initialization
   ! -----------------------------------------------------------------------

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0) then
      write(*,*) 'const_init failed'
      call mesa_error(__FILE__, __LINE__)
   end if

   call math_init()

   call colors_init(.false., '', ierr)
   if (ierr /= 0) then
      write(*,*) 'colors_init failed, ierr =', ierr
      stop 1
   end if

   ! -----------------------------------------------------------------------
   ! handle setup: empty inlist string -> defaults (Kurucz2003 + Johnson)
   ! -----------------------------------------------------------------------

   handle = alloc_colors_handle_using_inlist('', ierr)
   if (ierr /= 0) then
      write(*,*) 'alloc_colors_handle_using_inlist failed, ierr =', ierr
      stop 1
   end if

   call colors_ptr(handle, cs, ierr)
   if (ierr /= 0) then
      write(*,*) 'colors_ptr failed, ierr =', ierr
      stop 1
   end if

   ! enable photometry and explicitly load colors data, since the default
   ! namelist leaves use_colors = .false. and skips setup at handle.
   ! To do : This could be replaced by reading from an inlist instead in the future.
   cs%use_colors = .true.
   cs%mag_system = 'Vega'
   call colors_setup_tables(handle, ierr)
   if (ierr /= 0) then
      write(*,*) 'colors_setup_tables failed, ierr =', ierr
      stop 1
   end if

   ! probably unecessary since this is empty right now.
   call colors_setup_hooks(handle, ierr)
   if (ierr /= 0) then
      write(*,*) 'colors_setup_hooks failed, ierr =', ierr
      stop 1
   end if

   n_cols = how_many_colors_history_columns(handle)
   if (n_cols == 0) then
      write(*,*) 'how_many_colors_history_columns returned 0'
      stop 1
   end if

   allocate(col_names(n_cols), col_vals(n_cols))
   model_num = 0

   ! -----------------------------------------------------------------------
   ! group 1: representative stellar types
   ! -----------------------------------------------------------------------

   call write_section_header('# Group1  system=Vega  grid=Kurucz2003  filters=Johnson')
   do j = 1, n_cases
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         test_teff(j), test_logg(j), test_R(j), test_meta(j), model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group1 case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, a)') '# case: ', trim(adjustl(labels(j)))
      do k = 1, n_cols
         write(*,'(a40, 1pe23.13)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! group 2a: vary [M/H]
   ! -----------------------------------------------------------------------

   call write_section_header('# Group2a  vary_MH  Teff=5778  logg=4.44')
   do j = 1, n_meta
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         5778d0, 4.44d0, rsun, sweep_meta(j), model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group2a case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, f6.2)') '# MH= ', sweep_meta(j)
      do k = 1, n_cols
         write(*,'(a40, 1pe23.13)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! group 2b: vary log g
   ! -----------------------------------------------------------------------

   call write_section_header('# Group2b  vary_logg  Teff=5778  MH=0.0')
   do j = 1, n_logg
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         5778d0, sweep_logg(j), rsun, 0.0d0, model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group2b case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, f6.2)') '# logg= ', sweep_logg(j)
      do k = 1, n_cols
         write(*,'(a40, 1pe23.13)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! group 2c: vary Teff
   ! -----------------------------------------------------------------------

   call write_section_header('# Group2c  vary_Teff  logg=4.0  MH=0.0')
   do j = 1, n_teff
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         sweep_teff(j), 4.0d0, rsun, 0.0d0, model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group2c case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, f10.1)') '# Teff= ', sweep_teff(j)
      do k = 1, n_cols
         write(*,'(a40, 1pe23.13)') trim(col_names(k)), col_vals(k)
      end do
   end do

   sed_filepath = trim(resolve_path(cs%stellar_atm))

   ! -----------------------------------------------------------------------
   ! group 3: solar SED sample plus wavelength-coverage sanity checks
   ! -----------------------------------------------------------------------

   call write_section_header('# Group3  SED sample + wavelength_coverage_sanity')
   write(*,'(a)') '# SED sample  case=solar  Teff=5778  logg=4.44  FeH=0.0  columns=wavelength_AA  flux_erg_s_cm2_AA'
   call calculate_bolometric( &
      cs, test_teff(1), test_logg(1), test_meta(1), test_R(1), d_10pc, &
      bol_mag, bol_flux, wavelengths, fluxes, sed_filepath, interp_rad)

   n_wav = size(wavelengths)
   stride = max(1, n_wav / n_sed_samples)

   do i = 1, n_wav, stride
      write(*,'(1pe23.13, 1x, 1pe23.13)') wavelengths(i), fluxes(i)
   end do

   write(*,'(a)') ''
   write(*,'(a, 1pe10.2)') '# wavelength_coverage_sanity  logg=4.0  FeH=0.0  rel_tol=', bol_flux_rel_tol
   do j = 1, n_teff
      call check_bolometric_coverage( &
         cs, sed_filepath, 'Teff', sweep_teff(j), 4.0d0, 0.0d0, rsun, bol_flux_rel_tol)
   end do

   ! -----------------------------------------------------------------------
   ! cleanup
   ! -----------------------------------------------------------------------

   deallocate(col_names, col_vals)
   if (allocated(wavelengths)) deallocate(wavelengths)
   if (allocated(fluxes))      deallocate(fluxes)

   call free_colors_handle(handle)
   call colors_shutdown()

   write(*,*) 'test_colors: passed'

contains

   subroutine write_section_header(title)
      character(len=*), intent(in) :: title

      write(*,'(a)') ''
      write(*,'(a)') trim(title)
      write(*,'(a)') ''
   end subroutine write_section_header

   subroutine check_bolometric_coverage(rq, sed_path, label, teff, log_g, metallicity, radius, rel_tol)
      type(Colors_General_Info), intent(inout) :: rq
      character(len=*), intent(in) :: sed_path, label
      real(dp), intent(in) :: teff, log_g, metallicity, radius, rel_tol

      real(dp), allocatable :: local_wavelengths(:), local_fluxes(:)
      real(dp) :: local_mag, local_flux, local_interp_rad
      real(dp) :: expected_flux, abs_err, rel_err

      call calculate_bolometric( &
         rq, teff, log_g, metallicity, radius, d_10pc, &
         local_mag, local_flux, local_wavelengths, local_fluxes, sed_path, local_interp_rad)

      expected_flux = boltz_sigma*teff**4*(radius/d_10pc)**2
      abs_err = abs(local_flux - expected_flux)
      rel_err = abs(local_flux - expected_flux)/expected_flux

      write(*,'(a, a, a, f10.1)') '# case: ', trim(label), '=', teff
      write(*,'(a40, 1pe23.13)') 'Wav_min_AA', minval(local_wavelengths)
      write(*,'(a40, 1pe23.13)') 'Wav_max_AA', maxval(local_wavelengths)
      write(*,'(a40, 1pe23.13)') 'Flux_actual', local_flux
      write(*,'(a40, 1pe23.13)') 'Flux_expected', expected_flux
      write(*,'(a40, 1pe23.13)') 'Flux_abserr', abs_err
      write(*,'(a40, 1pe23.13)') 'Flux_relerr', rel_err

      if (rel_err > rel_tol) then
         write(*,'(a, a, a, 1pe11.3, a, 1pe11.3)') &
            'wavelength coverage sanity check failed for ', trim(label), &
            ': rel_err=', rel_err, ' > rel_tol=', rel_tol
         stop 1
      end if
   end subroutine check_bolometric_coverage
end program test_colors
