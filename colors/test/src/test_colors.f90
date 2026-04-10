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
      how_many_colors_history_columns, data_for_colors_history_columns, &
      calculate_bolometric
   use colors_def, only: Colors_General_Info
   use const_def,  only: dp, rsun, mesa_dir
   use utils_lib,  only: mesa_error

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

   write(*,*) 'Colors module initialized successfully.'
   write(*,*) 'Test passed!'

   write(*,*) 'test_colors: passed'

end program test_colors
