! ***********************************************************************
!
!   Copyright (C) 2026  Niall Miller & The MESA Team
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

! colors post-processing for an existing MESA history file.
!
! reads Teff, log g, R and the surface composition from each row of a
! history file, maps the composition onto the atmosphere-table metallicity
! axis the same way star does for the inline colors history columns, and
! writes the colors history columns for every row to an output file.
!
! controls come from ./inlist: the usual &colors namelist plus &post_proc.

module colors_post_controls

   use const_def, only: strlen

   implicit none

   private
   public :: read_post_proc_controls

   character(len=strlen) :: history_file, output_file

   namelist /post_proc/ history_file, output_file

contains

   subroutine read_post_proc_controls(inlist, history_file_out, output_file_out, ierr)
      use utils_namelist, only: read_namelist, missing_namelist_warning
      character(len=*), intent(in) :: inlist
      character(len=*), intent(out) :: history_file_out, output_file_out
      integer, intent(out) :: ierr

      history_file = 'LOGS/history.data'
      output_file = 'LOGS/post_proc_output.data'

      call read_namelist(inlist, read_post_proc_file, 'post_proc', ierr, missing_namelist_warning)

      history_file_out = history_file
      output_file_out = output_file
   end subroutine read_post_proc_controls

   subroutine read_post_proc_file(unit, iostat, iomsg, extra_inlists, extra_inlists_mask)
      use utils_namelist, only: max_extra_inlists
      integer, intent(in) :: unit
      integer, intent(out) :: iostat
      character(len=strlen), intent(out) :: iomsg
      character(len=strlen), dimension(max_extra_inlists), intent(out) :: extra_inlists
      logical, dimension(max_extra_inlists), intent(out) :: extra_inlists_mask

      read(unit, nml=post_proc, iostat=iostat, iomsg=iomsg)
      extra_inlists(:) = ''
      extra_inlists_mask(:) = .false.
   end subroutine read_post_proc_file

end module colors_post_controls


program colors_post

   use const_def, only: dp, strlen, rsun
   use const_lib, only: const_init
   use math_lib, only: math_init, exp10
   use utils_lib, only: mesa_error, mkdir, folder_exists
   use colors_def, only: Colors_General_Info
   use colors_lib, only: colors_init, colors_shutdown, alloc_colors_handle_using_inlist, &
                         free_colors_handle, colors_ptr, colors_setup_tables, &
                         how_many_colors_history_columns, data_for_colors_history_columns
   use colors_post_controls, only: read_post_proc_controls

   implicit none

   character(len=*), parameter :: inlist = 'inlist'

   ! columns copied from the history file to the front of each output row
   integer, parameter :: num_lead_cols = 4

   character(len=strlen) :: history_file, output_file
   integer :: handle, ierr, n_cols, n_models, n_failed, i
   type(Colors_General_Info), pointer :: rq

   integer, allocatable :: model_number(:)
   real(dp), allocatable :: age(:), teff(:), log_g(:), radius(:), x_surf(:), z_surf(:), m_div_h(:)

   character(len=80), allocatable :: col_names(:)
   real(dp), allocatable :: vals(:, :)
   real(dp) :: min_m_div_h, max_m_div_h

   call const_init('', ierr)
   if (ierr /= 0) then
      write(*, '(a)') 'colors_post: const_init failed'
      call mesa_error(__FILE__, __LINE__)
   end if
   call math_init()

   call read_post_proc_controls(inlist, history_file, output_file, ierr)
   if (ierr /= 0) then
      write(*, '(2a)') 'colors_post: failed to read &post_proc from ', inlist
      call mesa_error(__FILE__, __LINE__)
   end if

   call colors_init(.true., '', ierr)
   if (ierr /= 0) then
      write(*, '(a)') 'colors_post: colors_init failed'
      call mesa_error(__FILE__, __LINE__)
   end if

   handle = alloc_colors_handle_using_inlist(inlist, ierr)
   if (ierr /= 0) then
      write(*, '(2a)') 'colors_post: failed to read &colors from ', inlist
      call mesa_error(__FILE__, __LINE__)
   end if

   call colors_ptr(handle, rq, ierr)
   if (ierr /= 0) then
      write(*, '(a)') 'colors_post: colors_ptr failed'
      call mesa_error(__FILE__, __LINE__)
   end if

   ! a run made without inline colors usually has use_colors = .false.,
   ! and alloc_colors_handle_using_inlist then skips loading the tables
   if (.not. rq% use_colors) then
      write(*, '(a)') 'colors_post: use_colors = .false. in &colors, turning it on for post-processing'
      rq% use_colors = .true.
      call colors_setup_tables(handle, ierr)
      if (ierr /= 0) then
         write(*, '(a)') 'colors_post: colors_setup_tables failed'
         call mesa_error(__FILE__, __LINE__)
      end if
   end if

   if (rq% z_over_x_ref <= 0d0) then
      write(*, '(a)') 'colors_post: z_over_x_ref must be positive'
      call mesa_error(__FILE__, __LINE__)
   end if

   call read_history_file(history_file, n_models, model_number, age, teff, log_g, radius, &
                          x_surf, z_surf, ierr)
   if (ierr /= 0) then
      write(*, '(2a)') 'colors_post: failed to read history file ', trim(history_file)
      call mesa_error(__FILE__, __LINE__)
   end if
   write(*, '(a,i0,2a)') 'colors_post: read ', n_models, ' models from ', trim(history_file)

   ! same mapping of the surface Z/X onto the atmosphere-table metallicity axis as in
   ! star/private/history.f90, so that a run and its post-processing agree on [M/H]
   min_m_div_h = minval(rq% lu_meta)
   max_m_div_h = maxval(rq% lu_meta)
   allocate(m_div_h(n_models))
   do i = 1, n_models
      if (x_surf(i) > 0d0 .and. z_surf(i) > 0d0) then
         m_div_h(i) = log10((z_surf(i)/x_surf(i))/rq% z_over_x_ref)
         m_div_h(i) = max(min_m_div_h, min(max_m_div_h, m_div_h(i)))
      else
         m_div_h(i) = min_m_div_h
      end if
   end do

   n_cols = how_many_colors_history_columns(handle)
   allocate(col_names(n_cols), vals(n_cols, n_models))
   col_names(:) = 'unknown'
   vals(:, :) = -1d99

   n_failed = 0
   do i = 1, n_models
      call data_for_colors_history_columns( &
         teff(i), log_g(i), radius(i), m_div_h(i), model_number(i), &
         handle, n_cols, col_names, vals(:, i), ierr)
      if (ierr /= 0) then
         n_failed = n_failed + 1
         write(*, '(a,i0,a,i0)') 'colors_post: colors calculation failed for model ', &
            model_number(i), ', ierr = ', ierr
      end if
   end do

   call write_output_file(output_file, ierr)
   if (ierr /= 0) then
      write(*, '(2a)') 'colors_post: failed to write ', trim(output_file)
      call mesa_error(__FILE__, __LINE__)
   end if
   write(*, '(2a)') 'colors_post: wrote ', trim(output_file)

   call free_colors_handle(handle)
   call colors_shutdown()

   if (n_failed == 0) then
      write(*, '(a,i0,a)') 'colors_post: successfully post-processed all ', n_models, ' models'
   else
      write(*, '(a,i0,a,i0,a)') 'colors_post: colors calculation failed for ', n_failed, &
         ' of ', n_models, ' models'
   end if

contains

   subroutine read_history_file(filename, n_models, model_number, age, teff, log_g, radius, &
                                x_surf, z_surf, ierr)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: n_models, ierr
      integer, allocatable, intent(out) :: model_number(:)
      real(dp), allocatable, intent(out) :: age(:), teff(:), log_g(:), radius(:), x_surf(:), z_surf(:)

      ! a MESA history file has 5 header lines, then the column names on line 6
      integer, parameter :: num_header_lines = 5

      integer :: io, ios, i, n_cols
      integer :: idx_model, idx_age, idx_teff, idx_logg, idx_radius, idx_h1, idx_he4
      logical :: teff_is_log, radius_is_log
      character(len=:), allocatable :: line
      character(len=80), allocatable :: headers(:)
      real(dp), allocatable :: row(:)

      ierr = 0
      n_models = 0

      open(newunit=io, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*, '(2a)') 'colors_post: cannot open history file ', trim(filename)
         ierr = -1
         return
      end if

      do i = 1, num_header_lines + 1
         call read_line(io, line, ios)
         if (ios /= 0) then
            write(*, '(2a)') 'colors_post: not a MESA history file (too few header lines): ', trim(filename)
            ierr = -1
            close(io)
            return
         end if
      end do

      call split_words(line, headers, n_cols)

      idx_model = find_col(headers, 'model_number')
      idx_age = find_col(headers, 'star_age')
      idx_logg = find_col(headers, 'log_g')
      idx_h1 = find_col(headers, 'surface_h1')
      idx_he4 = find_col(headers, 'surface_he4')

      idx_teff = find_col(headers, 'log_Teff')
      teff_is_log = (idx_teff > 0)
      if (.not. teff_is_log) idx_teff = find_col(headers, 'Teff')

      idx_radius = find_col(headers, 'log_R')
      radius_is_log = (idx_radius > 0)
      if (.not. radius_is_log) idx_radius = find_col(headers, 'radius')

      call require_col(idx_age, 'star_age', ierr)
      call require_col(idx_teff, 'log_Teff or Teff', ierr)
      call require_col(idx_logg, 'log_g', ierr)
      call require_col(idx_radius, 'log_R or radius', ierr)
      call require_col(idx_h1, 'surface_h1', ierr)
      call require_col(idx_he4, 'surface_he4', ierr)
      if (ierr /= 0) then
         write(*, '(a)') 'colors_post: add the missing columns to the history_columns.list of the run ' // &
            'that wrote this history file'
         close(io)
         return
      end if

      do
         call read_line(io, line, ios)
         if (ios /= 0) exit
         if (len_trim(line) > 0) n_models = n_models + 1
      end do

      allocate(model_number(n_models), age(n_models), teff(n_models), log_g(n_models), &
               radius(n_models), x_surf(n_models), z_surf(n_models), row(n_cols))

      rewind(io)
      do i = 1, num_header_lines + 1
         call read_line(io, line, ios)
      end do

      do i = 1, n_models
         read(io, *, iostat=ios) row
         if (ios /= 0) then
            write(*, '(a,i0,2a)') 'colors_post: failed to read data row ', i, ' of ', trim(filename)
            ierr = -1
            close(io)
            return
         end if

         if (idx_model > 0) then
            model_number(i) = nint(row(idx_model))
         else
            model_number(i) = i
         end if
         age(i) = row(idx_age)
         log_g(i) = row(idx_logg)

         if (teff_is_log) then
            teff(i) = exp10(row(idx_teff))
         else
            teff(i) = row(idx_teff)
         end if

         if (radius_is_log) then
            radius(i) = exp10(row(idx_radius))*rsun
         else
            radius(i) = row(idx_radius)*rsun
         end if

         x_surf(i) = row(idx_h1)
         z_surf(i) = 1d0 - row(idx_h1) - row(idx_he4)
      end do

      close(io)
   end subroutine read_history_file


   subroutine require_col(idx, name, ierr)
      integer, intent(in) :: idx
      character(len=*), intent(in) :: name
      integer, intent(inout) :: ierr
      if (idx < 1) then
         write(*, '(3a)') 'colors_post: history file has no ', name, ' column'
         ierr = -1
      end if
   end subroutine require_col


   subroutine write_output_file(filename, ierr)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ierr

      character(len=*), parameter :: lead_names(num_lead_cols) = &
         [character(len=8) :: 'star_age', 'Teff', 'log_g', 'm_div_h']
      integer :: unit, i, j, slash

      ierr = 0

      slash = index(filename, '/', back=.true.)
      if (slash > 1) then
         if (.not. folder_exists(filename(1:slash-1))) call mkdir(filename(1:slash-1))
      end if

      open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ierr)
      if (ierr /= 0) return

      do j = 1, num_lead_cols
         write(unit, '(a32)', advance='no') trim(lead_names(j))
      end do
      do j = 1, n_cols
         write(unit, '(a32)', advance='no') trim(col_names(j))
      end do
      write(unit, '(a)') ''

      do i = 1, n_models
         write(unit, '(4(1x,es31.16e3))', advance='no') age(i), teff(i), log_g(i), m_div_h(i)
         do j = 1, n_cols
            write(unit, '(1x,es31.16e3)', advance='no') vals(j, i)
         end do
         write(unit, '(a)') ''
      end do

      close(unit)
   end subroutine write_output_file


   ! read one record of any length
   subroutine read_line(unit, line, ios)
      integer, intent(in) :: unit
      character(len=:), allocatable, intent(out) :: line
      integer, intent(out) :: ios

      character(len=1024) :: chunk
      integer :: n

      line = ''
      do
         read(unit, '(a)', advance='no', iostat=ios, size=n) chunk
         line = line // chunk(1:n)
         if (ios /= 0) exit
      end do
      if (is_iostat_eor(ios)) ios = 0
   end subroutine read_line


   subroutine split_words(str, words, n)
      character(len=*), intent(in) :: str
      character(len=80), allocatable, intent(out) :: words(:)
      integer, intent(out) :: n

      integer :: i, start, pass
      logical :: in_word

      do pass = 1, 2
         n = 0
         in_word = .false.
         start = 1
         do i = 1, len(str) + 1
            if (i <= len(str)) then
               if (str(i:i) /= ' ') then
                  if (.not. in_word) start = i
                  in_word = .true.
                  cycle
               end if
            end if
            if (in_word) then
               n = n + 1
               if (pass == 2) words(n) = str(start:i-1)
            end if
            in_word = .false.
         end do
         if (pass == 1) allocate(words(n))
      end do
   end subroutine split_words


   integer function find_col(headers, name) result(idx)
      character(len=*), intent(in) :: headers(:)
      character(len=*), intent(in) :: name
      integer :: j

      idx = -1
      do j = 1, size(headers)
         if (trim(headers(j)) == name) then
            idx = j
            return
         end if
      end do
   end function find_col

end program colors_post
