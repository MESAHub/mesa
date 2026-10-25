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
! copies the history file row by row, adding the colors history columns
! computed from each row's Teff, log g, R and surface composition. columns
! that already exist (e.g. from a run with inline colors) are replaced.
!
! controls come from ./inlist: the usual &colors namelist plus &post_proc.

module colors_post_controls

   use const_def, only: dp, strlen

   implicit none

   private
   public :: read_post_proc_controls, unset_m_div_h

   ! fallback_m_div_h value meaning "not set"
   real(dp), parameter :: unset_m_div_h = -1d99

   character(len=strlen) :: history_file, output_file
   real(dp) :: fallback_m_div_h

   namelist /post_proc/ history_file, output_file, fallback_m_div_h

contains

   subroutine read_post_proc_controls(inlist, history_file_out, output_file_out, fallback_m_div_h_out, ierr)
      use utils_namelist, only: read_namelist, missing_namelist_error
      character(len=*), intent(in) :: inlist
      character(len=*), intent(out) :: history_file_out, output_file_out
      real(dp), intent(out) :: fallback_m_div_h_out
      integer, intent(out) :: ierr

      history_file = ''
      output_file = 'LOGS/history.data'
      fallback_m_div_h = unset_m_div_h

      call read_namelist(inlist, read_post_proc_file, 'post_proc', ierr, missing_namelist_error)
      if (ierr /= 0) return

      if (len_trim(history_file) == 0) then
         write(*, '(2a)') 'colors_post: set history_file in &post_proc of ', inlist
         ierr = -1
         return
      end if

      history_file_out = history_file
      output_file_out = output_file
      fallback_m_div_h_out = fallback_m_div_h
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
   use utils_lib, only: mesa_error, mkdir, folder_exists, mv
   use colors_def, only: Colors_General_Info
   use colors_lib, only: colors_init, colors_shutdown, alloc_colors_handle_using_inlist, &
                         free_colors_handle, colors_ptr, colors_setup_tables, &
                         how_many_colors_history_columns, data_for_colors_history_columns
   use colors_post_controls, only: read_post_proc_controls, unset_m_div_h

   implicit none

   character(len=*), parameter :: inlist = 'inlist'

   ! a MESA history file starts with 3 header lines and a blank line, then
   ! the column numbers and the column names. columns are this wide.
   integer, parameter :: num_header_lines = 4, col_width = 40

   character(len=strlen) :: history_file, output_file
   character(len=:), allocatable :: line, temp_file
   character(len=80), allocatable :: in_names(:), out_names(:), words(:), col_names(:)
   integer, allocatable :: dest(:)
   real(dp), allocatable :: vals(:)
   real(dp) :: fallback_m_div_h, min_m_div_h, max_m_div_h
   integer :: handle, ierr, ios, in_unit, out_unit, i
   integer :: n_in, n_out, n_words, n_colors, n_models, n_failed
   integer :: idx_model, idx_teff, idx_logg, idx_radius, idx_h1, idx_he4
   logical :: teff_is_log, radius_is_log, have_composition
   type(Colors_General_Info), pointer :: rq

   call const_init('', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: const_init failed')
   call math_init()

   call read_post_proc_controls(inlist, history_file, output_file, fallback_m_div_h, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: failed to read &post_proc')

   call colors_init(.true., '', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: colors_init failed')
   handle = alloc_colors_handle_using_inlist(inlist, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: failed to read &colors')
   call colors_ptr(handle, rq, ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: colors_ptr failed')

   ! a run made without inline colors usually has use_colors = .false.,
   ! and alloc_colors_handle_using_inlist then skips loading the tables
   if (.not. rq% use_colors) then
      write(*, '(a)') 'colors_post: use_colors = .false. in &colors, turning it on for post-processing'
      rq% use_colors = .true.
      call colors_setup_tables(handle, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: colors_setup_tables failed')
   end if

   if (rq% z_over_x_ref <= 0d0) call mesa_error(__FILE__, __LINE__, 'colors_post: z_over_x_ref must be positive')
   min_m_div_h = minval(rq% lu_meta)
   max_m_div_h = maxval(rq% lu_meta)

   n_colors = how_many_colors_history_columns(handle)
   allocate(col_names(n_colors), vals(n_colors), dest(n_colors))
   col_names(:) = 'unknown'

   open(newunit=in_unit, file=trim(history_file), status='old', action='read', iostat=ios)
   if (ios /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: cannot open '//trim(history_file))

   do i = 1, num_header_lines + 2
      call read_line(in_unit, line, ios)
      if (ios /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: not a MESA history file: '//trim(history_file))
   end do
   call split_words(line, in_names, n_in)

   idx_model = find_col(in_names, 'model_number')
   idx_logg = find_col(in_names, 'log_g')
   idx_h1 = find_col(in_names, 'surface_h1')
   idx_he4 = find_col(in_names, 'surface_he4')

   idx_teff = find_col(in_names, 'log_Teff')
   teff_is_log = (idx_teff > 0)
   if (.not. teff_is_log) idx_teff = find_col(in_names, 'Teff')

   idx_radius = find_col(in_names, 'log_R')
   radius_is_log = (idx_radius > 0)
   if (.not. radius_is_log) idx_radius = find_col(in_names, 'radius')

   if (idx_teff < 1 .or. idx_logg < 1 .or. idx_radius < 1) then
      write(*, '(a)') 'colors_post: the history file needs log_Teff (or Teff), log_g and log_R (or radius) columns'
      call mesa_error(__FILE__, __LINE__, 'colors_post: missing history columns')
   end if

   have_composition = (idx_h1 > 0 .and. idx_he4 > 0)
   if (.not. have_composition) then
      if (fallback_m_div_h == unset_m_div_h) then
         write(*, '(a)') 'colors_post: the history file has no surface_h1 and surface_he4 columns to set [M/H].'
         write(*, '(a)') 'colors_post: add surface h1 and surface he4 to the history columns of the run,'
         write(*, '(a)') 'colors_post: or set fallback_m_div_h in &post_proc.'
         call mesa_error(__FILE__, __LINE__, 'colors_post: missing history columns')
      end if
      write(*, '(a,f0.3)') 'colors_post: no surface composition in the history file, using [M/H] = ', &
         max(min_m_div_h, min(max_m_div_h, fallback_m_div_h))
   end if

   ! write to a temporary file and move it into place at the end, so that an
   ! error part way through, or output_file = history_file, cannot lose data
   call make_parent_dir(output_file)
   temp_file = trim(output_file)//'.tmp'
   open(newunit=out_unit, file=temp_file, status='replace', action='write', iostat=ios)
   if (ios /= 0) call mesa_error(__FILE__, __LINE__, 'colors_post: cannot write '//temp_file)

   rewind(in_unit)
   do i = 1, num_header_lines + 2
      call read_line(in_unit, line, ios)
      if (i <= num_header_lines) write(out_unit, '(a)') line
   end do

   n_models = 0
   n_failed = 0
   do
      call read_line(in_unit, line, ios)
      if (ios /= 0) exit
      call split_words(line, words, n_words)
      if (n_words == 0) cycle
      n_models = n_models + 1
      if (n_words /= n_in) then
         write(*, '(a,i0,a,i0,a,i0)') 'colors_post: data row ', n_models, ' has ', n_words, &
            ' values, expected ', n_in
         close(out_unit, status='delete')
         call mesa_error(__FILE__, __LINE__, 'colors_post: bad history file')
      end if

      call colors_for_row(words, ierr)
      if (ierr /= 0) n_failed = n_failed + 1

      if (n_models == 1) call write_column_names()
      call write_row(words)
   end do

   close(in_unit)
   if (n_models == 0) then
      close(out_unit, status='delete')
      call mesa_error(__FILE__, __LINE__, 'colors_post: no data rows in '//trim(history_file))
   end if
   close(out_unit)
   call mv(temp_file, output_file)

   call free_colors_handle(handle)
   call colors_shutdown()

   write(*, '(a,i0,4a)') 'colors_post: wrote ', n_models, ' models from ', trim(history_file), &
      ' to ', trim(output_file)
   if (n_failed == 0) then
      write(*, '(a,i0,a)') 'colors_post: successfully post-processed all ', n_models, ' models'
   else
      write(*, '(a,i0,a,i0,a)') 'colors_post: colors calculation failed for ', n_failed, &
         ' of ', n_models, ' models'
   end if

contains

   subroutine colors_for_row(words, ierr)
      character(len=*), intent(in) :: words(:)
      integer, intent(out) :: ierr

      real(dp) :: teff, log_g, radius, x, z, m_div_h
      integer :: model_number

      teff = value_of(words, idx_teff)
      if (teff_is_log) teff = exp10(teff)
      log_g = value_of(words, idx_logg)
      radius = value_of(words, idx_radius)
      if (radius_is_log) radius = exp10(radius)
      radius = radius*rsun

      model_number = n_models
      if (idx_model > 0) model_number = nint(value_of(words, idx_model))

      ! same mapping of the surface Z/X onto the atmosphere-table metallicity
      ! axis as star/private/history.f90 uses for the inline colors columns
      if (have_composition) then
         x = value_of(words, idx_h1)
         z = 1d0 - x - value_of(words, idx_he4)
         if (x > 0d0 .and. z > 0d0) then
            m_div_h = log10((z/x)/rq% z_over_x_ref)
         else
            m_div_h = min_m_div_h
         end if
      else
         m_div_h = fallback_m_div_h
      end if
      m_div_h = max(min_m_div_h, min(max_m_div_h, m_div_h))

      vals(:) = -1d99
      call data_for_colors_history_columns(teff, log_g, radius, m_div_h, model_number, &
         handle, n_colors, col_names, vals, ierr)
      if (ierr /= 0) write(*, '(a,i0,a,i0)') 'colors_post: colors calculation failed for model ', &
         model_number, ', ierr = ', ierr
   end subroutine colors_for_row


   ! colors columns replace input columns of the same name, the rest are appended
   subroutine write_column_names()
      character(len=80) :: num
      integer :: j

      allocate(out_names(n_in + n_colors))
      out_names(1:n_in) = in_names
      n_out = n_in
      do j = 1, n_colors
         dest(j) = find_col(in_names, trim(col_names(j)))
         if (dest(j) < 1) then
            n_out = n_out + 1
            dest(j) = n_out
            out_names(n_out) = col_names(j)
         end if
      end do

      do j = 1, n_out
         write(num, '(i0)') j
         call write_word(num)
      end do
      write(out_unit, '(a)') ''
      do j = 1, n_out
         call write_word(out_names(j))
      end do
      write(out_unit, '(a)') ''
   end subroutine write_column_names


   subroutine write_row(words)
      character(len=*), intent(in) :: words(:)
      character(len=80) :: out_words(n_out)
      integer :: j

      out_words(1:n_in) = words
      do j = 1, n_colors
         write(out_words(dest(j)), '(es26.16e3)') vals(j)
      end do
      do j = 1, n_out
         call write_word(adjustl(out_words(j)))
      end do
      write(out_unit, '(a)') ''
   end subroutine write_row


   subroutine write_word(word)
      character(len=*), intent(in) :: word
      write(out_unit, '(a)', advance='no') repeat(' ', max(1, col_width - len_trim(word)))//trim(word)
   end subroutine write_word


   real(dp) function value_of(words, idx)
      character(len=*), intent(in) :: words(:)
      integer, intent(in) :: idx
      integer :: ios
      read(words(idx), *, iostat=ios) value_of
      if (ios /= 0) then
         close(out_unit, status='delete')
         call mesa_error(__FILE__, __LINE__, 'colors_post: cannot read '//trim(in_names(idx))//' value '//trim(words(idx)))
      end if
   end function value_of


   subroutine make_parent_dir(path)
      character(len=*), intent(in) :: path
      integer :: slash
      slash = index(path, '/', back=.true.)
      if (slash > 1) then
         if (.not. folder_exists(path(1:slash-1))) call mkdir(path(1:slash-1))
      end if
   end subroutine make_parent_dir


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


   integer function find_col(names, name) result(idx)
      character(len=*), intent(in) :: names(:)
      character(len=*), intent(in) :: name
      integer :: j

      idx = -1
      do j = 1, size(names)
         if (trim(names(j)) == name) then
            idx = j
            return
         end if
      end do
   end function find_col

end program colors_post
