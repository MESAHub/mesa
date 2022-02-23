! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
! ***********************************************************************

      module history_specs

      use star_private_def
      use star_history_def
      use const_def
      use chem_def
      use num_lib, only: linear_interp, find0

      implicit none

      public ! history.f90 uses most of this

      ! spacing between these must be larger than max number of nuclides
      integer, parameter :: idel = 10000

      integer, parameter :: center_xa_offset = idel
      integer, parameter :: surface_xa_offset = center_xa_offset + idel
      integer, parameter :: average_xa_offset = surface_xa_offset + idel
      integer, parameter :: category_offset = average_xa_offset + idel
      integer, parameter :: total_mass_offset = category_offset + idel
      integer, parameter :: log_total_mass_offset = total_mass_offset + idel
      integer, parameter :: log_average_xa_offset = log_total_mass_offset + idel
      integer, parameter :: log_center_xa_offset = log_average_xa_offset + idel
      integer, parameter :: log_surface_xa_offset = log_center_xa_offset + idel
      integer, parameter :: cz_max_offset = log_surface_xa_offset + idel
      integer, parameter :: cz_top_max_offset = cz_max_offset + idel
      integer, parameter :: max_eps_nuc_offset = cz_top_max_offset + idel
      integer, parameter :: c_log_eps_burn_offset = max_eps_nuc_offset + idel

      integer, parameter :: bc_offset = c_log_eps_burn_offset + idel
      integer, parameter :: abs_mag_offset = bc_offset + idel
      integer, parameter :: lum_band_offset = abs_mag_offset + idel
      integer, parameter :: log_lum_band_offset = lum_band_offset + idel
      
      integer, parameter :: raw_rate_offset = log_lum_band_offset + idel
      integer, parameter :: screened_rate_offset = raw_rate_offset + idel
      integer, parameter :: eps_nuc_rate_offset = screened_rate_offset + idel
      integer, parameter :: eps_neu_rate_offset = eps_nuc_rate_offset + idel

      integer, parameter :: start_of_special_cases = eps_neu_rate_offset + idel
      ! mixing and burning regions must be given the largest offsets
      ! so they can be distinguished from the other ones

      integer, parameter :: mixing_offset = start_of_special_cases
      integer, parameter :: mix_relr_offset = mixing_offset + idel
      integer, parameter :: burning_offset = mix_relr_offset + idel
      integer, parameter :: burn_relr_offset = burning_offset + idel

      !integer, parameter :: next_available_offset = burning_offset + idel


      contains



      recursive subroutine add_history_columns( &
            s, level, capacity, spec, history_columns_file, report, ierr)
         use utils_lib
         use utils_def
         use chem_def
         use chem_lib
         use const_def, only: mesa_dir
         use colors_def, only: bc_total_num_colors
         use colors_lib
         type (star_info), pointer :: s
         integer, intent(in) :: level
         integer, intent(inout) :: capacity
         integer, pointer :: spec(:)
         character (len=*), intent(in) :: history_columns_file
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         integer :: iounit, n, i, t, id, j, k, cnt, ii, nxt_spec, spec_err
         character (len=strlen) :: buffer, string, filename
         integer, parameter :: max_level = 20
         logical :: bad_item, special_case
         logical, parameter :: dbg = .false.

         include 'formats'

         if (level > max_level) then
            write(*,*) 'too many levels of nesting for log column files', level
            ierr = -1
            return
         end if

         ierr = 0
         spec_err = 0

         ! first try local directory
         filename = history_columns_file
         if (len_trim(filename) == 0) filename = 'history_columns.list'
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in star/defaults
            filename = trim(mesa_dir) // '/star/defaults/' // trim(filename)
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then ! look for log_columns.list
               filename = 'log_columns.list'
               ierr = 0
               open(newunit=iounit, &
                     file=trim(filename), action='read', status='old', iostat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to open ' // trim(history_columns_file)
                  return
               end if
               write(*,*) 'please rename log_columns.list to history_columns.list'
            end if
         end if

         if (dbg) then
            write(*,'(A)')
            write(*,*) 'history_columns_file <' // trim(filename) // '>'
            write(*,'(A)')
         end if

         call count_specs

         n = 0
         i = 0

         do

            t = token(iounit, n, i, buffer, string)
            if (t == eof_token) exit
            if (t /= name_token) then
               call error; return
            end if

            if (string == 'include') then
               t = token(iounit, n, i, buffer, string)
               if (t /= string_token) then
                  call error; return
               end if
               call add_history_columns(s, level+1, capacity, spec, string, report, spec_err)
               if (spec_err /= 0) then
                  write(*,*) 'failed for included log columns list ' // trim(string)
                  call error; return
               end if
               call count_specs
               cycle
            end if

            spec_err = 0

            if (string == 'add_center_abundances') then
               call do_abundances(center_xa_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_surface_abundances') then
               call do_abundances(surface_xa_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_average_abundances') then
               call do_abundances(average_xa_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_log_center_abundances') then
               call do_abundances(log_center_xa_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_log_surface_abundances') then
               call do_abundances(log_surface_xa_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_log_average_abundances') then
               call do_abundances(log_average_xa_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if
            
            if (string == 'add_total_mass') then
               call do_abundances(total_mass_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_log_total_mass') then
               call do_abundances(log_total_mass_offset, spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_bc') then
               call do_colors(bc_offset,'bc_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_abs_mag') then
               call do_colors(abs_mag_offset,'abs_mag_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if
            
            if (string == 'add_lum_band') then
               call do_colors(lum_band_offset,'lum_band_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if
            
            if (string == 'add_log_lum_band') then
               call do_colors(log_lum_band_offset,'log_lum_band_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_raw_rates') then
               call do_rate(raw_rate_offset,'raw_rate_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_screened_rates') then
               call do_rate(screened_rate_offset,'screened_rate_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_eps_nuc_rates') then
               call do_rate(eps_nuc_rate_offset,'eps_nuc_rate_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            if (string == 'add_eps_neu_rates') then
               call do_rate(eps_neu_rate_offset,'eps_neu_rate_', spec_err)
               if (spec_err /= 0) then
                  call error; return
               end if
               call count_specs
               cycle
            end if

            nxt_spec = do1_history_spec( &
               iounit, t, n, i, string, buffer, special_case, report, spec_err)
            if (spec_err /= 0) then
               ierr = spec_err
            else
               if (.not. special_case) then
                  call insert_spec(nxt_spec, string, spec_err)
                  if (spec_err /= 0) then
                     write(*,*) 'failed for history item ' // trim(string)
                     ierr = -1; cycle
                  end if
               else if (nxt_spec == h_mixing_regions) then
                  t = token(iounit, n, i, buffer, string)
                  if (t /= name_token) then
                     ierr = -1; cycle
                  end if
                  read(string,fmt=*,iostat=spec_err) cnt
                  if (spec_err /= 0 .or. cnt <= 0 .or. cnt > 1000) then
                     write(*,*) 'bad integer count for mixing regions: ' // trim(string)
                     ierr = -1; cycle
                  end if
                  do ii=1,2*cnt
                     call insert_spec(mixing_offset + ii, string, spec_err)
                     if (spec_err /= 0) then
                        call error; return
                     end if
                  end do
               else if (nxt_spec == h_mix_relr_regions) then
                  t = token(iounit, n, i, buffer, string)
                  if (t /= name_token) then
                     ierr = -1; cycle
                  end if
                  read(string,fmt=*,iostat=spec_err) cnt
                  if (spec_err /= 0 .or. cnt <= 0 .or. cnt > 1000) then
                     write(*,*) 'bad integer count for mix_relr regions: ' // trim(string)
                     ierr = -1; cycle
                  end if
                  do ii=1,2*cnt
                     call insert_spec(mix_relr_offset + ii, string, spec_err)
                     if (spec_err /= 0) then
                        call error; return
                     end if
                  end do
               else if (nxt_spec == h_burning_regions) then
                  t = token(iounit, n, i, buffer, string)
                  if (t /= name_token) then
                     ierr = -1; cycle
                  end if
                  read(string,fmt=*,iostat=spec_err) cnt
                  if (spec_err /= 0 .or. cnt <= 0 .or. cnt > 1000) then
                     write(*,*) 'bad integer count following burning regions: ' // trim(string)
                     ierr = -1
                  end if
                  do ii=1,2*cnt
                     call insert_spec(burning_offset + ii, string, spec_err)
                     if (spec_err /= 0) then
                        call error; return
                     end if
                  end do
               else if (nxt_spec == h_burn_relr_regions) then
                  t = token(iounit, n, i, buffer, string)
                  if (t /= name_token) then
                     ierr = -1; cycle
                  end if
                  read(string,fmt=*,iostat=spec_err) cnt
                  if (spec_err /= 0 .or. cnt <= 0 .or. cnt > 1000) then
                     write(*,*) 'bad integer count following burn_relr regions: ' // trim(string)
                     ierr = -1
                  end if
                  do ii=1,2*cnt
                     call insert_spec(burn_relr_offset + ii, string, spec_err)
                     if (spec_err /= 0) then
                        call error; return
                     end if
                  end do
               else
                  write(*,*) 'confusion in history specs'
                  ierr = -1
               end if
            end if

         end do

         if (dbg) write(*,*) 'finished ' // trim(filename)

         close(iounit)

         if (dbg) then
            write(*,'(A)')
            write(*,*) 'done add_history_columns ' // trim(filename)
            write(*,'(A)')
         end if


         contains


         subroutine count_specs
            integer :: i
            j = 1
            do i=1, capacity
               if (spec(i) == 0) then
                  j = i; exit
               end if
            end do
         end subroutine count_specs


         subroutine make_room(ierr)
            integer, intent(out) :: ierr
            if (j < capacity) return
            capacity = 50 + (3*capacity)/2
            call realloc_integer(spec,capacity,ierr)
            spec(j+1:capacity) = 0
         end subroutine make_room


         subroutine insert_spec(c, name, ierr)
            integer, intent(in) :: c
            character (len=*) :: name
            integer, intent(out) :: ierr
            integer :: i
            include 'formats'
            do i=1,j-1
               if (spec(i) == c) return
            end do
            call make_room(ierr)
            if (ierr /= 0) return
            spec(j) = c
            if (dbg) write(*,2) trim(name), spec(j)
            j = j+1
         end subroutine insert_spec


         subroutine do_abundances(offset, ierr)
            integer, intent(in) :: offset
            integer, intent(out) :: ierr
            integer :: k
            ierr = 0
            do k=1,s% species
               call insert_spec( &
                  offset + s% chem_id(k), &
                  chem_isos% name(s% chem_id(k)), ierr)
               if (ierr /= 0) return
            end do
         end subroutine do_abundances

         subroutine do_colors(offset,prefix,ierr)
            integer, intent(in) :: offset
            character(len=*) :: prefix
            integer, intent(out) :: ierr
            integer :: k
            ierr = 0
            do k=1,bc_total_num_colors
               call insert_spec( &
                  offset + k,trim(prefix)//trim(get_bc_name_by_id(k,ierr)), ierr)
               if (ierr /= 0) return
            end do
         end subroutine do_colors

         subroutine do_rate(offset,prefix,ierr) ! raw_rate, screened_rate, eps_nuc_rate, eps_neu_rate
            use rates_def, only: reaction_name
            integer, intent(in) :: offset
            character(len=*) :: prefix
            integer, intent(out) :: ierr
            integer :: k
            ierr = 0
            do k=1,s% num_reactions
               call insert_spec( &
                  offset + k,trim(prefix)//trim(reaction_name(k)), ierr)
               if (ierr /= 0) return
            end do
         end subroutine do_rate


         subroutine error
            ierr = -1
            close(iounit)
         end subroutine error


      end subroutine add_history_columns


      integer function do1_history_spec( &
            iounit, t, n, i, string, buffer, special_case, report, ierr) result(spec)

         use utils_lib
         use utils_def
         use chem_def
         use chem_lib
         integer :: iounit, t, n, i

         character (len=*) :: string, buffer
         logical, intent(out) :: special_case
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         integer :: id

         ierr = 0
         spec = -1
         special_case = .false.

         if (string=='raw_rate') then
            call do1_rate(raw_rate_offset)

         else if (string=='screened_rate') then
            call do1_rate(screened_rate_offset)

         else if (string=='eps_nuc_rate') then
            call do1_rate(eps_nuc_rate_offset)

         else if (string=='eps_neu_rate') then
            call do1_rate(eps_neu_rate_offset)

         else if (string=='abs_mag') then
            call do1_color(abs_mag_offset)

         else if (string=='bc') then
            call do1_color(bc_offset)
            
         else if (string=='lum_band') then
            call do1_color(lum_band_offset)

        else if (string=='log_lum_band') then
            call do1_color(log_lum_band_offset)

         else if (string == 'center') then
            call do1_nuclide(center_xa_offset)

         else if (string == 'surface') then
            call do1_nuclide(surface_xa_offset)

         else if (string == 'average') then
            call do1_nuclide(average_xa_offset)

         else if (string == 'total_mass') then
            call do1_nuclide(total_mass_offset)

         else if (string == 'log_total_mass') then
            call do1_nuclide(log_total_mass_offset)

         else if (string == 'log_average') then
            call do1_nuclide(log_average_xa_offset)

         else if (string == 'log_center') then
            call do1_nuclide(log_center_xa_offset)

         else if (string == 'log_surface') then
            call do1_nuclide(log_surface_xa_offset)

         else if (string == 'max_eps_nuc_log_xa') then
            call do1_nuclide(max_eps_nuc_offset)

         else if (string == 'cz_top_log_xa') then
            call do1_nuclide(cz_top_max_offset)

         else if (string == 'cz_log_xa') then
            call do1_nuclide(cz_max_offset)

         else if (string == 'c_log_eps_burn') then
            call do1_rates_category(c_log_eps_burn_offset)

         else
            id = do_get_history_id(string)
            if (id > 0) then
               spec = id
               if (id == h_mixing_regions .or. &
                   id == h_mix_relr_regions .or. &
                   id == h_burning_regions .or. &
                   id == h_burn_relr_regions) then
                  special_case = .true.
               end if
               return
            end if
            id = rates_category_id(string)
            if (id > 0) then
               spec = category_offset + id
               return
            end if
            if (report) write(*,*) 'bad history list name: ' // trim(string)
            ierr = -1

         end if


         contains


         subroutine do1_nuclide(offset)
            integer, intent(in) :: offset
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token) then
               ierr = -1; return
            end if
            id = chem_get_iso_id(string)
            if (id > 0) then
               spec = offset + id
               return
            end if
            write(*,*) 'bad iso name: ' // trim(string)
            ierr = -1
         end subroutine do1_nuclide


         subroutine do1_rates_category(offset)
            integer, intent(in) :: offset
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token) then
               ierr = -1; return
            end if
            id = rates_category_id(string)
            if (id > 0) then
               spec = offset + id
               return
            end if
            write(*,*) 'bad rates category name: ' // trim(string)
            ierr = -1
         end subroutine do1_rates_category

         subroutine do1_color(offset)
            use colors_lib, only : get_bc_id_by_name
            integer, intent(in) :: offset
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token) then
               ierr = -1; return
            end if
            id = get_bc_id_by_name(string,ierr)
            if (ierr/=0) return
            if (id > 0) then
               spec = offset + id
               return
            end if
            write(*,*) 'bad filter name: ' // trim(string)
            ierr = -1
         end subroutine do1_color

         subroutine do1_rate(offset)
            use rates_lib, only: rates_reaction_id
            integer, intent(in) :: offset
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token) then
               ierr = -1; return
            end if
            id = rates_reaction_id(string)
            if (ierr/=0) return
            if (id > 0) then
               spec = offset + id
               return
            end if
            write(*,*) 'bad filter name: ' // trim(string)
            ierr = -1
         end subroutine do1_rate


      end function do1_history_spec


      subroutine set_history_columns(id, history_columns_file, report, ierr)
         use utils_lib, only: realloc_integer
         integer, intent(in) :: id
         character (len=*), intent(in) :: history_columns_file
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: capacity, cnt, i
         logical, parameter :: dbg = .false.
         integer, pointer :: old_history_column_spec(:) => null()
         character (len=strlen) :: fname
         logical :: history_file_exists
         if (dbg) write(*,*) 'set_history_columns'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         old_history_column_spec => null()
         if (associated(s% history_column_spec)) old_history_column_spec => s% history_column_spec
         nullify(s% history_column_spec)
         capacity = 100 ! will increase if needed
         allocate(s% history_column_spec(capacity), stat=ierr)
         if (ierr /= 0) return
         s% history_column_spec(:) = 0
         call add_history_columns( &
               s, 1, capacity, s% history_column_spec, history_columns_file, report, ierr)
         if (ierr /= 0) then
            if (associated(old_history_column_spec)) deallocate(old_history_column_spec)
            return
         end if
         ! delete trailing 0's
         cnt = capacity+1
         do i=1, capacity
            if (s% history_column_spec(i) == 0) then
               cnt = i; exit
            end if
         end do
         capacity = cnt-1
         call realloc_integer(s% history_column_spec, capacity, ierr)
         if (ierr /= 0) return
         if (associated(old_history_column_spec)) then
            ! check that haven't changed the cols specs for an existing log file
            fname = trim(s% log_directory) // '/' // trim(s% star_history_name)
            inquire(file=trim(fname), exist=history_file_exists)
            if (history_file_exists) then
               if (capacity /= size(old_history_column_spec)) then
                  ierr = -1
                  write(*,*) 'new size of history col specs', capacity
                  write(*,*) 'old size of history col specs', size(old_history_column_spec)
               else
                  do i=1,capacity
                     if (old_history_column_spec(i) /= s% history_column_spec(i)) then
                        write(*,*) 'change in history col spec', &
                           i, old_history_column_spec(i), s% history_column_spec(i)
                        ierr = -1
                     end if
                  end do
               end if
               if (ierr /= 0) then
                  write(*,*) 'ERROR: cannot change history columns when have an existing history file'
                  write(*,*) 'please delete the history file or go back to previous history columns list'
               end if
            end if
            deallocate(old_history_column_spec)
            if (ierr /= 0) return
         end if

      end subroutine set_history_columns


      end module history_specs

