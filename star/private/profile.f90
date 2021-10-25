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
! ***********************************************************************

      module profile

      use star_private_def
      use star_profile_def
      use const_def
      use profile_getval

      implicit none

      private
      public :: write_profile_info, set_profile_columns, &
         get_profile_val, do_save_profiles, &
         do_get_data_for_profile_columns, do_get_num_standard_profile_columns

      ! model log priorities
      integer, parameter :: delta_priority = 1
      integer, parameter :: phase_priority = 2


      contains


      recursive subroutine add_profile_columns( &
            s, level, capacity, spec, profile_columns_file, report, ierr)
         use utils_lib
         use utils_def
         use chem_def
         use chem_lib
         use const_def, only: mesa_dir
         type (star_info), pointer :: s
         integer, intent(in) :: level
         integer, intent(inout) :: capacity
         integer, pointer :: spec(:)
         character (len=*), intent(in) :: profile_columns_file
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         integer :: iounit, n, i, t, id, j, k, num, nxt_spec, spec_err
         character (len=strlen) :: buffer, string, filename
         integer, parameter :: max_level = 20

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
         filename = profile_columns_file
         if (len_trim(filename) == 0) filename = 'profile_columns.list'
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in star/defaults
            filename = trim(mesa_dir) // '/star/defaults/' // trim(filename)
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(filename)
               return
            end if
         end if

         if (dbg) then
            write(*,*)
            write(*,*) 'profile_columns_file <' // trim(filename) // '>'
         end if

         call count_specs

         n = 0
         i = 0

         do
            t = token(iounit, n, i, buffer, string)
            if (dbg) write(*,*) 'token', t
            if (t == eof_token) then
               if (dbg) write(*,*) 'eof token'
               exit
            end if
            if (t /= name_token) then
               call error; return
            end if

            select case(string)

            case ('include')
               t = token(iounit, n, i, buffer, string)
               if (dbg) write(*,*) 'include file token', t
               if (dbg) write(*,*) 'include file string len', len_trim(string)
               if (t /= string_token) then
                  call error; return
               end if
               if (dbg) write(*,*) 'include file <' // trim(string) // '>'
               call add_profile_columns(s, level+1, capacity, spec, string, report, spec_err)
               if (spec_err /= 0) then
                  write(*,*) 'failed for included profile columns list: ' // trim(string)
                  ierr = -1; call error; return
               end if
               call count_specs

            case ('add_abundances')
               ! add all of the isos that are in the current net
               call insert_spec(add_abundances, 'add_abundances', spec_err)
               if (spec_err /= 0) then
                  ierr = -1; call error; return
               end if

            case ('add_log_abundances')
               ! add logs of all of the isos that are in the current net
               call insert_spec(add_log_abundances, 'add_log_abundances', spec_err)
               if (spec_err /= 0) then
                  ierr = -1; call error; return
               end if

            case ('add_reaction_categories') ! add all the reaction categories
               do k = 1, num_categories
                  call insert_spec(category_offset + k, category_name(k), spec_err)
                  if (spec_err /= 0) then
                     ierr = -1; call error; return
                  end if
               end do

            case default
               spec_err = 0
               nxt_spec = do1_profile_spec(iounit, n, i, string, buffer, report, spec_err)
               if (spec_err /= 0) then
                  ierr = spec_err
               else
                  if (nxt_spec > 0) then
                     call insert_spec(nxt_spec, string, spec_err)
                     if (spec_err /= 0) ierr = spec_err
                  else
                     if (report) &
                        write(*,*) 'failed to recognize item for profile columns: ' // trim(string)
                     ierr = -1
                  end if
               end if
            end select

         end do

         if (dbg) write(*,*) 'finished ' // trim(filename)

         close(iounit)

         if (dbg) then
            write(*,*)
            write(*,*) 'done add_profile_columns ' // trim(filename)
            write(*,*)
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


         subroutine error
            ierr = -1
            close(iounit)
         end subroutine error


      end subroutine add_profile_columns


      subroutine set_profile_columns(id, profile_columns_file, report, ierr)
         use utils_lib, only: realloc_integer
         integer, intent(in) :: id
         character (len=*), intent(in) :: profile_columns_file
         logical, intent(in) :: report
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: capacity, cnt, i
         logical, parameter :: dbg = .false.
         if (dbg) write(*,*) 'set_profile_columns'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (associated(s% profile_column_spec)) deallocate(s% profile_column_spec)
         capacity = 100 ! will increase if needed
         allocate(s% profile_column_spec(capacity), stat=ierr)
         if (ierr /= 0) return
         s% profile_column_spec(:) = 0
         call add_profile_columns( &
            s, 1, capacity, s% profile_column_spec, profile_columns_file, report, ierr)
         if (ierr /= 0) return
         ! delete trailing 0's
         cnt = capacity+1
         do i=1, capacity
            if (s% profile_column_spec(i) == 0) then
               cnt = i; exit
            end if
            if (dbg) write(*,*) 'profile col', i, s% profile_column_spec(i)
         end do
         call realloc_integer(s% profile_column_spec, cnt-1, ierr)
         if (dbg) write(*,*) 'num profile columns', cnt-1
         if (dbg) stop 'debug: set_profile_columns'
      end subroutine set_profile_columns
      
      
      integer function do_get_num_standard_profile_columns(s) ! not inluding extra profile columns
         use star_def, only: star_info
         type (star_info), pointer :: s
         integer :: numcols, j, num_specs
         numcols = 0
         if (.not. associated(s% profile_column_spec)) then
            num_specs = 0
         else
            num_specs = size(s% profile_column_spec, dim=1)
         end if
         do j = 1, num_specs
            if (s% profile_column_spec(j) == add_abundances .or. &
                s% profile_column_spec(j) == add_log_abundances) then
               numcols = numcols + s% species
            else
               numcols = numcols + 1
            end if
         end do
         do_get_num_standard_profile_columns = numcols
      end function do_get_num_standard_profile_columns


      subroutine do_get_data_for_profile_columns(s, nz, &
            names, vals, is_int, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         character (len=maxlen_profile_column_name), pointer :: names(:) ! (num_profile_columns)
         real(dp), pointer :: vals(:,:) ! (nz,num_profile_columns)
         logical, pointer :: is_int(:) ! (num_profile_columns) true iff the values in the column are integers
         integer, intent(out) :: ierr
         character (len=0) :: fname
         logical, parameter :: write_flag = .false.
         call do_profile_info(s, fname, &
            write_flag, names, vals, is_int, ierr)
      end subroutine do_get_data_for_profile_columns


      subroutine write_profile_info(s, fname, ierr)
         use chem_def
         use net_def ! categories
         use rates_def, only: i_rate
         type (star_info), pointer :: s
         character (len=*) :: fname
         integer, intent(out) :: ierr
         character (len=maxlen_profile_column_name), pointer :: names(:) ! (num_profile_columns)
         real(dp), pointer :: vals(:,:) ! (nz,num_profile_columns)
         logical, pointer :: is_int(:)
         logical, parameter :: write_flag = .true.
         names => null()
         vals => null()
         is_int => null()
         call do_profile_info(s, fname, &
            write_flag, names, vals, is_int, ierr)
      end subroutine write_profile_info


      subroutine do_profile_info(s, fname, &
            write_flag, names, vals, is_int, ierr)
         use chem_def
         use net_def ! categories
         use rates_def, only: i_rate
         use ctrls_io, only: write_controls
         use write_model, only: do_write_model
         use pulse, only: export_pulse_data
         use math_lib, only: math_backend
         
         type (star_info), pointer :: s
         character (len=*) :: fname
         logical, intent(in) :: write_flag
         character (len=maxlen_profile_column_name), pointer :: names(:) ! (num_profile_columns)
         real(dp), pointer :: vals(:,:) ! (nz,num_profile_columns)
         logical, pointer :: is_int(:)
         integer, intent(out) :: ierr

         real(dp) :: msum, mstar, dt, Lnuc, frac
         integer :: io, i, j, jj, nz, col, k, kk, n, species, &
            h1, he4, num_specs, numcols, num_extra_cols, num_extra_header_items, num_digits
         integer, pointer :: chem_id(:)
         logical, parameter :: dbg = .false.
         character (len=strlen) :: fname1, dbl_fmt, int_fmt, txt_fmt, fname_out, fstring, str
         character (len=maxlen_profile_column_name), pointer :: &
            extra_col_names(:), extra_header_item_names(:)
         real(dp), pointer :: extra_col_vals(:,:), extra_header_item_vals(:)

         include "formats"

         dbl_fmt = s% profile_dbl_format
         int_fmt = s% profile_int_format
         txt_fmt = s% profile_txt_format

         ierr = 0
         nullify(extra_col_names, extra_col_vals)

         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         mstar = s% mstar
         dt = s% dt
         if (.not. associated(s% profile_column_spec)) then
            num_specs = 0
         else
            num_specs = size(s% profile_column_spec, dim=1)
         end if

         if (num_specs == 0) then
            write(*,*) 'WARNING: do not have any output specified for profiles.'
            return
         end if
      
         numcols = do_get_num_standard_profile_columns(s)

         num_extra_cols = s% how_many_extra_profile_columns(s% id)
         if (num_extra_cols > 0) then
            allocate( &
               extra_col_names(num_extra_cols), extra_col_vals(nz,num_extra_cols), stat=ierr)
            if (ierr /= 0) return
            extra_col_names(1:num_extra_cols) = 'unknown'
            extra_col_vals(1:nz,1:num_extra_cols)  = -1d99
            call s% data_for_extra_profile_columns( &
               s% id, num_extra_cols, nz, extra_col_names, extra_col_vals, ierr)
            if (ierr /= 0) then
               deallocate(extra_col_names, extra_col_vals)
               return
            end if
            do i=1,num_extra_cols
               if(trim(extra_col_names(i))=='unknown') then
                  write(*,*) "Warning empty profile name for extra_profile_column ",i
               end if
            end do
         end if

         if (.not. write_flag) then

            if (associated(names)) then
               if (size(names,dim=1) < numcols+num_extra_cols) then
                  write(*,2) 'size(names,dim=1)', size(names,dim=1)
                  write(*,2) 'numcols+num_extra_cols', numcols+num_extra_cols
                  write(*,2) 'numcols', numcols
                  write(*,2) 'num_extra_cols', num_extra_cols
                  write(*,*) 'bad size for names in do_profile_info'
                  ierr = -1
                  return
               end if
            else
               write(*,*) 'failed to provide names array for do_profile_info'
               ierr = -1
               return
            end if

            if (associated(vals)) then
               if (size(vals,dim=1) < nz) then
                  write(*,2) 'size(vals,dim=1)', size(vals,dim=1)
                  write(*,2) 'nz', nz
                  write(*,*) 'bad size dim=1 for vals in do_profile_info'
                  ierr = -1
                  return
               end if
               if (size(vals,dim=2) < numcols+num_extra_cols) then
                  write(*,2) 'size(vals,dim=2)', size(vals,dim=2)
                  write(*,2) 'numcols+num_extra_cols', numcols+num_extra_cols
                  write(*,2) 'numcols', numcols
                  write(*,2) 'num_extra_cols', num_extra_cols
                  write(*,*) 'bad size dim=1 for names in do_profile_info'
                  ierr = -1
                  return
               end if
            else
               write(*,*) 'failed to provide vals array for do_profile_info'
               ierr = -1
               return
            end if

         end if

         if (write_flag) then

            if (len_trim(s% profile_data_header_suffix) == 0) then
               fname1 = fname
            else
               fname1 = trim(fname) // s% profile_data_header_suffix
            end if
            open(newunit=io, file=trim(fname1), action='write', status='replace', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(fname1)
               return
            end if

            num_extra_header_items = s% how_many_extra_profile_header_items(s% id)
            if (num_extra_header_items > 0) then
               allocate( &
                  extra_header_item_names(num_extra_header_items), &
                  extra_header_item_vals(num_extra_header_items), stat=ierr)
               if (ierr /= 0) then
                  return
               end if
               extra_header_item_names(1:num_extra_header_items) = 'unknown'
               extra_header_item_vals(1:num_extra_header_items)  = -1d99
               call s% data_for_extra_profile_header_items( &
                  s% id, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
               if (ierr /= 0) then
                  deallocate(extra_header_item_names, extra_header_item_vals)
                  return
               end if
               do i=1,num_extra_header_items
                  if(trim(extra_header_item_names(i))=='unknown') then
                     write(*,*) "Warning empty profile name for extra_profile_header ",i
                  end if
               end do
            end if

            do i=1, 3
               col = 0
               call do_integer(i, 'model_number', s% model_number)
               call do_integer(i, 'num_zones', s% nz)
               call do_val(i, 'initial_mass', s% initial_mass)
               call do_val(i, 'initial_z', s% initial_z)
               call do_val(i, 'star_age', s% star_age)
               call do_val(i, 'time_step', s% time_step)

               if (s% M_center /= 0d0) &
                  call do_val(i, 'M_center', s% M_center)
               if (s% v_center /= 0d0) &
                  call do_val(i, 'v_center', s% v_center)
               if (s% R_center /= 0d0) &
                  call do_val(i, 'R_center', s% R_center)
               if (s% L_center /= 0d0) &
                  call do_val(i, 'L_center', s% L_center)

               call do_val(i, 'Teff', s% Teff)
               call do_val(i, 'photosphere_L', s% photosphere_L)
               call do_val(i, 'photosphere_r', s% photosphere_r)

               call do_val(i, 'center_eta', s% center_degeneracy)
               call do_val(i, 'center_h1', s% center_h1)
               call do_val(i, 'center_he3', s% center_he3)
               call do_val(i, 'center_he4', s% center_he4)
               call do_val(i, 'center_c12', s% center_c12)
               call do_val(i, 'center_n14', s% center_n14)
               call do_val(i, 'center_o16', s% center_o16)
               call do_val(i, 'center_ne20', s% center_ne20)
               call do_val(i, 'star_mass', s% star_mass)
               call do_val(i, 'star_mdot', s% star_mdot)
               call do_val(i, 'star_mass_h1', s% star_mass_h1)
               call do_val(i, 'star_mass_he3', s% star_mass_he3)
               call do_val(i, 'star_mass_he4', s% star_mass_he4)
               call do_val(i, 'star_mass_c12', s% star_mass_c12)
               call do_val(i, 'star_mass_n14', s% star_mass_n14)
               call do_val(i, 'star_mass_o16', s% star_mass_o16)
               call do_val(i, 'star_mass_ne20', s% star_mass_ne20)
               call do_val(i, 'he_core_mass', s% he_core_mass)
               call do_val(i, 'co_core_mass', s% co_core_mass)
               call do_val(i, 'fe_core_mass', s% fe_core_mass)
               call do_val(i, 'neutron_rich_core_mass', s% neutron_rich_core_mass)
               call do_val(i, 'dynamic_time', s% dynamic_timescale)
               call do_val(i, 'kh_timescale', s% kh_timescale)
               call do_val(i, 'nuc_timescale', s% nuc_timescale)

               call do_val(i, 'power_nuc_burn', s% power_nuc_burn)
               call do_val(i, 'power_h_burn', s% power_h_burn)
               call do_val(i, 'power_he_burn', s% power_he_burn)
               call do_val(i, 'power_neu', s% power_neutrinos)

               call do_val(i, 'burn_min1', s% burn_min1)
               call do_val(i, 'burn_min2', s% burn_min2)

               call do_val(i, 'time_seconds', s% time)

               call do_string(i, 'version_number', version_number)
               
               if (s% profile_header_include_sys_details) then ! make this optional
                  call do_string(i, 'compiler', compiler_name)
                  call do_string(i, 'build', compiler_version_name)
                  call do_string(i, 'MESA_SDK_version', mesasdk_version_name)
                  call do_string(i, 'math_backend', math_backend)
                  call do_string(i, 'date', date)
               end if

               call do_val(i, 'msun', msun)
               call do_val(i, 'rsun', rsun)
               call do_val(i, 'lsun', lsun)

               do j=1,num_extra_header_items
                 call do_val(i, extra_header_item_names(j), extra_header_item_vals(j))
               end do

               write(io, *)

            end do

            write(io, *)

            if (num_extra_header_items > 0) &
               deallocate(extra_header_item_names, extra_header_item_vals)

         end if

         do i = 1, 3

            if (i==3) then
               if (s% max_num_profile_zones > 1) then
                  n = min(nz, s% max_num_profile_zones)
               else
                  n = nz
               end if
               if (write_flag .and. len_trim(s% profile_data_header_suffix) > 0) then
                  close(io)
                  open(newunit=io, file=trim(fname), &
                     action='write', status='replace', iostat=ierr)
                  if (ierr /= 0) then
                     write(*,*) 'failed to open ' // trim(fname)
                     return
                  end if
               end if
            else
               n = 1
            end if

            do k=1, n
               col = 0
               if (n > 1 .and. n < nz) then
                  kk = floor(1.5d0 + dble(nz-1)*dble(k-1)/dble(n-1))
               else
                  kk = k
               end if
               do j = 1, num_specs
                  if (s% profile_column_spec(j) == add_abundances .or. &
                      s% profile_column_spec(j) == add_log_abundances) then
                     do jj = 1, species
                        col = col+1
                        call do_abundance_col(i, j, jj, kk)
                     end do
                  else
                     col = col+1
                     call do_col(i, j, kk)
                  end if
               end do
               do j=1,num_extra_cols
                  col = col+1
                  call do_extra_col(i, j, kk)
               end do
               if (write_flag) write(io, *)
            end do

         end do

         if (associated(extra_col_vals)) deallocate(extra_col_vals)
         if (associated(extra_col_names)) deallocate(extra_col_names)

         if (write_flag) then

            close(io)

            s% most_recent_profile_filename = trim(fname1)

            write(str,'(a)') 'save ' // trim(fname1)
            write(*,'(a)', advance='no') trim(str)
            call write_to_extra_terminal_output_file(s, str, .false.)

            if (s% write_pulse_data_with_profile) then
               fname_out = trim(fname) // '.' // trim(s% pulse_data_format)
               call export_pulse_data(s%id, s%pulse_data_format, fname_out, &
                    s%add_center_point_to_pulse_data, s%keep_surface_point_for_pulse_data, &
                    s%add_atmosphere_to_pulse_data, ierr)
               if (ierr /= 0) then
                  write(*,*) 'save_pulsation_info failed to open ' // trim(fname_out)
                  ierr = 0
               else
                  write(*,'(a)', advance='no') ' ' // trim(fname_out)
               end if
            end if

            if (s% write_model_with_profile) then
               fname_out = s% model_data_filename
               call do_write_model(s% id, fname_out, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to open ' // trim(fname_out)
                  ierr = 0
               else
                  write(*,'(a)', advance='no') ' ' // trim(fname_out)
               end if
            end if

            if (s% write_controls_info_with_profile) then
               fname_out = s% model_controls_filename
               call write_controls(s, fname_out, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to write ' // trim(fname_out)
                  ierr = 0
               else
                  s% most_recent_controls_filename = trim(fname_out)
               end if
                  write(*,'(a)', advance='no') ' ' // trim(fname_out)
            end if ! write_controls_info_with_profile

            num_digits = 1 + log10(dble(max(1,s% model_number)))
            write(fstring,'( "(a,i",i2.2,".",i2.2,")" )') num_digits, num_digits

            write(str,fstring) ' for model ', s% model_number
            write(*,'(a)') trim(str)
            call write_to_extra_terminal_output_file(s, str, .false.)

         end if


         contains

           
           subroutine do_string(pass, col_name, val)
             integer, intent(in) :: pass
             character (len=*), intent(in) :: col_name, val
             character(len=strlen) :: my_val
             col = col+1
             my_val = '"'//trim(val)//'"'
             if (pass == 1) then
                write(io, fmt=int_fmt, advance='no') col
             else if (pass == 2) then
                write(io, fmt=txt_fmt, advance='no') trim(col_name)
             else if (pass == 3) then
                write(io, fmt=txt_fmt, advance='no') adjustr(trim('"'//trim(val)//'"'))
             end if
           end subroutine do_string


         subroutine do_integer(pass, col_name, val)
            integer, intent(in) :: pass
            character (len=*), intent(in) :: col_name
            integer, intent(in) :: val
            col = col+1
            if (pass == 1) then
               write(io, fmt=int_fmt, advance='no') col
            else if (pass == 2) then
               write(io, fmt=txt_fmt, advance='no') trim(col_name)
            else if (pass == 3) then
               write(io, fmt=int_fmt, advance='no') val
            end if
         end subroutine do_integer


         subroutine do_val(pass, col_name, val)
            integer, intent(in) :: pass
            character (len=*), intent(in) :: col_name
            real(dp), intent(in) :: val
            real(dp) :: v
            include 'formats'
            col = col+1
            if (pass == 1) then
               write(io, fmt=int_fmt, advance='no') col
            else if (pass == 2) then
               write(io, fmt=txt_fmt, advance='no') trim(col_name)
            else if (pass == 3) then
               v = val
               if (is_bad_num(v)) then
                  write(*,1) 'bad value for ' // trim(col_name), v
                  if (s% stop_for_bad_nums) stop 'profile do_val'
                  v = 0
               end if
               write(io, fmt=dbl_fmt, advance='no') v
            end if
         end subroutine do_val


         subroutine do_extra_col(pass, j, k)
            use rates_def
            integer, intent(in) :: pass, j, k
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') col
            else if (pass == 2) then
               if (write_flag) then
                  write(io, fmt=txt_fmt, advance='no') trim(extra_col_names(j))
               else
                  names(col) = trim(extra_col_names(j))
               end if
            else if (pass == 3) then
               if (write_flag) then
                  write(io, fmt=dbl_fmt, advance='no') extra_col_vals(k,j)
               else
                  vals(k,col) = extra_col_vals(k,j)
                  is_int(col) = .false.
               end if
            end if
         end subroutine do_extra_col


         subroutine do_col(pass, j, k)
            use rates_def
            use profile_getval, only: getval_for_profile
            integer, intent(in) :: pass, j, k
            integer :: i, c, ii, int_val
            real(dp) :: val, cno, z, dr, eps, eps_alt
            logical :: int_flag
            character (len=128) :: col_name
            logical, parameter :: dbg = .false.
            include 'formats'
            c = s% profile_column_spec(j)
            val = 0; int_val = 0
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') col
            else if (pass == 2) then
               if (c > extra_offset) then
                  i = c - extra_offset
                  col_name = trim(s% profile_extra_name(i))
               else if (c > diffusion_D_offset) then
                  i = c - diffusion_D_offset
                  col_name = 'diffusion_D_' // trim(chem_isos% name(i))
               else if (c > diffusion_dX_offset) then
                  i = c - diffusion_dX_offset
                  col_name = 'diffusion_dX_' // trim(chem_isos% name(i))
               else if (c > log_concentration_offset) then
                  i = c - log_concentration_offset
                  col_name = 'log_concentration_' // trim(chem_isos% name(i))
               else if (c > log_g_rad_offset) then
                  i = c - log_g_rad_offset
                  col_name = 'log_g_rad_' // trim(chem_isos% name(i))
               else if (c > v_rad_offset) then
                  i = c - v_rad_offset
                  col_name = 'v_rad_' // trim(chem_isos% name(i))
               else if (c > extra_diffusion_factor_offset) then
                  i = c - extra_diffusion_factor_offset
                  col_name = 'extra_diffusion_factor_' // trim(chem_isos% name(i))
               else if (c > edv_offset) then
                  i = c - edv_offset
                  col_name = 'edv_' // trim(chem_isos% name(i))
               else if (c > typical_charge_offset) then
                  i = c - typical_charge_offset
                  col_name = 'typical_charge_' // trim(chem_isos% name(i))
               else if (c > ionization_offset) then
                  i = c - ionization_offset
                  col_name = 'ionization_' // trim(chem_isos% name(i))
               else if (c > xaprev_offset) then
                  i = c - xaprev_offset
                  col_name = 'xaprev_' // trim(chem_isos% name(i))
               else if (c > xadot_offset) then
                  i = c - xadot_offset
                  col_name = 'xadot_' // trim(chem_isos% name(i))
               else if (c > log_abundance_offset) then
                  i = c - log_abundance_offset
                  col_name = 'log_' // trim(chem_isos% name(i))
               else if (c > abundance_offset) then
                  i = c - abundance_offset
                  col_name = trim(chem_isos% name(i))
               else if (c > category_offset) then
                  i = c - category_offset
                  col_name = trim(category_name(i))
               else
                  col_name = trim(profile_column_name(c))
               end if
               if (write_flag) then
                  write(io, fmt=txt_fmt, advance='no') trim(col_name)
               else
                  names(col) = trim(col_name)
               end if
            else if (pass == 3) then
               call getval_for_profile(s, c, k, val, int_flag, int_val)
               if (write_flag) then
                  if (int_flag) then
                     write(io, fmt=int_fmt, advance='no') int_val
                  else
                     if (is_bad(val)) val = 1d99
                     write(io, fmt=dbl_fmt, advance='no') val
                  end if
               else
                  if (int_flag) then
                     vals(k,col) = dble(int_val)
                     is_int(col) = .true.
                  else
                     vals(k,col) = val
                     is_int(col) = .false.
                  end if
               end if
            end if
         end subroutine do_col


         subroutine do_abundance_col(pass, j, jj, k)
            integer, intent(in) :: pass, j, jj, k
            integer :: i, c, ii
            real(dp) :: val
            logical :: int_flag, log_abundance
            character (len=128) :: col_name
            logical, parameter :: dbg = .false.
            include 'formats'
            log_abundance = (s% profile_column_spec(j) == add_log_abundances)
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') col
            else if (pass == 2) then
               if (log_abundance) then
                  col_name = 'log_'
               else
                  col_name = ''
               end if
               col_name = trim(col_name) // trim(chem_isos% name(s% chem_id(jj)))
               if (write_flag) then
                  write(io, fmt=txt_fmt, advance='no') trim(col_name)
               else
                  names(col) = trim(col_name)
               end if
            else if (pass == 3) then
               val = s% xa(jj,k)
               if (log_abundance) val = safe_log10(val)
               if (write_flag) then
                  write(io, fmt=dbl_fmt, advance='no') val
               else
                  vals(k,col) = val
                  is_int(col) = .false.
               end if
            end if
         end subroutine do_abundance_col


      end subroutine do_profile_info


      subroutine do_save_profiles( &
            s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer, pointer, dimension(:) :: model_numbers, model_priorities, model_logs
         integer :: nz, max_num_mods, num_models, model_profile_number, k
         character (len=strlen) :: fname
         integer :: model_priority

         include 'formats'

         ierr = 0
         nz = s% nz

         if (.not. s% write_profiles_flag) return
         if (.not. s% v_flag) s% v(1:nz) = 0
         if (.not. s% u_flag) s% u(1:nz) = 0
         if (.not. s% rotation_flag) s% omega(1:nz) = 0

         max_num_mods = s% max_num_profile_models
         if (max_num_mods < 0) max_num_mods = s% model_number
         model_priority = s% save_profiles_model_priority

         allocate(model_numbers(max_num_mods), model_priorities(max_num_mods), &
            model_logs(max_num_mods), stat=ierr)
         if (ierr /= 0) return

         write(fname, '(3a)') trim(s% log_directory), '/', trim(s% profiles_index_name)

         call read_profiles_info( &
            fname, max_num_mods, num_models, model_numbers, model_priorities, model_logs)

         call make_room_for_profile_info( &
            s% model_number, max_num_mods, num_models, model_numbers, model_priorities, model_logs, ierr)
         if (ierr /= 0) then
            call dealloc; return
         end if

         call pick_model_profile_number( &
            max_num_mods, num_models, model_logs, model_profile_number, ierr)
         if (ierr /= 0) then
            call dealloc; return
         end if

         call get_model_profile_filename(s, model_profile_number)

         ! add the new model to the list at the end
         num_models = num_models+1
         model_numbers(num_models) = s% model_number
         model_priorities(num_models) = model_priority
         model_logs(num_models) = model_profile_number

         s% save_profiles_model_priority = delta_priority ! reset it to the default value

         ! write the profiles before adding them to the list
         ! so if user interrupts during write, the index is still okay.
         call write_profile_info(s, s% model_profile_filename, ierr)
         if (ierr /= 0) then
            call dealloc; return
         end if

         call write_profiles_list( &
            fname, num_models, model_numbers, model_priorities, model_logs, ierr)
         if (ierr /= 0) then
            call dealloc; return
         end if

         call dealloc


         contains

         subroutine dealloc
            deallocate(model_numbers, model_priorities, model_logs)
         end subroutine dealloc

      end subroutine do_save_profiles


      subroutine write_profiles_list( &
            fname, num_models, model_numbers, model_priorities, model_logs, ierr)
         character (len=*), intent(in) :: fname
         integer, intent(in) :: num_models
         integer, pointer, dimension(:) :: model_numbers, model_priorities, model_logs
         integer, intent(out) :: ierr
         integer :: iounit, i
         ierr = 0
         ! write the new list
         open(newunit=iounit, file=trim(fname), action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open ' // trim(fname)
         else
            if (num_models == 1) then
               write(iounit, *) num_models, &
                  'model.    lines hold model number, priority, and profile number.'
            else
               write(iounit, *) num_models, &
                  'models.    lines hold model number, priority, and profile number.'
            end if
            do i=1, num_models
               write(iounit, *) model_numbers(i), model_priorities(i), model_logs(i)
            end do
            close(iounit)
         end if
      end subroutine write_profiles_list


      subroutine pick_model_profile_number( &
            max_num_mods, num_models, model_logs, model_profile_number, ierr)
         integer, intent(in) :: max_num_mods
         integer, intent(inout) :: num_models
         integer, pointer, dimension(:) :: model_logs
         integer, intent(out) :: model_profile_number, ierr
         logical :: in_use(max_num_mods)
         integer :: i
         ! pick log number for the new model
         ierr = 0
         in_use = .false.
         do i=1, num_models
            in_use(model_logs(i)) = .true.
         end do
         model_profile_number = 0
         do i=1, max_num_mods
            if (.not. in_use(i)) then
               model_profile_number = i; exit
            end if
         end do
         if (model_profile_number == 0) then
            write(*, *) 'model_profile_number == 0, cannot happen?'
            ierr = -1
            return
         end if
      end subroutine pick_model_profile_number


      subroutine make_room_for_profile_info( &
            model_number, max_num_mods, num_models, model_numbers, model_priorities, model_logs, ierr)
         integer, intent(in) :: model_number, max_num_mods
         integer, intent(inout) :: num_models
         integer, pointer, dimension(:) :: model_numbers, model_priorities, model_logs
         integer, intent(out) :: ierr
         integer :: i, j, nm
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         ! delete models with model number greater or equal to current model number
         nm = num_models; j = 0
         do i=1, nm
            if (model_numbers(i) < model_number .and. model_logs(i) <= max_num_mods) then
               ! keep this one
               j = j+1
               if (j < i) then
                  model_numbers(j) = model_numbers(i)
                  model_priorities(j) = model_priorities(i)
                  model_logs(j) = model_logs(i)
               end if
            end if
         end do
         num_models = j
         if (num_models == max_num_mods) then ! pick one to delete
            j = 1
            do i=2, num_models
               if (dbg) then
                  write(*,3) 'model_priorities(i)', i, model_priorities(i)
                  write(*,3) 'model_priorities(j)', j, model_priorities(j)
                  write(*,3) 'model_numbers(i)', i, model_numbers(i)
                  write(*,3) 'model_numbers(j)', j, model_numbers(j)
                  write(*,*) 'model_priorities(i) < model_priorities(j)', model_priorities(i) < model_priorities(j)
                  write(*,*) 'model_numbers(i) < model_numbers(j)', model_numbers(i) < model_numbers(j)
               end if
               if (model_priorities(i) < model_priorities(j)) then
                  if (dbg) write(*,3) '1 change j'
                  j = i
               else if (model_priorities(i) == model_priorities(j) .and. &
                        model_numbers(i) < model_numbers(j)) then
                  if (dbg) write(*,3) '2 change j'
                  j = i
               end if
               if (dbg) write(*,3) 'new j', j
               if (dbg) write(*,*)
            end do
            ! delete j
            if (dbg) write(*,*) 'delete j', j
            do i=j+1, num_models
               model_numbers(i-1) = model_numbers(i)
               model_priorities(i-1) = model_priorities(i)
               model_logs(i-1) = model_logs(i)
            end do
            num_models = num_models-1
         end if
      end subroutine make_room_for_profile_info


      subroutine read_profiles_info( &
            fname, max_num_mods, num_models, model_numbers, model_priorities, model_logs)
         character (len=*), intent(in) :: fname
         integer, intent(in) :: max_num_mods
         integer, intent(out) :: num_models
         integer, pointer, dimension(:) :: model_numbers, model_priorities, model_logs
         integer :: iounit, i, ierr
         num_models = 0
         ierr = 0
         open(newunit=iounit, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr == 0) then ! file exists
            read(iounit, *, iostat=ierr) num_models
            if (ierr == 0) then
               if (num_models > max_num_mods) num_models = max_num_mods
               do i=1, num_models
                  read(iounit, *, iostat=ierr) model_numbers(i), model_priorities(i), model_logs(i)
                  if (ierr /= 0) exit
               end do
            end if
            close(iounit)
            if (ierr /= 0) num_models = 0
         end if
      end subroutine read_profiles_info


      subroutine get_model_profile_filename(s, model_profile_number)
         ! sets s% model_profile_filename and s% model_controls_filename
         type (star_info), pointer :: s
         integer, intent(in) :: model_profile_number
         character (len=strlen) :: &
            profile_prefix, controls_prefix, model_prefix, num_str, fstring
         integer :: num_digits

         profile_prefix = trim(s% log_directory) // '/' // trim(s% profile_data_prefix)
         controls_prefix = trim(s% log_directory) // '/' // trim(s% controls_data_prefix)
         model_prefix = trim(s% log_directory) // '/' // trim(s% model_data_prefix)

         num_digits = 1 + log10(dble(max(1,model_profile_number)))
         write(fstring,'( "(a,i",i2.2,".",i2.2,",a)" )') num_digits, num_digits

         write(s% model_profile_filename, fmt=fstring) &
            trim(profile_prefix), model_profile_number, trim(s% profile_data_suffix)
         write(s% model_controls_filename, fmt=fstring) &
            trim(controls_prefix), model_profile_number, trim(s% controls_data_suffix)
         write(s% model_data_filename, fmt=fstring) &
            trim(model_prefix), model_profile_number, trim(s% model_data_suffix)

      end subroutine get_model_profile_filename

      end module profile

