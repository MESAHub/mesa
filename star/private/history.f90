! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module history

      use star_private_def
      use star_history_def
      use const_def
      use chem_def
      use history_specs
      use star_utils

      implicit none

      private
      public :: get_history_specs, get_history_values, get1_hist_value, &
         do_get_data_for_history_columns, write_history_info


      logical, parameter :: open_close_log = .true.



      contains


      subroutine do_get_data_for_history_columns( &
            s, &
            ierr)
         type (star_info), pointer :: s

         integer, intent(out) :: ierr
         logical, parameter :: write_flag = .false.
         call do_history_info( &
            s, &
            write_flag, ierr)
      end subroutine do_get_data_for_history_columns


      subroutine write_history_info(s, ierr)
         type (star_info), pointer :: s

         integer, intent(out) :: ierr
         logical, parameter :: write_flag = .true.
         call do_history_info( &
            s, &
            write_flag, ierr)
      end subroutine write_history_info


      subroutine do_history_info(s, write_flag, ierr)
         use utils_lib, only: integer_dict_create_hash, integer_dict_free
         use chem_def, only: category_name
         use math_lib, only: math_backend
         use rates_def, only: rates_reaction_id_max, i_rate
         type (star_info), pointer :: s

         logical, intent(in) :: write_flag
         integer, intent(out) :: ierr

         character (len=strlen) :: fname, dbl_fmt, int_fmt, txt_fmt
         integer :: numcols, io, i, nz, col, j, i0, &
            num_extra_cols, num_binary_cols, num_extra_binary_cols, num_extra_header_items, n
         integer, parameter :: num_epsnuc_out = 12
         real(dp) :: &
            epsnuc_out(num_epsnuc_out), csound_surf, v_surf, envelope_fraction_left
         integer :: mixing_regions, mix_relr_regions, burning_regions, burn_relr_regions
         integer, pointer :: mixing_type(:), mix_relr_type(:), burning_type(:), burn_relr_type(:)
         character (len=maxlen_history_column_name), pointer, dimension(:) :: &
            extra_col_names, binary_col_names, extra_binary_col_names, extra_header_item_names
         real(dp), pointer, dimension(:) :: &
            extra_col_vals, binary_col_vals, extra_binary_col_vals, extra_header_item_vals

         character (len=maxlen_history_column_name), pointer :: &
            names(:) ! (num_history_columns)
         real(dp), pointer :: vals(:) ! (num_history_columns)
         logical, pointer :: is_int(:) ! (num_history_columns)

         logical :: history_file_exists

         include 'formats'

         dbl_fmt = s% star_history_dbl_format
         int_fmt = s% star_history_int_format
         txt_fmt = s% star_history_txt_format

         ierr = 0

         if (.not. associated(s% history_column_spec)) then
            numcols = 0
         else
            numcols = size(s% history_column_spec, dim=1)
         end if
         num_extra_cols = s% how_many_extra_history_columns(s% id)
         if (s% include_binary_history_in_log_file) then
            num_binary_cols = s% how_many_binary_history_columns(s% binary_id)
            num_extra_binary_cols = s% how_many_extra_binary_history_columns(s% binary_id)
         else
            num_binary_cols = 0
            num_extra_binary_cols = 0
         end if
         n = numcols + num_extra_cols + num_binary_cols + num_extra_binary_cols
         if (n == 0) then
            write(*,*) 'WARNING: do not have any output specified for logs.'
            ierr = -1
            return
         end if

         if (s% number_of_history_columns < 0) then
            s% number_of_history_columns = n
         else if (s% number_of_history_columns /= n) then
            if (associated(s% history_values)) then
               deallocate(s% history_values)
               nullify(s% history_values)
            end if
            if (associated(s% history_names)) then
               deallocate(s% history_names)
               nullify(s% history_names)
            end if
            if (associated(s% history_value_is_integer)) then
               deallocate(s% history_value_is_integer)
               nullify(s% history_value_is_integer)
            end if
            if (associated(s% history_names_dict)) then
               call integer_dict_free(s% history_names_dict)
               nullify(s% history_names_dict)
            end if
            s% need_to_set_history_names_etc = .true.
            s% number_of_history_columns = n
         end if

         if (.not. associated(s% history_values)) then
            allocate(s% history_values(n))
         else if (size(s% history_values,dim=1) /= n) then
            ierr = -1
            write(*,2) 'bad size s% history_values', &
               size(s% history_values,dim=1), n
         end if
         vals => s% history_values

         if (.not. associated(s% history_names)) then
            allocate(s% history_names(n))
         else if (size(s% history_names,dim=1) /= n) then
            ierr = -1
            write(*,2) 'bad size s% history_names', &
               size(s% history_names,dim=1), n
         end if
         
         if (s% need_to_set_history_names_etc) then
            names => s% history_names
         else
            nullify(names)
         end if

         if (.not. associated(s% history_value_is_integer)) then
            allocate(s% history_value_is_integer(n))
         else if (size(s% history_value_is_integer,dim=1) /= n) then
            ierr = -1
            write(*,2) 'bad size s% history_value_is_integer', &
               size(s% history_value_is_integer,dim=1), n
         end if
         if (s% need_to_set_history_names_etc) then
            is_int => s% history_value_is_integer
         else
            nullify(is_int)
         end if

         nz = s% nz
         nullify(mixing_type)
         nullify(mix_relr_type)
         nullify(burning_type)
         nullify(burn_relr_type)
         nullify(extra_col_names)
         nullify(extra_col_vals)
         nullify(binary_col_names)
         nullify(binary_col_vals)
         nullify(extra_binary_col_names)
         nullify(extra_binary_col_vals)

         if (num_extra_cols > 0) then
            allocate( &
               extra_col_names(num_extra_cols), extra_col_vals(num_extra_cols), stat=ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            extra_col_names(1:num_extra_cols) = 'unknown'
            extra_col_vals(1:num_extra_cols)  = -1d99
            call s% data_for_extra_history_columns( &
               s% id, num_extra_cols, extra_col_names, extra_col_vals, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            do i=1,num_extra_cols
               if(trim(extra_col_names(i))=='unknown') then
                  write(*,*) "Warning empty history name for extra_history_column ",i
               end if
            end do
         end if

         if (num_binary_cols > 0) then
            allocate( &
               binary_col_names(num_binary_cols), binary_col_vals(num_binary_cols), stat=ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            call s% data_for_binary_history_columns( &
               s% binary_id, num_binary_cols, binary_col_names, binary_col_vals, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
         end if

         if (num_extra_binary_cols > 0) then
            allocate( &
               extra_binary_col_names(num_extra_binary_cols), extra_binary_col_vals(num_extra_binary_cols), stat=ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            extra_binary_col_names(1:num_extra_binary_cols) = 'unknown'
            extra_binary_col_vals(1:num_extra_binary_cols)  = -1d99
            call s% data_for_extra_binary_history_columns( &
               s% binary_id, num_extra_binary_cols, extra_binary_col_names, extra_binary_col_vals, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            do i=1,num_extra_binary_cols
               if(trim(extra_binary_col_names(i))=='unknown') then
                  write(*,*) "Warning empty history name for extra_binary_history_column ",i
               end if
            end do
         end if

         i0 = 1
         if (write_flag .and. (open_close_log .or. s% model_number == -100)) then
            fname = trim(s% log_directory) // '/' // trim(s% star_history_name)
            inquire(file=trim(fname), exist=history_file_exists)
            if ((.not. history_file_exists) .or. &
                  s% doing_first_model_of_run .or. s% need_to_set_history_names_etc) then
               ierr = 0
               if (len_trim(s% star_history_header_name) > 0) then
                  fname = trim(s% log_directory) // '/' // trim(s% star_history_header_name)
               end if
               open(newunit=io, file=trim(fname), action='write', iostat=ierr)
            else
               i0 = 3
               open(newunit=io, file=trim(fname), action='write', position='append', iostat=ierr)
            end if
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(fname)
               call dealloc
               return
            end if
         end if

         csound_surf = eval_csound(s,1,ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if

         if (s% u_flag) then
            v_surf = s% u(1)
         else if (s% v_flag) then
            v_surf = s% v(1)
         else
            v_surf = s% r(1)*s% dlnR_dt(1)
         end if

         if (s% initial_mass > s% he_core_mass) then
            envelope_fraction_left = &
               (s% star_mass - s% he_core_mass)/(s% initial_mass - s% he_core_mass)
         else
            envelope_fraction_left = 1
         end if

         if (write_flag .and. i0 == 1) then ! write values at start of log

            num_extra_header_items = s% how_many_extra_history_header_items(s% id)
            if (num_extra_header_items > 0) then
               allocate( &
                  extra_header_item_names(num_extra_header_items), &
                  extra_header_item_vals(num_extra_header_items), stat=ierr)
               if (ierr /= 0) then
                  call dealloc
                  return
               end if
               extra_header_item_names(1:num_extra_header_items) = 'unknown'
               extra_header_item_vals(1:num_extra_header_items)  = -1d99
               call s% data_for_extra_history_header_items( &
                  s% id, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
               if (ierr /= 0) then
                  call dealloc
                  return
               end if
               do i=1,num_extra_header_items
                  if(trim(extra_header_item_names(i))=='unknown') then
                     write(*,*) "Warning empty history name for extra_history_header ",i
                  end if
               end do
            end if

            do i=1,3
               col = 0
               call write_string(io, col, i, 'version_number', version_number)
               call write_string(io, col, i, 'compiler', compiler_name)
               call write_string(io, col, i, 'build', compiler_version_name)
               call write_string(io, col, i, 'MESA_SDK_version', mesasdk_version_name)
               call write_string(io, col, i, 'math_backend',math_backend)
               call write_string(io, col, i, 'date', date)
               !call write_val(io, col, i, 'initial_mass', s% initial_mass)
               !call write_val(io, col, i, 'initial_z', s% initial_z)
               call write_val(io, col, i, 'burn_min1', s% burn_min1)
               call write_val(io, col, i, 'burn_min2', s% burn_min2)
               
               do j=1,num_extra_header_items
                 call write_val(io, col, i, &
                    extra_header_item_names(j), extra_header_item_vals(j))
               end do

               write(io,*)

            end do

            write(io,*)

            if (num_extra_header_items > 0) &
               deallocate(extra_header_item_names, extra_header_item_vals)

         end if
         
         ! count the number of declared log regions in pass2. (= param in history_columns.list)
         ! write values for the current regions in pass3.
         burning_regions = 0
         mixing_regions = 0
         mix_relr_regions = 0
         burn_relr_regions = 0

         do i=i0,3 ! add a row to the log

            col = 0
            if (i==3) then

               if (write_flag .and. i0 == 1 .and. len_trim(s% star_history_header_name) > 0) then
                  close(io)
                  fname = trim(s% log_directory) // '/' // trim(s% star_history_name)
                  open(newunit=io, file=trim(fname), action='write', status='replace', iostat=ierr)
                  if (ierr /= 0) then
                     call dealloc
                     return
                  end if
               end if

               epsnuc_out(1:4) = s% burn_zone_mass(1:4,1)
               epsnuc_out(5:8) = s% burn_zone_mass(1:4,2)
               epsnuc_out(9:12) = s% burn_zone_mass(1:4,3)

               mixing_regions = count_output_mix_regions(mixing_offset)
               if (mixing_regions > 0) then
                  allocate(mixing_type(nz), stat=ierr)
                  if (ierr /= 0) exit
                  call set_mix_types(mixing_type, mixing_regions)
                  call prune_weak_mixing_regions(mixing_type, mixing_regions)
               end if

               mix_relr_regions = count_output_mix_regions(mix_relr_offset)
               if (mix_relr_regions > 0) then
                  allocate(mix_relr_type(nz), stat=ierr)
                  if (ierr /= 0) exit
                  call set_mix_types(mix_relr_type, mix_relr_regions)
                  call prune_weak_mixing_regions(mix_relr_type, mix_relr_regions)
               end if

               burning_regions = count_output_burn_regions(burning_offset)
               if (burning_regions > 0) then
                  allocate(burning_type(nz), stat=ierr)
                  if (ierr /= 0) exit
                  call set_burn_types(burning_type, burning_offset)
               end if
               
               burn_relr_regions = count_output_burn_regions(burn_relr_offset)
               if (burn_relr_regions > 0) then
                  allocate(burn_relr_type(nz), stat=ierr)
                  if (ierr /= 0) exit
                  call set_burn_types(burn_relr_type, burn_relr_offset)
               end if

            end if

            do j=1,numcols
               call do_col(i,j)
            end do
            do j=1,num_extra_cols
               call do_extra_col(i,j,numcols)
            end do
            do j=1,num_binary_cols
               call do_binary_col(i,j,numcols+num_extra_cols)
            end do
            do j=1,num_extra_binary_cols
               call do_extra_binary_col(i,j,numcols+num_extra_cols+num_binary_cols)
            end do

            if (write_flag) write(io,*)

         end do

         if (open_close_log .and. write_flag) close(io)

         call dealloc

         s% model_number_of_history_values = s% model_number

         if (s% need_to_set_history_names_etc) then
            call integer_dict_create_hash(s% history_names_dict, ierr)
         end if

         s% need_to_set_history_names_etc = .false.


         contains


         subroutine dealloc
            if (associated(mixing_type)) deallocate(mixing_type)
            if (associated(mix_relr_type)) deallocate(mix_relr_type)
            if (associated(burning_type)) deallocate(burning_type)
            if (associated(burn_relr_type)) deallocate(burn_relr_type)
            if (associated(extra_col_names)) deallocate(extra_col_names)
            if (associated(extra_col_vals)) deallocate(extra_col_vals)
            if (associated(binary_col_names)) deallocate(binary_col_names)
            if (associated(binary_col_vals)) deallocate(binary_col_vals)
            if (associated(extra_binary_col_names)) deallocate(extra_binary_col_names)
            if (associated(extra_binary_col_vals)) deallocate(extra_binary_col_vals)
            
            nullify(mixing_type)
            nullify(mix_relr_type)
            nullify(burning_type)
            nullify(burn_relr_type)
            nullify(extra_col_names)
            nullify(extra_col_vals)
            nullify(binary_col_names)
            nullify(binary_col_vals)
            nullify(extra_binary_col_names)
            nullify(extra_binary_col_vals)

            
         end subroutine dealloc


         subroutine do_extra_col(pass, j, col_offset)
            integer, intent(in) :: pass, j, col_offset
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') j + col_offset
            else if (pass == 2) then
               call do_name(j + col_offset, extra_col_names(j))
            else if (pass == 3) then
               call do_val(j + col_offset, extra_col_vals(j))
            end if
         end subroutine do_extra_col


         subroutine do_binary_col(pass, j, col_offset)
            integer, intent(in) :: pass, j, col_offset
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') j + col_offset
            else if (pass == 2) then
               call do_name(j + col_offset, binary_col_names(j))
            else if (pass == 3) then
               call do_val(j + col_offset, binary_col_vals(j))
            end if
         end subroutine do_binary_col


         subroutine do_extra_binary_col(pass, j, col_offset)
            integer, intent(in) :: pass, j, col_offset
            if (pass == 1) then
               if (write_flag) write(io, fmt=int_fmt, advance='no') j + col_offset
            else if (pass == 2) then
               call do_name(j + col_offset, extra_binary_col_names(j))
            else if (pass == 3) then
               call do_val(j + col_offset, extra_binary_col_vals(j))
            end if
         end subroutine do_extra_binary_col


         subroutine do_name(j, col_name)
            use utils_lib, only: integer_dict_define
            integer, intent(in) :: j
            integer :: ierr
            character (len=*), intent(in) :: col_name
            if (write_flag) write(io, fmt=txt_fmt, advance='no') trim(col_name)
            if (associated(names)) names(j) = trim(col_name)
            if (s% need_to_set_history_names_etc) then
               call integer_dict_define(s% history_names_dict, col_name, j, ierr)
               if (ierr /= 0) write(*,*) 'failed in integer_dict_define ' // trim(col_name)
            end if
         end subroutine do_name


         integer function count_output_burn_regions(b_regions)
            integer, intent(in) :: b_regions
            integer :: j, cnt, c, i, ii
            cnt = 0
            do j=1,numcols
               c = s% history_column_spec(j)
               if (c > b_regions .and. c <= b_regions + idel) then
                  i = c - b_regions
                  ii = (i+1)/2
                  if (ii > cnt) cnt = ii
               end if
            end do
            count_output_burn_regions = cnt
         end function count_output_burn_regions


         subroutine set_burn_types(b_type, b_regions)
            integer :: b_type(:), b_regions
            integer :: cnt, min_ktop, min_kbot, imax, i, prev_cnt, k
            real(dp) :: val
            logical, parameter :: dbg = .false.
            do k=1,nz
               val = s% eps_nuc(k) - s% non_nuc_neu(k)
               b_type(k) = int(sign(1d0,val)*max(0d0,1d0+safe_log10(abs(val))))
            end do
            ! remove smallest regions until <= b_regions remain
            imax = nz
            do i=1,imax
               call count_regions(b_type, cnt, min_ktop, min_kbot)
               if (dbg) write(*,*) 'count_regions', cnt
               if (cnt <= b_regions) return
               if (i > 1 .and. cnt >= prev_cnt) then
                  write(*,*) 'bug in set_burn_types: cnt, prev_cnt', cnt, prev_cnt
                  if (dbg) stop 'debug: set_burn_types'
                  return
               end if
               prev_cnt = cnt
               if (dbg) write(*,*) 'remove_region', min_ktop, min_kbot, cnt
               call remove_region(b_type, min_ktop, min_kbot)
            end do
            if (dbg) stop 'debug: set_burn_types'
         end subroutine set_burn_types


         integer function count_output_mix_regions(mx_offset)
            integer, intent(in) :: mx_offset
            integer :: j, cnt, c, i, ii
            cnt = 0
            do j=1,numcols
               c = s% history_column_spec(j)
               if (c > mx_offset .and. c < mx_offset+idel) then
                  i = c - mx_offset
                  ii = (i+1)/2
                  if (ii > cnt) cnt = ii
               end if
            end do
            count_output_mix_regions = cnt
         end function count_output_mix_regions


         subroutine prune_weak_mixing_regions(mx_type, mx_regions)
            integer :: mx_type(:), mx_regions
            real(dp) :: D_max_in_region, D_cutoff
            integer :: min_ktop, min_kbot
            integer :: i
            logical, parameter :: dbg = .false.
            include 'formats'
            D_cutoff = s% mixing_D_limit_for_log
            do i = 1, 1000
               call find_weakest_mixing_region( &
                  mx_type, mx_regions, D_max_in_region, min_ktop, min_kbot)
               if (D_max_in_region > D_cutoff) then
                  if (dbg) write(*,3) 'done', min_ktop, min_kbot, D_max_in_region, D_cutoff
                  return
               end if
               if (mx_type(min_ktop) == no_mixing) then
                  if (dbg) write(*,3) 'no mixing regions left'
                  return
               end if
               if (dbg) write(*,3) 'prune_weak_mixing_regions', min_ktop, min_kbot, D_max_in_region
               if (dbg) write(*,3) 'i model mass', &
                  i, s% model_number, s% m(min_kbot)/Msun, s% m(min_ktop)/Msun
               if (dbg) write(*,*)
               mx_type(min_ktop:min_kbot) = no_mixing
            end do
         end subroutine prune_weak_mixing_regions


         subroutine find_weakest_mixing_region( &
               mx_type, mx_regions, D_max_in_region, min_ktop, min_kbot)
            integer :: mx_type(:), mx_regions
            real(dp), intent(out) :: D_max_in_region
            integer, intent(out) :: min_ktop, min_kbot

            integer :: k, kbot, cur_type
            real(dp) :: D_max
            logical, parameter :: dbg = .false.
            include 'formats'
            min_ktop = 1
            min_kbot = nz
            D_max_in_region = 1d99
            kbot = nz
            cur_type = mx_type(nz)
            do k=nz-1,1,-1
               if (cur_type == mx_type(k)) cycle
               ! k is bottom of new region; k+1 is top of region cnt
               if (cur_type /= no_mixing &
                     .and. cur_type /= thermohaline_mixing &
                     .and. cur_type /= semiconvective_mixing) then
                  D_max = maxval(s% D_mix_non_rotation(k+1:kbot))
                  if (D_max < D_max_in_region) then
                     D_max_in_region = D_max; min_ktop = k+1; min_kbot = kbot
                  end if
               end if
               kbot = k
               cur_type = mx_type(k)
            end do
            if (cur_type == no_mixing) then
               D_max = maxval(s% D_mix_non_rotation(1:kbot))
               if (D_max < D_max_in_region) then
                  D_max_in_region = D_max; min_ktop = 1; min_kbot = kbot
               end if
            end if

         end subroutine find_weakest_mixing_region


         subroutine set_mix_types(mx_type, mx_regions)
            integer :: mx_type(:), mx_regions

            integer :: cnt, min_ktop, min_kbot, imax, i, prev_cnt, k
            logical, parameter :: dbg = .false.

            do k=1,nz
               mx_type(k) = s% mixing_type(k)
            end do
            ! remove smallest regions until <= mixing_regions remain
            imax = nz
            do i=1,imax
               call count_regions(mx_type, cnt, min_ktop, min_kbot)
               if (dbg) write(*,*) 'count_regions', cnt
               if (cnt <= mx_regions) exit
               if (i > 1 .and. cnt >= prev_cnt) then
                  write(*,*) 'bug in set_mix_types: cnt, prev_cnt', cnt, prev_cnt
                  if (dbg) stop 'set_mix_types'
                  return
               end if
               prev_cnt = cnt
               if (dbg) write(*,*) 'remove_region', min_ktop, min_kbot, cnt
               call remove_region(mx_type, min_ktop, min_kbot)
            end do
            if (dbg) call show_regions(mx_type)
         end subroutine set_mix_types


         subroutine count_regions(mx_type, cnt, min_ktop, min_kbot)
            integer, intent(in) :: mx_type(:)
            integer, intent(out) :: cnt, min_ktop, min_kbot
            integer :: k, kbot, cur_type
            real(dp) :: prev_qtop, qtop, dq, min_dq
            logical, parameter :: dbg = .false.
            include 'formats'
            cnt = 1
            min_ktop = 1
            min_kbot = nz
            prev_qtop = 0
            min_dq = 1
            kbot = nz
            cur_type = mx_type(nz)
            do k=nz-1,1,-1
               if (cur_type == mx_type(k)) cycle
               ! k is bottom of new region; k+1 is top of region cnt
               qtop = s% q(k+1)
               dq = qtop - prev_qtop
               if (dq < min_dq) then
                  min_dq = dq; min_ktop = k+1; min_kbot = kbot
                  if (dbg) write(*,3) &
                     'loop min_ktop min_kbot min_dq', min_ktop, min_kbot, min_dq
               end if
               kbot = k
               cnt = cnt+1
               cur_type = mx_type(k)
               prev_qtop = qtop
            end do
            qtop = 1
            dq = qtop - prev_qtop
            if (dq < min_dq) then
               min_dq = dq; min_ktop = 1; min_kbot = kbot
               if (dbg) write(*,3) &
                  'final min_ktop min_kbot min_dq', min_ktop, min_kbot, min_dq
            end if
         end subroutine count_regions


         subroutine show_regions(mx_type)
            integer, intent(in) :: mx_type(:)
            integer :: k, kbot, cur_type
            include 'formats'
            kbot = nz
            cur_type = mx_type(nz)
            do k=nz-1,1,-1
               if (cur_type == mx_type(k)) cycle
               ! k is bottom of new region; k+1 is top of region cnt
               write(*,4) 'mix region', cur_type, k+1, kbot, s% m(k+1)/Msun, s% m(kbot)/Msun
               kbot = k
               cur_type = mx_type(k)
            end do
            write(*,4) 'mix region', cur_type, 1, kbot, s% m(1)/Msun, s% m(kbot)/Msun
         end subroutine show_regions


         subroutine remove_region(mx_type, min_ktop, min_kbot)
            integer, intent(inout) :: mx_type(:)
            integer, intent(in) :: min_ktop, min_kbot
            integer :: new_type
            logical, parameter :: dbg = .false.
            include 'formats'
            if (dbg) then
               write(*,2) 'q top', min_ktop, s% q(min_ktop)
               write(*,2) 'q bot', min_kbot, s% q(min_kbot)
               write(*,2) 'dq', min_kbot, s% q(min_ktop) - s% q(min_kbot)
            end if
            if (min_ktop > 1) then ! merge with above
               new_type = mx_type(min_ktop-1)
            else if (min_kbot < nz) then ! merge with below
               new_type = mx_type(min_kbot+1)
            else
               write(*,*) 'confusion in args for remove_region', min_ktop, min_kbot
               return
            end if
            mx_type(min_ktop:min_kbot) = new_type
            if (dbg) write(*,2) 'new type', min_kbot, new_type
         end subroutine remove_region


         integer function region_top(mx_type,j)
            integer, intent(in) :: mx_type(:), j
            integer :: k, cnt, cur_type
            cnt = 1
            cur_type = mx_type(nz)
            do k=nz-1,1,-1
               if (cur_type == mx_type(k)) cycle
               ! k is start of new region; k+1 is top of region cnt
               if (cnt == j) then
                  region_top = k+1
                  return
               end if
               cnt = cnt+1
               cur_type = mx_type(k)
            end do
            if (cnt == j) then
               region_top = 1
            else
               region_top = -1
            end if
         end function region_top


         real(dp) function interpolate_burn_bdy_q(k) result(val)
            use num_lib, only: find0
            integer, intent(in) :: k
            integer :: bv, bv0, bv1
            real(dp) :: eps, eps0, eps1, d0, d1, q0, q1
            include 'formats'
            if (k <= 0) then
               val = -1; return
            end if
            val = s% q(k)
            if (k == 1) return
            eps1 = s% eps_nuc(k) - s% eps_nuc_neu_total(k) - s% non_nuc_neu(k)
            bv1 = sign(1d0,eps1)*log10(max(1d0,abs(eps1)))
            eps0 = s% eps_nuc(k-1) - s% eps_nuc_neu_total(k-1) - s% non_nuc_neu(k-1)
            bv0 = sign(1d0,eps0)*log10(max(1d0,abs(eps0)))
            bv = max(bv0,bv1)
            eps = pow(10d0,bv)
            d0 = eps0 - eps
            d1 = eps1 - eps
            if (d0*d1 >= 0) return
            q1 = s% q(k) - s% dq(k)*0.5d0
            q0 = s% q(k-1) - s% dq(k-1)*0.5d0
            val = find0(q1,d1,q0,d0)
         end function interpolate_burn_bdy_q
         

         real(dp) function interpolate_burn_bdy_r(k) result(val)
            use num_lib, only: find0
            integer, intent(in) :: k
            integer :: bv, bv0, bv1
            real(dp) :: eps, eps0, eps1, d0, d1, q0, q1
            include 'formats'
            if (k <= 0) then
               val = -1; return
            end if
            val = s% r(k)/s%r(1)
            if (k == 1) return
            eps1 = s% eps_nuc(k) - s% eps_nuc_neu_total(k) - s% non_nuc_neu(k)
            bv1 = sign(1d0,eps1)*log10(max(1d0,abs(eps1)))
            eps0 = s% eps_nuc(k-1) - s% eps_nuc_neu_total(k-1) - s% non_nuc_neu(k-1)
            bv0 = sign(1d0,eps0)*log10(max(1d0,abs(eps0)))
            bv = max(bv0,bv1)
            eps = pow(10d0,bv)
            d0 = eps0 - eps
            d1 = eps1 - eps
            if (d0*d1 >= 0) return
            val = find0(s%rmid(k),d1,s%rmid(k-1),d0)/s%r(1)
         end function interpolate_burn_bdy_r

         subroutine do_col(pass, j)
            integer, intent(in) :: pass, j
            if (pass == 1) then
               call do_col_pass1
            else if (pass == 2) then
               call do_col_pass2(j)
            else if (pass == 3) then
               call do_col_pass3(s% history_column_spec(j))
            end if
         end subroutine do_col


         subroutine do_col_pass1 ! write the column number
            col = col+1
            if (write_flag) write(io, fmt=int_fmt, advance='no') col
         end subroutine do_col_pass1


         ! The order of if statements matter, they should be in reverse order
         ! to the order in history_specs
         subroutine do_col_pass2(j) ! get the column name
            use colors_lib, only : get_bc_name_by_id
            integer, intent(in) :: j
            character (len=100) :: col_name
            character (len=10) :: str
            integer :: c, i, ii
            c = s% history_column_spec(j)
            if (c > burn_relr_offset) then
               i = c - burn_relr_offset
               ii = (i+1)/2
               if (ii > burn_relr_regions) burn_relr_regions = ii ! count the regions in pass2
               if (ii < 10) then
                  write(str,'(i1)') ii
               else if (ii < 100) then
                  write(str,'(i2)') ii
               else
                  write(str,'(i3)') ii
               end if
               if (mod(i,2)==1) then ! burning type
                  col_name = 'burn_relr_type_' // trim(str)
               else ! location of top
                  col_name = 'burn_relr_top_' // trim(str)
               end if
            else if (c > burning_offset) then
               i = c - burning_offset
               ii = (i+1)/2
               if (ii > burning_regions) burning_regions = ii ! count the regions in pass2
               if (ii < 10) then
                  write(str,'(i1)') ii
               else if (ii < 100) then
                  write(str,'(i2)') ii
               else
                  write(str,'(i3)') ii
               end if
               if (mod(i,2)==1) then ! burning type
                  col_name = 'burn_type_' // trim(str)
               else ! location of top
                  col_name = 'burn_qtop_' // trim(str)
               end if
            else if (c > mix_relr_offset) then
               i = c - mix_relr_offset
               ii = (i+1)/2
               if (ii > mix_relr_regions) mix_relr_regions = ii ! count the regions in pass2
               if (ii < 10) then
                  write(str,'(i1)') ii
               else if (ii < 100) then
                  write(str,'(i2)') ii
               else
                  write(str,'(i3)') ii
               end if
               if (mod(i,2)==1) then
                  col_name = 'mix_relr_type_' // trim(str)
               else ! location of top
                  col_name = 'mix_relr_top_' // trim(str)
               end if
            else if (c > mixing_offset) then
               i = c - mixing_offset
               ii = (i+1)/2
               if (ii > mixing_regions) mixing_regions = ii ! count the regions in pass2
               if (ii < 10) then
                  write(str,'(i1)') ii
               else if (ii < 100) then
                  write(str,'(i2)') ii
               else
                  write(str,'(i3)') ii
               end if
               if (mod(i,2)==1) then ! mixing type
                  col_name = 'mix_type_' // trim(str)
               else ! location of top
                  col_name = 'mix_qtop_' // trim(str)
               end if
            else if (c > log_lum_band_offset) then
               i = c - log_lum_band_offset
               col_name = 'log_lum_band_' // trim(get_bc_name_by_id(i,ierr))
            else if (c > lum_band_offset) then
               i = c - lum_band_offset
               col_name = 'lum_band_' // trim(get_bc_name_by_id(i,ierr))
            else if (c > abs_mag_offset) then
               i = c - abs_mag_offset
               col_name = 'abs_mag_' // trim(get_bc_name_by_id(i,ierr))
            else if (c > bc_offset) then
               i = c - bc_offset
               col_name = 'bc_' // trim(get_bc_name_by_id(i,ierr))
            else if (c > c_log_eps_burn_offset) then
               i = c - c_log_eps_burn_offset
               col_name = 'c_log_eps_burn_' // trim(category_name(i))
            else if (c > max_eps_nuc_offset) then
               i = c - max_eps_nuc_offset
               col_name = 'max_eps_nuc_log_' // trim(chem_isos% name(i))
            else if (c > cz_top_max_offset) then
               i = c - cz_top_max_offset
               col_name = 'cz_top_log_' // trim(chem_isos% name(i))
            else if (c > cz_max_offset) then
               i = c - cz_max_offset
               col_name = 'cz_log_' // trim(chem_isos% name(i))
            else if (c > log_surface_xa_offset) then
               i = c - log_surface_xa_offset
               col_name = 'log_surface_' // trim(chem_isos% name(i))
            else if (c > log_center_xa_offset) then
               i = c - log_center_xa_offset
               col_name = 'log_center_' // trim(chem_isos% name(i))
            else if (c > log_average_xa_offset) then
               i = c - log_average_xa_offset
               col_name = 'log_average_' // trim(chem_isos% name(i))
            else if (c > log_total_mass_offset) then
               i = c - log_total_mass_offset
               col_name = 'log_total_mass_' // trim(chem_isos% name(i))
            else if (c > total_mass_offset) then
               i = c - total_mass_offset
               col_name = 'total_mass_' // trim(chem_isos% name(i))
            else if (c > category_offset) then
               i = c - category_offset
               col_name = category_name(i)
            else if (c > average_xa_offset) then
               i = c - average_xa_offset
               col_name = 'average_' // trim(chem_isos% name(i))
            else if (c > surface_xa_offset) then
               i = c - surface_xa_offset
               col_name = 'surface_' // trim(chem_isos% name(i))
            else if (c > center_xa_offset) then
               i = c - center_xa_offset
               col_name = 'center_' // trim(chem_isos% name(i))
            else
               col_name = trim(history_column_name(c))
            end if
            call do_name(j, col_name)
         end subroutine do_col_pass2


         subroutine do_col_pass3(c) ! get the column value
            use rates_def
            integer, intent(in) :: c
            integer :: i, ii, k, int_val
            logical :: is_int_val
            real(dp) :: val, val1, Ledd, power_photo, frac
            int_val = 0; val = 0; is_int_val = .false.
            
            if (c > burn_relr_offset) then
               i = c - burn_relr_offset
               ii = (i+1)/2
               k = region_top(burn_relr_type, ii)
               if (mod(i,2)==1) then ! burning type
                  is_int_val = .true.
                  if (k > 0) then
                     int_val = burn_relr_type(k)
                  else
                     int_val = -9999
                  end if
               else ! location of top
                  val = interpolate_burn_bdy_r(k)
               end if
            else if (c > burning_offset) then
               i = c - burning_offset
               ii = (i+1)/2
               k = region_top(burning_type, ii)
               if (mod(i,2)==1) then ! burning type
                  is_int_val = .true.
                  if (k > 0) then
                     int_val = burning_type(k)
                  else
                     int_val = -9999
                  end if
               else ! location of top
                  val = interpolate_burn_bdy_q(k)
               end if
            else if (c > mix_relr_offset) then
               i = c - mix_relr_offset
               ii = (i+1)/2
               k = region_top(mix_relr_type, ii)
               if (mod(i,2)==1) then ! mixing type
                  is_int_val = .true.
                  if (k > 0) then
                     int_val = mix_relr_type(k)
                  else
                     int_val = -1
                  end if
               else ! r/rstar location of boundary
                  if (k <= 1) then
                     val = 1d0
                  else
                     frac = s% cz_bdy_dq(k-1)/s% dq(k-1)
                     val = (1d0 - frac)*pow3(s% r(k-1)) + frac*pow3(s% r(k))
                     val = pow(val,one_third)/s% r(1)
                  end if
               end if
            else if (c > mixing_offset) then
               i = c - mixing_offset
               ii = (i+1)/2
               k = region_top(mixing_type, ii)
               if (mod(i,2)==1) then ! mixing type
                  is_int_val = .true.
                  if (k > 0) then
                     int_val = mixing_type(k)
                  else
                     int_val = -1
                  end if
               else ! q location of boundary
                  if (k <= 1) then
                     val = 1d0
                  else
                     val = s% q(k-1) - s% cz_bdy_dq(k-1)
                  end if
               end if
            else
               call history_getval( &
                  s, c, val, int_val, is_int_val, &
                  nz, v_surf, csound_surf, envelope_fraction_left, epsnuc_out, ierr)
               if (ierr /= 0) then
                  write(*,*) 'missing log info for ' // trim(history_column_name(c)), j, k
                  int_val = -99999999
                  is_int_val = .true.
                  ierr = -1
               end if
            end if
            if (is_int_val) then
               call do_int_val(j,int_val)
            else
               call do_val(j,val)
            end if
         end subroutine do_col_pass3


         subroutine do_val(j, val)
            integer, intent(in) :: j
            real(dp), intent(in) :: val
            if (write_flag) then
               if (is_bad_num(val)) then
                  write(io, fmt=dbl_fmt, advance='no') -1d99
               else
                  write(io, fmt=dbl_fmt, advance='no') val
               end if
            end if
            if (associated(vals)) vals(j) = val
            if (associated(is_int)) is_int(j) = .false.
         end subroutine do_val


         subroutine do_int_val(j, val)
            integer, intent(in) :: j
            integer, intent(in) :: val
            if (write_flag) write(io, fmt=int_fmt, advance='no') val
            if (associated(vals)) vals(j) = dble(val)
            if (associated(is_int)) is_int(j) = .true.
         end subroutine do_int_val
         
         subroutine write_string(io, col, pass, name, val) !for header items only
           integer, intent(in) :: io, pass
           integer, intent(inout) :: col
           character(len=*), intent(in) :: name, val
           character(len=strlen) :: my_val
           
           my_val = '"'//trim(val)//'"'
           if (pass == 1) then
              col = col+1
              write(io, fmt=int_fmt, advance='no') col
           else if (pass == 2) then
              write(io, fmt=txt_fmt, advance='no') trim(name)
           else if (pass == 3) then
              write(io, fmt=txt_fmt, advance='no') trim(my_val)
           end if
         end subroutine write_string

      
         subroutine write_integer(io, col, pass, name, val) ! for header items only
           integer, intent(in) :: io, pass
           integer, intent(inout) :: col
           character (len=*), intent(in) :: name
           integer, intent(in) :: val
           if (pass == 1) then
              col = col+1
              write(io, fmt=int_fmt, advance='no') col
           else if (pass == 2) then
              write(io, fmt=txt_fmt, advance='no') trim(name)
           else if (pass == 3) then
              write(io, fmt=int_fmt, advance='no') val
           end if
         end subroutine write_integer


         subroutine write_val(io, col, pass, name, val) ! for header items only
           integer, intent(in) :: io, pass
           integer, intent(inout) :: col
           character (len=*), intent(in) :: name
           real(dp), intent(in) :: val
           if (pass == 1) then
              col = col+1
              write(io, fmt=int_fmt, advance='no') col
           else if (pass == 2) then
              write(io, fmt=txt_fmt, advance='no') trim(name)
           else if (pass == 3) then
              write(io, fmt=dbl_fmt, advance='no') val
           end if
         end subroutine write_val
         

      end subroutine do_history_info

      subroutine history_getval( &
            s, c, val, int_val, is_int_val, &
            nz, v_surf, csound_surf, envelope_fraction_left, epsnuc_out, ierr)
         use rates_def, only: i_rate
         use colors_lib, only: get_abs_mag_by_id, get_bc_by_id, get_lum_band_by_id
         use chem_lib, only: chem_M_div_h
         use rsp_def, only: rsp_phase_time0
         use gravity_darkening
         type (star_info), pointer :: s
         integer, intent(in) :: c, nz
         real(dp), intent(in) :: &
            v_surf, csound_surf, envelope_fraction_left, epsnuc_out(:)
         real(dp), intent(out) :: val
         integer, intent(out) :: int_val
         logical, intent(out) :: is_int_val
         integer, intent(out) :: ierr
         integer :: k, i, min_k
         real(dp) :: Ledd, L_rad, phi_Joss, power_photo, tmp, r, m_div_h, w_div_w_Kep, &
            min_gamma1
         real(dp), pointer :: v(:)
         logical :: v_flag

         include 'formats'

         ierr = 0
         is_int_val = .false.
         int_val = 0
         val = 0
         
         v_flag = .true.
         v(1:nz) => s% v(1:nz)  ! need this outside of conditional to keep compiler happy
         if (s% u_flag) then
            v(1:nz) => s% u(1:nz)
         else if (.not. s% v_flag) then
            v_flag = .false.
         end if
         
         if (c > log_lum_band_offset) THEN
            ! We want log Teff, Log g, M/H, Lum/lsun at the photosphere
            k = s% photosphere_cell_k
            m_div_h = chem_M_div_h(s% X(k),s% Z(k),s% job% initial_zfracs)
            if (k > 0) then
               i = c - log_lum_band_offset
               val = get_lum_band_by_id(i, safe_log10(s% photosphere_T), &
                  s% photosphere_logg, m_div_h, s% photosphere_L, ierr)
               if (ierr /= 0) return
            end if
            val=safe_log10(val*lsun)
         else if (c > lum_band_offset) THEN
            ! We want log Teff, Log g, M/H, Lum/lsun at the photosphere
            k = s% photosphere_cell_k
            m_div_h = chem_M_div_h(s% X(k),s% Z(k),s% job% initial_zfracs)
            if (k > 0) then
               i = c - lum_band_offset
               !val = get_lum_band_by_id(i,safe_log10(s% T(k)),safe_log10(s% grav(k)),m_div_h,s% L(k)/lsun, ierr)
               val = get_lum_band_by_id(i, safe_log10(s% photosphere_T), &
                  s% photosphere_logg, m_div_h, s% photosphere_L, ierr)
               val=val*lsun
               if (ierr /= 0) return
            end if
         else if (c > abs_mag_offset) THEN
            ! We want log Teff, Log g, M/H, Lum/lsun at the photosphere
            k = s% photosphere_cell_k
            m_div_h = chem_M_div_h(s% X(k),s% Z(k),s% job% initial_zfracs)
            if (k > 0) then
               i = c - abs_mag_offset
               !val = get_abs_mag_by_id(i,safe_log10(s% T(k)),safe_log10(s% grav(k)),m_div_h,s% L(k)/lsun, ierr)
               val = get_abs_mag_by_id(i, safe_log10(s% photosphere_T), &
                  s% photosphere_logg, m_div_h, s% photosphere_L, ierr)
               if (ierr /= 0) return
            end if
         else if (c > bc_offset) THEN
            ! We want log Teff, Log g, M/H at the photosphere
            k = s% photosphere_cell_k
            m_div_h = chem_M_div_h(s% X(k),s% Z(k),s% job% initial_zfracs)
            if (k > 0) then
               i = c - bc_offset
               !val = get_bc_by_id(i,safe_log10(s% T(k)),safe_log10(s% grav(k)),m_div_h, ierr)
               val = get_bc_by_id(i, safe_log10(s% photosphere_T), &
                  s% photosphere_logg, m_div_h, ierr)
               if (ierr /= 0) return
            end if
         else if (c > c_log_eps_burn_offset) then
            i = c - c_log_eps_burn_offset
            val = safe_log10(abs(s% center_eps_burn(i))) ! abs is for photo
         else if (c > max_eps_nuc_offset) then
            i = c - max_eps_nuc_offset
            val = safe_log10(max_eps_nuc_log_x(s% net_iso(i)))
         else if (c > cz_top_max_offset) then
            i = c - cz_top_max_offset
            val = safe_log10(cz_top_max_log_x(s% net_iso(i)))
         else if (c > cz_max_offset) then
            i = c - cz_max_offset
            val = safe_log10(cz_max_log_x(s% net_iso(i)))
         else if (c > log_surface_xa_offset) then
            i = c - log_surface_xa_offset
            val = safe_log10(surface_avg_x(s,s% net_iso(i)))
         else if (c > log_center_xa_offset) then
            i = c - log_center_xa_offset
            val = safe_log10(center_avg_x(s,s% net_iso(i)))
         else if (c > log_average_xa_offset) then
            i = c - log_average_xa_offset
            val = safe_log10(star_avg_x(s,s% net_iso(i)))
         else if (c > log_total_mass_offset) then
            i = c - log_total_mass_offset
            val = safe_log10(star_avg_x(s,s% net_iso(i))*s% xmstar/Msun)
         else if (c > total_mass_offset) then
            i = c - total_mass_offset
            val = star_avg_x(s,s% net_iso(i))*s% xmstar/Msun
         else if (c > category_offset) then
            i = c - category_offset
            val = category_L(i)
         else if (c > average_xa_offset) then
            i = c - average_xa_offset
            val = star_avg_x(s,s% net_iso(i))
         elseif (c > surface_xa_offset) then
            i = c - surface_xa_offset
            val = surface_avg_x(s,s% net_iso(i))
         else if (c > center_xa_offset) then
            i = c - center_xa_offset
            val = center_avg_x(s,s% net_iso(i))
         else

            select case(c)

            case(h_model_number)
               is_int_val = .true.
               int_val = s% model_number
               
            case(h_log_star_age)
               val = safe_log10(s% star_age)
            case(h_star_age)
               val = s% star_age
            case(h_log_star_age_sec)
               val = safe_log10(s% star_age*secyer)
            case(h_star_age_sec)
               val = s% star_age*secyer
            case(h_star_age_min)
               val = s% star_age*secyer/60
            case(h_star_age_hr)
               val = s% star_age*secyer/60/60
            case(h_star_age_day)
               val = s% star_age*secyer/60/60/24
            case(h_day)
               val = s% star_age*secyer/60/60/24

            case(h_time_step)
               val = s% time_step
            case(h_log_dt)
               val = safe_log10(s% time_step)
            case(h_time_step_sec)
               val = s% time_step*secyer
            case(h_log_dt_sec)
               val = safe_log10(s% time_step*secyer)
            case(h_time_step_days)
               val = s% time_step*secyer/60/60/24
            case(h_log_dt_days)
               val = safe_log10(s% time_step*secyer/60/60/24)

            case(h_log_star_mass)
               val = safe_log10(s% star_mass)
            case(h_star_mass)
               val = s% star_mass
            case(h_log_xmstar)
               val = safe_log10(s% xmstar)
            case(h_delta_mass)
               val = s% star_mass - s% initial_mass
            case(h_star_mdot)
               val = s% star_mdot
            case(h_log_abs_mdot)
               val = safe_log10(abs(s% star_mdot))

            case(h_m_center)
               val = s% M_center/Msun
            case(h_r_center)
               val = s% R_center/Rsun
            case(h_m_center_gm)
               val = s% M_center
            case(h_r_center_cm)
               val = s% R_center
            case(h_r_center_km)
               val = s% R_center*1d-5
            case(h_L_center)
               val = s% L_center/Lsun
            case(h_log_L_center_ergs_s)
               val = safe_log10(s% L_center)
            case(h_log_L_center)
               val = safe_log10(s% L_center/Lsun)
            case(h_v_center)
               val = s% v_center
            case(h_v_center_kms)
               val = s% v_center*1d-5
            case(h_infall_div_cs)
               if (s% v_center < 0d0) val = -s% v_center/s% csound(s% nz)

            case(h_mdot_timescale)
               val = s% star_mass/max(1d-99,abs(s% star_mdot))

            case(h_kh_div_mdot_timescales)
               val = s% kh_timescale/ &
                  (s% star_mass/max(1d-99,abs(s% star_mdot)))
            case(h_dlnR_dlnM)
               if (abs(s% star_mdot) > 1d-99) &
                  val = (s% lnR(1) - s% lnR_start(1)) / &
                     ((s% mstar - s% mstar_old)/(0.5d0*(s% mstar + s% mstar_old)))
            case(h_star_gravitational_mass)
               val = s% m_grav(1)/Msun
            case(h_star_mass_grav_div_mass)
               val = s% m_grav(1)/s% m(1)

            case(h_e_thermal)
               val = sum(s% dm(1:nz)*s% T(1:nz)*s% cp(1:nz))
            case(h_total_angular_momentum)
               val = s% total_angular_momentum
            case(h_log_total_angular_momentum)
               val = safe_log10(s% total_angular_momentum)
            case(h_species)
               int_val = s% species
               is_int_val = .true.
            case(h_Tsurf_factor)
               val = s% Tsurf_factor
            case(h_tau_factor)
               val = s% tau_factor
            case(h_tau_surface)
               val = s% tau_factor*s% tau_base
            case(h_log_tau_center)
               val = safe_log10(s% tau(s% nz))

            case(h_logT_max)
               val = s% log_max_temperature
            case(h_gamma1_min)
               val = s% min_gamma1

            case(h_logQ_max)
               val = maxval(s% lnd(1:nz)/ln10 - 2*s% lnT(1:nz)/ln10 + 12)
            case(h_logQ_min)
               val = minval(s% lnd(1:nz)/ln10 - 2*s% lnT(1:nz)/ln10 + 12)

            case(h_num_zones)
               int_val = nz
               is_int_val = .true.
            case(h_num_retries)
               int_val = s% num_retries
               is_int_val = .true.
               
            case(h_avg_skipped_setvars_per_step)
               val = dble(s% num_skipped_setvars)/max(1,s% model_number)
            case(h_avg_setvars_per_step)
               val = dble(s% num_setvars)/max(1,s% model_number)
            case(h_avg_solver_setvars_per_step)
               val = dble(s% num_solver_setvars)/max(1,s% model_number)

            case(h_total_num_solver_iterations)
               int_val = s% total_num_solver_iterations
               is_int_val = .true.
            case(h_total_num_solver_calls_made)
               int_val = s% total_num_solver_calls_made
               is_int_val = .true.
            case(h_total_num_solver_calls_converged)
               int_val = s% total_num_solver_calls_converged
               is_int_val = .true.
            case(h_total_num_solver_calls_failed)
               int_val = s% total_num_solver_calls_made - &
                  s% total_num_solver_calls_converged
               is_int_val = .true.

            case(h_total_num_solver_relax_iterations)
               int_val = s% total_num_solver_relax_iterations
               is_int_val = .true.
            case(h_total_num_solver_relax_calls_made)
               int_val = s% total_num_solver_relax_calls_made
               is_int_val = .true.
            case(h_total_num_solver_relax_calls_converged)
               int_val = s% total_num_solver_relax_calls_converged
               is_int_val = .true.
            case(h_total_num_solver_relax_calls_failed)
               int_val = s% total_num_solver_relax_calls_made - &
                  s% total_num_solver_relax_calls_converged
               is_int_val = .true.

            case(h_total_step_attempts)
               int_val = s% total_step_attempts
               is_int_val = .true.
            case(h_total_step_retries)
               int_val = s% total_step_retries
               is_int_val = .true.
            case(h_total_step_redos)
               int_val = s% total_step_redos
               is_int_val = .true.
            case(h_total_steps_taken)
               int_val = s% total_step_attempts - &
                  s% total_step_retries - s% total_step_redos
               is_int_val = .true.
            case(h_total_steps_finished)
               int_val = s% total_steps_finished
               is_int_val = .true.

            case(h_total_relax_step_attempts)
               int_val = s% total_relax_step_attempts
               is_int_val = .true.
            case(h_total_relax_step_retries)
               int_val = s% total_relax_step_retries
               is_int_val = .true.
            case(h_total_relax_step_redos)
               int_val = s% total_relax_step_redos
               is_int_val = .true.
            case(h_total_relax_steps_taken)
               int_val = s% total_relax_step_attempts - &
                  s% total_relax_step_retries - s% total_relax_step_redos
               is_int_val = .true.
            case(h_total_relax_steps_finished)
               int_val = s% total_relax_steps_finished
               is_int_val = .true.
               
            case(h_avg_num_solver_iters)
               val = dble(s% total_num_solver_iterations)/ &
                     dble(s% total_num_solver_calls_made)

            case(h_num_solver_iterations)
               int_val = s% num_solver_iterations
               is_int_val = .true.
            case(h_num_iters)
               int_val = s% num_solver_iterations
               is_int_val = .true.
               
            case(h_h1_czb_mass)
               val = s% h1_czb_mass
            case(h_surf_c12_minus_o16)
               val = s% surface_c12 - s% surface_o16
            case(h_surf_num_c12_div_num_o16)
               val = (16d0/12d0)*s% surface_c12/max(1d-99,s% surface_o16)
            case(h_conv_mx1_top)
               val = s% conv_mx1_top
            case(h_conv_mx1_bot)
               val = s% conv_mx1_bot
            case(h_conv_mx2_top)
               val = s% conv_mx2_top
            case(h_conv_mx2_bot)
               val = s% conv_mx2_bot
            case(h_mx1_top)
               val = s% mx1_top
            case(h_mx1_bot)
               val = s% mx1_bot
            case(h_mx2_top)
               val = s% mx2_top
            case(h_mx2_bot)
               val = s% mx2_bot
            case(h_conv_mx1_top_r)
               val = s% conv_mx1_top_r
            case(h_conv_mx1_bot_r)
               val = s% conv_mx1_bot_r
            case(h_conv_mx2_top_r)
               val = s% conv_mx2_top_r
            case(h_conv_mx2_bot_r)
               val = s% conv_mx2_bot_r
            case(h_mx1_top_r)
               val = s% mx1_top_r
            case(h_mx1_bot_r)
               val = s% mx1_bot_r
            case(h_mx2_top_r)
               val = s% mx2_top_r
            case(h_mx2_bot_r)
               val = s% mx2_bot_r
            case(h_epsnuc_M_1)
               val = epsnuc_out(1)
            case(h_epsnuc_M_2)
               val = epsnuc_out(2)
            case(h_epsnuc_M_3)
               val = epsnuc_out(3)
            case(h_epsnuc_M_4)
               val = epsnuc_out(4)
            case(h_epsnuc_M_5)
               val = epsnuc_out(5)
            case(h_epsnuc_M_6)
               val = epsnuc_out(6)
            case(h_epsnuc_M_7)
               val = epsnuc_out(7)
            case(h_epsnuc_M_8)
               val = epsnuc_out(8)
               
            case(h_power_h_burn)
               val = s% power_h_burn
            case(h_log_LH)
               val = safe_log10(s% power_h_burn)
               
            case(h_power_he_burn)
               val = s% power_he_burn
            case(h_log_LHe)
               val = safe_log10(s% power_he_burn)
               
            case(h_power_photo)
               val = s% power_photo
            case(h_Lnuc_photo)
               val = safe_log10(abs(s% power_photo))
               
            case(h_power_z_burn)
               val = s% power_z_burn
            case(h_log_LZ)
               val = safe_log10(s% power_z_burn)
               
            case(h_log_Lneu)
               val = safe_log10(s% power_neutrinos)
            case(h_log_Lneu_nuc)
               val = safe_log10(s% power_nuc_neutrinos)
            case(h_log_Lneu_nonnuc)
               val = safe_log10(s% power_nonnuc_neutrinos)

           case(h_Lsurf_m)
               val = s% m(1)/Msun
           case(h_luminosity)
               val = s% L_surf
           case(h_log_L)
               val = safe_log10(s% L_surf)
           case(h_luminosity_ergs_s)
               val = s% L_surf*Lsun
           case(h_log_L_ergs_s)
               val = safe_log10(s% L_surf*Lsun)

            case(h_photosphere_cell_log_density)
               val = s% lnd(s% photosphere_cell_k)/ln10
            case(h_photosphere_cell_density)
               val = s% rho(s% photosphere_cell_k)
            case(h_photosphere_cell_log_opacity)
               val = safe_log10(s% opacity(s% photosphere_cell_k))
            case(h_photosphere_cell_opacity)
               val = s% opacity(s% photosphere_cell_k)
            case(h_photosphere_cell_log_free_e)
               val = s% lnfree_e(s% photosphere_cell_k)/ln10
            case(h_photosphere_cell_free_e)
               val = exp(s% lnfree_e(s% photosphere_cell_k))

            case(h_log_Teff)
               val = safe_log10(s% Teff)
            case(h_Teff)
               val = s% Teff
            case(h_effective_T)
               val = s% Teff
               
            case(h_photosphere_cell_k)
               int_val = s% photosphere_cell_k
               is_int_val = .true.
            case(h_photosphere_cell_log_T)
               val = s% lnT(s% photosphere_cell_k)/ln10
            case(h_photosphere_cell_T)
               val = s% T(s% photosphere_cell_k)
            case(h_photosphere_T)
               val = s% photosphere_T
            case(h_photosphere_black_body_T)
               val = s% photosphere_black_body_T
            case(h_photosphere_logg)
               val = s% photosphere_logg
            case(h_photosphere_m)
               val = s% photosphere_m
            case(h_photosphere_xm)
               val = s% star_mass - s% photosphere_m
            case(h_photosphere_L)
               val = s% photosphere_L
            case(h_photosphere_r)
               val = s% photosphere_r
            case(h_photosphere_log_L)
               val = safe_log10(s% photosphere_L)
            case(h_photosphere_log_r)
               val = safe_log10(s% photosphere_r)
            case(h_photosphere_csound)
               val = s% photosphere_csound
            case(h_photosphere_opacity)
               val = s% photosphere_opacity
            case(h_photosphere_column_density)
               val = s% photosphere_column_density
            case(h_photosphere_log_column_density)
               val = safe_log10(s% photosphere_column_density)
            case(h_photosphere_v_km_s)
               val = s% photosphere_v/1d5
            case(h_v_phot_km_s)
               val = s% photosphere_v/1d5
            case(h_photosphere_v_div_cs)
               val = s% photosphere_v/s% photosphere_csound
               
            case(h_one_div_yphot)
               val = 1d0/s% photosphere_column_density
            case(h_log_one_div_yphot)
               val = safe_log10(1d0/s% photosphere_column_density)
            case(h_min_opacity)
               val = minval(s% opacity(1:s% nz))
            case(h_log_min_opacity)
               val = safe_log10(minval(s% opacity(1:s% nz)))

            case(h_radius_cm)
               val = s% R(1)
            case(h_log_R_cm)
               val = safe_log10(s% R(1))
            case(h_radius)
               val = s% R(1)/Rsun
            case(h_log_R)
               val = safe_log10(s% R(1)/Rsun)

            case(h_gravity)
               val = s% grav(1)
            case(h_log_g)
               val = safe_log10(s% grav(1))

            case(h_log_cntr_dr_cm)
               val = safe_log10(s% r(s% nz) - s% R_center)
            case(h_log_max_T)
               val = s% log_max_temperature
            case(h_log_cntr_T)
               val = s% log_center_temperature
            case(h_log_cntr_Rho)
               val = s% log_center_density
            case(h_log_cntr_P)
               val = s% log_center_pressure
            case(h_log_center_T)
               val = s% log_center_temperature
            case(h_log_center_Rho)
               val = s% log_center_density
            case(h_log_center_P)
               val = s% log_center_pressure

            case(h_max_T)
               val = exp10(s% log_max_temperature)
            case(h_center_T)
               val = exp10(s% log_center_temperature)
            case(h_center_Rho)
               val = exp10(s% log_center_density)
            case(h_center_P)
               val = exp10(s% log_center_pressure)

            case(h_log_mesh_adjust_IE_conservation)
               val = safe_log10(s% mesh_adjust_IE_conservation)
            case(h_log_mesh_adjust_PE_conservation)
               val = safe_log10(s% mesh_adjust_PE_conservation)
            case(h_log_mesh_adjust_KE_conservation)
               val = safe_log10(s% mesh_adjust_KE_conservation)

            case(h_rms_dvdt_div_v)
               val = eval_rms_dvdt_div_v(s, 1, s% nz)
           case(h_total_IE_div_IE_plus_KE)
               val = s% total_internal_energy_end / &
                        (s% total_internal_energy_end + s% total_radial_kinetic_energy_end)
           case(h_total_entropy)
               val = dot_product(s% dm(1:nz), s% entropy(1:nz))

           case(h_total_internal_energy_after_adjust_mass)
               val = s% total_internal_energy_after_adjust_mass
           case(h_total_gravitational_energy_after_adjust_mass)
               val = s% total_gravitational_energy_after_adjust_mass
           case(h_total_radial_kinetic_energy_after_adjust_mass)
               val = s% total_radial_kinetic_energy_after_adjust_mass
           case(h_total_rotational_kinetic_energy_after_adjust_mass)
               val = s% total_rotational_kinetic_energy_after_adjust_mass
           case(h_total_turbulent_energy_after_adjust_mass)
               val = s% total_turbulent_energy_after_adjust_mass
           case(h_total_energy_after_adjust_mass)
               val = s% total_energy_after_adjust_mass
               
           case(h_total_internal_energy)
               val = s% total_internal_energy_end
           case(h_total_gravitational_energy)
               val = s% total_gravitational_energy_end
           case(h_total_radial_kinetic_energy)
               val = s% total_radial_kinetic_energy_end
           case(h_total_rotational_kinetic_energy)
               val = s% total_rotational_kinetic_energy_end
           case(h_total_turbulent_energy)
               val = s% total_turbulent_energy_end
           case(h_total_energy)
               val = s% total_energy_end
           case(h_total_energy_foe)
               val = s% total_energy_end*1d-51
               
           case(h_log_total_internal_energy)
               val = safe_log10(s% total_internal_energy_end)
           case(h_log_total_gravitational_energy)
               val = safe_log10(abs(s% total_gravitational_energy_end))
           case(h_log_total_radial_kinetic_energy)
               val = safe_log10(s% total_radial_kinetic_energy_end)
           case(h_log_total_rotational_kinetic_energy)
               val = safe_log10(s% total_rotational_kinetic_energy_end)
           case(h_log_total_turbulent_energy)
               val = safe_log10(s% total_turbulent_energy_end)
           case(h_log_total_energy)
               val = safe_log10(abs(s% total_energy_end))

           case(h_avg_abs_v_div_cs)
               if (v_flag) &
                  val = sum(abs(v(1:nz))/s% csound(1:nz))/nz
           case(h_log_avg_abs_v_div_cs)
               if (v_flag) &
                  val = safe_log10(sum(abs(v(1:nz))/s% csound(1:nz))/nz)
           case(h_max_abs_v_div_cs)
               if (v_flag) &
                  val = maxval(abs(v(1:nz))/s% csound(1:nz))
           case(h_log_max_abs_v_div_cs)
               if (v_flag) &
                  val = safe_log10(maxval(abs(v(1:nz))/s% csound(1:nz)))

           case(h_avg_abs_v)
               if (v_flag) &
                  val = sum(abs(v(1:nz)))/nz
           case(h_log_avg_abs_v)
               if (v_flag) &
                  val = safe_log10(sum(abs(v(1:nz)))/nz)
           case(h_max_abs_v)
               if (v_flag) &
                  val = maxval(abs(v(1:nz)))
           case(h_log_max_abs_v)
               if (v_flag) &
                  val = safe_log10(maxval(abs(v(1:nz))))

           case(h_virial_thm_P_avg)
               val = s% virial_thm_P_avg
           case(h_virial_thm_rel_err)
               val = s% virial_thm_P_avg
               val = (val + s% total_gravitational_energy_end)/val
               
           case(h_total_eps_grav)
               val = s% total_eps_grav
           case(h_work_outward_at_surface)
               val = s% work_outward_at_surface
           case(h_work_inward_at_center)
               val = s% work_inward_at_center
           case(h_total_nuclear_heating)
               val = s% total_nuclear_heating
           case(h_total_non_nuc_neu_cooling)
               val = s% total_non_nuc_neu_cooling
           case(h_total_irradiation_heating)
               val = s% total_irradiation_heating
           case(h_total_WD_sedimentation_heating)
               if (s% do_element_diffusion) val = s% total_WD_sedimentation_heating
           case(h_total_extra_heating)
               val = s% total_extra_heating

           case(h_total_energy_sources_and_sinks)
               val = s% total_energy_sources_and_sinks
               
           case(h_error_in_energy_conservation)
               val = s% error_in_energy_conservation
           case(h_rel_error_in_energy_conservation)
               val = s% error_in_energy_conservation/abs(s% total_energy_end)
           case(h_log_rel_error_in_energy_conservation)
               val = safe_log10(abs(s% error_in_energy_conservation/s% total_energy_end))
               
           case(h_tot_E_equ_err)
               val = sum(s% E_residual(1:nz)*s% dm(1:nz))
           case(h_tot_E_err)
               val = s% error_in_energy_conservation
           case(h_rel_E_err)
               if (s% total_energy_end /= 0d0) &
                  val = s% error_in_energy_conservation/abs(s% total_energy_end)
           case(h_abs_rel_E_err)
               if (s% total_energy_end /= 0d0) &
                  val = abs(s% error_in_energy_conservation/s% total_energy_end)
           case(h_log_rel_E_err)
               if (s% total_energy_end /= 0d0) &
                  val = safe_log10(abs(s% error_in_energy_conservation/s% total_energy_end))

           case(h_cumulative_energy_error)
               val = s% cumulative_energy_error
           case(h_rel_cumulative_energy_error)
               if (s% total_energy_end /= 0d0) &
                  val = s% cumulative_energy_error/abs(s% total_energy_end)
           case(h_log_rel_cumulative_energy_error)
               if (s% total_energy_end /= 0d0) &
                  val = safe_log10(abs(s% cumulative_energy_error/s% total_energy_end))
           case(h_rel_run_E_err)
               if (s% total_energy_end /= 0d0) &
                  val = s% cumulative_energy_error/s% total_energy_end
           case(h_log_rel_run_E_err)
               if (s% total_energy_end /= 0d0) &
                  val = safe_log10(abs(s% cumulative_energy_error/s% total_energy_end))

           case(h_log_residual_norm)
               val = safe_log10(s% residual_norm)
           case(h_log_max_residual)
               val = safe_log10(s% max_residual)

           case(h_log_max_dvdt_residual)
               val = safe_log10(maxval(abs(s% v_residual(1:nz))))
           case(h_log_max_lnd_residual)
               val = safe_log10(maxval(abs(s% lnd_residual(1:nz))))
           case(h_log_max_dEdt_residual)
               val = safe_log10(maxval(abs(s% E_residual(1:nz))))
           case(h_log_max_drdt_residual)
               val = safe_log10(maxval(abs(s% lnR_residual(1:nz))))

           case(h_avg_v_residual)
               val = dot_product(s% dq(1:nz),s% v_residual(1:nz))
           case(h_log_avg_v_residual)
               val = safe_log10(abs(dot_product(s% dq(1:nz),s% v_residual(1:nz))))

           case(h_max_abs_v_residual)
               k = maxloc(abs(s% v_residual(1:nz)),dim=1)
               val = s% v_residual(k)
           case(h_log_max_abs_v_residual)
               val = safe_log10(maxval(abs(s% v_residual(1:nz))))

           case(h_avg_E_residual)
               val = dot_product(s% dq(1:nz),s% E_residual(1:nz))/ln10
           case(h_log_avg_E_residual)
               val = safe_log10(abs(dot_product(s% dq(1:nz),s% E_residual(1:nz)))/ln10)

           case(h_max_abs_E_residual)
               val = maxval(abs(s% E_residual(1:nz)))/ln10
           case(h_log_max_abs_E_residual)
               val = safe_log10(maxval(abs(s% E_residual(1:nz)))/ln10)

            case(h_u_surf_km_s)
               if (s% u_flag) val = s% u_face_18(1)%val*1d-5
            case(h_u_surf)
               if (s% u_flag) val = s% u_face_18(1)%val
            case(h_u_div_csound_max)
               if (s% u_flag) val = maxval(abs(s% u(1:nz))/s% csound(1:nz))
            case(h_u_div_csound_surf)
               if (s% u_flag) val = s% u_face_18(1)%val/s% csound_face(1)

            case(h_surf_escape_v)
               val = sqrt(2*s% cgrav(1)*s% m(1)/(s% r(1)))
            case(h_v_surf_div_escape_v)
               val = v_surf/sqrt(2*s% cgrav(1)*s% m(1)/(s% r(1)))
            case(h_v_surf_km_s)
               val = v_surf*1d-5
            case(h_v_surf)
               val = v_surf
            case(h_v_surf_div_v_kh)
               val = v_surf/(s% photosphere_r/s% kh_timescale)
            case(h_v_div_csound_max)
               if (v_flag) val = maxval(abs(v(1:nz))/s% csound_face(1:nz))
            case(h_v_div_csound_surf)
               val = v_surf/csound_surf
            case(h_v_div_cs)
               if (s% u_flag) then
                  val = s% u(1) / s% csound(1)
               else if (s% v_flag) then
                  val = s% v(1) / s% csound(1)
               else
                  val = 0d0 ! s% r(1)*s% dlnR_dt(1)
               end if
            case(h_remnant_M)
               val = get_remnant_mass(s)/Msun
            case(h_ejecta_M)
               val = get_ejecta_mass(s)/Msun

            case(h_log_L_div_Ledd)
               Ledd = eval_Ledd(s, ierr)
               if (ierr /= 0 .or. Ledd == 0d0) then
                  ierr = 0
               else
                  val = safe_log10(s% L_surf*Lsun/Ledd)
               end if

            case(h_lum_div_Ledd)
               Ledd = eval_Ledd(s, ierr)
               if (ierr /= 0 .or. Ledd == 0d0) then
                  ierr = 0
               else
                  val = s% L_surf*Lsun/Ledd
               end if

            case(h_gradT_excess_alpha)
               val = s% gradT_excess_alpha
            case(h_gradT_excess_min_beta)
               val = s% gradT_excess_min_beta
            case(h_gradT_excess_max_lambda)
               val = s% gradT_excess_max_lambda

            case(h_max_L_rad_div_Ledd)
               do k=1,nz
                  tmp = get_Lrad_div_Ledd(s,k)
                  if (tmp > val) val = tmp
               end do

            case(h_max_L_rad_div_Ledd_div_phi_Joss)
               do k=1,nz
                  tmp = get_Lrad_div_Ledd(s,k)
                  phi_Joss = get_phi_Joss(s,k)
                  if (tmp/phi_Joss > val) val = tmp/phi_Joss
               end do

            case(h_i_rot_total)
               val = dot_product(s% dm_bar(1:nz), s%i_rot(1:nz))
            case(h_surf_avg_j_rot)
               val = if_rot(s% j_rot_avg_surf)
            case(h_surf_avg_omega)
               val = if_rot(s% omega_avg_surf)
            case(h_surf_avg_omega_crit)
               val = if_rot(s% omega_crit_avg_surf)
            case(h_surf_avg_omega_div_omega_crit)
               val = if_rot(s% w_div_w_crit_avg_surf)

            case(h_surf_avg_v_rot)
               val = if_rot(s% v_rot_avg_surf)*1d-5 ! km/sec
            case(h_surf_avg_v_crit)
               val = if_rot(s% v_crit_avg_surf)*1d-5 ! km/sec
            case(h_surf_avg_v_div_v_crit)
               val = if_rot(s% v_div_v_crit_avg_surf)

            case(h_surf_avg_Lrad_div_Ledd)
               val = s% Lrad_div_Ledd_avg_surf
            case(h_surf_avg_opacity)
               val = s% opacity_avg_surf
            case(h_surf_avg_logT)
               val = s% logT_avg_surf
            case(h_surf_avg_logRho)
               val = s% logRho_avg_surf

            case(h_luminosity_for_BB_outer_BC)
               if (s% tau_for_L_BB >= 0) then
                  val = s% L_for_BB_outer_BC/Lsun
               else
                  val = s% L(1)/Lsun
               end if
            case(h_logL_for_BB_outer_BC)
               if (s% tau_for_L_BB >= 0) then
                  val = s% L_for_BB_outer_BC/Lsun
               else
                  val = s% L(1)/Lsun
               end if
               val = safe_log10(val)

            case(h_v_wind_Km_per_s)
               val = 1d-5*s% opacity(1)*max(0d0,-s% mstar_dot)/ &
                        (pi4*s% photosphere_r*Rsun*s% tau_base)

            case (h_kh_mdot_limit)
               if(s% rotation_flag) then
                  val = s% rotational_mdot_kh_fac*s% star_mass/s% kh_timescale
               else
                  val = 0d0
               end if
            case (h_log_rotational_mdot_boost)
               val = safe_log10(if_rot(s% rotational_mdot_boost))
            case (h_rotational_mdot_boost)
               val = if_rot(s% rotational_mdot_boost)

            case(h_min_Pgas_div_P)
               val = minval(s% Pgas(1:nz)/s% P(1:nz))

            case(h_center_degeneracy)
               val = s% center_degeneracy

            case(h_log_center_eps_nuc)
               val = safe_log10(s% center_eps_nuc)
            case(h_center_eps_nuc)
               val = s% center_eps_nuc
            case(h_d_center_eps_nuc_dlnT)
               val = s% d_center_eps_nuc_dlnT
            case(h_d_center_eps_nuc_dlnd)
               val = s% d_center_eps_nuc_dlnd

            case(h_center_dlogT)
               val = s% dt*center_value(s, s% dlnT_dt)/ln10
            case(h_center_dlogRho)
               val = s% dt*center_value(s, s% dlnd_dt)/ln10

            case(h_center_dlnT_dt)
               val = center_value(s, s% dlnT_dt)
            case(h_center_dlnd_dt)
               val = center_value(s, s% dlnd_dt)

            case(h_center_dL_dm)
               val = center_value(s, s% dL_dm_expected)
            case(h_center_eps_grav)
               val = center_value(s, s% eps_grav)

            case(h_center_non_nuc_neu)
               val = s% center_non_nuc_neu
            case(h_center_gamma)
               val = s% center_gamma
            case(h_center_zbar)
               val = s% center_zbar
            case(h_center_abar)
               val = s% center_abar
            case(h_center_mu)
               val = s% center_mu
            case(h_center_ye)
               val = s% center_ye
            case(h_center_entropy)
               val = s% center_entropy
            case(h_max_entropy)
               val = s% max_entropy
            case(h_compactness)
               if (s% m(1) > 2.5d0*Msun) then
                  do k=nz-1, 1, -1
                     if (s% m(k) > 2.5d0*Msun) exit
                  end do
                  r = s% r(k+1) + (s% r(k)-s% r(k+1))*(2.5d0*Msun - s% m(k+1))/s% dm(k)
                  val = 2.5d0/(r/1d8)
               end if
            case(h_compactness_parameter)
               if (s% m(1) > 2.5d0*Msun) then
                  do k=nz-1, 1, -1
                     if (s% m(k) > 2.5d0*Msun) exit
                  end do
                  r = s% r(k+1) + (s% r(k)-s% r(k+1))*(2.5d0*Msun - s% m(k+1))/s% dm(k)
                  val = 2.5d0/(r/1d8)
               end if
            case(h_max_infall_speed)
               if (s% u_flag) then
                  val = -minval(s% u(1:s% nz))*1d-5 ! convert to km/sec
               else if (s% v_flag) then
                  val = -minval(s% v(1:s% nz))*1d-5 ! convert to km/sec
               end if
            case(h_fe_core_infall)
               val = s% fe_core_infall*1d-5 ! convert to km/sec
            case(h_non_fe_core_infall)
               val = s% non_fe_core_infall*1d-5 ! convert to km/sec
            case(h_non_fe_core_rebound)
               val = s% non_fe_core_rebound*1d-5 ! convert to km/sec
            case(h_center_omega)
               val = if_rot(s% center_omega)
            case(h_center_omega_div_omega_crit)
               val = if_rot(s% center_omega_div_omega_crit)

            case(h_surf_r_equatorial_div_r_polar)
               if(s%rotation_flag) then
                  val = s% r_equatorial(1)/s% r_polar(1)
               else
                  val = 1.0d0
               end if
            case(h_surf_r_equatorial_div_r)
               if(s%rotation_flag) then
                  val = s% r_equatorial(1)/s% r(1)
               else
                  val = 1.0d0
               end if
            case(h_surf_r_polar_div_r)
               if(s%rotation_flag) then
                  val = s% r_polar(1)/s% r(1)
               else
                  val = 1.0d0
               end if
            case(h_h_rich_layer_mass)
               val = s% star_mass - s% he_core_mass
            case(h_he_rich_layer_mass)
               val = max(0d0, s% he_core_mass - s% co_core_mass)
            case(h_co_rich_layer_mass)
               val = max(0d0, s% co_core_mass - s% he_core_mass)

            case(h_he_core_mass)
               val = s% he_core_mass
            case(h_he_core_radius)
               val = s% he_core_radius
            case(h_he_core_lgT)
               val = s% he_core_lgT
            case(h_he_core_lgRho)
               val = s% he_core_lgRho
            case(h_he_core_L)
               val = s% he_core_L
            case(h_he_core_v)
               val = s% he_core_v
            case(h_he_core_omega)
               val = if_rot(s% he_core_omega)
            case(h_he_core_omega_div_omega_crit)
               val = if_rot(s% he_core_omega_div_omega_crit)
            case(h_he_core_k)
               int_val = s% he_core_k
               is_int_val = .true.

            case(h_co_core_mass)
               val = s% co_core_mass
            case(h_co_core_radius)
               val = s% co_core_radius
            case(h_co_core_lgT)
               val = s% co_core_lgT
            case(h_co_core_lgRho)
               val = s% co_core_lgRho
            case(h_co_core_L)
               val = s% co_core_L
            case(h_co_core_v)
               val = s% co_core_v
            case(h_co_core_omega)
               val = if_rot(s% co_core_omega)
            case(h_co_core_omega_div_omega_crit)
               val = if_rot(s% co_core_omega_div_omega_crit)
            case(h_co_core_k)
               int_val = s% co_core_k
               is_int_val = .true.

            case(h_fe_core_mass)
               val = s% fe_core_mass
            case(h_fe_core_radius)
               val = s% fe_core_radius
            case(h_fe_core_lgT)
               val = s% fe_core_lgT
            case(h_fe_core_lgRho)
               val = s% fe_core_lgRho
            case(h_fe_core_L)
               val = s% fe_core_L
            case(h_fe_core_v)
               val = s% fe_core_v
            case(h_fe_core_omega)
               val = if_rot(s% fe_core_omega)
            case(h_fe_core_omega_div_omega_crit)
               val = if_rot(s% fe_core_omega_div_omega_crit)
            case(h_fe_core_k)
               int_val = s% fe_core_k
               is_int_val = .true.

            case(h_neutron_rich_core_mass)
               val = s% neutron_rich_core_mass
            case(h_neutron_rich_core_radius)
               val = s% neutron_rich_core_radius
            case(h_neutron_rich_core_lgT)
               val = s% neutron_rich_core_lgT
            case(h_neutron_rich_core_lgRho)
               val = s% neutron_rich_core_lgRho
            case(h_neutron_rich_core_L)
               val = s% neutron_rich_core_L
            case(h_neutron_rich_core_v)
               val = s% neutron_rich_core_v
            case(h_neutron_rich_core_omega)
               val = if_rot(s% neutron_rich_core_omega)
            case(h_neutron_rich_core_omega_div_omega_crit)
               val = if_rot(s% neutron_rich_core_omega_div_omega_crit)
            case(h_neutron_rich_core_k)
               int_val = s% neutron_rich_core_k
               is_int_val = .true.

            case(h_envelope_mass)
               val = s% star_mass - s% he_core_mass
            case(h_envelope_fraction_left)
               val = envelope_fraction_left
            case(h_tau10_mass)
               val = s% tau10_mass
            case(h_tau10_radius)
               val = s% tau10_radius
            case(h_tau10_lgP)
               val = s% tau10_lgP
            case(h_tau10_T)
               val = exp10(s% tau10_lgT)
            case(h_tau10_lgT)
               val = s% tau10_lgT
            case(h_tau10_lgRho)
               val = s% tau10_lgRho
            case(h_tau10_L)
               val = s% tau10_L
            case(h_tau100_mass)
               val = s% tau100_mass
            case(h_tau100_radius)
               val = s% tau100_radius
            case(h_tau100_lgP)
               val = s% tau100_lgP
            case(h_tau100_T)
               val = exp10(s% tau100_lgT)
            case(h_tau100_lgT)
               val = s% tau100_lgT
            case(h_tau100_lgRho)
               val = s% tau100_lgRho
            case(h_tau100_L)
               val = s% tau100_L
            case(h_dynamic_timescale)
               val = s% dynamic_timescale
            case(h_kh_timescale)
               val = s% kh_timescale
            case(h_nuc_timescale)
               val = s% nuc_timescale
            case(h_eps_grav_integral)
               val = dot_product(s% dm(1:nz), s% eps_grav(1:nz))/Lsun
            case(h_extra_L)
               val = dot_product(s% dm(1:nz), s% extra_heat(1:nz))/Lsun
            case(h_log_extra_L)
               val = safe_log10(dot_product(s% dm(1:nz), s% extra_heat(1:nz))/Lsun)
            case(h_log_abs_Lgrav)
               val = safe_log10(abs(dot_product(s% dm(1:nz), s% eps_grav(1:nz))/Lsun))
            case(h_log_Lnuc)
               power_photo = dot_product(s% dm(1:nz), s% eps_nuc_categories(iphoto,1:nz))/Lsun
               val = safe_log10(s% power_nuc_burn - power_photo)
            case(h_log_Lnuc_ergs_s)
               power_photo = dot_product(s% dm(1:nz), s% eps_nuc_categories(iphoto,1:nz))
               val = safe_log10(s% power_nuc_burn*Lsun - power_photo)
            case(h_log_power_nuc_burn)
               val = safe_log10(s% power_nuc_burn)
            case(h_log_Lnuc_sub_log_L)
               power_photo = dot_product(s% dm(1:nz), s% eps_nuc_categories(iphoto,1:nz))/Lsun
               val = safe_log10(s% power_nuc_burn - power_photo)
               val = val - safe_log10(s% L_surf*Lsun/Lsun)
            case(h_mass_loc_of_max_eps_nuc)
               k = maxloc(s% eps_nuc(1:nz), dim=1)
               val = (s% m(k) - s% dm(k)/2)/Msun
            case(h_mass_ext_to_max_eps_nuc)
               k = maxloc(s% eps_nuc(1:nz), dim=1)
               val = (1d0 - s% q(k) + 0.5d0*s% dq(k))*s% xmstar/Msun

            case(h_trace_mass_location)
               val = s% trace_mass_location
            case(h_trace_mass_radius)
               val = s% trace_mass_radius
            case(h_trace_mass_lgT)
               val = s% trace_mass_lgT
            case(h_trace_mass_lgRho)
               val = s% trace_mass_lgRho
            case(h_trace_mass_L)
               val = s% trace_mass_L
            case(h_trace_mass_v)
               val = s% trace_mass_v
            case(h_trace_mass_omega)
               val = if_rot(s% trace_mass_omega)
            case(h_trace_mass_omega_div_omega_crit)
               val = if_rot(s% trace_mass_omega_div_omega_crit)

            case(h_trace_mass_lgP)
               val = s% trace_mass_lgP
            case(h_trace_mass_g)
               val = s% trace_mass_g
            case(h_trace_mass_X)
               val = s% trace_mass_X
            case(h_trace_mass_Y)
               val = s% trace_mass_Y
            case(h_trace_mass_edv_H)
               val = s% trace_mass_edv_H
            case(h_trace_mass_edv_He)
               val = s% trace_mass_edv_He
            case(h_trace_mass_scale_height)
               val = s% trace_mass_scale_height
            case(h_trace_mass_dlnX_dr)
               val = s% trace_mass_dlnX_dr
            case(h_trace_mass_dlnY_dr)
               val = s% trace_mass_dlnY_dr
            case(h_trace_mass_dlnRho_dr)
               val = s% trace_mass_dlnRho_dr

            case(h_diffusion_time_H_He_bdy)
               if (s% he_core_k > 0) then
                  val = (s% tau(s% he_core_k) - s% tau_factor*s% tau_base)* &
                     s% r(s% he_core_k)/clight
               end if
            case(h_temperature_H_He_bdy)
               if (s% he_core_k > 0) val = s% T(s% he_core_k)

            case(h_max_abs_v_velocity)
               val = s% max_abs_v_velocity
            case(h_max_abs_v_csound)
               val = s% max_abs_v_csound
            case(h_max_abs_v_v_div_cs)
               val = s% max_abs_v_v_div_cs
            case(h_max_abs_v_lgT)
               val = s% max_abs_v_lgT
            case(h_max_abs_v_lgRho)
               val = s% max_abs_v_lgRho
            case(h_max_abs_v_lgP)
               val = s% max_abs_v_lgP
            case(h_max_abs_v_mass)
               val = s% max_abs_v_mass
            case(h_max_abs_v_radius)
               val = s% max_abs_v_radius
            case(h_max_abs_v_radius_cm)
               val = s% max_abs_v_radius*Rsun
            case(h_max_abs_v_lgR)
               val = safe_log10(s% max_abs_v_radius)
            case(h_max_abs_v_lgR_cm)
               val = safe_log10(s% max_abs_v_radius*Rsun)
            case(h_max_abs_v_L)
               val = s% max_abs_v_L
            case(h_max_abs_v_gamma1)
               val = s% max_abs_v_gamma1
            case(h_max_abs_v_entropy)
               val = s% max_abs_v_entropy
            case(h_max_abs_v_eps_nuc)
               val = s% max_abs_v_eps_nuc
            case(h_max_abs_v_E0) ! 4/3 pi R^3 crad T^4
               val = s% max_abs_v_radius*Rsun
               val = four_thirds_pi*val*val*val*crad*exp10(4*s% max_abs_v_lgT)
               
            case(h_total_ni_co_56)
               if (s% net_iso(ico56) > 0 .and. s% net_iso(ini56) > 0) &
                  val = dot_product(s% dm(1:nz), &
                     s% xa(s% net_iso(ico56),1:nz) + &
                     s% xa(s% net_iso(ini56),1:nz))/Msun
               
            case(h_inner_mach1_velocity)
               val = s% inner_mach1_velocity
            case(h_inner_mach1_csound)
               val = s% inner_mach1_csound
            case(h_inner_mach1_v_div_cs)
               if (s% inner_mach1_csound > 0) &
                  val = s% inner_mach1_velocity/s% inner_mach1_csound
            case(h_inner_mach1_lgT)
               val = s% inner_mach1_lgT
            case(h_inner_mach1_lgRho)
               val = s% inner_mach1_lgRho
            case(h_inner_mach1_lgP)
               val = s% inner_mach1_lgP
            case(h_inner_mach1_q)
               val = s% inner_mach1_q
            case(h_inner_mach1_tau)
               val = s% inner_mach1_tau
            case(h_inner_mach1_mass)
               val = s% inner_mach1_mass
            case(h_inner_mach1_radius)
               val = s% inner_mach1_radius
            case(h_inner_mach1_gamma1)
               val = s% inner_mach1_gamma1
            case(h_inner_mach1_entropy)
               val = s% inner_mach1_entropy
            case(h_inner_mach1_k)
               int_val = s% inner_mach1_k
               is_int_val = .true.

            case(h_outer_mach1_velocity)
               val = s% outer_mach1_velocity
            case(h_outer_mach1_csound)
               val = s% outer_mach1_csound
            case(h_outer_mach1_v_div_cs)
               if (s% outer_mach1_csound > 0) &
                  val = s% outer_mach1_velocity/s% outer_mach1_csound
            case(h_outer_mach1_lgT)
               val = s% outer_mach1_lgT
            case(h_outer_mach1_lgRho)
               val = s% outer_mach1_lgRho
            case(h_outer_mach1_lgP)
               val = s% outer_mach1_lgP
            case(h_outer_mach1_q)
               val = s% outer_mach1_q
            case(h_outer_mach1_tau)
               val = s% outer_mach1_tau
            case(h_outer_mach1_mass)
               val = s% outer_mach1_mass
            case(h_outer_mach1_radius)
               val = s% outer_mach1_radius
            case(h_outer_mach1_gamma1)
               val = s% outer_mach1_gamma1
            case(h_outer_mach1_entropy)
               val = s% outer_mach1_entropy
            case(h_outer_mach1_k)
               int_val = s% outer_mach1_k
               is_int_val = .true.

            case(h_shock_velocity)
               if (s% shock_k > 0) val = s% shock_velocity
            case(h_shock_csound)
               if (s% shock_k > 0) val = s% shock_csound
            case(h_shock_v_div_cs)
               if (s% shock_csound > 0) &
                  val = s% shock_velocity/s% shock_csound
            case(h_shock_lgT)
               if (s% shock_k > 0) val = s% shock_lgT
            case(h_shock_lgRho)
               if (s% shock_k > 0) val = s% shock_lgRho
            case(h_shock_lgP)
               if (s% shock_k > 0) val = s% shock_lgP
            case(h_shock_q)
               if (s% shock_k > 0) val = s% shock_q
            case(h_shock_tau)
               if (s% shock_k > 0) val = s% shock_tau
            case(h_shock_mass)
               if (s% shock_k > 0) val = s% shock_mass
            case(h_shock_mass_gm)
               if (s% shock_k > 0) val = s% shock_mass*Msun
            case(h_shock_radius)
               if (s% shock_k > 0) val = s% shock_radius
            case(h_shock_radius_cm)
               if (s% shock_k > 0) val = s% shock_radius*Rsun
            case(h_shock_gamma1)
               if (s% shock_k > 0) val = s% shock_gamma1
            case(h_shock_entropy)
               if (s% shock_k > 0) val = s% shock_entropy
            case(h_shock_pre_lgRho)
               if (s% shock_k > 0) val = s% shock_pre_lgRho
            case(h_shock_k)
               if (s% shock_k > 0) int_val = s% shock_k
               is_int_val = .true.

            case(h_max_T_shell_binding_energy)
               val = s% max_T_shell_binding_energy
            case(h_max_T_lgP_thin_shell)
               val = s% max_T_lgP_thin_shell
            case(h_max_T_lgT)
               val = s% max_T_lgT
            case(h_max_T_lgP)
               val = s% max_T_lgP
            case(h_max_T_mass)
               val = s% max_T_mass
            case(h_max_T_radius)
               val = s% max_T_radius
            case(h_max_T_lgRho)
               val = s% max_T_lgRho
            case(h_max_T_L)
               val = s% max_T_L
            case(h_max_T_entropy)
               val = s% max_T_entropy
            case(h_max_T_eps_nuc)
               val = s% max_T_eps_nuc

            case(h_surface_optical_depth)
               val = s% tau_base*s% tau_factor
            case(h_log_surf_optical_depth)
               val = safe_log10(s% tau_base*s% tau_factor)

            case(h_log_surf_cell_opacity)
               val = safe_log10(s% opacity(1))
            case(h_log_surf_cell_density)
               val = s% lnd(1)/ln10
            case(h_surface_cell_temperature)
               val = s% T(1)
            case(h_log_surf_cell_temperature)
               val = s% lnT(1)/ln10
            case(h_surface_cell_entropy)
               val = s% entropy(1)
            case(h_log_surf_cell_P)
               val = s% lnP(1)/ln10
            case(h_log_surf_cell_pressure)
               val = s% lnP(1)/ln10
            case(h_log_surf_cell_z)
               val = 0
               if (s% net_iso(ih1) /= 0) val = val + s% xa(s% net_iso(ih1),1)
               if (s% net_iso(ih2) /= 0) val = val + s% xa(s% net_iso(ih2),1)
               if (s% net_iso(ihe3) /= 0) val = val + s% xa(s% net_iso(ihe3),1)
               if (s% net_iso(ihe4) /= 0) val = val + s% xa(s% net_iso(ihe4),1)
               val = safe_log10(1d0 - val)

            case(h_log_Ledd)
               Ledd = eval_Ledd(s, ierr)
               if (ierr /= 0 .or. Ledd <= 0d0) then
                  ierr = 0
               else
                  val = safe_log10(Ledd/Lsun)
               end if

            case(h_dt_div_dt_cell_collapse)
               val = s% dt/eval_min_cell_collapse_time(s, 2, nz, min_k, ierr)
            case(h_dt_cell_collapse)
               val = eval_min_cell_collapse_time(s, 2, nz, min_k, ierr)

            case(h_min_dr_div_cs_k)
               val = min_dr_div_cs(s,int_val)
               val = dble(int_val)
               is_int_val = .true.
            case(h_min_dr_div_cs)
               val = min_dr_div_cs(s,min_k)
            case(h_log_min_dr_div_cs)
               val = safe_log10(min_dr_div_cs(s,min_k))
            case(h_min_dr_div_cs_yr)
               val = min_dr_div_cs(s,min_k)/secyer
            case(h_log_min_dr_div_cs_yr)
               val = safe_log10(min_dr_div_cs(s,min_k)/secyer)
            case(h_dt_div_min_dr_div_cs)
               val = min_dr_div_cs(s,min_k)
               if (min_k <= 0) then
                  val = 1d99
               else
                  val = s% dt/val
               end if
            case(h_log_dt_div_min_dr_div_cs)
               val = min_dr_div_cs(s,min_k)
               if (min_k <= 0) then
                  val = 1d99
               else
                  val = safe_log10(s% dt/val)
               end if

            case(h_cz_bot_mass)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  if (k == nz) then
                     val = s% M_center/Msun
                  else
                     val = s% m(k)/Msun
                  end if
               end if
            case(h_cz_mass)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% m(k)/Msun
               end if
            case(h_cz_log_xmass)
               if (s% largest_conv_mixing_region == 0) then
                  val = -99
               else
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = safe_log10(s% xmstar*sum(s% dq(1:k-1)))
               end if
            case(h_cz_log_xmsun)
               if (s% largest_conv_mixing_region == 0) then
                  val = -99
               else
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = safe_log10(s% xmstar*sum(s% dq(1:k-1))/Msun)
               end if
            case(h_cz_xm)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% xmstar*sum(s% dq(1:k-1))/Msun
               end if
            case(h_cz_logT)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% lnT(k)/ln10
               end if
            case(h_cz_logRho)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% lnd(k)/ln10
               end if
            case(h_cz_logP)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% lnP(k)/ln10
               end if
            case(h_cz_log_column_depth)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = safe_log10(s% xmstar*sum(s% dq(1:k-1))/(pi4*s% r(k)*s% r(k)))
               end if
            case(h_cz_log_radial_depth)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = safe_log10(s% r(1) - s% r(k))
               end if
            case(h_cz_luminosity)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% L(k)/Lsun
               end if
            case(h_cz_log_tau)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = safe_log10(s% tau(k))
               end if
            case(h_cz_opacity)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% opacity(k)
               end if
            case(h_cz_log_eps_nuc)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = safe_log10(s% eps_nuc(k))
               end if
            case(h_cz_t_heat)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  if (s% eps_nuc(k) <= 0) then
                     val = 1d99
                  else
                     val = s% Cp(k)*s% T(k)/s% eps_nuc(k)
                  end if
               end if
            case(h_cz_eta)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% eta(k)
               end if
            case(h_cz_csound)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% csound(k)
               end if
            case(h_cz_scale_height)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% scale_height(k)
               end if
            case(h_cz_grav)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% grav(k)
               end if
            case(h_cz_bot_radius)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  if (k == nz) then
                     val = s% R_center/Rsun
                  else
                     val = s% R(k)/Rsun
                  end if
               end if
            case(h_cz_zone)
               if (s% largest_conv_mixing_region == 0) then
                  k = 0
               else
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
               end if
               int_val = k
               is_int_val = .true.
            case(h_cz_omega)
               if (s% largest_conv_mixing_region /= 0 .and. s% rotation_flag) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% omega(k)
               end if
            case(h_cz_omega_div_omega_crit)
               if (s% largest_conv_mixing_region /= 0 .and. s% rotation_flag) then
                  k = s% mixing_region_bottom(s% largest_conv_mixing_region)
                  val = s% omega(k)/omega_crit(s,k)
               end if

            case(h_cz_top_mass)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% m(k)/Msun
               end if
            case(h_cz_top_log_xmass)
               if (s% largest_conv_mixing_region == 0) then
                  val = -99
               else
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% xmstar*sum(s% dq(1:k-1)))
               end if
            case(h_cz_top_log_xmsun)
               if (s% largest_conv_mixing_region == 0) then
                  val = -99
               else
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% xmstar*sum(s% dq(1:k-1))/Msun)
               end if
            case(h_cz_top_xm)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% xmstar*sum(s% dq(1:k-1))/Msun
               end if
            case(h_cz_top_logT)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% lnT(k)/ln10
               end if
            case(h_cz_top_logRho)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% lnd(k)/ln10
               end if
            case(h_cz_top_logP)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% lnP(k)/ln10
               end if
            case(h_cz_top_log_column_depth)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% xmstar*sum(s% dq(1:k-1))/(pi4*s% r(k)*s% r(k)))
               end if
            case(h_cz_top_log_radial_depth)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% r(1) - s% r(k))
               end if
            case(h_cz_top_luminosity)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% L(k)/Lsun
               end if
            case(h_cz_top_log_tau)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% tau(k))
               end if
            case(h_cz_top_opacity)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% opacity(k)
               end if
            case(h_cz_top_log_eps_nuc)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% eps_nuc(k))
               end if
            case(h_cz_top_t_heat)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  if (s% eps_nuc(k) <= 0) then
                     val = 1d99
                  else
                     val = s% Cp(k)*s% T(k)/s% eps_nuc(k)
                  end if
               end if
            case(h_cz_top_eta)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% eta(k)
               end if
            case(h_cz_top_csound)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% csound(k)
               end if
            case(h_cz_top_scale_height)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% scale_height(k)
               end if
            case(h_cz_top_grav)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% grav(k)
               end if
            case(h_cz_top_radius)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% R(k)/Rsun
               end if

            case(h_mass_conv_core)
               val = s% mass_conv_core
            case(h_mass_semiconv_core)
               val = s% mass_semiconv_core

            case(h_cz_top_zone_logdq)
               if (s% largest_conv_mixing_region /= 0) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = safe_log10(s% dq(k))
               end if
            case(h_cz_top_zone)
               if (s% largest_conv_mixing_region == 0) then
                  k = 0
               else
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
               end if
               int_val = k
               is_int_val = .true.
            case(h_cz_top_omega)
               if (s% largest_conv_mixing_region /= 0 .and. s% rotation_flag) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% omega(k)
               end if
            case(h_cz_top_omega_div_omega_crit)
               if (s% largest_conv_mixing_region /= 0 .and. s% rotation_flag) then
                  k = s% mixing_region_top(s% largest_conv_mixing_region)
                  val = s% omega(k)/omega_crit(s,k)
               end if

            case(h_max_gradT_div_grada)
               val = 0
               do k = 2, nz
                  if (s% grada_face(k) == 0) cycle
                  if (s% gradT(k)/s% grada_face(k) > val) &
                     val = s% gradT(k)/s% grada_face(k)
               end do
            case(h_max_gradT_sub_grada)
               val = 0
               do k = 2, nz
                  if (s% gradT(k) - s% grada_face(k) > val) &
                     val = s% gradT(k) - s% grada_face(k)
               end do
            case(h_min_log_mlt_Gamma)
               val = 1d99
               do k = 2, nz
                  if (s% mlt_Gamma(k) > 0 .and. s% mlt_Gamma(k) < val) val = s% mlt_Gamma(k)
               end do
               val = safe_log10(val)

            case(h_max_conv_vel_div_csound)
               val = 0
               do k = 2, nz
                  if (s% q(k) > s% max_conv_vel_div_csound_maxq .or. s% csound(k) == 0) cycle 
                  if (s% conv_vel(k)/s% csound(k) > val) val = s% conv_vel(k)/s% csound(k)
               end do

            case(h_min_t_eddy)
               val = 1d99
               do k = 2, nz
                  if (s% conv_vel(k) <= 0) cycle
                  if (s% scale_height(k)/s% conv_vel(k) < val) &
                     val = s% scale_height(k)/s% conv_vel(k)
               end do

            case(h_elapsed_time)
               val = s% total_elapsed_time

            case(h_delta_nu)
               if (.not. s% get_delta_nu_from_scaled_solar) then
                  val = 1d6/(2*s% photosphere_acoustic_r) ! microHz
               else
                  val = &
                     s% delta_nu_sun*sqrt(s% star_mass)*pow3(s% Teff/s% Teff_sun) / &
                        pow(s% L_phot,0.75d0)
               end if
            case(h_delta_Pg)
               if (s% calculate_Brunt_N2) val = s% delta_Pg
            case(h_log_delta_Pg)
               if (s% calculate_Brunt_N2) val = safe_log10(s% delta_Pg)
            case(h_nu_max)
               val = s% nu_max
            case(h_nu_max_3_4th_div_delta_nu)
               val = pow(s% nu_max,0.75d0)/(1d6/(2*s% photosphere_acoustic_r))
            case(h_acoustic_cutoff)
               val = s% acoustic_cutoff
            case(h_acoustic_radius)
               val = s% photosphere_acoustic_r
            case(h_gs_per_delta_nu)
               if (s% calculate_Brunt_N2 .and. s% nu_max > 0 .and. s% delta_Pg >= 0) then
                  val = 1d6/(2*s% photosphere_acoustic_r) ! delta_nu
                  val = 1d6*val/(s% nu_max*s% nu_max*s% delta_Pg)
               end if
            case(h_ng_for_nu_max)
               if (s% calculate_Brunt_N2 .and. s% nu_max > 0 .and. s% delta_Pg >= 0) then
                  val = 1d6/(s% nu_max*s% delta_Pg)
               end if

            case(h_int_k_r_dr_nu_max_Sl1)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,1,1.0d0)
            case(h_int_k_r_dr_2pt0_nu_max_Sl1)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,1,2d0)
            case(h_int_k_r_dr_0pt5_nu_max_Sl1)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,1,0.5d0)
            case(h_int_k_r_dr_nu_max_Sl2)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,2,1.0d0)
            case(h_int_k_r_dr_2pt0_nu_max_Sl2)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,2,2d0)
            case(h_int_k_r_dr_0pt5_nu_max_Sl2)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,2,0.5d0)
            case(h_int_k_r_dr_nu_max_Sl3)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,3,1.0d0)
            case(h_int_k_r_dr_2pt0_nu_max_Sl3)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,3,2d0)
            case(h_int_k_r_dr_0pt5_nu_max_Sl3)
               if (s% calculate_Brunt_N2) val = get_int_k_r_dr(s,3,0.5d0)

            case (h_surface_extra_Pgas)
               val = s% surface_extra_Pgas

            case (h_min_L)
               val = minval(s% L(1:nz))/Lsun
            case (h_min_dL_dm)
               val = minval(s% dL_dm_expected(1:nz))
            case (h_min_dL_dm_m)
               val = s% m(minloc(s% dL_dm_expected(1:nz),dim=1))/Msun

            case (h_k_below_const_q)
               int_val = s% k_below_const_q
               is_int_val = .true.
            case (h_q_below_const_q)
               val = s% q(s% k_below_const_q)
            case (h_logxq_below_const_q)
               val = safe_log10(sum(s% dq(1:s% k_below_const_q-1)))

            case (h_k_const_mass)
               int_val = s% k_const_mass
               is_int_val = .true.
            case (h_q_const_mass)
               val = s% q(s% k_const_mass)
            case (h_logxq_const_mass)
               val = safe_log10(sum(s% dq(1:s% k_const_mass-1)))

            case (h_k_below_just_added)
               int_val = s% k_below_just_added
               is_int_val = .true.
            case (h_q_below_just_added)
               val = s% q(s% k_below_just_added)
            case (h_logxq_below_just_added)
               val = safe_log10(sum(s% dq(1:s% k_below_just_added-1)))

            case (h_k_for_test_CpT_absMdot_div_L)
               int_val = s% k_for_test_CpT_absMdot_div_L
               is_int_val = .true.
            case (h_q_for_test_CpT_absMdot_div_L)
               if (s% k_for_test_CpT_absMdot_div_L == nz) then
                  val = 0d0
               else
                  val = s% q(s% k_for_test_CpT_absMdot_div_L)
               end if
            case (h_logxq_for_test_CpT_absMdot_div_L)
               if (s% k_for_test_CpT_absMdot_div_L == nz) then
                  val = 0d0
               else
                  val = safe_log10(sum(s% dq(1:s% k_for_test_CpT_absMdot_div_L-1)))
               end if

            case (h_rotation_solver_steps)
               int_val = s% num_rotation_solver_steps
               is_int_val = .true.

            case (h_burn_solver_maxsteps)
               if (s% op_split_burn) &
                  int_val = maxval(s% burn_num_iters(1:s% nz))
               is_int_val = .true.

            case (h_diffusion_solver_steps)
               int_val = s% num_diffusion_solver_steps
               is_int_val = .true.

            case (h_diffusion_solver_iters)
               int_val = s% num_diffusion_solver_iters
               is_int_val = .true.

           case(h_tot_IE_div_IE_plus_KE)
               val = s% total_internal_energy_end / &
                     (s% total_internal_energy_end + s% total_radial_kinetic_energy_end)

           case(h_tot_E)
               val = s% total_energy_end
           case(h_log_tot_E)
               val = safe_log10(abs(s% total_energy_end))

           case(h_tot_KE)
               val = s% total_radial_kinetic_energy_end
           case(h_log_tot_KE)
               val = safe_log10(s% total_radial_kinetic_energy_end)

           case(h_tot_PE)
               val = s% total_gravitational_energy_end
           case(h_log_tot_PE)
               val = safe_log10(abs(s% total_gravitational_energy_end))

           case(h_tot_IE)
               val = s% total_internal_energy_end
           case(h_log_tot_IE)
               val = safe_log10(s% total_internal_energy_end)

           case(h_tot_Et)
               val = s% total_turbulent_energy_end
           case(h_log_tot_Et)
               val = safe_log10(s% total_turbulent_energy_end)
            
            case(h_num_hydro_merges)
               int_val = s% num_hydro_merges
               is_int_val = .true.
            case(h_num_hydro_splits)
               int_val = s% num_hydro_splits
               is_int_val = .true.

            case(h_RSP_DeltaR)
               if (s% RSP_flag) val = s% rsp_DeltaR
            case(h_RSP_DeltaMag)
               if (s% RSP_flag) val = s% rsp_DeltaMag
            case(h_RSP_GRPDV)
               if (s% RSP_flag) val = s% rsp_GRPDV
            case(h_RSP_GREKM)
               if (s% RSP_flag) val = s% rsp_GREKM
            case(h_RSP_GREKM_avg_abs)
               if (s% RSP_flag) val = s% rsp_GREKM_avg_abs

            case(h_RSP_phase)
               if (s% RSP_flag) val = (s% time - rsp_phase_time0())/s% RSP_period
            case(h_RSP_period_in_days)
               if (s% RSP_flag) val = s% RSP_period/(24*60*60) ! days
            case(h_RSP_num_periods)
               if (s% RSP_flag) int_val = s% RSP_num_periods
               is_int_val = .true.

            case(h_RSP_LINA_period_F_days)
               if (s% RSP_flag) val = s% rsp_LINA_periods(1)/86400.d0
            case(h_RSP_LINA_period_O1_days)
               if (s% RSP_flag) val = s% rsp_LINA_periods(2)/86400.d0
            case(h_RSP_LINA_period_O2_days)
               if (s% RSP_flag) val = s% rsp_LINA_periods(3)/86400.d0
            case(h_RSP_LINA_growth_rate_F)
               if (s% RSP_flag) val = s% rsp_LINA_growth_rates(1)
            case(h_RSP_LINA_growth_rate_O1)
               if (s% RSP_flag) val = s% rsp_LINA_growth_rates(2)
            case(h_RSP_LINA_growth_rate_O2)
               if (s% RSP_flag) val = s% rsp_LINA_growth_rates(3)

            case(h_grav_dark_L_polar) ! pole is at inclination = 0
               if(s% rotation_flag) then
                  w_div_w_Kep = if_rot(s% omega(1)*sqrt(pow3(s% r_equatorial(1))/(s% cgrav(1)*s% m(1))))
                  val = gravity_darkening_L_coeff(w_div_w_Kep, 0.0d0)*s% L_surf
               else
                  val = 0d0
               end if
            case(h_grav_dark_Teff_polar)
               if(s% rotation_flag) then
                  w_div_w_Kep = if_rot(s% omega(1)*sqrt(pow3(s% r_equatorial(1))/(s% cgrav(1)*s% m(1))))
                  val = gravity_darkening_Teff_coeff(w_div_w_Kep, 0.0d0)*s% Teff
               else
                  val = 0d0
               end if
            case(h_grav_dark_L_equatorial) ! equator is at inclination = pi/2
               if(s% rotation_flag) then
                  w_div_w_Kep = if_rot(s% omega(1)*sqrt(pow3(s% r_equatorial(1))/(s% cgrav(1)*s% m(1))))
                  val = gravity_darkening_L_coeff(w_div_w_Kep, 0.5d0*pi) * s% L_surf
               else
                  val = 0d0
               end if
            case(h_grav_dark_Teff_equatorial)
               if(s% rotation_flag) then
                  w_div_w_Kep = if_rot(s% omega(1)*sqrt(pow3(s% r_equatorial(1))/(s% cgrav(1)*s% m(1))))
                  val = gravity_darkening_Teff_coeff(w_div_w_Kep, 0.5d0*pi) * s% Teff
               else
                  val = 0d0
               end if

            case(h_apsidal_constant_k2)
               val = apsidal_constant(s, 2)

            ! following items correspond to names on terminal output lines

            case(h_lg_Lnuc)
               val = safe_log10(s% power_nuc_burn)
               
            case(h_H_rich)
               val = s% star_mass - max(s% he_core_mass, s% co_core_mass)

            case(h_N_cntr)
               val = s% center_n14

            case(h_lg_Lneu)
               val = safe_log10(abs(s% power_neutrinos))

            case(h_He_core)
               val = s% he_core_mass

            case(h_O_cntr)
               val = s% center_o16

            case(h_lg_Lphoto)
               val = safe_log10(abs(s% power_photo))

            case(h_CO_core)
               val = s% co_core_mass

            case(h_Fe_core)
               val = s% fe_core_mass

            case(h_Ne_cntr)
               val = s% center_ne20

            case(h_Mass)
               val = s% star_mass

            case(h_H_cntr)
               val = s% center_h1

            case(h_Si_cntr)
               val = s% center_si28

            case(h_lg_Mdot)
               val = safe_log10(abs(s% star_mdot))

            case(h_He_cntr)
               val = s% center_he3 + s% center_he4

            case(h_eta_cntr)
               val = s% eta(nz)

            case(h_gam_cntr)
               val = s% gam(nz)

            case(h_lg_Dsurf)
               val = s% lnd(1)/ln10

            case(h_C_cntr)
               val = s% center_c12


            case(h_zones)
               int_val = nz
               is_int_val = .true.

            case(h_retries)
               int_val = s% num_retries
               is_int_val = .true.
               
            case default
               ierr = -1

            end select

         end if


         contains


         real(dp) function max_eps_nuc_log_x(j)
            integer, intent(in) :: j
            real(dp) :: sum_x, sum_dq
            integer :: k
            max_eps_nuc_log_x = 0
            if (j == 0) return
            k = maxloc(s% eps_nuc(1:nz), dim=1)
            if (k < 1 .or. k > nz) return
            max_eps_nuc_log_x = s% xa(j,k)
         end function max_eps_nuc_log_x


         real(dp) function cz_top_max_log_x(j)
            integer, intent(in) :: j
            real(dp) :: sum_x, sum_dq
            integer :: k
            cz_top_max_log_x = 0
            if (s% largest_conv_mixing_region == 0) return
            k = s% mixing_region_top(s% largest_conv_mixing_region)
            if (j == 0 .or. k <= 1) return
            cz_top_max_log_x = s% xa(j,k-1)
         end function cz_top_max_log_x


         real(dp) function cz_max_log_x(j)
            integer, intent(in) :: j
            real(dp) :: sum_x, sum_dq
            integer :: k
            cz_max_log_x = 0
            if (s% largest_conv_mixing_region == 0) return
            k = s% mixing_region_bottom(s% largest_conv_mixing_region)
            if (j == 0 .or. k <= 1) return
            cz_max_log_x = s% xa(j,k-1)
         end function cz_max_log_x


         real(dp) function category_L(i)
            integer, intent(in) :: i
            if (i == 0) then
               category_L = 0
               return
            end if
            category_L = &
               safe_log10(dot_product(s% eps_nuc_categories(i,1:nz),s% dm(1:nz))/Lsun)
         end function category_L


         real(dp) function if_rot(v, alt)
            real(dp), intent(in) :: v
            real(dp), optional, intent(in) :: alt
            if (s% rotation_flag) then
               if_rot = v
            else
               if (present(alt)) then
                  if_rot = alt
               else
                  if_rot = 0
               end if
            endif
         end function if_rot

          
      end subroutine history_getval


      real(dp) function get_int_k_r_dr(s, el, nu_factor)
         use utils_lib, only: is_bad_num
         use chem_def, only: ih1
         type (star_info), pointer :: s
         integer, intent(in) :: el
         real(dp), intent(in) :: nu_factor

         real(dp) :: integral, cs2, r2, n2, sl2, omega2, &
            L2, kr2, dr, r0_outer, r0_inner, sl2_next, xh1
         integer :: k, k1, k_inner, k_outer, h1

         logical :: dbg

         include 'formats'

         dbg = .false. !(el == 1 .and. nu_factor == 1d0)


         get_int_k_r_dr = 0
         L2 = el*(el+1)
         omega2 = pow2(1d-6*2*pi*s% nu_max*nu_factor)

         ! k_inner and k_outer are bounds of evanescent region

         ! k_outer is outermost k where Sl2 <= omega2 at k-1 and Sl2 > omega2 at k
         ! 1st find outermost where Sl2 <= omega2

         h1 = s% net_iso(ih1)

         k1 = 0
         do k = 1, s% nz
            r2 = s% r(k)*s% r(k)
            cs2 = s% csound_face(k)*s% csound_face(k)
            sl2 = L2*cs2/r2
            if (sl2 <= omega2) then
               k1 = k; exit
            end if
         end do
         if (k1 == 0) return
         ! then find next k where Sl2 >= omega2
         k_outer = 0
         do k = k1+1, s% nz
            r2 = s% r(k)*s% r(k)
            cs2 = s% csound_face(k)*s% csound_face(k)
            sl2 = L2*cs2/r2
            if (sl2 > omega2) then
               k_outer = k; exit
            end if
         end do
         if (k_outer == 0) return

         ! k_inner is next k where N2 >= omega2 at k+1 and N2 < omega2 at k
         k_inner = 0
         do k = k_outer+1, s% nz
            if ( (s% brunt_N2(k) - s% brunt_N2_composition_term(k)) >= omega2) then
               ! Use the thermal component of the Brunt as starting point
               k_inner= k; exit
            end if
         end do
         if (k_inner == 0) return

         integral = 0
         do k = k_inner-1, k_outer, -1
            r2 = s% r(k)*s% r(k)
            cs2 = s% csound_face(k)*s% csound_face(k)
            n2 = s% brunt_N2(k)
            sl2 = L2*cs2/r2
            xh1 = s% xa(h1,k)
            dr = s% rmid(k-1) - s% rmid(k)
            kr2 = (1 - n2/omega2)*(1 - Sl2/omega2)*omega2/cs2
            if (kr2 < 0 .and. n2 < omega2 .and. omega2 < Sl2 ) &
               integral = integral + sqrt(-kr2)*dr
         end do


         if (integral == 0) return

         get_int_k_r_dr = integral

         if (dbg) write(*,3) 'r0 inner outer', &
            k_inner, k_outer, r0_inner/Rsun, r0_outer/Rsun, get_int_k_r_dr

         if (.not. is_bad_num(get_int_k_r_dr)) return

         write(*,2) 'el', el
         write(*,1) 'nu_factor', nu_factor
         write(*,1) 's% nu_max*nu_factor', s% nu_max*nu_factor
         write(*,1) 'log10 nu_max*nu_factor', log10(s% nu_max*nu_factor)
         write(*,1) 'Radius at k_inner', s% r(k_inner)
         write(*,1) 'Radius at k_outer', s% r(k_outer)


         write(*,1) 'get_int_k_r_dr', get_int_k_r_dr
         write(*,1) 'integral', integral
         write(*,2) 'k_inner', k_inner
         write(*,2) 'k_outer', k_outer

         stop 'get_int_k_r_dr'

      end function get_int_k_r_dr


      real(dp) function apsidal_constant(s,j)
        type (star_info), pointer :: s
        integer, intent(in) :: j
        real(dp) :: y, dy, rho_bar, dr_div_r, fprmid3
        integer :: k

        y = j - 2 ! value at r = 0
        do k = s% nz, 1, -1
           fprmid3 = pi4 * pow(s% rmid(k), 3)
           rho_bar = 3d0 * s% m(k) / fprmid3
           if (k == s% nz) then
              dr_div_r = s% r(k) / s% rmid(k) ! r(nz+1) would be zero
           else
              dr_div_r = (s% r(k) - s% r(k+1)) / s% rmid(k)
           end if
           dy = dr_div_r * (j*(j+1) - (6d0 * s% rho(k) / rho_bar) * (y + 1d0) - y * (y - 1d0))
           y = y + dy
        end do

        apsidal_constant = (j + 1d0 - y) / (2d0 * (j + y))

      end function apsidal_constant


      real(dp) function star_avg_x(s,j)
         type (star_info), pointer :: s
         integer, intent(in) :: j
         if (j == 0) then
            star_avg_x = 0
            return
         end if
         star_avg_x = dot_product(s% xa(j,1:s% nz),s% dq(1:s% nz))/sum(s% dq(1:s% nz))
      end function star_avg_x


      subroutine get_history_specs(s, num, names, specs, report)

         use utils_lib
         use utils_def

         type (star_info), pointer :: s
         integer, intent(in) :: num
         character (len=*), intent(in) :: names(:)
         integer, intent(out) :: specs(:)
         logical, intent(in) :: report

         integer :: i, ierr, n, j, iounit, t
         logical :: special_case
         character (len=strlen) :: buffer, string

         include 'formats'
         ierr = 0
         if (num <= 0) return
         iounit = -1
         specs(1:num) = 0
         do i = 1, num
            buffer = names(i)
            n = len_trim(buffer) + 1
            buffer(n:n) = ' '
            j = 0
            t = token(iounit, n, j, buffer, string)
            if (t /= name_token) then
               if (len_trim(names(i)) > 0 .and. report) &
                  write(*,*) 'bad value for name of history item ' // trim(names(i))
               specs(i) = -1
               ierr = 0
               cycle
            end if
            special_case = .false.
            specs(i) = do1_history_spec( &
               iounit, t, n, j, string, buffer, special_case, report, ierr)
            if (ierr /= 0 .or. special_case) then
               if (report) write(*,*) 'get_history_specs failed for ' // trim(names(i))
               specs(i) = -1
               ierr = 0
            end if
         end do

      end subroutine get_history_specs


      logical function get1_hist_value(s, name, val)
         ! includes other_history_columns from run_star_extras
         use utils_lib, only: integer_dict_lookup
         type (star_info), pointer :: s
         character (len=*) :: name
         real(dp), intent(out) :: val
         integer :: i, ierr, num_extra_cols, num_binary_cols
         character (len=80), pointer, dimension(:) :: &
            extra_col_names, binary_col_names
         real(dp), pointer, dimension(:) :: &
            extra_col_vals, binary_col_vals
         include 'formats'

         get1_hist_value = .false.
         call integer_dict_lookup(s% history_names_dict, name, i, ierr)
         if (ierr /= 0 .or. i <= 0) return ! didn't find it
         if (associated(s% pgstar_hist)) then
            if (associated(s% pgstar_hist% vals)) then
               if (size(s% pgstar_hist% vals,dim=1) >= i) then
                  val = s% pgstar_hist% vals(i)
                  get1_hist_value = .true.
                  return
               end if
            end if
         end if

         ! try extras
         if (associated(s% how_many_extra_history_columns) .and. &
             associated(s% data_for_extra_history_columns)) then
            num_extra_cols = s% how_many_extra_history_columns(s% id)
            if (num_extra_cols > 0) then
               allocate( &
                  extra_col_names(num_extra_cols), &
                  extra_col_vals(num_extra_cols), stat=ierr)
               call s% data_for_extra_history_columns( &
                  s% id, num_extra_cols, extra_col_names, extra_col_vals, ierr)
               do i=1,num_extra_cols
                  if (extra_col_names(i) == name) then
                     val = extra_col_vals(i)
                     get1_hist_value = .true.
                     exit
                  end if
               end do
               deallocate(extra_col_names, extra_col_vals)
               if (get1_hist_value) return
            end if
         end if

         ! try binary history
         num_binary_cols = s% how_many_binary_history_columns(s% binary_id)
         if (num_binary_cols > 0) then
            allocate( &
               binary_col_names(num_binary_cols), &
               binary_col_vals(num_binary_cols))
            call s% data_for_binary_history_columns( &
               s% binary_id, num_binary_cols, binary_col_names, binary_col_vals, ierr)
            if (ierr == 0) then
               do i=1,num_binary_cols
                  if (binary_col_names(i) == name) then
                     val = binary_col_vals(i)
                     get1_hist_value = .true.
                     exit
                  end if
               end do
            end if
            deallocate(binary_col_names, binary_col_vals)
            if (get1_hist_value) return
         end if

      end function get1_hist_value


      subroutine get_history_values(s, num, specs, &
            is_int_value, int_values, values, failed_to_find_value)
         ! note: this doesn't handle user-defined extra columns

         use utils_lib
         use utils_def

         type (star_info), pointer :: s
         integer, intent(in) :: num
         integer, intent(in) :: specs(:)
         logical, intent(out) :: is_int_value(:)
         integer, intent(out) :: int_values(:)
         real(dp), intent(inout) :: values(:)
         logical, intent(out) :: failed_to_find_value(:)

         integer :: i, c, int_val, ierr, n, t, j, iounit
         real(dp) :: val, epsnuc_out(12), v_surf, csound_surf, envelope_fraction_left
         logical :: is_int_val, special_case
         character (len=strlen) :: buffer, string

         include 'formats'
         ierr = 0
         if (num <= 0) return

         epsnuc_out(1:4) = s% burn_zone_mass(1:4,1)
         epsnuc_out(5:8) = s% burn_zone_mass(1:4,2)
         epsnuc_out(9:12) = s% burn_zone_mass(1:4,3)
         csound_surf = eval_csound(s,1,ierr)

         if (s% u_flag) then
            v_surf = s% u(1)
         else if (s% v_flag) then
            v_surf = s% v(1)
         else
            v_surf = s% r(1)*s% dlnR_dt(1)
         end if

         if (s% initial_mass > s% he_core_mass) then
            envelope_fraction_left = &
               (s% star_mass - s% he_core_mass)/(s% initial_mass - s% he_core_mass)
         else
            envelope_fraction_left = 1
         end if

         do i = 1, num
            failed_to_find_value(i) = .false.
            c = specs(i)
            if (c <= 0) then
               failed_to_find_value(i) = .true.
            else
               call history_getval( &
                  s, c, values(i), int_values(i), is_int_value(i), &
                  s% nz, v_surf, csound_surf, envelope_fraction_left, epsnuc_out, ierr)
               if (ierr /= 0) then
                  failed_to_find_value(i) = .true.
                  ierr = 0
               end if
            end if
         end do

      end subroutine get_history_values


      subroutine get_iso_val(s,str,val,ierr)
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         character (len=*), intent(in) :: str
         real(dp), intent(out) :: val
         integer, intent(out) :: ierr
         integer :: n, split, id, i
         ierr = 0
         val = 0
         n = len_trim(str)
         split = 0
         do i=1,n
            if (str(i:i) == ' ') then
               split = i
               exit
            end if
         end do
         if (split <= 1 .or. split >= n) then ! no interior space to split str
            ierr = -1
            return
         end if
         id = chem_get_iso_id(str(split+1:n))
         if (id <= 0) then ! not a valid iso name
            ierr = -1
            return
         end if
         i = s% net_iso(id)
         select case (str(1:split-1))
            case('center')
               val = center_avg_x(s,i)
            case('surface')
               val = surface_avg_x(s,i)
            case('average')
               val = star_avg_x(s,i)
            case('total_mass')
               val = star_avg_x(s,i)*s% xmstar/Msun
            case('log_total_mass')
               val = safe_log10(star_avg_x(s,i)*s% xmstar/Msun)
            case('log_average')
               val = safe_log10(star_avg_x(s,i))
            case('log_center')
               val = safe_log10(center_avg_x(s,i))
            case('log_surface')
               val = safe_log10(surface_avg_x(s,i))
            case default
               ierr = -1
         end select
      end subroutine get_iso_val

      end module history

