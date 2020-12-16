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

      module rates_initialize
      use const_def
      use math_lib
      use rates_def
      
      
      implicit none
            
      
      contains
      
      
      subroutine finish_rates_def_init
         use utils_lib, only: integer_dict_define, integer_dict_create_hash, integer_dict_size
         use reaclib_input, only: do_extract_rates
         use chem_def, only: nuclide_set, chem_isos, num_chem_isos, iso_name_length
         use chem_lib, only: generate_nuclide_set
         integer :: i, ierr
         character(len=iso_name_length), pointer :: names(:)
         type(nuclide_set), pointer :: set(:)
         integer, pointer :: chem_id(:) 
            ! will be allocated by extract_nuclides_from_chem_isos
         logical :: use_weaklib
         include 'formats.dek'
         
         ierr = 0
         call integer_dict_create_hash(reaction_names_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: rates_def_init failed in integer_dict_create_hash'
            return
         end if
         
         ! set up reaclib info
         allocate(set(num_chem_isos))
         
         call generate_nuclide_set(chem_isos% name, set)
         
         use_weaklib = .true.
         call do_extract_rates(set, chem_isos, reaclib_rates, use_weaklib, ierr)
         deallocate(set)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: extract_reaclib_rates failed in rates_def_init'
            return
         end if
         
      end subroutine finish_rates_def_init
         
         
      subroutine do_add_reaction_for_handle(reaction_handle, ierr)
         use reaclib_support, only: do_parse_reaction_handle
         character (len=*), intent(in) :: reaction_handle ! to be added
         integer, intent(out) :: ierr
         
         integer :: ir, num_in, num_out
         integer :: j, particles_in, particles_out
         logical :: already_defined
         integer :: iso_ids(max_num_reaction_inputs+max_num_reaction_outputs)
         integer :: cin(max_num_reaction_inputs), cout(max_num_reaction_outputs)
         character (len=16) :: op ! e.g., 'pg', 'wk', 'to', or ...
         
         logical, parameter :: weak = .false.
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         ierr = 0
         
         if (dbg) write(*,*) 'do_add_reaction_for_handle ' // trim(reaction_handle)
         
         call do_parse_reaction_handle( &
            reaction_handle, particles_in, particles_out, iso_ids, op, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'add_reaction_for_handle failed in reaclib_parse_handle ' //  &
               trim(reaction_handle)       
            return
         end if
         
         call alloc_reaction_ir(reaction_handle, ir, already_defined, ierr)
         if (already_defined) return
         if (ierr /= 0) return
      
         reaction_inputs(:,ir) = 0
         reaction_outputs(:,ir) = 0
         
         cin(:) = 1
         cout(:) = 1

         call setup(reaction_inputs(:,ir), num_in, cin, 0, particles_in)
         call setup(reaction_outputs(:,ir), num_out, cout, particles_in, particles_out)
         
         call set_reaction_info( &
               ir, num_in, num_out, particles_in, particles_out, weak, reaction_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in set_reaction_info for ' // trim(reaction_handle), ir
         end if

         if (dbg) write(*,*) 'done do_add_reaction_for_handle ' // trim(reaction_handle)
         
            
         contains
         

         subroutine setup(reaction_io, num, cnt, k, num_particles)
            integer :: reaction_io(:), num, cnt(:), k, num_particles
            integer :: j
            reaction_io(1) = 1
            reaction_io(2) = iso_ids(k+1)
            num = 1
            cnt(num) = 1
            do j=2,num_particles
               if (iso_ids(k+j) == iso_ids(k+j-1)) then
                  cnt(num) = cnt(num) + 1
               else
                  num = num + 1
                  reaction_io(2*num) = iso_ids(k+j)
                  cnt(num) = 1
               end if
               reaction_io(2*num-1) = cnt(num)
            end do
         end subroutine setup


      end subroutine do_add_reaction_for_handle
         
         
      subroutine do_add_reaction_from_reaclib(reaction_handle, reverse_handle, indx, ierr)
         
         character (len=*), intent(in) :: reaction_handle ! to be added
         character (len=*), intent(in) :: reverse_handle ! = '' if not a reverse
         integer, intent(in) :: indx ! index in reaclib rates
         integer, intent(out) :: ierr
         
         integer :: i, ir, chapter, num_in, num_out, iso_in, iso_out
         integer :: j, weak_j, particles_in, particles_out
         logical :: weak, reverse, already_defined
         integer :: cin(max_num_reaction_inputs), cout(max_num_reaction_outputs)
         type (reaction_data), pointer :: r
         
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         ierr = 0
         i = indx
         reverse = (len_trim(reverse_handle) > 0)
         r => reaclib_rates
         
         cin(:) = 1
         cout(:) = 1
         
         if (dbg) write(*,'(a, 2x, i5)') 'do_add_reaction_from_reaclib ' // trim(reaction_handle), i
         
         if (reverse) then
            if (reverse_handle /= r% reaction_handle(i)) then
               write(*,'(a)') trim(reverse_handle) // ' ' // trim(r% reaction_handle(i))
               write(*,'(a,3x,i8)') 'bad reverse_handle for add_reaction_from_reaclib', indx
               ierr = -1
               return
            end if
         else
            if (reaction_handle /= r% reaction_handle(i)) then
               write(*,'(a)') trim(reaction_handle) // ' ' // trim(r% reaction_handle(i))
               write(*,'(a,3x,i8)') 'bad reaction_handle for add_reaction_from_reaclib', indx
               ierr = -1
               return
            end if
         end if
         
         chapter = r% chapter(i)
         weak = (adjustl(r% reaction_flag(i)) == 'w')
         
         call alloc_reaction_ir(reaction_handle, ir, already_defined, ierr)
         if (already_defined) return
         if (ierr /= 0) return
      
         reaction_inputs(:,ir) = 0
         reaction_outputs(:,ir) = 0

         if (.not. reverse) then
            particles_in = Nin(chapter)
            particles_out = Nout(chapter)
            call setup(reaction_inputs(:,ir), num_in, cin, 0, particles_in)
            call setup(reaction_outputs(:,ir), num_out, cout, particles_in, particles_out)
         else
            particles_out = Nin(chapter)
            particles_in = Nout(chapter)
            call setup(reaction_inputs(:,ir), num_in, cin, particles_out, particles_in)
            call setup(reaction_outputs(:,ir), num_out, cout, 0, particles_out)
         end if
         
         call set_reaction_info( &
               ir, num_in, num_out, particles_in, particles_out, weak, reaction_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in set_reaction_info for ' // trim(reaction_handle), ir
         end if

         if (dbg) write(*,*) 'done do_add_reaction_from_reaclib'
         
            
         contains
         

         subroutine setup(reaction_test, num, cnt, k, num_particles)
            integer :: reaction_test(:), num, cnt(:), k, num_particles
            integer :: j
            reaction_test(1) = 1
            reaction_test(2) = r% pspecies(k+1,i)
            num = 1
            cnt(num) = 1
            do j=2,num_particles
               if (r% pspecies(k+j,i) == r% pspecies(k+j-1,i)) then
                  cnt(num) = cnt(num) + 1
               else
                  num = num + 1
                  reaction_test(2*num) = r% pspecies(k+j,i)
                  cnt(num) = 1
               end if
               reaction_test(2*num-1) = cnt(num)
            end do
         end subroutine setup


      end subroutine do_add_reaction_from_reaclib
      
      
      subroutine alloc_reaction_ir(reaction_handle, ir, already_defined, ierr)
         character (len=*), intent(in) :: reaction_handle
         integer, intent(out) :: ir
         logical, intent(out) :: already_defined
         integer, intent(out) :: ierr         
         logical, parameter :: dbg = .false.         
         include 'formats.dek'
         ierr = 0        
         ir = get_rates_reaction_id(reaction_handle)
         if (ir > 0) then
            already_defined = .true.
            return
         end if
!$omp critical (lock_alloc_reaction_ir)
         ir = get_rates_reaction_id(reaction_handle)
         if (ir > 0) then
            already_defined = .true.
         else
            already_defined = .false.
            if (dbg) write(*,2) 'call increase_num_reactions', rates_reaction_id_max
            call increase_num_reactions(ierr)
            if (ierr == 0) then
               ir = rates_reaction_id_max
               if (dbg) write(*,2) 'after increase_num_reactions', rates_reaction_id_max
            end if
         end if
!$omp end critical (lock_alloc_reaction_ir)
      end subroutine alloc_reaction_ir
      
      
      subroutine set_reaction_info( &
            ir, num_in, num_out, particles_in, particles_out, weak, reaction_handle, ierr)
         use chem_def
         use utils_lib, only: integer_dict_define
         integer, intent(in) :: ir, num_in, num_out, particles_in, particles_out
         logical, intent(in) :: weak
         character (len=*), intent(in) :: reaction_handle
         integer, intent(out) :: ierr
         
         integer :: j, iso_in, iso_out, weak_j, cin1
         
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         ierr = 0
         
         reaction_Name(ir) = reaction_handle
        
         if (dbg) write(*,*) 'call get_Qtotal'
         std_reaction_Qs(ir) = get_Qtotal(ir)
         std_reaction_neuQs(ir) = 0d0
         weak_lowT_rate(ir) = -1d99
         
         iso_in = reaction_inputs(num_in*2,ir)
         cin1 = reaction_inputs(num_in*2-1,ir)
         iso_out = reaction_outputs(num_out*2,ir)

         if (iso_in < 0 .or. iso_in > num_chem_isos) then
            write(*,3) 'bad iso_in', iso_in, num_in
            write(*,*) trim(reaction_handle)
            do j = 1, num_in
               write(*,3) 'reaction input', j, reaction_inputs(2*j,ir)
            end do
            stop 'set_reaction_info'
         end if
         
         if (iso_out < 0 .or. iso_out > num_chem_isos) then
            write(*,2) 'bad iso_out', iso_out
            write(*,2) 'num_out', num_out
            write(*,2) 'num_in', num_in
            write(*,*) trim(reaction_handle)
            do j = 1, num_out
               write(*,3) 'reaction output', j, reaction_outputs(2*j,ir)
            end do
            stop 'set_reaction_info'
         end if
      
         if (weak) then
            weak_reaction_info(1,ir) = iso_in
            weak_reaction_info(2,ir) = iso_out
            weak_j = do_get_weak_info_list_id(chem_isos% name(iso_in),chem_isos% name(iso_out))
            if (weak_j > 0) std_reaction_neuQs(ir) = weak_info_list_Qneu(weak_j)
            call set_weak_lowT_rate(ir, ierr)
            if (ierr /= 0) return
         else
            weak_reaction_info(1:2,ir) = 0
         end if

         reaction_ye_rho_exponents(1,ir) = 0 ! 1 for electron captures, 0 for rest.
      
         if (particles_in > 1) then
            reaction_screening_info(1,ir) = reaction_inputs(2,ir)
            if (cin1 > 1) then
               reaction_screening_info(2,ir) = reaction_screening_info(1,ir)
            else
               reaction_screening_info(2,ir) = reaction_inputs(4,ir)
            end if
            reaction_ye_rho_exponents(2,ir) = particles_in - 1
         else
            reaction_screening_info(1:2,ir) = 0
            reaction_ye_rho_exponents(2,ir) = 0
         end if
         reaction_screening_info(3,ir) = 0
         reaction_categories(ir) = iother
         if (is_pp_reaction(reaction_handle)) then
            reaction_categories(ir) = ipp
         else if (is_cno_reaction(reaction_handle)) then
            reaction_categories(ir) = icno
         else if (num_in == 1 .and. cin1 == 2) then
            if (iso_in == ic12) then
               reaction_categories(ir) = icc
            else if (iso_in == io16) then
               reaction_categories(ir) = ioo
            end if
         else if (num_in > 1) then
            select case(chem_isos% Z(iso_in))
               case (6)
                  reaction_categories(ir) = i_burn_c
               case (7)
                  reaction_categories(ir) = i_burn_n
               case (8)
                  reaction_categories(ir) = i_burn_o
               case (10)
                  reaction_categories(ir) = i_burn_ne
               case (11)
                  reaction_categories(ir) = i_burn_na
               case (12)
                  reaction_categories(ir) = i_burn_mg
               case (14)
                  reaction_categories(ir) = i_burn_si
               case (16)
                  reaction_categories(ir) = i_burn_s
               case (18)
                  reaction_categories(ir) = i_burn_ar
               case (20)
                  reaction_categories(ir) = i_burn_ca
               case (22)
                  reaction_categories(ir) = i_burn_ti
               case (24)
                  reaction_categories(ir) = i_burn_cr
               case (26)
                  reaction_categories(ir) = i_burn_fe
            end select
         else if (particles_in == 1 .and. particles_out == 2 .and. .not. weak) then
            reaction_categories(ir) = iphoto
         end if 

         reaction_Info(ir) = reaction_handle
         
         if (dbg) write(*,'(a,3x,i5)') 'call integer_dict_define ' // trim(reaction_handle), ir

         call integer_dict_define(reaction_names_dict, reaction_handle, ir, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: set_reaction_info failed in integer_dict_define'
            return
         end if         
      
      end subroutine set_reaction_info

      
      subroutine set_weak_lowT_rate(ir, ierr)
         use chem_def
         use utils_lib, only: integer_dict_define
         use reaclib_eval, only: do_reaclib_indices_for_reaction
         use ratelib, only: reaclib_rate_and_dlnT
         use load_weak, only: do_get_weak_rate_id

         integer, intent(in) :: ir
         integer, intent(out) :: ierr
         integer :: i, lo, hi, weak_reaclib_id_i
         real(dp) :: half_life, lambda, dlambda_dlnT, rlambda, drlambda_dlnT
         
         include 'formats'
         
         weak_reaclib_id_i = 0
         
         if (weak_reaction_info(1,ir) == 0 .or. weak_reaction_info(2,ir) == 0) return

         call do_reaclib_indices_for_reaction( &
            reaction_Name(ir), reaclib_rates, lo, hi, ierr)
         if (ierr /= 0) then ! not in reaclib
            ierr = 0
            i = do_get_weak_info_list_id( &
               chem_isos% name(weak_reaction_info(1,ir)), &
               chem_isos% name(weak_reaction_info(2,ir)))
            if (i > 0) then 
               half_life = weak_info_list_halflife(i)
               if (half_life > 0d0) then
                  weak_lowT_rate(ir) = ln2/half_life
                  return
               end if
            end if
            weak_lowT_rate(ir) = -1d-99
            return
         end if

         ! evaluate rates at T= 10^7 K (i.e. T9 = 0.01)
         call reaclib_rate_and_dlnT( &
            lo, hi, reaction_Name(ir), 0.01_dp, &
            lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
         if (ierr /= 0) then
            write(*,2) 'set_reaction_info failed in reaclib_rate_and_dlnT ' // &
               trim(reaction_Name(ir)), ir
            call mesa_error(__FILE__,__LINE__)
         end if
         if (lo <= reaclib_rates% num_from_weaklib) then
            if (.not. reaclib_rates% also_in_reaclib(lo)) lambda = -1
         end if
         weak_lowT_rate(ir) = lambda

         weak_reaclib_id_i = do_get_weak_rate_id( &
            chem_isos% name(weak_reaction_info(1,ir)), &
            chem_isos% name(weak_reaction_info(2,ir)))
         if (weak_reaclib_id_i == 0) return
         weak_reaclib_id(weak_reaclib_id_i) = ir
             
      end subroutine set_weak_lowT_rate
      
      
      logical function is_pp_reaction(reaction_handle)
         character (len=*), intent(in) :: reaction_handle
         is_pp_reaction =  &
            reaction_handle == 'r_h1_h1_wk_h2' .or.  &
            reaction_handle == 'r_h1_h1_ec_h2' .or.  &
            reaction_handle == 'r_h2_pg_he3' .or.  &
            reaction_handle == 'r_he3_he3_to_h1_h1_he4' .or.  &
            reaction_handle == 'r_he3_ag_be7' .or.  &
            reaction_handle == 'r_be7_wk_li7' .or.  &
            reaction_handle == 'r_h1_li7_to_he4_he4' .or.  &
            reaction_handle == 'r_li7_pa_he4' .or.  &
            reaction_handle == 'r_be7_pg_b8' .or.  &
            reaction_handle == 'r_b8_wk_he4_he4'                                             
      end function is_pp_reaction
      
      
      logical function is_cno_reaction(reaction_handle)
         character (len=*), intent(in) :: reaction_handle
         is_cno_reaction =  &
            reaction_handle == 'r_c12_pg_n13' .or.  &
            reaction_handle == 'r_c13_pg_n14' .or.  &
            reaction_handle == 'r_n13_wk_c13' .or.  &
            reaction_handle == 'r_n13_pg_o14' .or.  &
            reaction_handle == 'r_n14_pg_o15' .or.  &
            reaction_handle == 'r_n15_pg_o16' .or.  &
            reaction_handle == 'r_n15_pa_c12' .or.  &
            reaction_handle == 'r_o14_wk_n14' .or.  &
            reaction_handle == 'r_o14_ap_f17' .or.  &
            reaction_handle == 'r_o15_wk_n15' .or.  &
            reaction_handle == 'r_o16_pg_f17' .or.  &
            reaction_handle == 'r_o17_pg_f18' .or.  &
            reaction_handle == 'r_o17_pa_n14' .or.  &
            reaction_handle == 'r_o18_pg_f19' .or.  &
            reaction_handle == 'r_o18_pa_n15' .or.  &
            reaction_handle == 'r_f17_wk_o17' .or.  &
            reaction_handle == 'r_f17_pg_ne18' .or.  &
            reaction_handle == 'r_f18_pa_o15' .or.  &
            reaction_handle == 'r_f18_wk_o18' .or.  &
            reaction_handle == 'r_f19_pa_o16' .or.  &
            reaction_handle == 'r_ne18_wk_f18'
      end function is_cno_reaction
      
      
      subroutine free_raw_rates_records
         type (rate_table_info), pointer :: ri
         integer :: i
         if (ASSOCIATED(raw_rates_records)) then
            do i = 1, rates_reaction_id_max
               ri => raw_rates_records(i)
               if (ASSOCIATED(ri% T8s)) deallocate(ri% T8s)
               if (ASSOCIATED(ri% f1)) deallocate(ri% f1)
            end do
            deallocate(raw_rates_records)
         end if
      end subroutine free_raw_rates_records
      
      
      subroutine init_raw_rates_records(ierr)
         use utils_lib
         use utils_def
         integer, intent(out) :: ierr
         
         type (rate_table_info), pointer :: ri
         integer :: i, iounit, n, t, ir
         character (len=256) :: line, dir, rate_name, rate_fname, filename, message
         character (len=256) :: buffer, string
         logical :: okay
         
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         if (dbg) write(*,*) 'init_raw_rates_records'
         
         ierr = 0
         
         ! first try local rate_tables_dir
         dir = rates_table_dir
         filename = trim(dir) // '/rate_list.txt'
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in rates_dir
            dir = trim(rates_dir) // '/rate_tables'
            filename = trim(dir) // '/rate_list.txt'
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open rates list file ' // trim(filename)
               return
            end if
         end if
         rates_table_dir = dir
         
         n = 0
         i = 0

         if (dbg) write(*,*) 'read rate list file ' // trim(filename)

         rate_loop: do
            t = token(iounit, n, i, buffer, rate_name)
            if (t == eof_token) exit
            if (t /= name_token) then
               call error; return
            end if
            if (dbg) write(*,*) 'use rate table from file for ', trim(rate_name)
            ir = lookup_rate_name(rate_name)
            if (ir <= 0) then
               call do_add_reaction_for_handle(rate_name, ierr)
               if (ierr == 0) ir = lookup_rate_name(rate_name)
               if (ierr /= 0 .or. ir <= 0) then
                  write(*,*) 'invalid rate file ' // trim(rate_name)
                  call error; return
               end if
            end if
            t = token(iounit, n, i, buffer, rate_fname)
            if (t /= string_token) then
               call error; return
            end if
            if (dbg) write(*,*) 'rate_fname ', trim(rate_fname)
            ri => raw_rates_records(ir)
            ri% use_rate_table = .true.
            ri% need_to_read = .true.
            ri% rate_fname = trim(rate_fname)
            if (dbg) write(*,*) 'done'
         end do rate_loop
         
         close(iounit)
         
         if (dbg) stop 'read rates'
         !call check
         
         
         
         contains
         
         
         subroutine check
            ! check that there are cases for all of the rates
            use ratelib, only: tfactors
            use raw_rates, only: set_raw_rate
            real(dp) :: logT, temp, raw_rate
            integer :: which_rate, i, ierr
            type (T_Factors), pointer :: tf
         
            allocate(tf)
            logT = 8
            temp = exp10(logT)
            call tfactors(tf, logT, temp)
            ierr = 0
            which_rate = 1
            okay = .true.
            do i = 1, rates_reaction_id_max
               call set_raw_rate(i, which_rate, temp, tf, raw_rate, ierr)
               if (ierr /= 0) then
                  write(*,'(a)') 'set_raw_rate failed for ' // reaction_Name(i)
                  okay = .false.
                  ierr = 0
               end if
            end do
            deallocate(tf)
            if (.not. okay) stop 'init_raw_rates_records'
         end subroutine check
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = .false.
            if (ierr == 0) return
            if (len_trim(str) > 0) then
               write(*,*) trim(str) // ' failed in reading ' // trim(filename)
            else ! non-fatal error, so just quietly stop reading
               ierr = 0
            end if
            failed = .true.
            return
         end function failed
         
         
         subroutine error
            ierr = -1
            close(iounit)
         end subroutine error
      
      
      end subroutine init_raw_rates_records

      
      subroutine read_reaction_parameters(reactionlist_filename, ierr)
         use utils_lib
         use chem_lib
         use chem_def, only: chem_isos, category_id, category_Name
         use const_def, only: mesa_data_dir
         character (len=*), intent(in) :: reactionlist_filename
         integer, intent(out) :: ierr
         
         character (len=256) :: line, filename, rname, cname, iname, str, message
         integer :: iounit, len, i, j, jj, k, cnt, ir, ic, ii, ye, rho, n, num_reactions, new_max
         character (len=maxlen_reaction_Info) :: info
         real(dp) :: Q
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'

         ierr = 0
         
         call alloc_and_init_reaction_parameters(ierr)
         if (ierr /= 0) return
         
         ! first try the reaction_filename alone
         filename = trim(reactionlist_filename)
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in rates_data
            filename = trim(mesa_data_dir) // '/rates_data/' // trim(reactionlist_filename)
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*)  &
                     'failed to open reaction_parameters file ' // trim(reactionlist_filename)
               return
            end if
         end if

         num_reactions = 0
         do cnt = 1, rates_reaction_id_max*10 ! will stop when reach end of file
            if (dbg) write(*, *) 'cnt', cnt

            read(iounit,'(a)',iostat=ierr) line
            if (ierr /= 0) then
               if (dbg) write(*,*) 'reached end of file'
               exit
            end if
            len = len_trim(line)
            if (len == 0) then
               if (dbg) write(*,*) '(len == 0)'
               cycle
            end if
            if (line(1:1) == '!') then
               if (dbg) write(*,*) '(line(1:1) == !)'
               cycle
            end if
            
            i = 1; j = 35
            rname = line(i:j)
            
            if (dbg) write(*,*) trim(rname)
            
            if (line(i:i+1) == 'r ') then
               call increase_num_reactions(ierr)
               if (ierr /= 0) then
                  write(*,*) 'FATAL ERROR: rates_def_init failed in increase_num_reactions'
                  return
               end if
               ir = rates_reaction_id_max
            else
               ir = get_rates_reaction_id(rname)
            end if
            if (dbg) write(*,*)   'name: ' // trim(rname), ir
            if (ir == 0) then
               call increase_num_reactions(ierr)
               if (ierr /= 0) then
                  write(*,*) 'FATAL ERROR: rates_def_init failed in increase_num_reactions'
                  return
               end if
               ir = rates_reaction_id_max
               if (dbg) write(*,*) 'size(reaction_Name,dim=1), ir', size(reaction_Name,dim=1), ir
               reaction_Name(ir) = rname               
               call integer_dict_define(reaction_names_dict, reaction_Name(ir), ir, ierr)
               if (ierr /= 0) then
                  write(*,*) 'FATAL ERROR: rates_def_init failed in integer_dict_define'
                  call mesa_error(__FILE__,__LINE__)
               end if
            end if

            i = 36; j = 70           
            call read_inputs
            if (ierr /= 0) return

            i = 74; j = 108
            call read_outputs
            if (ierr /= 0) return

            i = 110; j = 127
            call read_Q
            if (ierr /= 0) return
      
            i = 128; j = 143
            call read_Qneu
            if (ierr /= 0) return
         
            i = 144; j = 149
            call read_ye_rho_exponent1
            if (ierr /= 0) return
         
            i = 150; j = 155
            call read_ye_rho_exponent2
            if (ierr /= 0) return
         
            i = 156; j = 160
            call read_screening_info(1)
            if (ierr /= 0) return
         
            i = 162; j = 166
            call read_screening_info(2)
            if (ierr /= 0) return
         
            i = 168; j = 172
            call read_screening_info(3)
            if (ierr /= 0) return
         
            i = 174; j = 178
            call read_weak_info(1)
            if (ierr /= 0) return
         
            i = 180; j = 184
            call read_weak_info(2)
            if (ierr /= 0) return
               
            if (std_reaction_neuQs(ir) > 0) then
               weak_reaction_info(1,ir) = reaction_inputs(2,ir)
               weak_reaction_info(2,ir) = reaction_outputs(2,ir)
            end if
               
            if (std_reaction_neuQs(ir) == 0 .and. &
                  weak_reaction_info(1,ir) > 0 .and. weak_reaction_info(2,ir) > 0) then
               j = do_get_weak_info_list_id(  &
                     chem_isos% name(weak_reaction_info(1,ir)),  &
                     chem_isos% name(weak_reaction_info(2,ir)))
               if (j > 0) std_reaction_neuQs(ir) = weak_info_list_Qneu(j)
               if (ierr /= 0) return
            end if
            
            if (std_reaction_neuQs(ir) > 0) then
               call set_weak_lowT_rate(ir, ierr)
               if (ierr /= 0) return
            end if

            i = 190; j = 203
            call read_category_id
            if (ierr /= 0) return
         
            i = 208
            call read_reaction_Info
            if (ierr /= 0) return
            
            if (dbg) write(*,*)
            
            num_reactions = cnt
            
         end do
         
         if (dbg) write(*,*) 'num_reactions', num_reactions
         
         ierr = 0
         close(iounit)
         
         num_reactions = rates_reaction_id_max
         
         call check_std_reaction_Qs
         call check_std_reaction_neuQs
         call check_reaction_categories
         call check_reaction_info
         
         if (dbg) stop 'read_reaction_parameters'
         
         
         contains
            
            
         subroutine read_inputs
            if (dbg) write(*,*) '  inputs <' // line(i:j) // '>',i,j

            jj = j
            do k = 1, 2*max_num_reaction_inputs-1, 2
               j = i; j = i+2
               n = read_int()
               if (n == 0) exit
               if (dbg) write(*,*) 'n <' // line(i:j) // '>', i, j
               reaction_inputs(k,ir) = n
 
! fxt
              i = j+1; j = i+7

!               i = j+1; j = i+4

               if (dbg) write(*,*) 'iso <' // line(i:j) // '>', i, j
               ii = read_iso()
               if (ii <= 0) then
                  write(*,'(a)') 'bad input iso in reaction_parameters file <' // line(i:j) // '>'
                  write(*,'(a)') trim(line)
                  stop 'read_reaction_parameters'
                  ierr = -1
                  return
               end if
               i = j+2
               reaction_inputs(k+1,ir) = ii
               if (dbg) write(*,*) 'in', n, trim(chem_isos% name(ii)), i
            end do
         end subroutine read_inputs

            
         subroutine read_outputs
            if (dbg) write(*,*) ' outputs: ' // line(i:j)

            jj = j
            do k = 1, 2*max_num_reaction_outputs-1, 2
               j = i; j = i+2
               n = read_int()
               if (n == 0) exit
               if (dbg) write(*,*) 'n <' // line(i:j) // '>'
               reaction_outputs(k,ir) = n

! fxt
               i = j+1; j = i+7

!               i = j+1; j = i+4

               if (dbg) write(*,*) 'iso <' // line(i:j) // '>'
               ii = read_iso()
               if (ii <= 0) then
                  write(*,'(a)') 'bad output iso in reaction_parameters file <' // line(i:j) // '>'
                  write(*,'(a)') trim(line)
                  stop 'read_reaction_parameters'
                  ierr = -1
                  return
               end if
               i = j+2
               reaction_outputs(k+1,ir) = ii
               if (dbg) write(*,*) 'out', n, trim(chem_isos% name(ii))
            end do
         end subroutine read_outputs


         subroutine read_Q
            if (missing_dbl()) then ! use standard Q from chem_lib
               std_reaction_Qs(ir) = get_Qtotal(ir)
            else
               std_reaction_Qs(ir) = read_dbl()
               if (dbg) write(*,*) 'std_reaction_Qs(ir)', std_reaction_Qs(ir)
            end if
         end subroutine read_Q


         subroutine read_Qneu
            include 'formats'
            std_reaction_neuQs(ir) = read_dbl()
            if (dbg) write(*,1) 'std_reaction_neuQs(ir)', std_reaction_neuQs(ir)
         end subroutine read_Qneu


         subroutine read_ye_rho_exponent1
            reaction_ye_rho_exponents(1,ir) = read_int()
            if (dbg) write(*,*)  'ye', reaction_ye_rho_exponents(1,ir)
            if (reaction_ye_rho_exponents(1,ir) > 2) then
               write(*,'(a)') 'ERROR: must revise rates_support for large ye exponent'
               write(*,'(a)') trim(line)
               stop 'read_reaction_parameters'
            end if
         end subroutine read_ye_rho_exponent1


         subroutine read_ye_rho_exponent2
            ii = read_int()
            if (ii > 0) ii = ii-1
            reaction_ye_rho_exponents(2,ir) = ii
            if (ii > 4) then
               write(*,'(a)') 'ERROR: must revise rates_support for large rho exponent'
               write(*,'(a)') trim(line)
               stop 'read_reaction_parameters'
            end if
            if (dbg) write(*,*)  'rho', reaction_ye_rho_exponents(2,ir)+1
         end subroutine read_ye_rho_exponent2


         subroutine read_screening_info(which)
            integer, intent(in) :: which
            logical :: empty
            integer :: jj
            ii = read_iso()
            reaction_screening_info(which,ir) = ii
            
          !Hack to get around a bug in ifort 17,18, which returns len_trim(line(i:j)) < 0
            empty=.true.
            do jj=i,j
               if(len_trim(line(jj:jj)) /= 0) empty=.false.
            end do
            
            if (ii <= 0 .and. .not. empty) then
               write(*,'(a)') 'bad iso name for screening in reaction_parameters file <' // line(i:j) // '>'
               write(*,'(a)') trim(line)
               call mesa_error(__FILE__,__LINE__,'read_screening_info')
               ierr = -1
               return
            end if
            if (dbg .and. ii > 0)  &
                  write(*,*) 'screen: ' // trim(chem_isos% name(ii)) // ' ' // line(i:j), ii
         end subroutine read_screening_info


         subroutine read_weak_info(which)
            integer, intent(in) :: which
            logical :: empty
            integer :: jj
            
            ii = read_iso()
            weak_reaction_info(which,ir) = ii
            
            !Hack to get around a bug in ifort 7, 18 which returns len_trim(line(i:j)) < 0
            empty=.true.
            do jj=i,j
               if(len_trim(line(jj:jj)) /= 0) empty=.false.
            end do
            
            if (ii <= 0 .and. .not. empty) then
               write(*,'(a)') 'bad iso name for weak in reaction_parameters file <' // line(i:j) // '>'
               write(*,'(a)') trim(line)
               call mesa_error(__FILE__,__LINE__,'read_reaction_parameters len')
               ierr = -1
               return
            end if
            if (ii > 0) then
               write(*,*)  'weak: ' // trim(chem_isos% name(ii)) // ' ' // line(i:j), ii
               write(*,'(a)') trim(line)
               write(*,'(a)') 'DO NOT USE reactions.list FOR WEAK ISOS; just give Qneu'
               call mesa_error(__FILE__,__LINE__,'read_reaction_parameters weak')
               ierr = -1
               return
            end if
         end subroutine read_weak_info


         subroutine read_category_id
            cname = line(i:j)
            ic = category_id(cname)
            if (ic == 0) then
               write(*,*) 'bad category name in reaction_parameters file <' // line(i:j) // '>'
               write(*,'(a)') trim(line)
               call mesa_error(__FILE__,__LINE__,'read_category_id')
               ierr = -1
               return
            end if
            reaction_categories(ir) = ic
            if (dbg) write(*,*) 'category: ' // trim(cname), ic, trim(category_name(ic))
         end subroutine read_category_id


         subroutine read_reaction_Info
            ii = min(maxlen_reaction_Info, len - i + 1)
            do j = 1, maxlen_reaction_Info
               if (j <= ii) then
                  info(j:j) = line(i:i)
                  i = i+1
               else
                  info(j:j) = ' ' 
               end if
            end do
            reaction_Info(ir) = info
            if (dbg) write(*,*)  'info: ' // trim(reaction_Info(ir))
         end subroutine read_reaction_Info
         
         
         integer function read_iso()
            use chem_def, only: iso_name_length
            character (len=64) :: str
            integer :: ierr
            str = line(i:j)
            read_iso = chem_get_iso_id(str)
         end function read_iso
         
         
         integer function read_int()
            character (len=64) :: str
            integer :: ierr
            str = line(i:j)
            read(str,fmt=*,iostat=ierr) read_int
            if (ierr /= 0) read_int = 0
         end function read_int
         
         
         real(dp) function read_dbl()
            use math_lib, only: str_to_double
            integer :: ierr
            ierr = 0
            if (len_trim(line(i:j)) > 0) then
               call str_to_double(line(i:j),read_dbl,ierr)
               if (ierr /= 0) read_dbl = 0d0
            else
               read_dbl = 0d0
            end if
         end function read_dbl
         
         
         logical function missing_dbl()
            character (len=64) :: str
            integer :: ierr
            str = line(i:j)
            missing_dbl = (len_trim(str) == 0)
         end function missing_dbl


         subroutine check_std_reaction_Qs
            integer :: i, cnt
            cnt = 0
            do i=1,num_reactions
               if (std_reaction_Qs(i) < -1d50) then
                  write(*,*) 'missing reaction_Q for reaction ' // trim(reaction_Name(i)), i, cnt+1
                  write(*,*) 
                  cnt = cnt+1
               end if
            end do
            if (cnt > 0) stop 'check_std_reaction_Qs'
         end subroutine check_std_reaction_Qs
      
      
         subroutine check_std_reaction_neuQs
            integer :: i, cnt
            cnt = 0
            do i=1,num_reactions
               if (std_reaction_neuQs(i) < -1d50) then
                  write(*,*) 'missing std_reaction_neuQs for reaction ' // trim(reaction_Name(i))
                  write(*,*) 
                  cnt = cnt+1
               end if
            end do
            if (cnt > 0) stop 'check_std_reaction_neuQs'
         end subroutine check_std_reaction_neuQs
      
      
         subroutine check_reaction_categories
            integer :: cnt, i
            cnt = 0
            do i=1,num_reactions
               if (reaction_categories(i) < 0) then
                  write(*,*) 'missing reaction_category for reaction ' // trim(reaction_Name(i))
                  write(*,*) 
                  cnt = cnt+1
               end if
            end do
            if (cnt > 0) stop 'check_reaction_categories'      
         end subroutine check_reaction_categories


         subroutine check_reaction_info
            integer :: i, cnt
            cnt = 0
            do i=1,num_reactions
               if (len_trim(reaction_Info(i)) == 0) then
                  write(*,*) 'missing info for reaction', i
                  if (i > 1) write(*,*) 'following ' // trim(reaction_Info(i-1))
                  write(*,*) 
                  cnt = cnt+1
               end if
            end do
            if (cnt > 0) stop 'check_reaction_info'
         end subroutine check_reaction_info


      end subroutine read_reaction_parameters


      subroutine increase_num_reactions(ierr)
         integer, intent(out) :: ierr
         
         integer :: old_max, new_max, i
         type (rate_table_info), pointer :: old_raw_rates_records(:)
         type (rate_table_info), pointer :: ri
         
         include 'formats'

         old_max = rates_reaction_id_max
         rates_reaction_id_max = rates_reaction_id_max + 1
      
         if (rates_reaction_id_max > size(std_reaction_Qs,dim=1)) then
         
            new_max = rates_reaction_id_max*2 + 1000
            !write(*,3) 'increase size', rates_reaction_id_max, new_max
            
            old_raw_rates_records => raw_rates_records
            allocate(raw_rates_records(new_max))
            do i=1,old_max
               raw_rates_records(i) = old_raw_rates_records(i)
            end do
            deallocate(old_raw_rates_records)
            
            call grow_reactions_arrays(old_max, new_max, ierr)
            if (ierr /= 0) return
         end if

         ri => raw_rates_records(rates_reaction_id_max)
         ri% nT8s = 0
         ri% use_rate_table = .false.
         ri% need_to_read = .false.
         nullify(ri% T8s)
         nullify(ri% f1)

      end subroutine increase_num_reactions


      subroutine alloc_and_init_reaction_parameters(ierr)
         integer, intent(out) :: ierr
      
         allocate( &
               reaction_Info(rates_reaction_id_max),  &
               reaction_categories(rates_reaction_id_max),  &
               reaction_is_reverse(rates_reaction_id_max),  &
               reaction_reaclib_lo(rates_reaction_id_max),  &
               reaction_reaclib_hi(rates_reaction_id_max),  &
               reverse_reaction_id(rates_reaction_id_max),  &
               reaction_screening_info(3,rates_reaction_id_max),  &
               weak_reaction_info(2,rates_reaction_id_max),  &
               reaction_ye_rho_exponents(2,rates_reaction_id_max),  &
               reaction_inputs(2*max_num_reaction_inputs,rates_reaction_id_max),  &
               reaction_outputs(2*max_num_reaction_outputs,rates_reaction_id_max),  &
               std_reaction_Qs(rates_reaction_id_max),  &
               std_reaction_neuQs(rates_reaction_id_max),   &
               weak_lowT_rate(rates_reaction_id_max),   &
               stat=ierr)
         if (ierr /= 0) return
         
         reaction_Info(:) = ''
         reaction_categories(:) = -1
         reaction_is_reverse(:) = 0
         reaction_reaclib_lo(:) = 0
         reaction_reaclib_hi(:) = 0
         reverse_reaction_id(:) = 0
         reaction_screening_info(:,:) = 0
         weak_reaction_info(:,:) = 0
         reaction_ye_rho_exponents(:,:) = 0
         reaction_inputs(:,:) = 0
         reaction_outputs(:,:) = 0
         std_reaction_Qs(:) = -1d99
         std_reaction_neuQs(:) = -1d99
         weak_lowT_rate(:) = -1d99
      
      end subroutine alloc_and_init_reaction_parameters
      
      
      real(dp) function get_Qtotal(ir)
         use chem_lib, only: reaction_Qtotal
         use chem_def, only: chem_isos
         integer, intent(in) :: ir
         
         integer :: num_in, num_out, reactants(100), k, n, i, ii, j
         include 'formats.dek'
         i = 0
         do k = 1, 2*max_num_reaction_inputs-1, 2
            n = reaction_inputs(k,ir)
            if (n == 0) exit
            ii = reaction_inputs(k+1,ir)
            do j = 1, n
               i = i+1
               reactants(i) = ii
            end do
         end do
         num_in = i
         do k = 1, 2*max_num_reaction_outputs-1, 2
            n = reaction_outputs(k,ir)
            if (n == 0) exit
            ii = reaction_outputs(k+1,ir)
            do j = 1, n
               i = i+1
               reactants(i) = ii
            end do
         end do
         num_out = i - num_in
         get_Qtotal = reaction_Qtotal(num_in,num_out,reactants,chem_isos)
      end function get_Qtotal


   
      
      subroutine init_rates_info(reactionlist_filename, ierr)
         character (len=*), intent(in) :: reactionlist_filename
         integer, intent(out) :: ierr ! 0 means AOK.           
         include 'formats'
         
         ierr = 0
         call init1_rates_info

         call start_rates_def_init(ierr)
         if (ierr /= 0) then
            write(*,*) 'start_rates_def_init failed'
            return
         end if
         
         call read_reaction_parameters(reactionlist_filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'rates_init failed in read_reaction_parameters'
            return
         end if

         call finish_rates_def_init
         call do_rates_init(ierr)
         if (ierr /= 0) then
            write(*,*) 'rates_init failed in do_rates_init'
            return
         end if        
      
      end subroutine init_rates_info
      
      
      subroutine init1_rates_info
         use rates_names, only: set_reaction_names
         type (rate_table_info), pointer :: ri
         integer :: i
         rates_reaction_id_max = num_predefined_reactions
         allocate( &
            raw_rates_records(rates_reaction_id_max), &
            reaction_Name(rates_reaction_id_max))
         do i = 1, rates_reaction_id_max
            ri => raw_rates_records(i)
            ri% nT8s = 0
            ri% use_rate_table = .false.
            ri% need_to_read = .false.
            nullify(ri% T8s)
            nullify(ri% f1)
         end do
         call set_reaction_names
      end subroutine init1_rates_info
      
         
      integer function lookup_rate_name(str) ! -1 if not found
         use rates_def
         character (len=*), intent(in) :: str
         integer :: i
         lookup_rate_name = -1
         do i = 1, rates_reaction_id_max
            if (trim(str) == trim(reaction_Name(i))) then
               lookup_rate_name = i
               return
            end if
         end do
      end function lookup_rate_name

      subroutine grow_reactions_arrays(old_max, new_max, ierr)
         use utils_lib
         integer, intent(in) :: old_max, new_max
         integer, intent(out) :: ierr

         character (len=maxlen_reaction_Name), pointer :: new_reaction_Name(:)
         character (len=maxlen_reaction_Info), pointer :: new_reaction_Info(:)
         integer :: i
         
         include 'formats'
         
         allocate(new_reaction_Name(new_max))
         do i=1,old_max
            new_reaction_Name(i) = reaction_Name(i)
         end do
         deallocate(reaction_Name)
         reaction_Name => new_reaction_Name

         allocate(new_reaction_Info(new_max))
         do i=1,old_max
            new_reaction_Info(i) = reaction_Info(i)
         end do
         deallocate(reaction_Info)
         reaction_Info => new_reaction_Info

         call realloc_integer(reaction_categories,new_max,ierr); if (ierr /= 0) return
         call realloc_integer(reaction_is_reverse,new_max,ierr); if (ierr /= 0) return
         call realloc_integer(reaction_reaclib_lo,new_max,ierr); if (ierr /= 0) return
         call realloc_integer(reaction_reaclib_hi,new_max,ierr); if (ierr /= 0) return
         call realloc_integer(reverse_reaction_id,new_max,ierr); if (ierr /= 0) return

         call realloc_double(std_reaction_Qs,new_max,ierr); if (ierr /= 0) return
         call realloc_double(std_reaction_neuQs,new_max,ierr); if (ierr /= 0) return
         
         call realloc_double(weak_lowT_rate,new_max,ierr); if (ierr /= 0) return
         
         call realloc_integer2( &
            reaction_screening_info,size( &
            reaction_screening_info,dim=1),new_max,ierr); if (ierr /= 0) return
         call realloc_integer2( &
            weak_reaction_info,2,new_max,ierr); if (ierr /= 0) return
         call realloc_integer2( &
            reaction_ye_rho_exponents,2,new_max,ierr); if (ierr /= 0) return
         call realloc_integer2( &
            reaction_inputs,2*max_num_reaction_inputs,new_max,ierr); if (ierr /= 0) return
         call realloc_integer2( &
            reaction_outputs,2*max_num_reaction_outputs,new_max,ierr); if (ierr /= 0) return
         
         reaction_Info(rates_reaction_id_max:new_max) = ''
         reaction_categories(rates_reaction_id_max:new_max) = -1

         reaction_is_reverse(rates_reaction_id_max:new_max) = 0
         reaction_reaclib_lo(rates_reaction_id_max:new_max) = 0
         reaction_reaclib_hi(rates_reaction_id_max:new_max) = 0
         reverse_reaction_id(rates_reaction_id_max:new_max) = 0

         reaction_screening_info(:,rates_reaction_id_max:new_max) = 0
         weak_reaction_info(:,rates_reaction_id_max:new_max) = 0
         reaction_ye_rho_exponents(:,rates_reaction_id_max:new_max) = 0
         reaction_inputs(:,rates_reaction_id_max:new_max) = 0
         reaction_outputs(:,rates_reaction_id_max:new_max) = 0
         std_reaction_Qs(rates_reaction_id_max:new_max) = -1d99
         std_reaction_neuQs(rates_reaction_id_max:new_max) = -1d99
         weak_lowT_rate(rates_reaction_id_max:new_max) = -1d99
      
      end subroutine grow_reactions_arrays

      
      subroutine free_reaction_arrays()

        if (ASSOCIATED(reaction_Name)) deallocate(reaction_Name)
        if (ASSOCIATED(reaction_Info)) deallocate(reaction_Info)

        if (ASSOCIATED(reaction_categories)) deallocate(reaction_categories)
        if (ASSOCIATED(reaction_is_reverse)) deallocate(reaction_is_reverse)
        if (ASSOCIATED(reaction_reaclib_lo)) deallocate(reaction_reaclib_lo)
        if (ASSOCIATED(reaction_reaclib_hi)) deallocate(reaction_reaclib_hi)
        if (ASSOCIATED(reverse_reaction_id)) deallocate(reverse_reaction_id)

        if (ASSOCIATED(std_reaction_Qs)) deallocate(std_reaction_Qs)
        if (ASSOCIATED(std_reaction_neuQs)) deallocate(std_reaction_neuQs)

        if (ASSOCIATED(weak_lowT_rate)) deallocate(weak_lowT_rate)

        if (ASSOCIATED(reaction_screening_info)) deallocate(reaction_screening_info)
        if (ASSOCIATED(weak_reaction_info)) deallocate(weak_reaction_info)
        if (ASSOCIATED(reaction_ye_rho_exponents)) deallocate(reaction_ye_rho_exponents)
        if (ASSOCIATED(reaction_inputs)) deallocate(reaction_inputs)
        if (ASSOCIATED(reaction_outputs)) deallocate(reaction_outputs)

        nullify(reaction_Name)
        nullify(reaction_Info)
        nullify(reaction_is_reverse)
        nullify(reaction_categories)
        nullify(reaction_reaclib_lo)
        nullify(reaction_reaclib_hi)
        nullify(reverse_reaction_id)
        nullify(std_reaction_Qs)
        nullify(std_reaction_neuQs)
        nullify(weak_lowT_rate)
        nullify(reaction_screening_info)
        nullify(reaction_ye_rho_exponents)
        nullify(reaction_inputs)
        nullify(reaction_outputs)

      end subroutine free_reaction_arrays

      
      subroutine do_rates_init(ierr)
         use ratelib, only: mazurek_init
         integer, intent(out) :: ierr
         ierr = 0
         ! setup interpolation info for mazurek's 1973 fits for the ni56 electron capture rate
         call mazurek_init(ierr)
         if (ierr /= 0) return
      end subroutine do_rates_init


      end module rates_initialize


