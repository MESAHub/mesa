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

      module rates_lib
      
      use const_def, only: dp
      use utils_lib, only: mesa_error
      
      implicit none


      contains
      
      
      ! call this routine to initialize the rates module. 
      ! only needs to be done once at start of run.
            
      subroutine rates_init( &
           reactionlist_filename, jina_reaclib_filename, &
           rates_table_dir_in, &
           use_suzuki_weak_rates, &
           use_special_weak_rates, &
           special_weak_states_file, &
           special_weak_transitions_file, &
           cache_dir, &
           ierr)
         use rates_def
         use reaclib_input, only: do_read_reaclib
         use load_weak, only: load_weak_data
         use load_ecapture, only: load_ecapture_data
         use rates_initialize, only: init_rates_info
         
         character (len=*), intent(in) :: reactionlist_filename, jina_reaclib_filename, rates_table_dir_in
         logical, intent(in) :: use_special_weak_rates, use_suzuki_weak_rates
         character (len=*), intent(in) :: special_weak_states_file, special_weak_transitions_file
         character (len=*), intent(in) :: cache_dir ! '' means use default
         integer, intent(out) :: ierr ! 0 means AOK.  
         
         logical, parameter :: dbg = .false.
         include 'formats'
         
         ierr = 0

         rates_table_dir = rates_table_dir_in

         call set_rates_cache_dir(cache_dir, ierr)
         if (ierr /= 0) return

         use_suzuki_tables = use_suzuki_weak_rates
         
         if (dbg) write(*,*) 'call weaklib_init'
         call weaklib_init(ierr)
         if (ierr /= 0) return
         call load_weak_data(ierr)
         if (ierr /= 0) return

         ! map special weak rates controls to old names
         do_ecapture = use_special_weak_rates
         ecapture_states_file = special_weak_states_file
         ecapture_transitions_file = special_weak_transitions_file
         
         if (dbg) write(*,*) 'call ecapture_init'
         call ecapture_init(ierr)
         if (ierr /= 0) return
         if (do_ecapture) call load_ecapture_data(ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call reaclib_init'
         call reaclib_init(jina_reaclib_filename)
         if (dbg) write(*,*) 'call do_read_reaclib'
         call do_read_reaclib(ierr)
         if (ierr /= 0) return
         
         if (dbg) write(*,*) 'call init_rates_info'
         call init_rates_info(reactionlist_filename, ierr)
         if (ierr /= 0) return
        
         have_finished_initialization = .true.
         
      end subroutine rates_init
      
      
      subroutine rates_warning_init( &
           warn_rates_for_high_temp_in, max_safe_logT_for_rates_in)
         use rates_def
         logical, intent(in) :: warn_rates_for_high_temp_in
         real(dp), intent(in) :: max_safe_logT_for_rates_in
         ! Setup warnings
         warn_rates_for_high_temp = warn_rates_for_high_temp_in
         max_safe_logT_for_rates = max_safe_logT_for_rates_in           
      end subroutine rates_warning_init

      
      subroutine read_raw_rates_records(ierr)
         use rates_initialize, only: init_raw_rates_records
         integer, intent(out) :: ierr ! 0 means AOK.  
         ierr = 0
         call init_raw_rates_records(ierr)
         if (ierr /= 0) then
            write(*,*) 'rates_init failed in init_raw_rates_records'
            return
         end if
      end subroutine read_raw_rates_records


      subroutine read_rate_from_file(rate_name, filename, ierr)
         use rates_initialize
         character(len=*), intent(in) :: rate_name, filename
         integer, intent(out) :: ierr

         call rate_from_file(rate_name, filename, ierr)
         if(ierr/=0) then
            write(*,*) 'read_rate_from_file for',trim(rate_name), trim(filename)
            return
         end if

      end subroutine read_rate_from_file

      subroutine read_rates_from_files(rate_names, filenames, ierr)
         use rates_initialize
         character(len=*), intent(in) :: rate_names(:), filenames(:)
         integer, intent(out) :: ierr
         integer :: i

         ierr = 0
         do i=1, ubound(rate_names,dim=1)
            call rate_from_file(rate_names(i), filenames(i), ierr)
            if(ierr/=0) then
               write(*,*) 'read_rates_from_files for num=',i,trim(rate_names(i)), trim(filenames(i))
               return
            end if
         end do

      end subroutine read_rates_from_files      

      
      subroutine rates_shutdown

         use rates_def
         use rates_initialize, only: free_reaction_arrays, free_raw_rates_records
         use reaclib_input, only: reaclib
         use utils_lib

         call integer_dict_free(skip_warnings_dict)
         call integer_dict_free(reaction_names_dict)
         
         call free_weak_info()
         call free_ecapture_info()
         call free_reaclib_data(reaclib)
         call free_reaction_data(reaclib_rates)
         call free_reaction_arrays()
         call free_raw_rates_records()
         
      end subroutine rates_shutdown
         
         
      subroutine add_reaction_from_reaclib(reaction_handle, reverse_handle, indx, ierr)
         use rates_initialize, only: do_add_reaction_from_reaclib
         character (len=*), intent(in) :: reaction_handle ! to be added
         character (len=*), intent(in) :: reverse_handle ! = '' if not a reverse
         integer, intent(in) :: indx ! index in reaclib rates
         integer, intent(out) :: ierr
         call do_add_reaction_from_reaclib(reaction_handle, reverse_handle, indx, ierr)
      end subroutine add_reaction_from_reaclib
         
         
      subroutine add_reaction_for_handle(handle, ierr)
         use rates_initialize, only: do_add_reaction_for_handle
         character (len=*), intent(in) :: handle ! to be added
         integer, intent(out) :: ierr
         call do_add_reaction_for_handle(handle, ierr)
      end subroutine add_reaction_for_handle


      subroutine make_rate_tables( &
            num_reactions, cache_suffix, net_reaction_id, &
            rattab, rattab_f1, nT8s, ttab, logttab, ierr)  
         use rates_support, only : do_make_rate_tables
         integer, intent(in) :: num_reactions, nT8s, net_reaction_id(:)
         character (len=*), intent(in) :: cache_suffix
         real(dp) :: rattab(:,:), ttab(:), logttab(:)
         real(dp), pointer :: rattab_f1(:)
         integer, intent(out) :: ierr
         call do_make_rate_tables( &
               num_reactions, cache_suffix, net_reaction_id,  &
               rattab, rattab_f1, nT8s, ttab, logttab, ierr)
      end subroutine make_rate_tables
      
           
      subroutine show_reaction_rates_from_cache(cache_filename, ierr) 
         use rates_support, only: do_show_reaction_from_cache
         character (len=*) :: cache_filename
         integer, intent(out) :: ierr
         call do_show_reaction_from_cache(cache_filename, ierr) 
      end subroutine show_reaction_rates_from_cache
      
            
      subroutine extract_reaclib_rates(set,nuclides,rates,use_weaklib,ierr)
         use rates_def
         use chem_def, only: nuclide_set, nuclide_data
         use reaclib_input, only: do_extract_rates
         type(nuclide_set), dimension(:), intent(in) :: set
         type(nuclide_data), intent(in), target :: nuclides
         logical, intent(in) :: use_weaklib
         type(reaction_data), intent(out) :: rates
         integer, intent(out) :: ierr
         call do_extract_rates(set,nuclides,rates,use_weaklib,ierr)
      end subroutine extract_reaclib_rates
      

      subroutine output_reaclib_rates(unitno,rates,nuclides,format)
         use reaclib_print
         use rates_def
         integer, intent(in) :: unitno
         type(reaction_data),intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         integer, intent(in) :: format
         integer :: err = 0
         select case (format)
            case(internal_format)
               call write_reaction_data(unitno,rates,err)
            case(pretty_print_format)
               call pretty_print_reactions(unitno,rates,nuclides,err)
            case(short_format)
               call print_short_format_reactions(unitno,rates,nuclides,err)
         end select
      end subroutine output_reaclib_rates
      
      
      subroutine reaclib_pretty_print_reaction(unitno, i, rates, nuclides, reverse, str, ierr)
         use reaclib_print, only: do_pretty_print_reaction
         use rates_def
         integer, intent(in) :: unitno, i
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: reverse
         character (len=100), intent(inout) :: str
         integer, intent(out) :: ierr
         call do_pretty_print_reaction(unitno, i, rates, nuclides, reverse, str, ierr)
      end subroutine reaclib_pretty_print_reaction
      
               
      subroutine eval_tfactors(tf, logT, temp)
         use rates_def, only : T_Factors
         use ratelib, only: tfactors
         type (T_Factors), pointer :: tf ! allocate this before calling
         real(dp), intent(in) :: logT, temp
         call tfactors(tf, logT, temp)
      end subroutine eval_tfactors
      

      subroutine get_raw_rate(ir, temp, tf, raw_rate, ierr)
         use rates_def, only : T_Factors
         use raw_rates
         integer, intent(in) :: ir
         real(dp), intent(in) :: temp
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: raw_rate
         integer, intent(out) :: ierr
         call set_raw_rate(ir,  temp, tf, raw_rate, ierr)
      end subroutine get_raw_rate


      subroutine get_raw_rates(n, irs, temp, tf, rates, ierr)
         use rates_def, only : T_Factors
         use raw_rates, only: set_raw_rates
         integer, intent(in) :: n
         integer, intent(in) :: irs(:) ! (n) maps 1..n to reaction id
         real(dp), intent(in) :: temp
         type (T_Factors), pointer :: tf
         real(dp), intent(inout) :: rates(:)
         integer, intent(out) :: ierr
         call set_raw_rates(n, irs, temp, tf, rates, ierr)
      end subroutine get_raw_rates
      
      integer function rates_reaction_id(rname)
         use rates_def, only: get_rates_reaction_id
         character (len=*), intent(in)  :: rname ! reaction name such as 'rpp' 
         ! returns id for the reaction if there is a matching entry in reaction_Name
         ! returns 0 otherwise.
         rates_reaction_id = get_rates_reaction_id(rname)
      end function rates_reaction_id
      
      integer function eval_num_reaction_inputs(ir)
         use rates_def, only: get_num_reaction_inputs
         integer, intent(in) :: ir
         eval_num_reaction_inputs = get_num_reaction_inputs(ir)
      end function eval_num_reaction_inputs
      
      integer function eval_num_reaction_outputs(ir)
         use rates_def, only: get_num_reaction_outputs
         integer, intent(in) :: ir
         eval_num_reaction_outputs = get_num_reaction_outputs(ir)
      end function eval_num_reaction_outputs

      
      subroutine rates_eval_reaclib_21( &
            ir, temp, den, rate_raw, reverse_rate_raw, ierr)
         use rates_support, only: do_eval_reaclib_21
         integer, intent(in) :: ir ! reaction_id
         real(dp), intent(in) :: temp, den
         real(dp), intent(inout) :: rate_raw(:), reverse_rate_raw(:)
         integer, intent(out) :: ierr
         call do_eval_reaclib_21( &
            ir, temp, den, rate_raw, reverse_rate_raw, ierr)
      end subroutine rates_eval_reaclib_21

      
      subroutine rates_eval_reaclib_22( &
            ir, temp, den, rate_raw, reverse_rate_raw, ierr)
         use rates_support, only: do_eval_reaclib_22
         integer, intent(in) :: ir ! reaction_id
         real(dp), intent(in) :: temp, den
         real(dp), intent(inout) :: rate_raw(:), reverse_rate_raw(:)
         integer, intent(out) :: ierr
         call do_eval_reaclib_22( &
            ir, temp, den, rate_raw, reverse_rate_raw, ierr)
      end subroutine rates_eval_reaclib_22
      
      
      subroutine rates_two_to_one_coeffs_for_reverse_factor( &
            Q, iso_A, iso_B, iso_C, a, b, ierr)
         use chem_def, only: chem_isos
         use math_lib
         real(dp), intent(in) :: Q
         integer, intent(in) :: iso_A, iso_B, iso_C
         real(dp), intent(out) :: a, b
         integer, intent(out) :: ierr
         real(dp) :: W_A, W_B, W_C, g_A, g_B, g_C         
         if (Q < 0) then
            ierr = -1
            return
         end if
         ierr = 0         
         W_A = chem_isos% W(iso_A)
         W_B = chem_isos% W(iso_B)
         W_C = chem_isos% W(iso_C)
         g_A = 2d0*chem_isos% spin(iso_A) + 1d0
         g_B = 2d0*chem_isos% spin(iso_B) + 1d0
         g_C = 2d0*chem_isos% spin(iso_C) + 1d0         
         ! Arnett, Supernovae and Nucleosynthesis, eqn 3.136
         a = 9.8678d9*(g_A*g_B/g_C)*pow(W_A*W_B/W_C,1.5d0)
         b = -11.605d0*Q         
      end subroutine rates_two_to_one_coeffs_for_reverse_factor
         
      
      ! note: assumes ground state spins and requires Q > 0.
      ! i.e., A + B -> C exothermic
      subroutine rates_two_to_one_reverse_factor( &
            Q, T9, T932, iso_A, iso_B, iso_C, rev, d_rev_dT, ierr) ! A + B <-> C
         use chem_def, only: chem_isos
         use math_lib
         real(dp), intent(in) :: Q, T9, T932
         integer, intent(in) :: iso_A, iso_B, iso_C
         real(dp), intent(out) :: rev, d_rev_dT
         integer, intent(out) :: ierr         
         real(dp) :: a, b      
         call rates_two_to_one_coeffs_for_reverse_factor( &
            Q, iso_A, iso_B, iso_C, a, b, ierr)
         if (ierr /= 0) return
         rev = a*T932*exp(b/T9)
         d_rev_dT = rev*(1.5d0*T9 - b)/(T9*T9*1d9)         
      end subroutine rates_two_to_one_reverse_factor


      subroutine rates_two_to_two_coeffs_for_reverse_factor( &
            Q, iso_A, iso_B, iso_C, iso_D, a, b, ierr)
         use chem_def, only: chem_isos
         real(dp), intent(in) :: Q
         integer, intent(in) :: iso_A, iso_B, iso_C, iso_D
         real(dp), intent(out) :: a, b
         integer, intent(out) :: ierr
         real(dp) :: W_A, W_B, W_C, W_D, g_A, g_B, g_C, g_D, a1
         if (Q < 0) then
            ierr = -1
            return
         end if
         ierr = 0         
         W_A = chem_isos% W(iso_A)
         W_B = chem_isos% W(iso_B)
         W_C = chem_isos% W(iso_C)
         W_D = chem_isos% W(iso_D)
         g_A = 2d0*chem_isos% spin(iso_A) + 1d0
         g_B = 2d0*chem_isos% spin(iso_B) + 1d0
         g_C = 2d0*chem_isos% spin(iso_C) + 1d0
         g_D = 2d0*chem_isos% spin(iso_D) + 0.5d0         
         ! Arnett, Supernovae and Nucleosynthesis, eqn 3.137
         a1 = ((g_A*g_B)/(g_C*g_D))*((W_A*W_B)/(W_C*W_D))
         a = a1*sqrt(a1)
         b = -11.605d0*Q         
      end subroutine rates_two_to_two_coeffs_for_reverse_factor
      
      
      ! note: assumes ground state spins and requires Q > 0.
      ! i.e., A + B -> C + D exothermic
      subroutine rates_two_to_two_reverse_factor( &
            Q, T9, iso_A, iso_B, iso_C, iso_D, rev, d_rev_dT, ierr) ! A + B <-> C + D
         use chem_def, only: chem_isos
         use math_lib
         real(dp), intent(in) :: Q, T9
         integer, intent(in) :: iso_A, iso_B, iso_C, iso_D
         real(dp), intent(out) :: rev, d_rev_dT
         integer, intent(out) :: ierr         
         real(dp) :: a, b
         call rates_two_to_two_coeffs_for_reverse_factor( &
            Q, iso_A, iso_B, iso_C, iso_D, a, b, ierr)
         if (ierr /= 0) return
         rev = a*exp(b/T9)
         d_rev_dT = -rev*b/(T9*T9*1d9)
      end subroutine rates_two_to_two_reverse_factor
            
      
      logical function is_weak_reaction(ir) ! not just weaklib.  any weak reaction.
         use rates_def, only: weak_reaction_info, std_reaction_neuQs
         integer, intent(in) :: ir ! reaction index
         is_weak_reaction = &
            (weak_reaction_info(1,ir) > 0 .and. weak_reaction_info(2,ir) > 0) .or. &
            (std_reaction_neuQs(ir) > 0)
      end function is_weak_reaction
      
      
      ! weaklib sources
      !   FFN: G.M. Fuller, W.A. Fowler, M.J. Newman, Ap. J. 293 (1985)
      !   OHMT: Oda, Hino, Muto, Takahara, and Sato. Atomic Data and Nuclear Data Tables, 1994.
      !   LMP: K. Langanke, G. Martínez-Pinedo / Nuclear Physics A 673 (2000) 481–508

      integer function get_weak_rate_id(lhs, rhs) ! returns 0 if reaction not found
         use rates_def, only: do_get_weak_rate_id
         character (len=*), intent(in)  :: lhs, rhs 
         get_weak_rate_id = do_get_weak_rate_id(lhs, rhs)
      end function get_weak_rate_id
      
      integer function get_weak_info_list_id(lhs, rhs) ! returns 0 if reaction not found
         ! value can be used to index weak_info_list_halflife and weak_info_list_Qneu
         use rates_def, only: do_get_weak_info_list_id
         character (len=*), intent(in)  :: lhs, rhs ! names as in weak_info.list file
         get_weak_info_list_id = do_get_weak_info_list_id(lhs, rhs)
      end function get_weak_info_list_id


      ! ecapture

      integer function get_ecapture_rate_id(lhs, rhs) ! returns 0 if reaction not found
         use rates_def
         use utils_lib
         character (len=*), intent(in)   :: lhs, rhs 
         ! names of the nuclides as given in ecapturereactions.tables (e.g. 'p', 'n', 'ca42', etc.)
         integer :: ierr, i
         character (len=2*iso_name_length+1) :: key
         character (len=iso_name_length) :: lhs_name, rhs_name
         ierr = 0
         get_ecapture_rate_id = 0
         lhs_name = adjustl(lhs)
         rhs_name = adjustl(rhs)
         call create_ecapture_dict_key(lhs_name, rhs_name, key)
         call integer_dict_lookup(ecapture_reactions_dict, key, i, ierr)
         if (ierr /= 0) then
             !write(*,*) 'failed in integer_dict_lookup for key ' // trim(key)
             return
         end if
         get_ecapture_rate_id = i
      end function get_ecapture_rate_id

      integer function get_ecapture_info_list_id(lhs, rhs) ! returns 0 if reaction not found
         ! value can be used to index ecapture_info_life_halflife and ecapture_info_list_Qneu
         use rates_def
         use utils_lib
         character (len=*), intent(in)   :: lhs, rhs ! names as in ecapture_info.list file
         integer :: ierr, i
         character (len=2*iso_name_length+1) :: key
         character (len=iso_name_length) :: lhs_name, rhs_name
         ierr = 0
         get_ecapture_info_list_id = 0
         lhs_name = adjustl(lhs)
         rhs_name = adjustl(rhs)
         call create_ecapture_dict_key(lhs_name, rhs_name, key)
         call integer_dict_lookup(ecapture_reactions_dict, key, i, ierr)
         if (ierr /= 0) then
             !write(*,'(a)') 'get_ecapture_info_list_id failed for ' // trim(key)
             return
         end if
         get_ecapture_info_list_id = i
      end function get_ecapture_info_list_id

      ! reaclib
      
      subroutine reaclib_parse_handle(handle, num_in, num_out, iso_ids, op, ierr)
         use reaclib_support, only: do_parse_reaction_handle
         character (len=*), intent(in) :: handle
         integer, intent(out) :: num_in, num_out
         integer, intent(out) :: iso_ids(:) ! holds chem_ids for input and output species
         character (len=*), intent(out) :: op ! e.g., 'pg', 'wk', 'to', or ...
         integer, intent(out) :: ierr
         call do_parse_reaction_handle(handle, num_in, num_out, iso_ids, op, ierr)  
      end subroutine reaclib_parse_handle
      
      subroutine reaclib_create_handle(num_in, num_out, iso_ids, handle)
         use reaclib_support, only: reaction_handle
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: iso_ids(:) ! holds chem_ids for input and output species
         character (len=*), intent(out) :: handle
         character (len=1) :: reaction_flag = '-'
         call reaction_handle(num_in, num_out, iso_ids, reaction_flag, handle)   
      end subroutine reaclib_create_handle
      
      subroutine reaclib_create_ec_handle(num_in, num_out, iso_ids, handle)
         use reaclib_support, only: reaction_handle
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: iso_ids(:) ! holds chem_ids for input and output species
         character (len=*), intent(out) :: handle
         character (len=1) :: reaction_flag = 'e'
         call reaction_handle(num_in, num_out, iso_ids, reaction_flag, handle)   
      end subroutine reaclib_create_ec_handle
      
      subroutine reaclib_create_wk_handle(num_in, num_out, iso_ids, handle)
         use reaclib_support, only: reaction_handle
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: iso_ids(:) ! holds chem_ids for input and output species
         character (len=*), intent(out) :: handle
         character (len=1) :: reaction_flag = 'w'
         call reaction_handle(num_in, num_out, iso_ids, reaction_flag, handle)   
      end subroutine reaclib_create_wk_handle
      
      subroutine reaclib_create_reverse_handle(num_in, num_out, iso_ids, handle)
         use reaclib_support, only: reverse_reaction_handle
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: iso_ids(:) ! holds chem_ids for input and output species
         character (len=*), intent(out) :: handle
         call reverse_reaction_handle(num_in, num_out, iso_ids, handle) 
      end subroutine reaclib_create_reverse_handle
      
      
      integer function reaclib_lookup(handle, rates_dict) result(indx)
         ! returns first reaction index that matches handle. 
         ! there may be several following that one having the same handle.
         ! returns 0 if handle doesn't match any of the reactions
         use rates_def
         use reaclib_eval, only: do_reaclib_lookup
         character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         type (integer_dict), pointer :: rates_dict ! from create_reaclib_rates_dict
         indx = do_reaclib_lookup(handle, rates_dict)
      end function reaclib_lookup
      
      subroutine create_reaction_handle( &
            num_in, num_out, pspecies, nuclides, reverse, reaction_flag, handle)
         use reaclib_support, only: get1_reaction_handle
         use rates_def
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: pspecies(:)
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: reverse
         character (len=*), intent(in) :: reaction_flag
         character (len=*), intent(out) :: handle
         call get1_reaction_handle( &
            num_in, num_out, pspecies, nuclides, reverse, reaction_flag, handle)
      end subroutine create_reaction_handle
      
      
      subroutine reaclib_indices_for_reaction(handle, rates, lo, hi, ierr)
         use reaclib_eval, only: do_reaclib_indices_for_reaction
         use rates_def
         character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         type(reaction_data), intent(in) :: rates
         integer, intent(out) :: lo, hi
         integer, intent(out) :: ierr
         call do_reaclib_indices_for_reaction(handle, rates, lo, hi, ierr)
      end subroutine reaclib_indices_for_reaction
      
      
      subroutine reaclib_reaction_rates( &
            lo, hi, T9, rates, nuclides, forward_only, &
            lambda, dlambda_dlnT, &
            rlambda, drlambda_dlnT, &
            ierr)
         use rates_def
         use reaclib_eval, only: do_reaclib_reaction_rates
         integer, intent(in) :: lo, hi ! from reaclib_indices_for_reaction
         real(dp), intent(in) :: T9
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: forward_only
         real(dp), intent(out) :: lambda, dlambda_dlnT
         real(dp), intent(out) :: rlambda, drlambda_dlnT
         integer, intent(out) :: ierr
         call do_reaclib_reaction_rates( &
            lo, hi, T9, rates, nuclides, forward_only, &
            lambda, dlambda_dlnT, &
            rlambda, drlambda_dlnT, &
            ierr)
      end subroutine reaclib_reaction_rates
      
      
      ! screen
      
      subroutine screen_init_AZ_info( &
               a1, z1, a2, z2, &
               zs13, zhat, zhat2, lzav, aznut, zs13inv, &
               ierr)
         use screen5, only: screen5_init_AZ_info
         real(dp), intent(in) :: a1, z1, a2, z2
         real(dp), intent(out) :: zs13, zhat, zhat2, lzav, aznut, zs13inv
         integer, intent(out) :: ierr
         if (ierr /= 0) return
         call screen5_init_AZ_info( &
               zs13, zhat, zhat2, lzav, aznut, zs13inv, a1, z1, a2, z2, ierr)
      end subroutine screen_init_AZ_info
      
      integer function screening_option(which_screening_option, ierr)
         use rates_def
         use utils_lib, only: StrLowCase
         character (len=*), intent(in) :: which_screening_option
         integer, intent(out) :: ierr
         character (len=64) :: option
         ierr = 0

         option = StrLowCase(which_screening_option)  
         if (associated(rates_other_screening)) then
            screening_option = other_screening
         else if (option == 'no_screening' .or. len_trim(option) == 0) then
            screening_option = no_screening                   
         else if (option == 'extended') then
            screening_option = extended_screening           
         else if (option == 'salpeter') then
            screening_option = salpeter_screening    
         else if (option == 'chugunov') then
            screening_option = chugunov_screening  
         else
            ierr = -1
            screening_option = -1
         end if        
      end function screening_option
      
      subroutine screening_option_str(which_screening_option, screening_option, ierr)
         use rates_def
         integer, intent(in) :: which_screening_option
         character (len=*), intent(out) :: screening_option
         integer, intent(out) :: ierr
         ierr = 0    
         
         if (which_screening_option == other_screening) then
            screening_option = 'other_screening'
         else if (which_screening_option == no_screening) then
            screening_option = 'no_screening'                    
         else if (which_screening_option == extended_screening) then
            screening_option = 'extended'            
         else if (which_screening_option == salpeter_screening) then
            screening_option = 'salpeter'        
         else if (which_screening_option == chugunov_screening) then
            screening_option = 'chugunov'  
         else
            ierr = -1
            screening_option = ''
         end if          
      end subroutine screening_option_str
      
      
      ! note: if do_ecapture is true, then eval_weak_reaction_info calls this.
      subroutine eval_ecapture_reaction_info( &
             n, ids, cc, T9, YeRho, &
             eta, d_eta_dlnT, d_eta_dlnRho, &
             lambda, dlambda_dlnT, dlambda_dlnRho, &
             Q, dQ_dlnT, dQ_dlnRho, &
             Qneu, dQneu_dlnT, dQneu_dlnRho, &
             ierr)
         use rates_def, only: Coulomb_Info
         use eval_ecapture, only: do_eval_ecapture_reaction_info
         integer, intent(in) :: n, ids(:) 
         type(Coulomb_Info), intent(in) :: cc
         real(dp), intent(in) :: T9, YeRho, eta, d_eta_dlnT, d_eta_dlnRho
         ! lambda = combined rate (capture and decay)
         ! Q and Qneu are for combined rate of beta decay and electron capture.
         ! Q is total, so Q-Qneu is the actual thermal energy.
         ! note: lambdas include Ye Rho factors for electron captures.
         ! so treat the rates as if just beta decays
         real(dp), dimension(:), pointer, intent(inout) :: &
                lambda, dlambda_dlnT, dlambda_dlnRho, &
                Q, dQ_dlnT, dQ_dlnRho, &
                Qneu, dQneu_dlnT, dQneu_dlnRho
         integer, intent(out) :: ierr
         call do_eval_ecapture_reaction_info( &
                n, ids, cc, T9, YeRho, &
                eta, d_eta_dlnT, d_eta_dlnRho, &
                lambda, dlambda_dlnT, dlambda_dlnRho, &
                Q, dQ_dlnT, dQ_dlnRho, &
                Qneu, dQneu_dlnT, dQneu_dlnRho, &
                ierr)
      end subroutine eval_ecapture_reaction_info
      
      
      subroutine eval_salpeter_screening(sc, z1, z2, scor, scordt, scordd, ierr)
         ! weak screening only.  following Salpeter (1954),
         ! with equations (4-215) and (4-221) of Clayton (1968).
         use rates_def
         use math_lib
         type (Screen_Info) :: sc ! previously setup 
         real(dp), intent(in) :: z1, z2
         real(dp), intent(out) :: scor ! screening factor
         real(dp), intent(out) :: scordt ! partial wrt temperature
         real(dp), intent(out) :: scordd ! partial wrt density
         integer, intent(out) :: ierr
         real(dp) :: zeta, lnf, rho, T, dlnf_dd, dlnf_dt
         ierr = 0
         rho = sc% den
         T = sc% temp
         zeta = (sc% z2bar + sc% zbar) / sc% abar
         lnf = 1.88d8*z1*z2*sqrt(rho*zeta/(T*T*T))
         dlnf_dd = lnf/(2*rho)
         dlnf_dt = -lnf*3/(2*T)
         scor = exp(lnf)
         scordd = scor*dlnf_dd
         scordt = scor*dlnf_dt
      end subroutine eval_salpeter_screening

      subroutine eval_weak_reaction_info( &
            n, ids, reaction_ids, cc, T9, YeRho, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)
           use rates_def, only: Coulomb_Info
         use eval_weak, only: do_eval_weak_reaction_info
           use rates_def, only : do_ecapture
         integer, intent(in) :: n, ids(:), reaction_ids(:)
           type(Coulomb_Info), pointer :: cc
         real(dp), intent(in) :: T9, YeRho, eta, d_eta_dlnT, d_eta_dlnRho
         ! lambda = combined rate (capture and decay)
         ! Q and Qneu are for combined rate of beta decay and electron capture.
         ! Q is total, so Q-Qneu is the actual thermal energy.
         ! note: lambdas include Ye Rho factors for electron captures.
         ! so treat the rates as if just beta decays
         real(dp), dimension(:), pointer :: &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho
         integer, intent(out) :: ierr
         call do_eval_weak_reaction_info( &
            n, ids, reaction_ids, T9, YeRho, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)
         if (ierr /= 0) return
         if (.not. do_ecapture) return 
         call eval_ecapture_reaction_info( &
             n, ids, cc, T9, YeRho, &
             eta, d_eta_dlnT, d_eta_dlnRho, &
             lambda, dlambda_dlnT, dlambda_dlnRho, &
             Q, dQ_dlnT, dQ_dlnRho, &
             Qneu, dQneu_dlnT, dQneu_dlnRho, &
             ierr)
      end subroutine eval_weak_reaction_info
      
      subroutine eval_using_rate_tables( &
            num_reactions, reaction_id, rattab, rattab_f1, nT8s, &
            ye, logtemp, btemp, bden, raw_rate_factor, logttab, &
            rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
         use rates_support, only : do_get_raw_rates
         integer, intent(in) :: num_reactions, reaction_id(:), nT8s
         real(dp), intent(in) ::  &
            ye, logtemp, btemp, bden, raw_rate_factor(:),  &
            rattab(:,:), logttab(:)
         real(dp), pointer :: rattab_f1(:)
         real(dp), intent(out), dimension(:) :: rate_raw, rate_raw_dT, rate_raw_dRho
         integer, intent(out) :: ierr
         call do_get_raw_rates(num_reactions, reaction_id, rattab, rattab_f1, nT8s, &
               ye, logtemp, btemp, bden, raw_rate_factor, logttab, &
               rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
      end subroutine eval_using_rate_tables
      
      ! call this once before calling screen_pair for each reaction
      ! sets info that depends only on temp, den, and overall composition         
      subroutine screen_set_context( &
            sc, temp, den, logT, logRho, zbar, abar, z2bar,  &
            screening_mode, num_isos, y, iso_z158)
         use rates_def
         use screen, only: do_screen_set_context
         type (Screen_Info) :: sc
         integer, intent(in) :: num_isos
         real(dp), intent(in) ::  &
               temp, den, logT, logRho, zbar, abar, z2bar, y(:), iso_z158(:)
            ! y(:) = x(:)/chem_A(chem_id(:))
            ! iso_z(:) = chem_Z(chem_id(:))**1.58
         integer, intent(in) :: screening_mode 
         call do_screen_set_context( &
            sc, temp, den, logT, logRho, zbar, abar, z2bar, &
            screening_mode, num_isos, y, iso_z158)

      end subroutine screen_set_context


      ! set jscr = 0 before 1st call.
      ! make calls in exactly the same order as for screen_init_AZ_info   
      subroutine screen_pair( &
               sc, a1, z1, a2, z2, screening_mode, &
               zs13, zhat, zhat2, lzav, aznut, zs13inv, low_logT_lim, &
               scor, scordt, scordd, ierr)
         use rates_def
         use screen5, only: fxt_screen5
         use screening_chugunov, only: eval_screen_chugunov
         
         type (Screen_Info) :: sc ! previously setup 
         real(dp), intent(in) :: a1, z1, a2, z2
         integer, intent(in) :: screening_mode ! see screen_def.
         ! cached info
         real(dp), intent(in) :: zs13, zhat, zhat2, lzav, aznut, zs13inv
         real(dp), intent(in) :: low_logT_lim ! scor==0 for T < low_logT_lim
         ! outputs
         real(dp), intent(out) :: scor ! screening factor
         real(dp), intent(out) :: scordt ! partial wrt temperature
         real(dp), intent(out) :: scordd ! partial wrt density
         integer, intent(out) :: ierr
         
         if(sc% logT < low_logT_lim ) then
            scor = 0d0
            scordt=0d0
            scordd=0d0
            return
         end if

         if(screening_mode == other_screening) then
            call rates_other_screening(sc, z1, z2, a1, a2, scor, scordt, scordd, ierr)
         else if (screening_mode == extended_screening) then
            call fxt_screen5( &
               sc, zs13, zhat, zhat2, lzav, aznut, zs13inv,  &
               a1, z1, a2, z2, scor, scordt, scordd, ierr)
         else if (screening_mode == salpeter_screening) then
            call eval_salpeter_screening(sc, z1, z2, scor, scordt, scordd, ierr)
         else if (screening_mode == chugunov_screening) then
            call eval_screen_chugunov(sc, z1, z2, a1, a2, scor, scordt, scordd, ierr)
         else if (screening_mode == no_screening) then
            scor = 1; scordt = 0; scordd = 0
         else
            ierr = -1
            write(*,*) 'screen_pair: unknown value for screening_mode', screening_mode
         end if
      end subroutine screen_pair

      subroutine eval_ecapnuc_rate(etakep,temp,rho,rpen,rnep,spen,snep)
         use ratelib, only: ecapnuc
         real(dp), intent(in) :: etakep,temp,rho
         real(dp), intent(out) :: rpen,rnep,spen,snep
         !  given the electron degeneracy parameter etakep (chemical potential
         !  without the electron's rest mass divided by kt) and the temperature temp,
         !  this routine calculates rates for 
         !  electron capture on protons rpen (captures/sec/proton),
         !  positron capture on neutrons rnep (captures/sec/neutron), 
         !  and their associated neutrino energy loss rates 
         !  spen (ergs/sec/proton) and snep (ergs/sec/neutron)
         call ecapnuc(etakep,temp,rho,rpen,rnep,spen,snep)
      end subroutine eval_ecapnuc_rate

      subroutine eval_mazurek_rate(btemp,bden,y56,ye,rn56ec,sn56ec)       
         use ratelib, only: mazurek
         real(dp), intent(in) :: btemp,bden,y56,ye
         real(dp), intent(out) :: rn56ec,sn56ec
         call mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)
      end subroutine eval_mazurek_rate
      
      subroutine eval_FL_epsnuc_3alf(T, Rho, Y, UE, eps_nuc, deps_nuc_dT, deps_nuc_dRho)       
         ! based on analytic expressions in Fushiki and Lamb, Apj, 317, 368-388, 1987.
         
         ! Note: if you plot the results of this, you'll see abrupt changes in rate at
         ! logRho about 9.74 and 10.25 -- these aren't bugs in the code. 
         ! They are discussed in F&L, and show up as step functions in their expressions.
         
         ! They provide expressions for both pyconuclear regime and strong screening regime.
         ! The transition between the regimes happens at U = 1, where U is defined below.
         ! Unfortunately, at U = 1, their expressions for pycnonuclear rate and
         ! strong screening rate disagree!
         ! Bummer.  For example, at logRho = 8.0, U = 1 for logT = 7.1955. 
         ! For these values, and pure He,
         ! their strong screening expression is larger than their pycno expression
         ! by a factor of about 25.
         
         ! need to add transition region in U instead of having an abrupt change at U = 1
         
         use pycno, only: FL_epsnuc_3alf
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: Y ! helium mass fraction
         real(dp), intent(in) :: UE ! electron molecular weight
         real(dp), intent(out) :: eps_nuc ! eps_nuc in ergs/g/sec
         real(dp), intent(out) :: deps_nuc_dT ! partial wrt temperature
         real(dp), intent(out) :: deps_nuc_dRho ! partial wrt density
         call FL_epsnuc_3alf(T, Rho, Y, UE, eps_nuc, deps_nuc_dT, deps_nuc_dRho)
      end subroutine eval_FL_epsnuc_3alf
      
      subroutine eval_n14_electron_capture_rate(T,Rho,UE,rate)
         use ratelib, only: n14_electron_capture_rate
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: UE ! electron molecular weight
         real(dp), intent(out) :: rate ! (s^-1)
         call n14_electron_capture_rate(T,Rho,UE,rate)
      end subroutine eval_n14_electron_capture_rate      


      ! call this once before calling eval_ecapture
      ! sets info that depends only on temp, den, and overall composition         
      subroutine coulomb_set_context( &
            cc, temp, den, logT, logRho, zbar, abar, z2bar)
         use rates_def, only: Coulomb_Info
         use coulomb, only: do_coulomb_set_context
         type (Coulomb_Info), pointer :: cc
         real(dp), intent(in) ::  &
               temp, den, logT, logRho, zbar, abar, z2bar
         call do_coulomb_set_context( &
            cc, temp, den, logT, logRho, zbar, abar, z2bar)
      end subroutine coulomb_set_context


      ! translate from option string to integer
      integer function get_mui_value(option_string)
        use eval_coulomb, only : None, DGC1973, I1993, PCR2009

        character(len=*), intent(in) :: option_string

        select case (trim(option_string))
        case('none')
           get_mui_value = None
        case('DGC1973')
           get_mui_value = DGC1973
        case('I1993')
           get_mui_value = I1993
        case('PCR2009')
           get_mui_value = PCR2009
        case DEFAULT
           call mesa_error(__FILE__,__LINE__,'Incorrect option for ion_coulomb_corrections')
        end select

        return

      end function get_mui_value


      integer function get_vs_value(option_string)
        use eval_coulomb, only : None, ThomasFermi, Itoh2002

        character(len=*), intent(in) :: option_string

        select case (trim(option_string))
        case('none')
           get_vs_value = None
        case('ThomasFermi')
           get_vs_value = ThomasFermi
        case('Itoh2002')
           get_vs_value = Itoh2002
        case DEFAULT
           call mesa_error(__FILE__,__LINE__,'Incorrect option for electron_coulomb_corrections')
        end select

        return

      end function get_vs_value

      subroutine rates_get_density_factors(ir, ye, rho, factor, factor_drho)
         use rates_support
         integer,intent(in) :: ir
         real(dp), intent(in) :: ye,rho
         real(dp),intent(out) :: factor, factor_drho

        call get_density_factors(ir, ye, rho, factor, factor_drho)

      end subroutine rates_get_density_factors


      end module rates_lib
