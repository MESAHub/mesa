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
 
      module net_lib
      ! library for calculating nuclear reaction rates and energy production
      ! the data interface for the library is defined in net_def
      
      use chem_def
      use const_def, only: dp
      
      implicit none


      contains ! the procedure interface for the library
      ! client programs should only call these routines.
      
      
      ! call this routine to initialize the net module. 
      ! only needs to be done once at start of run.
      
      subroutine net_init(ierr)      
         use net_def, only : do_net_def_init
         use net_initialize, only : init_special_case_reaction_info
         integer, intent(out) :: ierr ! 0 means AOK.      
         ierr = 0   
         call do_net_def_init
         call init_special_case_reaction_info
         
      end subroutine net_init
      
      
      subroutine net_shutdown
      end subroutine net_shutdown

      
      ! after net_init has finished, you can allocate a "handle".
      
      integer function alloc_net_handle(ierr)
         use net_def, only: do_alloc_net, init_net_handle_data
         integer, intent(out) :: ierr
         alloc_net_handle = do_alloc_net(ierr)
         if (ierr /= 0) return
      end function alloc_net_handle      
      
      subroutine free_net_handle(handle)
         ! frees the handle and all associated data
         use net_def, only: do_free_net
         integer, intent(in) :: handle
         call do_free_net(handle)
      end subroutine free_net_handle      
      
      ! if you want to access the Net_General_Info record directly, 
      ! you'll need a pointer to it.
      subroutine net_ptr(handle, g, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle ! from alloc_net_handle
         type (Net_General_Info), pointer :: g
         integer, intent(out):: ierr
         call get_net_ptr(handle, g, ierr)
      end subroutine net_ptr
      
      
      ! routines for defining the net isos and reactions

      
      ! call this before starting to define
      ! the set of isotopes and reactions for the net.
      subroutine net_start_def(handle, ierr)
         use net_initialize, only:start_net_def
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         call start_net_def(handle, ierr)
      end subroutine net_start_def
      
      
      ! call this after you've finished defining
      ! the set of isotopes and reactions for the net.
      subroutine net_finish_def(handle, ierr)
         use net_initialize, only: finish_net_def
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         call finish_net_def(handle, ierr)
      end subroutine net_finish_def
      
      ! note: after net_finish_def returns,
      ! you have the option of reordering the isotopes
      ! before you set up the full set of tables for the net.
      ! use get_chem_id_table_ptr and get_net_iso_table_ptr
      ! and change both tables to permute the set of isotopes.
      ! the default isotope ordering is by increasing chem_id number.
      
      
      ! read_net_file first tries opening the filename in the current directory.
      ! if doesn't find that file, then tries the data_dir from the call on net_init.
      ! i.e., looks for <data_dir>/net_data/nets/<filename>
      ! check net_data/nets/README for info about net files.
      subroutine read_net_file(filename, handle, ierr)
         use net_initialize, only: do_read_net_file
         character (len=*), intent(in) :: filename
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         call do_read_net_file(filename, handle, ierr)
      end subroutine read_net_file
      
      
      subroutine net_add_iso(handle, iso_id, ierr)
         use net_initialize, only:add_net_iso
         integer, intent(in) :: handle
         integer, intent(in) :: iso_id
         integer, intent(out) :: ierr
         call add_net_iso(handle, iso_id, ierr)
      end subroutine net_add_iso
      
      
      subroutine net_add_isos(handle, num_isos, iso_ids, ierr)
         use net_initialize, only:add_net_isos
         integer, intent(in) :: handle
         integer, intent(in) :: num_isos, iso_ids(num_isos)
         integer, intent(out) :: ierr
         call add_net_isos(handle, num_isos, iso_ids, ierr)
      end subroutine net_add_isos
      
      
      subroutine net_remove_iso(handle, iso_id, ierr)
         use net_initialize, only:remove_net_iso
         integer, intent(in) :: handle
         integer, intent(in) :: iso_id
         integer, intent(out) :: ierr
         call remove_net_iso(handle, iso_id, ierr)
      end subroutine net_remove_iso
      
      
      subroutine net_remove_isos(handle, num_isos, iso_ids, ierr)
         use net_initialize, only:remove_net_isos
         integer, intent(in) :: handle
         integer, intent(in) :: num_isos, iso_ids(num_isos)
         integer, intent(out) :: ierr
         call remove_net_isos(handle, num_isos, iso_ids, ierr)
      end subroutine net_remove_isos
      
      
      subroutine net_add_reaction(handle, reaction_id, ierr)
         use net_initialize, only:add_net_reaction
         integer, intent(in) :: handle
         integer, intent(in) :: reaction_id
         integer, intent(out) :: ierr
         call add_net_reaction(handle, reaction_id, ierr)
      end subroutine net_add_reaction
      
      
      subroutine net_add_reactions(handle, num_reactions, reaction_ids, ierr)
         use net_initialize, only:add_net_reactions
         integer, intent(in) :: handle
         integer, intent(in) :: num_reactions, reaction_ids(num_reactions)
         integer, intent(out) :: ierr
         call add_net_reactions(handle, num_reactions, reaction_ids, ierr)
      end subroutine net_add_reactions
      
      
      subroutine net_remove_reaction(handle, reaction_id, ierr)
         use net_initialize, only:remove_net_reaction
         integer, intent(in) :: handle
         integer, intent(in) :: reaction_id
         integer, intent(out) :: ierr
         call remove_net_reaction(handle, reaction_id, ierr)
      end subroutine net_remove_reaction
      
      
      subroutine net_remove_reactions(handle, num_reactions, reaction_ids, ierr)
         use net_initialize, only:remove_net_reactions
         integer, intent(in) :: handle
         integer, intent(in) :: num_reactions, reaction_ids(num_reactions)
         integer, intent(out) :: ierr
         call remove_net_reactions(handle, num_reactions, reaction_ids, ierr)
      end subroutine net_remove_reactions
      
      
      subroutine show_net_reactions(handle, iounit, ierr)
         use net_def
         use rates_def, only: reaction_Name
         integer, intent(in) :: handle
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer :: i, id
         type (Net_General_Info), pointer  :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         do i = 1, g% num_reactions
            id = g% reaction_id(i)
            if (id > 0) write(iounit,'(a)') trim(reaction_Name(id))
         end do
      end subroutine show_net_reactions

      
      subroutine show_net_reactions_and_info(handle, iounit, ierr)
         use net_def
         use rates_def, only:  &
               reaction_Name, reaction_Info, maxlen_reaction_Info, &
               std_reaction_Qs, std_reaction_neuQs, reaction_categories
         integer, intent(in) :: handle
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer :: i, id, weak_id, icat
         real(dp) :: Q, Qneu
         logical :: weaklib_reaction
         character (len=maxlen_reaction_Info) :: info
         type (Net_General_Info), pointer  :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         write(iounit,*)
         write(iounit,'(4x,a30,3a16,4x,a4,4x,a20)') 'name', 'Qtotal', 'Qneu', 'category', 'info', 'source'
         do i = 1, g% num_reactions
            id = g% reaction_id(i)
            weaklib_reaction = .false.
            weak_id = g% weak_reaction_index(i)
            if (weak_id > 0) then
               if (g% weaklib_ids(weak_id) > 0) weaklib_reaction = .true.
            end if
            if (id > 0) then
               info = reaction_Info(id)
               if (reaction_Name(id) == info) info = ''
               Q = std_reaction_Qs(id)
               Qneu = std_reaction_neuQs(id)
               icat = reaction_categories(id)
               if (weaklib_reaction) then
                  write(iounit,'(i4,a30,16x,a7,9x,4x,a10,4x,a66)') i, &
                        trim(reaction_Name(id)), 'weaklib', trim(category_name(icat)), info
               else if (Qneu /= 0) then
                  write(iounit,'(i4,a30,2f16.6,4x,a10,4x,a66)') i, trim(reaction_Name(id)),  &
                        Q, Qneu, trim(category_name(icat)), info
               else if (Q /= 0) then
                  write(iounit,'(i4,a30,f16.6,16x,4x,a10,4x,a66)') i, trim(reaction_Name(id)),  &
                        Q, trim(category_name(icat)), info
               else 
                  write(iounit,'(i4,a30,16x,16x,4x,a10,4x,a66)') i, trim(reaction_Name(id)),  &
                        trim(category_name(icat)), info
               end if
            end if
         end do
         write(iounit,*)
      end subroutine show_net_reactions_and_info
      
      
      subroutine show_net_species(handle, iounit, ierr)
         use net_def
         use chem_def
         integer, intent(in) :: handle
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer :: i, id
         type (Net_General_Info), pointer  :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         do i = 1, g% num_isos
            id = g% chem_id(i)
            if (id > 0) write(iounit,'(i4,a10)') i, trim(chem_isos% name(id))
         end do
      end subroutine show_net_species
      
      
      subroutine show_net_params(handle, iounit, ierr)
         use net_def
         integer, intent(in) :: handle
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer :: i, id
         type (Net_General_Info), pointer  :: g
         include 'formats'
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return         
         write(iounit,2) 'logTcut_lo =', g% logTcut_lo
         write(iounit,2) 'logTcut_lim =', g% logTcut_lim
      end subroutine show_net_params
      
      
      ! before calling net_setup_tables, you can specify options for various rates
      ! NOTE: these choices are used in building the rates tables, so if you want
      ! to change your choice, you'll need to rebuild the rate tables.
      
      
      subroutine net_set_which_rates(handle, which_rates, ierr)
         use net_def, only: do_net_set_which_rates
         integer, intent(in) :: handle, which_rates(:)
         integer, intent(out) :: ierr
         ierr = 0
         call do_net_set_which_rates(handle, which_rates, ierr)
         if (ierr /= 0) return
      end subroutine net_set_which_rates
      
      
      subroutine net_set_fe56ec_fake_factor( &
            handle, fe56ec_fake_factor, min_T_for_fe56ec_fake_factor, ierr)
         use net_def, only: do_net_set_fe56ec_fake_factor
         integer, intent(in) :: handle
         real(dp), intent(in) :: fe56ec_fake_factor, min_T_for_fe56ec_fake_factor
         integer, intent(out) :: ierr
         ierr = 0
         call do_net_set_fe56ec_fake_factor( &
            handle, fe56ec_fake_factor, min_T_for_fe56ec_fake_factor, ierr)
         if (ierr /= 0) return
      end subroutine net_set_fe56ec_fake_factor
      
      
      subroutine net_set_logTcut(handle, logTcut_lo, logTcut_lim, ierr)
         use net_def, only: do_net_set_logTcut
         integer, intent(in) :: handle
         real(dp), intent(in) :: logTcut_lo 
         real(dp), intent(in) :: logTcut_lim 
         integer, intent(out) :: ierr
         ierr = 0
         call do_net_set_logTcut(handle, logTcut_lo, logTcut_lim, ierr)
         if (ierr /= 0) return
      end subroutine net_set_logTcut
      
      
      subroutine net_set_eps_nuc_cancel(handle, logT_lo_eps_nuc_cancel, logT_hi_eps_nuc_cancel, ierr)
         use net_def, only: do_net_set_eps_nuc_cancel
         integer, intent(in) :: handle
         real(dp), intent(in) :: logT_lo_eps_nuc_cancel, logT_hi_eps_nuc_cancel 
         integer, intent(out) :: ierr
         ierr = 0
         call do_net_set_eps_nuc_cancel(handle, logT_lo_eps_nuc_cancel, logT_hi_eps_nuc_cancel, ierr)
         if (ierr /= 0) return
      end subroutine net_set_eps_nuc_cancel


      ! call this after you have finished defining the net
      subroutine net_setup_tables(handle, cache_suffix, ierr)
         ! This routine fills in a Net_General_Info structure.
         use net_initialize, only : alloc_net_general_info
         use rates_lib, only: read_raw_rates_records
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         character (len=*), intent(in) :: cache_suffix
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call read_raw_rates_records(ierr)
         if (ierr /= 0) then
            write(*,*) 'read_raw_rates_records failed in net_setup_tables'
            return
         end if
         call alloc_net_general_info(handle, cache_suffix, ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_general_info failed in net_setup_tables'
            return
         end if
      end subroutine net_setup_tables

      
      ! general info about the net
      
      
      integer function net_num_isos(handle, ierr) ! total number in current net
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         net_num_isos = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_num_isos -- did you call alloc_net_handle?'
            return
         end if
         net_num_isos = g% num_isos
      end function net_num_isos
      
      integer function net_num_reactions(handle, ierr) ! total number of rates for net
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         net_num_reactions = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_num_reactions -- did you call alloc_net_handle?'
            return
         end if
         net_num_reactions = g% num_reactions
      end function net_num_reactions
      
      ! this routine supplies the required size for the work array needed by net_get.
      integer function net_work_size(handle, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         use net_initialize, only: work_size
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         net_work_size = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_work_size -- did you call alloc_net_handle?'
            return
         end if
         net_work_size = work_size(g)
      end function net_work_size
      
      subroutine get_chem_id_table(handle, num_isos, chem_id, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle, num_isos ! num_isos must be number of isos in current net
         integer, intent(out) :: chem_id(num_isos) ! maps net iso number to chem id
            ! index from 1 to num_isos in current net
            ! value is between 1 and num_chem_isos
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_chem_id_table -- did you call alloc_net_handle?'
            return
         end if
         if (num_isos /= g% num_isos) then
            ierr = -1
            return
         end if
         chem_id(1:num_isos) = g% chem_id(1:num_isos)
         ierr = 0
      end subroutine get_chem_id_table

      subroutine get_chem_id_table_ptr(handle, chem_id_ptr, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, pointer :: chem_id_ptr(:) ! maps net iso number to chem id
            ! index from 1 to num_isos in current net
            ! value is between 1 and num_chem_isos
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_chem_id_table_ptr -- did you call alloc_net_handle?'
            return
         end if
         chem_id_ptr => g% chem_id
      end subroutine get_chem_id_table_ptr
      
      subroutine get_net_iso_table(handle, net_iso_table, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, intent(out) :: net_iso_table(num_chem_isos) ! maps chem id to net iso number
            ! index from 1 to num_chem_isos
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and num_isos in current net
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_net_iso_table -- did you call alloc_net_handle?'
            return
         end if
         net_iso_table(1:num_chem_isos) = g% net_iso(1:num_chem_isos)
         ierr = 0
      end subroutine get_net_iso_table
      
      subroutine get_net_iso_table_ptr(handle, net_iso_ptr, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, pointer :: net_iso_ptr(:) ! maps chem id to net iso number
            ! index from 1 to num_chem_isos
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and num_isos in current net
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_net_iso_table_ptr -- did you call alloc_net_handle?'
            return
         end if
         net_iso_ptr => g% net_iso
         ierr = 0
      end subroutine get_net_iso_table_ptr

      subroutine get_reaction_id_table(handle, num_reactions, reaction_id, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle, num_reactions 
            ! num_reactions must be number of reactions in current net
         integer, intent(out) :: reaction_id(num_reactions) 
            ! maps net reaction number to reaction id
            ! index from 1 to num_reactions in current net
            ! value is between 1 and num_reactions
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_reaction_id_table -- did you call alloc_net_handle?'
            return
         end if
         if (num_reactions /= g% num_reactions) then
            ierr = -1
            return
         end if
         reaction_id(1:num_reactions) = g% reaction_id(1:num_reactions)
         ierr = 0
      end subroutine get_reaction_id_table
      
      subroutine get_reaction_id_table_ptr(handle, reaction_id_ptr, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, pointer :: reaction_id_ptr(:) ! maps net reaction number to reaction id
            ! index from 1 to num_reactions in current net
            ! value is between 1 and rates_reaction_id_max
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_reaction_id_table_ptr -- did you call alloc_net_handle?'
            return
         end if
         reaction_id_ptr => g% reaction_id
      end subroutine get_reaction_id_table_ptr
      
      subroutine get_net_reaction_table(handle, net_reaction_table, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         use rates_def, only: rates_reaction_id_max
         integer, intent(in) :: handle
         integer, intent(out) :: net_reaction_table(rates_reaction_id_max) 
            ! maps reaction id to net reaction number
            ! index from 1 to rates_reaction_id_max
            ! value is 0 if the reaction is not in the current net
            ! else is value between 1 and rates_reaction_id_max in current net
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) &
               'invalid handle for get_net_reaction_table -- did you call alloc_net_handle?'
            return
         end if
         net_reaction_table(1:rates_reaction_id_max) = g% net_reaction(1:rates_reaction_id_max)
         ierr = 0
      end subroutine get_net_reaction_table
      
      subroutine get_net_reaction_table_ptr(handle, net_reaction_ptr, ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, pointer :: net_reaction_ptr(:) 
            ! maps reaction id to net reaction number
            ! index from 1 to num_reactions
            ! value is 0 if the reaction is not in the current net
            ! else is value between 1 and num_reactions in current net
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) &
               'invalid handle for get_net_reaction_table_ptr -- did you call alloc_net_handle?'
            return
         end if
         net_reaction_ptr => g% net_reaction
         ierr = 0
      end subroutine get_net_reaction_table_ptr


      ! net evaluation routines
      
      subroutine net_get( &
            handle, just_dxdt, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)
         use chem_def, only: num_categories
         use net_eval, only: eval_net
         use net_def, only: Net_General_Info, Net_Info, get_net_ptr
            
         use rates_def, only: num_rvs
      
         ! provide T or logT or both (the code needs both, so pass 'em if you've got 'em!)
         ! same for Rho and logRho
      
         integer, intent(in) :: handle
         logical, intent(in) :: just_dxdt
         type (Net_Info), pointer:: n
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
            ! mean number free electrons per nucleon, assuming complete ionization
            ! d_dxdt_dx(i, j) is d_dxdt(i)_dx(j), 
            ! i.e., partial derivative of rate for i'th isotope wrt j'th isotope abundance
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
            ! this arg is only used for prot(e-nu)neut and neut(e+nu)prot.
            ! if your net doesn't include those, you can safely ignore this arg.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
            ! when rates are calculated, they are multiplied by the
            ! corresponding values in this array.
            ! rate_factors array is indexed by reaction number.
            ! use net_reaction_table to map reaction id to reaction number.
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         logical, intent(in) :: reuse_rate_raw, reuse_rate_screened ! if true. use given rate_screened

         real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning after subtract reaction neutrinos
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos) 
            ! partial derivatives wrt mass fractions
      
         real(dp), intent(inout) :: dxdt(:) ! (num_isos)
            ! rate of change of mass fractions caused by nuclear reactions
         real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)
            ! partial derivatives of rates wrt mass fractions
            
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
            ! eps_nuc subtotals for each reaction category

         real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions

         integer, intent(in) :: screening_mode 
            
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         
         integer, intent(out) :: ierr ! 0 means okay
                  
         integer(8) :: time0, time1
         type (Net_General_Info), pointer :: g
         real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs
         logical, pointer :: from_weaklib(:) ! ignore if null
         logical, parameter :: symbolic = .false.
         logical, parameter :: rates_only = .false.
         actual_Qs => null()
         actual_neuQs => null()
         from_weaklib => null()

         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_get -- did you call alloc_net_handle?'
            return
         end if
         
         if (g% doing_timing) then
            call system_clock(time0)
         else
            time0 = 0
         endif
         
         call eval_net( &
               n, g, rates_only, just_dxdt, num_isos, num_reactions, g% num_wk_reactions, &
               x, temp, log10temp, rho, log10rho,  &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
               screening_mode,  &
               eps_nuc_categories, eps_neu_total, &
               lwork, work, actual_Qs, actual_neuQs, from_weaklib, symbolic, &
               ierr)
         if (g% doing_timing) then
            call system_clock(time1)
            g% clock_net_get = g% clock_net_get + (time1 - time0)
         end if

      end subroutine net_get
      
      subroutine get_net_rate_ptrs(handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            ierr)
         use net_def, only: Net_General_Info, get_net_ptr
         use net_initialize, only: set_rate_ptrs
         integer, intent(in) :: handle
         real(dp), pointer, dimension(:) :: &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho
         integer, intent(in) :: lwork
         real(dp), pointer :: work(:)
         integer, intent(out) :: ierr
         integer :: i
         type (Net_General_Info), pointer  :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for get_net_rate_ptrs -- did you call alloc_net_handle?'
            return
         end if
         call set_rate_ptrs(g, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            i, ierr)
      end subroutine get_net_rate_ptrs
      
      
      subroutine net_get_rates_only( &
            handle, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            screening_mode,  &
            lwork, work, ierr)
         use chem_def, only: num_categories
         use net_eval, only: eval_net
         use net_def, only: Net_General_Info, Net_Info, get_net_ptr
         use rates_def, only: num_rvs
      
         ! provide T or logT or both (the code needs both, so pass 'em if you've got 'em!)
         ! same for Rho and logRho
      
         integer, intent(in) :: handle
         type (Net_Info), pointer:: n
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
            ! mean number free electrons per nucleon, assuming complete ionization
            ! d_dxdt_dx(i, j) is d_dxdt(i)_dx(j), 
            ! i.e., partial derivative of rate for i'th isotope wrt j'th isotope abundance
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
            ! this arg is only used for prot(e-nu)neut and neut(e+nu)prot.
            ! if your net doesn't include those, you can safely ignore this arg.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
            ! when rates are calculated, they are multiplied by the
            ! corresponding values in this array.
            ! rate_factors array is indexed by reaction number.
            ! use net_reaction_table to map reaction id to reaction number.
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
            
         ! rate_raw and rate_screened are described in the declaration of the Net_Info derived type

         integer, intent(in) :: screening_mode
            
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         
         integer, intent(out) :: ierr ! 0 means okay
                  
         integer(8) :: time0, time1
         type (Net_General_Info), pointer :: g
         real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs
         logical, pointer :: from_weaklib(:) ! ignore if null
         logical, parameter :: symbolic = .false., &
            reuse_rate_raw = .false., reuse_rate_screened = .false.

         real(dp) :: eps_nuc, d_eps_nuc_dT, d_eps_nuc_dRho, eps_neu_total
         
         real(dp), target :: empty_array1(0), empty_array2(0,0)
         real(dp), pointer, dimension(:) :: &
            d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT, eps_nuc_categories
         real(dp), pointer, dimension(:,:) :: d_dxdt_dx
         logical, parameter :: rates_only = .true.
         logical, parameter :: just_dxdt = .true.

         d_eps_nuc_dx => empty_array1
         dxdt => empty_array1
         d_dxdt_dRho => empty_array1
         d_dxdt_dT => empty_array1
         
         d_dxdt_dx => empty_array2
         eps_nuc_categories => empty_array1

         actual_Qs => null()
         actual_neuQs => null()
         from_weaklib => null()

         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_get -- did you call alloc_net_handle?'
            return
         end if
         
         if (g% doing_timing) then
            call system_clock(time0)
         else
            time0 = 0
         endif
         
         call eval_net( &
               n, g, rates_only, just_dxdt, num_isos, num_reactions, g% num_wk_reactions, &
               x, temp, log10temp, rho, log10rho,  &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
               screening_mode,  &
               eps_nuc_categories, eps_neu_total, &
               lwork, work, actual_Qs, actual_neuQs, from_weaklib, symbolic, &
               ierr)
         if (g% doing_timing) then
            call system_clock(time1)
            g% clock_net_get = g% clock_net_get + (time1 - time0)
         end if
         
      end subroutine net_get_rates_only
      
      
      ! this sets d_dxdt_dx to 1 in locations where can have a nonzero partial
      ! it doesn't set other things such as eps_nuc or rates.
      ! takes the same set of args as net_get even though doesn't use them all.
      subroutine net_get_symbolic_d_dxdt_dx( &
            handle, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho, &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)
         use chem_def, only: num_categories
         use net_eval, only: eval_net
         use net_def, only: Net_General_Info, Net_Info, get_net_ptr
         use rates_def, only: num_rvs
      
         integer, intent(in) :: handle
         type (Net_Info), pointer :: n
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
            ! mean number free electrons per nucleon, assuming complete ionization
            ! d_dxdt_dx(i, j) is d_dxdt(i)_dx(j), 
            ! i.e., partial derivative of rate for i'th isotope wrt j'th isotope abundance
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
            ! this arg is only used for prot(e-nu)neut and neut(e+nu)prot.
            ! if your net doesn't include those, you can safely ignore this arg.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
            ! when rates are calculated, they are multiplied by the
            ! corresponding values in this array.
            ! rate_factors array is indexed by reaction number.
            ! use net_reaction_table to map reaction id to reaction number.
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)

         real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning after subtract reaction neutrinos
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos) 
            ! partial derivatives wrt mass fractions
      
         real(dp), intent(inout) :: dxdt(:) ! (num_isos)
            ! rate of change of mass fractions caused by nuclear reactions
         real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)
            ! partial derivatives of rates wrt mass fractions
            
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
            ! eps_nuc subtotals for each reaction category

         real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions
            
         ! rate_raw and rate_screened are described in the declaration of the Net_Info derived type

         integer, intent(in) :: screening_mode ! Selects which screening mode to use, see rates_def for definition 
            
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         
         integer, intent(out) :: ierr
            ! ierr = 0 means AOK
            ! ierr = -1 means mass fractions don't add to something very close to 1.0
            ! ierr = -2 means neither T nor logT were provided
            ! ierr = -3 means neither Rho nor logRho were provided
                  
         type (Net_General_Info), pointer :: g
         real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs
         logical, pointer :: from_weaklib(:) ! ignore if null
         logical, parameter :: symbolic = .true.
         integer :: num_rates_reduced
         real(dp) :: max_old_rate_div_new_rate
         logical, parameter :: rates_only = .false.
         logical, parameter :: just_dxdt = .false.
         
         actual_Qs => null()
         actual_neuQs => null()
         from_weaklib => null()

         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,'(a)')  &
               'invalid handle for net_get_symbolic_d_dxdt_dx -- did you call alloc_net_handle?'
            return
         end if
         
         call eval_net( &
            n, g, rates_only, just_dxdt, num_isos, num_reactions, g% num_wk_reactions, &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, .false., .false., &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, actual_Qs, actual_neuQs, from_weaklib, symbolic, &
            ierr)
         
      end subroutine net_get_symbolic_d_dxdt_dx
      
      subroutine net_get_with_Qs( &
            handle, just_dxdt, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, actual_Qs, actual_neuQs, from_weaklib, &
            ierr)
         use chem_def, only: num_categories
         use net_eval, only: eval_net
         use net_def, only: Net_General_Info, Net_Info, get_net_ptr
         use rates_def, only: num_rvs
         integer, intent(in) :: handle
         logical, intent(in) :: just_dxdt
         type (Net_Info), pointer:: n
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         logical, intent(in) :: reuse_rate_raw, reuse_rate_screened
         real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning after subtract reaction neutrinos
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos) 
         real(dp), intent(inout) :: dxdt(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)            
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions
         integer, intent(in) :: screening_mode
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs ! ignore if null  (num_reactions)
         logical, pointer :: from_weaklib(:) ! ignore if null
         integer, intent(out) :: ierr
                  
         logical, parameter :: rates_only = .false.
         logical, parameter :: symbolic = .false.
         integer(8) :: time0, time1
         type (Net_General_Info), pointer :: g

         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_get_with_Qs -- did you call alloc_net_handle?'
            return
         end if

         if (g% doing_timing) then
            call system_clock(time0)
         else
            time0 = 0
         endif
         
         call eval_net( &
            n, g, rates_only, just_dxdt, num_isos, num_reactions, g% num_wk_reactions, &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, actual_Qs, actual_neuQs, from_weaklib, symbolic, &
            ierr)
         if (g% doing_timing) then
            call system_clock(time1)
            g% clock_net_get = g% clock_net_get + (time1 - time0)
         end if
         
      end subroutine net_get_with_Qs

      integer function net_1_zone_burn_work_size(handle,ierr) result(sz)
         use net_burn, only: burn_1_zone_work_size
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_1_zone_burn_work_size -- did you call alloc_net_handle?'
            sz = 0
            return
         end if
         sz = burn_1_zone_work_size(g)
      end function net_1_zone_burn_work_size

      ! a 1-zone integrator for nets -- for given temperature and density as functions of time
      subroutine net_1_zone_burn( &
            net_handle, eos_handle, num_isos, num_reactions, t_start, t_end, starting_x, &
            num_times_for_interpolation, times, log10Ts_f1, log10Rhos_f1, etas_f1, &
            dxdt_source_term, rate_factors, &
            weak_rate_factor, reaction_Qs, reaction_neuQs, &
            screening_mode, &
            stptry, max_steps, eps, odescal, &
            okay_to_reuse_rate_screened, &
            use_pivoting, trace, dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, &
            net_lwork, net_work_array, &
            ! results
            ending_x, eps_nuc_categories, avg_eps_nuc, eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
         use net_burn, only: burn_1_zone
         use net_def
         use chem_def, only: num_categories
         
         integer, intent(in) :: net_handle, eos_handle
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in) :: t_start, t_end, starting_x(:) ! (num_isos)
         
         integer, intent(in) :: num_times_for_interpolation 
            ! ending time is times(num_times); starting time is 0
         real(dp), pointer, intent(in) :: times(:) ! (num_times) 
         real(dp), pointer, intent(in) :: log10Ts_f1(:) ! =(4,numtimes) interpolant for log10T(time)
         real(dp), pointer, intent(in) :: log10Rhos_f1(:) ! =(4,numtimes) interpolant for log10Rho(time)
         real(dp), pointer, intent(in) :: etas_f1(:) ! =(4,numtimes) interpolant for eta(time)
         real(dp), pointer, intent(in) :: dxdt_source_term(:) ! (num_isos)  or null if no source term.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode ! see screen_def
         real(dp), intent(in) :: stptry ! try this for 1st step.  0 means try in 1 step.
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         real(dp), intent(in) :: eps, odescal ! tolerances.  e.g., set both to 1d-6
         logical, intent(in) :: okay_to_reuse_rate_screened
            ! this flag should be false if there will be large changes in composition.
            ! if the composition changes will be small, then can gain efficiency
            ! by only evaluating the screening factors for the starting composition.
            ! reuse_rate_raw should be false unless reuse_rate_screened is true, and
            ! there is no change in temperature or density since the last call.
            ! in both cases, the burn_work_array must be the same as used previously.
         logical, intent(in) :: use_pivoting ! for matrix solves
         logical, intent(in) :: trace, dbg
         interface
            include 'burner_finish_substep.inc'
         end interface
         integer, intent(in) :: net_lwork, burn_lwork
         real(dp), intent(inout), pointer :: burn_work_array(:) ! (burn_lwork)
         real(dp), intent(inout), pointer :: net_work_array(:) ! (net_lwork)
         real(dp), intent(inout) :: ending_x(:) ! (num_isos)
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: avg_eps_nuc, eps_neu_total
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         integer, intent(out) :: ierr
         
         call burn_1_zone( &
            net_handle, eos_handle, num_isos, num_reactions, t_start, t_end, starting_x, &
            num_times_for_interpolation, times, log10Ts_f1, log10Rhos_f1, etas_f1, &
            dxdt_source_term, rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, screening_mode, &
            stptry, max_steps, eps, odescal, &
            okay_to_reuse_rate_screened, &
            use_pivoting, trace, dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, &
            net_lwork, net_work_array, &
            ending_x, eps_nuc_categories, avg_eps_nuc, eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
         
      end subroutine net_1_zone_burn
      
      subroutine get_burn_work_array_pointers(handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            burn_lwork, burn_work_array, ierr)
         integer, intent(in) :: handle
         real(dp), pointer :: burn_work_array(:)
         integer, intent(in) :: burn_lwork
         real(dp), pointer, dimension(:) :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         integer, intent(out) :: ierr
         call get_net_rate_ptrs(handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, burn_lwork, burn_work_array, &
            ierr)
      end subroutine get_burn_work_array_pointers



      integer function net_1_zone_burn_const_density_work_size(handle,ierr) result(sz)
         use net_burn_const_density, only: burn_const_density_1_zone_work_size
         use net_def, only: Net_General_Info, get_net_ptr
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_1_zone_burn_const_density_work_size'
            sz = 0
            return
         end if
         sz = burn_const_density_1_zone_work_size(g)
      end function net_1_zone_burn_const_density_work_size


      ! a 1-zone integrator for nets -- for given density
      ! evolve lnT according to dlnT/dt = eps_nuc/(Cv*T)
      subroutine net_1_zone_burn_const_density( &
            net_handle, eos_handle, num_isos, num_reactions, t_start, t_end, &
            starting_x, starting_log10T, log10Rho, &
            get_eos_info_for_burn_at_const_density, &
            rate_factors, weak_rate_factor, reaction_Qs, reaction_neuQs, &
            screening_mode,  &
            stptry, max_steps, eps, odescal, &
            use_pivoting, trace, dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, net_lwork, net_work_array, &
            ! results
            ending_x, eps_nuc_categories, ending_log10T, &
            avg_eps_nuc, ending_eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
         use net_burn_const_density, only: burn_const_density_1_zone
         use net_def
         use chem_def, only: num_categories
         
         integer, intent(in) :: net_handle, eos_handle, num_isos, num_reactions
         real(dp), intent(in) :: t_start, t_end, starting_x(:) ! (num_isos)
         real(dp), intent(in) :: starting_log10T, log10Rho
         interface
            include 'burner_const_density_get_eos_info.inc'
         end interface
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode ! see screen_def
         real(dp), intent(in) :: stptry ! try this for 1st step.  0 means try in 1 step.
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         real(dp), intent(in) :: eps, odescal ! tolerances.  e.g., set both to 1d-6
         logical, intent(in) :: use_pivoting ! for matrix solves
         logical, intent(in) :: trace, dbg
         interface
            include 'burner_finish_substep.inc'
         end interface
         integer, intent(in) :: net_lwork, burn_lwork
         real(dp), intent(inout), pointer :: burn_work_array(:) ! (burn_lwork)
         real(dp), intent(inout), pointer :: net_work_array(:) ! (net_lwork)
         real(dp), intent(inout) :: ending_x(:) ! (num_isos)
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: ending_log10T, avg_eps_nuc, ending_eps_neu_total
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         integer, intent(out) :: ierr
         
         call burn_const_density_1_zone( &
            net_handle, eos_handle, num_isos, num_isos+1, num_reactions, t_start, t_end, &
            starting_x, starting_log10T, log10Rho, &
            get_eos_info_for_burn_at_const_density, &
            rate_factors, weak_rate_factor, reaction_Qs, reaction_neuQs, &
            screening_mode,  &
            stptry, max_steps, eps, odescal, &
            use_pivoting, trace, dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, net_lwork, net_work_array, &
            ending_x, eps_nuc_categories, ending_log10T, avg_eps_nuc, ending_eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
         
      end subroutine net_1_zone_burn_const_density
      
      subroutine get_burn_const_density_work_array_pointers(handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            burn_lwork, burn_work_array, ierr)
         integer, intent(in) :: handle
         real(dp), pointer :: burn_work_array(:)
         integer, intent(in) :: burn_lwork
         real(dp), pointer, dimension(:) :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         integer, intent(out) :: ierr
         call get_net_rate_ptrs(handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, burn_lwork, burn_work_array, &
            ierr)
      end subroutine get_burn_const_density_work_array_pointers

      
      ! evolve T according to dT/dt = eps_nuc/Cp while using given P.
      !  then find new Rho and Cp to match P and new T.
      subroutine net_1_zone_burn_const_P( &
            net_handle, eos_handle, num_isos, num_reactions,  &
            which_solver, starting_temp, starting_x, clip, &
            ! for interpolating log10P wrt time &
            num_times_for_interpolation, times, log10Ps_f1, &
            ! other args for net_get &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, screening_mode,  &
            ! args to control the solver &
            h, max_step_size, max_steps, rtol, atol, itol, x_min, x_max, which_decsol,  &
            ! results &
            caller_id, solout, iout,  &
            ending_x, ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS, &
            nfcn, njac, nstep, naccpt, nrejct, time_doing_net, time_doing_eos, ierr)
         use net_burn_const_P, only: burn_1_zone_const_P
         use chem_def, only: num_categories
         use rates_def, only: num_rvs
         
         integer, intent(in) :: net_handle, eos_handle
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), pointer, intent(in) :: starting_x(:) ! (num_isos)
         real(dp), intent(in) :: starting_temp
         logical, intent(in) :: clip ! if true, set negative x's to zero during burn.
         
         integer, intent(in) :: which_solver ! as defined in num_def.f
         integer, intent(in) :: num_times_for_interpolation ! ending time is times(num_times); starting time is 0
         real(dp), pointer, intent(in) :: times(:) ! (num_times) 
         real(dp), pointer, intent(in) :: log10Ps_f1(:) ! =(4,numtimes) interpolant for log10P(time)

         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode
         
         ! args to control the solver -- see num/public/num_isolve.dek
         real(dp), intent(inout) :: h 
         real(dp), intent(in) :: max_step_size ! maximal step size.
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         ! absolute and relative error tolerances
         real(dp), pointer :: rtol(:) ! relative error tolerance (num_isos)
         real(dp), pointer :: atol(:) ! absolute error tolerance (num_isos)
         integer, intent(in) :: itol ! switch for rtol and atol
         real(dp), intent(in) :: x_min, x_max ! bounds on allowed values
         integer, intent(in) :: which_decsol ! from mtx_def
         integer, intent(in) :: caller_id
         interface ! subroutine called after each successful step
            include "num_solout.dek"
         end interface
         integer, intent(in)  :: iout
         real(dp), intent(inout), pointer :: ending_x(:)
         real(dp), intent(out) :: ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         real(dp), intent(inout) :: time_doing_net
            ! if < 0, then ignore
            ! else on return has input value plus time spent doing eval_net
         real(dp), intent(inout) :: time_doing_eos
            ! if < 0, then ignore
            ! else on return has input value plus time spent doing eos
         integer, intent(out) :: ierr
         
         call burn_1_zone_const_P( &
            net_handle, eos_handle, num_isos, num_reactions,  &
            which_solver, starting_temp, starting_x, clip, &
            num_times_for_interpolation, times, log10Ps_f1, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, screening_mode,  &
            h, max_step_size, max_steps, rtol, atol, itol, x_min, x_max, which_decsol,  &
            caller_id, solout, iout, &
            ending_x, ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS, &
            nfcn, njac, nstep, naccpt, nrejct, time_doing_net, time_doing_eos, ierr)
         
      end subroutine net_1_zone_burn_const_P
      
      ! approximate beta decay neutrino energies (in MeV)
      ! Fowler, Caughlan, Zimmerman, Annual Review Astro. Astrophys., 1975.12:69-112. eqn (1).
      real(dp) function eval_neutrino_Q(i1, i2)
         use net_initialize, only:neutrino_Q
         integer, intent(in) :: i1, i2 ! i1 decays to i2.  e.g., i1=in13 and i2=ic13
         eval_neutrino_Q = neutrino_Q(i1, i2)
      end function eval_neutrino_Q      
      
      
      ! for calculating reaction Q
      real(dp) function isoB(ci)
         use chem_def, only: del_Mp, del_Mn
         integer, intent(in) :: ci
         isoB = chem_isos% binding_energy(ci) - chem_isos% Z(ci)*del_Mp - chem_isos% N(ci)*del_Mn
      end function isoB
      
      
      subroutine clean_up_fractions(nzlo, nzhi, species, nz, xa, max_sum_abs, xsum_tol, ierr)
         ! make sure all fractions are okay and sum to 1.0
         use net_eval, only: do_clean_up_fractions
         integer, intent(in) :: nzlo, nzhi, species, nz
         real(dp), intent(inout) :: xa(:,:) ! (species, nz) ! mass fractions
         real(dp), intent(in) :: max_sum_abs
         ! if any k has sum(abs(xa(:,k))) greater than this, set ierr = -1 and return.
         ! else clip each element abundance to the range 0 to 1 and continue.
         real(dp), intent(in) :: xsum_tol
         ! if any sum of abundances is now different from 1 by more than this, set ierr = -1
         ! otherwise, rescale the abundances so that they sum to 1.
         integer, intent(out) :: ierr
         call do_clean_up_fractions(nzlo, nzhi, species, nz, xa, max_sum_abs, xsum_tol, ierr)
      end subroutine clean_up_fractions
      

      subroutine clean1(species, xa, max_sum_abs, xsum_tol, ierr)
         use net_eval, only: do_clean1
         integer, intent(in) :: species
         real(dp), intent(inout) :: xa(:) ! (species)
         real(dp), intent(in) :: max_sum_abs, xsum_tol
         integer, intent(out) :: ierr
         call do_clean1(species, xa, 1, max_sum_abs, xsum_tol, ierr)
      end subroutine clean1
      
      
      end module net_lib

