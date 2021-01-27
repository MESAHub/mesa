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

      module net_initialize
      use net_def
      use net_screen
      use chem_def
      use math_lib

      use net_approx21, only : num_reactions_func => num_reactions, &
                               num_mesa_reactions_func => num_mesa_reactions
      
      implicit none


         
         
      integer, parameter :: max_num_special_case_reactants = 5
      integer, parameter :: max_num_special_case_reactions = 70
      integer :: num_special_case_reactions
      integer :: special_case_reactants( &
         max_num_special_case_reactants, max_num_special_case_reactions)
      character (len=maxlen_reaction_Name) :: &
         special_case_reactions(2,max_num_special_case_reactions)


      contains

      integer function work_size(g)
         type (Net_General_Info), pointer  :: g         
         integer :: num_isos, num_reactions, num_wk_reactions
         num_reactions = g% num_reactions
         if (g% doing_approx21) then
            num_reactions  = num_reactions_func(g% add_co56_to_approx21)
         end if
         num_isos = g% num_isos
         num_wk_reactions = g% num_wk_reactions
         work_size = &
            num_isos*(2+num_isos) + &
            num_weak_info_arrays_in_Net_Info*num_wk_reactions + &
            6*num_reactions
               ! ratraw, dratrawdt, dratrawdd
               ! ratdum, dratdumdt, dratdumdd
         if (g% doing_approx21) &
            work_size = work_size + num_isos*num_isos + & ! dfdy
               2*num_reactions + & ! dratdumdy1, dratdumdy2
               5*num_isos ! d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho
      end function work_size

      subroutine set_ptrs_for_approx21( &
            add_co56, i, work, dfdy, dratdumdy1, dratdumdy2, &
            d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho)

         logical, intent(in) :: add_co56
         integer, intent(inout) :: i
         real(dp), pointer :: work(:) ! (lwork)
         real(dp), pointer :: dfdy(:,:)
         real(dp), dimension(:), pointer :: &
            dratdumdy1, dratdumdy2, &
            d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho
            
         integer :: num_isos, num_reactions

         num_reactions = num_reactions_func(add_co56)
         if (add_co56) then
            num_isos = 22
         else
            num_isos = 21
         end if
         
         dfdy(1:num_isos,1:num_isos) => work(i:i+num_isos*num_isos-1)
         i=i+num_isos*num_isos
      
         dratdumdy1(1:num_reactions) => work(i:i+num_reactions-1)
         i=i+num_reactions
      
         dratdumdy2(1:num_reactions) => work(i:i+num_reactions-1)
         i=i+num_reactions
      
         d_epsnuc_dy(1:num_isos) => work(i:i+num_isos-1)
         i=i+num_isos
      
         d_epsneu_dy(1:num_isos) => work(i:i+num_isos-1)
         i=i+num_isos
      
         dydt1(1:num_isos) => work(i:i+num_isos-1)
         i=i+num_isos
      
         dfdT(1:num_isos) => work(i:i+num_isos-1)
         i=i+num_isos
      
         dfdRho(1:num_isos) => work(i:i+num_isos-1)
         i=i+num_isos
         
      end subroutine set_ptrs_for_approx21

      subroutine set_rate_ptrs(g, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            i, ierr)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: lwork
         integer, intent(inout) :: i
         real(dp), pointer, dimension(:) :: work, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho
         integer, intent(out) :: ierr
         integer :: sz
         ierr = 0
         if (g% doing_approx21) then
            sz = num_reactions_func(g% add_co56_to_approx21)
         else
            sz = g% num_reactions
         end if
         if (6*sz > lwork) then
            ierr = -1
            return
         end if
         i = 0         
         rate_screened(1:sz) => work(i+1:i+sz); i=i+sz      
         rate_screened_dT(1:sz) => work(i+1:i+sz); i=i+sz      
         rate_screened_dRho(1:sz) => work(i+1:i+sz); i=i+sz      
         rate_raw(1:sz) => work(i+1:i+sz); i=i+sz      
         rate_raw_dT(1:sz) => work(i+1:i+sz); i=i+sz      
         rate_raw_dRho(1:sz) => work(i+1:i+sz); i=i+sz      
      end subroutine set_rate_ptrs
      
      
      subroutine setup_net_info( &
            g, n, eps_nuc_categories,  &
            screening_mode,  &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            reuse_rate_raw, reuse_rate_screened, &
            i, ierr)
         use chem_def
         type (Net_General_Info), pointer  :: g
         type (Net_Info), pointer :: n
         integer, intent(in) :: lwork
         real(dp), intent(inout), target :: eps_nuc_categories(:) ! (num_categories)
         integer, intent(in) :: screening_mode
         real(dp), intent(inout), dimension(:), target :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         real(dp), pointer :: work(:) ! (lwork)
         logical, intent(in) :: reuse_rate_raw, reuse_rate_screened
         integer, intent(inout) :: i
         integer, intent(out) :: ierr
         
         integer :: num_reactions, num_isos, num_wk_reactions
         character (len=256) :: err_msg    
         real(dp) :: qneu
         
         include 'formats'
         
         ierr = 0
         n% g => g
         num_isos = g% num_isos
         num_wk_reactions = g% num_wk_reactions
         if (g% doing_approx21) then
            num_reactions = num_reactions_func(g% add_co56_to_approx21)
         else
            num_reactions = g% num_reactions
         end if
         
         n% screening_mode = screening_mode
 
         n% eps_nuc_categories => eps_nuc_categories
         
         if (.not. reuse_rate_screened) then
            rate_screened(:) = 0
            rate_screened_dT(:) = 0
            rate_screened_dRho(:) = 0
         end if
         n% rate_screened => rate_screened
         n% rate_screened_dT => rate_screened_dT
         n% rate_screened_dRho => rate_screened_dRho
         
         if (.not. reuse_rate_raw) then
            rate_raw(:) = 0
            rate_raw_dT(:) = 0
            rate_raw_dRho(:) = 0
         end if
         
         n% rate_raw => rate_raw
         n% rate_raw_dT => rate_raw_dT
         n% rate_raw_dRho => rate_raw_dRho
         
         ! input value of i is amount already used
         n% y => work(i+1:i+num_isos); i=i+num_isos
         
         n% d_eps_nuc_dy => work(i+1:i+num_isos); i=i+num_isos
         
         n% d_dydt_dy(1:num_isos,1:num_isos) => work(i+1:i+num_isos*num_isos)
         i = i+num_isos*num_isos 
         
         n% lambda(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         n% dlambda_dlnT(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         n% dlambda_dlnRho(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         
         n% Q(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         n% dQ_dlnT(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         n% dQ_dlnRho(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         
         n% Qneu(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         n% dQneu_dlnT(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions
         n% dQneu_dlnRho(1:num_wk_reactions) => work(i+1:i+num_wk_reactions)
         i=i+num_wk_reactions

      end subroutine setup_net_info
      

      subroutine alloc_net_general_info(handle, cache_suffix, ierr)
         use rates_def, only: extended_screening
         use rates_lib, only: make_rate_tables
         
         integer, intent(in) :: handle
         character (len=*), intent(in) :: cache_suffix
         integer, intent(out) :: ierr
         
         type (Net_Info), target :: netinfo 
            ! just used during initialization and then discarded
         type (Net_Info), pointer :: n
         integer :: ios, status, lwork, num_reactions, &
            num_isos, num_wk_reactions, i, iwork
         real(dp), dimension(:), pointer :: eps_nuc_categories
         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         real(dp), pointer :: work(:)
         type (Net_General_Info), pointer  :: g
         
         ! the following not used during initialization, but need as args.
         logical, parameter :: reuse_rate_raw = .false.
         logical, parameter :: reuse_rate_screened = .false.
         integer, parameter :: screening_mode = extended_screening
         
         include 'formats'
         
         ierr = 0

         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return

         netinfo % g => g
         n => netinfo
         
         n% reaction_Qs => std_reaction_Qs
         n% reaction_neuQs => std_reaction_neuQs
                  
         num_reactions = g% num_reactions
         num_isos = g% num_isos
         num_wk_reactions = g% num_wk_reactions
         lwork = work_size(g)
         
         allocate(work(lwork), eps_nuc_categories(num_categories), stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_general_info failed in allocate'
            return
         end if
         
         eps_nuc_categories = 0

         call set_rate_ptrs(g, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            iwork, ierr) ! iwork is number of entries in work used for rates
         if (ierr /= 0) then
            write(*,*) 'failed in set_ptrs_in_work'
            return
         end if
         
         call setup_net_info( &
            g, n, eps_nuc_categories,  &
            screening_mode, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            reuse_rate_raw, reuse_rate_screened, &
            iwork, ierr) ! iwork updated for amount now used in work
         if (ierr /= 0) then
            write(*,*) 'alloc_net_general_info failed in setup_net_info'
            return
         end if
         
         g% cache_suffix = trim(cache_suffix)
         
         call make_rate_tables( &
            g% num_reactions, g% cache_suffix, g% reaction_id, g% which_rates,  &
            g% rate_table, g% rattab_f1, nrattab, g% ttab, g% logttab, ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_general_info failed in call on make_rate_tables'
            return
         end if
         
         call make_screening_tables(n, ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_general_info failed in make_screening_tables'
            return
         end if
         
         deallocate(work, eps_nuc_categories)
         
      end subroutine alloc_net_general_info
      
         
      subroutine set_reaction_max_Z(g)
         type (Net_General_Info), pointer  :: g
         integer :: i, ir, max_Z, max_Z_plus_N
         include 'formats'
         do i=1, g% num_reactions
            ir = g% reaction_id(i)
            max_Z = 0
            max_Z_plus_N = 0
            call update_max_Z(reaction_inputs(2,ir))
            call update_max_Z(reaction_inputs(4,ir))
            call update_max_Z(reaction_inputs(6,ir))
            call update_max_Z(reaction_outputs(2,ir))
            call update_max_Z(reaction_outputs(4,ir))
            call update_max_Z(reaction_outputs(6,ir))
            g% reaction_max_Z(i) = max_Z
            g% reaction_max_Z_plus_N_for_max_Z(i) = max_Z_plus_N
            !write(*,3) trim(reaction_name(ir)), max_Z, max_Z_plus_N
         end do
         
         contains
         
         subroutine update_max_Z(iso)
            use chem_def, only: chem_isos
            integer, intent(in) :: iso
            integer :: Z, Z_plus_N
            if (iso == 0) return
            Z = chem_isos% Z(iso)
            if (Z < max_Z) return
            Z_plus_N = chem_isos% Z_plus_N(iso)
            if (Z > max_Z) then
               max_Z = Z
               max_Z_plus_N = Z_plus_N
            else if (Z_plus_N > max_Z_plus_N) then
               max_Z_plus_N = Z_plus_N
            end if
         end subroutine update_max_Z         
         
      end subroutine set_reaction_max_Z         
      
            
      recursive subroutine do_read_net_file(net_filename, handle, ierr)
         use utils_def
         use utils_lib
         use chem_lib, only: chem_get_iso_id
         use chem_def, only: chem_isos, ih1, ihe4, ineut
         character (len=*), intent(in) :: net_filename
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         
         integer :: iounit, n, i, j, k, t, id, h1, he4, neut
         character (len=256) :: buffer, string, filename
         logical, parameter :: dbg = .false.
         type (Net_General_Info), pointer  :: g
         
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (len_trim(g% net_filename) == 0) g% net_filename = trim(net_filename)
         
         ! first look in local directory
         filename = trim(net_filename)
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! look in local nets directory
            filename = 'nets/' // trim(net_filename)
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then ! look in global nets directory           
               filename = trim(net_dir) // '/nets/' // trim(net_filename)
               ierr = 0
               open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed to open net file ' // trim(filename)
                  return
               end if
            end if
         end if
         
         if (dbg) then
            write(*,*)
            write(*,*) 'read_net_file <' // trim(filename) // '>'
            write(*,*) 'net_filename <' // trim(net_filename) // '>'
            write(*,*)
         end if

         n = 0
         i = 0
         
         do
            t = token(iounit, n, i, buffer, string)
            select case(t)
               case(name_token)
                  select case(string)
                  
                     case ('add_isos')
                        call do_isos(.true., ierr)
                        if (ierr /= 0) return
                  
                     case ('add_iso')
                        call do_isos(.true., ierr)
                        if (ierr /= 0) return
                        
                     case ('remove_isos')
                        call do_isos(.false., ierr)
                        if (ierr /= 0) return
                        
                     case ('remove_iso')
                        call do_isos(.false., ierr)
                        if (ierr /= 0) return
                        
                     case ('add_reaction')
                        call do_reactions(.true., ierr)
                        if (ierr /= 0) return
                        
                     case ('add_reactions')
                        call do_reactions(.true., ierr)
                        if (ierr /= 0) return
                        
                     case ('remove_reaction')
                        call do_reactions(.false., ierr)
                        if (ierr /= 0) return
                        
                     case ('remove_reactions')
                        call do_reactions(.false., ierr)
                        if (ierr /= 0) return
                        
                     case ('add_iso_and_reactions')
                        call do_basic_reactions_for_isos(ierr)
                        if (ierr /= 0) return
                        
                     case ('add_isos_and_reactions')
                        call do_basic_reactions_for_isos(ierr)
                        if (ierr /= 0) return
                        
                     case ('approx21')
                        call do_approx21(1,ierr)
                        if (ierr /= 0) return
                        
                     case ('approx21_plus_co56')
                        call do_approx21(2,ierr)
                        if (ierr /= 0) return
                        
                     case ('include')
                        t = token(iounit, n, i, buffer, string)
                        if (t /= string_token) then
                           call error; return
                        end if
                        call do_read_net_file(string, handle, ierr)
                        if (ierr /= 0) return
                        
                     case default
                        call error; return
                        
                  end select
               case(eof_token)
                  exit
               case default
                  call error; return
            end select
            
         end do
         
         close(iounit)


! network veracity checks go here

! fxt check that al26 and al26-1 or al26-2 are not both specified
           i = chem_get_iso_id('al26')
!           if (dbg) write(6,*) i, trim(chem_isos% name(i)), g% net_iso(i)
           j = chem_get_iso_id('al26-1')
!           if (dbg) write(6,*) j, trim(chem_isos% name(j)), g% net_iso(j)
           k = chem_get_iso_id('al26-2')
!           if (dbg) write(6,*) k, trim(chem_isos% name(k)), g% net_iso(k)

           if ( (g% net_iso(i) == 1 .and. g% net_iso(j) == 1) .or. & 
                (g% net_iso(i) == 1 .and. g% net_iso(k) == 1)) then
            string = 'cannot specify al26 and al26-1 or al26-2'
            call error ; return
           end if

! fxt check that both al26-1 and al26-2 is specified if one or the other is given
           if (g% net_iso(j) == 1  .and. g% net_iso(k) == 0) then
            string = 'must specify al26-2 if al26-1 is set'
            call error ; return
           end if
           if (g% net_iso(k) == 1  .and. g% net_iso(j) == 0) then
            string = 'must specify al26-1 if al26-2 is set'
            call error ; return
           end if

! done with network veracity checks
         
         if (dbg) then
            write(*,*)
            write(*,*) 'done read_net_file ' // trim(filename)
            write(*,*)
         end if
         
         
         contains
         
         
         subroutine error
            character (len=256) :: message
            ierr = -1
            write(*,'(a)') ' problem reading net file ' // trim(filename) &
                  // ' error somewhere around here <' // trim(string) // &
                  '>.  please check for missing comma or other typo.'
            close(iounit)
         end subroutine error     
             
         
         subroutine do_approx21(which_case, ierr) ! e.g. approx21(cr56)
            use chem_lib, only: chem_get_iso_id
            use chem_def, only: chem_isos
            integer, intent(in) :: which_case
            integer, intent(out) :: ierr
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
            t = token(iounit, n, i, buffer, string)
            if (t /= name_token) then
               call error; return
            end if
            id = chem_get_iso_id(string)
            if (id <= 0) then
               call error; return
            end if
            t = token(iounit, n, i, buffer, string)
            if (t /= right_paren_token) then
               call error; return
            end if
            g% approx21_ye_iso = id
            g% fe56ec_n_neut = &
               (chem_isos% N(id) - chem_isos% N(ife56)) - &
               (chem_isos% Z(ife56) - chem_isos% Z(id))
            g% doing_approx21 = .true.
            g% add_co56_to_approx21 = (which_case == 2) 
         end subroutine do_approx21
             
         
         subroutine do_isos(add_flag, ierr)
            use chem_lib, only: chem_get_element_id
            use chem_lib, only: lookup_ZN
            logical, intent(in) :: add_flag
            integer, intent(out) :: ierr
            logical :: have_next_token
            integer :: A, A1, A2, Z
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
            have_next_token = .false.
         iso_loop: do
               if (.not. have_next_token) t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  call error; return
               end if
               id = chem_get_iso_id(string)
               if (id > 0) then
                  if (add_flag) then
                     if (dbg) write(*,'(a30,i5,a30)') 'add_net_iso', id, trim(chem_isos% name(id))
                     call add_net_iso(handle, id, ierr)
                     if (ierr /= 0) then
                        write(*,*) 'failed in add_net_iso for ' // trim(chem_isos% name(id))
                        call error; return
                     end if
                  else
                     if (dbg) write(*,'(a30,i5,a30)') 'remove_net_iso', &
                        id, trim(chem_isos% name(id))
                     call remove_net_iso(handle, id, ierr)
                     if (ierr /= 0) then
                        write(*,*) 'failed in add_net_iso for ' // trim(chem_isos% name(id))
                        call error; return
                     end if
                  end if
               else
                  Z = chem_get_element_id(string)
                  if (Z < 0) then
                     read(string,fmt=*,iostat=ierr) Z
                     if (Z < 0) then
                        write(*,*) 'unknown name ' // trim(string)
                        call error; return
                     end if
                  end if
                  t = token(iounit, n, i, buffer, string)
                  !write(*,*) 'iso_loop 2nd token ' // trim(string)
                  if (t /= name_token) then
                     !write(*,*) 'e 2'
                     call error; return
                  end if
                  read(string,fmt=*,iostat=ierr) A1
                  if (ierr /= 0) then
                     !write(*,*) 'string <' // trim(string) // '>'
                     call error; return
                  end if
                  t = token(iounit, n, i, buffer, string)
                  !write(*,*) 'iso_loop 3rd token ' // trim(string)
                  if (t /= name_token) then
                     !write(*,*) 'e 3'
                     call error; return
                  end if
                  read(string,fmt=*,iostat=ierr) A2
                  if (ierr /= 0) then
                     !write(*,*) 'string <' // trim(string) // '>'
                     call error; return
                  end if
                  do A = A1, A2
                     id = lookup_ZN(Z, A-Z)
                     if (id <= 0) then
                        write(*,'(a30,5i5,a30)') 'bad iso for ' // trim(string), &
                           Z, A-Z, A, A1, A2
                        ierr = -1; call error; return
                        stop
                     end if
                     if (g% net_iso(id) == 0) then
                        if (add_flag) then
                           if (dbg) write(*,'(a30,i5,a30)') 'add_net_iso', &
                              id, trim(chem_isos% name(id))
                           call add_net_iso(handle, id, ierr)
                           if (ierr /= 0) then
                              write(*,*) 'failed in add_net_iso for ' // trim(chem_isos% name(id))
                              call error; return
                           end if
                        else
                           if (dbg) write(*,'(a30,i5,a30)') 'remove_net_iso', &
                              id, trim(chem_isos% name(id))
                           call remove_net_iso(handle, id, ierr)
                           if (ierr /= 0) then
                              write(*,*) 'failed in add_net_iso for ' // trim(chem_isos% name(id))
                              call error; return
                           end if
                        end if
                     
                     end if            
                  end do
               end if
               t = token(iounit, n, i, buffer, string)
               if (t == right_paren_token) exit iso_loop
               if (t /= comma_token) then
                  have_next_token = .true.
               else
                  have_next_token = .false.
               end if
            end do iso_loop
         end subroutine do_isos

         
         subroutine do_reactions(add_flag, ierr)
            use rates_lib, only: rates_reaction_id
            logical, intent(in) :: add_flag
            integer, intent(out) :: ierr
            integer :: cnt, ir
            character (len=16) :: str2
            logical :: have_next_token
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
            cnt = 0
            have_next_token = .false.
         reaction_loop: do
               if (.not. have_next_token) t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  call error; return
               end if
               id = rates_reaction_id(string)
               if (id == 0) then
                  if (add_flag) then
                     call add_this_reaction(string, ir, ierr)
                     if (ierr /= 0) then
                        ierr = 0
                        write(*,*) 'add_reaction failed: unknown reaction name ' // trim(string)
                        cnt = cnt+1
                     end if
                  else
                     write(*,*) 'remove: unknown reaction name ' // trim(string)
                     cnt = cnt+1
                  end if
               else
                  if (dbg) write(*,'(a30,i5,a30)') 'reaction ' // trim(string), &
                     id, trim(reaction_Name(id))
                  if (add_flag) then
                     call add_net_reaction(handle, id, ierr)
                  else
                     call remove_net_reaction(handle, id, ierr)
                  end if
                  if (ierr /= 0) then
                     call error; return
                  end if
               end if
               t = token(iounit, n, i, buffer, string)
               if (t == right_paren_token) exit reaction_loop
               if (t /= comma_token) then
                  have_next_token = .true.
               else
                  have_next_token = .false.
               end if
            end do reaction_loop
            if (cnt > 0) ierr = -1
         end subroutine do_reactions
         
         
         subroutine do_basic_reactions_for_isos(ierr)
            use chem_lib, only: chem_get_element_id
            use chem_lib, only: lookup_ZN
            integer, intent(out) :: ierr
            integer :: A, A1, A2, Z
            logical :: have_next_token, dbg
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
            have_next_token = .false.
            dbg = .false.
         iso_loop: do
               if (.not. have_next_token) t = token(iounit, n, i, buffer, string)
               !write(*,*) 'iso_loop 1st token ' // trim(string)
               !dbg = string == 'al26-1'
               if (t /= name_token) then
                  !write(*,*) 'e 1'
                  call error; return
               end if
               id = chem_get_iso_id(string)
               if (dbg) write(*,*) trim(string) // ' id', id
               if (id > 0) then
                  Z = chem_isos% Z(id)
                  A1 = chem_isos% N(id) + Z
                  A2 = A1
                  A = A1
                  call do1_iso(id,Z,A,A1,A2,ierr)
                  if (ierr /= 0) then
                     call error; return
                  end if
               else
                  Z = chem_get_element_id(string)
                  if (Z < 0) then
                     read(string,fmt=*,iostat=ierr) Z
                     if (Z < 0) then
                        write(*,*) 'unknown name ' // trim(string)
                        call error; return
                     end if
                  end if
                  t = token(iounit, n, i, buffer, string)
                  !write(*,*) 'iso_loop 2nd token ' // trim(string)
                  if (t /= name_token) then
                     !write(*,*) 'e 2'
                     call error; return
                  end if
                  read(string,fmt=*,iostat=ierr) A1
                  if (ierr /= 0) then
                     !write(*,*) 'string <' // trim(string) // '>'
                     call error; return
                  end if
                  t = token(iounit, n, i, buffer, string)
                  !write(*,*) 'iso_loop 3rd token ' // trim(string)
                  if (t /= name_token) then
                     !write(*,*) 'e 3'
                     call error; return
                  end if
                  read(string,fmt=*,iostat=ierr) A2
                  if (ierr /= 0) then
                     !write(*,*) 'string <' // trim(string) // '>'
                     call error; return
                  end if
                  do A = A1, A2
                     id = lookup_ZN(Z, A-Z)
                     call do1_iso(id,Z,A,A1,A2,ierr)
                     if (ierr /= 0) then
                        call error; return
                     end if
                  end do
               end if
               t = token(iounit, n, i, buffer, string)
               !write(*,*) 'iso_loop final token ' // trim(string)
               if (t == right_paren_token) exit iso_loop
               if (t /= comma_token) then
                  have_next_token = .true.
               else
                  have_next_token = .false.
               end if
            end do iso_loop
            
         end subroutine do_basic_reactions_for_isos
         
         
         subroutine do1_iso(id,Z,A,A1,A2,ierr)
            integer, intent(in) :: id,Z,A,A1,A2
            integer, intent(out) :: ierr
            ierr = 0
            if (id <= 0) then
               write(*,'(a30,5i5,a30)') 'bad iso for ' // trim(string), Z, A-Z, A, A1, A2
               ierr = -1; call error; return
               stop
            end if
            if (g% net_iso(id) == 0) then
               if (dbg) write(*,'(a30,i5,a30)') 'add_net_iso', id, trim(chem_isos% name(id))
               call add_net_iso(handle, id, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in add_net_iso for ' // trim(chem_isos% name(id))
                  call error; return
               end if
            end if            
            call do_basic_reactions_for_iso(id, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do_basic_reactions_for_iso for ' // &
                  trim(chem_isos% name(id))
               call error; return
            end if
         end subroutine do1_iso
         
         
         subroutine do_basic_reactions_for_iso(id, ierr)
            ! pg X, X gp, ag X, X ag, ap X, X pa, ng X, X gn,
            ! X wk, X wk_minus, X wk_h1, X wk_he4
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            integer :: Z, N, Z2, N2, other_id
            integer :: i, j, k, ir, r_ir
            logical :: matches
            
            include 'formats'

            ierr = 0
            special_loop: do j = 1, num_special_case_reactions
               matches = .false.
               do i = 1, max_num_special_case_reactants
                  k = special_case_reactants(i,j)
                  if (k == 0) exit
                  if (k == id) then
                     matches = .true.; exit
                  end if
               end do
               if (.not. matches) cycle special_loop
               ! check that all are present
               do i = 1, max_num_special_case_reactants
                  k = special_case_reactants(i,j)
                  if (k /= 0) then
                     if (g% net_iso(k) == 0) cycle special_loop
                  end if
               end do

               if (len_trim(special_case_reactions(1,j)) > 0) then
                  call add_this_reaction(special_case_reactions(1,j), ir, ierr)
                  if (ierr /= 0) return
               else
                  ir = 0
               end if

               if (len_trim(special_case_reactions(2,j)) > 0) then
                  call add_this_reaction(special_case_reactions(2,j), r_ir, ierr)
                  if (ierr /= 0) return
               else
                  r_ir = 0
               end if

               if (ir > 0) then
                  reverse_reaction_id(ir) = r_ir
               else if (r_ir > 0) then
                  !write(*,*) 'failed to find reverse for special ' // trim(reaction_name(r_ir))
               end if
               if (r_ir > 0) then
                  reverse_reaction_id(r_ir) = ir 
               else if (ir > 0) then
                  !write(*,*) 'failed to find reverse for special ' // trim(reaction_name(ir))
               end if

            end do special_loop
            
            Z = chem_isos% Z(id)
            N = chem_isos% N(id)

            if (g% net_iso(ih1) > 0) then              
               other_id = get_iso_id(Z-1, N)
               call add_basic_reaction(other_id, id, 'pg', 'gp', .true., ir, ierr)
               if (ierr /= 0) return            
               call add_basic_reaction(id, other_id, 'gp', 'pg', .true., r_ir, ierr)
               if (ierr /= 0) return
               if (ir > 0) reverse_reaction_id(ir) = r_ir
               if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
               other_id = get_iso_id(Z+1, N)
               call add_basic_reaction(id, other_id, 'pg', 'gp', .true., ir, ierr)
               if (ierr /= 0) return            
               call add_basic_reaction(other_id, id, 'gp', 'pg', .true., r_ir, ierr)
               if (ierr /= 0) return        
               if (ir > 0) reverse_reaction_id(ir) = r_ir
               if (r_ir > 0) reverse_reaction_id(r_ir) = ir                   
            end if
            
            if (g% net_iso(ihe4) > 0) then
               if (id == ib8) then
                  call add_this_reaction('r_b8_wk_he4_he4', ir, ierr)
                  if (ierr /= 0) return
               end if
               other_id = get_iso_id(Z-2, N-2)
               call add_basic_reaction(other_id, id, 'ag', 'ga', .true., ir, ierr)
               if (ierr /= 0) return
               call add_basic_reaction(id, other_id, 'ga', 'ag', .true., r_ir, ierr)
               if (ierr /= 0) return            
               if (ir > 0) reverse_reaction_id(ir) = r_ir
               if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
               other_id = get_iso_id(Z+2, N+2)
               call add_basic_reaction(id, other_id, 'ag', 'ga', .true., ir, ierr)
               if (ierr /= 0) return            
               call add_basic_reaction(other_id, id, 'ga', 'ag', .true., r_ir, ierr)
               if (ierr /= 0) return            
               if (ir > 0) reverse_reaction_id(ir) = r_ir
               if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
            end if
            
            if (g% net_iso(ih1) > 0 .and. g% net_iso(ihe4) > 0) then
               if (id /= ihe4) then
                  other_id = get_iso_id(Z-1, N-2)
                  call add_basic_reaction(other_id, id, 'ap', 'pa', .true., ir, ierr)
                  if (ierr /= 0) return            
                  call add_basic_reaction(id, other_id, 'pa', 'ap', .true., r_ir, ierr)
                  if (ierr /= 0) return
                  if (ir > 0) reverse_reaction_id(ir) = r_ir
                  if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
                  other_id = get_iso_id(Z+1, N+2)
                  call add_basic_reaction(id, other_id, 'ap', 'pa', .true., ir, ierr)
                  if (ierr /= 0) return            
                  call add_basic_reaction(other_id, id, 'pa', 'ap', .true., r_ir, ierr)
                  if (ierr /= 0) return
                  if (ir > 0) reverse_reaction_id(ir) = r_ir
                  if (r_ir > 0) reverse_reaction_id(r_ir) = ir
               end if          
            end if
            
            if (id /= ih1 .and. g% net_iso(ih1) > 0 .and. g% net_iso(ineut) > 0) then
               other_id = get_iso_id(Z-1, N+1)
               call add_basic_reaction(id, other_id, 'np', 'pn', .true., ir, ierr)
               if (ierr /= 0) return
               call add_basic_reaction(other_id, id, 'pn', 'np', .true., r_ir, ierr)
               if (ierr /= 0) return          
               if (ir > 0) reverse_reaction_id(ir) = r_ir
               if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
               if (id /= ineut) then
                  other_id = get_iso_id(Z+1, N-1)
                  call add_basic_reaction(other_id, id, 'np', 'pn', .true., ir, ierr)
                  if (ierr /= 0) return
                  call add_basic_reaction(id, other_id, 'pn', 'np', .true., r_ir, ierr)
                  if (ierr /= 0) return 
                  if (ir > 0) reverse_reaction_id(ir) = r_ir
                  if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
               end if
            end if
            
            if (id /= ihe4 .and. g% net_iso(ihe4) > 0 .and. g% net_iso(ineut) > 0) then
               other_id = get_iso_id(Z-2, N-1)
               if (other_id /= ihe4) then
                  call add_basic_reaction(id, other_id, 'na', 'an', .true., ir, ierr)
                  if (ierr /= 0) return
                  call add_basic_reaction(other_id, id, 'an', 'na', .true., r_ir, ierr)
                  if (ierr /= 0) return
                  if (ir > 0) reverse_reaction_id(ir) = r_ir
                  if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
               end if   
               other_id = get_iso_id(Z+2, N+1)
               if (other_id /= ihe4) then
                  call add_basic_reaction(other_id, id, 'na', 'an', .true., ir, ierr)
                  if (ierr /= 0) return
                  call add_basic_reaction(id, other_id, 'an', 'na', .true., r_ir, ierr)
                  if (ierr /= 0) return
                  if (ir > 0) reverse_reaction_id(ir) = r_ir
                  if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
               end if   
            end if
            
            if (g% net_iso(ineut) > 0 .and. id /= ih2) then            
               other_id = get_iso_id(Z, N-1)
               call add_basic_reaction(other_id, id, 'ng', 'gn', .true., ir, ierr)
               if (ierr /= 0) return            
               call add_basic_reaction(id, other_id, 'gn', 'ng', .true., r_ir, ierr)
               if (ierr /= 0) return
                  if (ir > 0) reverse_reaction_id(ir) = r_ir
                  if (r_ir > 0) reverse_reaction_id(r_ir) = ir 
            end if
            
            other_id = get_iso_id(Z-1, N+1)
            call add_basic_reaction(id, other_id, 'wk', '', .true., ir, ierr)
            if (ierr /= 0) return
            call add_basic_reaction(other_id, id, 'wk-minus', '', .false., r_ir, ierr)
            if (ierr /= 0) then
               write(*,3) 'failed in add_basic_reaction wk-minus', other_id, id
               return
            end if
            if (ir > 0) reverse_reaction_id(ir) = r_ir
            if (r_ir > 0) reverse_reaction_id(r_ir) = ir
            
            other_id = get_iso_id(Z+1, N-1)
            call add_basic_reaction(other_id, id, 'wk', '', .true., ir, ierr)
            if (ierr /= 0) return
            call add_basic_reaction(id, other_id, 'wk-minus', '', .false., r_ir, ierr)
            if (ierr /= 0) then
               write(*,3) 'failed in add_basic_reaction wk-minus', other_id, id
               return
            end if
            if (ir > 0) reverse_reaction_id(ir) = r_ir
            if (r_ir > 0) reverse_reaction_id(r_ir) = ir
            
            if (g% net_iso(ih1) > 0) then  
               other_id = get_iso_id(Z-2, N+1)
               call add_basic_reaction(id, other_id, 'wk-h1', '', .false., ir, ierr)
               if (ierr /= 0) return
               other_id = get_iso_id(Z+2, N-1)
               call add_basic_reaction(other_id, id, 'wk-h1', '', .false., ir, ierr)
               if (ierr /= 0) return
            end if
            
            if (g% net_iso(ihe4) > 0) then  
               other_id = get_iso_id(Z-3, N-1)
               call add_basic_reaction(id, other_id, 'wk-he4', '', .false., ir, ierr)
               if (ierr /= 0) return
               if (ir > 0) reverse_reaction_id(ir) = 0
               other_id = get_iso_id(Z+3, N+1)
               call add_basic_reaction(other_id, id, 'wk-he4', '', .false., ir, ierr)
               if (ierr /= 0) return
               if (ir > 0) reverse_reaction_id(ir) = 0
            end if


! fxt for al26 isomers
            if (trim(chem_isos% name(id)) == 'al26-1' .and. g% net_iso(id) > 0) then
             call add_this_reaction('r_al26-1_to_al26-2', ir, ierr)
            end if
            if (trim(chem_isos% name(id)) == 'al26-2' .and. g% net_iso(id) > 0) then
             call add_this_reaction('r_al26-2_to_al26-1', ir, ierr)
            end if

         end subroutine do_basic_reactions_for_iso
         
         
         integer function get_iso_id(Z, N)
            use chem_lib, only: lookup_ZN
            integer, intent(in) :: Z, N
            get_iso_id = lookup_ZN(Z, N)
            if (get_iso_id <= 0) return
            if (g% net_iso(get_iso_id) <= 0) get_iso_id = 0
         end function get_iso_id
         
         
         subroutine add_this_reaction(string, ir, ierr)
            use rates_lib, only: rates_reaction_id
            character (len=*), intent(in) :: string
            integer, intent(out) :: ir, ierr
            logical :: okay, dbg
            dbg = .false. ! (string == 'r_h1_li7_to_he4_he4')
            okay = add_this_reaclib_forward_reaction(string, ir, ierr)
            if (ierr /= 0) then
               write(*,*) 'ierr from add_this_reaclib_forward_reaction ' // trim(string)
               !stop
               return
            end if
            if (okay) then
               if (dbg) write(*,*) 'add_this_reaclib_forward_reaction ' // &
                  trim(string), rates_reaction_id(string)
               return
            end if
            okay = add_this_reaclib_reverse_reaction(string, ir, ierr)
            if (ierr /= 0) then
               write(*,*) 'ierr from add_this_reaclib_reverse_reaction ' // trim(string)
               !stop
               return
            end if
            if (okay) then
               if (dbg) write(*,*) 'add_this_reaclib_reverse_reaction ' // trim(string)
               return
            end if
            okay = add_reaction_for_this_handle(string, ir, ierr)
            if (ierr /= 0) then
               write(*,*) 'ierr from add_reaction_for_this_handle ' // trim(string)
               !stop
               return
            end if
            if (.not. okay) then
               write(*,*) 'WARNING: failed to add reaction ' // trim(string)
               stop
            end if
            if (dbg) write(*,*) 'add_reaction_for_this_handle ' // trim(string)
            if (dbg) stop 'add_this_reaction'
         end subroutine add_this_reaction
         

         logical function add_reaction_for_this_handle(string, ir, ierr)
            use rates_lib, only: rates_reaction_id, add_reaction_for_handle
            character (len=*), intent(in) :: string
            integer, intent(out) :: ir, ierr
            include 'formats'
            ierr = 0
            add_reaction_for_this_handle = .false.
            ir = rates_reaction_id(string)
            if (ir == 0) then ! check if reaction is defined in reaclib
               call add_reaction_for_handle(string, ierr)
               if (ierr /= 0) then  
                  write(*,*) 'failed in add_reaction_for_handle for ' // trim(string)
                  return
               end if
               ir = rates_reaction_id(string)
            end if            
            if (ir > 0) then
               call add_net_reaction(handle, ir, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in add_reaction_for_this_handle ' // trim(string)
                  call error
               end if
               add_reaction_for_this_handle = .true.
            end if            
         end function add_reaction_for_this_handle
         

         logical function add_this_reaclib_forward_reaction(string, ir, ierr)
            use rates_lib, only: rates_reaction_id, add_reaction_from_reaclib, reaclib_lookup
            character (len=*), intent(in) :: string
            integer, intent(out) :: ir, ierr
            integer :: indx
            include 'formats'
            ierr = 0
            add_this_reaclib_forward_reaction = .false.
            ir = rates_reaction_id(string)
            if (ir == 0) then ! check if reaction is defined in reaclib
               indx = reaclib_lookup(string, reaclib_rates% reaction_dict)
               if (indx > 0) then ! add a definition for it
                  call add_reaction_from_reaclib(string, '', indx, ierr)
                  if (ierr /= 0) then  
                     write(*,*) 'failed in add_this_reaclib_forward_reaction for ' // trim(string)
                     return
                  end if
                  ir = rates_reaction_id(string)
               end if
            end if            
            if (ir > 0) then
               call add_net_reaction(handle, ir, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in add_this_reaclib_forward_reaction ' // trim(string)
                  call error
               end if
               add_this_reaclib_forward_reaction = .true.
            end if            
         end function add_this_reaclib_forward_reaction
         
         
         logical function add_this_reaclib_reverse_reaction(string, ir, ierr)
            use rates_lib, only: rates_reaction_id, add_reaction_from_reaclib, reaclib_lookup
            character (len=*), intent(in) :: string
            integer, intent(out) :: ir, ierr
            integer :: indx
            character (len=256) :: forward
            ierr = 0
            add_this_reaclib_reverse_reaction = .false.
            ir = rates_reaction_id(string)
            if (ir == 0) then ! check if reaction is defined in reaclib
               indx = reaclib_lookup(string, reaclib_rates% reverse_dict)              
               if (indx > 0) then ! add a definition for it
                  call add_reaction_from_reaclib( &
                           string, reaclib_rates% reaction_handle(indx), indx, ierr)
                  if (ierr /= 0) return
                  ir = rates_reaction_id(string)
               end if
            end if            
            if (ir > 0) then
               call add_net_reaction(handle, ir, ierr)
               if (ierr /= 0) then
                  call error
               end if
               add_this_reaclib_reverse_reaction = .true.
               return
            end if            
         end function add_this_reaclib_reverse_reaction
         
         
         subroutine add_basic_reaction(iso_in, iso_out, op, reverse_op, warn, ir, ierr)
            use rates_lib, only: rates_reaction_id, add_reaction_from_reaclib, reaclib_lookup
            use rates_def, only: reaction_inputs
            integer, intent(in) :: iso_in, iso_out
            character (len=*), intent(in) :: op, reverse_op
            logical, intent(in) :: warn
            integer, intent(out) :: ir, ierr
            character (len=100) :: string, reverse
            integer :: indx, i
            logical, parameter :: dbg = .false.
            include 'formats'
            ierr = 0
            ir = 0
            
            if (iso_in <= 0 .or. iso_out <= 0) return
            
            string = 'r_' // trim(chem_isos% name(iso_in)) // '_' // &
                  trim(op) // '_' // trim(chem_isos% name(iso_out))
            
            ir = rates_reaction_id(string)
                  
            if (ir < 0 .or. ir > rates_reaction_id_max) then
               write(*,*) 'failed in rates_reaction_id for ' // trim(string)
               stop 'net_init'
            end if
                  
            if (ir == 0) then ! check if reaction or reverse is defined in reaclib
               indx = reaclib_lookup(string, reaclib_rates% reaction_dict)
               if (indx > 0) then ! add a definition for it
                  call add_reaction_from_reaclib(string, '', indx, ierr)
                  if (ierr /= 0) then  
                     write(*,*) 'failed in add_reaction_from_reaclib for ' // trim(string)
                     return
                  end if
                  ir = rates_reaction_id(string)
                  
                  if (ir < 0 .or. ir > rates_reaction_id_max) then
                     write(*,*) 'failed in rates_reaction_id for ' // trim(string)
                     stop 'net_init'
                  end if
                  
               else
                  ierr = 0
                  if (len_trim(reverse_op) > 0) then ! check for the reverse reaction in reaclib
                     reverse = 'r_' // trim(chem_isos% name(iso_out)) // '_' // &
                        trim(reverse_op) // '_' // trim(chem_isos% name(iso_in))
                     indx = reaclib_lookup(reverse, reaclib_rates% reaction_dict)
                     if (indx > 0) then ! add a definition for it
                        call add_reaction_from_reaclib(string, reverse, indx, ierr)
                        if (ierr /= 0) then  
                           write(*,'(a)') 'failed in add_reaction_from_reaclib for '  &
                              // trim(string) // ' reverse ' // trim(reverse)
                           return
                        end if
                        ir = rates_reaction_id(string)
                     else
                        ierr = 0
                     end if
                  end if
               end if
            end if
            
            if (ir > 0) then
               call add_net_reaction(handle, ir, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in add_net_reaction ' // trim(string)
                  call error
               end if
               if (dbg) write(*,*) 'add_basic_reaction ' // trim(string)
               return
            end if
            
            if (warn) then
               call integer_dict_lookup(skip_warnings_dict, string, i, ierr)
               if (ierr /= 0 .or. i <= 0) then
                  ierr = 0
               end if
            end if
            
         end subroutine add_basic_reaction
         

         
      end subroutine do_read_net_file

            
      subroutine start_net_def(handle, ierr)
         use rates_def, only: rates_reaction_id_max
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer  :: g
         
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (.not. associated(g% net_iso)) then
            allocate(g% net_iso(num_chem_isos), stat=ierr)
            if (ierr /= 0) return
         end if
         g% net_iso(:) = 0
         
         if (.not. associated(g% net_reaction)) then
            allocate(g% net_reaction(rates_reaction_id_max), stat=ierr)
            if (ierr /= 0) return
         end if
         g% net_reaction(:) = 0    
         
         g% net_has_been_defined = .false. 
         g% net_filename = ''  
         
         
         !call show_scr3('start_net_def')
         
      end subroutine start_net_def
      
      
      subroutine show_scr3(str)
         use rates_def
         use chem_def
         character (len=*), intent(in) :: str
         integer :: ir
         include 'formats'
         write(*,*) trim(str)
         do ir = 1, rates_reaction_id_max
            if (reaction_screening_info(3,ir) <= 0) cycle
            write(*,2) 'scr 3 ' // trim(reaction_Name(ir))  &
                  // ' ' // trim(chem_isos% name(reaction_screening_info(1,ir)))  &
                  // ' ' // trim(chem_isos% name(reaction_screening_info(2,ir)))  &
                  // ' ' // trim(chem_isos% name(reaction_screening_info(3,ir)))
         end do
         write(*,*)
      end subroutine show_scr3
      
      
      subroutine finish_net_def(handle, ierr)
         use rates_def, only: reaction_names_dict
         use utils_lib, only: integer_dict_create_hash, realloc_integer
         integer, intent(in) :: handle
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer  :: g
         integer :: i
         
         logical, parameter :: dbg = .false.
         !logical, parameter :: dbg = .true.
         
         include 'formats'
         
         ierr = 0
         
         if (dbg) write(*,*) 'finish_net_def'
         
         ! may have defined some new reactions as part of loading the net
         ! so recreate the dict hash 
         call integer_dict_create_hash(reaction_names_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: finish_net_def failed in integer_dict_create_hash'
            return
         end if

         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (.not. associated(g% net_iso)) then
            ierr = -1
            return
         end if
         
         if (.not. associated(g% net_reaction)) then
            ierr = -1
            return
         end if
         
         if (size(g% net_reaction, dim=1) > rates_reaction_id_max) then
            if (dbg) write(*,*) 'call realloc_integer'
            call realloc_integer(g% net_reaction, rates_reaction_id_max, ierr)
            if (ierr /= 0) then
               write(*,*) 'finish_net_def failed in realloc_integer'
               return
            end if
         end if
         
         if (g% doing_approx21) then
            if (dbg) write(*,*) 'call mark_approx21'
            call do_mark_approx21(handle, ierr)
            if (ierr /= 0) return
         end if
         
         if (dbg) write(*,*) 'call setup_iso_info'
         call setup_iso_info(g, ierr)
         if (ierr /= 0) return

         if (dbg) write(*,*) 'call setup_reaction_info'
         call setup_reaction_info(g, ierr)
         if (ierr /= 0) return         
         
         allocate( &
            g% which_rates(rates_reaction_id_max),  &
            g% reaction_max_Z(g% num_reactions),  &
            g% reaction_max_Z_plus_N_for_max_Z(g% num_reactions),  &
            stat=ierr)
         if (ierr /= 0) return
         
         g% which_rates(:) = 1
         
         allocate( &
            g% z158(g% num_isos),  &
            g% z52(g% num_isos),  &
            stat=ierr)
         if (ierr /= 0) return
         ! Precompute some powers of Z for screening and coulomb
         do i=1, g% num_isos
            g% z158(i) = pow(real(chem_isos% Z(g% chem_id(i)),kind=dp),1.58d0)
            g% z52(i)  = pow(real(chem_isos% Z(g% chem_id(i)),kind=dp),2.50d0) ! 5.d0/2.d0)
         end do
         
         call set_reaction_max_Z(g)

         g% net_has_been_defined = .true.   
         
         if (g% doing_approx21) then
            if (dbg) write(*,*) 'call set_approx21'
            call do_set_approx21(handle, ierr)
            if (ierr /= 0) return
         end if

         if (dbg) write(*,*) 'done finish_net_def'
         
         contains
         
         subroutine do_mark_approx21(handle, ierr)
            use net_approx21, only: mark_approx21
            integer, intent(in) :: handle
            integer, intent(out) :: ierr
            call mark_approx21(handle, ierr)
         end subroutine do_mark_approx21
      
         
         subroutine do_set_approx21(handle, ierr)
            use net_approx21, only: set_approx21
            integer, intent(in) :: handle
            integer, intent(out) :: ierr
            call set_approx21(handle, ierr)
         end subroutine do_set_approx21
         
      
      end subroutine finish_net_def
      
         
      subroutine check_for_hardwired_pairs ! especially for approx21
         call check_pair('r_he4_he4_he4_to_c12', 'r_c12_to_he4_he4_he4')
         call check_pair('rhe4_breakup', 'rhe4_rebuild') 
                              ! << should treat this as he4 + n + p -> 3 n + 3 p ???
         call check_pair('rprot_to_neut', 'rneut_to_prot')
         call check_pair('r_c12_ag_o16', 'r_o16_ga_c12')
         call check_pair('rc12ap_to_o16', 'ro16gp_to_c12')
         call check_pair('r_o16_ag_ne20', 'r_ne20_ga_o16')
         call check_pair('ro16ap_to_ne20', 'rne20gp_to_o16')      
         call check_pair('r_ne20_ag_mg24', 'r_mg24_ga_ne20')
         call check_pair('rne20ap_to_mg24', 'rmg24gp_to_ne20')
         call check_pair('r_mg24_ag_si28', 'r_si28_ga_mg24')
         call check_pair('rmg24ap_to_si28', 'rsi28gp_to_mg24')
         call check_pair('r_si28_ag_s32', 'r_s32_ga_si28')
         call check_pair('rsi28ap_to_s32', 'rs32gp_to_si28')
         call check_pair('r_s32_ag_ar36', 'r_ar36_ga_s32')
         call check_pair('rs32ap_to_ar36', 'rar36gp_to_s32')
         call check_pair('r_ar36_ag_ca40', 'r_ca40_ga_ar36')
         call check_pair('rar36ap_to_ca40', 'rca40gp_to_ar36')
         call check_pair('r_ca40_ag_ti44', 'r_ti44_ga_ca40')
         call check_pair('rca40ap_to_ti44', 'rti44gp_to_ca40')
         call check_pair('r_ti44_ag_cr48', 'r_cr48_ga_ti44')
         call check_pair('rti44ap_to_cr48', 'rcr48gp_to_ti44')
         call check_pair('r_cr48_ag_fe52', 'r_fe52_ga_cr48')
         call check_pair('rcr48ap_to_fe52', 'rfe52gp_to_cr48')              
         call check_pair('rfe52aprot_to_fe54', 'rfe54prot_to_fe52')     
         call check_pair('rfe52neut_to_fe54', 'rfe54g_to_fe52')     
         call check_pair('rfe54ng_to_fe56', 'rfe56gn_to_fe54')     
         call check_pair('r_fe52_ag_ni56', 'r_ni56_ga_fe52')     
         call check_pair('rfe52aprot_to_ni56', 'rni56gprot_to_fe52')     
         call check_pair('rfe54prot_to_ni56', 'rni56gprot_to_fe54')     
         call check_pair('rfe54prot_to_ni56', 'rni56gprot_to_fe54')     
      end subroutine check_for_hardwired_pairs
         
         
      subroutine check_pair(s1,s2)
         use rates_def, only: get_rates_reaction_id, reverse_reaction_id
         character (len=*), intent(in) :: s1, s2
         integer :: ir, r_ir
         ir = get_rates_reaction_id(s1)
         if (ir <= 0) then
            write(*,*) 'missing ' // trim(reaction_name(ir))
            return
         end if
         r_ir = get_rates_reaction_id(s2)
         if (r_ir <= 0) then
            write(*,*) 'missing ' // trim(reaction_name(r_ir))
            return
         end if
         reverse_reaction_id(ir) = r_ir
         reverse_reaction_id(r_ir) = ir
      end subroutine check_pair
         
      
         
      subroutine set_reaction_kinds(g)
         use chem_def, only: chem_isos
         use rates_lib, only: &
            reaclib_create_handle, reaclib_parse_handle, is_weak_reaction
         type (Net_General_Info), pointer  :: g
         
         integer :: iso_in1, iso_in2, iso_out1, iso_out2, &
            num_in1, num_in2, num_out1, num_out2, in_a, in_b, out_c, out_d, &
            num_basic_ng_kind, num_basic_pn_kind, num_basic_pg_kind, &
            num_basic_ap_kind, num_basic_an_kind, num_basic_ag_kind, &
            num_general_one_one_kind, &
            num_general_two_one_kind, num_general_two_two_kind
         integer :: i, ir, num_in, num_out, iso_ids(10), r_ir, r_i, ierr
         character (len=100) :: reverse_handle, op
         character (len=4) :: rstr
         integer, pointer :: kind(:), reverse_id(:)
         
         include 'formats'
         
         kind => g% reaction_reaclib_kind
         reverse_id => g% reverse_id_for_kind_ne_other
         
         num_basic_ng_kind = 0
         num_basic_pn_kind = 0
         num_basic_pg_kind = 0
         num_basic_ap_kind = 0
         num_basic_an_kind = 0
         num_basic_ag_kind = 0
         num_general_one_one_kind = 0
         num_general_two_one_kind = 0
         num_general_two_two_kind = 0
         
         do i=1, g% num_reactions
         
            kind(i) = other_kind
            reverse_id(i) = 0
            ir = g% reaction_id(i)         
            
            if (is_weak_reaction(ir)) cycle ! don't do weak reactions
            if (reaction_name(ir)(1:2) /= 'r_') cycle ! don't mess with special ones.

            r_ir = reverse_reaction_id(ir)
            if (r_ir <= 0) then
               cycle ! no reverse reaction
            end if
            reverse_id(i) = r_ir
            r_i = g% net_reaction(r_ir)
            if (r_i <= 0) then
               cycle ! reverse reaction not in this net 
            end if
            
            if (reaction_inputs(6,ir) /= 0) then
               cycle ! more than 2 input species
            end if
            if (reaction_outputs(6,ir) /= 0) then
               cycle ! more than 2 output species
            end if

            num_in1 = reaction_inputs(1,ir)
            iso_in1 = reaction_inputs(2,ir)
            
            num_in2 = reaction_inputs(3,ir)
            iso_in2 = reaction_inputs(4,ir)
            
            num_out1 = reaction_outputs(1,ir)
            iso_out1 = reaction_outputs(2,ir)
            
            num_out2 = reaction_outputs(3,ir)
            iso_out2 = reaction_outputs(4,ir)

            if (iso_in1 == 0 .or. iso_out1 == 0) then
               cycle ! non-standard reaction
            end if
            if (num_in1 > 3 .or. num_out1 > 3) then
               cycle ! non-standard reaction
            end if
            
            if (iso_in2 == 0) then ! only 1 species on lhs
               if (iso_out2 > 0) cycle ! 1 to many is treated as reverse of many to 1
               ! 1 species to 1 species reaction
               if (num_in1 == 1 .and. num_out1 == 1) cycle
               kind(i) = general_one_one_kind
               if (r_ir < ir) then
                  kind(i) = other_kind
               else
                  num_general_one_one_kind = num_general_one_one_kind + 1
               end if
               cycle
            end if
            
            if (is_weak_reaction(ir)) cycle
            
            if (iso_out2 == 0) then ! 2 to 1
               if (num_in1 > 1 .or. num_in2 > 1 .or. num_out1 > 1) then
                  kind(i) = general_two_one_kind
                  if (r_ir <= 0) then
                     kind(i) = other_kind
                  else
                     num_general_two_one_kind = num_general_two_one_kind + 1
                  end if
                  cycle
               end if
            else ! 2 to 2
               if (num_in1 > 1 .or. num_in2 > 1 .or. num_out1 > 1 .or. num_out2 > 1) then
                  kind(i) = general_two_two_kind
                  if (r_ir <= 0) then
                     kind(i) = other_kind
                  else
                     num_general_two_two_kind = num_general_two_two_kind + 1
                  end if
                  cycle
               end if
            end if

            ! if get here, have a + b on left with different species a and b
            ! and only have 1 of each species on left and right
            
            if (chem_isos% Z(iso_in1) <= chem_isos% Z(iso_in2)) then
               in_a = iso_in1
               in_b = iso_in2
            else
               in_a = iso_in2
               in_b = iso_in1
            end if
            ! a + b => c; Z(a) <= Z(b)

            if (iso_out2 == 0) then ! 2-to-1
               
               if (in_a == ineut) then
                  kind(i) = ng_kind
                  num_basic_ng_kind = num_basic_ng_kind + 1
               else if (in_a == ihe4) then
                  kind(i) = ag_kind
                  num_basic_ag_kind = num_basic_ag_kind + 1                     
               else if (in_a == ih1) then
                  kind(i) = pg_kind
                  num_basic_pg_kind = num_basic_pg_kind + 1
               else
                  kind(i) = general_two_one_kind
                  if (r_ir <= 0) then
                     kind(i) = other_kind
                     cycle
                  else
                     num_general_two_one_kind = num_general_two_one_kind + 1
                     reverse_id(i) = r_ir
                     !write(*,*) 'general 2-1 ' // trim(reaction_name(ir))
                  end if
               end if
               
            else ! 2-to-2. only do ap,an,pn.  skip pa,na,np.
            
               if (chem_isos% Z(iso_out1) <= chem_isos% Z(iso_out2)) then
                  out_c = iso_out1
                  out_d = iso_out2
               else
                  out_c = iso_out2
                  out_d = iso_out1
               end if
               ! a + b => c + d; Z(a) <= Z(b) and Z(c) <= Z(d)
               if (in_a == ihe4) then
                  if (out_c == ih1) then
                     kind(i) = ap_kind
                     num_basic_ap_kind = num_basic_ap_kind + 1
                  else if (out_c == ineut) then
                     kind(i) = an_kind
                     num_basic_an_kind = num_basic_an_kind + 1
                  else
                     kind(i) = general_two_two_kind
                  end if
               else if (in_a == ih1) then
                  if (out_c == ineut) then
                     kind(i) = pn_kind
                     num_basic_pn_kind = num_basic_pn_kind + 1
                  else if (out_c == ihe4) then
                     kind(i) = other_kind
                     cycle
                  else
                     kind(i) = general_two_two_kind
                  end if
               else if (in_a == ineut) then
                  if (out_c == ih1 .or. out_c == ihe4) then
                     kind(i) = other_kind
                     cycle
                  else
                     kind(i) = general_two_two_kind
                  end if
               else
                  kind(i) = general_two_two_kind
               end if
                              
            end if
            
            if (kind(i) == general_two_two_kind) then
               if (r_ir > ir) then
                  kind(i) = other_kind
                  cycle
               else
                  num_general_two_two_kind = num_general_two_two_kind + 1
               end if
            end if
                     
         end do
         
         return
         
         if (.true.) then
            do i=1, g% num_reactions
               ir = g% reaction_id(i)  
               write(*,3) trim(reaction_name(ir)), i, kind(i)
               !r_ir = reverse_reaction_id(ir)
               !if (kind(i) == ag_kind) write(*,'(a60,i5)') trim(reaction_name(ir)) // &
               !   ' ' // trim(reaction_name(r_ir)), kind(i)
               !if (kind(i) /= other_kind) write(*,'(a60,i5)') trim(reaction_name(ir)) // &
               !   ' ' // trim(reaction_name(r_ir)), kind(i)
            end do
         end if
         
         stop
         
         !return
         
         write(*,2) 'num_basic_ag_kind', num_basic_ag_kind
         write(*,2) 'num_basic_pg_kind', num_basic_pg_kind
         write(*,2) 'num_basic_ng_kind', num_basic_ng_kind
         write(*,2) 'num_basic_pn_kind', num_basic_pn_kind
         write(*,2) 'num_basic_ap_kind', num_basic_ap_kind
         write(*,2) 'num_basic_an_kind', num_basic_an_kind
         write(*,2) 'num_general_one_one_kind', num_general_one_one_kind
         write(*,2) 'num_general_two_one_kind', num_general_two_one_kind
         write(*,2) 'num_general_two_two_kind', num_general_two_two_kind
         write(*,*)
         
         !stop

      end subroutine set_reaction_kinds      
      
      
      subroutine add_net_iso(handle, iso_id, ierr)
         use chem_def, only: chem_isos
         integer, intent(in) :: handle
         integer, intent(in) :: iso_id
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer  :: g
         
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (.not. associated(g% net_iso)) then
            ierr = -1
            return
         end if
         
         if (g% net_has_been_defined) then
            ierr = -1
            return
         end if
         
         if (iso_id < 0 .or. iso_id > num_chem_isos) then
            ierr = -1
            return
         end if
         
         g% net_iso(iso_id) = 1 ! mark as added

      end subroutine add_net_iso
      
      
      subroutine add_net_isos(handle, num_isos, iso_ids, ierr)
         integer, intent(in) :: handle
         integer, intent(in) :: num_isos, iso_ids(num_isos)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do i = 1, num_isos
            call add_net_iso(handle, iso_ids(i), ierr)
            if (ierr /= 0) return
         end do
      end subroutine add_net_isos
      
      
      subroutine remove_net_iso(handle, iso_id, ierr)
         integer, intent(in) :: handle
         integer, intent(in) :: iso_id
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer  :: g
         
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (.not. associated(g% net_iso)) then
            ierr = -1
            return
         end if
         
         if (g% net_has_been_defined) then
            ierr = -1
            return
         end if
         
         if (iso_id < 0 .or. iso_id > num_chem_isos) then
            ierr = -1
            return
         end if
         
         g% net_iso(iso_id) = 0 ! mark as removed
         
         call remove_reactions_for_iso(g, iso_id, ierr)
         
      end subroutine remove_net_iso
      
      
      subroutine remove_net_isos(handle, num_isos, iso_ids, ierr)
         integer, intent(in) :: handle
         integer, intent(in) :: num_isos, iso_ids(:)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do i = 1, num_isos
            call remove_net_iso(handle, iso_ids(i), ierr)
            if (ierr /= 0) return
         end do
      end subroutine remove_net_isos
      
      
      subroutine add_net_reaction(handle, reaction_id, ierr)
         use utils_lib, only: realloc_integer
         integer, intent(in) :: handle
         integer, intent(in) :: reaction_id
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer  :: g
         integer :: old_sz, new_sz
         
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (.not. associated(g% net_reaction)) then
            ierr = -1
            write(*,*) 'must call net_start_def before calling net_add_reaction'
            return
         end if
         
         if (g% net_has_been_defined) then
            ierr = -1
            write(*,*) 'must call net_start_def before calling net_add_reaction'
            return
         end if
         
         if (reaction_id < 0 .or. reaction_id > rates_reaction_id_max) then
            ierr = -1
            write(*,*) 'invalid reaction_id for net_add_reaction'
            return
         end if
         
         old_sz = size(g% net_reaction, dim=1)
         if (reaction_id > old_sz) then
            new_sz = (reaction_id*11)/10 + 100
            call realloc_integer(g% net_reaction,new_sz,ierr)
            if (ierr /= 0) then
               write(*,*) 'add_net_reaction failed in realloc_integer'
               return
            end if
            g% net_reaction(old_sz+1:new_sz) = 0
         end if
         
         g% net_reaction(reaction_id) = 1 ! mark as added
         
      end subroutine add_net_reaction
      
      
      subroutine add_net_reactions(handle, num_reactions, reaction_ids, ierr)
         integer, intent(in) :: handle
         integer, intent(in) :: num_reactions, reaction_ids(:)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do i = 1, num_reactions
            call add_net_reaction(handle, reaction_ids(i), ierr)
            if (ierr /= 0) return
         end do
      end subroutine add_net_reactions
      
      
      subroutine remove_net_reaction(handle, reaction_id, ierr)
         integer, intent(in) :: handle
         integer, intent(in) :: reaction_id
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer  :: g
         
         ierr = 0
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) return
         
         if (.not. associated(g% net_reaction)) then
            ierr = -1
            return
         end if
         
         if (g% net_has_been_defined) then
            ierr = -1
            return
         end if
         
         if (reaction_id <= 0 .or. reaction_id > rates_reaction_id_max) then
            ierr = -1
            return
         end if
         
         g% net_reaction(reaction_id) = 0 ! mark as removed
         
      end subroutine remove_net_reaction
      
      
      subroutine remove_net_reactions(handle, num_reactions, reaction_ids, ierr)
         integer, intent(in) :: handle
         integer, intent(in) :: num_reactions, reaction_ids(:)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do i = 1, num_reactions
            call remove_net_reaction(handle, reaction_ids(i), ierr)
            if (ierr /= 0) return
         end do
      end subroutine remove_net_reactions
      

      subroutine setup_iso_info(g, ierr)
         type (Net_General_Info), pointer  :: g
         integer, intent(out) :: ierr
         
         integer :: i, iso_num, num_isos
         integer, pointer :: itab(:), chem_id(:)
         
         ierr = 0
         itab => g% net_iso
         
         num_isos = sum(itab(:))
         if (num_isos <= 0) then
            ierr = -1
            return
         end if

         g% num_isos = num_isos
         allocate(g% chem_id(num_isos), stat=ierr)
         if (ierr /= 0) return
         chem_id => g% chem_id
         
         iso_num = 0
         do i = 1, num_chem_isos
            if (itab(i) == 0) cycle
            iso_num = iso_num + 1
            chem_id(iso_num) = i
            itab(i) = iso_num
         end do

      end subroutine setup_iso_info


      subroutine remove_reactions_for_iso(g, target_iso, ierr)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: target_iso
         integer, intent(out) :: ierr

         integer, pointer, dimension(:) :: itab, rtab
         integer :: i, j, ir, &
            num_reaction_inputs, num_reaction_outputs

         ierr = 0
         rtab => g% net_reaction
         itab => g% net_iso
         
         !if (target_iso /= 0) &
         !   write(*,*) 'remove_reactions_for_iso ' // trim(chem_isos% name(target_iso))
         
         do ir = 1, size(rtab,dim=1)
            if (rtab(ir) == 0) cycle ! not used
            num_reaction_inputs = get_num_reaction_inputs(ir)
            do i = 1, num_reaction_inputs
               j = reaction_inputs(2*i,ir)
               if (j == target_iso) then
                  rtab(ir) = 0 ! mark as removed
                  !write(*,*) 'remove reaction ' // trim(reaction_Name(ir))
                  exit
               end if
            end do
            num_reaction_outputs = get_num_reaction_outputs(ir)
            do i = 1, num_reaction_outputs
               j = reaction_outputs(2*i,ir)
               if (j == target_iso) then
                  rtab(ir) = 0 ! mark as removed
                  !write(*,*) 'remove reaction ' // trim(reaction_Name(ir))
                  exit
               end if
            end do
         end do
         
      end subroutine remove_reactions_for_iso


      subroutine setup_reaction_info(g, ierr) ! assumes have already called setup_iso_info
         use chem_def, only: chem_isos
         use rates_lib, only: get_weak_rate_id, is_weak_reaction
         use num_lib, only: qsort_string_index
         type (Net_General_Info), pointer  :: g
         integer, intent(out) :: ierr
         
         integer :: i, j, k, kind, reaction_num, num_reactions, &
            r1, r2, r3, &
            num_wk_reactions, ir, cid_lhs, cid_rhs, r_ir, r_i, &
            num_reaction_inputs, num_reaction_outputs, cid
         logical :: has_neut, has_prot, found_one
         integer, pointer, dimension(:) :: &
            rtab, index, ids, reaction_kind, &
            reaction_reaclib_kind, reaction_id, reverse_id_for_kind_ne_other
         
         include 'formats'
         
         ierr = 0
         rtab => g% net_reaction
         
         num_reactions = sum(rtab(:))
         g% num_reactions = num_reactions
         allocate( &
            g% reaction_kind(num_reactions), &
            g% reaction_reaclib_kind(num_reactions), &
            g% reverse_id_for_kind_ne_other(num_reactions), &
            g% reaction_id(num_reactions), &
            g% zs13(num_reactions), &
            g% zhat(num_reactions), &
            g% zhat2(num_reactions), &
            g% lzav(num_reactions), &
            g% aznut(num_reactions), &
            g% zs13inv(num_reactions), &
            g% rate_table(num_reactions, nrattab), &
            g% ttab(nrattab), &
            g% logttab(nrattab), &
            g% rattab_f1(4*nrattab*num_reactions), &
            stat=ierr)
         if (ierr /= 0) return
         reaction_reaclib_kind => g% reaction_reaclib_kind
         reaction_id => g% reaction_id
         reverse_id_for_kind_ne_other => g% reverse_id_for_kind_ne_other
         
         reaction_num = 0
         num_wk_reactions = 0    
         do ir = 1, rates_reaction_id_max
            if (rtab(ir) == 0) cycle ! reaction not in this net
            reaction_num = reaction_num + 1
            reaction_id(reaction_num) = ir
            if (is_weak_reaction(ir)) &
               num_wk_reactions = num_wk_reactions + 1
         end do

         if (reaction_num /= num_reactions) then
            write(*,*) 'reaction_num /= num_reactions', reaction_num, num_reactions
            stop 'setup_reaction_info'
         end if
         
         call check_for_hardwired_pairs
                  
         ! need to order reactions to ensure bit-for-bit results
         ! so sort reaction_id by reaction_Name
         index(1:num_reactions) => g% reaction_kind(1:num_reactions)
         ids(1:num_reactions) => g% reverse_id_for_kind_ne_other(1:num_reactions)
         do i=1, num_reactions
            ids(i) = reaction_id(i)
         end do
         call qsort_string_index(index,num_reactions,ids,reaction_Name)
         ! use reverse_id_for_kind_ne_other as temp for reordering reaction_id
         do i=1, num_reactions
            reaction_id(i) = ids(index(i))
            !write(*,*) trim(reaction_Name(reaction_id(i)))
         end do

         call set_reaction_kinds(g)
         
         do i = 1, num_reactions
            ir = reaction_id(i)
            if (ir < 1 .or. ir > rates_reaction_id_max) then
               write(*,*) '(ir < 1 .or. ir > rates_reaction_id_max)', &
                  ir, i, rates_reaction_id_max
               stop 'setup_reaction_info'
            end if
            rtab(ir) = i
         end do
         
         i = 1 
kind_loop: do kind = 1, max_kind ! reorder by kind of reaction; other_kind goes last
            do while (reaction_reaclib_kind(i) == kind)
               i = i+1
               if (i > num_reactions) exit kind_loop
            end do
            do ! move all of the instances of current kind
               found_one = .false.
               do j=i+1,num_reactions ! locate the next instance of current kind of reaction
                  if (reaction_reaclib_kind(j) /= kind) cycle
                  
                  ! exchange
                  
                  r1 = reaction_reaclib_kind(j)
                  r2 = reverse_id_for_kind_ne_other(j)
                  r3 = reaction_id(j)

                  reaction_reaclib_kind(j) = reaction_reaclib_kind(i)
                  reverse_id_for_kind_ne_other(j) = reverse_id_for_kind_ne_other(i)
                  reaction_id(j) = reaction_id(i)
                  
                  reaction_reaclib_kind(i) = r1
                  reverse_id_for_kind_ne_other(i) = r2
                  reaction_id(i) = r3

                  rtab(reaction_id(i)) = i
                  rtab(reaction_id(j)) = j
                  
                  found_one = .true.
                  exit
               end do
               
               if (.not. found_one) exit ! look for another kind
               i = i+1 ! next destination
               if (i > num_reactions) exit
            end do   
         end do kind_loop
         
         i = 0
         !g% have_all_reverses = .true.
         do ! reorder so that forward + reverse pairs are adjacent.
            i = i+1
            if (i >= num_reactions) exit
            if (reaction_reaclib_kind(i) == other_kind) cycle
            r_ir = reverse_id_for_kind_ne_other(i)
            if (r_ir <= 0) then
               write(*,*) trim(reaction_name(reaction_id(i)))
               write(*,2) 'reaction kind', reaction_reaclib_kind(i)
               stop 'setup_reaction_info: missing reverse id'
            end if
            r_i = rtab(r_ir)
            if (r_i == 0) then ! reverse reaction not in net
               !g% have_all_reverses = .false.
               cycle
            end if
            if (reaction_id(r_i) /= r_ir) stop 'setup_reaction_info: bad reverse'
            ir = reaction_id(i)
            if (r_i <= i) then
               write(*,2) 'r_i', r_i
               write(*,2) 'i', i
               write(*,2) 'reverse_id_for_kind_ne_other(r_i)', reverse_id_for_kind_ne_other(r_i)
               write(*,2) trim(reaction_Name(ir)), i
               write(*,2) trim(reaction_Name(r_ir)), r_i
               write(*,*) 'r_i <= i'
               stop
            end if
            
            i = i+1
            
            do k=r_i-1,i,-1
               reaction_reaclib_kind(k+1) = reaction_reaclib_kind(k)
               reverse_id_for_kind_ne_other(k+1) = reverse_id_for_kind_ne_other(k)
               reaction_id(k+1) = reaction_id(k)
               rtab(reaction_id(k+1)) = k+1
            end do
            
            reaction_reaclib_kind(i) = other_kind
            reverse_id_for_kind_ne_other(i) = ir
            reaction_id(i) = r_ir
            rtab(r_ir) = i
            
         end do
         
         g% num_wk_reactions = num_wk_reactions
         allocate( &
               g% weak_reaction_num(num_wk_reactions), &
               g% reaction_id_for_weak_reactions(num_wk_reactions), &
               g% weaklib_ids(num_wk_reactions), &
               g% weak_reaction_index(num_reactions), &
               stat=ierr)
         if (ierr /= 0) return

         j = 0
         do i = 1, reaction_num
            ir = reaction_id(i)
            if (.not. is_weak_reaction(ir)) then
               g% weak_reaction_index(i) = 0
               num_reaction_inputs = get_num_reaction_inputs(ir)
               num_reaction_outputs = get_num_reaction_outputs(ir)
               has_neut = .false.
               has_prot = .false.
               do k=1, num_reaction_inputs
                  cid = reaction_inputs(2*k,ir)
                  if (cid == ineut) has_neut = .true.
                  if (cid == iprot .or. cid == ih1) has_prot = .true.
               end do
               do k=1, num_reaction_outputs
                  cid = reaction_outputs(2*k,ir)
                  if (cid == ineut) has_neut = .true.
                  if (cid == iprot .or. cid == ih1) has_prot = .true.
               end do
               if (has_neut .and. .not. has_prot) then
                  g% reaction_kind(i) = neut_kind
               else if (has_prot) then
                  g% reaction_kind(i) = prot_kind
               else
                  g% reaction_kind(i) = other_strong_kind
               end if
               cycle
            end if
            g% reaction_kind(i) = weak_kind
            j = j+1
            g% weak_reaction_index(i) = j
            g% weak_reaction_num(j) = i
            g% reaction_id_for_weak_reactions(j) = ir
            cid_lhs = weak_reaction_info(1,ir)
            cid_rhs = weak_reaction_info(2,ir)
            g% weaklib_ids(j) =  &
               get_weak_rate_id(chem_isos% name(cid_lhs), chem_isos% name(cid_rhs))
         end do
         if (j /= num_wk_reactions) then
            write(*,3) 'problem with num_wk_reactions in setup_reaction_info', j, num_wk_reactions
            stop 'setup_reaction_info'
         end if

      end subroutine setup_reaction_info
      

      ! Fowler, Caughlan, Zimmerman, Annual Review Astro. Astrophys., 1975.12:69-112. eqn (1).
      real(dp) function neutrino_Q(i1, i2)
         integer, intent(in) :: i1, i2 ! i1 decays to i2
         real(dp) :: sum, sum2
         sum    = chem_isos% binding_energy(i2) - chem_isos% binding_energy(i1) - 0.782d0 - 1.022d0
         sum    = 1.0d0 + sum/0.511d0
         sum2   = sum*sum
         neutrino_Q = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
                     * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      end function neutrino_Q


      subroutine init_special_case_reaction_info
         integer :: i
         include 'formats'
         
         i = 0
         
         call set(ih1, ih2, 0, 0, 0, 'r_h1_h1_wk_h2', 'r_h1_h1_ec_h2')
         
         call set(ineut, ih1, ih2, 0, 0, 'r_neut_h1_h1_to_h1_h2', 'r_h1_h2_to_neut_h1_h1')
         
         call set(ih1, ih2, ih3, 0, 0, 'r_h2_h2_to_h1_h3', 'r_h1_h3_to_h2_h2')
         
         call set(ih3, ihe3, 0, 0, 0, 'r_he3_ec_h3', '')
         
         call set(ineut, ih2, ihe3, 0, 0, 'r_h2_h2_to_neut_he3', 'r_neut_he3_to_h2_h2')
         
         call set(ih1, ihe3, ihe4, 0, 0, 'r_h1_he3_wk_he4', '')
         
         call set(ih2, ihe4, 0, 0, 0, 'r_h2_h2_to_he4', 'r_he4_to_h2_h2')
         
         call set(ineut, ih2, ih3, ihe4, 0, 'r_h2_h3_to_neut_he4', 'r_neut_he4_to_h2_h3')
         
         call set(ih1, ih2, ihe3, ihe4, 0, 'r_h2_he3_to_h1_he4', 'r_h1_he4_to_h2_he3')
         
         call set(ih2, ih3, ihe3, ihe4, 0, 'r_h3_he3_to_h2_he4', 'r_h2_he4_to_h3_he3')
         
         call set(ineut, ih3, ihe4, 0, 0, 'r_h3_h3_to_neut_neut_he4', 'r_neut_neut_he4_to_h3_h3')
         
         call set(ih1, ihe3, ihe4, 0, 0, 'r_he3_he3_to_h1_h1_he4', 'r_h1_h1_he4_to_he3_he3')
         
         call set(ineut, ih1, ihe4, ili6, 0, 'r_li6_to_neut_h1_he4', 'r_neut_h1_he4_to_li6')
         
         call set(ineut, ih2, ihe4, ili7, 0, 'r_h2_li7_to_neut_he4_he4', &
            'r_neut_he4_he4_to_h2_li7')
         
         call set(ineut, ih3, ihe4, ili7, 0, &
               'r_h3_li7_to_neut_neut_he4_he4', 'r_neut_neut_he4_he4_to_h3_li7')
         
         call set(ih1, ih2, ili6, ili7, 0, 'r_h1_li7_to_h2_li6', 'r_h2_li6_to_h1_li7')
         
         call set(ibe7, ili7, 0, 0, 0, 'r_be7_wk_li7', '')
         
         ! might also need reverse for r_n13_wk_c13, r_o15_wk_n15
         
         call set(ih1, ih2, ihe4, ibe7, 0, 'r_h2_be7_to_h1_he4_he4', 'r_h1_he4_he4_to_h2_be7')
         
         call set(ineut, ihe4, ibe7, 0, 0, 'r_neut_be7_to_he4_he4', 'r_he4_he4_to_neut_be7')
         
         call set(ih1, ihe3, ihe4, ibe7, 0, 'r_he3_be7_to_h1_h1_he4_he4', &
            'r_h1_h1_he4_he4_to_he3_be7')
         
         call set(ineut, ih2, ili6, ibe7, 0, 'r_neut_be7_to_h2_li6', 'r_h2_li6_to_neut_be7')
         
         call set(ineut, ih3, ili7, ibe9, 0, 'r_neut_be9_to_h3_li7', 'r_h3_li7_to_neut_be9')
         
         call set(ineut, ihe4, ibe9, 0, 0, 'r_be9_to_neut_he4_he4', 'r_neut_he4_he4_to_be9')
         
         call set(ih1, ih2, ihe4, ibe9, 0, 'r_h1_be9_to_h2_he4_he4', 'r_h2_he4_he4_to_h1_be9')
         
         call set(ineut, ih1, ihe4, ib8, 0, 'r_neut_b8_to_h1_he4_he4', 'r_h1_he4_he4_to_neut_b8')
         
         call set(ih1, ihe4, ib11, 0, 0, 'r_h1_b11_to_he4_he4_he4', 'r_he4_he4_he4_to_h1_b11')
         
         call set(ineut, ih3, ibe9, ib11, 0, 'r_neut_b11_to_h3_be9', 'r_h3_be9_to_neut_b11')
         
         call set(ih1, ihe4, ic9, 0, 0, 'r_c9_wk_h1_he4_he4', '')
         
         call set(ineut, ihe4, ic11, 0, 0, 'r_neut_c11_to_he4_he4_he4', &
            'r_he4_he4_he4_to_neut_c11')
         
         call set(ihe4, ic12, 0, 0, 0, 'r_he4_he4_he4_to_c12', 'r_c12_to_he4_he4_he4')
         
         call set(ineut, ih2, ic13, in14, 0, 'r_neut_n14_to_h2_c13', 'r_h2_c13_to_neut_n14')
         
         call set(ineut, ih2, ic14, in15, 0, 'r_neut_n15_to_h2_c14', 'r_h2_c14_to_neut_n15')
         
         call set(ih1, ihe4, io13, io15, 0, 'r_he4_o13_to_h1_h1_o15', 'r_h1_h1_o15_to_he4_o13')
         
         call set(ihe4, ic12, ine20, 0, 0, 'r_c12_c12_to_he4_ne20', 'r_he4_ne20_to_c12_c12')
         
         call set(ih1, ic12, ina23, 0, 0, 'r_c12_c12_to_h1_na23', 'r_h1_na23_to_c12_c12')
         
         call set(ineut, ic12, img23, 0, 0, 'r_neut_mg23_to_c12_c12', 'r_c12_c12_to_neut_mg23')
         
         call set(ihe4, ic12, io16, img24, 0, 'r_c12_o16_to_he4_mg24', 'r_he4_mg24_to_c12_o16')
         
         call set(ih1, ic12, io16, ial27, 0, 'r_c12_o16_to_h1_al27', 'r_h1_al27_to_c12_o16')
         
         call set(ineut, ic12, io16, isi27, 0, 'r_neut_si27_to_c12_o16', 'r_c12_o16_to_neut_si27')
         
         call set(ihe4, io16, isi28, 0, 0, 'r_o16_o16_to_he4_si28', 'r_he4_si28_to_o16_o16')
         
         call set(ihe4, ic12, ine20, isi28, 0, 'r_c12_ne20_to_he4_si28', 'r_he4_si28_to_c12_ne20')
         
         call set(ih1, io16, ip31, 0, 0, 'r_o16_o16_to_h1_p31', 'r_h1_p31_to_o16_o16')
         
!         call set(ih2, io16, ip30, 0, 0, 'r_o16_o16_to_h2_p30', 'r_h2_p30_to_o16_o16')
         
         call set(ih1, ic12, ine20, ip31, 0, 'r_c12_ne20_to_h1_p31', 'r_h1_p31_to_c12_ne20')
         
         call set(ih1, ial25, is27, 0, 0, 'r_s27_wk_h1_h1_al25', '')
         
         call set(ineut, io16, is31, 0, 0, 'r_o16_o16_to_neut_s31', 'r_neut_s31_to_o16_o16')
         
         call set(ineut, ic12, ine20, is31, 0, 'r_c12_ne20_to_neut_s31', 'r_neut_s31_to_c12_ne20')
         
         call set(ih1, ip29, iar31, 0, 0, 'r_ar31_wk_h1_h1_p29', '')
         
         call set(ih1, ihe4, ica36, ica38, 0, 'r_he4_ca36_to_h1_h1_ca38', &
            'r_h1_h1_ca38_to_he4_ca36')
         
         call set(ineut, ih1, ih3, ihe3, ihe4, 'r_h3_he3_to_neut_h1_he4', '')
         
         call set(ineut, ih1, ihe3, ihe4, ili7, 'r_he3_li7_to_neut_h1_he4_he4', '')
         
         call set(ineut, ih1, ih3, ihe4, ibe7, 'r_h3_be7_to_neut_h1_he4_he4', '')
         
         num_special_case_reactions = i
         
         !write(*,2) 'num_special_case_reactions', num_special_case_reactions
         
         
         contains
         
         subroutine set(i1, i2, i3, i4, i5, s1_in, s2_in)
            integer, intent(in) :: i1, i2, i3, i4, i5
            character (len=*), intent(in) :: s1_in, s2_in
            character (len=maxlen_reaction_Name) :: s1, s2
            i = i+1
            special_case_reactants(1,i) = i1
            special_case_reactants(2,i) = i2
            special_case_reactants(3,i) = i3
            special_case_reactants(4,i) = i4
            special_case_reactants(5,i) = i5
            s1 = s1_in
            s2 = s2_in
            special_case_reactions(1,i) = s1      
            special_case_reactions(2,i) = s2    
         end subroutine set
      
      end subroutine init_special_case_reaction_info

      

      end module net_initialize

