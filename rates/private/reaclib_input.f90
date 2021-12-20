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
 
      module reaclib_input
      use rates_def
      use utils_lib, only: mesa_error      

      implicit none

      type(reaclib_data) :: reaclib
      integer :: nreaclib

      contains
      
      
      subroutine do_extract_rates(set,nuclides,rates,use_weaklib,ierr)
         use reaclib_support
         use chem_def, only: nuclide_set, nuclide_data
         type(nuclide_set), dimension(:), intent(in) :: set
         type(nuclide_data), intent(in), target :: nuclides
         logical, intent(in) :: use_weaklib
         type(reaction_data), intent(out) :: rates
         integer, intent(out) :: ierr
         logical, parameter :: dbg = .false.

         if (dbg) write(*,*) 'call extract_rates_from_reaclib', nreaclib
         call extract_rates_from_reaclib(reaclib,nreaclib,nuclides,rates,set,use_weaklib,ierr)
         if (failed('extract_rates_from_reaclib')) return
         
         if (dbg) write(*,*) 'call set_up_network_information'
         call set_up_network_information(rates)
         
         if (dbg) write(*,*) 'call assign_weights'
         call assign_weights(rates)
         
         if (dbg) write(*,*) 'call compute_rev_ratio'
         call compute_rev_ratio(rates,nuclides)
         
         if (dbg) write(*,*) 'return from do_extract_rates'

         contains
         
         logical function failed(msg)
            character (len=*), intent(in) :: msg
            if (ierr /= 0) then
               failed = .true.
               write(*,*) 'do_extract_rates failed in ' // trim(msg)
            else
               failed = .false.
            end if
         end function failed
         
      end subroutine do_extract_rates
      
      
      subroutine do_read_reaclib(ierr)
         use utils_lib, only: integer_dict_define, integer_dict_create_hash
         use math_lib, only: str_to_double
         integer, intent(out) :: ierr
         integer :: i, j, count, reaclib_unitno, cache_io_unit, ios
         ! file and table format
         character(len=*), parameter :: line0 = '(i2)'
         character(len=*), parameter :: line1 = '(5x,6a5,8x,a4,a1,a1,3x,a12)'
         !character(len=*), parameter :: line2 = '(4es13.6)'
         !character(len=*), parameter :: line3 = '(3es13.6)'
         character(len=iso_name_length) :: species
         character(len=256) :: filename, cache_filename, buf
         character(len=12) :: Qvalue_str
         
         logical, parameter :: use_cache = .true.
         
         include 'formats'

         ierr = 0
         
         cache_filename = trim(rates_cache_dir) // '/jina_reaclib.bin'
         
         filename = trim(reaclib_dir) // '/' //trim(reaclib_filename)
         open(newunit=reaclib_unitno, file=filename, iostat=ierr, status="old", action="read")
         if ( ierr /= 0 ) then
            write(*,*) 'Error opening ' // filename
            stop
         end if

         ! allocate the library
         call allocate_reaclib_data(reaclib,max_nreaclib,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'unable to allocate storage for reaclib')

         if (use_cache) then
            ios = 0
            open(newunit=cache_io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
            if (ios == 0) then ! opened it okay
               call read_reaclib_cache(cache_io_unit,ios)
               close(cache_io_unit)
               if (ios == 0) then ! read it okay
                  close(reaclib_unitno)
                  return
               end if
            end if
         end if

         count = 0
         do i = 1, max_nreaclib
            read(unit=reaclib_unitno, fmt=line0, iostat=ierr) reaclib% chapter(i)
            if (ierr /= 0 ) then ! assume end of file
               ierr = 0; exit 
            end if
            read(unit=reaclib_unitno,fmt=line1,iostat=ierr,err=100) &
               reaclib% species(1,i), reaclib% species(2,i), reaclib% species(3,i), &
               reaclib% species(4,i), reaclib% species(5,i), reaclib% species(6,i), &
               reaclib% label(i), reaclib% reaction_flag(i), reaclib% reverse_flag(i), Qvalue_str
            call str_to_double(Qvalue_str, reaclib% Qvalue(i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            do j=1,6
               species = adjustl(reaclib% species(j,i))
               if (species == 'p') then
                  reaclib% species(j,i) = 'h1'
               else if (species == 'n') then
                  reaclib% species(j,i) = 'neut'
               else if (species == 'd') then
                  reaclib% species(j,i) = 'h2'
               else if (species== 't') then
                  reaclib% species(j,i) = 'h3'
               else if (species== 'al-6') then
                  reaclib% species(j,i) = 'al26-1'
               else if (species== 'al*6') then
                  reaclib% species(j,i) = 'al26-2'
               end if
            end do
            read(unit=reaclib_unitno,fmt='(a)',iostat=ierr,err=100) buf
            call str_to_double(buf(1:13), reaclib% coefficients(1,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call str_to_double(buf(14:26), reaclib% coefficients(2,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call str_to_double(buf(27:39), reaclib% coefficients(3,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call str_to_double(buf(40:52), reaclib% coefficients(4,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            read(unit=reaclib_unitno,fmt='(a)',iostat=ierr,err=100) buf
            call str_to_double(buf(1:13), reaclib% coefficients(5,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call str_to_double(buf(14:26), reaclib% coefficients(6,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            call str_to_double(buf(27:39), reaclib% coefficients(7,i), ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            count = count + 1
         end do
         
         nreaclib = count
         close(reaclib_unitno)

         if (use_cache) then
            open(newunit=cache_io_unit, file=trim(cache_filename), iostat=ios, &
                  action='write', form='unformatted')
            if (ios == 0) then
               !write(*,'(a)') 'write ' // trim(cache_filename)
               call write_reaclib_cache(cache_io_unit)
               close(cache_io_unit)
            end if
         end if

         return
      100 call mesa_error(__FILE__,__LINE__,'error in do_read_reaclib')
      
      
         contains
         
         
         subroutine read_reaclib_cache(io,ios)
            integer, intent(in) :: io
            integer, intent(out) :: ios
            integer :: i,n            
            ios = 0
            read(io,iostat=ios) nreaclib
            if (ios /= 0) return
            n = nreaclib
            read(io,iostat=ios) &
               reaclib% chapter(1:n), &
               reaclib% species(1:max_species_per_reaction,1:n), &
               reaclib% label(1:n), &
               reaclib% reaction_flag(1:n), &
               reaclib% reverse_flag(1:n), &
               reaclib% Qvalue(1:n), &
               reaclib% coefficients(1:ncoefficients,1:n) 
            if (ios /= 0) return          
         end subroutine read_reaclib_cache
         
         
         subroutine write_reaclib_cache(io)
            integer, intent(in) :: io
            integer :: i, n
            write(io) nreaclib
            n = nreaclib
            write(io) &
               reaclib% chapter(1:n), &
               reaclib% species(1:max_species_per_reaction,1:n), &
               reaclib% label(1:n), &
               reaclib% reaction_flag(1:n), &
               reaclib% reverse_flag(1:n), &
               reaclib% Qvalue(1:n), &
               reaclib% coefficients(1:ncoefficients,1:n) 
         end subroutine write_reaclib_cache         
         
         
      end subroutine do_read_reaclib
      

      subroutine extract_rates_from_reaclib(reaclib,nreaclib,nuclides,rates,set,use_weaklib,ierr)
         use chem_def, only: &
            nuclide_set, nuclide_not_found, del_Mp, del_Mn, &
            ipp, icno, i3alf, i_burn_c, i_burn_n, i_burn_o, i_burn_ne, i_burn_na, i_burn_mg, &
            i_burn_si, i_burn_s, i_burn_ar, i_burn_ca, i_burn_ti, i_burn_cr, i_burn_fe, icc, &
            ico, ioo, iphoto, iother
         use chem_lib, only : get_nuclide_index_in_set, get_Q
         use utils_lib, only: integer_dict_define, integer_dict_create_hash, integer_dict_lookup
         use reaclib_support, only: get_reaction_handle, get_reverse_reaction_handle
         type(reaclib_data), intent(in) :: reaclib
         integer, intent(in) :: nreaclib
         type(nuclide_data), intent(in), target :: nuclides
         type(reaction_data), intent(out) :: rates
         type(nuclide_set), dimension(:), intent(in) :: set
         logical, intent(in) :: use_weaklib
         integer, intent(out) :: ierr
         type(reaction_data) :: r ! temporary storage
         integer :: i,j,l,count,loc_count,nt,indx,cat, &
            weaklib_count,chapter,num_in,num_out,num_skip_for_weaklib,num_from_reaclib, &
            max_lhs_Z, min_Z, i1, i2, cid_in, cid_out
         logical :: include_this_rate, already_included_from_weaklib, found_it
         integer, dimension(max_species_per_reaction) :: pspecies
         character(len=max_id_length) :: handle
         character(len=iso_name_length) :: name_i, name_j, name1, name2, label
         integer, parameter :: max_weaklib_rates = 1000
         real(dp) :: Q
         integer :: rate_ipp, rate_ipep

         logical, parameter :: dbg = .false.
         
         include 'formats'

         ierr = 0
         
         if (dbg) write(*,*) 'call allocate_reaction_data'
         call allocate_reaction_data(r,max_nreaclib,max_weaklib_rates,ierr)
         if (ierr /= 0) then
            print *,'unable to allocate temporary storage for rates'
            return
         end if

         count = 0         
         if (use_weaklib) then ! add weaklib rates first
            if (dbg) write(*,*) 'add weaklib rates'
            do i = 1, nuclides% nnuclides
               do j = 1, nuclides% nnuclides
                  if (nuclides% Z(i)+nuclides% N(i) /= nuclides% Z(j)+nuclides% N(j)) cycle
                  if (abs(nuclides% Z(i) - nuclides% Z(j)) /= 1) cycle
                  ! check if this one is in weaklib
                  call get_weaklib_name(i, name_i)
                  call get_weaklib_name(j, name_j)
                  indx = do_get_weak_rate_id(name_i, name_j)
                  if (indx == 0) cycle
                  count = count + 1
                  r% weaklib_ids(count) = indx
                  r% also_in_reaclib(count) = .false.
                  r% chapter(count) = 1
                  r% pspecies(1,count) = get_nuclide_index_in_set(name_i,set) 
                  r% pspecies(2,count) = get_nuclide_index_in_set(name_j,set)
                  r% pspecies(3:max_species_per_reaction,count) = 0
                  r% coefficients(:,count) = 0
                  r% coefficients(1,count) = -99
                  r% reaction_flag(count) = ''
                  call get_reaction_handle( &
                     1, 1, r% pspecies(:,count), nuclides, &
                     r% reaction_flag(count), r% reaction_handle(count))
                  r% reverse_handle(count) = '' 
                  r% Q(count) = 0
                  r% Qneu(count) = 0
                  r% reaction_flag(count) = 'w'
                  if (.false.) write(*,*) 'found weaklib rate for    ' // name_i &
                        // name_j // ' ' // trim(r% reaction_handle(count)), count
               end do
            end do
         end if
         
         weaklib_count = count         
         num_skip_for_weaklib = 0
         num_from_reaclib = 0
         loc_count = 0
         rate_ipp = 0; rate_ipep = 0
                  
         if (dbg) write(*,*) 'loop_over_rates'
         loop_over_rates: do i = 1,nreaclib
            include_this_rate = .true.
            pspecies(:) = 0
            chapter = reaclib% chapter(i)
            num_in = Nin(chapter)
            num_out = Nout(chapter)
            nt = num_in+num_out
            loop_over_nuclides: do j = 1, nt
               l = get_nuclide_index_in_set(reaclib% species(j,i),set)
               if (l == nuclide_not_found) then
                  include_this_rate = .false.
                  exit
               else
                  pspecies(j) = l
               end if

            end do loop_over_nuclides

            ! only include forward rates
            ! Define the reverse rate as being the endothermic reaction, always
            ! Some photo disintegrations can be exothermic, so dont trust what REACLIB calls a reverse rate
            if(include_this_rate) then
               if(sum(nuclides% binding_energy(pspecies(1:num_in))) - &
                  sum(nuclides% binding_energy(pspecies(num_in+1:Nt)))  > 0 ) cycle loop_over_rates
            end if
            if (include_this_rate) then
               if (use_weaklib .and. num_in == 1 .and. num_out == 1) then
                  already_included_from_weaklib = .false.
                  name1 = adjustl(reaclib% species(1,i))
                  name2 = adjustl(reaclib% species(2,i))
                  indx = do_get_weak_rate_id(name1, name2)
                  if (indx /= 0) then
                     already_included_from_weaklib = .true.
                  else
                     indx = do_get_weak_rate_id(name2, name1)
                     if (indx /= 0) then
                        already_included_from_weaklib = .true.
                     end if
                  end if
                  if (already_included_from_weaklib) then ! find it and store coefficients
                     found_it = .false.
                     do j=1,weaklib_count
                        if (r% weaklib_ids(j) == indx) then
                           r% coefficients(:,j) = reaclib% coefficients(:,i)
                           r% also_in_reaclib(j) = .true.
                           found_it = .true.
                           exit
                        end if
                     end do
                     if (.not. found_it) then
                        write(*,*) 'problem in reaclib_input'
                        write(*,*) 'failed to find', indx
                        call mesa_error(__FILE__,__LINE__)
                     end if
                     cycle
                  end if
               end if
               count = count + 1
               num_from_reaclib = num_from_reaclib + 1
               r% chapter(count) = chapter
               r% pspecies(:,count) = pspecies(:)
               r% reaction_flag(count) = reaclib% reaction_flag(i)
               r% coefficients(:,count) = reaclib% coefficients(:,i)
               ! mark the ec rates so they'll get an extra factor of Ye*Rho
               label = adjustl(reaclib% label(i))
               if (label(1:2) == 'ec') r% reaction_flag(count) = 'e'
               call get_reaction_handle( &
                  num_in, num_out, r% pspecies(:,count), nuclides, &
                  r% reaction_flag(count), handle)
               r% reaction_handle(count) = handle
               if (r% reaction_flag(count) == 'e' .or. r% reaction_flag(count) == 'w') then
                  r% reverse_handle(count) = ''
               else
                  call get_reverse_reaction_handle( &
                     num_in, num_out, r% pspecies(:,count), nuclides, r% reverse_handle(count))
               end if

               r% Q(count) = 0
               do j=1,num_in+num_out
                  l = pspecies(j)
                  Q = get_Q(nuclides,l)
                  if (j <= num_in) then
                     r% Q(count) = r% Q(count) - Q
                  else
                     r% Q(count) = r% Q(count) + Q
                  end if
               end do
               
               r% Qneu(count) = 0
               
               if (handle == 'r_b8_to_he4_he4') then
                  r% Qneu(count) = 0.6735D+01
                  
               else if (handle == 'r_h1_he3_to_he4') then
                  r% Qneu(count) = 9.628D0
                  
               else if (handle == 'r_h1_h1_ec_h2') then
                  rate_ipep = count
                  r% Qneu(count) = 1.445D0
                  
               else if (handle == 'r_h1_h1_wk_h2') then
                  rate_ipp = count
                  r% Qneu(count) = 0.2668D0
                  
               else if (handle == 'r_he3_ec_h3') then
                  r% Qneu(count) = 10D0 ! who knows?  who cares?
                  
               else if (adjustl(reaclib% reaction_flag(i)) == 'w') then ! check weak_info list
                  name1 = reaclib% species(1,i)
                  if (num_out == 1) then
                     name2 = reaclib% species(2,i)
                  else if (num_out == 2) then
                     name2 = reaclib% species(3,i)
                  else
                     name2 = ''
                  end if
                  j = do_get_weak_info_list_id(name1, name2)
                  if (j > 0) r% Qneu(count) = weak_info_list_Qneu(j)
                  if (.false. .and. num_out == 2) &
                     write(*,2) trim(name1) // ' ' // trim(name2) // &
                        ' ' // trim(handle) // ' Qneu', count, r% Qneu(count)
               end if
               
               if (reaclib% reaction_flag(i) == 'w' .and. .false.) then
                  write(*,2) 'reaclib weak ' // trim(handle) // ' Qneu', count, r% Qneu(count)
               end if

            end if
         end do loop_over_rates
         
         if (.false.) then
            write(*,2) 'num_skip_for_weaklib', num_skip_for_weaklib
            write(*,2) 'weaklib_count', weaklib_count
            write(*,2) 'num_from_reaclib', num_from_reaclib
            write(*,2) 'reaclib+weaklib', num_from_reaclib + weaklib_count
            write(*,2) 'total num reactions', count
            call mesa_error(__FILE__,__LINE__,'extract_rates_from_reaclib')
         end if
         

         ! we can now stuff our temporary file into the output and discard the temporary
         
         if (dbg) write(*,*) 'call allocate_reaction_data'
         call allocate_reaction_data(rates,count,weaklib_count,ierr)
         if (ierr /= 0) then
            print *,'unable to allocate storage for rates'
            return
         end if
         
         rates% nreactions = count
         rates% nuclides => nuclides

         rates% num_from_weaklib = weaklib_count
         do i=1,weaklib_count
            rates% weaklib_ids(i) = r% weaklib_ids(i)
            rates% also_in_reaclib(i) = r% also_in_reaclib(i)
         end do
         
         do i=1,count
            rates% reaction_handle(i) = r% reaction_handle(i)
            rates% reverse_handle(i) = r% reverse_handle(i)
            rates% chapter(i) = r% chapter(i)
            rates% reaction_flag(i) = r% reaction_flag(i)
            rates% Q(i) = r% Q(i)
            rates% Qneu(i) = r% Qneu(i)
            do j=1,max_species_per_reaction
               rates% pspecies(j,i) = r% pspecies(j,i)
            end do
            do j=1,ncoefficients
               rates% coefficients(j,i) = r% coefficients(j,i)
            end do
         end do
         
         nullify(rates% reaction_dict)
         nullify(rates% reverse_dict)
         do i=1,count
            chapter = rates% chapter(i)
            num_in = Nin(chapter)
            num_out = Nout(chapter)
            max_lhs_Z = -1
            do j=1,num_in
               max_lhs_Z = max(max_lhs_Z, nuclides% Z(rates% pspecies(j,i)))
            end do
            min_Z = 999999
            do j=1,num_in+num_out
               min_Z = min(min_Z, nuclides% Z(rates% pspecies(j,i)))
            end do
            cat = -1   
            ! NOTE: reaction categories that are used by net are set in rates_initialize
            if (chapter == r_one_one) then
               if (min_Z > 0) then
                  if (max_lhs_Z <= 5) then
                     cat = ipp
                  else if (max_lhs_Z <= 10) then
                     cat = icno
                  end if
               end if
            else if (chapter == r_two_one .or. chapter == r_two_two) then
               if (rates% reaction_handle(i) == 'r_he4_c12_to_o16') then
                  cat = i_burn_c
               else if (rates% reaction_handle(i) == 'r_c12_ag_o16') then
                  cat = i_burn_c
               else if (rates% reaction_handle(i) == 'r_o16_ag_ne20') then
                  cat = i_burn_o
               else if (rates% reaction_handle(i) == 'r_he4_o16_to_ne20') then
                  cat = i_burn_o
               else if (rates% reaction_handle(i) == 'r_c12_c12_to_mg24') then
                  cat = icc
               else if (rates% reaction_handle(i) == 'r_c12_c12_to_he4_ne20') then
                  cat = icc
               else if (rates% reaction_handle(i) == 'r_c12_c12_to_p_na23') then
                  cat = icc
               else if (rates% reaction_handle(i) == 'r_c12_o16_to_si28') then
                  cat = ico
               else if (rates% reaction_handle(i) == 'r_c12_o16_to_he4_mg24') then
                  cat = ico
               else if (rates% reaction_handle(i) == 'r_c12_o16_to_p_al27') then
                  cat = ico
               else if (rates% reaction_handle(i) == 'r_o16_o16_to_s32') then
                  cat = ioo
               else if (rates% reaction_handle(i) == 'r_o16_o16_to_he4_si28') then
                  cat = ioo
               else if (rates% reaction_handle(i) == 'r_o16_o16_to_p_p31') then
                  cat = ioo
               else if (min_Z > 0) then
                  if (max_lhs_Z <= 5) then
                     cat = ipp
                  else if (max_lhs_Z <= 9) then
                     cat = icno
                  end if
               end if
            else if (chapter == r_two_three) then
               if (min_Z > 0 .and. max_lhs_Z <= 5) then
                  cat = ipp
               end if
            else if (chapter == r_three_one) then
               if (rates% reaction_handle(i) == 'r_he4_he4_he4_to_c12') then
                  cat = i3alf
               end if
            else if (chapter == r_three_two) then
               if (rates% reaction_handle(i) == 'r_h1_h1_he4_to_he3_he3') then
                  cat = ipp
               end if
            else if (chapter == r_one_two .or. chapter == r_one_three .or. chapter == r_one_four) then
               if (rates% reaction_handle(i) == 'r_b8_to_he4_he4') then
                  cat = ipp
               else
                  cat = iphoto
               end if
            end if
            if (cat == -1) then
               if (max_lhs_Z <= 5) then
                  if (min_Z > 0) then
                     cat = ipp
                  else
                     cat = iother
                  end if
               else if (max_lhs_Z <= 6) then
                  cat = i_burn_c
               else if (max_lhs_Z <= 7) then
                  cat = i_burn_n
               else if (max_lhs_Z <= 8) then
                  cat = i_burn_o
               else if (max_lhs_Z <= 10) then
                  cat = i_burn_ne
               else if (max_lhs_Z <= 11) then
                  cat = i_burn_na
               else if (max_lhs_Z <= 12) then
                  cat = i_burn_mg
               else if (max_lhs_Z <= 14) then
                  cat = i_burn_si
               else if (max_lhs_Z <= 16) then
                  cat = i_burn_s
               else if (max_lhs_Z <= 18) then
                  cat = i_burn_ar
               else if (max_lhs_Z <= 20) then
                  cat = i_burn_ca
               else if (max_lhs_Z <= 22) then
                  cat = i_burn_ti
               else if (max_lhs_Z <= 24) then
                  cat = i_burn_cr
               else if (max_lhs_Z <= 28) then
                  cat = i_burn_fe
               else 
                  cat = iother
               end if
            end if
            rates% category(i) = cat
            
            call integer_dict_define(rates% reaction_dict, rates% reaction_handle(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: extract_rates_from_reaclib failed in integer_dict_define'
               call mesa_error(__FILE__,__LINE__)
            end if
            if (len_trim(rates% reverse_handle(i)) > 0) then
               call integer_dict_define(rates% reverse_dict, rates% reverse_handle(i), i, ierr)
               if (ierr /= 0) then
                  write(*,*) 'FATAL ERROR: extract_rates_from_reaclib failed in integer_dict_define'
                  call mesa_error(__FILE__,__LINE__)
               end if
            end if
         end do
         
         if (dbg) write(*,*) 'call integer_dict_create_hash reaction_dict'
         call integer_dict_create_hash(rates% reaction_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: extract_rates_from_reaclib failed in integer_dict_create_hash'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (dbg) write(*,*) 'call integer_dict_create_hash reverse_dict'
         call integer_dict_create_hash(rates% reverse_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: extract_rates_from_reaclib failed in integer_dict_create_hash'
            call mesa_error(__FILE__,__LINE__)
         end if

         if (dbg) write(*,*) 'call free_reaction_data'
         call free_reaction_data(r)

         ! don't allow blanks in reaction flag
         where (rates% reaction_flag == ' ') rates% reaction_flag = '-'

         if (dbg) write(*,*) 'done extract_rates_from_reaclib'
         
         
         contains
         
         
         subroutine get_weaklib_name(i,name)
            integer, intent(in) :: i
            character (len=iso_name_length), intent(out) :: name
            if (nuclides% Z(i) == 0) then
               name = 'neut'
            else if (nuclides% Z(i) == 1 .and. nuclides% N(i) == 0) then
               name = 'h1'
            else
               name = nuclides% name(i)
            end if
         end subroutine get_weaklib_name


      end subroutine extract_rates_from_reaclib


      ! Fowler, Caughlan, Zimmerman, Annual Review Astro. Astrophys., 1975.12:69-112. eqn (1).
      real(dp) function neutrino_Q(b1, b2)
         use chem_def, only: chem_isos
         real(dp), intent(in) :: b1, b2
         real(dp) :: sum, sum2
         sum    = b2 - b1 - 0.782d0 - 1.022d0
         sum    = 1.0d0 + sum/0.511d0
         sum2   = sum*sum
         neutrino_Q = abs(0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
                     * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2)))
      end function neutrino_Q
      
      

      end module reaclib_input
