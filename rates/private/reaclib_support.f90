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
 
      module reaclib_support
      use rates_def
      use math_lib
      use chem_lib
      
      implicit none

      contains
      

      subroutine set_up_network_information(rates)
         type(reaction_data), intent(inout) :: rates
         integer :: current_chapter, i

         ! set up bookmarks
         current_chapter = 0
         rates% nchapters_included = 0
         rates% chapters_present = 0
         rates% bookmarks = 0
         do i = 1, rates% nreactions
            new_chapter : if (rates% chapter(i) /= current_chapter) then
               ! close out the chapter we just left.
               if (current_chapter /= 0) rates% bookmarks(2,current_chapter) = i-1
               ! set up information on the new chapter
               current_chapter = rates% chapter(i)
               rates% nchapters_included = rates% nchapters_included + 1
               rates% chapters_present(rates% nchapters_included) = current_chapter
               rates% bookmarks(1,current_chapter) = i
            end if new_chapter
         end do
         ! mark the end of the last chapter
         rates% bookmarks(2,current_chapter) = rates% nreactions
      end subroutine set_up_network_information
      

      subroutine assign_weights(rates)
         type(reaction_data), intent(inout) :: rates
         integer :: i, i1, i2, i3, i4
         
         include 'formats.dek'

         ! check for allocation
         if (.not.associated(rates% weight)) then
            return
         end if
         
         do i = 1, rates% nreactions
            i1 = -1; i2 = -2; i3 = -3; i4 = -4
            select case (rates% chapter(i))
               case (r_one_one)
               case (r_one_two)
               case (r_one_three)
               case (r_two_one)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
               case (r_two_two)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
               case (r_two_three)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
               case (r_two_four)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
               case (r_three_one)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
                  i3 = rates% pspecies(3,i)
               case (r_three_two)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
                  i3 = rates% pspecies(3,i)
               case (r_four_two)
                  i1 = rates% pspecies(1,i)
                  i2 = rates% pspecies(2,i)
                  i3 = rates% pspecies(3,i)
                  i4 = rates% pspecies(4,i)
               case (r_one_four)
            end select
            call set_weight(rates% weight(i))
         end do
         
         do i = 1, rates% nreactions
            i1 = -1; i2 = -2; i3 = -3; i4 = -4
            select case (rates% chapter(i))
               case (r_one_one)
               case (r_one_two)
                  i1 = rates% pspecies(2,i)
                  i2 = rates% pspecies(3,i)
               case (r_one_three)
                  i1 = rates% pspecies(2,i)
                  i2 = rates% pspecies(3,i)
                  i3 = rates% pspecies(4,i)
               case (r_two_one)
               case (r_two_two)
                  i1 = rates% pspecies(3,i)
                  i2 = rates% pspecies(4,i)
               case (r_two_three)
                  i1 = rates% pspecies(3,i)
                  i2 = rates% pspecies(4,i)
                  i3 = rates% pspecies(5,i)
               case (r_two_four)
                  i1 = rates% pspecies(3,i)
                  i2 = rates% pspecies(4,i)
                  i3 = rates% pspecies(5,i)
                  i4 = rates% pspecies(6,i)
               case (r_three_one)
               case (r_three_two)
                  i1 = rates% pspecies(4,i)
                  i2 = rates% pspecies(5,i)
               case (r_four_two)
                  i1 = rates% pspecies(5,i)
                  i2 = rates% pspecies(6,i)
               case (r_one_four)
                  i1 = rates% pspecies(2,i)
                  i2 = rates% pspecies(3,i)
                  i3 = rates% pspecies(4,i)
                  i4 = rates% pspecies(5,i)
            end select

            call set_weight(rates% weight_reverse(i))

         end do
         
         
         contains
         
         
         subroutine set_weight(w)
            ! nuclei are sorted, so if identical, then are adjacent in list
            real(dp), intent(out) :: w
            if (i1 == i2 .and. i2 == i3 .and. i3 == i4) then
               w = 1d0/24d0
            else if (i2 == i3 .and. (i1 == i2 .or. i3 == i4)) then
               w = 1d0/6d0
            else if (i1 == i2) then
               if (i3 == i4) then
                  w = 1d0/4d0
               else
                  w = 1d0/2d0
               end if
            else if (i2 == i3 .or. i3 == i4) then
               w = 1d0/2d0
            else
               w = 1d0
            end if
         end subroutine set_weight
         

      end subroutine assign_weights
      

      subroutine compute_rev_ratio(rates,winvn)
         use const_def, only : pi, kB=>boltzm, NA=>avo, hbar, &
            c=>clight, conv=>mev_to_ergs
         type(reaction_data), intent(inout) :: rates
         type(nuclide_data), intent(in) :: winvn
         real(dp) :: mp, mn
         real(dp),  dimension(max_species_per_reaction) :: g   ! statistical weights of nuclides
         real(dp), dimension(max_species_per_reaction) :: mass ! mass no's of nuclides
         integer, dimension(max_species_per_reaction) :: ps
         integer :: Ni,No,Nt,i
         real(dp) :: fac, massfac, sum1, sum2, tmp
      
         
         include 'formats.dek'
         
         ! Get these consistently from the isotopes.data file
         mp=winvn%W(chem_get_iso_id('prot'))
         mn=winvn%W(chem_get_iso_id('neut'))
         
         fac = pow(1d9*kB/(2d0*pi*hbar*hbar*NA),1.5d0)/NA
         massfac = conv*NA/(c*c)

         rates% weak_mask = 1d0
         loop_over_rates: do i = 1,rates% nreactions
            ! weak rates don't have inverses, so set to innocuous values and mask out
            if (rates% reaction_flag(i) == 'w' .or. rates% reaction_flag(i) == 'e') then
               rates% weak_mask(i) = 0d0
               rates% inverse_coefficients(:,i) = (/-huge(1d0), 0d0/)
               rates% inverse_exp(i) = 0d0
               rates% inverse_part(:,i) = 1d0
               cycle
            end if
            Ni = Nin(rates% chapter(i))
            No = Nout(rates% chapter(i))
            Nt = Ni+No
            ps(1:Nt) = rates% pspecies(1:Nt,i)
            g(1:Nt) = 2d0*winvn% spin(ps(1:Nt)) + 1d0
!           mass(1:Nt) = winvn% Z(ps(1:Nt))*mp + winvn% N(ps(1:Nt))*mn - winvn% binding_energy(ps(1:Nt))*massfac
            mass(1:Nt) = winvn% W(ps(1:Nt))

            ! log(prefactor of reverse_ratio)
            tmp = product(mass(1:Ni))/product(mass(Ni+1:Nt))
            rates% inverse_coefficients(1,i) = pow(tmp,1.5d0)*(product(g(1:Ni))/product(g(Ni+1:Nt)))

            ! -Q/(kB*10**9)
            sum1 = sum(winvn% binding_energy(ps(1:Ni)))
            sum2 = sum(winvn% binding_energy(ps(Ni+1:Nt)))
            rates% inverse_coefficients(2,i) = (sum1-sum2)*conv/kB/1d9

            ! This should be 0 for non-photo-disintegration reverse rates and 1 for photos's in the reverse channel
            if (No==1) then
               rates% inverse_exp(i) = 1
               if(rates% inverse_coefficients(2,i)<0) then
                  ! negative values denote endothermic photodisintegrations
                  ! We want rate_photo/rate_forward
                  rates% inverse_coefficients(1,i) = rates% inverse_coefficients(1,i) * fac
               else
                  ! positive values denote exothermic photodisintegrations
                  ! We divide by fac and invert the T^3/2 as we want to compute
                  ! rate_reverse/rate_photo
                  rates% inverse_coefficients(1,i) = rates% inverse_coefficients(1,i) / fac
                  rates% inverse_exp(i) = -1 ! We us this term in a log() expression
               end if
            else
               rates% inverse_exp(i) = 0
            end if
            rates% inverse_coefficients(1,i) = log(rates% inverse_coefficients(1,i))


            rates% inverse_part(:,i) = product(winvn% pfcn(1:npart,ps(1:Ni)),dim=2)/ &
               & product(winvn% pfcn(1:npart,ps(Ni+1:Nt)),dim=2)
         end   do loop_over_rates
      end subroutine compute_rev_ratio

      subroutine do_parse_reaction_handle(handle, num_in, num_out, iso_ids, op, ierr)
         use chem_def
         use chem_lib
         character (len=*), intent(in) :: handle
         integer, intent(out) :: num_in, num_out
         integer, intent(out) :: iso_ids(:) ! holds chem_ids for input and output species
         character (len=*), intent(out) :: op ! e.g., 'pg', 'wk', 'to', or ...
         integer, intent(out) :: ierr
         
         integer :: len, i, j, cnt, cid, extra_in, extra_out
         logical :: doing_inputs
         
         num_in = 0; num_out = 0; op = ''
         ierr = -1
         len = len_trim(handle)
         if (handle(1:2) /= 'r_') return
         i = 3
         cnt = 0
         doing_inputs = .true.
         do while (i <= len)
            call nxt ! set j to last char of token
            cid = chem_get_iso_id(handle(i:j))
            if (cid == nuclide_not_found) then
               if (doing_inputs) then
                  op = handle(i:j)
                  extra_in = -1
                  extra_out = -1
                  if (j == i+1) then ! check 2 character ops
                     select case (op(1:1))
                        case ('p')
                           extra_in = ih1
                        case ('a')
                           extra_in = ihe4
                        case ('n')
                           extra_in = ineut
                        case ('g')
                           extra_in = 0
                        case default
                     end select
                     if (extra_in >= 0) then
                        if (extra_in > 0) then
                           cnt = cnt+1
                           if (cnt /= 2) then
                              !write(*,*) 'failed to parse ' // &
                              !   trim(handle) // ' -- problem with ' // handle(i:j)
                              return
                           end if
                           if (chem_isos% Z(iso_ids(1)) >= chem_isos% Z(extra_in)) then
                              iso_ids(2) = iso_ids(1)
                              iso_ids(1) = extra_in
                           else
                              iso_ids(2) = extra_in
                           end if
                        end if
                        select case (op(2:2))
                           case ('p')
                              extra_out = ih1
                           case ('a')
                              extra_out = ihe4
                           case ('n')
                              extra_out = ineut
                           case ('g')
                              extra_out = 0
                           case default
                              !write(*,*) 'failed to parse ' // &
                              !   trim(handle) // ' -- problem with ' // handle(i:j)
                              return
                        end select
                     end if               
                  end if
                  num_in = cnt
                  doing_inputs = .false.
                  if (extra_out > 0) then
                     cnt = cnt+1
                     iso_ids(cnt) = extra_out
                  end if
               else
                  !write(*,*) 'failed to parse ' // &
                  !   trim(handle) // ' -- problem with ' // handle(i:j)
                  return
               end if
            else
               cnt = cnt+1
               iso_ids(cnt) = cid
            end if
            i = j+2
         end do
         num_out = cnt - num_in
         ierr = 0
         
         contains
         
         subroutine nxt
            j = i
            do
               if (j >= len) return
               j = j+1
               if (handle(j:j) == '_') then
                  j = j-1; return
               end if
            end do
         end subroutine nxt
         

      end subroutine do_parse_reaction_handle

      subroutine reaction_handle(num_in, num_out, iso_ids, reaction_flag, handle)
         use chem_def, only: chem_isos
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: iso_ids(:)
         character (len=*), intent(in) :: reaction_flag
         character (len=*), intent(out) :: handle
         logical, parameter :: reverse = .false.
         call get1_reaction_handle(num_in, num_out, iso_ids, chem_isos, reverse, reaction_flag, handle)
      end subroutine reaction_handle
      
      subroutine reverse_reaction_handle(num_in, num_out, iso_ids, handle)
         use chem_def, only: chem_isos
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: iso_ids(:)
         character (len=*), intent(out) :: handle
         logical, parameter :: reverse = .true.
         character (len=1) :: reaction_flag = '-'
         call get1_reaction_handle(num_in, num_out, iso_ids, chem_isos, reverse, reaction_flag, handle)
      end subroutine reverse_reaction_handle         
      
      subroutine get_reaction_handle(num_in, num_out, pspecies, nuclides, reaction_flag, handle)
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: pspecies(:)
         type(nuclide_data), intent(in) :: nuclides
         character (len=*), intent(in) :: reaction_flag
         character (len=*), intent(out) :: handle
         logical, parameter :: reverse = .false.
         call get1_reaction_handle(num_in, num_out, pspecies, nuclides, reverse, reaction_flag, handle)
      end subroutine get_reaction_handle
      
      subroutine get_reverse_reaction_handle(num_in, num_out, pspecies, nuclides, handle)
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: pspecies(:)
         type(nuclide_data), intent(in) :: nuclides
         character (len=*), intent(out) :: handle
         logical, parameter :: reverse = .true.
         character (len=1) :: reaction_flag = '-'
         call get1_reaction_handle(num_in, num_out, pspecies, nuclides, reverse, reaction_flag, handle)
      end subroutine get_reverse_reaction_handle
      
      subroutine get1_reaction_handle( &
            num_in, num_out, pspecies_in, nuclides, reverse, reaction_flag, handle)
         use chem_def, only: ih1, ih2, ih3, ihe3, ihe4, ibe7, ili7, chem_isos
         integer, intent(in) :: num_in, num_out
         integer, intent(in) :: pspecies_in(:)
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: reverse
         character (len=*), intent(in) :: reaction_flag
         character (len=*), intent(out) :: handle

         integer :: in1, in2, out1, out2, num, pspecies(num_in + num_out)
         logical :: do_long_form, ec_flag, wk_flag
         
         include 'formats.dek'
         
         num = num_in + num_out
         pspecies(1:num) = pspecies_in(1:num)
         call sort(num_in, pspecies(1:num_in))
         call sort(num_out, pspecies(num_in+1:num))
         ec_flag = (reaction_flag == 'e')
         wk_flag = (reaction_flag == 'w')
         
         if (ec_flag) then ! special cases
            if (reverse) then
               handle = ''
               return
            end if
            if (num_in == 2 .and. num_out == 1) then
               if (nuclides% chem_id(pspecies(1)) == ih1 .and. &
                   nuclides% chem_id(pspecies(2)) == ih1 .and. &
                   nuclides% chem_id(pspecies(3)) == ih2) then
                  handle = 'r_h1_h1_ec_h2'
                  return
               end if
            else if (num_in == 1 .and. num_out == 1) then
               if (nuclides% chem_id(pspecies(1)) == ihe3 .and. &
                   nuclides% chem_id(pspecies(2)) == ih3) then
                  handle = 'r_he3_ec_h3'
                  return
               end if
               if (nuclides% chem_id(pspecies(1)) == ibe7 .and. &
                   nuclides% chem_id(pspecies(2)) == ili7) then
                  handle = 'r_be7_wk_li7'
                  return 
               end if
            end if
         else if (wk_flag) then
            if (reverse) then
               handle = ''
               return
            end if
            if (num_in == 2 .and. num_out == 1) then
               if (nuclides% chem_id(pspecies(1)) == ih1 .and. &
                   nuclides% chem_id(pspecies(2)) == ih1 .and. &
                   nuclides% chem_id(pspecies(3)) == ih2) then
                  handle = 'r_h1_h1_wk_h2'
                  return
               end if
            else if (num_in == 2 .and. num_out == 1) then
               if (nuclides% chem_id(pspecies(1)) == ih1 .and. &
                   nuclides% chem_id(pspecies(2)) == ihe3 .and. &
                   nuclides% chem_id(pspecies(3)) == ihe4) then
                  handle = 'r_h1_he3_wk_he4'
                  return
               end if
            end if
         end if

         in1 = 0; in2 = 0; out1 = 0; out2 = 0
         do_long_form = .true.
         if (num_in == 1 .and. num_out == 1) then
            call do_n_to_m(1,1)
            do_long_form = one_one()
         else if (num_in == 1 .and. num_out == 2) then
            call do_n_to_m(1,2)
            if (reverse) then
               do_long_form = two_one()
            else
               do_long_form = one_two()
            end if
         else if (num_in == 2 .and. num_out == 1) then
            call do_n_to_m(2,1)
            if (reverse) then
               do_long_form = one_two()
            else
               do_long_form = two_one()
            end if
         else if (num_in == 2 .and. num_out == 2) then
            call do_n_to_m(2,2)
            do_long_form = two_two()            
         end if

         if (do_long_form) then
            call long_form
         else if (out2 /= 0) then
            handle = trim(handle) // '_' // nuclides% name(out2)
         else
            handle = trim(handle) // '_' // nuclides% name(out1)
         end if

         
         contains
         
         subroutine sort(n, species)
            integer :: n
            integer :: species(n)
            integer :: i, j, Zi, Ni, Zj, Nj, isomer_j, isomer_i, cid
            include 'formats.dek'
            do i=1,n-1
               cid = species(i)
               if (cid <= 0) cycle
               Zi = chem_isos% Z(cid)
               Ni = chem_isos% N(cid)
               isomer_i = chem_isos% isomeric_state(cid)
               do j=i+1,n
                  cid = species(j)
                  if (cid <= 0) cycle
                  Zj = chem_isos% Z(cid)
                  Nj = chem_isos% N(cid)
                  isomer_j = chem_isos% isomeric_state(cid)
                  if (Zj > Zi) cycle
                  if (Zj == Zi) then
                     if (Nj > Ni) cycle
                     if (Nj == Ni) then
                        if (isomer_j >= isomer_i) cycle
                     end if
                  end if
                  ! exchange i and j
                  species(j) = species(i)
                  species(i) = cid
                  Zi = Zj
                  Ni = Nj
                  isomer_i = isomer_j
               end do
            end do
         end subroutine sort
         
         subroutine long_form
            integer :: i, cid
            character (len=3) :: op
            handle = 'r_'
            if (wk_flag) then
               op = 'wk_'
            else if (ec_flag) then
               op = 'ec_'
            else
               op = 'to_'
            end if
            if (reverse) then
               do i = num_in+1,num_in+num_out
                  cid = pspecies(i)
                  handle = trim(handle) // trim(nuclides% name(cid)) // '_'
               end do
               handle = trim(handle) // op
               do i = 1,num_in
                  cid = pspecies(i)
                  handle = trim(handle) // trim(nuclides% name(cid))
                  if (i < num_in) handle = trim(handle) // '_'
               end do
            else
               do i = 1,num_in
                  cid = pspecies(i)
                  handle = trim(handle) // trim(nuclides% name(cid)) // '_'
               end do
               handle = trim(handle) // op
               do i = num_in+1,num_in+num_out
                  cid = pspecies(i)
                  handle = trim(handle) // trim(nuclides% name(cid))
                  if (i < num_in+num_out) handle = trim(handle) // '_'
               end do
            end if
         end subroutine long_form
         
         logical function one_one()
            one_one = .true.
            if (in1 == 0 .or. out1 == 0) return
            one_one = .false.
            if (nuclides% Z(out1) == nuclides% Z(in1) - 1 .and. &
                nuclides% N(out1) == nuclides% N(in1) + 1) then
               handle = trim(handle) // '_wk'
            else if (nuclides% Z(out1) == nuclides% Z(in1) + 1 .and. &
                     nuclides% N(out1) == nuclides% N(in1) - 1) then
               handle = trim(handle) // '_wk-minus'
            else
               one_one = .true.
            end if
         end function one_one
         
         logical function one_two()
            one_two = .true.
            if (in1 == 0 .or. out1 == 0 .or. out2 == 0 .or. out1 == out2) return
            one_two = .false.
            if (nuclides% Z(out1) == 0 .and. nuclides% N(out1) == 1 .and. &
                     nuclides% Z(out2) == nuclides% Z(in1) .and. &
                     nuclides% N(out2) == nuclides% N(in1) - 1) then
               handle = trim(handle) // '_gn'
            else if (nuclides% Z(out1) == 1 .and. nuclides% N(out1) == 0 .and. &
                     nuclides% Z(out2) == nuclides% Z(in1) - 1 .and. &
                     nuclides% N(out2) == nuclides% N(in1)) then
               handle = trim(handle) // '_gp'
            else if (nuclides% Z(out1) == 2 .and. nuclides% N(out1) == 2 .and. &
                     nuclides% Z(out2) == nuclides% Z(in1) - 2 .and. &
                     nuclides% N(out2) == nuclides% N(in1) - 2) then
               handle = trim(handle) // '_ga'
            else if (nuclides% Z(out1) == 1 .and. nuclides% N(out1) == 0 .and. &
                     nuclides% Z(out2) == nuclides% Z(in1) - 2 .and. &
                     nuclides% N(out2) == nuclides% N(in1) + 1) then
               handle = trim(handle) // '_wk_h1'
            else if (nuclides% Z(out1) == 2 .and. nuclides% N(out1) == 2 .and. &
                     nuclides% Z(out2) == nuclides% Z(in1) - 3 .and. &
                     nuclides% N(out2) == nuclides% N(in1) - 1) then
               handle = trim(handle) // '_wk_he4'
            else
               one_two = .true.
            end if
         end function one_two
         
         logical function two_one()
            include 'formats.dek'
            two_one = .true.
            if (in1 == 0 .or. in2 == 0 .or. out1 == 0 .or. in1 == in2) return
            two_one = .false.
            if (nuclides% Z(in1) == 0 .and. nuclides% N(in1) == 1 .and. &
                     nuclides% Z(out1) == nuclides% Z(in2) .and. &
                     nuclides% N(out1) == nuclides% N(in2) + 1) then
               handle = trim(handle) // '_ng'
            else if (nuclides% Z(in1) == 1 .and. nuclides% N(in1) == 0 .and. &
                     nuclides% Z(out1) == nuclides% Z(in2) + 1 .and. &
                     nuclides% N(out1) == nuclides% N(in2)) then
               handle = trim(handle) // '_pg'
            else if (nuclides% Z(in1) == 2 .and. nuclides% N(in1) == 2 .and. &
                     nuclides% Z(out1) == nuclides% Z(in2) + 2 .and. &
                     nuclides% N(out1) == nuclides% N(in2) + 2) then
               handle = trim(handle) // '_ag'
            else
               two_one = .true.
            end if
         end function two_one
         
         logical function two_two()
            two_two = .true.
            if (in1 == 0 .or. in2 == 0 .or. out1 == 0 .or. out2 == 0 .or. &
                in1 == in2 .or. out1 == out2) return
            two_two = .false.
            if (nuclides% Z(in1) == 2 .and. nuclides% N(in1) == 2 .and. &
                     nuclides% Z(out1) == 1 .and. nuclides% N(out1) == 0 .and. &
                     nuclides% Z(out2) == nuclides% Z(in2) + 1 .and. &
                     nuclides% N(out2) == nuclides% N(in2) + 2) then
               handle = trim(handle) // '_ap'
            else if (nuclides% Z(in1) == 1 .and. nuclides% N(in1) == 0 .and. &
                     nuclides% Z(out1) == 2 .and. nuclides% N(out1) == 2 .and. &
                     nuclides% Z(out2) == nuclides% Z(in2) - 1 .and. &
                     nuclides% N(out2) == nuclides% N(in2) - 2) then
                handle = trim(handle) // '_pa'
            else if (nuclides% Z(in1) == 2 .and. nuclides% N(in1) == 2 .and. &
                     nuclides% Z(out1) == 0 .and. nuclides% N(out1) == 1 .and. &
                     nuclides% Z(out2) == nuclides% Z(in2) + 2 .and. &
                     nuclides% N(out2) == nuclides% N(in2) + 1) then
                handle = trim(handle) // '_an'
            else if (nuclides% Z(in1) == 0 .and. nuclides% N(in1) == 1 .and. &
                     nuclides% Z(out1) == 2 .and. nuclides% N(out1) == 2 .and. &
                     nuclides% Z(out2) == nuclides% Z(in2) - 2 .and. &
                     nuclides% N(out2) == nuclides% N(in2) - 1) then
                handle = trim(handle) // '_na'
            else if (nuclides% Z(in1) == 1 .and. nuclides% N(in1) == 0 .and. &
                     nuclides% Z(out1) == 0 .and. nuclides% N(out1) == 1 .and. &
                     nuclides% Z(out2) == nuclides% Z(in2) + 1 .and. &
                     nuclides% N(out2) == nuclides% N(in2) - 1) then
                handle = trim(handle) // '_pn'
            else if (nuclides% Z(in1) == 0 .and. nuclides% N(in1) == 1 .and. &
                     nuclides% Z(out1) == 1 .and. nuclides% N(out1) == 0 .and. &
                     nuclides% Z(out2) == nuclides% Z(in2) - 1 .and. &
                     nuclides% N(out2) == nuclides% N(in2) + 1) then
                handle = trim(handle) // '_np'
            else
               two_two = .true.
            end if
         end function two_two
         
         subroutine do_n_to_m(n,m)
            integer, intent(in) :: n, m ! each is either 1 or 2
            integer :: j
            in1 = 0; in2 = 0; out1 = 0; out2 = 0
            if (.not. reverse) then
               in1 = pspecies(1)
               if (n == 2) then
                  in2 = pspecies(2)
                  if (in2 == 0) then
                     in1 = 0; return
                  end if
                  call switch_if_necessary(in1,in2)
                  handle = 'r_' // nuclides% name(in2)
               else if (n == 1) then
                  handle = 'r_' // nuclides% name(in1)
               else
                  in1 = 0
                  return
               end if
               out1 = pspecies(n+1)
               if (m == 2) then
                  out2 = pspecies(n+2)
                  call switch_if_necessary(out1,out2)
               end if
            else
               in1 = pspecies(n+1)
               if (m == 2) then
                  in2 = pspecies(n+2)
                  if (in2 == 0) then
                     in1 = 0; return
                  end if
                  call switch_if_necessary(in1,in2)
                  handle = 'r_' // nuclides% name(in2)
               else if (m == 1) then
                  handle = 'r_' // nuclides% name(in1)
               else
                  in1 = 0
                  return
               end if
               out1 = pspecies(1)
               if (n == 2) then
                  out2 = pspecies(2)
                  call switch_if_necessary(out1,out2)
               end if
            end if
         end subroutine do_n_to_m
         
         subroutine switch_if_necessary(iso1,iso2)
            integer, intent(inout) :: iso1, iso2
            integer :: j
            if (nuclides% Z(iso2) == 1 .and. nuclides% N(iso2) == 0) then ! iso2 is ih1
               j = iso1; iso1 = iso2; iso2 = j; return
            end if
            if (nuclides% Z(iso2) == 2 .and. nuclides% N(iso2) == 2) then ! iso2 is ihe4
               if (nuclides% Z(iso1) == 1 .and. nuclides% N(iso1) == 0) return ! iso1 is ih1
               j = iso1; iso1 = iso2; iso2 = j; return
            end if
         end subroutine switch_if_necessary
         
         
      end subroutine get1_reaction_handle
      

      end module reaclib_support
