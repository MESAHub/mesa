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
 
      module reaclib_print
      use rates_def
      
      implicit none

      contains


      subroutine write_reaction_data(unitno, rates, ierr)
         integer, intent(in) :: unitno
         type(reaction_data), intent(in) :: rates
         integer, intent(out) :: ierr
         integer :: i
         character(len=*), parameter :: fmt = &
         '(a36, i3, tr1, 6(i5, tr1), a1, tr1, f5.1, tr1, /, 4es14.6, /, 3es14.6, /, 2es14.6, tr1, i2, tr1, 3(/, 8(f12.9, tr1)))'
         write(unitno, *) rates% nreactions
         write(unitno, *) rates% nchapters_included
         write(unitno, '(i6)') rates% chapters_present(:)
         write(unitno, '(2i6, tr1)') (rates% bookmarks(:, i), i=1, nchapters)
         do i = 1, rates% nreactions
            write(unitno, fmt=fmt, iostat=ierr) &
               rates% reaction_handle(i), rates% chapter(i), rates% pspecies(:, i), rates% reaction_flag(i), &
               1.0d0/rates% weight(i), rates% coefficients(:, i), rates% inverse_coefficients(:, i), rates% inverse_exp(i), &
               rates% inverse_part(:, i)
         end do
      end subroutine write_reaction_data
      

      subroutine pretty_print_reactions(unitno, rates, nuclides, ierr)
         integer, intent(in) :: unitno
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         integer, intent(out) :: ierr
         integer :: i, pass
         logical :: reverse
         character (len=100) :: str
         write(unitno,'(/,a,/)') 'forward reactions'
         do pass = 1, 2
            str = ''
            reverse = (pass == 2)
            do i = 1, rates% nreactions
               call do_pretty_print_reaction(unitno, i, rates, nuclides, reverse, str, ierr)
               if (ierr /= 0) exit
            end do
            if (pass == 1) write(unitno,'(/,a,/)') 'reverse reactions'
         end do
      end subroutine pretty_print_reactions
      
      
      subroutine do_pretty_print_reaction(unitno, i, rates, nuclides, reverse, str, ierr)
         integer, intent(in) :: unitno, i
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: reverse
         character (len=100), intent(inout) :: str
         integer, intent(out) :: ierr
         character (len=100) :: str_nxt
         ierr = 0
         str_nxt = ''
         if (reverse .and. i <= rates% num_from_weaklib) return
         select case(rates% chapter(i))
            case(r_one_one)
               call write_n_to_m(1,1)
            case(r_one_two)
               call write_n_to_m(1,2)
            case(r_one_three)
               call write_n_to_m(1,3)
            case(r_two_one)
               call write_n_to_m(2,1)
            case(r_two_two)
               call write_n_to_m(2,2)
            case(r_two_three)
               call write_n_to_m(2,3)
            case(r_two_four)
               call write_n_to_m(2,4)
            case(r_three_one)
               call write_n_to_m(3,1)
            case(r_three_two)
               call write_n_to_m(3,2)
            case(r_four_two)
               call write_n_to_m(4,2)
            case(r_one_four)
               call write_n_to_m(1,4)
         end select
         if (trim(str_nxt) /= trim(str) .and. len_trim(str_nxt) > 0) then
            if (i <= rates% num_from_weaklib) str_nxt = trim(str_nxt) // '    weaklib'
            write(unitno,fmt='(a)') trim(str_nxt)
            str = str_nxt
         end if
         
         contains
         
         subroutine write_n_to_m(n,m)
            integer, intent(in) :: n, m
            integer :: j
            if (.not. reverse) then
               do j=1,n
                  str_nxt = trim(str_nxt) // ' ' // nuclides% name(rates% pspecies(j, i))
               end do
               str_nxt = trim(str_nxt) //  " -> "
               do j=1,m
                  str_nxt = trim(str_nxt) // ' ' // nuclides% name(rates% pspecies(n+j, i))
               end do
            else if (rates% reaction_flag(i) /= 'w' .and. rates% reaction_flag(i) /= 'e') then
               do j=1,m
                  str_nxt = trim(str_nxt) // ' ' // nuclides% name(rates% pspecies(n+j, i))
               end do
               str_nxt = trim(str_nxt) // " -> "
               do j=1,n
                  str_nxt = trim(str_nxt) // ' ' // nuclides% name(rates% pspecies(j, i))
               end do
            end if
         end subroutine write_n_to_m
         
         
      end subroutine do_pretty_print_reaction


      subroutine print_short_format_reactions(unitno, rates, nuclides, ierr)
         integer, intent(in) :: unitno
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         integer, intent(out) :: ierr
         integer :: i, pass
         logical :: reverse
         character (len=100) :: str
         write(unitno,'(/,a,/)') 'forward reactions'
         do pass = 1, 2
            str = ''
            reverse = (pass == 2)
            do i = 1, rates% nreactions
               call do_print_short_format_reaction(unitno, i, rates, nuclides, reverse, str, ierr)
               if (ierr /= 0) exit
            end do
            if (pass == 1) write(unitno,'(/,a,/)') 'reverse reactions'
         end do
      end subroutine print_short_format_reactions
      
      
      subroutine do_print_short_format_reaction(unitno, i, rates, nuclides, reverse, str, ierr)
         use reaclib_support, only: get1_reaction_handle
         integer, intent(in) :: unitno, i
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: reverse
         character (len=100), intent(inout) :: str
         integer, intent(out) :: ierr
         
         character (len=100) :: str_nxt
         integer :: chapter, num_in, num_out
         
         ierr = 0
         str = ''
         
         if (reverse .and. &
            (rates% reaction_flag(i) == 'w' .or. rates% reaction_flag(i) == 'e')) return
         
         chapter = rates% chapter(i)
         num_in = Nin(chapter)
         num_out = Nout(chapter)
         
         call get1_reaction_handle( &
            num_in, num_out, rates% pspecies(:,i), nuclides, reverse, &
            rates% reaction_flag(i), str_nxt)         
         if (trim(str_nxt) /= trim(str) .and. len_trim(str_nxt) > 0) then
            write(unitno,fmt='(a)') trim(str_nxt)
            str = str_nxt
         end if
         
      end subroutine do_print_short_format_reaction


      end module reaclib_print
