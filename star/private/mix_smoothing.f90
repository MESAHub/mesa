! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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


      module mix_smoothing

      use const_def
      use num_lib
      use utils_lib
      use star_private_def
      use star_utils, only: find_cell_for_mass

      implicit none

      private
      public :: set_newly_non_conv


      contains


      subroutine set_newly_non_conv(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: num, num_old, j, nz
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         nz = s% nz
         s% newly_nonconvective(1:nz) = .false.
         num = s% n_conv_regions
         num_old = s% n_conv_regions_old
         if (num == 0 .or. num_old == 0) return
         if (dbg) then
            write(*,*)
            do j=1,s% n_conv_regions_old
               write(*,2) 'old conv region', j, s% cz_bot_mass_old(j)/Msun, s% cz_top_mass_old(j)/Msun
            end do
            write(*,*)
            do j=1,s% n_conv_regions
               write(*,2) 'conv region', j, s% cz_bot_mass(j)/Msun, s% cz_top_mass(j)/Msun
            end do
            write(*,*)
         end if
         if (dbg) write(*,*) 'set_newly_non_conv: call do_all_regions'
         call do_all_regions( &
            s, set_top_moved_down, set_bottom_moved_up, set_for_departed_region, &
            num_old, s% cz_top_mass_old, s% cz_bot_mass_old, &
            num, s% cz_top_mass, s% cz_bot_mass, ierr)
      end subroutine set_newly_non_conv


      subroutine set_top_moved_down(s, nz, species, top_old, top_new)
         type (star_info), pointer :: s
         integer, intent(in) :: nz, species
         real(dp), intent(in) :: top_old, top_new
         integer :: ktop_old, ktop_new
         ktop_old = find_cell_for_mass(s,top_old)
         ktop_new = find_cell_for_mass(s,top_new)
         s% newly_nonconvective(ktop_old:ktop_new) = .true.
      end subroutine set_top_moved_down


      subroutine set_bottom_moved_up(s, nz, species, bot_old, bot_new)
         type (star_info), pointer :: s
         integer, intent(in) :: nz, species
         real(dp), intent(in) :: bot_old, bot_new
         integer :: kbot_old, kbot_new
         kbot_old = find_cell_for_mass(s,bot_old)
         kbot_new = find_cell_for_mass(s,bot_new)
         s% newly_nonconvective(kbot_new:kbot_old) = .true.
      end subroutine set_bottom_moved_up


      subroutine set_for_departed_region(s, nz, species, top, bot)
         type (star_info), pointer :: s
         integer, intent(in) :: nz, species
         real(dp), intent(in) :: top, bot
         integer :: ktop, kbot
         ktop = find_cell_for_mass(s,top)
         kbot = find_cell_for_mass(s,bot)
         s% newly_nonconvective(ktop:kbot) = .true.
      end subroutine set_for_departed_region


      subroutine do_all_regions( &
            s, top_moved_down, bottom_moved_up, departed_region, &
            num_old, cz_top_mass_old, cz_bot_mass_old, &
            num, cz_top_mass, cz_bot_mass, ierr)
         type (star_info), pointer :: s
         interface
            subroutine top_moved_down(s, nz, species, top, bot)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: nz, species
               real(dp), intent(in) :: top, bot
            end subroutine top_moved_down
            subroutine bottom_moved_up(s, nz, species, top, bot)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: nz, species
               real(dp), intent(in) :: top, bot
            end subroutine bottom_moved_up
            subroutine departed_region(s, nz, species, top, bot)
               use star_def, only: star_info
               use const_def, only: dp
               type (star_info), pointer :: s
               integer, intent(in) :: nz, species
               real(dp), intent(in) :: top, bot
            end subroutine departed_region
         end interface
         integer, intent(in) :: num, num_old
         real(dp), dimension(:) :: cz_top_mass_old, cz_bot_mass_old, &
            cz_top_mass, cz_bot_mass
         integer, intent(out) :: ierr

         real(dp), dimension(max_num_mixing_regions) :: region_dm, region_dm_old
         integer :: nz, species, i, j, k, j_top, j_bot
         real(dp) :: top, bot, top_old, bot_old
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         nz = s% nz
         species = s% species

         if (dbg) write(*,*)
         if (dbg) write(*,*) 'do_all_regions'

         do i=1,num
            region_dm(i) = cz_top_mass(i) - cz_bot_mass(i)
            if (dbg) write(*,2) 'cz_top_mass cz_bot_mass', i, &
               cz_top_mass(i)/Msun, cz_bot_mass(i)/Msun
         end do

         do i=1,num_old
            region_dm_old(i) = cz_top_mass_old(i) - cz_bot_mass_old(i)
            if (dbg) write(*,2) 'cz_top_mass_old cz_bot_mass_old', i, &
               cz_top_mass_old(i)/Msun, cz_bot_mass_old(i)/Msun
         end do

         do i=1,num
            j = maxloc(region_dm(1:num),dim=1)
            region_dm(j) = -1 ! mark as done
            top = cz_top_mass(j)
            bot = cz_bot_mass(j)
            j_top = minloc(abs(top - cz_top_mass_old(1:num_old)),dim=1)
            top_old = cz_top_mass_old(j_top)
            j_bot = minloc(abs(bot - cz_bot_mass_old(1:num_old)),dim=1)
            bot_old = cz_bot_mass_old(j_bot)
            if (min(top,top_old) > max(bot,bot_old)) then ! overlap
               if (top < top_old) then
                  if (dbg) write(*,2) 'top_moved_down: top_old, top', i, &
                     top_old/Msun, top/Msun
                  call top_moved_down(s, nz, species, top_old, top) ! top moved down
               end if
               if (bot > bot_old) then
                  if (dbg) write(*,2) 'bottom_moved_up: bot_old, bot', i, &
                     bot_old/Msun, bot/Msun
                  call bottom_moved_up(s, nz, species, bot_old, bot) ! bottom moved up
               end if
            end if
            do j=1,num_old
               if (cz_top_mass_old(j) <= top_old .and. cz_bot_mass_old(j) >= bot_old) &
                  region_dm_old(j) = -1 ! mark as used
            end do
         end do

         do j=1,num_old
            if (region_dm_old(j) < 0) cycle ! was used
            ! wasn't used, so region is no longer convective
            call departed_region(s, nz, species, cz_top_mass_old(j), cz_bot_mass_old(j))
            if (dbg) write(*,1) 'departed_region top bot', &
               cz_top_mass_old(j)/Msun, cz_bot_mass_old(j)/Msun
         end do

      end subroutine do_all_regions


      end module mix_smoothing

