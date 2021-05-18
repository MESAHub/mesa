! ***********************************************************************
!
!   Copyright (C) 2014-2019  The MESA Team
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


      module adjust_net

      use star_private_def

      implicit none

      contains


      subroutine check_adjust_net(s, species, &
            min_x_for_keep, min_x_for_n, min_x_for_add, &
            max_Z_for_add, max_N_for_add, max_A_for_add, ierr)
         use adjust_xyz, only: change_net
         use utils_lib, only: mkdir
         use chem_def, only: &
            chem_isos, ineut, max_el_z, el_name, element_min_N, element_max_N
         use chem_lib, only: chem_get_iso_id
         use star_utils, only: get_string_for_model_number
         use net_lib, only: show_net_reactions

         type (star_info), pointer :: s
         integer, intent(in) :: species
         real(dp), intent(in) :: &
            min_x_for_keep, min_x_for_n, min_x_for_add, max_Z_for_add, max_N_for_add, max_A_for_add
         integer, intent(out) :: ierr

         character (len=strlen) :: &
            net_name, fname, temp_fname, z_plus_n_str, cname, line_buf, species_str
         integer :: num_digits, i, j, io, io_new, nz, cid, Z, N, A, &
            new_species, next_model_number
         integer, dimension(0:max_el_z) :: min_N, max_N
         real(dp) :: max_x_for_species(species), max_x
         logical :: have_new_isos, still_in_net(species)
         logical, parameter :: adjust_abundances_for_new_isos = .false.

         include 'formats'

         ierr = 0
         nz = s% nz
         min_N = -1
         max_N = -1

         next_model_number = s% model_number + 1

         call include_iso(0,1)   ! neut
         call include_iso(1,0)   ! h1
         call include_iso(2,2)   ! he4
         call include_iso(6,6)   ! c12
         call include_iso(6,7)   ! c13
         call include_iso(7,7)   ! n14
         call include_iso(8,8)   ! o16
         call include_iso(10,10) ! ne20
         call include_iso(10,12) ! ne22
         call include_iso(12,12) ! mg24
         call include_iso(12,14) ! mg26
         call include_iso(14,14) ! si28
         call include_iso(16,16) ! s32
         call include_iso(26,30) ! fe56

         do j=1,species
            max_x = maxval(s% xa(j,1:nz))
            max_x_for_species(j) = max_x
            if (max_x < min_x_for_keep) cycle
            cid = s% chem_id(j)
            Z = chem_isos% Z(cid)
            N = chem_isos% N(cid)
            A=Z+N
            call include_iso(Z,N)
            if (max_x >= min_x_for_n .and. N+1 <= max_N_for_add .and. A+1 <= max_A_for_add) then
               call include_iso(Z,N+1)   ! (n,g)
               call include_iso(Z,N-1)   ! (g,n)
            end if
            if (max_x >= min_x_for_add) then
               if (Z+1 <= max_Z_for_add.and. A+1 <= max_A_for_add) then
                  call include_iso(Z+1,N)   ! (p,g)
                  call include_iso(Z-1,N)   ! (g,p)
               end if
               if (Z+2 <= max_Z_for_add .and. N+2 <= max_N_for_add .and. A+4 <= max_A_for_add) then
                  call include_iso(Z+2,N+2) ! (a,g)
                  call include_iso(Z-2,N-2) ! (g,a)
               end if
               if (Z+2 <= max_Z_for_add .and. N+1 <= max_N_for_add .and. A+3 <= max_A_for_add) then
                  call include_iso(Z+2,N+1) ! (a,n)
                  call include_iso(Z-2,N-1) ! (n,a)
               end if
               if (Z+1 <= max_Z_for_add .and. N+2 <= max_N_for_add .and. A+3 <= max_A_for_add) then
                  call include_iso(Z+1,N+2) ! (a,p)
                  call include_iso(Z-1,N-2) ! (p,a)
               end if
               if (Z+1 <= max_Z_for_add .and. N+1 <= max_N_for_add .and. A+2 <= max_A_for_add) then
                  call include_iso(Z+1,N-1) ! (p,n)
                  call include_iso(Z-1,N+1) ! (n,p)
               end if
               if (Z+4 <= max_Z_for_add .and. N+4 <= max_N_for_add .and. A+8 <= max_A_for_add) then
                  call include_iso(Z+4,N+4) ! (2a,g) ! extend alpha chain by 2
                  call include_iso(Z+3,N+4) ! (2a,p) ! extend alpha chain by 2
               end if
            end if
         end do

         temp_fname = '.temp_net'
         open(newunit=io, file=trim(temp_fname), action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(temp_fname)
            return
         end if

         write(io,'(a)') 'add_isos_and_reactions('
         write(io,'(a)') '!       iso    log10 max x'

         new_species = 0
         have_new_isos = .false.
         still_in_net(:) = .false.
         do Z=0,max_el_z
            if (min_N(Z) < 0 .or. max_N(Z) < min_N(Z)) cycle
            if (Z == 0) then
               i = s% net_iso(ineut)
               if (i == 0) then
                  write(*,*) '  add ' // trim(el_name(Z))
                  have_new_isos = .true.
                  write(io,'(3x,a8,a)') trim(el_name(Z)), '   ! newly added'
               else
                  still_in_net(i) = .true.
                  write(io,'(3x,a8,a,f10.5)') trim(el_name(Z)), '   ! ', &
                     log10(max(1d-199,max_x_for_species(i)))
               end if
               new_species = new_species + 1
               cycle
            end if
            write(io,*)
            do N = min_N(Z), max_N(Z)
               write(z_plus_n_str,'(i4)') Z+N
               write(cname,'(a)') trim(el_name(Z)) // trim(adjustl(z_plus_n_str))
               cid = chem_get_iso_id(cname)
               if (cid <= 0) cycle
               i = s% net_iso(cid)
               if (i == 0) then
                  write(*,*) '  add ' // trim(cname)
                  have_new_isos = .true.
                  write(io,'(3x,a8,a)') trim(cname), '   ! newly added'
               else
                  still_in_net(i) = .true.
                  write(io,'(3x,a8,a,f10.5)') trim(cname), '   ! ', &
                     log10(max(1d-199,max_x_for_species(i)))
               end if
               new_species = new_species + 1
            end do
         end do

         do i=1,species
            if (still_in_net(i)) cycle
            write(*,*) ' drop ' // trim(chem_isos% name(s% chem_id(i)))
         end do

         write(io,'(a)') ')'
         write(io,*)

         close(io)

         if (new_species == species .and. .not. have_new_isos) return

         open(newunit=io, file=trim(temp_fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(temp_fname)
            return
         end if

         num_digits = 5
         call get_string_for_model_number( &
               '', next_model_number, num_digits, net_name)
         write(species_str,'(i4)') new_species
         net_name = trim(net_name) // '_' // trim(adjustl(species_str)) // '.net'
         fname = 'nets/' // trim(net_name)
         call mkdir('nets/')
         open(newunit=io_new, file=trim(fname), action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            return
         end if

         write(io_new,'(a,i8,a)') '! ' //  trim(net_name), new_species, ' species'
         write(io_new,*)

         do i=1,10000
            read(io, fmt='(a)', iostat=ierr) line_buf
            if (ierr /= 0) exit
            write(io_new, fmt='(a)') trim(line_buf)
         end do
         ierr = 0

         close(io)
         close(io_new)

         write(*,'(i11,a,i8,a)') next_model_number, &
            '   change to ' //  trim(fname), new_species, ' species'

         call change_net( &
            s% id, adjust_abundances_for_new_isos, net_name, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in change_net ' // trim(net_name)
            stop 'check_adjust_net'
            return
         end if

         if (net_name /= s% net_name) then
            write(*,*) '   new net_name ', trim(net_name)
            write(*,*) 'old s% net_name ', trim(s% net_name)
            write(*,*) 'failed to change'
            stop 'check_adjust_net'
         end if

         s% using_revised_net_name = .true.
         s% revised_net_name = s% net_name
         s% need_to_setvars = .true.


         contains


         subroutine include_iso(Z,N)
            integer, intent(in) :: Z,N
            if (Z < 0 .or. Z > max_el_z .or. N < 0) return
            if (N < element_min_N(Z) .or. N > element_max_N(Z)) return
            if (N < min_N(Z) .or. min_N(Z) < 0) min_N(Z) = N
            if (N > max_N(Z)) max_N(Z) = N
         end subroutine include_iso


      end subroutine check_adjust_net



      end module adjust_net












