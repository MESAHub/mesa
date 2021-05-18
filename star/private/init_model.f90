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

      module init_model

      use star_private_def
      use const_def

      implicit none


      integer :: min_when_created = 20081212

      contains


      subroutine get_revised_mass(s, fullname, ierr)
         use read_model, only: do_read_saved_model
         use relax, only:do_relax_mass
         use relax, only: do_relax_Z
         type (star_info), pointer :: s
         character (len=*), intent(in) :: fullname
         integer, intent(out) :: ierr

         real(dp), parameter :: lg_max_abs_mdot = -3.5d0
         real(dp) :: init_mass, init_z
         logical :: want_RSP_model, is_rsp_model, want_RSP2_model, is_rsp2_model

         init_mass = s% initial_mass
         init_z = s% initial_z

         want_RSP_model = .false.
         want_RSP2_model = .false.
         call do_read_saved_model(s, fullname, &
            want_RSP_model, is_rsp_model, want_RSP2_model, is_rsp2_model, ierr)
         if (ierr /= 0) return

         if (abs(s% initial_z - init_z) > 1d-3*init_z .or. is_rsp_model .or. is_rsp2_model) then
            ierr = -1
            return
         end if

         call do_relax_mass(s% id, init_mass, lg_max_abs_mdot, ierr)
         if (ierr /= 0) return

         s% initial_mass = init_mass
         s% dt_next = min(s% dt_next, 1d2*secyer)

      end subroutine get_revised_mass


      subroutine get_zams_model(s, zams_filename, ierr)
         use alloc, only: allocate_star_info_arrays
         use utils_lib, only: is_bad
         use relax, only: do_relax_mass
         type (star_info), pointer :: s
         character (len=*) :: zams_filename
         integer, intent(out) :: ierr

         integer :: nz
         real(dp), dimension(:,:), pointer :: xh, xa
         real(dp), dimension(:), pointer :: q, dq, omega
         real(dp) :: init_mass
         logical :: in_range
         real(dp), parameter :: lg_max_abs_mdot = -3.5d0

         include 'formats'

         ierr = 0
         if (is_bad(s% initial_mass)) then
            write(*,1) 's% initial_mass', s% initial_mass
            stop 'get_zams_model'
         end if

         init_mass = s% initial_mass
         s% mstar = s% initial_mass*Msun
         s% xmstar = s% mstar
         s% M_center = 0

         call get1_zams_model( &
            s, zams_filename, nz, xh, xa, q, dq, &
            omega, in_range, ierr)
         if (ierr /= 0) then
            write(*,1) 'failed in get1_zams_model'
            stop 'get_zams_model'
         end if


         s% nz = nz
         call allocate_star_info_arrays(s, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if

         ! copy, then deallocate
         s% xh(:,1:nz) = xh(:,1:nz)
         s% xa(:,1:nz) = xa(:,1:nz)
         s% q(1:nz) = q(1:nz)
         s% dq(1:nz) = dq(1:nz)
         s% omega(1:nz) = omega(1:nz)

         call dealloc

         if (.not. in_range) then ! have revised s% initial_mass
            s% mstar = s% initial_mass*Msun
            s% xmstar = s% mstar
            s% M_center = 0
            s% dt_next = 1d2*secyer
            call do_relax_mass(s% id, init_mass, lg_max_abs_mdot, ierr)
            if (ierr /= 0) return
            s% initial_mass = init_mass
            s% dt_next = min(s% dt_next, 1d2*secyer)
         end if

         contains

         subroutine dealloc
            deallocate(xh, xa, q, dq, omega)
         end subroutine dealloc

      end subroutine get_zams_model


      subroutine get1_zams_model( &
            s, zams_filename, nz, xh, xa, q, dq, &
            omega, in_range, ierr)
         use utils_lib
         use const_def, only: mesa_data_dir
         use net, only: set_net
         type (star_info), pointer :: s
         character (len=*), intent(in) :: zams_filename
         integer, intent(out) :: nz
         real(dp), dimension(:,:), pointer :: xh, xa
         real(dp), dimension(:), pointer :: q, dq, omega, j_rot
         logical, intent(out) :: in_range
         integer, intent(out) :: ierr

         integer :: iounit, nz1, nz2, file_type, nvar_hydro, species
         character (len=250) :: fname, line
         real(dp) :: m1, m2, initial_mass
         logical :: okay

         include 'formats'

         ierr = 0
         nz = 0

         fname = zams_filename
         open(newunit=iounit, file=trim(fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            ierr = 0
            fname = trim(mesa_data_dir) // '/star_data/zams_models/' // trim(zams_filename)
            open(newunit=iounit, file=trim(fname), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'failed to open ' // trim(zams_filename)
               write(*, *) 'failed to open ' // trim(fname)
               return
            end if
         end if

         read(iounit, *, iostat=ierr) file_type
         if (ierr /= 0) then
            write(*, *) 'ERROR: first line needs to contains file type number: ' // trim(fname)
            close(iounit)
            return
         end if

         initial_mass = s% initial_mass
         nvar_hydro = s% nvar_hydro

         call read_zams_header ! sets net_name
         if (ierr /= 0) then
            close(iounit)
            write(*,*) 'failed in read_zams_header'
            stop 'get1_zams_model'
            return
         end if

         call set_net(s, s% net_name, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in set_net'
            return
         end if

         species = s% species

         ! find mass that is closest to init_m
         m1 = 0
         nz1 = 0
         okay = .false.
         in_range = .true.

         read(iounit, *, iostat=ierr) ! this is the 'M/Msun n_shells' header line
         if (ierr == 0) then
            index_loop: do
               read(iounit, *, iostat=ierr) m2, nz2
               if (ierr /= 0) exit index_loop
               if (m2 <= 0) then ! end of list
                  read(iounit, *, iostat=ierr) ! blank line
                  m2 = m1
                  nz2 = nz1
                  in_range = .false.
                  s% initial_mass = m2
                  initial_mass = m2
                  exit index_loop
               end if
               if (m2 >= initial_mass) then
                  if (m1 == 0) then
                     m1 = m2
                     nz1 = nz2
                     in_range = .false.
                     s% initial_mass = m2
                     initial_mass = m2
                  end if
                  ! skip to end of index
                  skip_loop: do
                     read(iounit, fmt='(a)', iostat=ierr) line
                     if (len_trim(line) == 0) then ! blank line indicates end of index
                        okay = .true.
                        exit index_loop
                     end if
                     if (ierr /= 0) exit index_loop
                  end do skip_loop
                  exit index_loop
               end if
               m1 = m2
               nz1 = nz2
            end do index_loop
         end if

         nz = nz1

         allocate(xh(nvar_hydro,nz), xa(species,nz), q(nz), dq(nz), &
            omega(nz), j_rot(nz), stat=ierr)
         if (ierr /= 0) then
            close(iounit)
            return
         end if

         call get1_mass( &
               s, iounit, m1, nz1, m2, nz2, initial_mass, &
               nvar_hydro, species, xh, xa, q, dq, &
               omega, j_rot, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get1_mass'
            stop 'get_zams_model'
         end if
         close(iounit)
         deallocate(j_rot)

         contains

         subroutine read_zams_header
            use read_model, only: read_properties
            integer :: year_month_day_when_created, iprop
            real(dp) :: dprop, initial_z, initial_y
            character (len=net_name_len) :: net_name
            read(iounit, *, iostat=ierr) ! skip blank line before property list
            include 'formats'
            if (ierr /= 0) return
            year_month_day_when_created = -1


            call read_properties(iounit, &
               s% net_name, iprop, iprop, year_month_day_when_created, &
               dprop, initial_z, initial_y, &
               dprop, iprop, dprop, dprop, &
               dprop, dprop, dprop, dprop, &
               dprop, dprop, dprop, dprop, dprop, &
               dprop, dprop, dprop, dprop, iprop, ierr)
            if (ierr /= 0) then
               write(*,2) 'year_month_day_when_created', year_month_day_when_created
               write(*,*) 'net_name' // trim(s% net_name)
               stop 'read_zams_header'
               return
            end if

            if (year_month_day_when_created < min_when_created) then
               ierr = -1
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, *) 'sorry: you need to update the zams data file.'
               write(*, *) 'found "year_month_day_when_created" =', year_month_day_when_created
               write(*, *) 'but need at least', min_when_created
               write(*, *)
               write(*, *)
               write(*, *)
               return
            end if
            if (abs(initial_z - s% initial_z) > 1d-3*s% initial_z) then
               ierr = -1
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, *) 'WARNING: requested initial_z does not match zams file initial_z.'
               write(*, 1) 'zams file initial_z', initial_z
               write(*, 1) 'requested initial_z', s% initial_z
               write(*, *)
               write(*, *)
               write(*, *)
               return
            end if
            if (s% initial_y > 0 .and. &
                  abs(initial_y - s% initial_y) > 1d-3*s% initial_y) then
               ierr = -1
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, *) 'WARNING: requested initial_y does not match zams file initial_y.'
               write(*, 1) 'zams file initial_y', initial_y
               write(*, 1) 'requested initial_y', s% initial_y
               write(*, *)
               write(*, *)
               write(*, *)
               return
            end if
         end subroutine read_zams_header

      end subroutine get1_zams_model


      subroutine get1_mass( &
            s, iounit, m1, nz1, m2, nz2, initial_mass, &
            nvar_hydro, species, xh, xa, q, dq, &
            omega, j_rot, ierr)
         use read_model, only: read_properties, read1_model
         use chem_def, only: iso_name_length
         use read_model, only: get_chem_col_names
         use star_utils, only: interp_q
         type (star_info), pointer :: s
         integer, intent(in) :: iounit, nz1, nz2, nvar_hydro, species
         real(dp), intent(in) :: m1, m2, initial_mass
         real(dp), intent(inout) :: xh(:,:) ! (nvar_hydro,nz1)
         real(dp), intent(inout) :: xa(:,:) ! (species,nz1)
         real(dp), intent(inout), dimension(:) :: &
            q, dq, omega, j_rot ! (nz1)
         integer, intent(out) :: ierr

         integer :: i, k, nz, nz_in, iprop
         real(dp) :: m_in, m_read, dprop, lnm1, lnm2
         real(dp), dimension(:, :), pointer :: xh2, xa2
         real(dp), dimension(:), pointer :: &
            q2, dq2, omega2, j_rot2
         real(dp) :: alfa, struct(nvar_hydro), comp(species)
         logical :: okay
         character (len=net_name_len) :: net_name
         character(len=iso_name_length), pointer :: names(:) ! (species)
         integer, pointer :: perm(:) ! (species)

         include 'formats'

         nz = nz1
         m_read = m1

         allocate( &
            xh2(nvar_hydro, nz2), xa2(species, nz2), q2(nz2), dq2(nz2), &
            omega2(nz2), j_rot2(nz2), &
            names(species), perm(species), stat=ierr)
         if (ierr /= 0) return
         okay = .false.
         mass_loop: do ! loop until find desired mass

            m_in = -1; nz_in = -1; net_name = ''
            call read_properties(iounit, &
               net_name, iprop, nz_in, iprop, m_in, &
               dprop, dprop, dprop, iprop, &
               dprop, dprop, dprop, dprop, dprop, &
               dprop, dprop, dprop, dprop, dprop, &
               dprop, dprop, dprop, dprop, dprop, iprop, ierr)
            if (ierr /= 0 .or. m_in < 0 .or. nz_in < 0) then
               write(*,*) 'missing required properties'
               write(*,*) 'ierr', ierr
               write(*,*) 'm_in', m_in
               write(*,*) 'nz_in', nz_in
               ierr = -1
               exit
            end if

            call get_chem_col_names(s, iounit, species, names, perm, ierr)
            if (ierr /= 0) exit

            if (abs(m_in-m_read) > 1d-4) then

               do i = 1, nz_in ! skip this one
                  read(iounit, *, iostat=ierr)
                  if (ierr /= 0) exit mass_loop
               end do

            else ! store this one

               if (nz /= nz_in) then
                  write(*, '(a, 2i6)') &
                     'nz /= nz_in: logic error in routine for reading model file', nz, nz_in
                  ierr = -1
                  return
               end if

               if (m_read == m1) then
                  call read1_model( &
                     s, species, nvar_hydro, nz, iounit, &
                     .false., .false., .false., .false., &
                     xh, xa, q, dq, omega, j_rot, perm, ierr)
                  if (ierr /= 0) exit mass_loop
                  okay = .true.
                  if (m2 == m1) exit mass_loop
                  m_read = m2
                  nz = nz2
               else
                  call read1_model( &
                     s, species, nvar_hydro, nz, iounit, &
                     .false., .false., .false., .false., &
                     xh2, xa2, q2, dq2, omega2, j_rot2, perm, ierr)
                  if (ierr /= 0) exit mass_loop
                  okay = .true.
                  nz = nz1
                  exit mass_loop
               end if

            end if
            read(iounit, *, iostat=ierr) ! skip line following the last zone
            if (ierr /= 0) exit mass_loop

         end do mass_loop

         if (.not. okay) then
            ierr = -1
            call dealloc
            return
         end if

         if (m1 /= m2) then ! interpolate linearly in log(m)
            lnm1 = log(m1)
            lnm2 = log(m2)
            alfa = (log(initial_mass) - lnm2) / (lnm1 - lnm2)
            do k=1,nz1
               call interp_q( &
                  nz2, nvar_hydro, species, q(k), xh2, xa2, q2, dq2, struct, comp, ierr)
               if (ierr /= 0) then
                  call dealloc
                  return
               end if
               xh(1:nvar_hydro,k) = alfa*xh(1:nvar_hydro,k) + (1-alfa)*struct(1:nvar_hydro)
               xa(1:species,k) = alfa*xa(1:species,k) + (1-alfa)*comp(1:species)
            end do
            call dealloc
         end if

         contains

         subroutine dealloc
            deallocate(xh2, xa2, q2, dq2, omega2, j_rot2, names, perm)
         end subroutine dealloc

      end subroutine get1_mass


      end module init_model

