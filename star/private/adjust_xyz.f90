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


      module adjust_xyz

      use star_private_def
      use const_def
      use chem_def
      use utils_lib

      implicit none

      logical, parameter :: dbg = .false.


      contains


      subroutine change_net( &
            id, adjust_abundances_for_new_isos, new_net_name, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: adjust_abundances_for_new_isos
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr
         call do_composition_fixup( &
            id, adjust_abundances_for_new_isos, new_net_name, ierr)
      end subroutine change_net


      subroutine change_small_net( &
            id, adjust_abundances_for_new_isos, new_net_name, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: adjust_abundances_for_new_isos
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr
         call do_composition_fixup( &
            id, adjust_abundances_for_new_isos, new_net_name, ierr)
      end subroutine change_small_net


      subroutine do_composition_fixup( &
            id, adjust_abundances_for_new_isos, new_net_name, ierr)
         use net, only: do_micro_change_net
         use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results
         integer, intent(in) :: id
         logical, intent(in) :: adjust_abundances_for_new_isos
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: old_num_isos, old_num_reactions
         integer :: old_net_iso(num_chem_isos)
         integer :: old_chem_id(num_chem_isos)
         integer :: species, nz, num_reactions
         integer, pointer :: chem_id(:)
         character (len=net_name_len) :: old_net_name

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (new_net_name == s% net_name) return

         old_net_name = s% net_name
         old_num_isos = s% species
         old_num_reactions = s% num_reactions
         old_net_iso = s% net_iso
         old_chem_id(1:old_num_isos) = s% chem_id(1:old_num_isos)

         call do_micro_change_net(s, new_net_name, ierr)
         if (ierr /= 0) then
            s% retry_message = 'micro_change_net failed'
            if (s% report_ierr) write(*,*) s% retry_message
            return
         end if

         species = s% species
         nz = max(s% nz, s% prev_mesh_nz)
         chem_id => s% chem_id
         num_reactions = s% num_reactions
         
         write(*,*) 'change to "' // trim(new_net_name)//'"'
         write(*,*) 'number of species', s% species

         call realloc(s% dxdt_nuc); if (ierr /= 0) return
         call realloc(s% dxdt_nuc_start); if (ierr /= 0) return
         call realloc(s% d_epsnuc_dX); if (ierr /= 0) return
         call realloc(s% d_dxdt_nuc_dRho); if (ierr /= 0) return
         call realloc(s% d_dxdt_nuc_dT); if (ierr /= 0) return

         call realloc(s% dxdt_mix); if (ierr /= 0) return

         call realloc(s% d_eps_grav_dX); if (ierr /= 0) return
         call realloc(s% dlnE_dxa_for_partials); if (ierr /= 0) return
         call realloc(s% dlnP_dxa_for_partials); if (ierr /= 0) return

         call realloc(s% extra_diffusion_factor); if (ierr /= 0) return
         call realloc(s% edv); if (ierr /= 0) return
         call realloc(s% v_rad); if (ierr /= 0) return
         call realloc(s% g_rad); if (ierr /= 0) return
         call realloc(s% typical_charge); if (ierr /= 0) return
         call realloc(s% diffusion_dX); if (ierr /= 0) return
         call realloc(s% diffusion_D_self); if (ierr /= 0) return

         if (associated(s% d_dXdt_nuc_dX)) deallocate(s% d_dXdt_nuc_dX)
         allocate(s% d_dXdt_nuc_dX(species, species, nz + nz_alloc_extra), stat=ierr)
         if (ierr /= 0) return

         if (associated(s% d_eos_dxa)) deallocate(s% d_eos_dxa)
         allocate(s% d_eos_dxa(num_eos_d_dxa_results, species, nz + nz_alloc_extra), stat=ierr)
         if (ierr /= 0) return
         
         call realloc(s% xa_sub_xa_start); if (ierr /= 0) return
         call realloc(s% xa_start); if (ierr /= 0) return
         call realloc(s% prev_mesh_xa); if (ierr /= 0) return
         call do_xa(s% nz, s% xh, s% xa)
         if (s% generations > 1) call do_xa(s% nz_old, s% xh_old, s% xa_old)

         s% need_to_setvars = .true.
         s% prev_mesh_species_or_nvar_hydro_changed = .true.

         contains

         subroutine do_xa(nz, xh, xa_startv)
            use net_lib, only: clean_up_fractions
            integer, intent(in) :: nz
            real(dp), pointer :: xh(:,:)
            real(dp), pointer :: xa_startv(:,:)
            real(dp), pointer :: xa_new(:,:)
            real(dp), parameter :: max_sum_abs = 10d0
            real(dp), parameter :: xsum_tol = 1d-2
            integer :: k, i
            allocate(xa_new(species, nz + nz_alloc_extra), stat=ierr)
            if (ierr /= 0) return
            call set_x_new( &
               s, adjust_abundances_for_new_isos, &
               nz, old_num_isos, species, xh, xa_startv, xa_new, &
               old_chem_id, old_net_iso, chem_id, ierr)
            if (associated(xa_startv)) deallocate(xa_startv)
            xa_startv => xa_new
         end subroutine do_xa

         subroutine realloc(ptr)
            real(dp), pointer :: ptr(:, :)
            if (associated(ptr)) deallocate(ptr)
            allocate(ptr(species, nz + nz_alloc_extra), stat=ierr)
         end subroutine realloc

         subroutine realloc_integer(ptr)
            integer, pointer :: ptr(:, :)
            if (associated(ptr)) deallocate(ptr)
            allocate(ptr(species, nz + nz_alloc_extra), stat=ierr)
         end subroutine realloc_integer

      end subroutine do_composition_fixup


      subroutine set_x_new( &
            s, adjust_abundances_for_new_isos, &
            nz, old_num_isos, species, xh, xa_startv, xa_new, &
            old_chem_id, old_net_iso, chem_id, ierr)
         use net_lib, only: clean1
         type (star_info), pointer :: s
         logical, intent(in) :: adjust_abundances_for_new_isos
         integer, intent(in) :: nz, old_num_isos, species
         real(dp), pointer :: xh(:,:)
         real(dp), pointer :: xa_startv(:,:)
         real(dp), pointer :: xa_new(:,:)
         integer, pointer :: chem_id(:)
         integer, intent(in) :: old_chem_id(old_num_isos), old_net_iso(num_chem_isos)
         integer, intent(out) :: ierr

         integer :: k, op_err
         real(dp), parameter :: max_sum_abs = 10
         real(dp), parameter :: xsum_tol = 1d-2
         include 'formats'

         ierr = 0
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k=1, nz
            call set_new_abundances( &
               s, adjust_abundances_for_new_isos, &
               nz, old_num_isos, species, xh, xa_startv, xa_new, k, &
               old_chem_id, old_net_iso, chem_id, op_err)
            if (op_err /= 0) ierr = op_err
            call clean1(species, xa_new(1:species,k), max_sum_abs, xsum_tol, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
         if (ierr /= 0) then
            s% retry_message = 'set_new_abundances failed in set_x_new'
            if (s% report_ierr) write(*, *) s% retry_message
            return
         end if

         s% need_to_setvars = .true.

      end subroutine set_x_new


      subroutine set_new_abundances( &
            s, adjust_abundances_for_new_isos, nz, &
            old_num_isos, species, xh, xa_startv, xa_new, &
            k, old_chem_id, old_net_iso, chem_id, ierr)
         use chem_lib, only: chem_Xsol
         type (star_info), pointer :: s
         logical, intent(in) :: adjust_abundances_for_new_isos
         integer, intent(in) :: nz, old_num_isos, species, k
         real(dp), pointer :: xh(:,:)
         real(dp), pointer :: xa_startv(:,:)
         real(dp), pointer :: xa_new(:,:)
         integer, pointer :: chem_id(:)
         integer, intent(in) :: &
            old_chem_id(old_num_isos), old_net_iso(num_chem_isos)
         integer, intent(out) :: ierr

         real(dp) :: &
            old_total_neut, old_total_h, old_total_he, &
            old_total_c, old_total_n, old_total_o, old_other, &
            total_neut, total_h, total_he, total_c, total_n, total_o, &
            other, lgT, zfrac, xsol, lgT_lo, lgT_hi
         integer :: i, j, cid, Z
         character(len=solnamelen) :: sol_name

         logical :: dbg, did_total_neut, did_total_h, did_total_he, &
            did_total_c, did_total_n, did_total_o, did_total_other

         include 'formats'

         dbg = .false. !(k == nz) .or. (k == nz-1)

         ierr = 0
         zfrac = s% initial_z / zsol
         lgT_lo = s% lgT_lo_for_set_new_abundances
         lgT_hi = s% lgT_hi_for_set_new_abundances

         if (dbg) then
            write(*,*)
            write(*,2) 'set_new_abundances', k
            write(*,2) 'old_num_isos', old_num_isos
            do j=1,old_num_isos
               write(*,2) 'old ' // chem_isos% name(old_chem_id(j)), k, xa_startv(j,k)
            end do
            write(*,*)
         end if

         old_total_neut = 0
         old_total_h = 0
         old_total_he = 0
         old_total_c = 0
         old_total_n = 0
         old_total_o = 0
         old_other = 0
         do j=1, old_num_isos
            Z = chem_isos% Z(old_chem_id(j))
            if (Z == 0) then
               old_total_neut = old_total_neut + xa_startv(j,k)
            else if (Z == 1) then
               old_total_h = old_total_h + xa_startv(j,k)
            else if (Z == 2) then
               old_total_he = old_total_he + xa_startv(j,k)
            else if (Z == 6) then
               old_total_c = old_total_c + xa_startv(j,k)
            else if (Z == 7) then
               old_total_n = old_total_n + xa_startv(j,k)
            else if (Z == 8) then
               old_total_o = old_total_o + xa_startv(j,k)
            else
               old_other = old_other + xa_startv(j,k)
            end if
         end do

         lgT = get_test_lgT(k)
         ! copy old isos and set abundances for new ones
         do j=1, species
            cid = chem_id(j)
            i = old_net_iso(cid) ! old index number for this species
            if (i /= 0) then
               xa_new(j,k) = xa_startv(i,k)
            else if (.not. adjust_abundances_for_new_isos) then
               xa_new(j,k) = 0d0
            else
               if (chem_isos% Z(cid) > 5) then
                  sol_name = chem_isos% name(cid)
                  xsol = chem_Xsol(sol_name)
                  xa_new(j,k) = zfrac*xsol
                  if (dbg) then
                     write(*,2) 'xa new ' // trim(sol_name), j, xa_new(j,k), xsol, zfrac
                  end if
               else
                  ! cannot simply add light elements to high temperature region
                  ! or will create an explosive situation that won't converge
                  if (lgT >= lgT_hi) then ! don't add any
                     xa_new(j,k) = 0
                  else
                     sol_name = chem_isos% name(cid)
                     xsol = chem_Xsol(sol_name)
                     if (lgT > lgT_lo) then ! interpolate
                        xa_new(j, k) = &
                           xsol*0.5d0*(1+cospi((lgT-lgT_lo)/(lgT_hi-lgT_lo)))
                     else ! add solar
                        xa_new(j,k) = xsol
                     end if
                  end if
               end if
            end if
         end do

         total_neut = 0
         total_h = 0
         total_he = 0
         total_c = 0
         total_n = 0
         total_o = 0
         other = 0
         do j=1, species
            select case(int(chem_isos% Z(chem_id(j))))
               case (0)
                  total_neut = total_neut + xa_new(j,k)
               case (1)
                  total_h = total_h + xa_new(j,k)
               case (2)
                  total_he = total_he + xa_new(j,k)
               case (6)
                  total_c = total_c + xa_new(j,k)
               case (7)
                  total_n = total_n + xa_new(j,k)
               case (8)
                  total_o = total_o + xa_new(j,k)
               case default
                  other = other + xa_new(j,k)
            end select
         end do

         if (dbg) then
            write(*,1) 'old_total_n', old_total_n
            write(*,1) 'total_n', total_n
            write(*,1) 'old_total_n/total_n', old_total_n/total_n
         end if

         did_total_neut = .false.
         did_total_h = .false.
         did_total_he = .false.
         did_total_c = .false.
         did_total_n = .false.
         did_total_o = .false.
         did_total_other = .false.
         do j=1, species
            select case(int(chem_isos% Z(chem_id(j))))
               case (0)
                  if (total_neut > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_total_neut/total_neut
                     did_total_neut = .true.
                  end if
               case (1)
                  if (total_h > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_total_h/total_h
                     did_total_h = .true.
                  end if
               case (2)
                  if (total_he > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_total_he/total_he
                     did_total_he = .true.
                  end if
               case (6)
                  if (total_c > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_total_c/total_c
                     did_total_c = .true.
                  end if
               case (7)
                  if (total_n > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_total_n/total_n
                     did_total_n = .true.
                  end if
               case (8)
                  if (total_o > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_total_o/total_o
                     did_total_o = .true.
                  end if
               case default
                  if (other > 0) then
                     xa_new(j,k) = xa_new(j,k)*old_other/other
                     did_total_other = .true.
                  end if
            end select
         end do

         ! check for leftovers and dump them into the last iso
         j = species
         if (old_total_neut > 0 .and. .not. did_total_neut) &
            xa_new(j,k) = xa_new(j,k) + old_total_neut
         if (old_total_h > 0 .and. .not. did_total_h) &
            xa_new(j,k) = xa_new(j,k) + old_total_h
         if (old_total_he > 0 .and. .not. did_total_he) &
            xa_new(j,k) = xa_new(j,k) + old_total_he
         if (old_total_c > 0 .and. .not. did_total_c) &
            xa_new(j,k) = xa_new(j,k) + old_total_c
         if (old_total_n > 0 .and. .not. did_total_n) &
            xa_new(j,k) = xa_new(j,k) + old_total_n
         if (old_total_o > 0 .and. .not. did_total_o) &
            xa_new(j,k) = xa_new(j,k) + old_total_o
         if (old_other > 0 .and. .not. did_total_other) &
            xa_new(j,k) = xa_new(j,k) + old_other

         if (abs(sum(xa_new(:,k)) - 1d0) > 0.1d0) then
!$omp critical (new_abund)
            write(*,*)
            write(*,2) 'bad sum: set_new_abundances', k, sum(xa_new(:,k))
            write(*,2) 'species', species
            do j=1,species
               write(*,2) trim(chem_isos% name(chem_id(j))), k, xa_new(j,k)
            end do
            write(*,*)
            write(*,2) 'old_total_neut', k, old_total_neut
            write(*,2) 'old_total_h', k, old_total_h
            write(*,2) 'old_total_he', k, old_total_he
            write(*,2) 'old_total_c', k, old_total_c
            write(*,2) 'old_total_n', k, old_total_n
            write(*,2) 'old_total_o', k, old_total_o
            write(*,2) 'old_other', k, old_other
            write(*,*)
            stop 'debug: set_new_abundances'
!$omp end critical (new_abund)

            return
         end if

         if (dbg) then
            write(*,*)
            write(*,2) 'new num species', species
            do j=1,species
               write(*,2) trim(chem_isos% name(chem_id(j))), j, xa_new(j,k)
            end do
            write(*,*)
            stop 'debug: set_new_abundances'
         end if

         contains

         real(dp) function get_test_lgT(loc)
            integer, intent(in) :: loc
            integer :: k, kmax, i_lnT
            i_lnT = s% i_lnT
            if (i_lnT == 0) then
               get_test_lgT = 0d0
            else
               get_test_lgT = xh(i_lnT,loc)/ln10
            end if
            kmax = min(min(s%nz,nz),size(s% mlt_D,dim=1))
            ! may be using mlt_D for different generation
            ! okay for this purpose, but don't want to go beyond size.
            ! Also make sure we would have the array set to something
            do k=loc+1, kmax
               if (s% mlt_D(k) < 1d2) return
               get_test_lgT = xh(i_lnT,k)/ln10
            end do
         end function get_test_lgT

      end subroutine set_new_abundances

      subroutine get_ag89_composition(s, species, xa, ierr)
         use chem_lib, only: chem_Xsol
         type (star_info), pointer :: s
         integer, intent(in) :: species
         real(dp) :: xa(species)
         integer, intent(out) :: ierr
         integer :: i, mg24
         ierr = 0
         do i=1, species
            xa(i) = chem_Xsol(chem_isos% name(s% chem_id(i)))
         end do
         mg24 = s% net_iso(img24)
         xa(mg24) = xa(mg24) + (1 - sum(xa(1:species)))
      end subroutine get_ag89_composition


      subroutine set_y(s, y, nzlo, nzhi, ierr)
         use star_utils, only: eval_current_abundance
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         real(dp), intent(in) :: y
         integer, intent(out) :: ierr

         real(dp) :: xh1, xhe3, xhe4, z, ratio, desired_xh1, desired_xhe4, &
            new_xh1, new_xhe3, new_xhe4, new_z, new_ratio
         include 'formats'

         ierr = 0

         if (y == 0) then ! convert all he4 to h1
            call set_abundance_ratio(s% id, ihe4, ih1, 0d0, nzlo, nzhi, ierr)
         end if

         xh1 = eval_current_abundance(s, s% net_iso(ih1), nzlo, nzhi, ierr)
         xhe3 = eval_current_abundance(s, s% net_iso(ihe3), nzlo, nzhi, ierr)
         xhe4 = eval_current_abundance(s, s% net_iso(ihe4), nzlo, nzhi, ierr)
         z = 1d0 - (xh1 + xhe3 + xhe4)
         ! keep xhe3 and z constant; change ratio of xh1 and xhe4
         desired_xhe4 = y - xhe3
         desired_xh1 = 1d0 - z - y

         if (desired_xh1 <= 0d0) then ! convert all h1 to he4
            ratio = 0d0
         else if (desired_xhe4 <= 0d0) then ! convert all he4 to h1
            ratio = 1d0
         else
            ratio = desired_xh1/desired_xhe4
         end if

         call set_abundance_ratio(s% id, ih1, ihe4, ratio, nzlo, nzhi, ierr)
         
      end subroutine set_y


      subroutine set_abundance_ratio(id, i1, i2, ratio, nzlo, nzhi, ierr)
         integer, intent(in) :: id, i1, i2, nzlo, nzhi
         real(dp), intent(in) :: ratio ! new x(i1)/x(i2)
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: k, j1, j2
         real(dp) :: xsum

         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         j1 = s% net_iso(i1)
         j2 = s% net_iso(i2)
         if (j1 == 0 .or. j2 == 0) then
            ierr = -1
            return
         end if

         do k=nzlo, nzhi
            xsum = s% xa(j1, k) + s% xa(j2, k)
            s% xa(j2, k) = xsum/(1+ratio)
            s% xa(j1, k) = xsum - s% xa(j2, k)
         end do

         s% need_to_setvars = .true.

      end subroutine set_abundance_ratio


      subroutine read_xa(s, species, xa, filename, ierr)
         use chem_def
         use chem_lib
         use utils_def
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(in) :: species
         real(dp) :: xa(species)
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         integer :: iounit, t, n, i, cid, j
         character (len=strlen) :: buffer, string
         logical, parameter :: dbg = .false.

         ierr = 0

         xa(1:species) = 0 ! default

         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open abundance specification file "' // trim(filename)//'"'
            return
         end if

         n = 0
         i = 0

         do
            t = token(iounit, n, i, buffer, string)
            select case(t)
               case(name_token) ! iso name
                  cid = chem_get_iso_id(string)
                  if (cid <= 0) then
                     write(*,*) 'reading ' // trim(filename)
                     write(*,*) 'unknown chem name ' // trim(string)
                     ierr = -1
                     call cleanup
                     return
                  end if
                  j = s% net_iso(cid)
                  if (j == 0) then
                     write(*,*) 'reading ' // trim(filename)
                     write(*,*) 'iso not in current net ' // trim(string)
                     ierr = -1
                     call cleanup
                     return
                  end if
                  t = token(iounit, n, i, buffer, string)
                  if (t /= name_token) then
                     write(*,*) 'reading ' // trim(filename)
                     write(*,*) 'failed in reading abundance value for ' // &
                        trim(chem_isos% name(cid))
                     ierr = -1
                     call cleanup
                     return
                  end if

                  !read(string,fmt=*,iostat=ierr) xa(j)
                  call str_to_double(string, xa(j), ierr)
                  if (ierr /= 0 .or. xa(j) > 1 .or. xa(j) < 0) then
                     write(*,*) 'reading ' // trim(filename)
                     write(*,*) 'invalid number given as abundance value for ' // &
                        trim(chem_isos% name(cid))
                     ierr = -1
                     call cleanup
                     return
                  end if
               case(eof_token)
                  exit
               case default
                  write(*,*) 'error in reading abundance file ' // trim(filename)
                  ierr = -1
                  call cleanup
                  return
            end select

         end do

         call cleanup

         if (abs(sum(xa(1:species)) - 1d0) > 1d-6) then
            write(*,*) 'reading ' // trim(filename)
            write(*,*) 'abundances fail to add to 1.0:  sum - 1 = ', sum(xa(:)) - 1d0
            ierr = -1
            return
         end if

         xa(1:species) = xa(1:species)/sum(xa)

         contains

         subroutine cleanup
            close(iounit)
         end subroutine cleanup

      end subroutine read_xa


      subroutine set_uniform_xa_from_file(id, file_for_uniform_xa, ierr)
         integer, intent(in) :: id
         character (len=*), intent(in) :: file_for_uniform_xa
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call doit(s% species)

         contains

         subroutine doit(species)
            integer, intent(in) :: species
            real(dp) :: xa(species)
            call read_xa(s, species, xa, file_for_uniform_xa, ierr)
            if (ierr /= 0) return
            call set_uniform_composition(id, species, xa, ierr)
         end subroutine doit

      end subroutine set_uniform_xa_from_file


      subroutine set_uniform_composition(id, species, xa, ierr)
         integer, intent(in) :: species
         real(dp), intent(in) :: xa(species)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call set_composition(id, 1, s% nz, species, xa, ierr)
      end subroutine set_uniform_composition


      subroutine set_composition(id, nzlo, nzhi, num_species, xa_new, ierr)
         integer, intent(in) :: nzlo, nzhi, num_species
         real(dp), intent(in) :: xa_new(num_species)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer :: j, k, nz, species
         type (star_info), pointer :: s
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s% nz
         species = s% species
         if (num_species /= species) then
            ierr = -1
            s% retry_message = 'set_composition requires number of species to match current model'
            if (s% report_ierr) write(*, *) s% retry_message
            return
         end if
         if (nzlo < 1 .or. nzhi > nz) then
            ierr = -1
            s% retry_message = 'set_composition requires nzlo and nzhi to be within 1 to current num zones'
            if (s% report_ierr) write(*, *) s% retry_message
            return
         end if
         if (abs(1d0-sum(xa_new(1:species))) > 1d-6) then
            ierr = -1
            s% retry_message = 'set_composition requires new mass fractions to add to 1.'
            if (s% report_ierr) write(*, *) s% retry_message
            return
         end if
         do k=nzlo,nzhi
            do j=1,species
               s% xa(j,k) = xa_new(j)
            end do
         end do
         ierr = 0

         s% need_to_setvars = .true.

      end subroutine set_composition


      subroutine set_standard_composition( &
            s, species, h1, h2, he3, he4, which_zfracs, &
            dump_missing_into_heaviest, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: species
         real(dp), intent(in) :: h1, h2, he3, he4 ! mass fractions
         integer, intent(in) :: which_zfracs ! defined in chem_def. e.g., GS98_zfracs
         logical, intent(in) :: dump_missing_into_heaviest
         integer, intent(out) :: ierr
         real(dp) :: xa(species)
         ierr = 0
         call get_xa_for_standard_metals(s, &
            species, s% chem_id, s% net_iso, h1, h2, he3, he4, which_zfracs, &
            dump_missing_into_heaviest, xa, ierr)
         if (ierr /= 0) return
         call set_uniform_composition(s% id, species, xa, ierr)
      end subroutine set_standard_composition


      subroutine get_xa_for_accretion(s, xa, ierr)
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         real(dp) :: xa(:) ! (species)
         integer, intent(out) :: ierr
         integer :: j, i, species, cid
         include 'formats'
         ierr = 0
         species = s% species
         if (s% accrete_given_mass_fractions) then
            xa(1:species) = 0
            do j=1,s% num_accretion_species
               if (len_trim(s% accretion_species_id(j)) == 0) cycle
               cid = chem_get_iso_id(s% accretion_species_id(j))
               if (cid <= 0) cycle
               i = s% net_iso(cid)
               if (i == 0) cycle
               xa(i) = s% accretion_species_xa(j)
            end do
            if (abs(1d0 - sum(xa(1:species))) > 1d-2) then
               write(*,'(a)') &
                  'get_xa_for_accretion: accretion species mass fractions do not add to 1.0'
               write(*,1) 'sum(xa(1:species))', sum(xa(1:species))
               do j=1,s% num_accretion_species
                  write(*,2) trim(s% accretion_species_id(j)), j, xa(j)
               end do
               ierr = -1
               return
            end if
            xa(1:species) = xa(1:species)/sum(xa(1:species))
            return
         end if
         call get_xa_for_standard_metals(s, &
            s% species, s% chem_id, s% net_iso, &
            s% accretion_h1, s% accretion_h2, s% accretion_he3, s% accretion_he4, &
            s% accretion_zfracs, &
            s% accretion_dump_missing_metals_into_heaviest, xa, ierr)
      end subroutine get_xa_for_accretion


      subroutine do_change_to_xa_for_accretion(id, nzlo_in, nzhi_in, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nzlo_in, nzhi_in
         integer, intent(out) :: ierr
         integer :: j, k, species
         real(dp), pointer :: xa(:) ! (species)
         type (star_info), pointer :: s
         integer :: nzlo, nzhi
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         species = s% species
         nzlo = nzlo_in
         nzhi = nzhi_in
         if (nzlo < 1) nzlo = 1
         if (nzhi < 1 ) nzhi = s% nz
         allocate(xa(species))
         if (s% accrete_same_as_surface) then
            do j=1,species
               xa(j) = s% xa(j,1)
            end do
         else
            call get_xa_for_accretion(s, xa, ierr)
            if (ierr /= 0) then
               s% retry_message = 'get_xa_for_accretion failed in change_to_xa_for_accretion'
               if (s% report_ierr) write(*, *) s% retry_message
               deallocate(xa)
               return
            end if
         end if
         do k=nzlo,min(s% nz,nzhi)
            do j=1,species
               s% xa(j,k) = xa(j)
            end do
         end do
         deallocate(xa)
         s% need_to_setvars = .true.
      end subroutine do_change_to_xa_for_accretion


      subroutine get_xa_for_standard_metals( &
            s, species, chem_id, net_iso, h1_in, h2_in, he3_in, he4_in, which_zfracs, &
            dump_missing_into_heaviest, xa, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: species, chem_id(:), net_iso(:), which_zfracs
         real(dp), intent(in) :: h1_in, h2_in, he3_in, he4_in
         logical, intent(in) :: dump_missing_into_heaviest
         real(dp) :: xa(:) ! (species)
         integer, intent(out) :: ierr
         real(dp) :: zfrac(num_chem_elements), Z, h1, h2, he3, he4
         integer :: i
         include 'formats'
         ierr = 0
         h1 = h1_in; h2 = h2_in; he3 = he3_in; he4 = he4_in
         select case(which_zfracs)
            case (AG89_zfracs)
               zfrac(:) = AG89_element_zfrac(:)
            case (GN93_zfracs)
               zfrac(:) = GN93_element_zfrac(:)
            case (GS98_zfracs)
               zfrac(:) = GS98_element_zfrac(:)
            case (L03_zfracs)
               zfrac(:) = L03_element_zfrac(:)
            case (AGS05_zfracs)
               zfrac(:) = AGS05_element_zfrac(:)
            case (AGSS09_zfracs)
               zfrac(:) = AGSS09_element_zfrac(:)
            case (L09_zfracs)
               zfrac(:) = L09_element_zfrac(:)
            case (A09_Prz_zfracs)
               zfrac(:) = A09_Prz_zfrac(:)
            case (0) ! use non-standard values given in controls
               zfrac(:) = 0
               zfrac(e_li) = s% z_fraction_li
               zfrac(e_be) = s% z_fraction_be
               zfrac(e_b)  = s% z_fraction_b
               zfrac(e_c)  = s% z_fraction_c
               zfrac(e_n)  = s% z_fraction_n
               zfrac(e_o)  = s% z_fraction_o
               zfrac(e_f)  = s% z_fraction_f
               zfrac(e_ne) = s% z_fraction_ne
               zfrac(e_na) = s% z_fraction_na
               zfrac(e_mg) = s% z_fraction_mg
               zfrac(e_al) = s% z_fraction_al
               zfrac(e_si) = s% z_fraction_si
               zfrac(e_p)  = s% z_fraction_p
               zfrac(e_s)  = s% z_fraction_s
               zfrac(e_cl) = s% z_fraction_cl
               zfrac(e_ar) = s% z_fraction_ar
               zfrac(e_k)  = s% z_fraction_k
               zfrac(e_ca) = s% z_fraction_ca
               zfrac(e_sc) = s% z_fraction_sc
               zfrac(e_ti) = s% z_fraction_ti
               zfrac(e_v)  = s% z_fraction_v
               zfrac(e_cr) = s% z_fraction_cr
               zfrac(e_mn) = s% z_fraction_mn
               zfrac(e_fe) = s% z_fraction_fe
               zfrac(e_co) = s% z_fraction_co
               zfrac(e_ni) = s% z_fraction_ni
               zfrac(e_cu) = s% z_fraction_cu
               zfrac(e_zn) = s% z_fraction_zn
               if (abs(sum(zfrac)-1) > 1d-6) then
                  write(*,*) 'bad sum(zfrac) for specified z fractions', sum(zfrac)
                  ierr = -1
               end if
               do i = 1, size(zfrac, dim=1)
                  if (zfrac(i) < 0) then
                     write(*,2) 'zfrac(i)', i, zfrac(i)
                     ierr = -1
                     return
                  end if
               end do
               zfrac(:) = zfrac(:) / sum(zfrac(:))
            case default
               if (abs(1d0 - (h1 + h2 + he3 + he4)) < 1d-8) then ! okay -- no metals
                  if (h1 > he4) then
                     h1 = max(0d0, min(1d0, 1d0 - (h2 + he3 + he4)))
                  else
                     he4 = max(0d0, min(1d0, 1d0 - (h1 + h2 + he3)))
                  end if
                  zfrac(:) = 0
               else
                  ierr = -1
               end if
         end select
         if (ierr /= 0) return
         Z = max(0d0, min(1d0, 1d0 - (h1 + h2 + he3 + he4)))
         call get_xa( &
            species, chem_id, net_iso, h1, h2, he3, he4, zfrac, Z, &
            dump_missing_into_heaviest, xa, ierr)
         do i=1,species
            if (xa(i) < 0 .or. is_bad(xa(i))) then
               write(*,2) 'get_xa_for_standard_metals xa(i)', i, xa(i)
               ierr = -1
               return
            end if
         end do
      end subroutine get_xa_for_standard_metals


      subroutine get_xa( &
            species, chem_id, net_iso, h1, h2, he3, he4, zfrac, Zinit, &
            dump_missing_into_heaviest, xa, ierr)
         use chem_def
         integer, intent(in) :: species, chem_id(:), net_iso(:)
         real(dp), intent(in) :: h1, h2, he3, he4, Zinit
         real(dp), intent(in) :: zfrac(:) ! (num_chem_elements)
         logical, intent(in) :: dump_missing_into_heaviest
         real(dp) :: xa(:) ! (species)
         real(dp), parameter :: tiny = 1d-15
         integer, intent(out) :: ierr
         include 'formats'
         if (h1 < 0 .and. abs(h1)>tiny &
             .or. h2 < 0 .and. abs(h2)>tiny &
             .or. he3 < 0.and. abs(he3)>tiny &
             .or. he4 < 0.and. abs(he4)>tiny &
             .or.Zinit < 0.and.abs(Zinit)>tiny ) then
            ierr = -1
            write(*,*) 'when setting composition, need to provide values for h1, h2, he3, and he4'
            write(*,*) 'H1=', h1
            write(*,*) 'H2=', h2
            write(*,*) 'He3=', he3
            write(*,*) 'He4=', he4
            write(*,*) 'Z  =', Zinit
            return
         else if ( abs(1d0-(h1+h2+he3+he4+Zinit)) > tiny ) then
            ierr = -2
            write(*,1) 'h1', h1
            write(*,1) 'h2', h2
            write(*,1) 'he3', he3
            write(*,1) 'he4', he4
            write(*,1) 'Zinit', Zinit
            write(*,*) 'sum of (H1+H2+He3+He4+Z) does not sum to 1: (1-sum)= ', &
               1d0-(h1+h2+he3+he4+Zinit)
            stop 'get_xa'
         end if
         ierr = 0
         xa(:) = 0
         if (net_iso(ih2) /= 0) then
            xa(net_iso(ih2)) = h2
            xa(net_iso(ih1)) = h1
         else if (net_iso(ih1) /= 0) then
            xa(net_iso(ih1)) = h1 + h2
         else
            ierr = -1
            write(*,*) 'require h1 to be in net'
            return
         end if
         if (net_iso(ihe3) /= 0) xa(net_iso(ihe3)) = he3
         if (net_iso(ihe4) /= 0) xa(net_iso(ihe4)) = he4
         call adjust_z_fractions( &
            species, chem_id, 1d0-(h1+h2+he3+he4), zfrac, &
            dump_missing_into_heaviest, xa, ierr)
         xa(:) = xa(:)/sum(xa)
      end subroutine get_xa


      subroutine adjust_z_fractions( &
            species, chem_id, ztotal, zfrac, dump_missing_into_heaviest, xa, ierr)
         use chem_def
         use chem_lib, only: lodders03_element_atom_percent
         integer, intent(in) :: species, chem_id(:) ! (species)
         real(dp), intent(in) :: ztotal
         real(dp), intent(in) :: zfrac(:) ! (num_chem_elements)
         logical, intent(in) :: dump_missing_into_heaviest
         real(dp) :: xa(:) ! (species) h & he not touched
         integer, intent(out) :: ierr

         integer :: i, j, cid, cid2, element_id, jskip, Z(species)
         real(dp) :: new_xa(species), frac, frac_sum, iso_frac, zsum

         include 'formats'

         ierr = 0
         new_xa(:) = 0

         if (dump_missing_into_heaviest) then
            jskip = maxloc(chem_isos% W(chem_id(:)),dim=1) ! find the heaviest
         else
            jskip = 0
         end if

         do j=1,species
            Z(j) = chem_isos% Z(chem_id(j))
         end do

         do j=1,species
            if (j == jskip) cycle
            if (Z(j) <= 2) cycle
            cid = chem_id(j)
            ! multiply lodders' number fractions by chem_isos% A since want fractions by mass
            frac = 1d-2*lodders03_element_atom_percent(chem_isos% name(cid))*chem_isos% W(cid)
            if (frac == 0) cycle
            frac_sum = 0
            do i=1,species
               if (Z(i) /= Z(j)) cycle
               cid2 = chem_id(i)
               frac_sum = frac_sum + &
                  1d-2*lodders03_element_atom_percent(chem_isos% name(cid2))*chem_isos% W(cid2)
            end do
            element_id = Z(j) ! element index equals number of protons
            iso_frac = frac/frac_sum ! this iso is iso_frac of the element abundance by mass
            new_xa(j) = ztotal*zfrac(element_id)*iso_frac
         end do

         if (jskip > 0) then
            ! put the missing metals into jskip
            new_xa(jskip) = max(0d0, ztotal-sum(new_xa(:)))
            do j=1,species
               if (j == jskip) cycle
               if (Z(j) <= 2) cycle
               xa(j) = new_xa(j)
            end do
            xa(jskip) = new_xa(jskip)
         else ! renormalize all of the metals
            zsum = sum(new_xa(:)) ! only have metals in new_xa
            if (zsum > 0d0) then
               do j=1,species
                  if (Z(j) <= 2) cycle
                  xa(j) = new_xa(j)*ztotal/zsum
               end do
            end if
         end if

      end subroutine adjust_z_fractions


      subroutine set_z(s, new_z, nzlo, nzhi, ierr)
         use net_lib, only:clean_up_fractions
         use star_utils, only: eval_current_z
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_z
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         integer :: nz, h1, he3, he4, species
         real(dp) :: frac_z, current_z, initial_z
         real(dp), parameter :: max_sum_abs = 10d0
         real(dp), parameter :: xsum_tol = 1d-2

         include 'formats'

         ierr = 0

         if (nzlo > nzhi) then
            ierr = -1; return
         end if

         nz = s% nz
         species = s% species
         h1 = s% net_iso(ih1)
         he3 = s% net_iso(ihe3)
         he4 = s% net_iso(ihe4)

         call clean_up_fractions(1, nz, species, nz, s% xa, max_sum_abs, xsum_tol, ierr)
         if (ierr /= 0) return

         current_z = eval_current_z(s, nzlo, nzhi, ierr)
         initial_z = current_z

         if (current_z <= 0 .and. new_z > 0) then
            ierr = -1
            return
         end if
         frac_z = new_z/current_z

         if (abs(1-frac_z) < 1d-20) return

         call convert
         if (ierr /= 0) return

         ! check
         current_z = eval_current_z(s, nzlo, nzhi, ierr)

         if (abs(current_z-new_z) > 1d-8 .or. is_bad(current_z)) then
            ierr = -1
            s% retry_message = 'set_z failed'
            if (s% report_ierr) then
               write(*, *) 'set_z failed'
               write(*, 1) 'initial_z', initial_z
               write(*, 1) 'requested new_z', new_z
               write(*, 1) 'actual current_z', current_z
               write(*, *)
               write(*, 1) 'requested/initial', new_z/initial_z
               write(*, 1) 'requested/final', new_z/current_z
               write(*, *)
            end if
            if (s% stop_for_bad_nums) then
               write(*, 1) 'initial_z', initial_z
               write(*, 1) 'requested new_z', new_z
               write(*, 1) 'actual current_z', current_z
               stop 'set_z'
            end if
         end if

         if (s% doing_first_model_of_run) then
            s% initial_z = new_z
            write(*,1) 's% initial_z =', new_z
         end if

         s% need_to_setvars = .true.


         contains

         subroutine convert
            integer :: i, k
            real(dp) :: old_z, old_xy, sumfracs, sum_z, frac_xy
            include 'formats'
            do k=nzlo, nzhi
               old_z = 0; old_xy = 0
               do i=species, 1, -1
                  if (i==h1 .or. i==he3 .or. i==he4) then
                     old_xy = old_xy + s% xa(i, k)
                  end if
               end do
               old_z = 1 - old_xy
               sum_z = 0
               do i=species, 1, -1
                  if (i==h1 .or. i==he3 .or. i==he4) cycle
                  s% xa(i, k) = s% xa(i, k)*frac_z
                  sum_z = sum_z + s% xa(i, k)
               end do
               ! sum_z is the new z for this cell
               ! modify x and y to make mass fractions sum to 1
               if (old_z < 1) then
                  frac_xy = (1 - old_z*frac_z)/old_xy
                  s% xa(h1, k) = s% xa(h1, k)*frac_xy
                  if (he3 /= 0) s% xa(he3, k) = s% xa(he3, k)*frac_xy
                  s% xa(he4, k) = s% xa(he4, k)*frac_xy
               else
                  s% xa(h1, k) = (1-sum_z)/2
                  if (he3 /= 0) s% xa(he3, k) = 0
                  s% xa(he4, k) = (1-sum_z)/2
               end if
               sumfracs = 0
               do i=species, 1, -1
                  sumfracs = sumfracs + s% xa(i, k)
               end do
               s% xa(1:species, k) = s% xa(1:species, k)/sumfracs
               do i=1, species
                  if (is_bad(s% xa(i, k))) then
                     ierr = -1
                     s% retry_message = 'set_z failed - bad mass fraction'
                     if (s% report_ierr) then
                        write(*,2) 'set_z s% xa(i, k)', k, s% xa(i, k)
                     end if
                     if (s% stop_for_bad_nums) then
                        write(*,2) 'set_z s% xa(i, k)', k, s% xa(i, k)
                        stop 'set_z'
                     end if
                     return
                  end if
               end do
            end do
         end subroutine convert

      end subroutine set_z


      subroutine do_replace(s, id1, id2, nzlo, nzhi, ierr)
         ! replaces species1 by species2
         type (star_info), pointer :: s
         integer, intent(in) :: id1, id2 ! values are chem_id's such as ihe4
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         integer :: k, nz, species, species1, species2

         ierr = 0
         nz = s% nz
         species = s% species
         species1 = s% net_iso(id1)
         species2 = s% net_iso(id2)

         do k=nzlo, nzhi
            s% xa(species2, k) = s% xa(species1, k) + s% xa(species2, k)
            s% xa(species1, k) = 0
         end do

         s% need_to_setvars = .true.

      end subroutine do_replace


      subroutine do_uniform_mix_section(s, species, nzlo_in, nzhi_in, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: species, nzlo_in, nzhi_in
         integer, intent(out) :: ierr
         integer :: j, k, nzlo, nzhi
         real(dp) :: dqsum, newxa(species), xsum
         include 'formats'
         ierr = 0
         nzlo = max(1, nzlo_in)
         nzhi = min(s% nz, nzhi_in)
         dqsum = sum(s% dq(nzlo:nzhi))
         do j=1,species
            newxa(j) = dot_product(s% dq(nzlo:nzhi), s% xa(j,nzlo:nzhi))/dqsum
         end do
         xsum = sum(newxa(:))
         do j=1,species
            newxa(j) = newxa(j)/xsum
         end do
         do k=nzlo,nzhi
            do j=1,species
               s% xa(j,k) = newxa(j)
            end do
         end do
         s% need_to_setvars = .true.
      end subroutine do_uniform_mix_section


      subroutine do_uniform_mix_envelope_down_to_T(s, T, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: T
         integer, intent(out) :: ierr
         integer :: k, kmix
         ierr = 0
         if (s% T(1) > T) return
         kmix = s% nz
         do k=2,s% nz
            if (s% T(k) > T) then
               kmix = k-1
               exit
            end if
         end do
         call do_uniform_mix_section(s, s% species, 1, kmix, ierr)
      end subroutine do_uniform_mix_envelope_down_to_T


      subroutine do_set_abundance(s, chem_id, new_frac, nzlo, nzhi, ierr)
         ! set mass fraction of species to new_frac uniformly in cells nzlo to nzhi
         type (star_info), pointer :: s
         integer, intent(in) :: chem_id
         real(dp), intent(in) :: new_frac
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         integer :: k, j, i, nz, species
         real(dp) :: old_frac, rescale

         ierr = 0
         nz = s% nz
         species = s% species
         j = s% net_iso(chem_id)
         if (j == 0) then
            write(*,*) 'do_set_abundance: failed to find requested iso in current net'
            ierr = -1
            return
         end if

         do k=nzlo, nzhi
            old_frac = s% xa(j, k)
            s% xa(j, k) = new_frac
            if (1d0-old_frac > 1d-10) then
               rescale = (1-new_frac)/(1-old_frac)
               do i=1, species
                  if (i==j) cycle
                  s% xa(i, k) = rescale*s% xa(i, k)
               end do
            else
               rescale = (1-new_frac)/(species-1)
               do i=1, species
                  if (i==j) cycle
                  s% xa(i, k) = rescale
               end do
            end if
         end do
         s% need_to_setvars = .true.

      end subroutine do_set_abundance


      end module adjust_xyz












