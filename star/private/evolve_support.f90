 ! ***********************************************************************
!
!   Copyright (C) 2012-2019  The MESA Team
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

      module evolve_support

      use star_private_def
      use const_def

      implicit none

      private
      public :: set_current_to_old, new_generation, output, output_to_file


      contains


      subroutine new_generation(s, ierr)
         use utils_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: nz, j, k

         include 'formats'

         ierr = 0
         nz = s% nz
         
         if (.not. s% rsp_flag) then

            call copy_to_old(s% dq, s% dq_old, ierr)
            if (ierr /= 0) return

            call copy_to_old(s% omega, s% omega_old, ierr)
            if (ierr /= 0) return

            call copy_to_old(s% j_rot, s% j_rot_old, ierr)
            if (ierr /= 0) return
            
            call copy_to_old(s% mlt_vc, s% mlt_vc_old, ierr)
            if (ierr /= 0) return

            call enlarge_if_needed_2(s% xh_old,s% nvar_hydro,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(s% xh_old)

            call enlarge_if_needed_2(s% xa_old,s% species,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(s% xh_old)

            do k = 1, s% nz
               do j=1, s% nvar_hydro
                  s% xh_old(j,k) = s% xh(j,k)
               end do
               do j=1,s% species
                  s% xa_old(j,k) = s% xa(j,k)
               end do
            end do
         
         end if
         
         s% model_number_old = s% model_number
         s% nz_old = s% nz
         s% time_old = s% time
         s% dt_old = s% dt
         s% mstar_old = s% mstar
         s% xmstar_old = s% xmstar
         s% M_center_old = s% M_center
         s% R_center_old = s% R_center
         s% L_center_old = s% L_center
         s% v_center_old = s% v_center
         s% cumulative_energy_error_old = s% cumulative_energy_error
         s% total_energy_old = s% total_energy
         s% total_internal_energy_old = s% total_internal_energy
         s% total_angular_momentum_old = s% total_angular_momentum
         s% Teff_old = s% Teff
         s% power_nuc_burn_old = s% power_nuc_burn
         s% power_h_burn_old = s% power_h_burn
         s% power_he_burn_old = s% power_he_burn
         s% power_z_burn_old = s% power_z_burn
         s% power_photo_old = s% power_photo         
         s% mstar_dot_old = s% mstar_dot
         s% L_phot_old = s% L_phot
         s% L_surf_old = s% L_surf
         s% dt_limit_ratio_old = s% dt_limit_ratio
         s% gradT_excess_alpha_old = s% gradT_excess_alpha

         do j = 1, s% len_extra_work
            s% extra_work_old(j) = s% extra_work(j)
         end do

         do j = 1, s% len_extra_iwork
            s% extra_iwork_old(j) = s% extra_iwork(j)
         end do

         s% ixtra_old = s% ixtra
         s% xtra_old = s% xtra
         s% lxtra_old = s% lxtra

         call s% other_new_generation(s% id, ierr)
         
         s% need_to_setvars = .true.

         contains

         subroutine copy_to_old(ptr, ptr_old, ierr)
            real(dp), pointer, dimension(:) :: ptr, ptr_old
            integer, intent(out) :: ierr
            logical :: first_time
            ierr = 0
            first_time = (.not. associated(ptr_old))
            call realloc_if_needed_1(ptr_old,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) then
               call fill_with_NaNs(ptr_old)
            else if (s% zero_when_allocate) then
               ptr_old(:) = 0
            else if (first_time) then
               ptr_old(1:nz) = -9d99
            end if
            do j = 1, s% nz
               ptr_old(j) = ptr(j)
            end do
         end subroutine copy_to_old

      end subroutine new_generation


      subroutine set_current_to_old(s)
         use star_utils, only: total_angular_momentum, set_m_and_dm, set_dm_bar, set_qs
         use hydro_rotation, only: use_xh_to_update_i_rot
         use utils_lib
         type (star_info), pointer :: s
         real(dp), pointer :: p1(:)
         integer :: j, k, ierr

         include 'formats'

         s% nz = s% nz_old
         s% mstar = s% mstar_old
         s% xmstar = s% xmstar_old
         s% M_center = s% M_center_old
         s% R_center = s% R_center_old
         s% L_center = s% L_center_old
         s% v_center = s% v_center_old
         s% cumulative_energy_error = s% cumulative_energy_error_old
         s% total_energy = s% total_energy_old
         s% total_internal_energy = s% total_internal_energy_old
         s% total_angular_momentum = s% total_angular_momentum_old
         s% Teff = s% Teff_old
         s% power_nuc_burn = s% power_nuc_burn_old
         s% power_h_burn = s% power_h_burn_old
         s% power_he_burn = s% power_he_burn_old
         s% power_z_burn = s% power_z_burn_old
         s% power_photo = s% power_photo_old
         s% mstar_dot = s% mstar_dot_old
         s% L_phot = s% L_phot_old
         s% L_surf = s% L_surf_old
         s% dt_limit_ratio = s% dt_limit_ratio_old
         s% gradT_excess_alpha = s% gradT_excess_alpha_old
         
         if (.not. s% RSP_flag) then
            if (s% fill_arrays_with_NaNs) then
               call fill_with_NaNs_2d(s% xh)
               call fill_with_NaNs_2d(s% xa)
            end if
            do k = 1, s% nz
               do j=1, s% nvar_hydro
                  s% xh(j,k) = s% xh_old(j,k)
               end do
               do j=1,s% species
                  s% xa(j,k) = s% xa_old(j,k)
               end do
               s% dq(k) = s% dq_old(k)
               s% mlt_vc(k) = s% mlt_vc_old(k)
            end do
            s% okay_to_set_mlt_vc = .true.
            
            call set_qs(s, s% nz, s% q, s% dq, ierr)
            if (ierr /= 0) then
               write(*,*) 'set_current_to_old failed in set_qs'
               stop
            end if
            call set_m_and_dm(s)
            call set_dm_bar(s, s% nz, s% dm, s% dm_bar)

            if (s% rotation_flag) then
               do k=1,s% nz
                  s% j_rot(k) = s% j_rot_old(k)
               end do
               call use_xh_to_update_i_rot(s)
               do k=1,s% nz
                  s% omega(k) = s% j_rot(k)/s% i_rot(k)% val
                  if (is_bad_num(s% omega(k)) .or. abs(s% omega(k)) > 1d50) then
                     if (s% stop_for_bad_nums) then
                        write(*,2) 's% omega(k)', k, s% omega(k)
                        stop 'set_current_to_old'
                     end if
                  end if
               end do
               s% total_angular_momentum = total_angular_momentum(s)
            end if

         end if

         do j = 1, s% len_extra_work
            s% extra_work(j) = s% extra_work_old(j)
         end do

         do j = 1, s% len_extra_iwork
            s% extra_iwork(j) = s% extra_iwork_old(j)
         end do

         s% ixtra = s% ixtra_old
         s% xtra = s% xtra_old
         s% lxtra = s% lxtra_old

         call s% other_set_current_to_old(s% id)

         s% need_to_setvars = .true.

      end subroutine set_current_to_old


      subroutine output(id, ierr)
         use star_utils, only: get_name_for_restart_file
         interface
            subroutine save_restart_info(iounit, id, ierr)
               integer, intent(in) :: iounit
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine save_restart_info
         end interface
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         character (len=strlen) :: filename, num_str, fstring
         type (star_info), pointer :: s
         integer :: num_digits

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call get_name_for_restart_file(s% model_number, s% photo_digits, num_str)
         filename = trim(s% photo_directory) // '/' // trim(num_str)
         call output_to_file(filename, id, ierr)
         if (ierr /= 0) return

         write(*, '(a)', advance='no') 'save ' // trim(filename)
         num_digits = 1 + log10(dble(max(1,s% model_number)))
         write(fstring,'( "(a,i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(*,fstring) ' for model ', s% model_number

      end subroutine output


      subroutine output_to_file(filename, id, ierr)
         use photo_out, only: output_star_photo
         use utils_lib, only: folder_exists, mkdir
         character (len=*) :: filename
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         integer :: iounit, k
         type (star_info), pointer :: s
         character(len=strlen) :: iomsg

         include 'formats'

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if(.not. folder_exists(trim(s% photo_directory))) call mkdir(trim(s% photo_directory))

         open(newunit=iounit, file=trim(filename), action='write', &
            status='replace', iostat=ierr, iomsg=iomsg, form='unformatted')
         if (ierr == 0) then
            s% most_recent_photo_name = trim(filename)
            call output_star_photo(s, iounit, ierr)
            close(iounit)
         else
            write(*,*) trim(iomsg)
         endif

      end subroutine output_to_file


      end module evolve_support


