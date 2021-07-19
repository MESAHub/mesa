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

         real(dp), pointer :: tmp(:,:), tmp1(:)
         integer, pointer :: itmp1(:)
         integer :: nz, i, k, i_u

         include 'formats'

         ierr = 0
         nz = s% nz
         
         if (.not. s% rsp_flag) then

            call flip(s% dq, s% dq_old, ierr)
            if (ierr /= 0) return

            call flip(s% omega, s% omega_old, ierr)
            if (ierr /= 0) return

            call flip(s% j_rot, s% j_rot_old, ierr)
            if (ierr /= 0) return
            
            call flip(s% mlt_vc, s% mlt_vc_old, ierr)
            if (ierr /= 0) return

            tmp => s% xh_old
            s% xh_old => s% xh
            call enlarge_if_needed_2(tmp,s% nvar_hydro,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(tmp)
            s% xh => tmp

            tmp => s% xa_old
            s% xa_old => s% xa
            call enlarge_if_needed_2(tmp,s% species,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(tmp)
            s% xa => tmp
         
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
         s% crystal_core_boundary_mass_old = s% crystal_core_boundary_mass

         do i = 1, s% len_extra_work
            s% extra_work_old(i) = s% extra_work(i)
         end do

         do i = 1, s% len_extra_iwork
            s% extra_iwork_old(i) = s% extra_iwork(i)
         end do

         s% ixtra_old = s% ixtra
         s% xtra_old = s% xtra
         s% lxtra_old = s% lxtra

         call s% other_new_generation(s% id, ierr)
         
         s% need_to_setvars = .true.

         contains

         subroutine flip(ptr, ptr_old, ierr)
            real(dp), pointer, dimension(:) :: ptr, ptr_old
            integer, intent(out) :: ierr
            logical :: first_time
            ierr = 0
            tmp1 => ptr_old
            ptr_old => ptr
            first_time = (.not. associated(tmp1))
            call realloc_if_needed_1(tmp1,nz,nz_alloc_extra,ierr)
            if (ierr /= 0) return
            if (s% fill_arrays_with_NaNs) then
               call fill_with_NaNs(tmp1)
            else if (s% zero_when_allocate) then
               tmp1(:) = 0
            else if (first_time) then
               tmp1(1:nz) = -9d99
            end if
            ptr => tmp1
         end subroutine flip

      end subroutine new_generation


      subroutine set_current_to_old(s)
         type (star_info), pointer :: s
         real(dp), pointer :: p1(:)
         integer :: i, k

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
         s% crystal_core_boundary_mass = s% crystal_core_boundary_mass_old

         do i = 1, s% len_extra_work
            s% extra_work(i) = s% extra_work_old(i)
         end do

         do i = 1, s% len_extra_iwork
            s% extra_iwork(i) = s% extra_iwork_old(i)
         end do

         s% ixtra = s% ixtra_old
         s% xtra = s% xtra_old
         s% lxtra = s% lxtra_old

         call s% other_set_current_to_old(s% id)

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
         character (len=*) :: filename
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         integer :: iounit, k
         type (star_info), pointer :: s
         character(len=strlen) :: iomsg

         include 'formats'

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

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


