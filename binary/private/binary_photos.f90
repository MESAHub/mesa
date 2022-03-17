! ***********************************************************************
!
!   Copyright (C) 2010-2019  Pablo Marchant & The MESA Team
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


      module binary_photos

      use const_def
      use math_lib
      use star_lib
      use star_def
      use binary_def

      implicit none

      contains
         
      subroutine do_saves_for_binary(b, ierr)
         type(binary_info), pointer :: b
         integer, intent(out) :: ierr
         integer :: iounit, id
         character (len=strlen) :: str_photo, filename, iomsg, report_str

         call string_for_model_number('x', b% model_number, b% photo_digits, str_photo)

         filename = trim(trim(b% photo_directory) // '/b_' // str_photo)
         report_str = trim('save ' // filename)
         open(newunit=iounit, file=trim(filename), action='write', &
            status='replace', iostat=ierr, iomsg=iomsg, form='unformatted')
         if (ierr /= 0) then
            write(*,*) 'failed in do_saves_for_binary', trim(filename)
            return
         end if
         call binary_photo_write(b% binary_id, iounit)
         close(iounit)

         if (b% have_star_1) then
            filename = trim(trim(b% s1% ctrl% photo_directory) // '/1_' // str_photo)
            call star_save_for_restart(b% s1% id, filename, ierr)
            report_str = trim(trim(report_str) // ', ' // filename)
         end if
         if (b% have_star_2) then
            filename = trim(trim(b% s2% ctrl% photo_directory) // '/2_' // str_photo)
            call star_save_for_restart(b% s2% id, filename, ierr)
            report_str = trim(trim(report_str) // ', ' // filename)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed in do_saves_for_binary'
            return
         end if

         write(*,*) trim(trim(report_str) // ' for model'), b% model_number
         
      end subroutine do_saves_for_binary

      subroutine binary_photo_write(binary_id, iounit)
         integer, intent(in) :: binary_id, iounit
         type(binary_info), pointer :: b

         integer :: ierr

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         write(iounit) star_def_version

         write(iounit, iostat=ierr) &
             b% binary_age, b% binary_age_old, &
             b% model_number, b% model_number_old, &
             b% mtransfer_rate, b% mtransfer_rate_old, &
             b% angular_momentum_j, b% angular_momentum_j_old, & 
             b% separation, b% separation_old, &
             b% eccentricity, b% eccentricity_old, &
             b% rl_relative_gap(1), b% rl_relative_gap_old(1), &
             b% rl_relative_gap(2), b% rl_relative_gap_old(2), &
             b% r(1), b% r_old(1), &
             b% r(2), b% r_old(2), &
             b% rl(1), b% rl_old(1), &
             b% rl(2), b% rl_old(2), &
             b% m(1), b% m_old(1), &
             b% m(2), b% m_old(2), &
             b% dt, b% dt_old, &
             b% env, b% env_old, &
             b% eq_initial_bh_mass, &
             b% period, b% period_old, & 
             b% max_timestep, b% max_timestep_old, &
             b% change_factor, b% change_factor_old, &
             b% min_binary_separation, &
             b% using_jdot_mb(1), b% using_jdot_mb_old(1), &
             b% using_jdot_mb(2), b% using_jdot_mb_old(2), &
             b% d_i, b% d_i_old, b% a_i, b% a_i_old, &
             b% point_mass_i, b% point_mass_i_old, &
             b% ignore_rlof_flag, b% ignore_rlof_flag_old, &
             b% model_twins_flag, b% model_twins_flag_old, &
             b% dt_why_reason, b% dt_why_reason_old, &
             b% have_star_1, b% have_star_2, &
             b% CE_flag, b% CE_flag_old, &
             b% CE_init, b% CE_init_old, &
             b% CE_nz, b% CE_initial_radius, b% CE_initial_separation, b% CE_initial_Mdonor, &
             b% CE_initial_Maccretor, b% CE_initial_age, b% CE_initial_model_number, &
             b% CE_b_initial_age, b% CE_b_initial_model_number, &
             b% CE_num1, b% CE_num1_old, &
             b% CE_num2, b% CE_num2_old, &
             b% CE_lambda1, b% CE_lambda1_old, &
             b% CE_lambda2, b% CE_lambda2_old, &
             b% CE_Ebind1, b% CE_Ebind1_old, &
             b% CE_Ebind2, b% CE_Ebind2_old, &
             b% ixtra(:), b% ixtra_old(:), &
             b% xtra(:), b% xtra_old(:), &
             b% lxtra(:), b% lxtra_old(:)

         if (b% CE_init) then
            write(iounit, iostat=ierr) &
               b% CE_m(:), b% CE_entropy(:), b% CE_U_in(:), b% CE_U_out(:), b% CE_Omega_in(:), b% CE_Omega_out(:)
         end if

         call b% other_binary_photo_write(binary_id, iounit)

         if (ierr /= 0) stop "error in binary_photo_write"

      end subroutine binary_photo_write

      subroutine binary_load_photo(b, photo_filename, ierr)
         type(binary_info), pointer :: b
         character (len=strlen) :: photo_filename
         integer, intent(out) :: ierr
         integer :: iounit, version

         open(newunit=iounit, file=trim(photo_filename), action='read', &
            status='old', iostat=ierr, form='unformatted')
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(photo_filename)
            return
         end if

         read(iounit, iostat=ierr) version
         if (ierr /= 0) then
            write(*,*) 'failed to read version number'
            return
         end if
         if (version /= star_def_version) then
            write(*,'(/,a,/)') ' FAILURE: the restart data' // &
               ' is from a previous version of the code and is no longer usable.'
            ierr = -1
            return
         end if

         call binary_photo_read(b% binary_id, iounit, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_photo_read'
            return
         end if

         close(iounit)
         
      end subroutine binary_load_photo

      subroutine binary_photo_read(binary_id, iounit, ierr)
         integer, intent(in) :: binary_id, iounit
         integer, intent(out) :: ierr
         type(binary_info), pointer :: b
         integer :: nz

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         read(iounit, iostat=ierr) &
             b% binary_age, b% binary_age_old, &
             b% model_number, b% model_number_old, &
             b% mtransfer_rate, b% mtransfer_rate_old, &
             b% angular_momentum_j, b% angular_momentum_j_old, & 
             b% separation, b% separation_old, &
             b% eccentricity, b% eccentricity_old, &
             b% rl_relative_gap(1), b% rl_relative_gap_old(1), &
             b% rl_relative_gap(2), b% rl_relative_gap_old(2), &
             b% r(1), b% r_old(1), &
             b% r(2), b% r_old(2), &
             b% rl(1), b% rl_old(1), &
             b% rl(2), b% rl_old(2), &
             b% m(1), b% m_old(1), &
             b% m(2), b% m_old(2), &
             b% dt, b% dt_old, &
             b% env, b% env_old, &
             b% eq_initial_bh_mass, &
             b% period, b% period_old, & 
             b% max_timestep, b% max_timestep_old, &
             b% change_factor, b% change_factor_old, &
             b% min_binary_separation, &
             b% using_jdot_mb(1), b% using_jdot_mb_old(1), &
             b% using_jdot_mb(2), b% using_jdot_mb_old(2), &
             b% d_i, b% d_i_old, b% a_i, b% a_i_old, &
             b% point_mass_i, b% point_mass_i_old, &
             b% ignore_rlof_flag, b% ignore_rlof_flag_old, &
             b% model_twins_flag, b% model_twins_flag_old, &
             b% dt_why_reason, b% dt_why_reason_old, &
             b% have_star_1, b% have_star_2, &
             b% CE_flag, b% CE_flag_old, &
             b% CE_init, b% CE_init_old, &
             b% CE_nz, b% CE_initial_radius, b% CE_initial_separation, b% CE_initial_Mdonor, &
             b% CE_initial_Maccretor, b% CE_initial_age, b% CE_initial_model_number, &
             b% CE_b_initial_age, b% CE_b_initial_model_number, &
             b% CE_num1, b% CE_num1_old, &
             b% CE_num2, b% CE_num2_old, &
             b% CE_lambda1, b% CE_lambda1_old, &
             b% CE_lambda2, b% CE_lambda2_old, &
             b% CE_Ebind1, b% CE_Ebind1_old, &
             b% CE_Ebind2, b% CE_Ebind2_old, &
             b% ixtra(:), b% ixtra_old(:), &
             b% xtra(:), b% xtra_old(:), &
             b% lxtra(:), b% lxtra_old(:)

         if (b% CE_flag .and. b% CE_init) then
            nz = b% CE_nz
            allocate(b% CE_m(nz), b% CE_entropy(4*nz), &
               b% CE_U_in(4*nz), b% CE_U_out(4*nz), b% CE_Omega_in(4*nz), b% CE_Omega_out(4*nz), stat=ierr)
            if (ierr /= 0) stop "error during allocation in binary_photo_read"
            read(iounit, iostat=ierr) &
               b% CE_m(:), b% CE_entropy(:), b% CE_U_in(:), b% CE_U_out(:), b% CE_Omega_in(:), b% CE_Omega_out(:)
         end if

         call b% other_binary_photo_read(binary_id, iounit, ierr)

         if (ierr /= 0) stop "error in binary_photo_read"
      end subroutine binary_photo_read

      end module binary_photos

