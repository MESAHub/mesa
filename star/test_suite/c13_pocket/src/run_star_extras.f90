! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff

      implicit none

      real(dp) :: ms_t0, cheb_t0, ms_t1, cheb_t1, m_1DUP, mcore_TACHeB
      real(dp) :: mcore_1TP, age_1TP
      real(dp) :: mcore_at_TP, age_at_TP, mcore_min_after_TP
      real(dp) :: mcore_2TP_with_3DUP, age_2TP_with_3DUP
      integer :: TP_count, TP_with_3DUP
      logical :: in_LHe_peak
      real(dp) :: initial_surface_c12

      include "test_suite_extras_def.inc"

      contains

      include "test_suite_extras.inc"

      subroutine extras_photo_read(id, iounit, ierr)
        integer, intent(in) :: id, iounit
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        select case (s% x_integer_ctrl(1))
        case(1)
           read(iounit,iostat=ierr) ms_t0, cheb_t0, ms_t1, cheb_t1, m_1DUP, mcore_TACHeB
        case(2)
           read(iounit,iostat=ierr) mcore_at_TP, age_at_TP, mcore_min_after_TP, &
              mcore_1TP, age_1TP, TP_count, in_LHe_peak, &
              mcore_2TP_with_3DUP, age_2TP_with_3DUP, TP_with_3DUP
        case(3)
        case(4)
           read(iounit,iostat=ierr) TP_count, in_LHe_peak, initial_surface_c12
        end select

      end subroutine extras_photo_read

      subroutine extras_photo_write(id, iounit)
        integer, intent(in) :: id, iounit
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        select case (s% x_integer_ctrl(1))
        case(1)
           write(iounit) ms_t0, cheb_t0, ms_t1, cheb_t1, m_1DUP, mcore_TACHeB
        case(2)
           write(iounit) mcore_at_TP, age_at_TP, mcore_min_after_TP,&
              mcore_1TP, age_1TP, TP_count, in_LHe_peak, &
              mcore_2TP_with_3DUP, age_2TP_with_3DUP, TP_with_3DUP
        case(3)
        case(4)
           write(iounit) TP_count, in_LHe_peak, initial_surface_c12
        end select

      end subroutine extras_photo_write

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call test_suite_startup(s, restart, ierr)

         if (.not. restart) then
            select case(s% x_integer_ctrl(1))
            case(1)
               ms_t0 = 0
               ms_t1 = 0
               cheb_t0 = 0
               cheb_t1 = 0
               m_1DUP = 1d99
               mcore_TACHeB = 0
            case(2)
               mcore_1TP = 0
               age_1TP = 0
               mcore_at_TP = 0
               mcore_min_after_TP = 0
               TP_count = 0
               mcore_2TP_with_3DUP = 0
               age_2TP_with_3DUP = 0
               TP_with_3DUP = 0
               in_LHe_peak = .false.
            case(3)
            case(4)
               TP_count = 0
               in_LHe_peak = .false.
               initial_surface_c12 = s% surface_c12
            end select
         end if

      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) call mesa_error(__FILE__,__LINE__,'data_for_extra_profile_columns')
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         use chem_def, only: ic13
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         integer :: c13, k_max_c13

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         select case (s% x_integer_ctrl(1))
         case (1)

            ! measure MS lifetime
            if (ms_t0 .eq. 0) then
               if (s% power_h_burn .gt. 0.99 * s% L_surf) then
                  ms_t0 = s% star_age
                  write(*,*) 'started MS', ms_t0
               end if
            else
               if (ms_t1 .eq. 0) then
                  if (s% center_h1 .lt. 1d-4) then
                     ms_t1 = s% star_age
                     write(*,*) 'finished MS', ms_t1
                  end if
               end if
            end if

            ! measure CHeB lifetime
            if ((ms_t1 .ne. 0) .and. (cheb_t0 .eq. 0)) then
               if (s% power_he_burn .gt. 0.01 * s% L_surf) then
                  cheb_t0 = s% star_age
                  write(*,*) 'started CHeB', cheb_t0
               end if
            else
               if (cheb_t1 .eq. 0) then
                  if (s% center_he4 .lt. 1d-4) then
                     cheb_t1 = s% star_age
                     write(*,*) 'finished CHeB', cheb_t1
                     mcore_TACHeB = s% he_core_mass
                  end if
               end if
            end if

            ! measure extent of 1DUP
            if ((ms_t1 .gt. 0) .and. (cheb_t0 .eq. 0)) then
               if ((s% conv_mx1_top - s% conv_mx1_bot) .gt. 0.1) then
                  m_1DUP = min(s% conv_mx1_bot * s% star_mass, m_1DUP)
               end if
            end if

         case(2)

            ! record thermal pulses
            if (.not. in_LHe_peak) then
               ! check for peak
               if (s% power_he_burn .gt. 1e4) then
                  in_LHe_peak = .true.
                  TP_count = TP_count + 1
                  write(*,*) 'starting thermal pulse', TP_count
                  mcore_at_TP = s% he_core_mass
                  age_at_TP = s% star_age
                  mcore_min_after_TP = mcore_at_TP
                  if (TP_count == 1) then
                     mcore_1TP = s% he_core_mass
                     age_1TP = s% star_age
                  end if
                  if ((TP_with_3DUP > 0) .and. (TP_count - TP_with_3DUP == 1)) then
                     mcore_2TP_with_3DUP = s% he_core_mass
                     age_2TP_with_3DUP = s% star_age
                  end if
               end if
            else
               if (s% power_h_burn/s% power_he_burn .gt. 10) in_LHe_peak = .false. ! pulse over
            end if

            ! checking for 3DUP
            if (TP_count > 0) then
               mcore_min_after_TP = min(mcore_min_after_TP, s% he_core_mass)
            end if

            ! mark when signifcant 3DUP has first occured
            if (TP_with_3DUP == 0) then
               if ((mcore_min_after_TP - mcore_at_TP) < -1d-4) then
                  TP_with_3DUP = TP_count
                  write(*,'(A, I8)')    '>> 3DUP occurred at pulse: ', TP_with_3DUP
               end if
            end if

            ! stop after 3rd TP after dredge up starts
            if (.not. in_LHe_peak) then ! pulse is over
               if ((TP_with_3DUP > 0) .and. (TP_count - TP_with_3DUP == 2)) then
                  termination_code_str(t_xtra1) = 'third pulse with 3DUP has occurred'
                  s% termination_code = t_xtra1
                  extras_finish_step = terminate
               end if
            end if

         case(3)

         case(4)

            ! record thermal pulses
            if (.not. in_LHe_peak) then
               ! check for peak
               if (s% power_he_burn .gt. 1e4) then
                  in_LHe_peak = .true.
                  TP_count = TP_count + 1
                  write(*,*) 'starting thermal pulse'
                  mcore_min_after_TP = mcore_at_TP
               end if
            else
               if (s% power_h_burn/s% power_he_burn .gt. 10) in_LHe_peak = .false. ! pulse over
            end if

            ! stop after one TP
            if (.not. in_LHe_peak) then ! pulse is over
               if (TP_count == 1) then
                  termination_code_str(t_xtra1) = 'one thermal pulse cycle complete'
                  s% termination_code = t_xtra1
                  extras_finish_step = terminate
               end if
            end if

         end select

         ! Dynamic pgstar axis limits
         if (s% x_integer_ctrl(1) == 4) then
            c13 = s% net_iso(ic13)
            k_max_c13 = maxloc(s% xa(c13,1:s% nz),dim=1)
            if (s% xa(c13,k_max_c13) .gt. 0.01) then
               s% Abundance_xmin = s% m(k_max_c13)/Msun - 0.0001
               s% Abundance_xmax = s% m(k_max_c13)/Msun + 0.0001
            else
               s% Abundance_xmin = s% he_core_mass - 0.0125
               s% Abundance_xmax = s% he_core_mass + 0.0025
            end if
         end if

      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         use chem_def, only: ic13
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         integer :: c13, k_max_c13, k
         real(dp) :: max_c13, mass_max_c13, pocket_mass_c13
         real(dp) :: max_c13_expected, mass_max_c13_expected, pocket_mass_c13_expected
         real(dp) :: delta_surface_c12, delta_surface_c12_expected

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         include 'formats'

         write(*,'(A)')
         select case (s% x_integer_ctrl(1))
         case (1)
            write(*,'(A70, F8.3)') '[TestHub] Main-sequence lifetime (Gyr): ', (ms_t1 - ms_t0) / 1d9
            write(*,'(A70, F8.3)') '[TestHub] Deepest penetration of first dredge-up (Msun): ', m_1DUP
            write(*,'(A70, F8.3)') '[TestHub] H-free core mass at the end of He-core burning (Msun): ', mcore_TACHeB

            ! put comparison info in TestHub output
            testhub_extras_names(1) = 'MS_lifetime'; testhub_extras_vals(1) = (ms_t1 - ms_t0) / 1d9
            testhub_extras_names(2) = 'm_1DUP' ; testhub_extras_vals(2) = m_1DUP
            testhub_extras_names(3) = 'mcore_TACHeB'; testhub_extras_vals(3) = mcore_TACHeB

         case (2)
            write(*,'(A70, F8.3)') '[TestHub] Core mass at first thermal pulse (Msun): ', mcore_1TP
            write(*,'(A70, F8.3)') '[TestHub] Age at first thermal pulse (Gyr): ', age_1TP / 1d9
            write(*,'(A70, F8.3)') '[TestHub] Core mass at second thermal pulse with 3DUP (Msun): ', mcore_2TP_with_3DUP
            write(*,'(A70, F8.3)') '[TestHub] Following interpulse time (kyr): ', (age_at_TP - age_2TP_with_3DUP) / 1e3
            write(*,'(A70, F8.3)') '[TestHub] Following pulse-to-pulse core growth (1e-3 Msun): ', (mcore_at_TP - mcore_2TP_with_3DUP) / 1d-3
            write(*,'(A70, F8.3)') '[TestHub] Dredge up mass at following pulse (1e-3 Msun): ', (mcore_at_TP - mcore_min_after_TP) / 1d-3

            ! put comparison info in TestHub output
            testhub_extras_names(1) = 'mcore_1TP'; testhub_extras_vals(1) = mcore_1TP
            testhub_extras_names(2) = 'age_1TP' ; testhub_extras_vals(2) = age_1TP / 1d9
            testhub_extras_names(3) = 'mcore_2TP_with_3DUP'; testhub_extras_vals(3) = mcore_2TP_with_3DUP
            testhub_extras_names(4) = 'interpulse_time'; testhub_extras_vals(4) = (age_at_TP - age_2TP_with_3DUP) / 1e3
            testhub_extras_names(5) = 'delta_mcore_TP'; testhub_extras_vals(5) = (mcore_at_TP - mcore_2TP_with_3DUP) / 1d-3
            testhub_extras_names(6) = 'mass_DUP'; testhub_extras_vals(6) = (mcore_at_TP - mcore_min_after_TP) / 1d-3

         case(3)

         case(4)

            ! characterize c13 pocket location and mass
            c13 = s% net_iso(ic13)

            ! mass coordinate of peak c13
            k_max_c13 = maxloc(s% xa(c13,1:s% nz),dim=1)
            max_c13 = s% xa(c13,k_max_c13)
            mass_max_c13 = s% mstar*(s% q(k_max_c13) + s% dq(k_max_c13)/2)/Msun

            ! extent of region with significant c13
            pocket_mass_c13 = 0
            do k = 1, s% nz
               if (s% xa(c13,k) > 5d-3) then
                  pocket_mass_c13 = pocket_mass_c13 + s% dq(k)
               end if
            end do
            pocket_mass_c13 = pocket_mass_c13*s% star_mass ! mass in Msun units

            delta_surface_c12 = s% surface_c12 - initial_surface_c12

            write(*,1) 'max_c13', max_c13
            write(*,1) 'mass_max_c13', mass_max_c13
            write(*,1) 'pocket_mass_c13', pocket_mass_c13
            write(*,1) 'change in surface_c12', delta_surface_c12

            ! put target info in TestHub output
            testhub_extras_names(1) = 'max_c13'; testhub_extras_vals(1) = max_c13
            testhub_extras_names(2) = 'mass_max_c13' ; testhub_extras_vals(2) = mass_max_c13
            testhub_extras_names(3) = 'pocket_mass_c13'; testhub_extras_vals(3) = pocket_mass_c13
            testhub_extras_names(4) = 'delta_surface_c12'; testhub_extras_vals(4) = delta_surface_c12

            ! get targets from inlist
            max_c13_expected = s% x_ctrl(1)
            mass_max_c13_expected = s% x_ctrl(2)
            pocket_mass_c13_expected = s% x_ctrl(3)
            delta_surface_c12_expected = s% x_ctrl(4)

            if (abs(max_c13 - max_c13_expected) > 1d-2) then
               write(*,*) 'bad value for max_c13'
               write(*,1) 'max_c13', max_c13
               write(*,1) 'expected', max_c13_expected
               write(*,1) 'max_c13-expected', max_c13-max_c13_expected
            else if (abs(mass_max_c13 - mass_max_c13_expected) > 1d-2) then
               write(*,*) 'bad value for mass_max_c13'
               write(*,1) 'mass_max_c13', mass_max_c13
               write(*,1) 'expected', mass_max_c13_expected
               write(*,1) 'mass_max_c13-expected', mass_max_c13-mass_max_c13_expected
            else if (abs(pocket_mass_c13 - pocket_mass_c13_expected) > 1d-5) then
               write(*,*) 'bad value for pocket_mass_c13'
               write(*,1) 'pocket_mass_c13', pocket_mass_c13
               write(*,1) 'expected', pocket_mass_c13_expected
               write(*,1) 'pocket_mass_c13-expected', pocket_mass_c13-pocket_mass_c13_expected
            else if (abs(delta_surface_c12 - delta_surface_c12_expected) > 2d-4) then
               write(*,*) 'bad value for delta_surface_c12'
               write(*,1) 'delta_surface_c12', delta_surface_c12
               write(*,1) 'expected', delta_surface_c12_expected
               write(*,1) 'delta_surface_c12-expected', delta_surface_c12-delta_surface_c12_expected
            else
               write(*,'(a)') 'all values are within tolerance'
            end if

         end select
         write(*,'(A)')

         call test_suite_after_evolve(s, ierr)

      end subroutine extras_after_evolve


      end module run_star_extras
