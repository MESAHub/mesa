! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      implicit none

      include 'test_suite_extras_def.inc'

      real(dp) :: e_res_err_step, e_res_err_run, e_eos_err_step, e_eos_err_run, &
         e_eos_err_step_blend, e_eos_err_run_blend

      contains

      include 'test_suite_extras.inc'


      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit, iostat=ierr) e_res_err_step, e_res_err_run, &
            e_eos_err_step, e_eos_err_run, &
            e_eos_err_step_blend, e_eos_err_run_blend
      end subroutine extras_photo_read


      subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) e_res_err_step, e_res_err_run, &
            e_eos_err_step, e_eos_err_run, &
            e_eos_err_step_blend, e_eos_err_run_blend
      end subroutine extras_photo_write


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (.not. restart) then
            e_res_err_step = 0
            e_res_err_run = 0
            e_eos_err_step = 0
            e_eos_err_run = 0
            e_eos_err_step_blend = 0
            e_eos_err_run_blend = 0
         end if

         call test_suite_startup(s, restart, ierr)

      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


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
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 6
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

         names(1) = 'rel_e_res_err_step'
         vals(1) = e_res_err_step / s% total_energy_end

         names(2) = 'rel_e_res_err_run'
         vals(2) = e_res_err_run / s% total_energy_end

         names(3) = 'rel_e_eos_err_step'
         vals(3) = e_eos_err_step / s% total_energy_end

         names(4) = 'rel_e_eos_err_run'
         vals(4) = e_eos_err_run / s% total_energy_end

         names(5) = 'rel_e_eos_err_step_blend'
         vals(5) = e_eos_err_step_blend / s% total_energy_end

         names(6) = 'rel_e_eos_err_run_blend'
         vals(6) = e_eos_err_run_blend / s% total_energy_end

         
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 4
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'eos_frac_HELM'
         names(2) = 'eos_frac_FreeEOS'
         names(3) = 'eps_eos'
         names(4) = 'rse_eps_grav'
         do k = 1, nz
            vals(k,1) = s% eos_frac_HELM(k)
            vals(k,2) = s% eos_frac_FreeEOS(k)
            vals(k,3) = ((s% energy(k) - s% energy_start(k)) + &
               -0.5d0 * (s% T_start(k)*s% cv_start(k) + s% T(k)*s% cv(k)) * s% dxh_lnT(k) + &
               -0.5d0 * (s% Rho_start(k)*s% dE_dRho_start(k) + s% Rho(k)*s% dE_dRho(k)) * s% dxh_lnd(k)) * s% dt
            vals(k,4) = (&
               -0.5d0 * (s% T_start(k)*s% cv_start(k) + s% T(k)*s% cv(k)) * s% dxh_lnT(k) + &
               -0.5d0 * (s% Rho_start(k)*s% dE_dRho_start(k) + s% Rho(k)*s% dE_dRho(k)) * s% dxh_lnd(k)) * s% dt
         end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr, k, nz
         real(dp) :: e_eos_err_cell
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         nz = s% nz

         ! total energy residual
         e_res_err_step = sum(s% ergs_error(1:nz))

         ! difference between dedt and eps_grav - style estimate
         e_eos_err_step = 0
         e_eos_err_step_blend = 0
         do k = 1, nz
            e_eos_err_cell = s% dm(k) * &
               ((s% energy(k) - s% energy_start(k)) + &
               -0.5d0 * (s% T_start(k)*s% cv_start(k) + s% T(k)*s% cv(k)) * s% dxh_lnT(k) + &
               -0.5d0 * (s% Rho_start(k)*s% dE_dRho_start(k) + s% Rho(k)*s% dE_dRho(k)) * s% dxh_lnd(k))
            e_eos_err_step = e_eos_err_step + e_eos_err_cell
            if (in_eos_blend(s,k)) e_eos_err_step_blend = e_eos_err_step_blend + e_eos_err_cell
         end do

         ! accumulate cumulative quantities
         e_res_err_run = e_res_err_run + e_res_err_step
         e_eos_err_run = e_eos_err_run + e_eos_err_step
         e_eos_err_run_blend = e_eos_err_run_blend + e_eos_err_step_blend

      end function extras_finish_step


      logical function in_eos_blend(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k

         in_eos_blend = &
            ((s% eos_frac_OPAL_SCVH(k) .gt. 0) .and. (s% eos_frac_OPAL_SCVH(k) .lt. 1)) .or. &
            ((s% eos_frac_HELM(k) .gt. 0) .and. (s% eos_frac_HELM(k) .lt. 1)) .or. &
            ((s% eos_frac_Skye(k) .gt. 0) .and. (s% eos_frac_Skye(k) .lt. 1)) .or. &
            ((s% eos_frac_PC(k) .gt. 0) .and. (s% eos_frac_PC(k) .lt. 1)) .or. &
            ((s% eos_frac_CMS(k) .gt. 0) .and. (s% eos_frac_CMS(k) .lt. 1)) .or. &
            ((s% eos_frac_FreeEOS(k) .gt. 0) .and. (s% eos_frac_FreeEOS(k) .lt. 1))
         
      end function in_eos_blend
      

      end module run_star_extras
