! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module eos_timing

   use const_def, only: dp
   use eos_def, only: num_eos_frac_results

   implicit none

   private

   integer, parameter, public :: num_skye_dxa_timing_parts = 3
   integer, parameter, public :: i_skye_dxa_ideal = 1
   integer, parameter, public :: i_skye_dxa_coul = 2
   integer, parameter, public :: i_skye_dxa_pack = 3

   public :: eos_timing_set_dxa_enabled
   public :: eos_timing_reset_dxa
   public :: eos_timing_start
   public :: eos_timing_record_component
   public :: eos_timing_record_table_eval
   public :: eos_timing_record_table_expand
   public :: eos_timing_record_moments
   public :: eos_timing_record_blend
   public :: eos_timing_record_check
   public :: eos_timing_record_skye
   public :: eos_timing_get_dxa

   ! Coarse thread-time buckets for expensive EOS composition-partial calls;
   ! these are diagnostics only and do not change the solver path.
   logical :: dxa_timing_enabled = .false.

   real(dp) :: dxa_component_time(num_eos_frac_results) = 0d0
   real(dp) :: dxa_table_eval_time(num_eos_frac_results) = 0d0
   real(dp) :: dxa_table_expand_time(num_eos_frac_results) = 0d0
   real(dp) :: dxa_skye_time(num_skye_dxa_timing_parts) = 0d0
   real(dp) :: dxa_moments_time = 0d0
   real(dp) :: dxa_blend_time = 0d0
   real(dp) :: dxa_check_time = 0d0

   integer :: dxa_component_count(num_eos_frac_results) = 0
   integer :: dxa_table_eval_count(num_eos_frac_results) = 0
   integer :: dxa_table_expand_count(num_eos_frac_results) = 0
   integer :: dxa_skye_count(num_skye_dxa_timing_parts) = 0
   integer :: dxa_moments_count = 0
   integer :: dxa_blend_count = 0
   integer :: dxa_check_count = 0

contains

   subroutine eos_timing_set_dxa_enabled(enabled)
      logical, intent(in) :: enabled
      dxa_timing_enabled = enabled
   end subroutine eos_timing_set_dxa_enabled


   subroutine eos_timing_reset_dxa
      dxa_timing_enabled = .false.
      dxa_component_time = 0d0
      dxa_table_eval_time = 0d0
      dxa_table_expand_time = 0d0
      dxa_skye_time = 0d0
      dxa_moments_time = 0d0
      dxa_blend_time = 0d0
      dxa_check_time = 0d0
      dxa_component_count = 0
      dxa_table_eval_count = 0
      dxa_table_expand_count = 0
      dxa_skye_count = 0
      dxa_moments_count = 0
      dxa_blend_count = 0
      dxa_check_count = 0
   end subroutine eos_timing_reset_dxa


   subroutine eos_timing_start(time0, clock_rate)
      integer, intent(out) :: time0, clock_rate

      time0 = 0
      clock_rate = 1
      if (.not. dxa_timing_enabled) return
      call system_clock(time0, clock_rate)
   end subroutine eos_timing_start


   subroutine eos_timing_record_component(which_eos, time0, clock_rate)
      integer, intent(in) :: which_eos, time0, clock_rate
      call record_array(dxa_component_time, dxa_component_count, which_eos, time0, clock_rate)
   end subroutine eos_timing_record_component


   subroutine eos_timing_record_table_eval(which_eos, time0, clock_rate)
      integer, intent(in) :: which_eos, time0, clock_rate
      call record_array(dxa_table_eval_time, dxa_table_eval_count, which_eos, time0, clock_rate)
   end subroutine eos_timing_record_table_eval


   subroutine eos_timing_record_table_expand(which_eos, time0, clock_rate)
      integer, intent(in) :: which_eos, time0, clock_rate
      call record_array(dxa_table_expand_time, dxa_table_expand_count, which_eos, time0, clock_rate)
   end subroutine eos_timing_record_table_expand


   subroutine eos_timing_record_skye(which_part, time0, clock_rate)
      integer, intent(in) :: which_part, time0, clock_rate
      call record_array(dxa_skye_time, dxa_skye_count, which_part, time0, clock_rate)
   end subroutine eos_timing_record_skye


   subroutine eos_timing_record_moments(time0, clock_rate)
      integer, intent(in) :: time0, clock_rate
      call record_scalar(dxa_moments_time, dxa_moments_count, time0, clock_rate)
   end subroutine eos_timing_record_moments


   subroutine eos_timing_record_blend(time0, clock_rate)
      integer, intent(in) :: time0, clock_rate
      call record_scalar(dxa_blend_time, dxa_blend_count, time0, clock_rate)
   end subroutine eos_timing_record_blend


   subroutine eos_timing_record_check(time0, clock_rate)
      integer, intent(in) :: time0, clock_rate
      call record_scalar(dxa_check_time, dxa_check_count, time0, clock_rate)
   end subroutine eos_timing_record_check


   subroutine eos_timing_get_dxa( &
         component_time, component_count, table_eval_time, table_eval_count, &
         table_expand_time, table_expand_count, skye_time, skye_count, &
         moments_time, moments_count, blend_time, blend_count, &
         check_time, check_count)
      real(dp), intent(out) :: component_time(num_eos_frac_results)
      real(dp), intent(out) :: table_eval_time(num_eos_frac_results)
      real(dp), intent(out) :: table_expand_time(num_eos_frac_results)
      real(dp), intent(out) :: skye_time(num_skye_dxa_timing_parts)
      real(dp), intent(out) :: moments_time, blend_time, check_time
      integer, intent(out) :: component_count(num_eos_frac_results)
      integer, intent(out) :: table_eval_count(num_eos_frac_results)
      integer, intent(out) :: table_expand_count(num_eos_frac_results)
      integer, intent(out) :: skye_count(num_skye_dxa_timing_parts)
      integer, intent(out) :: moments_count, blend_count, check_count

      component_time = dxa_component_time
      table_eval_time = dxa_table_eval_time
      table_expand_time = dxa_table_expand_time
      skye_time = dxa_skye_time
      moments_time = dxa_moments_time
      blend_time = dxa_blend_time
      check_time = dxa_check_time

      component_count = dxa_component_count
      table_eval_count = dxa_table_eval_count
      table_expand_count = dxa_table_expand_count
      skye_count = dxa_skye_count
      moments_count = dxa_moments_count
      blend_count = dxa_blend_count
      check_count = dxa_check_count
   end subroutine eos_timing_get_dxa


   subroutine record_scalar(time_bucket, count_bucket, time0, clock_rate)
      real(dp), intent(inout) :: time_bucket
      integer, intent(inout) :: count_bucket
      integer, intent(in) :: time0, clock_rate

      integer :: time1
      real(dp) :: dt

      if (.not. dxa_timing_enabled) return
      if (time0 == 0 .or. clock_rate <= 0) return

      call system_clock(time1)
      dt = dble(time1 - time0)/dble(clock_rate)
!$OMP critical (eos_timing_critical)
      time_bucket = time_bucket + dt
      count_bucket = count_bucket + 1
!$OMP end critical (eos_timing_critical)
   end subroutine record_scalar


   subroutine record_array(time_bucket, count_bucket, index, time0, clock_rate)
      real(dp), intent(inout) :: time_bucket(:)
      integer, intent(inout) :: count_bucket(:)
      integer, intent(in) :: index, time0, clock_rate

      integer :: time1
      real(dp) :: dt

      if (.not. dxa_timing_enabled) return
      if (time0 == 0 .or. clock_rate <= 0) return
      if (index < 1 .or. index > size(time_bucket)) return

      call system_clock(time1)
      dt = dble(time1 - time0)/dble(clock_rate)
!$OMP critical (eos_timing_critical)
      time_bucket(index) = time_bucket(index) + dt
      count_bucket(index) = count_bucket(index) + 1
!$OMP end critical (eos_timing_critical)
   end subroutine record_array

end module eos_timing
