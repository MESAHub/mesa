! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
   
   include "test_suite_extras_def.inc"
   
   real(dp) :: constant_lnP, constant_lnT
   real(dp) :: flame_position, flame_width, flame_r0, flame_t0

contains

include "test_suite_extras.inc"
   
   subroutine build_isothermal_ball(id, ierr)
      
      use chem_def
      use chem_lib, only : basic_composition_info
      use eos_def
      use eos_lib, only : eosDT_get
      
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k, nz
      
      real(dp), allocatable, dimension(:) :: xa
      real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
      real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
      real(dp), dimension(:, :), allocatable :: d_dxa
      
      real(dp) :: rho_c, T_c
      real(dp) :: dmk, r0, dr, rp, rm
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      ! array for composition derivatives
      allocate(d_dxa(num_eos_d_dxa_results, s% species))
      
      ! choose a uniform composition
      allocate(xa(s% species))
      xa = 0
      xa(s% net_iso(ic12)) = s% x_ctrl(6)
      xa(s% net_iso(io16)) = s% x_ctrl(7)
      xa(s% net_iso(ine20)) = s% x_ctrl(8)
      
      call basic_composition_info(&
         s% species, s% chem_id, xa, &
         X, Y, Z, abar, zbar, z2bar, z53bar, ye, &
         mass_correction, sumx)
      
      ! get arrays for model
      nz = s% x_integer_ctrl(1)
      s% nz = nz
      call star_allocate_arrays(id, ierr)
      if (ierr /= 0) then
         return
      end if
      
      s% M_center = 0
      s% L_center = 0
      s% R_center = 0
      s% v_center = 0
      
      s% mstar = s% x_ctrl(1)
      s% star_mass = s% mstar / Msun
      s% xmstar = s% mstar
      
      Rho_c = s% x_ctrl(2)
      T_c = s% x_ctrl(3)
      
      ! initialize at constant density and temperature
      s% xh(s% i_lnd, 1:nz) = log(Rho_c)
      s% xh(s% i_lnT, 1:nz) = log(T_c)
      
      ! initialize uniform composition
      do k = 1, nz
         s% xa(1:s% species, k) = xa(1:s% species)
      end do
      
      ! call eos to get pressure
      call eosDT_get(&
         s% eos_handle, &
         s% species, s% chem_id, s% net_iso, xa, &
         Rho_c, log10(Rho_c), T_c, log10(T_c), &
         res, d_dlnd, d_dlnT, d_dxa, ierr)
      
      constant_lnT = log(T_c)
      constant_lnP = res(i_lnPgas)
      
      ! set radius
      r0 = pow(3d0 / (4d0 * pi) * s% mstar / rho_c, 1d0 / 3d0)
      dr = r0 / nz
      
      do k = 1, nz
         s% xh(s% i_lnR, k) = log((nz - k + 1) * dr)
      end do
      
      ! set mass coordinate
      s% dq(nz) = 4d0 * pi * pow(dr, 3d0) / 3d0 * rho_c / s% mstar
      s% q(nz) = s% dq(nz)
      do k = nz - 1, 1, -1
         rp = exp(s% xh(s% i_lnR, k))
         rm = exp(s% xh(s% i_lnR, k + 1))
         s% dq(k) = 4d0 * pi * (pow(rp, 3d0) - pow(rm, 3d0)) / 3d0 * rho_c / s% mstar
         s% q(k) = s% q(k + 1) + s% dq(k)
      end do
      
      ! set luminosity
      do k = 1, nz
         s% xh(s% i_lum, k) = 0
      end do
      
      ! bump the center
      do k = nz, 1, -1
         if (s% q(k) < s% x_ctrl(4)) then
            s% xh(s% i_lnT, k) = log(s% x_ctrl(5))
         end if
      end do
      
      deallocate(xa, d_dxa)
   
   end subroutine build_isothermal_ball
   
   
   subroutine constant_surface_PT(id, &
      skip_partials, &
      lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
      lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
      use const_def, only : dp
      integer, intent(in) :: id
      logical, intent(in) :: skip_partials
      real(dp), intent(out) :: &
         lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
      integer, intent(out) :: ierr
      
      type (star_info), pointer :: s
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      lnT_surf = constant_lnT
      dlnT_dL = 0
      dlnT_dlnR = 0
      dlnT_dlnM = 0
      dlnT_dlnkap = 0
      
      lnP_surf = constant_lnP
      dlnP_dL = 0
      dlnP_dlnR = 0
      dlnP_dlnM = 0
      dlnP_dlnkap = 0
   
   end subroutine constant_surface_PT
   
   
   subroutine extras_photo_read(id, iounit, ierr)
      integer, intent(in) :: id, iounit
      integer, intent(out) :: ierr
      ierr = 0
      read(iounit, iostat = ierr) constant_lnP, constant_lnT, &
         flame_r0, flame_t0
   end subroutine extras_photo_read
   
   
   subroutine extras_photo_write(id, iounit)
      integer, intent(in) :: id, iounit
      write(iounit) constant_lnP, constant_lnT, &
         flame_r0, flame_t0
   end subroutine extras_photo_write
   
   
   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      ! this is the place to set any procedure pointers you want to change
      ! e.g., other_wind, other_mixing, other_energy   (see star_data.inc)
      s% other_surface_PT => constant_surface_PT
      s% other_build_initial_model => build_isothermal_ball
      
      s% other_photo_read => extras_photo_read
      s% other_photo_write => extras_photo_write
      
      ! Uncomment these lines if you wish to use the functions in this file,
      ! otherwise we use a null_ version which does nothing.
      s% extras_startup => extras_startup
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
      
      ! Once you have set the function pointers you want,
      ! then uncomment this (or set it in your star_job inlist)
      ! to disable the printed warning message,
      s% job% warn_run_star_extras = .false.
   
   end subroutine extras_controls
   
   ! None of the following functions are called unless you set their
   ! function point in extras_control.
   
   
   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      if (.not. restart) then
         flame_r0 = -1
         flame_t0 = -1
      end if
      call test_suite_startup(s, restart, ierr)
   end subroutine extras_startup
   
   
   subroutine flame_properties(s, position, width)
      type (star_info), pointer :: s
      real(dp), intent(out) :: position, width
      integer :: k, k_burn, kmax
      real(dp) :: rtop, rbot
      
      kmax = maxloc(s% eps_nuc(1:s% nz), dim = 1)
      
      rtop = 0d0
      do k = kmax, 1, -1
         if (s% eps_nuc(k) .lt. 0.1d0 * maxval(s% eps_nuc(1:s% nz))) then
            rtop = s% r(k)
            exit
         end if
      end do
      
      rbot = 0d0
      do k = kmax, s% nz
         if (s% eps_nuc(k) .lt. 0.1d0 * maxval(s% eps_nuc(1:s% nz))) then
            rbot = s% r(k)
            exit
         end if
      end do
      
      position = s% r(kmax)
      width = max(0d0, rtop - rbot)
   
   end subroutine flame_properties
   
   
   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr, kmax
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going
      
      ! if you want to check multiple conditions, it can be useful
      ! to set a different termination code depending on which
      ! condition was triggered.   MESA provides 9 customizeable
      ! termination codes, named t_xtra1 .. t_xtra9.   You can
      ! customize the messages that will be printed upon exit by
      ! setting the corresponding termination_code_str value.
      ! termination_code_str(t_xtra1) = 'my termination condition'
      
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
      how_many_extra_history_columns = 2
   end function how_many_extra_history_columns
   
   
   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      ! note: do NOT add the extras names to history_columns.list
      ! the history_columns.list is only for the built-in log column options.
      ! it must not include the new column names you are adding here.
      
      call flame_properties(s, flame_position, flame_width)
      
      names(1) = 'r_flame'
      vals(1) = flame_position
      
      names(2) = 'l_flame'
      vals(2) = flame_width
   
   end subroutine data_for_extra_history_columns
   
   
   integer function how_many_extra_profile_columns(id)
      use star_def, only : star_info
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 0
   end function how_many_extra_profile_columns
   
   
   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      use star_def, only : star_info, maxlen_profile_column_name
      use const_def, only : dp
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      !note: do NOT add the extra names to profile_columns.list
      ! the profile_columns.list is only for the built-in profile column options.
      ! it must not include the new column names you are adding here.
      
      ! here is an example for adding a profile column
      !if (n /= 1) call mesa_error(__FILE__,__LINE__,'data_for_extra_profile_columns')
      !names(1) = 'beta'
      !do k = 1, nz
      !    vals(k,1) = s% Pgas(k)/s% P(k)
      !end do
   
   end subroutine data_for_extra_profile_columns
   
   
   integer function how_many_extra_history_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_header_items = 5
   end function how_many_extra_history_header_items
   
   
   subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if(ierr/=0) return
      
      names(1) = 'x_ctrl(1)'
      vals(1) = s% x_ctrl(1)
      
      names(2) = 'x_ctrl(2)'
      vals(2) = s% x_ctrl(2)
      
      names(3) = 'x_ctrl(3)'
      vals(3) = s% x_ctrl(3)
      
      names(4) = 'x_ctrl(4)'
      vals(4) = s% x_ctrl(4)
      
      names(5) = 'x_ctrl(5)'
      vals(5) = s% x_ctrl(5)
   
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
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if(ierr/=0) return
      
      ! here is an example for adding an extra profile header item
      ! also set how_many_extra_history_profile_items
      ! names(1) = 'mixing_length_alpha'
      ! vals(1) = s% mixing_length_alpha
   
   end subroutine data_for_extra_profile_header_items
   
   
   ! returns either keep_going or terminate.
   ! note: cannot request retry; extras_check_model can do that.
   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going
      
      ! to save a profile,
      ! s% need_to_save_profiles_now = .true.
      ! to update the star log,
      ! s% need_to_update_history_now = .true.
      
      ! get flame properties
      call flame_properties(s, flame_position, flame_width)
      
      ! when flame first is 30% through domain, record properties
      ! will be used to calculate flame speed
      if (flame_r0 .lt. 0) then
         if (flame_position .gt. 0.3d0 * s%r (1)) then
            flame_r0 = flame_position
            flame_t0 = s% star_age
         end if
      end if
      
      ! stop once flame is halfway through domain
      if (flame_position .gt. 0.5d0 * s%r (1)) then
         extras_finish_step = terminate
         s% termination_code = t_xtra1
         termination_code_str(t_xtra1) = 'flame reached halfway point'
         return
      end if
      
      ! see extras_check_model for information about custom termination codes
      ! by default, indicate where (in the code) MESA terminated
      if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
   end function extras_finish_step
   
   
   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: flame_speed, flame_speed_expected, flame_width_expected
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      include 'formats'
      
      call flame_properties(s, flame_position, flame_width)
      
      flame_speed = (flame_position - flame_r0) / ((s% star_age - flame_t0) * secyer)
      
      ! put target info in TestHub output
      testhub_extras_names(1) = 'flame_speed'; testhub_extras_vals(1) = flame_speed
      testhub_extras_names(2) = 'flame_width'; testhub_extras_vals(2) = flame_width
      
      write(*, '(A)')
      write(*, 1) testhub_extras_names(1), testhub_extras_vals(1)
      write(*, 1) testhub_extras_names(2), testhub_extras_vals(2)
      write(*, '(A)')
      
      ! get targets from inlist
      flame_speed_expected = s% x_ctrl(9)
      flame_width_expected = s% x_ctrl(10)
      
      if (abs(flame_speed - flame_speed_expected) > 0.1 * flame_speed_expected) then
         write(*, *) 'bad value for flame_speed'
         write(*, 1) 'flame_speed', flame_speed
         write(*, 1) 'expected', flame_speed_expected
         write(*, 1) 'flame_speed-expected', flame_speed - flame_speed_expected
      else if (abs(flame_width - flame_width_expected) > 0.1 * flame_width_expected) then
         write(*, *) 'bad value for flame_width'
         write(*, 1) 'flame_width', flame_width
         write(*, 1) 'expected', flame_width_expected
         write(*, 1) 'flame_width-expected', flame_width - flame_width_expected
      else
         write(*, '(a)') 'all values are within tolerance'
      end if
      
      call test_suite_after_evolve(s, ierr)
   end subroutine extras_after_evolve

end module run_star_extras
