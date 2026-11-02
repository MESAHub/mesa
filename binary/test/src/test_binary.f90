module test_binary_mod

   use binary_lib
   use binary_timestep, only: binary_pick_next_timestep
   use binary_wind, only: Tout_enhance_wind
   use binary_def
   use binary_private_def
   use star_def
   use star_lib
   use const_def
   use math_lib
   use utils_lib, only: mesa_error
   implicit none

contains

   subroutine do_test
      integer :: binary_id, id, ierr, res, wind_id, wind_case
      real(dp) :: fL2, expected_wind
      type(binary_info), pointer :: b
      type(star_info), pointer :: s, wind_secondary, wind_star

      include 'formats'

      call math_init  ! we need this for exp10

      write (*, *) 'check time_delta_coeff behavior'
      binary_id = alloc_binary(ierr)
      call binary_ptr(binary_id, b, ierr)
      call alloc_star(id, ierr)
      call star_ptr(id, s, ierr)
      b%s_donor => s  ! set donor star
      b%d_i = 1
      b%a_i = 2
      b%point_mass_i = 1  ! no second star
      b%rl_relative_gap_old(2) = 0d0
      b%rl_relative_gap(2) = 0d0

      b%max_timestep = 10d0  ! current max timestep
      s%time_step = 10d0  ! current time step

      ! these are dummy vals so function doesn't give nans
      s%star_mass = 1d0
      s%he_core_mass = 1d0
      b%mtransfer_rate = 0d0
      b%env_old(b%d_i) = 0d0
      b%angular_momentum_j_old = 0d0
      b%separation_old = 0d0
      b%eccentricity_old = 0d0
      b%fa = -1d0
      b%fe = -1d0
      b%fj = -1d0
      b%fm = -1d0
      b%fdm = -1d0
      b%mtransfer_rate = 0d0
      b%s_donor%center_h1 = 0d0
      b%varcontrol_ms = 0d0
      b%varcontrol_post_ms = 0d0

      ! test for rl
      b%rl_relative_gap_old(b%d_i) = -0.102d0
      b%rl_relative_gap(b%d_i) = -0.100d0  ! so change = -2d-3
      b%fr = 1d-2
      b%fr_limit = 1d-2
      b%fr_dt_limit = 0d0
      b%fr_hard = -1d0
      b%ignore_hard_limits_this_step = .true.

      b%dt_softening_factor = 0d0  ! no softening
      b%time_delta_coeff = 1d0
      ! rel change is 2d-2 so timestep should be cut in half
      ! max(rel_change, fr_limit) / fr = 2d0 = cut factor
      res = binary_pick_next_timestep(b)

      write (*, 1) 'b% max_timestep / secyer', b%max_timestep/secyer  ! should be 5

      b%max_timestep = 10d0
      b%time_delta_coeff = 0.5d0
      ! now fr is effectively 0.5d-2, so timestep cut in 4 this time
      res = binary_pick_next_timestep(b)

      write (*, 1) 'b% max_timestep / secyer', b%max_timestep/secyer  ! should be 2.5

      ! Test accretion disk functions
      write (*, *) 'check L2_mass_loss_fraction behavior'

      fL2 = binary_L2_mass_loss_fraction(20.0_dp, 10.0_dp, 0.001_dp, 2.0_dp, 0.1_dp, 1.3_dp/2.4_dp, ierr)

      write (*, 1) 'fL2', fL2  ! should be ~0.9

      fL2 = binary_L2_mass_loss_fraction(20.0_dp, 10.0_dp, 0.0001_dp, 2.0_dp, 0.1_dp, 1.3_dp/2.4_dp, ierr)

      write (*, 1) 'fL2', fL2  ! should be ~0.2

      fL2 = binary_L2_mass_loss_fraction(20.0_dp, 18.0_dp, 0.0_dp, 304.0_dp, 0.1_dp, 0.64_dp, ierr)

      write (*, 1) 'fL2', fL2  ! should be 0.0

      ! Check component selection in the actual wind enhancement routine.
      call alloc_star(wind_id, ierr)
      call star_ptr(wind_id, wind_secondary, ierr)
      b% anomaly_steps = 2
      allocate(b% theta_co(2), b% time_co(2))
      b% theta_co = [0d0, pi]
      b% time_co = [0d0, 1d0]
      b% eccentricity = 0d0
      b% r = 1d0
      b% rl = 2d0
      b% do_enhance_wind_1 = .true.
      b% do_enhance_wind_2 = .true.
      b% tout_B_wind_1 = 128d0
      b% tout_B_wind_2 = 64d0
      do wind_case = 1, 4
         b% point_mass_i = 0
         if (wind_case == 3) b% point_mass_i = 1
         if (wind_case == 4) b% point_mass_i = 2
         b% have_star_1 = wind_case /= 3
         b% have_star_2 = wind_case /= 4
         b% star_ids = [id, wind_id]
         b% s1 => s
         b% s2 => wind_secondary
         if (.not. b% have_star_1) then
            ! An absent primary may retain the surviving star's ID.
            b% star_ids(1) = wind_id
            nullify(b% s1)
         end if
         if (.not. b% have_star_2) nullify(b% s2)
         wind_star => s
         if (wind_case == 2 .or. wind_case == 3) wind_star => wind_secondary
         expected_wind = merge(-3d0, -2d0, wind_case == 1 .or. wind_case == 4)
         wind_star% mstar_dot = -1d0
         call Tout_enhance_wind(b, wind_star)
         if (abs(wind_star% mstar_dot - expected_wind) > 1d-12) &
            call mesa_error(__FILE__, __LINE__, 'Incorrect wind enhancement component')
      end do
      deallocate(b% theta_co, b% time_co)
      nullify(b% theta_co, b% time_co, b% s1, b% s2, wind_star)
      call free_star(wind_id, ierr)
      nullify(wind_secondary)

      write (*, '(a)') 'done'

   end subroutine do_test

end module test_binary_mod

program test_binary
   use test_binary_mod
   implicit none
   call do_test
end program test_binary
