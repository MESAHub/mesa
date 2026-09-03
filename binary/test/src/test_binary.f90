module test_binary_mod

   use binary_lib
   use binary_timestep, only: binary_pick_next_timestep
   use binary_roche_deformation, only: build_roche_interpolators, eval_psi
   use binary_def
   use binary_private_def
   use star_def
   use star_lib
   use const_lib, only: const_init
   use const_def
   use math_lib
   use utils_lib, only: mesa_error
   use ieee_arithmetic, only: ieee_is_finite
   implicit none

contains

   subroutine do_test
      integer :: binary_id, id, ierr, res
      real(dp) :: fL2
      type(binary_info), pointer :: b
      type(star_info), pointer :: s
      character(len=32) :: my_mesa_dir

      include 'formats'

      my_mesa_dir = '../..'
      call const_init(my_mesa_dir, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
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

      call test_roche_psi

      write (*, '(a)') 'done'

   end subroutine do_test

   subroutine test_roche_psi
      integer :: i
      real(dp) :: psi_inner, psi_cutoff, psi_below_cutoff
      real(dp), parameter :: u_cutoff = 0.0193458d0
      real(dp), parameter :: du = 1d-5 * u_cutoff
      real(dp), parameter :: lq_values(5) = [-2.9d0, -1.5d0, 0d0, 1.5d0, 2.9d0]

      call build_roche_interpolators

      do i = 1, size(lq_values)
         psi_inner = eval_finite_psi(lq_values(i), 0.5d0 * u_cutoff)
         psi_cutoff = eval_finite_psi(lq_values(i), u_cutoff)
         psi_below_cutoff = eval_finite_psi(lq_values(i), u_cutoff - du)

         if (psi_inner >= psi_below_cutoff) call mesa_error(__FILE__, __LINE__)
         if (abs(psi_cutoff - psi_below_cutoff) > &
               1d-3 * max(1d0, abs(psi_cutoff), abs(psi_below_cutoff))) &
            call mesa_error(__FILE__, __LINE__)
      end do
   end subroutine test_roche_psi

   real(dp) function eval_finite_psi(lq, u) result(psi)
      real(dp), intent(in) :: lq, u
      integer :: ierr

      psi = eval_psi(lq, u, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (.not. ieee_is_finite(psi)) call mesa_error(__FILE__, __LINE__)
   end function eval_finite_psi

end module test_binary_mod

program test_binary
   use test_binary_mod
   implicit none
   call do_test
end program test_binary
