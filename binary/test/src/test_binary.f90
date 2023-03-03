module test_binary_mod
   use binary_lib
   use binary_timestep, only: binary_pick_next_timestep
   use binary_def
   use binary_private_def
   use star_def
   use star_lib
   use const_def
   use math_lib
   implicit none

contains

   subroutine do_test
      integer :: binary_id, id, ierr, res
      type (binary_info), pointer :: b
      type (star_info), pointer :: s

      include 'formats'

      call math_init  ! we need this for exp10

      write(*, *) 'check time_delta_coeff behavior'
      binary_id = alloc_binary(ierr)
      call binary_ptr(binary_id, b, ierr)
      call alloc_star(id, ierr)
      call star_ptr(id, s, ierr)
      b% s_donor => s  ! set donor star
      b% d_i = 1
      b% point_mass_i = 1  ! no second star
      b% rl_relative_gap_old(2) = 0d0
      b% rl_relative_gap(2) = 0d0

      b% max_timestep = 10d0  ! current max timestep
      s% time_step = 10d0  ! current time step

      ! these are dummy vals so function doesn't give nans
      s% star_mass = 1d0
      s% he_core_mass = 1d0
      b% env_old(b% d_i) = 0d0
      b% angular_momentum_j_old = 0d0
      b% separation_old = 0d0
      b% eccentricity_old = 0d0
      b% fa = -1d0
      b% fe = -1d0
      b% fj = -1d0
      b% fm = -1d0
      b% fdm = -1d0

      ! test for rl
      b% rl_relative_gap_old(b% d_i) = -0.102d0
      b% rl_relative_gap(b% d_i) = -0.100d0  ! so change = -2d-3
      b% fr = 1d-2
      b% fr_limit = 1d-2
      b% fr_dt_limit = 0d0
      b% fr_hard = -1d0
      b% ignore_hard_limits_this_step = .true.

      b% dt_softening_factor = 0d0  ! no softening
      b% time_delta_coeff = 1d0
      ! rel change is 2d-2 so timestep should be cut in half
      ! max(rel_change, fr_limit) / fr = 2d0 = cut factor
      res = binary_pick_next_timestep(b)

      write(*,1) 'b% max_timestep / secyer', b% max_timestep / secyer  ! should be 5

      b% max_timestep = 10d0
      b% time_delta_coeff = 0.5d0
      ! now fr is effectively 0.5d-2, so timestep cut in 4 this time
      res = binary_pick_next_timestep(b)

      write(*,1) 'b% max_timestep / secyer', b% max_timestep / secyer  ! should be 2.5

      write(*,'(a)') 'done'

   end subroutine do_test


end module test_binary_mod


program test_binary
   use test_binary_mod
   implicit none
   call do_test
end program
