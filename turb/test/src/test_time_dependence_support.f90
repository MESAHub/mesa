module test_time_dependence_support

   use math_lib
   use auto_diff
   use const_def, only: dp, pi, rsun, lsun, msun, kerg, mp, boltz_sigma, standard_cgrav, sqrt_2_div_3
   use turb, only: set_MLT, set_TDC, TDC_arnett_growth_target_mlt, &
      TDC_arnett_growth_target_tdc_no_mlt_corr, TDC_arnett_growth_target_tdc_with_mlt_corr

   implicit none

   type :: tdc_mode_data
      character(len=48) :: name
      logical :: use_arnett
      logical :: use_acceleration_limit
      logical :: use_Af_split
      logical :: include_mlt_correction
      integer :: growth_target
   end type tdc_mode_data

   type :: tdc_case_data
      real(dp) :: mixing_length_alpha, TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt
      real(dp) :: cgrav, m, scale, max_conv_vel, TDC_alpha_C, TDC_alpha_S
      real(dp) :: conv_vel_start_ref
      logical :: report, use_TDC_enthalpy_flux_limiter
      type(auto_diff_real_star_order1) :: &
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, &
         gradr, grada, scale_height, gradL, Eq_div_w, grav, energy
   end type tdc_case_data

   integer, parameter :: num_tdc_modes = 7
   type(tdc_mode_data), parameter :: tdc_modes(num_tdc_modes) = [ &
      tdc_mode_data('tdc', .false., .false., .false., .false., TDC_arnett_growth_target_mlt), &
      tdc_mode_data('tdc_with_Af_split', .false., .false., .true., .false., TDC_arnett_growth_target_mlt), &
      tdc_mode_data('tdc_with_mlt_corr', .false., .false., .false., .true., TDC_arnett_growth_target_mlt), &
      tdc_mode_data('tdc_with_arnett_closure', .true., .false., .false., .true., TDC_arnett_growth_target_mlt), &
      tdc_mode_data('tdc_with_arnett_closure_tdc_ss', .true., .false., .false., .false., TDC_arnett_growth_target_tdc_no_mlt_corr), &
      tdc_mode_data('tdc_with_arnett_closure_tdc_ss_mlt_corr', .true., .false., .false., .true., TDC_arnett_growth_target_tdc_with_mlt_corr), &
      tdc_mode_data('tdc_with_acceleration_limit', .false., .true., .false., .false., TDC_arnett_growth_target_mlt) ]

   private
   public :: write_time_dependence_csv
   public :: check_time_dependence

contains

   subroutine write_time_dependence_csv(output_path, ierr)
      character(len=*), intent(in) :: output_path
      integer, intent(out) :: ierr
      integer :: io, mode_i

      ierr = 0
      open(newunit=io, file=trim(output_path), status='replace', action='write', iostat=ierr)
      if (ierr /= 0) return

      write(io, '(a)') 'mode,group,scenario,step,time,dt,conv_vel_start,L_conv_start,L,gradr,gradL,grada,gradT,Y_face,conv_vel,L_conv_end,L_conv_div_L,tau_phys,ierr'

      do mode_i = 1, num_tdc_modes
         call run_mode(io, tdc_modes(mode_i))
         call run_cburn_growth_evolution(io, tdc_modes(mode_i))
      end do

      call run_mlt_evolution(io)
      close(io)
   end subroutine write_time_dependence_csv


   subroutine check_time_dependence()
      integer :: mode_i
      type(tdc_case_data) :: decay_case, growth_case, envelope_decay_case, cburn_case
      real(dp), parameter :: decay_gradT_start = 2.5204370043250246d-01
      real(dp), parameter :: cburn_gradT_start = 2.8060983127196415d-01
      real(dp) :: envelope_decay_gradT_start
      write(*, '(a)') ''
      write(*, '(a)') 'Decay'
      write(*, '(a)') ''
      call setup_decay_case(decay_case)
      call print_case_setup(decay_case, decay_case%conv_vel_start_ref, decay_gradT_start, &
         'Idealized deep-interior decay test with a large inward luminosity; a generic strongly contracting negative-luminosity state.')

      do mode_i = 1, num_tdc_modes
         call print_decay_track(tdc_modes(mode_i))
      end do

      write(*, '(a)') ''
      write(*, '(a)') 'Growth'
      write(*, '(a)') ''
      call setup_growth_case(growth_case)
      call print_case_setup(growth_case, 0d0, growth_case%gradr%val, &
         'Idealized generic onset test in a dilute, marginally unstable layer initialized on the radiative branch.')

      do mode_i = 1, num_tdc_modes
         call print_growth_track(tdc_modes(mode_i))
      end do

      write(*, '(a)') ''
      write(*, '(a)') 'Envelope Decay'
      write(*, '(a)') ''
      call setup_envelope_decay_case(envelope_decay_case, envelope_decay_gradT_start)
      call print_case_setup(envelope_decay_case, envelope_decay_case%conv_vel_start_ref, envelope_decay_gradT_start, &
         'Idealized generic envelope decay test in a dilute, slightly stable layer initialized from a nearby unstable convective state.')

      do mode_i = 1, num_tdc_modes
         call print_envelope_decay_track(tdc_modes(mode_i))
      end do

      write(*, '(a)') ''
      write(*, '(a)') 'Carbon Burning Growth'
      write(*, '(a)') ''
      call setup_cburn_case(cburn_case)
      call print_case_setup(cburn_case, cburn_case%conv_vel_start_ref, cburn_gradT_start, &
         '12 Msun pre-core-collapse model at core-carbon depletion; a hot convective carbon-burning zone taken directly from the saved model.')

      do mode_i = 1, num_tdc_modes
         call print_cburn_growth_track(tdc_modes(mode_i), mode_i < num_tdc_modes)
      end do
   end subroutine check_time_dependence


   subroutine run_mode(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode

      call run_decay_memory_sweep(io, mode)
      call run_decay_velocity_sweep(io, mode)
      call run_growth_sweep(io, mode)
      call run_growth_flux_memory_sweep(io, mode)
      call run_growth_velocity_sweep(io, mode)
      call run_decay_evolution(io, mode)
      call run_growth_evolution(io, mode)
      call run_envelope_decay_evolution(io, mode)
   end subroutine run_mode


   subroutine print_decay_track(mode)
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: decay_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time
      integer :: ierr, j
      real(dp), parameter :: gradT_start_old = 2.5204370043250246d-01
      integer, parameter :: num_steps = 12

      call setup_decay_case(decay_case)
      conv_vel_start = decay_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT_start_old)
      tau_phys = arnett_relaxation_time(decay_case%mixing_length_alpha * decay_case%scale_height%val, decay_case%grav%val, &
         decay_case%scale_height%val, decay_case%gradL%val, decay_case%gradr%val, conv_vel_start, conv_vel_start)
      dt = pow(10d0, 2.3d0)

      write(*, '(a)') 'Convection model: ' // trim(mode%name)
      call print_mode_controls(mode)
      write(*, '(a)') '      step                 time[s]                   dt[s]                   gradT' // &
         '                  Y_face      conv_vel_start[cm/s]           conv_vel[cm/s]                L_conv/L'

      gradT = gradT_start_old
      Y_face = gradT_start_old - decay_case%gradL
      conv_vel = conv_vel_start
      call print_track_row(0, 0d0, dt, conv_vel_start, decay_case, gradT, Y_face, conv_vel)

      do j = 1, num_steps
         call call_tdc(decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call print_track_row(j, time, dt, conv_vel_start, decay_case, gradT, Y_face, conv_vel)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT%val)
      end do

      write(*, '(a,1x,es24.16)') '      tau_phys[s]', tau_phys
      write(*, '(a)') ''
   end subroutine print_decay_track


   subroutine print_growth_track(mode)
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: growth_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys
      real(dp) :: time
      integer :: ierr, j
      integer, parameter :: num_steps = 12

      call setup_growth_case(growth_case)
      call call_tdc(growth_case, mode, 0d0, 0d0, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      if (ierr /= 0) return

      tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
         growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, 0d0, eq_conv_vel%val)
      call reference_evolution_dt(growth_case, 0d0, 0d0, dt, ierr)
      if (ierr /= 0) return
      conv_vel_start = 0d0
      L_conv_start = 0d0

      write(*, '(a)') 'Convection model: ' // trim(mode%name)
      call print_mode_controls(mode)
      write(*, '(a)') '      step                 time[s]                   dt[s]                   gradT' // &
         '                  Y_face      conv_vel_start[cm/s]           conv_vel[cm/s]           conv_vel/conv_vel_eq                L_conv/L'

      gradT = growth_case%gradr
      Y_face = growth_case%gradr - growth_case%gradL
      conv_vel = 0d0
      call print_track_row(0, 0d0, dt, conv_vel_start, growth_case, gradT, Y_face, conv_vel, eq_conv_vel%val)

      do j = 1, num_steps
         call call_tdc(growth_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call print_track_row(j, time, dt, conv_vel_start, growth_case, gradT, Y_face, conv_vel, eq_conv_vel%val)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(growth_case%L%val, growth_case%gradr%val, gradT%val)
      end do

      write(*, '(a,1x,es24.16)') '      tau_phys[s]', tau_phys
      write(*, '(a,1x,es24.16)') '      conv_vel_eq[cm/s]', eq_conv_vel%val
      write(*, '(a)') ''
   end subroutine print_growth_track


   subroutine print_envelope_decay_track(mode)
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: envelope_decay_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time, gradT_start
      integer :: ierr, j
      integer, parameter :: num_steps = 12

      call setup_envelope_decay_case(envelope_decay_case, gradT_start)
      conv_vel_start = envelope_decay_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(envelope_decay_case%L%val, envelope_decay_case%gradr%val, gradT_start)
      tau_phys = arnett_relaxation_time(envelope_decay_case%mixing_length_alpha * envelope_decay_case%scale_height%val, &
         envelope_decay_case%grav%val, envelope_decay_case%scale_height%val, envelope_decay_case%gradL%val, &
         envelope_decay_case%gradr%val, conv_vel_start, 0d0)
      call reference_evolution_dt(envelope_decay_case, conv_vel_start, L_conv_start, dt, ierr)
      if (ierr /= 0) return

      write(*, '(a)') 'Convection model: ' // trim(mode%name)
      call print_mode_controls(mode)
      write(*, '(a)') '      step                 time[s]                   dt[s]                   gradT' // &
         '                  Y_face      conv_vel_start[cm/s]           conv_vel[cm/s]                L_conv/L'

      gradT = gradT_start
      Y_face = gradT_start - envelope_decay_case%gradL
      conv_vel = conv_vel_start
      call print_track_row(0, 0d0, dt, conv_vel_start, envelope_decay_case, gradT, Y_face, conv_vel)

      do j = 1, num_steps
         call call_tdc(envelope_decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call print_track_row(j, time, dt, conv_vel_start, envelope_decay_case, gradT, Y_face, conv_vel)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(envelope_decay_case%L%val, envelope_decay_case%gradr%val, gradT%val)
      end do

      write(*, '(a,1x,es24.16)') '      tau_phys[s]', tau_phys
      write(*, '(a)') ''
   end subroutine print_envelope_decay_track


   subroutine print_cburn_growth_track(mode, emit_blank_line)
      type(tdc_mode_data), intent(in) :: mode
      logical, intent(in), optional :: emit_blank_line
      type(tdc_case_data) :: cburn_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time
      logical :: do_blank_line
      integer :: ierr, j
      real(dp), parameter :: gradT_start = 2.8060983127196415d-01
      integer, parameter :: num_steps = 12

      do_blank_line = .true.
      if (present(emit_blank_line)) do_blank_line = emit_blank_line

      call setup_cburn_case(cburn_case)
      conv_vel_start = cburn_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(cburn_case%L%val, cburn_case%gradr%val, gradT_start)
      call call_tdc(cburn_case, mode, conv_vel_start, L_conv_start, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      if (ierr /= 0) return

      tau_phys = arnett_relaxation_time(cburn_case%mixing_length_alpha * cburn_case%scale_height%val, cburn_case%grav%val, &
         cburn_case%scale_height%val, cburn_case%gradL%val, cburn_case%gradr%val, conv_vel_start, eq_conv_vel%val)
      call reference_evolution_dt(cburn_case, conv_vel_start, L_conv_start, dt, ierr)
      if (ierr /= 0) return

      write(*, '(a)') 'Convection model: ' // trim(mode%name)
      call print_mode_controls(mode)
      write(*, '(a)') '      step                 time[s]                   dt[s]                   gradT' // &
         '                  Y_face      conv_vel_start[cm/s]           conv_vel[cm/s]           conv_vel/conv_vel_eq                L_conv/L'

      gradT = gradT_start
      Y_face = gradT_start - cburn_case%gradL
      conv_vel = conv_vel_start
      call print_track_row(0, 0d0, dt, conv_vel_start, cburn_case, gradT, Y_face, conv_vel, eq_conv_vel%val)

      do j = 1, num_steps
         call call_tdc(cburn_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call print_track_row(j, time, dt, conv_vel_start, cburn_case, gradT, Y_face, conv_vel, eq_conv_vel%val)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(cburn_case%L%val, cburn_case%gradr%val, gradT%val)
      end do

      write(*, '(a,1x,es24.16)') '      tau_phys[s]', tau_phys
      write(*, '(a,1x,es24.16)') '      conv_vel_eq[cm/s]', eq_conv_vel%val
      if (do_blank_line) write(*, '(a)') ''
   end subroutine print_cburn_growth_track


   subroutine print_track_row(step, time, dt, conv_vel_start, case_data, gradT, Y_face, conv_vel, conv_vel_eq)
      integer, intent(in) :: step
      real(dp), intent(in) :: time, dt, conv_vel_start
      type(tdc_case_data), intent(in) :: case_data
      type(auto_diff_real_star_order1), intent(in) :: gradT, Y_face, conv_vel
      real(dp), intent(in), optional :: conv_vel_eq
      real(dp) :: L_conv_div_L, conv_vel_div_eq
      real(dp), parameter :: tiny = 1d-99

      L_conv_div_L = Lconv_from_gradT(case_data%L%val, case_data%gradr%val, gradT%val) / case_data%L%val
      if (present(conv_vel_eq)) then
         conv_vel_div_eq = conv_vel%val / max(abs(conv_vel_eq), tiny)
         write(*, '(i10,8(1x,es24.16))') step, time, dt, gradT%val, Y_face%val, conv_vel_start, conv_vel%val, &
            conv_vel_div_eq, L_conv_div_L
      else
         write(*, '(i10,7(1x,es24.16))') step, time, dt, gradT%val, Y_face%val, conv_vel_start, conv_vel%val, L_conv_div_L
      end if
   end subroutine print_track_row


   subroutine run_decay_memory_sweep(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: decay_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, Lambda, L_conv_ref
      integer :: i, j, ierr
      real(dp), parameter :: gradT_start_old = 2.5204370043250246d-01
      real(dp), parameter :: memory_factors(5) = [0d0, 0.25d0, 0.5d0, 0.75d0, 1d0]
      character(len=24), parameter :: names(5) = [character(len=24) :: &
         'zero_memory', 'quarter_memory', 'half_memory', 'three_quarter_memory', 'full_memory']

      call setup_decay_case(decay_case)
      conv_vel_start = decay_case%conv_vel_start_ref
      L_conv_ref = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT_start_old)
      Lambda = decay_case%mixing_length_alpha * decay_case%scale_height%val

      do i = 1, size(memory_factors)
         L_conv_start = memory_factors(i) * L_conv_ref
         do j = 0, 40
            dt = pow(10d0, 2.3d0 + 0.02d0*j)
            tau_phys = arnett_relaxation_time(Lambda, decay_case%grav%val, decay_case%scale_height%val, &
               decay_case%gradL%val, decay_case%gradr%val, conv_vel_start, conv_vel_start)
            call call_tdc(decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
            call write_row(io, mode%name, 'decay_memory_sweep', names(i), -1, -1d0, dt, conv_vel_start, L_conv_start, &
               decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         end do
      end do
   end subroutine run_decay_memory_sweep


   subroutine run_decay_velocity_sweep(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: decay_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, Lambda, L_conv_ref
      integer :: i, j, ierr
      real(dp), parameter :: gradT_start_old = 2.5204370043250246d-01
      real(dp), parameter :: velocity_factors(4) = [0.25d0, 0.5d0, 0.75d0, 1d0]
      character(len=24), parameter :: names(4) = [character(len=24) :: &
         'quarter_a0', 'half_a0', 'three_quarter_a0', 'full_a0']

      call setup_decay_case(decay_case)
      L_conv_ref = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT_start_old)
      Lambda = decay_case%mixing_length_alpha * decay_case%scale_height%val
      L_conv_start = L_conv_ref

      do i = 1, size(velocity_factors)
         conv_vel_start = velocity_factors(i) * decay_case%conv_vel_start_ref
         do j = 0, 40
            dt = pow(10d0, 2.3d0 + 0.02d0*j)
            tau_phys = arnett_relaxation_time(Lambda, decay_case%grav%val, decay_case%scale_height%val, &
               decay_case%gradL%val, decay_case%gradr%val, conv_vel_start, decay_case%conv_vel_start_ref)
            call call_tdc(decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
            call write_row(io, mode%name, 'decay_velocity_sweep', names(i), -1, -1d0, dt, conv_vel_start, L_conv_start, &
               decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         end do
      end do
   end subroutine run_decay_velocity_sweep


   subroutine run_growth_sweep(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: growth_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Gamma
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, L_conv_eq
      integer :: i, j, ierr
      character(len=24), parameter :: names(2) = [character(len=24) :: 'radiative_start', 'equilibrium_start']

      call setup_growth_case(growth_case)
      call call_tdc(growth_case, mode, 0d0, 0d0, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      L_conv_eq = Lconv_from_gradT(growth_case%L%val, growth_case%gradr%val, eq_gradT%val)
      tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
         growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, 0d0, eq_conv_vel%val)

      call write_row(io, mode%name, 'reference', 'tdc_equilibrium', -1, -1d0, 1d40, 0d0, 0d0, &
         growth_case, eq_gradT, eq_Y_face, eq_conv_vel, tau_phys, ierr)

      do i = 1, size(names)
         do j = 0, 40
            dt = pow(10d0, -2d0 + 0.25d0*j)
            if (i == 1) then
               conv_vel_start = 0d0
               L_conv_start = 0d0
            else
               conv_vel_start = eq_conv_vel%val
               L_conv_start = L_conv_eq
            end if
            tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
               growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, conv_vel_start, eq_conv_vel%val)
            call call_tdc(growth_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
            call write_row(io, mode%name, 'growth_sweep', names(i), -1, -1d0, dt, conv_vel_start, L_conv_start, &
               growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         end do
      end do

      call call_cox_mlt(growth_case, gradT, Y_face, conv_vel, D, Gamma, ierr)
      tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
         growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, 0d0, conv_vel%val)
      call write_row(io, mode%name, 'reference', 'cox_mlt_equilibrium', -1, -1d0, 1d40, 0d0, 0d0, &
         growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
   end subroutine run_growth_sweep


   subroutine run_growth_flux_memory_sweep(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: growth_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, L_conv_eq
      integer :: i, j, ierr
      real(dp), parameter :: memory_factors(5) = [0d0, 0.25d0, 0.5d0, 0.75d0, 1d0]
      character(len=32), parameter :: names(5) = [character(len=32) :: &
         'zero_flux_memory', 'quarter_flux_memory', 'half_flux_memory', 'three_quarter_flux_memory', 'full_flux_memory']

      call setup_growth_case(growth_case)
      call call_tdc(growth_case, mode, 0d0, 0d0, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      L_conv_eq = Lconv_from_gradT(growth_case%L%val, growth_case%gradr%val, eq_gradT%val)
      conv_vel_start = 0d0

      do i = 1, size(memory_factors)
         L_conv_start = memory_factors(i) * L_conv_eq
         do j = 0, 40
            dt = pow(10d0, -2d0 + 0.25d0*j)
            tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
               growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, conv_vel_start, eq_conv_vel%val)
            call call_tdc(growth_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
            call write_row(io, mode%name, 'growth_flux_memory_sweep', names(i), -1, -1d0, dt, conv_vel_start, L_conv_start, &
               growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         end do
      end do
   end subroutine run_growth_flux_memory_sweep


   subroutine run_growth_velocity_sweep(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: growth_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys
      integer :: i, j, ierr
      real(dp), parameter :: velocity_factors(5) = [0d0, 0.25d0, 0.5d0, 0.75d0, 1d0]
      character(len=32), parameter :: names(5) = [character(len=32) :: &
         'zero_velocity_memory', 'quarter_velocity_memory', 'half_velocity_memory', 'three_quarter_velocity_memory', 'full_velocity_memory']

      call setup_growth_case(growth_case)
      call call_tdc(growth_case, mode, 0d0, 0d0, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      L_conv_start = 0d0

      do i = 1, size(velocity_factors)
         conv_vel_start = velocity_factors(i) * eq_conv_vel%val
         do j = 0, 40
            dt = pow(10d0, -2d0 + 0.25d0*j)
            tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
               growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, conv_vel_start, eq_conv_vel%val)
            call call_tdc(growth_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
            call write_row(io, mode%name, 'growth_velocity_sweep', names(i), -1, -1d0, dt, conv_vel_start, L_conv_start, &
               growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         end do
      end do
   end subroutine run_growth_velocity_sweep


   subroutine run_decay_evolution(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: decay_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time
      integer :: j, ierr
      real(dp), parameter :: gradT_start_old = 2.5204370043250246d-01
      integer, parameter :: num_steps = 12
      integer, parameter :: num_steps_long = 400

      call setup_decay_case(decay_case)
      conv_vel_start = decay_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT_start_old)
      tau_phys = arnett_relaxation_time(decay_case%mixing_length_alpha * decay_case%scale_height%val, decay_case%grav%val, &
         decay_case%scale_height%val, decay_case%gradL%val, decay_case%gradr%val, conv_vel_start, conv_vel_start)
      dt = pow(10d0, 2.3d0)

      gradT = gradT_start_old
      Y_face = gradT_start_old - decay_case%gradL
      conv_vel = conv_vel_start
      call write_row(io, mode%name, 'decay_evolution', 'track', 0, 0d0, dt, conv_vel_start, L_conv_start, &
         decay_case, gradT, Y_face, conv_vel, tau_phys, 0)

      do j = 1, num_steps
         call call_tdc(decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call write_row(io, mode%name, 'decay_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT%val)
      end do

      conv_vel_start = decay_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT_start_old)
      gradT = gradT_start_old
      Y_face = gradT_start_old - decay_case%gradL
      conv_vel = conv_vel_start
      call write_row(io, mode%name, 'decay_long_evolution', 'track', 0, 0d0, dt, conv_vel_start, L_conv_start, &
         decay_case, gradT, Y_face, conv_vel, tau_phys, 0)

      do j = 1, num_steps_long
         call call_tdc(decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call write_row(io, mode%name, 'decay_long_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT%val)
      end do
   end subroutine run_decay_evolution


   subroutine run_growth_evolution(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: growth_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time
      integer :: j, ierr
      integer, parameter :: num_steps = 120

      call setup_growth_case(growth_case)
      call call_tdc(growth_case, mode, 0d0, 0d0, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      if (ierr /= 0) return

      tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
         growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, 0d0, eq_conv_vel%val)
      call reference_evolution_dt(growth_case, 0d0, 0d0, dt, ierr)
      if (ierr /= 0) return
      conv_vel_start = 0d0
      L_conv_start = 0d0
      gradT = growth_case%gradr
      Y_face = growth_case%gradr - growth_case%gradL
      conv_vel = 0d0
      call write_row(io, mode%name, 'growth_evolution', 'track', 0, 0d0, dt, conv_vel_start, L_conv_start, &
         growth_case, gradT, Y_face, conv_vel, tau_phys, 0)

      do j = 1, num_steps
         call call_tdc(growth_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call write_row(io, mode%name, 'growth_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(growth_case%L%val, growth_case%gradr%val, gradT%val)
      end do
   end subroutine run_growth_evolution


   subroutine run_envelope_decay_evolution(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: envelope_decay_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time, gradT_start
      integer :: j, ierr
      integer, parameter :: num_steps = 120

      call setup_envelope_decay_case(envelope_decay_case, gradT_start)
      conv_vel_start = envelope_decay_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(envelope_decay_case%L%val, envelope_decay_case%gradr%val, gradT_start)
      tau_phys = arnett_relaxation_time(envelope_decay_case%mixing_length_alpha * envelope_decay_case%scale_height%val, &
         envelope_decay_case%grav%val, envelope_decay_case%scale_height%val, envelope_decay_case%gradL%val, &
         envelope_decay_case%gradr%val, conv_vel_start, 0d0)
      call reference_evolution_dt(envelope_decay_case, conv_vel_start, L_conv_start, dt, ierr)
      if (ierr /= 0) return
      gradT = gradT_start
      Y_face = gradT_start - envelope_decay_case%gradL
      conv_vel = conv_vel_start
      call write_row(io, mode%name, 'envelope_decay_evolution', 'track', 0, 0d0, dt, conv_vel_start, L_conv_start, &
         envelope_decay_case, gradT, Y_face, conv_vel, tau_phys, 0)

      do j = 1, num_steps
         call call_tdc(envelope_decay_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call write_row(io, mode%name, 'envelope_decay_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            envelope_decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(envelope_decay_case%L%val, envelope_decay_case%gradr%val, gradT%val)
      end do
   end subroutine run_envelope_decay_evolution


   subroutine run_cburn_growth_evolution(io, mode)
      integer, intent(in) :: io
      type(tdc_mode_data), intent(in) :: mode
      type(tdc_case_data) :: cburn_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D
      type(auto_diff_real_star_order1) :: eq_gradT, eq_Y_face, eq_conv_vel, eq_D
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time
      integer :: j, ierr
      real(dp), parameter :: gradT_start = 2.8060983127196415d-01
      integer, parameter :: num_steps = 120

      call setup_cburn_case(cburn_case)
      conv_vel_start = cburn_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(cburn_case%L%val, cburn_case%gradr%val, gradT_start)
      call call_tdc(cburn_case, mode, conv_vel_start, L_conv_start, 1d40, eq_gradT, eq_Y_face, eq_conv_vel, eq_D, ierr)
      if (ierr /= 0) return

      tau_phys = arnett_relaxation_time(cburn_case%mixing_length_alpha * cburn_case%scale_height%val, cburn_case%grav%val, &
         cburn_case%scale_height%val, cburn_case%gradL%val, cburn_case%gradr%val, conv_vel_start, eq_conv_vel%val)
      call reference_evolution_dt(cburn_case, conv_vel_start, L_conv_start, dt, ierr)
      if (ierr /= 0) return
      gradT = gradT_start
      Y_face = gradT_start - cburn_case%gradL
      conv_vel = conv_vel_start
      call write_row(io, mode%name, 'cburn_growth_evolution', 'track', 0, 0d0, dt, conv_vel_start, L_conv_start, &
         cburn_case, gradT, Y_face, conv_vel, tau_phys, 0)

      do j = 1, num_steps
         call call_tdc(cburn_case, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
         time = j * dt
         call write_row(io, mode%name, 'cburn_growth_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            cburn_case, gradT, Y_face, conv_vel, tau_phys, ierr)
         conv_vel_start = conv_vel%val
         L_conv_start = Lconv_from_gradT(cburn_case%L%val, cburn_case%gradr%val, gradT%val)
      end do
   end subroutine run_cburn_growth_evolution


   subroutine run_mlt_evolution(io)
      integer, intent(in) :: io
      type(tdc_case_data) :: decay_case, growth_case, cburn_case
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Gamma
      real(dp) :: conv_vel_start, L_conv_start, dt, tau_phys, time
      integer :: ierr, j
      real(dp), parameter :: gradT_start_old = 2.5204370043250246d-01
      integer, parameter :: decay_steps = 12
      integer, parameter :: decay_long_steps = 400
      integer, parameter :: growth_steps = 120
      real(dp) :: envelope_decay_gradT_start

      call setup_decay_case(decay_case)
      conv_vel_start = decay_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(decay_case%L%val, decay_case%gradr%val, gradT_start_old)
      tau_phys = arnett_relaxation_time(decay_case%mixing_length_alpha * decay_case%scale_height%val, decay_case%grav%val, &
         decay_case%scale_height%val, decay_case%gradL%val, decay_case%gradr%val, conv_vel_start, 0d0)
      dt = pow(10d0, 2.3d0)

      call call_cox_mlt(decay_case, gradT, Y_face, conv_vel, D, Gamma, ierr)

      do j = 0, decay_steps
         time = j * dt
         call write_row(io, 'mlt', 'decay_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
      end do

      do j = 0, decay_long_steps
         time = j * dt
         call write_row(io, 'mlt', 'decay_long_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            decay_case, gradT, Y_face, conv_vel, tau_phys, ierr)
      end do

      call setup_growth_case(growth_case)
      call call_cox_mlt(growth_case, gradT, Y_face, conv_vel, D, Gamma, ierr)
      tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
         growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, 0d0, conv_vel%val)
      call reference_evolution_dt(growth_case, 0d0, 0d0, dt, ierr)
      if (ierr /= 0) return

      do j = 0, growth_steps
         time = j * dt
         call write_row(io, 'mlt', 'growth_evolution', 'track', j, time, dt, 0d0, 0d0, &
            growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
      end do

      call setup_envelope_decay_case(growth_case, envelope_decay_gradT_start)
      conv_vel_start = growth_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(growth_case%L%val, growth_case%gradr%val, envelope_decay_gradT_start)
      call call_cox_mlt(growth_case, gradT, Y_face, conv_vel, D, Gamma, ierr)
      tau_phys = arnett_relaxation_time(growth_case%mixing_length_alpha * growth_case%scale_height%val, growth_case%grav%val, &
         growth_case%scale_height%val, growth_case%gradL%val, growth_case%gradr%val, conv_vel_start, 0d0)
      call reference_evolution_dt(growth_case, conv_vel_start, L_conv_start, dt, ierr)
      if (ierr /= 0) return

      do j = 0, growth_steps
         time = j * dt
         call write_row(io, 'mlt', 'envelope_decay_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            growth_case, gradT, Y_face, conv_vel, tau_phys, ierr)
      end do

      call setup_cburn_case(cburn_case)
      conv_vel_start = cburn_case%conv_vel_start_ref
      L_conv_start = Lconv_from_gradT(cburn_case%L%val, cburn_case%gradr%val, 2.8060983127196415d-01)
      call call_cox_mlt(cburn_case, gradT, Y_face, conv_vel, D, Gamma, ierr)
      tau_phys = arnett_relaxation_time(cburn_case%mixing_length_alpha * cburn_case%scale_height%val, cburn_case%grav%val, &
         cburn_case%scale_height%val, cburn_case%gradL%val, cburn_case%gradr%val, conv_vel_start, conv_vel%val)
      call reference_evolution_dt(cburn_case, conv_vel_start, L_conv_start, dt, ierr)
      if (ierr /= 0) return

      do j = 0, growth_steps
         time = j * dt
         call write_row(io, 'mlt', 'cburn_growth_evolution', 'track', j, time, dt, conv_vel_start, L_conv_start, &
            cburn_case, gradT, Y_face, conv_vel, tau_phys, ierr)
      end do
   end subroutine run_mlt_evolution


   subroutine reference_evolution_dt(case_data, conv_vel_start, L_conv_start, dt, ierr)
      type(tdc_case_data), intent(in) :: case_data
      real(dp), intent(in) :: conv_vel_start, L_conv_start
      real(dp), intent(out) :: dt
      integer, intent(out) :: ierr
      type(auto_diff_real_star_order1) :: ref_gradT, ref_Y_face, ref_conv_vel, ref_D
      real(dp) :: tau_ref

      call call_tdc(case_data, tdc_modes(1), conv_vel_start, L_conv_start, 1d40, ref_gradT, ref_Y_face, ref_conv_vel, ref_D, ierr)
      if (ierr /= 0) return

      tau_ref = arnett_relaxation_time(case_data%mixing_length_alpha * case_data%scale_height%val, case_data%grav%val, &
         case_data%scale_height%val, case_data%gradL%val, case_data%gradr%val, conv_vel_start, ref_conv_vel%val)
      dt = max(1d-6, 5d-2 * tau_ref)
   end subroutine reference_evolution_dt


   subroutine call_tdc(case_data, mode, conv_vel_start, L_conv_start, dt, gradT, Y_face, conv_vel, D, ierr)
      type(tdc_case_data), intent(in) :: case_data
      type(tdc_mode_data), intent(in) :: mode
      real(dp), intent(in) :: conv_vel_start, L_conv_start, dt
      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D
      integer, intent(out) :: ierr
      integer :: mixing_type, tdc_num_iters
      real(dp) :: gradT_start, Y_face_start
      real(dp), parameter :: tiny = 1d-30

      if (abs(case_data%L%val) <= tiny .or. abs(case_data%gradr%val) <= tiny) then
         gradT_start = case_data%gradL%val
      else
         gradT_start = case_data%gradr%val * (1d0 - L_conv_start/case_data%L%val)
      end if
      Y_face_start = gradT_start - case_data%gradL%val

      call set_TDC( &
         conv_vel_start, Y_face_start, case_data%mixing_length_alpha, case_data%TDC_alpha_D, case_data%TDC_alpha_R, &
         case_data%TDC_alpha_Pt, dt, case_data%cgrav, case_data%m, case_data%report, &
         mixing_type, case_data%scale, case_data%chiT, case_data%chiRho, case_data%gradr, case_data%r, case_data%P, &
         case_data%T, case_data%rho, case_data%dV, case_data%Cp, case_data%opacity, case_data%scale_height, &
         case_data%gradL, case_data%grada, conv_vel, D, Y_face, gradT, tdc_num_iters, case_data%max_conv_vel, &
         case_data%Eq_div_w, case_data%grav, mode%include_mlt_correction, case_data%TDC_alpha_C, case_data%TDC_alpha_S, &
         case_data%use_TDC_enthalpy_flux_limiter, mode%use_arnett, mode%use_acceleration_limit, mode%use_Af_split, mode%growth_target, &
         case_data%energy, ierr)
   end subroutine call_tdc


   subroutine call_cox_mlt(case_data, gradT, Y_face, conv_vel, D, Gamma, ierr)
      type(tdc_case_data), intent(in) :: case_data
      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D, Gamma
      integer, intent(out) :: ierr
      integer :: mixing_type
      real(dp) :: Henyey_MLT_nu_param, Henyey_MLT_y_param

      Henyey_MLT_nu_param = 0d0
      Henyey_MLT_y_param = 0d0

      call set_MLT('Cox', case_data%mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
         case_data%chiT, case_data%chiRho, case_data%Cp, case_data%grav, &
         case_data%mixing_length_alpha*case_data%scale_height, case_data%rho, case_data%P, case_data%T, &
         case_data%opacity, case_data%gradr, case_data%grada, case_data%gradL, &
         Gamma, gradT, Y_face, conv_vel, D, mixing_type, case_data%max_conv_vel, ierr)
   end subroutine call_cox_mlt


   subroutine setup_decay_case(case_data)
      type(tdc_case_data), intent(out) :: case_data

      case_data%max_conv_vel = 1d99
      case_data%conv_vel_start_ref = 52320587.415154047d0
      case_data%mixing_length_alpha = 2.0d0
      case_data%TDC_alpha_D = 1.0d0
      case_data%TDC_alpha_R = 0.0d0
      case_data%TDC_alpha_Pt = 0.0d0
      case_data%TDC_alpha_C = 1.0d0
      case_data%TDC_alpha_S = 1.0d0
      case_data%cgrav = 6.6743000000000004d-8
      case_data%m = 5.8707400456875664d34
      case_data%scale = 5.0386519362246294d45
      case_data%L = 1.0941528815883500015d0*(-5.0386519362246294d45)
      case_data%r = 10314294541.567163d0
      case_data%P = 5.0581587249808894d20
      case_data%T = 613193666.51783681d0
      case_data%rho = 5204.5732574745753d0
      case_data%dV = 3.8256494463482604d-7
      case_data%Cp = 6628075118.4606590d0
      case_data%opacity = 9.0750171231469945d-2
      case_data%scale_height = 2638686602.0063782d0
      case_data%gradL = 0.25207587267343501d0
      case_data%grada = 0.25204697256872738d0
      case_data%report = .false.
      case_data%chiT = 1d0
      case_data%chiRho = 1d0
      case_data%grav = case_data%m*case_data%cgrav/pow2(case_data%r)
      case_data%gradr = 3d0*case_data%P*case_data%opacity*case_data%L / &
         (64d0*pi*boltz_sigma*pow4(case_data%T)*case_data%cgrav*case_data%m)
      case_data%Eq_div_w = 0d0
      case_data%energy = 0d0
      case_data%use_TDC_enthalpy_flux_limiter = .false.
   end subroutine setup_decay_case


   subroutine setup_growth_case(case_data)
      type(tdc_case_data), intent(out) :: case_data

      case_data%max_conv_vel = 1d99
      case_data%conv_vel_start_ref = 0d0
      case_data%mixing_length_alpha = 2.0d0
      case_data%chiT = 1d0
      case_data%chiRho = 1d0
      case_data%T = 1d5
      case_data%rho = 1d-5
      case_data%r = Rsun
      case_data%m = Msun
      case_data%cgrav = standard_cgrav
      case_data%grav = case_data%m*case_data%cgrav/pow2(case_data%r)
      case_data%Cp = 2.5d0*kerg/mp
      case_data%P = case_data%rho*case_data%T*kerg/mp
      case_data%scale_height = case_data%P/(case_data%rho*case_data%grav)
      case_data%opacity = 1d0
      case_data%grada = 0.4d0
      case_data%gradL = case_data%grada
      case_data%L = 70d0*Lsun
      case_data%gradr = 3d0*case_data%P*case_data%opacity*case_data%L / &
         (64d0*pi*boltz_sigma*pow4(case_data%T)*case_data%grav*pow2(case_data%r))
      case_data%L = case_data%L*(1d0 + 1d-5)*(case_data%grada/case_data%gradr)
      case_data%gradr = 3d0*case_data%P*case_data%opacity*case_data%L / &
         (64d0*pi*boltz_sigma*pow4(case_data%T)*case_data%grav*pow2(case_data%r))
      case_data%TDC_alpha_D = 1.0d0
      case_data%TDC_alpha_R = 0.0d0
      case_data%TDC_alpha_Pt = 0.0d0
      case_data%TDC_alpha_C = 1.0d0
      case_data%TDC_alpha_S = 1.0d0
      case_data%dV = 0d0
      case_data%energy = 0d0
      case_data%scale = case_data%L%val*1d-3
      case_data%report = .false.
      case_data%Eq_div_w = 0d0
      case_data%use_TDC_enthalpy_flux_limiter = .false.
   end subroutine setup_growth_case


   subroutine setup_envelope_decay_case(case_data, gradT_start_ref)
      type(tdc_case_data), intent(out) :: case_data
      real(dp), intent(out) :: gradT_start_ref
      type(tdc_case_data) :: unstable_case
      type(auto_diff_real_star_order1) :: gradT_start, Y_face_start, conv_vel_start, D_start, Gamma_start
      integer :: ierr

      call setup_growth_case(unstable_case)
      call call_cox_mlt(unstable_case, gradT_start, Y_face_start, conv_vel_start, D_start, Gamma_start, ierr)

      case_data = unstable_case
      case_data%conv_vel_start_ref = conv_vel_start%val
      gradT_start_ref = gradT_start%val
      case_data%L = unstable_case%L * ((1d0 - 1d-5)/(1d0 + 1d-5))
      case_data%gradr = 3d0*case_data%P*case_data%opacity*case_data%L / &
         (64d0*pi*boltz_sigma*pow4(case_data%T)*case_data%grav*pow2(case_data%r))
      case_data%scale = case_data%L%val*1d-3
   end subroutine setup_envelope_decay_case


   subroutine setup_cburn_case(case_data)
      type(tdc_case_data), intent(out) :: case_data

      case_data%max_conv_vel = 1d99
      case_data%conv_vel_start_ref = 5.0383392852019933d4
      case_data%mixing_length_alpha = 2.0d0
      case_data%TDC_alpha_D = 1.0d0
      case_data%TDC_alpha_R = 0.0d0
      case_data%TDC_alpha_Pt = 0.0d0
      case_data%TDC_alpha_C = 1.0d0
      case_data%TDC_alpha_S = 1.0d0
      case_data%cgrav = standard_cgrav
      case_data%m = 1.1259123990149387d0*Msun
      case_data%scale = 1.0372316545177369d6*Lsun*1d-3
      case_data%L = 1.0372316545177369d6*Lsun
      case_data%r = 1.0818528002666345d9
      case_data%P = 8.7371150320499226d21
      case_data%T = 8.9448446840612292d8
      case_data%rho = 1.6585817351580659d5
      case_data%dV = 0d0
      case_data%Cp = 3.7545969237063849d8
      case_data%opacity = 6.8854059514891264d-2
      case_data%scale_height = 4.1262096121030027d8
      case_data%grada = 2.7977735453489061d-1
      ! The saved model provides gradT - grada, so use gradL = grada here.
      case_data%gradL = case_data%grada
      case_data%report = .false.
      case_data%chiT = 1.5046551755548503d0
      case_data%chiRho = 8.4356679250011213d-1
      case_data%grav = case_data%m*case_data%cgrav/pow2(case_data%r)
      case_data%gradr = 6.5970726825239474d0
      case_data%Eq_div_w = 0d0
      case_data%energy = 0d0
      case_data%use_TDC_enthalpy_flux_limiter = .false.
   end subroutine setup_cburn_case


   real(dp) function arnett_relaxation_time(Lambda, grav, Hp, gradL, gradr, conv_vel_start, conv_vel_eq) result(tau_phys)
      real(dp), intent(in) :: Lambda, grav, Hp, gradL, gradr, conv_vel_start, conv_vel_eq
      real(dp), parameter :: tiny = 1d-99
      real(dp) :: N, rate

      if (gradr > gradL) then
         tau_phys = max(Lambda, tiny) / max(abs(conv_vel_eq), tiny)
      else
         N = sqrt(max(0d0, grav/max(Hp, tiny) * (gradL - gradr)))
         rate = N + abs(conv_vel_start)/max(Lambda, tiny)
         tau_phys = 1d0 / max(rate, tiny)
      end if
   end function arnett_relaxation_time


   real(dp) function Lconv_from_gradT(L, gradr, gradT) result(L_conv)
      real(dp), intent(in) :: L, gradr, gradT

      if (gradr == 0d0) then
         L_conv = 0d0
      else
         L_conv = L * (1d0 - gradT/gradr)
      end if
   end function Lconv_from_gradT


   subroutine print_mode_controls(mode)
      type(tdc_mode_data), intent(in) :: mode

      if (mode%include_mlt_correction) then
         write(*, '(a)') '      include_mlt_corr_to_tdc = .true.'
      else
         write(*, '(a)') '      include_mlt_corr_to_tdc = .false.'
      end if

      if (mode%use_arnett) then
         write(*, '(a)') '      use_TDC_arnett_velocity_closure = .true.'
      else
         write(*, '(a)') '      use_TDC_arnett_velocity_closure = .false.'
      end if

      if (mode%use_acceleration_limit) then
         write(*, '(a)') '      use_TDC_acceleration_limit = .true.'
      else
         write(*, '(a)') '      use_TDC_acceleration_limit = .false.'
      end if
      if (mode%use_Af_split) then
         write(*, '(a)') '      use_TDC_Af_split = .true.'
      else
         write(*, '(a)') '      use_TDC_Af_split = .false.'
      end if
      write(*, '(a,a)') '      TDC_arnett_growth_target = ', trim(growth_target_name(mode%growth_target))
   end subroutine print_mode_controls


   subroutine print_case_setup(case_data, conv_vel_start, gradT_start, interpretation)
      type(tdc_case_data), intent(in) :: case_data
      real(dp), intent(in) :: conv_vel_start, gradT_start
      character(len=*), intent(in) :: interpretation
      real(dp) :: L_conv_start, L_conv_div_L

      L_conv_start = Lconv_from_gradT(case_data%L%val, case_data%gradr%val, gradT_start)
      if (case_data%L%val == 0d0) then
         L_conv_div_L = 0d0
      else
         L_conv_div_L = L_conv_start/case_data%L%val
      end if

      write(*, '(a,1x,es24.16,1x,a,1x,es24.16)') '      T[K] =', case_data%T%val, 'rho[g/cm^3] =', case_data%rho%val
      write(*, '(a,1x,es24.16,1x,a,1x,es24.16)') '      P[dyn/cm^2] =', case_data%P%val, 'L[erg/s] =', case_data%L%val
      write(*, '(a,1x,es24.16,1x,a,1x,es24.16)') '      r[cm] =', case_data%r%val, 'Hp[cm] =', case_data%scale_height%val
      write(*, '(a,1x,es24.16,1x,a,1x,es24.16,1x,a,1x,es24.16)') '      gradL =', case_data%gradL%val, &
         'gradr =', case_data%gradr%val, 'grada =', case_data%grada%val
      write(*, '(a,1x,es24.16,1x,a,1x,es24.16,1x,a,1x,es24.16)') '      conv_vel_start[cm/s] =', conv_vel_start, &
         'gradT_start =', gradT_start, 'L_conv_start/L =', L_conv_div_L
      write(*, '(a,a)') '      Interpretation: ', trim(interpretation)
      write(*, '(a)') ''
   end subroutine print_case_setup

   character(len=32) function growth_target_name(growth_target) result(name)
      integer, intent(in) :: growth_target

      select case (growth_target)
      case (TDC_arnett_growth_target_mlt)
         name = 'mlt'
      case (TDC_arnett_growth_target_tdc_no_mlt_corr)
         name = 'tdc_no_mlt_corr'
      case (TDC_arnett_growth_target_tdc_with_mlt_corr)
         name = 'tdc_with_mlt_corr'
      case default
         name = 'unknown'
      end select
   end function growth_target_name


   subroutine write_row(io, mode_name, group_name, scenario_name, step, time, dt, conv_vel_start, L_conv_start, &
         case_data, gradT, Y_face, conv_vel, tau_phys, ierr)
      integer, intent(in) :: io, step, ierr
      character(len=*), intent(in) :: mode_name, group_name, scenario_name
      real(dp), intent(in) :: time, dt, conv_vel_start, L_conv_start, tau_phys
      type(tdc_case_data), intent(in) :: case_data
      type(auto_diff_real_star_order1), intent(in) :: gradT, Y_face, conv_vel
      real(dp) :: L_conv_end, L_conv_div_L

      L_conv_end = Lconv_from_gradT(case_data%L%val, case_data%gradr%val, gradT%val)
      if (case_data%L%val == 0d0) then
         L_conv_div_L = 0d0
      else
         L_conv_div_L = L_conv_end / case_data%L%val
      end if

      write(io, '(a,",",a,",",a,",",*(g0,:,","))') trim(mode_name), trim(group_name), trim(scenario_name), &
         step, time, dt, conv_vel_start, L_conv_start, case_data%L%val, case_data%gradr%val, case_data%gradL%val, &
         case_data%grada%val, gradT%val, Y_face%val, conv_vel%val, L_conv_end, L_conv_div_L, tau_phys, ierr
   end subroutine write_row

end module test_time_dependence_support
