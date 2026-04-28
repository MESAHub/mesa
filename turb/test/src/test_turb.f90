program test_turb

   use math_lib
   use auto_diff
   use const_def, only: dp, pi, rsun, lsun, msun, kerg, mp, boltz_sigma, standard_cgrav
   use turb
   use test_time_dependence_support, only: check_time_dependence, write_time_dependence_csv

   implicit none

   integer, parameter :: num_tdc_modes = 4
   character(len=32), parameter :: tdc_mode_names(num_tdc_modes) = [character(len=32) :: &
      'plain TDC', 'TDC + Af split', 'TDC + Arnett closure', 'TDC + acceleration limit']
   logical, parameter :: tdc_mode_use_arnett(num_tdc_modes) = [.false., .false., .true., .false.]
   logical, parameter :: tdc_mode_use_accel(num_tdc_modes) = [.false., .false., .false., .true.]
   logical, parameter :: tdc_mode_use_split(num_tdc_modes) = [.false., .true., .false., .false.]

   call check_efficient_MLT_scaling()
   call check_TDC()
   call compare_TDC_and_Cox_MLT()
   call header('Test Time Dependence')
   call check_time_dependence()
   call write_test_time_dependence_csv()

contains

   subroutine header(text)
      character(len=*), intent(in) :: text

      write (*, '(a)') ' ----------------------------------------------------------------'
      write (*, '(a)') ' '
      write (*, '(a)') ' '//text
      write (*, '(a)') ' '
      write (*, '(a)') ' ----------------------------------------------------------------'

   end subroutine header

   subroutine check_efficient_MLT_scaling()
      type(auto_diff_real_star_order1) :: chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=3) :: MLT_option
      real(dp) :: mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, max_conv_vel
      type(auto_diff_real_star_order1) :: Gamma, gradT, Y_face, conv_vel, conv_vel2, D, r, L
      integer :: mixing_type, ierr

      include 'formats'

      call header('Test MLT')

      MLT_option = 'Cox'
      mixing_length_alpha = 1d0
      Henyey_MLT_nu_param = 0d0
      Henyey_MLT_y_param = 0d0
      chiT = 1d0
      chiRho = 1d0

      T = 1d5
      rho = 1d-5
      grav = 1d4
      r = Rsun
      Cp = 2.5d0*kerg/mp
      P = rho*T*kerg/mp
      Lambda = P/(rho*grav)
      opacity = 1d0
      grada = 0.4d0
      gradL = grada

      L = 1d5*Lsun
      gradr = 3d0*P*opacity*L/(64*pi*boltz_sigma*pow4(T)*grav*pow2(r))

      max_conv_vel = 1d99

      call set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                   chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                   gradr, grada, gradL, &
                   Gamma, gradT, Y_face, conv_vel, D, mixing_type, max_conv_vel, ierr)

      write (*, 1) 'vc at 1d5 Lsun', conv_vel%val

      L = 1d8*Lsun
      gradr = 3d0*P*opacity*L/(64*pi*boltz_sigma*pow4(T)*grav*pow2(r))

      call set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                   chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                   gradr, grada, gradL, &
                   Gamma, gradT, Y_face, conv_vel2, D, mixing_type, max_conv_vel, ierr)

      write (*, 1) 'vc at 1d8 Lsun', conv_vel2%val

      write (*, 1) 'Ratio:', conv_vel2%val/conv_vel%val
      write (*, '(a)') 'Expected ~10 because in the efficient limit vc ~ L^{1/3}'
   end subroutine check_efficient_MLT_scaling

   subroutine compare_TDC_and_Cox_MLT()
      real(dp) :: mixing_length_alpha, conv_vel_start, &
         TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, scale, TDC_alpha_C, TDC_alpha_S
      type(auto_diff_real_star_order1) :: &
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, gradL, grav, Lambda
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Gamma, Eq_div_w, energy
      real(dp) :: Henyey_MLT_nu_param, Henyey_MLT_y_param, max_conv_vel

      character(len=3) :: MLT_option
      integer :: mixing_type, ierr, tdc_num_iters, mode_i
      logical :: report, include_mlt_corr_to_TDC, use_TDC_enthalpy_flux_limiter

      include 'formats'

      call header('Compare MLT and TDC Modes')

      ! For limiting the conv_vel coming out of mlt/TDC with Csound.
      max_conv_vel = 1d99 ! we don't limit the conv_vel for testing.

      ! General
      mixing_length_alpha = 2.0d0
      chiT = 1d0
      chiRho = 1d0
      T = 1d5
      rho = 1d-5
      r = Rsun
      m = Msun
      cgrav = standard_cgrav
      grav = m*cgrav/pow2(r)
      Cp = 2.5d0*kerg/mp
      P = rho*T*kerg/mp
      scale_height = P/(rho*grav)
      Lambda = mixing_length_alpha*scale_height
      opacity = 1d0
      grada = 0.4d0
      gradL = grada

      L = 70*Lsun
      gradr = 3d0*P*opacity*L/(64*pi*boltz_sigma*pow4(T)*grav*pow2(r))

      ! Adjust L down to get just slightly superadiabatic gradR)
      L = L*(1d0 + 1d-5)*(grada/gradr)
      gradr = 3d0*P*opacity*L/(64*pi*boltz_sigma*pow4(T)*grav*pow2(r))

      ! TDC
      TDC_alpha_D = 1.0d0
      TDC_alpha_R = 0.0d0
      TDC_alpha_Pt = 0.0d0
      TDC_alpha_C = 1.0d0
      TDC_alpha_S = 1.0d0
      dV = 0d0
      energy = 0d0
      conv_vel_start = 0d0  !1d10
      scale = L%val*1d-3
      report = .false.
      dt = 1d40 ! Long time-step so we get into equilibrium
      Eq_div_w = 0d0 ! TDC_alpha_M is implicit in this term
      include_mlt_corr_to_TDC = .true.
      use_TDC_enthalpy_flux_limiter = .false.
      ! MLT
      MLT_option = 'Cox'
      Henyey_MLT_nu_param = 0d0
      Henyey_MLT_y_param = 0d0

      write (*, 1) 'gradR - gradA', gradr%val - grada%val

      call set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                   chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                   gradr, grada, gradL, &
                   Gamma, gradT, Y_face, conv_vel, D, mixing_type, max_conv_vel, ierr)

      write (*, '(a)') 'Mode: MLT'
      write (*, 1) 'Y, conv_vel_start, conv_vel, Gamma', Y_face%val, conv_vel_start, conv_vel%val, Gamma%val

      do mode_i = 1, num_tdc_modes
         conv_vel_start = 0d0
         call set_TDC( &
            conv_vel_start, 0d0, mixing_length_alpha, TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, report, &
            mixing_type, scale, chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, max_conv_vel, &
            Eq_div_w, grav, include_mlt_corr_to_TDC, TDC_alpha_C, TDC_alpha_S, use_TDC_enthalpy_flux_limiter, &
            tdc_mode_use_arnett(mode_i), tdc_mode_use_accel(mode_i), tdc_mode_use_split(mode_i), &
            .false., 0d0, 0d0, TDC_arnett_growth_target_mlt, energy, ierr)
         write (*, '(a)') 'Mode: ' // trim(tdc_mode_names(mode_i))
         write (*, 1) 'Y, conv_vel_start, conv_vel, dt', Y_face%val, conv_vel_start, conv_vel%val, dt
      end do

   end subroutine compare_TDC_and_Cox_MLT

   subroutine check_TDC()
      real(dp) :: mixing_length_alpha, conv_vel_start
      real(dp) :: TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, scale, max_conv_vel, TDC_alpha_C, TDC_alpha_S
      type(auto_diff_real_star_order1) :: &
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, gradL
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Eq_div_w, grav, energy
      integer :: mixing_type, ierr, tdc_num_iters, mode_i
      logical :: report, include_mlt_corr_to_TDC, use_TDC_enthalpy_flux_limiter
      integer :: j
      real(dp), parameter :: gradT_start_old = 2.5204370043250246d-01

      include 'formats'

      call header('Test TDC Modes')

      ! For limiting the conv_vel coming out of mlt/TDC with Csound.
      max_conv_vel = 1d99 ! we don't limit the conv_vel for testing.

      conv_vel_start = 52320587.415154047d0

      mixing_length_alpha = 2.0d0
      TDC_alpha_D = 1.0d0
      TDC_alpha_R = 0.0d0
      TDC_alpha_Pt = 0.0d0
      TDC_alpha_C = 1.0d0
      TDC_alpha_S = 1.0d0
      cgrav = 6.6743000000000004d-8
      m = 5.8707400456875664d34
      scale = 5.0386519362246294d45
      L = 1.0941528815883500015d0*(-5.0386519362246294d45)
      r = 10314294541.567163d0
      P = 5.0581587249808894d20
      T = 613193666.51783681d0
      rho = 5204.5732574745753d0
      dV = 3.8256494463482604d-7
      Cp = 6628075118.4606590d0
      opacity = 9.0750171231469945d-2
      scale_height = 2638686602.0063782d0
      gradL = 0.25207587267343501d0
      grada = 0.25204697256872738d0
      report = .false.
      chiT = 1d0
      chiRho = 1d0

      gradr = 3d0 * P * opacity * L / (64 * pi * boltz_sigma * pow4(T) * cgrav * m)
      grav = m * cgrav / pow2(r)
      Eq_div_w = 0d0 ! TDC_alpha_M is implicit in this term
      energy = 0d0
      include_mlt_corr_to_TDC = .true.
      use_TDC_enthalpy_flux_limiter = .false.
      do mode_i = 1, num_tdc_modes
         write (*, '(a)') '####################################'
         write (*, '(a)') 'Mode: ' // trim(tdc_mode_names(mode_i))
         write (*, '(a)') 'Running dt test'

         do j = 0, 30
            dt = 500d0*pow(1.02d0, j)
            call set_TDC( &
               conv_vel_start, gradT_start_old - gradL%val, mixing_length_alpha, TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, report, &
               mixing_type, scale, chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
               scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, max_conv_vel, &
               Eq_div_w, grav, include_mlt_corr_to_TDC, TDC_alpha_C, TDC_alpha_S, use_TDC_enthalpy_flux_limiter, &
               tdc_mode_use_arnett(mode_i), tdc_mode_use_accel(mode_i), tdc_mode_use_split(mode_i), &
               .false., 0d0, 0d0, TDC_arnett_growth_target_mlt, energy, ierr)

            write (*, 1) 'dt, gradT, conv_vel_start, conv_vel', dt, gradT%val, conv_vel_start, conv_vel%val
            if (report) stop
         end do
      end do

   end subroutine check_TDC

   subroutine write_test_time_dependence_csv()
      integer :: ierr

      call write_time_dependence_csv('plotter/time_dependence.csv', ierr)
      if (ierr /= 0) stop 1
   end subroutine write_test_time_dependence_csv

end program test_turb
