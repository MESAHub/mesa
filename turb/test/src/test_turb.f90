program test_turb

   use math_lib
   use auto_diff
   use const_def, only: dp, pi, rsun, lsun, msun, kerg, mp, boltz_sigma, standard_cgrav, sqrt_2_div_3
   use turb
   use const_def, only: crad, clight
   use utils_lib, only: mesa_error

   implicit none

   call check_efficient_MLT_scaling()
   call check_TDC()
   call check_TDC_zero_start()
   call check_TDC_degenerate()
   call check_TDC_radiative_root()
   call check_TDC_structure_guard()
   call check_TDC_far_active_branch()
   call compare_TDC_and_Cox_MLT()

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
      type(auto_diff_real_star_order1) :: &
         mixing_length_alpha, chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=3) :: MLT_option
      real(dp) :: Henyey_MLT_nu_param, Henyey_MLT_y_param, max_conv_vel
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
      real(dp) :: conv_vel_start, &
         TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, scale, L_start, TDC_alpha_C, TDC_alpha_S
      type(auto_diff_real_star_order1) :: &
         mixing_length_alpha, r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, &
         scale_height, gradL, grav, Lambda
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Gamma, Eq_div_w, energy
      real(dp) :: Henyey_MLT_nu_param, Henyey_MLT_y_param, max_conv_vel

      character(len=3) :: MLT_option
      integer :: mixing_type, ierr, tdc_num_iters
      logical :: report, include_mlt_corr_to_TDC, use_TDC_enthalpy_flux_limiter

      include 'formats'

      call header('Compare TDC with MLT Cox')

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

      call set_TDC( &
         conv_vel_start, mixing_length_alpha, Lambda, TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, report, &
         mixing_type, scale, chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
         scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, max_conv_vel, &
         Eq_div_w, grav, include_mlt_corr_to_TDC, TDC_alpha_C, TDC_alpha_S, use_TDC_enthalpy_flux_limiter, &
         energy, 0d0, ierr)


      write (*, 1) 'TDC: Y, conv_vel_start, conv_vel, dt   ', Y_face%val, conv_vel_start, conv_vel%val, dt

      call set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                   chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                   gradr, grada, gradL, &
                   Gamma, gradT, Y_face, conv_vel, D, mixing_type, max_conv_vel, ierr)

      write (*, 1) 'MLT: Y, conv_vel_start, conv_vel, Gamma', Y_face%val, conv_vel_start, conv_vel%val, Gamma%val

   end subroutine compare_TDC_and_Cox_MLT

   subroutine check_TDC_radiative_root()
      ! Stable face at the helium-flash remesh failure (model 10011).
      ! Long steps admit both active and radiative roots. Short steps must
      ! retain the starting turbulent velocity rather than force radiation.
      type(auto_diff_real_star_order1) :: &
         mixing_length_alpha, Lambda, r, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, &
         scale_height, gradL, grav, gradT, Y_face, conv_vel, D, Eq_div_w, energy
      real(dp) :: cgrav, m, conv_vel_start, L0, c0, expected_Y, dt
      integer :: j, mixing_type, ierr, tdc_num_iters

      cgrav = standard_cgrav
      m = 8.8836249890374269d32
      r = 1.7294445833353324d9
      T = 3.2863823906154651d7
      P = 1.3538674774076137d19
      rho = 5.6379903213572588d3
      opacity = 8.3687119532310605d-1
      Cp = 1.2298073267445153d8
      scale_height = 1.2113487345139594d8
      mixing_length_alpha = 2d0
      Lambda = mixing_length_alpha*scale_height
      grav = cgrav*m/pow2(r)
      chiRho = 1d0
      chiT = 1d0
      grada = 3.9264489938751634d-1
      gradL = grada
      gradr = 3.3550225088181102d-1
      gradr%d1Array(1) = 1d0
      dV = 0d0
      Eq_div_w = 0d0
      energy = 0d0
      conv_vel_start = 8.5914510924155285d4
      L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T%val)/(P%val*opacity%val)
      c0 = 0.5d0*sqrt_2_div_3*mixing_length_alpha%val*rho%val*T%val*Cp%val*4d0*pi*pow2(r%val)
      expected_Y = L0*(gradr%val - gradL%val)/(L0 + c0*conv_vel_start/sqrt_2_div_3)

      do j = 1, 2
         dt = 1d-20
         if (j == 2) dt = 2477.8759383372562d0
         call set_TDC( &
            conv_vel_start, mixing_length_alpha, Lambda, 1d0, 0d0, 0d0, dt, cgrav, m, .false., &
            mixing_type, L0, chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, 1d99, &
            Eq_div_w, grav, .false., 1d0, 1d0, .false., energy, 0d0, ierr, gradr%val - gradL%val)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'TDC stable-face solve failed')
         if (j == 1) then
            if (abs(conv_vel%val/conv_vel_start - 1d0) > 1d-10 .or. &
                  abs(Y_face%val/expected_Y - 1d0) > 1d-7) &
               call mesa_error(__FILE__, __LINE__, 'TDC short step lost convective memory')
         else
            if (conv_vel%val /= 0d0 .or. D%val /= 0d0 .or. &
                  abs(gradT%val - gradr%val) > 1d-12 .or. &
                  abs(gradT%d1Array(1) - 1d0) > 1d-12) &
               call mesa_error(__FILE__, __LINE__, 'TDC did not select the valid radiative root')
         end if
      end do
   end subroutine check_TDC_radiative_root

   subroutine check_TDC_structure_guard()
      ! Recorded millisecond trial: active and radiative roots coexist, but
      ! only the active root matches the stellar temperature gradient.
      type(auto_diff_real_star_order1) :: mixing_length_alpha, Lambda, r, T, P, opacity, rho, dV, &
         chiRho, chiT, Cp, gradr, grada, scale_height, gradL, grav, gradT, Y_face, conv_vel, D, Eq_div_w, energy
      real(dp) :: cgrav, m, v_start, L0, c0, L, residual, Y_structure
      integer :: j, mixing_type, ierr, tdc_num_iters

      cgrav = standard_cgrav
      m = 8.883624989037365d32
      r = 1.726807947376931d9
      T = 3.284365522706588d7
      P = 1.3617767634326184d19
      rho = 5.658725701446729d3
      opacity = 8.320572356602799d-1
      Cp = 1.2315419542635994d8
      scale_height = 1.2102620447360665d8
      mixing_length_alpha = 2d0
      Lambda = mixing_length_alpha*scale_height
      grav = cgrav*m/pow2(r)
      chiRho = 1d0
      chiT = 1d0
      grada = 3.926951774730907d-1
      gradL = grada
      gradr = -4.1122551917529997d5
      dV = 0d0
      Eq_div_w = 0d0
      energy = 0d0
      v_start = 8.591451092415529d4
      Y_structure = 3.9256603764530656d-1 - gradL%val
      L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T%val)/(P%val*opacity%val)
      c0 = 0.5d0*sqrt_2_div_3*mixing_length_alpha%val*rho%val*T%val*Cp%val*4d0*pi*pow2(r%val)
      L = L0*gradr%val

      do j = 1, 2
         call set_TDC( &
            v_start, mixing_length_alpha, Lambda, 1d0, 0d0, 0d0, 1d-3, cgrav, m, .false., &
            mixing_type, abs(L), chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, 1d99, &
            Eq_div_w, grav, .false., 1d0, 1d0, .false., energy, 0d0, ierr, Y_structure)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'TDC structure guard solve failed')
         if (j == 1) then
            if (abs(conv_vel%val - 8.591390530398137d4) > 1d-4 .or. &
                  abs(Y_face%val + 1.291406056182205d-4) > 1d-10) &
               call mesa_error(__FILE__, __LINE__, 'TDC structure guard discarded a valid active root')
            residual = L - L0*(gradL%val + Y_face%val) - c0*conv_vel%val/sqrt_2_div_3*Y_face%val
            if (abs(residual)/abs(L) > 1d-8) &
               call mesa_error(__FILE__, __LINE__, 'TDC guarded active-root luminosity closure failed')
         else
            if (conv_vel%val /= 0d0 .or. abs(gradT%val/gradr%val - 1d0) > 1d-12) &
               call mesa_error(__FILE__, __LINE__, 'TDC structure guard failed to permit quenching')
         end if
         Y_structure = gradr%val - gradL%val
      end do
   end subroutine check_TDC_structure_guard

   subroutine check_TDC_far_active_branch()
      ! Late decay at the recorded equL failure in cell 1183. Construct the
      ! luminosity from the temperature gradient, which lies beyond the fold.
      type(auto_diff_real_star_order1) :: mixing_length_alpha, Lambda, r, T, P, opacity, rho, dV, &
         chiRho, chiT, Cp, gradr, grada, scale_height, gradL, grav, gradT, Y_face, conv_vel, D, Eq_div_w, energy
      real(dp) :: cgrav, m, v_start, dt, Y_structure, L0, c0, S0, D0, A0, J0, expected_Af, L, selection_Y
      integer :: j, mixing_type, ierr, tdc_num_iters

      cgrav = standard_cgrav
      m = 9.218977444390872d32
      r = 2.6113190958022785d9
      T = 3.545580611754425d7
      P = 4.7170213391744184d16
      rho = 1.1765569185679698d1
      opacity = 2.769445570492667d-1
      Cp = 4.6573048478812575d8
      scale_height = 4.443102293673566d8
      mixing_length_alpha = 2d0
      Lambda = mixing_length_alpha*scale_height
      grav = cgrav*m/pow2(r)
      chiRho = 1d0
      chiT = 1d0
      grada = 3.326410377083683d-1
      gradL = grada
      dV = 0d0
      Eq_div_w = 0d0
      energy = 0d0
      v_start = 4.5848688619071245d0
      dt = 1.0935d-6
      Y_structure = -3.947322779852037d-1 - gradL%val
      A0 = v_start/sqrt_2_div_3
      S0 = 0.5d0*sqrt_2_div_3*mixing_length_alpha%val*Cp%val*T%val/scale_height%val*grada%val
      D0 = (8d0/3d0)*sqrt_2_div_3/Lambda%val
      J0 = sqrt(-4d0*S0*Y_structure*D0)
      expected_Af = -J0*tan(J0*dt/4d0 - atan(2d0*A0*D0/J0))/(2d0*D0)
      L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T%val)/(P%val*opacity%val)
      c0 = 0.5d0*sqrt_2_div_3*mixing_length_alpha%val*rho%val*T%val*Cp%val*4d0*pi*pow2(r%val)
      L = L0*(gradL%val + Y_structure) + c0*expected_Af*Y_structure
      gradr = L/L0
      gradr%d1Array(1) = 1d0

      do j = 1, 3
         selection_Y = Y_structure
         if (j == 2) selection_Y = 0d0 ! No structural preference: retain the original closest root.
         ! An unavailable constraint must also retain the original branch,
         ! even when the numerical value would select the far active root.
         call set_TDC( &
            v_start, mixing_length_alpha, Lambda, 1d0, 0d0, 0d0, dt, cgrav, m, .false., &
            mixing_type, abs(L), chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, 1d99, &
            Eq_div_w, grav, .false., 1d0, 1d0, .false., energy, 0d0, ierr, selection_Y, j /= 3)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'TDC far active-branch solve failed')
         if (j == 1) then
            if (abs(Y_face%val/Y_structure - 1d0) > 1d-6 .or. &
                  abs(conv_vel%val/(sqrt_2_div_3*expected_Af) - 1d0) > 1d-6 .or. Y_face%d1Array(1) >= 0d0) &
               call mesa_error(__FILE__, __LINE__, 'TDC discarded the active branch required by the temperature equation')
         else
            if (Y_face%val <= Y_structure .or. conv_vel%val <= 0d0 .or. Y_face%d1Array(1) <= 0d0) &
               call mesa_error(__FILE__, __LINE__, 'TDC original closest-root fallback changed')
         end if
      end do
   end subroutine check_TDC_far_active_branch

   subroutine check_TDC()
      real(dp) :: conv_vel_start
      real(dp) :: TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, scale, max_conv_vel, L_start, TDC_alpha_C, TDC_alpha_S
      type(auto_diff_real_star_order1) :: &
         mixing_length_alpha, Lambda, r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, gradL
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Eq_div_w, grav, energy
      integer :: mixing_type, ierr, tdc_num_iters
      logical :: report, include_mlt_corr_to_TDC, use_TDC_enthalpy_flux_limiter
      integer :: j

      include 'formats'

      call header('Test TDC')

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
      Lambda = mixing_length_alpha*scale_height
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

      write (*, *) "####################################"
      write (*, *) "Running dt test"

      do j = 0, 30
         dt = 500d0*pow(1.02d0, j)
         call set_TDC( &
            conv_vel_start, mixing_length_alpha, Lambda, TDC_alpha_D, TDC_alpha_R, TDC_alpha_Pt, dt, cgrav, m, report, &
            mixing_type, scale, chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, max_conv_vel, &
            Eq_div_w, grav, include_mlt_corr_to_TDC, TDC_alpha_C, TDC_alpha_S, use_TDC_enthalpy_flux_limiter, &
            energy, 0d0, ierr)


         write (*, 1) 'dt, gradT, conv_vel_start, conv_vel', dt, gradT%val, conv_vel_start, conv_vel%val
         if (report) stop
      end do

   end subroutine check_TDC

   subroutine check_TDC_zero_start()
      use tdc, only: get_TDC_solution
      use tdc_support, only: tdc_info, compute_Q, convert
      type(tdc_info) :: info
      type(auto_diff_real_tdc) :: Zlb, Zub, Y, Q, Af
      type(auto_diff_real_star_order1) :: conv_vel, Y_face, shear_Y, shear_v
      real(dp) :: expected_Y, expected_Af
      integer :: j, ierr, tdc_num_iters

      info%report = .false.
      info%include_mlt_corr_to_TDC = .false.
      info%use_TDC_enthalpy_flux_limiter = .false.
      info%mixing_length_alpha = 1d0
      info%Lambda = 1d0
      info%TDC_alpha_C = 1d0
      info%TDC_alpha_S = 2d0/sqrt_2_div_3
      info%TDC_alpha_D = 0d0
      info%TDC_alpha_R = 0d0
      info%TDC_alpha_Pt = 0d0
      info%dt = 1d0
      info%e = 0d0
      info%c0 = 1d0
      info%L = -1d0
      info%L%d1Array(1) = 1d0
      info%L0 = 1d0
      info%gradL = 0d0
      info%grada = 1d0
      info%T = 1d0
      info%rho = 1d0
      info%dV = 0d0
      info%Cp = 1d0
      info%kap = 1d0
      info%Hp = 1d0
      info%Gamma = 0d0
      info%P = 1d0
      info%h = 1d0
      Zlb = -50d0
      Zub = 50d0
      Zlb%d1val1 = 1d0
      Zub%d1val1 = 1d0

      do j = 1, 7
         info%A0 = 0d0
         info%Eq_div_w = 4d0
         info%L = -1d0
         info%L%d1Array(1) = 1d0
         if (j == 2) info%A0 = 5d-11
         if (j == 3) info%Eq_div_w = 0d0
         if (j == 4) info%Eq_div_w = 0.5d0
         if (j >= 5) info%L%val = 0d0
         if (j == 6) info%Eq_div_w = 0d0
         if (j == 7) info%A0 = 1d0
         call get_TDC_solution(info, 1d0, Zlb, Zub, conv_vel, Y_face, tdc_num_iters, 0d0, ierr)
         if (ierr /= 0) stop 1
         if (j == 1) then
            shear_Y = Y_face
            shear_v = conv_vel
            call get_TDC_solution(info, 1d0, Zlb, Zub, conv_vel, Y_face, tdc_num_iters, 0d0, ierr, -1d0)
            if (ierr /= 0) stop 1
            if (Y_face%val /= shear_Y%val .or. any(Y_face%d1Array /= shear_Y%d1Array) .or. &
                  conv_vel%val /= shear_v%val .or. any(conv_vel%d1Array /= shear_v%d1Array)) stop 1
         end if
         Y = convert(Y_face)
         call compute_Q(info, Y, Q, Af)
         expected_Af = max(0d0, info%A0%val + 0.5d0*(Y_face%val + info%Eq_div_w%val))
         if (abs(Af%val - expected_Af) > 1d-12 .or. &
               abs(conv_vel%val/sqrt_2_div_3 - expected_Af) > 1d-12) stop 1
         if (abs(Q%val) > 1d-7 .or. maxval(abs(Q%d1Array)) > 1d-6) then
            write(*,*) 'Zero-start TDC luminosity closure failed', j, Q%val
            stop 1
         end if
         if (j <= 2) then
            ! Af = A0 + (4 + Y)/2, so Q = -1 - (3 + A0)*Y - Y**2/2.
            expected_Y = -(3d0 + info%A0%val) + sqrt(pow2(3d0 + info%A0%val) - 2d0)
            if (abs(Y_face%val - expected_Y) > 1d-7 .or. conv_vel%val <= 0d0) stop 1
         else if (j <= 4) then
            if (Y_face%val /= -1d0 .or. conv_vel%val /= 0d0 .or. tdc_num_iters /= 0) stop 1
         else
            if (Y_face%val /= 0d0 .or. tdc_num_iters /= 0) stop 1
            if (j == 6) then
               if (conv_vel%val /= 0d0) stop 1
            else
               if (conv_vel%val <= 0d0) stop 1
            end if
         end if
      end do
   end subroutine check_TDC_zero_start


   subroutine check_TDC_degenerate()
      use tdc_support, only: eval_Af
      type(auto_diff_real_tdc) :: A0, xi0, xi1, xi2, Af

      xi0 = -1d0
      xi1 = 2d0
      xi2 = -1d0
      A0 = 0.5d0
      ! This solution reaches zero before the pole of its continuation.
      Af = eval_Af(10d0, A0, xi0, xi1, xi2)
      if (Af%val /= 0d0 .or. maxval(abs(Af%d1Array)) /= 0d0) stop 1
      A0 = 2d0
      Af = eval_Af(10d0, A0, xi0, xi1, xi2)
      if (abs(Af%val - 7d0/6d0) > 1d-14) stop 1
      xi0 = 1d0
      xi1 = 0d0
      xi2 = 0d0
      A0 = 0.5d0
      Af = eval_Af(2d0, A0, xi0, xi1, xi2)
      if (Af%val /= 1.5d0) stop 1
   end subroutine check_TDC_degenerate

end program test_turb
