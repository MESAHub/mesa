program test_turb
   use math_lib
   use auto_diff
   use const_def
   use turb

   call check_efficient_MLT_scaling()
   call check_TDC()

   contains

   subroutine header(text)
      character(len=*), intent(in) :: text

      write(*,'(a)') ' ----------------------------------------------------------------'
      write(*,'(a)') ' '
      write(*,'(a)') ' '//text
      write(*,'(a)') ' '
      write(*,'(a)') ' ----------------------------------------------------------------'

   end subroutine header


   subroutine check_efficient_MLT_scaling()
      type(auto_diff_real_star_order1) :: chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=3) :: MLT_option = 'Cox'
      real(dp) :: mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param
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
      Cp = 2.5d0 * kerg / mp
      P = rho * T * kerg / mp
      Lambda = P / (rho * grav)
      opacity = 1d0
      grada = 0.4d0
      gradL = grada

      L = 1d5*Lsun
      gradr = 3d0 * P * opacity * L / (64 * pi * boltz_sigma * pow4(T) * grav * pow2(r))

      call set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                        chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                        gradr, grada, gradL, &
                        Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)

      write(*,1) 'vc at 1d5 Lsun',conv_vel%val

      L = 1d8*Lsun
      gradr = 3d0 * P * opacity * L / (64 * pi * boltz_sigma * pow4(T) * grav * pow2(r))

      call set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                        chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                        gradr, grada, gradL, &
                        Gamma, gradT, Y_face, conv_vel2, D, mixing_type, ierr)

      write(*,1) 'vc at 1d8 Lsun',conv_vel2%val

      write(*,1) 'Ratio:',conv_vel2%val/conv_vel%val
      write(*,'(a)') 'Expected ~10 because in the efficient limit vc ~ L^{1/3}'
   end subroutine check_efficient_MLT_scaling



   subroutine check_TDC()
      real(dp) :: mixing_length_alpha, conv_vel_start, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, scale
      type(auto_diff_real_star_order1) :: &
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, gradL
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Gamma
      integer :: mixing_type, ierr, tdc_num_iters
      logical :: report
      integer :: j

      include 'formats'

      call header('Test TDC')

      conv_vel_start = 52320587.415154047

      mixing_length_alpha=2.0000000000000000
      alpha_TDC_DAMP=1.0000000000000000
      alpha_TDC_DAMPR=0.0000000000000000
      alpha_TDC_PtdVdt=0.0000000000000000
      cgrav=6.6743000000000004d-8
      m=5.8707400456875664d34
      scale=5.0386519362246294d45
      L=1.0941528815883500015d0*(-5.0386519362246294d45)
      r=10314294541.567163d0
      P=5.0581587249808894d20
      T=613193666.51783681d0
      rho=5204.5732574745753d0
      dV=3.8256494463482604d-7
      Cp=6628075118.4606590d0
      opacity=9.0750171231469945d-2
      scale_height=2638686602.0063782d0
      gradL=0.25207587267343501d0
      grada=0.25204697256872738d0
      report = .false.


      write(*,*) "####################################"
      write(*,*) "Running dt test"

      do j=0,30
         dt = 500d0*pow(1.02d0,j)
         call set_TDC(&
            conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, report, &
            mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, ierr)

         write(*,1) 'dt, gradT, conv_vel_start, conv_vel', dt, gradT%val, conv_vel_start, conv_vel% val
         if (report) stop
      end do

   end subroutine check_TDC

end program test_turb
