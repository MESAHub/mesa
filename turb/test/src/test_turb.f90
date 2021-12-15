program test_turb
   use math_lib
   use auto_diff
   use const_def
   use turb

   call check_efficient_MLT_scaling()

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

end program test_turb
