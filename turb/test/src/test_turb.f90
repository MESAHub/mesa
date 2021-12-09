program test_turb
   use math_lib
   use auto_diff
   use const_def
   use turb

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

   subroutine check_TDC()
      real(dp) :: mixing_length_alpha, conv_vel_start, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, scale
      type(auto_diff_real_star_order1) :: &
         r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, gradL
      type(auto_diff_real_star_order1) :: gradT, Y_face, conv_vel, D, Gamma
      integer :: mixing_type, ierr, tdc_num_iters
      real(dp) :: dt1, dt2
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


      write(*,*) "####################################"
      write(*,*) "Running first dt test"

      do j=-20,0
         dt = pow(1.5d0, j)
         call set_TDC(&
            conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, .false., &
            mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, ierr)

         write(*,1) 'dt', dt
         write(*,1) 'gradT', gradT%val
         write(*,1) 'conv_vels', conv_vel_start, conv_vel% val
      end do
stop
      write(*,*) "####################################"
      write(*,*) "Running second dt test"

      dt=9.9999999998979316d-3
      dt2 = dt

      call set_TDC(&
         conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, .true., &
         mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
         scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, ierr)

      write(*,1) 'dt', dt
      write(*,1) 'gradT', gradT%val
      write(*,1) 'conv_vels', conv_vel_start, conv_vel% val

      write(*,*) "check relative difference in dt", (dt1-dt2)/dt1

   end subroutine check_TDC

end program test_turb
