! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module turb_lib

   use const_def, only: dp
   use num_lib
   use math_lib
   use utils_lib
   use auto_diff
   use turb_def
   use thermohaline, only: set_info_HG19, set_info_FRG24 ! Used by plotter routine
   
   implicit none

   private
   public :: set_thermohaline
   public :: set_info_HG19
   public :: set_info_FRG24
   public :: set_mlt
   public :: set_tdc
   public :: set_semiconvection

contains

   !> Computes the diffusivity of thermohaline mixing when the
   !! thermal gradient is stable and the composition gradient is unstable.
   !!
   !! @param thermohaline_option A string specifying which thermohaline prescription to use.
   !! @param grada Adiabatic gradient dlnT/dlnP
   !! @param gradr Radiative temperature gradient dlnT/dlnP, equals the actual gradient because there's no convection
   !! @param N2_T Structure part of brunt squared (excludes composition term)
   !! @param T Temperature
   !! @param rho Density
   !! @param Cp Heat capacity at constant pressure
   !! @param opacity opacity
   !! @param gradL_composition_term dlnMu/dlnP where Mu is the mean molecular weight.
   !! @param XH1 Mass fraction of H1.
   !! @param eta Magnetic diffusivity.
   !! @param iso The index of the species that drives thermohaline mixing.
   !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
   !! @param thermohaline_mag_B Magnetic field strength (guass) for HG19 and FRG24 prescriptions.
   !! @param thermohaline_FRG24_safety Safety parameter for choosing approximations in FRG24 prescription.
   !! @param thermohaline_FRG24_nks Number of vertical wavenumbers to search over in FRG24 prescription.
   !! @param thermohaline_FRG24_N Maximal mode index in FRG24 prescription.
   !! @param D Output, thermohaline diffusivity.
   !! @param gradT Output, temperature gradient.
   !! @param Y_face Output, superadiabaticity at face.
   !! @param conv_vel Output, convective velocity.
   !! @param mixing_type Output, mixing type.
   !! @param ierr Output, error index.
   !! @param th_results Output, detailed thermohaline results (see thermohaline.f90)
   subroutine set_thermohaline(thermohaline_option, Lambda, grada, gradr, N2_T, T, rho, Cp, opacity, &
      gradL_composition_term, XH1, eta, iso, &
      thermohaline_coeff, thermohaline_mag_B, thermohaline_r_ext, &
      thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_N, &
      D, gradT, Y_face, conv_vel, mixing_type, ierr, th_info)

      use thermohaline

      character(*), intent(in)                      :: thermohaline_option
      type(auto_diff_real_star_order1), intent(in)  :: Lambda
      type(auto_diff_real_star_order1), intent(in)  :: grada
      type(auto_diff_real_star_order1), intent(in)  :: gradr
      type(auto_diff_real_star_order1), intent(in)  :: N2_T
      type(auto_diff_real_star_order1), intent(in)  :: T
      type(auto_diff_real_star_order1), intent(in)  :: rho
      type(auto_diff_real_star_order1), intent(in)  :: Cp
      type(auto_diff_real_star_order1), intent(in)  :: opacity
      real(dp), intent(in)                          :: gradL_composition_term
      real(dp), intent(in)                          :: XH1
      real(dp), intent(in)                          :: eta
      integer, intent(in)                           :: iso
      real(dp),intent(in)                           :: thermohaline_coeff
      real(dp), intent(in)                          :: thermohaline_mag_B
      real(dp), intent(in)                          :: thermohaline_r_ext
      integer, intent(in)                           :: thermohaline_FRG24_safety
      integer, intent(in)                           :: thermohaline_FRG24_nks
      integer, intent(in)                           :: thermohaline_FRG24_N
      type(auto_diff_real_star_order1), intent(out) :: D
      type(auto_diff_real_star_order1), intent(out) :: gradT
      type(auto_diff_real_star_order1), intent(out) :: Y_face
      type(auto_diff_real_star_order1), intent(out) :: conv_vel
      integer, intent(out)                          :: mixing_type
      integer, intent(out)                          :: ierr
      type(th_info_t), intent(out), optional        :: th_info

      type(th_info_t) :: th_info_

      call get_thermohaline_info( &
         thermohaline_option, grada%val, gradr%val, N2_T%val, T%val, rho%val, Cp%val, opacity%val, &
         gradL_composition_term, XH1, eta, iso, &
         thermohaline_coeff, thermohaline_mag_B, thermohaline_r_ext, &
         thermohaline_FRG24_safety, thermohaline_FRG24_nks, thermohaline_FRG24_N, &
         th_info_, ierr)

      D = th_info_%D_thrm
      gradT = gradr
      Y_face = gradT - grada
      conv_vel = 3._dp*D/Lambda
      mixing_type = th_info_%mixing_type

      if (PRESENT(th_info)) th_info = th_info_

   end subroutine set_thermohaline

   !> Computes the outputs of time-dependent convection theory following the model specified in
   !! Radek Smolec's thesis [https://users.camk.edu.pl/smolec/phd_smolec.pdf], which in turn
   !! follows the model of Kuhfuss 1986.
   !!
   !! Internally this solves the equation L = L_conv + L_rad.
   !!
   !! @param conv_vel_start The convection speed at the start of the step.
   !! @param mixing_length_alpha The mixing length parameter.
   !! @param alpha_TDC_DAMP TDC turbulent damping parameter
   !! @param alpha_TDC_DAMPR TDC radiative damping parameter
   !! @param alpha_TDC_PtdVdt TDC coefficient on P_turb*dV/dt. Physically should probably be 1.
   !! @param The time-step (s).
   !! @param cgrav gravitational constant (erg*cm/g^2).
   !! @param m Mass inside the face (g).
   !! @param report Write debug output if true, not if false.
   !! @param mixing_type Set to semiconvective if convection operates (output).
   !! @param scale The scale for computing residuals to the luminosity equation (erg/s).
   !! @param chiT dlnP/dlnT|rho
   !! @param chiRho dlnP/dlnRho|T
   !! @param gradr Radiative temperature gradient.
   !! @param r radial coordinate of the face (cm).
   !! @param P pressure (erg/cm^3).
   !! @param T temperature (K).
   !! @param rho density (g/cm^3).
   !! @param dV The change in specific volume of the face (cm^3/g) since the start of the step.
   !! @param Cp Specific heat at constant pressure (erg/g/K).
   !! @param opacity opacity (cm^2/g).
   !! @param scale_height The pressure scale-height (cm).
   !! @param gradL The Ledoux temperature gradient dlnT/dlnP
   !! @param grada The adiabatic temperature gradient dlnT/dlnP|s
   !! @param conv_vel The convection speed (cm/s).
   !! @param D The chemical diffusion coefficient (cm^2/s).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param gradT The temperature gradient dlnT/dlnP (output).
   !! @param tdc_num_iters Number of iterations taken in the TDC solver.
   !! @param ierr Tracks errors (output).
   subroutine set_TDC( &
            conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, report, &
            mixing_type, scale, chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, ierr)
      use tdc
      use tdc_support
      real(dp), intent(in) :: conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, scale
      type(auto_diff_real_star_order1), intent(in) :: &
         chiT, chiRho, gradr, r, P, T, rho, dV, Cp, opacity, scale_height, gradL, grada
      logical, intent(in) :: report
      type(auto_diff_real_star_order1),intent(out) :: conv_vel, Y_face, gradT, D
      integer, intent(out) :: tdc_num_iters, mixing_type, ierr
      type(tdc_info) :: info
      type(auto_diff_real_star_order1) :: L, grav, Lambda, Gamma
      real(dp), parameter :: alpha_c = (1d0/2d0)*sqrt_2_div_3
      real(dp), parameter :: lower_bound_Z = -1d2
      real(dp), parameter :: upper_bound_Z = 1d2
      real(dp), parameter :: eps = 1d-2 ! Threshold in logY for separating multiple solutions.
      type(auto_diff_real_tdc) :: Zub, Zlb
      include 'formats'

      ! Do a call to MLT
      grav = cgrav * m / pow2(r)
      L = 64 * pi * boltz_sigma * pow4(T) * grav * pow2(r) * gradr / (3d0 * P * opacity)
      Lambda = mixing_length_alpha * scale_height
      call set_MLT('Cox', mixing_length_alpha, 0d0, 0d0, &
                     chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                     gradr, grada, gradL, &
                     Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)

      ! Pack TDC info
      info%report = report
      info%mixing_length_alpha = mixing_length_alpha
      info%alpha_TDC_DAMP = alpha_TDC_DAMP
      info%alpha_TDC_DAMPR = alpha_TDC_DAMPR
      info%alpha_TDC_PtdVdt = alpha_TDC_PtdVdt
      info%dt = dt
      info%L = convert(L)
      info%gradL = convert(gradL)
      info%grada = convert(grada)
      info%c0 = convert(mixing_length_alpha*alpha_c*rho*T*Cp*4d0*pi*pow2(r))
      info%L0 = convert((16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T)/(P*opacity)) ! assumes QHSE for dP/dm
      info%A0 = conv_vel_start/sqrt_2_div_3
      info%T = T
      info%rho = rho
      info%dV = dV
      info%Cp = Cp
      info%kap = opacity
      info%Hp = scale_height
      info%Gamma = Gamma

      ! Get solution
      Zub = upper_bound_Z
      Zlb = lower_bound_Z
      call get_TDC_solution(info, scale, Zlb, Zub, conv_vel, Y_face, tdc_num_iters, ierr)

      ! Unpack output
      gradT = Y_face + gradL
      D = conv_vel*scale_height*mixing_length_alpha/3d0     ! diffusion coefficient [cm^2/sec]
      if (conv_vel > 0d0) then
         mixing_type = convective_mixing
      else
         mixing_type = no_mixing
      end if
   end subroutine set_TDC

   !> Calculates the outputs of semiconvective mixing theory.
   !!
   !! @param L Luminosity across a face (erg/s).
   !! @param Lambda The mixing length (cm).
   !! @param m Mass inside the face (g).
   !! @param T temperature (K).
   !! @param P pressure (erg/cm^3).
   !! @param Pr radiation pressure (erg/cm^3).
   !! @param beta ratio of gas pressure to radiation pressure.
   !! @param opacity opacity (cm^2/g).
   !! @param rho density (g/cm^3).
   !! @param alpha_semiconvection The semiconvective alpha parameter.
   !! @param semiconvection_option A string specifying which semiconvection theory to use. Currently supported are 'Langer_85 mixing; gradT = gradr' and 'Langer_85'.
   !! @param cgrav gravitational constant (erg*cm/g^2).
   !! @param Cp Specific heat at constant pressure (erg/g/K).
   !! @param gradr The radiative temperature gradient dlnT/dlnP_{rad}
   !! @param grada The adiabatic temperature gradient dlnT/dlnP|s
   !! @param gradL The Ledoux temperature gradient dlnT/dlnP
   !! @param gradL_composition_term The contribution of composition gradients to the Ledoux temperature gradient.
   !! @param gradT The temperature gradient dlnT/dlnP (output).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param conv_vel The convection speed (cm/s).
   !! @param D The chemical diffusion coefficient (cm^2/s).
   !! @param mixing_type Set to semiconvective if convection operates (output).
   !! @param ierr Tracks errors (output).
   subroutine set_semiconvection(L, Lambda, m, T, P, Pr, beta, opacity, rho, alpha_semiconvection, &
                                 semiconvection_option, cgrav, Cp, gradr, grada, gradL, &
                                 gradL_composition_term, &
                                 gradT, Y_face, conv_vel, D, mixing_type, ierr)
      use semiconvection
      type(auto_diff_real_star_order1), intent(in) :: L, Lambda, T, P, Pr, beta, opacity, rho
      type(auto_diff_real_star_order1), intent(in) :: Cp, gradr, grada, gradL
      character(len=*), intent(in) :: semiconvection_option
      real(dp), intent(in) :: alpha_semiconvection, cgrav, gradL_composition_term, m
      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      call calc_semiconvection(L, Lambda, m, T, P, Pr, beta, opacity, rho, alpha_semiconvection, &
                                 semiconvection_option, cgrav, Cp, gradr, grada, gradL, &
                                 gradL_composition_term, &
                                 gradT, Y_face, conv_vel, D, mixing_type, ierr)
   end subroutine set_semiconvection

   !> Calculates the outputs of convective mixing length theory.
   !!
   !! @param MLT_option A string specifying which MLT option to use. Currently supported are Cox, Henyey, ML1, ML2, Mihalas. Note that 'TDC' is also a valid input and will return the Cox result. This is for use when falling back from TDC -> MLT, as Cox is the most-similar prescription to TDC.
   !! @param mixing_length_alpha The mixing length parameter.
   !! @param Henyey_MLT_nu_param The nu parameter in Henyey's MLT prescription.
   !! @param Henyey_MLT_y_param The y parameter in Henyey's MLT prescription.
   !! @param chiT dlnP/dlnT|rho
   !! @param chiRho dlnP/dlnRho|T
   !! @param Cp Specific heat at constant pressure (erg/g/K).
   !! @param grav The acceleration due to gravity (cm/s^2).
   !! @param Lambda The mixing length (cm).
   !! @param rho density (g/cm^3).
   !! @param T temperature (K).
   !! @param opacity opacity (cm^2/g)
   !! @param gradr The radiative temperature gradient dlnT/dlnP_{rad}
   !! @param grada The adiabatic temperature gradient dlnT/dlnP|s
   !! @param gradL The Ledoux temperature gradient dlnT/dlnP
   !! @param Gamma The convective efficiency parameter (output).
   !! @param gradT The temperature gradient dlnT/dlnP (output).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param conv_vel The convection speed (cm/s).
   !! @param D The chemical diffusion coefficient (cm^2/s).
   !! @param mixing_type Set to convective if convection operates (output).
   !! @param ierr Tracks errors (output).
   subroutine set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                     chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                     gradr, grada, gradL, &
                     Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)
      use mlt
      type(auto_diff_real_star_order1), intent(in) :: chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=*), intent(in) :: MLT_option
      real(dp), intent(in) :: mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param
      type(auto_diff_real_star_order1), intent(out) :: Gamma, gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      call calc_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                     chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                     gradr, grada, gradL, &
                     Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)
   end subroutine set_MLT

end module turb_lib
