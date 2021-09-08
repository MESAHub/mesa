! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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


module MLT

use star_private_def
use const_def
use num_lib
use utils_lib
use auto_diff_support
use star_utils

implicit none

private
public :: set_MLT

contains

   subroutine set_MLT(MLT_option, mixing_length_alpha, report, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                     chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                     gradr, grada, gradL, k, &
                     Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)
      type(auto_diff_real_star_order1), intent(in) :: chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=*), intent(in) :: MLT_option
      real(dp), intent(in) :: mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param
      integer, intent(in) :: k
      logical, intent(in) :: report
      type(auto_diff_real_star_order1), intent(out) :: Gamma, gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      real(dp) :: ff1, ff2, ff3, ff4
      type(auto_diff_real_star_order1) :: &
         Q, omega, a0, ff4_omega2_plus_1, A_1, A_2, &
         A_numerator, A_denom, A, Bcubed, delta, Zeta, &
         f, f0, f1, f2, radiative_conductivity, convective_conductivity
      include 'formats' 
      if (gradr > gradL) then
         ! Convection zone

         Q = chiT/chiRho ! 'Q' param  C&G 14.24
         if (MLT_option == 'Cox' .or. MLT_option == 'TDC') then ! this assumes optically thick
            a0 = 9d0/4d0
            convective_conductivity = &
               Cp*grav*pow2(Lambda)*rho*(sqrt(Q*rho/(2d0*P)))/9d0 ! erg / (K cm sec)
            radiative_conductivity = &
               (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
            A = convective_conductivity / radiative_conductivity !  unitless.
         else
            select case(trim(MLT_option))
            case ('Henyey')
               ff1=1.0d0/Henyey_MLT_nu_param
               ff2=0.5d0 
               ff3=8.0d0/Henyey_MLT_y_param
               ff4=1.0d0/Henyey_MLT_y_param
            case ('ML1')
               ff1=0.125d0 
               ff2=0.5d0 
               ff3=24.0d0
               ff4=0.0d0
            case ('ML2')
               ff1=1.0d0
               ff2=2.0d0
               ff3=16.0d0
               ff4=0.0d0
            case ('Mihalas')
               ff1=0.125d0 
               ff2=0.5d0 
               ff3=16.0d0
               ff4=2.0d0
            case default
               write(*,'(3a)') 'Error: ', trim(MLT_option), ' is not an allowed MLT option'
               call mesa_error(__FILE__,__LINE__)
            end select

            omega = Lambda*rho*opacity
            ff4_omega2_plus_1 = ff4/pow2(omega) + 1d0
            a0 = (3d0/16d0)*ff2*ff3/ff4_omega2_plus_1
            A_1 = 4d0*Cp*sqrt(ff1*P*Q*rho)
            A_2 = mixing_length_alpha*omega*ff4_omega2_plus_1
            A_numerator = A_1*A_2
            A_denom = ff3*crad*clight*pow3(T)
            A = A_numerator/A_denom   
         end if 

         ! 'B' param  C&G 14.81
         Bcubed = (pow2(A)/a0)*(gradr - gradL)   

         ! now solve cubic equation for convective efficiency, Gamma
         ! a0*Gamma^3 + Gamma^2 + Gamma - a0*Bcubed == 0   C&G 14.82,
         ! leave it to Mathematica to find an expression for the root we want      
         delta = a0*Bcubed               
         f = -2d0 + 9d0*a0 + 27d0*a0*a0*delta
         if (f > 1d100) then
            f0 = f
         else
            f0 = pow2(f) + 4d0*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)
            if (f0 <= 0d0) then
               f0 = f
            else
               f0 = sqrt(f0)         
            end if
         end if
         f1 = -2d0 + 9d0*a0 + 27d0*a0*a0*delta + f0  
         if (f1 <= 0d0) return
         f1 = pow(f1,one_third)     
         f2 = 2d0*two_13*(1d0 - 3d0*a0) / f1       
         Gamma = (four_13*f1 + f2 - 2d0) / (6d0*a0)

         if (Gamma <= 0d0) return

         ! average convection velocity   C&G 14.86b
         conv_vel = mixing_length_alpha*sqrt(Q*P/(8d0*rho))*Gamma / A
         D = conv_vel*Lambda/3d0     ! diffusion coefficient [cm^2/sec]

         !Zeta = pow3(Gamma)/Bcubed  ! C&G 14.80
         Zeta = exp(3d0*log(Gamma) - log(Bcubed)) ! write it this way to avoid overflow problems
      else
         ! Radiative zone, because this means that gradr < gradL
         Zeta = 0d0
         conv_vel = 0d0
         D = 0d0
      end if

      ! Zeta must be >= 0 and <= 1.
      ! By construction (above) it cannot be less than zero,
      ! so we just check that it is a valid number and is not greater
      ! than one.
      if (is_bad(Zeta%val)) return
      if (Zeta > 1d0) then
         Zeta = 1d0
      end if            
      
      gradT = (1d0 - Zeta)*gradr + Zeta*gradL ! C&G 14.79      
      Y_face = gradT - gradL
      
      if (Y_face > 0d0) then
         mixing_type = convective_mixing
      end if

      if (report) then
         write(*,2) 'set_MLT val for Zeta gradr grada gradT Y_face', k, &
            Zeta%val, gradr%val, grada%val, gradT%val, Y_face%val
         write(*,2) 'set_MLT d_dlnd_00 for Zeta gradr grada gradT Y_face', k, &
            Zeta%d1Array(i_lnd_00), gradr%d1Array(i_lnd_00), &
            grada%d1Array(i_lnd_00), gradT%d1Array(i_lnd_00), &
            Y_face%d1Array(i_lnd_00)
      end if

   end subroutine set_MLT   

end module MLT