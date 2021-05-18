! ***********************************************************************
!
!   Copyright (C) 2010-2019  Joris Vos & The MESA Team
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

    module binary_edot
    
    use const_def
    use star_lib
    use star_def
    use binary_def
    use binary_tides
    use math_lib

    implicit none

    contains

    real(dp) function get_edot(b) result(edot)
       type (binary_info), pointer :: b
       integer :: ierr
       ierr = 0

       if (b% eccentricity == 0d0) then
           b% edot_tidal = 0d0
           b% edot_enhance = 0d0
           b% extra_edot = 0d0
           b% edot = 0d0
       else
           ! tidal circularisation
           if (b% do_tidal_circ) then
              if (.not. b% use_other_edot_tidal) then
                 call edot_tidal(b% binary_id, ierr)
              else
                 call b% other_edot_tidal(b% binary_id, ierr)
              end if
              if (ierr /= 0) then
                 write(*,*) 'ERROR while calculating edot_tidal'
                 return
              end if
           else
              b% edot_tidal = 0d0
           end if
           
           if (b% edot_tidal < -b% max_abs_edot_tidal) then
              b% edot_tidal = -b% max_abs_edot_tidal
           end if
           
           ! eccentricity enhancement
           if (b% use_eccentricity_enhancement) then
              if (.not. b% use_other_edot_enhance) then
                 call edot_enhancement_Isotropic(b% binary_id, ierr)
              else
                 call b% other_edot_enhance(b% binary_id, ierr)
              end if
              if (ierr /= 0) then
                 write(*,*) 'ERROR while calculating edot_enhance'
                 return
              end if
           else
              b% edot_enhance = 0d0
           end if
           
           if (b% edot_enhance > b% max_abs_edot_enhance) then
              b% edot_enhance = b% max_abs_edot_enhance
           end if
           
           ! user defined eccentricity changes
           if (b% use_other_extra_edot) then
              call b% other_extra_edot(b% binary_id, ierr)
              if (ierr /= 0) then
                 write(*,*) 'ERROR while calculating extra_edot'
                 return
              end if
           else
              b% extra_edot = 0d0
           end if
           
           b% edot = b% edot_tidal + b% edot_enhance + b% extra_edot
           
       end if

       edot = b% edot

    end function get_edot

    ! ==========================================
    ! edot TIDAL
    ! ==========================================
    subroutine edot_tidal(binary_id, ierr)
       integer, intent(in) :: binary_id
       integer, intent(out) :: ierr
       type (binary_info), pointer :: b

       ierr = 0
       call binary_ptr(binary_id, b, ierr)
       if (ierr /= 0) then
          write(*,*) 'failed in binary_ptr'
          return
       end if
       
       b% edot_tidal = 0d0

       if (b% point_mass_i /= 1) then
          if (b% circ_type_1 == "Hut_conv") then
             b% edot_tidal = edot_tidal_Hut(b, b% s1, .true.)
          else if (b% circ_type_1 == "Hut_rad") then
             b% edot_tidal = edot_tidal_Hut(b, b% s1, .false.)
          else
             write(*,*) "Unrecognized circ_type_1", b% circ_type_1
          end if
       end if
       if (b% point_mass_i /= 2) then
          if (b% circ_type_2 == "Hut_conv") then
             b% edot_tidal = b% edot_tidal + edot_tidal_Hut(b, b% s2, .true.)
          else if (b% circ_type_2 == "Hut_rad") then
             b% edot_tidal = b% edot_tidal + edot_tidal_Hut(b, b% s2, .false.)
          else
             write(*,*) "Unrecognized circ_type_2", b% circ_type_2
          end if
       end if

       if (b% model_twins_flag) then
          b% edot_tidal = b% edot_tidal + b% edot_tidal
       end if
    end subroutine edot_tidal

    real(dp) function edot_tidal_Hut(b, s , has_convective_envelope) result(edot_tidal)
       type (binary_info), pointer :: b
       type (star_info), pointer :: s
       logical, intent(in) :: has_convective_envelope
       real(dp) :: m, porb, r_phot, osep, qratio, omega_s, omega_sync

       edot_tidal = 0d0

       porb = b% period
       omega_sync = 2d0*pi/b% period
       omega_s = s% omega_avg_surf
       osep = b% separation

       qratio = b% m(b% a_i) / b% m(b% d_i)
       if (is_donor(b, s)) then
          m = b% m(b% d_i)
          r_phot = b% r(b% d_i)
       else
          qratio = 1.0d0/qratio
          m = b% m(b% a_i)
          r_phot = b% r(b% a_i)
       end if

       ! eq. (10) of Hut, P. 1981, A&A, 99, 126
       edot_tidal = -27.0d0*qratio*(1+qratio)*pow8(r_phot/osep) &
           * b% eccentricity*pow(1-pow2(b% eccentricity),-6.5d0)*b% Ftid_1
       ! add multiplication by (k/T), eq. (29) of Hurley et al. 2002
       edot_tidal = edot_tidal*k_div_T(b, s, has_convective_envelope)
       ! add terms dependant on omega
       edot_tidal = edot_tidal*(f3(b% eccentricity) - &
           11d0/18d0 * omega_s / omega_sync * f4(b% eccentricity) * &
           pow(1-pow2(b% eccentricity),1.5d0))
    
    end function edot_tidal_Hut
    
    ! ==========================================
    ! Edot MASS LOSS
    ! ==========================================
    
    subroutine edot_enhancement_Isotropic(binary_id, ierr)
       integer, intent(in) :: binary_id
       integer, intent(out) :: ierr
       type (binary_info), pointer :: b
       integer :: i
       real(dp) :: de, Mtot, xfer, costh

       ierr = 0
       call binary_ptr(binary_id, b, ierr)
       if (ierr /= 0) then
          write(*,*) 'failed in binary_ptr'
          return
       end if

       b% edot_enhance = 0d0

       ! cos_cr isn't vectorised, so we have to do this in a loop
       do i = 1, b% anomaly_steps
          costh = cos(b% theta_co(i))
       
          b% e1(i) = b% eccentricity + costh
          b% e2(i) = 2d0*costh + b% eccentricity*(1d0 + costh*costh)
          b% e3(i) = b% eccentricity*(1d0-costh*costh)  ! = b% eccentricity*sin(b% theta_co)**2
       end do
       
!        xfer = min(b% wind_xfer_fraction, b% xfer_fraction)
       Mtot = b% m(1) + b% m(2) ! total mass in gr
       
       b% edot_theta = - b% mdot_donor_theta / Mtot * b% e1 !-&
!               b% mdot_donor_theta * xfer / b% m(b% a_i) * (b% m(b% d_i) / Mtot *&
!               ((b% m(b% a_i)**2 / b% m(b% d_i)**2 - 1 ) * e2 - e3 ))
       
       !integrate to get total eccentricity enhancement
       de = 0d0
       do i = 2,b% anomaly_steps ! trapezoidal integration
          de = de + 0.5d0 * (b% edot_theta(i-1) + b% edot_theta(i)) * (b% time_co(i) - b% time_co(i-1)) 
       end do
       
       b% edot_enhance = de
    
    end subroutine edot_enhancement_Isotropic
    

    end module binary_edot

