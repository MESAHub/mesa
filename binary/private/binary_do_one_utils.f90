! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton, Pablo Marchant & The MESA Team
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

      module binary_do_one_utils
      
      use binary_def
      use const_def
      use math_lib

      implicit none
      
      contains
      
      subroutine write_binary_terminal_header(b)
         type (binary_info), pointer :: b
         if (b% model_number <= b% recent_binary_log_header) return
         if (b% just_wrote_binary_terminal_header) return
         b% recent_binary_log_header = b% model_number
         call do_show_binary_terminal_header(b)
         b% just_wrote_binary_terminal_header = .true.
      end subroutine write_binary_terminal_header
      
      
      subroutine do_show_binary_log_description(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call get_binary_ptr(id, b, ierr)
         if (ierr /= 0) return
         write(*,*)
         write(*,'(a)') " The binary terminal output contains the following information"
         write(*,*)
         write(*,'(a)') "      'step' is the number of steps since the start of the run,"
         write(*,'(a)') "      'lg_dt' is log10 timestep in years,"
         write(*,'(a)') "      'age_yr' is the simulated years since the start run,"
         write(*,'(a)') "      'M1+M2' is the total mass of the system (Msun),"
         write(*,'(a)') "      'M1' is the mass of the primary (Msun)"
         write(*,'(a)') "      'M2' is the mass of the secondary (Msun)"
         write(*,'(a)') "      'separ' is the semi-major axis of the orbit (Rsun),"
         write(*,'(a)') "      'R1' is the radius of the primary (Rsun)"
         write(*,'(a)') "      'R2' is the radius of the secondary (Rsun)"
         write(*,'(a)') "      'Porb' is the orbital period (days),"
         write(*,'(a)') "      'P1' is the rotation period of star 1 (days, zero if not modeling rotation),"
         write(*,'(a)') "      'P2' is the rotation period of star 2 (days, zero if not modeling rotation),"
         write(*,'(a)') "      'e' orbital eccentricity,"
         write(*,'(a)') "      'dot_e' time derivative of e (1/yr),"
         write(*,'(a)') "      'Eorb' orbital energy G*M1*M2/2*separation (ergs),"
         write(*,'(a)') "      'M2/M1' mass ratio,"
         write(*,'(a)') "      'vorb1' orbital velocity of star 1 (km/s),"
         write(*,'(a)') "      'vorb2' orbital velocity of star 2 (km/s),"
         write(*,'(a)') "      'pm_i' index of star evolved as point mass, zero if both stars are modeled,"
         write(*,'(a)') "      'RL1' Roche lobe radius of star 1 (Rsun),"
         write(*,'(a)') "      'Rl2' Roche lobe radius of star 2 (Rsun),"
         write(*,'(a)') "      'donor_i' index of star taken as donor,"
         write(*,'(a)') "      'RL_gap1' (R1-Rl1)/Rl1,"
         write(*,'(a)') "      'RL_gap2' (R2-Rl2)/Rl2,"
         write(*,'(a)') "      'dot_Mmt', mass transfer rate (Msun/yr),"
         write(*,'(a)') "      'dot_M1', time derivative for the mass of star 1 (Msun/yr),"
         write(*,'(a)') "      'dot_M2', time derivative for the mass of star 2 (Msun/yr),"
         write(*,'(a)') "      'eff', mass transfer efficiency, computed as -dot_M2/dot_M1 (zero if dot_M1=0),"
         write(*,'(a)') "      'dot_Medd', Eddington accretion rate (Msun/yr),"
         write(*,'(a)') "      'L_acc', accretion luminosity when accreting to a point mass (ergs/s),"
         write(*,'(a)') "      'Jorb', orbital angular momentum (g*cm^2/s)"
         write(*,'(a)') "      'spin1', spin angular momentum of star 1 (g*cm^2/s),"
         write(*,'(a)') "      'spin2', spin angular momentum of star 2 (g*cm^2/s),"
         write(*,'(a)') "      'dot_J', time derivative of Jorb (g*cm^2/s^2),"
         write(*,'(a)') "      'dot_Jgr', time derivative of Jorb due to gravitational waves (g*cm^2/s^2),"
         write(*,'(a)') "      'dot_Jml', time derivative of Jorb due to mass loss (g*cm^2/s^2),"
         write(*,'(a)') "      'dot_Jmb', time derivative of Jorb due to magnetic braking (g*cm^2/s^2),"
         write(*,'(a)') "      'dot_Jls', time derivative of Jorb due to spin-orbit coupling (g*cm^2/s^2),"
         write(*,'(a)') "      'rlo_iters', number of iterations for implicit calculation of mass transfer,"
         write(*,*)
         write(*,'(a)') " All this and more can be saved in binary_history.data during the run."
      end subroutine do_show_binary_log_description
      
      
      subroutine do_show_binary_terminal_header(b)
         type (binary_info), pointer :: b
         call output_binary_terminal_header(b,terminal_iounit)
         if (b% extra_binary_terminal_iounit > 0) &
            call output_binary_terminal_header(b,b% extra_binary_terminal_iounit)
      end subroutine do_show_binary_terminal_header
      
      
      subroutine output_binary_terminal_header(b,io)
         type (binary_info), pointer :: b
         integer, intent(in) :: io
         write(io,'(a)') &
            '_______________________________________________________________________' // &
            '___________________________________________________________________________'
         write(io,*)
         write(io,'(a)') &
            'binary_step      M1+M2      separ       Porb          e      M2/M1' // &
            '       pm_i    donor_i    dot_Mmt        eff       Jorb      dot_J    dot_Jmb'
         write(io,'(a)') &
            '      lg_dt         M1         R1         P1      dot_e      vorb1' // &
            '        RL1    Rl_gap1     dot_M1   dot_Medd      spin1    dot_Jgr    dot_Jls'
         write(io,'(a)') &
            '     age_yr         M2         R2         P2       Eorb      vorb2' // &
            '        RL2    Rl_gap2     dot_M2      L_acc      spin2    dot_Jml  rlo_iters'
         write(io,'(a)') &
            '_______________________________________________________________________' // &
            '___________________________________________________________________________'
         write(io,*)
         
      end subroutine output_binary_terminal_header
      
      
      subroutine do_binary_terminal_summary(b)
         type (binary_info), pointer :: b
         call output_binary_terminal_summary(b,terminal_iounit)
         if (b% extra_binary_terminal_iounit > 0) then
            call output_binary_terminal_summary(b,b% extra_binary_terminal_iounit)
            flush(b% extra_binary_terminal_iounit)
         end if
      end subroutine do_binary_terminal_summary
      
      
      subroutine output_binary_terminal_summary(b,io)
         type (binary_info), pointer :: b
         integer, intent(in) :: io
         
         real(dp) :: age, time_step, total_mass
         real(dp) :: Eorb, vorb1, vorb2, dot_M1, dot_M2, eff, dot_Medd, spin1, spin2, P1, P2
         integer :: model, ierr, rlo_iters
         character (len=90) :: fmt, fmt1, fmt2, fmt3, fmt4, fmt5, fmt6
         
         include 'formats'
         
         age = b% binary_age
         time_step = b% time_step
         model = b% model_number
         rlo_iters = b% num_tries

         total_mass = (b% m(1) + b% m(2))/Msun

         vorb1 = 2.0d0 * pi * b% m(2)/(b% m(1) + b% m(2)) * b% separation / b% period / 1.0d5
         vorb2 = 2.0d0 * pi * b% m(1)/(b% m(1) + b% m(2)) * b% separation / b% period / 1.0d5
         Eorb = -standard_cgrav * b% m(1) * b% m(2) / (2*b% separation)
         dot_M1 = b% component_mdot(1)/Msun*secyer
         if (abs(dot_M1) < 1d-99) dot_M1 = 0d0
         dot_M2 = b% component_mdot(2)/Msun*secyer
         if (abs(dot_M2) < 1d-99) dot_M2 = 0d0
         if (b% component_mdot(b% d_i) == 0d0) then
            eff = 1d0
         else
            eff = -b% component_mdot(b% a_i)/b% component_mdot(b% d_i)
         end if
         if (b% limit_retention_by_mdot_edd) then
            dot_Medd = b% mdot_edd/Msun*secyer
         else
            dot_Medd = 1d99
         end if
         if (b% point_mass_i /= 1) then
            spin1 = b% s1% total_angular_momentum
         else
            spin1 = 0d0
         end if

         if (b% point_mass_i /= 2) then
            spin2 = b% s2% total_angular_momentum
         else
            if (.not. b% model_twins_flag) then
               spin2 = 0d0
            else
               spin2 = b% s1% total_angular_momentum
            end if
         end if

         P1 = 0d0
         if (b% point_mass_i /= 1) then
            if (b% s1% rotation_flag) then
               if (abs(b% s1% omega_avg_surf) > 0) then
                  P1 = 2*pi/(b% s1% omega_avg_surf*24d0*3600d0)
               end if
            end if
         end if

         P2 = 0d0
         if (b% point_mass_i /= 2)then 
            if (b% s2% rotation_flag) then
               if (abs(b% s2% omega_avg_surf) > 0) then
                  P2 = 2*pi/(b% s2% omega_avg_surf*24d0*3600d0)
               end if
            end if
         else if (b% model_twins_flag) then
            if (b% s1% rotation_flag) then
               if (abs(b% s1% omega_avg_surf) > 0) then
                  P2 = 2*pi/(b% s1% omega_avg_surf*24d0*3600d0)
               end if
            end if
         end if

         ierr = 0         

         !make format strings for first line of output
         !step and m1+m2
         if (total_mass >= 1d2) then
            fmt1 = '(a3,i8,1pe11.3,0p,'
         else
            fmt1 = '(a3,i8,f11.6,'
         end if
         !separation
         if (b% separation/Rsun >= 1d2) then
            fmt2 = '1pe11.3,0p,'
         else
            fmt2 = 'f11.6,'
         end if
         !Porb
         if (b% period/(3600d0*24d0) >= 1d2) then
            fmt3 = '1pe11.3,0p,'
         else
            fmt3 = 'f11.6,'
         end if
         !eccentricity
         if (b% eccentricity >= 1d-2) then
            fmt4 = 'f11.6,'
         else
            fmt4 = '1pe11.3,0p,'
         end if
         !mass ratio, pm_i, donor_index and dot_Mmt
         if (b% m(2)/b% m(1) >= 1d-2) then
            fmt5 = 'f11.6,2i11,1pe11.3,0p,'
         else
            fmt5 = '1pe11.3,0p,2i11,1pe11.3,0p,'
         end if
         !eff
         if (eff >= 1d-2) then
            fmt6 = 'f11.6,3(1pe11.3))'
         else
            fmt6 = '4(1pe11.3))'
         end if

         
         fmt = trim(fmt1) // trim(fmt2) // trim(fmt3) // trim(fmt4) // trim(fmt5) // trim(fmt6)
         write(io,fmt=fmt) &
            'bin', &
            model, &
            total_mass, &
            b% separation / Rsun, &
            b% period / (3600d0*24d0), &
            b% eccentricity, &
            b% m(2)/b% m(1), &            
            b% point_mass_i, &
            b% d_i, &
            b% step_mtransfer_rate/Msun*secyer, &
            eff, &
            b% angular_momentum_j, &
            b% jdot, &
            b% jdot_mb

         !make format strings for second line of output
         !lg_dt and M1
         if (b% m(1)/Msun >= 1d2) then
            fmt1 = '(f11.6,1pe11.3,0p,'
         else
            fmt1 = '(2f11.6,'
         end if
         !R1
         if (b% r(1)/Rsun >= 1d2) then
            fmt2 = '1pe11.3,0p,'
         else
            fmt2 = 'f11.6,'
         end if
         !P1, dot_e
         if (P1 >= 1d2) then
            fmt3 = '2(1pe11.3),0p'
         else
            fmt3 = 'f11.6,1pe11.3,0p'
         end if
         !vorb1
         if (vorb1 >= 1d-2 .and. vorb1 < 1d4) then
            fmt4 = 'f11.6,'
         else
            fmt4 = '1pe11.3,0p,'
         end if
         !RL1 and the rest
         if (b% rl(1)/Rsun >= 1d2) then
            fmt5 = '7(1pe11.3))'
         else
            fmt5 = 'f11.6,6(1pe11.3))'
         end if
         
         fmt = trim(fmt1) // trim(fmt2) // trim(fmt3) // trim(fmt4) // trim(fmt5)
         write(io,fmt=fmt) &
            safe_log10(time_step),  &
            b% m(1)/Msun, &
            b% r(1)/Rsun, &
            P1, &
            b% edot, &
            vorb1, &
            b% rl(1)/Rsun, &
            b% rl_relative_gap(1), &
            dot_M1, &
            dot_Medd, &
            spin1, &
            b% jdot_gr, &
            b% jdot_ls

         !make format strings for third line of output
         !age_yr and M2
         if (b% m(2)/Msun >= 1d2) then
            fmt1 = '(1pe11.4,0p,1pe11.3,0p,'
         else
            fmt1 = '(1pe11.4,0p,f11.6,'
         end if
         !R2
         if (b% r(2)/Rsun >= 1d2) then
            fmt2 = '1pe11.3,0p,'
         else
            fmt2 = 'f11.6,'
         end if
         !P2, Eorb
         if (P2 >= 1d2) then
            fmt3 = '2(1pe11.3),0p'
         else
            fmt3 = 'f11.6,1pe11.3,0p'
         end if
         !vorb2
         if (vorb2 >= 1d-2 .and. vorb2 < 1d4) then
            fmt4 = 'f11.6,'
         else
            fmt4 = '1pe11.3,0p,'
         end if
         !RL2 and the rest
         if (b% rl(2)/Rsun >= 1d2) then
            fmt5 = '6(1pe11.3),i11)'
         else
            fmt5 = 'f11.6,5(1pe11.3),0p,i11)'
         end if
         
         fmt = trim(fmt1) // trim(fmt2) // trim(fmt3) // trim(fmt4) // trim(fmt5)
         write(io,fmt=fmt) &
            age,  &
            b% m(2)/Msun, &
            b% r(2)/Rsun, &
            P2, &
            Eorb, &
            vorb2, &
            b% rl(2)/Rsun, &
            b% rl_relative_gap(2), &
            dot_M2, &
            b% accretion_luminosity, &
            spin2, &
            b% jdot_ml, &
            rlo_iters

         write(io,*)
         
         b% just_wrote_binary_terminal_header = .false.
         
      end subroutine output_binary_terminal_summary

      
      end module binary_do_one_utils
      
