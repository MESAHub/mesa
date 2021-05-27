! ***********************************************************************
!
!   Copyright (C) 2011-2019  The MESA Team
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

      module write_model

      use star_private_def
      use const_def
      use read_model

      implicit none

      private
      public :: do_write_model


      contains


      subroutine do_write_model(id, filename, ierr)
         use utils_lib
         use chem_def
         use star_utils, only: get_scale_height_face_val
         integer, intent(in) :: id
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         integer :: iounit, i, k, nvar_hydro, nz, species, file_type
         integer, pointer :: chem_id(:)
         type (star_info), pointer :: s
         logical :: v_flag, RTI_flag, &
            RSP2_flag, u_flag, prev_flag, rotation_flag, &
            write_mlt_vc, RSP_flag

         1 format(a32, 2x, 1pd26.16)
         11 format(a32, 2x, 1pd26.16, 2x, a, 2x, 99(1pd26.16))
         2 format(a32, 2x, i30)
         3 format(a32, 3x, a8)
         4 format(a32, 3x, a)

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         chem_id => s% chem_id
         nvar_hydro = s% nvar_hydro
         nz = s% nz
         v_flag = s% v_flag
         u_flag = s% u_flag
         RTI_flag = s% RTI_flag
         rotation_flag = s% rotation_flag
         RSP_flag = s% RSP_flag
         RSP2_flag = s% RSP2_flag 
         write_mlt_vc = s% have_mlt_vc
         
         species = s% species
         
         open(newunit=iounit, file=trim(filename), action='write', status='replace')
         write(iounit,'(a)') '! note: initial lines of file can contain comments'
         write(iounit,'(a)') '!'
         prev_flag = (s% nz_old == s% nz .and. s% generations > 1)
         file_type = 0
         if (RTI_flag) file_type = file_type + 2**bit_for_RTI
         if (prev_flag) file_type = file_type + 2**bit_for_2models
         if (v_flag) file_type = file_type + 2**bit_for_velocity
         if (u_flag) file_type = file_type + 2**bit_for_u
         if (rotation_flag) file_type = file_type + 2**bit_for_rotation
         if (rotation_flag) file_type = file_type + 2**bit_for_j_rot
         if (RSP_flag) file_type = file_type + 2**bit_for_RSP
         if (RSP2_flag) file_type = file_type + 2**bit_for_RSP2
         if (write_mlt_vc) file_type = file_type + 2**bit_for_mlt_vc
                  
         write(iounit, '(i14)', advance='no') file_type
         write(iounit,'(a)',advance='no') ' -- model for mesa/star'
         if (BTEST(file_type, bit_for_velocity)) &
            write(iounit,'(a)',advance='no') ', cell boundary velocities (v)'
         if (BTEST(file_type, bit_for_rotation)) &
            write(iounit,'(a)',advance='no') ', angular velocities (omega)'
         if (BTEST(file_type, bit_for_j_rot)) &
            write(iounit,'(a)',advance='no') ', specific angular momentum (j_rot)'
         if (BTEST(file_type, bit_for_D_omega)) &
            write(iounit,'(a)',advance='no') ', omega diffusion coefficients (D_omega)'
         if (BTEST(file_type, bit_for_am_nu_rot)) &
            write(iounit,'(a)',advance='no') ', am_nu_rot diffusion coefficients'
         if (BTEST(file_type, bit_for_u)) &
            write(iounit,'(a)',advance='no') ', cell center Riemann velocities (u)'
         if (BTEST(file_type, bit_for_RTI)) &
            write(iounit,'(a)',advance='no') ', Rayleigh-Taylor instabilities (alpha_RTI)'
         if (BTEST(file_type, bit_for_mlt_vc)) &
            write(iounit,'(a)',advance='no') ', mlt convection velocity (mlt_vc)'
         if (BTEST(file_type, bit_for_RSP)) &
            write(iounit,'(a)',advance='no') ', RSP luminosity (L), turbulent energy (et_RSP), and radiative flux (erad_RSP)'
         if (BTEST(file_type, bit_for_RSP2)) &
            write(iounit,'(a)',advance='no') ', RSP2 turbulent energy (w^2) and pressure scale height (Hp)'
         write(iounit,'(a)',advance='no') &
            '. cgs units. lnd=ln(density), lnT=ln(temperature), lnR=ln(radius)'
         if (s% i_lum /= 0) then
            write(iounit,'(a)',advance='no') ', L=luminosity'
         end if
         if (s% M_center /= 0) then
            write(iounit,'(a)',advance='no') &
               ', dq=fraction of xmstar=(mstar-mcenter) in cell; remaining cols are mass fractions.'
         else
            write(iounit,'(a)',advance='no') &
               ', dq=fraction of total mstar in cell; remaining cols are mass fractions.'
         end if
         write(iounit,*)
         ! write property list
         write(iounit, '(a)') ! blank line before start of property list
         write(iounit, 4) 'version_number', "'" // trim(version_number) // "'"
         write(iounit, 1) 'M/Msun', s% star_mass
         write(iounit, 2) 'model_number', s% model_number
         write(iounit, 1) 'star_age', s% star_age
         write(iounit, 1) 'initial_z', s% initial_z
         write(iounit, 2) 'n_shells', nz
         write(iounit, 4) 'net_name', "'" // trim(s% net_name) // "'"
         write(iounit, 2) 'species', species
         if (s% M_center /= 0) then
            write(iounit, 11) 'xmstar', s% xmstar, &
               '! above core (g).  core mass: Msun, grams:', &
               s% M_center/Msun, s% M_center
         end if
         if (s% R_center /= 0) then
            write(iounit, 11) 'R_center', s% R_center, &
               '! radius of core (cm).  R/Rsun, avg core density (g/cm^3):', &
                  s% R_center/Rsun, s% M_center/(four_thirds_pi*pow3(s% R_center))
         end if
         if (s% v_center /= 0) then
            write(iounit, 11) 'v_center', s% v_center, &
               '! velocity of outer edge of core (cm/s)'
         end if
         if (s% L_center /= 0) then
            write(iounit, 11) 'L_center', s% L_center, &
               '! luminosity of core (erg/s). L/Lsun, avg core eps (erg/g/s):', &
                  s% L_center/Lsun, s% L_center/max(1d0,s% M_center)
         end if
         if (s% tau_factor /= 1) then
            write(iounit, 1) 'tau_factor', s% tau_factor
         end if
         if (s% Tsurf_factor /= 1) then
            write(iounit, 1) 'Tsurf_factor', s% Tsurf_factor
         end if
         if (s% opacity_factor /= 1) then
            write(iounit, 1) 'opacity_factor', s% opacity_factor
         end if
         write(iounit, 1) 'Teff', s% Teff
         write(iounit, 1) 'power_nuc_burn', s% power_nuc_burn
         write(iounit, 1) 'power_h_burn', s% power_h_burn
         write(iounit, 1) 'power_he_burn', s% power_he_burn
         write(iounit, 1) 'power_z_burn', s% power_z_burn
         write(iounit, 1) 'power_photo', s% power_photo
         write(iounit, 1) 'total_energy', s% total_energy
         write(iounit, 1) 'cumulative_energy_error', s% cumulative_energy_error
         if (abs(s% total_energy) > 1d0) &
            write(iounit, 11) 'cumulative_error/total_energy', &
               s% cumulative_energy_error/s% total_energy, &
               'log_rel_run_E_err', &
               safe_log10(abs(s% cumulative_energy_error/s% total_energy))
         write(iounit, 2) 'num_retries', s% num_retries
         write(iounit, '(a)') ! blank line for end of property list      

         call header
         do k=1, nz
            if (ierr /= 0) exit
            write(iounit, fmt='(i5, 1x)', advance='no') k
            call write1(s% lnd(k),ierr); if (ierr /= 0) exit
            call write1(s% lnT(k),ierr); if (ierr /= 0) exit
            call write1(s% lnR(k),ierr); if (ierr /= 0) exit            
            if (RSP_flag) then
               call write1(s% RSP_Et(k),ierr); if (ierr /= 0) exit
               call write1(s% erad(k),ierr); if (ierr /= 0) exit
               call write1(s% Fr(k),ierr); if (ierr /= 0) exit
            else if (RSP2_flag) then
               call write1(s% w(k),ierr); if (ierr /= 0) exit
               call write1(s% Hp_face(k),ierr)
               if (ierr /= 0) exit
            end if            
            call write1(s% L(k),ierr); if (ierr /= 0) exit
            call write1(s% dq(k),ierr); if (ierr /= 0) exit
            if (v_flag) then
               call write1(s% v(k),ierr); if (ierr /= 0) exit
            end if
            if (rotation_flag) then
               call write1(s% omega(k),ierr); if (ierr /= 0) exit
               call write1(s% j_rot(k),ierr); if (ierr /= 0) exit
            end if
            if (u_flag) then
               call write1(s% u(k),ierr); if (ierr /= 0) exit
            end if
            if (RTI_flag) then
               call write1(s% alpha_RTI(k),ierr); if (ierr /= 0) exit
            end if
            if (write_mlt_vc) then
               call write1(s% mlt_vc(k),ierr); if (ierr /= 0) exit
            end if
            do i=1, species
               call write1(s% xa(i,k),ierr); if (ierr /= 0) exit
            end do
            write(iounit, *)
         end do

         if (prev_flag) then
            write(iounit, '(a)')
            write(iounit, '(a)') '        previous model'
            ! write prev properties
            write(iounit, '(a)')
            write(iounit, 2) 'previous n_shells', s% nz_old
            write(iounit, 1, advance='no') 'previous mass (grams)'
            call write1_eol(s% mstar_old,ierr)
            write(iounit, 1, advance='no') 'timestep (seconds)'
            call write1_eol(s% dt,ierr)
            write(iounit, 1, advance='no') 'dt_next (seconds)'
            call write1_eol(s% dt_next_unclipped,ierr)
            write(iounit, '(a)')
         end if
         close(iounit)

         contains

         subroutine write1_eol(val,ierr)
            real(dp), intent(in) :: val
            integer, intent(out) :: ierr
            call write1(val,ierr)
            write(iounit,*)
         end subroutine write1_eol

         subroutine write1(val,ierr)
            real(dp), intent(in) :: val
            integer, intent(out) :: ierr
            integer, parameter :: str_len = 26
            character (len=str_len) :: string
            ierr = 0
            call double_to_str(val,string)
            write(iounit, fmt='(a26,1x)', advance='no') string
         end subroutine write1

         subroutine header
            write(iounit, fmt='(10x, a9, 1x, 99(a26, 1x))', advance='no') 'lnd', 'lnT', 'lnR'
            if (RSP_flag) then
               write(iounit, fmt='(a26, 1x)', advance='no') 'et_RSP'
               write(iounit, fmt='(a26, 1x)', advance='no') 'erad_RSP'
               write(iounit, fmt='(a26, 1x)', advance='no') 'Fr_RSP'
            else if (RSP2_flag) then
               write(iounit, fmt='(a26, 1x)', advance='no') 'w'
               write(iounit, fmt='(a26, 1x)', advance='no') 'Hp'
            end if
            write(iounit, fmt='(a26, 1x)', advance='no') 'L'
            write(iounit, fmt='(a26, 1x)', advance='no') 'dq'
            if (v_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'v'
            if (rotation_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'omega'
            if (rotation_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'j_rot'
            if (u_flag) write(iounit, fmt='(a26, 1x)', advance='no') 'u'
            if (RTI_flag) &
               write(iounit, fmt='(a26, 1x)', advance='no') 'alpha_RTI'
            if (write_mlt_vc) write(iounit, fmt='(a26, 1x)', advance='no') 'mlt_vc'
            do i=1, species
               write(iounit, fmt='(a26, 1x)', advance='no') chem_isos% name(chem_id(i))
            end do
            write(iounit, *)
         end subroutine header

      end subroutine do_write_model


      end module write_model
