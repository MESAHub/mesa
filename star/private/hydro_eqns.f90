! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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


      module hydro_eqns

      use star_private_def
      use const_def
      use star_utils, only: em1, e00, ep1
      use utils_lib, only: mesa_error, is_bad
      use auto_diff
      use auto_diff_support

      implicit none

      private
      public :: eval_equ

      logical, parameter :: dbg = .false.
      integer, parameter :: dbg_cell = -1


      contains


      subroutine eval_equ(s, nvar, skip_partials, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         call eval_equ_for_solver(s, nvar, 1, s% nz, skip_partials, ierr)
      end subroutine eval_equ


      subroutine eval_equ_for_solver(s, nvar, nzlo, nzhi, skip_partials, ierr)
         use chem_def
         use utils_lib, only: set_nan
         use mesh_functions
         use hydro_riemann, only: do_uface_and_Pface, do1_Riemann_momentum_eqn
         use hydro_momentum, only: do1_momentum_eqn, do1_radius_eqn
         use hydro_chem_eqns, only: do_chem_eqns, do1_chem_eqns
         use hydro_energy, only: do1_energy_eqn
         use hydro_tdc, only: do1_turbulent_energy_eqn, do1_tdc_L_eqn, do1_tdc_Hp_eqn
         use hydro_alpha_rti_eqns, only: do1_dalpha_RTI_dt_eqn
         use eps_grav, only: zero_eps_grav_and_partials
         use profile, only: do_save_profiles
         use star_utils, only: show_matrix, &
            no_extra_profile_columns, no_data_for_extra_profile_columns

         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nzlo, nzhi
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         integer :: &
            i_dv_dt, i_du_dt, i_du_dk, i_equL, i_dlnd_dt, i_dlnE_dt, i_dlnR_dt, &
            i_dalpha_RTI_dt, i_dln_cvpv0_dt, i_equ_w_div_wc, i_dj_rot_dt, i_detrb_dt, &
            i, k, j, nvar_hydro, nz, op_err
         integer :: &
            i_lnd, i_lnR, i_lnT, i_lum, i_v, i_u, i_du, i_ln_cvpv0, i_w_div_wc, i_j_rot, &
            i_alpha_RTI, i_chem1, i_xh1, i_xhe4, kmax_equ(nvar), species
         real(dp) :: max_equ(nvar), L_phot_old
         real(dp), dimension(:), pointer :: &
            L, lnR, lnP, lnT, energy
         logical :: v_flag, u_flag, conv_vel_flag, cv_flag, w_div_wc_flag, j_rot_flag, dump_for_debug, &
            do_chem, do_mix, do_dlnd_dt, do_dv_dt, do_du_dt, do_dlnR_dt, &
            do_alpha_RTI, do_conv_vel, do_w_div_wc, do_j_rot, do_dlnE_dt, do_equL, do_detrb_dt

         include 'formats'

         ierr = 0

         if (s% u_flag .and. s% use_mass_corrections) &
            stop 'use_mass_corrections dP not supported with u_flag true'

         if (s% u_flag) then
            call do_uface_and_Pface(s,ierr)
            if (ierr /= 0) then
               if (len_trim(s% retry_message) == 0) s% retry_message = 'do_uface_and_Pface failed'
               if (s% report_ierr) write(*,*) 'ierr from do_uface_and_Pface'
               return
            end if
         end if

         dump_for_debug = .false.
         !dump_for_debug = .true.

         do_mix = s% do_mix
         do_chem = (do_mix .or. s% do_burn)

         call unpack
         
         ! set flags indicating the variables currently in use
         do_dlnd_dt = (i_dlnd_dt > 0 .and. i_dlnd_dt <= nvar)
         do_dv_dt = (i_dv_dt > 0 .and. i_dv_dt <= nvar)
         do_du_dt = (i_du_dt > 0 .and. i_du_dt <= nvar)
         do_dlnR_dt = (i_dlnR_dt > 0 .and. i_dlnR_dt <= nvar)
         do_alpha_RTI = (i_alpha_RTI > 0 .and. i_alpha_RTI <= nvar)
         do_conv_vel = (i_ln_cvpv0 > 0 .and. i_ln_cvpv0 <= nvar)
         do_w_div_wc = (i_w_div_wc > 0 .and. i_w_div_wc <= nvar)
         do_j_rot = (i_j_rot > 0 .and. i_j_rot <= nvar)
         do_dlnE_dt = (i_dlnE_dt > 0 .and. i_dlnE_dt <= nvar)
         do_equL = (i_equL > 0 .and. i_equL <= nvar)
         do_detrb_dt = (i_detrb_dt > 0 .and. i_detrb_dt <= nvar)

         if (s% fill_arrays_with_NaNs) call set_nan(s% equ1)

      ! solving structure equations

!$OMP PARALLEL DO PRIVATE(op_err,k) SCHEDULE(dynamic,2)
         do k = nzlo, nzhi
            op_err = 0
            ! hack for composition partials
            if (s% use_d_eos_dxa .and. .not. skip_partials) then
               call fix_d_eos_dxa_partials(s, k, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr from fix_d_eos_dxa_partials', k
                  ierr = op_err
               end if
            end if
         end do
!$OMP END PARALLEL DO

         if (s% use_other_energy_implicit) then
            call s% other_energy_implicit(s% id, ierr)
            if (ierr /= 0) then
               if (len_trim(s% retry_message) == 0) s% retry_message = 'other_energy_implicit failed'
               if (s% report_ierr) &
                  write(*,*) 'ierr from other_energy_implicit'
               return
            end if
         end if

         if (s% use_other_momentum_implicit) then
            call s% other_momentum_implicit(s% id, ierr)
            if (ierr /= 0) then
               if (len_trim(s% retry_message) == 0) s% retry_message = 'other_momentum_implicit failed'
               if (s% report_ierr) &
                  write(*,*) 'ierr from other_momentum_implicit'
               return
            end if
         end if
            
!$OMP PARALLEL DO PRIVATE(op_err,k) SCHEDULE(dynamic,2)
         do k = nzlo, nzhi
            if (.not. skip_partials) then
               s% dblk(:,:,k) = 0
               s% ublk(:,:,k) = 0
               s% lblk(:,:,k) = 0
            end if

            op_err = 0
            if (do_dlnd_dt) then
               call do1_density_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_density_eqn'
                  if (s% report_ierr) write(*,2) 'ierr in do1_density_eqn', k
                  ierr = op_err
               end if
            end if
            if (k > 1) then ! k=1 is surf P BC
               if (do_du_dt) then
                  call do1_Riemann_momentum_eqn(s, k, skip_partials, nvar, op_err)
                  if (op_err /= 0) then
                     if (s% report_ierr) write(*,2) 'ierr in do1_Riemann_momentum_eqn', k
                     if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_Riemann_momentum_eqn'
                     ierr = op_err
                  end if
               end if
               if (do_dv_dt) then
                  call do1_momentum_eqn(s, k, skip_partials, nvar, op_err)
                  if (op_err /= 0) then
                     if (s% report_ierr) write(*,2) 'ierr in do1_momentum_eqn', k
                     if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_momentum_eqn'
                     ierr = op_err
                  end if
               end if
            end if
            if (do_dlnR_dt) then
               call do1_radius_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_radius_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_radius_eqn'
                  ierr = op_err
               end if
            end if
            if (do_alpha_RTI) then
               call do1_dalpha_RTI_dt_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_dalpha_RTI_dt_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_dalpha_RTI_dt_eqn'
                  ierr = op_err
               end if
            end if
            if (do_conv_vel) then
               call do1_dln_cvpv0_dt_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_dln_cvpv0_dt_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_dln_cvpv0_dt_eqn'
                  ierr = op_err
               end if
            end if
            if (do_w_div_wc) then
               call do1_w_div_wc_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_w_div_wc_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_w_div_wc_eqn'
                  ierr = op_err
               end if
            end if
            if (do_j_rot) then
               call do1_dj_rot_dt_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_dj_rot_dt_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_dj_rot_dt_eqn'
                  ierr = op_err
               end if
            end if

            if (do_dlnE_dt) then
               call zero_eps_grav_and_partials(s, k)
               call do1_energy_eqn(s, k, skip_partials, do_chem, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_energy_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_energy_eqn'
                  ierr = op_err
               end if
            end if
            if (do_detrb_dt) then
               call do1_turbulent_energy_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_turbulent_energy_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_turbulent_energy_eqn'
                  ierr = op_err
               end if
               call do1_tdc_Hp_eqn(s, k, skip_partials, nvar, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_tdc_Hp_eqn', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_tdc_Hp_eqn'
                  ierr = op_err
               end if
            end if
            if (do_equL) then
               if (s% TDC_flag) then
                  call do1_tdc_L_eqn(s, k, skip_partials, nvar, op_err)
                  if (op_err /= 0) then
                     if (s% report_ierr) write(*,2) 'ierr in do1_tdc_L_eqn', k
                     if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_tdc_L_eqn'
                     ierr = op_err
                  end if
               else if (k > 1) then ! k==1 is done by T_surf BC
                  call do1_dlnT_dm_eqn(s, k, skip_partials, nvar, op_err)
                  if (op_err /= 0) then
                     if (s% report_ierr) write(*,2) 'ierr in do1_dlnT_dm_eqn', k
                     if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_dlnT_dm_eqn'
                     ierr = op_err
                  end if
               end if
            end if

            if (do_chem) then
               call do1_chem_eqns(s, k, nvar, skip_partials, op_err)
               if (op_err /= 0) then
                  if (s% report_ierr) write(*,2) 'ierr in do1_chem_eqns', k
                  if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_chem_eqns'
                  ierr = op_err
               end if
            end if
         end do
!$OMP END PARALLEL DO
         
         if (ierr == 0 .and. nzlo == 1) then
            if (dbg) write(*,*) 'call PT_eqns_surf'
            call PT_eqns_surf( &
               s, skip_partials, nvar, &
               do_du_dt, do_dv_dt, do_equL, ierr)
            if (dbg) write(*,*) 'done PT_eqns_surf'
            if (ierr /= 0) then
               if (s% report_ierr) write(*,2) 'ierr in PT_eqns_surf', ierr
               if (len_trim(s% retry_message) == 0) s% retry_message = 'error in PT_eqns_surf'
            end if
         end if

         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'ierr in eval_equ_for_solver'
            return
         end if
         
         !if (.not. skip_partials) stop 'after eval_equ_for_solver'
         
         if (.false. .and. .not. skip_partials .and. s% model_number == 2) then !  .and. .not. s% doing_relax) then
            if (.false.) then
               i = s% i_dv_dt
               k = 1
               do j=1,5
                  write(*,5) 'dblk(i,j,k) ' // &
                     trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), &
                     s% solver_iter, i, j, k, s% dblk(i,j,k)
               end do
            end if
            !write(*,*) 'call show_matrix'
            !call show_matrix(s, s% dblk(1:s% nvar_hydro,1:s% nvar_hydro,1), s% nvar_hydro)
            call dump_equ ! debugging
            stop 'after dump_equ'
         end if

         
         contains

         subroutine dump_equ
            integer :: k, j, k0, k1
            include 'formats'
            do k=1,s% nz
               do j=1,nvar
                  write(*,3) 'equ ' // trim(s% nameofequ(j)), &
                     k, s% solver_iter, s% equ(j, k)
               end do
               write(*,3) 'eps_mdot', k, s% solver_iter, s% eps_mdot(k)
               write(*,3) 'energy_start', k, s% solver_iter, s% energy_start(k)
               write(*,3) 'energy', k, s% solver_iter, s% energy(k)
               write(*,3) 'L', k, s% solver_iter, s% L(k)
               write(*,3) 'L', k+1, s% solver_iter, s% L(k+1)
               write(*,*)
               !if (k == 6) exit
            end do
         end subroutine dump_equ

         subroutine unpack

            include 'formats'

            nz = s% nz
            species = s% species
            nvar_hydro = s% nvar_hydro

            lnT => s% lnT
            lnR => s% lnR
            L => s% L
            lnP => s% lnP
            energy => s% energy

            i_dv_dt = s% i_dv_dt
            i_du_dt = s% i_du_dt
            i_equL = s% i_equL
            i_dlnd_dt = s% i_dlnd_dt
            i_dlnE_dt = s% i_dlnE_dt
            i_dlnR_dt = s% i_dlnR_dt
            i_dalpha_RTI_dt = s% i_dalpha_RTI_dt
            i_dln_cvpv0_dt = s% i_dln_cvpv0_dt
            i_equ_w_div_wc = s% i_equ_w_div_wc
            i_dj_rot_dt = s% i_dj_rot_dt
            i_detrb_dt = s% i_detrb_dt

            i_lnd = s% i_lnd
            i_lnT = s% i_lnT
            i_lnR = s% i_lnR
            i_lum = s% i_lum
            i_v = s% i_v
            i_u = s% i_u
            i_alpha_RTI = s% i_alpha_RTI
            i_ln_cvpv0 = s% i_ln_cvpv0
            i_w_div_wc = s% i_w_div_wc
            i_j_rot = s% i_j_rot

            i_chem1 = s% i_chem1
            i_xh1 = i_chem1-1+s% net_iso(ih1)
            i_xhe4 = i_chem1-1+s% net_iso(ihe4)

            L_phot_old = s% L_phot_old
            v_flag = s% v_flag
            u_flag = s% u_flag

         end subroutine unpack

         subroutine fix_d_eos_dxa_partials(s, k, ierr)

            ! revise composition partials
            ! subroutine can be removed when EOS fully provides composition partials

            use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_lnE, i_lnPgas
            use eos_support, only: get_eos
            use star_utils, only: lookup_nameofvar
            type (star_info), pointer :: s
            integer, intent(in) :: k
            integer, intent(out) :: ierr

            integer :: j

            logical, parameter :: debug = .false.

            ! these vars are for faking composition derivatives
            real(dp), dimension(num_eos_basic_results) :: &
               res, dres_dlnd, dres_dlnT, dres_dabar, dres_dzbar
            real(dp) :: dres_dxa(num_eos_d_dxa_results, s% species)
            real(dp) :: dxa
            real(dp) :: xa_start_1(s% species)
            real(dp) :: frac_without_dxa
            real(dp) :: lnE_with_xa_start, lnPgas_with_xa_start

            integer :: i_var, i_var_sink

            real(dp) :: dxa_threshold = 1d-4

            logical, parameter :: checking = .true.

            include 'formats'

            ierr = 0

            ! some EOSes have composition partials and some do not
            ! those currently without dx partials are PC & Skye
            frac_without_dxa = s% eos_frac_PC(k) + s% eos_frac_Skye(k)
            
            if (debug .and. k == s% solver_test_partials_k) then
              write(*,2) 's% eos_frac_PC(k)', k, s% eos_frac_PC(k)
              write(*,2) 's% eos_frac_Skye(k)', k, s% eos_frac_Skye(k)
              write(*,2) 'frac_without_dxa', k, frac_without_dxa
            end if

            if (k == s% solver_test_partials_k .and. s% solver_iter == s% solver_test_partials_iter_number) then
               i_var = lookup_nameofvar(s, s% solver_test_partials_var_name)
               if (i_var .gt. s% nvar_hydro) then
                  i_var_sink = lookup_nameofvar(s, s% solver_test_partials_sink_name)
               end if
            end if

            ! if we're on an EOS where there aren't composition partials,
            ! approximate derivatives with finite differences
            if (frac_without_dxa .gt. 0) then

               do j=1, s% species
                  dxa = s% xa(j,k) - s% xa_start(j,k)

                  if (debug .and. k == s% solver_test_partials_k .and. &
                        s% solver_iter == s% solver_test_partials_iter_number) &
                     write(*,2) 'dxa', j, dxa

                  if (abs(dxa) .ge. dxa_threshold) then

                     ! first, get eos with xa_start

                     call get_eos( &
                        s, k, s% xa_start(:,k), &
                        s% rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
                        res, dres_dlnd, dres_dlnT, dres_dxa, ierr)
                     if (ierr /= 0) then
                        if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start', k
                        return
                     end if

                     lnE_with_xa_start = res(i_lnE)
                     lnPgas_with_xa_start = res(i_lnPgas)

                     ! now, get eos with 1 iso perturbed

                     xa_start_1 = s% xa_start(:,k)
                     xa_start_1(j) = s% xa_start(j,k) + dxa

                     call get_eos( &
                        s, k, xa_start_1, &
                        s% Rho(k), s% lnd(k)/ln10, s% T(k), s% lnT(k)/ln10, &
                        res, dres_dlnd, dres_dlnT, dres_dxa, ierr)
                     if (is_bad(res(i_lnPgas))) ierr = -1
                     if (ierr /= 0) then
                        
                        ! punt silently for now
                        s% dlnE_dxa_for_partials(:,k) = 0d0
                        s% dlnP_dxa_for_partials(:,k) = 0d0
                        ierr = 0
                        return
                        
                        if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start_1', k
                        return
                     end if

                     ! fix up derivatives

                     if (debug .and. k == s% solver_test_partials_k) & ! .and. s% solver_iter == s% solver_test_partials_iter_number) &
                        write(*,2) 'res(i_lnE) - lnE_with_xa_start', j, res(i_lnE) - lnE_with_xa_start

                     s% dlnE_dxa_for_partials(j,k) = dres_dxa(i_lnE, j) + &
                        frac_without_dxa * (res(i_lnE) - lnE_with_xa_start) / dxa
                     if (checking) call check_dxa(j,k,s% dlnE_dxa_for_partials(j,k),'dlnE_dxa_for_partials')

                     s% dlnP_dxa_for_partials(j,k) = s% Pgas(k)*dres_dxa(i_lnPgas,j)/s% P(k) + &
                        frac_without_dxa * (s% Pgas(k)/s% P(k)) * (res(i_lnPgas) - lnPgas_with_xa_start) / dxa
                     if (.false. .and. s% model_number == 1100 .and. k == 151 .and. j == 1 .and. is_bad(s% dlnP_dxa_for_partials(j,k))) then
                        write(*,2) 's% Pgas(k)', k, s% Pgas(k)
                        write(*,2) 'dres_dxa(i_lnPgas,j)', k, dres_dxa(i_lnPgas,j)
                        write(*,2) 's% P(k)', k, s% P(k)
                        write(*,2) 'frac_without_dxa', k, frac_without_dxa
                        write(*,2) 'res(i_lnPgas)', k, res(i_lnPgas)
                        write(*,2) 'lnPgas_with_xa_start', k, lnPgas_with_xa_start
                        write(*,2) 'dxa', k, dxa
                        write(*,2) 'dxa_threshold', k, dxa_threshold
                        write(*,2) 's% xa_start(j,k)', j, s% xa_start(j,k)
                        write(*,2) 'xa_start_1(j)', j, xa_start_1(j)
                        write(*,2) 's% eos_frac_PC(k)', k, s% eos_frac_PC(k)
                        write(*,2) 's% eos_frac_Skye(k)', k, s% eos_frac_Skye(k)
                     end if
                     if (checking) call check_dxa(j,k,s% dlnP_dxa_for_partials(j,k),'dlnP_dxa_for_partials')

                  else

                     if (k == s% solver_test_partials_k .and. s% solver_iter == s% solver_test_partials_iter_number) then
                        if (i_var_sink > 0 .and. i_var > s% nvar_hydro) then
                           if (dxa < dxa_threshold) then
                              if (j .eq. i_var - s% nvar_hydro) then
                                 write(*,*) 'fix_d_eos_dxa_partials: skipping dxa derivative fix for ', trim (s% solver_test_partials_var_name), &
                                    ' (dxa < dxa_threshold): ', abs(dxa), ' < ', dxa_threshold
                              endif
                              if (j .eq. i_var_sink - s% nvar_hydro) then
                                 write(*,*) 'fix_d_eos_dxa_partials: skipping dxa derivative fix for ', trim (s% solver_test_partials_sink_name), &
                                    ' (dxa < dxa_threshold): ', abs(dxa), ' < ', dxa_threshold
                              end if
                           end if
                        end if
                     end if

                  end if
               end do

            end if

         end subroutine fix_d_eos_dxa_partials

         subroutine check_dxa(j, k, dxa, str)
            integer, intent(in) :: j, k
            real(dp), intent(in) :: dxa
            character (len=*), intent(in) :: str
            include 'formats'
            if (is_bad(dxa)) then
!$omp critical (eval_equ_for_solver_crit1)
               ierr = -1
               if (s% report_ierr) then
                  write(*,3) 'eval_equ_for_solver: bad ' // trim(str), j, k, dxa
               end if
               if (s% stop_for_bad_nums) stop 'eval_equ_for_solver'
!$omp end critical (eval_equ_for_solver_crit1)
               return
            end if
         end subroutine check_dxa


      end subroutine eval_equ_for_solver


      subroutine do1_density_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: &
            rho, dr3, r_p1, rp13, lnR_expected, lnR_actual, resid_ad
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         rho = wrap_d_00(s,k)
         dr3 = (s% dm(k)/rho)/(4d0*pi/3d0)
         r_p1 = wrap_r_p1(s,k)
         rp13 = pow3(r_p1)

         lnR_expected = log(rp13 + dr3)/3d0
         lnR_actual = wrap_lnR_00(s,k)
         
         resid_ad = lnR_actual - lnR_expected
         s% equ(s% i_dlnd_dt, k) = resid_ad%val

         if (test_partials) then
            s% solver_test_partials_val = 0
         end if

         if (skip_partials) return
         call save_eqn_residual_info( &
            s, k, nvar, s% i_dlnd_dt, resid_ad, 'do1_density_eqn', ierr)           

         if (test_partials) then   
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_density_eqn', s% solver_test_partials_var
         end if

      end subroutine do1_density_eqn


      subroutine do1_dln_cvpv0_dt_eqn(s, k, skip_partials, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         integer :: i_dln_cvpv0_dt, i_ln_cvpv0
         real(dp) :: residual, dln_cvpv0_m1, dln_cvpv0_00, dln_cvpv0_p1
         real(dp) :: d_actual, d_d_actual_dvc, d_expected, d_d_expected_dvc, &
            dt_div_tau, scale, tau, dt, cs, N2_start, grav_start, &
            dt_div_lambda, lambda, D, siglimit, sigmid_00, sigmid_m1, rc, f, &
            mix, dmix_dcvm1, dmix_dcv00, dmix_dcvp1, a, da_dvc, inv, d_inv_dvc, &
            c, dc_dvc, alfa, d_ln_cvpv0_dq
         logical :: test_partials
         include 'formats'
         
         ierr = 0         
         !test_partials = (k == s% solver_test_partials_k-1)
         test_partials = .false.
         
         stop 'PABLO needs to rewrite using auto_diff'
                  
      end subroutine do1_dln_cvpv0_dt_eqn
      

      subroutine do1_w_div_wc_eqn(s, k, skip_partials, nvar, ierr)
         use hydro_rotation
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         integer :: i_equ_w_div_wc, i_w_div_wc
         real(dp) :: scale, wwc, dimless_rphi, dimless_rphi_given_wwc, w1, w2
         real(dp) :: jrot_ratio, sigmoid_jrot_ratio, d_sigmoid_jrot_ratio_dx, &
            jr_lim1, jr_lim2, A, C, dA_dw, dC_dw
         logical :: test_partials
         include 'formats'
         
         ierr = 0

         
         stop 'PABLO needs to rewrite using auto_diff'
         
      end subroutine do1_w_div_wc_eqn
      
      
      subroutine do1_dj_rot_dt_eqn(s, k, skip_partials, nvar, ierr)
         use hydro_rotation
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         integer :: i_dj_rot_dt, i_j_rot
         real(dp) :: scale
         real(dp) :: F00, dF00_dw00, dF00_dwp1, dF00_dj00, dF00_djp1, dF00_dlnr00, dF00_dlnrp1, dF00_dlnd00
         real(dp) :: Fm1, dFm1_dwm1, dFm1_dw00, dFm1_djm1, dFm1_dj00, dFm1_dlnrm1, dFm1_dlnr00, dFm1_dlndm1
         logical :: test_partials
         include 'formats'
         
         ierr = 0
         
         stop 'PABLO needs to rewrite using auto_diff'
         
      end subroutine do1_dj_rot_dt_eqn


      ! just relate L_rad to T gradient.
      ! d_P_rad/dm = -<opacity_face>*L_rad/(clight*area^2) -- see, e.g., K&W (5.12)
      ! P_rad = (1/3)*crad*T^4
      ! d_P_rad/dm = (crad/3)*(T(k-1)^4 - T(k)^4)/dm_bar
      ! L_rad = L - L_non_rad, L_non_rad = L_start - L_rad_start
      ! L_rad_start = (-d_P_rad/dm_bar*clight*area^2/<opacity_face>)_start
      subroutine do1_alt_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         real(dp) :: alfa, beta, scale, dm_bar
         type(auto_diff_real_star_order1) :: L_ad, r_00, area, area2, Lrad_ad, &
            kap_00, kap_m1, kap_face, d_P_rad_expected_ad, T_m1, T4_m1, T_00, T4_00, &
            P_rad_m1, P_rad_00, d_P_rad_actual_ad, resid
         
         integer :: i_equL, i
         logical :: dbg
         logical :: test_partials

         include 'formats'
         ierr = 0
         i_equL = s% i_equL
         if (i_equL == 0) return

         if (.not. s% use_dPrad_dm_form_of_T_gradient_eqn) then
            ierr = -1
            return
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         dbg = .false.

         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         beta = 1d0 - alfa

         scale = s% energy_start(k)*s% rho_start(k)
         dm_bar = s% dm_bar(k)
         L_ad = wrap_L_00(s,k)
         r_00 = wrap_r_00(s,k)
         area = pi4*pow2(r_00); area2 = pow2(area)

         if ((check_flag_and_val(s% conv_vel_flag, s% conv_vel, k)) .or. &
               (.not. s% conv_vel_flag .and. s% lnT(k)/ln10 <= s% max_logT_for_mlt &
               .and. s% mixing_type(k) == convective_mixing .and. s% gradr(k) > 0d0 &
               .and. abs(s% gradr(k) - s% gradT(k)) > abs(s% gradr(k))*1d-5)) then
            Lrad_ad = L_ad*s% gradT_ad(k)/s% gradr_ad(k) ! C&G 14.109
         else
            Lrad_ad = L_ad
         end if

         kap_00 = wrap_kap_00(s,k)
         kap_m1 = wrap_kap_m1(s,k)
         kap_face = alfa*kap_00 + beta*kap_m1
         if (kap_face%val < s% min_kap_for_dPrad_dm_eqn) &
            kap_face = s% min_kap_for_dPrad_dm_eqn
                  
         ! calculate expected d_P_rad from current L_rad
         d_P_rad_expected_ad = -dm_bar*kap_face*Lrad_ad/(clight*area2)
         
         ! calculate actual d_P_rad in current model
         T_m1 = wrap_T_m1(s,k); T4_m1 = pow4(T_m1)
         T_00 = wrap_T_00(s,k); T4_00 = pow4(T_00)

         !d_P_rad_expected = d_P_rad_expected*s% gradr_factor(k) !TODO(Pablo): check this

         P_rad_m1 = (crad/3d0)*T4_m1
         P_rad_00 = (crad/3d0)*T4_00
         d_P_rad_actual_ad = P_rad_m1 - P_rad_00
         
         ! residual
         resid = (d_P_rad_expected_ad - d_P_rad_actual_ad)/scale 
         s% equ(i_equL, k) = resid%val
         
         if (is_bad(resid%val)) then
!$OMP critical (star_alt_dlntdm_bad_num)
            write(*,2) 'resid%val', k, resid%val
            if (s% stop_for_bad_nums) stop 'do1_alt_dlnT_dm_eqn'
!$OMP end critical (star_alt_dlntdm_bad_num)
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% gradT(k)
         end if

         if (skip_partials) return
         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_alt_dlnT_dm_eqn', ierr)

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'do1_alt_dlnT_dm_eqn', s% solver_test_partials_var
         end if

         contains 
         
         logical function check_flag_and_val(flag, array, index)
            logical,intent(in) :: flag
            real(dp), dimension(:),intent(in) :: array
            integer, intent(in) :: index

            check_flag_and_val = .false.
            if(flag) then
               if(array(index)>0d0) check_flag_and_val = .true.
            end if
         end function check_flag_and_val

      end subroutine do1_alt_dlnT_dm_eqn


      subroutine do1_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)
         use eos_def
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: resid, &
            dlnPdm, Ppoint, gradT, dlnTdm, T00, Tm1, dT, Tpoint, lnTdiff
         real(dp) :: delm, alfa
         integer :: i_equL
         logical :: test_partials

         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_equL = s% i_equL
         if (i_equL == 0) return

         if (s% use_dPrad_dm_form_of_T_gradient_eqn .or. s% conv_vel_flag) then
            call do1_alt_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)            
            return
         end if
         
         ! dT/dm = dP/dm * T/P * grad_T, grad_T = dlnT/dlnP from MLT.
         ! but use hydrostatic value for dP/dm in this.
         ! this is because of limitations of MLT for calculating grad_T.
         ! (MLT assumes hydrostatic equilibrium)
         ! see comment in K&W chpt 9.1.
         
         call eval_dlnPdm_qhse(s, k, dlnPdm, Ppoint, ierr)
         if (ierr /= 0) return

         gradT = s% gradT_ad(k)
         dlnTdm = dlnPdm*gradT

         Tm1 = wrap_T_m1(s,k)
         T00 = wrap_T_00(s,k)
         dT = Tm1 - T00
         alfa = s% dm(k-1)/(s% dm(k-1) + s% dm(k))
         Tpoint = alfa*T00 + (1d0 - alfa)*Tm1
         lnTdiff = dT/Tpoint ! use this in place of lnT(k-1)-lnT(k)
         delm = (s% dm(k) + s% dm(k-1))/2
         
         resid = delm*dlnTdm - lnTdiff
         s% equ(i_equL, k) = resid%val

         if (is_bad(s% equ(i_equL, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            if (s% stop_for_bad_nums) stop 'hydro eqns'
            return
            write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            write(*,2) 'lnTdiff', k, lnTdiff
            write(*,2) 'delm', k, delm
            write(*,2) 'dlnPdm', k, dlnPdm
            write(*,2) 'gradT', k, gradT
            stop 'i_equL'
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% equ(i_equL,k)
         end if

         if (skip_partials) return
         call save_eqn_residual_info( &
            s, k, nvar, i_equL, resid, 'do1_dlnT_dm_eqn', ierr)

      end subroutine do1_dlnT_dm_eqn


      subroutine PT_eqns_surf( &
            s, skip_partials, nvar, do_du_dt, do_dv_dt, do_equL, ierr)

         use star_utils, only: save_eqn_residual_info
         use eos_lib, only: Radiation_Pressure

         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         logical, intent(in) :: do_du_dt, do_dv_dt, do_equL
         integer, intent(out) :: ierr

         type(auto_diff_real_star_order1) :: &
            P_bc_ad, T_bc_ad, lnP_bc_ad, lnT_bc_ad, resid_ad
         integer :: i_P_eqn, i_T_eqn
         logical :: need_P_surf, need_T_surf, test_partials

         include 'formats'

         !test_partials = (s% solver_iter == s% solver_test_partials_iter_number)
         test_partials = .false.         
         ierr = 0
         if (s% u_flag) then
            i_P_eqn = s% i_du_dt
         else
            i_P_eqn = s% i_dv_dt
         end if
         i_T_eqn = s% i_equL

         need_P_surf = .false.
         if (s% use_compression_outer_BC) then
            call set_compression_BC(ierr)
         else if (s% use_zero_Pgas_outer_BC) then
            P_bc_ad = Radiation_Pressure(s% T_start(1))
            call set_momentum_BC(ierr)
         else if (s% use_fixed_Psurf_outer_BC) then
            P_bc_ad = s% fixed_Psurf
            call set_momentum_BC(ierr)
         else if (s% use_fixed_vsurf_outer_BC) then
            call set_fixed_vsurf_outer_BC(ierr)
         else 
            need_P_surf = .true.
         end if
         if (ierr /= 0) return

         need_T_surf = .false.
         if (s% TDC_flag .or. .not. do_equL) then 
            ! no Tsurf BC
         else if (s% use_zero_dLdm_outer_BC) then
            call set_zero_dL_dm_BC(ierr)
         else
            need_T_surf = .true.
         end if
         if (ierr /= 0) return
         
         if (need_P_surf .or. need_T_surf) then
            call get_PT_bc_ad(ierr) 
            if (ierr /= 0) return
         else
            s% P_surf = s% P(1)
            s% T_surf = s% T(1)
            return
         end if
         
         ! use P_bc_ad and T_bc_ad
         
         if (need_P_surf) then
            if (s% use_momentum_outer_BC) then
               call set_momentum_BC(ierr)
            else
               call set_Psurf_BC(ierr)
            end if
            if (ierr /= 0) return
         end if

         if (need_T_surf) then         
            call set_Tsurf_BC(ierr)
            if (ierr /= 0) return
         end if

         if (test_partials) then
            s% solver_test_partials_val = 0
         end if

         if (skip_partials) return

         if (test_partials) then
            s% solver_test_partials_var = 0
            s% solver_test_partials_dval_dx = 0
            write(*,*) 'PT_eqns_surf', s% solver_test_partials_var
         end if

         contains
         
         subroutine get_PT_bc_ad(ierr) ! set P_bc_ad and T_bc_ad
            use hydro_vars, only: set_Teff_info_for_eqns
            use chem_def
            use atm_def
            integer, intent(out) :: ierr
            real(dp) :: r, L, Teff, &
               lnT_surf, dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
               lnP_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap
            real(dp) :: &
               dlnT_bc_dlnd, dlnT_bc_dlnT, dlnT_bc_dlnR, &
               dlnT_bc_dL, dlnP_bc_dlnd, dlnP_bc_dlnT, dlnP_bc_dL, dlnP_bc_dlnR, &
               dlnkap_dlnd, dlnkap_dlnT, dPinv_dlnd, dPinv_dlnT, dP0, dT0, &
               P_surf, T_surf, dlnP_bc_dlnPsurf, P_rad, &
               dlnT_bc_dlnTsurf, P_bc, T_bc, lnT_bc, lnP_bc, &
               dP0_dlnR, dT0_dlnR, dT0_dlnT, dT0_dlnd, dT0_dL, dlnP_bc_dP0, dlnT_bc_dT0, &
               dlnT_bc_dlnE_const_Rho, dlnT_dlnE_const_Rho, dlnP_dlnE_c_Rho, &
               dlnP_bc_dlnE_c_Rho, dlnT_bc_dlnd_c_E, dlnP_bc_dlnd_c_E, &
               d_gradT_dlnR, d_gradT_dlnT00, d_gradT_dlnd00, d_gradT_dL, &
               dlnR00, dlnT00, dlnd00
            include 'formats'
            ierr = 0
         
            call set_Teff_info_for_eqns(s, skip_partials, r, L, Teff, &
               lnT_surf, dlnTsurf_dL, dlnTsurf_dlnR, dlnTsurf_dlnM, dlnTsurf_dlnkap, &
               lnP_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnM, dlnPsurf_dlnkap, &
               ierr)
            if (ierr /= 0) then
               if (s% report_ierr) then
                  write(*,*) 'PT_eqns_surf: ierr from set_Teff_info_for_eqns'
               end if
               return
            end if

            ! P_surf and T_surf are at outer boundary of cell 1
            P_surf = exp(lnP_surf)
            T_surf = exp(lnT_surf)
            
            s% P_surf = P_surf
            s% T_surf = T_surf

            if (s% use_atm_PT_at_center_of_surface_cell .or. r <= 0d0) then
               dP0 = 0
               dT0 = 0
            else
               if (s% use_momentum_outer_BC) then
                  dP0 = 0
               else
                  dP0 = s% cgrav(1)*s% m_grav(1)*s% dm(1)/(8*pi*r*r*r*r)
               end if
               dT0 = dP0*s% gradT(1)*s% T(1)/s% P(1)
            end if
         
            P_bc = P_surf + dP0
            T_bc = T_surf + dT0

            lnP_bc = log(P_bc)
            lnT_bc = log(T_bc)
         
            if (is_bad(P_bc)) then
               write(*,1) 'lnP_bc', lnP_bc
               write(*,1) 'P_bc', P_bc
               write(*,1) 'P_surf', P_surf
               write(*,1) 'dP0', dP0
               write(*,1) 'lnP_surf', lnP_surf
               stop 'bc'
            end if

            if (s% use_atm_PT_at_center_of_surface_cell .or. r <= 0d0) then

               dP0_dlnR = 0
               dT0_dlnR = 0
               dT0_dlnT = 0
               dT0_dlnd = 0
               dT0_dL = 0

            else
         
               d_gradT_dlnR = s% gradT_ad(1)%d1Array(i_lnR_00)
               d_gradT_dlnT00 = s% gradT_ad(1)%d1Array(i_lnT_00)
               d_gradT_dlnd00 = s% gradT_ad(1)%d1Array(i_lnd_00)
               d_gradT_dL = s% gradT_ad(1)%d1Array(i_L_00)

               dP0_dlnR = -4*dP0
               dT0_dlnR = -4*dT0 + dP0*d_gradT_dlnR*s% T(1)/s% P(1)

               dPinv_dlnT = -s% chiT_for_partials(1)/s% P(1)
               dT0_dlnT = &
                    dT0 + &
                    dP0*d_gradT_dlnT00*s% T(1)/s% P(1) + &
                    dP0*s% gradT(1)*s% T(1)*dPinv_dlnT
               dPinv_dlnd = -s% chiRho_for_partials(1)/s% P(1)
               dT0_dlnd = &
                    dP0*d_gradT_dlnd00*s% T(1)/s% P(1) + &
                    dP0*s% gradT(1)*s% T(1)*dPinv_dlnd

               dT0_dL = dP0*d_gradT_dL*s% T(1)/s% P(1)

            endif

            dlnP_bc_dP0 = 1/P_bc
            dlnT_bc_dT0 = 1/T_bc

            dlnP_bc_dlnPsurf = P_surf/P_bc
            dlnT_bc_dlnTsurf = T_surf/T_bc

            dlnkap_dlnd = s% d_opacity_dlnd(1)/s% opacity(1)
            dlnkap_dlnT = s% d_opacity_dlnT(1)/s% opacity(1)

            dlnP_bc_dlnd = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnd
            dlnP_bc_dlnT = dlnP_bc_dlnPsurf*dlnPsurf_dlnkap*dlnkap_dlnT
            dlnP_bc_dL = dlnP_bc_dlnPsurf*dlnPsurf_dL
            dlnP_bc_dlnR = dlnP_bc_dlnPsurf*dlnPsurf_dlnR + dlnP_bc_dP0*dP0_dlnR

            dlnR00 = P_bc*dlnP_bc_dlnR
            if (.not. s% solver_use_lnR) dlnR00 = dlnR00/s% r(1)
            dlnT00 = P_bc*dlnP_bc_dlnT
            if (.not. s% solver_use_lnT) dlnT00 = dlnT00/s% T(1)
            dlnd00 = P_bc*dlnP_bc_dlnd
            if (.not. s% solver_use_lnd) dlnd00 = dlnd00/s% rho(1)
            call wrap(P_bc_ad, P_bc, &
               0d0, dlnd00, 0d0, &
               0d0, dlnT00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, dlnR00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, P_bc*dlnP_bc_dL, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0)
            dlnR00 = dlnP_bc_dlnR
            if (.not. s% solver_use_lnR) dlnR00 = dlnR00/s% r(1)
            dlnT00 = dlnP_bc_dlnT
            if (.not. s% solver_use_lnT) dlnT00 = dlnT00/s% T(1)
            dlnd00 = dlnP_bc_dlnd
            if (.not. s% solver_use_lnd) dlnd00 = dlnd00/s% rho(1)
            call wrap(lnP_bc_ad, lnP_bc, &
               0d0, dlnd00, 0d0, &
               0d0, dlnT00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, dlnR00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, dlnP_bc_dL, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0)

            dlnT_bc_dlnT = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnT &
                  + dlnT_bc_dT0*dT0_dlnT
            dlnT_bc_dlnd = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnd &
                  + dlnT_bc_dT0*dT0_dlnd
            dlnT_bc_dL = dlnT_bc_dlnTsurf*dlnTsurf_dL + dlnT_bc_dT0*dT0_dL
            dlnT_bc_dlnR = dlnT_bc_dlnTsurf*dlnTsurf_dlnR + dlnT_bc_dT0*dT0_dlnR

            dlnR00 = T_bc*dlnT_bc_dlnR
            if (.not. s% solver_use_lnR) dlnR00 = dlnR00/s% r(1)
            dlnT00 = T_bc*dlnT_bc_dlnT
            if (.not. s% solver_use_lnT) dlnT00 = dlnT00/s% T(1)
            dlnd00 = T_bc*dlnT_bc_dlnd
            if (.not. s% solver_use_lnd) dlnd00 = dlnd00/s% rho(1)
            call wrap(T_bc_ad, T_bc, &
               0d0, dlnd00, 0d0, &
               0d0, dlnT00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, dlnR00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, T_bc*dlnT_bc_dL, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0)
            dlnR00 = dlnT_bc_dlnR
            if (.not. s% solver_use_lnR) dlnR00 = dlnR00/s% r(1)
            dlnT00 = dlnT_bc_dlnT
            if (.not. s% solver_use_lnT) dlnT00 = dlnT00/s% T(1)
            dlnd00 = dlnT_bc_dlnd
            if (.not. s% solver_use_lnd) dlnd00 = dlnd00/s% rho(1)
            call wrap(lnT_bc_ad, lnT_bc, &
               0d0, dlnd00, 0d0, &
               0d0, dlnT00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, dlnR00, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, dlnT_bc_dL, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0, &
               0d0, 0d0, 0d0)
            
         end subroutine get_PT_bc_ad


         subroutine set_Tsurf_BC(ierr)
            integer, intent(out) :: ierr
            logical :: test_partials
            type(auto_diff_real_star_order1) :: lnT1_ad
            include 'formats'
            !test_partials = (1 == s% solver_test_partials_k)
            test_partials = .false.
            ierr = 0           
            lnT1_ad = wrap_lnT_00(s,1)            
            resid_ad = lnT_bc_ad/lnT1_ad - 1d0
            s% equ(i_T_eqn, 1) = resid_ad%val
            if (test_partials) then
               s% solver_test_partials_val = 0
            end if
            if (skip_partials) return
            call save_eqn_residual_info( &
               s, 1, nvar, i_T_eqn, resid_ad, 'set_Tsurf_BC', ierr)           
            if (test_partials) then
               s% solver_test_partials_var = 0
               s% solver_test_partials_dval_dx = 0
               write(*,*) 'set_Tsurf_BC', s% solver_test_partials_var
            end if
         end subroutine set_Tsurf_BC

 
         subroutine set_Psurf_BC(ierr)
            integer, intent(out) :: ierr
            logical :: test_partials
            type(auto_diff_real_star_order1) :: lnP1_ad
            include 'formats'
            !test_partials = (1 == s% solver_test_partials_k)
            test_partials = .false.
            ierr = 0           
            lnP1_ad = wrap_lnP_00(s,1)            
            resid_ad = lnP_bc_ad/lnP1_ad - 1d0
            s% equ(i_P_eqn, 1) = resid_ad%val
            if (test_partials) then
               s% solver_test_partials_val = 0
            end if
            if (skip_partials) return
            call save_eqn_residual_info( &
               s, 1, nvar, i_P_eqn, resid_ad, 'set_Psurf_BC', ierr)           
            if (test_partials) then
               s% solver_test_partials_var = 0
               s% solver_test_partials_dval_dx = 0
               write(*,*) 'set_Psurf_BC', s% solver_test_partials_var
            end if
         end subroutine set_Psurf_BC
         

         subroutine set_momentum_BC(ierr)
            use hydro_riemann, only: do_surf_Riemann_dudt_eqn
            use hydro_momentum, only: do_surf_momentum_eqn
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0            
            if (s% u_flag) then
               call do_surf_Riemann_dudt_eqn( &
                  s, P_bc_ad, skip_partials, nvar, ierr)
            else
               call do_surf_momentum_eqn( &
                  s, P_bc_ad, skip_partials, nvar, ierr)
            end if            
         end subroutine set_momentum_BC


         subroutine set_compression_BC(ierr)
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: &
               rho1, rho2, lnd1, lnd2, dlnd1, dlnd2
            include 'formats'
            ! gradient of compression vanishes fixes density for cell 1
               ! d_dt(1/rho(1)) = d_dt(1/rho(2))  e.g., Grott, Chernigovski, Glatzel, 2005.
               ! -dlnd_dt(1)/rho(1) = -dlnd_dt(2)/rho(2)
               ! rho(2)*dlnd_dt(1) = rho(1)*dlnd_dt(2)
            ierr = 0            
            rho1 = wrap_d_00(s,1)
            rho2 = wrap_d_p1(s,1)
            lnd1 = wrap_lnd_00(s,1)
            lnd2 = wrap_lnd_p1(s,1)
            dlnd1 = lnd1 - s% lnd_start(1)
            dlnd2 = lnd2 - s% lnd_start(2)
            resid_ad = (rho2*dlnd1 - rho1*dlnd2)/s% dt            
            s% equ(i_P_eqn, 1) = resid_ad%val     
            if (skip_partials) return            
            call save_eqn_residual_info( &
               s, 1, nvar, i_P_eqn, resid_ad, 'set_compression_BC', ierr)           
         end subroutine set_compression_BC


         subroutine set_fixed_vsurf_outer_BC(ierr)
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: vsurf
            include 'formats'
            ierr = 0
            if (do_du_dt) then
               vsurf = wrap_u_00(s,1)
            else if (do_dv_dt) then
               vsurf = wrap_v_00(s,1)
            else
               ierr = -1
               write(*,*) 'set_fixed_vsurf_outer_BC requires u_flag or v_flag true'
               return
            end if 
            resid_ad = (vsurf - s% fixed_vsurf)/s% csound_start(1)
            s% equ(i_P_eqn,1) = resid_ad%val
            if (skip_partials) return
            call save_eqn_residual_info( &
               s, 1, nvar, i_P_eqn, resid_ad, 'set_fixed_vsurf_outer_BC', ierr)           
         end subroutine set_fixed_vsurf_outer_BC


         subroutine set_zero_dL_dm_BC(ierr)
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: L1, L2
            include 'formats'
            ierr = 0
            L1 = wrap_L_00(s,1)
            L2 = wrap_L_p1(s,1)
            resid_ad = L2/L1 - 1d0
            s% equ(i_T_eqn,1) = resid_ad%val
            if (skip_partials) return
            call save_eqn_residual_info( &
               s, 1, nvar, i_T_eqn, resid_ad, 'set_zero_dL_dm_BC', ierr)           
         end subroutine set_zero_dL_dm_BC


      end subroutine PT_eqns_surf


      ! only used for dlnT_dm equation
      subroutine eval_dlnPdm_qhse(s, k, & ! calculate the expected dlnPdm for HSE
            dlnPdm_qhse, Ppoint, ierr)
         use hydro_momentum, only: expected_HSE_grav_term
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: dlnPdm_qhse, Ppoint
         integer, intent(out) :: ierr

         real(dp) :: alfa
         type(auto_diff_real_star_order1) :: grav, area, P00, Pm1
         include 'formats'

         ierr = 0

         ! basic eqn is dP/dm = -G m / (4 pi r^4)
         ! divide by Ppoint to make it unitless

         ! for rotation, multiply gravity by factor fp.  MESA 2, eqn 22.
         
         call expected_HSE_grav_term(s, k, grav, area, ierr)
         if (ierr /= 0) return
         
         P00 = wrap_p_00(s,k)
         if (s% using_velocity_time_centering) P00 = 0.5d0*(P00 + s% P_start(k))
         
         if (k == 1) then
            Pm1 = 0d0
            Ppoint = P00
         else
            Pm1 = wrap_p_m1(s,k)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*P00 + (1d0-alfa)*Pm1
         end if
         
         dlnPdm_qhse = grav/(area*Ppoint) ! note that expected_HSE_grav_term is negative
                  
         if (is_bad(dlnPdm_qhse%val)) then
            ierr = -1
            s% retry_message = 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
            if (s% report_ierr) then
!$OMP critical (hydro_vars_crit1)
               write(*,*) 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
               stop
!$OMP end critical (hydro_vars_crit1)
            end if
            if (s% stop_for_bad_nums) then
               write(*,2) 'dlnPdm_qhse', k, dlnPdm_qhse
               stop 'eval_dlnPdm_qhse'
            end if
            return
         end if

      end subroutine eval_dlnPdm_qhse      


      integer function equ_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         equ_extra_profile_columns = s% nvar_hydro
      end function equ_extra_profile_columns


      subroutine equ_data_for_extra_profile_columns( &
            id, n, nz, names, vals, ierr)
         use star_def, only: maxlen_profile_column_name, star_info
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: i, k
         type (star_info), pointer :: s
         real(dp), dimension(:, :), pointer :: equ
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% nvar_hydro /= n) then
            write(*,3) 'nvar_hydro /= n', s% nvar_hydro, n
            stop 'equ_data_for_extra_profile_columns'
         end if
         do i=1,n
            do k=1,nz
               vals(k,i) = s% equ(i,k)
            end do
            names(i) = s% nameofequ(i)
         end do
      end subroutine equ_data_for_extra_profile_columns


      end module hydro_eqns

