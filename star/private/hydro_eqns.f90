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
         use hydro_tdc, only: do1_turbulent_energy_eqn, do1_tdc_L_eqn
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
            i_dalpha_RTI_dt, i_dln_cvpv0_dt, i_equ_w_div_wc, i_dj_rot_dt, i_dw_dt, &
            i, k, j, nvar_hydro, nz, op_err
         integer :: &
            i_lnd, i_lnR, i_lnT, i_lum, i_v, i_u, i_du, i_ln_cvpv0, i_w_div_wc, i_j_rot, &
            i_alpha_RTI, i_chem1, i_xh1, i_xhe4, kmax_equ(nvar), species
         real(dp) :: max_equ(nvar), dVARdot_dVAR, L_phot_old
         real(dp), dimension(:), pointer :: &
            L, lnR, lnP, lnT, energy
         logical :: v_flag, u_flag, conv_vel_flag, cv_flag, w_div_wc_flag, j_rot_flag, dump_for_debug, &
            do_chem, do_mix, do_struct_hydro, do_struct_thermo, &
            do_dlnd_dt, do_dv_dt, do_du_dt, do_dlnR_dt, &
            do_alpha_RTI, do_conv_vel, do_w_div_wc, do_j_rot, do_dlnE_dt, do_equL, do_dw_dt

         include 'formats'

         ierr = 0

         if (s% do_struct_hydro .and. s% u_flag .and. s% use_mass_corrections) &
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
         
         do_struct_hydro = s% do_struct_hydro
         do_struct_thermo = s% do_struct_thermo

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
         do_dw_dt = (i_dw_dt > 0 .and. i_dw_dt <= nvar)

         if (s% fill_arrays_with_NaNs) call set_nan(s% equ1)

         if (s% v_flag) s% v_residual(1) = 0

         if (.not. (do_struct_hydro .or. do_struct_thermo)) then

            stop '(.not. (do_struct_hydro .or. do_struct_thermo)) is probably not working'
            if (.not. do_chem) then
               write(*,*) 'bug: do_chem, do_struct_hydro, and do_struct_thermo'
               call mesa_error(__FILE__,__LINE__)
            end if
            ! hold structure constant while solve burn and/or mix
            do j=1,nvar_hydro
               call null_eqn(j)
            end do
            if (ierr == 0) &
               call do_chem_eqns(s, nvar, skip_partials, ierr)

         else ! solving structure equations

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
               if (do_struct_hydro) then
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
                        call do1_Riemann_momentum_eqn(s, k, -1d0, skip_partials, nvar, op_err)
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
               end if
               if (do_struct_thermo) then
                  if (do_dlnE_dt) then
                     call zero_eps_grav_and_partials(s, k)
                     call do1_energy_eqn(s, k, skip_partials, do_chem, nvar, op_err)
                     if (op_err /= 0) then
                        if (s% report_ierr) write(*,2) 'ierr in do1_energy_eqn', k
                        if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_energy_eqn'
                        ierr = op_err
                     end if
                  end if
                  if (do_dw_dt) then
                     call do1_turbulent_energy_eqn(s, k, skip_partials, nvar, op_err)
                     if (op_err /= 0) then
                        if (s% report_ierr) write(*,2) 'ierr in do1_turbulent_energy_eqn', k
                        if (len_trim(s% retry_message) == 0) s% retry_message = 'error in do1_turbulent_energy_eqn'
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
            
            if (ierr == 0 .and. nzlo == 1 .and. &
                  (do_struct_hydro .or. do_struct_thermo)) then
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

            if (.not. do_struct_hydro) then
               call dummy_eqn(i_dlnd_dt,i_lnR,nzlo,nzhi)
               call dummy_eqn(i_dv_dt,i_lnd,max(2,nzlo),nzhi)
               if (i_dlnR_dt /= 0) call dummy_eqn(i_dlnR_dt,i_v,nzlo,nzhi)
            end if

            if (.not. do_struct_thermo) then
               call dummy_eqn(i_dlnE_dt,i_lum,nzlo,nzhi)
               call dummy_eqn(i_equL,i_lnT,max(2,nzlo),nzhi)
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
            dVARdot_dVAR = s% dVARdot_dVAR

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
            i_dw_dt = s% i_dw_dt

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


         subroutine null_eqn(j)
            integer, intent(in) :: j
            integer :: k
            do k=nzlo,nzhi
               s% equ(j,k) = 0 ! s% xs(j,k) - s% xs_pre_pass(j,k)
               if (.not. skip_partials) call e00(s,j,j,k,nvar,1d0)
            end do
         end subroutine null_eqn


         subroutine dummy_eqn(j,i,nzlo,nzhi)
            integer, intent(in) :: j, i, nzlo, nzhi
            integer :: k
            do k=nzlo,nzhi
               s% equ(j,k) = 0
               if (.not. skip_partials) call e00(s,j,i,k,nvar,1d0)
            end do
         end subroutine dummy_eqn


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
                     if (ierr /= 0) then
                        if (s% report_ierr) write(*,2) 'failed in get_eos with xa_start_1', k
                        return
                     end if

                     ! fix up derivatives

                     if (debug .and. k == s% solver_test_partials_k) & ! .and. s% solver_iter == s% solver_test_partials_iter_number) &
                        write(*,2) 'res(i_lnE) - lnE_with_xa_start', j, res(i_lnE) - lnE_with_xa_start

                     s% dlnE_dxa_for_partials(j,k) = dres_dxa(i_lnE, j) + &
                        frac_without_dxa * (res(i_lnE) - lnE_with_xa_start) / dxa

                     s% dlnP_dxa_for_partials(j,k) = s% Pgas(k)*dres_dxa(i_lnPgas,j)/s% P(k) + &
                        frac_without_dxa * (s% Pgas(k)/s% P(k)) * (res(i_lnPgas) - lnPgas_with_xa_start) / dxa

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


      end subroutine eval_equ_for_solver


      subroutine do1_density_eqn(s, k, skip_partials, nvar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         real(dp) :: rp13, dm, rho, dlnd_dt, dr3, dequ_ddr3, ddr3_dlnd, &
            ddr3_dlnPgas_const_T, ddr3_dlnT_const_Pgas, res, &
            rho_div_rho_start, rho_dvA_dm, R200, R2p1, v00, vp1, &
            scale, dlnd_dt_base, expected, actual
         integer :: nz, i_dlnd_dt, i_lnd, i_lnR
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         nz = s% nz
         i_dlnd_dt = s% i_dlnd_dt
         i_lnd = s% i_lnd
         i_lnR = s% i_lnR
         dm = s% dm(k)
         rho = s% rho(k)

         dr3 = (dm/rho)/(pi4/3d0)

         if (k < nz) then
            rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
         else
            rp13 = s% R_center*s% R_center*s% R_center
         end if
         ! dm = four_thirds_pi*(r(k)**3 - rp13)*rho
         ! r(k)**3 = rp13 + (dm/rho)/four_thirds_pi = rp13 + dr3
         res = log(rp13 + dr3)
         s% equ(i_dlnd_dt, k) = s% lnR(k) - res/3d0

         s% lnd_residual(k) = s% equ(i_dlnd_dt, k)

         if (test_partials) then
            s% solver_test_partials_val = s% equ(i_dlnd_dt, k)
         end if

         if (skip_partials) return

         call e00(s, i_dlnd_dt, i_lnR, k, nvar, 1d0)
         if (k < nz) call ep1(s, i_dlnd_dt, i_lnR, k, nvar, -rp13/(rp13 + dr3))

         dequ_ddr3 = -one_third/(rp13 + dr3)
         ddr3_dlnd = -dr3
         call e00(s, i_dlnd_dt, s% i_lnd, k, nvar, dequ_ddr3*ddr3_dlnd)

         if (test_partials) then   
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = dequ_ddr3*ddr3_dlnd
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

         dt = s% dt

         i_dln_cvpv0_dt = s% i_dln_cvpv0_dt
         i_ln_cvpv0 = s% i_ln_cvpv0
         
         if (.false.) then ! don't change conv_vel
            s% equ(i_dln_cvpv0_dt, k) = 0d0
            if (skip_partials) return
            call e00(s,i_dln_cvpv0_dt,i_ln_cvpv0,k,nvar,1d0)
            return
         end if
         
         if (.false.) then ! force conv_vel == 0
            s% equ(i_dln_cvpv0_dt, k) = s% conv_vel(k)
            if (skip_partials) return
            call e00(s,i_dln_cvpv0_dt,i_ln_cvpv0,k,nvar,1d0)
            return
         end if

         D = s% conv_vel_D
         siglimit = s% conv_vel_siglimit
         if (k == s% nz) then
            sigmid_00 = 0d0
         else
            rc = 0.5d0*(s% r_start(k) + s% r_start(k+1))
            f = pi4*rc*rc*s% rho_start(k) ! gm/cm
            sigmid_00 = min(D*f*f, siglimit*dt)/s% dm(k) ! gm/s
         end if

         if (k == 1) then
            sigmid_m1 = 0d0
         else
            rc = 0.5d0*(s% r_start(k-1) + s% r_start(k))
            f = pi4*rc*rc*s% rho_start(k-1)
            sigmid_m1 = min(D*f*f, siglimit*dt)/s% dm(k-1) ! gm/s
         end if
         
         f = s% conv_vel_mix_factor/s% dm_bar(k)
         mix = -f*(sigmid_00 + sigmid_m1)*s% conv_vel(k)
         if (k > 1) mix = mix + f*sigmid_m1*s% conv_vel(k-1)
         if (k < s% nz) mix = mix + f*sigmid_00*s% conv_vel(k+1)
         dmix_dcvm1 = f*sigmid_m1
         dmix_dcv00 = -f*(sigmid_00 + sigmid_m1)
         dmix_dcvp1 = f*sigmid_00

         ! reduce tolerance for either tiny conv vel or when very near to mlt value
         if (s% conv_vel(k) < 1d0) then
            scale = 1d0/(max(s% conv_vel(k),1d-10))
         else if (s% conv_vel_use_mlt_vc_start .and. s% mlt_vc_start(k) > 0d0) then
            scale = 1d0/(max(abs(s% mlt_vc_start(k) - s% conv_vel(k))/s% mlt_vc_start(k),1d-10))
         else if (.not. s% conv_vel_use_mlt_vc_start .and. s% mlt_vc(k) > 0d0) then
            scale = 1d0/(max(abs(s% mlt_vc(k) - s% conv_vel(k))/s% mlt_vc(k),1d-10))
         else
            scale = 1d0 ! s% csound_start(k)
         end if
         
         if (k==1) then
            s% equ(i_dln_cvpv0_dt, 1) = &
               (log(s% conv_vel(1)+s% conv_vel_v0) - log(s% conv_vel(2)+s% conv_vel_v0))/scale
            if (test_partials) then
               write(*,*) "test partials!", s% equ(i_dln_cvpv0_dt, 1)
               s% solver_test_partials_val = s% equ(i_dln_cvpv0_dt, 1)
            end if
            if (skip_partials) return
            call e00(s, i_dln_cvpv0_dt, i_ln_cvpv0, 1, nvar, &
               1d0/scale)
            call ep1(s, i_dln_cvpv0_dt, i_ln_cvpv0, 1, nvar, &
               -1d0/scale)

            if (test_partials) then
               s% solver_test_partials_var = s% i_ln_cvpv0
               s% solver_test_partials_dval_dx = &
                     -1d0/scale
               write(*,*) 'do1_dln_cvpv0_dt_eqn', s% solver_test_partials_var
            end if

            return
         else if (k==s% nz) then
            s% equ(i_dln_cvpv0_dt, s% nz) = &
               (log(s% conv_vel(s% nz)+s% conv_vel_v0) - log(s% conv_vel(s% nz-1)+s% conv_vel_v0))/scale
            if (skip_partials) return
            call e00(s, i_dln_cvpv0_dt, i_ln_cvpv0, s% nz, nvar, &
               1d0/scale)
            call em1(s, i_dln_cvpv0_dt, i_ln_cvpv0, s% nz, nvar, &
               -1d0/scale)
            return
         end if

         tau = 1d50
         if(s% gradr_start(k) - s% gradL_start(k) < 0d0) then
            grav_start = s% cgrav(k)*s% m(k)/pow2(s% r_start(k))
            N2_start = -s% chiT_start(k)/s% chiRho_start(k) &
               *(s% gradT_start(k) - s% gradL_start(k))*grav_start/s% scale_height_start(k)
            if (N2_start > 0d0) then ! .and. s% mlt_vc_start(k) <= 0d0) then
               tau = max(1d0/sqrt(N2_start),1d-5)
            end if
            tau = min(tau, 1d50)
         end if
         dt_div_tau = dt/tau
         lambda = (s% scale_height_start(k)*s% mixing_length_alpha)
         dt_div_lambda = dt/lambda

         if (s% conv_vel_fully_lagrangian .or. s% q(k) <= s% q(s% k_const_mass)) then
            alfa = 0d0
         else if (s% q(k) >= s% q(s% k_below_const_q)) then
            alfa = 1d0
         else
            alfa = (s% q(k)-s% q(s% k_const_mass)) &
               /(s% q(s% k_below_const_q) - s% q(s% k_const_mass))
         end if 

         if (.not. s% conv_vel_fully_lagrangian .and. s% conv_vel_include_homologous_term) then
            if (alfa > 0d0 .and. s% mstar_dot <= 0d0) then
               d_ln_cvpv0_dq = log(s% conv_vel(k)+s% conv_vel_v0) - log(s% conv_vel(k+1)+s% conv_vel_v0)
               d_ln_cvpv0_dq = d_ln_cvpv0_dq/(s% q(k) - s% q(k+1))
            else if (alfa > 0d0 .and. s% mstar_dot > 0d0) then
               d_ln_cvpv0_dq = log(s% conv_vel(k-1)+s% conv_vel_v0) - log(s% conv_vel(k)+s% conv_vel_v0)
               d_ln_cvpv0_dq = d_ln_cvpv0_dq/(s% q(k-1) - s% q(k))
            else
               d_ln_cvpv0_dq = 0d0
            end if
         else  
            d_ln_cvpv0_dq = 0d0
         end if

         d_actual = dt*((1-alfa)*s% dln_cvpv0_dt(k) &
            + alfa*(s% dln_cvpv0_dt_const_q(k) - s% mstar_dot/s% m(1)*d_ln_cvpv0_dq))

         if (s% conv_vel_use_mlt_vc_start) then
            a = (pow2(s% mlt_vc_start(k)) - s% conv_vel(k)*s% conv_vel(k))*dt_div_lambda
         else
            a = (pow2(s% mlt_vc(k)) - s% conv_vel(k)*s% conv_vel(k))*dt_div_lambda
         end if
         inv = 1d0/(s% conv_vel(k)+s% conv_vel_v0)
         c = a - s% conv_vel(k)*dt_div_tau + mix*dt
         d_expected = c*inv

         s% equ(i_dln_cvpv0_dt, k) = (d_expected - d_actual)/scale


         if (is_bad(s% equ(i_dln_cvpv0_dt, k))) then
            ierr = -1
            return
!$omp critical (hydro_eqns_crit2)
            write(*,2) 'equ(i_dln_cvpv0_dt, k)', k, s% equ(i_dln_cvpv0_dt, k)
            write(*,2) 's% dln_cvpv0_dt(k)', k, s% dln_cvpv0_dt(k)
            write(*,2) 's% mlt_vc(k)', k, s% mlt_vc(k)
            write(*,2) 's% conv_vel(k)', k, s% conv_vel(k)
            write(*,2) 'scale', k, scale
            write(*,2) 'inv', k, inv
            write(*,2) 'd_expected', k, d_expected
            write(*,2) 'd_actual', k, d_actual
            write(*,2) 'c', k, c
            write(*,2) 'a', k, a
            write(*,2) 's% mlt_vc_start(k)', k, s% mlt_vc_start(k)
            write(*,2) 's% mlt_vc_start(k)', k, s% mlt_vc_start(k)
            write(*,2) 's% conv_vel(k)', k, s% conv_vel(k)
            write(*,2) 'dt_div_lambda', k, dt_div_lambda
            write(*,2) 'lambda', k, lambda
            write(*,2) 's% scale_height_start(k)', k, s% scale_height_start(k)
            write(*,2) 's% mixing_length_alpha', k, s% mixing_length_alpha
            write(*,2) 'tau', k, tau
            stop 'do1_dln_cvpv0_dt_eqn'
!$omp end critical (hydro_eqns_crit2)
         end if

         if (test_partials) then
            write(*,*) "test partials!", (d_expected - d_actual)/scale, s% mstar_dot/msun*secyer
            s% solver_test_partials_val = (d_expected - d_actual)/scale
         end if

         !if (k==669) then
         !   write(*,*) "check more", k, s% conv_vel(k), s% mlt_vc(k), s% mlt_vc_start(k),&
         !      (d_expected - d_actual)/scale, d_expected, d_actual, &
         !     (pow2(s% mlt_vc_start(k)) - pow2(s% conv_vel(k)))*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0), &
         !     - s% conv_vel(k)*dt_div_tau/(s% conv_vel(k)+s% conv_vel_v0), &
         !     mix*dt/(s% conv_vel(k)+s% conv_vel_v0), scale, tau, N2_start
         !end if

         if (skip_partials) return

         da_dvc = -2*s% conv_vel(k)*dt_div_lambda
         d_inv_dvc = -inv*inv
         dc_dvc = da_dvc - dt_div_tau + dmix_dcv00*dt
         d_d_expected_dvc = (dc_dvc*inv + c*d_inv_dvc)/inv ! /inv = 1/(d(log(cv+v0))_dvc) = cv+v0
         d_d_actual_dvc = dt*s% dVARdot_dVAR

         if (.not. s% conv_vel_fully_lagrangian .and. s% conv_vel_include_homologous_term) then
            if (s% mstar_dot <= 0d0) then
               call e00(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
                  (d_d_expected_dvc - d_d_actual_dvc)/scale &
                     +alfa*dt*s% mstar_dot/s% m(1)/(s% q(k) - s% q(k+1))/scale)
            else
               call e00(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
                  (d_d_expected_dvc - d_d_actual_dvc)/scale &
                     -alfa*dt*s% mstar_dot/s% m(1)/(s% q(k-1) - s% q(k))/scale)
            end if
         else
            call e00(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
               (d_d_expected_dvc - d_d_actual_dvc)/scale)
         end if

         if (.not. s% conv_vel_use_mlt_vc_start) then
            ! partials of mlt_vc*dt_div_tau/scale
            call e00(s, i_dln_cvpv0_dt, s% i_lnR, k, nvar, &
               2*s% mlt_vc(k)*s% d_mlt_vc_dlnR(k)*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0)/scale)

            call e00(s, i_dln_cvpv0_dt, s% i_lum, k, nvar, &
               2*s% mlt_vc(k)*s% d_mlt_vc_dL(k)*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0)/scale)

            call e00(s, i_dln_cvpv0_dt, s% i_lnd, k, nvar, &
               2*s% mlt_vc(k)*s% d_mlt_vc_dlnd00(k)*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0)/scale)

            call e00(s, i_dln_cvpv0_dt, s% i_lnT, k, nvar, &
               2*s% mlt_vc(k)*s% d_mlt_vc_dlnT00(k)*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0)/scale)
         end if
            
         if (k > 1) then
            call em1(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
               dmix_dcvm1*dt/(s% conv_vel(k)+s% conv_vel_v0)*(s% conv_vel(k-1)+s% conv_vel_v0)/scale)
            if (.not. s% conv_vel_fully_lagrangian .and. s% conv_vel_include_homologous_term) then
               if (s% mstar_dot > 0d0) then
                  call em1(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
                     alfa*dt*s% mstar_dot/s% m(1)/(s% q(k-1) - s% q(k))/scale)
               end if
            end if
            if (.not. s% conv_vel_use_mlt_vc_start) then
               ! partials of mlt_vc*dt_div_tau/scale
               call e00(s, i_dln_cvpv0_dt, s% i_lnd, k, nvar, &
                  2*s% mlt_vc(k)*s% d_mlt_vc_dlndm1(k)*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0)/scale)
               call e00(s, i_dln_cvpv0_dt, s% i_lnT, k, nvar, &
                  2*s% mlt_vc(k)*s% d_mlt_vc_dlnTm1(k)*dt_div_lambda/(s% conv_vel(k)+s% conv_vel_v0)/scale)
            end if
         end if
         
         if (k < s% nz) then
            call ep1(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
               dmix_dcvp1*dt/(s% conv_vel(k)+s% conv_vel_v0)*(s% conv_vel(k+1)+s% conv_vel_v0)/scale)
            if (.not. s% conv_vel_fully_lagrangian .and. s% conv_vel_include_homologous_term) then
               if (s% mstar_dot <= 0d0) then
                  call ep1(s, i_dln_cvpv0_dt, i_ln_cvpv0, k, nvar, &
                     -alfa*dt*s% mstar_dot/s% m(1)/(s% q(k) - s% q(k+1))/scale)
               end if
            end if
         end if

         if (test_partials) then
            s% solver_test_partials_var = s% i_ln_cvpv0
            s% solver_test_partials_dval_dx = &
                  (d_d_expected_dvc - d_d_actual_dvc)/scale
            write(*,*) 'do1_dln_cvpv0_dt_eqn', s% solver_test_partials_var
         end if         
                  
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

         scale = 1d0
         
         !test_partials = (k == s% solver_test_partials_k-1)
         test_partials = .false.

         i_equ_w_div_wc = s% i_equ_w_div_wc
         i_w_div_wc = s% i_w_div_wc
         wwc = s% w_div_wcrit_max
         A = 1d0-0.1076d0*pow4(wwc)-0.2336d0*pow6(wwc)-0.5583d0*log(1d0-pow4(wwc))
         C = 1d0+17d0/60d0*pow2(wwc)-0.3436d0*pow4(wwc)-0.4055d0*pow6(wwc)-0.9277d0*log(1d0-pow4(wwc))
         jr_lim1 = two_thirds*wwc*C/A

         wwc = s% w_div_wcrit_max2
         A = 1d0-0.1076d0*pow4(wwc)-0.2336d0*pow6(wwc)-0.5583d0*log(1d0-pow4(wwc))
         C = 1d0+17d0/60d0*pow2(wwc)-0.3436d0*pow4(wwc)-0.4055d0*pow6(wwc)-0.9277d0*log(1d0-pow4(wwc))
         jr_lim2 = two_thirds*wwc*C/A

         wwc = s% w_div_w_crit_roche(k)
         A = 1d0-0.1076d0*pow4(wwc)-0.2336d0*pow6(wwc)-0.5583d0*log(1d0-pow4(wwc))
         C = 1d0+17d0/60d0*pow2(wwc)-0.3436d0*pow4(wwc)-0.4055d0*pow6(wwc)-0.9277d0*log(1d0-pow4(wwc))

         jrot_ratio = s% j_rot(k)/sqrt(s% cgrav(k)*s% m(k)*s% r(k))
         if (abs(jrot_ratio) > jr_lim1) then
            sigmoid_jrot_ratio = 2*(jr_lim2-jr_lim1)/(1+exp(-2*(abs(jrot_ratio)-jr_lim1)/(jr_lim2-jr_lim1)))-jr_lim2+2*jr_lim1
            if (jrot_ratio < 0d0) then
               sigmoid_jrot_ratio = -sigmoid_jrot_ratio
            end if
            !if (k==79) write(*,*) "check k sigmoid",k,sigmoid_jrot_ratio,jr_lim1, jr_lim2, two_thirds*wwc*C/A, wwc
            s% equ(i_equ_w_div_wc, k) = (sigmoid_jrot_ratio - two_thirds*wwc*C/A)/scale
         else
            !if (k==79) write(*,*) "check k normal",k,jrot_ratio,jr_lim1,jr_lim2,two_thirds*wwc*C/A, wwc
            s% equ(i_equ_w_div_wc, k) = (jrot_ratio - two_thirds*wwc*C/A)/scale
         end if

         if (is_bad(s% equ(i_equ_w_div_wc, k))) then
            ierr = -1
            return
         end if

         if (test_partials) then
            s% solver_test_partials_val = (jrot_ratio-2d0/3d0*wwc*C/A)/scale
         end if

         if (skip_partials) return

         dA_dw = -4d0*0.1076d0*pow3(wwc)-6d0*0.2336d0*pow5(wwc) &
            -0.5583d0/(1d0-pow4(wwc))*(-4d0*pow3(wwc))
         dC_dw = 17d0/30d0*wwc-4d0*0.3436d0*pow3(wwc)-6d0*0.4055d0*pow5(wwc)-0.9277d0/(1d0-pow4(wwc))*(-4d0*pow3(wwc))

         call e00(s, i_equ_w_div_wc, i_w_div_wc, k, nvar, &
            -two_thirds*(C/A+wwc*(dC_dw/A - C*dA_dw/pow2(A))))

         if (jrot_ratio > jr_lim1) then
            d_sigmoid_jrot_ratio_dx = -2d0*(jr_lim2-jr_lim1)/pow2(1+exp(-2*(abs(jrot_ratio)-jr_lim1)/(jr_lim2-jr_lim1)))&
               *(-2d0/(jr_lim2-jr_lim1))*exp(-2d0*(abs(jrot_ratio)-jr_lim1)/(jr_lim2-jr_lim1))
            if (jrot_ratio < 0d0) then
               d_sigmoid_jrot_ratio_dx = -d_sigmoid_jrot_ratio_dx
            end if
            call e00(s, i_equ_w_div_wc, s% i_lnR, k, nvar, &
               d_sigmoid_jrot_ratio_dx*(-0.5d0*jrot_ratio))
            if (s% j_rot_flag) then
               call e00(s, i_equ_w_div_wc, s% i_j_rot, k, nvar, &
                  d_sigmoid_jrot_ratio_dx/sqrt(s% cgrav(k)*s% m(k)*s% r(k)))
            end if
         else
            call e00(s, i_equ_w_div_wc, s% i_lnR, k, nvar, &
               -0.5d0*jrot_ratio)
            if (s% j_rot_flag) then
               call e00(s, i_equ_w_div_wc, s% i_j_rot, k, nvar, &
                  1/sqrt(s% cgrav(k)*s% m(k)*s% r(k)))
            end if
         end if

         if (test_partials) then
            s% solver_test_partials_var = s% i_w_div_wc
            s% solver_test_partials_dval_dx = 1d0
            write(*,*) 'do1_w_div_wc_eqn', s% solver_test_partials_var
         end if
         
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

         scale = 1d4*max(1d2*sqrt(s% cgrav(k)*s% m(k)*s%r_start(k))/s%dt,&
            s% total_abs_angular_momentum/(s% dm(k)*s% dt*s% nz))
         
         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         i_dj_rot_dt = s% i_dj_rot_dt
         i_j_rot = s% i_j_rot

         if (k>1) then
            Fm1 = s% j_flux(k-1)
            dFm1_dwm1 = s% dj_flux_dw00(k-1)
            dFm1_dw00 = s% dj_flux_dwp1(k-1)
            dFm1_djm1 = s% dj_flux_dj00(k-1)
            dFm1_dj00 = s% dj_flux_djp1(k-1)
            dFm1_dlnrm1 = s% dj_flux_dlnr00(k-1)
            dFm1_dlnr00 = s% dj_flux_dlnrp1(k-1)
            dFm1_dlndm1 = s% dj_flux_dlnd(k-1)
         else
            Fm1 = 0d0
            dFm1_dwm1 = 0d0 
            dFm1_dw00 = 0d0
            dFm1_djm1 = 0d0
            dFm1_dj00 = 0d0
            dFm1_dlnrm1 = 0d0
            dFm1_dlnr00 = 0d0
            dFm1_dlndm1 = 0d0
         end if

         F00 = s% j_flux(k)
         dF00_dw00 = s% dj_flux_dw00(k)
         dF00_dwp1 = s% dj_flux_dwp1(k)
         dF00_dj00 = s% dj_flux_dj00(k)
         dF00_djp1 = s% dj_flux_djp1(k)
         dF00_dlnr00 = s% dj_flux_dlnr00(k)
         dF00_dlnrp1 = s% dj_flux_dlnrp1(k)
         dF00_dlnd00 = s% dj_flux_dlnd(k)

         s% dj_rot_dt(k) = (s% j_rot(k)-s% j_rot_start(k))/s% dt ! for some reason not working from hydro_mtx
         !equ(i_dj_rot_dt, k) = s% dj_rot_dt(k)/scale!(s% dj_rot_dt(k)-(Fplus-Fminus)/s% dm_bar(k))/scale
         s% equ(i_dj_rot_dt, k) = (s% dj_rot_dt(k)+(Fm1-F00)/s% dm_bar(k)-s% extra_jdot(k))/scale

         !if (k==171) then
         !   write(*,*) "check eqn", k, s% equ(i_dj_rot_dt, k), s%dj_rot_dt(k), &
         !           (Fm1-F00)/s% dm_bar(k), s% extra_jdot(k)
         !end if

         if (is_bad(s% equ(i_dj_rot_dt, k))) then
            ierr = -1
            return
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% j_flux(19)
         end if

         if (skip_partials) return

         call e00(s, i_dj_rot_dt, s% i_w_div_wc, k, nvar, &
            +(dFm1_dw00-dF00_dw00)/s% dm_bar(k)/scale)
         call e00(s, i_dj_rot_dt, s% i_j_rot, k, nvar, &
            (1/s% dt+(dFm1_dj00-dF00_dj00)/s% dm_bar(k))/scale)
         call e00(s, i_dj_rot_dt, s% i_lnr, k, nvar, &
            +(dFm1_dlnr00-dF00_dlnr00)/s% dm_bar(k)/scale)
         call e00(s, i_dj_rot_dt, s% i_lnd, k, nvar, &
            -dF00_dlnd00/s% dm_bar(k)/scale)
         
         if (k > 1) then
            call em1(s, i_dj_rot_dt, s% i_w_div_wc, k, nvar, &
               +dFm1_dwm1/s% dm_bar(k)/scale)
            call em1(s, i_dj_rot_dt, s% i_j_rot, k, nvar, &
               +dFm1_djm1/s% dm_bar(k)/scale)
            call em1(s, i_dj_rot_dt, s% i_lnr, k, nvar, &
               +dFm1_dlnrm1/s% dm_bar(k)/scale)
            call em1(s, i_dj_rot_dt, s% i_lnd, k, nvar, &
               +dFm1_dlndm1/s% dm_bar(k)/scale)
         end if

         if (k < s% nz) then
            call ep1(s, i_dj_rot_dt, s% i_w_div_wc, k, nvar, &
               -dF00_dwp1/s% dm_bar(k)/scale)
            call ep1(s, i_dj_rot_dt, s% i_j_rot, k, nvar, &
               -dF00_djp1/s% dm_bar(k)/scale)
            call ep1(s, i_dj_rot_dt, s% i_lnr, k, nvar, &
               -dF00_dlnrp1/s% dm_bar(k)/scale)
         end if

         if (test_partials) then
            s% solver_test_partials_var = s% i_j_rot
            s% solver_test_partials_dval_dx =  s% dj_flux_dj00(19)
            write(*,*) 'do1_dj_rot_dt_eqn', s% solver_test_partials_var
         end if
         
      end subroutine do1_dj_rot_dt_eqn


      ! just relate L_rad to T gradient.
      ! d_P_rad/dm = -<opacity_face>*L_rad/(clight*area^2) -- see, e.g., K&W (5.12)
      ! P_rad = (1/3)*crad*T^4
      ! d_P_rad/dm = (crad/3)*(T(k-1)^4 - T(k)^4)/dm_bar
      ! L_rad = L - L_non_rad, L_non_rad = L_start - L_rad_start
      ! L_rad_start = (-d_P_rad/dm_bar*clight*area^2/<opacity_face>)_start
      subroutine do1_alt_dlnT_dm_eqn(s, k, skip_partials, nvar, ierr)
         use eos_def
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         real(dp) :: alfa, beta, opacity_face, area, area2, T4_m1, T4_00, &
            L, P_rad_start, P_rad_m1, P_rad_00, scale, dm_bar, L_rad, L_non_rad, &
            d_P_rad_expected, d_P_rad_actual, r, dr_dL00, dr_dw00, &
            d_kap_dlnT00, d_expected_dlnT00, d_actual_dlnT00, dr_dlnT00, &
            d_kap_dlnTm1, d_expected_dlnTm1, d_actual_dlnTm1, dr_dlnTm1, &
            d_kap_dlnd00, d_expected_dlnd00, dr_dlnd00, &
            d_kap_dlndm1, d_expected_dlndm1, dr_dlndm1, &
            dr_dlnT00_const_Pgas, dr_dlnTm1_const_Pgas, dr_dlnR00, dr_dln_cvpv0, &
            dr_dlnPgas00_const_T, dr_dlnPgasm1_const_T, d_Lrad_dL, d_Lrad_dw, &
            d_Lrad_dlnT00, d_Lrad_dlnTm1, d_Lrad_dlnd00, d_Lrad_dlndm1, &
            d_Lrad_dlnR, d_Lrad_dln_cvpv0, gradr2, dlnd, d_dlnd, d_dlnT, &
            d_dlnd_const_E, d_dlnE_const_Rho, d_T400_dlnT00, d_T4m1_dlnTm1, &
            d_expected_dL, d_expected_dw, d_expected_dlnR, d_expected_dln_cvpv0, cs
         integer :: i_equL, i
         logical :: dbg
         logical :: test_partials

         include 'formats'
         ierr = 0
         i_equL = s% i_equL
         if (i_equL == 0) return

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         dbg = .false.

         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         beta = 1d0 - alfa

         scale = s% energy_start(k)*s% rho_start(k)
         dm_bar = s% dm_bar(k)
         L = s% L(k)
         area = pi4*s% r(k)*s% r(k); area2 = area*area

         ! set Lrad partials for usual case; revise below if necessary.
         d_Lrad_dL = 1d0
         d_Lrad_dlnT00 = 0d0
         d_Lrad_dlnTm1 = 0d0
         d_Lrad_dlnR = 0d0
         d_Lrad_dlnd00 = 0d0
         d_Lrad_dlndm1 = 0d0
         d_Lrad_dln_cvpv0 = 0d0
         d_Lrad_dw = 0d0

         if (s% use_dPrad_dm_form_of_T_gradient_eqn) then
            if ((check_flag_and_val(s% conv_vel_flag, s% conv_vel, k)) .or. &
                  (.not. s% conv_vel_flag .and. s% lnT(k)/ln10 <= s% max_logT_for_mlt &
                  .and. s% mixing_type(k) == convective_mixing .and. s% gradr(k) > 0d0 &
                  .and. abs(s% gradr(k) - s% gradT(k)) > abs(s% gradr(k))*1d-5)) then
               L_rad = L*s% gradT(k)/s% gradr(k) ! C&G 14.109
               gradr2 = s% gradr(k)*s% gradr(k)
               d_Lrad_dL = s% gradT(k)/s% gradr(k) + L*( &
                  s% d_gradT_dL(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dL(k)/gradr2)
               d_Lrad_dw = s% gradT(k)/s% gradr(k) + L*( &
                  s% d_gradT_dw_div_wc(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dw_div_wc(k)/gradr2)
               d_Lrad_dlnT00 = L*( &
                  s% d_gradT_dlnT00(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dlnT00(k)/gradr2)
               d_Lrad_dlnTm1 = L*( &
                  s% d_gradT_dlnTm1(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dlnTm1(k)/gradr2)
               d_Lrad_dlnR = L*( &
                  s% d_gradT_dlnR(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dlnR(k)/gradr2)
               d_Lrad_dlnd00 = L*( &
                  s% d_gradT_dlnd00(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dlnd00(k)/gradr2)
               d_Lrad_dlndm1 = L*( &
                  s% d_gradT_dlndm1(k)/s% gradr(k) &
                  - s% gradT(k)*s% d_gradr_dlndm1(k)/gradr2)
               if (s% conv_vel_flag) &
                  d_Lrad_dln_cvpv0 = L*s% d_gradT_dln_cvpv0(k)/s% gradr(k)
            else
               L_rad = L
            end if
            L_non_rad = L - L_rad
         else
            ierr = -1
            return
         end if

         ! calculate expected d_P_rad from current L_rad
         opacity_face = alfa*s% opacity(k) + beta*s% opacity(k-1)
         
         if (opacity_face < s% min_kap_for_dPrad_dm_eqn) then
            opacity_face = s% min_kap_for_dPrad_dm_eqn
            d_kap_dlnd00 = 0d0
            d_kap_dlndm1 = 0d0
            d_kap_dlnT00 = 0d0
            d_kap_dlnTm1 = 0d0
         else
            d_kap_dlnd00 = alfa*s% d_opacity_dlnd(k)
            d_kap_dlndm1 = beta*s% d_opacity_dlnd(k-1)
            d_kap_dlnT00 = alfa*s% d_opacity_dlnT(k)
            d_kap_dlnTm1 = beta*s% d_opacity_dlnT(k-1)
         end if
                  
         d_P_rad_expected = -dm_bar*opacity_face*L_rad/(clight*area2)
         
         d_expected_dlnR = -d_P_rad_expected*(4 + d_Lrad_dlnR)
         d_expected_dL = -dm_bar*opacity_face*d_Lrad_dL/(clight*area2)
         d_expected_dln_cvpv0 = -dm_bar*opacity_face*d_Lrad_dln_cvpv0/(clight*area2)
         d_expected_dw = -dm_bar*opacity_face*d_Lrad_dw/(clight*area2)
         
         d_expected_dlnT00 = -dm_bar/(clight*area2)*( &
            d_kap_dlnT00*L_rad + opacity_face*d_Lrad_dlnT00)
         d_expected_dlnTm1 = -dm_bar/(clight*area2)*( &
            d_kap_dlnTm1*L_rad + opacity_face*d_Lrad_dlnTm1)
            
         d_expected_dlnd00 = -dm_bar/(clight*area2)*( &
            d_kap_dlnd00*L_rad + opacity_face*d_Lrad_dlnd00)
         d_expected_dlndm1 = -dm_bar/(clight*area2)*( &
            d_kap_dlndm1*L_rad + opacity_face*d_Lrad_dlndm1)
         
         T4_m1 = s% T(k-1)*s% T(k-1)*s% T(k-1)*s% T(k-1)
         d_T4m1_dlnTm1 = 4d0*T4_m1
         
         T4_00 = s% T(k)*s% T(k)*s% T(k)*s% T(k)
         d_T400_dlnT00 = 4d0*T4_00

         !d_P_rad_expected = d_P_rad_expected*s% gradr_factor(k) !TODO(Pablo): check this
         

         ! calculate actual d_P_rad in current model
         P_rad_m1 = (crad/3d0)*T4_m1
         P_rad_00 = (crad/3d0)*T4_00
         d_P_rad_actual = P_rad_m1 - P_rad_00
         d_actual_dlnT00 = -4d0*P_rad_00
         d_actual_dlnTm1 = 4d0*P_rad_m1
         
         ! residual
         r = (d_P_rad_expected - d_P_rad_actual)/scale 
         s% equ(i_equL, k) = r
         s% equL_residual(k) = r

         dr_dlnTm1 = (d_expected_dlnTm1 - d_actual_dlnTm1)/scale
         dr_dlnT00 = (d_expected_dlnT00 - d_actual_dlnT00)/scale
         dr_dlndm1 = d_expected_dlndm1/scale
         dr_dlnd00 = d_expected_dlnd00/scale
         
         if (is_bad(r)) then
!$OMP critical (star_alt_dlntdm_bad_num)
            write(*,2) 'r', k, r
            write(*,2) 'd_P_rad_expected', k, d_P_rad_expected
            write(*,2) 'd_P_rad_actual', k, d_P_rad_actual
            write(*,2) 'scale', k, scale
            if (s% stop_for_bad_nums) stop 'do1_alt_dlnT_dm_eqn'
!$OMP end critical (star_alt_dlntdm_bad_num)
         end if

         if (test_partials) then
            s% solver_test_partials_val = s% gradT(k)
         end if
       
         if (k == -2) then
            write(*,*) 'skip_partials', skip_partials
            write(*,2) 'r', k, r
            stop 'do1_alt_dlnT_dm_eqn'
         end if

         if (skip_partials) return

         dr_dL00 = d_expected_dL/scale
         if (is_bad(dr_dL00)) then
!$OMP critical (star_alt_dlntdm_dl00_bad_num)
            write(*,2) 'dr_dL00', k, dr_dL00
            write(*,2) 'dm_bar', k, dm_bar
            write(*,2) 'opacity_face', k, opacity_face
            write(*,2) 'd_Lrad_dL', k, d_Lrad_dL
            write(*,2) 'area2', k, area2
            write(*,2) 'scale', k, scale
            write(*,2) 's% gradT(k)', k, s% gradT(k)
            write(*,2) 's% gradr(k)', k, s% gradr(k)
            write(*,2) 's% d_gradT_dL(k)', k, s% d_gradT_dL(k)
            write(*,2) 's% d_gradr_dL(k)', k, s% d_gradr_dL(k)
            if (s% stop_for_bad_nums) stop 'do1_alt_dlnT_dm_eqn'
!$OMP end critical (star_alt_dlntdm_dl00_bad_num)
         end if

         if(s% i_w_div_wc /= 0) then
            dr_dw00 = d_expected_dw/scale
            if (is_bad(dr_dw00)) then
!$OMP critical (star_alt_dlntdm_dw00_bad_num)
               write(*,2) 'dr_dw00', k, dr_dw00
               write(*,2) 'dm_bar', k, dm_bar
               write(*,2) 'opacity_face', k, opacity_face
               write(*,2) 'd_Lrad_dw', k, d_Lrad_dw
               write(*,2) 'area2', k, area2
               write(*,2) 'scale', k, scale
               write(*,2) 's% gradT(k)', k, s% gradT(k)
               write(*,2) 's% gradr(k)', k, s% gradr(k)
               write(*,2) 's% d_gradT_dw_div_wc(k)', k, s% d_gradT_dw_div_wc(k)
               if (s% stop_for_bad_nums) stop 'do1_alt_dlnT_dm_eqn'
!$OMP end critical (star_alt_dlntdm_dw00_bad_num)
            end if
         end if

         if (s% i_lum /= 0) & ! d/d_L00
            call e00(s, i_equL, s% i_lum, k, nvar, dr_dL00)

         if (s% i_w_div_wc /= 0) & ! d/d_w00
            call e00(s, i_equL, s% i_w_div_wc, k, nvar, dr_dw00)
            
         if (s% conv_vel_flag) then
            dr_dln_cvpv0 = -dm_bar*opacity_face*d_Lrad_dln_cvpv0/(clight*area2)/scale
            call e00(s, i_equL, s% i_ln_cvpv0, k, nvar, dr_dln_cvpv0)
         end if

         call e00(s, i_equL, s% i_lnT, k, nvar, dr_dlnT00)
         call em1(s, i_equL, s% i_lnT, k, nvar, dr_dlnTm1)

         if (s% do_struct_hydro) then ! partials wrt lnR and lnd
            dr_dlnR00 = &
               -4*d_P_rad_expected/scale &
               -dm_bar*opacity_face*d_Lrad_dlnR/(clight*area2)/scale
            call e00(s, i_equL, s% i_lnR, k, nvar, dr_dlnR00) ! d/dlnR
            call e00(s, i_equL, s% i_lnd, k, nvar, dr_dlnd00)
            call em1(s, i_equL, s% i_lnd, k, nvar, dr_dlndm1)
         end if

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = s% d_gradT_dlnT00(k)
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
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         real(dp) :: alfa, beta, r, dr_dL00, dr_dlnR00, dr_dlnd00, dr_dlnT00, &
            dr_dlndm1, dr_dlnTm1, dr_dw00, m, dlnPdm, &
            d_dlnPdm_dlnR, d_dlnPdm_dL, d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            dP_dlnPgas00_const_T, dP_dlnPgasm1_const_T, &
            dP_dlnT00_const_Pgas, dP_dlnTm1_const_Pgas, &
            d_dlnPdm_dw, &
            gradT, d_gradT_dL, d_gradT_dw, d_gradT_dlnR, d_gradT_du, &
            d_gradT_dlnd00, d_gradT_dlndm1, &
            d_gradT_dlnT00, d_gradT_dlnTm1, d_gradT_dln_cvpv0, &
            Ppoint, scale, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            d_dlnTdm_dLum, d_dlnTdm_dlnR, d_dlnTdm_dlnd00, &
            d_dlnTdm_dlnT00, d_dlnTdm_dlndm1, d_dlnTdm_dlnTm1, &
            d_dlnTdm_dw, &
            delm, T00, Tm1, dT, Tpoint, &
            lnTdiff, d_lnTdiff_dlnT00, d_lnTdiff_dlnTm1, diff, diff_old, &
            d_gradT_dlnT00_const_Pgas, d_gradT_dlnTm1_const_Pgas, &
            d_gradT_dlnPgas00_const_T, d_gradT_dlnPgasm1_const_T, &
            d_dlnTdm_dlnPgas00_const_T, d_dlnTdm_dlnPgasm1_const_T, &
            d_dlnTdm_dlnT00_const_Pgas, d_dlnTdm_dlnTm1_const_Pgas, &
            dr_dlnPgas00_const_T, dr_dlnPgasm1_const_T, dr_dln_cvpv0, &
            dr_dlnT00_const_Pgas, dr_dlnTm1_const_Pgas, dlnTdm
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
         
         call eval_dlnPdm_qhse(s, k, m, &
            dlnPdm, d_dlnPdm_dlnR, d_dlnPdm_dL, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            d_dlnPdm_dw, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            ierr)
         if (ierr /= 0) return

         call eval_gradT_info( &
            s, k, gradT, d_gradT_dL, d_gradT_dw, d_gradT_dlnR, d_gradT_du, &
            d_gradT_dlnd00, d_gradT_dlndm1, &
            d_gradT_dlnT00, d_gradT_dlnTm1, &
            d_gradT_dlnT00_const_Pgas, d_gradT_dlnTm1_const_Pgas, &
            d_gradT_dlnPgas00_const_T, d_gradT_dlnPgasm1_const_T, &
            d_gradT_dln_cvpv0, ierr)
         if (ierr /= 0) return

         dlnTdm = dlnPdm*gradT
         s% dlnT_dm_expected(k) = dlnTdm

         delm = (s% dm(k) + s% dm(k-1))/2
         Tm1 = s% T(k-1)
         alfa = s% dm(k-1)/(s% dm(k-1) + s% dm(k))
         beta = 1 - alfa

         T00 = s% T(k)
         dT = Tm1 - T00
         Tpoint = alfa*T00 + beta*Tm1
         lnTdiff = dT/Tpoint ! use this in place of lnT(k-1)-lnT(k)
         
         scale = 1d0
         r = (delm*s% dlnT_dm_expected(k) - lnTdiff)*scale
         s% equ(i_equL, k) = r
         s% equL_residual(k) = s% equ(i_equL,k)

         if (k == s% trace_k) then
            write(*,5) 'i_equL', k, s% solver_iter, s% solver_adjust_iter, &
               s% model_number, s% equ(i_equL, k)
         end if

         if (is_bad(s% equ(i_equL, k))) then
            ierr = -1
            if (s% report_ierr) write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            if (s% stop_for_bad_nums) stop 'hydro eqns'
            return
            write(*,2) 'equ(i_equL, k)', k, s% equ(i_equL, k)
            write(*,2) 'lnTdiff', k, lnTdiff
            write(*,2) 'delm', k, delm
            write(*,2) 'dlnT_dm_expected(k)', k, s% dlnT_dm_expected(k)
            write(*,2) 'dlnPdm', k, dlnPdm
            write(*,2) 'gradT', k, gradT
            stop 'i_equL'
         end if

         if (test_partials) then
            s% solver_test_partials_val = gradT ! s% equ(i_equL,k)
         end if

         if (skip_partials) return

         d_lnTdiff_dlnTm1 = T00*Tm1/(Tpoint*Tpoint)
         d_lnTdiff_dlnT00 = -d_lnTdiff_dlnTm1

         d_dlnTdm_dLum = dlnPdm*d_gradT_dL + d_dlnPdm_dL*gradT
         dr_dL00 = delm*d_dlnTdm_dLum*scale

         d_dlnTdm_dlnR = dlnPdm*d_gradT_dlnR + d_dlnPdm_dlnR*gradT
         dr_dlnR00 = delm*d_dlnTdm_dlnR*scale

         if (s% w_div_wc_flag) then
            !d_dlnTdm_dlnw = dlnPdm*d_gradT_dlnR + d_dlnPdm_dlnR*gradT !TODO gradT rotation correction
            d_dlnTdm_dw = d_dlnPdm_dw*gradT + dlnPdm*d_gradT_dw
            dr_dw00 = delm*d_dlnTdm_dw*scale
            call e00(s, i_equL, s% i_w_div_wc, k, nvar, dr_dw00)
         end if

         if (s% conv_vel_flag) then
            dr_dln_cvpv0 = delm*dlnPdm*d_gradT_dln_cvpv0*scale
            call e00(s, i_equL, s% i_ln_cvpv0, k, nvar, dr_dln_cvpv0)
         end if

         if (s% i_lum /= 0) &
            call e00(s, i_equL, s% i_lum, k, nvar, dr_dL00) ! d/d_L00

         d_dlnTdm_dlnT00 = dlnPdm*d_gradT_dlnT00 + d_dlnPdm_dlnT00*gradT
         dr_dlnT00 = (delm*d_dlnTdm_dlnT00 - d_lnTdiff_dlnT00)*scale
         d_dlnTdm_dlnTm1 = dlnPdm*d_gradT_dlnTm1 + d_dlnPdm_dlnTm1*gradT
         dr_dlnTm1 = (delm*d_dlnTdm_dlnTm1 - d_lnTdiff_dlnTm1)*scale
         
         d_dlnTdm_dlnd00 = dlnPdm*d_gradT_dlnd00 + d_dlnPdm_dlnd00*gradT
         dr_dlnd00 = delm*d_dlnTdm_dlnd00*scale
         d_dlnTdm_dlndm1 = dlnPdm*d_gradT_dlndm1 + d_dlnPdm_dlndm1*gradT
         dr_dlndm1 = delm*d_dlnTdm_dlndm1*scale
      
         call e00(s, i_equL, s% i_lnT, k, nvar, dr_dlnT00)
         call em1(s, i_equL, s% i_lnT, k, nvar, dr_dlnTm1)

         if (s% do_struct_hydro) then
            call e00(s, i_equL, s% i_lnR, k, nvar, dr_dlnR00) ! d/dlnR
            call e00(s, i_equL, s% i_lnd, k, nvar, dr_dlnd00)
            call em1(s, i_equL, s% i_lnd, k, nvar, dr_dlndm1)
         end if

         !dr_dL00       1e-12
         !dr_dlnR00     4e-12
         !dr_dln_cvpv0
         !dr_dlnT00     3e-12
         !dr_dlnd00     1e-12
         !dr_dlnTm1
         !dr_dlndm1
         
         ! gradT    d_gradT_dlnd00
         ! dlnPdm   d_dlnPdm_dlnd00  1d-12
         ! dr_dL00 = delm*d_dlnTdm_dLum*scale

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = d_gradT_dlnR ! dr_dlnR00
            write(*,*) 'do1_dlnT_dm_eqn', s% solver_test_partials_var
            !write(*,*) 's% conv_vel_flag', s% conv_vel_flag
            !write(*,3) 'mixing_type gradT gradr', k, s% mixing_type(k), gradT, s% gradr(k)
         end if

      end subroutine do1_dlnT_dm_eqn


      subroutine PT_eqns_surf( &
            s, skip_partials, nvar, do_du_dt, do_dv_dt, do_equL, ierr)

         use hydro_vars, only: set_Teff_info_for_eqns
         use chem_def
         use atm_def
         use eos_lib, only: Radiation_Pressure

         type (star_info), pointer :: s
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         logical, intent(in) :: do_du_dt, do_dv_dt, do_equL
         integer, intent(out) :: ierr

         integer :: i_lnd, i_lnT, i_lnR, i_lum, &
            i_du_dt, i_u, i_v, i_dv_dt, i_equL, k
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
            dlnP_bc_dlnPgas_const_T, dlnP_bc_dlnT_const_Pgas, dlnP_dlnPgas_const_T, &
            dlnP_dlnT_const_Pgas, dlnT_bc_dlnPgas_const_T, dlnT_bc_dlnT_const_Pgas, &
            dlnT_bc_dlnE_const_Rho, dlnT_dlnE_const_Rho, dlnP_dlnE_c_Rho, &
            dlnP_bc_dlnE_c_Rho, dlnT_bc_dlnd_c_E, dlnP_bc_dlnd_c_E
         logical :: test_partials

         include 'formats'

         !test_partials = (s% solver_iter == s% solver_test_partials_iter_number)
         test_partials = .false.
         
         ierr = 0

         i_dv_dt = s% i_dv_dt
         i_du_dt = s% i_du_dt
         i_equL = s% i_equL

         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_v = s% i_v
         i_u = s% i_u

         s% dlnP_dm_expected(1) = s% dlnP_dm_expected(2) ! not used except for output
         s% dlnT_dm_expected(1) = s% dlnT_dm_expected(2) ! not used except for output
         
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
            ! r == 0 means no photosphere
            dP0 = 0
            dT0 = 0
         else ! offset P and T from outer edge of cell 1 to center of cell 1
            dP0 = s% cgrav(1)*s% m_grav(1)*s% dm(1)/(8*pi*r*r*r*r)
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

            dP0_dlnR = -4*dP0
            dT0_dlnR = -4*dT0 + dP0*s% d_gradT_dlnR(1)*s% T(1)/s% P(1)

            dPinv_dlnT = -s% chiT_for_partials(1)/s% P(1)
            dT0_dlnT = &
                 dT0 + &
                 dP0*s% d_gradT_dlnT00(1)*s% T(1)/s% P(1) + &
                 dP0*s% gradT(1)*s% T(1)*dPinv_dlnT
            dPinv_dlnd = -s% chiRho_for_partials(1)/s% P(1)
            dT0_dlnd = &
                 dP0*s% d_gradT_dlnd00(1)*s% T(1)/s% P(1) + &
                 dP0*s% gradT(1)*s% T(1)*dPinv_dlnd

            dT0_dL = dP0*s% d_gradT_dL(1)*s% T(1)/s% P(1)

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

         dlnT_bc_dlnT = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnT &
               + dlnT_bc_dT0*dT0_dlnT
         dlnT_bc_dlnd = dlnT_bc_dlnTsurf*dlnTsurf_dlnkap*dlnkap_dlnd &
               + dlnT_bc_dT0*dT0_dlnd
         dlnT_bc_dL = dlnT_bc_dlnTsurf*dlnTsurf_dL + dlnT_bc_dT0*dT0_dL
         dlnT_bc_dlnR = dlnT_bc_dlnTsurf*dlnTsurf_dlnR + dlnT_bc_dT0*dT0_dlnR
         
         if (s% use_compression_outer_BC) then
            call set_compression_BC(ierr)
         else if (s% use_zero_Pgas_outer_BC) then
            P_rad = Radiation_Pressure(s% T_start(1))
            call set_momentum_BC(P_rad, 0d0, 0d0, 0d0, 0d0, ierr)
         else if (s% use_fixed_Psurf_outer_BC) then
            call set_momentum_BC(s% fixed_Psurf, 0d0, 0d0, 0d0, 0d0, ierr)
         else if (s% use_momentum_outer_BC) then
            call set_momentum_BC( &
               P_surf, dlnPsurf_dL, dlnPsurf_dlnR, &
               dlnPsurf_dlnkap*dlnkap_dlnd, dlnPsurf_dlnkap*dlnkap_dlnT, ierr)
         else if (s% use_fixed_vsurf_outer_BC) then
            call set_fixed_vsurf_outer_BC(ierr)
         else
            call set_Psurf_BC(ierr)
         end if
         if (ierr /= 0) return
         
         if (s% TDC_flag .or. .not. do_equL) return ! no Tsurf BC
         
         if (s% use_zero_dLdm_outer_BC) then
            call set_zero_dL_dm_BC(ierr)
         else if (s% use_T_black_body_outer_BC) then
            call set_T_black_body_BC(ierr)
         else
            call set_Tsurf_BC(ierr)
         end if

         if (test_partials) then
            s% solver_test_partials_val = lnP_bc
         end if

         if (skip_partials) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = dlnP_bc_dlnd
            write(*,*) 'PT_eqns_surf', s% solver_test_partials_var
         end if

         contains


         subroutine set_fixed_vsurf_outer_BC(ierr)
            integer, intent(out) :: ierr
            integer :: i_eqn
            include 'formats'
            ierr = 0
            if ((.not. do_du_dt) .and. (.not. do_dv_dt)) then
               ierr = -1
               write(*,*) 'set_fixed_vsurf_outer_BC requires u_flag or v_flag true'
               return
            end if
            if (do_du_dt) then
               s% equ(i_du_dt, 1) = (s% u(1) - s% fixed_vsurf)/s% csound_start(1)
               if (is_bad(s% equ(i_du_dt, 1))) then
                  write(*,1) 'equ(i_du_dt, 1)', s% equ(i_du_dt, 1)
                  stop 'set_fixed_vsurf_outer_BC'
               end if
               s% u_residual(1) = s% equ(i_du_dt, 1)
               if (.not. skip_partials) &
                  call e00(s, i_du_dt, i_u, 1, nvar, 1d0/s% csound_start(1))
            end if
            if (do_dv_dt) then
               s% equ(i_dv_dt, 1) = (s% v(1) - s% fixed_vsurf)/s% csound_start(1)
               if (is_bad(s% equ(i_dv_dt, 1))) then
                  write(*,1) 'equ(i_dv_dt, 1)', s% equ(i_dv_dt, 1)
                  stop 'set_fixed_vsurf_outer_BC'
               end if
               s% v_residual(1) = s% equ(i_dv_dt, 1)
               if (.not. skip_partials) &
                  call e00(s, i_dv_dt, i_v, 1, nvar, 1d0/s% csound_start(1))
            end if 
         end subroutine set_fixed_vsurf_outer_BC

         subroutine set_zero_dL_dm_BC(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            if (s% L(1) <= 0d0) then
               s% equ(i_equL,1) = s% L(1) - s% L(2)
               s% equL_residual(1) = s% equ(i_equL,1)
               if (.not. skip_partials) then
                  call e00(s,i_equL,i_lum,1,nvar,1d0)
                  call ep1(s,i_equL,i_lum,1,nvar,-1d0)
               end if
               return
            end if
            s% equ(i_equL,1) = s% L(2)/s% L(1) - 1d0
            s% equL_residual(1) = s% equ(i_equL,1)
            if (skip_partials) return
            call e00(s,i_equL,i_lum,1,nvar,-s% L(2)/(s% L(1)*s% L(1)))
            call ep1(s,i_equL,i_lum,1,nvar,1d0/s% L(1))
         end subroutine set_zero_dL_dm_BC


         subroutine set_T_black_body_BC(ierr)
            integer, intent(out) :: ierr
            call set_BB_BC(4d0,ierr)
         end subroutine set_T_black_body_BC


         subroutine set_BB_BC(factor,ierr)
            real(dp), intent(in) :: factor
            integer, intent(out) :: ierr

            real(dp) :: rmid, Lmid, Tscale, T1, dT1_dL, &
               dT1_dL00, dT1_dLp1, dT1_dlnR, dT1_dlnR00, dT1_dlnRp1
            logical :: test_partials

            include 'formats'
            ierr = 0
            
            !test_partials = (s% solver_test_partials_k == 1) 
            test_partials = .false.

            if (.not. s% do_struct_thermo) then ! dummy eqn
               s% equ(i_equL,1) = 0
               s% equL_residual(1) = s% equ(i_equL,1)
               if (skip_partials) return
               call e00(s,i_equL,i_lnT,1,nvar,1d0)
               return
            end if

            rmid = s% rmid(1)
            if (s% use_fixed_L_for_BB_outer_BC) then
               Lmid = s% fixed_L_for_BB_outer_BC
               T1 = s% Tsurf_factor* &
                  pow(Lmid/(factor*pi*rmid*rmid*boltz_sigma), 0.25d0)
               dT1_dL = 0
            else
               if (s% tau_for_L_BB > 0d0) then
                  Lmid = s% L_for_BB_outer_BC
               else
                  Lmid = (s% L(1) + s% L(2))/2
               end if
               T1 = s% Tsurf_factor* &
                  pow(Lmid/(factor*pi*rmid*rmid*boltz_sigma), 0.25d0)
               dT1_dL = T1/(4*Lmid)
            end if
            if (Lmid <= 0) then
               if (s% report_ierr) then
                  write(*,2) 'Lmid <= 0 for set_BB_BC', s% model_number
                  write(*,*) 's% use_fixed_L_for_BB_outer_BC', s% use_fixed_L_for_BB_outer_BC
                  write(*,2) 's% tau_for_L_BB', s% model_number, s% tau_for_L_BB
                  write(*,2) 's% L_for_BB_outer_BC', s% model_number, s% L_for_BB_outer_BC
                  write(*,2) 's% L(1)', s% model_number, s% L(1)
                  write(*,2) 'L(2)', s% model_number, s% L(2)
                  stop 'set_BB_BC'
               end if
               ierr = -1
               return
            end if

            dT1_dL00 = dT1_dL/2
            dT1_dLp1 = dT1_dL/2
            dT1_dlnR = -T1/2
            dT1_dlnR00 = dT1_dlnR*s% drmid_dlnR00(1)/rmid
            dT1_dlnRp1 = dT1_dlnR*s% drmid_dlnRp1(1)/rmid

            Tscale = 1d6*s% T_start(1) ! 1d6 to reduce equ size compared to k > 1 cases
            s% equ(i_equL, 1) = (T1 - s% T(1))/Tscale
            s% equL_residual(1) = s% equ(i_equL,1)
            if (test_partials) then
               s% solver_test_partials_val = s% equ(i_equL, 1)
            end if

            if (is_bad(s% equ(i_equL, 1))) then
               write(*,1) 'equ(i_equL, 1)', s% equ(i_equL, 1)
               write(*,1) 's% r(1)', s% r(1)
               write(*,1) 'T1', T1
               write(*,1) 'Lmid', Lmid
               write(*,1) 'rmid', rmid
               write(*,1) 's% L(1)', s% L(1)
               write(*,1) 's% L(2)', s% L(2)
               write(*,1) 's% r(1)', s% r(1)
               write(*,1) 's% r(2)', s% r(2)
               write(*,1) 's% T(1)', s% T(1)
               write(*,1) 's% T_start(1)', s% T_start(1)
               write(*,1) 'Tscale', Tscale
               stop 'set_BB_BC'
            end if

            if (skip_partials) return

            call e00(s, i_equL, i_lnT, 1, nvar, -s% T(1)/Tscale)
            
            if (s% use_fixed_L_for_BB_outer_BC .or. s% tau_for_L_BB > 0d0) then
            else if (i_lum /= 0) then
               call e00(s, i_equL, i_lum, 1, nvar, dT1_dL00/Tscale)
               call ep1(s, i_equL, i_lum, 1, nvar, dT1_dLp1/Tscale)
            end if

            if (s% do_struct_hydro) then
               call e00(s, i_equL, i_lnR, 1, nvar, dT1_dlnR00/Tscale)
               call ep1(s, i_equL, i_lnR, 1, nvar, dT1_dlnRp1/Tscale)
            end if
            
            if (test_partials) then
               s% solver_test_partials_var = i_lnT
               s% solver_test_partials_dval_dx = -s% T(1)/Tscale
               write(*,*) 'set_BB_BC', s% solver_test_partials_var
            end if

         end subroutine set_BB_BC


         subroutine set_Tsurf_BC(ierr)
            integer, intent(out) :: ierr

            logical :: test_partials
            integer :: i_T_BC
            real(dp) :: scale
            include 'formats'

            !test_partials = (1 == s% solver_test_partials_k)
            test_partials = .false.

            ierr = 0
            
            if (s% TDC_flag) return
            i_T_BC = i_equL            
            if (i_T_BC == 0) then
               write(*,2) 'i_T_BC', s% model_number, i_T_BC
               stop 'set_Tsurf_BC'
            end if

            if (.not. s% do_struct_thermo) then ! dummy eqn
               s% equ(i_T_BC,1) = 0
               s% equL_residual(1) = s% equ(i_T_BC,1)
               if (skip_partials) return
               call e00(s,i_T_BC,i_lnT,1,nvar,1d0)
               return
            end if
            
            scale = 1d0 ! 1d-3
            
            s% equ(i_T_BC, 1) = (lnT_bc - s% lnT(1))*scale
            s% equL_residual(1) = s% equ(i_T_BC,1)
            
            if (is_bad(s% equ(i_T_BC,1))) then
               write(*,1) 'lnT_bc', lnT_bc
               write(*,1) 's% lnT(1)', s% lnT(1)
               write(*,1) 's% L(1)', s% L(1)
               stop 'set_Tsurf_BC'
            end if

            if (test_partials) then
               s% solver_test_partials_val = lnT_bc
               write(*,1) 'lnT_bc', lnT_bc
               write(*,1) 'logT', s% lnT(1)/ln10
               write(*,1) 'logPgas', s% lnPgas(1)/ln10
            end if

            if (skip_partials) return

            call e00(s, i_T_BC, i_lnT, 1, nvar, (dlnT_bc_dlnT - 1d0)*scale)

            if (i_lum /= 0) call e00(s, i_T_BC, i_lum, 1, nvar, dlnT_bc_dL*scale)

            if (.not. s% do_struct_hydro) return

            ! partial of temperature eqn wrt (lnR or v) and (lnd or lnPgas)

            call e00(s, i_T_BC, i_lnR, 1, nvar, dlnT_bc_dlnR*scale)

            call e00(s, i_T_BC, i_lnd, 1, nvar, dlnT_bc_dlnd*scale)
            if (test_partials) then
               s% solver_test_partials_var = i_lnT
               s% solver_test_partials_dval_dx = dlnT_bc_dlnT_const_Pgas
               write(*,*) 'set_Tsurf_BC', s% solver_test_partials_var
            end if

         end subroutine set_Tsurf_BC

 
         subroutine set_Psurf_BC(ierr) ! use Psurf from atm
            integer, intent(out) :: ierr
            logical :: test_partials

            integer :: i_eqn

            include 'formats'
            ierr = 0
            
            !test_partials = (s% solver_iter == s% solver_test_partials_iter_number)
            test_partials = .false.

            if (s% u_flag) then
               i_eqn = i_du_dt
            else
               i_eqn = i_dv_dt
            end if

            s% equ(i_eqn,1) = lnP_bc - s% lnP(1)
            
            if (test_partials) then
               s% solver_test_partials_val = lnP_bc ! s% equ(i_eqn,1)
            end if

            if (skip_partials) return

            call e00(s, i_eqn, i_lnR, 1, nvar, dlnP_bc_dlnR)

            call e00(s, i_eqn, i_lnd, 1, nvar, &
               dlnP_bc_dlnd - s% chiRho_for_partials(1))

            if (.not. s% do_struct_thermo) return

            call e00(s, i_eqn, i_lnT, 1, nvar,&
                dlnP_bc_dlnT - s% chiT_for_partials(1))

            dlnP_bc_dL = dlnP_bc_dlnPsurf*dlnPsurf_dL
            if (i_lum /= 0) call e00(s, i_eqn, i_lum, 1, nvar, dlnP_bc_dL)
            
            if (test_partials) then
               s% solver_test_partials_var = s% i_lnT
               s% solver_test_partials_dval_dx = dlnP_bc_dlnT ! - s% chiT_for_partials(1)
               write(*,*) 'set_Psurf_BC', s% solver_test_partials_var
            end if

         end subroutine set_Psurf_BC
         

         subroutine set_momentum_BC( & ! momentum eqn using P(1) and P
               P, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, ierr)
            use hydro_riemann, only: do_surf_Riemann_dudt_eqn

            use hydro_momentum, only: do_surf_momentum_eqn
            use auto_diff_support
            real(dp), intent(in) :: P, &
               dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
            integer, intent(out) :: ierr

            type(auto_diff_real_star_order1) :: P_surf_ad
            logical :: test_partials
            include 'formats'
            ierr = 0
            
            !test_partials = (s% solver_iter == s% solver_test_partials_iter_number)
            test_partials = .false.

            if (s% u_flag) then
               call do_surf_Riemann_dudt_eqn( &
                  s, P, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
                  skip_partials, nvar, ierr)
            else
               call wrap(P_surf_ad, P, &
                  0d0, P*dlnPsurf_dlnd, 0d0, &
                  0d0, P*dlnPsurf_dlnT, 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, P*dlnPsurf_dlnR, 0d0, &
                  0d0, 0d0, 0d0, &
                  0d0, P*dlnPsurf_dL, 0d0)
               call do_surf_momentum_eqn( &
                  s, P_surf_ad, skip_partials, nvar, ierr)
            end if
            
            if (test_partials) then
               s% solver_test_partials_val = P
            end if
            
            if (skip_partials) return
            if (test_partials) then
               s% solver_test_partials_var = s% i_lum
               s% solver_test_partials_dval_dx = P*dlnPsurf_dL
               write(*,*) 'set_momentum_BC', s% solver_test_partials_var
            end if

         end subroutine set_momentum_BC


         subroutine set_compression_BC(ierr)
            integer, intent(out) :: ierr

            include 'formats'
            ierr = 0

            if (i_lnd == 0) then
               write(*,*) 'set_compression_BC: no support for lnPgas'
               ierr = -1
               !return
               stop 'set_compression_BC'
            end if
            
            if (.not. (s% u_flag .or. s% v_flag)) then
               write(*,*) 'set_compression_BC: must have velocities on'
               ierr = -1
               !return
               stop 'set_compression_BC'
            end if

            if (s% u_flag) call set_compression_eqn(i_du_dt, s% u_residual)
            
            if (s% v_flag) call set_compression_eqn(i_dv_dt, s% v_residual)

         end subroutine set_compression_BC
         
         
         subroutine set_compression_eqn(i_eqn, residual)
            integer, intent(in) :: i_eqn
            real(dp), pointer :: residual(:)

            real(dp) :: d
            include 'formats'

            ! gradient of compression vanishes fixes density for cell 1
               ! d_dt(1/rho(1)) = d_dt(1/rho(2))
               ! widely used.  e.g., Grott, Chernigovski, Glatzel, 2005.
            
            s% equ(i_eqn, 1) = s% rho(2)*s% dlnd_dt(1) - s% rho(1)*s% dlnd_dt(2)
            residual(1) = s% equ(i_eqn, 1)
            
            if (is_bad(s% equ(i_eqn, 1))) then
               write(*,1) 'equ(i_eqn, 1)', s% equ(i_eqn, 1)
               write(*,1) 's% rho(1)', s% rho(1)
               write(*,1) 's% rho_start(1)', s% rho_start(1)
               write(*,1) 's% rho(2)', s% rho(2)
               write(*,1) 's% rho_start(2)', s% rho_start(2)
               write(*,*)
               stop 'set_compression_BC'
            end if

            if (skip_partials) return
            
            d = s% rho(2)*s% dVARDOT_dVAR - s% rho(1)*s% dlnd_dt(2)
            call e00(s, i_eqn, i_lnd, 1, nvar, d)
            d = s% rho(2)*s% dlnd_dt(1) - s% rho(1)*s% dVARDOT_dVAR
            call ep1(s, i_eqn, i_lnd, 1, nvar, d)
         
         end subroutine set_compression_eqn


      end subroutine PT_eqns_surf


      subroutine eval_gradT_info( &
            s, k, gradT, d_gradT_dL, d_gradT_dw, d_gradT_dlnR, d_gradT_du, &
            d_gradT_dlnd00, d_gradT_dlndm1, &
            d_gradT_dlnT00, d_gradT_dlnTm1, &
            d_gradT_dlnT00_const_Pgas, d_gradT_dlnTm1_const_Pgas, &
            d_gradT_dlnPgas00_const_T, d_gradT_dlnPgasm1_const_T, &
            d_gradT_dln_cvpv0, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: &
            gradT, d_gradT_dL, d_gradT_dw, d_gradT_dlnR, d_gradT_du, &
            d_gradT_dlnd00, d_gradT_dlndm1, &
            d_gradT_dlnT00, d_gradT_dlnTm1, &
            d_gradT_dlnT00_const_Pgas, d_gradT_dlnTm1_const_Pgas, &
            d_gradT_dlnPgas00_const_T, d_gradT_dlnPgasm1_const_T, &
            d_gradT_dln_cvpv0
         integer, intent(out) :: ierr
         logical :: okay

         include 'formats'

         ierr = 0
         gradT = s% gradT(k)

         d_gradT_dL = s% d_gradT_dL(k)
         d_gradT_dw = s% d_gradT_dw_div_wc(k)
         d_gradT_dlnR = s% d_gradT_dlnR(k)
         d_gradT_du = 0d0
         if (s% conv_vel_flag) then
            d_gradT_dln_cvpv0 = s% d_gradT_dln_cvpv0(k)
         else
            d_gradT_dln_cvpv0 = 0d0
         end if

         d_gradT_dlnd00 = s% d_gradT_dlnd00(k)
         d_gradT_dlndm1 = s% d_gradT_dlndm1(k)
         d_gradT_dlnT00 = s% d_gradT_dlnT00(k)
         d_gradT_dlnTm1 = s% d_gradT_dlnTm1(k)

      end subroutine eval_gradT_info


      ! only used for dlnT_dm equation
      subroutine eval_dlnPdm_qhse(s, k, & ! calculate the expected dlnPdm for HSE
            m, dlnPdm_qhse, d_dlnPdm_dlnR, d_dlnPdm_dL, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            d_dlnPdm_dw, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas, &
            ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: &
            m, dlnPdm_qhse, d_dlnPdm_dlnR, d_dlnPdm_dL, &
            d_dlnPdm_dlnd00, d_dlnPdm_dlnT00, &
            d_dlnPdm_dlndm1, d_dlnPdm_dlnTm1, &
            d_dlnPdm_dlnPgas00_const_T, d_dlnPdm_dlnT00_const_Pgas, &
            d_dlnPdm_dlnPgasm1_const_T, d_dlnPdm_dlnTm1_const_Pgas, &
            d_dlnPdm_dw, &
            Ppoint, &
            dPpoint_dlnd00, dPpoint_dlndm1, dPpoint_dlnT00, dPpoint_dlnTm1, &
            dPpoint_dlnPgas00_const_T, dPpoint_dlnPgasm1_const_T, &
            dPpoint_dlnT00_const_Pgas, dPpoint_dlnTm1_const_Pgas
         integer, intent(out) :: ierr

         real(dp) :: alfa, r4, d_r4_dlnR, dlnq_dq, P00, Pm1, Ppoint_start, &
            theta, rtheta, d_rtheta_dlnR, R2, d_R2_dlnR, extra_dlnPdm, area_Ppoint

         include 'formats'

         ierr = 0

         ! basic eqn is dP/dm = -G m / (4 pi r^4)
         ! divide by <P> to make it unitless
         ! simple average is adequate for <P> since is only for normalizing the equation.
         ! however, be careful to use same <P> for both sides of equation.....

         ! for rotation, multiply gravity by factor fp.  MESA 2, eqn 22.

         ! for tau < 2/3, multiply by Paczynski factor for dilution of radiation
            ! B. Paczynski, 1969, Acta Astr., vol. 19.  eqn 13

         ! dlnPdm_qhse = -G m / (4 pi r^4 <P>)

         if (s% using_velocity_time_centering) then
            theta = 0.5d0
            rtheta = s% r_start(k)
            d_rtheta_dlnR = 0
         else
            theta = 1d0
            rtheta = s% r(k)
            d_rtheta_dlnR = s% r(k)
         end if
         
         R2 = s% R2(k)
         d_R2_dlnR = s% d_R2_dlnR(k)
         
         r4 = R2*s% r(k)*rtheta
         d_r4_dlnR = (d_R2_dlnR + R2)*s% r(k)*rtheta + &
            R2*s% r(k)*d_rtheta_dlnR

         P00 = theta*s% P(k) + (1d0-theta)*s% P_start(k)
         if (k == 1) then
            alfa = 1
            Pm1 = 0d0
            Ppoint = alfa*P00
            dPpoint_dlndm1 = 0
            dPpoint_dlnTm1 = 0
            dPpoint_dlnPgasm1_const_T = 0
            dPpoint_dlnTm1_const_Pgas = 0
         else
            Pm1 = theta*s% P(k-1) + (1d0-theta)*s% P_start(k-1)
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            Ppoint = alfa*P00 + (1-alfa)*Pm1
            dPpoint_dlndm1 = (1-alfa)*theta*s% P(k-1)*s% chiRho_for_partials(k-1)
            dPpoint_dlnTm1 = (1-alfa)*theta*s% P(k-1)*s% chiT_for_partials(k-1)
         end if

         dPpoint_dlnd00 = alfa*theta*s% P(k)*s% chiRho_for_partials(k)
         dPpoint_dlnT00 = alfa*theta*s% P(k)*s% chiT_for_partials(k)

         m = s% m_grav(k)
         dlnPdm_qhse = -s% cgrav(k)*m/(pi4*r4*Ppoint)

         area_Ppoint = pi4*R2*Ppoint
         if (s% use_other_momentum .or. s% use_other_momentum_implicit) then
            extra_dlnPdm = s% extra_grav(k)/area_Ppoint
            dlnPdm_qhse = dlnPdm_qhse + extra_dlnPdm
         end if

         if (s% rotation_flag .and. s% use_gravity_rotation_correction) then
            if (s% w_div_wc_flag) then
               d_dlnPdm_dw = dlnPdm_qhse*s% dfp_rot_dw_div_wc(k)
            end if
            dlnPdm_qhse = dlnPdm_qhse*s% fp_rot(k)
         end if

         d_dlnPdm_dlnR = -(d_r4_dlnR/r4)*dlnPdm_qhse
         dlnq_dq = 1/s% q(k)

         d_dlnPdm_dlnd00 = -dPpoint_dlnd00/Ppoint*dlnPdm_qhse
         d_dlnPdm_dlnT00 = -dPpoint_dlnT00/Ppoint*dlnPdm_qhse
         d_dlnPdm_dlndm1 = -dPpoint_dlndm1/Ppoint*dlnPdm_qhse
         d_dlnPdm_dlnTm1 = -dPpoint_dlnTm1/Ppoint*dlnPdm_qhse

         d_dlnPdm_dL = 0d0
         
         if (s% use_other_momentum_implicit) then
            d_dlnPdm_dlndm1 = d_dlnPdm_dlndm1 + s% d_extra_grav_dlndm1(k)/area_Ppoint &
               - (dPpoint_dlndm1/Ppoint)*extra_dlnPdm
            d_dlnPdm_dlnd00 = d_dlnPdm_dlnd00 + s% d_extra_grav_dlnd00(k)/area_Ppoint &
               - (dPpoint_dlnd00/Ppoint)*extra_dlnPdm
            d_dlnPdm_dlnTm1 = d_dlnPdm_dlnTm1 + s% d_extra_grav_dlnTm1(k)/area_Ppoint &
               - (dPpoint_dlnTm1/Ppoint)*extra_dlnPdm
            d_dlnPdm_dlnT00 = d_dlnPdm_dlnT00 + s% d_extra_grav_dlnT00(k)/area_Ppoint &
               - (dPpoint_dlnT00/Ppoint)*extra_dlnPdm
            d_dlnPdm_dlnR = d_dlnPdm_dlnR + s% d_extra_grav_dlnR(k)/area_Ppoint
            d_dlnPdm_dL = d_dlnPdm_dL + s% d_extra_grav_dL(k)/area_Ppoint
         end if

         if (is_bad(dlnPdm_qhse)) then
            ierr = -1
            s% retry_message = 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
            if (s% report_ierr) then
!$OMP critical (hydro_vars_crit1)
               write(*,*) 'eval_dlnPdm_qhse: is_bad(dlnPdm_qhse)'
               write(*,2) 'dlnPdm_qhse', k, dlnPdm_qhse
               write(*,2) 's% tau(k)', k, s% tau(k)
               write(*,2) 's% fp_rot(k)', k, s% fp_rot(k)
               write(*,2) 's% w_div_w_crit_roche(k)', k, s% w_div_w_crit_roche(k)
               write(*,2) 'r4', k, r4
               write(*,2) 'Ppoint', k, Ppoint
               write(*,2) 'm', k, m
               write(*,2) 's% cgrav(k)', k, s% cgrav(k)
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

