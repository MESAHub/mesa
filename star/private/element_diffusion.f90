! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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


      module element_diffusion

      use star_private_def
      use const_def

      implicit none

      private
      public :: do_element_diffusion, finish_element_diffusion

      logical, parameter :: dbg = .false.


      contains


      subroutine do_element_diffusion(s, dt_in, ierr)
         ! return ierr /= 0 if cannot satisfy accuracy requirements
         use chem_def, only: chem_isos, ihe4, ih1
         use chem_lib, only: chem_get_iso_id
         use star_utils, only: start_time, update_time
         use diffusion, only: &
            do_solve_diffusion, set_diffusion_classes, diffusion_min_nc
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt_in
         integer, intent(out) :: ierr

         integer :: i, j, k, kk, nc, m, nzlo, nzhi, nz, species, iounit, &
            steps_used, total_num_iters, total_num_retries, cid, he4
         integer(8) :: time0
         real(dp) :: s1, s2, dqsum, dist, r, Hp, dt, total, &
            gradT_mid, gradRho_mid, alfa, gradRho_face, chiRho_face, chiT_face
         real(dp) :: rho, Pgas, T, &
            logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
            e, de, e_with_xa, Amass, Zcharge, min_D_mix

         integer, dimension(:), allocatable :: &
            class, class_chem_id, mixing_type, mixing_type_arg
         real(dp), dimension(:), allocatable :: &
            gamma, free_e, &
            dlnPdm_face, dlnT_dm_face, dlnRho_dm_face, &
            dlnPdm, dlnT_dm, dlnPdm_mid, dlnT_dm_mid, dlnRho_dm_mid, dm_hat
         real(dp), dimension(:,:), allocatable :: &
            X_init, X_final, typical_charge, &
            D_self, v_advection, v_total, &
            vlnP, vlnT, v_rad, g_rad, GT, xa_save
         real(dp), dimension(:,:,:), allocatable :: CD
         character (len=8), allocatable :: class_name(:)


         logical :: dumping, okay

         include 'formats'

         ierr = 0
         dt = dt_in
         nz = s% nz

         s% num_diffusion_solver_steps = 0
         
         s% eps_diffusion(1:nz) = 0d0
         do k = 1, nz
            s% energy_start(k) = s% energy(k)
            ! This is just to track energy changes due to diffusion.
            ! s% energy_start gets reset later in do_struct_burn_mix before structure solve.
         end do

         if ((.not. s% do_element_diffusion) .or. dt < s% diffusion_dt_limit) then
            s% num_diffusion_solver_iters = 0 ! Flush diff iters to avoid crashing the timestep.
            s% edv(:,1:nz) = 0
            s% eps_WD_sedimentation(1:nz) = 0d0
            if (s% do_element_diffusion .and. s% report_ierr .and. dt < s% diffusion_dt_limit) &
               write(*,2) 'skip diffusion this step: dt < s% diffusion_dt_limit', &
                  s% model_number, dt, s% diffusion_dt_limit
            return
         end if
                         
         s% need_to_setvars = .true.

         if (s% use_other_diffusion_factor) then
            call s% other_diffusion_factor(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'do_element_diffusion failed in other_diffusion_factor'
               return
            end if
         end if

         if (s% use_other_diffusion) then
            call s% other_diffusion(s% id, dt_in, ierr)
            if (ierr /= 0 .and. s% report_ierr) &
               write(*,*) 'do_element_diffusion failed in other_diffusion_factor'
            return
         end if
         
         if (s% doing_timing) call start_time(s, time0, total)  

         nz = s% nz
         nzlo = 1
         nzhi = nz
         min_D_mix = s% D_mix_ignore_diffusion
         species = s% species
         if ( s% diffusion_use_full_net ) then
            if( species > 100 ) then
               stop "Net is too large for diffusion_use_full_net. max_num_diffusion_classes = 100"
            end if
            nc = species
         else
            nc = s% diffusion_num_classes
         end if
         m = nc+1

         allocate( &
            class(species), class_chem_id(nc), class_name(nc), CD(nc,nc,nz), mixing_type(nz), &
            gamma(nz), free_e(nz), dlnPdm_face(nz), dlnT_dm_face(nz), dlnRho_dm_face(nz), &
            dlnPdm(nz), dlnT_dm(nz), dlnPdm_mid(nz), dlnT_dm_mid(nz), dlnRho_dm_mid(nz), dm_hat(nz), &
            X_init(nc,nz), X_final(nc,nz), typical_charge(nc,nz), &
            D_self(nc,nz), v_advection(nc,nz), v_total(nc,nz), xa_save(species,nz), &
            vlnP(nc,nz), vlnT(nc,nz), v_rad(nc,nz), g_rad(nc,nz), GT(nc,nz))

         call set_extras(ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,*) 'do_element_diffusion failed in set_extras'
            return
         end if


         !write(*,1) 's% diffusion_min_T_at_surface', s% diffusion_min_T_at_surface
         !reset nzlo if necessary; nzlo=1 above
         nzlo = 1
         if (s% diffusion_min_T_at_surface > 0) then
            k = nzlo
            do while (s% T(k) < s% diffusion_min_T_at_surface)
               k = k+1
               if (k > nz) exit
               nzlo = k
            end do
         end if

         !write(*,2) 'diffusion_min_T_at_surface', nzlo

         if (s% diffusion_min_dq_ratio_at_surface > 0) then
            dqsum = sum(s% dq(1:nzlo))
            k = nzlo
            do while (dqsum < s% diffusion_min_dq_ratio_at_surface*s% dq(k+1))
               k = k+1
               if (k >= nz) exit
               dqsum = dqsum + s% dq(k)
               nzlo = k
            end do
         end if

         !write(*,2) 'before diffusion_min_dq_ratio_at_surface', nzlo, s% diffusion_min_dq_at_surface

         if (s% diffusion_min_dq_at_surface > 0) then
            dqsum = sum(s% dq(1:nzlo))
            k = nzlo
            do while (dqsum < s% diffusion_min_dq_at_surface)
               k = k+1
               if (k > nz) exit
               ! don't go across composition transition
               if (maxloc(s% xa(:,k),dim=1) /= maxloc(s% xa(:,k-1),dim=1)) exit
               dqsum = dqsum + s% dq(k)
               nzlo = k
            end do
            !write(*,2) 'after diffusion_min_dq_ratio_at_surface', nzlo, dqsum
         end if

         !write(*,3) 'nzlo mixing_type', nzlo, s% mixing_type(nzlo)

         if (s% D_mix(nzlo) > min_D_mix) then
            kk = nzlo
            do k = nzlo, nz-1
               if (maxloc(s% xa(:,k),dim=1) /= maxloc(s% xa(:,k-1),dim=1) .or. &
                    s% D_mix(k+1) <= min_D_mix) then
                  nzlo=k; exit
               endif
            end do
            nzlo = kk + (nzlo - kk)*4/5 ! back up into the convection zone
         end if

         !write(*,3) 'nzlo mixing_type', nzlo, s% mixing_type(nzlo)
         !write(*,3) 'nzlo-1 mixing_type', nzlo-1, s% mixing_type(nzlo-1)
         !write(*,*)
         !write(*,*)

         !reset nzhi if necessary; nzhi=nz above
         if (s% D_mix(nzhi) > min_D_mix) then
            do k = nz, 2, -1
               if (s% D_mix(k-1) <= min_D_mix) then
                  nzhi=k; exit
               end if
            end do
            nzhi = (3*nzhi+nz)/4 ! back up some into the convection zone
         end if

         if(s% diffusion_use_full_net) then
            do j=1,nc
               class_chem_id(j) = s% chem_id(j) ! Just a 1-1 map between classes and chem_ids.
            end do
         else
            do j=1,nc
               cid = chem_get_iso_id(s% diffusion_class_representative(j))
               if (cid <= 0) then
                  write(*,'(a,3x,i3)') 'bad entry for diffusion_class_representative: ' // &
                       trim(s% diffusion_class_representative(j)), j
                  return
               end if
               class_chem_id(j) = cid
            end do
         end if

         call set_diffusion_classes( &
            nc, species, s% chem_id, class_chem_id, s% diffusion_class_A_max, s% diffusion_use_full_net, &
            class, class_name)

         s% diffusion_call_number = s% diffusion_call_number + 1
         dumping = (s% diffusion_call_number == s% diffusion_dump_call_number)

         if ( s% diffusion_use_full_net ) then
            s% diffusion_calculates_ionization = .true. ! class_typical_charges can't be used, so make sure they aren't.
         end if
         
         if (.not. s% diffusion_calculates_ionization) then
            do j=1,nc
               typical_charge(j,1:nz) = s% diffusion_class_typical_charge(j)
            end do
         end if

         do k=1,nz
            do j=1,species
               xa_save(j,k) = s% xa(j,k)
            end do
         end do

         mixing_type(1:nz) = no_mixing

         do k=1,nz-1
            s1 = dlnPdm(k)
            s2 = dlnPdm(k+1)
            if (s1*s2 <= 0) then
               dlnPdm_mid(k) = 0
            else
               dlnPdm_mid(k) = 2*s1*s2/(s1+s2)
            end if
            gradT_mid = 0.5d0*(s% gradT(k) + s% gradT(k+1))
            dlnT_dm_mid(k) = gradT_mid*dlnPdm_mid(k)
            gradRho_mid = (1 - s% chiT(k)*gradT_mid)/s% chiRho(k)
            dlnRho_dm_mid(k) = gradRho_mid*dlnPdm_mid(k)
         end do
         dlnPdm_mid(nz) = dlnPdm(nz)
         dlnT_dm_mid(nz) = dlnT_dm(nz)
         gradRho_mid = (1 - s% chiT(nz)*s% gradT(nz))/s% chiRho(nz)
         dlnRho_dm_mid(nz) = gradRho_mid*dlnPdm_mid(nz)

         do k=2,nz
            dlnPdm_face(k) = dlnPdm(k)
            dlnT_dm_face(k) = dlnT_dm(k)
            alfa = s% dm(k-1)/(s% dm(k) + s% dm(k-1))
            chiT_face = alfa*s% chiT(k) + (1-alfa)*s% chiT(k-1)
            chiRho_face = alfa*s% chiRho(k) + (1-alfa)*s% chiRho(k-1)
            gradRho_face = (1 - chiT_face*s% gradT(k))/chiRho_face
            dlnRho_dm_face(k) = gradRho_face*dlnPdm(k)
         end do
         dlnPdm_face(1) = 0
         dlnT_dm_face(1) = 0
         dlnRho_dm_face(1) = 0

         if (dumping) call dump_diffusion_info

         ! print *, "About to call do_solve_diffusion. Here's the class structure"
         ! print *, "chem_id:        ", s% chem_id
         ! print *, "class_chem_id:  ", class_chem_id
         ! print *, "class:          ", class
         ! print *, "class name:     ", class_name
         
         ! args are at cell center points.
         !if (s% show_diffusion_info) write(*,*) 'call solve_diffusion'
         !write(*,4) 'call do_solve_diffusion nzlo nzhi nz', nzlo, nzhi, nz, &
         !   sum(s% xa(1,1:nzlo))
         call do_solve_diffusion( &
            s, nz, species, nc, m, class, class_chem_id, s% net_iso, s% chem_id, &
            s% abar, s% ye, free_e, s% dm_bar, s% dm, &
            s% T, s% lnT, s% rho, s% lnd, s% rmid, &
            dlnPdm_mid, dlnT_dm_mid, dlnRho_dm_mid, &
            s% L, s% r, dlnPdm_face, dlnT_dm_face, dlnRho_dm_face, &
            s% diffusion_use_iben_macdonald, dt, s% diffusion_dt_div_timescale, &
            s% diffusion_steps_hard_limit, s% diffusion_iters_hard_limit, &
            s% diffusion_max_iters_per_substep, &
            s% diffusion_calculates_ionization, typical_charge, &
            s% diffusion_nsmooth_typical_charge, &
            s% diffusion_min_T_for_radaccel, s% diffusion_max_T_for_radaccel, &
            s% diffusion_min_Z_for_radaccel, s% diffusion_max_Z_for_radaccel, &
            s% diffusion_screening_for_radaccel, &
            s% op_mono_data_path, s% op_mono_data_cache_filename, &
            s% diffusion_v_max, s% R_center, &
            gamma, s% diffusion_gamma_full_on, s% diffusion_gamma_full_off, &
            s% diffusion_T_full_on, s% diffusion_T_full_off, &
            s% diffusion_class_factor, s% xa, &
            steps_used, total_num_iters, total_num_retries, nzlo, nzhi, X_init, X_final, &
            D_self, v_advection, v_total, vlnP, vlnT, v_rad, g_rad, &
            s% E_field, s% g_field_element_diffusion, ierr )
         s% num_diffusion_solver_steps = steps_used
         s% num_diffusion_solver_iters = total_num_iters

         if (dbg .or. s% show_diffusion_info .or. ierr /= 0) then
            if (ierr == 0) then
               write(*,'(a,f6.3,3x,a,1pe10.3,3x,99(a,i5,3x))') &
                  'log_dt', log10(s% dt/secyer), 'age', s% star_age, 'model', s% model_number, &
                  'iters', total_num_iters, 'steps', steps_used, 'retries', total_num_retries, &
                  'nzlo', nzlo, 'nzhi', nzhi, 'n', nzhi-nzlo+1, 'nz', nz, &
                  'diffusion_call_number', s% diffusion_call_number
            else
               write(*,'(a,2x,f10.3,3x,99(a,i5,3x))') &
                  'do_solve_diffusion FAILED: log_dt', log10(s% dt/secyer), 'model', s% model_number, &
                  'iters', total_num_iters, 'steps', steps_used, 'retries', total_num_retries, &
                  'nzlo', nzlo, 'nzhi', nzhi, 'n', nzhi-nzlo+1, 'nz', nz, &
                  'diffusion_call_number', s% diffusion_call_number
            end if
         end if

         if (ierr /= 0) then
            do k=1,nz
               do j=1,species
                  s% xa(j,k) = xa_save(j,k)
               end do
            end do
            if (s% report_ierr) then
               write(*, *)
               write(*, *) 'solve_diffusion returned false'
               write(*, *) 's% model_number', s% model_number
               write(*, *) 's% diffusion_call_number', s% diffusion_call_number
               write(*, *)
            end if
         end if

         if (dumping) stop 'debug: dump_diffusion_info'
         
         do k=nzlo+1,nzhi
            do j=1,species
               i = class(j)
               s% diffusion_D_self(j,k) = D_self(i,k)
               s% edv(j,k) = v_total(i,k)
               s% v_rad(j,k) = v_rad(i,k)
               s% g_rad(j,k) = g_rad(i,k)
               s% typical_charge(j,k) = typical_charge(i,k)
               s% diffusion_dX(j,k) = s% xa(j,k) - xa_save(j,k)
            end do
         end do

         if(s% do_diffusion_heating .and. s% do_WD_sedimentation_heating) then
            write(*,*) "do_diffusion_heating is incompatible with do_WD_sedimentation_heating"
            write(*,*) "at least one of these options must be set to .false."
            stop 'do_element_diffusion'
         end if         

         s% eps_WD_sedimentation(1:nz) = 0d0

         if(s% do_WD_sedimentation_heating) then
            do k=nzlo+1,nzhi
               ! loop over all ion species
               do i = 1,s% species
                  ! limit heating to species with significant mass fractions
                  if(s% xa(i,k) > s% min_xa_for_WD_sedimentation_heating) then
                     Amass = chem_isos% Z_plus_N(s% chem_id(i))
                     Zcharge = chem_isos% Z(s% chem_id(i))
                     s% eps_WD_sedimentation(k) = s% eps_WD_sedimentation(k) + &
                          s% eps_WD_sedimentation_factor * &
                          ( Amass - Zcharge * qe * s% E_field(k)/(amu * s% g_field_element_diffusion(k)) ) * &
                          sign(1d0,-1d0*s% edv(i,k)) * min( s% diffusion_v_max, abs(s% edv(i,k)) ) * s% xa(i,k) * &
                          s% g_field_element_diffusion(k) / Amass
                  end if
               end do
               ! For diagnostics:
               if( .false. .and. mod(k,100) == 0) then
                  write(*,*) "Zone: ", k
                  write(*,*) "eE/mg = ", (qe * s% E_field(k)/(amu * s% g_field_element_diffusion(k)))
                  write(*,*) "Net Force on Ne22 (units of mp*g) = ", &
                       (22d0 - 10d0 * qe * s% E_field(k)/(amu * s% g_field_element_diffusion(k)))
                  write(*,*) "g_field_element_diffusion/grav = ", s% g_field_element_diffusion(k) / s% grav(k)
               end if
            end do
         end if

         
         do k=1,nzlo
            do j=1,species
               s% diffusion_D_self(j,k) = s% diffusion_D_self(j,nzlo+1)
               s% edv(j,k) = 0d0 ! s% edv(j,nzlo+1)
               s% v_rad(j,k) = s% v_rad(j,nzlo+1)
               s% g_rad(j,k) = s% g_rad(j,nzlo+1)
               s% typical_charge(j,k) = s% typical_charge(j,nzlo+1)
               s% diffusion_dX(j,k) = s% xa(j,k) - xa_save(j,k)
            end do
            s% E_field(k) = 0d0 ! s% E_field(nzlo+1)
            s% g_field_element_diffusion(k) = 0d0 ! s% g_field_element_diffusion(nzlo+1)
         end do

         do k=nzhi+1,nz
            do j=1,species
               s% diffusion_D_self(j,k) = s% diffusion_D_self(j,nzhi)
               s% edv(j,k) = 0d0 ! s% edv(j,nzhi)
               s% v_rad(j,k) = s% v_rad(j,nzhi)
               s% g_rad(j,k) = s% g_rad(j,nzhi)
               s% typical_charge(j,k) = s% typical_charge(j,nzhi)
               s% diffusion_dX(j,k) = s% xa(j,k) - xa_save(j,k)
            end do
            s% E_field(k) = 0d0 ! s% E_field(nzhi)
            s% g_field_element_diffusion(k) = 0d0 ! s% g_field_element_diffusion(nzhi)
         end do

         if (s% doing_timing) call update_time(s, time0, total, s% time_element_diffusion)

         contains

         subroutine check_xa_sums(ierr)
            integer, intent(out) :: ierr
            integer :: k
            include 'formats'
            do k=1, nz
               if (abs(sum(s% xa(1:species, k)) - 1d0) > 1d-3) then
                  write(*,*) 'k', k
                  write(*,1) 'sum', sum(s% xa(1:species, k))
                  write(*,1) 'sum-1d0', sum(s% xa(1:species, k))-1d0
                  write(*,1) 'abs(sum-1d0)', abs(sum(s% xa(1:species, k))-1d0)
                  do j=1,species
                     write(*,2) 's% xa(j,k)', j, s% xa(j, k)
                  end do
                  ierr = -1
               end if
            end do
         end subroutine check_xa_sums

         subroutine dump_diffusion_info
            use utils_lib
            use chem_def, only: chem_isos
            integer :: i, k, ierr
            real(dp) :: alfa, rho_face, chiT_face, chiRho_face, &
               dm_dr, gradRho, dlnRho_dm

            ierr = 0
            write(*, *)
            write(*, *) 'write diffusion.data'
            write(*, *)
            open(newunit=iounit, file='diffusion.data', action='write', status='replace', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'failed to open diffusion dump file'
               return
            end if

            ! args
            write(iounit, '(99i20)') nz, nzlo, nzhi, species, nc, &
               s% diffusion_steps_hard_limit, &
               s% diffusion_nsmooth_typical_charge

            do i=1,species
               write(iounit,*) trim(chem_isos% name(s% chem_id(i)))
            end do

            write(iounit, '(99i20)') class(1:species)

            write(iounit, '(99i20)') chem_isos% Z(s% chem_id(1:species))

            write(iounit, '(99i20)') chem_isos% Z_plus_N(s% chem_id(1:species))

            write(iounit, '(99e22.10)') chem_isos% W(s% chem_id(1:species))

            do i=1,nc
               write(iounit, '(a)') trim(chem_isos% name(class_chem_id(i)))
            end do

            do i=1,nc
               write(iounit, '(a)') trim(class_name(i))
            end do

            if (s% diffusion_calculates_ionization) then
               write(iounit,*) 1
            else
               write(iounit,*) 0
            end if

            if (s% diffusion_use_iben_macdonald) then
               write(iounit,*) 1
            else
               write(iounit,*) 0
            end if

            if (s% diffusion_screening_for_radaccel) then
               write(iounit,*) 1
            else
               write(iounit,*) 0
            end if

            write(iounit, '(99(1pd26.16))') &
               s% mstar, s% dt, &
               s% diffusion_v_max, s% diffusion_gamma_full_on, s% diffusion_gamma_full_off, &
               s% diffusion_T_full_on, s% diffusion_T_full_off, &
               s% diffusion_max_T_for_radaccel, s% diffusion_dt_div_timescale, &
               s% R_center

            do k=1, nz
               write(iounit, '(99(1pd26.16))') &
                  gamma(k), s% abar(k), s% ye(k), free_e(k), s% dm_bar(k), s% dm(k), &
                  s% T(k), s% lnT(k), s% Rho(k), s% lnd(k), s% rmid(k), &
                  dlnPdm_mid(k), dlnT_dm_mid(k), dlnRho_dm_mid(k), &
                  s% L(k), s% r(k), dlnPdm_face(k), dlnT_dm_face(k), dlnRho_dm_face(k)
               write(iounit, '(99(1pd26.16))') s% xa(1:species, k)
            end do

            close(iounit)

         end subroutine dump_diffusion_info


         subroutine set_extras(ierr)
            integer, intent(out) :: ierr
            integer :: k, op_err
            op_err = 0
            ierr = 0
!$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(dynamic,2)
            do k=1,nz
               call set1_extras(k, op_err)
               if (op_err /= 0) ierr = op_err
            end do
!$OMP END PARALLEL DO
         end subroutine set_extras


         subroutine set1_extras(k,ierr)
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            real(dp) :: grav, area, P_face
            ierr = 0
            free_e(k) = exp(s% lnfree_e(k))
            gamma(k) = s% gam(k)
            if (k==1) then
               dlnPdm(k) = 0; dlnT_dm(k) = 0; return
            end if
            grav = -s% cgrav(k)*s% m(k)/s% r(k)**2
            area = pi4*s% r(k)**2
            P_face = 0.5d0*(s% Peos(k) + s% Peos(k-1))
            dlnPdm(k) = grav/(area*P_face) ! estimate based on QHSE
            dlnT_dm(k) = s% gradT(k)*dlnPdm(k)
         end subroutine set1_extras


      end subroutine do_element_diffusion
      
      subroutine finish_element_diffusion(s,dt)
        type (star_info), pointer :: s
        real(dp), intent(in) :: dt
        integer :: k
        
        do k=1,s% nz
           s% eps_diffusion(k) = (s% energy_start(k) - s% energy(k))/dt
        end do
        
      end subroutine finish_element_diffusion

      end module element_diffusion












