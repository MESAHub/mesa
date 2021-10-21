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


      module adjust_mass

      use star_private_def
      use const_def
      use utils_lib

      implicit none

      private
      public :: do_adjust_mass, compute_prev_mesh_dm


      logical, parameter :: dbg_adjm = .false.


      contains

      real(dp) function compute_delta_m(s) result(delta_m)
         type (star_info), pointer :: s
         include 'formats'
         delta_m = s% dt*s% mstar_dot

         if (s% super_eddington_wind_mdot /= 0 .and. s% super_eddington_wind_mdot > -s% mstar_dot) then
            ! NOTE: super_eddington should not be greater than max_wind
            if (s% max_wind > 0 .and. &
                  s% max_wind*Msun/secyer < s% super_eddington_wind_mdot) then
               s% mstar_dot = -s% max_wind*Msun/secyer
            else
               s% mstar_dot = -s% super_eddington_wind_mdot
            end if
         else if (delta_m == 0 &
            .or. (delta_m < 0 .and. s% star_mass <= s% min_star_mass_for_loss) &
            .or. (delta_m > 0 .and. s% max_star_mass_for_gain > 0 &
                  .and. s% star_mass >= s% max_star_mass_for_gain)) then
            s% mstar_dot = 0d0
         end if

         delta_m = s% dt*s% mstar_dot
         
         if (is_bad(s% dt)) then
            write(*,1) 's% dt', s% dt
            stop 'compute_delta_m'
         end if
         
         if (is_bad(s% mstar_dot)) then
            write(*,1) 's% mstar_dot', s% mstar_dot
            stop 'compute_delta_m'
         end if
         
         if (is_bad(delta_m)) then
            write(*,1) 'delta_m', delta_m
            stop 'compute_delta_m'
         end if

      end function compute_delta_m

      subroutine save_for_eps_mdot(s)
         use star_utils, only: eval_deltaM_total_energy_integrals
         type (star_info), pointer :: s

         integer :: j
            
         call eval_deltaM_total_energy_integrals( &
            s, 1, s% nz, s% mstar, .true., &
            s% total_energy_profile_before_adjust_mass, &
            s% total_internal_energy_before_adjust_mass, &
            s% total_gravitational_energy_before_adjust_mass, &
            s% total_radial_kinetic_energy_before_adjust_mass, &
            s% total_rotational_kinetic_energy_before_adjust_mass, &
            s% total_turbulent_energy_before_adjust_mass, &
            s% total_energy_before_adjust_mass)

         do j=1,s%nz
            s% dm_before_adjust_mass(j) = s% dm(j)
         end do

      end subroutine save_for_eps_mdot


      ! Reconstruct previous mesh based on what the mass differences should have been.
      ! This is to deal with the changes in cell mass being small (relatively speaking) at depth.
      ! That results in differences at the limit of double precision, so we need to go to quad precision.
      ! However dm and prev_mesh_dm are only known to double precision, whereas the difference between them
      ! is known to quad precision. So we take dm as given and reconstruct prev_mesh_dm using the difference.
      subroutine compute_prev_mesh_dm(s, prev_mesh_dm, dm, change_in_dm)
         ! Inputs
         type (star_info), pointer :: s

         ! Intermediates
         integer :: j, nz
         real(dp) :: change_sum

         ! Outputs
         real(qp), dimension(:), intent(out) :: prev_mesh_dm, dm, change_in_dm


         nz = s% nz

         change_sum = 0
         do j=1,nz
            dm(j) = s%dm(j)
            prev_mesh_dm(j) = s%dm_before_adjust_mass(j)
            if (j < s%k_below_const_q) then
               change_in_dm(j) = s%adjust_mass_outer_frac_sub1 * prev_mesh_dm(j)
            else if (j < s%k_const_mass) then
               change_in_dm(j) = s%adjust_mass_mid_frac_sub1 * prev_mesh_dm(j)
            else
               change_in_dm(j) = s%adjust_mass_inner_frac_sub1 * prev_mesh_dm(j)
            end if
            change_sum = change_sum + change_in_dm(j)
            prev_mesh_dm(j) = dm(j) - change_in_dm(j)
         end do

      end subroutine compute_prev_mesh_dm

      ! Updates the radii of cells after adjust mass. Note that the radius appears in many
      ! places in s. Currently those are: r, r_start, rmid, rmid_start, R2 (r^2), and lnR (log(r)).
      ! If more of these are added they should also be updated here.
      subroutine update_radius(s)
         use star_utils, only: store_r_in_xh, get_r_and_lnR_from_xh
         ! Inputs
         type (star_info), pointer :: s

         ! Intermediates
         integer j, nz
         real(dp) r_new, vol00, volp1, cell_vol         

         nz = s%nz

         vol00 = four_thirds_pi*s% R_center*s% R_center*s% R_center
         do j=nz,1,-1

            volp1 = vol00
            cell_vol = s% dm(j)/s% rho(j)
            vol00 = volp1 + cell_vol

            r_new = exp(log(vol00/four_thirds_pi)*one_third)

            s%r(j) = r_new
            s%r_start(j) = s%r(j)
            s%rmid(j) = r_new
            s%rmid_start(j) = r_new
            call store_r_in_xh(s, j, r_new)
            call get_r_and_lnR_from_xh(s, j, s% r(j), s% lnR(j))
            s% xh_start(s% i_lnR,j) = s% xh(s% i_lnR,j)
            
         end do

      end subroutine update_radius

      subroutine do_adjust_mass(s, species, ierr)
         use adjust_xyz, only: get_xa_for_accretion
         use star_utils, only: report_xa_bad_nums, &
            start_time, update_time
         use hydro_rotation, only: set_omega
         use chem_def

         type (star_info), pointer :: s
         integer, intent(in) :: species
         integer, intent(out) :: ierr

         real(dp) :: &
            dt, delta_m, old_mstar, new_mstar, old_J, new_J, factor, total, &
            frac, env_mass, mmax, alfa, new_xmstar, old_xmstar, removed, &
            q_for_just_added, xq_for_CpT_absMdot_div_L, sum_dq, dm, sumx

         real(dp), dimension(species) :: &
            xaccrete, mtot_init, mtot_final
         real(dp), dimension(:), allocatable :: &
            rxm_old, rxm_new, old_cell_mass, new_cell_mass, &
            oldloc, newloc, oldval, newval, xm_old, xm_new, &
            xq_center_old, xq_center_new
         real(dp), dimension(:,:), allocatable :: xa_old
         real(dp), pointer :: work(:)

         integer :: j, k, l, m, k_const_mass, nz, k_below_just_added, k_newval
         integer(8) :: time0
         real(dp) :: partial_xa_mass, region_total_mass
         real(dp) :: starting_j_rot_surf

         logical, parameter :: dbg = .false.

         include 'formats'

         if (dbg) write(*,*) 'do_adjust_mass'

         call save_for_eps_mdot(s)

         ierr = 0

         nz = s% nz
         dt = s% dt
         

         s% k_const_mass = 1

         s% k_for_test_CpT_absMdot_div_L = 1

         q_for_just_added = 1d0
         s% k_below_just_added = 1

         ! NOTE: don't assume that vars are all set at this point.
         ! exceptions: omega, i_rot, j_rot, and total_angular_momentum have been set.
         ! store surface j_rot, to later adjust angular momentum lost
         if (s% rotation_flag) then
            starting_j_rot_surf = s% j_rot(1)
         end if

         delta_m = compute_delta_m(s)
         
         if (delta_m == 0) then
            return
         end if

         if (s% doing_timing) call start_time(s, time0, total)

         s% need_to_setvars = .true.

         old_mstar = s% mstar
         old_xmstar = s% xmstar

         new_mstar = old_mstar + delta_m
         new_xmstar = old_xmstar + delta_m
         
         if (is_bad(new_xmstar)) then
            write(*,1) 'new_xmstar', new_xmstar
            stop 'do_adjust_mass'
         end if

         if (delta_m > 0 .and. s% max_star_mass_for_gain > 0 &
               .and. new_mstar > Msun*s% max_star_mass_for_gain) then
            new_mstar = Msun*s% max_star_mass_for_gain
            delta_m = new_mstar - old_mstar
            s% mstar_dot = delta_m/dt
         else if (delta_m < 0 .and. new_mstar < Msun*s% min_star_mass_for_loss) then
            new_mstar = Msun*s% min_star_mass_for_loss
            delta_m = new_mstar - old_mstar
            s% mstar_dot = delta_m/dt
         end if

         frac = old_xmstar/new_xmstar
         new_xmstar = old_xmstar/frac
         if (new_xmstar <= 0) then
            ierr = -1
            return
         end if
         s% xmstar = new_xmstar
         s% mstar = s% xmstar + s% M_center

         if (dbg_adjm) then
            env_mass = old_mstar - s% he_core_mass*Msun
            write(*,'(a40,f26.16)') 'env_mass/old_mstar', env_mass/old_mstar
            write(*,*)
            write(*,1) 'delta_m/old_mstar', delta_m/old_mstar
            write(*,1) 's% he_core_mass*Msun', s% he_core_mass*Msun
            write(*,1) 'env_mass', env_mass
            write(*,1) 'delta_m/env_mass', delta_m/env_mass
            write(*,1) 'log10(abs(delta_m/env_mass))', safe_log10(abs(delta_m/env_mass))
            write(*,*)
         end if

         call do_alloc(ierr)
         if (ierr /= 0) return

         do k=1,nz
            old_cell_mass(k) = old_xmstar*s% dq(k)
         end do
         xm_old(1) = 0
         xq_center_old(1) = 0.5d0*s% dq(1)
         do k=2,nz
            xm_old(k) = xm_old(k-1) + old_cell_mass(k-1)
            xq_center_old(k) = xq_center_old(k-1) + 0.5d0*(s% dq(k-1) + s% dq(k))
         end do

         s% xa_removed(1:species) = 0
         s% angular_momentum_removed = 0
         s% angular_momentum_added = 0
         if (delta_m < 0) then
            s% angular_momentum_removed = angular_momentum_removed(ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            removed = -delta_m
            do k=1,nz
               if (xm_old(k) >= removed) then
                  do j=1,species
                     s% xa_removed(j) = s% xa_removed(j)/removed
                  end do
                  exit
               end if
               do j=1,species
                  dm = min(old_cell_mass(k), removed-xm_old(k))
                  s% xa_removed(j) = s% xa_removed(j) + s% xa(j,k)*dm
               end do
            end do
         end if
         !s% angular_momentum_added is set after j_rot is adjusted in outer layers

         call revise_q_and_dq( &
            s, nz, old_xmstar, new_xmstar, delta_m, k_const_mass, ierr)
         if (ierr /= 0) then
            s% retry_message = 'revise_q_and_dq failed in adjust mass'
            if (s% report_ierr) write(*,*) s% retry_message
            call dealloc
            return
         end if

         mtot_init(1:species) = 0
         do k=1,k_const_mass ! save total mass by species down to where constant
            sumx = 0
            do j=1,species
               s% xa(j,k) = max(0d0, min(1d0, s% xa(j,k)))
               sumx = sumx + s% xa(j,k)
            end do
            do j=1,species
               s% xa(j,k) = s% xa(j,k)/sumx
               mtot_init(j) = mtot_init(j) + old_cell_mass(k)*s% xa(j,k)
            end do
         end do

         do k=1,nz
            new_cell_mass(k) = new_xmstar*s% dq(k)
         end do

         xm_new(1) = 0
         xq_center_new(1) = 0.5d0*s% dq(1)
         do k=2,nz
            xm_new(k) = xm_new(k-1) + new_cell_mass(k-1)
            xq_center_new(k) = xq_center_new(k-1) + 0.5d0*(s% dq(k-1) + s% dq(k))
         end do

         xq_for_CpT_absMdot_div_L = abs(delta_m)/new_xmstar
         s% k_for_test_CpT_absMdot_div_L = nz
         sum_dq = 0
         do k = 1, nz
            if (sum_dq >= xq_for_CpT_absMdot_div_L) then
               s% k_for_test_CpT_absMdot_div_L = k; exit
            end if
            sum_dq = sum_dq + s% dq(k)
         end do

         if (delta_m > 0) then
            q_for_just_added = old_xmstar/new_xmstar
            s% k_below_just_added = nz
            do k = 1, nz
               if (s% q(k) <= q_for_just_added) then
                  s% k_below_just_added = k; exit
               end if
            end do
         end if
         
         k_below_just_added = s% k_below_just_added
         k_newval = 1

         do k=1,nz
            do j=1,species
               xa_old(j,k) = s% xa(j,k)
            end do
         end do

         if (delta_m < 0) then
            xaccrete(1:species) = 0 ! xaccrete not used when removing mass
         else ! set xaccrete for composition of added material
            if (s% accrete_same_as_surface) then
               do j=1,species
                  xaccrete(j) = xa_old(j,1)
               end do
            else
               call get_xa_for_accretion(s, xaccrete, ierr)
               if (ierr /= 0) then
                  s% retry_message = 'get_xa_for_accretion failed in adjust mass'
                  if (s% report_ierr) write(*, *) s% retry_message
                  call dealloc
                  return
               end if
            end if
         end if

         mmax = max(old_mstar, new_mstar)

         ! rxm_old and rxm_new are for interpolating by mass
         ! but instead of using mass coord, we want to use external mass
         ! in order to get better accuracy near the surface.
         ! to simplify this, the zero point is the same for both rxm_old and rxm_new.
         ! that makes rxm_new = rxm_old for k >= k_const_mass.
         if (delta_m < 0) then
            rxm_old(1) = 0
            rxm_new(1) = -delta_m ! note that rxm_new(1) > 0 since delta_m < 0
         else
            rxm_old(1) = delta_m
            rxm_new(1) = 0
         end if
         do k = 2, nz
            rxm_old(k) = rxm_old(k-1) + old_cell_mass(k-1)
            if (k >= k_const_mass) then
               rxm_new(k) = rxm_old(k)
               new_cell_mass(k) = old_cell_mass(k)
            else
               rxm_new(k) = rxm_new(k-1) + new_cell_mass(k-1)
            end if
         end do

         call set_xa(s, nz, k_const_mass, species, xa_old, xaccrete, &
            rxm_old, rxm_new, mmax, old_cell_mass, new_cell_mass, ierr)
         if (ierr /= 0) then
            s% retry_message = 'set_xa failed in adjust mass'
            if (s% report_ierr) write(*,*) s% retry_message
            call dealloc
            return
         end if

         if (s% rotation_flag) then
            ! update j_rot if have extra angular momentum change
            call set_omega_adjust_mass( &
               s, nz, k_const_mass, k_below_just_added, delta_m, &
               rxm_old, rxm_new, mmax, old_cell_mass, new_cell_mass, &
               ! reuse these work arrays
               oldloc, newloc, oldval, newval, xm_old, xm_new, &
               ierr)
            if (ierr /= 0) then
               s% retry_message = 'set_omega_adjust_mass failed in adjust mass'
               if (s% report_ierr) write(*,*) s% retry_message
               call dealloc
               return
            end if
            if (s% do_adjust_J_lost) then
               call adjust_J_lost(s, k_below_just_added, starting_j_rot_surf, ierr)
               if (ierr /= 0) then
                  s% retry_message = 'do_adjust_J_lost failed in adjust mass'
                  if (s% report_ierr) write(*,*) s% retry_message
                  call dealloc
                  return
               end if
            end if
            if (s% D_omega_flag) then
               call set_D_omega( &
                  s, nz, k_const_mass, k_newval, &
                  rxm_old, rxm_new, delta_m, old_xmstar, new_xmstar, &
                  s% D_omega, oldloc, newloc, oldval, newval, work, ierr)
               if (ierr /= 0) then
                  s% retry_message = 'set_D_omega failed in adjust mass'
                  if (s% report_ierr) write(*,*) s% retry_message
                  call dealloc
                  return
               end if
            end if
            if (s% am_nu_rot_flag) then
               call set_D_omega( & ! reuse the set_D_omega routine to set am_nu_rot
                  s, nz, k_const_mass, k_newval, &
                  rxm_old, rxm_new, delta_m, old_xmstar, new_xmstar, &
                  s% am_nu_rot, oldloc, newloc, oldval, newval, work, ierr)
               if (ierr /= 0) then
                  s% retry_message = 'set_am_nu_rot failed in adjust mass'
                  if (s% report_ierr) write(*,*) s% retry_message
                  call dealloc
                  return
               end if
            end if
         end if

         ! soften outer xa    (Pablo -- this should be the very last step)
         if (s% smooth_outer_xa_big > 0.0d0 .and. s% smooth_outer_xa_small > 0.0d0) then
            m = 1
            do k = 1, s% nz
               if (s% m(1) - s% m(k) > s% smooth_outer_xa_big * s% m(1)) exit
               region_total_mass = 0
               m = k
               do l = k, s% nz
                  if (s% m(k) - s% m(l) > s% smooth_outer_xa_small * s% m(1) * &
                     (1-(s% m(1) - s% m(k))/(s% smooth_outer_xa_big * s% m(1)))) exit
                  m = l
                  region_total_mass = region_total_mass + s% dm(l)
               end do
               if (m == k) cycle
               do j=1,s% species
                  partial_xa_mass = 0
                  do l = k, m
                     partial_xa_mass = partial_xa_mass + s% dm(l) * s% xa(j,l)
                  end do
                  do l = k, m
                     s% xa(j,l) = partial_xa_mass / region_total_mass
                  end do
               end do
            end do
            ! fix small errors to ensure xa's add up to unity
            do k=1,m
               frac = 1d0/sum(s% xa(1:s% species,k))
               do j=1,s% species
                  s% xa(j,k) = s% xa(j,k)*frac
               end do
            end do
         end if

         
         call update_radius(s)
         

         call dealloc
         

         if (s% doing_timing) call update_time(s, time0, total, s% time_adjust_mass)

         if (dbg_adjm) stop 'debugging: do_adjust_mass'
         if (dbg) write(*,*) 'do_adjust_mass return'

         contains


         real(dp) function angular_momentum_removed(ierr) result(J)
            ! when call this, s% j_rot is still for old mass
            integer, intent(out) :: ierr
            integer :: k
            real(dp) :: r2, dmm1, dm00, dm, dm_sum, dm_lost
            include 'formats'
            ierr = 0
            J = 0
            if (.not. s% rotation_flag) return
            dm00 = 0
            dm_sum = 0
            dm_lost = -delta_m
            do k = 1, nz
               dmm1 = dm00
               dm00 = old_cell_mass(k)
               dm = 0.5d0*(dmm1+dm00)
               if (dm_sum + dm > dm_lost) then
                  dm = dm_lost - dm_sum
                  dm_sum = dm_lost
               else
                  dm_sum = dm_sum + dm
               end if
               J = J + dm*s% j_rot(k)
               if (dm_sum == dm_lost) exit
            end do
         end function angular_momentum_removed


         real(dp) function eval_total_angular_momentum(s,cell_mass,nz_last) result(J)
            type (star_info), pointer :: s
            real(dp) :: cell_mass(:)
            integer, intent(in) :: nz_last
            integer :: k
            real(dp) :: dmm1, dm00, dm
            include 'formats'
            J = 0
            if (.not. s% rotation_flag) return
            dm00 = 0
            do k = 1, nz_last
               dmm1 = dm00
               dm00 = cell_mass(k)
               if (k == s% nz) then
                  dm = 0.5d0*dmm1+dm00
               else if (k == nz_last) then
                  dm = 0.5d0*dmm1
               else
                  dm = 0.5d0*(dmm1+dm00)
               end if
               J = J + dm*s% j_rot(k)
            end do
         end function eval_total_angular_momentum
            
         subroutine do_alloc(ierr)
            integer, intent(out) :: ierr
            allocate(rxm_old(nz), rxm_new(nz), old_cell_mass(nz), new_cell_mass(nz), &
               xa_old(species,nz), oldloc(nz), newloc(nz), oldval(nz), newval(nz), &
               xm_old(nz), xm_new(nz), xq_center_old(nz), xq_center_new(nz))
            call do_work_arrays(.true.,ierr)
         end subroutine do_alloc

         subroutine dealloc
            integer :: ierr
            call do_work_arrays(.false.,ierr)
         end subroutine dealloc

         subroutine do_work_arrays(alloc_flag, ierr)
            use interp_1d_def
            use alloc, only: work_array
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
                work, nz*pm_work_size, nz_alloc_extra, 'adjust_mass work', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays            

      end subroutine do_adjust_mass


      subroutine revise_q_and_dq( &
            s, nz, old_xmstar, new_xmstar, delta_m, k_const_mass, ierr)
         use star_utils, only: get_lnT_from_xh
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         real(dp), intent(in) :: old_xmstar, new_xmstar, delta_m
         integer, intent(out) :: k_const_mass, ierr

         integer :: k, kA, kB, j00, jp1, k_check
         real(dp) :: lnTlim_A, lnTlim_B, sumdq, sumdq1, sumdq2, sumdq3, &
            min_xq_const_mass, min_q_for_kB, mold_o_mnew, lnTmax, lnT_A, lnT_B, &
            xqA, xqB_old, xqB_new, qfrac, frac, dqacc
         real(dp) :: xq(nz)
         real(qp) :: qfrac_qp, frac_qp, mold_o_mnew_qp, q1, q2, q3, q4
         real(qp) :: adjust_mass_outer_frac, adjust_mass_mid_frac, adjust_mass_inner_frac
         integer, parameter :: min_kA = 5
         logical :: dbg, flag
         logical :: okay_to_move_kB_inward

         include 'formats'

         ierr = 0
         dbg = .false.
         flag = .false.
         
         if (is_bad(old_xmstar)) then
            write(*,1) 'old_xmstar', old_xmstar
            stop 'revise_q_and_dq'
         end if

         if (is_bad(new_xmstar)) then
            write(*,1) 'new_xmstar', new_xmstar
            stop 'revise_q_and_dq'
         end if

         if (is_bad(delta_m)) then
            write(*,1) 'delta_m', delta_m
            stop 'revise_q_and_dq'
         end if


         

         okay_to_move_kB_inward = .false.

         lnTlim_A = ln10*s% max_logT_for_k_below_const_q
         lnTlim_B = ln10*s% max_logT_for_k_const_mass
         
         q1 = old_xmstar
         q2 = new_xmstar
         mold_o_mnew_qp = q1/q2
         mold_o_mnew = mold_o_mnew_qp
         s% max_q_for_k_below_const_q = mold_o_mnew
         dqacc = delta_m / new_xmstar

         ! might have drifted away from summing to 1 in mesh adjustment, fix up
         sumdq = 0
         do k=1,nz
            sumdq = sumdq + s% dq(k)
         end do
         frac = 1.0d0/sumdq
         if (is_bad(frac)) then
            write(*,1) 'frac for initial renorm', frac
            stop 'revise_q_and_dq'
         end if
         do k = 1, nz
            if (is_bad(s% dq(k))) then
               write(*,2) 'bad dq input', s% dq(k)
               stop 'revise_q_and_dq'
            end if
            s% dq(k) = s% dq(k) * frac
         end do

         ! will work in xq := 1-q to reduce truncation errors
         xq(1)=0
         do k = 2, nz
            xq(k) = xq(k-1) + s% dq(k-1)
         end do
         
         k = maxloc(s% xh(s% i_lnT,1:nz),dim=1)
         lnTmax = get_lnT_from_xh(s, k)
         lnT_A = min(lnTmax, lnTlim_A)
         
         if (is_bad(s% max_q_for_k_below_const_q)) then
            write(*,*) 's% max_q_for_k_below_const_q', s% max_q_for_k_below_const_q
            stop 'revise_q_and_dq'
         end if
         if (is_bad(s% min_q_for_k_below_const_q)) then
            write(*,*) 's% min_q_for_k_below_const_q', s% min_q_for_k_below_const_q
            stop 'revise_q_and_dq'
         end if
         
         kA = min_kA
         do k = min_kA, nz-1
            kA = k
            if ( (1-xq(k)) > s% max_q_for_k_below_const_q) cycle
            if ( 1.0d0-xq(k) <= s% min_q_for_k_below_const_q) exit
            if (get_lnT_from_xh(s, k) >= lnT_A) exit
         end do
         xqA = xq(kA)

         lnT_B = min(lnTmax, lnTlim_B)
         
         if (is_bad(s% max_q_for_k_const_mass)) then
            write(*,*) 's% max_q_for_k_const_mass', s% max_q_for_k_const_mass
            stop 'revise_q_and_dq'
         end if
         if (is_bad(s% min_q_for_k_const_mass)) then
            write(*,*) 's% min_q_for_k_const_mass', s% min_q_for_k_const_mass
            stop 'revise_q_and_dq'
         end if

         kB = kA+1
         do k = kB, nz
            kB = k
            if ( (1-xq(k)) > s% max_q_for_k_const_mass) cycle
            if ( (1-xq(k)) <= s% min_q_for_k_const_mass) exit
            if (get_lnT_from_xh(s, k) >= lnT_B) exit
         end do

         xqB_old = xq(kB)

         if (dbg_adjm) then
            write(*,*) 'kA', kA
            write(*,1) 'xqA', xqA
            write(*,*) 'kB', kB
            write(*,1) 'xqB_old', xqB_old
            write(*,*)
            write(*,1) 'xqA-xqB_old', xqA-xqB_old
            write(*,*)
         end if

         xqB_new = dqacc + xqB_old*mold_o_mnew  ! in order to keep m interior to kB constant
         ! same as  qB_new = qB_old*mold/mnew, but without the truncation problems

         do ! make sure qfrac is not too far from 1 by moving kB inward
            q1 = xqB_new - xqA
            q2 = max(1d-99,xqB_old - xqA)
            qfrac_qp = q1/q2
            qfrac = qfrac_qp
            if (kB == nz) exit
            if (kB-kA > 10) then
               if (qfrac > 0.67d0 .and. qfrac < 1.5d0) exit
               if (qfrac > 0.50d0 .and. qfrac < 2.0d0) then
                  j00 = maxloc(s% xa(:,kB),dim=1) ! most abundant species at kB
                  jp1 = maxloc(s% xa(:,kB+1),dim=1) ! most abundant species at kB+1
                  if (j00 /= jp1) then ! change in composition.
                     if (dbg) write(*,*) 'change in composition.  back up kB.'
                     kB = max(1,kB-5)
                     exit
                  end if
               end if
            end if
            kB = kB+1
            xqB_old = xq(kB)
            xqB_new = dqacc + xqB_old*mold_o_mnew
         end do

         ! set new dq's
         ! s% dq(1:kA-1) unchanged
         s% dq(kA:kB-1) = s% dq(kA:kB-1)*qfrac
         frac_qp = mold_o_mnew_qp
         frac = frac_qp
         if (is_bad(frac)) then
            write(*,1) 'frac for kA:kB-1', frac
            stop 'revise_q_and_dq'
         end if
         s% dq(kB:nz) = s% dq(kB:nz)*frac
         
         adjust_mass_outer_frac = 1d0
         adjust_mass_mid_frac = qfrac_qp
         adjust_mass_inner_frac = frac_qp

         if (dbg_adjm) then
            write(*,1) 'frac for kb:nz', frac
            write(*,1) 'qfrac for kA:kB-1', qfrac
            write(*,1) 'revise_q_and_dq sum dqs', sum(s% dq(1:nz))
            write(*,2) 'qfrac region', kB, qfrac, s% q(kB), s% lnT(kB)/ln10
            write(*,2) 'frac region', kA, frac, s% q(kA), s% lnT(kA)/ln10
            write(*,2) 'kA', kA
            write(*,2) 'kB', kB
            write(*,2) 'nz', nz
            write(*,*)
            stop 'adjust_mass'
         end if

         ! renorm
         sumdq = 0
         do k=1,nz
            sumdq = sumdq + s% dq(k)
         end do
         q1 = 1d0
         q2 = sumdq
         frac_qp = q1/q2
         frac = frac_qp
         if (is_bad(frac)) then
            write(*,1) 'frac for renorm', frac
            stop 'revise_q_and_dq'
         end if
         do k = 1, nz
            s% dq(k) = s% dq(k) * frac
         end do
         
         s% adjust_mass_outer_frac_sub1 = frac_qp*adjust_mass_outer_frac*new_xmstar / old_xmstar - 1.0_qp
         s% adjust_mass_mid_frac_sub1 = frac_qp*adjust_mass_mid_frac*new_xmstar / old_xmstar - 1.0_qp
         s% adjust_mass_inner_frac_sub1 = frac_qp*adjust_mass_inner_frac*new_xmstar / old_xmstar - 1.0_qp

         ! set q's to match new dq's
         s% q(1) = 1d0
         sumdq = s% dq(1)
         do k = 2, nz-1
            s% q(k) = 1d0 - sumdq
            sumdq = sumdq + s% dq(k)
         end do
         s% q(nz) = 1d0 - sumdq

         s% dq(nz) = s% q(nz)

         ! update dm's for new dq's
         do k=1,nz
            s% dm(k) = s% dq(k) * new_xmstar
         end do

         s% k_below_const_q = kA
         s% k_const_mass = kB
         k_const_mass = kB

      end subroutine revise_q_and_dq


      subroutine set_xa( &
            s, nz, k_const_mass, species, xa_old, xaccrete, &
            old_cell_xbdy, new_cell_xbdy, mmax, old_cell_mass, new_cell_mass, ierr)
         ! set new values for s% xa(:,:)
         type (star_info), pointer :: s
         integer, intent(in) :: nz, k_const_mass, species
         real(dp), intent(in) :: mmax
         real(dp), intent(in) :: xa_old(:, :), xaccrete(:)
         real(dp), dimension(:), intent(in) :: &
            old_cell_xbdy, new_cell_xbdy, old_cell_mass, new_cell_mass ! (nz)
         integer, intent(out) :: ierr
         integer :: k, j, op_err
         real(dp), parameter :: max_sum_abs = 10d0
         real(dp), parameter :: xsum_tol = 1d-2
         include 'formats'
         ierr = 0
         if (dbg_adjm) &
            write(*,2) 'set_xa: k_const_mass', k_const_mass
         if (k_const_mass < nz) then
            ! for k >= k_const_mass have m_new(k) = m_old(k),
            ! so no change in xa_new(:,k) for k > k_const_mass
            do k=k_const_mass+1,nz
               do j=1,species
                  s% xa(j,k) = xa_old(j,k)
               end do
            end do
         end if
!$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(dynamic,2)
         do k = 1, k_const_mass
            op_err = 0
            call set1_xa(s, k, nz, species, xa_old, xaccrete, &
               old_cell_xbdy, new_cell_xbdy, mmax, old_cell_mass, new_cell_mass, op_err)
            if (op_err /= 0) then
               if (s% report_ierr) write(*,2) 'set1_xa error', k
               ierr = op_err
            end if
         end do
!$OMP END PARALLEL DO

      end subroutine set_xa


      subroutine set1_xa(s, k, nz, species, xa_old, xaccrete, &
            old_cell_xbdy, new_cell_xbdy, mmax, old_cell_mass, new_cell_mass, ierr)
         ! set new values for s% xa(:,k)
         use num_lib, only: binary_search
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: k, nz, species
         real(dp), intent(in) :: mmax
         real(dp), intent(in) :: xa_old(:,:), xaccrete(:)
         real(dp), dimension(:), intent(in) :: &
            old_cell_xbdy, new_cell_xbdy, old_cell_mass, new_cell_mass
         integer, intent(out) :: ierr

         real(dp) :: msum(species), xm_inner, sumx, &
            xm0, xm1, new_cell_dm, dm_sum, dm, xm_outer
         integer :: kk, k_outer, j

         integer, parameter :: k_dbg = -1
         logical, parameter :: xa_dbg = .false.

         logical, parameter :: do_not_mix_accretion = .false.

         logical :: partially_accreted_cell

         include 'formats'

         ierr = 0
         msum(:) = -1 ! for testing

         xm_outer = new_cell_xbdy(k)
         if (k == nz) then
            new_cell_dm = mmax - xm_outer - s% M_center
         else
            new_cell_dm = new_cell_mass(k)
         end if
         xm_inner = xm_outer + new_cell_dm

         dm_sum = 0d0

         partially_accreted_cell = .false.
         if (xm_outer < old_cell_xbdy(1)) then ! there is some accreted material in new cell
            if (do_not_mix_accretion .or. xm_inner <= old_cell_xbdy(1)) then
               ! new cell is entirely accreted material
               !write(*,2) 'new cell is entirely accreted material', k, new_cell_dm
               do j=1,species
                  s% xa(j,k) = xaccrete(j)
               end do
               return
            end if
            dm = min(new_cell_dm, old_cell_xbdy(1) - xm_outer)
            dm_sum = dm
            do j=1,species
               msum(j) = xaccrete(j)*dm
            end do
            partially_accreted_cell = .true.
            xm_outer = old_cell_xbdy(1)
            k_outer = 1
         else ! new cell entirely composed of old material
            msum(:) = 0
            if (xm_outer >= old_cell_xbdy(nz)) then
               ! new cell contained entirely in old cell nz
               k_outer = nz
            else
               ! binary search for k_outer such that
               ! xm_outer >= old_cell_xbdy(k_outer)
               ! and old_cell_xbdy(k_outer+1) > xm_outer
               k_outer = binary_search(nz, old_cell_xbdy, 0, xm_outer)

               ! check
               if (k_outer <= 0 .or. k_outer > nz) then

                  ierr = -1
                  if (.not. xa_dbg) return

                  write(*,2) 'k', k
                  write(*,2) 'k_outer', k_outer
                  write(*,1) 'xm_outer', xm_outer
                  write(*,2) 'old_cell_xbdy(1)', 1, old_cell_xbdy(1)
                  write(*,2) 'old_cell_xbdy(nz)', nz, old_cell_xbdy(nz)
                  stop 'debugging: set1_xa'
               end if

               if (xm_outer < old_cell_xbdy(k_outer)) then

                  ierr = -1
                  if (.not. xa_dbg) return

                  write(*,*) 'k', k
                  write(*,*) 'k_outer', k_outer
                  write(*,1) 'xm_outer', xm_outer
                  write(*,1) 'old_cell_xbdy(k_outer)', old_cell_xbdy(k_outer)
                  write(*,*) '(xm_outer < old_cell_xbdy(k_outer))'
                  stop 'debugging: set1_xa'
               end if

               if (k_outer < nz) then
                  if (old_cell_xbdy(k_outer+1) <= xm_outer) then

                     ierr = -1
                     if (.not. xa_dbg) return

                     write(*,*) 'k', k
                     write(*,*) 'k_outer', k_outer
                     write(*,1) 'xm_outer', xm_outer
                     write(*,1) 'old_cell_xbdy(k_outer+1)', old_cell_xbdy(k_outer+1)
                     write(*,*) '(old_cell_xbdy(k_outer+1) <= xm_outer)'
                     stop 'debugging: set1_xa'
                  end if
               end if

            end if
         end if

         if (k == -1) then
            ierr = -1
            if (.not. xa_dbg) return

            write(*,2) 'nz', nz
            write(*,2) 'k_outer', k_outer
            write(*,1) 'xm_outer', xm_outer
            write(*,1) 'xm_inner', xm_inner
         end if

         do kk = k_outer, nz ! loop until reach m_inner
            xm0 = old_cell_xbdy(kk)

            if (xm0 >= xm_inner) then
               if (dm_sum < new_cell_dm .and. kk > 1) then
                  ! need to add a bit more from the previous source cell
                  dm = new_cell_dm - dm_sum
                  dm_sum = new_cell_dm
                  do j=1,species
                     msum(j) = msum(j) + xa_old(j,kk-1)*dm
                  end do
               end if
               exit
            end if

            if (kk == nz) then
               xm1 = mmax - s% M_center
            else
               xm1 = old_cell_xbdy(kk+1)
            end if

            if (xm1 < xm_outer) then
               ierr = -1
               if (.not. xa_dbg) return
               write(*,*)
               write(*,*) 'k', k
               write(*,*) 'kk', kk
               write(*,1) 'xm1', xm1
               write(*,1) 'xm_outer', xm_outer
               write(*,*) 'xm1 < xm_outer'
               stop 'debugging: set1_xa'
            end if

            if (xm0 >= xm_outer .and. xm1 <= xm_inner) then
               ! entire old cell kk is in new cell k

               dm = old_cell_mass(kk)
               dm_sum = dm_sum + dm

               if (dm_sum > new_cell_dm) then
                  ! dm too large -- numerical roundoff problems
                  dm = dm - (new_cell_dm - dm_sum)
                  dm_sum = new_cell_dm
               end if

               do j=1,species
                  msum(j) = msum(j) + xa_old(j,kk)*dm
               end do

            else if (xm0 <= xm_outer .and. xm1 >= xm_inner) then
               ! entire new cell k is in old cell kk
               !  except if it is a partially accreted cell for which xm_outer has been reset
               if (partially_accreted_cell) then
                  dm = xm_inner-xm_outer
               else
                  dm = new_cell_mass(k)
               end if
               dm_sum = dm_sum + dm
               do j=1,species
                  msum(j) = msum(j) + xa_old(j,kk)*dm
               end do

            else ! only use the part of old cell kk that is in new cell k

               if (xm_inner <= xm1) then ! this is the last part of new cell k

                  dm = new_cell_dm - dm_sum
                  dm_sum = new_cell_dm

               else ! notice that we avoid this case if possible because of numerical roundoff

                  dm = max(0d0, xm1 - xm_outer)
                  if (dm_sum + dm > new_cell_dm) dm = new_cell_dm - dm_sum
                  dm_sum = dm_sum + dm

               end if

               do j=1,species
                  msum(j) = msum(j) + xa_old(j,kk)*dm
               end do

               if (dm <= 0) then
                  ierr = -1
                  if (.not. xa_dbg) return
                  write(*,*) 'dm <= 0', dm
                  stop 'debugging: set1_xa'
               end if

            end if

            if (dm_sum >= new_cell_dm) then
               exit
            end if

         end do

         ! revise and renormalize
         do j=1,species
            s% xa(j,k) = msum(j) / new_cell_mass(k)
         end do
         sumx = sum(s% xa(:,k))

         sumx = sum(s% xa(:,k))
         if (sumx <= 0d0) then
            ierr = -1
            return
         end if

         do j=1,species
            s% xa(j,k) = s% xa(j,k)/sumx
         end do

      end subroutine set1_xa


      subroutine set_omega_adjust_mass( &
            s, nz, k_const_mass, k_below_just_added, delta_m, &
            old_cell_xbdy, new_cell_xbdy, mmax, old_cell_mass, new_cell_mass, &
            ! work arrays (nz)
            old_xout, new_xout, old_dmbar, new_dmbar, old_j_rot, extra_work, &
            ierr)
         use star_utils, only: total_angular_momentum
         type (star_info), pointer :: s
         integer, intent(in) :: nz, k_const_mass, k_below_just_added
         real(dp), intent(in) :: mmax, delta_m
         real(dp), dimension(:), intent(in) :: &
            old_cell_xbdy, new_cell_xbdy, old_cell_mass, new_cell_mass ! (nz)
         real(dp), dimension(:) :: &
            old_xout, new_xout, old_dmbar, new_dmbar, old_j_rot, extra_work
         integer, intent(out) :: ierr

         integer :: k, k0, op_err, old_k, new_k, k_uniform
         logical :: okay
         real(dp) :: old_j_tot, new_j_tot, goal_total_added, actual_total_added, &
            f, jtot_bdy, goal_total, bdy_j, bdy_total, inner_total, outer_total, &
            msum, isum, jsum, omega_uniform

         include 'formats'

         ierr = 0

         old_xout(1) = old_cell_xbdy(1)
         new_xout(1) = new_cell_xbdy(1)
         old_dmbar(1) = old_cell_mass(1)/2
         new_dmbar(1) = new_cell_mass(1)/2
         old_j_rot(1) = s% j_rot(1)
         do k=2,nz
            old_xout(k) = old_xout(k-1) + old_dmbar(k-1)
            new_xout(k) = new_xout(k-1) + new_dmbar(k-1)
            old_dmbar(k) = (old_cell_mass(k-1) + old_cell_mass(k))/2
            new_dmbar(k) = (new_cell_mass(k-1) + new_cell_mass(k))/2
            old_j_rot(k) = s% j_rot(k)
         end do
         old_dmbar(nz) = old_cell_mass(nz-1)/2 + old_cell_mass(nz)
         new_dmbar(nz) = new_cell_mass(nz-1)/2 + new_cell_mass(nz)

         old_j_tot = dot_product(s% j_rot(1:nz),old_dmbar(1:nz))

         if (k_below_just_added == 1) then
            k0 = 1
         else
            k0 = k_below_just_added + 1
         end if
!$OMP PARALLEL DO PRIVATE(k, op_err) SCHEDULE(dynamic,2)
         do k = 1, k_const_mass
            if (k < k0) then
               cycle
            end if
            op_err = 0
            call set1_omega(s, k, k_below_just_added, nz, &
               old_xout, new_xout, mmax, old_dmbar, new_dmbar, old_j_rot, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO

         if (k_below_just_added > 1) then
            ! set omega in cells with newly added material
            if (s% use_accreted_material_j) then
               actual_total_added = 0d0
               do k=1,k_below_just_added-2 ! remaining 2 done below
                  s% j_rot(k) = s% accreted_material_j
                  call set1_irot(s, k, k_below_just_added, .true.)
                  s% omega(k) = s% j_rot(k)/s% i_rot(k)
                  actual_total_added = actual_total_added + s% j_rot(k)*new_dmbar(k)
               end do
               k = k_below_just_added
               goal_total_added = delta_m*s% accreted_material_j
               goal_total = old_j_tot + goal_total_added
               inner_total = dot_product(s% j_rot(k+1:nz),new_dmbar(k+1:nz))
               outer_total = dot_product(s% j_rot(1:k-2),new_dmbar(1:k-2))
               bdy_total = goal_total - (inner_total + outer_total)
               bdy_j = bdy_total/sum(new_dmbar(k-1:k))
               if (bdy_j <= 0) then
                  ierr = -1
                  return
               end if
               do k=k_below_just_added-1,k_below_just_added
                  s% j_rot(k) = bdy_j
                  call set1_irot(s, k, k_below_just_added, .true.)
                  s% omega(k) = s% j_rot(k)/s% i_rot(k)
               end do
            else ! use old surface omega in all the new material
               do k=1,k_below_just_added-1
                  s% omega(k) = s% omega(k_below_just_added)
                  call set1_irot(s, k, k_below_just_added, .false.)
                  s% j_rot(k) = s% omega(k)*s% i_rot(k)
               end do
            end if
         end if

         okay = .true.
         do k=1,nz
            if (is_bad(s% omega(k)) .or. &
              abs(s% omega(k)) > 1d50) then
               if (s% report_ierr) then
!$OMP critical (star_adjust_mass_omega)
                  write(*,2) 's% omega(k)', k, s% omega(k)
!$OMP end critical (star_adjust_mass_omega)
               end if
               if (s% stop_for_bad_nums) then
                  write(*,2) 's% omega(k)', k, s% omega(k)
                  stop 'set_omega_adjust_mass'
               end if
               okay = .false.
            end if
         end do
         if (.not. okay) then
            write(*,2) 'model_number', s% model_number
            stop 'set_omega_adjust_mass'
         end if

      end subroutine set_omega_adjust_mass


      subroutine set1_irot(s, k, k_below_just_added, jrot_known) ! using lnR_for_d_dt_const_m
         use hydro_rotation, only: eval_i_rot, w_div_w_roche_jrot, w_div_w_roche_omega
         use star_utils, only: get_r_from_xh
         type (star_info), pointer :: s
         integer, intent(in) :: k, k_below_just_added
         logical, intent(in) :: jrot_known

         real(dp) :: r00, r003, ri, ro, rp13, rm13
         real(dp) :: w_div_wcrit_roche

         r00 = get_r_from_xh(s,k)
         ! The moment of inertia depends
         ! on the ratio of rotational frequency to its critical value.
         ! This ratio is computed in two different ways depending on whether
         ! omega or j_rot is known.
         if (jrot_known) then
            w_div_wcrit_roche = w_div_w_roche_jrot(r00,s% m(k),s% j_rot(k),s% cgrav(k), &
               s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
         else
            w_div_wcrit_roche = w_div_w_roche_omega(r00,s% m(k),s% omega(k),s% cgrav(k), &
               s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
         end if

         call eval_i_rot(s, k, r00, w_div_wcrit_roche,&
            s% i_rot(k), s% di_rot_dlnr(k), s% di_rot_dw_div_wc(k))

      end subroutine set1_irot


      ! this works like set1_xa except shifted to cell edge instead of cell center
      subroutine set1_omega(s, k, k_below_just_added, nz, &
            old_xout, new_xout, mmax, old_dmbar, new_dmbar, old_j_rot, ierr)
         ! set new value for s% omega(k)
         use num_lib, only: binary_search
         type (star_info), pointer :: s
         integer, intent(in) :: k, k_below_just_added, nz
         real(dp), intent(in) :: mmax
         real(dp), dimension(:), intent(in) :: &
            old_xout, new_xout, old_dmbar, new_dmbar, old_j_rot
         integer, intent(out) :: ierr

         real(dp) :: xm_outer, xm_inner, j_tot, xm0, xm1, new_point_dmbar, &
            dm_sum, dm
         integer :: kk, k_outer, j

         integer, parameter :: k_dbg = -1

         include 'formats'

         ierr = 0

         xm_outer = new_xout(k)
         if (k == nz) then
            new_point_dmbar = mmax - xm_outer - s% M_center
         else
            new_point_dmbar = new_dmbar(k)
         end if
         xm_inner = xm_outer + new_point_dmbar

         if (k == k_dbg) then
            write(*,2) 'xm_outer', k, xm_outer
            write(*,2) 'xm_inner', k, xm_inner
            write(*,2) 'new_point_dmbar', k, new_point_dmbar
         end if

         !write(*,*)
         !write(*,2) 'xm_outer', k, xm_outer

         dm_sum = 0d0

         if (xm_outer < old_xout(1)) then ! there is some accreted material in new
            if (xm_inner <= old_xout(1)) then
               ! new is entirely accreted material
               !write(*,2) 'new is entirely accreted material', k, new_point_dmbar
               s% omega(k) = 0
               return
            end if
            dm = min(new_point_dmbar, old_xout(1) - xm_outer)
            dm_sum = dm
            j_tot = 0
            xm_outer = old_xout(1)
            k_outer = 1
         else ! new entirely composed of old material
            if (k == k_dbg) write(*,*) 'new entirely composed of old material'
            j_tot = 0
            if (xm_outer >= old_xout(nz)) then
               ! new contained entirely in old nz
               k_outer = nz
            else
               ! binary search for k_outer such that
               ! xm_outer >= old_xout(k_outer)
               ! and old_xout(k_outer+1) > xm_outer
               k_outer = binary_search(nz, old_xout, 0, xm_outer)

               if (k == k_dbg) &
                  write(*,2) 'k_outer', k_outer, old_xout(k_outer), old_xout(k_outer+1)

               ! check
               if (k_outer <= 0 .or. k_outer > nz) then
                  ierr = -1
                  return
               end if

               if (xm_outer < old_xout(k_outer)) then
                  ierr = -1
                  return
               end if

               if (k_outer < nz) then
                  if (old_xout(k_outer+1) <= xm_outer) then
                     ierr = -1
                     return
                  end if
               end if

            end if
         end if

         if (k == -1) then
            ierr = -1
            return
         end if

         do kk = k_outer, nz ! loop until reach xm_inner
            xm0 = old_xout(kk)

            if (k == k_dbg) write(*,2) 'kk', kk, old_xout(kk), old_xout(kk+1)

            if (xm0 >= xm_inner) then
               if (dm_sum < new_point_dmbar .and. kk > 1) then
                  ! need to add a bit more from the previous source
                  dm = new_point_dmbar - dm_sum
                  dm_sum = new_point_dmbar
                  j_tot = j_tot + old_j_rot(kk-1)*dm
               end if
               exit
            end if

            if (kk == nz) then
               xm1 = mmax - s% M_center
            else
               xm1 = old_xout(kk+1)
            end if

            if (xm1 < xm_outer) then
               ierr = -1
               return
            end if

            if (xm0 >= xm_outer .and. xm1 <= xm_inner) then ! entire old kk is in new k

               dm = old_dmbar(kk)
               dm_sum = dm_sum + dm

               if (dm_sum > new_point_dmbar) then
                  ! dm too large -- numerical roundoff problems
                  dm = dm - (new_point_dmbar - dm_sum)
                  dm_sum = new_point_dmbar
               end if

               j_tot = j_tot + old_j_rot(kk)*dm

            else if (xm0 <= xm_outer .and. xm1 >= xm_inner) then ! entire new k is in old kk

               dm = new_dmbar(k)
               dm_sum = dm_sum + dm
               j_tot = j_tot + old_j_rot(kk)*dm

            else ! only use the part of old kk that is in new k

               if (k == k_dbg) then
                  write(*,*) 'only use the part of old kk that is in new k', xm_inner <= xm1
                  write(*,1) 'xm_outer', xm_outer
                  write(*,1) 'xm_inner', xm_inner
                  write(*,1) 'xm0', xm0
                  write(*,1) 'xm1', xm1
                  write(*,1) 'dm_sum', dm_sum
                  write(*,1) 'new_point_dmbar', new_point_dmbar
                  write(*,1) 'new_point_dmbar - dm_sum', new_point_dmbar - dm_sum
               end if

               if (xm_inner <= xm1) then ! this is the last part of new k

                  if (k == k_dbg) write(*,3) 'this is the last part of new k', k, kk

                  dm = new_point_dmbar - dm_sum
                  dm_sum = new_point_dmbar

               else
                  ! notice that we avoid this case if possible because of numerical roundoff

                  if (k == k_dbg) write(*,3) 'we avoid this case if possible', k, kk

                  dm = max(0d0, xm1 - xm_outer)
                  if (dm_sum + dm > new_point_dmbar) dm = new_point_dmbar - dm_sum
                  dm_sum = dm_sum + dm

               end if

               j_tot = j_tot + old_j_rot(kk)*dm

               if (dm <= 0) then
                  ierr = -1
                  return
               end if

            end if

            if (dm_sum >= new_point_dmbar) then
               if (k == k_dbg) then
                  write(*,2) 'exit for k', k
                  write(*,2) 'dm_sum', kk, dm_sum
                  write(*,2) 'new_point_dmbar', kk, new_point_dmbar
               end if
               dm_sum = new_point_dmbar
               exit
            end if

         end do

         if (dm_sum /= new_point_dmbar) then
            if (k == nz) then
               j_tot = j_tot + old_j_rot(nz)*(new_point_dmbar - dm_sum)
            else
               ierr = -1
               return
            end if
         end  if

         s% j_rot(k) = j_tot/new_point_dmbar
         call set1_irot(s, k, k_below_just_added, .true.)
         s% omega(k) = s% j_rot(k)/s% i_rot(k)

         if (k_dbg == k) then
            write(*,2) 's% omega(k)', k, s% omega(k)
            write(*,2) 's% j_rot(k)', k, s% j_rot(k)
            write(*,2) 's% i_rot(k)', k, s% i_rot(k)
            stop 'debugging: set1_omega'
         end if

      end subroutine set1_omega

      ! before mix, remove  actual_J_lost - s% angular_momentum_removed
      ! then set s% angular_momentum_removed to actual_J_lost

      subroutine adjust_J_lost(s, k_below_just_added, starting_j_rot_surf, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k_below_just_added
         real(dp), intent(in) :: starting_j_rot_surf
         integer, intent(out) :: ierr

         real(dp) :: dmm1, dm00, dm, J, mass_lost, j_for_mass_loss, &
            actual_J_lost, delta_J, frac, qlast, jnew, J_removed

         integer :: k, last_k

         include 'formats'

         ierr = 0

         mass_lost = s% mstar_old - s% mstar
         if (mass_lost <= 0) return

         ! can use a different j to account for things like wind braking
         if (s% use_other_j_for_adjust_J_lost) then
            call s% other_j_for_adjust_J_lost(s% id, starting_j_rot_surf, j_for_mass_loss, ierr)
            if (ierr /= 0) then
               write(*,*) "Error in other_j_for_adjust_J_lost"
               return
            end if
         else
            j_for_mass_loss = starting_j_rot_surf
         end if

         actual_J_lost = &
            s% adjust_J_fraction*mass_lost*j_for_mass_loss + &
            (1d0 - s% adjust_J_fraction)*s% angular_momentum_removed
         delta_J = actual_J_lost - s% angular_momentum_removed

         if (delta_J == 0d0) then
            ! this means the model is not rotating, but with rotation flag on, nothing to do
            s% adjust_J_q = 1d0
            return
         end if

         dm00 = 0
         J = 0
         last_k = s% nz
         do k = 1, s% nz
            dmm1 = dm00
            dm00 = s% dm(k)
            dm = 0.5d0*(dmm1+dm00)
            J = J + dm*s% j_rot(k)
            if (J < s% min_J_div_delta_J*abs(delta_J) &
               .or. s% q(k) > s% min_q_for_adjust_J_lost) cycle
            last_k = k
            exit
         end do
         if (last_k == s% nz) then
            ierr = 1
            write(*,*) "Error in adjust_J_lost: last_q = nz"
            return
         end if

         ! adjust in a manner that distributes the shear
         ! at each layer above last_k the specific angular momentum is scaled by a factor
         ! (frac(q-qlast)+1), where frac is computed such that the total change in angular momentum
         ! matches delta_J
         ! NOTE: this could in principle cause counter-rotation depending on the choice of
         ! parameters, but since j_rot normally rises steeply towards the surface its unlikely.
         dm00 = 0
         frac = 0
         qlast = s% q(last_k)
         s% adjust_J_q = qlast
         do k = 1, last_k-1
            dmm1 = dm00
            dm00 = s% dm(k)
            dm = 0.5d0*(dmm1+dm00)
            frac = frac + (s% q(k) - qlast)*s% j_rot(k)*dm
         end do
         frac = -delta_J/frac

         dm00 = 0
         J_removed = 0d0
         do k = 1, last_k-1
            dmm1 = dm00
            dm00 = s% dm(k)
            dm = 0.5d0*(dmm1+dm00)
            jnew = (frac*(s% q(k)-qlast)+1)*s% j_rot(k)
            J_removed = J_removed + dm*(s% j_rot(k) - jnew)
            s% j_rot(k) = jnew
            call set1_irot(s, k, k_below_just_added, .true.)
            s% omega(k) = s% j_rot(k)/s% i_rot(k)
         end do

         s% angular_momentum_removed = actual_J_lost

      end subroutine adjust_J_lost


      subroutine set_D_omega( &
            s, nz, k_const_mass, k_newval, &
            rxm_old, rxm_new, delta_m, old_xmstar, new_xmstar, &
            D_omega, oldloc, newloc, oldval, newval, work, ierr)
         use interp_1d_lib
         use interp_1d_def
         type (star_info), pointer :: s
         integer, intent(in) :: nz, k_const_mass, k_newval
         real(dp), dimension(:), intent(in) :: rxm_old, rxm_new ! (nz)
         real(dp), intent(in) :: delta_m, old_xmstar, new_xmstar
         real(dp), dimension(:) :: &
            D_omega, oldloc, newloc, oldval, newval

         real(dp), pointer :: work(:)

         integer, intent(out) :: ierr

         integer :: n, nwork, k
         logical :: dbg

         include 'formats'

         ierr = 0

         dbg = .false.
         n = nz! k_const_mass
         nwork = pm_work_size

         oldloc(1) = 0
         do k=2,n
            oldloc(k) = rxm_old(k)
         end do
         do k=1,n
            newloc(k) = rxm_new(k)
            oldval(k) = D_omega(k)
         end do
         
         call interpolate_vector( &
            n, oldloc, n, newloc, oldval, newval, interp_pm, nwork, work, &
            'adjust_mass set_D_omega', ierr)
         if (ierr /= 0) return
         do k=1,k_newval-1
            D_omega(k) = 0d0
         end do
         do k=k_newval,n
            D_omega(k) = newval(k)
         end do

      end subroutine set_D_omega


      end module adjust_mass












