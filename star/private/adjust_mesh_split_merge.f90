! ***********************************************************************
!
!   Copyright (C) 2016-2019  The MESA Team
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

      module adjust_mesh_split_merge

      use star_private_def
      use const_def
      use chem_def, only: ih1, ihe3, ihe4
      use utils_lib
      use auto_diff_support

      implicit none

      private
      public :: remesh_split_merge

      contains

      ! makes new mesh and sets new values for xh and xa.
      integer function remesh_split_merge(s)
         use star_utils, only: total_angular_momentum, set_dm_bar

         type (star_info), pointer :: s
         real(dp) :: old_J, new_J
         integer :: ierr

         include 'formats'
         
         if (s% RSP2_flag) then
            stop 'need to add mlt_vc and Hp_face to remesh_split_merge'
         end if

         s% amr_split_merge_has_undergone_remesh(:) = .false.

         remesh_split_merge = keep_going
         if (.not. s% okay_to_remesh) return

         if (s% rotation_flag) old_J = total_angular_momentum(s)

         call amr(s,ierr)
         if (ierr /= 0) then
            s% retry_message = 'remesh_split_merge failed'
            if (s% report_ierr) write(*, *) s% retry_message
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh_split_merge = terminate
            return
         end if

         call set_dm_bar(s, s% nz, s% dm, s% dm_bar)

         if (s% rotation_flag) then
            new_J = total_angular_momentum(s)
            if (abs((old_J-new_J)/old_J)>1d-14) then
               write(*,*) "Error in angular momemtum conservation from amr split merge"
               s% result_reason = adjust_mesh_failed
               s% termination_code = t_adjust_mesh_failed
               remesh_split_merge = terminate
               return
            end if
         end if

      end function remesh_split_merge


      subroutine amr(s,ierr)
         use chem_def, only: ih1
         use hydro_rotation, only: w_div_w_roche_jrot, update1_i_rot_from_xh
         use star_utils, only: get_r_from_xh
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: TooBig, TooSmall, MaxTooBig, MaxTooSmall, dr, minE
         real(dp) :: grad_xa(s% species), cell_time, test_dr, new_xa(s% species), &
            tau_center, r00
         integer :: iTooBig, iTooSmall, iter, k, k0, species, &
            nz, i_h1, num_split, num_merge, nz_old
         include 'formats'
         species = s% species
         nz_old = s% nz
         ierr = 0         
         num_split = 0
         num_merge = 0
         MaxTooSmall = s% split_merge_amr_MaxShort
         MaxTooBig = s% split_merge_amr_MaxLong
         k = nz_old
         tau_center = s% tau(k) + &
            s% dm(k)*s% opacity(k)/(pi4*s% rmid(k)*s% rmid(k))
         do iter = 1, s% split_merge_amr_max_iters
            call biggest_smallest(s, tau_center, TooBig, TooSmall, iTooBig, iTooSmall)
            if (mod(iter,2) == 0) then
               if (iTooSmall > 0 .and. TooSmall > MaxTooSmall) then
                  call split1
               else if (iTooBig > 0 .and. TooBig > MaxTooBig) then
                  call merge1
               else
                  exit
               end if
            else
               if (iTooBig > 0 .and. TooBig > MaxTooBig) then
                  call merge1
               else if (iTooSmall > 0 .and. TooSmall > MaxTooSmall) then
                  call split1
               else
                  exit
               end if
            end if
            if (ierr /= 0) return
         end do
         do iter = 1, 100
            call emergency_merge(s, iTooSmall)
            if (iTooSmall == 0) exit
            if (s% split_merge_amr_avoid_repeated_remesh .and. &
                  s% amr_split_merge_has_undergone_remesh(iTooSmall)) exit
            if (s% trace_split_merge_amr) &
               write(*,2) 'emergency_merge', iTooSmall, s% dq(iTooSmall)
            call do_merge(s, iTooSmall, species, new_xa, ierr)
            if (ierr /= 0) return
            num_merge = num_merge + 1
            if (s% trace_split_merge_amr) then
               call report_energies(s,'after emergency_merge')
               !write(*,*)
            end if
         end do
         do iter = 1, 100
            call emergency_split(s, iTooBig)
            if (iTooBig == 0) exit
            if (s% split_merge_amr_avoid_repeated_remesh .and. &
                  s% amr_split_merge_has_undergone_remesh(iTooBig)) exit
            if (s% trace_split_merge_amr) &
               write(*,2) 'emergency_split', iTooBig, s% dq(iTooBig)
            call do_split(s, iTooBig, species, tau_center, grad_xa, new_xa, ierr)
            if (ierr /= 0) return
            num_split = num_split + 1
            if (s% trace_split_merge_amr) then
               call report_energies(s,'after emergency_split')
               !write(*,*)
            end if
         end do
         s% mesh_call_number = s% mesh_call_number + 1
         nz = s% nz
         s% q(nz) = s% dq(nz)
         s% m(nz) = s% dm(nz) + s% m_center
         if (s% show_mesh_changes) then
            write(*,*) 'doing mesh_call_number', s% mesh_call_number
            write(*,*) '                nz_old', nz_old
            write(*,*) '                nz_new', nz
            write(*,*) '                 split', num_split
            write(*,*) '                merged', num_merge
         end if
         if (s% rotation_flag) then !update moments of inertia and omega
            do k=1, s% nz
               r00 = get_r_from_xh(s,k)
               s% w_div_w_crit_roche(k) = &
                  w_div_w_roche_jrot(r00,s% m(k),s% j_rot(k),s% cgrav(k), &
                  s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
               call update1_i_rot_from_xh(s, k)
               s% omega(k) = s% j_rot(k)/s% i_rot(k)
            end do
         end if
         if (s% model_number == -6918) stop 'amr'
         
         contains
         
         subroutine split1 ! ratio of desired/actual is too large
            include 'formats'
            if (s% trace_split_merge_amr) then
               write(*,4) 'do_merge TooSmall', &
                  iTooSmall, s% nz, s% model_number, &
                  TooSmall, MaxTooSmall, s% dq(iTooSmall)
               call report_energies(s,'before merge TooSmall')
            end if
            call do_merge(s, iTooSmall, species, new_xa, ierr)
            if (ierr /= 0) return
            num_merge = num_merge + 1
            if (s% trace_split_merge_amr) then
               call report_energies(s,'after merge')
               !write(*,*)
            end if
         end subroutine split1
         
         subroutine merge1  ! ratio of actual/desired is too large
            include 'formats'
            if (s% trace_split_merge_amr) then
               write(*,4) 'do_split TooBig', iTooBig, &
                  s% nz, s% model_number, TooBig, MaxTooBig, s% dq(iTooBig)
               call report_energies(s,'before split')
            end if
            call do_split(s, iTooBig, species, tau_center, grad_xa, new_xa, ierr)
            if (ierr /= 0) return
            num_split = num_split + 1
            if (s% trace_split_merge_amr) then
               call report_energies(s,'after split')
               !write(*,*)
            end if
         end subroutine merge1
         
      end subroutine amr


      subroutine emergency_merge(s, iTooSmall)
         type (star_info), pointer :: s
         integer, intent(out) :: iTooSmall
         integer :: k_min_dq
         include 'formats'
         k_min_dq = minloc(s% dq(1:s% nz),dim=1)
         if (s% dq(k_min_dq) < s% split_merge_amr_dq_min) then
            iTooSmall = k_min_dq
         else
            iTooSmall = 0
         end if
      end subroutine emergency_merge


      subroutine emergency_split(s, iTooBig)
         type (star_info), pointer :: s
         integer, intent(out) :: iTooBig
         integer :: k_max_dq
         include 'formats'
         k_max_dq = maxloc(s% dq(1:s% nz),dim=1)
         if (s% dq(k_max_dq) > s% split_merge_amr_dq_max) then
            iTooBig = k_max_dq
         else
            iTooBig = 0
         end if
      end subroutine emergency_split


      subroutine biggest_smallest( &
            s, tau_center, TooBig, TooSmall, iTooBig, iTooSmall)
         type (star_info), pointer :: s
         real(dp), intent(in) :: tau_center
         real(dp), intent(out) :: TooBig, TooSmall
         integer, intent(out) :: iTooBig, iTooSmall
         real(dp) :: &
            oversize_ratio, undersize_ratio, abs_du_div_cs, outer_fraction, &
            xmin, xmax, dx_actual, xR, xL, dq_min, dq_max, dx_baseline, &
            outer_dx_baseline, inner_dx_baseline, inner_outer_q, r_core_cm, &
            target_dr_core, target_dlnR_envelope, target_dlnR_core, target_dr_envelope
         logical :: hydrid_zoning, flipped_hydrid_zoning, log_zoning, logtau_zoning, &
            du_div_cs_limit_flag
         integer :: nz, nz_baseline, k, kmin, nz_r_core
         real(dp), pointer :: v(:), r_for_v(:)

         include 'formats'
         
         nz = s% nz
         hydrid_zoning = s% split_merge_amr_hybrid_zoning
         flipped_hydrid_zoning = s% split_merge_amr_flipped_hybrid_zoning
         log_zoning = s% split_merge_amr_log_zoning
         logtau_zoning = s% split_merge_amr_logtau_zoning
         nz_baseline = s% split_merge_amr_nz_baseline         
         nz_r_core = s% split_merge_amr_nz_r_core
         if (s% split_merge_amr_mesh_delta_coeff /= 1d0) then
            nz_baseline = int(dble(nz_baseline)/s% split_merge_amr_mesh_delta_coeff)
            nz_r_core = int(dble(nz_r_core)/s% split_merge_amr_mesh_delta_coeff)
         end if
         if (s% split_merge_amr_r_core_cm > 0d0) then
            r_core_cm = s% split_merge_amr_r_core_cm
         else if (s% split_merge_amr_nz_r_core_fraction > 0d0) then
            r_core_cm = s% R_center + &
               s% split_merge_amr_nz_r_core_fraction*(s% r(1) - s% R_center)
         else
            r_core_cm = s% R_center
         end if
         dq_min = s% split_merge_amr_dq_min
         dq_max = s% split_merge_amr_dq_max
         inner_outer_q = 0d0
         if (s% u_flag) then
            v => s% u
            r_for_v => s% rmid
         else if (s% v_flag) then
            v => s% v
            r_for_v => s% r
         else
            nullify(v,r_for_v)
         end if
         
         if (hydrid_zoning) then
            target_dr_core = (r_core_cm - s% R_center)/nz_r_core
            target_dlnR_envelope = &
               (s% lnR(1) - log(max(1d0,r_core_cm)))/(nz_baseline - nz_r_core)
         else if (flipped_hydrid_zoning) then
            target_dlnR_core = (log(max(1d0,r_core_cm)) - s% R_center)/(nz_baseline - nz_r_core)
            target_dr_envelope = (s% r(1) - r_core_cm)/nz_r_core
         else if (logtau_zoning) then
            k = nz
            xmin = log(tau_center)
            xmax = log(s% tau(1))
            inner_dx_baseline = (xmin - xmax)/nz_baseline ! keep it > 0
            outer_dx_baseline = inner_dx_baseline 
         else if (log_zoning) then
            xmin = log(max(1d0,s% R_center))
            xmax = s% lnR(1)
            inner_dx_baseline = (xmax - xmin)/nz_baseline
            outer_dx_baseline = inner_dx_baseline
         else
            xmin = s% R_center
            xmax = s% r(1)
            inner_dx_baseline = (xmax - xmin)/nz_baseline
            outer_dx_baseline = inner_dx_baseline
         end if
         dx_baseline = outer_dx_baseline
         
         TooBig = 0d0
         TooSmall = 0d0
         iTooBig = -1
         iTooSmall = -1
         xR = xmin ! start at center
         do k = nz, 1, -1
         
            xL = xR
            dx_baseline = inner_dx_baseline
            if (hydrid_zoning) then
               if (s% r(k) < r_core_cm) then
                  xR = s% r(k)
                  if (k == nz) then
                     xL = s% R_center
                  else
                     xL = s% r(k+1)
                  end if
                  dx_baseline = target_dr_core
               else
                  xR = log(s% r(k))
                  if (k == nz) then
                     xL = log(max(1d0,s% R_center))
                  else
                     xL = log(s% r(k+1))
                  end if
                  dx_baseline = target_dlnR_envelope
               end if
            else if (flipped_hydrid_zoning) then
               if (s% r(k) <= r_core_cm) then
                  xR = log(s% r(k))
                  if (k == nz) then
                     xL = log(max(1d0,s% R_center))
                  else
                     xL = log(s% r(k+1))
                  end if
                  dx_baseline = target_dlnR_core
               else
                  xR = s% r(k)
                  if (k == nz) then
                     xL = s% R_center
                  else
                     xL = s% r(k+1)
                  end if
                  dx_baseline = target_dr_core
               end if
            else if (logtau_zoning) then
               xR = log(s% tau(k))
            else if (log_zoning) then
               xR = log(s% r(k)) ! s% lnR(k) may not be set since making many changes
            else
               xR = s% r(k)
            end if
            
            if (s% split_merge_amr_avoid_repeated_remesh .and. &
                  (s% split_merge_amr_avoid_repeated_remesh .and. &
                     s% amr_split_merge_has_undergone_remesh(k))) cycle
            dx_actual = xR - xL
            if (logtau_zoning) dx_actual = -dx_actual ! make dx_actual > 0
            
            ! first check for cells that are too big and need to be split
            oversize_ratio = dx_actual/dx_baseline
            if (TooBig < oversize_ratio .and. s% dq(k) > 5d0*dq_min) then
               if (k < nz .or. s% split_merge_amr_okay_to_split_nz) then
                  if (k > 1 .or. s% split_merge_amr_okay_to_split_1) then
                     TooBig = oversize_ratio; iTooBig = k
                  end if
               end if
            end if
            
            ! next check for cells that are too small and need to be merged

            if (s% merge_amr_ignore_surface_cells .and. &
                  k<=s% merge_amr_k_for_ignore_surface_cells) cycle

            if (abs(dx_actual)>0d0) then
               undersize_ratio = max(dx_baseline/dx_actual, dq_min/s% dq(k))
            else
               undersize_ratio = dq_min/s% dq(k)
            end if
            
            if (s% merge_amr_max_abs_du_div_cs >= 0d0) then
               call check_merge_limits
            else if (TooSmall < undersize_ratio .and. s% dq(k) < dq_max/5d0) then
               TooSmall = undersize_ratio; iTooSmall = k
            end if
            
         end do
         
         
         contains
         
         subroutine check_merge_limits
            ! Pablo's additions to modify when merge
            ! merge_amr_max_abs_du_div_cs
            ! merge_amr_du_div_cs_limit_only_for_compression
            ! merge_amr_inhibit_at_jumps

            du_div_cs_limit_flag = .false.

            if (.not. s% merge_amr_du_div_cs_limit_only_for_compression) then
               du_div_cs_limit_flag = .true.
            else if (associated(v)) then
               if (k < nz) then
                  if (v(k+1)*pow2(r_for_v(k+1)) > v(k)*pow2(r_for_v(k))) then
                     du_div_cs_limit_flag = .true.
                  end if
               end if
               if (.not. du_div_cs_limit_flag .and. k > 1) then
                  if (v(k)*pow2(r_for_v(k)) > v(k-1)*pow2(r_for_v(k-1))) then
                     du_div_cs_limit_flag = .true.
                  end if
               end if
            end if

            if (du_div_cs_limit_flag .and. associated(v)) then
               if (k == 1) then 
                  abs_du_div_cs = abs(v(k) - v(k+1))/s% csound(k)
               else if (k == nz) then
                  abs_du_div_cs = abs(v(nz-1) - v(nz))/s% csound(nz)
               else
                  abs_du_div_cs = max(abs(v(k) - v(k+1)), &
                            abs(v(k) - v(k-1)))/s% csound(k)
               end if
            else
               abs_du_div_cs = 0d0
            end if
         
            if (du_div_cs_limit_flag) then
               if (s% merge_amr_inhibit_at_jumps) then 
                  ! reduce undersize_ratio for large jumps
                  ! i.e. large jumps inhibit merges but don't prohibit completely
                  if (abs_du_div_cs > s% merge_amr_max_abs_du_div_cs) &
                     undersize_ratio = undersize_ratio * &
                        s% merge_amr_max_abs_du_div_cs/abs_du_div_cs
                  if (TooSmall < undersize_ratio .and. s% dq(k) < dq_max/5d0) then ! switch
                     TooSmall = undersize_ratio; iTooSmall = k
                  end if
               else if (TooSmall < undersize_ratio .and. &
                        abs_du_div_cs <= s% merge_amr_max_abs_du_div_cs .and. &
                        s% dq(k) < dq_max/5d0) then
                  TooSmall = undersize_ratio; iTooSmall = k
               end if
            else
               if (TooSmall < undersize_ratio .and. s% dq(k) < dq_max/5d0) then
                  TooSmall = undersize_ratio; iTooSmall = k
               end if
            end if
         end subroutine check_merge_limits
         
         
      end subroutine biggest_smallest


      subroutine do_merge(s, i_merge, species, new_xa, ierr)
         use mesh_adjust, only: set_lnT_for_energy
         use star_utils, only: set_rmid
         type (star_info), pointer :: s
         integer, intent(in) :: i_merge, species
         real(dp), intent(inout) :: new_xa(species)
         integer, intent(out) :: ierr
         logical :: merge_center
         integer :: i, ip, i0, im, k, q, nz, qi_max, qim_max, op_err
         real(dp) :: &
            rR, rL, drR, drL, rC, rho, P, v, &
            dm, dm_i, dm_ip, m_old, star_PE0, star_PE1, &
            cell_mom, cell_ie, cell_etrb, min_IE, d_IE, d_KE, d_Esum, &
            Esum_i, KE_i, PE_i, IE_i, Etrb_i, &
            Esum_ip, KE_ip, PE_ip, IE_ip, Etrb_ip, &
            Esum, KE, PE, IE, Esum1, KE1, PE1, IE1, &
            Etot0, KEtot0, PEtot0, IEtot0, &
            Etot1, KEtot1, PEtot1, IEtot1, &
            vt_i, vt_ip, j_rot_new, j_rot_p1_new, J_old, &
            dmbar_old, dmbar_p1_old, dmbar_p2_old, &
            dmbar_new, dmbar_p1_new
         include 'formats'

         ierr = 0
         s% need_to_setvars = .true.         
         star_PE0 = get_star_PE(s)
         nz = s% nz

         s% num_hydro_merges = s% num_hydro_merges+1
         i = i_merge
         if (i > 1 .and. i < s% nz) then
            ! don't merge across change in most abundance species
            qi_max = maxloc(s% xa(1:species,i), dim=1)
            qim_max = maxloc(s% xa(1:species,i-1), dim=1)
            if (qi_max == qim_max) then ! merge with smaller neighbor
               if (i+1 == nz) then
                  drL = s% r(nz) - s% R_center
               else
                  drL = s% r(i) - s% r(i+1)
               end if
               drR = s% r(i-1) - s% r(i)
               if (drR < drL) i = i-1
            ! else i-1 has different most abundant species,
               ! so don't consider merging with it
            end if
         end if

         merge_center = (i == nz)         
         if (merge_center) i = i-1
         ip = i+1
         if (s% split_merge_amr_avoid_repeated_remesh .and. &
               (s% amr_split_merge_has_undergone_remesh(i) .or. &
                  s% amr_split_merge_has_undergone_remesh(ip))) then
            s% amr_split_merge_has_undergone_remesh(i) = .true.
            s% amr_split_merge_has_undergone_remesh(ip) = .true.
            
            return
         end if

         ! merge contents of zone ip into zone i; remove zone ip
         ! get energies for i and ip before merge
         call get_cell_energies(s, i, Esum_i, KE_i, PE_i, IE_i, Etrb_i)
         call get_cell_energies(s, ip, Esum_ip, KE_ip, PE_ip, IE_ip, Etrb_ip)

         if (s% rotation_flag) then
            ! WARNING! this is designed to conserve angular momentum, but not to explicitly conserve energy
            j_rot_new = 0d0
            j_rot_p1_new = 0d0
            if (i==1) then
               dmbar_old = 0.5d0*s% dm(i)
            else
               dmbar_old = 0.5d0*(s% dm(i)+s% dm(i-1))
            end if
            if (ip /= nz) then
               dmbar_p1_old = 0.5d0*(s% dm(ip)+s% dm(ip-1))
               if (ip+1 /= nz) then
                  dmbar_p2_old = 0.5d0*(s% dm(ip+1)+s% dm(ip))
               else
                  dmbar_p2_old = s% dm(nz) + 0.5d0*s% dm(nz-1)
               end if
               ! before merge we have (for i/=nz)
               ! dmbar_old + dmbar_p1_old = 0.5*(m(i-1)+m(i))-0.5*(m(i+1)+m(i+2))
               ! after merge we have the merged dm_bar_new
               ! dmbar_new = 0.5*(m(i-1)+m(i)) - 0.5*(m(i)+m(i+2)) = 0.5*(m(i-1)-m(i+2))
               ! where m is the old mass coordinate. From this we have dmbar_new < dmbar_old,
               ! dmbar_old + dmbar_p1_old = dmbar_new + 0.5*(m(i) - m(i+1)) = dmbar_new + 0.5*dm(i)
               ! this last expression also holds if i=nz
               dmbar_new = dmbar_old + dmbar_p1_old - 0.5d0*s% dm(i)
               ! this difference in mass goes to the dm_bar below
               dmbar_p1_new = dmbar_p2_old + 0.5d0*s% dm(i)
               ! for the merged cell we take the specific angular momentum of both merged cells
               J_old = dmbar_old*s% j_rot(i) + dmbar_p1_old*s% j_rot(ip)
               j_rot_new = J_old/(dmbar_old + dmbar_p1_old)
               ! and the j_rot of the cell downwards is adjusted to preserve angular momentum
               J_old = J_old + dmbar_p2_old*s% j_rot(ip+1)
               j_rot_p1_new = (J_old - j_rot_new*dmbar_new)/dmbar_p1_new
               ! which satisfies:
               ! J_old = j_rot_p1_new*dmbar_p1_new + j_rot_new*dmbar_new
            else
               dmbar_p1_old = s% dm(nz) + 0.5d0*s% dm(nz-1)
               ! simple case, dmbar_new is equal to dmbar_old plus dmbar_p1_old
               ! no need to adjust j_rot of another cell downwards
               dmbar_new = dmbar_old + dmbar_p1_old
               j_rot_new = (dmbar_old*s% j_rot(i) + dmbar_p1_old*s% j_rot(ip))/dmbar_new
               j_rot_p1_new = 0d0
               !write(*,*) "check center", i, ip, j_rot_new, j_rot_p1_new
            end if
         end if

         dm_i = s% dm(i)
         dm_ip = s% dm(ip)
         dm = dm_i + dm_ip
         s% dm(i) = dm
         s% dq(i) = dm/s% xmstar

         if (s% RTI_flag) &
            s% alpha_RTI(i) = (dm_i*s% alpha_RTI(i) + dm_ip*s% alpha_RTI(ip))/dm
         do q=1,species
            s% xa(q,i) = (dm_i*s% xa(q,i) + dm_ip*s% xa(q,ip))/dm
         end do

         KE = KE_i + KE_ip
         if (s% u_flag) then
            v = sqrt(KE/(0.5d0*dm))
            if (s% u(i) < 0d0) v = -v
            s% u(i) = v
         else if (s% v_flag) then
            ! there's no good solution for this.
            ! so just leave s% v(i) unchanged.
         end if

         cell_ie = IE_i + IE_ip
         s% energy(i) = cell_ie/dm
         
         if (s% RSP2_flag) then
            cell_etrb = Etrb_i + Etrb_ip
            s% w(i) = sqrt(cell_etrb/dm)
         end if

         if (s% rotation_flag) then
            s% j_rot(i) = j_rot_new
            if (ip/=nz) then
               s% j_rot(ip+1) = j_rot_p1_new
            end if
         end if

         ! shift ip+1...nz down by 1 to ip...nz-1
         do i0 = ip+1, nz
            im = i0-1
            !write(*,3) 'shift to im from i0', im, i0
            do q=1,species
               s% xa(q,im) = s% xa(q,i0)
            end do
            do q=1,s% nvar_hydro
               s% xh(q,im) = s% xh(q,i0)
            end do
            s% r(im) = s% r(i0)
            s% rmid(im) = s% rmid(i0)
            s% dm(im) = s% dm(i0)
            s% m(im) = s% m(i0)
            s% dq(im) = s% dq(i0)
            s% q(im) = s% q(i0)
            if (s% u_flag) then
               s% u(im) = s% u(i0)
            else if (s% v_flag) then
               s% v(im) = s% v(i0)
            end if
            if (s% RSP2_flag) then
               s% w(im) = s% w(i0)
               s% Hp_face(im) = s% Hp_face(i0)
            end if
            s% energy(im) = s% energy(i0)
            s% dPdr_dRhodr_info(im) = s% dPdr_dRhodr_info(i0)
            s% cgrav(im) = s% cgrav(i0)
            s% alpha_mlt(im) = s% alpha_mlt(i0)
            s% lnT(im) = s% lnT(i0)
            s% D_mix(im) = s% D_mix(i0)
            s% mlt_vc(im) = s% mlt_vc(i0)
            s% csound(im) = s% csound(i0)
            s% tau(im) = s% tau(i0)
            s% opacity(im) = s% opacity(i0)
            s% amr_split_merge_has_undergone_remesh(im) = s% amr_split_merge_has_undergone_remesh(i0)
            if (s% RTI_flag) s% alpha_RTI(im) = s% alpha_RTI(i0)
            if (s% rotation_flag) s% j_rot(im) = s% j_rot(i0)
         end do
         s% amr_split_merge_has_undergone_remesh(i) = .true.

         nz = nz - 1
         s% nz = nz
         
         if (s% u_flag) then
            s% xh(s% i_u,i) = s% u(i)
         else if (s% v_flag) then
            s% xh(s% i_v,i) = s% v(i)
         end if
         
         if (s% RTI_flag) s% xh(s% i_alpha_RTI,i) = s% alpha_RTI(i)
         
         if (s% RSP2_flag) then
            s% xh(s% i_w,i) = s% w(i)
            s% xh(s% i_Hp,i) = s% Hp_face(i)
         end if

         ! do this after move cells since need new r(ip) to calc new rho(i).
         call update_xh_eos_and_kap(s,i,species,new_xa,ierr)
         if (ierr /= 0) return ! stop 'update_xh_eos_and_kap failed in do_merge'
         
         s% rmid_start(i) = -1
         call set_rmid(s, i, i, ierr)
         if (ierr /= 0) return ! stop 'update_xh_eos_and_kap failed in do_merge'

         star_PE1 = get_star_PE(s)
         call revise_star_radius(s, star_PE0, star_PE1)

      end subroutine do_merge
      
      
      subroutine revise_star_radius(s, star_PE0, star_PE1)
         use star_utils, only: store_r_in_xh, get_lnR_from_xh
         type (star_info), pointer :: s
         real(dp), intent(in) :: star_PE0, star_PE1
         integer :: k
         real(dp) :: frac, r, star_PE, new_frac
         include 'formats'
         if (star_PE1 == 0d0 .or. star_PE0 == star_PE1) return
         frac = star_PE1/star_PE0
         if (s% model_number == -6918) write(*,1) 'frac', frac
         if (s% model_number == -6918) write(*,1) 'star_PE0', star_PE0
         if (s% model_number == -6918) write(*,1) 'star_PE1', star_PE1
         do k=1,s% nz
            s% r(k) = s% r(k)*frac
            if (s% model_number == -6918) write(*,2) 's% r(k)', k, s% r(k)
            call store_r_in_xh(s, k, s% r(k))
            s% lnR(k) = get_lnR_from_xh(s,k)
         end do
         s% r_center = frac*s% r_center
      end subroutine revise_star_radius

      
      real(dp) function get_star_PE(s) result(totPE)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: PE, rL, rC, dm, mC
         totPE = 0d0
         do k=1,s% nz
            if (k == s% nz) then
               rL = s% R_center
            else
               rL = s% r(k+1)
            end if
            rC = 0.5d0*(rL + s% r(k))
            dm = s% dm(k)
            mC = s% m(k) - 0.5d0*dm
            PE = -s% cgrav(k)*mC*dm/rC
            totPE = totPE + PE
         end do
      end function get_star_PE


      subroutine get_cell_energies(s, k, Etot, KE, PE, IE, Etrb)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: Etot, KE, PE, IE, Etrb
         real(dp) :: dm, mC, v0, v1, rL, rC
         include 'formats'
         dm = s% dm(k)
         if (s% u_flag) then
            KE = 0.5d0*dm*s% u(k)**2
         else if (s% v_flag) then
            v0 = s% v(k)
            if (k < s% nz) then
               v1 = s% v(k+1)
            else
               v1 = s% v_center
            end if
            KE = 0.25d0*dm*(v0**2 + v1**2)
         else
            KE = 0d0
         end if
         IE = s% energy(k)*dm
         if (s% RSP2_flag) then
            Etrb = pow2(s% w(k))*dm
         else
            Etrb = 0d0
         end if
         if (s% cgrav(k) == 0) then
            PE = 0
         else
            if (k == s% nz) then
               rL = s% R_center
            else
               rL = s% r(k+1)
            end if
            rC = 0.5d0*(rL + s% r(k))
            mC = s% m(k) - 0.5d0*dm
            PE = -s% cgrav(k)*mC*dm/rC
         end if
         Etot = IE + KE + PE
         if (is_bad(Etot + IE + KE + PE) .or. &
             IE <= 0 .or. KE < 0) then
            write(*,2) 'nz', s% nz
            write(*,2) 'Etot', k, Etot
            write(*,2) 'IE', k, IE
            write(*,2) 'KE', k, KE
            write(*,2) 'PE', k, PE
            stop 'get_cell_energies'
         end if
      end subroutine get_cell_energies


      subroutine split1_non_negative( &
            val, grad, dr, dV, dVR, dVL, new_valL, new_valR)
         real(dp), intent(in) :: val, grad, dr, dV, dVR, dVL
         real(dp), intent(out) :: new_valL, new_valR
         real(dp) :: xL, xR, xML, xMR, xM, f
         include 'formats'
         new_valR = val
         new_valL = val
         if (val < 1d-99 .or. grad == 0d0) return
         xL = max(0d0, min(1d0, val - grad*dr/4))
         xR = max(0d0, min(1d0, val + grad*dr/4))
         xML = xL*dVL
         xMR = xR*dVR
         xM = val*dV
         if (xM < 1d-99 .or. xML + xMR < 1d-99) then
            new_valR = val
            new_valL = val
         else
            f = xM/(xML + xMR)
            new_valR = f*xR
            new_valL = f*xL
         end if
      end subroutine split1_non_negative


      subroutine do_split(s, i_split, species, tau_center, grad_xa, new_xa, ierr)
         use alloc, only: reallocate_star_info_arrays
         use star_utils, only: set_rmid, store_r_in_xh
         type (star_info), pointer :: s
         integer, intent(in) :: i_split, species
         real(dp) :: tau_center, grad_xa(species), new_xa(species)
         integer, intent(out) :: ierr
         integer :: i, ip, j, jp, q, nz, nz_old, &
            iR, iC, iL, imin, imax, op_err
         real(dp) :: &
            cell_Esum_old, cell_KE_old, cell_PE_old, cell_IE_old, cell_Etrb_old, &
            rho_RR, rho_iR, rR, rL, dr, dr_old, rC, dV, dVR, dVL, dM, dML, dMR, rho, &
            v, v2, energy, v2_R, energy_R, rho_R, v2_C, energy_C, rho_C, v2_L, energy_L, rho_L, &
            dLeft, dRght, dCntr, grad_rho, grad_energy, grad_v2, &
            sumx, sumxp, new_xaL, new_xaR, star_PE0, star_PE1, got_cell_Esum, &
            got_cell_Esum_R, got_cell_KE_R, got_cell_PE_R, got_cell_IE_R, &
            got_cell_Esum_L, got_cell_KE_L, got_cell_PE_L, got_cell_IE_L, &
            grad_alpha, f, new_alphaL, new_alphaR, v_R, v_C, v_L, min_dm, &
            mlt_vcL, mlt_vcR, tauL, tauR, etrb, etrb_L, etrb_C, etrb_R, grad_etrb, &
            j_rot_new, dmbar_old, dmbar_p1_old, dmbar_new, dmbar_p1_new, dmbar_p2_new, J_old
         logical :: okay, done, use_new_grad_rho
         include 'formats'

         ierr = 0
         star_PE0 = get_star_PE(s)
         s% need_to_setvars = .true.         
         nz = s% nz
         s% num_hydro_splits = s% num_hydro_splits + 1
         done = .false.
         nz_old = nz

         i = i_split
         ip = i+1

         call get_cell_energies( &
            s, i, cell_Esum_old, cell_KE_old, cell_PE_old, cell_IE_old, cell_Etrb_old)

         if (s% rotation_flag .and. i<nz) then
            ! WARNING! this is designed to conserve angular momentum, but not to explicitly conserve energy
            if (i==1) then
               dmbar_old = 0.5d0*s% dm(i)
            else
               dmbar_old = 0.5d0*(s% dm(i)+s% dm(i-1))
            end if
            if (ip /= nz) then
               dmbar_p1_old = 0.5d0*(s% dm(ip)+s% dm(ip-1))
            else
               dmbar_p1_old = s% dm(ip)+0.5d0*s% dm(ip-1)
            end if

            ! We need to conserve the angular momentum in these two cells
            J_old = dmbar_old*s% j_rot(i) + dmbar_p1_old*s% j_rot(ip)
         end if

         rR = s% r(i)
         mlt_vcR = s% mlt_vc(i)
         if (i == nz) then
            rL = s% R_center
            mlt_vcL = 0d0
            tauL = tau_center
         else
            rL = s% r(ip)
            mlt_vcL = s% mlt_vc(ip)
            tauL = s% tau(ip)
         end if
         
         tauR = s% tau(i)
         if (i == nz) then
            tauL = tau_center
         else
            tauL = s% tau(ip)
         end if
         if (is_bad(tauL+tauR)) then
            !$omp critical (adjust_mesh_split_merge_crit1)
            write(*,2) 'tauL', ip, tauL
            write(*,2) 'tauR', i, tauR
            write(*,2) 'nz', nz
            stop 'do_split'
            !$omp end critical (adjust_mesh_split_merge_crit1)         
         end if
                 
         dr = rR - rL
         dr_old = dr
         rC = 0.5d0*(rR + rL) ! split at center by radius

         dV = four_thirds_pi*(rR*rR*rR - rL*rL*rL)
         dVR = four_thirds_pi*(rR*rR*rR - rC*rC*rC)
         dVL = dV - dVR

         dm = s% dm(i)
         rho = dm/dV
         
         if (s% u_flag) then
            v = s% u(i)
            v2 = v*v
         else ! not used
            v = 0
            v2 = 0
         end if
         
         energy = s% energy(i)
         if (s% RSP2_flag) etrb = pow2(s% w(i))

         ! use iR, iC, and iL for getting values to determine slopes
         if (i > 1 .and. i < nz_old) then
            iR = i-1
            iC = i
            iL = i+1
         else if (i == 1) then
            iR = 1
            iC = 2
            iL = 3
         else ! i == nz_old
            iR = nz_old-2
            iC = nz_old-1
            iL = nz_old
         end if
         
         energy_R = s% energy(iR)
         rho_R = s% dm(iR)/get_dV(s,iR)
         
         energy_C = s% energy(iC)
         rho_C = s% dm(iC)/get_dV(s,iC)
         
         energy_L = s% energy(iL)
         rho_L = s% dm(iL)/get_dV(s,iL)
            
         ! get gradients before move cell contents

         if (iL == nz_old) then
            dLeft = 0.5d0*(s% r(iC) - s% R_center)
         else
            dLeft = 0.5d0*(s% r(iC) - s% r(iL+1))
         end if
         dRght = 0.5d0*(s% r(iR) - s% r(iL))
         dCntr = dLeft + dRght

         if (s% equal_split_density_amr) then ! same density in both parts
            rho_RR = 0
            grad_rho = 0d0
            use_new_grad_rho = .false.
         else if (.false.) then
            rho_RR = s% dm(iR-1)/get_dV(s,iR-1)
            grad_rho = 2*(rho_RR - rho_R)/(s% r(iR-1) - s% r(iR+1))
            use_new_grad_rho = .true.
         else
            rho_RR = 0
            grad_rho = get1_grad(rho_L, rho_C, rho_R, dLeft, dCntr, dRght)
            use_new_grad_rho = .false.
         end if

         grad_energy = get1_grad(energy_L, energy_C, energy_R, dLeft, dCntr, dRght)
            
         if (s% RTI_flag) grad_alpha = get1_grad( &
            s% alpha_RTI(iL), s% alpha_RTI(iC), s% alpha_RTI(iR), dLeft, dCntr, dRght)
         
         if (s% RSP2_flag) then
            etrb_R = pow2(s% w(iR))
            etrb_C = pow2(s% w(iC))
            etrb_L = pow2(s% w(iL))
            grad_etrb = get1_grad(etrb_L, etrb_C, etrb_R, dLeft, dCntr, dRght)
         end if
         
         if (s% u_flag) then
            v_R = s% u(iR)
            v2_R = v_R*v_R
            v_C = s% u(iC)
            v2_C = v_C*v_C         
            v_L = s% u(iL)
            v2_L = v_L*v_L
            if ((v_L - v_C)*(v_C - v_R) <= 0) then ! not strictly monotonic velocities
               grad_v2 = 0d0
            else
               grad_v2 = get1_grad(v2_L, v2_C, v2_R, dLeft, dCntr, dRght)
            end if
         else if (s% v_flag) then
            if (iL == s% nz) then
               v_L = s% v_center
            else
               v_L = s% v(ip)
            end if
            v2_L = v_L*v_L
            v_R = s% v(i)
            v2_R = v_R*v_R
         end if

         do q = 1, species
            grad_xa(q) = get1_grad( &
               s% xa(q,iL), s% xa(q,iC), s% xa(q,iR), dLeft, dCntr, dRght)
         end do

         nz = nz + 1
         s% nz = nz
         call reallocate_star_info_arrays(s, ierr)
         if (ierr /= 0) then
            write(*,2) 'reallocate_star_info_arrays ierr', ierr
            stop 'split failed'
         end if

         if (i_split < nz_old) then ! move i_split..nz-1 to i_split+1..nz
            do j=nz-1,i_split,-1
               jp = j+1
               do q=1,species
                  s% xa(q,jp) = s% xa(q,j)
               end do
               do q=1,s% nvar_hydro
                  s% xh(q,jp) = s% xh(q,j)
               end do
               s% r(jp) = s% r(j)
               s% rmid(jp) = s% rmid(j)
               s% dm(jp) = s% dm(j)
               s% m(jp) = s% m(j)
               s% dq(jp) = s% dq(j)
               s% q(jp) = s% q(j)
               if (s% u_flag) then
                  s% u(jp) = s% u(j)
               else if (s% v_flag) then
                  s% v(jp) = s% v(j)
               end if
               s% energy(jp) = s% energy(j)
               s% dPdr_dRhodr_info(jp) = s% dPdr_dRhodr_info(j)
               s% cgrav(jp) = s% cgrav(j)
               s% alpha_mlt(jp) = s% alpha_mlt(j)
               s% lnT(jp) = s% lnT(j)
               s% D_mix(jp) = s% D_mix(j)
               s% mlt_vc(jp) = s% mlt_vc(j)
               s% csound(jp) = s% csound(j)
               s% tau(jp) = s% tau(j)
               s% opacity(jp) = s% opacity(j)
               s% amr_split_merge_has_undergone_remesh(jp) = s% amr_split_merge_has_undergone_remesh(j)
               if (s% RTI_flag) s% alpha_RTI(jp) = s% alpha_RTI(j)
               if (s% rotation_flag) s% j_rot(jp) = s% j_rot(j)
            end do
         end if
         s% amr_split_merge_has_undergone_remesh(i) = .true.
         s% amr_split_merge_has_undergone_remesh(ip) = .true.
         
         dM = rho*dV
         if (.not. use_new_grad_rho) then ! do it the old way

            rho_L = rho - grad_rho*dr/4
            rho_R = rho + grad_rho*dr/4
            if (grad_rho == 0d0 .or. &
                  rho_L <= 0d0 .or. &
                  rho_R <= 0d0) then
               rho_R = rho
               rho_L = rho
               dML = rho_L*dVL
               dMR = dM - dML ! should = rhoR*dVR
            else
               dML = rho_L*dVL
               dMR = rho_R*dVR
               f = dM/(dML + dMR)
               rho_L = f*rho_L
               rho_R = f*rho_R
               dML = rho_L*dVL
               dMR = rho_R*dVR
               if (abs(dML + dMR - dM) > 1d-14*dM) then
                  write(*,2) '(dML + dMR - dM)/dM', i, (dML + dMR - dM)/dM
                  stop 'split'
               end if
               dMR = dM - dML
            end if
            
         else
            
            ! at this point, rho_R is the density of the cell iR
            ! we are about to redefine it as the density of the right 1/2 of the split
            ! similarly for rho_L 
            rho_iR = rho_R
            dR = -(dRght/4 + (s% r(iR) - s% r(iC))/2)
            rho_R = rho_iR + grad_rho*dR
            dMR = rho_R*dVR
            dML = dM - dMR
            rho_L = dML/dVL
            if (rho_R <= 1d-16 .or. rho_L <= 1d-16) then
   !$omp critical (adjust_mesh_split_merge_crit2)
               write(*,*)
               write(*,2) 'nz', nz
               write(*,2) 'rho_RR', iR-1, rho_RR
               write(*,2) 'rho_iR', iR, rho_iR
               write(*,2) 'rho', iC, rho
               write(*,2) 's% r(iR-1)', iR-1, s% r(iR-1)
               write(*,2) 's% r(iR)', iR, s% r(iR)
               write(*,2) 's% r(iR+1)', iR+1, s% r(iR+1)
               write(*,2) 'rho_RR', iR-1, rho_RR
               write(*,2) 'rho_iR', iR, rho_iR
               write(*,2) 'rho_RR - rho_iR', iR, rho_RR - rho_iR
               write(*,2) 'dr for right', iR, (s% r(iR-1) - s% r(iR+1))/2
               write(*,2) 'grad_rho', iR, grad_rho
               write(*,*)
               write(*,2) 's% r(iL)', iL, s% r(iL)
               write(*,2) 'dR', iC, dR
               write(*,2) 'dRho', iC, grad_rho*dR
               write(*,2) 'rho_R', iC, rho_R
               write(*,2) 'rho_L', iC, rho_L
               write(*,*)
               stop 'failed in do_split extrapolation of density from above'
   !$omp end critical  (adjust_mesh_split_merge_crit2)        
            end if
         
         end if
         
         min_dm = s% split_merge_amr_dq_min*s% xmstar
         if (dML < min_dm .or. dMR < min_dm) then
            rho_R = rho
            rho_L = rho
            dM = rho*dV
            dML = rho_L*dVL
            dMR = dM - dML ! should = rhoR*dVR
         end if

         energy_R = energy + grad_energy*dr/4
         energy_L = (dm*energy - dmR*energy_R)/dmL
         if (energy_R < 0d0 .or. energy_L < 0d0) then
            energy_R = energy
            energy_L = energy
         end if
         
         s% energy(i) = energy_R
         s% energy(ip) = energy_L
         
         if (s% RSP2_flag) then
            etrb_R = etrb + grad_etrb*dr/4
            etrb_L = (dm*etrb - dmR*etrb_R)/dmL
            if (etrb_R < 0d0 .or. etrb_L < 0d0) then
               etrb_R = etrb
               etrb_L = etrb
            end if
            s% w(i) = sqrt(max(0d0,etrb_R))
            s% w(ip) = sqrt(max(0d0,etrb_L))
         end if

         if (s% u_flag) then
            v2_R = v2 + grad_v2*dr/4
            v2_L = (dm*v2 - dmR*v2_R)/dmL
            if (v2_R < 0d0 .or. v2_L < 0d0) then
               v2_R = v2
               v2_L = v2
            end if
            s% u(i) = sqrt(v2_R)
            s% u(ip) = sqrt(v2_L)
            if (v < 0d0) then
               s% u(i) = -s% u(i)
               s% u(ip) = -s% u(ip)
            end if
         else if (s% v_flag) then ! just make a rough approximation. 
            s% v(ip) = sqrt(0.5d0*(v2_L + v2_R))
         end if
         
         if (s% RTI_flag) then ! set new alpha
            if (i == 1) then
               s% alpha_RTI(ip) = s% alpha_RTI(i)
            else
               call split1_non_negative( &
                  s% alpha_RTI(i), grad_alpha, &
                  dr, dV, dVR, dVL, new_alphaL, new_alphaR)
               s% alpha_RTI(i) = new_alphaR
               s% alpha_RTI(ip) = new_alphaL
            end if
            s% dPdr_dRhodr_info(ip) = s% dPdr_dRhodr_info(i)
         end if
         
         if (i == 1) then
            s% mlt_vc(ip) = s% mlt_vc(i)
         else
            s% mlt_vc(ip) = (mlt_vcL*dML + mlt_vcR*dMR)/dM
         end if
         
         s% tau(ip) = tauR + (tauL - tauR)*dMR/dM
         if (is_bad(s% tau(ip))) then
            write(*,2) 'tau', ip, s% tau(ip), tauL, tauR, dMR/dM
            stop 'split'
         end if

         if (i == 1) then
            do q=1,species
               s% xa(q,ip) = s% xa(q,i)
            end do
         else ! split mass fractions
            ! check that sum of grads for mass fractions is 0
            sumx = sum(grad_xa)
            if (sumx > 0) then ! reduce largest grad
               j = maxloc(grad_xa, dim=1)
            else ! increase smallest grad
               j = minloc(grad_xa, dim=1)
            end if
            grad_xa(j) = 0
            grad_xa(j) = -sum(grad_xa)
            ! set new mass fractions
            do q = 1, species
               call split1_non_negative( &
                  s% xa(q,i), grad_xa(q), &
                  dr, dV, dVR, dVL, new_xaL, new_xaR)
               s% xa(q,i) = new_xaR
               s% xa(q,ip) = new_xaL
            end do
            !check mass fractions >= 0 and <= 1 and sum to 1.0
            do q = 1, species
               s% xa(q,i) = min(1d0, max(0d0, s% xa(q,i)))
               s% xa(q,ip) = min(1d0, max(0d0, s% xa(q,ip)))
            end do
            sumx = sum(s% xa(1:species,i))
            sumxp = sum(s% xa(1:species,ip))
            do q = 1, species
               s% xa(q,i) = s% xa(q,i)/sumx
               s% xa(q,ip) = s% xa(q,ip)/sumxp
            end do
         end if
         
         if (s% u_flag) s% u_face_ad(ip)%val = 0.5d0*(s% u(i) + s% u(ip))
            ! just for setting u_face_start so don't need partials

         ! r, q, m, u_face unchanged for face i
         dVR = dV - dVL
         dmR = s% dm(i) - dmL
         s% dm(i) = dmR
         s% dq(i) = s% dm(i)/s% xmstar

         s% r(ip) = rC
         s% m(ip) = s% m(i) - s% dm(i) ! <<< using new value dm(i)
         s% q(ip) = s% q(i) - s% dq(i) ! <<< using new value dq(i)
         s% dm(ip) = dmL
         s% dq(ip) = s% dm(ip)/s% xmstar

         s% cgrav(ip) = s% cgrav(i)
         s% alpha_mlt(ip) = s% alpha_mlt(i)
         s% lnT(ip) = s% lnT(i)
         s% T(ip) = s% T(i)
         s% D_mix(ip) = s% D_mix(i)
         s% mlt_vc(ip) = s% mlt_vc(i)
         s% csound(ip) = s% csound(i)
         s% opacity(ip) = s% opacity(i)
         if (s% RTI_flag) s% alpha_RTI(ip) = s% alpha_RTI(i)

         if (s% rotation_flag) then
            ! WARNING! this is designed to conserve angular momentum, but not to explicitly conserve energy
            j_rot_new = 0d0 ! specific angular momentum for the newly created cell
            if (i==nz_old) then
               ! simple case, dm_bar contained in cells i and i+1 after merge is conserved,
               ! so use same j_rot
               j_rot_new = s% j_rot(i)
            else
               if (i==1) then
                  dmbar_new = 0.5d0*s% dm(i)
               else
                  dmbar_new = 0.5d0*(s% dm(i)+s% dm(i-1))
               end if
               if (ip /= nz) then
                  dmbar_p1_new = 0.5d0*(s% dm(ip)+s% dm(ip-1))
               else
                  dmbar_p1_new = s% dm(ip)+0.5d0*s% dm(ip-1)
               end if
               if (ip+1 /= nz) then
                  dmbar_p2_new = 0.5d0*(s% dm(ip+1)+s% dm(ip))
               else
                  dmbar_p2_new = s% dm(ip+1)+0.5d0*s% dm(ip)
               end if

               ! we keep the same j_rot for cells i and ip (cells i and ip+1 after merge),
               ! we compute the j_rot needed for the new cell to preserve angular momentum.
               j_rot_new = (J_old - dmbar_new*s% j_rot(i) - dmbar_p2_new*s% j_rot(ip+1))/(dmbar_p1_new)
               ! this satisfies:
               ! dmbar_old*s% j_rot(i) + dmbar_p1_old*s% j_rot(ip) = dmbar_new*s% j_rot(i) + dmbar_p2_new*s% j_rot(ip) + dmbar_p1_new*j_rot_new
            end if
            s% j_rot(ip) = j_rot_new
         end if

         if (s% i_lum /= 0) then
            if (ip < nz_old) then
               s% xh(s% i_lum,ip) = &
                  0.5d0*(s% xh(s% i_lum,i) + s% xh(s% i_lum,ip+1))
            else
               s% xh(s% i_lum,ip) = &
                  0.5d0*(s% xh(s% i_lum,i) + s% L_center)
            end if
         end if
         
         call store_r_in_xh(s, ip, s% r(ip))
         if (s% u_flag) then
            s% xh(s% i_u,i) = s% u(i)
            s% xh(s% i_u,ip) = s% u(ip)
         else if (s% v_flag) then
            s% xh(s% i_v,i) = s% v(i)
            s% xh(s% i_v,ip) = s% v(ip)
         end if
         if (s% RTI_flag) then
            s% xh(s% i_alpha_RTI,i) = s% alpha_RTI(i)
            s% xh(s% i_alpha_RTI,ip) = s% alpha_RTI(ip)
         end if

         call update_xh_eos_and_kap(s,i,species,new_xa,ierr)
         if (ierr /= 0) return ! stop 'update_xh_eos_and_kap failed in do_split'

         call update_xh_eos_and_kap(s,ip,species,new_xa,ierr)
         if (ierr /= 0) return ! stop 'update_xh_eos_and_kap failed in do_split'
         
         s% rmid_start(i) = -1
         s% rmid_start(ip) = -1
         call set_rmid(s, i, ip, ierr)
         if (ierr /= 0) return ! stop 'update_xh_eos_and_kap failed in do_split'

         star_PE1 = get_star_PE(s)
         call revise_star_radius(s, star_PE0, star_PE1)

      end subroutine do_split


      subroutine update_xh_eos_and_kap(s,i,species,new_xa,ierr)
         use mesh_adjust, only: set_lnT_for_energy
         use micro, only: do_kap_for_cell
         use eos_lib, only: eos_gamma_DE_get_PT
         use chem_lib, only: basic_composition_info
         use star_utils, only: store_lnT_in_xh, get_T_and_lnT_from_xh, &
            store_rho_in_xh, get_rho_and_lnd_from_xh
         type (star_info), pointer :: s
         integer, intent(in) :: i, species
         real(dp) :: new_xa(species)
         integer, intent(out) :: ierr
         real(dp) :: rho, logRho, new_lnT, revised_energy, xsum
         integer :: q
         include 'formats'
         ierr = 0
         rho = s% dm(i)/get_dV(s,i)
         call store_rho_in_xh(s, i, rho)
         call get_rho_and_lnd_from_xh(s, i, s% rho(i), s% lnd(i))
         logRho = s% lnd(i)/ln10
         do q=1,species
            new_xa(q) = s% xa(q,i)
         end do
         call set_lnT_for_energy(s, i, &
            s% net_iso(ih1), s% net_iso(ihe3), s% net_iso(ihe4), &
            species, new_xa, rho, logRho, s% energy(i), s% lnT(i), &
            new_lnT, revised_energy, ierr)
         if (ierr /= 0) return
         call store_lnT_in_xh(s, i, new_lnT)
         call get_T_and_lnT_from_xh(s, i, s% T(i), s% lnT(i))
      end subroutine update_xh_eos_and_kap


      real(dp) function get_dV(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         if (k == s% nz) then
            get_dV = four_thirds_pi*(s% r(k)*s% r(k)*s% r(k) - s% R_center*s% R_center*s% R_center)
         else
            get_dV = four_thirds_pi*(s% r(k)*s% r(k)*s% r(k) - s% r(k+1)*s% r(k+1)*s% r(k+1))
         end if
      end function get_dV


      real(dp) function minmod(a, b, c) result(m)
         real(dp), intent(in) :: a, b, c
         m = a
         if (a*b < 0d0) then
            m = 0d0
            return
         end if
         if (abs(b) < abs(m)) m = b
         if (b*c < 0d0) then
            m = 0d0
            return
         end if
         if (abs(c) < abs(m)) m = c
      end function minmod


      real(dp) function get1_grad(vL, vC, vR, dLeft, dCntr, dRght) &
            result(grad)
         real(dp), intent(in) :: vL, vC, vR, dLeft, dCntr, dRght
         real(dp) :: sL, sR, sC
         include 'formats'
         if (dLeft <= 0d0 .or. dCntr <= 0d0 .or. dRght <= 0d0) then
            grad = 0d0
            return
         end if
         sL = (vC - vL)/dLeft
         sR = (vR - vC)/dRght
         sC = (vR - vL)/dCntr
         grad = minmod(sL, sC, sR)
      end function get1_grad


      real(dp) function total_KE(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: v0, v1
         include 'formats'
         total_KE = 0
         if (s% u_flag) then
            do k=1,s% nz
               total_KE = total_KE + 0.5d0*s% dm(k)*s% u(k)**2
            end do
         else if (s% v_flag) then
            do k=1,s% nz-1
               total_KE = total_KE + &
                  0.25d0*s% dm(k)*(s% v(k)**2 + s% v(k+1)**2)
            end do
            k = s% nz
            total_KE = total_KE + &
               0.25d0*s% dm(k)*(s% v(k)**2 + s% v_center**2)
         end if
      end function total_KE


      real(dp) function total_PE(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: rL, rR, rC, mC
         total_PE = 0d0
         rR = s% R_center
         do k = s% nz,1,-1
            rL = rR
            rR = s% r(k)
            rC = 0.5d0*(rL + rR)
            mC = s% m(k) - 0.5d0*s% dm(k)
            total_PE = total_PE - s% cgrav(k)*mC*s% dm(k)/rC
         end do
      end function total_PE


      real(dp) function total_IE(s)
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: specific_ie, egas
         total_IE = 0
         do k=1,s% nz
            specific_ie = s% energy(k)
            total_IE = total_IE + specific_ie*s% dm(k)
         end do
      end function total_IE


      subroutine report_energies(s, str)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: str
         real(dp) :: KE, IE, PE, Etot
         include 'formats'
         
         return
         
         
         KE = total_KE(s)
         IE = total_IE(s)
         PE = total_PE(s)
         Etot = KE + IE + PE
         write(*,1) 'Etot, KE, IE, PE ' // trim(str), Etot, KE, IE, PE
      end subroutine report_energies


      end module adjust_mesh_split_merge


