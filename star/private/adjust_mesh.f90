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

      module adjust_mesh

      use star_private_def
      use const_def
      use adjust_mesh_support
      use adjust_mesh_plot

      implicit none

      private
      public :: remesh

      logical, parameter :: dbg_remesh = .false.


      contains


      ! makes new mesh and sets new values for xh and xa.
      integer function remesh(s)
         use alloc
         use mesh_functions, only: max_allowed_gvals
         use mesh_plan, only: do_mesh_plan
         use mesh_adjust, only: do_mesh_adjust
         use rates_def, only: i_rate
         use net_lib, only: clean_up_fractions
         use utils_lib
         use star_utils
         use chem_def
         use chem_lib, only: chem_get_iso_id

         type (star_info), pointer :: s

         integer :: alloc_level, i, k, ierr, species, nvar, nz, nz_new, nz_old, &
            num_gvals, j, cid, cid_max, unchanged, split, merged, revised
         type (star_info), target :: copy_info
         type (star_info), pointer :: c, prv
         real(dp) :: delta_coeff, LH, sum_L_other, sum_L_other_limit, A_max, &
            mesh_max_allowed_ratio, tmp, J_tot1, J_tot2, center_logT, alfa, beta, &
            d_dlnR00, d_dlnRp1, d_dv00, d_dvp1
         real(dp), pointer, dimension(:) :: &
            delta_gval_max, specific_PE, specific_KE, &
            xq_old, xq_new, dq_new
         real(dp), pointer :: gvals(:,:), gvals1(:)
         integer, pointer :: which_gval(:) , comes_from(:), cell_type(:)
         logical, pointer :: do_not_split(:)
         character (len=32) :: gval_names(max_allowed_gvals)
         logical, dimension(max_allowed_gvals) :: &
            gval_is_xa_function, gval_is_logT_function
         logical :: changed_mesh, dumping
         logical, parameter :: dbg = .false.

         real(dp), parameter :: max_sum_abs = 10d0
         real(dp), parameter :: xsum_tol = 1d-2
         real(dp), parameter :: h_cntr_limit = 0.5d0 ! for pre-MS decision

         include 'formats'

         ierr = 0
         remesh = keep_going
         alloc_level = 0

         s% mesh_adjust_KE_conservation = 0
         s% mesh_adjust_IE_conservation = 0
         s% mesh_adjust_PE_conservation = 0

         if (dbg_remesh) write(*,*) 'enter adjust mesh'

         if (.not. s% okay_to_remesh) then
            if (dbg_remesh) write(*,*) 's% okay_to_remesh', s% okay_to_remesh
            return
         end if

         if (s% remesh_max_allowed_logT < 1d2) then ! check it
            if (maxval(s% lnT(1:s% nz))/ln10 > s% remesh_max_allowed_logT) then
               if (dbg_remesh) write(*,2) &
                  's% remesh_max_allowed_logT', s% model_number, s% remesh_max_allowed_logT
               return
            end if
         end if

         call nullify_work_ptrs

         if (s% dt < s% remesh_dt_limit) then
            if (dbg_remesh) write(*,*) 's% dt < s% remesh_dt_limit'
            return
         end if

         if (s% L_nuc_burn_total > 0 .and. s% M_center == 0 .and. &
               s% chem_id(maxloc(s% xa(:,s% nz),dim=1)) == ih1 .and. &
               safe_log10(s% L_nuc_burn_total) < s% remesh_log_L_nuc_burn_min) then
            if (dbg_remesh) write(*,*) 'remesh_log_L_nuc_burn_min'
            return
         end if

         species = s% species
         nz_old = s% nz
         nz = nz_old

         J_tot1 = total_angular_momentum(s)
         call clean_up_fractions(1, nz, species, nz, s% xa, max_sum_abs, xsum_tol, ierr)
         if (ierr /= 0) then
            if (dbg_remesh) write(*,*) 'clean_up_fractions failed during adjust mesh'
            s% result_reason = adjust_mesh_failed
            s% retry_message = 'clean_up_fractions failed during adjust mesh'
            remesh = retry
            return
         end if

         call do_alloc1(ierr)
         if (ierr /= 0) then
            if (dbg_remesh) write(*,*) 'do_alloc1 failed for adjust mesh'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if

         do k=1,nz
            specific_PE(k) = cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
            specific_KE(k) = cell_specific_KE(s,k,d_dv00,d_dvp1)
         end do

         s% mesh_call_number = s% mesh_call_number + 1
         dumping = (s% mesh_call_number == s% mesh_dump_call_number)

         if (dumping) then
            write(*,*) 'write mesh plot info'
            write(*,*) 's% model_number', s% model_number
            write(*,*) 's% mesh_dump_call_number', s% mesh_dump_call_number
            write(*,*)
         end if

         if (.not. associated(s% other_star_info)) then
            allocate(s% other_star_info)
            prv => s% other_star_info
            c => null()
         else
            prv => s% other_star_info
            c => copy_info
            c = prv
         end if

         prv = s ! this makes copies of pointers and scalars

         sum_L_other = 0
         sum_L_other_limit = Lsun
         do j=1,num_categories
            if (j == ipp .or. j == icno .or. j == i3alf .or. j == iphoto) cycle
            sum_L_other = sum_L_other + dot_product(s% dq(1:nz), s% eps_nuc_categories(j,1:nz))
            if (sum_L_other > sum_L_other_limit) exit
         end do

         center_logT = s% lnT(s% nz)/ln10
         if (center_logT <= s% logT_max_for_standard_mesh_delta_coeff) then
            delta_coeff = s% mesh_delta_coeff
         else if (center_logT >= s% logT_min_for_highT_mesh_delta_coeff) then
            delta_coeff = s% mesh_delta_coeff_for_highT
         else
            alfa = (center_logT - s% logT_max_for_standard_mesh_delta_coeff)/ &
               (s% logT_min_for_highT_mesh_delta_coeff - s% logT_max_for_standard_mesh_delta_coeff)
            beta = 1d0 - alfa
            delta_coeff = alfa*s% mesh_delta_coeff_for_highT + beta*s% mesh_delta_coeff
         end if

         mesh_max_allowed_ratio = s% mesh_max_allowed_ratio
         if (mesh_max_allowed_ratio < 2.5d0) then
            write(*,*) 'WARNING: increasing mesh_max_allowed_ratio to 2.5'
            s% mesh_max_allowed_ratio = 2.5d0
         end if

         ! find the heaviest species in the current net
         cid_max = 0; A_max = 0

         do i = 1, s% species
            cid = s% chem_id(i)
            if (chem_isos% W(cid) > A_max) then
               A_max = chem_isos% W(cid); cid_max = cid
            end if
         end do

         do j = 1, num_xa_function
            if (s% xa_mesh_delta_coeff(j) == 1) cycle
            if (len_trim(s% xa_function_species(j)) == 0) cycle
            cid = chem_get_iso_id(s% xa_function_species(j))
            if (cid <= 0) cycle
            if (s% net_iso(cid) == 0) cycle
            if (chem_isos% W(cid) /= A_max) cycle
            delta_coeff = delta_coeff*s% xa_mesh_delta_coeff(j)
         end do

         call check_validity(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'check_validity failed at entry to adjust mesh'
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            s% result_reason = adjust_mesh_failed
            call dealloc
            return
         end if

         call get_num_gvals
         if (num_gvals == 0) then
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            write(*,*) 'must have at least 1 mesh gradient function with nonzero weight'
            s% result_reason = adjust_mesh_failed
            call dealloc
            return
         end if

         call do_alloc2(ierr)
         if (ierr /= 0) then
            if (dbg_remesh) write(*,*) 'do_alloc2 failed for adjust mesh'
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            s% result_reason = adjust_mesh_failed
            call dealloc
            return
         end if

         if (dbg_remesh) write(*,*) 'call mesh_plan'

         if (s% show_mesh_changes) then
            write(*,*) 'doing mesh_call_number', s% mesh_call_number
            write(*,*) '      mesh_plan nz_old', nz_old
         end if

         call get_gval_info( &
               s, delta_gval_max, gvals1, nz, &
               num_gvals, gval_names, &
               gval_is_xa_function, gval_is_logT_function, ierr)
         if (ierr /= 0) then
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            if (dbg_remesh) write(*,*) 'get_gval_info failed for adjust mesh'
            s% result_reason = adjust_mesh_failed
            call dealloc
            return
         end if

         do k=1,nz-1
            do_not_split(k) = &
               (s% lnR(k) - s% lnR(k+1) < 2*s% mesh_min_dlnR) .or. &
               (s% r(k) - s% r(k+1) < 2*s% mesh_min_dr_div_cs*s% csound(k))
         end do
         do_not_split(nz) = &
            (s% r(nz) - s% R_center < 2*s% mesh_min_dr_div_cs*s% csound(nz))

         if (dbg_remesh) write(*,*) 'call do_mesh_plan'

         if (dbg) write(*,1) 'sum(s% dq(1:nz))', sum(s% dq(1:nz))
         if (dbg) write(*,1) 's% dq(nz)', s% dq(nz)
         if (dbg) write(*,1) 'sum(s% dq(1:nz-1))', sum(s% dq(1:nz-1))
         tmp = s% dq(nz)
         s% dq(nz) = 1d0 - sum(s% dq(1:nz-1))
         if (s% dq(nz) <= 0) then
            if (dbg) write(*,1) 'fix starting s% dq(nz) <= 0', s% dq(nz), sum(s% dq(1:nz))
            s% dq(nz) = tmp
            s% dq(1:nz) = s% dq(1:nz)/sum(s% dq(1:nz))
            s% dq(nz) = 1d0 - sum(s% dq(1:nz-1))
            if (s% dq(nz) <= 0) then
               ierr = -1
               if (dbg) write(*,*) 's% dq(nz) <= 0'
               if (dbg) stop 'debug adjust mesh'
               s% result_reason = adjust_mesh_failed
               s% termination_code = t_adjust_mesh_failed
               remesh = terminate
               call dealloc
               return
            end if
         end if

         call set_xqs(nz, xq_old, s% dq, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'set_xqs ierr for xq_old in adjust mesh'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if

         do k=1,nz
            delta_gval_max(k) = delta_gval_max(k)*delta_coeff
         end do

         call do_mesh_plan( &
            s, nz, s% max_allowed_nz, s% mesh_ok_to_merge, s% D_mix, &
            s% mesh_max_k_old_for_split, s% mesh_min_k_old_for_split, xq_old, s% dq, &
            s% min_dq, s% max_dq*s% mesh_delta_coeff, s% min_dq_for_split, mesh_max_allowed_ratio, &
            do_not_split, num_gvals, gval_names, &
            gval_is_xa_function, gval_is_logT_function, gvals, delta_gval_max, &
            s% max_center_cell_dq*s% mesh_delta_coeff, s% max_surface_cell_dq*s% mesh_delta_coeff, &
            s% max_num_subcells, s% max_num_merge_cells, &
            nz_new, xq_new, dq_new, which_gval, comes_from, ierr)

         if (dbg_remesh .or. dbg) write(*,*) 'back from mesh_plan'

         if (ierr /= 0) then
            write(*,*)
            write(*,*) 'mesh_plan problem'
            write(*,*) 'doing mesh_call_number', s% mesh_call_number
            write(*,*) 's% model_number', s% model_number
            write(*,*)
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            s% result_reason = adjust_mesh_failed
            call dealloc
            return
         end if
         
         nz = nz_new
         s% nz = nz
         nvar = s% nvar

         if (dbg_remesh .or. dbg) write(*,*) 'call resize_star_info_arrays'
         call resize_star_info_arrays(s, c, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'resize_star_info_arrays ierr'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if

         if (dbg_remesh .or. dbg) write(*,*) 'call do_alloc3'
         call do_alloc3(ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'do_alloc3 failed for adjust mesh'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if

         if (dbg_remesh .or. dbg) write(*,*) 'call set_types_of_new_cells'
         call set_types_of_new_cells(cell_type, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'set_types_of_new_cells ierr'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if

         do k=1,nz
            s% dq(k) = dq_new(k)
         end do
         s% dq(nz) = 1d0 - sum(s% dq(1:nz-1))
         if (s% dq(nz) <= 0) then
            write(*,1) 'fix s% dq(nz) <= 0', s% dq(nz), sum(s% dq(1:nz))
            s% dq(nz) = dq_new(nz)
            s% dq(1:nz) = s% dq(1:nz)/sum(s% dq(1:nz))
            s% dq(nz) = 1d0 - sum(s% dq(1:nz-1))
            if (s% dq(nz) <= 0) then
               ierr = -1
               if (dbg) write(*,*) 's% dq(nz) <= 0 in adjust mesh'
               s% result_reason = adjust_mesh_failed
               s% termination_code = t_adjust_mesh_failed
               remesh = terminate
               call dealloc
               return
            end if
         end if

         call set_qs(s, nz, s% q, s% dq, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'set_qs ierr in adjust mesh'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if
         call set_xqs(nz, xq_new, s% dq, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'set_qs ierr in adjust mesh'
            s% result_reason = adjust_mesh_failed
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if

         ! testing -- check for q strictly decreasing
         do k = 2, nz
            if (xq_new(k) <= xq_new(k-1)) then
               write(*,3) 'bad xq_new before call do_mesh_adjust', &
                  k, nz, xq_new(k), xq_new(k-1), dq_new(k-1), xq_new(k-1) + dq_new(k-1)
               ierr = -1
               call dealloc
               return
               stop 'debug adjust mesh'
            end if
         end do

         call set_m_and_dm(s)
         call set_dm_bar(s, s% nz, s% dm, s% dm_bar)

         if (dumping) then
            call write_plot_data_for_mesh_plan( &
               s, nz_old, nz, prv% xh, prv% xa, &
               prv% lnd, prv% lnT, prv% lnPgas, prv% lnE, prv% eturb, &
               prv% D_mix, prv% mixing_type, &
               prv% dq, prv% q, xq_old, prv% q, &
               s% species, s% i_lnR, s% i_lum, s% i_v, s% i_u, comes_from, &
               num_gvals, gval_names, gvals, delta_gval_max, &
               prv% xmstar, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            write(*,*)
         end if

         if (dbg_remesh) write(*,*) 'call do_mesh_adjust'
         
         call do_mesh_adjust( &
            s, nz, nz_old, prv% xh, prv% xa, &
            prv% energy, prv% eta, prv% lnd, prv% lnPgas, &
            prv% j_rot, prv% i_rot, prv% omega, prv% D_omega, &
            prv% conv_vel, prv% lnT, prv% eturb, specific_PE, specific_KE, &
            prv% m, prv% r, prv% rho, prv% dPdr_dRhodr_info, prv% D_mix, &
            cell_type, comes_from, prv% dq, xq_old, s% xh, s% xa, s% dq, xq_new, ierr)
         if (ierr /= 0) then
            s% retry_message = 'do_mesh_adjust failed in mesh_adjust'
            if (s% report_ierr) write(*, *) s% retry_message
            s% termination_code = t_adjust_mesh_failed
            remesh = terminate
            call dealloc
            return
         end if
         if (dbg_remesh) write(*,*) 'back from do_mesh_adjust'

         ! restore prev_mesh info
         do k=1, s% prev_mesh_nz
            s% prev_mesh_xa(:,k) = prv% prev_mesh_xa(:,k)
            s% prev_mesh_xh(:,k) = prv% prev_mesh_xh(:,k)
            s% prev_mesh_j_rot(k) = prv% prev_mesh_j_rot(k)
            s% prev_mesh_omega(k) = prv% prev_mesh_omega(k)
            s% prev_mesh_dq(k) = prv% prev_mesh_dq(k)
         end do

         if (s% show_mesh_changes) then
            ! note: do_mesh_adjust can change cell_type from unchanged to revised
            ! so need to recount
            unchanged=0; split=0; merged=0; revised=0
            do k=1,nz
               select case(cell_type(k))
               case (unchanged_type)
                  unchanged = unchanged + 1
               case (split_type)
                  split = split + 1
               case (merged_type)
                  merged = merged + 1
               case (revised_type)
                  revised = revised + 1
               case default
                  write(*,3) 'bad value for cell_type(k)', k, cell_type(k)
                  stop 'adjust_mesh'
               end select
            end do
            write(*,*) '      mesh_plan nz_new', nz_new
            write(*,*) '             unchanged', unchanged
            write(*,*) '                 split', split
            write(*,*) '                merged', merged
            write(*,*) '               revised', revised
            write(*,*)

         end if

         ! testing
         do k = 2, nz
            if (xq_new(k) <= xq_new(k-1)) then
               write(*,3) 'bad xq_new after call do_mesh_adjust', k, nz, xq_new(k), xq_new(k-1)
               ierr = -1
               call dealloc
               return
               stop 'debug: adjust mesh'
            end if
         end do

         if (ierr /= 0 .and. s% report_ierr) then
            write(*,*) 'mesh_adjust problem'
            write(*,*) 'doing mesh_call_number', s% mesh_call_number
            write(*,*) 's% model_number', s% model_number
            write(*,*) 's% nz', s% nz
            write(*,*) 's% num_retries', s% num_retries
            write(*,*)
         end if

         if (remesh /= keep_going) then
            s% result_reason = adjust_mesh_failed
            s% retry_message = 'mesh_adjust failed'
            if (s% report_ierr) write(*, *) s% retry_message
            call dealloc
            return
         end if

         if (dumping) then
            call write_plot_data_for_new_mesh( &
               s, nz, nz_old, prv% xh, prv% xa, &
               prv% D_mix, prv% q, &
               s% xh, s% xa, s% dq, s% q, xq_new, s% species, s% net_iso, &
               num_gvals, gval_names, gvals, &
               which_gval, comes_from, cell_type, delta_gval_max, &
               prv% xmstar, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if
            write(*,*)
         end if

         if (dumping) call end_dump

         if (s% rotation_flag) then
            J_tot2 = total_angular_momentum(s)
            if (abs(J_tot2 - J_tot1) > 1d-2*max(J_tot1,J_tot2)) then
               s% retry_message = 'adjust_mesh failed'
               if (s% report_ierr) then
                  write(*,2) 'adjust_mesh J_tot', &
                     s% model_number, (J_tot2 - J_tot1)/J_tot2, J_tot2, J_tot1
                  write(*,*) 'failure to conserve angular momentum in adjust_mesh'
                  write(*,*)
               end if
               ierr = -1
               call dealloc
               return
               !stop 'adjust_mesh J_tot conservation error'
            end if
         end if

         call dealloc

         if (dbg_remesh .or. dbg) write(*,*) 'finished adjust mesh'

         contains

         subroutine set_types_of_new_cells(cell_type, ierr)
            integer, pointer :: cell_type(:)
            integer, intent(out) :: ierr
            integer :: k, k_old, k_old_prev, new_type

            include 'formats'
            ierr = 0
            unchanged=0; split=0; merged=0; revised=0

            do k=1,nz_new
               k_old = comes_from(k)
               new_type = -111
               if (xq_new(k) < xq_old(k_old)) then
                  write(*,*) 'xq_new(k) < xq_old(k_old)'
                  write(*,2) 'xq_new(k)', k, xq_new(k)
                  write(*,2) 'xq_old(k_old)', k_old, xq_old(k_old)
                  write(*,*) 'adjust mesh set_types_of_new_cells'
                  ierr = -1
                  return
               else if (xq_new(k) > xq_old(k_old)) then
                  new_type = split_type
               else if (k_old == nz_old) then
                  if (k == nz_new) then
                     new_type = unchanged_type
                  else
                     new_type = split_type
                  end if
               else if (k == nz_new) then
                  new_type = split_type
               else ! k_old < nz_old .and. k < nz .and. xq_new(k) == xq_old(k_old)
                  if (xq_new(k+1) == xq_old(k_old+1)) then
                     new_type = unchanged_type
                  else if (xq_new(k+1) > xq_old(k_old+1)) then
                     new_type = merged_type
                  else
                     new_type = split_type
                  end if
               end if
               cell_type(k) = new_type
               select case (new_type)
                  case (split_type)
                     split = split + 1
                  case (unchanged_type)
                     unchanged = unchanged + 1
                  case (merged_type)
                     merged = merged + 1
                  case default
                     write(*,*) 'failed to set new_type in adjust mesh set_types_of_new_cells'
                     ierr = -1
                     return
               end select
            end do

            if (unchanged + split + merged /= nz_new) then
               write(*,2) 'unchanged + split + merged', unchanged + split + merged
               write(*,2) 'nz_new', nz_new
               ierr = -1
               return
            end if

         end subroutine set_types_of_new_cells


         subroutine nullify_work_ptrs
            nullify( &
               which_gval, xq_old, xq_new, dq_new, comes_from, &
               delta_gval_max, do_not_split, gvals, cell_type)
         end subroutine nullify_work_ptrs


         subroutine do_alloc1(ierr)
            integer, intent(out) :: ierr
            ierr = 0
            call get_integer_work_array(s, comes_from, nz, nz_alloc_extra, ierr)
            if (ierr /= 0) return
            call get_integer_work_array(s, which_gval, nz, nz_alloc_extra, ierr)
            if (ierr /= 0) return
            call do_work_arrays1(.true., ierr)
            alloc_level = 1
         end subroutine do_alloc1


         subroutine do_dealloc1
            call return_integer_work_array(s, comes_from)
            call return_integer_work_array(s, which_gval)
            call do_work_arrays1(.false., ierr)
         end subroutine do_dealloc1
            
            
         subroutine do_work_arrays1(alloc_flag, ierr)
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               specific_PE, nz, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               specific_KE, nz, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               xq_old, nz, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               xq_new, nz, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               dq_new, nz, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays1


         subroutine do_alloc2(ierr)
            integer, intent(out) :: ierr
            ierr = 0
            call do_work_arrays2(.true., ierr)
            if (ierr /= 0) return
            call get_logical_work_array(s, do_not_split, nz, nz_alloc_extra, ierr)
            if (ierr /= 0) return
            gvals(1:nz,1:num_gvals) => gvals1(1:nz*num_gvals)
            alloc_level = 2
         end subroutine do_alloc2


         subroutine do_dealloc2
            call do_work_arrays2(.false., ierr)
            if (ierr /= 0) return
            call return_logical_work_array(s, do_not_split)
         end subroutine do_dealloc2
            
            
         subroutine do_work_arrays2(alloc_flag, ierr)
            logical, intent(in) :: alloc_flag
            integer, intent(out) :: ierr
            logical, parameter :: crit = .false.
            ierr = 0
            call work_array(s, alloc_flag, crit, &
               delta_gval_max, nz, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
            call work_array(s, alloc_flag, crit, &
               gvals1, nz*num_gvals, nz_alloc_extra, 'adjust_mesh', ierr)
            if (ierr /= 0) return
         end subroutine do_work_arrays2


         subroutine do_alloc3(ierr)
            integer, intent(out) :: ierr
            ierr = 0
            call get_integer_work_array(s, cell_type, nz, nz_alloc_extra, ierr)
            alloc_level = 3
         end subroutine do_alloc3


         subroutine do_dealloc3
            call return_integer_work_array(s, cell_type)
         end subroutine do_dealloc3


         subroutine dealloc
            if (alloc_level >= 3) call do_dealloc3
            if (alloc_level >= 2) call do_dealloc2
            if (alloc_level >= 1) call do_dealloc1
         end subroutine dealloc


         subroutine get_num_gvals
            use mesh_functions, only: num_mesh_functions
            num_gvals = num_mesh_functions(s)
         end subroutine get_num_gvals


         subroutine end_dump
            write(*,*) '      model_num', s% model_number
            write(*,*) '         nz_old', nz_old
            write(*,*) '         nz_new', nz_new
            write(*,*) 'finished dump_mesh'
            stop 'debugging: end_dump remesh'
         end subroutine end_dump


      end function remesh


      end module adjust_mesh


