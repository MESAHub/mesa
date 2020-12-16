! ***********************************************************************
!
!   Copyright (C) 2015-2019  Bill Paxton & The MESA Team
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

      module remove_shells

      use star_private_def
      use const_def
      use utils_lib

      implicit none


      contains


      subroutine do_remove_center_at_cell_k(id, k, ierr)
         integer, intent(in) :: id, k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_at_cell_k: get_star_ptr ierr', ierr
            return
         end if
         call do_remove_inner_fraction_q(id, s% q(k), ierr)
      end subroutine do_remove_center_at_cell_k


      subroutine do_remove_center_by_temperature(id, temperature, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: temperature
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_temperature: get_star_ptr ierr', ierr
            return
         end if
         do k=1,s% nz
            if (s% T(k) >= temperature) then
               call do_remove_inner_fraction_q(id, s% q(k), ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_center_by_temperature


      subroutine do_remove_center_by_he4(id, x, ierr)
         use chem_def, only: ihe4
         integer, intent(in) :: id
         real(dp), intent(in) :: x
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, he4
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_he4: get_star_ptr ierr', ierr
            return
         end if
         he4 = s% net_iso(ihe4)
         if (he4 <= 0) then
            ierr = -1
            write(*,*) 'do_remove_center_by_he4: no he4 in current net'
            return
         end if
         do k=1,s% nz
            if (s% xa(he4,k) >= x) then
               call do_remove_inner_fraction_q(id, s% q(k), ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_center_by_he4


      subroutine do_remove_center_by_si28(id, x, ierr)
         use chem_def, only: isi28
         integer, intent(in) :: id
         real(dp), intent(in) :: x
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, si28
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_si28: get_star_ptr ierr', ierr
            return
         end if
         si28 = s% net_iso(isi28)
         if (si28 <= 0) then
            ierr = -1
            write(*,*) 'do_remove_center_by_si28: no si28 in current net'
            return
         end if
         do k=1,s% nz
            if (s% xa(si28,k) >= x) then
               call do_remove_inner_fraction_q(id, s% q(k), ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_center_by_si28


      subroutine do_remove_center_to_reduce_co56_ni56(id, x, ierr)
         use chem_def, only: ico56, ini56
         integer, intent(in) :: id
         real(dp), intent(in) :: x
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, jj, k, species, nz, co56, ni56
         real(dp) :: mtotal, dm56, alfa_co56, dm56_new, dm56_old, x56_new, x56_old, xsum
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_to_reduce_co56_ni56: get_star_ptr ierr', ierr
            return
         end if
         nz = s% nz
         species = s% species
         co56 = s% net_iso(ico56)
         ni56 = s% net_iso(ini56)
         if (co56 <= 0 .or. ni56 <= 0) then
            ierr = -1
            write(*,*) 'do_remove_center_to_reduce_co56_ni56: must have both in current net'
            return
         end if
         mtotal = 0d0
         do k=1,nz
            dm56 = s% dm(k)*(s% xa(co56,k) + s% xa(ni56,k))/Msun
            if (mtotal+dm56 >= x) then
               write(*,2) 'mtotal+dm56 >= x', k, mtotal+dm56, x, mtotal
               if (mtotal+dm56 - x < x - mtotal) then
                  if (k < nz) call do_remove_inner_fraction_q(id, s% q(k+1), ierr)
               else
                  call do_remove_inner_fraction_q(id, s% q(k), ierr)
               end if
               if (ierr == 0) then ! adjust ni+co in center zone
                  nz = s% nz
                  mtotal = dot_product(s% dm(1:nz), &
                     s% xa(co56,1:nz) + s% xa(ni56,1:nz))/Msun
                  write(*,1) 'mtotal after remove', mtotal
                  write(*,2) 'nz after remove', nz
                  dm56 = x - mtotal ! change Ni+Co in nz by this much
                  x56_old = s% xa(co56,nz) + s% xa(ni56,nz)
                  dm56_old = s% dm(nz)*x56_old/Msun
                  if (x56_old <= 0d0) then
                     alfa_co56 = 0d0
                  else
                     alfa_co56 = s% xa(co56,nz)/x56_old
                  end if
                  dm56_new = dm56_old + dm56
                  x56_new = dm56_new/(s% dm(nz)/Msun)
                  write(*,3) 'old x co56', co56, species, s% xa(co56,nz)
                  write(*,3) 'old x ni56', ni56, species, s% xa(ni56,nz)

                  s% xa(co56,nz) = 0d0
                  s% xa(ni56,nz) = 0d0
                  j = 1
                  do jj=2,species
                     if (s% xa(jj,nz) > s% xa(j,nz)) j = jj
                  end do
                  
                  s% xa(co56,nz) = x56_new*alfa_co56
                  s% xa(ni56,nz) = x56_new*(1d0 - alfa_co56)
                  s% xa(j,nz) = s% xa(j,nz) - (x56_new - x56_old)
                  
                  write(*,1) 'new s% xa(co56,nz)', s% xa(co56,nz)
                  write(*,1) 'new s% xa(ni56,nz)', s% xa(ni56,nz)
                  mtotal = dot_product(s% dm(1:nz), &
                     s% xa(co56,1:nz) + s% xa(ni56,1:nz))/Msun
                  write(*,1) 'final mass for Ni56+Co56', mtotal
                  write(*,*)
                  stop 'do_remove_center_to_reduce_co56_ni56'
               end if
               return
            end if
            mtotal = mtotal + dm56
            write(*,2) 'mtotal', k, mtotal
         end do
         write(*,1) 'failed in do_remove_center_to_reduce_co56_ni56'
         write(*,1) 'not enough Ni56+Co56 in model', mtotal, x
         ierr = -1
      end subroutine do_remove_center_to_reduce_co56_ni56


      subroutine do_remove_fallback(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k0, k1, nz
         real(dp) :: ie, ke, pe, rR, rL, rC, m_cntr, &
            sum_total_energy, speed_limit
         
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_fallback: get_star_ptr ierr', ierr
            return
         end if
         
         if (.not. s% u_flag) then
            write(*,*) 'only remove fallback with u_flag'
            ierr = -1
            return
         end if
         
         nz = s% nz
         
         ! check to see how far extend fallback above innermost cell
         k0 = nz
         if (s% job% fallback_check_total_energy) then ! remove_bound_inner_region
            ! integrate total energy outward looking for sign going negative.
            ! if find, then continue until reach minimum integral and cut there.
            sum_total_energy = 0d0
            do k = nz,1,-1
               ie = s% energy(k)*s% dm(k)
               ke = 0.5d0*s% u(k)*s% u(k)*s% dm(k)
               if (k == s% nz) then   
                  rL = s% R_center
               else
                  rL = s% r(k+1)
               end if
               rR = s% r(k)
               rC = 0.5d0*(rR + rL)
               m_cntr = s% m(k) - 0.5d0*s% dm(k)
               pe = -s% cgrav(k)*m_cntr*s% dm(k)/rC
               if (is_bad(ie + ke + pe)) then
                  write(*,2) 'ie', k, ie
                  write(*,2) 'ke', k, ke
                  write(*,2) 'pe', k, pe
                  stop 'do_remove_fallback'
               end if
               sum_total_energy = sum_total_energy + ie + ke + pe
               if (is_bad(sum_total_energy)) then
                  write(*,2) 'sum_total_energy', k, sum_total_energy
                  write(*,2) 'ie', k, ie
                  write(*,2) 'ke', k, ke
                  write(*,2) 'pe', k, pe
                  stop 'do_remove_fallback'
               end if
               !write(*,3) 'sum_total_energy', k, nz, sum_total_energy, ie + ke + pe
               if (sum_total_energy < 0d0) exit
               k0 = k
            end do
            if (sum_total_energy >= 0d0) then
               !write(*,1) 'no bound inner region', sum_total_energy
               return ! no bound inner region
            end if
            do k=k0-1,1,-1
               ie = s% energy(k)*s% dm(k)
               ke = 0.5d0*s% u(k)*s% u(k)*s% dm(k)
               rL = s% r(k+1)
               rR = s% r(k)
               rC = 0.5d0*(rR + rL)
               m_cntr = s% m(k) - 0.5d0*s% dm(k)
               pe = -s% cgrav(k)*m_cntr*s% dm(k)/rC
               if (ie + ke + pe > 0d0) then ! now back to unbound cell
                  k0 = k+1
                  !write(*,2) 'top', k0, ie + ke + pe
                  exit
               end if
            end do
         else
            speed_limit = s% job% remove_fallback_speed_limit
            !write(*,3) '-u/limit', s% model_number, nz, -s% u(nz)/(speed_limit*s% csound(nz))
            if (-s% u(nz) <= speed_limit*s% csound(nz)) return
            do k = nz-1,1,-1
               k0 = k+1
               !write(*,3) '-u/limit', s% model_number, k, -s% u(k)/(speed_limit*s% csound(k))
               if (-s% u(k) < speed_limit*s% csound(k)) exit
            end do
         end if
         
         !stop 'do_remove_fallback'
         
         !write(*,3) 'k0 old nz', k0, s% nz, s% m(k0)/Msun
         
         ! remove cells k0..nz
         call do_remove_inner_fraction_q(id, s% q(k0), ierr)

         if (s% job% remove_center_set_zero_v_center) s% v_center = 0d0

      end subroutine do_remove_fallback


      subroutine do_remove_center_by_logRho(id, logRho_limit, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: logRho_limit
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k0
         real(dp) :: lnd_limit
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_logRho: get_star_ptr ierr', ierr
            return
         end if
         
         lnd_limit = logRho_limit*ln10
         k0 = 0
         do k = s% nz,1,-1
            if (s% lnd(k) < lnd_limit) then
               k0 = k
               exit
            end if
            if (s% q(k) > 0.01d0) return
         end do
         
         ! k0 is innermost cell with density below limit
         ! search out from there for outermost with density too low
         
         do k = k0,1,-1
            if (s% lnd(k) < lnd_limit) cycle
            call do_remove_inner_fraction_q(id, s% q(k), ierr)
            return
         end do

      end subroutine do_remove_center_by_logRho


      subroutine do_limit_center_logP(id, logP_limit, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: logP_limit
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k0
         real(dp) :: lnP_limit
         include 'formats'        
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_limit_center_logP: get_star_ptr ierr', ierr
            return
         end if         
         lnP_limit = logP_limit*ln10
         k = s% nz
         if (s% lnP(k) > lnP_limit) then
            call do_remove_inner_fraction_q(id, s% q(k), ierr)
            if (ierr == 0) &
               write(*,3)' remove center for logP_limit', &
                  k, s% nz, s% m(k)/Msun, s% lnP(k)/ln10, logP_limit
         end if
      end subroutine do_limit_center_logP


      subroutine do_remove_center_by_ye(id, ye, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: ye
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_ye: get_star_ptr ierr', ierr
            return
         end if
         do k=1,s% nz
            if (s% ye(k) <= ye) then
               call do_remove_inner_fraction_q(id, s% q(k), ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_center_by_ye


      subroutine do_remove_center_by_entropy(id, entropy, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: entropy
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_entropy: get_star_ptr ierr', ierr
            return
         end if
         do k=1,s% nz
            if (s% entropy(k) <= entropy) then
               call do_remove_inner_fraction_q(id, s% q(k), ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_center_by_entropy


      subroutine do_remove_center_by_infall_kms(id, infall_kms, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: infall_kms
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: v
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_infall_kms: get_star_ptr ierr', ierr
            return
         end if
         v = -abs(infall_kms)*1d5
         do k=1,s% nz
            if (s% v(k) <= v) then
               call do_remove_inner_fraction_q(id, s% q(k), ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_center_by_infall_kms


      subroutine do_remove_center_by_radius_cm(id, r, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: r
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: q_r, rp1, r00, qp1
         integer :: k

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_radius_cm: get_star_ptr ierr', ierr
            return
         end if
         rp1 = s% R_center
         if (rp1 > r) then
            ierr = -1
            s% retry_message = 'error in remove center by radius: r < R_center'
            if (s% report_ierr) write(*, *) s% retry_message
         end if
         if (s% r(1) <= r) then
            ierr = -1
            s% retry_message = 'error in remove center by radius: r >= R_surface'
            if (s% report_ierr) write(*, *) s% retry_message
         end if
         if (rp1 == r) return
         qp1 = 0d0
         do k=s% nz, 1, -1
            r00 = s% r(k)
            if (r00 > r .and. r >= rp1) then
               q_r = qp1 + s% dq(k)* &
                  (r*r*r - rp1*rp1*rp1)/(r00*r00*r00 - rp1*rp1*rp1)
               exit
            end if
            rp1 = r00
            qp1 = s% q(k)
         end do
         call do_remove_inner_fraction_q(id, q_r, ierr)

      end subroutine do_remove_center_by_radius_cm


      subroutine do_remove_center_at_inner_max_abs_v(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: q_max, abs_v_max
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_at_inner_max_abs_v: get_star_ptr ierr', ierr
            return
         end if
         q_max = 0
         abs_v_max = 0
         do k=s% nz-1, 1, -1
            if (abs(s% v(k)) > abs_v_max) then
               q_max = s% q(k)
               abs_v_max = abs(s% v(k))
               exit
            end if
         end do
         call do_remove_inner_fraction_q(id, q_max, ierr)
      end subroutine do_remove_center_at_inner_max_abs_v


      subroutine do_remove_fe_core(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_fe_core: get_star_ptr ierr', ierr
            return
         end if
         if (s% fe_core_k <= 0 .or. s% fe_core_k > s% nz) return
         call do_remove_inner_fraction_q(id, s% q(s% fe_core_k), ierr)
      end subroutine do_remove_fe_core


      subroutine do_remove_center_by_mass_gm(id, m, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: m
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: q_m
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_mass_gm: get_star_ptr ierr', ierr
            return
         end if
         q_m = (m - s% M_center)/s% xmstar
         if (q_m <= 0d0) return
         call do_remove_inner_fraction_q(id, q_m, ierr)
      end subroutine do_remove_center_by_mass_gm


      subroutine do_remove_inner_fraction_q(id, q, ierr)
         ! note: we have not implemented fractional cell removal.
         ! it adds a lot of complexity and, so far, we don't need it.
         use alloc, only: prune_star_info_arrays
         use star_utils, only: set_qs
         integer, intent(in) :: id
         real(dp), intent(in) :: q
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_inner_fraction_q: get_star_ptr ierr', ierr
            return
         end if
         if (q < 0d0 .or. q > 1d0) then
            ierr = -1
            s% retry_message = 'error in remove center: invalid location q'
            if (s% report_ierr) write(*, *) s% retry_message
            return
         end if
         do k = 1, s% nz
            if (s% q(k) <= q) exit
         end do
         call do_remove_center(id, k, ierr)
      end subroutine do_remove_inner_fraction_q


      subroutine do_remove_center(id, k, ierr)
         use read_model, only: finish_load_model
         use alloc, only: prune_star_info_arrays
         use star_utils, only: set_qs
         integer, intent(in) :: id, k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: old_xmstar, new_xmstar
         integer :: kk
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center: get_star_ptr ierr', ierr
            return
         end if
         if (k <= 1 .or. k > s% nz) return
         old_xmstar = s% xmstar
         s% M_center = s% m(k)
         new_xmstar = s% m(1) - s% M_center
         s% xmstar = new_xmstar
         s% R_center = s% r(k)
         if (s% job% remove_center_adjust_L_center) s% L_center = s% L(k)
         if (s% u_flag) then
            s% v_center = s% u(k)
            if (is_bad(s% v_center)) then
               write(*,2) 's% u(k)', k, s% u(k)
               stop 'do_remove_center'
            end if
         else if (s% v_flag) then
            s% v_center = s% v(k)
            if (is_bad(s% v_center)) then
               write(*,2) 's% v(k)', k, s% v(k)
               stop 'do_remove_center'
            end if
         else
            s% v_center = 0d0
         end if
         if (s% job% remove_center_set_zero_v_center) s% v_center = 0d0
         s% nz = k-1
         do kk=1,k-1
            s% dq(kk) = s% dm(kk)/new_xmstar
         end do
         call set_qs(s% nz, s% q, s% dq, ierr)
         if (ierr /= 0) return
         s% generations = 1 ! memory leak, but hopefully not necessary to fix
            ! assuming remove center is a rare operation
         call prune_star_info_arrays(s, ierr)
         if (ierr /= 0) return
         s% need_to_setvars = .true.
         call finish_load_model(s, .false., .false., .false., ierr)
      end subroutine do_remove_center


      subroutine do_zero_inner_v_by_mass_gm(id, m, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: m
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, i_u, i_v
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_zero_inner_v_by_mass_gm: get_star_ptr ierr', ierr
            return
         end if
         i_u = s% i_u
         i_v = s% i_v
         do k=s% nz, 1, -1
            if (s% m(k) > m) exit
            if (i_u > 0) then
               s% xh(i_u, k) = 0d0
               s% u(k) = 0d0
            end if
            if (i_v > 0) then
               s% xh(i_v, k) = 0d0
               s% v(k) = 0d0
            end if
         end do
      end subroutine do_zero_inner_v_by_mass_gm


      subroutine do_remove_surface_at_cell_k(id, k, ierr)
         integer, intent(in) :: id, k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_at_cell_k: get_star_ptr ierr', ierr
            return
         end if
         if (k <= 1) return
         call do_remove_surface(id, k, ierr)
      end subroutine do_remove_surface_at_cell_k


      subroutine do_remove_surface_by_optical_depth(id, optical_depth, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: optical_depth
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_optical_depth: get_star_ptr ierr', ierr
            return
         end if
         if (optical_depth <= s% tau(1)) return
         do k=1,s% nz
            if (s% tau(k) >= optical_depth) then
               call do_remove_surface(id, k, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_by_optical_depth


      subroutine do_remove_surface_by_v_surf_km_s(id, v_surf_km_s, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: v_surf_km_s
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp), dimension(:), pointer :: v
         real(dp) :: v_max
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_v_surf_km_s: get_star_ptr ierr', ierr
            return
         end if
         if (s% u_flag) then
            v => s% u
         else if (s% v_flag) then
            v => s% v
         else
            return
         end if
         v_max = 1d5*v_surf_km_s
         if (v(1) < v_max) return
         do k=2,3 ! s% nz
            if (v(k) < v_max) exit
            write(*,2) 'v', k-1, v(k-1)/1d5, v_surf_km_s
         end do
         write(*,2) 'do_remove_surface_by_v_surf_km_s', k-1, v(k-1)/1d5, v_surf_km_s
         call do_remove_surface(id, k-1, ierr)
         return
      end subroutine do_remove_surface_by_v_surf_km_s


      subroutine do_remove_surface_by_v_surf_div_cs(id, v_surf_div_cs, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: v_surf_div_cs
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp), dimension(:), pointer :: v
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_v_surf_div_cs: get_star_ptr ierr', ierr
            return
         end if
         if (s% u_flag) then
            v => s% u
         else if (s% v_flag) then
            v => s% v
         else
            return
         end if
         !write(*,1) 'v(1)/cs', v(1)/s% csound(1), v_surf_div_cs
         if (v(1) < s% csound(1)*v_surf_div_cs) return
         do k=2,30 ! s% nz
            if (v(k) < s% csound(k)*v_surf_div_cs) exit
            write(*,2) 'v/cs', k-1, v(k-1)/s% csound(k-1)
         end do
         write(*,2) 'do_remove_surface_by_v_surf_div_cs', k-1, v(k-1)/s% csound(k-1)
         call do_remove_surface(id, k-1, ierr)
         !stop 'do_remove_surface_by_v_surf_div_cs'
         return
      end subroutine do_remove_surface_by_v_surf_div_cs


      subroutine do_remove_surface_by_v_surf_div_v_escape(id, v_surf_div_v_escape, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: v_surf_div_v_escape
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: vesc, vesc_m1
         real(dp), dimension(:), pointer :: v
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_v_surf_div_v_escape: get_star_ptr ierr', ierr
            return
         end if
         if (s% u_flag) then
            v => s% u
         else if (s% v_flag) then
            v => s% v
         else
            return
         end if
         vesc = sqrt(2*s% cgrav(1)*s% m(1)/(s% r(1)))
         if (v(1) < vesc*v_surf_div_v_escape) return
         do k=2,3 ! s% nz
            vesc_m1 = vesc
            vesc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (v(k) < vesc*v_surf_div_v_escape) exit
            write(*,2) 'v/vesc', k-1, v(k-1)/vesc_m1
         end do
         write(*,2) 'do_remove_surface_by_v_surf_div_v_escape', k-1, v(k-1)/vesc_m1
         call do_remove_surface(id, k-1, ierr)
         return
      end subroutine do_remove_surface_by_v_surf_div_v_escape


      subroutine do_remove_surface_by_pressure(id, pressure, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: pressure
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_pressure: get_star_ptr ierr', ierr
            return
         end if
         if (pressure <= s% P(1)) return
         do k=1,s% nz
            if (s% P(k) >= pressure) then
               call do_remove_surface(id, k, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_by_pressure


      subroutine do_remove_surface_by_density(id, density, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: density
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: avg_rho
         ierr = 0
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_density: get_star_ptr ierr', ierr
            return
         end if
         if (s% rho(1) < density) then
            write(*,2) 'do_remove_surface', 1, s% rho(1), density
            call do_remove_surface(id, 2, ierr)
            return
         end if
      end subroutine do_remove_surface_by_density


      subroutine do_remove_surface_by_radius_cm(id, r, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: r
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_radius_cm: get_star_ptr ierr', ierr
            return
         end if
         if (r >= s% r(1)) return
         do k=1,s% nz
            if (s% r(k) <= r) then
               call do_remove_surface(id, k, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_by_radius_cm


      subroutine do_remove_surface_by_mass_gm(id, m, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: m
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_mass_gm: get_star_ptr ierr', ierr
            return
         end if
         if (m >= s% m(1)) return
         do k=1,s% nz
            if (s% m(k) <= m) then
               call do_remove_surface(id, k, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_by_mass_gm


      subroutine do_remove_surface_by_q(id, q, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: q
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_by_q: get_star_ptr ierr', ierr
            return
         end if
         if (q >= 1d0) return
         do k=1,s% nz
            if (s% q(k) <= q) then
               call do_remove_surface(id, k, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_by_q


      subroutine do_remove_surface(id, surface_k, ierr)
         use read_model, only: finish_load_model
         use mesh_adjust, only: do_prune_mesh_surface
         use alloc, only: resize_star_info_arrays
         use star_utils, only: tau_eff
         integer, intent(in) :: id, surface_k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s, c, prv
         type (star_info), target :: copy_info
         real(dp) :: tau_surf_new, tau_factor_new, Lmid, Rmid, T, P, T_black_body
         integer :: k, k_old, nz, nz_old, skip

         logical, parameter :: dbg = .false.

         include 'formats'

         if (surface_k == 1) then
            ierr = 0
            return
         end if

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'do_remove_surface: get_star_ptr ierr'
            return
         end if

         nz_old = s% nz
         skip = surface_k - 1

         if (dbg) write(*,2) 'do remove surface skip', skip
         if (skip < 1 .or. skip >= nz_old) return

         tau_surf_new = tau_eff(s,1+skip)
         tau_factor_new = s% tau_factor*tau_surf_new/s% tau(1)

         if (dbg) write(*,1) 'tau_surf_old', s% tau(1)
         if (dbg) write(*,1) 'tau_factor_old', s% tau_factor
         if (dbg) write(*,1) 'tau_surf_new', tau_surf_new
         if (dbg) write(*,1) 'tau_factor_new', tau_factor_new

         rmid = s% rmid(1+skip)
         Lmid = (s% L(1+skip) + s% L(2+skip))/2
         T = s% T(1+skip)
         P = s% P(1+skip)

         if (.not. associated(s% other_star_info)) then
            allocate(s% other_star_info)
            prv => s% other_star_info
            c => null()
            if (dbg) write(*,1) 'c is null'
         else
            prv => s% other_star_info
            c => copy_info
            c = prv
         end if

         prv = s ! this makes copies of pointers and scalars

         nz = nz_old - skip
         s% nz = nz
         if (dbg) write(*,2) 'nz_old', nz_old
         if (dbg) write(*,2) 'nz', nz

         if (dbg) write(*,1) 'call resize_star_info_arrays'
         call resize_star_info_arrays(s, c, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'resize_star_info_arrays failed in do_remove_surface'
            return
         end if

         if (dbg) write(*,2) 'prv% dm(1)/Msun', 1, prv% dm(1)/Msun
         if (dbg) write(*,2) 'prv% m(1)/msun', 1, prv% m(1)/Msun
         if (dbg) write(*,2) 'prv% dm(nz_old)/Msun', nz_old, prv% dm(nz_old)/Msun
         if (dbg) write(*,2) 'prv% m(nz_old)/msun', nz_old, prv% m(nz_old)/Msun

         s% mstar = prv% m(1 + skip)
         s% xmstar = s% mstar - prv% M_center
         s% q(1) = 1d0
         do k = 1, nz
            k_old = k + skip
            s% dq(k) = prv% dm(k_old)/s% xmstar
            if (k > 1) s% q(k) = s% q(k-1) - s% dq(k-1)
         end do
         s% dq(nz) = s% q(nz)
         do k = 1, nz
            s% dm(k) = s% dq(k)*s% xmstar
            s% m(k) = s% q(k)*s% xmstar + s% M_center
         end do
         if (dbg) write(*,2) 's% dq(1)', 1, s% dq(1)
         if (dbg) write(*,2) 's% dm(1)/Msun', 1, s% dm(1)/Msun
         if (dbg) write(*,2) 's% m(1)/msun', 1, s% m(1)/Msun
         if (dbg) write(*,2) 's% dm(nz)/Msun', nz, s% dm(nz)/Msun
         if (dbg) write(*,2) 's% m(nz)/msun', nz, s% m(nz)/Msun

         if (dbg) write(*,1) 'call do_prune_mesh_surface'
         call do_prune_mesh_surface( &
            s, nz, nz_old, prv% xh, prv% xa, &
            prv% j_rot, prv% i_rot, &
            prv% omega, prv% D_omega, prv% am_nu_rot, &
            prv% D_smooth, prv% conv_vel, prv% lnT, &
            prv% dPdr_dRhodr_info, prv% nu_ST, prv% D_ST, prv% D_DSI, prv% D_SH, &
            prv% D_SSI, prv% D_ES, prv% D_GSF, prv% D_mix, &
            s% xh, s% xa, ierr)
         if (ierr /= 0) then
            return
         end if
         if (dbg) write(*,2) 's% dq(1)', 1, s% dq(1)
         if (dbg) write(*,2) 's% dm(1)/Msun', 1, s% dm(1)/Msun
         if (dbg) write(*,2) 's% m(1)/msun', 1, s% m(1)/Msun
         if (dbg) write(*,2) 's% dm(nz)/Msun', nz, s% dm(nz)/Msun
         if (dbg) write(*,2) 's% m(nz)/msun', nz, s% m(nz)/Msun

         if (Lmid > 0d0) then
            T_black_body = pow(Lmid/(4*pi*rmid*rmid*boltz_sigma), 0.25d0)
            s% Tsurf_factor = T/T_black_body
         else
            s% Tsurf_factor = 1d0
         end if
         s% force_Tsurf_factor = s% Tsurf_factor

         if (s% use_momentum_outer_bc) then
            s% tau_factor = tau_factor_new
            s% force_tau_factor = s% tau_factor
         end if

         s% need_to_setvars = .true.

         if (dbg) write(*,1) 'call finish_load_model'
         call finish_load_model(s, .false., .false., .false., ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'finish_load_model failed in do_remove_surface'
            return
         end if

         if (dbg) write(*,1) 'do_remove_surface tau_factor, Tsurf_factor', &
            s% tau_factor, s% Tsurf_factor
            
         if (dbg) stop 'do_remove_surface'

      end subroutine do_remove_surface


      end module remove_shells
