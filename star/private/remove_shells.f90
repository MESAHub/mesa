! ***********************************************************************
!
!   Copyright (C) 2015-2019  The MESA Team
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
         real(dp), pointer :: v(:)
         
         include 'formats'
         
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_fallback: get_star_ptr ierr', ierr
            return
         end if
         
         if (s% u_flag) then
            v => s% u
         else if (s% v_flag) then
            v => s% v
         else
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
               ke = 0.5d0*v(k)*v(k)*s% dm(k)
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
               ke = 0.5d0*v(k)*v(k)*s% dm(k)
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
            if (-v(nz) <= speed_limit*s% csound(nz)) return
            do k = nz-1,1,-1
               k0 = k+1
               if (-v(k) < speed_limit*s% csound(k)) exit
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
         if (s% lnPeos(k) > lnP_limit) then
            call do_remove_inner_fraction_q(id, s% q(k), ierr)
            if (ierr == 0) &
               write(*,3)' remove center for logP_limit', &
                  k, s% nz, s% m(k)/Msun, s% lnPeos(k)/ln10, logP_limit
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
         integer :: k, k_infall
         real(dp) :: v_infall, v
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_by_infall_kms: get_star_ptr ierr', ierr
            return
         end if
         v_infall = -abs(infall_kms)*1d5
         k_infall = 0
         do k=s% nz,1,-1
            if (s% v_flag) then
               v = s% v(k)
            else if (s% u_flag) then
               v = s% u(k)
            else
               write(*,*) 'must have v or u for do_remove_center_by_infall_kms'
               ierr = -1
               return
            end if
            if (v > v_infall) exit ! not falling fast enough
            k_infall = k
         end do
         if (k_infall == 0) return ! no infall
         call do_remove_inner_fraction_q(id, s% q(k_infall), ierr)
         write(*,1) 'new inner boundary mass', s% m_center/Msun
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
         integer :: k_max
         real(dp), pointer :: v(:)
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_center_at_inner_max_abs_v: get_star_ptr ierr', ierr
            return
         end if
         if (s% u_flag) then
            v => s% u
         else if (s% v_flag) then
            v => s% v
         else
            stop 'no u or v for do_remove_center_at_inner_max_abs_v?'
            return
         end if
         k_max = minloc(v(1:s% nz),dim=1)
         q_max = s% q(k_max)
         abs_v_max = abs(v(k_max))
         !write(*,2) 'v abs_v_max q_max', k_max, v(k_max), abs_v_max, q_max
         call do_remove_center(id, k_max, ierr)
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
         use star_utils, only: normalize_dqs, set_qs
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
            kk = minloc(s% u(1:s% nz),dim=1)
            s% v_center = s% u(k)
            if (is_bad(s% v_center)) then
               write(*,2) 'center s% u(k)', k, s% u(k)
               stop 'do_remove_center'
            end if
         else if (s% v_flag) then
            s% v_center = s% v(k)
            if (is_bad(s% v_center)) then
               write(*,2) 'center s% v(k)', k, s% v(k)
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
         if (.not. s% do_normalize_dqs_as_part_of_set_qs) then
            call normalize_dqs(s, s% nz, s% dq, ierr)
            if (ierr /= 0) then
               if (s% report_ierr) write(*,*) 'normalize_dqs failed in do_remove_center'
               return
            end if
         end if
         call set_qs(s, s% nz, s% q, s% dq, ierr)
         if (ierr /= 0) return
         s% generations = 1 ! memory leak, but hopefully not necessary to fix
            ! assuming remove center is a rare operation
         call prune_star_info_arrays(s, ierr)
         if (ierr /= 0) return
         s% need_to_setvars = .true.
         call finish_load_model(s, .false., .false., .false., .false., .false., ierr)
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


      subroutine do_remove_surface_at_he_core_boundary(id, h1_fraction, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: h1_fraction
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'do_remove_surface_at_he_core_boundary: get_star_ptr ierr', ierr
            return
         end if
         if (s% X(1) <= h1_fraction) return
         do k=2,s% nz
            if (s% X(k) <= h1_fraction) then
               call do_remove_surface(id, k-1, ierr)
               return
            end if
         end do
         ierr = -1
      end subroutine do_remove_surface_at_he_core_boundary


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
         integer :: k, k_vesc
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
         k_vesc = 0
         do k=2, s% nz
            if (s% q(k) > s% job% max_q_for_remove_surface_by_v_surf_div_v_escape) cycle
            if (s% q(k) < s% job% min_q_for_remove_surface_by_v_surf_div_v_escape) exit
            vesc = sqrt(2*s% cgrav(k)*s% m(k)/(s% r(k)))
            if (v(k) >= vesc*v_surf_div_v_escape) k_vesc = k
         end do
         if (k_vesc == 0) return
         write(*,2) 'do_remove_surface_by_v_surf_div_v_escape q', k_vesc, s% q(k_vesc)
         call do_remove_surface(id, k_vesc, ierr)
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
         if (pressure <= s% Peos(1)) return
         do k=1,s% nz
            if (s% Peos(k) >= pressure) then
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

         if (s% job% remove_surface_by_relax_to_star_cut) then
            if (s% R_center /= 0d0) then
               write(*,*) 'remove surface currently requires model with inner boundary at true center of star'
               ierr = -1
               stop 'do_remove_surface'
            end if         
            call do_relax_to_star_cut( &
               id, surface_k, s% job% remove_surface_do_jrot, &
               s% job% remove_surface_do_entropy, &
               s% job% remove_surface_turn_off_energy_sources_and_sinks, ierr)
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
         P = s% Peos(1+skip)

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
            prv% conv_vel, prv% lnT, &
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
            T_black_body = pow(Lmid/(pi4*rmid*rmid*boltz_sigma), 0.25d0)
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
         call finish_load_model(s, .false., .false., .false., .false., .false., ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*,*) 'finish_load_model failed in do_remove_surface'
            return
         end if

         if (dbg) write(*,1) 'do_remove_surface tau_factor, Tsurf_factor', &
            s% tau_factor, s% Tsurf_factor
            
         if (dbg) stop 'do_remove_surface'
            
      end subroutine do_remove_surface


      ! Relax to a trimmed stellar model with surface cells removed down to k_remove
      ! (the cell k_remove will be the outermost in the new model).
      subroutine do_relax_to_star_cut( &
            id, k_remove, do_jrot, do_entropy, turn_off_energy_sources_and_sinks, ierr)

         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interp_pm, interp_values, interp_value
         use adjust_xyz, only: change_net
         use alloc, only: set_conv_vel_flag, set_v_flag, set_u_flag, set_rotation_flag
         use rotation_mix_info, only: set_rotation_mixing_info
         use hydro_rotation, only: set_i_rot, set_rotation_info
         use relax, only: do_relax_composition, do_relax_angular_momentum, do_relax_entropy
         use init, only: load_zams_model

         integer, intent(in) :: id, k_remove
         logical, intent(in) :: do_jrot, do_entropy
         logical, intent(in) :: turn_off_energy_sources_and_sinks 
            ! determines if we turn off non_nuc_neu and eps_nuc for entropy relax
         integer, intent(out) :: ierr

         logical :: conv_vel_flag, v_flag, u_flag, rotation_flag
         type (star_info), pointer :: s
         character (len=net_name_len) :: net_name
         integer :: model_number, num_trace_history_values, photo_interval
         real(dp) :: eps_nuc_factor, non_nuc_neu_factor, &
            initial_z, initial_y, initial_mass, &
            cumulative_energy_error, cumulative_extra_heating

         real(dp), pointer :: interp_work(:), conv_vel_interp(:)
         real(dp), pointer :: q(:), xq(:), xa(:,:), j_rot(:), entropy(:)
         real(dp) :: conv_vel_temp, time
         integer :: num_pts, k, k0, species
         logical :: save_have_mlt_vc
         logical :: dbg = .false.

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         eps_nuc_factor = s% eps_nuc_factor
         non_nuc_neu_factor = s% non_nuc_neu_factor
         net_name = s% net_name
         num_trace_history_values = s% num_trace_history_values

         time = s% time
         model_number = s% model_number
         num_trace_history_values = s% num_trace_history_values
         cumulative_energy_error = s% cumulative_energy_error
         cumulative_extra_heating = s% cumulative_extra_heating

         ! zero model_number and time (will restore later)
         s% model_number = 0
         s% time = 0

         species = s% species
         num_pts = s% nz - k_remove + 1
         allocate(q(num_pts), xq(num_pts), xa(species, num_pts))
         rotation_flag = .false.
         if (do_jrot .and. s% rotation_flag) then
            allocate(j_rot(num_pts))
            rotation_flag = .true.
         end if
         if (do_entropy) then
            allocate(entropy(num_pts))
         end if
         !need to compute cell-centered q for remaining mass
         xq(1) = s% dq(k_remove)/2/s% q(k_remove)
         do k0 = 1, num_pts-1
            xq(1+k0) = xq(1+k0-1) + (s% dq(k_remove+k0) + s% dq(k_remove+k0-1))/s% q(k_remove)/2
         end do

         !create interpolant for convective velocities
         conv_vel_flag = .false.
         if (s% conv_vel_flag) then
            conv_vel_flag = .true.
            allocate(interp_work((num_pts)*pm_work_size), &
               conv_vel_interp(4*(num_pts)), stat=ierr)
            do k0 = 1, num_pts
               conv_vel_interp(4*k0-3) = s% conv_vel(k0+k_remove-1)
               q(k0) = s% q(k0+k_remove-1)/s% q(k_remove)
            end do
            call interp_pm(q, num_pts, conv_vel_interp,&
               pm_work_size, interp_work, 'conv_vel interpolant', ierr)

            ! turn off conv vel flag to load model
            call set_conv_vel_flag(id, .false., ierr)
            if (dbg) write(*,*) "set_conv_vel_flag ierr", ierr
         end if

         ! save have_mlt_vc and set to false (to load ZAMS model)
         save_have_mlt_vc = s% have_mlt_vc
         s% have_mlt_vc = .false.

         !save composition and entropy profiles
         xa(:,:) = s% xa(:,k_remove:s% nz)
         if (rotation_flag) then
            j_rot(:) = s% j_rot(k_remove:s% nz)
         end if
         if (do_entropy) then
            entropy(:) = s% entropy(k_remove:s% nz)*avo*kerg
         end if

         ! various flags need to be turned off for the ZAMS model to load
         v_flag = .false.
         if (s% v_flag) then
            call set_v_flag(id, .false., ierr)
            if (dbg) write(*,*) "set_v_flag ierr", ierr
            v_flag = .true.
         end if
         u_flag = .false.
         if (s% u_flag) then
            call set_u_flag(id, .false., ierr)
            if (dbg) write(*,*) "set_u_flag ierr", ierr
            u_flag = .true.
         end if

         if (s% rotation_flag) then
            call set_rotation_flag(id, .false., ierr)
            if (dbg) write(*,*) "set_rotation_flag ierr", ierr
         end if

         ! avoid making photos
         photo_interval = s% photo_interval
         s% photo_interval = 10000000
         s% have_previous_conv_vel = .false.
         s% have_j_rot = .false.
         ! WARNING, might need to add stuff here to actually get the ZAMS model to load.
         ! otherwise can get an error of the form "error in reading model data  j+species > nvec"
         ! if you happen to run into these problem, check for flags being checked in read1_model in read_model.f90
         ! and be sure they're turned off.

         ! set values used to load the starting model that will be relaxed
         initial_z = s% initial_z
         initial_y = s% initial_y
         initial_mass = s% initial_mass
         s% initial_z = 0.02d0
         s% initial_y = 0.28d0
         s% initial_mass = s% m(k_remove)/Msun

         s% prev_mesh_nz = 0

         call change_net(id, .true., 'basic.net', ierr) ! TODO:need to allow specification of different net
         if (dbg) write(*,*) "check change_net ierr", ierr
         if (ierr /= 0) return
         call load_zams_model(id, ierr)
         if (dbg) write(*,*) "check load_zams ierr", ierr
         if (ierr /= 0) return
         call change_net(id, .true., net_name, ierr)
         if (dbg) write(*,*) "check ierr", ierr
         if (ierr /= 0) return

         if (conv_vel_flag) then
            call set_conv_vel_flag(id, .true., ierr)
            if (dbg) write(*,*) "check set_conv_vel_flag ierr", ierr
            if (ierr /= 0) return
         end if

         ! restore have_mlt_vc
         s% have_mlt_vc = save_have_mlt_vc

         if (turn_off_energy_sources_and_sinks) then
            s% non_nuc_neu_factor = 0d0
            s% eps_nuc_factor = 0d0
         end if

         s% num_trace_history_values = 0
         call do_relax_composition( &
            id, s% job% num_steps_to_relax_composition, num_pts, species, xa, xq, ierr)
         if (dbg) write(*,*) "check ierr", ierr
         if (ierr /= 0) return
         deallocate(xa)

         if (rotation_flag) then
            call set_rotation_flag(id, .true., ierr)
            if (dbg) write(*,*) "set_rotation_flag true ierr", ierr
            if (ierr /= 0) return
            call set_rotation_info(s, .false., ierr)
            if (dbg) write(*,*) "set_rotation_info ierr", ierr
            if (ierr /= 0) return
            call set_rotation_mixing_info(s, ierr)
            if (dbg) write(*,*) "set_rotation_mixing_info ierr", ierr
            if (ierr /= 0) return
            call do_relax_angular_momentum( &
               id, s% job% max_steps_to_relax_angular_momentum, num_pts, j_rot, xq, ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
            deallocate(j_rot)
         end if

         if (do_entropy) then
            call do_relax_entropy( &
               id, s% job% max_steps_to_relax_entropy, num_pts, entropy, xq, ierr)
            if (dbg) write(*,*) "check ierr", ierr
            if (ierr /= 0) return
            deallocate(entropy)
         end if

         !take care of convective velocities
         if (s% conv_vel_flag) then
            do k0=1, s% nz
               call interp_value(q, num_pts, conv_vel_interp, s% q(k0), s% conv_vel(k0), ierr)
               !avoid extending regions with non-zero conv vel
               do k=2, num_pts-1
                  if(s% q(k0) < q(k) .and. s% q(k0) > q(k+1) &
                     .and. (conv_vel_interp(4*k-3)<1d-5 .or. conv_vel_interp(4*(k+1)-3)<1d-5)) then
                     s% conv_vel(k0) = 0d0
                     exit
                  end if
               end do
               s% xh(s% i_ln_cvpv0, k0) = log(s% conv_vel(k0)+s% conv_vel_v0)
            end do
            write(*,*) 'need to rewrite some things here in do_relax_to_star_cut'
            stop 'do_relax_to_star_cut'
            deallocate(conv_vel_interp, interp_work)
         end if

         s% generations = 1

         ! restore v_flag and u_flag
         if (v_flag) then
            call set_v_flag(id, .true., ierr)
         end if
         if (u_flag) then
            call set_u_flag(id, .true., ierr)
         end if

         ! this avoids the history file from being rewritten
         s% doing_first_model_of_run = .false.

         s% time = time
         s% model_number = model_number
         s% num_trace_history_values = num_trace_history_values
         s% cumulative_energy_error = cumulative_energy_error
         s% cumulative_extra_heating = cumulative_extra_heating

         s% non_nuc_neu_factor = non_nuc_neu_factor
         s% eps_nuc_factor = eps_nuc_factor

         s% initial_z = initial_z
         s% initial_y = initial_y
         s% initial_mass = initial_mass
         s% photo_interval = photo_interval

         deallocate(q, xq)
         
         s% need_to_setvars = .true.

      end subroutine do_relax_to_star_cut


      end module remove_shells
