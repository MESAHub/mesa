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

      module star_utils

      use star_private_def
      use const_def
      use num_lib
      use utils_lib

      implicit none

      private :: mdb
      logical, parameter :: mdb = .false.

      contains


      subroutine foreach_cell(s,nzlo,nzhi,use_omp,do1,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         logical, intent(in) :: use_omp
         interface
            subroutine do1(s,k,ierr)
               use star_private_def
               type (star_info), pointer :: s
               integer, intent(in) :: k
               integer, intent(out) :: ierr
            end subroutine do1
         end interface
         integer, intent(out) :: ierr

         integer :: k, op_err
         logical :: okay
         ierr = 0

         if (nzlo == nzhi) then
            call do1(s,nzlo,ierr)
            return
         end if

         if (use_omp) then
            okay = .true.
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
            do k = nzlo, nzhi
               if (.not. okay) cycle
               op_err = 0
               call do1(s,k,op_err)
               if (op_err /= 0) okay = .false. ! cannot just exit from a parallel loop
            end do
!$OMP END PARALLEL DO
            if (.not. okay) ierr = -1
         else
            do k = nzlo, nzhi
               call do1(s,k,ierr)
               if (ierr /= 0) exit
            end do
         end if

      end subroutine foreach_cell


      subroutine get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         real(dp), intent(out) :: y_avg, z_avg
         integer, intent(out) :: ierr
         integer :: k, nz,  h1, h2, he3, he4
         real(dp) :: total_mass_h, total_mass_he, total_mass_z, &
            cell_mass, total_mass
         ierr = 0
         nz = s% nz
         h1 = s% net_iso(ih1)
         h2 = s% net_iso(ih2)
         he3 = s% net_iso(ihe3)
         he4 = s% net_iso(ihe4)
         total_mass=0; total_mass_h=0; total_mass_he=0; total_mass_z=0
         do k=nzlo, nzhi
            cell_mass = s% dm(k)
            total_mass = total_mass + cell_mass
            total_mass_h = total_mass_h + cell_mass*s% xa(h1, k)
            if (h2 /= 0) total_mass_h = total_mass_h + cell_mass*s% xa(h2, k)
            total_mass_he = total_mass_he + cell_mass*s% xa(he4, k)
            if (he3 /= 0) total_mass_he = total_mass_he + cell_mass*s% xa(he3, k)
         end do
         total_mass_z = total_mass - (total_mass_h + total_mass_he)
         z_avg = total_mass_z / total_mass
         y_avg = total_mass_he / total_mass
      end subroutine get_average_Y_and_Z


      real(dp) function eval_current_y(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         real(dp) :: y_avg, z_avg
         call get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         eval_current_y = y_avg
      end function eval_current_y


      real(dp) function eval_current_z(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         real(dp) :: y_avg, z_avg
         call get_average_Y_and_Z(s, nzlo, nzhi, y_avg, z_avg, ierr)
         eval_current_z = z_avg
      end function eval_current_z


      real(dp) function eval_current_abundance(s, j, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: j, nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: cell_mass, jmass, total_mass

         ierr = 0

         if (j == 0) then
            eval_current_abundance = 0
            return
         end if

         nz = s% nz
         total_mass=0; jmass=0
         do k=nzlo, nzhi
            cell_mass = s% dm(k)
            total_mass = total_mass + cell_mass
            jmass = jmass + cell_mass*s% xa(j, k)
         end do
         eval_current_abundance = jmass / total_mass

      end function eval_current_abundance


      subroutine smooth_abundances(s, cnt, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: cnt ! make this many passes
         integer, intent(in) :: nzlo, nzhi ! only smooth zones nzlo to nzhi inclusive
         integer, intent(out) :: ierr
         integer :: k, j, nz
         ierr = 0
         nz = s% nz
         do j = 1, cnt
            do k = max(nzlo,2), min(nzhi, nz)
               s% xa(:,k) = (s% xa(:,k-1) + s% xa(:,k) + s% xa(:,k+1))/3
            end do
            if (nzhi == nz) s% xa(:,nz) = (s% xa(:,nz-1) + s% xa(:,nz) + s% xa(:,nz))/3
            if (nzlo == 1) s% xa(:,1) = (s% xa(:,2) + s% xa(:,1) + s% xa(:,1))/3
         end do
         s% need_to_setvars = .true.
      end subroutine smooth_abundances


      integer function k_for_q(s, q)
         ! return k s.t. q(k) >= q > q(k)-dq(k)
         type (star_info), pointer :: s
         real(dp), intent(in) :: q
         integer :: k, nz
         nz = s% nz
         if (q >= 1) then
            k_for_q = 1; return
         else if (q <= s% q(nz)) then
            k_for_q = nz; return
         end if
         do k = 1, nz-1
            if (q > s% q(k+1)) then
               k_for_q = k; return
            end if
         end do
         k_for_q = nz
      end function k_for_q


      subroutine get_name_for_restart_file(n, num_digits, num)
         integer, intent(in) :: n, num_digits
         character (len=*), intent(out) :: num
         call get_string_for_model_number('x', n, num_digits, num)
      end subroutine get_name_for_restart_file


      subroutine get_string_for_model_number(prefix, n, num_digits, num)
         character (len=*), intent(in) :: prefix
         integer, intent(in) :: n, num_digits
         character (len=*), intent(out) :: num
         integer :: val
         character (len=32) :: fstring
         include 'formats'
         val = mod(n, 10**num_digits) ! wrap around
         if (val == 0) then
            write(num,*) n
            num = adjustl(num)
            return
         end if
        write(fstring,'( "(a,i",i2.2,".",i2.2,")" )') num_digits, num_digits
        write(num,fstring) trim(prefix), val
      end subroutine get_string_for_model_number


      subroutine report_xa_bad_nums(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, j
         ierr = 0
         do k=1,s% nz
            do j=1,s% species
               if (is_bad(s% xa(j,k))) then
                  ierr = -1
                  write(*,*) j, k, s% xa(j,k)
               end if
            end do
         end do
      end subroutine report_xa_bad_nums


      real(dp) function eval_csound(s,k,ierr) result(cs)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         real(dp) :: cs2
         include 'formats'
         ierr = 0
         cs2 = s% gamma1(k)*s% P(k)/s% rho(k)
         if (cs2 < 0d0) then
            cs = 0d0
            ierr = -1
            return
         end if
         cs = sqrt(cs2)
      end function eval_csound


      subroutine set_m_grav_and_grav(s) ! using mass_corrections
         type (star_info), pointer :: s
         integer :: k, nz
         include 'formats'
         nz = s% nz
         if (.not. s% use_mass_corrections) then
            do k=1,nz
               s% m_grav(k) = s% m(k)
            end do
         else
            s% m_grav(nz) = &
               s% M_center + s% dm(nz)*s% mass_correction(nz)
            do k=nz-1,1,-1
               s% m_grav(k) = &
                  s% m_grav(k+1) + s% dm(k)*s% mass_correction(k)
            end do
         end if
         do k=1,nz
            s% grav(k) = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
         end do
      end subroutine set_m_grav_and_grav


      subroutine use_xh_to_set_rho_to_dm_div_dV(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k, nz, i_lnR, i_lnd
         real(dp) :: rL, rR, dm, dV, rho, old_lnd, new_lnd
         include 'formats'
         ierr = 0
         i_lnR = s% i_lnR
         i_lnd = s% i_lnd
         if (i_lnd == 0 .or. i_lnR == 0) return
         nz = s% nz
         rR = s% R_center
         do k = nz, 1, -1
            rL = rR
            rR = exp(s% xh(i_lnR,k))
            dm = s% dm(k)
            dV = 4d0*pi/3d0*(rR*rR*rR - rL*rL*rL)
            rho = dm/dV
            if (rho <= 0d0) then
               write(*,3) 'set_rho_to_dm_div_dV: rho <= 0', &
                  k, nz, rho, dm, dV, rR, rL
               ierr = -1
               return
            end if
            new_lnd = log(rho)
            s% xh(i_lnd,k) = new_lnd
         end do
      end subroutine use_xh_to_set_rho_to_dm_div_dV


      subroutine set_m_and_dm(s)
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         do k = 1, s% nz
            s% m(k) = s% M_center + s% q(k)*s% xmstar
            s% dm(k) = s% dq(k)*s% xmstar
            if (s% dm(k) <= 0d0 .or. is_bad(s% m(k) + s% dm(k))) then
               write(*,2) 'dm m dq q M_center', k, &
                  s% dm(k), s% m(k), s% dq(k), s% q(k), s% M_center
               stop 'set_m_and_dm'
            end if
         end do
      end subroutine set_m_and_dm


      subroutine set_dm_bar(s, nz, dm, dm_bar)
         type (star_info), pointer :: s
         integer, intent(in) :: nz
         real(dp), intent(in) :: dm(:) ! (nz)
         real(dp), intent(inout) :: dm_bar(:) ! (nz)
         integer :: k
         do k=2,nz-1
            dm_bar(k) = 0.5d0*(dm(k-1) + dm(k))
         end do
         dm_bar(1) = 0.5d0*dm(1)
         if (s% rsp_flag) then ! rsp uses this definition
            dm_bar(nz) = 0.5d0*(dm(nz-1) + dm(nz))
         else
            dm_bar(nz) = 0.5d0*dm(nz-1) + dm(nz)
         end if
      end subroutine set_dm_bar


      subroutine normalize_dqs(nz, dq, ierr)
         ! rescale dq's so that add to 1.000
         ! work in from boundaries to meet at largest dq
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k, midq
         real(dp) :: dqsum1, dqsum2, dq_min
         include 'formats'
         k = minloc(dq(1:nz),dim=1)
         dq_min = dq(k)
         if (dq_min <= 0d0) then
            write(*,2) 'bad dq', k, dq(k)
            ierr = -1
            stop
            return
         end if
         midq = maxloc(dq(1:nz),dim=1)
         ! surface inward
         dqsum1 = 0
         do k=1, midq
            dqsum1 = dqsum1 + dq(k)
            if (dq(k) <= 0) then
               ierr = -1
               return
            end if
         end do
         ! center outward
         dqsum2 = 0
         do k=nz, midq+1, -1
            dqsum2 = dqsum2 + dq(k)
            if (dq(k) <= 0) then
               ierr = -1
               return
            end if
         end do
         do k=1,nz
            dq(k) = dq(k)/(dqsum1 + dqsum2)
         end do
      end subroutine normalize_dqs


      subroutine set_qs(nz, q, dq, ierr) ! set q's using normalized dq's
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         real(dp), intent(inout) :: q(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k, midq
         real(dp) :: dqsum1, dqsum2
         logical :: okay
         include 'formats'
         ierr = 0
         call normalize_dqs(nz, dq, ierr)
         if (ierr /= 0) return
         q(1) = 1d0
         okay = .true.
         do k=2,nz
            q(k) = q(k-1) - dq(k-1)
            if (q(k) < 0d0 .or. q(k) > 1d0) then
               okay = .false.
               exit
            end if
         end do
         if (okay) return
         midq = maxloc(dq(1:nz),dim=1)
         ! surface inward
         dqsum1 = 0
         do k=1, midq
            q(k) = 1d0 - dqsum1
            dqsum1 = dqsum1 + dq(k)
         end do
         ! center outward
         dqsum2 = 0
         do k=nz, midq+1, -1
            dqsum2 = dqsum2 + dq(k)
            q(k) = dqsum2
         end do         
      end subroutine set_qs


      subroutine set_xqs(nz, xq, dq, ierr) ! set xq's using dq's
         integer, intent(in) :: nz
         real(dp), intent(inout) :: dq(:) ! (nz)
         real(dp), intent(inout) :: xq(:) ! (nz)
         integer, intent(out) :: ierr
         integer :: k
         include 'formats'
         ierr = 0
         xq(1) = 0
         do k=2,nz-1
            xq(k) = xq(k-1) + dq(k-1)
         end do
         xq(nz) = 1 - dq(nz)
         if (xq(nz) < xq(nz-1)) then
            xq(nz) = xq(nz-1) + dq(nz-1)
            dq(nz) = 1 - xq(nz)
            if (dq(nz) <= 0) then
               ierr = -1
               return
            end if
         end if
      end subroutine set_xqs


      real(dp) function interp_val_to_pt(v,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         integer, intent(in) :: k, sz
         real(dp), intent(in) :: v(:), dq(:)
         character (len=*), intent(in) :: str
         integer :: ierr
         include 'formats'
         if (k == 1) then
            interp_val_to_pt = v(k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               v(k-2), v(k-1), v(k), v(k+1), &
               interp_val_to_pt, str, ierr)
            if (ierr == 0) return
            write(*,1) '0.5d0*(dq(k-2)+dq(k-1))', 0.5d0*(dq(k-2)+dq(k-1))
            write(*,1) '0.5d0*(dq(k-1)+dq(k))', 0.5d0*(dq(k-1)+dq(k))
            write(*,1) '0.5d0*(dq(k)+dq(k+1))', 0.5d0*(dq(k)+dq(k+1))
            write(*,2) 'dq(k-2)', k-2, dq(k-2)
            write(*,2) 'dq(k-1)', k-1, dq(k-1)
            write(*,2) 'dq(k)', k, dq(k)
            write(*,2) 'dq(k+1)', k+1, dq(k+1)

            stop 'interp_val_to_pt'
         endif
         interp_val_to_pt = (v(k)*dq(k-1) + v(k-1)*dq(k))/(dq(k-1) + dq(k))
      end function interp_val_to_pt


      real(dp) function interp_xa_to_pt(xa,j,k,sz,dq,str)
         use interp_1d_lib, only: interp_4_to_1
         real(dp), intent(in) :: xa(:,:), dq(:)
         character (len=*), intent(in) :: str
         integer, intent(in) :: j, k, sz
         integer :: ierr
         include 'formats'
         if (j == 0) then
            interp_xa_to_pt = 0
            return
         end if
         if (k == 1) then
            interp_xa_to_pt = xa(j,k)
            return
         end if
         if (k > 2 .and. k < sz) then
            ierr = 0
            call interp_4_to_1( &
               0.5d0*(dq(k-2)+dq(k-1)), &
               0.5d0*(dq(k-1)+dq(k)), &
               0.5d0*(dq(k)+dq(k+1)), &
               0.5d0*dq(k-2)+dq(k-1), &
               xa(j,k-2), xa(j,k-1), xa(j,k), xa(j,k+1), &
               interp_xa_to_pt, str, ierr)
            interp_xa_to_pt = min(1d0,max(0d0,interp_xa_to_pt))
            if (ierr == 0) return
         endif
         interp_xa_to_pt = (xa(j,k)*dq(k-1) + xa(j,k-1)*dq(k))/(dq(k-1) + dq(k))
         interp_xa_to_pt = min(1d0,max(0d0,interp_xa_to_pt))
      end function interp_xa_to_pt


      real(dp) function get_dtau1(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: kap, kap_min
         include 'formats'
         ierr = 0
         k = 1
         kap = s% opacity(1)
         get_dtau1 = s% dm(1)*kap/(4*pi*s% rmid(1)*s% rmid(1))
         if (is_bad(get_dtau1)) then
            ierr = -1
            if (.not. s% report_ierr) return
            k = 1
            write(*,2) 'get_dtau1', k, get_dtau1
            write(*,2) 's% dm(1)', k, s% dm(k)
            write(*,2) 's% opacity(1)', k, s% opacity(k)
            write(*,2) 's% rmid(1)', k, s% rmid(k)
            write(*,2) 's% r(1)', k, s% r(k)
            write(*,2) 's% r(2)', 2, s% r(2)
            stop 'get_dtau1'
         end if
      end function get_dtau1


      subroutine get_tau(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ! tau(k) is optical depth at outer boundary of cell k
         real(dp) :: dtau, dr, Z, kap_min, kap, dm_sum, L_sum
         integer :: k
         logical, parameter :: dbg = .true.
         include 'formats'
         ierr = 0
         dtau = get_dtau1(s, ierr)
         if (ierr /= 0) return
         s% tau(1) = s% tau_factor*s% tau_base
         s% lntau(1) = safe_log(s% tau(1))
         s% tau_start(1) = s% tau(1)
         dm_sum = 0
         L_sum = 0
         do k = 2, s% nz
            s% tau(k) = s% tau(k-1) + dtau
            s% lntau(k) = log(s% tau(k))
            if (s% tau_for_L_BB > 0d0 .and. s% L_for_BB_outer_BC <= 0d0) then
               dm_sum = dm_sum + s% dm(k-1)
               L_sum = L_sum + s% L(k-1)*s% dm(k-1)
               if (s% tau(k) >= s% tau_for_L_BB) then
                  s% L_for_BB_outer_BC = L_sum/dm_sum
                  !write(*,2) 's% L_for_BB_outer_BC', k, s% L_for_BB_outer_BC
               end if
            end if
            if (s% tau_start(k) < 0) s% tau_start(k) = s% tau(k)
            kap = s% opacity(k)
            dtau = s% dm(k)*kap/(4*pi*s% rmid(k)*s% rmid(k))
            if (is_bad(dtau)) then
               ierr = -1
               if (.not. s% report_ierr) return
               write(*,2) 'dtau', k, dtau
               write(*,2) 's% dm(k)', k, s% dm(k)
               write(*,2) 's% opacity(k)', k, s% opacity(k)
               write(*,2) 's% rmid(k)', k, s% rmid(k)
               stop 'get_tau'
            end if
            if (k == s% nz) s% tau_center = s% tau(k) + dtau
            !write(*,*) 'dtau, dlogtau', k, tau(k) - tau(k-1), &
            !   log10(tau(k)/tau(k-1))
         end do
         if (s% tau_for_L_BB > 0d0 .and. s% L_for_BB_outer_BC < 0d0) then
            write(*,1) 'failed to set s% L_for_BB_outer_BC', s% L_for_BB_outer_BC
            stop 'get_tau'
         end if
      end subroutine get_tau


      integer function find_cell_for_mass(s, m)
         type (star_info), pointer :: s
         real(dp), intent(in) :: m
         integer :: k
         find_cell_for_mass = s% nz
         do k = 1, s% nz-1
            if (s% m(k) >= m .and. m > s% m(k+1)) then
               find_cell_for_mass = k
               return
            end if
         end do
      end function find_cell_for_mass


      subroutine get_delta_Pg(s, nu_max, delta_Pg)
         type (star_info), pointer :: s
         real(dp), intent(in) :: nu_max ! microHz
         real(dp), intent(out) :: delta_Pg ! seconds
         ! g-mode period spacing for l=1
         real(dp) :: integral, N2, omega2, kr2, L2, el, &
            dr, r, r2, cs2, sl2, I_integral, I_integral_limit
         integer :: k, k_sl2
         logical, parameter :: dbg = .false.
         include 'formats'
         if (dbg) then
            write(*,2) 'nu_max', s% model_number, nu_max
            write(*,2) 's% star_mass', s% model_number, s% star_mass
            write(*,2) 's% photosphere_r', s% model_number, s% photosphere_r
            write(*,2) 's% Teff', s% model_number, s% Teff
         end if
         delta_Pg = 0
         if (.not. s% calculate_Brunt_N2) return
         integral = 0
         I_integral = 0
         I_integral_limit = 0.5d0
         omega2 = pow2(2*pi*nu_max/1d6)
         if (dbg) write(*,1) 'log omega2', log10(omega2)
         el = 1
         L2 = el*(el+1)
         k_sl2 = 0
         do k = 2, s% nz
            N2 = s% brunt_N2(k)
            r = s% r(k)
            r2 = r*r
            cs2 = s% csound_face(k)*s% csound_face(k)
            sl2 = L2*cs2/r2
            dr = s% rmid(k-1) - s% rmid(k)
            if (omega2 >= sl2) then
               cycle
            end if
            if (k_sl2 == 0) then
               k_sl2 = k
               if (dbg) write(*,2) 'k_sl2', k
            end if
            if (N2 > omega2) then ! in g-cavity
               if (dbg .and. integral == 0) write(*,2) 'enter g-cavity', k
               integral = integral + sqrt(N2)*dr/r
            else ! in decay region
               if (integral == 0) cycle ! ! haven't been in g-cavity yet
               if (dbg .and. I_integral == 0) write(*,2) 'enter decay', k
               ! in decay region below g-cavity; I_integral estimates decay
               kr2 = (1 - n2/omega2)*(1 - Sl2/omega2)*omega2/cs2
               I_integral = I_integral + sqrt(-kr2)*dr
               if (I_integral > I_integral_limit) exit
            end if
         end do

         if (dbg) write(*,2) 'omega2 nu_max integral I_integral', &
            s% model_number, omega2, nu_max, integral, I_integral

         if (integral == 0) return
         delta_Pg = sqrt(2d0)*pi*pi/integral
         if (is_bad(delta_Pg)) delta_Pg = 0

         if (dbg) write(*,2) 'delta_Pg', s% model_number, delta_Pg

      end subroutine get_delta_Pg


      subroutine set_rmid(s, nzlo, nzhi, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr
         integer :: k, nz
         real(dp) :: r003, rp13, rmid, rmid2
         logical :: dbg
         include 'formats'
         ierr = 0
         nz = s% nz
         dbg = .false.
         if (s% RSP_flag) then ! mid = avg r
            !dbg = s% model_number >= s% max_model_number - 1
            do k=nzlo, nzhi
               if (k < nz) then
                  rmid = 0.5d0*(s% r(k) + s% r(k+1))
               else
                  rmid = 0.5d0*(s% r(k) + s% R_center)
               end if
               s% rmid(k) = rmid
               if (dbg) write(*,3) 'set_rmid s% r(k)', k, s% model_number, s% r(k)
               if (dbg) write(*,3) 'set_rmid s% rmid(k)', k, s% model_number, s% rmid(k)
               if (s% rmid_start(k) < 0) s% rmid_start(k) = s% rmid(k)
               rmid2 = rmid*rmid
               s% drmid_dlnR00(k) = 0.5d0*s% r(k)
               s% drmid2_dlnR00(k) = 2d0*rmid*s% drmid_dlnR00(k)
               if (k < nz) then
                  s% drmid_dlnRp1(k) = 0.5d0*s% r(k+1)
                  s% drmid2_dlnRp1(k) = 2d0*rmid*s% drmid_dlnRp1(k)
               else
                  s% drmid_dlnRp1(k) = 0d0
                  s% drmid2_dlnRp1(k) = 0d0
               end if
            end do
            return
         end if
         ! mid = middle by volume
         do k=nzlo, nzhi
            r003 = s% r(k)*s% r(k)*s% r(k)
            if (k < nz) then
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
            else
               rp13 = s% R_center*s% R_center*s% R_center
            end if
            rmid = pow(0.5d0*(r003 + rp13),1d0/3d0)
            s% rmid(k) = rmid
            if (s% rmid_start(k) < 0) s% rmid_start(k) = s% rmid(k)
            rmid2 = rmid*rmid
            s% drmid_dlnR00(k) = 0.5d0*r003/rmid2
            s% drmid2_dlnR00(k) = r003/rmid
            if (k < nz) then
               s% drmid_dlnRp1(k) = 0.5d0*rp13/rmid2
               s% drmid2_dlnRp1(k) = rp13/rmid
            else
               s% drmid_dlnRp1(k) = 0d0
               s% drmid2_dlnRp1(k) = 0d0
            end if
         end do
      end subroutine set_rmid


      real(dp) function get_tau_at_r(s, r, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: r
         integer, intent(out) :: ierr
         real(dp) :: dtau, dr, tau_m1, tau_00
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         dtau = get_dtau1(s, ierr)
         if (ierr /= 0) return
         tau_00 = s% tau_factor*s% tau_base
         get_tau_at_r = tau_00
         if (r >= s% r(1)) return
         do k = 2, s% nz
            tau_m1 = tau_00
            tau_00 = tau_m1 + dtau
            if (r < s% r(k-1) .and. r >= s% r(k)) then
               get_tau_at_r = &
                  (tau_00*(s% r(k-1)-r) + tau_m1*(r-s% r(k)))/(s% r(k-1)-s% r(k))
               return
            end if
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
         end do
      end function get_tau_at_r


      integer function find_tau_phot(s, tau00, taup1, ierr)
         ! return k for the cell containing optical depth = tau_base
         type (star_info), pointer :: s
         real(dp), intent(out) :: tau00, taup1
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: dtau, tau_phot

         include 'formats'
         ierr = 0
         tau00 = 0
         taup1 = 0
         find_tau_phot = 1
         if (s% tau_factor >= 1 .and. .not. s% RSP_flag) return
         tau_phot = s% tau_base
         tau00 = s% tau_factor*s% tau_base
         do k = 1, s% nz
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            if (taup1 >= tau_phot) then
               find_tau_phot = k
               return
            end if
            tau00 = taup1
         end do
         ierr = -1
      end function find_tau_phot

      
      real(dp) function get_r_phot(s)
         type (star_info), pointer :: s  
         real(dp) :: r, m, v, L, T_phot, cs, kap, logg, ysum
         integer :: k_phot
         call get_phot_info(s,r,m,v,L,T_phot,cs,kap,logg,ysum,k_phot)
         get_r_phot = r
      end function get_r_phot


      subroutine get_phot_info(s,r,m,v,L,T_phot,cs,kap,logg,ysum,k_phot)
         type (star_info), pointer :: s
         real(dp), intent(out) :: r, m, v, L, T_phot, cs, kap, logg, ysum
         integer, intent(out) :: k_phot

         integer :: k
         real(dp) :: tau00, taup1, dtau, r003, rp13, r3, tau_phot, &
            Tface_0, Tface_1

         include 'formats'

         tau00 = 0
         taup1 = 0
         ysum = 0
         r = s% r(1)
         m = s% m(1)
         if (s% u_flag) then
            v = s% u(1)
         else if (s% v_flag) then
            v = s% v(1)
         else
            v = 0d0
         end if
         L = max(1d0, s% L(1)) ! don't use negative L(1)
         T_phot = s% T(1)
         cs = s% csound(1)
         kap = s% opacity(1)
         logg = safe_log10(s% cgrav(1)*m/(r*r))
         k_phot = 1
         tau_phot = s% tau_base
         tau00 = s% tau_factor*s% tau_base
         if (tau00 >= tau_phot) return
         do k = 1, s% nz-1
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            ysum = ysum + s% rho(k)*(s% r(k) - s% r(k+1))
            if (taup1 >= tau_phot .and. dtau > 0d0) then
               if (k == 1) then
                  Tface_0 = s% T(k)
               else
                  Tface_0 = 0.5d0*(s% T(k) + s% T(k-1))
               end if
               Tface_1 = 0.5d0*(s% T(k) + s% T(k+1))
               T_phot = Tface_0 + (Tface_1 - Tface_0)*(tau_phot - tau00)/dtau
               r003 = s% r(k)*s% r(k)*s% r(k)
               rp13 = s% r(k+1)*s% r(k+1)*s% r(k+1)
               r3 = r003 + (rp13 - r003)*(tau_phot - tau00)/dtau
               r = pow(r3,1d0/3d0)
               m = s% m(k) - s% dm(k)*(tau_phot - tau00)/dtau
               if (s% u_flag) then
                  v = s% u_face(k) + &
                     (s% u_face(k+1) - s% u_face(k))*(tau_phot - tau00)/dtau
                  if (is_bad(v)) then
                     write(*,2) 'v', k, v
                     write(*,2) 's% u_face(k)', k, s% u_face(k)
                     write(*,2) 's% u_face(k+1)', k, s% u_face(k+1)
                     stop 'get_phot_info'
                  end if
               else if (s% v_flag) then
                  v = s% v(k) + (s% v(k+1) - s% v(k))*(tau_phot - tau00)/dtau
               end if
               L = s% L(k) + (s% L(k+1) - s% L(k))*(tau_phot - tau00)/dtau
               k_phot = k
               cs = s% csound(k_phot)
               kap = s% opacity(k_phot)
               logg = safe_log10(s% cgrav(k_phot)*m/(r*r))
               return
            end if
            tau00 = taup1
         end do
         !write(*,*) 'get_phot_info failed to find photosphere'
         k_phot = s% nz
         r = s% R_center
         m = s% m_center
         v = s% v_center
         T_phot = s% T(k_phot)
         L = max(1d0, s% L_center)
         cs = s% csound(k_phot)
         kap = s% opacity(k_phot)
         logg = safe_log10(s% cgrav(k_phot)*m/(r*r))
      end subroutine get_phot_info


      real(dp) function get_phot_kap(s)
         type (star_info), pointer :: s

         integer :: k
         real(dp) :: tau00, taup1, dtau, tau_phot

         include 'formats'

         get_phot_kap = s% opacity(1)
         tau_phot = s% tau_base
         tau00 = s% tau_factor*s% tau_base
         if (tau00 >= tau_phot) return
         do k = 1, s% nz-1
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            taup1 = tau00 + dtau
            if (taup1 >= tau_phot .and. dtau > 0d0) then
               get_phot_kap = s% opacity(k)
               return
            end if
            tau00 = taup1
         end do
         get_phot_kap = s% opacity(s% nz)
      end function get_phot_kap


      real(dp) function center_value(s, p)
         type (star_info), pointer :: s
         real(dp), intent(in) :: p(:)
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = p(k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x + dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_value = sum_x/sum_dq
      end function center_value


      subroutine interp_q( &
            nz2, nvar_hydro, species, qval, xh, xa, q, dq, struct, comp, ierr)
         use num_lib, only: binary_search
         integer, intent(in) :: nz2, nvar_hydro, species
         real(dp), intent(in) :: qval
         real(dp), intent(in) :: xh(:,:), xa(:,:), q(:), dq(:)
         real(dp), intent(inout) :: struct(:), comp(:)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: alfa
         ierr = 0
         if (qval <= q(nz2)) then
            if (nvar_hydro > 0) &
               struct(1:nvar_hydro) = xh(1:nvar_hydro,nz2)
            if (species > 0) &
               comp(1:species) = xa(1:species,nz2)
            return
         end if
         k = binary_search(nz2, q, 0, qval)
         if (k < 1 .or. k >= nz2) then
            ierr = -1
            return
         end if
         if (qval <= q(k) .and. qval > q(k+1)) then
            alfa = (qval - q(k+1)) / dq(k)
            if (nvar_hydro > 0) &
               struct(1:nvar_hydro) = &
                  alfa*xh(1:nvar_hydro,k) + (1-alfa)*xh(1:nvar_hydro,k+1)
            if (species > 0) &
               comp(1:species) = alfa*xa(1:species,k) + (1-alfa)*xa(1:species,k+1)
            return
         end if
         ierr = -1
      end subroutine interp_q


      subroutine std_write_internals_to_file(id, num)
         use utils_lib, only : mkdir
         integer, intent(in) :: num, id
         character (len=strlen) :: fname
         integer :: ierr
         ierr = 0
         call mkdir('plot_data')
         write(fname, '(a, i1, a)') 'plot_data/internals', mod(abs(num), 10), '.data'
         write(*,*) 'call write_internals_to_file ' // trim(fname)
         call write_internals_to_file(id, fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in write_internals_to_file ' // trim(fname)
         end if
      end subroutine std_write_internals_to_file


      subroutine write_internals_to_file(id, filename, ierr)
         use utils_lib
         character (len=*), intent(in) :: filename
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer :: iounit
         ierr = 0
         open(newunit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr == 0) then
            call write_internals(id, iounit, ierr)
            close(iounit)
         else
            write(*, *) 'failed to open internals file ' // trim(filename)
         end if
      end subroutine write_internals_to_file


      subroutine write_internals(id, iounit, ierr)
         use chem_def
         integer, intent(in) :: iounit, id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'write_internals: get_star_ptr ierr', ierr
            return
         end if
         call write_model_info(s, iounit, ierr)
      end subroutine write_internals


      subroutine write_model_info(s, iounit, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer, pointer :: chem_id(:)
         integer :: k, i, nz, species
         integer :: he4

         ierr = 0
         
         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         he4 = s% net_iso(ihe4)

         write(iounit,'(a)') '            mass         initial_z       n_shells'
         write(iounit,'(2x,2f15.4,i15)') s% star_mass, s% initial_z, nz
         write(iounit,fmt='(i5)',advance='no') 1
         do i=2,88
            write(iounit,fmt='(i12,15x)',advance='no') i
         end do
         write(iounit,*)
         write(iounit,fmt='(a5,1x,99(a26,1x))',advance='no') &
               'grid', 'r', 'm', 'log_dq', &
               'log10d', 'log10T', 'log10m', 'log10r', 'L', 'r_div_rstar', &
               'log10P', 'log10Pgas', 'chiT', 'chiRho', &
               'dlnRho_dlnPgas_const_T', 'dlnRho_dlnT_const_Pgas', &
               'extra5', 'extra6', 'extra7', 'extra8'
         ! composition info
         do i=1, species
            write(iounit, fmt='(a26, 1x)', advance='no') trim(chem_isos% name(chem_id(i)))
         end do
         do i=1, species
            write(iounit, fmt='(a26, 1x)', advance='no') 'lg_' // trim(chem_isos% name(chem_id(i)))
         end do
         write(iounit,fmt=*)

         do k=1, nz
            write(iounit,'(i5,1x,99(1pe26.16,1x))',advance='no') k,  &
               s% r(k)/Rsun, s% m(k)/Msun, safe_log10(s% dq(k)), &
               s% lnd(k)/ln10, s% lnT(k)/ln10, log10(s% m(k)),  &
               s% lnR(k)/ln10, s% L(k)/Lsun, s% r(k)/s% r(1), &
               s% lnP(k)/ln10, s% lnPgas(k)/ln10, s% chiT(k), s% chiRho(k), &
               s% dlnRho_dlnPgas_const_T(k), s% dlnRho_dlnT_const_Pgas(k), &
               s% profile_extra(k,5), s% profile_extra(k,6), &
               s% profile_extra(k,7), s% profile_extra(k,8)
            do i=1, species
               write(iounit, fmt='(1pe26.16, 1x)', advance='no') s% xa(i, k)
            end do
            do i=1, species
               write(iounit, fmt='(1pe26.16, 1x)', advance='no') safe_log10(s% xa(i, k))
            end do
            write(iounit,*)
         end do

      end subroutine write_model_info


      subroutine dump_model_info_for_ndiff(s, iounit, ierr)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: iounit
         integer, intent(out) :: ierr
         integer, pointer :: chem_id(:)
         integer :: k, j, nz, species
         include 'formats'
         ierr = 0
         nz = s% nz
         species = s% species
         chem_id => s% chem_id
         write(iounit,*) 'nz', nz
         write(iounit,1) 'star_mass', s% star_mass
         write(iounit,1) 'initial_z', s% initial_z
         do k=1, nz
            do j=1, s% nvar_hydro
               write(iounit,2) trim(s% nameofvar(j)), k, s% xh(j,k)
            end do
            do j=1,species
               write(iounit,2) trim(chem_isos% name(chem_id(j))), k, clip(s% xa(j,k))
            end do
         end do

         contains

         real(dp) function clip(x)
            real(dp), intent(in) :: x
            if (.true. .or. x > 1d-30) then
               clip = x
            else
               clip = 0d0
            end if
         end function clip

      end subroutine dump_model_info_for_ndiff

      
      subroutine set_abs_du_div_cs(s)
         type (star_info), pointer :: s
         
         integer :: k, nz, j
         real(dp) :: abs_du, cs
         include 'formats'
         nz = s% nz

         if (s% u_flag) then
            do k=2,nz-1
               abs_du = &
                  max(abs(s% u_start(k) - s% u_start(k+1)), &
                      abs(s% u_start(k) - s% u_start(k-1)))
               cs = maxval(s% csound(max(1,k-5):min(nz,k+5)))
               s% abs_du_plus_cs(k) = abs_du + cs
               s% abs_du_div_cs(k) = abs_du/cs
            end do
            k = nz
            s% abs_du_plus_cs(k) = &
               abs(s% u_start(k) - s% u_start(k-1)) + s% csound_start(k)
            s% abs_du_div_cs(k) = &
               abs(s% u_start(k) - s% u_start(k-1))/s% csound_start(k)
            k = 2
            s% abs_du_plus_cs(k) = &
               abs(s% u_start(k) - s% u_start(k+1)) + s% csound_start(k)
            s% abs_du_div_cs(k) = &
               abs(s% u_start(k) - s% u_start(k+1))/s% csound_start(k)
            k = 1
            s% abs_du_plus_cs(k) = s% abs_du_plus_cs(k+1)
            s% abs_du_div_cs(k) = s% abs_du_div_cs(k+1)
            do j = 1,3
               do k=2,nz-1
                  s% abs_du_div_cs(k) = sum(s% abs_du_div_cs(k-1:k+1))/3d0
               end do
            end do
         else
            do k=1,nz
               s% abs_du_plus_cs(k) = 1d99
               s% abs_du_div_cs(k) = 1d99
            end do
         end if
      
      end subroutine set_abs_du_div_cs


      real(dp) function rsi_div_rsimelt(s,k,species) 
         ! rsi = ion density parameter for cell k
         ! rsimelt = ion density parameter of quantum melting
         ! rsi < rsimelt => liquid, independent of T
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: k, species
         
         integer :: IX, j
         real(dp), dimension(species) :: AZion, ACMI, AY
         real(dp) :: Y, CMImean, Z73, RS, RSI
         real(dp), parameter :: RSIMELT=140d0, TINY=1d-7, &
            AUM=1822.888d0 ! a.m.u./m_e
         
         include 'formats'
         
         ! details from eos/private/pc_eos.f
         
         AZion(1:species) = chem_isos% Z(s% chem_id(1:species))
         ACMI(1:species) = chem_isos% W(s% chem_id(1:species))
         do j=1,species
            if (s% xa(j,k) < s% eos_rq% mass_fraction_limit_for_PC) then
               AY(j) = 0
            else
               AY(j) = s% xa(j,k)/ACMI(j)
            end if
         end do
         
         Y=0.d0
         do IX=1,species
            Y=Y+AY(IX)
         end do
         if (dabs(Y-1.d0).gt.TINY) then
           do IX=1,species
              AY(IX)=AY(IX)/Y
           end do
         end if

         CMImean=0.d0
         Z73=0.d0
         do IX=1,species
            if (AY(IX) < TINY) cycle
            Z73 = Z73 + AY(IX)*pow(AZion(IX),7d0/3d0)
            CMImean = CMImean + AY(IX)*ACMI(IX)
         end do

         RS=pow(0.75d0/PI/s% rho(k),1d0/3d0)
         RSI=RS*CMImean*Z73*AUM
         
         if (is_bad(RSI)) then
            write(*,2) 'RSI', k, RSI
            write(*,2) 'Z73', k, Z73
            write(*,2) 'CMImean', k, CMImean
            write(*,2) 'RS', k, RS
            write(*,2) 's% rho(k)', k, s% rho(k)
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            stop 'rsi_div_rsimelt'
         end if
         
         rsi_div_rsimelt = RSI/RSIMELT
         
      end function rsi_div_rsimelt


      subroutine get_shock_info(s)
         type (star_info), pointer :: s
         integer :: k, nz, kk, kmin, k0, k1, ierr
         real(dp) :: v_div_cs_00, v_div_cs_m1, rmax, vmax, alfa
         real(dp), pointer :: v(:)

         include 'formats'
         
         if (s% u_flag) then
            v => s% u ! may not have u_face
         else if (s% v_flag) then
            v => s% v
         else
            return
         end if

         nz = s% nz
         kmin = nz
         do k=1,nz-1
            if (s% q(k) <= s% max_q_for_outer_mach1_location) then
               kmin = k
               exit
            end if
         end do

         s% shock_velocity = 0
         s% shock_csound = 0
         s% shock_lgT = 0
         s% shock_lgRho = 0
         s% shock_lgP = 0
         s% shock_mass = 0
         s% shock_q = 0
         s% shock_radius = 0
         s% shock_gamma1 = 0
         s% shock_entropy = 0
         s% shock_tau = 0
         s% shock_k = 0
         s% shock_pre_lgRho = 0
         
         if (kmin < nz) then ! search inward for 1st shock moving outward
            v_div_cs_00 = v(kmin)/s% csound_face(kmin)
            do k = kmin+1,nz
               v_div_cs_m1 = v_div_cs_00
               v_div_cs_00 = v(k)/s% csound_face(k)
               if (v_div_cs_00 > 0 .and. v_div_cs_m1 > 0 .and. &
                     v_div_cs_00 >= 1d0 .and. v_div_cs_m1 < 1d0) then
                  do kk = k+1, nz ! search inward for local max v
                     if (v(kk) <= v(kk-1)) then
                        s% shock_k = kk - 1
                        exit
                     end if
                  end do
                  exit
               end if
            end do
         end if

         k = s% shock_k
         if (k < nz .and. k > 1) then
            if (v(k) > max(v(k-1),v(k+1))) then
               rmax = s% r(k)
               vmax = v(k)
               s% shock_velocity = vmax
               s% shock_radius = rmax/Rsun
               if (rmax <= s% r(k)) then
                  k0 = k-1
                  k1 = k
               else if (rmax <= s% r(k+1)) then
                  k0 = k
                  k1 = k+1
               else
                  k0 = k+1
                  k1 = k+2
               end if
               alfa = (s% r(k0) - rmax)/(s% r(k0) - s% r(k1))
               s% shock_q = s% q(k0) - alfa*s% dq(k0)
               s% shock_mass = (s% m(k0) - alfa*s% dm(k0))/Msun
               s% shock_csound = s% csound(k0)
               s% shock_lgT = s% lnT(k0)/ln10
               s% shock_lgRho = s% lnd(k0)/ln10
               s% shock_lgP = s% lnP(k0)/ln10
               s% shock_gamma1 = s% gamma1(k0)
               s% shock_entropy = s% entropy(k0)
               s% shock_tau = s% tau(k0)
               s% shock_k = k0
               do kk=k0-1,1,-1
                  if (v(kk) < 0.1d0*s% csound_face(kk)) then
                     s% shock_pre_lgRho = s% lnd(kk)/ln10
                     exit
                  end if
               end do
            end if
         end if

      end subroutine get_shock_info


      real(dp) function min_dr_div_cs(s,min_k) ! seconds
         type (star_info), pointer :: s
         integer, intent(out) :: min_k
         integer :: k, nz, j, k_min
         real(dp) :: dr, dt, D, abs_du, cs, min_q, max_q, &
            min_abs_du_div_cs, r00, rp1, dr_div_cs
         include 'formats'
         nz = s% nz
         min_k = nz
         min_dr_div_cs = 1d99
         if (s% v_flag) then
            r00 = s% R_center
            do k=nz,1,-1
               rp1 = r00
               r00 = s% r(k)
               dr_div_cs = (r00 - rp1)/s% csound(k)
               if (dr_div_cs < min_dr_div_cs) then
                  min_dr_div_cs = dr_div_cs
                  min_k = k
               end if
            end do
            return
         end if
         if (.not. s% u_flag) return
         min_abs_du_div_cs = &
            s% min_abs_du_div_cs_for_dt_div_min_dr_div_cs_limit
         min_q = s% min_q_for_dt_div_min_dr_div_cs_limit
         max_q = s% max_q_for_dt_div_min_dr_div_cs_limit
         k_min = max(1, s% min_k_for_dt_div_min_dr_div_cs_limit)
         do k = k_min, nz-1
            if (s% q(k) > max_q) cycle
            if (s% q(k) < min_q) exit
            if (s% abs_du_div_cs(k) < min_abs_du_div_cs) cycle
            dr = s% r(k) - s% r(k+1)
            dt = dr/s% abs_du_plus_cs(k)
            if (dt < min_dr_div_cs) then
               min_dr_div_cs = dt
               min_k = k
            end if
         end do
      end function min_dr_div_cs


      ! largest k s.t. for all k' < k, cell k' has Cp(k')*T(k')*mstar_dot < L(k).
      subroutine set_k_CpTMdot_lt_L(s)
         type (star_info), pointer :: s
         integer :: k, nz
         if (s% mstar_dot <= 0d0 .or. s% gamma_law_hydro > 0d0) then
            s% k_CpTMdot_lt_L = 1
            return
         end if
         nz = s% nz
         do k = 2, nz
            if (s% Cp(k)*s% T(k)*s% mstar_dot >= max(1d-99,s% L(k))) then
               s% k_CpTMdot_lt_L = k-1
               return
            end if
         end do
         s% k_CpTMdot_lt_L = nz
      end subroutine set_k_CpTMdot_lt_L


      subroutine reset_epsnuc_vectors(s)
         type (star_info), pointer :: s
         integer :: k, nz
         nz = s% nz
         do k=1,s% nz
            s% eps_nuc(k) = 0d0
            s% d_epsnuc_dlnd(k) = 0d0
            s% d_epsnuc_dlnT(k) = 0d0
            s% d_epsnuc_dx(:,k) = 0d0
            s% eps_nuc_categories(:,k) = 0d0
            s% dxdt_nuc(:,k) =  0d0
            s% d_dxdt_nuc_dRho(:,k) =  0d0
            s% d_dxdt_nuc_dT(:,k) =  0d0
            s% d_dxdt_nuc_dx(:,:,k) =  0d0
            s% eps_nuc_neu_total(k) = 0d0
         end do
      end subroutine reset_epsnuc_vectors


      subroutine reset_starting_vectors(s)
         type (star_info), pointer :: s
         integer :: k, nz
         nz = s% nz
         do k=1,s% nz
            s% T_start(k) = -1d99
            s% r_start(k) = -1d99
            s% rmid_start(k) = -1d99
            s% v_start(k) = -1d99
            s% u_start(k) = -1d99
            s% lnd_start(k) = -1d99
            s% lnT_start(k) = -1d99
            s% csound_start(k) = -1d99
            s% eta_visc_start(k) = -1d99
            s% rho_start(k) = -1d99
            s% tau_start(k) = -1d99
            s% erad_start(k) = -1d99
            s% alpha_RTI_start(k) = -1d99
            s% w_start(k) = -1d99
            s% dPdr_dRhodr_info(k) = -1d99
         end do
      end subroutine reset_starting_vectors
      
      
      subroutine store_partials(s, k, xscale, i_eqn, nvar, d_dm1, d_d00, d_dp1)
         type (star_info), pointer :: s
         integer, intent(in) :: k, i_eqn, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), intent(in) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         integer :: nz, j
         nz = s% nz
         do j=1,nvar
            if (k > 1) call em1(s, xscale, i_eqn, j, k, nvar, d_dm1(j))
            call e00(s, xscale, i_eqn, j, k, nvar, d_d00(j))
            if (k < nz) call ep1(s, xscale, i_eqn, j, k, nvar, d_dp1(j))
         end do            
      end subroutine store_partials


      subroutine unpack_res18_partials(s, k, nvar, xscale, i_eqn, &
            res18, d_dm1, d_d00, d_dp1)
         use auto_diff
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar, i_eqn
         real(dp), pointer :: xscale(:,:)
         type(auto_diff_real_18var_order1) :: res18
         real(dp) :: d_dm1(nvar), d_d00(nvar), d_dp1(nvar)
         
         real(dp) :: val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
            deturb_m1, deturb_00, deturb_p1, dlnR_m1, dlnR_00, dlnR_p1, &
            dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1
         integer :: j

         include 'formats'

         call unwrap(res18, val, dlnd_m1, dlnd_00, dlnd_p1, dlnT_m1, dlnT_00, dlnT_p1, &
                     deturb_m1, deturb_00, deturb_p1, dlnR_m1, dlnR_00, dlnR_p1, &
                     dv_m1, dv_00, dv_p1, dL_m1, dL_00, dL_p1) 
                     
         call unpack1(s% i_lnd, dlnd_m1, dlnd_00, dlnd_p1)
         call unpack1(s% i_lnT, dlnT_m1, dlnT_00, dlnT_p1)
         call unpack1(s% i_lnR, dlnR_m1, dlnR_00, dlnR_p1)
         if (s% i_v /= 0) call unpack1(s% i_v, dv_m1, dv_00, dv_p1)
         if (s% i_lum /= 0) call unpack1(s% i_lum, dL_m1, dL_00, dL_p1)
         if (s% i_eturb /= 0) call unpack1(s% i_eturb, deturb_m1, deturb_00, deturb_p1)
         
         contains
         
         subroutine unpack1(j, dvar_m1, dvar_00, dvar_p1)
            integer, intent(in) :: j
            real(dp), intent(in) :: dvar_m1, dvar_00, dvar_p1
            d_dm1(j) = dvar_m1
            d_d00(j) = dvar_00
            d_dp1(j) = dvar_p1
         end subroutine unpack1         
         
      end subroutine unpack_res18_partials


      subroutine set_scale_height(s)
         type (star_info), pointer :: s
         real(dp) :: Hp, alt_Hp, alfa, beta, rho_face, P_face
         integer :: k
         include 'formats'
         !if (s% gamma_law_hydro > 0d0) return
         do k=1,s% nz
            if (s% cgrav(k) == 0) then
               s% scale_height(k) = 0
               cycle
            end if
            if (k == 1) then
               alfa = 1
            else
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            end if
            beta = 1 - alfa
            if (alfa == 1) then
               rho_face = s% rho(k)
               P_face = s% P(k)
            else
               rho_face = alfa*s% rho(k) + beta*s% rho(k-1)
               P_face = alfa*s% P(k) + beta*s% P(k-1)
            end if
            Hp = P_face/(rho_face*s% grav(k))
            if (s% cgrav(k) <= 0) then
               alt_Hp = s% r(k)
            else
               alt_Hp = sqrt(P_face / s% cgrav(k)) / rho_face
            end if
            s% scale_height(k) = min(Hp, alt_Hp)
         end do
      end subroutine set_scale_height


      real(dp) function tau_eff(s,k)
         ! tau_eff = tau that gives the local P == P_atm if this location at surface
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: P, g, Pextra_factor
         if (k == 1) then
            tau_eff = s% tau(1)
            return
         end if
         if (s% cgrav(k) <= 0d0) then
            tau_eff = 0d0
            return
         end if
         P = (s% dq(k-1)*s% P(k) + s% dq(k)*s% P(k-1))/(s% dq(k-1) + s% dq(k))
         g = s% cgrav(k)*s% m_grav(k)/(s% r(k)*s% r(k))
         Pextra_factor = s% Pextra_factor
         tau_eff = s% opacity(k)*(P/g - &
               Pextra_factor*(s% L(k)/s% m_grav(k))/(6d0*pi*clight*s% cgrav(k)))
      end function tau_eff


      real(dp) function eval_Ledd(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp) :: dtau1, dtau, dr, tau, dqsum, Ledd_sum
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         eval_Ledd = 0d0
         if (s% cgrav(1) <= 0d0) return
         dtau1 = get_dtau1(s, ierr)
         if (ierr /= 0) return
         dtau = dtau1
         tau = s% tau_factor*s% tau_base
         dqsum = s% dq(1)
         Ledd_sum = s% dq(1)*4*pi*clight*s% cgrav(1)*s% m_grav(1)/s% opacity(1)
         do k = 2, s% nz
            tau = tau + dtau
            if (tau > s% surf_avg_tau) exit
            dtau = s% dm(k)*s% opacity(k)/(4*pi*s% rmid(k)*s% rmid(k))
            dqsum = dqsum + s% dq(k)
            Ledd_sum = Ledd_sum + &
               s% dq(k)*4*pi*clight*s% cgrav(1)*s% m_grav(1)/s% opacity(k)
         end do
         eval_Ledd = Ledd_sum/dqsum
      end function eval_Ledd


      real(dp) function eval_min_cell_collapse_time(s,k_lo,k_hi,min_collapse_k,ierr) &
            result(min_collapse_time)
         type (star_info), pointer :: s
         integer, intent(in) :: k_lo, k_hi
         integer, intent(out) :: min_collapse_k, ierr
         real(dp) :: rp1, vp1, r00, v00, time
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         rp1 = s% R_center
         vp1 = s% v_center
         min_collapse_time = 1d99
         min_collapse_k = -1
         if (.not. s% v_flag) return
         do k = k_hi, k_lo, -1
            v00 = s% v(k)
            r00 = s% r(k)
            if (r00 <= rp1) then
               !ierr = -1
               min_collapse_time = -1
               min_collapse_k = -1
               return
               write(*,2) 'bad radii', k, r00, rp1
               stop 'eval_min_cell_collapse_time'
            end if
            if (vp1 > v00) then
               time = (r00 - rp1)/(vp1 - v00)
               if (time < min_collapse_time) then
                  min_collapse_time = time
                  min_collapse_k = k
               end if
            end if
            rp1 = r00
            vp1 = v00
         end do
      end function eval_min_cell_collapse_time


      real(dp) function total_angular_momentum(s) result(J)
         type (star_info), pointer :: s
         include 'formats'
         if (.not. s% rotation_flag) then
            J = 0
         else
            J = dot_product(s% dm_bar(1:s% nz), s% j_rot(1:s% nz))
         end if
      end function total_angular_momentum


      real(dp) function eval_irradiation_heat(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: irradiation_dq, xq, eps
         eval_irradiation_heat = 0
         if (s% irradiation_flux /= 0) then
            irradiation_dq = 4*pi*s% r(1)*s% r(1)*s% column_depth_for_irradiation/s% xmstar
            xq = 1 - s% q(k)
            if (irradiation_dq > xq) then ! add irradiation heat for cell k
               eps = 0.25d0 * s% irradiation_flux / s% column_depth_for_irradiation
               if (irradiation_dq < xq + s% dq(k)) then ! only part of cell gets heated
                  eval_irradiation_heat = eps*(irradiation_dq - xq)/s% dq(k)
               else ! all of cell gets heated
                  eval_irradiation_heat = eps
               end if
            end if
         end if
      end function eval_irradiation_heat


      subroutine start_time(s, time0, total_all_before)
         type (star_info), pointer :: s
         integer(8), intent(out) :: time0
         real(dp), intent(out) :: total_all_before
         integer(8) :: clock_rate
         if (.not. s% doing_timing) return
         total_all_before = total_times(s)
         call system_clock(time0,clock_rate)
      end subroutine start_time


      subroutine update_time(s, time0, total_all_before, total)
         type (star_info), pointer :: s
         integer(8), intent(in) :: time0
         real(dp), intent(in) :: total_all_before
         real(dp), intent(inout) :: total
         real(dp) :: total_all_after, other_stuff
         integer(8) :: time1, clock_rate
         if (.not. s% doing_timing) return
         call system_clock(time1,clock_rate)
         total_all_after = total_times(s)
         other_stuff = total_all_after - total_all_before
            ! don't double count any other stuff
         total = total + (dble(time1-time0)/clock_rate - other_stuff)
      end subroutine update_time


      real(dp) function total_times(s)
         type (star_info), pointer :: s
         total_times = &
            s% time_evolve_step + &
            s% time_remesh + &
            s% time_adjust_mass + &
            s% time_conv_premix + &
            s% time_element_diffusion + &
            s% time_struct_burn_mix + &
            s% time_solver_matrix + &
            s% time_solve_mix + &
            s% time_solve_burn + &
            s% time_solve_omega_mix + &
            s% time_eos + &
            s% time_neu_kap + &
            s% time_nonburn_net + &
            s% time_mlt + &
            s% time_set_hydro_vars + &
            s% time_set_mixing_info
      end function total_times



      subroutine smooth(dc, sz)
         real(dp), intent(inout) :: dc(:)
         integer, intent(in) :: sz
         integer :: k
         k = 1
         dc(k) = (3*dc(k) + dc(k+1))/4
         do k=2,sz-1
            dc(k) = (dc(k-1) + 2*dc(k) + dc(k+1))/4
         end do
         k = sz
         dc(k) = (dc(k-1) + 3*dc(k))/4
      end subroutine smooth


      subroutine do_boxcar_mixing( &
            s, min_mass, max_mass, boxcar_nominal_mass, number_iterations, ierr)
         ! code based on SNEC
         type (star_info), pointer :: s
         real(dp), intent(in) :: min_mass, max_mass, boxcar_nominal_mass
         integer, intent(in) :: number_iterations
         integer, intent(out) :: ierr
         integer :: nz, iter, k_inner, k_outer, j, k
         real(dp) :: dm_sum, target_mass, actual_mass, xa, min_m, max_m
         include 'formats'
         ierr = 0
         nz = s% nz
         target_mass = boxcar_nominal_mass*Msun
         min_m = min_mass*Msun
         max_m = max_mass*Msun
         write(*,2) 'iters min max box', number_iterations, min_mass, max_mass, boxcar_nominal_mass
         do iter = 1, number_iterations
            do k_inner = nz, 1, -1
               if (s% m(k_inner) < min_m) cycle
               dm_sum = 0
               k_outer = 1
               do k = k_inner, 1, -1
                  dm_sum = dm_sum + s% dm(k)
                  if (dm_sum > target_mass) then
                     k_outer = k
                     exit
                  end if
               end do
               if (s% m(k_outer) > max_m) cycle
               actual_mass = dm_sum
               do j=1,s% species
                  xa = 0d0
                  do k=k_outer,k_inner
                     xa = xa + s% xa(j,k)*s% dm(k)
                  end do
                  do k=k_outer,k_inner
                     s% xa(j,k) = xa/dm_sum
                  end do
               end do
               if (actual_mass < target_mass) exit
            end do
         end do
         s% need_to_setvars = .true.
      end subroutine do_boxcar_mixing


      subroutine get_XYZ(s, xa, X, Y, Z)
         use chem_def, only: ih1, ih2, ihe3, ihe4
         type (star_info), pointer :: s
         real(dp), intent(in) :: xa(:)
         real(dp), intent(out) :: X, Y, Z
         X = 0d0
         if (s% net_iso(ih1) /= 0) X = X + xa(s% net_iso(ih1))
         if (s% net_iso(ih2) /= 0) X = X + xa(s% net_iso(ih2))
         X = min(1d0, max(0d0, X))
         Y = 0d0
         if (s% net_iso(ihe3) /= 0) Y = Y + xa(s% net_iso(ihe3))
         if (s% net_iso(ihe4) /= 0) Y = Y + xa(s% net_iso(ihe4))
         Y = min(1d0, max(0d0, Y))
         Z = min(1d0, max(0d0, 1d0 - (X + Y)))
      end subroutine get_XYZ


      subroutine get_face_values(s, v_mid, v_face, ierr)
         ! simple interpolation by mass
         type (star_info), pointer :: s
         real(dp), intent(in) :: v_mid(:)
         real(dp), intent(inout) :: v_face(:)
         integer, intent(out) :: ierr
         integer :: k
         real(dp) :: dq_sum
         ierr = 0
         v_face(1) = v_mid(1)
         do k=2, s% nz
            dq_sum = s% dq(k-1) + s% dq(k)
            v_face(k) = (v_mid(k)*s% dq(k-1) + v_mid(k-1)*s% dq(k))/dq_sum
         end do
      end subroutine get_face_values


      real(dp) function get_Ledd(s,k) result(Ledd)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: kap_face
         integer :: j
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         kap_face = interp_val_to_pt(s% opacity,j,s% nz,s% dq,'get_Ledd')
         Ledd = pi4*clight*s% cgrav(j)*s% m_grav(j)/kap_face
      end function get_Ledd


      real(dp) function get_Lrad(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         if (k == 1) then
            get_Lrad = s% L(k)
            return
         end if
         get_Lrad = s% L(k) - s% L_conv(k) ! L_conv set by last call on mlt
      end function get_Lrad


      real(dp) function get_Ladv(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: T, Erad, v, r
         T = s% T(k)
         Erad = crad*T*T*T*T
         if (s% u_flag) then
            v = s% u(k)
         else if (s% v_flag) then
            v = s% v(k)
         else
            v = s% r(1)*s% dlnR_dt(1)
         end if
         r = s% rmid(k)
         get_Ladv = 4*pi*r*r*v*Erad
      end function get_Ladv


      real(dp) function get_Lrad_div_Ledd(s,k) result(L_rad_div_Ledd)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer :: j
         real(dp) :: del_m, del_T4, area
         if (s% cgrav(k) <= 0) then
            L_rad_div_Ledd = 0d0
            return
         end if
         if (k == 1) then
            j = 2
         else
            j = k
         end if
         del_m = s% dm_bar(j)
         del_T4 = pow4(s% T(j-1)) - pow4(s% T(j))
         area = 4*pi*s% r(j)*s% r(j)
         L_rad_div_Ledd = &
            -(area*area*crad*(del_T4/del_m)/3)/(pi4*s% cgrav(j)*s% m_grav(j))
      end function get_Lrad_div_Ledd


      real(dp) function eval_rms_dvdt_div_v(s, klo, khi)
         type (star_info), pointer :: s
         integer, intent(in) :: klo, khi ! sum from klo to khi
         integer :: k
         real(dp) :: term, sum
         if (khi <= klo) then
            eval_rms_dvdt_div_v = 0d0
            return
         end if
         sum = 0
         do k=klo, khi
            term = s% dv_dt(k)/max(1d-50,abs(s% v(k)))
            sum = sum + term*term
         end do
         eval_rms_dvdt_div_v = sqrt(sum/(khi - klo + 1))
      end function eval_rms_dvdt_div_v
      
      
      real(dp) function cell_start_specific_KE(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: ke
         cell_start_specific_KE = cell_start_specific_KE_qp(s,k)
      end function cell_start_specific_KE
      
      
      real(qp) function cell_start_specific_KE_qp(s,k)
         ! for consistency with dual cells at faces, use <v**2> instead of <v>**2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(qp) :: qhalf, v0, v1, v2
         qhalf = 0.5d0
         if (s% u_flag) then
            v0 = s% u_start(k)
            cell_start_specific_KE_qp = qhalf*v0**2
         else if (s% v_flag) then
            v0 = s% v_start(k)
            if (k < s% nz) then
               v1 = s% v_start(k+1)
            else
               v1 = s% v_center
            end if
            v2 = qhalf*(v0**2 + v1**2)
            cell_start_specific_KE_qp = qhalf*v2
         else ! ignore kinetic energy if no velocity variables
            cell_start_specific_KE_qp = 0d0
         end if
      end function cell_start_specific_KE_qp
      
      
      real(dp) function cell_specific_KE(s,k,d_dv00,d_dvp1)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: d_dv00, d_dvp1
         real(dp) :: v2, dv2_dv00, dv2_dvp1
         cell_specific_KE = cell_specific_KE_qp(s,k,d_dv00,d_dvp1)
      end function cell_specific_KE
      
      
      real(qp) function cell_specific_KE_qp(s,k,d_dv00,d_dvp1)
         ! for consistency with dual cells at faces, use <v**2> instead of <v>**2
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: d_dv00, d_dvp1
         real(dp) :: dv2_dv00, dv2_dvp1
         real(qp) :: qhalf, v0, v1, v2
         qhalf = 0.5d0
         if (s% u_flag) then
            v0 = s% u(k)
            cell_specific_KE_qp = qhalf*v0**2
            d_dv00 = s% u(k)
            d_dvp1 = 0d0
         else if (s% v_flag) then
            v0 = s% v(k)
            if (k < s% nz) then
               v1 = s% v(k+1)
               dv2_dvp1 = s% v(k+1)
            else
               v1 = s% v_center
               dv2_dvp1 = 0d0
            end if
            v2 = qhalf*(v0**2 + v1**2)
            dv2_dv00 = s% v(k)
            cell_specific_KE_qp = qhalf*v2
            d_dv00 = 0.5d0*dv2_dv00
            d_dvp1 = 0.5d0*dv2_dvp1
         else ! ignore kinetic energy if no velocity variables
            cell_specific_KE_qp = 0d0
            d_dv00 = 0d0
            d_dvp1 = 0d0
         end if
      end function cell_specific_KE_qp
      
      
      real(dp) function cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: d_dlnR00,d_dlnRp1
         cell_specific_PE = cell_specific_PE_qp(s,k,d_dlnR00,d_dlnRp1)
      end function cell_specific_PE
      
      
      real(qp) function cell_specific_PE_qp(s,k,d_dlnR00,d_dlnRp1)
         ! for consistency with dual cells at faces, <m/r>_cntr => (m(k)/r(k) + m(k+1)/r(k+1))/2 /= m_cntr/r_cntr
         ! i.e., use avg of m/r at faces of cell rather than ratio of cell center mass over cell center r.
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: d_dlnR00,d_dlnRp1
         real(qp) :: qhalf, rp1, r00, mp1, m00, Gp1, G00, gravp1, grav00
         real(dp) :: d_grav00_dlnR00, d_gravp1_dlnRp1
         include 'formats'
         qhalf = 0.5d0
         if (k == s% nz) then
            rp1 = s% R_center
            mp1 = s% m_center
            Gp1 = s% cgrav(s% nz)
         else
            rp1 = s% r(k+1)
            mp1 = s% m(k+1)
            Gp1 = s% cgrav(k+1)
         end if
         if (rp1 <= 0d0) then
            gravp1 = 0d0
            d_gravp1_dlnRp1 = 0d0
         else
            gravp1 = -Gp1*mp1/rp1
            d_gravp1_dlnRp1 = -gravp1
         end if
         r00 = s% r(k)
         m00 = s% m(k)
         G00 = s% cgrav(k)
         grav00 = -G00*m00/r00
         d_grav00_dlnR00 = -grav00
         cell_specific_PE_qp = qhalf*(gravp1 + grav00)
         d_dlnR00 = 0.5d0*d_grav00_dlnR00
         d_dlnRp1 = 0.5d0*d_gravp1_dlnRp1
         if (is_bad(cell_specific_PE_qp)) then
            write(*,2) 'cell_specific_PE_qp', k, cell_specific_PE_qp
            write(*,2) 'gravp1', k, gravp1
            write(*,2) 'grav00', k, grav00
            stop 'cell_specific_PE'
         end if
      end function cell_specific_PE_qp
      
      
      real(dp) function cell_start_specific_PE(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         cell_start_specific_PE = cell_start_specific_PE_qp(s,k)
      end function cell_start_specific_PE
      
      
      real(dp) function cell_start_specific_PE_qp(s,k)
         ! for consistency with dual cells at faces, <m/r>_cntr => (m(k)/r(k) + m(k+1)/r(k+1))/2 /= m_cntr/r_cntr
         ! i.e., use avg of m/r at faces of cell rather than ratio of cell center mass over cell center r.
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(qp) :: qhalf, rp1, r00, mp1, m00, Gp1, G00, gravp1, grav00
         include 'formats'
         qhalf = 0.5d0
         if (k == s% nz) then
            rp1 = s% R_center
            mp1 = s% m_center
            Gp1 = s% cgrav(s% nz)
         else
            rp1 = s% r_start(k+1)
            mp1 = s% m(k+1)
            Gp1 = s% cgrav(k+1)
         end if
         if (rp1 <= 0d0) then
            gravp1 = 0d0
         else
            gravp1 = -Gp1*mp1/rp1
         end if
         r00 = s% r_start(k)
         m00 = s% m(k)
         G00 = s% cgrav(k)
         grav00 = -G00*m00/r00
         cell_start_specific_PE_qp = qhalf*(gravp1 + grav00)
         if (is_bad(cell_start_specific_PE_qp)) then
            write(*,2) 'cell_start_specific_PE_qp', k, cell_start_specific_PE_qp
            write(*,2) 'gravp1', k, gravp1
            write(*,2) 'grav00', k, grav00
            stop 'cell_start_specific_PE_qp'
         end if
      end function cell_start_specific_PE_qp
      
      
      real(dp) function cell_specific_rotational_energy(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: e_00, e_p1
         e_00 = s% i_rot(k)*s% omega(k)*s% omega(k)
         if (k < s% nz) then
            e_p1 = s% i_rot(k+1)*s% omega(k+1)*s% omega(k+1)
         else
            e_p1 = 0
         end if
         cell_specific_rotational_energy = 0.5d0*(e_p1 + e_00)
      end function cell_specific_rotational_energy
      
      
      real(dp) function eval_deltaM_total_from_profile( &
            deltaM, premesh_dm, profile)
         real(dp), intent(in) :: deltaM
         real(dp), intent(in) :: premesh_dm(:), profile(:)
         real(dp) :: sum_dm, dm, total, f
         integer :: k
         include 'formats'
         total = 0
         sum_dm = 0
         do k=1,size(premesh_dm,dim=1)
            if (sum_dm >= deltaM) exit
            dm = premesh_dm(k)
            if (sum_dm + dm > deltaM) then
               f = (deltaM - sum_dm)/dm
               total = total + f*profile(k)
               exit
            end if
            total = total + profile(k)
            sum_dm = sum_dm + dm
         end do
         eval_deltaM_total_from_profile = total
      end function eval_deltaM_total_from_profile
      
      
      real(dp) function cell_specific_total_energy(s, k) result(cell_total)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: d_dv00,d_dvp1,d_dlnR00,d_dlnRp1
         include 'formats'
         cell_total = s% energy(k)
         if (s% v_flag .or. s% u_flag) &
            cell_total = cell_total + cell_specific_KE(s,k,d_dv00,d_dvp1)
         if (.not. s% zero_gravity) &
            cell_total = cell_total + cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
         if (s% rotation_flag .and. s% include_rotation_in_total_energy) &
               cell_total = cell_total + cell_specific_rotational_energy(s,k)
         if (s% Eturb_flag) cell_total = cell_total + s% Eturb(k)
         if (s% rsp_flag) cell_total = cell_total + s% Et(k)
      end function cell_specific_total_energy
      
      
      subroutine eval_integrated_total_energy_profile(s, arr, direction, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: direction
         integer, intent(out) :: ierr
         real(dp), intent(out) :: arr(:)
         integer :: k, start, finish

         ierr = 0
         if (direction == 1) then
            start = 1
            finish = s%nz
         else if (direction == -1) then
            start = s%nz
            finish = 1
         else
            ierr = 1
            call mesa_error(__FILE__,__LINE__)
         end if

         arr(start) = cell_specific_total_energy(s,start) * s%dm(start)
         do k=start+direction, finish, direction
            arr(k) = arr(k-direction) + cell_specific_total_energy(s,k) * s%dm(k)
         end do
            
      end subroutine eval_integrated_total_energy_profile

      subroutine eval_deltaM_total_energy_integrals( &
            s, klo, khi, deltaM, save_profiles, &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
         type (star_info), pointer :: s
         integer, intent(in) :: klo, khi ! sum from klo to khi
         real(dp), intent(in) :: deltaM
         logical, intent(in) :: save_profiles
         real(dp), intent(out), dimension(:) :: total_energy_profile
         real(dp), intent(out) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total
         integer :: k
         real(dp) :: rR, rL, rC, dVR, dV, dm, m_cntr, rho, egas, v, sum_dm, &
            cell_total, cell1, d_dv00, d_dvp1, d_dlnR00, d_dlnRp1
         include 'formats'
         
         total_internal_energy = 0d0
         total_gravitational_energy = 0d0
         total_radial_kinetic_energy = 0d0
         total_rotational_kinetic_energy = 0d0
         total_turbulent_energy = 0d0
         sum_total = 0d0

         if (klo < 1 .or. khi > s% nz .or. klo > khi) return
         
         sum_dm = 0
         do k=klo,khi
            if (sum_dm >= deltaM) exit
            cell_total = 0
            dm = s% dm(k)
            if (sum_dm + dm > deltaM) dm = deltaM - sum_dm
            cell1 = dm*s% energy(k)
            cell_total = cell_total + cell1
            total_internal_energy = total_internal_energy + cell1
            if (s% v_flag .or. s% u_flag) then
               cell1 = dm*cell_specific_KE(s,k,d_dv00,d_dvp1)
               cell_total = cell_total + cell1
               total_radial_kinetic_energy = total_radial_kinetic_energy + cell1
            end if
            if (.not. s% zero_gravity) then
               cell1 = dm*cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
               cell_total = cell_total + cell1
               total_gravitational_energy = total_gravitational_energy + cell1
            end if
            if (s% rotation_flag) then
               cell1 = dm*cell_specific_rotational_energy(s,k)
               total_rotational_kinetic_energy = total_rotational_kinetic_energy + cell1
               if (s% include_rotation_in_total_energy) &
                  cell_total = cell_total + cell1
            end if
            if (s% Eturb_flag) then
               cell1 = dm*s% Eturb(k)
               cell_total = cell_total + cell1
               total_turbulent_energy = total_turbulent_energy + cell1
            end if
            if (s% rsp_flag) then
               cell1 = dm*s% Et(k)
               cell_total = cell_total + cell1
               total_turbulent_energy = total_turbulent_energy + cell1
            end if
            if (save_profiles) then
               total_energy_profile(k) = cell_total
            end if
         end do

         sum_total = total_internal_energy + total_gravitational_energy + &
            total_radial_kinetic_energy + total_turbulent_energy
            
         if (s% include_rotation_in_total_energy) &
            sum_total = sum_total + total_rotational_kinetic_energy

      end subroutine eval_deltaM_total_energy_integrals


      subroutine eval_total_energy_profile(s, total_energy_profile)
         type (star_info), pointer :: s
         real(dp), intent(out), dimension(:) :: total_energy_profile

         integer :: k
         real(dp) :: rR, rL, rC, dVR, dV, dm, m_cntr, rho, egas, v, &
            cell_total, cell1, d_dv00, d_dvp1, d_dlnR00, d_dlnRp1
         include 'formats'
         
         do k=1, s%nz
            cell_total = 0
            dm = s% dm(k)
            cell1 = dm*s% energy(k)
            cell_total = cell_total + cell1
            if (s% v_flag .or. s% u_flag) then
               cell1 = dm*cell_specific_KE(s,k,d_dv00,d_dvp1)
               cell_total = cell_total + cell1
            end if
            if (.not. s% zero_gravity) then
               cell1 = dm*cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
               cell_total = cell_total + cell1
            end if
            if (s% rotation_flag) then
               cell1 = dm*cell_specific_rotational_energy(s,k)
               if (s% include_rotation_in_total_energy) &
                  cell_total = cell_total + cell1
            end if
            if (s% Eturb_flag) then
               cell1 = dm*s% Eturb(k)
               cell_total = cell_total + cell1
            end if
            if (s% rsp_flag) then
               cell1 = dm*s% Et(k)
               cell_total = cell_total + cell1
            end if
            total_energy_profile(k) = cell_total
         end do
            
      end subroutine eval_total_energy_profile
      
      
      real(dp) function eval_cell_section_total_energy( &
            s, klo, khi) result(sum_total)
         type (star_info), pointer :: s
         integer, intent(in) :: klo, khi ! sum from klo to khi
         real(dp) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy
         real(dp), allocatable, dimension(:) :: total_energy_profile
         allocate(total_energy_profile(1:s% nz))
         call eval_deltaM_total_energy_integrals( &
            s, klo, khi, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end function eval_cell_section_total_energy


      subroutine eval_total_energy_integrals(s, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
         type (star_info), pointer :: s
         real(dp), intent(out) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total
         real(dp), allocatable, dimension(:) :: total_energy_profile
         allocate(total_energy_profile(1:s% nz))
         call eval_deltaM_total_energy_integrals( &
            s, 1, s% nz, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end subroutine eval_total_energy_integrals


      subroutine eval_total_energy_integrals_and_profile(s, &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
         type (star_info), pointer :: s
         real(dp), intent(out) :: total_energy_profile(:), &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total
         call eval_deltaM_total_energy_integrals( &
            s, 1, s% nz, s% mstar, .true., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end subroutine eval_total_energy_integrals_and_profile


      real(dp) function get_total_energy_integral(s,k) result(sum_total)
         ! from surface down to k
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy
         real(dp), allocatable, dimension(:) :: total_energy_profile
         allocate(total_energy_profile(1:s% nz))
         call eval_deltaM_total_energy_integrals( &
            s, 1, k, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end function get_total_energy_integral


      real(dp) function get_total_energy_integral_outward(s,k) result(sum_total)
         ! from surface down to k
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy
         real(dp), allocatable, dimension(:) :: total_energy_profile
         allocate(total_energy_profile(1:s% nz))
         call eval_deltaM_total_energy_integrals( &
            s, k, s% nz, s% mstar, .false., &
            total_energy_profile, &
            total_internal_energy, total_gravitational_energy, &
            total_radial_kinetic_energy, total_rotational_kinetic_energy, &
            total_turbulent_energy, sum_total)
      end function get_total_energy_integral_outward


      real(dp) function get_log_concentration(s,j,k) result(log_c)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: j, k
         ! concentration = number density / number density of electrons
         !  Ci = (Xi/Ai) / sum(Zi*Xi/Ai)   [see Thoul et al, ApJ 421:828-842, 1994]
         integer :: i, cid, species
         real(dp) :: tmp, c
         log_c = -1d99
         if (s% chem_id(j) == 0) return
         species = s% species
         tmp = 0d0
         do i=1,species
            cid = s% chem_id(i)
            tmp = tmp + chem_isos% Z(cid)*s% xa(i,k)/chem_isos% Z_plus_N(cid)
         end do
         cid = s% chem_id(j)
         c = (s% xa(j,k)/chem_isos% Z_plus_N(cid))/tmp
         log_c = safe_log10(c)
      end function get_log_concentration


      real(dp) function get_phi_Joss(s,k) result(phi)
         use eos_def, only: i_lnPgas
         ! Joss, Salpeter, Ostriker, 1973. density inversion when Lrad/Ledd > phi.
         type (star_info), pointer :: s
         integer, intent(in) :: k
         phi = 1d0/(1d0 + (s% Pgas(k)/(4* s% Prad(k)))*s% d_eos_dlnT(i_lnPgas,k))
      end function get_phi_Joss


      logical function after_He_burn(s, he4_limit)
         use chem_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: he4_limit
         integer :: nz, h1, he4
         real(dp) :: small = 1d-4
         after_He_burn = .false.
         nz = s% nz
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         if (h1 == 0 .or. he4 == 0) return
         if (s% xa(h1,nz) > small .or. s% xa(he4,nz) > he4_limit) return
         after_He_burn = .true.
      end function after_He_burn


      logical function after_C_burn(s, c12_limit)
         use chem_def
         type (star_info), pointer :: s
         real(dp), intent(in) :: c12_limit
         integer :: nz, h1, he4, c12
         real(dp) :: small = 1d-4
         after_C_burn = .false.
         nz = s% nz
         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         c12 = s% net_iso(ic12)
         if (h1 == 0 .or. he4 == 0 .or. c12 == 0) return
         if (s% xa(h1,nz) > small .or. s% xa(he4,nz) > small .or. &
             s% xa(c12,nz) > c12_limit) return
         after_C_burn = .true.
      end function after_C_burn
      
      
      integer function lookup_nameofvar(s, namestr)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: namestr
         integer :: i
         lookup_nameofvar = 0
         do i=1,s% nvar
            if (namestr == s% nameofvar(i)) then
               lookup_nameofvar = i
               return
            end if
         end do
      end function lookup_nameofvar
      
      
      integer function lookup_nameofequ(s, namestr)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: namestr
         integer :: i
         lookup_nameofequ = 0
         do i=1,s% nvar
            if (namestr == s% nameofequ(i)) then
               lookup_nameofequ = i
               return
            end if
         end do
      end function lookup_nameofequ


      real(dp) function omega_crit(s, k)  ! result always is > 0
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: Ledd, gamma_factor, Lrad_div_Ledd, rmid, r003, rp13
         include 'formats'
         if (s% fitted_fp_ft_i_rot .and. s% rotation_flag) then
            ! Use equatorial radius (at center of cell by volume)
            r003 = s% r_equatorial(k)*s% r_equatorial(k)*s% r_equatorial(k)
            if (k < s% nz) then
               rp13 = s% r_equatorial(k+1)*s% r_equatorial(k+1)*s% r_equatorial(k+1)
            else
               rp13 = s% R_center*s% R_center*s% R_center
            end if
            rmid = pow(0.5d0*(r003 + rp13),1d0/3d0)
         else
           rmid = s% rmid(k)
         end if

         Ledd = 4*pi*clight*s% cgrav(k)*s% m_grav(k)/s% opacity(k)
         Lrad_div_Ledd = get_Lrad_div_Ledd(s,k)
         gamma_factor = 1d0 - min(Lrad_div_Ledd, 0.9999d0)
         omega_crit = sqrt(gamma_factor*s% cgrav(k)*s% m_grav(k)/pow3(rmid))
      end function omega_crit


      subroutine set_surf_avg_rotation_info(s)
         use hydro_rotation, only: set_rotation_info
         type (star_info), pointer :: s
         real(dp) :: &
            dm, dmsum, omega_sum, omega_crit_sum, omega_div_omega_crit_sum, &
            v_rot_sum, v_crit_sum, v_div_v_crit_sum, Lrad_div_Ledd_sum, &
            kap_face, Ledd, gamma_factor, omega_crit, omega, kap_sum, &
            j_rot_sum, j_rot, v_rot, v_crit, Lrad_div_Ledd, dtau, tau, &
            cgrav, kap, mmid, Lmid, rmid, logT_sum, logRho_sum
         integer :: k
         logical, parameter :: dbg = .false.
         include 'formats'

         if (.not. s% rotation_flag) then
            s% omega_avg_surf = 0
            s% omega_crit_avg_surf = 0
            s% w_div_w_crit_avg_surf = 0
            s% j_rot_avg_surf = 0
            s% v_rot_avg_surf = 0
            s% v_crit_avg_surf = 0
            s% v_div_v_crit_avg_surf = 0
            s% Lrad_div_Ledd_avg_surf = 0
            s% opacity_avg_surf = 0
            s% logT_avg_surf = 0
            s% logRho_avg_surf = 0
            return
         end if

         call set_rotation_info(s,.true.,k)

         tau = s% tau_factor*s% tau_base
         dmsum = 0d0
         Lrad_div_Ledd_sum = 0d0
         rmid = 0d0

         do k = 1, s% nz - 1
            kap = s% opacity(k)
            rmid = s% rmid(k)
            mmid = 0.5d0*(s% m_grav(k) + s% m_grav(k+1))
            Lmid = 0.5d0*(s% L(k) + s% L(k+1))
            cgrav = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
            dm = s% dm(k)
            dtau = dm*kap/(4*pi*rmid*rmid)

            if (tau + dtau <= s% surf_avg_tau_min) then
               tau = tau + dtau
               cycle
            end if

            ! check for partial contribution from cell
            ! the tau < s% surf_avg_tau is meant for the case in which the surface tau is set
            ! equal or larger to surf_avg_tau. In that case we just use the values of the surface cell.
            if (tau < s% surf_avg_tau) then
               if (tau < s% surf_avg_tau_min) then ! only use part of this cell
                  dm = dm*(tau + dtau - s% surf_avg_tau_min)/dtau
               else if (tau + dtau > s% surf_avg_tau) then ! only use part of this cell
                  dm = dm*(s% surf_avg_tau - tau)/dtau
                  !write(*,2) 'tau limit', k, (s% surf_avg_tau - tau)/dtau
               end if
            end if
            dmsum = dmsum + dm
            Lrad_div_Ledd = get_Lrad_div_Ledd(s,k)
            Lrad_div_Ledd_sum = Lrad_div_Ledd_sum + dm*Lrad_div_Ledd
            tau = tau + dtau
            if (tau >= s% surf_avg_tau) exit
         end do

         s% Lrad_div_Ledd_avg_surf = Lrad_div_Ledd_sum/dmsum
         if (s% generations > 2) & ! time average
            s% Lrad_div_Ledd_avg_surf = &
               0.5d0*(s% Lrad_div_Ledd_avg_surf + s% Lrad_div_Ledd_avg_surf_old)

         gamma_factor = 1d0 - min(s% Lrad_div_Ledd_avg_surf, 0.9999d0)

         tau = s% tau_factor*s% tau_base
         dmsum = 0
         j_rot_sum = 0
         omega_sum = 0
         omega_crit_sum = 0
         omega_div_omega_crit_sum = 0
         v_rot_sum = 0
         v_crit_sum = 0
         v_div_v_crit_sum = 0
         kap_sum = 0
         logT_sum = 0
         logRho_sum = 0

         do k = 1, s% nz - 1

            kap = s% opacity(k)
            if (s% fitted_fp_ft_i_rot) then
               ! TODO: better explain
               ! Use equatorial radius
               rmid = 0.5d0*(s% r_equatorial(k) + s% r_equatorial(k+1))
            else
              rmid = s% rmid(k)
            end if
            dm = s% dm(k)
            dtau = dm*kap/(4*pi*rmid*rmid)

            if (tau + dtau <= s% surf_avg_tau_min) then
               tau = tau + dtau
               cycle
            end if

            ! check for partial contribution from cell
            ! the tau < s% surf_avg_tau is meant for the case in which the surface tau is set
            ! equal or larger to surf_avg_tau. In this case we just use the values of the surface cell.
            if (tau < s% surf_avg_tau) then
               if (tau < s% surf_avg_tau_min) then ! only use part of this cell
                  dm = dm*(tau + dtau - s% surf_avg_tau_min)/dtau
               else if (tau + dtau > s% surf_avg_tau) then ! only use part of this cell
                  dm = dm*(s% surf_avg_tau - tau)/dtau
               end if
            end if

            dmsum = dmsum + dm
            cgrav = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
            mmid = 0.5d0*(s% m_grav(k) + s% m_grav(k+1))
            omega = 0.5d0*(s% omega(k) + s% omega(k+1))
            j_rot = 0.5d0*(s% j_rot(k) + s% j_rot(k+1))

            kap_sum = kap_sum + dm*kap
            j_rot_sum = j_rot_sum + dm*j_rot

            omega_crit = sqrt(gamma_factor*cgrav*mmid/pow3(rmid))
            omega_div_omega_crit_sum = omega_div_omega_crit_sum + dm*abs(omega/omega_crit)

            v_rot = omega*rmid
            v_crit = omega_crit*rmid
            omega_sum = omega_sum + dm*omega
            omega_crit_sum = omega_crit_sum + dm*omega_crit
            v_rot_sum = v_rot_sum + dm*v_rot
            v_crit_sum = v_crit_sum + dm*v_crit
            v_div_v_crit_sum = v_div_v_crit_sum + dm*abs(v_rot/v_crit)
            logT_sum = logT_sum + dm*s% lnT(k)/ln10
            logRho_sum = logRho_sum + dm*s% lnd(k)/ln10
            kap_sum = kap_sum + dm*kap
            tau = tau + dtau
            
            
            if (tau >= s% surf_avg_tau) exit

         end do

         s% logT_avg_surf = logT_sum/dmsum
         s% logRho_avg_surf = logRho_sum/dmsum
         s% opacity_avg_surf = kap_sum/dmsum
         s% j_rot_avg_surf = j_rot_sum/dmsum
         s% omega_avg_surf = omega_sum/dmsum
         s% omega_crit_avg_surf = omega_crit_sum/dmsum
         s% w_div_w_crit_avg_surf = omega_div_omega_crit_sum/dmsum
         s% v_rot_avg_surf = v_rot_sum/dmsum
         s% v_crit_avg_surf = v_crit_sum/dmsum
         s% v_div_v_crit_avg_surf = v_div_v_crit_sum/dmsum

      end subroutine set_surf_avg_rotation_info


      subroutine median_smoothing(dd, n, ns, dmed)
         use num_lib, only: qsort
         real(dp), intent(inout) :: dd(:) ! (n)
         integer, intent(in) :: n, ns
         real(dp), intent(inout) :: dmed(:) ! (n) work array

         real(dp) :: x(2*ns+1)
         integer :: i, j, k, nmed, index(2*ns+1)

         nmed = 2*ns+1

         do i=1,n
            if ((i > 1+ns) .and. (i < n-ns)) then
               k = 1
               do j = i-ns, i+ns
                  x(k) = dd(j)
                  k = k+1
               end do
               call qsort(index,nmed,x)
               dmed(i) = x(index(ns+1))
            else
               dmed(i) = dd(i)
            end if
         end do

         do i=1,n
            if (dmed(i) /= 0) dd(i) = dmed(i)
         end do

      end subroutine median_smoothing


      subroutine weighed_smoothing(dd, n, ns, preserve_sign, ddold)
      !     based on routine written by S.-C. Yoon, 18 Sept. 2002
      !     for smoothing  any variable (dd) with size n over 2*ns+1 cells.
         real(dp), intent(inout) :: dd(:) ! (n)
         integer, intent(in) :: n, ns
         logical, intent(in) :: preserve_sign
         real(dp), intent(inout) :: ddold(:) ! (n) work array

         integer :: nweight, mweight, i, j, k
         real(dp) :: weight(2*ns+1), sweight, v0

         include 'formats'

         do i = 1,n
           ddold(i) = dd(i)
         end do

         !--preparation for smoothing --------
         nweight = ns
         mweight = 2*nweight+1
         do i = 1,mweight
            weight(i) = 0d0
         end do
         weight(1) = 1d0
         do i = 1,mweight-1
            do j = i+1,2,-1
               weight(j) = weight(j) + weight(j-1)
            end do
         end do

         !--smoothing ------------------------
         do i=2,n-1
            sweight=0d0
            dd(i)=0d0
            v0 = ddold(i)
            do j = i, max(1,i-nweight), -1
               k=j-i+nweight+1
               if (preserve_sign .and. v0*ddold(j) <= 0) exit
               sweight = sweight+weight(k)
               dd(i) = dd(i)+ddold(j)*weight(k)
            end do
            do j = i+1, min(n,i+nweight)
               k=j-i+nweight+1
               if (preserve_sign .and. v0*ddold(j) <= 0) exit
               sweight = sweight+weight(k)
               dd(i) = dd(i)+ddold(j)*weight(k)
            end do
            if (sweight > 0) then
               sweight = 1d0/sweight
               dd(i) = dd(i)*sweight
            end if
         end do

      end subroutine weighed_smoothing

      
      subroutine threshold_smoothing (dd, dd_thresh, n, ns, preserve_sign, ddold)

        ! Same as weighed_smoothing, but only smooth contiguous regions where |dd| >= dd_thresh

        real(dp), intent(inout) :: dd(:)    ! (n)
        real(dp), intent(in)    :: dd_thresh
        integer, intent(in)     :: n
        integer, intent(in)     :: ns
        logical, intent(in)     :: preserve_sign
        real(dp), intent(inout) :: ddold(:) ! (n) work array

        logical :: in_region
        integer :: i
        integer :: i_a
        integer :: i_b
        
        include 'formats'

        ! Process regions

        in_region = .FALSE.

        i_a = 1
        do i = 1, n

           if (in_region) then

              if (ABS(dd(i)) < dd_thresh) then
                 i_b = i-1
                 if (i_b > i_a) call weighed_smoothing(dd(i_a:i_b), i_b-i_a+1, ns, preserve_sign, ddold(i_a:i_b))
                 in_region = .FALSE.
              endif

           else
              if (ABS(dd(i)) >= dd_thresh) then
                 i_a = i
                 in_region = .TRUE.
              endif

           end if

        end do

        ! Handle the final region

        if (in_region) then

           i_b = n
           if (i_b > i_a) call weighed_smoothing(dd(i_a:i_b), i_b-i_a+1, ns, preserve_sign, ddold(i_a:i_b))

        endif

        ! Finish

        return

      end subroutine threshold_smoothing
      

      real(dp) function eval_kh_timescale(G,M,R,L) result(kh)
         real(dp), intent(in) :: G,M,R,L
         if (L <= 0) then
            kh = 0d0
         else
            kh = 0.75d0*G*M*M/(R*L) ! 0.75 is based on sun.  Hansen & Kawaler eqn 1.30
         end if
      end function eval_kh_timescale



      real(dp) function yrs_for_init_timestep(s)
         type (star_info), pointer :: s
         if (s% initial_mass <= 1) then
            yrs_for_init_timestep = 1d5
         else
            yrs_for_init_timestep = 1d5 / pow(s% initial_mass,2.5d0)
         end if
      end function yrs_for_init_timestep


      subroutine set_phase_of_evolution(s) ! from evolve after call do_report
         use rates_def, only: i_rate
         use chem_def
         type (star_info), pointer :: s
         real(dp) :: power_he_burn, power_c_burn, power_neutrinos, &
            center_h1, center_he4
         integer :: nz, j
         include 'formats'
         nz = s% nz
         
         j = s% net_iso(ih1)
         if (j > 0) then
            center_h1 = center_avg_x(s,j)
         else
            center_h1 = 1d99
         end if
         j = s% net_iso(ihe4)
         if (j > 0) then
            center_he4 = center_avg_x(s,j)
         else
            center_he4 = 1d99
         end if

         if (s% center_gamma > s% gamma_center_limit) then
            s% phase_of_evolution = phase_WDCS
         else if (s% L_by_category(i_burn_si) > 1d2) then
            s% phase_of_evolution = phase_Si_Burn
         else if (s% L_by_category(i_burn_o) > 1d2 .and. center_he4 < 1d-4) then
            s% phase_of_evolution = phase_O_Burn
         else if (s% L_by_category(i_burn_ne) > 1d2 .and. center_he4 < 1d-4) then
            s% phase_of_evolution = phase_Ne_Burn
         else if (s% L_by_category(i_burn_c) > 1d2 .and. center_he4 < 1d-4) then
            s% phase_of_evolution = phase_C_Burn
         else if (center_he4 < 1d-4 .and. &
            s% he_core_mass - s% c_core_mass <= 0.1d0 .and. &
            any(s% burn_he_conv_region(1:s% num_conv_boundaries))) then
            s% phase_of_evolution = phase_TP_AGB
         else if (center_he4 <= 1d-4) then
            s% phase_of_evolution = phase_TACHeB          
         else if (s% center_eps_burn(i3alf) > Lsun) then
            s% phase_of_evolution = phase_ZACHeB
         else if (s% L_by_category(i3alf) > 1d2) then
            s% phase_of_evolution = phase_He_Burn
         else if (center_h1 <= 1d-6) then
            s% phase_of_evolution = phase_TAMS            
         else if (center_h1 <= 0.3d0) then
            s% phase_of_evolution = phase_IAMS            
         else if (s% L_nuc_burn_total >= s% L_phot*s% Lnuc_div_L_zams_limit) then
            s% phase_of_evolution = phase_ZAMS
         else if (s% log_center_temperature > 5d0) then
            s% phase_of_evolution = phase_PreMS
         else
            s% phase_of_evolution = phase_starting
         end if 
         
         
         
         
         
         
!         if (.not. arrived_main_seq(s) .or. s% phase_of_evolution == phase_carbon_burning) return
!         power_he_burn = s% power_he_burn
!         power_c_burn = dot_product(s% dm(1:nz), s% eps_nuc_categories(i_burn_c,1:nz))/Lsun
!         power_neutrinos = s% power_neutrinos
!         if (s% phase_of_evolution == phase_helium_burning .and. power_c_burn > power_neutrinos) then
!            s% phase_of_evolution = phase_carbon_burning
!         else if (power_c_burn + power_he_burn > power_neutrinos) then
!            s% phase_of_evolution = phase_helium_burning
!         else if (s% center_he4 < center_he_going) then
!            s% phase_of_evolution = phase_helium_burning
!         else if (s% center_h1 < center_h_gone) then
!            s% phase_of_evolution = phase_wait_for_he
!         else if (s% center_h1 < center_h_going) then
!            s% phase_of_evolution = phase_mid_main_seq
!         else
!            s% phase_of_evolution = phase_early_main_seq
!         end if
      end subroutine set_phase_of_evolution



      logical function arrived_main_seq(s)
         type (star_info), pointer :: s
         include 'formats'
         arrived_main_seq = &
            (s% L_nuc_burn_total >= s% L_phot) .and. &
            (s% power_h_burn >= s% L_nuc_burn_total/2)
         return
         write(*,1) 's% L_nuc_burn_total', s% L_nuc_burn_total
         write(*,1) 's% L_phot', s% L_phot
         write(*,1) 's% power_h_burn', s% L_phot
         write(*,*) 'arrived_main_seq',  arrived_main_seq
         write(*,*)
      end function arrived_main_seq


      subroutine save_for_d_dt(s)
         ! these values will be modified as necessary by adjust mass
         type (star_info), pointer :: s
         integer :: k, nz, i_lnR, i_lnT, i_lnd, i_eturb, &
            i_v, i_u, i_alpha_RTI, i_ln_cvpv0
         include 'formats'
         nz = s% nz
         i_lnR = s% i_lnR
         i_lnT = s% i_lnT
         i_lnd = s% i_lnd
         i_eturb = s% i_eturb
         i_v = s% i_v
         i_u = s% i_u
         i_alpha_RTI = s% i_alpha_RTI
         i_ln_cvpv0 = s% i_ln_cvpv0
         do k=1, nz
            s% lnR_for_d_dt_const_m(k) = s% xh(i_lnR, k)
         end do
         if (i_lnT /= 0) then
            do k=1, nz
               s% lnT_for_d_dt_const_m(k) = s% xh(i_lnT, k)
               s% lnT_for_d_dt_const_q(k) = s% xh(i_lnT, k)
            end do
         end if
         if (i_lnd /= 0) then
            do k=1, nz
               s% lnd_for_d_dt_const_m(k) = s% xh(i_lnd, k)
               s% lnd_for_d_dt_const_q(k) = s% xh(i_lnd, k)
            end do
         end if
         if (i_ln_cvpv0 /= 0) then
            do k=1, nz
               s% ln_cvpv0_for_d_dt_const_m(k) = s% xh(i_ln_cvpv0, k)
               s% ln_cvpv0_for_d_dt_const_q(k) = s% xh(i_ln_cvpv0, k)
            end do
         end if
         if (i_v /= 0) then
            do k=1, nz
               s% v_for_d_dt_const_m(k) = s% xh(i_v, k)
            end do
         end if
         if (i_eturb /= 0) then
            do k=1, nz
               s% eturb_for_d_dt_const_m(k) = s% xh(i_eturb, k)
            end do
         end if
         if (i_u /= 0) then
            do k=1, nz
               s% u_for_d_dt_const_m(k) = s% xh(i_u, k)
            end do
         end if
         if (i_alpha_RTI /= 0) then
            do k=1, nz
               s% alpha_RTI_for_d_dt_const_m(k) = s% xh(i_alpha_RTI, k)
            end do
         end if
      end subroutine save_for_d_dt
      
      
      subroutine set_rv_info(s,k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: r2
         include 'formats'
         r2 = s% r(k)*s% r(k)
         if (s% using_Fraley_time_centering) then
            s% R2(k) = &
               (r2 + s% r_start(k)*s% r(k) + s% r_start(k)*s% r_start(k))/3d0
            s% d_R2_dlnR(k) = (2d0*r2 + s% r_start(k)*s% r(k))/3d0            
         else
            s% R2(k) = r2
            s% d_R2_dlnR(k) = 2d0*r2
         end if
         if (s% v_flag) then
            if (s% using_Fraley_time_centering) then
               s% vc(k) = 0.5d0*(s% v_start(k) + s% v(k))
            else
               s% vc(k) = s% v(k)
            end if
         else
            s% vc(k) = 0d0
         end if
      end subroutine set_rv_info
      
      
      subroutine show_matrix(s, dmat, nvar)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         real(dp) :: dmat(nvar,nvar)
         integer :: i, j
         write(*,*)
         write(*,'(18x)', advance = 'no') 
         do j = 1, nvar
            write(*,'(a15)', advance = 'no') s% nameofvar(j)
         end do
         write(*,*)
         do i = 1, nvar
            write(*,'(a15)', advance = 'no') s% nameofequ(i)
            do j = 1, nvar
               write(*,'(e13.4,2x)', advance = 'no') dmat(i,j)
            end do
            write(*,*)
         end do
         write(*,*)
      end subroutine show_matrix


      ! e00(i,j,k) is partial of equ(i,k) wrt var(j,k)
      subroutine e00(s,xscale,i,j,k,nvar,v)
         use num_def, only: &
            block_tridiag_dble_matrix_type, block_tridiag_quad_matrix_type
         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: i, j, k, nvar
         real(dp), intent(in) :: v
         integer :: b, q, v00
         real(qp) :: q1, q2
         !logical, parameter :: dbg = .true.
         logical, parameter :: dbg = .false.
         include 'formats'
         
         if (mdb .and. k==397 .and. v /= 0d0) &
            write(*,4) 'e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         
         !if (j == s% i_lnd .and. k /= s% nz) return ! this variable is being held constant
         
         if (v == 0d0) return
         
         if (.false. .and. j == s% i_lnT .and. k == 30) then
            write(*,4) 'e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v, xscale(j,k)
         end if
         
         if (is_bad(v)) then
!$omp critical (star_utils_e00_crit1)
            write(*,4) 'e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
            if (s% stop_for_bad_nums) stop 'e00'
!$omp end critical (star_utils_e00_crit1)
         end if
         
         if (i <= 0 .or. j <= 0 .or. k <= 0 .or. k > s% nz) then
            write(*,4) 'bad i,j,k e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
            stop 'e00'
         end if
         
         if (j > nvar) return ! hybrid
         
         if (i > nvar) then
!$omp critical (star_utils_e00_crit2)
            write(*,5) 'bad i e00(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), &
               s% solver_iter, i, j, k, v
            stop 'e00'
!$omp end critical (star_utils_e00_crit2)
         end if

         if (abs(v) < 1d-250) return

         if (associated(xscale)) then
            s% dblk(i,j,k) = s% dblk(i,j,k) + v*xscale(j,k)
         else
            s% dblk(i,j,k) = s% dblk(i,j,k) + v
         end if

      end subroutine e00


      ! em1(i,j,k) is partial of equ(i,k) wrt var(j,k-1)
      subroutine em1(s,xscale,i,j,k,nvar,v)
         use num_def, only: &
            block_tridiag_dble_matrix_type, block_tridiag_quad_matrix_type
         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: i, j, k, nvar
         real(dp), intent(in) :: v
         integer :: b, q, vm1
         real(qp) :: q1, q2
         !logical, parameter :: dbg = .true.
         logical, parameter :: dbg = .false.
         if (k == 1) return
         include 'formats'
         
         if (mdb .and. k==397 .and. v /= 0d0) &
            write(*,4) 'em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         
         if (v == 0d0) return
         
         if (.false. .and. j == s% i_lnT .and. k == 31) then
            write(*,4) 'em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v, xscale(j,k-1)
         end if
         
         if (s% Eturb_flag .and. j == s% i_lum) then ! assume j = 0 means partial wrt L
            write(*,2) 'cannot have Eturb_flag and partials wrt L(k-1)', k
            stop 'em1'
         end if
         
         if (is_bad(v)) then
!$omp critical (star_utils_em1_crit1)
            write(*,4) 'em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
            if (s% stop_for_bad_nums) stop 'em1'
!$omp end critical (star_utils_em1_crit1)
         end if
         
         if (i <= 0 .or. j <= 0 .or. k <= 0 .or. k > s% nz) then
            write(*,4) 'bad i,j,k em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
            stop 'em1'
         end if
         
         if (j > nvar) return ! hybrid
         
         if (i > nvar) then
            write(*,5) 'bad i em1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), &
               s% solver_iter, i, j, k, v
            stop 'em1'
         end if

         if (abs(v) < 1d-250) return

         if (associated(xscale)) then
            s% lblk(i,j,k) = s% lblk(i,j,k) + v*xscale(j,k-1)
         else
            s% lblk(i,j,k) = s% lblk(i,j,k) + v
         end if

      end subroutine em1


      ! ep1(i,j,k) is partial of equ(i,k) wrt var(j,k+1)
      subroutine ep1(s,xscale,i,j,k,nvar,v)
         use num_def, only: &
            block_tridiag_dble_matrix_type, block_tridiag_quad_matrix_type
         type (star_info), pointer :: s
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: i, j, k, nvar
         real(dp), intent(in) :: v
         integer :: b, q, vp1
         real(qp) :: q1, q2
         !logical, parameter :: dbg = .true.
         logical, parameter :: dbg = .false.
         include 'formats'
         
         if (mdb .and. k==397 .and. v /= 0d0) &
            write(*,4) 'ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
         
         if (v == 0d0) return
         
         
         if (.false. .and. j == s% i_lnT .and. k == 29) then
            write(*,4) 'ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v, xscale(j,k+1)
         end if
         
         if (is_bad(v)) then
!$omp critical (star_utils_ep1_crit1)
            write(*,4) 'ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
            if (s% stop_for_bad_nums) stop 'ep1'
!$omp end critical (star_utils_ep1_crit1)
         end if
         
         if (i <= 0 .or. j <= 0 .or. k <= 0 .or. k > s% nz) then
            write(*,4) 'bad i,j,k ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), i, j, k, v
            stop 'ep1'
         end if
         
         if (j > nvar) return ! hybrid
         
         if (i > nvar) then
            write(*,5) 'bad i ep1(i,j,k) ' // &
               trim(s% nameofequ(i)) // ' ' // trim(s% nameofvar(j)), &
               s% solver_iter, i, j, k, v
            stop 'ep1'
         end if

         if (abs(v) < 1d-250) return

         if (associated(xscale)) then
            s% ublk(i,j,k) = s% ublk(i,j,k) + v*xscale(j,k+1)
         else
            s% ublk(i,j,k) = s% ublk(i,j,k) + v
         end if

      end subroutine ep1


      real(dp) function current_min_xa_hard_limit(s) result(min_xa_hard_limit)
         type (star_info), pointer :: s
         real(dp) :: logTc, alfa
         logTc = s% lnT(s% nz)/ln10
         if (logTc <= s% logT_max_for_min_xa_hard_limit) then
            min_xa_hard_limit = s% min_xa_hard_limit
         else if (logTc >= s% logT_min_for_min_xa_hard_limit_for_highT) then
            min_xa_hard_limit = s% min_xa_hard_limit_for_highT
         else
            alfa = (logTc - s% logT_max_for_min_xa_hard_limit) / &
                   (s% logT_min_for_min_xa_hard_limit_for_highT - s% logT_max_for_min_xa_hard_limit)
            min_xa_hard_limit = &
               alfa*s% min_xa_hard_limit_for_highT + (1d0 - alfa)*s% min_xa_hard_limit
         end if
      end function current_min_xa_hard_limit


      real(dp) function current_sum_xa_hard_limit(s) result(sum_xa_hard_limit)
         type (star_info), pointer :: s
         real(dp) :: logTc, alfa
         include 'formats'
         logTc = s% lnT(s% nz)/ln10
         if (logTc <= s% logT_max_for_sum_xa_hard_limit) then
            sum_xa_hard_limit = s% sum_xa_hard_limit
         else if (logTc >= s% logT_min_for_sum_xa_hard_limit_for_highT) then
            sum_xa_hard_limit = s% sum_xa_hard_limit_for_highT
         else
            alfa = (logTc - s% logT_max_for_sum_xa_hard_limit) / &
                   (s% logT_min_for_sum_xa_hard_limit_for_highT - s% logT_max_for_sum_xa_hard_limit)
            sum_xa_hard_limit = &
               alfa*s% sum_xa_hard_limit_for_highT + (1d0 - alfa)*s% sum_xa_hard_limit
         end if
         !write(*,1) 'logTc hard_limit', logTc, sum_xa_hard_limit
      end function current_sum_xa_hard_limit


      integer function no_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         no_extra_profile_columns = 0
      end function no_extra_profile_columns


      subroutine no_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: maxlen_profile_column_name, star_info
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine no_data_for_extra_profile_columns


      ! Stiriba, Youssef, Appl, Numer. Math. 45, 499-511. 2003.

         ! LPP-HARMOD -- local piecewise parabolic reconstruction

         ! interpolant is derived to conserve integral of v in cell k
         ! interpolant slope at cell midpoint is harmonic mean of slopes between adjacent cells
         ! where these slopes between cells are defined as the difference between cell averages
         ! divided by the distance between cell midpoints.
         ! interpolant curvature based on difference between the midpoint slope
         ! and the smaller in magnitude of the slopes between adjacent cells.

         ! interpolant f(dq) = a + b*dq + (c/2)*dq^2, with dq = q - q_midpoint
         ! c0 holds a's, c1 holds b's, and c2 holds c's.



      subroutine get1_lpp(k, nz, &
            dm1, d00, dp1, vm1, v00, vp1, quadratic, monotonic, dbg, c0, c1, c2)
         integer, intent(in) :: k, nz
         real(dp), intent(in) :: dm1, d00, dp1, vm1, v00, vp1
         logical, intent(in) :: quadratic, monotonic
         logical, intent(in) :: dbg
         real(dp), intent(out) :: c0, c1, c2
         real(dp) :: vbdy1, vbdy2, dhalf, sm1, s00, sp1, sprod, dv1, dv2
         logical :: okay
         include 'formats'

         c0 = v00
         c2 = 0d0
         if (k == 1) then
            sm1 = 0d0
            sp1 = (v00 - vp1) / (0.5d0*(d00 + dp1))
            c1 = sp1
         else if (k == nz) then
            sp1 = 0d0
            sm1 = (vm1 - v00) / (0.5d0*(dm1 + d00))
            c1 = sm1
            if (dbg) then
               write(*,2) 'get1_lpp c1', k, c1
            end if
         else
            sm1 = (vm1 - v00) / (0.5d0*(dm1 + d00))
            sp1 = (v00 - vp1) / (0.5d0*(d00 + dp1))
            sprod = sm1*sp1
            if (sprod <= 0d0) then
               ! at local min or max, so set slope and curvature to 0.
               c1 = 0d0
               return
            end if
            if (.not. quadratic) then ! minmod
               s00 = (vm1 - vp1)/(d00 + (dm1 + dp1)/2)
               c1 = sm1
               if (sm1*s00 < 0d0) c1 = 0d0
               if (abs(s00) < abs(c1)) c1 = s00
               if (s00*sp1 < 0d0) c1 = 0d0
               if (abs(sp1) < abs(c1)) c1 = sp1
            else
               c1 = sprod*2/(sp1 + sm1) ! harmonic mean slope
               if (abs(sm1) <= abs(sp1)) then
                  c2 = (sm1 - c1)/(2d0*d00)
               else
                  c2 = (c1 - sp1)/(2d0*d00)
               end if
               c0 = v00 - c2*d00*d00/24d0
            end if
         end if

         if (.not. monotonic) return

         dhalf = 0.5d0*d00
         dv1 = c1*dhalf
         dv2 = 0.5d0*c2*dhalf*dhalf
         vbdy1 = c0 + dv1 + dv2 ! value at face(k)
         vbdy2 = c0 - dv1 + dv2 ! value at face(k+1)

         okay = .true.

         if (k > 1) then ! check v00 <= vbdy1 <= vm1 or v00 >= vbdy1 >= vm1
            if ((vm1 - vbdy1)*(vbdy1 - v00) < 0d0) okay = .false.
         end if

         if (k < nz) then ! check vp1 <= vdby2 <= v00 or vp1 >= vdby2 >= v00
            if ((v00 - vbdy2)*(vbdy2 - vp1) < 0d0) okay = .false.
         end if

         if (okay) return

         c0 = v00
         c2 = 0d0

         if (k == 1) then
            c1 = (v00 - vp1) / (0.5d0*(d00 + dp1))
         else if (k == nz) then
            c1 = (vm1 - v00) / (0.5d0*(dm1 + d00))
         else if (abs(sm1) <= abs(sp1)) then
            c1 = sm1
         else
            c1 = sp1
         end if

      end subroutine get1_lpp


      subroutine calc_Pt_18_tw(s, k, Pt, ierr)
         use auto_diff
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_18var_order1), intent(out) :: Pt
         integer, intent(out) :: ierr
         type(auto_diff_real_18var_order1) :: Eturb, rho
         real(dp) :: Pt_start
         logical :: time_center, test_partials
         include 'formats'
         ierr = 0
         if (s% Eturb_alfap == 0 .or. s% Eturb_alfa == 0) then
            Pt = 0d0
            return
         end if

         rho = wrap_d_00(s,k)
         Eturb = wrap_eturb_00(s,k)
         Pt = s% Eturb_alfap*Eturb*rho
         time_center = (s% using_Fraley_time_centering .and. &
                  s% include_P_in_Fraley_time_centering)
         if (time_center) then
            Pt_start = s% Eturb_alfap*s% Eturb_start(k)*s% rho_start(k)
            Pt = 0.5d0*(Pt + Pt_start)
         end if

         if (is_bad(Pt%val)) then
!$omp critical (hydro_eturb_crit2)
            write(*,2) 'Pt', k, Pt%val
            stop 'calc_Pt_tw'
!$omp end critical (hydro_eturb_crit2)
         end if

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         if (test_partials) then
            s% solver_test_partials_val = Pt%val
            !s% solver_test_partials_var = i_var_R
            !s% solver_test_partials_dval_dx = 0 ! d_residual_dr_00
            write(*,*) 'calc_Pt_18_tw', s% solver_test_partials_var
         end if
         
      end subroutine calc_Pt_18_tw


      ! XP_18 = P_18 + avQ_18 + Pt_18 + mlt_Pturb_18 with time weighting
      subroutine calc_XP_18_tw(s, k, skip_P, skip_mlt_Pturb, XP_18, d_XP_dxa, ierr)
         use auto_diff_support
          type (star_info), pointer :: s 
         integer, intent(in) :: k
         logical, intent(in) :: skip_P, skip_mlt_Pturb
         type(auto_diff_real_18var_order1), intent(out) :: XP_18
         real(dp), dimension(s% species), intent(out) :: d_XP_dxa
         integer, intent(out) :: ierr
         integer :: j
         real(dp) :: mlt_Pturb_start
         type(auto_diff_real_18var_order1) :: rho_m1, rho_00, &
            P_18, avQ_18, Pt_18, mlt_Pturb_18
         logical :: time_center
         
         ierr = 0
         d_XP_dxa = 0d0
         
         time_center = (s% using_Fraley_time_centering .and. &
                  s% include_P_in_Fraley_time_centering)
         
         P_18 = 0d0         
         if (.not. skip_P) then
            P_18 = wrap_p_00(s, k)
            if (time_center) P_18 = 0.5d0*(P_18 + s% P_start(k))
            do j=1,s% species
               d_XP_dxa(j) = s% P(k)*s% dlnP_dxa_for_partials(j,k)
               if (time_center) d_XP_dxa(j) = 0.5d0*d_XP_dxa(j)
            end do
         end if

         avQ_18 = 0d0
         if (s% use_avQ_art_visc) then
            call get_avQ_18(s, k, avQ_18, ierr)
            if (ierr /= 0) return
            if (time_center) avQ_18 = 0.5d0*(avQ_18 + s% avQ_start(k))
         end if
         
         Pt_18 = 0d0
         if (s% Eturb_flag) then
            call calc_Pt_18_tw(s, k, Pt_18, ierr) 
            if (ierr /= 0) return
            ! note that Pt_18 is already time weighted
         end if

         mlt_Pturb_18 = 0d0
         if ((.not. skip_mlt_Pturb) .and. &
             s% mlt_Pturb_factor > 0d0 .and. s% mlt_vc_start(k) > 0d0 .and. k > 1) then
            rho_m1 = wrap_d_m1(s,k)
            rho_00 = wrap_d_00(s,k)
            mlt_Pturb_18 = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*(rho_m1 + rho_00)/6d0
            if (time_center) then
               mlt_Pturb_start = s% mlt_Pturb_factor*s% mlt_vc_start(k)**2*(s% rho_start(k-1) + s% rho_start(k))/6d0
               mlt_Pturb_18 = 0.5d0*(mlt_Pturb_18 + mlt_Pturb_start)
            end if
         end if           
         
         XP_18 = P_18 + avQ_18 + Pt_18 + mlt_Pturb_18
         
         if (s% use_other_pressure) XP_18%val = XP_18%val + s% extra_pressure(k)

      end subroutine calc_XP_18_tw
      
      
      subroutine get_avQ_18(s, k, avQ_18, ierr)
         use auto_diff
         use auto_diff_support
         type (star_info), pointer :: s      
         integer, intent(in) :: k 
         type(auto_diff_real_18var_order1), intent(out) :: avQ_18
         integer, intent(out) :: ierr
         real(dp) :: avQ, d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1
         call get_avQ(s, k, avQ, &
            d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1, ierr)
         if (ierr /= 0) return
         avQ_18%val = avQ
         avQ_18%d1Array(i_lnd_00) = d_avQ_dlnd
         avQ_18%d1Array(i_lnT_00) = d_avQ_dlnT
         avQ_18%d1Array(i_v_00) = d_avQ_dv00
         avQ_18%d1Array(i_v_p1) = d_avQ_dvp1
      end subroutine get_avQ_18
      
      
      subroutine get_avQ(s, k, avQ, & ! artificial pressure for cell k
            d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1, ierr)
         ! avQ(k) = CQ*rho(k)*max(0d0,dv)**2
         ! dv = v(k+1) - v(k) - ZSH*sqrt(P(k)/rho(k))
         ! CQ = avQ_cq, ZSH = avQ_zsh
         type (star_info), pointer :: s      
         integer, intent(in) :: k 
         real(dp), intent(out) :: avQ
         real(dp), intent(out) :: &
            d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1
         integer, intent(out) :: ierr
         call get1_avQ(s, k, avQ, &
            d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1, ierr)
         if (ierr /= 0) return
         if (s% avQ_start(k) < -1d90) s% avQ_start(k) = avQ
      end subroutine get_avQ
      
      
      subroutine get1_avQ(s, k, avQ, &
            d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1, ierr)
         ! avQ(k) = CQ*rho(k)*max(0d0,dv)**2
         ! dv = v(k+1) - v(k) - ZSH*sqrt(P(k)/rho(k))
         ! CQ = avQ_cq, ZSH = avQ_zsh
         type (star_info), pointer :: s      
         integer, intent(in) :: k 
         real(dp), intent(out) :: avQ
         real(dp), intent(out) :: &
            d_avQ_dlnd, d_avQ_dlnT, d_avQ_dv00, d_avQ_dvp1
         integer, intent(out) :: ierr
         real(dp) :: d_P_div_rho_dlnd, &
            d_P_div_rho_dlnT, d_sqrt_P_div_rho_dlnd, &
            d_sqrt_P_div_rho_dlnT, d_dv2_dlnd, d_dv2_dlnT, &
            d_dv2_dv00, d_dv2_dvp1
         real(qp) :: cq, zsh, v00, vp1, P, rho, P_div_rho, &
            sqrt_P_div_rho, dv, dv2
         logical :: test_partials
         include 'formats'
         ierr = 0
         
         avQ = 0d0
         d_avQ_dlnd = 0d0
         d_avQ_dlnT = 0d0
         d_avQ_dv00 = 0d0
         d_avQ_dvp1 = 0d0
         
         if (.not. s% v_flag) return
         if (.not. s% use_avQ_art_visc) return

         !test_partials = (k+1 == s% solver_test_partials_k)
         test_partials = .false.
         
         if (test_partials) then
            s% solver_test_partials_val = avQ
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = d_avQ_dlnT
         end if
         
         cq = s% avQ_cq
         if (cq == 0d0) return
         
         zsh = s% avQ_zsh
         
         v00 = s% v(k)
         if (k < s% nz) then
            vp1 = s% v(k+1)
         else
            vp1 = s% v_center
         end if
            
         P = s% P(k)    ! Note: not time-centered.  time-centering of avQ done at higher level.
         rho = s% rho(k)
         P_div_rho = P/rho
         d_P_div_rho_dlnd = (s% chiRho_for_partials(k) - 1d0)*P_div_rho
         d_P_div_rho_dlnT = s% chiT_for_partials(k)*P_div_rho
         
         sqrt_P_div_rho = sqrt(P_div_rho)
         d_sqrt_P_div_rho_dlnd = 0.5d0*d_P_div_rho_dlnd/sqrt_P_div_rho
         d_sqrt_P_div_rho_dlnT = 0.5d0*d_P_div_rho_dlnT/sqrt_P_div_rho
         
         dv = (vp1 - v00) - zsh*sqrt_P_div_rho
         if (dv <= 0d0) return
         
         dv2 = dv**2
         d_dv2_dlnd = -2d0*dv*zsh*d_sqrt_P_div_rho_dlnd
         d_dv2_dlnT = -2d0*dv*zsh*d_sqrt_P_div_rho_dlnT
         d_dv2_dv00 = -2d0*dv
         d_dv2_dvp1 = 2d0*dv
         
         avQ = cq*rho*dv2
         
         d_avQ_dlnd = cq*rho*d_dv2_dlnd + avQ
         d_avQ_dlnT = cq*rho*d_dv2_dlnT
         d_avQ_dv00 = cq*rho*d_dv2_dv00
         d_avQ_dvp1 = cq*rho*d_dv2_dvp1
      
         if (test_partials) then
            s% solver_test_partials_val = avQ
            s% solver_test_partials_var = s% i_v
            s% solver_test_partials_dval_dx = d_avQ_dvp1
            write(*,*) 'get1_avQ', s% solver_test_partials_var
         end if
      
      end subroutine get1_avQ
      
      
      ! marsaglia and zaman random number generator. period is 2**43 with
      ! 900 million different sequences. the state of the generator (for restarts)
      subroutine init_random(s)
         type (star_info), pointer :: s      
         integer :: ijkl,ij,kl,i,j,k,l,ii,jj,m
         real(dp) :: x,t
         ijkl = 54217137
         ij   = ijkl/30082
         kl   = ijkl - 30082 * ij
         i    = mod(ij/177,177) + 2
         j    = mod(ij,177) + 2
         k    = mod(kl/169,178) + 1
         l    = mod(kl,169)
         do ii=1,97
            x = 0.0d0
            t = 0.5d0
            do jj =1,24
               m = mod(mod(i*j,179)*k,179)
               i = j
               j = k
               k = m
               l = mod(53*l + 1, 169)
               if (mod(l*m,64) .ge. 32) x = x + t
               t = 0.5d0 * t
            enddo
            s% rand_u(ii) = x
         enddo
         s% rand_c   = 362436.0d0/16777216.0d0
         s% rand_cd  = 7654321.0d0/16777216.0d0
         s% rand_cm  = 16777213.0d0/16777216.0d0
         s% rand_i97 = 97
         s% rand_j97 = 33
      end subroutine init_random
      
      
      real(dp) function rand(s)
         type (star_info), pointer :: s      
         real(dp) :: uni
         uni = s% rand_u(s% rand_i97) - s% rand_u(s% rand_j97)
         if (uni .lt. 0.0d0) uni = uni + 1.0d0
         s% rand_u(s% rand_i97) = uni
         s% rand_i97 = s% rand_i97 - 1
         if (s% rand_i97 .eq. 0) s% rand_i97 = 97
         s% rand_j97 = s% rand_j97 - 1
         if (s% rand_j97 .eq. 0) s% rand_j97 = 97
         s% rand_c = s% rand_c - s% rand_cd
         if (s% rand_c .lt. 0.0d0) s% rand_c = s% rand_c + s% rand_cm
         uni = uni - s% rand_c
         if (uni .lt. 0.0d0) uni = uni + 1.0d0
         rand = uni
      end function rand
      
      
      real(dp) function eval1_eps_num_visc(s, k, dt) result(eps) ! erg/g/s
         ! eps_num_visc = - dKE/dt - dPE/dt - d(AvP)/dm - P*dlnd_dt/rho
         type (star_info), pointer :: s      
         integer, intent(in) :: k 
         real(dp), intent(in) :: dt
         integer :: nz
         real(dp) :: rC_start, rC_new, mC, u_new, u_old, &
            source, u_face, P_face, alfa, beta, dm
         real(qp) :: q1, q2, dKE_dt, dPE_dt, AvP00, AvPp1, d_AvP_dm
         
         include 'formats'
         
         nz = s% nz
         dm = s% dm(k)

         rC_start = s% rmid_start(k)
         rC_new = s% rmid(k)
         mC = s% m(k) - 0.5d0*dm
         q1 = 1/rC_new
         q2 = 1/rC_start
         q1 = q1 - q2
         dPE_dt = -s% cgrav(k)*mC*q1/dt ! erg/g/s
         
         if (s% u_flag) then
            u_new = s% u(k)
            u_old = s% u_start(k)
         else if (s% v_flag) then
            if (k < nz) then
               u_new = 0.5d0*(s% v(k) + s% v(k+1))
               u_old = 0.5d0*(s% v_start(k) + s% v_start(k+1))
            else
               u_new = 0.5d0*(s% v(k) + s% v_center)
               u_old = 0.5d0*(s% v_start(k) + s% v_center)
            end if
         else ! ignore kinetic energy if no velocity variables
            u_new = 0d0
            u_old = 0d0
         end if
         q1 = u_new
         q2 = u_old
         q1 = q1*q1 - q2*q2
         dKE_dt = 0.5d0*q1/dt ! erg/g/s
         
         q1 = s% lnd(k)
         q2 = s% lnd_start(k)
         q1 = q1 - q2
         source = s% P(k)*q1/(dt*s% rho(k))
         
         if (s% u_flag) then
            u_face = s% u_face(k)
            P_face = s% P_face(k)
         else if (s% v_flag) then
            u_face = s% v(k)
            if (k > 1) then
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
               beta = 1d0 - alfa
               P_face = alfa*s% P(k) + beta*s% P(k-1)
            else
               P_face = s% P(k)
            end if
         else
            u_face = 0d0
            P_face = 0d0
         end if
         AvP00 = 4d0*pi*s% r(k)*s% r(k)*u_face*P_face

         if (k == nz) then
            AvPp1 = 0d0
         else
            if (s% u_flag) then
               u_face = s% u_face(k+1)
               P_face = s% P_face(k+1)
            else if (s% v_flag) then
               u_face = s% v(k+1)
               alfa = s% dq(k)/(s% dq(k) + s% dq(k+1))
               beta = 1d0 - alfa
               P_face = alfa*s% P(k+1) + beta*s% P(k)
            else
               u_face = 0d0
               P_face = 0d0
            end if
            AvPp1 = 4d0*pi*s% r(k+1)*s% r(k+1)*u_face*P_face
         end if
         d_AvP_dm = (AvP00 - AvPp1)/dm

         q1 = dKE_dt + dPE_dt + source + d_AvP_dm
         eps = -q1
         
      end function eval1_eps_num_visc


      subroutine write_to_extra_terminal_output_file(s, str, advance)
         type (star_info), pointer :: s
         character(len=*), intent(in) :: str
         logical, intent(in) :: advance
         integer :: id, ierr, io
         if (len_trim(s% extra_terminal_output_file) > 0) then
            ierr = 0
            open(newunit=io, file=trim(s% extra_terminal_output_file), &
               action='write', position='append', iostat=ierr)
            if (ierr == 0) then
               if (advance) then
                  write(io,'(a)') trim(str)
               else
                  write(io,'(a)',advance='no') trim(str)
               end if
               close(io)
            else
               write(*,*) 'failed to open extra_terminal_output_file ' // &
                  trim(s% extra_terminal_output_file)
            end if
         end if
      end subroutine write_to_extra_terminal_output_file
      
      
      subroutine write_eos_call_info(s,k)
         use chem_def
         type (star_info), pointer :: s
         integer, intent(in) :: k ! 0 means not being called for a particular cell
         integer :: j
         include 'formats'
         !$OMP critical (omp_write_eos_call_info)
         write(*,*)
         do j=1,s% species
            write(*,4) 'xa(j,k) ' // trim(chem_isos% name(s% chem_id(j))), j, j+s% nvar_hydro, k, s% xa(j,k)
         end do
         write(*,1) 'sum(xa) =', sum(s% xa(:,k))
         write(*,1) 'Pgas = ', s% Pgas(k)
         write(*,1) 'logPgas = ', s% lnPgas(k)/ln10
         write(*,1) 'rho = ', s% rho(k)
         write(*,1) 'T = ', s% T(k)
         write(*,1) 'logQ = ', s% lnd(k)/ln10 - 2*s% lnT(k)/ln10 + 12
         write(*,*)
         write(*,1) 'eos_frac_OPAL_SCVH',    s% eos_frac_OPAL_SCVH(k)
         write(*,1) 'eos_frac_HELM',    s% eos_frac_HELM(k)
         write(*,1) 'eos_frac_Skye',    s% eos_frac_Skye(k)
         write(*,1) 'eos_frac_PC',      s% eos_frac_PC(k)
         write(*,1) 'eos_frac_FreeEOS', s% eos_frac_FreeEOS(k)
         write(*,1) 'eos_frac_CMS',     s% eos_frac_CMS(k)
         write(*,*)
         write(*,1) 'P = ', s% P(k)
         write(*,1) 'Prad = ', s% Prad(k)
         write(*,1) 'logP = ', s% lnP(k)/ln10
         write(*,1) 'logS = ', s% lnS(k)/ln10
         write(*,1) 'logE = ', s% lnE(k)/ln10
         write(*,1) 'energy = ', s% energy(k)
         write(*,1) 'grada = ', s% grada(k)
         write(*,1) 'Cv = ', s% Cv(k)
         write(*,1) 'dE_dRho = ', s% dE_dRho(k)
         write(*,1) 'cp = ', s% cp(k)
         write(*,1) 'gamma1 = ', s% gamma1(k)
         write(*,1) 'gamma3 = ', s% gamma3(k)
         write(*,1) 'eta = ', s% eta(k)
         write(*,1) 'gam = ', s% gam(k)
         write(*,1) 'mu = ', s% mu(k)
         write(*,1) 'log_free_e = ', s% lnfree_e(k)/ln10
         write(*,1) 'chiRho = ', s% chiRho(k)
         write(*,1) 'chiT = ', s% chiT(k)
         write(*,*)
         write(*,*) 'do_eos_for_cell k, nz', k, s% nz
         write(*,1) 'logRho = ', s% lnd(k)/ln10
         write(*,1) 'logT = ', s% lnT(k)/ln10
         write(*,1) 'z = ', s% Z(k)
         write(*,1) 'x = ', s% X(k)
         write(*,1) 'abar = ', s% abar(k)
         write(*,1) 'zbar = ', s% zbar(k)
         write(*,*)
         write(*,1) 'tau = ', s% tau(k)
         write(*,*)
         write(*,*) 's% eos_rq% tiny_fuzz', s% eos_rq% tiny_fuzz
         write(*,*)
         !stop 'write_eos_call_info'
         !$OMP end critical (omp_write_eos_call_info)
      end subroutine write_eos_call_info


      real(dp) function surface_avg_x(s,j)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp) :: sum_x, sum_dq
         integer :: k
         include 'formats'
         sum_x = 0
         sum_dq = 0
         do k = 1, s% nz
            sum_x = sum_x + s% xa(j,k)*s% dq(k)
            sum_dq = sum_dq + s% dq(k)
            if (sum_dq >= s% surface_avg_abundance_dq) exit
         end do
         surface_avg_x = sum_x/sum_dq
      end function surface_avg_x


      real(dp) function center_avg_x(s,j)
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = s% xa(j,k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x+ dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_avg_x = sum_x/sum_dq
      end function center_avg_x


      end module star_utils
