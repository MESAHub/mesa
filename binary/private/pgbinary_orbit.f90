! ***********************************************************************
!
!   Copyright (C) 2022  Matthias Fabry
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

module pgbinary_orbit

   use binary_private_def
   use pgbinary_support

   implicit none


contains


   subroutine Orbit_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(id, b, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_Orbit_plot(b, id, device_id, &
         b% pg% Orbit_xleft, b% pg% Orbit_xright, &
         b% pg% Orbit_ybot, b% pg% Orbit_ytop, .false., &
         b% pg% Orbit_title, b% pg% Orbit_txt_scale, ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine Orbit_plot


   subroutine do_Orbit_plot(b, id, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: id, device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call orbit_panel(b, device_id, &
         winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
   end subroutine do_Orbit_plot


   subroutine orbit_panel(b, device_id, &
      winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)

      use num_lib, only : safe_root_with_guess
      use math_lib, only : pow
      
      type (binary_info), pointer :: b
      integer, intent(in) :: device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      integer :: i, icut
      integer, parameter :: num_points = 50
      ! cannot be dp's! PGPLOT doesn't know what to do with them...
      real, dimension(2 * num_points + 1) :: &
         thetas, r1s, r2s, x1s, x2s, y1s, y2s, rs, phis, &
         x1s_RL, x2s_RL, y1s_RL, y2s_RL
      real :: a1, a2, e, x1max, x2max, xmax
      integer, pointer :: ipar(:) ! (lipar)
      real(dp), pointer :: rpar(:) ! (lrpar)
      real(dp) :: cosp, q, this_psi, xl1

      include 'formats'

     ierr = 0
      call pgsave
      call pgsvp(winxmin, winxmax, winymin, winymax)
      if (.not. subplot) then
         call show_model_number_pgbinary(b)
         call show_age_pgbinary(b)
      end if
      call show_title_pgbinary(b, title)
      call pgunsa
      
      a1 = 1 / (1 + b% m(1) / b% m(2))
      a2 = 1 / (1 + b% m(2) / b% m(1))
      e = b% eccentricity
      
      !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(dynamic,2)
      do i = 1, num_points
         thetas(i) = (i - 0.5) * pi / num_points
         r1s(i) = a1 * (1 - e**2) / (1 + e * cos(thetas(i)))
         r2s(i) = a2 / a1 * r1s(i)
         
         x1s(i) = -r1s(i) * cos(thetas(i))  ! minus to flip orbit
         x1s(2 * num_points - i + 1) = x1s(i)
         y1s(i) = -r1s(i) * sin(thetas(i))
         y1s(2 * num_points - i + 1) = -y1s(i)  ! flip y for other side of orbit
         x2s(i) = r2s(i) * cos(thetas(i))
         x2s(2 * num_points - i + 1) = x2s(i)
         y2s(i) = r2s(i) * sin(thetas(i))
         y2s(2 * num_points - i + 1) = -y2s(i)
      end do
      !$OMP END PARALLEL DO
      x1s(2 * num_points + 1) = x1s(1)
      y1s(2 * num_points + 1) = y1s(1)
      x2s(2 * num_points + 1) = x2s(1)
      y2s(2 * num_points + 1) = y2s(1)
      
      x1max = maxval(abs(x1s))
      x2max = maxval(abs(x2s))
      xmax = max(x1max, x2max)
      
      if (b% Orbit_show_RL .and. abs(log10(q)) <= 2) then
         call pgsci(clr_Goldenrod)
         call pgpt1(x1s(1), y1s(1), -1)
         q = b% m(2) / b% m(1)
         if (b% point_mass_i /= 1) then
            this_psi = Psi_fit(b% r(1) / b% separation, q)
            xl1 = xl1_fit(q)
            do i = 1, num_points
               phis(i) = (i - 0.5) * pi / num_points
               cosp = cos(phis(i))
               rs(i) = safe_root_with_guess(f, 1d-1, 1d-3, &  ! function, guess, dx for bracket
                  1d-3, 1d0 - 1d-3, & ! left, right bracket
                  roche(1d-3, cosp), roche(1d0 - 1d-3, cosp), & ! f(left, right bracket)
                  50, 100, 1d-6, 1d-8, & ! i_next, imax, x_tol, y_tol
                  0, rpar, 0, ipar, & ! func_params
                  ierr)
               
               x1s_RL(i) = rs(i) * cosp
               y1s_RL(i) = rs(i) * sin(phis(i))
               y1s_RL(2 * num_points - i + 1) = -y1s_RL(i)
            end do
            icut = 1
            do while (x1s_RL(icut) > xl1)  ! find most extreme x within cut
               icut = icut + 1
            end do
            do i = 1, icut - 1  ! move points beyond l1 to point close to xl1
               x1s_RL(i) = x1s_RL(icut)
               y1s_RL(i) = y1s_RL(icut)
            end do
            do i = 1, num_points  ! displace the xs
               x1s_RL(i) = x1s_RL(i) - a1 * (1-e)
               x1s_RL(2 * num_points - i + 1) = x1s_RL(i)
            end do
            x1s_RL(2 * num_points + 1) = x1s_RL(1)  ! close contour
            y1s_RL(2 * num_points + 1) = y1s_RL(1)
            x1max = maxval(abs(x1s_RL))
            xmax = max(x1max, xmax)
         end if
         
         if (b% point_mass_i /= 2) then
            q = 1d0 / q  ! flip q for other star
            this_psi = Psi_fit(b% r(2) / b% separation, q)
            xl1 = xl1_fit(q)
            do i = 1, num_points
               phis(i) = (i - 0.5) * pi / num_points
               cosp = cos(phis(i))
               rs(i) = safe_root_with_guess(f, 1d-1, 1d-3, &  ! function, guess, dx for bracket
                  1d-3, 1d0 - 1d-3, & ! left, right bracket
                  roche(1d-3, cosp), roche(1d0 - 1d-3, cosp), & ! f(left, right bracket)
                  25, 50, 1d-4, 1d-6, & ! i_next, imax, x_tol, y_tol
                  0, rpar, 0, ipar, & ! func_params
                  ierr)
               
               x2s_RL(i) = rs(i) * cosp
               y2s_RL(i) = rs(i) * sin(phis(i))
               y2s_RL(2 * num_points - i + 1) = -y2s_RL(i)
            end do
            icut = 1
            do while (x2s_RL(icut) > xl1)  ! find most extreme x within cut
               icut = icut + 1
            end do
            do i = 1, icut - 1  ! move points beyond l1 to point close to xl1
               x2s_RL(i) = x2s_RL(icut)
               y2s_RL(i) = y2s_RL(icut)
            end do
            do i = 1, num_points  ! displace the xs
               x2s_RL(i) = -(x2s_RL(i) - a2 * (1-e))  ! flip x for 2nd star!
               x2s_RL(2 * num_points - i + 1) = x2s_RL(i)
            end do
            x2s_RL(2 * num_points + 1) = x2s_RL(1)
            y2s_RL(2 * num_points + 1) = y2s_RL(1)
            x2max = maxval(abs(x2s_RL))
            xmax = max(x2max, xmax)
         end if
      else if (b% Orbit_show_RL .and. abs(log10(q)) > 2) then
         write(*, 1) "pgbinary: Not plotting RL, q too extreme: abs(log(q)) = ", abs(log10(q))
      end if
      
      call pgsave
      call pgsci(1)
      call pgscf(1)
      call pgsch(txt_scale)
      call pgslw(b% pg% pgbinary_lw / 2)
      call pgwnad(-1.1 * xmax, 1.1 * xmax, -1.1 * xmax, 1.1 * xmax)
      call show_box_pgbinary(b, 'BCSTN', 'BCSTNMV')
      call show_xaxis_label_pgbinary(b, 'separation')
      call show_left_yaxis_label_pgbinary(b, 'separation')
      
      call pgsci(clr_Goldenrod)
      call pgline(2 * num_points + 1, x1s, y1s)
      call pgslw(1)
      call pgmtxt('T', -2.0, 0.05, 0.0, 'Star 1')
      call pgsci(clr_LightSkyBlue)
      call pgslw(b% pg% pgbinary_lw)
      call pgline(2 * num_points + 1, x2s, y2s)
      
      if (b% Orbit_show_RL .and. abs(log10(q)) <= 2) then
         call pgslw(int(2.0 * b% pg% pgbinary_lw / 3.0))
         call pgsfs(3)
         call pgshs(45.0, 0.33, 0.0)
         call pgsci(clr_Goldenrod)
         call pgline(2 * num_points + 1, x1s_RL, y1s_RL)
         call pgpoly(2 * num_points + 1, x1s_RL, y1s_RL)
         call pgsci(clr_LightSkyBlue)
         call pgline(2 * num_points + 1, x2s_RL, y2s_RL)
         call pgpoly(2 * num_points + 1, x2s_RL, y2s_RL)
      end if
      
      call pgslw(1)
      call pgmtxt('T', -2.0 - 1.3, 0.05, 0.0, 'Star 2')
      
      call pgsci(1)
      call pgpt1(0.0, 0.0, 5)
      call pgunsa
   
   contains
      
      real function xl1_fit(q)
         real(dp), intent(in) :: q
         real(dp) :: logq
         
         logq = log10(q)
         if (q > 1) logq = -logq
         xl1_fit = - 1.72452947 / pi * atan(logq * 0.21625699) + 0.5 &
            + 0.01559149 * logq &
            - 1.3924d-05 * logq * (logq + 1.5) * (logq + 4.0)
         if (q > 1) xl1_fit = 1 - xl1_fit
      
      end function xl1_fit
      
      real(dp) function Psi_fit(req, q)
         ! fit of Roche potential versus q = m_other / m_this and r_eq, the &
         ! dimensionless volume equivalent radius (== r / separation of the model)
         real(dp), intent(in) :: req, q
         Psi_fit = -1d0 / req - q &
            - 0.5533 * (1 + q) * req ** 2 &
            - 0.3642 * req ** 2 * (req ** 2 - 1) &
            - 1.8693 * req * (req - 0.1) * (req - 0.3) * (req - 0.7) * (req - 1.0414)
      end function Psi_fit
      
      real(dp) function roche(r, cosp)
         real(dp), intent(in) :: r, cosp
         
         roche = -1d0 / r &
            - q * (pow(1 - 2 * r * cosp + r**2, -0.5d0) - r * cosp) &
            - (1 + q) / 2 * r**2
      end function roche
      
      real(dp) function f(r, dfdx, lrpar, rpar, lipar, ipar, ierr)
         real(dp), intent(in) :: r
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(out) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         f = roche(r, cosp) - this_psi
         dfdx = 1d0 / r**2 &
            + q * (pow(1 - 2 * r * cosp + r**2, -1.5d0) * (r - cosp) + cosp) &
            - (1 + q) * r
         ierr = 0
      end function f
   
   end subroutine orbit_panel


end module pgbinary_orbit

