! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team & Matthias Fabry
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
      integer :: i, imax
      integer, parameter :: num_points = 50
      ! cannot be dp's! PGPLOT doesn't know what to do with them...
      real, dimension(2 * num_points + 1) :: &
         thetas, r1s, r2s, x1s, x2s, y1s, y2s, rs, &
         x1s_star, x2s_star, y1s_star, y2s_star, &
         x1s_RL, x2s_RL, y1s_RL, y2s_RL
      real :: a1, a2, e, x1max, x2max, xmax
      integer, pointer :: ipar(:) ! (lipar)
      real(dp), pointer :: rpar(:) ! (lrpar)
      real(dp), dimension(num_points) :: cosp
      real(dp) :: q, q_in, this_psi, xl1, psi_RL
      
      include 'formats'
      
      ierr = 0

      a1 = 1 / (1 + b% m(1) / b% m(2))
      a2 = 1 / (1 + b% m(2) / b% m(1))
      e = b% eccentricity
      
      !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static, 4)
      do i = 1, num_points
         thetas(i) = (i - 0.5) * pi / num_points
         cosp(i) = cos(thetas(i))
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
      
      q = b% m(2) / b% m(1)
      q_in = 1d0 / q  ! flip q for other star
      allocate(rpar(3))

      if (b% pg% Orbit_show_RL .and. abs(log10(q)) <= 2) then
         psi_RL = Psi_fit(b% rl(1) / b% separation, q)
         xl1 = xl1_fit(q)
         rpar(1) = psi_RL
         rpar(3) = q
         do i = 1, num_points
            rpar(2) = cosp(i)
            rs(i) = safe_root_with_guess(f_roche, 1d-1, 1d-3, &  ! function, guess, dx for bracket
               1d-3, 1d0 - 1d-3, & ! left, right bracket
               roche(1d-3, cosp(i), q) - psi_RL, roche(1d0 - 1d-3, cosp(i), q) - psi_RL, & ! f(left, right bracket)
               25, 50, 1d-4, 1d-6, & ! i_next, imax, x_tol, y_tol
               0, rpar, 0, ipar, & ! func_params
               ierr)
            x1s_RL(i) = rs(i) * cosp(i)
            y1s_RL(i) = rs(i) * sin(thetas(i))
            x1s_RL(i) = x1s_RL(i) - a1 * (1 - e)  ! displace xs
            x1s_RL(2 * num_points - i + 1) = x1s_RL(i)  ! fill out lower lobe
            y1s_RL(2 * num_points - i + 1) = -y1s_RL(i) ! flip y for lower lobe
         end do
         x1s_RL(2 * num_points + 1) = x1s_RL(1)  ! close contour
         y1s_RL(2 * num_points + 1) = y1s_RL(1)

         psi_RL = Psi_fit(b% rl(2) / b% separation, q_in)
         xl1 = xl1_fit(q_in)
         rpar(1) = psi_RL
         rpar(3) = q_in
         do i = 1, num_points
            rpar(2) = cosp(i)
            rs(i) = safe_root_with_guess(f_roche, 1d-1, 1d-3, &  ! function, guess, dx for bracket
               1d-3, 1d0 - 1d-3, & ! left, right bracket
               roche(1d-3, cosp(i), q) - psi_RL, roche(1d0 - 1d-3, cosp(i), q) - psi_RL, & ! f(left, right bracket)
               25, 50, 1d-4, 1d-6, & ! i_next, imax, x_tol, y_tol
               0, rpar, 0, ipar, & ! func_params
               ierr)
            x2s_RL(i) = rs(i) * cosp(i)
            y2s_RL(i) = rs(i) * sin(thetas(i))
            x2s_RL(i) = -(x2s_RL(i) - a2 * (1 - e))  ! displace xs
            x2s_RL(2 * num_points - i + 1) = x2s_RL(i)  ! fill out lower lobe
            y2s_RL(2 * num_points - i + 1) = -y2s_RL(i) ! flip y for lower lobe
         end do

         x2s_RL(2 * num_points + 1) = x2s_RL(1)  ! close contour
         y2s_RL(2 * num_points + 1) = y2s_RL(1)
      else if (b% pg% Orbit_show_RL .and. abs(log10(q)) > 2) then
         write(*, 1) "pgbinary: Not plotting RL, q too extreme: abs(log(q)) = ", abs(log10(q))
      end if

      if (b% pg% Orbit_show_stars .and. abs(log10(q)) <= 2) then
         if (b% point_mass_i /= 1) then
            this_psi = Psi_fit(b% r(1) / b% separation, q)
            rpar(1) = this_psi
            rpar(3) = q
            do i=1, num_points
               rpar(2) = cosp(i)
               rs(i) = safe_root_with_guess(f_roche, 1d-1, 1d-3, &  ! function, guess, dx for bracket
                  1d-3, 1d0 - 1d-3, & ! left, right bracket
                  roche(1d-3, cosp(i), q) - this_psi, roche(1d0 - 1d-3, cosp(i), q) - this_psi, & ! f(left, right bracket)
                  25, 50, 1d-4, 1d-6, & ! i_next, imax, x_tol, y_tol
                  0, rpar, 0, ipar, & ! func_params
                  ierr)
               x1s_star(i) = rs(i) * cosp(i)
               y1s_star(i) = rs(i) * sin(thetas(i))
            end do

            imax = 1  ! find highest i at which x was badly rooted
            do i = 1, num_points
               if (x1s_star(i) >= xl1) then
                  imax = i
               end if
            end do

            do i = 1, num_points
               if (i <= imax .and. x1s_star(i) > xl1) then
                  x1s_star(i) = xl1  ! fix bad xs
                  y1s_star(i) = y1s_star(imax + 1)  ! give them a good y as well
               end if
               x1s_star(i) = x1s_star(i) - a1 * (1 - e)  ! displace xs
               x1s_star(2 * num_points - i + 1) = x1s_star(i)  ! fill out lower lobe
               y1s_star(2 * num_points - i + 1) = -y1s_star(i) ! flip y for lower lobe
            end do
            x1s_star(2 * num_points + 1) = x1s_star(1)  ! close contour
            y1s_star(2 * num_points + 1) = y1s_star(1)
         else
            x1s_star(:) = 0d0
            y1s_star(:) = 0d0
         end if

         x1max = maxval(abs(x1s_star))
         xmax = max(x1max, xmax)
         x1max = maxval(abs(x1s_RL))
         xmax = max(x1max, xmax)

         if (b% point_mass_i /= 2) then
            this_psi = Psi_fit(b% r(2) / b% separation, q)
            rpar(1) = this_psi
            rpar(3) = q_in
            do i = 1, num_points
               rpar(2) = cosp(i)
               rs(i) = safe_root_with_guess(f_roche, 1d-1, 1d-3, &  ! function, guess, dx for bracket
                  1d-3, 1d0 - 1d-3, & ! left, right bracket
                  roche(1d-3, cosp(i), q) - this_psi, roche(1d0 - 1d-3, cosp(i), q) - this_psi, & ! f(left, right bracket)
                  25, 50, 1d-4, 1d-6, & ! i_next, imax, x_tol, y_tol
                  0, rpar, 0, ipar, & ! func_params
                  ierr)
               x2s_star(i) = rs(i) * cosp(i)
               y2s_star(i) = rs(i) * sin(thetas(i))
            end do

            imax = 1
            do i = 1, num_points  ! find highest i at whihc x was badly rooted
               if (x2s_star(i) >= xl1) then
                  imax = i
               end if
            end do

            do i = 1, num_points
               if (i <= imax .and. x2s_star(i) > xl1) then
                  x2s_star(i) = xl1  ! fix bad xs
                  y2s_star(i) = y2s_star(imax + 1)  ! give them a good y as well
               end if
               x2s_star(i) = -(x2s_star(i) - a2 * (1 - e))  ! displace xs
               x2s_star(2 * num_points - i + 1) = x2s_star(i)  ! fill out lower lobe
               y2s_star(2 * num_points - i + 1) = -y2s_star(i) ! flip y for lower lobe
            end do
            x2s_star(2 * num_points + 1) = x2s_star(1)  ! close contour
            y2s_star(2 * num_points + 1) = y2s_star(1)
         else
            x2s_star(:) = 0d0
            y2s_star(:) = 0d0
         end if
         x2max = maxval(abs(x2s_star))
         xmax = max(x2max, xmax)
         x2max = maxval(abs(x2s_RL))
         xmax = max(x2max, xmax)
      else if (b% pg% Orbit_show_stars .and. abs(log10(q)) > 2) then
         write(*, 1) "pgbinary: Not plotting stars, q too extreme: abs(log(q)) = ", abs(log10(q))
      end if
      deallocate(rpar)

      call pgsave  ! plot titles, axes and stuff
      call pgsch(txt_scale)
      call pgsvp(winxmin, winxmax, winymin, winymax)
      if (.not. subplot) then
         call show_model_number_pgbinary(b)
         call show_age_pgbinary(b)
      end if
      call show_title_pgbinary(b, title)

      call pgwnad(-1.1 * xmax, 1.1 * xmax, -1.1 * xmax, 1.1 * xmax)
      call pgscf(1)
      call pgsci(1)
      call show_box_pgbinary(b, 'BCSTN', 'BCSTNMV')
      call show_xaxis_label_pgbinary(b, 'separation')
      call show_left_yaxis_label_pgbinary(b, 'separation')

      call pgsci(b% pg% star_1_color)
      call pgmtxt('T', -2.0, 0.05, 0.0, 'Star 1')
      call pgslw(b% pg% pgbinary_lw / 2)
      call pgline(2 * num_points + 1, x1s, y1s)

      call pgsci(b% pg% star_2_color)
      call pgmtxt('T', -2.0 - 1.3, 0.05, 0.0, 'Star 2')
      call pgslw(b% pg% pgbinary_lw / 2)
      call pgline(2 * num_points + 1, x2s, y2s)
      call pgunsa

      if (b% pg% Orbit_show_stars .and. abs(log10(q)) <= 2) then
         call pgsave  ! plot stars
         call pgslw(int(2.0 * b% pg% pgbinary_lw / 3.0))
         call pgsfs(3)  ! select hatched
         call pgshs(45.0, 0.33, 0.0)  ! define hatching
         call pgsci(clr_Goldenrod)
         if (b% point_mass_i /= 1) then
            call pgline(2 * num_points + 1, x1s_star, y1s_star)
            call pgpoly(2 * num_points + 1, x1s_star, y1s_star)
         else
            call pgslw(int(2.0 * b% pg% pgbinary_lw))
            call pgpt1(- a1 * (1 - e), 0.0, 1)
         end if
         call pgsci(clr_LightSkyBlue)
         if (b% point_mass_i /= 2) then
            call pgline(2 * num_points + 1, x2s_star, y2s_star)
            call pgpoly(2 * num_points + 1, x2s_star, y2s_star)
         else
            call pgslw(int(2.0 * b% pg% pgbinary_lw))
            call pgpt1(a2 * (1 - e), 0.0, 1)
         end if
         call pgunsa
      end if
      if (b% pg% Orbit_show_RL .and. abs(log10(q)) <= 2) then
         call pgsave
         call pgslw(int(2.0 * b% pg% pgbinary_lw / 3.0))
         call pgline(2 * num_points + 1, x1s_RL, y1s_RL)
         call pgline(2 * num_points + 1, x2s_RL, y2s_RL)
         call pgunsa
      end if

      call pgpt1(0.0, 0.0, 5)  ! draw center of mass

   
   contains

      real function xl1_fit(qq)
         real(dp), intent(in) :: qq
         real(dp) :: logq
         
         logq = log10(qq)
         if (q > 1) logq = -logq
         xl1_fit = - 1.72452947 / pi * atan(logq * 0.21625699) + 0.5 &
            + 0.01559149 * logq &
            - 1.3924d-05 * logq * (logq + 1.5) * (logq + 4.0)
         if (q > 1) xl1_fit = 1 - xl1_fit
      
      end function xl1_fit
      
      real(dp) function Psi_fit(r_eq, qq)
         ! aprroximation of Roche potential versus q = m_other / m_this and r_eq, the &
         ! dimensionless volume equivalent radius (== r / separation of the model), with limit at 5%
         real(dp), intent(in) :: r_eq, qq
         real(dp) :: r_here
         if (req <= 0.05) then
            r_here = 0.05
         else
            r_here = r_eq
         end if
         Psi_fit = -1d0 / r_here - qq &
            - 0.5533 * (1 + qq) * r_here ** 2 &
            - 0.3642 * r_here ** 2 * (r_here ** 2 - 1) &
            - 1.8693 * r_here * (r_here - 0.1) * (r_here - 0.3) * (r_here - 0.7) * (r_here - 1.0414)
      end function Psi_fit
      
      real(dp) function roche(r, ccosp, qq)
         real(dp), intent(in) :: r, ccosp, qq
         
         roche = -1d0 / r &
            - qq * (pow(1 - 2 * r * ccosp + r**2, -0.5d0) - r * ccosp) &
            - (1 + qq) / 2 * r**2
      end function roche
      
      real(dp) function f_roche(r, dfdx, lrpar, rpar, lipar, ipar, ierr)
         real(dp), intent(in) :: r
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(out) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         f_roche = roche(r, rpar(2), rpar(3)) - rpar(1)
         dfdx = 1d0 / r**2 &
            + rpar(3) * (pow(1 - 2 * r * rpar(2) + r**2, -1.5d0) * (r - rpar(2)) + rpar(2)) &
            - (1 + rpar(3)) * r
         ierr = 0
      end function f_roche
   
   end subroutine orbit_panel


end module pgbinary_orbit

