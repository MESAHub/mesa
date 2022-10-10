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

      type (binary_info), pointer :: b
      integer, intent(in) :: device_id
      real, intent(in) :: winxmin, winxmax, winymin, winymax, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      integer :: i
      integer, parameter :: num_points = 50
      ! cannot be dp's! PGPLOT doesn't know what to do with them...
      real, dimension(num_points + 1) :: &
         thetas, r1s, r2s, x1s, x2s, y1s, y2s
      real :: a1, a2, e, a, x1max, x2max, xmax

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

      do i = 1, num_points
         thetas(i) = i * (2 * pi / num_points)
         r1s(i) = a1 * (1 - e**2) / (1 + e * cos(thetas(i)))
         r2s(i) = a2 / a1 * r1s(i)
         x1s(i) = r1s(i) * cos(thetas(i))
         y1s(i) = r1s(i) * sin(thetas(i))
         x2s(i) = -r2s(i) * cos(thetas(i))  ! minus to flip orbit
         y2s(i) = -r2s(i) * sin(thetas(i))
      end do

      x1s(num_points + 1) = x1s(1)
      x2s(num_points + 1) = x2s(1)
      y1s(num_points + 1) = y1s(1)
      y2s(num_points + 1) = y2s(1)

      x1max = maxval(abs(x1s))
      x2max = maxval(abs(x2s))
      xmax = max(x1max, x2max)

      call pgsave
      call pgsci(1)
      call pgscf(1)
      call pgsch(txt_scale)
      call pgslw(b% pg% pgbinary_lw)
      call pgwnad(-1.1 * xmax, 1.1 * xmax, -1.1 * xmax, 1.1 * xmax)
      call show_box_pgbinary(b, 'BCSTN', 'BCSTNMV')
      call show_xaxis_label_pgbinary(b, 'separation')
      call show_left_yaxis_label_pgbinary(b, 'separation')

      call pgsci(clr_Goldenrod)
      call pgline(num_points + 1, x1s, y1s)
      call pgslw(1)
      call pgmtxt('T', -2.0, 0.05, 0.0, 'Star 1')
      call pgsci(clr_LightSkyBlue)
      call pgslw(b% pg% pgbinary_lw)
      call pgline(num_points + 1, x2s, y2s)
      call pgslw(1)
      call pgmtxt('T', -2.0 - 1.3, 0.05, 0.0, 'Star 2')

      call pgsci(1)
      call pgpt(1, 0.0, 0.0, 5)

      call pgunsa

   end subroutine orbit_panel


end module pgbinary_orbit

