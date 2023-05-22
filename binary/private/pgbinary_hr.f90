! ***********************************************************************
!
!   Copyright (C) 2023  The MESA Team
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

module pgbinary_HR

   use const_def
   use pgbinary_support
   use binary_pgbinary

   implicit none


contains


   subroutine HR_Plot(binary_id, device_id, ierr)
      integer, intent(in) :: binary_id, device_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call get_binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return

      call pgslct(device_id)
      call pgbbuf()
      call pgeras()

      call do_HR_Plot(binary_id, device_id, &
         b% pg% HR_xleft, b% pg% HR_xright, &
         b% pg% HR_ybot, b% pg% HR_ytop, .false., &
         b% pg% HR_title, b% pg% HR_txt_scale, ierr)
      if (ierr /= 0) return

      call pgebuf()

   end subroutine HR_Plot

   subroutine do_HR_Plot(binary_id, device_id, &
      xleft, xright, ybot, ytop, subplot, &
      title, txt_scale, ierr)
      use pgbinary_star_hist_track, only: do_Star_hist_track
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, device_id
      real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      logical, parameter :: &
         reverse_xaxis = .true., reverse_yaxis = .false.
      ierr = 0

      call get_binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) return
      call do_Star_hist_track(b, id, device_id, &
         xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         'log_Teff', 'log_L', &
         'log Teff', 'log L/L\d\(2281)', &
         b% pg% HR_logT_min, b% pg% HR_logT_max, &
         b% pg% HR_logT_margin, b% pg% HR_dlogT_min, &
         b% pg% HR_logL_min, b% pg% HR_logL_max, &
         b% pg% HR_logL_margin, b% pg% HR_dlogL_min, &
         b% pg% HR_step_min, b% pg% HR_step_max, &
         reverse_xaxis, reverse_yaxis, .false., .false., &
         b% pg% show_HR_annotation1, &
         b% pg% show_HR_annotation2, &
         b% pg% show_HR_annotation3, &
         b% pg% HR_fname, &
         b% pg% HR_use_decorator, &
         b% pg% HR_pgbinary_decorator, &
         HR_decorate, ierr)


   end subroutine do_HR_Plot


end module pgbinary_HR

