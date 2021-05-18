! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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

      module pgstar_logg_logT

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine logg_logT_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_logg_logT_Plot(s, id, device_id, &
            s% logg_logT_xleft, s% logg_logT_xright, &
            s% logg_logT_ybot, s% logg_logT_ytop, .false., &
            s% logg_logT_title, s% logg_logT_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine logg_logT_Plot


      subroutine do_logg_logT_Plot(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, &
            title, txt_scale, ierr)
         use pgstar_hist_track, only: null_decorate, do_Hist_Track
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr
         logical, parameter :: &
            reverse_xaxis = .true., reverse_yaxis = .true.
         ierr = 0
         call do_Hist_Track(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, title, txt_scale, &
            'log_Teff', 'log_g', &
            'logT', 'log g', &
            s% logg_logT_logT_min, s% logg_logT_logT_max, &
            s% logg_logT_logT_margin, s% logg_logT_dlogT_min, &
            s% logg_logT_logg_min, s% logg_logT_logg_max, &
            s% logg_logT_logg_margin, s% logg_logT_dlogg_min, &
            s% logg_logT_step_min, s% logg_logT_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            s% show_logg_logT_target_box, s% logg_logT_target_n_sigma, &
            s% logg_logT_target_logT, s% logg_logT_target_logg, &
            s% logg_logT_target_logT_sigma, s% logg_logT_target_logg_sigma, &
            s% show_logg_logT_annotation1, &
            s% show_logg_logT_annotation2, &
            s% show_logg_logT_annotation3, &
            s% logg_logT_fname, &
            s% logg_logT_use_decorator, &
            s% logg_logT_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_logg_logT_Plot


      end module pgstar_logg_logT

