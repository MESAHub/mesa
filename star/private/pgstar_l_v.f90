! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
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

      module pgstar_l_v

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine L_v_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_L_v_Plot(s, id, device_id, &
            s% L_v_xleft, s% L_v_xright, &
            s% L_v_ybot, s% L_v_ytop, .false., &
            s% L_v_title, s% L_v_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine L_v_Plot


      subroutine do_L_v_Plot(s, id, device_id, &
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
            reverse_xaxis = .false., reverse_yaxis = .false.
         ierr = 0
         call do_Hist_Track(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, title, txt_scale, &
            'v_surf_km_s', 'luminosity', &
            'v km/s', 'L/L\d\(2281)', &
            s% L_v_v_min, s% L_v_v_max, &
            s% L_v_v_margin, s% L_v_dv_min, &
            s% L_v_L_min, s% L_v_L_max, &
            s% L_v_L_margin, s% L_v_dL_min, &
            s% L_v_step_min, s% L_v_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            s% show_L_v_target_box, s% L_v_target_n_sigma, &
            s% L_v_target_v, s% L_v_target_L, &
            s% L_v_target_v_sigma, s% L_v_target_L_sigma, &
            s% show_L_v_annotation1, &
            s% show_L_v_annotation2, &
            s% show_L_v_annotation3, &
            s% L_v_fname, &
            s% L_v_use_decorator, &
            s% L_v_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_L_v_Plot


      end module pgstar_l_v

